#!/usr/bin/env python3
"""
analyze_baseline.py

Baseline (pedestal) quality check for the AD9694 digitizer boards, from the
ROOT files written by raw2root_new (ad9694daqtrigchain), e.g.
run00001_file0000_20261028_172125.root. Meant for runs without signals
(software trigger), where every sample is pedestal plus noise.

Usage:
    python3 analyze_baseline.py              # asks for the file and output folder
    python3 analyze_baseline.py run.root
    python3 analyze_baseline.py run.root --out /path/to/outdir --label "run 1"
    python3 analyze_baseline.py run.root --max-triggers 100

Without a ROOT file on the command line the script asks for the file and
the output folder; Enter keeps the folder shown in [brackets] and Tab
completes paths. Every other setting keeps its default unless it is given
as an option below. Needs daqyy_common.py in the same folder.

Input: one TTree per channel, named ch<id> with id = 10000 + 10*board + channel.
Each entry is one trigger: chdata[nsamps] holds the ADC samples and ttt the
trigger time tag, which all channels of a board share. raw2root_new also
stores baseline and rms (mean and sample standard deviation of the first 20
samples), area (sum of baseline - chdata over samples 35..100), xmin and xmax
(positions of the lowest and highest sample), tbusy and t2prev.

Checks:
    pedestal and noise    Gaussian fit to the histogram of all samples of a channel
    stability             mean and standard deviation of every waveform
    noise spectrum        Welch power spectrum, frequency in cycles/sample
    correlated noise      correlation between the channels of each board
    data integrity        missing triggers, short waveforms, tbusy, saturation,
                          spikes, stuck lowest bit (flags explained at the end
                          of baseline_summary.txt)
All values stay in ADC counts and samples.

Outputs (written into --out, default: <run name>_baseline in the folder of the input file):
    baseline_summary.txt         the tables printed to the terminal
    baseline_channels.csv        one row per channel: fit, noise, flags
    baseline_events.csv          one row per channel and trigger
    baseline_boards.csv          one row per board: triggers, correlated noise
    noise_correlation.csv        correlation of every channel pair of a board
    noise_psd.csv                noise power spectrum of every channel
    summary.png                  pedestal, noise and event count of every channel
    pedestal_map.png             pedestal level and sigma of every channel, colored by sigma
    pedestal_hist_board<B>.png   ADC histograms with Gaussian fits, one page per board
    pedestal_stability.png       pedestal shift and noise of every waveform
    noise_spectra.png            noise power spectra
    noise_correlation.png        noise correlation matrix of every board
"""

import os

# One thread is plenty here; keep numpy's BLAS from starting one per core.
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse
import csv
import itertools
import math
import sys
from functools import partial

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import LogFormatter, MaxNLocator, NullLocator
from scipy.ndimage import median_filter
from scipy.optimize import curve_fit
from scipy.signal import welch

from daqyy_common import (
    ADC_CODES, CHANNELS_PER_BOARD, SURFACE, INK, INK2, MUTED, AXIS, BLUE, ORANGE,
    GOOD, WARNING, CRITICAL, SEQUENTIAL, DIVERGING,
    open_root, channel_trees, read_branches, waveform, run_name, default_output,
    make_output_folder, enable_path_completion, ask, existing_file, output_folder, Progress,
)

FLAG_HELP = {
    "NO_DATA": "no samples at all",
    "DEAD": "every sample has the same ADC code",
    "FIT_FAILED": "Gaussian fit failed; raw mean and standard deviation used instead",
    "NOISY": "noise sigma above --noisy-factor times the median of all channels",
    "QUIET": "noise sigma below the median of all channels divided by --noisy-factor",
    "SPIKES": "more samples beyond --spike-sigma than Gaussian noise gives (> max(3, 10x expected))",
    "SATURATED": "samples at ADC code 0 or 16383",
    "LSB": "share of odd ADC codes outside 0.4-0.6 (lowest bit stuck)",
    "PSD_LINE": "narrow line in the noise spectrum, more than 5 sigma above its running median",
    "MISSING": "no waveform for a trigger that other channels of the board recorded",
    "SHORT_WF": "waveform shorter than the full length",
    "LONG_WF": "waveform longer than the full length",
    "FEWER_TRIGGERS": "board recorded fewer triggers than the other boards",
}


def gaussian(x, amplitude, mean, sigma):
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def fit_gaussian(counts):
    """Fit a Gaussian to an ADC-code histogram, within 5 sigma of its median.

    Bins are weighted by their Poisson error sqrt(N) and empty bins are
    skipped, like ROOT's default chi-square fit. Returns a dict, or None if
    the fit could not be performed."""
    cdf = np.cumsum(counts) / counts.sum()
    median = float(np.searchsorted(cdf, 0.5))
    width = max((np.searchsorted(cdf, 0.8413) - np.searchsorted(cdf, 0.1587)) / 2, 0.5)
    lo = max(int(median - 5 * width), 0)
    hi = min(int(median + 5 * width) + 2, counts.size)
    x = np.arange(lo, hi, dtype=float)
    y = counts[lo:hi].astype(float)
    sel = y > 0
    if np.count_nonzero(sel) < 4:
        return None  # need more filled bins than fit parameters
    x, y = x[sel], y[sel]
    try:
        popt, pcov = curve_fit(gaussian, x, y, p0=[y.max(), median, width],
                               sigma=np.sqrt(y), absolute_sigma=True, maxfev=10000)
    except (RuntimeError, ValueError):
        return None
    perr = np.sqrt(np.diag(pcov))
    return {
        "amplitude": popt[0], "mean": popt[1], "sigma": abs(popt[2]),
        "mean_err": perr[1], "sigma_err": perr[2],
        "chi2": float((((y - gaussian(x, *popt)) / np.sqrt(y)) ** 2).sum()), "ndf": len(y) - 3,
    }


def noise_psd(wave, nfft):
    """Welch power spectrum of one waveform: the average over its consecutive
    nfft-sample segments, each with its mean removed and a Hann window. With
    fs = 1 the frequency is in cycles/sample and the PSD in ADC^2 per
    (cycle/sample), so its integral over 0..0.5 is the noise variance.
    Returns (psd, number of segments); (None, 0) if the waveform is too short."""
    nseg = len(wave) // nfft
    if nseg == 0:
        return None, 0
    _, psd = welch(wave[:nseg * nfft], fs=1.0, window="hann", nperseg=nfft,
                   noverlap=0, detrend="constant", scaling="density")
    return psd, nseg


class Channel:
    """Everything collected for one channel tree."""

    def __init__(self, name, board, chn, tree):
        self.name, self.board, self.chn, self.tree = name, board, chn, tree
        index = read_branches(tree, ("ttt", "nsamps"))
        self.ttt, self.nsamps = index["ttt"], index["nsamps"]
        self.counts = np.zeros(ADC_CODES, dtype=np.int64)  # histogram of all samples
        self.n_overflow = 0                                 # samples above code 16383
        self.psd_sum, self.psd_segments = 0.0, 0
        self.events = []                                    # one dict per trigger
        self.flags = []
        self._last_tbusy = 0

    def add_event(self, trigger, entry, wave, full_length, nfft):
        """Accumulate one waveform. The tree must hold this entry (GetEntry done)."""
        t = self.tree
        counts = np.bincount(wave, minlength=ADC_CODES)
        self.n_overflow += int(counts[ADC_CODES:].sum())
        self.counts += counts[:ADC_CODES]
        psd, nseg = noise_psd(wave, nfft)
        if nseg:
            self.psd_sum = self.psd_sum + psd * nseg
            self.psd_segments += nseg
        if len(wave) == full_length:
            status = "ok"
        else:
            status = "short" if len(wave) < full_length else "long"
        tbusy = int(t.tbusy)
        empty = len(wave) == 0
        self.events.append({
            "tree": self.name, "board": self.board, "channel": self.chn,
            "trigger": trigger, "status": status, "entry": entry,
            "ttt": int(t.ttt), "t2prev": t.t2prev, "tbusy": tbusy,
            "tbusy_step": tbusy - self._last_tbusy, "nsamps": len(wave),
            "mean": math.nan if empty else float(wave.mean()),
            "std": math.nan if empty else float(wave.std()),
            "min": "" if empty else int(wave.min()), "max": "" if empty else int(wave.max()),
            "daq_baseline": t.baseline, "daq_rms": t.rms, "daq_area": t.area,
            "daq_xmin": t.xmin, "daq_xmax": t.xmax,
        })
        self._last_tbusy = tbusy

    def add_missing(self, trigger, ttt):
        self.events.append({"tree": self.name, "board": self.board, "channel": self.chn,
                            "trigger": trigger, "status": "missing", "ttt": ttt})

    def finish(self, spike_sigma, nfft):
        """Fit and summary numbers, once all triggers are in."""
        nan = math.nan
        recorded = [e for e in self.events if e["status"] != "missing"]
        self.n_triggers = len(self.events)
        self.n_missing = self.n_triggers - len(recorded)
        self.n_short = sum(e["status"] == "short" for e in recorded)
        self.n_long = sum(e["status"] == "long" for e in recorded)
        self.n_samples = n = int(self.counts.sum()) + self.n_overflow
        self.fit = None
        self.raw_mean = self.raw_std = self.mu = self.sigma = nan
        self.adc_min = self.adc_max = self.odd_fraction = self.expected_outliers = nan
        self.n_low = self.n_high = self.n_saturated = 0
        self.spread = self.event_sigma = self.daq_rms = self.daq_baseline_scatter = nan
        self.area_mean = self.area_std = self.tbusy_last = nan
        self.psd = None
        self.psd_sigma = self.psd_line_freq = self.psd_line_ratio = self.psd_line_sigma = nan
        if n == 0:
            return

        # ---- pedestal and noise from the histogram of all samples ----
        codes = np.arange(ADC_CODES)
        filled = np.flatnonzero(self.counts)
        if len(filled):
            self.adc_min, self.adc_max = int(filled[0]), int(filled[-1])
            in_range = self.counts.sum()
            self.raw_mean = float((codes * self.counts).sum() / in_range)
            self.raw_std = float(np.sqrt(((codes - self.raw_mean) ** 2 * self.counts).sum() / in_range))
            self.odd_fraction = float(self.counts[1::2].sum() / in_range)
        if self.raw_std > 0:
            self.fit = fit_gaussian(self.counts)
        if self.fit:
            self.mu, self.sigma = self.fit["mean"], self.fit["sigma"]
        else:
            self.mu, self.sigma = self.raw_mean, self.raw_std
        cut = spike_sigma * self.sigma
        self.n_low = int(self.counts[codes < self.mu - cut].sum())
        self.n_high = int(self.counts[codes > self.mu + cut].sum()) + self.n_overflow
        self.expected_outliers = n * math.erfc(spike_sigma / math.sqrt(2))
        self.n_saturated = int(self.counts[0] + self.counts[-1]) + self.n_overflow

        # ---- per-waveform numbers ----
        recorded = [e for e in recorded if e["nsamps"]]
        means = np.array([e["mean"] for e in recorded])
        areas = np.array([e["daq_area"] for e in recorded])
        if len(recorded) > 1:
            self.spread = float(np.std(means, ddof=1))
            self.area_std = float(np.std(areas, ddof=1))
        self.event_sigma = float(np.mean([e["std"] for e in recorded]))
        self.daq_rms = float(np.mean([e["daq_rms"] for e in recorded]))
        self.daq_baseline_scatter = float(np.sqrt(np.mean(
            [(e["daq_baseline"] - e["mean"]) ** 2 for e in recorded])))
        self.area_mean = float(areas.mean())
        self.tbusy_last = recorded[-1]["tbusy"]

        # ---- noise spectrum ----
        if self.psd_segments:
            self.psd = self.psd_sum / self.psd_segments
            self.psd_sigma = math.sqrt(self.psd.sum() / nfft)  # integral, bin width 1/nfft
            # strongest narrow line: each bin over the running median of the spectrum,
            # which follows broad shapes (1/f rise, white floor) but not narrow lines
            inner = self.psd[1:-1]                             # without DC and Nyquist
            ratio = inner / median_filter(inner, size=31, mode="nearest")
            k = int(np.argmax(ratio))
            self.psd_line_freq = (k + 1) / nfft
            self.psd_line_ratio = float(ratio[k])
            # a bin averaged over K segments scatters by 1/sqrt(K) around the true PSD
            self.psd_line_sigma = (ratio[k] - 1) * math.sqrt(self.psd_segments)


class Correlation:
    """Noise correlation between the channels of one board, pooled over triggers."""

    def __init__(self):
        n = CHANNELS_PER_BOARD
        self.sxy, self.sxx, self.syy = np.zeros((n, n)), np.zeros((n, n)), np.zeros((n, n))
        self.samples = np.zeros((n, n))
        self.triggers = np.zeros((n, n), dtype=int)
        self.n_used = 0

    def add(self, waves):
        """waves: {channel: waveform} of one trigger, all of the same length."""
        if len(waves) < 2:
            return
        chs = sorted(waves)
        x = np.array([waves[c] - waves[c].mean() for c in chs])
        sxy = x @ x.T
        var = np.diag(sxy)
        idx = np.ix_(chs, chs)
        self.sxy[idx] += sxy
        self.sxx[idx] += var[:, None]
        self.syy[idx] += var[None, :]
        self.samples[idx] += x.shape[1]
        self.triggers[idx] += 1
        self.n_used += 1

    def result(self):
        """Correlation coefficient and covariance per sample (ADC^2); NaN without data."""
        with np.errstate(invalid="ignore", divide="ignore"):
            return self.sxy / np.sqrt(self.sxx * self.syy), self.sxy / self.samples


def open_channels(root_file):
    """One Channel per ch<id> tree, sorted by board and channel."""
    trees = [(board, chn, name) for board, chs in channel_trees(root_file).items()
             for chn, name in chs.items()]
    bar = Progress("Scanning channels", len(trees))
    chans = []
    for board, chn, name in trees:
        chans.append(Channel(name, board, chn, root_file.Get(name)))
        bar.update(note=name)
    bar.close()
    return chans


def analyze_board(chans, full_length, nfft, max_triggers, progress):
    """Loop over the triggers of one board. All channels of a board share the
    trigger time tag, so entries with equal ttt belong to one trigger.
    progress is told about every trigger. Returns (number of triggers, Correlation)."""
    by_chn = {c.chn: c for c in chans}
    corr = Correlation()
    entries = sorted((int(t), c.chn, i) for c in chans for i, t in enumerate(c.ttt))
    n_trig = 0
    for trigger, (ttt, group) in enumerate(itertools.groupby(entries, key=lambda e: e[0])):
        if max_triggers is not None and trigger >= max_triggers:
            break
        full, present = {}, set()
        for _, chn, entry in group:
            ch = by_chn[chn]
            wave = waveform(ch.tree, entry)
            ch.add_event(trigger, entry, wave, full_length, nfft)
            present.add(chn)
            if len(wave) == full_length:
                full.setdefault(chn, wave)
        for chn in by_chn.keys() - present:
            by_chn[chn].add_missing(trigger, ttt)
        corr.add(full)
        n_trig += 1
        progress.update(note=f"board {chans[0].board}")
    return n_trig, corr


def summarize_board(info, most_triggers):
    """Row for baseline_boards.csv; also stores r and cov in info."""
    chans, corr = info["chans"], info["corr"]
    r, cov = corr.result()
    info["r"], info["cov"] = r, cov
    row = {
        "board": info["board"], "n_channels": len(chans), "n_triggers": info["n_triggers"],
        "n_missing": sum(c.n_missing for c in chans),
        "n_short": sum(c.n_short + c.n_long for c in chans),
        "corr_triggers": corr.n_used, "mean_r": math.nan, "max_r": math.nan,
        "max_r_pair": "", "sigma_coh": math.nan,
        "flags": "FEWER_TRIGGERS" if info["n_triggers"] < most_triggers else "OK",
    }
    iu = np.triu_indices(CHANNELS_PER_BOARD, 1)
    have = corr.triggers[iu] > 0
    if have.any():
        rr, cc = r[iu][have], cov[iu][have]
        k = int(np.argmax(np.abs(rr)))
        row["mean_r"] = float(rr.mean())
        row["max_r"] = float(rr[k])
        row["max_r_pair"] = f"ch{iu[0][have][k]}-ch{iu[1][have][k]}"  # "4-5" would turn into a date in Excel
        # common-mode noise: for x_i = c + u_i, every pair has covariance var(c)
        row["sigma_coh"] = math.sqrt(max(float(cc.mean()), 0.0))
    info["row"] = row
    return row


def assign_flags(chans, noisy_factor):
    """Set c.flags for every channel (see FLAG_HELP). Returns the median sigma."""
    sigmas = [c.sigma for c in chans if c.sigma > 0]
    median_sigma = float(np.median(sigmas)) if sigmas else math.nan
    for c in chans:
        flags = []
        if c.n_samples == 0:
            flags.append("NO_DATA")
        elif not c.raw_std > 0:
            flags.append("DEAD")
        else:
            if c.fit is None:
                flags.append("FIT_FAILED")
            if c.sigma > noisy_factor * median_sigma:
                flags.append("NOISY")
            elif c.sigma < median_sigma / noisy_factor:
                flags.append("QUIET")
            if c.n_low + c.n_high > max(3, 10 * c.expected_outliers):
                flags.append("SPIKES")
            if abs(c.odd_fraction - 0.5) > 0.1:
                flags.append("LSB")
            if c.psd_line_sigma > 5:
                flags.append("PSD_LINE")
        if c.n_saturated:
            flags.append("SATURATED")
        if c.n_missing:
            flags.append("MISSING")
        if c.n_short:
            flags.append("SHORT_WF")
        if c.n_long:
            flags.append("LONG_WF")
        c.flags = flags
    return median_sigma


def write_csv(path, rows, columns):
    """Write dict rows; floats with 8 significant digits, NaN as empty."""
    def fmt(v):
        if isinstance(v, float):
            return "" if math.isnan(v) else f"{v:.8g}"
        return v
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, restval="", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: fmt(v) for k, v in row.items()})


def board_bounds(chans):
    """(board, first index, end index) of every board in a sorted channel list."""
    starts = [i for i, c in enumerate(chans) if i == 0 or c.board != chans[i - 1].board]
    return [(chans[s].board, s, e) for s, e in zip(starts, starts[1:] + [len(chans)])]


def colorbar(fig, im, ax, label):
    """Colorbar just right of ax, without shrinking ax."""
    cax = ax.inset_axes([1.01, 0.0, 0.015, 1.0])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(label)
    cb.outline.set_visible(False)
    return cb


def text_color_on(rgba):
    """White or ink, whichever contrasts more with a fill color (WCAG luminance)."""
    lin = [v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in rgba[:3]]
    lum = 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]
    return "#ffffff" if 1.05 / (lum + 0.05) > (lum + 0.05) / 0.0534 else INK


def plot_summary(chans, median_sigma, noisy_factor, path, label):
    x = np.arange(len(chans))
    fig, axes = plt.subplots(3, 1, figsize=(13, 9), sharex=True, gridspec_kw={"hspace": 0.45})

    ax = axes[0]
    ax.plot(x, [c.mu for c in chans], "o", ms=4, color=BLUE)
    ax.set_ylabel("ADC")
    ax.set_title("Pedestal mean (Gaussian fit to all samples)", loc="left")

    ax = axes[1]
    ax.plot(x, [c.sigma for c in chans], "o", ms=4, color=BLUE,
            label="Gaussian fit σ, all samples")
    ax.plot(x, [c.daq_rms for c in chans], "D", ms=3.5, color=ORANGE,
            label="DAQ rms, first 20 samples (mean over waveforms)")
    for level, text in ((median_sigma * noisy_factor, "NOISY above"),
                        (median_sigma / noisy_factor, "QUIET below")):
        ax.axhline(level, color=MUTED, lw=0.8, ls="--")
        ax.text(len(chans) - 0.4, level, f" {text}", va="center", ha="left", fontsize=7, color=INK2)
    ax.set_ylabel("ADC")
    ax.set_title("Noise", loc="left")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, borderaxespad=0.2)

    ax = axes[2]
    missing = np.array([c.n_missing for c in chans])
    short = np.array([c.n_short + c.n_long for c in chans])
    full = np.array([c.n_triggers for c in chans]) - missing - short
    bar = {"width": 0.72, "edgecolor": SURFACE, "linewidth": 1.0}
    ax.bar(x, full, color=GOOD, label="full-length waveform", **bar)
    ax.bar(x, short, bottom=full, color=WARNING, label="short waveform", **bar)
    ax.bar(x, missing, bottom=full + short, color=CRITICAL, label="missing", **bar)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylabel("triggers")
    ax.set_title("Waveforms per channel", loc="left")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=3, borderaxespad=0.2)

    bounds = board_bounds(chans)
    for ax in axes:
        ax.grid(axis="x", visible=False)
        for _, start, _ in bounds[1:]:
            ax.axvline(start - 0.5, color=AXIS, lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([str(c.chn) for c in chans], fontsize=6)
    ax.set_xlim(-0.6, len(chans) - 0.4)
    for board, start, end in bounds:
        ax.text((start + end - 1) / 2, -0.16, f"board {board}", transform=ax.get_xaxis_transform(),
                ha="center", va="top", color=INK2)
    ax.set_xlabel("channel number within each board", labelpad=26)
    fig.suptitle(f"Baseline summary: {label}", x=0.01, y=0.95, ha="left", fontsize=10)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_pedestal_map(chans, path, label):
    """Board x channel heatmap colored by the fit sigma (the width of the
    pedestal); each cell gives the pedestal level and, below it, the sigma."""
    boards = sorted({c.board for c in chans})
    level = np.full((len(boards), CHANNELS_PER_BOARD), np.nan)
    sigma = np.full_like(level, np.nan)
    for c in chans:
        level[boards.index(c.board), c.chn] = c.mu
        sigma[boards.index(c.board), c.chn] = c.sigma
    height = 0.62 * len(boards) + 1.6
    fig, ax = plt.subplots(figsize=(8.5, height))
    fig.subplots_adjust(top=1 - 0.7 / height)  # margin in inches, for any number of boards
    im = ax.imshow(sigma, cmap=SEQUENTIAL, aspect="auto", interpolation="nearest")
    colorbar(fig, im, ax, "fit σ (ADC)")
    for (i, j), v in np.ndenumerate(level):
        if np.isfinite(v):
            color = text_color_on(im.cmap(im.norm(sigma[i, j])))
            ax.text(j, i - 0.13, f"{v:.1f}", ha="center", va="center", fontsize=7.5, color=color)
            ax.text(j, i + 0.2, f"σ {sigma[i, j]:.2f}", ha="center", va="center", fontsize=6.5,
                    color=color)
        else:
            ax.text(j, i, "no data", ha="center", va="center", fontsize=6, color=MUTED)
    ax.set_xticks(range(CHANNELS_PER_BOARD))
    ax.set_yticks(range(len(boards)))
    ax.set_yticklabels([f"board {b}" for b in boards])
    ax.set_xlabel("channel")
    # 2 px surface-colored gaps between the cells instead of drawn borders
    ax.set_xticks(np.arange(-0.5, CHANNELS_PER_BOARD), minor=True)
    ax.set_yticks(np.arange(-0.5, len(boards)), minor=True)
    ax.grid(False)
    ax.grid(which="minor", color=SURFACE, linewidth=2)
    ax.tick_params(which="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("Gaussian fit to all samples: pedestal level (upper number) and σ (lower number, "
                 "color), ADC", loc="left")
    fig.suptitle(f"Pedestal per channel: {label}", x=0.125, y=1 - 0.1 / height, va="top",
                 ha="left", fontsize=10)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_histograms(board, chans, spike_sigma, path, label):
    by_chn = {c.chn: c for c in chans}
    fig, axes = plt.subplots(2, 4, figsize=(13, 6.6), gridspec_kw={"hspace": 0.55, "wspace": 0.3})
    for chn, ax in enumerate(axes.ravel()):
        c = by_chn.get(chn)
        if c is None or math.isnan(c.adc_max):
            ax.set_title(f"ch {chn}: no data", loc="left", color=MUTED)
            ax.axis("off")
            continue
        half = (spike_sigma + 2) * c.sigma if c.sigma > 0 else 3  # threshold lines stay in view
        lo = max(int(math.floor(c.mu - half)), 0)
        hi = min(int(math.ceil(c.mu + half)), ADC_CODES - 1)
        counts = c.counts[lo:hi + 1]
        edges = np.arange(lo, hi + 2) - 0.5
        ax.stairs(counts, edges, fill=True, color=BLUE, alpha=0.12, lw=0)
        ax.stairs(counts, edges, color=BLUE, lw=1.2)
        if c.fit:
            xx = np.linspace(edges[0], edges[-1], 400)
            ax.plot(xx, gaussian(xx, c.fit["amplitude"], c.mu, c.sigma), color=INK, lw=1)
        if c.sigma > 0:
            for side in (-1, 1):
                ax.axvline(c.mu + side * spike_sigma * c.sigma, color=MUTED, lw=0.8, ls="--")
        ax.set_yscale("log")
        ax.yaxis.set_minor_locator(NullLocator())  # minor log ticks cost seconds per page
        ax.set_ylim(0.5, max(counts.max(), 1) * 30)
        ax.set_xlim(edges[0], edges[-1])
        fit_text = f"χ²/ndf = {c.fit['chi2']:.0f}/{c.fit['ndf']}" if c.fit else "no fit"
        ax.set_title(f"{c.name} (ch {chn})\nμ = {c.mu:.2f}   σ = {c.sigma:.3f}   {fit_text}",
                     loc="left")
        if c.flags:
            ax.text(0.02, 0.97, "⚠ " + "\n   ".join(c.flags), transform=ax.transAxes,
                    va="top", fontsize=7, color=CRITICAL)
        outside = c.n_samples - int(counts.sum())
        if outside:
            ax.text(0.98, 0.97, f"{outside} samples\noutside view", transform=ax.transAxes,
                    va="top", ha="right", fontsize=7, color=INK2)
        ax.set_xlabel("ADC code")
        if chn % 4 == 0:
            ax.set_ylabel("samples")
    handles = [Patch(facecolor=to_rgba(BLUE, 0.12), edgecolor=BLUE, label="all samples, all waveforms"),
               Line2D([], [], color=INK, lw=1, label="Gaussian fit"),
               Line2D([], [], color=MUTED, lw=0.8, ls="--", label=f"±{spike_sigma:g}σ spike threshold")]
    fig.legend(handles=handles, loc="lower right", bbox_to_anchor=(0.9, 0.93), ncol=3)
    fig.suptitle(f"Pedestal histograms, board {board}: {label}", x=0.125, y=0.975, ha="left",
                 fontsize=10)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_stability(chans, n_triggers, path, label):
    shift = np.full((len(chans), n_triggers), np.nan)
    noise = np.full_like(shift, np.nan)
    missing, odd = [], []
    for row, c in enumerate(chans):
        means = [e["mean"] for e in c.events if e["status"] != "missing" and e["nsamps"]]
        ref = np.median(means) if means else np.nan
        for e in c.events:
            if e["status"] == "missing":
                missing.append((e["trigger"], row))
                continue
            shift[row, e["trigger"]] = e["mean"] - ref
            noise[row, e["trigger"]] = e["std"]
            if e["status"] != "ok":
                odd.append((e["trigger"], row))

    fig, axes = plt.subplots(1, 2, figsize=(12, 10), sharey=True, gridspec_kw={"wspace": 0.25})
    lim = np.nanmax(np.abs(shift)) if np.isfinite(shift).any() else 1.0
    panels = ((axes[0], shift, DIVERGING, -max(lim, 1e-3), max(lim, 1e-3),
               "Pedestal shift: waveform mean − channel median"),
              (axes[1], noise, SEQUENTIAL, None, None,
               "Noise: standard deviation of each waveform"))
    for ax, data, cmap, vmin, vmax, title in panels:
        im = ax.imshow(data, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                       interpolation="nearest")
        colorbar(fig, im, ax, "ADC")
        ax.grid(False)
        ax.tick_params(length=0)
        ax.set_title(title, loc="left")
        ax.set_xlabel("trigger number (per board)")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        for _, start, _ in board_bounds(chans)[1:]:
            ax.axhline(start - 0.5, color=SURFACE, lw=2)
        if n_triggers <= 100:  # markers would cover long runs; the CSV has them all
            if missing:
                ax.scatter(*zip(*missing), marker="x", s=16, color=INK2, lw=0.9)
            if odd:
                ax.scatter(*zip(*odd), marker="o", s=12, facecolors="none", edgecolors=INK, lw=0.8)
    axes[0].set_yticks(range(len(chans)))
    axes[0].set_yticklabels([c.name for c in chans], fontsize=5.5)
    if n_triggers <= 100:
        handles = [Line2D([], [], marker="x", ls="", color=INK2, label="missing"),
                   Line2D([], [], marker="o", ls="", mfc="none", mec=INK, label="short waveform")]
        fig.legend(handles=handles, loc="lower right", bbox_to_anchor=(0.9, 0.9), ncol=2)
    fig.suptitle(f"Pedestal stability: {label}", x=0.125, y=0.93, ha="left", fontsize=10)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_spectra(chans, nfft, path, label):
    rows = [c for c in chans if c.psd is not None]
    if not rows:
        return False
    # skip DC (each segment's mean is removed) and Nyquist (a one-sided spectrum
    # does not double that bin, so it plots at half height)
    freq = np.fft.rfftfreq(nfft)[1:-1]
    psd = np.array([c.psd[1:-1] for c in rows])
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True,
                             gridspec_kw={"height_ratios": [1, 2.2], "hspace": 0.12})
    ax = axes[0]
    p10, p50, p90 = np.percentile(psd, [10, 50, 90], axis=0)
    ax.fill_between(freq, p10, p90, color=BLUE, alpha=0.12, lw=0)
    ax.plot(freq, p50, color=BLUE, lw=1.2)
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(LogFormatter())
    ax.yaxis.set_minor_formatter(LogFormatter(minor_thresholds=(1, 0.5)))
    ax.set_ylabel("PSD (ADC² per cycle/sample)")
    ax.set_title(f"Median of {len(rows)} channels (line), 10–90 % of channels (band)", loc="left")

    ax = axes[1]
    df = freq[1] - freq[0]
    im = ax.imshow(np.log10(psd), aspect="auto", cmap=SEQUENTIAL, interpolation="nearest",
                   extent=(freq[0] - df / 2, freq[-1] + df / 2, len(rows) - 0.5, -0.5))
    colorbar(fig, im, ax, "log10 PSD (ADC² per cycle/sample)")
    ax.grid(False)
    ax.tick_params(axis="y", length=0)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([c.name for c in rows], fontsize=5.5)
    for _, start, _ in board_bounds(rows)[1:]:
        ax.axhline(start - 0.5, color=SURFACE, lw=2)
    ax.set_title("Every channel", loc="left")
    ax.set_xlabel("frequency (cycles/sample)")
    fig.suptitle(f"Noise power spectrum ({nfft}-sample segments): {label}",
                 x=0.125, y=0.93, ha="left", fontsize=10)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_correlation(boards, path, label):
    ncols = min(4, len(boards))
    nrows = math.ceil(len(boards) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.3 * ncols, 3.5 * nrows), squeeze=False,
                             gridspec_kw={"hspace": 0.5, "wspace": 0.35})
    offdiag = np.concatenate([info["r"][~np.eye(CHANNELS_PER_BOARD, dtype=bool)] for info in boards])
    offdiag = np.abs(offdiag[np.isfinite(offdiag)])
    lim = max(0.1, math.ceil((offdiag.max() if offdiag.size else 0) / 0.05) * 0.05)
    im = None
    for ax, info in zip(axes.ravel(), boards):
        r = info["r"].copy()
        np.fill_diagonal(r, np.nan)
        im = ax.imshow(r, cmap=DIVERGING, vmin=-lim, vmax=lim)
        for i, j in itertools.product(range(CHANNELS_PER_BOARD), repeat=2):
            if np.isfinite(r[i, j]) and abs(r[i, j]) >= 0.1:
                ax.text(j, i, f"{r[i, j]:.2f}", ha="center", va="center", fontsize=6,
                        color=text_color_on(im.cmap(im.norm(r[i, j]))))
        ax.set_xticks(range(CHANNELS_PER_BOARD))
        ax.set_yticks(range(CHANNELS_PER_BOARD))
        ax.grid(False)
        ax.tick_params(length=0)
        ax.set_xlabel("channel")
        ax.set_ylabel("channel")
        row = info["row"]
        if row["corr_triggers"]:
            ax.set_title(f"board {info['board']}: mean r = {row['mean_r']:.3f}\n"
                         f"σ_coh = {row['sigma_coh']:.2f} ADC, {row['corr_triggers']} triggers",
                         loc="left")
        else:
            ax.set_title(f"board {info['board']}: no full-length triggers", loc="left", color=MUTED)
    for ax in axes.ravel()[len(boards):]:
        ax.axis("off")
    cb = fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.02, pad=0.03)
    cb.set_label("correlation coefficient r (numbers shown where |r| ≥ 0.1)")
    cb.outline.set_visible(False)
    fig.suptitle(f"Noise correlation between channels, full-length waveforms only: {label}",
                 x=0.125, y=0.99, ha="left", fontsize=10)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def ask_settings(args):
    """Ask for the ROOT file and the output folder (when no ROOT file was
    given); every other setting keeps its default."""
    enable_path_completion()
    print("DAQ_YY baseline analysis. Enter keeps the value in [brackets]; Tab completes paths.")
    args.input = ask("ROOT file", convert=existing_file)
    args.out = ask("Output folder", default=args.out or default_output(args.input, "baseline"),
                   convert=output_folder)
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", nargs="?", default=None,
                    help="ROOT file written by raw2root_new (asked for if left out)")
    ap.add_argument("--out", default=None,
                    help="Output directory (default: <run name>_baseline in the folder of the "
                         "input file)")
    ap.add_argument("--label", default=None,
                    help="Label for the plot titles (default: the run name)")
    ap.add_argument("--max-triggers", type=int, default=None,
                    help="Analyze only the first N triggers of every board (default: all)")
    ap.add_argument("--nsamps", type=int, default=None,
                    help="Full waveform length in samples (default: the most common length)")
    ap.add_argument("--nfft", type=int, default=1024,
                    help="Segment length in samples for the noise spectrum (default: 1024)")
    ap.add_argument("--spike-sigma", type=float, default=5.0,
                    help="Samples further than this many sigma from the pedestal are "
                         "spikes (default: 5)")
    ap.add_argument("--noisy-factor", type=float, default=1.5,
                    help="Flag NOISY/QUIET when sigma is above/below this factor times the "
                         "median of all channels (default: 1.5)")
    args = ap.parse_args()
    if args.input is None:
        ask_settings(args)

    run = run_name(args.input)
    out_dir = args.out or default_output(args.input, "baseline")
    label = args.label or run

    root_file = open_root(args.input)
    if root_file is None:
        print(f"ERROR: cannot open {args.input}", file=sys.stderr)
        sys.exit(1)
    print(f"Reading {args.input}")
    chans = open_channels(root_file)
    if not chans:
        print("ERROR: no ch<id> trees in the file.", file=sys.stderr)
        sys.exit(1)
    make_output_folder(out_dir)

    # ---- full waveform length: the most common nsamps, unless given ----
    lengths, n_len = np.unique(np.concatenate([c.nsamps for c in chans]), return_counts=True)
    full_length = args.nsamps or int(lengths[np.argmax(n_len)])

    # ---- loop over the triggers of every board ----
    groups = [(board, list(group)) for board, group in itertools.groupby(chans, key=lambda c: c.board)]
    cap = args.max_triggers if args.max_triggers is not None else math.inf
    bar = Progress("Reading triggers", sum(
        min(len(np.unique(np.concatenate([c.ttt for c in group]))), cap) for _, group in groups))
    boards = []
    for board, group in groups:
        n_trig, corr = analyze_board(group, full_length, args.nfft, args.max_triggers, bar)
        boards.append({"board": board, "chans": group, "n_triggers": n_trig, "corr": corr})
    bar.close()
    for c in chans:
        c.finish(args.spike_sigma, args.nfft)
    median_sigma = assign_flags(chans, args.noisy_factor)
    most_triggers = max(info["n_triggers"] for info in boards)
    board_rows = [summarize_board(info, most_triggers) for info in boards]

    # ---- print tables, and keep a copy for baseline_summary.txt ----
    lines = []

    def report(text=""):
        print(text)
        lines.append(text)

    statuses = [e["status"] for c in chans for e in c.events if e["status"] != "missing"]
    # does tbusy rise right after a short waveform? (count over consecutive waveforms)
    after_short, after_full = [0, 0], [0, 0]
    for c in chans:
        rec = [e for e in c.events if e["status"] != "missing"]
        for prev, nxt in zip(rec, rec[1:]):
            tally = after_full if prev["status"] == "ok" else after_short
            tally[0] += nxt["tbusy_step"] > 0
            tally[1] += 1
    report()
    report(f"File       : {args.input}")
    report(f"Channels   : {len(chans)} on {len(boards)} boards, "
           f"{most_triggers} triggers on the fullest board")
    report(f"Waveforms  : {len(statuses)} analyzed, {statuses.count('ok')} at the full length of "
           f"{full_length} samples")
    report(f"Noise      : median sigma {median_sigma:.3f} ADC over all channels")
    report(f"tbusy      : rose after {after_short[0]} of {after_short[1]} short waveforms and "
           f"after {after_full[0]} of {after_full[1]} full-length ones (next waveform, same channel)")
    report()
    header = (f"{'tree':<8} {'brd':>3} {'ch':>2} {'trig':>4} {'miss':>4} {'short':>5} "
              f"{'FitMean':>9} {'FitSigma':>8} {'EvSigma':>7} {'Spread':>6} {'DAQrms':>6} "
              f"{'out<':>4} {'out>':>4} {'line':>6} {'lineSig':>7} {'tbusy':>11}  flags")
    report(header)
    report("-" * len(header))
    for c in chans:
        report(f"{c.name:<8} {c.board:>3} {c.chn:>2} {c.n_triggers:>4} {c.n_missing:>4} "
               f"{c.n_short + c.n_long:>5} {c.mu:>9.2f} {c.sigma:>8.3f} {c.event_sigma:>7.3f} "
               f"{c.spread:>6.3f} {c.daq_rms:>6.3f} {c.n_low:>4} {c.n_high:>4} "
               f"{c.psd_line_freq:>6.3f} {c.psd_line_sigma:>7.1f} {c.tbusy_last:>11}  "
               f"{';'.join(c.flags) or 'OK'}")
    report()
    header = (f"{'board':>5} {'chans':>5} {'trig':>4} {'missing':>7} {'short':>5} "
              f"{'corrTrig':>8} {'mean_r':>7} {'max_r':>6} {'pair':>7} {'sigCoh':>6}  flags")
    report(header)
    report("-" * len(header))
    for row in board_rows:
        report(f"{row['board']:>5} {row['n_channels']:>5} {row['n_triggers']:>4} "
               f"{row['n_missing']:>7} {row['n_short']:>5} {row['corr_triggers']:>8} "
               f"{row['mean_r']:>7.3f} {row['max_r']:>6.3f} {row['max_r_pair']:>7} "
               f"{row['sigma_coh']:>6.3f}  {row['flags']}")
    report()
    report("Columns: FitMean/FitSigma = Gaussian fit to all samples; EvSigma = mean standard "
           "deviation of single waveforms;")
    report("         Spread = standard deviation of the waveform means; DAQrms = mean of the DAQ "
           "rms (first 20 samples);")
    report(f"         out</out> = samples more than {args.spike_sigma:g} sigma below/above the "
           "pedestal;")
    report("         line/lineSig = frequency (cycles/sample) and significance (sigma) of the "
           "strongest narrow line in the noise spectrum;")
    report("         tbusy = DAQ tbusy at the last waveform (raw counter); sigCoh = common-mode "
           "noise from the mean pair covariance.")
    used = sorted({f for c in chans for f in c.flags} |
                  {row["flags"] for row in board_rows if row["flags"] != "OK"})
    if used:
        report("Flags:")
        for flag in used:
            report(f"  {flag:<15} {FLAG_HELP[flag]}")
        n_flagged = sum(bool(c.flags) for c in chans)
        report(f"{n_flagged} of {len(chans)} channels flagged.")
    else:
        report("No flags: every channel passed every check.")
    report()

    # ---- write tables ----
    with open(os.path.join(out_dir, "baseline_summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")

    channel_rows = []
    for c in chans:
        fit = c.fit or {}
        channel_rows.append({
            "tree": c.name, "board": c.board, "channel": c.chn,
            "n_triggers": c.n_triggers, "n_missing": c.n_missing,
            "n_short": c.n_short, "n_long": c.n_long, "n_samples": c.n_samples,
            "raw_mean": c.raw_mean, "raw_std": c.raw_std,
            "fit_mean": fit.get("mean", math.nan), "fit_mean_err": fit.get("mean_err", math.nan),
            "fit_sigma": fit.get("sigma", math.nan), "fit_sigma_err": fit.get("sigma_err", math.nan),
            "fit_chi2": fit.get("chi2", math.nan), "fit_ndf": fit.get("ndf", ""),
            "event_sigma_mean": c.event_sigma, "pedestal_spread": c.spread,
            "daq_rms_mean": c.daq_rms, "daq_baseline_scatter": c.daq_baseline_scatter,
            "daq_area_mean": c.area_mean, "daq_area_std": c.area_std,
            "adc_min": c.adc_min, "adc_max": c.adc_max,
            "n_low_outliers": c.n_low, "n_high_outliers": c.n_high,
            "expected_outliers": c.expected_outliers, "n_saturated": c.n_saturated,
            "odd_fraction": c.odd_fraction, "psd_sigma": c.psd_sigma,
            "psd_line_freq": c.psd_line_freq, "psd_line_ratio": c.psd_line_ratio,
            "psd_line_sigma": c.psd_line_sigma,
            "tbusy_last": c.tbusy_last, "flags": ";".join(c.flags) or "OK",
        })
    write_csv(os.path.join(out_dir, "baseline_channels.csv"), channel_rows,
              list(channel_rows[0].keys()))

    event_columns = ["tree", "board", "channel", "trigger", "status", "entry", "ttt", "t2prev",
                     "tbusy", "tbusy_step", "nsamps", "mean", "std", "min", "max",
                     "daq_baseline", "daq_rms", "daq_area", "daq_xmin", "daq_xmax"]
    event_rows = sorted((e for c in chans for e in c.events),
                        key=lambda e: (e["board"], e["trigger"], e["channel"]))
    write_csv(os.path.join(out_dir, "baseline_events.csv"), event_rows, event_columns)

    write_csv(os.path.join(out_dir, "baseline_boards.csv"), board_rows, list(board_rows[0].keys()))

    pair_rows = []
    for info in boards:
        names = {c.chn: c.name for c in info["chans"]}
        for i, j in itertools.combinations(range(CHANNELS_PER_BOARD), 2):
            if info["corr"].triggers[i, j]:
                pair_rows.append({"board": info["board"], "channel_a": i, "channel_b": j,
                                  "tree_a": names[i], "tree_b": names[j],
                                  "r": float(info["r"][i, j]), "covariance": float(info["cov"][i, j]),
                                  "n_triggers": int(info["corr"].triggers[i, j])})
    write_csv(os.path.join(out_dir, "noise_correlation.csv"), pair_rows,
              ["board", "channel_a", "channel_b", "tree_a", "tree_b", "r", "covariance",
               "n_triggers"])

    psd_chans = [c for c in chans if c.psd is not None]
    if psd_chans:
        with open(os.path.join(out_dir, "noise_psd.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["freq_cycles_per_sample"] + [c.name for c in psd_chans])
            for k, fr in enumerate(np.fft.rfftfreq(args.nfft)):
                w.writerow([f"{fr:.6g}"] + [f"{c.psd[k]:.6g}" for c in psd_chans])

    # ---- plots ----
    def out(name):
        return os.path.join(out_dir, name)

    jobs = [("summary", partial(plot_summary, chans, median_sigma, args.noisy_factor,
                                out("summary.png"), label)),
            ("pedestal map", partial(plot_pedestal_map, chans, out("pedestal_map.png"), label))]
    jobs += [(f"histograms board {info['board']}",
              partial(plot_histograms, info["board"], info["chans"], args.spike_sigma,
                      out(f"pedestal_hist_board{info['board']}.png"), label)) for info in boards]
    jobs += [("stability", partial(plot_stability, chans, most_triggers,
                                   out("pedestal_stability.png"), label)),
             ("noise spectra", partial(plot_spectra, chans, args.nfft, out("noise_spectra.png"), label)),
             ("correlation", partial(plot_correlation, boards, out("noise_correlation.png"), label))]
    bar = Progress("Making plots", len(jobs))
    made = {}
    for note, job in jobs:
        bar.update(0, note)     # show what is being drawn before it starts
        made[note] = job()
        bar.update()
    bar.close()
    if made["noise spectra"] is False:
        print(f"WARNING: no waveform has {args.nfft} samples; no noise spectrum.", file=sys.stderr)
    print(f"Saved tables and plots to {out_dir}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("\nCancelled.")
