#!/usr/bin/env python3
"""
analyze_sync.py

Timing synchronization of the AD9694 digitizer boards, from the ROOT files
written by raw2root_new (ad9694daqtrigchain). One signal split to several
channels (e.g. a square wave through a splitter, self-trigger mode) gives every
channel a waveform of the same pulse; the time of its rising edge on every
channel is compared with that on the reference channel.

Usage:
    python3 analyze_sync.py              # asks for the file(s) and output folder
    python3 analyze_sync.py run.root
    python3 analyze_sync.py run.root --out /path/to/outdir
    python3 analyze_sync.py run1.root run2.root run3.root    # also compares the runs
    python3 analyze_sync.py run.root --ref 10010 --fraction 0.5 --phase-correction

Without a ROOT file on the command line the script asks for the file(s) (several
separated by spaces, or a * pattern) and the output folder; Enter keeps the
folder shown in [brackets] and Tab completes paths. Every other setting keeps
its default unless it is given as an option below. Needs daqyy_common.py in the
same folder.

Input: one TTree per channel, named ch<id> with id = 10000 + 10*board + channel.
Each entry is one waveform: chdata[nsamps] holds the ADC samples and ttt the
trigger time tag. pulse_pre comes from <run>.conf next to the ROOT file (XIAOLAB
DAQ keeps the run's configuration there), or from --pulse-pre.

Method (times in ns, jitter in ps; a sample and a ttt tick are both 2 ns):
    edge time    t = ttt * 2 ns + (x - 2*pulse_pre) * 2 ns. ttt counts the board's
                 clock from the start of the run (the boards start together and
                 step ttt by 2 ticks, i.e. 4 ns); sample 2*pulse_pre of the record
                 is the place ttt refers to; x is where the rising edge crosses
                 --fraction (default 0.5) of the pulse height, interpolated
                 linearly between the two samples around it. Baseline: mean of the
                 20 samples from 40 to 21 before sample 2*pulse_pre; height: the
                 highest sample from 10 before to 29 after it.
    edges used   the record holds those samples, its baseline is quiet (standard
                 deviation below 15 ADC), its height is above --min-height, and the
                 crossing lies from 10 samples before to 1 after sample
                 2*pulse_pre (a record cut at its start has its edge elsewhere).
    same pulse   for each edge of the reference channel, the edge of the other
                 channel with the nearest ttt, if closer than 0.4 of the pulse
                 period seen on the reference channel (at most 100 ticks, 200 ns).
    per channel  dt = t(channel) - t(reference) over the common pulses. Offset and
                 jitter are the mean and sigma of a Gaussian fit to the dt
                 histogram within 5 robust sigma of its median (with fewer than 100
                 common pulses: the mean and standard deviation there). The median
                 and the robust sigma (1.4826 x median absolute deviation) are
                 listed too.
    phase bias   linear interpolation moves the crossing by up to about +-25 ps,
                 depending on where the edge falls between two samples.
                 --phase-correction maps the fractional sample position of each
                 channel through its own distribution, as the true positions are
                 uniform for a pulser that is not locked to the ADC clock. It needs
                 a few thousand edges per channel; with fewer it adds noise.

Outputs (written into --out, default: <run name>_sync in the folder of the
(first) input file):
    sync_summary.txt        the tables printed to the terminal
    sync_channels.csv       one row per run and channel: edges, common pulses,
                            offset, jitter, their errors, median, robust sigma
    sync_offsets.png        offset of every channel to the reference, error
                            bars = jitter (one run)
    sync_dt_ch<id>.png      dt histogram with the Gaussian fit, one file per
                            channel (one run)
With several runs every run gets its own folder <out>/<run name>/ for its
sync_offsets.png and sync_dt_ch<id>.png, and <out> also holds
    sync_runs.png           every channel's offset per run, minus its mean over
                            the runs (ps), error bars = statistical error
"""

import os

# One thread is plenty here; keep numpy's BLAS from starting one per core.
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse
import csv
import glob
import math
import re
import sys
import textwrap

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.optimize import curve_fit

from daqyy_common import (
    SURFACE, INK, INK2, MUTED, AXIS, BLUE,
    open_root, channel_trees, run_name, default_output,
    make_output_folder, enable_path_completion, ask, output_folder, Progress,
)

NS_PER_TICK = 2.0             # ttt ticks and samples: 2 ns (500 MS/s)
BASE_FROM, BASE_TO = -40, -20  # baseline samples, relative to sample 2*pulse_pre
EDGE_FROM, EDGE_TO = -10, 1    # where the crossing may lie, relative to it
HEIGHT_TO = 30                 # the height is the highest sample up to 29 after it
QUIET_ADC = 15.0               # baseline standard deviation limit
MATCH_MAX = 100                # matching window: at most 100 ticks (200 ns)
FIT_MIN = 100                  # common pulses needed for a Gaussian fit
HIST_BINS = 50
# categorical slots of the dataviz reference palette, in their fixed order (one per channel)
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]


def gaussian(x, amplitude, mean, sigma):
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def short_run(name):
    """run00050 for run00050_file0000_20261107_133636."""
    m = re.match(r"run\d+", name)
    return m.group(0) if m else name


def chan_label(chid):
    return f"board {(chid - 10000) // 10} ch {(chid - 10000) % 10}"


# ---- settings of the run ----

def pulse_pre_of(path, given):
    """pulse_pre from --pulse-pre, else from <run>.conf next to the ROOT file. The
    board library uses 20 for values outside 5..1020, and so does this."""
    if given is not None:
        value, source = given, "--pulse-pre"
    else:
        conf = os.path.splitext(path)[0] + ".conf"
        try:
            with open(conf, encoding="utf-8", errors="replace") as f:
                text = f.read()
        except OSError:
            sys.exit(f"ERROR: no {os.path.basename(conf)} next to {path} to read pulse_pre from; "
                     f"give it with --pulse-pre")
        m = re.search(r"^\s*pulse_pre\s*=\s*(\d+)", text, re.M)
        if not m:
            sys.exit(f"ERROR: no pulse_pre in {conf}; give it with --pulse-pre")
        value, source = int(m.group(1)), os.path.basename(conf)
    if not 5 <= value <= 1020:
        return 20, f"{source} says {value}, outside 5..1020, so the board library used 20"
    return value, source


def note_of(path):
    """The run's note from <run>.conf (the '# ---- note ----' block XIAOLAB DAQ writes), or ''."""
    try:
        with open(os.path.splitext(path)[0] + ".conf", encoding="utf-8", errors="replace") as f:
            text = f.read()
    except OSError:
        return ""
    m = re.search(r"^# ---- note ----\n(.*?)^# ---- end of note ----", text, re.M | re.S)
    if not m:
        return ""
    return " · ".join(line[1:].strip() for line in m.group(1).splitlines() if line.strip("# "))


# ---- reading the edges ----

def read_windows(tree, pre2, max_waveforms, progress, note):
    """ttt and the samples around sample pre2 of every record that holds them
    (sorted by ttt), and the number of records too short for that."""
    first = max(pre2 + BASE_FROM, 0)
    last = pre2 + HEIGHT_TO
    n = tree.GetEntries() if max_waveforms is None else min(tree.GetEntries(), max_waveforms)
    tree.SetBranchStatus("*", 0)
    for name in ("ttt", "nsamps", "chdata"):
        tree.SetBranchStatus(name, 1)
    ttt = np.zeros(n, dtype=np.int64)
    win = np.zeros((n, last - first), dtype=np.uint16)
    keep = np.zeros(n, dtype=bool)
    for i in range(n):
        tree.GetEntry(i)
        if tree.nsamps >= last:
            win[i] = np.frombuffer(tree.chdata, dtype=np.uint16, count=last)[first:]
            ttt[i] = tree.ttt
            keep[i] = True
        if i % 2000 == 1999:
            progress.update(2000, note)
    progress.update(n % 2000, note)
    tree.SetBranchStatus("*", 1)
    order = np.argsort(ttt[keep], kind="stable")
    return ttt[keep][order], win[keep][order], n, int(n - keep.sum()), pre2 - first


def find_edges(win, at, fraction, min_height):
    """Crossing of the rising edge, in samples relative to sample 2*pulse_pre
    (index `at` of the windows), and which waveforms have a usable edge."""
    w = win.astype(np.float64)
    nbase = at + BASE_TO
    base = w[:, :nbase].mean(axis=1)
    quiet = w[:, :nbase].std(axis=1) < QUIET_ADC
    top = w[:, at + EDGE_FROM:].max(axis=1)
    high = top - base > min_height
    level = base + fraction * (top - base)
    i = np.argmax(w > level[:, None], axis=1)
    placed = (i >= at + EDGE_FROM) & (i <= at + EDGE_TO)
    k = np.arange(len(w))
    i1 = np.clip(i, 1, w.shape[1] - 1)
    x = (i1 - 1) + (level - w[k, i1 - 1]) / np.maximum(w[k, i1] - w[k, i1 - 1], 1e-9) - at
    return x, quiet, high, placed


def phase_corrected(x):
    """Map the fractional part of x through its own distribution (uniform when right)."""
    whole = np.floor(x)
    frac = x - whole
    ranked = np.searchsorted(np.sort(frac), frac, side="right")
    return whole + (ranked - 0.5) / len(frac)


class Channel:
    def __init__(self, chid):
        self.chid = chid
        self.board, self.chn = (chid - 10000) // 10, (chid - 10000) % 10
        self.name = f"ch{chid}"
        self.records = self.short = self.noisy = self.low = self.misplaced = 0
        self.ttt = np.zeros(0, dtype=np.int64)
        self.x = np.zeros(0)

    @property
    def t_ns(self):
        return self.ttt * NS_PER_TICK + self.x * NS_PER_TICK


def read_run(path, args):
    """{chid: Channel} with the usable edges of every channel of one ROOT file."""
    pulse_pre, source = pulse_pre_of(path, args.pulse_pre)
    pre2 = 2 * pulse_pre
    if pre2 + BASE_TO < 8:
        sys.exit(f"ERROR: {path}: pulse_pre {pulse_pre} leaves too few samples before the edge "
                 f"for a baseline (needs pulse_pre >= 14)")
    root_file = open_root(path)
    if root_file is None:
        sys.exit(f"ERROR: cannot open {path}")
    trees = [name for chans in channel_trees(root_file).values() for name in chans.values()]
    if not trees:
        sys.exit(f"ERROR: no ch<id> trees in {path}")
    total = sum(min(root_file.Get(name).GetEntries(), args.max_waveforms or math.inf) for name in trees)
    bar = Progress(f"Reading {short_run(run_name(path))}", total)
    chans = {}
    for name in trees:
        tree = root_file.Get(name)
        ch = Channel(int(name[2:]))
        ttt, win, ch.records, ch.short, at = read_windows(tree, pre2, args.max_waveforms, bar, name)
        if len(ttt):
            x, quiet, high, placed = find_edges(win, at, args.fraction, args.min_height)
            good = quiet & high & placed
            ch.noisy = int((~quiet).sum())
            ch.low = int((quiet & ~high).sum())
            ch.misplaced = int((quiet & high & ~placed).sum())
            ch.ttt, ch.x = ttt[good], x[good]
            if args.phase_correction and len(ch.x):
                ch.x = phase_corrected(ch.x)
        chans[ch.chid] = ch
    bar.close()
    root_file.Close()
    return chans, pulse_pre, source


# ---- comparing with the reference ----

def match(ref_ttt, ttt, window):
    """Indices (into the reference, into the other channel) of the common pulses."""
    if len(ttt) == 0 or len(ref_ttt) == 0:
        return np.zeros(0, dtype=int), np.zeros(0, dtype=int)
    j = np.searchsorted(ttt, ref_ttt)
    lo, hi = np.clip(j - 1, 0, len(ttt) - 1), np.clip(j, 0, len(ttt) - 1)
    pick = np.where(np.abs(ttt[lo] - ref_ttt) <= np.abs(ttt[hi] - ref_ttt), lo, hi)
    ok = np.abs(ttt[pick] - ref_ttt) < window
    return np.flatnonzero(ok), pick[ok]


def matching_window(ref_ttt):
    """0.4 of the pulse period seen on the reference channel, at most MATCH_MAX ticks."""
    steps = np.diff(ref_ttt)
    steps = steps[steps > 0]
    period = float(np.median(steps)) if len(steps) else math.inf
    return min(MATCH_MAX, 0.4 * period), period


def dt_stats(dt):
    """Offset and jitter of dt (ns): Gaussian fit within 5 robust sigma of the
    median, or mean and standard deviation there with fewer than FIT_MIN values."""
    out = {"n": len(dt), "fit": False}
    if len(dt) == 0:
        return out
    med = float(np.median(dt))
    rsig = 1.4826 * float(np.median(np.abs(dt - med)))
    out.update(median=med, robust_sigma=rsig)
    near = dt[np.abs(dt - med) <= 5 * rsig] if rsig > 0 else dt
    out["outliers"] = len(dt) - len(near)
    mean = float(near.mean())
    std = float(near.std(ddof=1)) if len(near) > 1 else math.nan
    out.update(offset=mean, jitter=std, offset_err=std / math.sqrt(len(near)) if len(near) > 1 else math.nan,
               jitter_err=std / math.sqrt(2 * (len(near) - 1)) if len(near) > 1 else math.nan)
    if len(near) < FIT_MIN or not rsig > 0:
        return out
    counts, edges = np.histogram(near, bins=HIST_BINS, range=(med - 5 * rsig, med + 5 * rsig))
    centers = 0.5 * (edges[:-1] + edges[1:])
    sel = counts > 0
    out.update(counts=counts, edges=edges)
    if np.count_nonzero(sel) < 4:
        return out
    try:
        popt, pcov = curve_fit(gaussian, centers[sel], counts[sel], p0=[counts.max(), mean, std],
                               sigma=np.sqrt(counts[sel]), absolute_sigma=True, maxfev=10000)
    except (RuntimeError, ValueError):
        return out
    perr = np.sqrt(np.diag(pcov))
    model = gaussian(centers[sel], *popt)
    out.update(fit=True, offset=float(popt[1]), offset_err=float(perr[1]), jitter=abs(float(popt[2])),
               jitter_err=float(perr[2]), popt=popt,
               chi2=float((((counts[sel] - model) / np.sqrt(counts[sel])) ** 2).sum()), ndf=int(sel.sum()) - 3)
    return out


def compare(chans, ref_chid):
    ref = chans[ref_chid]
    window, period = matching_window(ref.ttt)
    results = {}
    for chid, ch in chans.items():
        if chid == ref_chid:
            continue
        i_ref, i_ch = match(ref.ttt, ch.ttt, window)
        dt = ch.t_ns[i_ch] - ref.t_ns[i_ref]
        results[chid] = dt_stats(dt)
        results[chid]["dt"] = dt
    return results, window, period


# ---- plots ----

def plot_offsets(chans, ref_chid, results, path, title):
    chids = sorted(chans)
    fig, ax = plt.subplots(figsize=(max(6.5, 1.0 * len(chids) + 2), 4.2))
    xs = np.arange(len(chids))
    offsets = [0.0 if c == ref_chid else results[c].get("offset", math.nan) for c in chids]
    jitters = [0.0 if c == ref_chid else results[c].get("jitter", math.nan) for c in chids]
    ax.axhline(0, color=AXIS, lw=1)
    for x, c, off, jit in zip(xs, chids, offsets, jitters):
        if c == ref_chid:
            ax.plot(x, 0, "o", ms=9, mfc=SURFACE, mec=INK2, mew=1.5)
            ax.annotate("reference", (x, 0), xytext=(0, 10), textcoords="offset points", ha="center",
                        fontsize=7, color=INK2)
            continue
        if math.isnan(off):
            ax.annotate("no common\npulses", (x, 0), xytext=(0, 10), textcoords="offset points",
                        ha="center", fontsize=7, color=MUTED)
            continue
        ax.errorbar(x, off, yerr=jit, fmt="o", ms=8, color=BLUE, ecolor=BLUE, elinewidth=1.5, capsize=4)
        ax.annotate(f"{off:+.3f} ns\nσ {jit * 1e3:.0f} ps", (x, off), xytext=(0, 9), textcoords="offset points",
                    ha="center", fontsize=7, color=INK)
    finite = [o for o in offsets if not math.isnan(o)]
    lo, hi = min(finite + [0.0]), max(finite + [0.0])
    pad = max(0.25 * (hi - lo), 0.5)
    ax.set_ylim(lo - pad, hi + 1.6 * pad)
    ax.set_xticks(xs, [f"{chan_label(c)}\nch{c}" for c in chids])
    ax.set_xlim(-0.6, len(chids) - 0.4)
    ax.set_ylabel(f"edge time offset to ch{ref_chid} (ns)")
    ax.grid(axis="x", visible=False)
    ax.set_title(title, loc="left", fontsize=9)
    ax.text(0.01, 0.02, "error bars: jitter σ (event by event), too small to see at this scale when it is tens of ps",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=7, color=INK2)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_dt(chid, ref_chid, res, path, run):
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    dt = res["dt"]
    med, rsig = res["median"], res["robust_sigma"]
    if "counts" in res:
        counts, edges = res["counts"], res["edges"]
    else:                                   # few pulses: a coarse histogram of all of them
        span = max(5 * rsig, 0.01)
        counts, edges = np.histogram(dt, bins=min(HIST_BINS, max(len(dt) // 3, 5)), range=(med - span, med + span))
    ax.stairs(counts, edges, fill=True, color=BLUE, alpha=0.12, lw=0)
    ax.stairs(counts, edges, color=BLUE, lw=1.2)
    handles = [Patch(facecolor=to_rgba(BLUE, 0.12), edgecolor=BLUE, label=f"{res['n']} common pulses")]
    if res["fit"]:
        xx = np.linspace(edges[0], edges[-1], 400)
        ax.plot(xx, gaussian(xx, *res["popt"]), color=INK, lw=1)
        handles.append(Line2D([], [], color=INK, lw=1, label=f"Gaussian fit, χ²/ndf = {res['chi2']:.0f}/{res['ndf']}"))
        what = "fit"
    else:
        what = "mean and standard deviation"
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0, max(counts.max(), 1) * 1.25)
    ax.set_xlabel(f"Δt = t(ch{chid}) − t(ch{ref_chid})  (ns)")
    ax.set_ylabel("common pulses per bin")
    ax.legend(handles=handles, loc="upper right")
    ax.set_title(f"{chan_label(chid)} (ch{chid}) vs {chan_label(ref_chid)} (ch{ref_chid}): {run}\n"
                 f"offset {res['offset']:+.4f} ns ± {res['offset_err'] * 1e3:.1f} ps    "
                 f"jitter σ = {res['jitter'] * 1e3:.1f} ± {res['jitter_err'] * 1e3:.1f} ps  ({what})",
                 loc="left", fontsize=8)
    outside = res["n"] - int(counts.sum())
    if outside:
        ax.text(0.02, 0.97, f"{outside} outside view", transform=ax.transAxes, va="top", fontsize=7, color=INK2)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def run_note_lines(note, width=14, lines=3):
    """A run's note in up to `lines` short lines for a tick label."""
    parts = textwrap.wrap(note, width)
    return parts[:lines - 1] + [parts[lines - 1] + ("…" if len(parts) > lines else "")] if len(parts) >= lines else parts


def plot_runs(runs, chids, ref_chid, table, means, notes, path):
    """Each channel's offset per run minus its mean over the runs, in ps."""
    fig, ax = plt.subplots(figsize=(max(6.5, 1.3 * len(runs) + 3), 4.8))
    xs = np.arange(len(runs))
    ax.axhline(0, color=AXIS, lw=1)
    k_total = len(chids)
    for k, chid in enumerate(chids):
        color = CATEGORICAL[k % len(CATEGORICAL)]
        dodge = (k - (k_total - 1) / 2) * min(0.08, 0.6 / max(k_total, 1))
        y = np.array([(table[r][chid]["offset"] - means[chid]) * 1e3 if chid in table[r] and "offset" in table[r][chid]
                      else math.nan for r in runs])
        e = np.array([table[r][chid]["offset_err"] * 1e3 if chid in table[r] and "offset" in table[r][chid]
                      else math.nan for r in runs])
        ax.plot(xs + dodge, y, color=color, lw=1.5, alpha=0.6)
        ax.errorbar(xs + dodge, y, yerr=e, fmt="o", ms=6, color=color, mec=SURFACE, mew=1.2, elinewidth=1.2,
                    capsize=2, label=f"{chan_label(chid)} (ch{chid}): mean {means[chid]:+.3f} ns")
    ax.set_xticks(xs, [r + "".join("\n" + line for line in run_note_lines(notes.get(r, ""))) for r in runs],
                  fontsize=7)
    ax.set_xlim(-0.5, len(runs) - 0.5)
    ax.grid(axis="x", visible=False)
    ax.set_ylabel(f"offset to ch{ref_chid} − its mean over the runs (ps)")
    ax.set_title(f"Edge time offset of every channel to ch{ref_chid}, per run", loc="left", fontsize=9)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---- output ----

def write_csv(path, rows, columns):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in row.items()})


def fmt(v, scale=1.0, digits=3, signed=False):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "-"
    return f"{v * scale:{'+' if signed else ''}.{digits}f}"


def files_from_text(text):
    """One or more ROOT files, separated by spaces; * patterns are expanded."""
    return files_from_list(text.split())


def files_from_list(parts):
    """The ROOT files named in parts (each a path or a * pattern)."""
    paths = []
    for part in parts:
        part = os.path.expanduser(part)
        found = sorted(glob.glob(part)) if any(c in part for c in "*?[") else [part]
        if not found:
            raise ValueError(f"nothing matches {part}")
        for p in found:
            p = os.path.abspath(p)
            if not os.path.isfile(p):
                raise ValueError(f"no such file: {p}")
            paths.append(p)
    if not paths:
        raise ValueError("give at least one ROOT file")
    return paths


def ask_settings(args):
    """Ask for the ROOT file(s) and the output folder (when no ROOT file was
    given); every other setting keeps its default."""
    enable_path_completion()
    print("DAQ_YY timing sync analysis. Enter keeps the value in [brackets]; Tab completes paths.")
    args.inputs = ask("ROOT file(s) (several: separate with spaces, or use a * pattern)", convert=files_from_text)
    args.out = ask("Output folder", default=args.out or default_output(args.inputs[0], "sync"),
                   convert=output_folder)
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="*", default=None,
                    help="ROOT file(s) written by raw2root_new (asked for if left out); several are compared")
    ap.add_argument("--out", default=None,
                    help="Output directory (default: <run name>_sync in the folder of the first input file)")
    ap.add_argument("--ref", type=int, default=None,
                    help="Reference channel id, e.g. 10010 (default: the lowest id, i.e. the first board)")
    ap.add_argument("--fraction", type=float, default=0.5,
                    help="Fraction of the pulse height where the edge time is taken (default: 0.5)")
    ap.add_argument("--min-height", type=float, default=500.0,
                    help="Smallest pulse height above the baseline, in ADC counts (default: 500)")
    ap.add_argument("--pulse-pre", type=int, default=None,
                    help="pulse_pre of the run (default: from <run>.conf next to the ROOT file)")
    ap.add_argument("--phase-correction", action="store_true",
                    help="Remove the phase bias of the linear interpolation (needs thousands of edges per channel)")
    ap.add_argument("--max-waveforms", type=int, default=None,
                    help="Read only the first N waveforms of every channel (default: all)")
    args = ap.parse_args()
    if not args.inputs:
        ask_settings(args)
    else:
        try:
            args.inputs = files_from_list(args.inputs)
        except ValueError as err:
            sys.exit(f"ERROR: {err}")
    if not 0.05 <= args.fraction <= 0.95:
        sys.exit("ERROR: --fraction must be from 0.05 to 0.95")

    several = len(args.inputs) > 1
    out_dir = args.out or default_output(args.inputs[0], "sync")
    make_output_folder(out_dir)

    lines = []

    def report(text=""):
        print(text)
        lines.append(text)

    rows, table, runs, notes, ref_chid = [], {}, [], {}, args.ref
    for path in args.inputs:
        name = run_name(path)
        print(f"Reading {path}")
        chans, pulse_pre, source = read_run(path, args)
        if ref_chid is None:
            ref_chid = min(chans)
        report()
        report(f"File       : {path}")
        report(f"pulse_pre  : {pulse_pre} ({source}): ttt refers to sample {2 * pulse_pre} of each record")
        note = note_of(path)
        if note:
            report(f"Note       : {note}")
        if ref_chid not in chans or len(chans[ref_chid].ttt) == 0:
            report(f"Reference  : ch{ref_chid} has no usable edges in this file: run skipped")
            continue
        results, window, period = compare(chans, ref_chid)
        report(f"Reference  : ch{ref_chid} ({chan_label(ref_chid)}), {len(chans[ref_chid].ttt)} usable edges")
        report(f"Same pulse : nearest ttt within {window:.0f} ticks ({window * NS_PER_TICK:.0f} ns); "
               f"pulse period on the reference {period * NS_PER_TICK:.0f} ns")
        report(f"Edge time  : {args.fraction:.0%} of the pulse height, linear interpolation"
               + (", phase bias corrected" if args.phase_correction else ""))
        report()
        header = (f"{'channel':<8} {'brd':>3} {'ch':>2} {'records':>8} {'usable':>8} {'short':>6} {'noisy':>6} "
                  f"{'low':>6} {'cut':>6} {'common':>7} {'offset(ns)':>11} {'±(ps)':>6} {'jitter(ps)':>10} "
                  f"{'±(ps)':>6} {'median(ns)':>10} {'robσ(ps)':>8} {'>5σ':>4}")
        report(header)
        report("-" * len(header))
        for chid in sorted(chans):
            ch = chans[chid]
            r = results.get(chid, {})
            line = (f"{ch.name:<8} {ch.board:>3} {ch.chn:>2} {ch.records:>8} {len(ch.ttt):>8} {ch.short:>6} "
                    f"{ch.noisy:>6} {ch.low:>6} {ch.misplaced:>6} ")
            if chid == ref_chid:
                line += f"{'-':>7} {'reference':>11}"
            else:
                line += (f"{r.get('n', 0):>7} {fmt(r.get('offset'), 1, 4, True):>11} {fmt(r.get('offset_err'), 1e3, 1):>6} "
                         f"{fmt(r.get('jitter'), 1e3, 1):>10} {fmt(r.get('jitter_err'), 1e3, 1):>6} "
                         f"{fmt(r.get('median'), 1, 4, True):>10} {fmt(r.get('robust_sigma'), 1e3, 1):>8} "
                         f"{r.get('outliers', 0):>4}")
            report(line)
            rows.append({"run": name, "reference": f"ch{ref_chid}", "channel": ch.name, "board": ch.board,
                         "ch": ch.chn, "records": ch.records, "usable_edges": len(ch.ttt), "too_short": ch.short,
                         "noisy_baseline": ch.noisy, "low_height": ch.low, "edge_out_of_place": ch.misplaced,
                         "common_pulses": r.get("n", "") if chid != ref_chid else "",
                         "offset_ns": r.get("offset", "") if chid != ref_chid else 0.0,
                         "offset_err_ps": r.get("offset_err", math.nan) * 1e3 if "offset_err" in r else "",
                         "jitter_ps": r.get("jitter", math.nan) * 1e3 if "jitter" in r else "",
                         "jitter_err_ps": r.get("jitter_err", math.nan) * 1e3 if "jitter_err" in r else "",
                         "median_ns": r.get("median", ""), "robust_sigma_ps": r.get("robust_sigma", math.nan) * 1e3
                         if "robust_sigma" in r else "", "outliers_5sigma": r.get("outliers", ""),
                         "gaussian_fit": ("yes" if r.get("fit") else "no") if chid != ref_chid else "",
                         "phase_correction": "yes" if args.phase_correction else "no"})
        report()
        report("records: waveforms read; usable: with a usable edge; short: too short to hold the samples "
               "around the edge; noisy: baseline not quiet; low: height below --min-height; cut: edge not "
               f"where it belongs (record cut at its start); common: pulses also seen on ch{ref_chid}; "
               ">5σ: common pulses further than 5 robust σ from the median (left out of offset and jitter)")

        plot_dir = os.path.join(out_dir, name) if several else out_dir
        make_output_folder(plot_dir)
        title = f"Edge time offset to ch{ref_chid}: {name}"
        plot_offsets(chans, ref_chid, results, os.path.join(plot_dir, "sync_offsets.png"), title)
        for chid, res in results.items():
            if res.get("n", 0) > 0:
                plot_dt(chid, ref_chid, res, os.path.join(plot_dir, f"sync_dt_ch{chid}.png"), name)
        runs.append(short_run(name))
        table[runs[-1]] = results
        notes[runs[-1]] = note

    if not runs:
        sys.exit("ERROR: no run had usable edges on the reference channel")

    if several and len(runs) > 1:
        chids = sorted({c for r in runs for c in table[r]})
        means = {}
        for chid in chids:
            vals = [table[r][chid]["offset"] for r in runs if chid in table[r] and "offset" in table[r][chid]]
            means[chid] = float(np.mean(vals)) if vals else math.nan
        report()
        report(f"Offset to ch{ref_chid} per run (ns); spread = largest minus smallest over the runs (ps)")
        header = f"{'channel':<8} {'brd':>3} " + " ".join(f"{r:>10}" for r in runs) + f" {'mean':>9} {'spread':>7}"
        report(header)
        report("-" * len(header))
        for chid in chids:
            vals = [table[r][chid].get("offset", math.nan) if chid in table[r] else math.nan for r in runs]
            good = [v for v in vals if not math.isnan(v)]
            spread = (max(good) - min(good)) * 1e3 if len(good) > 1 else math.nan
            report(f"ch{chid:<6} {(chid - 10000) // 10:>3} " + " ".join(f"{fmt(v, 1, 3, True):>10}" for v in vals)
                   + f" {fmt(means[chid], 1, 3, True):>9} {fmt(spread, 1, 0):>7}")
        report()
        report("Jitter σ per run (ps)")
        header = f"{'channel':<8} {'brd':>3} " + " ".join(f"{r:>10}" for r in runs)
        report(header)
        report("-" * len(header))
        for chid in chids:
            vals = [table[r][chid].get("jitter", math.nan) if chid in table[r] else math.nan for r in runs]
            report(f"ch{chid:<6} {(chid - 10000) // 10:>3} " + " ".join(f"{fmt(v, 1e3, 1):>10}" for v in vals))
        if any(notes.values()):
            report()
            report("Runs and their notes")
            for r in runs:
                report(f"  {r}: {notes[r] or '-'}")
        plot_runs(runs, chids, ref_chid, table, means, notes, os.path.join(out_dir, "sync_runs.png"))

    with open(os.path.join(out_dir, "sync_summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    columns = ["run", "reference", "channel", "board", "ch", "records", "usable_edges", "too_short",
               "noisy_baseline", "low_height", "edge_out_of_place", "common_pulses", "offset_ns", "offset_err_ps",
               "jitter_ps", "jitter_err_ps", "median_ns", "robust_sigma_ps", "outliers_5sigma", "gaussian_fit",
               "phase_correction"]
    write_csv(os.path.join(out_dir, "sync_channels.csv"), rows, columns)
    print()
    print(f"Outputs in {out_dir}")


if __name__ == "__main__":
    main()
