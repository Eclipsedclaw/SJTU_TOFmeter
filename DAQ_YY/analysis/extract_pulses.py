#!/usr/bin/env python3
"""
extract_pulses.py  -  step 1 of the paddle muon analysis

Measures the pulse in every waveform of a DAQ_YY run (the ROOT file written by
raw2root_new, one tree ch<id> per channel) and writes one entry per waveform to
pulses.root in the analysis folder, with a table per channel in
pulses_summary.txt.

Usage:
    python3 extract_pulses.py                  # asks for the ROOT file and the output folder
    python3 extract_pulses.py run.root         # output: <run name>_muon next to run.root
    python3 extract_pulses.py run.root --out /path/to/folder
    python3 extract_pulses.py run.root --setup setups/two_layer_2026-11-12.toml --max-waveforms 1000

Without a ROOT file on the command line the script asks for it and for the
output folder; Enter keeps the folder shown in [brackets], Tab completes paths.

Settings: the [pulse] section of the setup file (setups/*.toml; by default the
one valid for the run's date). The run's .conf next to the ROOT file gives
pulse_pre (the trigger's place in each record), the polarity of every channel
and the thresholds.

Method, for every waveform (samples are 2 ns; s = the trigger sample, 2 x pulse_pre):
    baseline     mean of samples 0 .. s + baseline_end (default s - 20); its rms
                 is the noise; 'slope' = mean of its last 30 samples - first 30
    height       highest running mean of `smooth` samples (default 5) within
                 s + search[0] .. s + search[1] (default -30 .. +50), above the
                 baseline in the pulse direction
    crossing     where the rising edge crosses cfd_fraction (default 0.2) of the
                 height, interpolated linearly between two samples ("constant
                 fraction"): the same point of every pulse, whatever its height
    time         timing = "fit" (default): the rising edge fitted with the channel's
                 own template. The template is the mean edge of the channel's first
                 template_pulses clean pulses (no flag, height >= template_min_height),
                 each divided by its height and placed at its crossing, in bins of 0.1
                 sample, smoothed by a Gaussian of 0.2 sample against the bins' noise,
                 then shifted to cross cfd_fraction exactly at 0. Every pulse is
                 fitted by least squares with height x template(t - t0) over the samples
                 fit_window (default -4 .. +4) around its crossing: the rising edge from
                 before it starts to near its top, not the flat top (its shape varies
                 from pulse to pulse). The height follows in closed form for every trial
                 t0 on a grid of 0.025 samples within fit_max_shift of the crossing; a
                 parabola through the best three refines t0. The time is t0, where the
                 fitted edge reaches cfd_fraction of the fitted height. Without a fit
                 (window past the end of the record, best t0 at the grid's limit, or
                 fewer than template_min_pulses clean pulses) the time is the crossing
                 and the flag NO_FIT is set. timing = "cfd": the time is the crossing.
                 t_ns = (ttt - first ttt of the run) x 2 ns + (time - s) x 2 ns
                 - the board and channel offsets of the setup file
    rise_ns      10 % to 90 % of the height on the rising edge
    charge       sum of the samples from the edge + charge_window, x 2 ns (ADC x ns)
    tests        flat top: a real pulse of this front end stays high after its
                 peak (mean over plateau_window / height >= min_plateau); noise:
                 no swing below the baseline deeper than max_undershoot x height;
                 pile-up: no later rise of 4 samples above pileup_fraction x height
                 and pileup_min ADC; saturation: no sample at the ADC limits (0, 16383)
                 or beyond max_raw_adc / min_raw_adc (where the amplifier clips)
    flags        bit mask, see pulses_summary.txt; SHORT, NO_EDGE, NOISE and LOW
                 pulses are left out by the later steps

Outputs (in the analysis folder):
    pulses.root          tree 'pulses' (one entry per waveform: edge = the time in samples,
                         edge_cfd = the crossing, fit_amp and fit_rms = the fit's height and
                         rms residual in ADC), tree 'templates' (each channel's edge
                         template), tree 'shapes' (average pulse shape per channel), and
                         'info' (settings)
    pulses_summary.txt   one line per channel: waveforms, flags, rate, noise; then the
                         templates and fits
"""

import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse
import math
import sys

import numpy as np

import muon_common as mc

CHUNK = 4096          # waveforms measured at a time
SHAPE_KEEP = 500      # waveforms per channel kept for the average pulse shape
SHAPE_FROM, SHAPE_TO = -20, 150   # samples around the edge in the shape
TEMPLATE_STEP = 0.1   # samples: bin width of the edge templates
TEMPLATE_SMOOTH = 0.2  # samples: Gaussian smoothing of a template against its bin noise (the edge spans 6)
FIT_STEP = 0.025      # samples: grid of trial times of the edge fit
FIT_BLOCK = 20000     # pulses fitted at a time (memory)


def window_length(ref, ps):
    """Samples of each record that the measurement uses."""
    s1 = int(ps["search"][1])
    k = int(ps["smooth"])
    return ref + s1 + k + max(int(ps["plateau_window"][1]), int(ps["charge_window"][1]),
                              SHAPE_TO) + 5


def lowpass(W, mhz, kind):
    """Zero-phase low-pass filter of every row of W (2 ns samples), NaN past the end of a
    record kept. "gauss": Gaussian smoothing whose response falls to -3 dB at `mhz` (no
    overshoot, no ringing of its own); "butter": 4th-order Butterworth with that design
    frequency, run forwards and backwards. Both leave a pulse's timing unshifted."""
    X = np.array(W, dtype=np.float64, copy=True)
    bad = ~np.isfinite(X)
    if bad.any():   # continue each record with its last sample, so the filter sees no NaN
        idx = np.where(~bad, np.arange(X.shape[1])[None, :], 0)
        np.maximum.accumulate(idx, axis=1, out=idx)
        X = X[np.arange(len(X))[:, None], idx]
        X[~np.isfinite(X)] = 0.0
    if kind == "gauss":
        from scipy.ndimage import gaussian_filter1d
        sigma = math.sqrt(math.log(2)) / (2 * math.pi * mhz * 1e6) / (mc.NS_PER_SAMPLE * 1e-9)
        Y = gaussian_filter1d(X, sigma, axis=1, mode="nearest")
    elif kind == "butter":
        from scipy.signal import butter, sosfiltfilt
        sos = butter(4, mhz * 1e6, fs=1e9 / mc.NS_PER_SAMPLE, output="sos")
        Y = sosfiltfilt(sos, X, axis=1)
    else:
        raise ValueError(f"lowpass_kind must be gauss or butter, not {kind}")
    Y[bad] = np.nan
    return Y


def measure(W, nsamps, ref, sign, ps):
    """Measure a chunk of waveforms. W: (n, L) float64, NaN past the end of a record.
    With [pulse] lowpass_mhz > 0 every measurement uses the filtered waveform, except the
    baseline noise (always the raw samples) and, with lowpass_timing = false, the time."""
    n, L = W.shape
    r = np.arange(n)
    nb = ref + int(ps["baseline_end"])
    k = int(ps["smooth"])
    lo, hi = ref + int(ps["search"][0]), ref + int(ps["search"][1])
    flags = np.zeros(n, dtype=np.int32)
    short = nsamps < hi + k
    flags[short] |= mc.FLAG_SHORT

    brms = W[:, :nb].std(axis=1)                   # noise of the raw samples, filter or not
    sat = ((W >= min(float(ps["max_raw_adc"]), mc.ADC_MAX)) | (W <= max(float(ps["min_raw_adc"]), 0))).any(axis=1)
    mhz = float(ps["lowpass_mhz"])
    raw = W
    if mhz > 0:
        W = lowpass(W, mhz, str(ps["lowpass_kind"]))
    base = W[:, :nb].mean(axis=1)
    slope = W[:, nb - 30:nb].mean(axis=1) - W[:, :30].mean(axis=1)

    y = sign * (W - base[:, None])
    y_time = y if (mhz <= 0 or ps["lowpass_timing"]) else sign * (raw - base[:, None])
    valid = np.isfinite(y)
    y0 = np.where(valid, y, 0.0)
    cs = np.concatenate([np.zeros((n, 1)), np.cumsum(y0, axis=1)], axis=1)       # cs[:, i] = sum y[:, :i]
    cn = np.concatenate([np.zeros((n, 1)), np.cumsum(valid, axis=1)], axis=1)
    sm = (cs[:, k:] - cs[:, :-k]) / k                                              # mean of y[i:i+k]
    sm = np.where(cn[:, k:] - cn[:, :-k] == k, sm, -np.inf)

    seg = sm[:, lo:hi + 1]
    ipk = lo + np.argmax(seg, axis=1)
    height = sm[r, ipk]
    height = np.where(np.isfinite(height), height, np.nan)
    ic = ipk + k // 2

    def crossing(level):
        a = max(lo - 15, 1)
        b = min(hi + k, L - 1)
        cols = np.arange(a, b)
        part = y_time[:, a:b]
        below = (part <= level[:, None]) & (cols[None, :] < ic[:, None])
        has = below.any(axis=1)
        j = a + (below.shape[1] - 1 - np.argmax(below[:, ::-1], axis=1))
        ya = y_time[r, j]
        yb = y_time[r, np.minimum(j + 1, L - 1)]
        with np.errstate(invalid="ignore", divide="ignore"):
            x = j + (level - ya) / (yb - ya)
        ok = has & (yb > ya) & np.isfinite(x)
        return np.where(ok, x, np.nan)

    frac = float(ps["cfd_fraction"])
    edge = crossing(frac * height)
    rise = (crossing(0.9 * height) - crossing(0.1 * height)) * mc.NS_PER_SAMPLE

    def window_mean(start, stop):
        start = np.clip(start, 0, L)
        stop = np.clip(stop, 0, L)
        cnt = cn[r, stop] - cn[r, start]
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(cnt > 0, (cs[r, stop] - cs[r, start]) / cnt, np.nan), cnt

    p0, p1 = int(ps["plateau_window"][0]), int(ps["plateau_window"][1])
    plat_mean, _ = window_mean(ipk + p0, ipk + p1)
    with np.errstate(invalid="ignore", divide="ignore"):
        plateau = plat_mean / height
        under = -np.nanmin(np.where(valid[:, lo:hi + k + 10], y[:, lo:hi + k + 10], np.inf), axis=1) / height

    e0 = np.where(np.isfinite(edge), np.floor(np.nan_to_num(edge)).astype(np.int64), ref)
    c0, c1 = int(ps["charge_window"][0]), int(ps["charge_window"][1])
    cmean, ccount = window_mean(e0 + c0, e0 + c1)
    charge = np.where(np.isfinite(edge), cmean * ccount * mc.NS_PER_SAMPLE, np.nan)

    with np.errstate(invalid="ignore"):       # -inf - -inf where a record has ended
        d4 = sm[:, 4:] - sm[:, :-4]
    cols = np.arange(d4.shape[1])
    d4 = np.where((cols[None, :] >= (ipk + 10)[:, None]) & np.isfinite(d4), d4, -np.inf)
    rise_later = d4.max(axis=1)
    pileup = rise_later > np.maximum(float(ps["pileup_fraction"]) * height, float(ps["pileup_min"]))

    ok_height = np.isfinite(height)
    flags[~short & (~ok_height | (height < float(ps["min_height"])))] |= mc.FLAG_LOW
    flags[~short & ok_height & ~np.isfinite(edge)] |= mc.FLAG_NO_EDGE
    noise = np.zeros(n, dtype=bool)
    if float(ps["min_plateau"]) > 0:
        noise |= np.isfinite(plateau) & (plateau < float(ps["min_plateau"]))
    noise |= np.isfinite(under) & (under > float(ps["max_undershoot"]))
    flags[~short & noise] |= mc.FLAG_NOISE
    flags[sat] |= mc.FLAG_SATURATED
    flags[~short & pileup] |= mc.FLAG_PILEUP
    flags[np.abs(slope) > float(ps["max_baseline_slope"])] |= mc.FLAG_SLOPE
    flags[brms > float(ps["max_baseline_rms"])] |= mc.FLAG_NOISY

    return {
        "edge": edge, "baseline": base, "baseline_rms": brms, "baseline_slope": slope,
        "height": height, "charge": charge, "rise_ns": rise, "plateau": plateau,
        "undershoot": under, "flags": flags, "y": y, "y_time": y_time,
        "noise_used": y_time[:, :nb].std(axis=1),   # baseline noise of the waveform the time comes from
    }


def samples_at(y, x0, k):
    """y[i, floor(x0[i]) + k] for every row i; NaN where x0 is NaN or the sample lies outside the record."""
    n, L = y.shape
    j = np.floor(np.nan_to_num(x0, nan=-1e6)).astype(np.int64)[:, None] + np.asarray(k)[None, :]
    inside = (j >= 0) & (j < L)
    return np.where(inside, y[np.arange(n)[:, None], np.clip(j, 0, L - 1)], np.nan)


def template_samples(ps):
    """Samples around a pulse's crossing kept for the template: enough for the fit window and its shifts."""
    w0, w1 = int(ps["fit_window"][0]), int(ps["fit_window"][1])
    m = math.ceil(float(ps["fit_max_shift"]))
    return np.arange(min(-15, w0 - m - 3), max(25, w1 + m + 3) + 1)


def build_template(win, fx, k, frac):
    """Mean rising edge of clean pulses. win[i, j] = y / height of pulse i at sample floor(x0) + k[j], where x0
    is its crossing of frac and fx = x0 - floor(x0). The values are binned (TEMPLATE_STEP samples) against the
    distance from x0 and smoothed (Gaussian, TEMPLATE_SMOOTH samples: with few pulses the bins are noisy, and
    the fit would follow their wiggles), then the axis is shifted so that the mean crosses frac exactly at 0
    (the last crossing before 90 % of its top). Returns a dict (x in samples from that point, y, spread = rms
    of the pulses in each bin, start = where the steepest tangent meets 0, t10 / t80 / t90 = crossings) or None."""
    from scipy.ndimage import gaussian_filter1d
    rel = k[None, :] - fx[:, None]
    lo, hi = float(k[0] + 1), float(k[-1] - 1)
    ok = np.isfinite(win) & (rel >= lo) & (rel <= hi)
    edges = np.arange(lo, hi + TEMPLATE_STEP / 2, TEMPLATE_STEP)
    cnt, _ = np.histogram(rel[ok], bins=edges)
    s1, _ = np.histogram(rel[ok], bins=edges, weights=win[ok].astype(np.float64))
    s2, _ = np.histogram(rel[ok], bins=edges, weights=win[ok].astype(np.float64) ** 2)
    x = 0.5 * (edges[1:] + edges[:-1])
    filled = cnt > 0
    if filled.sum() < 0.9 * len(x):
        return None
    y = np.interp(x, x[filled], s1[filled] / cnt[filled])
    spread = np.sqrt(np.maximum(np.interp(x, x[filled], s2[filled] / cnt[filled]) - y ** 2, 0.0))
    y = gaussian_filter1d(y, TEMPLATE_SMOOTH / TEMPLATE_STEP, mode="nearest")
    top = float(y.max())
    if not top > 1.5 * frac:
        return None
    i_top = int(np.argmax(y > 0.9 * top))
    i_max = int(np.argmax(y))
    below = np.where(y[:i_top] <= frac)[0]
    if not len(below):
        return None
    j = below[-1]
    x = x - float(x[j] + (frac - y[j]) / (y[j + 1] - y[j]) * TEMPLATE_STEP)   # 0 = the crossing of frac

    def crossing(level):
        under = np.where(y[:i_max + 1] <= level)[0]
        if not len(under) or under[-1] >= i_max:
            return float("nan")
        i = under[-1]
        return float(x[i] + (level - y[i]) / (y[i + 1] - y[i]) * TEMPLATE_STEP)

    slope = np.convolve(np.gradient(y, TEMPLATE_STEP), np.ones(5) / 5, mode="same")
    i_s = int(np.argmax(slope[:i_max + 1]))
    return {"x": x, "y": y, "spread": spread, "n": int(len(win)), "top": top,
            "start": float(x[i_s] - y[i_s] / slope[i_s]) if slope[i_s] > 0 else float("nan"),
            "t10": crossing(0.1), "t80": crossing(0.8), "t90": crossing(0.9)}


def fit_edges(Y, fx, k, tx, ty, max_shift, progress=None):
    """Fit every row of Y (samples at floor(x0) + k, fx = x0 - floor(x0)) with A x T(sample - x0 - dx) by least
    squares, T = the template (tx, ty). For every trial dx on a grid of FIT_STEP samples within +-max_shift the
    height A follows in closed form; a parabola through the three best chi2 refines dx. Returns dx, A, the rms
    residual (ADC) and ok (all samples there, the best dx inside the grid, A > 0)."""
    n = len(Y)
    grid = np.arange(-max_shift, max_shift + FIT_STEP / 2, FIT_STEP)
    out_dx, out_a, out_rms = np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan)
    out_ok = np.zeros(n, dtype=bool)
    right = float(ty[-1])
    for s in range(0, n, FIT_BLOCK):
        y = Y[s:s + FIT_BLOCK].astype(np.float64)
        f = fx[s:s + FIT_BLOCK]
        good = np.isfinite(y).all(axis=1) & np.isfinite(f)
        y = np.where(good[:, None], y, 0.0)
        rel = np.asarray(k)[None, :] - np.nan_to_num(f)[:, None]
        chi = np.empty((len(y), len(grid)))
        for g, d in enumerate(grid):
            T = np.interp(rel - d, tx, ty, left=0.0, right=right)
            a = (y * T).sum(axis=1) / np.maximum((T * T).sum(axis=1), 1e-12)
            chi[:, g] = ((y - a[:, None] * T) ** 2).sum(axis=1)
        best = np.argmin(chi, axis=1)
        r = np.arange(len(y))
        b = np.clip(best, 1, len(grid) - 2)
        c0, c1, c2 = chi[r, b - 1], chi[r, b], chi[r, b + 1]
        den = c0 - 2 * c1 + c2
        with np.errstate(invalid="ignore", divide="ignore"):
            frac_step = np.where(den > 0, np.clip(0.5 * (c0 - c2) / den, -1.0, 1.0), 0.0)
        d = grid[b] + frac_step * FIT_STEP
        T = np.interp(rel - d[:, None], tx, ty, left=0.0, right=right)
        a = (y * T).sum(axis=1) / np.maximum((T * T).sum(axis=1), 1e-12)
        out_dx[s:s + len(y)] = d
        out_a[s:s + len(y)] = a
        out_rms[s:s + len(y)] = np.sqrt(((y - a[:, None] * T) ** 2).mean(axis=1))
        out_ok[s:s + len(y)] = good & (best > 0) & (best < len(grid) - 1) & (a > 0)
        if progress is not None:
            progress(len(y))
    return out_dx, out_a, out_rms, out_ok


def process_channel(tree, chid, ref, sign, ps, max_waveforms, bar):
    """Measure every waveform of one channel and, with timing = "fit", fit its rising edges. Returns the columns,
    the average pulse shapes, the edge template (or None) and the fit's numbers for the summary."""
    n = tree.GetEntries() if max_waveforms is None else min(tree.GetEntries(), max_waveforms)
    L = window_length(ref, ps)
    tree.SetBranchStatus("*", 0)
    for name in ("ttt", "tbusy", "nsamps", "chdata"):
        tree.SetBranchStatus(name, 1)
    cols = {name: [] for name in ("entry", "ttt", "tbusy", "nsamps", "edge_cfd", "baseline", "baseline_rms",
                                  "baseline_slope", "height", "charge", "rise_ns", "plateau",
                                  "undershoot", "flags")}
    shapes = []
    fit = ps["timing"] == "fit"
    k_fit = np.arange(int(ps["fit_window"][0]), int(ps["fit_window"][1]) + 1)
    k_tmpl = template_samples(ps)
    n_tmpl, min_h = int(ps["template_pulses"]), float(ps["template_min_height"])
    fit_y, tmpl_y, tmpl_fx, noise = [], [], [], []
    for start in range(0, n, CHUNK):
        stop = min(start + CHUNK, n)
        W = np.full((stop - start, L), np.nan)
        ttt = np.zeros(stop - start, dtype=np.int64)
        tbusy = np.zeros(stop - start, dtype=np.int64)
        ns = np.zeros(stop - start, dtype=np.int64)
        for i in range(start, stop):
            tree.GetEntry(i)
            m = min(int(tree.nsamps), L)
            if m:
                W[i - start, :m] = np.frombuffer(tree.chdata, dtype=np.uint16, count=m)
            ttt[i - start] = tree.ttt
            tbusy[i - start] = tree.tbusy
            ns[i - start] = tree.nsamps
        res = measure(W, ns, ref, sign, ps)
        y, y_time = res.pop("y"), res.pop("y_time")
        noise.append(res.pop("noise_used"))
        res["edge_cfd"] = res.pop("edge")
        res.update({"entry": np.arange(start, stop, dtype=np.int64), "ttt": ttt, "tbusy": tbusy,
                    "nsamps": ns})
        for key in cols:
            cols[key].append(res[key])
        x0 = res["edge_cfd"]
        if len(shapes) < SHAPE_KEEP:
            good = np.where((res["flags"] == 0) & np.isfinite(x0))[0]
            for i in good[:SHAPE_KEEP - len(shapes)]:
                j = int(math.floor(x0[i]))
                if j + SHAPE_FROM >= 0 and j + SHAPE_TO <= L:
                    seg = y[i, j + SHAPE_FROM:j + SHAPE_TO] / res["height"][i]
                    if np.isfinite(seg).all():
                        shapes.append(seg)
        if fit:
            fit_y.append(samples_at(y_time, x0, k_fit).astype(np.float32))
            have = sum(len(a) for a in tmpl_fx)
            if have < n_tmpl:
                good = np.where((res["flags"] == 0) & np.isfinite(x0) & (res["height"] >= min_h))[0][:n_tmpl - have]
                if len(good):
                    tmpl_y.append((samples_at(y_time[good], x0[good], k_tmpl)
                                   / res["height"][good][:, None]).astype(np.float32))
                    tmpl_fx.append(x0[good] - np.floor(x0[good]))
        bar.update(stop - start, f"ch{chid}")
    tree.SetBranchStatus("*", 1)
    out = {k: np.concatenate(v) if v else np.zeros(0) for k, v in cols.items()}
    x0 = out["edge_cfd"]
    m = len(x0)
    out["edge"], out["fit_amp"], out["fit_rms"] = x0.copy(), np.full(m, np.nan), np.full(m, np.nan)
    noise = np.concatenate(noise) if noise else np.zeros(0)
    if not fit:
        return out, np.array(shapes), None, None
    has = np.isfinite(x0)
    n_t = sum(len(a) for a in tmpl_fx)
    stats = {"template_pulses": n_t, "with_crossing": int(has.sum()),
             "median_noise_adc": float(np.median(noise)) if len(noise) else float("nan")}
    tmpl = None
    if n_t >= int(ps["template_min_pulses"]):
        tmpl = build_template(np.concatenate(tmpl_y), np.concatenate(tmpl_fx), k_tmpl,
                              float(ps["cfd_fraction"]))
        if tmpl is None:
            stats["status"] = "no template: the mean of the clean pulses has no clear rising edge"
    else:
        stats["status"] = (f"no template: {n_t} clean pulses of >= {min_h:g} ADC "
                           f"(template_min_pulses {int(ps['template_min_pulses'])})")
    if tmpl is None:
        out["flags"][has] |= mc.FLAG_NO_FIT
        stats["fitted"] = 0
        bar.update(m, f"ch{chid}: no template")
        return out, np.array(shapes), None, stats
    dx, amp, rms, ok = fit_edges(np.concatenate(fit_y), x0 - np.floor(np.nan_to_num(x0)), k_fit,
                                 tmpl["x"], tmpl["y"], float(ps["fit_max_shift"]),
                                 progress=lambda k: bar.update(k, f"ch{chid} edge fit"))
    ok &= has
    out["edge"] = np.where(ok, x0 + dx, x0)
    out["fit_amp"] = np.where(ok, amp, np.nan)
    out["fit_rms"] = np.where(ok, rms, np.nan)
    out["flags"][has & ~ok] |= mc.FLAG_NO_FIT
    d = dx[ok] * mc.NS_PER_SAMPLE
    stats.update({
        "status": "fitted", "fitted": int(ok.sum()),
        "fit_minus_crossing_median_ns": float(np.median(d)) if len(d) else float("nan"),
        "fit_minus_crossing_rms_ns": float(1.4826 * np.median(np.abs(d - np.median(d)))) if len(d) else float("nan"),
        "median_fit_rms_adc": float(np.median(rms[ok])) if ok.any() else float("nan"),
        "template_start_ns": tmpl["start"] * mc.NS_PER_SAMPLE, "template_rise_ns":
            (tmpl["t90"] - tmpl["t10"]) * mc.NS_PER_SAMPLE, "template_t80_ns": tmpl["t80"] * mc.NS_PER_SAMPLE,
    })
    return out, np.array(shapes), tmpl, stats


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file", nargs="?", default=None, help="ROOT file of the run (asked for if left out)")
    ap.add_argument("--out", default=None, help="Analysis folder (default: <run name>_muon next to the ROOT file)")
    ap.add_argument("--setup", default=None,
                    help="Setup file (default: the one in setups/ valid for the run's date)")
    ap.add_argument("--pulse-pre", type=int, default=None,
                    help="pulse_pre of the run (default: from <run>.conf next to the ROOT file)")
    ap.add_argument("--polarity", choices=("positive", "negative"), default=None,
                    help="Pulse polarity of every channel (default: the setup's [pulse] polarity, "
                         "normally pulsePolarity of the .conf)")
    ap.add_argument("--max-waveforms", type=int, default=None,
                    help="Measure only the first N waveforms of every channel (default: all)")
    ap.add_argument("--lowpass", type=float, default=None,
                    help="Low-pass filter at this frequency (MHz, -3 dB) before measuring; 0 = off "
                         "(default: the setup's [pulse] lowpass_mhz)")
    ap.add_argument("--lowpass-kind", choices=("gauss", "butter"), default=None,
                    help="Filter shape (default: the setup's lowpass_kind)")
    ap.add_argument("--raw-timing", action="store_true",
                    help="Take the pulse time from the raw samples even with the filter on")
    ap.add_argument("--timing", choices=("fit", "cfd"), default=None,
                    help="Pulse time: fit = template fit of the rising edge, cfd = crossing between two samples "
                         "(default: the setup's [pulse] timing, normally fit)")
    args = ap.parse_args()
    asked = args.file is None
    path = mc.ask_data_file(args.file)
    out = mc.ask_out_folder(args.out, mc.default_out(path)) if asked else \
        mc.ask_out_folder(args.out or mc.default_out(path), None)
    run(path, out, args.setup, args.pulse_pre, args.polarity, args.max_waveforms,
        args.lowpass, args.lowpass_kind, args.raw_timing, args.timing)


def run(path, out, setup_path=None, pulse_pre=None, polarity=None, max_waveforms=None,
        lowpass_mhz=None, lowpass_kind=None, raw_timing=False, timing=None):
    """Step 1 for one ROOT file; returns the path of pulses.root."""
    setup = mc.choose_setup(setup_path, path)
    ps = dict(setup.pulse)
    if lowpass_mhz is not None:
        ps["lowpass_mhz"] = float(lowpass_mhz)
    if lowpass_kind is not None:
        ps["lowpass_kind"] = lowpass_kind
    if raw_timing:
        ps["lowpass_timing"] = False
    if timing is not None:
        ps["timing"] = timing
    fit = ps["timing"] == "fit"
    conf = mc.read_conf(path)
    if pulse_pre is None:
        if conf is None or conf.pulse_pre is None:
            sys.exit(f"ERROR: no <run>.conf with pulse_pre next to {path}; give it with --pulse-pre")
        pulse_pre = conf.pulse_pre
    ref = 2 * pulse_pre
    if ref + int(ps["baseline_end"]) < 40:
        sys.exit(f"ERROR: pulse_pre {pulse_pre} leaves fewer than 40 samples before the trigger for the "
                 f"baseline; take runs with pulse_pre >= 30")
    root_file = mc.open_root(path)
    if root_file is None:
        sys.exit(f"ERROR: cannot open {path}")
    trees = {}
    for board, chans in mc.channel_trees(root_file).items():
        for chn, name in chans.items():
            trees[int(name[2:])] = name
    if not trees:
        sys.exit(f"ERROR: no ch<id> trees in {path}")
    os.makedirs(out, exist_ok=True)

    lines = []

    def report(text=""):
        print(text)
        lines.append(text)

    report(f"File        : {path}")
    report(f"Setup       : {setup.name} ({setup.path})")
    if conf is not None:
        report(f"Run conf    : {os.path.basename(conf.path)}: trigmode {conf.trigmode}, ctrlmode {conf.ctrlmode}, "
               f"pulse_pre {conf.pulse_pre_typed}, pulse_post {conf.pulse_post}, daqruntime {conf.daqruntime} s")
        if conf.note:
            report(f"Note        : {conf.note}")
        if conf.trigmode not in (0, None):
            report(f"WARNING     : trigmode {conf.trigmode}: this analysis was checked on self-trigger "
                   f"runs (trigmode 0) only; the trigger is assumed at sample 2 x pulse_pre")
    report(f"Trigger     : ttt refers to sample {ref} of each record (2 x pulse_pre {pulse_pre})")
    if float(ps["lowpass_mhz"]) > 0:
        report(f"Filter      : {ps['lowpass_kind']} low-pass at {float(ps['lowpass_mhz']):g} MHz, zero phase; pulse time "
               f"from the {'filtered' if ps['lowpass_timing'] else 'raw'} samples; baseline noise from the raw samples")
    else:
        report("Filter      : none (raw samples)")
    frac = float(ps["cfd_fraction"])
    w0, w1 = int(ps["fit_window"][0]), int(ps["fit_window"][1])
    if fit:
        report(f"Pulse time  : rising-edge fit (timing = \"fit\"): each channel's template = the mean edge of its first "
               f"{int(ps['template_pulses'])} clean pulses")
        report(f"              of >= {float(ps['template_min_height']):g} ADC; every pulse fitted with height x template "
               f"over samples {w0:+d} .. {w1:+d} around its {frac:.0%} crossing;")
        report(f"              the time is where the fitted edge reaches {frac:.0%} of the fitted height")
    else:
        report(f"Pulse time  : {frac:.0%} crossing between two samples (timing = \"cfd\")")
    missing = sorted(set(setup.channels) - set(trees))
    unmapped = sorted(set(trees) - set(setup.channels))
    if missing:
        report(f"WARNING     : in the setup but not in the file: " +
               ", ".join(f"ch{c} ({setup.channel_label(c)})" for c in missing))
    if unmapped:
        report(f"Note        : in the file but not in the setup (measured, not used later): " +
               ", ".join(f"ch{c}" for c in unmapped))
    boards = sorted({mc.channel_board(c) for c in trees})
    if len(boards) > 1 and conf is not None and conf.ctrlmode != 3:
        report(f"WARNING     : {len(boards)} boards but ctrlmode {conf.ctrlmode}: their clocks did not start "
               f"together, so times cannot be compared across boards")

    total = sum(root_file.Get(name).GetEntries() if max_waveforms is None else
                min(root_file.Get(name).GetEntries(), max_waveforms) for name in trees.values())
    bar = mc.Progress("Measuring pulses" + (" and fitting edges" if fit else ""), total * (2 if fit else 1))
    columns, shapes, per_channel, templates = [], {}, {}, {}
    for chid in sorted(trees):
        tree = root_file.Get(trees[chid])
        if polarity is not None:
            positive = polarity == "positive"
        elif ps["polarity"] in ("positive", "negative"):
            positive = ps["polarity"] == "positive"
        else:
            positive = conf.positive(chid) if conf is not None else None
            if positive is None:
                positive = True
                report(f"WARNING     : ch{chid}: polarity unknown (no .conf); positive assumed")
        res, shape, tmpl, fit_stats = process_channel(tree, chid, ref, 1.0 if positive else -1.0, ps,
                                                      max_waveforms, bar)
        res["channel"] = np.full(len(res["entry"]), chid, dtype=np.int32)
        columns.append(res)
        shapes[chid] = shape
        if tmpl is not None:
            templates[chid] = tmpl
        per_channel[chid] = {"positive": positive, "threshold": conf.threshold(chid) if conf else None,
                             "records_in_file": int(tree.GetEntries()), "timing": fit_stats}
    bar.close()
    root_file.Close()

    data = {k: np.concatenate([c[k] for c in columns]) for k in columns[0]}
    origin = int(data["ttt"].min()) if len(data["ttt"]) else 0
    offsets = np.array([setup.time_offset(int(c)) for c in data["channel"]])
    edge_rel = np.where(np.isfinite(data["edge"]), data["edge"] - ref, 0.0)
    data["t_ns"] = (data["ttt"] - origin).astype(np.float64) * mc.NS_PER_TICK + edge_rel * mc.NS_PER_SAMPLE - offsets
    data["board"] = ((data["channel"] - 10000) // 10).astype(np.int32)
    span_ns = float((data["ttt"].max() - origin) * mc.NS_PER_TICK) if len(data["ttt"]) else 0.0

    out_cols = {
        "channel": data["channel"].astype(np.int32), "board": data["board"],
        "entry": data["entry"].astype(np.int64), "ttt": data["ttt"].astype(np.int64),
        "tbusy": data["tbusy"].astype(np.int64), "nsamps": data["nsamps"].astype(np.int32),
        "t_ns": data["t_ns"].astype(np.float64), "edge": data["edge"].astype(np.float64),
        "edge_cfd": data["edge_cfd"].astype(np.float64), "fit_amp": data["fit_amp"].astype(np.float32),
        "fit_rms": data["fit_rms"].astype(np.float32),
        "baseline": data["baseline"].astype(np.float32), "baseline_rms": data["baseline_rms"].astype(np.float32),
        "baseline_slope": data["baseline_slope"].astype(np.float32), "height": data["height"].astype(np.float32),
        "charge": data["charge"].astype(np.float32), "rise_ns": data["rise_ns"].astype(np.float32),
        "plateau": data["plateau"].astype(np.float32), "undershoot": data["undershoot"].astype(np.float32),
        "flags": data["flags"].astype(np.int32),
    }
    pulses_path = os.path.join(out, mc.PULSES_FILE)
    mc.write_tree(pulses_path, "pulses", out_cols)

    t_axis = (np.arange(SHAPE_FROM, SHAPE_TO)) * mc.NS_PER_SAMPLE
    sh = {"channel": [], "t_ns": [], "mean": [], "p16": [], "p84": [], "n": []}
    for chid, s in shapes.items():
        if len(s) == 0:
            continue
        sh["channel"].append(np.full(len(t_axis), chid, dtype=np.int32))
        sh["t_ns"].append(t_axis)
        sh["mean"].append(s.mean(axis=0))
        sh["p16"].append(np.percentile(s, 16, axis=0))
        sh["p84"].append(np.percentile(s, 84, axis=0))
        sh["n"].append(np.full(len(t_axis), len(s), dtype=np.int32))
    sh = {k: (np.concatenate(v) if v else np.zeros(0, dtype=np.int32 if k in ("channel", "n") else np.float64))
          for k, v in sh.items()}
    mc.write_tree(pulses_path, "shapes", sh, update=True)
    tm = {"channel": [], "x_ns": [], "y": [], "spread": [], "n": []}
    for chid, t in sorted(templates.items()):
        tm["channel"].append(np.full(len(t["x"]), chid, dtype=np.int32))
        tm["x_ns"].append(t["x"] * mc.NS_PER_SAMPLE)
        tm["y"].append(t["y"])
        tm["spread"].append(t["spread"])
        tm["n"].append(np.full(len(t["x"]), t["n"], dtype=np.int32))
    tm = {k: (np.concatenate(v) if v else np.zeros(0, dtype=np.int32 if k in ("channel", "n") else np.float64))
          for k, v in tm.items()}
    mc.write_tree(pulses_path, "templates", tm, update=True)

    # ---- summary table ----
    report()
    report(f"Run time    : {span_ns * 1e-9:.3f} s from the first to the last waveform "
           f"(daqruntime {conf.daqruntime if conf else '?'} s)")
    head = (f"{'channel':<9}{'paddle':<10}{'thr':>5}{'records':>9}{'used':>7}{'rate Hz':>9}{'noise':>7}"
            f"{'low':>5}{'short':>6}{'edge?':>6}{'sat':>5}{'pileup':>7}{'slope':>6}{'noisy':>6}{'nofit':>6}"
            f"{'height':>8}{'rms':>6}{'dead %':>7}")
    report(head)
    report("-" * len(head))
    summary = {}
    for chid in sorted(trees):
        sel = data["channel"] == chid
        fl = data["flags"][sel]
        used = (fl & mc.REJECT_MASK) == 0
        cnt = lambda bit: int(((fl & bit) != 0).sum())
        rate = used.sum() / (span_ns * 1e-9) if span_ns > 0 else float("nan")
        h = data["height"][sel][used]
        rms = data["baseline_rms"][sel]
        dead = data["tbusy"][sel].max() * mc.NS_PER_BUSY_TICK / span_ns * 100 if span_ns > 0 and sel.any() else 0.0
        thr = per_channel[chid]["threshold"]
        summary[str(chid)] = {
            "label": setup.channel_label(chid), "records": int(sel.sum()), "used": int(used.sum()),
            "rate_hz": rate, "threshold": thr, "positive": per_channel[chid]["positive"],
            "flags": {name: cnt(bit) for bit, name, _ in mc.PULSE_FLAGS},
            "median_height": float(np.median(h)) if len(h) else None,
            "median_baseline_rms": float(np.median(rms)) if len(rms) else None, "dead_percent": dead,
            "timing": per_channel[chid]["timing"],
        }
        report(f"ch{chid:<7}{setup.channel_label(chid):<10}{thr if thr is not None else '-':>5}{sel.sum():>9}"
               f"{used.sum():>7}{rate:>9.2f}{cnt(mc.FLAG_NOISE):>7}{cnt(mc.FLAG_LOW):>5}{cnt(mc.FLAG_SHORT):>6}"
               f"{cnt(mc.FLAG_NO_EDGE):>6}{cnt(mc.FLAG_SATURATED):>5}{cnt(mc.FLAG_PILEUP):>7}"
               f"{cnt(mc.FLAG_SLOPE):>6}{cnt(mc.FLAG_NOISY):>6}{cnt(mc.FLAG_NO_FIT):>6}"
               f"{(np.median(h) if len(h) else float('nan')):>8.0f}{(np.median(rms) if len(rms) else float('nan')):>6.1f}"
               f"{dead:>7.2f}")
    report()
    report("Columns: thr = DAQ threshold (ADC above the baseline); used = pulses the later steps use (not SHORT,")
    report("NO_EDGE, NOISE or LOW); rate = used per second; noise ... nofit = waveforms with that flag; height =")
    report("median height of the used pulses (ADC); rms = median baseline noise (ADC); dead = tbusy x 4 ns / run time.")
    report("Flags:")
    for bit, name, text in mc.PULSE_FLAGS:
        report(f"  {bit:>3} {name:<10}{text}")
    if fit:
        report()
        report(f"Rising-edge fit: samples {w0:+d} .. {w1:+d} around each crossing ({(w0 - 1) * mc.NS_PER_SAMPLE:+.0f} .. "
               f"{w1 * mc.NS_PER_SAMPLE:+.0f} ns from it); times in ns from the template's {frac:.0%} point")
        head = (f"{'channel':<9}{'paddle':<10}{'template':>9}{'rise':>7}{'start':>7}{'80 %':>7}{'fitted %':>10}"
                f"{'fit-cross':>11}{'spread':>8}{'fit rms':>9}{'noise':>7}")
        report(head)
        report("-" * len(head))
        for chid in sorted(trees):
            st = per_channel[chid]["timing"] or {}
            label = setup.channel_label(chid)
            if st.get("status") != "fitted":
                report(f"ch{chid:<7}{label:<10}{st.get('template_pulses', 0):>9}  {st.get('status', '-')}")
                continue
            pct = 100.0 * st["fitted"] / max(st["with_crossing"], 1)
            report(f"ch{chid:<7}{label:<10}{st['template_pulses']:>9}{st['template_rise_ns']:>7.1f}"
                   f"{st['template_start_ns']:>+7.1f}{st['template_t80_ns']:>+7.1f}{pct:>10.2f}"
                   f"{st['fit_minus_crossing_median_ns']:>+11.3f}{st['fit_minus_crossing_rms_ns']:>8.3f}"
                   f"{st['median_fit_rms_adc']:>9.1f}{st['median_noise_adc']:>7.1f}")
        report()
        report(f"template = clean pulses averaged; rise = 10-90 % of the template (ns); start = where its steepest")
        report(f"tangent meets the baseline, 80 % = where it reaches 80 % of the height (ns from its {frac:.0%} point);")
        report("fitted = pulses with a crossing whose fit worked (the others keep the crossing, flag NO_FIT);")
        report("fit-cross = median of (fit time - crossing), spread = its robust rms (ns); fit rms = median rms")
        report("residual of the fits, noise = median baseline noise of the waveform the fit uses (ADC): a fit rms")
        report("well above the noise means that the edges differ from pulse to pulse more than the noise does.")
        for chid in sorted(trees):
            st = per_channel[chid]["timing"] or {}
            if st.get("status") != "fitted":
                continue
            if w1 < st["template_t80_ns"] / mc.NS_PER_SAMPLE:
                report(f"WARNING     : ch{chid}: fit_window ends ({w1:+d}) before the edge reaches 80 % "
                       f"({st['template_t80_ns'] / mc.NS_PER_SAMPLE:+.1f} samples): widen fit_window")
            if np.isfinite(st["template_start_ns"]) and w0 > st["template_start_ns"] / mc.NS_PER_SAMPLE - 1:
                report(f"WARNING     : ch{chid}: fit_window starts ({w0:+d}) less than a sample before the edge "
                       f"({st['template_start_ns'] / mc.NS_PER_SAMPLE:+.1f} samples): lower fit_window[0]")

    info = mc.stamp("extract_pulses")
    info.update({
        "raw_file": path, "run": mc.run_name(path), "trigger_sample": ref, "pulse_pre": pulse_pre,
        "ttt_origin": origin, "span_ns": span_ns, "max_waveforms": max_waveforms,
        "setup": {"name": setup.name, "path": setup.path}, "pulse_settings": ps,
        "board_offsets_ns": setup.board_offsets, "channel_offsets_ns": setup.channel_offsets,
        "conf": conf.summary() if conf else None, "channels": summary,
        "missing_channels": missing, "unmapped_channels": unmapped,
    })
    mc.write_info(pulses_path, info)
    with open(os.path.join(out, "pulses_summary.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Saved {pulses_path} and pulses_summary.txt")
    return pulses_path


if __name__ == "__main__":
    main()
