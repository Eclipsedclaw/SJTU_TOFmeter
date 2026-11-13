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
    time         where the rising edge crosses cfd_fraction (default 0.2) of the
                 height, interpolated linearly between two samples ("constant
                 fraction"): the same point of every pulse, whatever its height.
                 t_ns = (ttt - first ttt of the run) x 2 ns + (crossing - s) x 2 ns
                 - the board and channel offsets of the setup file
    rise_ns      10 % to 90 % of the height on the rising edge
    charge       sum of the samples from the edge + charge_window, x 2 ns (ADC x ns)
    tests        flat top: a real pulse of this front end stays high after its
                 peak (mean over plateau_window / height >= min_plateau); noise:
                 no swing below the baseline deeper than max_undershoot x height;
                 pile-up: no later rise of 4 samples above pileup_fraction x height
                 and pileup_min ADC; saturation: no sample at code 0 or 16383
    flags        bit mask, see pulses_summary.txt; SHORT, NO_EDGE, NOISE and LOW
                 pulses are left out by the later steps

Outputs (in the analysis folder):
    pulses.root          tree 'pulses' (one entry per waveform), tree 'shapes'
                         (average pulse shape per channel), and 'info' (settings)
    pulses_summary.txt   one line per channel: waveforms, flags, rate, noise
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


def window_length(ref, ps):
    """Samples of each record that the measurement uses."""
    s1 = int(ps["search"][1])
    k = int(ps["smooth"])
    return ref + s1 + k + max(int(ps["plateau_window"][1]), int(ps["charge_window"][1]),
                              SHAPE_TO) + 5


def measure(W, nsamps, ref, sign, ps):
    """Measure a chunk of waveforms. W: (n, L) float64, NaN past the end of a record."""
    n, L = W.shape
    r = np.arange(n)
    nb = ref + int(ps["baseline_end"])
    k = int(ps["smooth"])
    lo, hi = ref + int(ps["search"][0]), ref + int(ps["search"][1])
    flags = np.zeros(n, dtype=np.int32)
    short = nsamps < hi + k
    flags[short] |= mc.FLAG_SHORT

    base = W[:, :nb].mean(axis=1)
    brms = W[:, :nb].std(axis=1)
    slope = W[:, nb - 30:nb].mean(axis=1) - W[:, :30].mean(axis=1)
    sat = ((W >= mc.ADC_MAX) | (W <= 0)).any(axis=1)

    y = sign * (W - base[:, None])
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
        part = y[:, a:b]
        below = (part <= level[:, None]) & (cols[None, :] < ic[:, None])
        has = below.any(axis=1)
        j = a + (below.shape[1] - 1 - np.argmax(below[:, ::-1], axis=1))
        ya = y[r, j]
        yb = y[r, np.minimum(j + 1, L - 1)]
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
        "undershoot": under, "flags": flags, "y": y,
    }


def process_channel(tree, chid, ref, sign, ps, max_waveforms, bar):
    n = tree.GetEntries() if max_waveforms is None else min(tree.GetEntries(), max_waveforms)
    L = window_length(ref, ps)
    tree.SetBranchStatus("*", 0)
    for name in ("ttt", "tbusy", "nsamps", "chdata"):
        tree.SetBranchStatus(name, 1)
    cols = {name: [] for name in ("entry", "ttt", "tbusy", "nsamps", "edge", "baseline", "baseline_rms",
                                  "baseline_slope", "height", "charge", "rise_ns", "plateau",
                                  "undershoot", "flags")}
    shapes = []
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
        y = res.pop("y")
        res.update({"entry": np.arange(start, stop, dtype=np.int64), "ttt": ttt, "tbusy": tbusy,
                    "nsamps": ns})
        for key in cols:
            cols[key].append(res[key])
        if len(shapes) < SHAPE_KEEP:
            good = np.where((res["flags"] == 0) & np.isfinite(res["edge"]))[0]
            for i in good[:SHAPE_KEEP - len(shapes)]:
                j = int(math.floor(res["edge"][i]))
                if j + SHAPE_FROM >= 0 and j + SHAPE_TO <= L:
                    seg = y[i, j + SHAPE_FROM:j + SHAPE_TO] / res["height"][i]
                    if np.isfinite(seg).all():
                        shapes.append(seg)
        bar.update(stop - start, f"ch{chid}")
    tree.SetBranchStatus("*", 1)
    out = {k: np.concatenate(v) if v else np.zeros(0) for k, v in cols.items()}
    return out, np.array(shapes)


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
    args = ap.parse_args()
    asked = args.file is None
    path = mc.ask_data_file(args.file)
    out = mc.ask_out_folder(args.out, mc.default_out(path)) if asked else \
        mc.ask_out_folder(args.out or mc.default_out(path), None)
    run(path, out, args.setup, args.pulse_pre, args.polarity, args.max_waveforms)


def run(path, out, setup_path=None, pulse_pre=None, polarity=None, max_waveforms=None):
    """Step 1 for one ROOT file; returns the path of pulses.root."""
    setup = mc.choose_setup(setup_path, path)
    ps = setup.pulse
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
    bar = mc.Progress("Measuring pulses", total)
    columns, shapes, per_channel = [], {}, {}
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
        res, shape = process_channel(tree, chid, ref, 1.0 if positive else -1.0, ps, max_waveforms, bar)
        res["channel"] = np.full(len(res["entry"]), chid, dtype=np.int32)
        columns.append(res)
        shapes[chid] = shape
        per_channel[chid] = {"positive": positive, "threshold": conf.threshold(chid) if conf else None,
                             "records_in_file": int(tree.GetEntries())}
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

    # ---- summary table ----
    report()
    report(f"Run time    : {span_ns * 1e-9:.3f} s from the first to the last waveform "
           f"(daqruntime {conf.daqruntime if conf else '?'} s)")
    head = (f"{'channel':<9}{'paddle':<10}{'thr':>5}{'records':>9}{'used':>7}{'rate Hz':>9}{'noise':>7}"
            f"{'low':>5}{'short':>6}{'edge?':>6}{'sat':>5}{'pileup':>7}{'slope':>6}{'noisy':>6}"
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
        }
        report(f"ch{chid:<7}{setup.channel_label(chid):<10}{thr if thr is not None else '-':>5}{sel.sum():>9}"
               f"{used.sum():>7}{rate:>9.2f}{cnt(mc.FLAG_NOISE):>7}{cnt(mc.FLAG_LOW):>5}{cnt(mc.FLAG_SHORT):>6}"
               f"{cnt(mc.FLAG_NO_EDGE):>6}{cnt(mc.FLAG_SATURATED):>5}{cnt(mc.FLAG_PILEUP):>7}"
               f"{cnt(mc.FLAG_SLOPE):>6}{cnt(mc.FLAG_NOISY):>6}"
               f"{(np.median(h) if len(h) else float('nan')):>8.0f}{(np.median(rms) if len(rms) else float('nan')):>6.1f}"
               f"{dead:>7.2f}")
    report()
    report("Columns: thr = DAQ threshold (ADC above the baseline); used = pulses the later steps use (not SHORT,")
    report("NO_EDGE, NOISE or LOW); rate = used per second; noise ... noisy = waveforms with that flag; height =")
    report("median height of the used pulses (ADC); rms = median baseline noise (ADC); dead = tbusy x 4 ns / run time.")
    report("Flags:")
    for bit, name, text in mc.PULSE_FLAGS:
        report(f"  {bit:>3} {name:<10}{text}")

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
