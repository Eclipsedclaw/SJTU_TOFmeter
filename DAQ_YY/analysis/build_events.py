#!/usr/bin/env python3
"""
build_events.py  -  step 2 of the paddle muon analysis

Pairs the pulses of the two ends of every paddle into hits, then groups the hits
of all paddles that come within a short time of each other into events. Reads
pulses.root (step 1) and writes events.root and events_summary.txt into the same
analysis folder.

Usage:
    python3 build_events.py                    # asks for the analysis folder
    python3 build_events.py /path/to/run00056_..._muon
    python3 build_events.py FOLDER --setup setups/two_layer_2026-11-12.toml

Settings: the [events] section of the setup file (by default the setup that
step 1 used).

Method (times in ns):
    pulses   only the pulses that step 1 did not reject (no SHORT, NO_EDGE, NOISE or LOW flag)
    hit      the A and B pulse of a paddle, paired one to one, closest pairs first,
             if |tA - tB| <= length / v_eff + pair_margin_ns (v_eff from the
             [calibration] section, default 15 cm/ns; margin default 5 ns). A muon
             at u cm from the paddle's centre (u > 0 towards the A end) gives
             tA - tB = -2 u / v_eff + offset, and (tA + tB) / 2 does not depend on u.
             A pulse without a partner is a one-end hit (kind A or B).
    hit time (tA + tB) / 2 for a two-end hit, the pulse time for a one-end hit
    event    hits sorted by time; an event starts with the first hit not yet used
             and takes every later hit within window_ns (default 30 ns) of it
    chance   for every pair of layers, the coincidences expected by chance:
             counted after shifting one layer's hits by 10, 20, ... 500 us (mean
             of 50 shifts), and computed as N1 x N2 x 2 x window / run time

Outputs (in the analysis folder):
    events.root          tree 'hits' (one entry per hit, with its event number),
                         tree 'events' (one entry per event), and 'info'
    events_summary.txt   hits per paddle, events by layers hit, coincidences
"""

import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse
import sys

import numpy as np

import muon_common as mc

SHIFTS_US = np.arange(10, 501, 10)    # chance coincidences: layer shifts of 10 .. 500 us


def pair_ends(tA, tB, max_dt):
    """One-to-one pairs (indices) of the times tA and tB closer than max_dt,
    closest first; and which of tA / tB stayed without a partner."""
    order_b = np.argsort(tB, kind="stable")
    tbs = tB[order_b]
    lo = np.searchsorted(tbs, tA - max_dt, "left")
    hi = np.searchsorted(tbs, tA + max_dt, "right")
    counts = hi - lo
    ia = np.repeat(np.arange(len(tA)), counts)
    first = np.repeat(np.cumsum(counts) - counts, counts)
    ib = order_b[np.repeat(lo, counts) + (np.arange(counts.sum()) - first)]
    order = np.argsort(np.abs(tA[ia] - tB[ib]), kind="stable")
    used_a = np.zeros(len(tA), dtype=bool)
    used_b = np.zeros(len(tB), dtype=bool)
    pairs = []
    for k in order:
        a, b = ia[k], ib[k]
        if not used_a[a] and not used_b[b]:
            used_a[a] = used_b[b] = True
            pairs.append((a, b))
    pairs = np.array(pairs, dtype=np.int64).reshape(-1, 2)
    return pairs, ~used_a, ~used_b


def count_near(t1, t2, window, shift=0.0):
    """How many times in t1 have a time of (sorted) t2 within window of t1 + shift."""
    lo = np.searchsorted(t2, t1 + shift - window, "left")
    hi = np.searchsorted(t2, t1 + shift + window, "right")
    return int((hi > lo).sum())


def group_events(t, window):
    """Event number of each time: an event starts at the first time more than
    `window` after the start of the previous event."""
    order = np.argsort(t, kind="stable")
    event = np.empty(len(t), dtype=np.int64)
    current, start = -1, -np.inf
    for i in order:
        if t[i] > start + window:
            current += 1
            start = t[i]
        event[i] = current
    return event


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", default=None,
                    help="Analysis folder with pulses.root (asked for if left out)")
    ap.add_argument("--setup", default=None, help="Setup file (default: the one step 1 used)")
    args = ap.parse_args()
    folder = mc.ask_analysis_folder(args.folder, mc.PULSES_FILE, "extract_pulses.py")
    run(folder, args.setup)


def run(folder, setup_path=None):
    """Step 2 for one analysis folder; returns the path of events.root."""
    pulses_path = mc.need_file(folder, mc.PULSES_FILE, "extract_pulses.py")
    pinfo = mc.read_info(pulses_path)
    setup = mc.setup_of(folder, setup_path)
    span_s = pinfo["span_ns"] * 1e-9
    P = mc.read_tree(pulses_path, "pulses", ["channel", "t_ns", "height", "charge", "flags"])
    good = (P["flags"] & mc.REJECT_MASK) == 0
    v_nom = float(setup.calibration["v_eff_cm_per_ns"])
    margin = float(setup.events["pair_margin_ns"])
    window = float(setup.events["window_ns"])

    lines = []

    def report(text=""):
        print(text)
        lines.append(text)

    report(f"Analysis folder : {folder}")
    report(f"Run             : {pinfo['run']} ({span_s:.3f} s)")
    report(f"Setup           : {setup.name}")
    report(f"Hits            : A and B pulse within length / {v_nom:g} cm/ns + {margin:g} ns")
    report(f"Events          : hits within {window:g} ns of the event's first hit")

    cols = {k: [] for k in ("paddle", "layer", "kind", "t_ns", "tA_ns", "tB_ns", "tdiff_ns", "hA", "hB",
                            "qA", "qB", "flagsA", "flagsB", "pulseA", "pulseB", "chA", "chB")}
    per_paddle = {}
    bar = mc.Progress("Pairing ends", len(setup.paddles))
    for p in setup.paddles:
        ends = {}
        for end, chid in p.channels.items():
            idx = np.where(good & (P["channel"] == chid))[0]
            ends[end] = idx[np.argsort(P["t_ns"][idx], kind="stable")]
        ia = ends.get("A", np.zeros(0, dtype=np.int64))
        ib = ends.get("B", np.zeros(0, dtype=np.int64))
        if p.both_ends:
            max_dt = p.length / v_nom + margin
            pairs, lone_a, lone_b = pair_ends(P["t_ns"][ia], P["t_ns"][ib], max_dt)
        else:
            pairs = np.zeros((0, 2), dtype=np.int64)
            lone_a, lone_b = np.ones(len(ia), dtype=bool), np.ones(len(ib), dtype=bool)
        groups = [(ia[pairs[:, 0]], ib[pairs[:, 1]], mc.HIT_AB),
                  (ia[lone_a], None, mc.HIT_A_ONLY), (None, ib[lone_b], mc.HIT_B_ONLY)]
        for a, b, kind in groups:
            n = len(a) if a is not None else len(b)
            if n == 0:
                continue
            nan = np.full(n, np.nan)
            none = np.full(n, -1, dtype=np.int64)
            ta = P["t_ns"][a] if a is not None else nan
            tb = P["t_ns"][b] if b is not None else nan
            cols["paddle"].append(np.full(n, p.index, dtype=np.int32))
            cols["layer"].append(np.full(n, p.layer.index, dtype=np.int32))
            cols["kind"].append(np.full(n, kind, dtype=np.int32))
            cols["t_ns"].append(0.5 * (ta + tb) if kind == mc.HIT_AB else (ta if a is not None else tb))
            cols["tA_ns"].append(ta)
            cols["tB_ns"].append(tb)
            cols["tdiff_ns"].append(ta - tb)
            for end, sel in (("A", a), ("B", b)):
                cols[f"h{end}"].append(P["height"][sel].astype(np.float64) if sel is not None else nan)
                cols[f"q{end}"].append(P["charge"][sel].astype(np.float64) if sel is not None else nan)
                cols[f"flags{end}"].append(P["flags"][sel].astype(np.int32) if sel is not None else
                                           np.full(n, -1, dtype=np.int32))
                cols[f"pulse{end}"].append(sel.astype(np.int64) if sel is not None else none)
                cols[f"ch{end}"].append(np.full(n, p.channels.get(end, -1), dtype=np.int32))
        per_paddle[p.name] = {"pulses_A": int(len(ia)), "pulses_B": int(len(ib)), "AB": int(len(pairs)),
                              "A_only": int(lone_a.sum()) if p.both_ends else int(len(ia)),
                              "B_only": int(lone_b.sum()) if p.both_ends else int(len(ib))}
        bar.update(1, p.name)
    bar.close()

    if not cols["paddle"]:
        sys.exit("ERROR: no usable pulses in any paddle of the setup")
    H = {k: np.concatenate(v) for k, v in cols.items()}
    H["event"] = group_events(H["t_ns"], window)
    order = np.lexsort((H["paddle"], H["t_ns"]))
    H = {k: v[order] for k, v in H.items()}
    H["hit"] = np.arange(len(H["t_ns"]), dtype=np.int64)
    n_ev = int(H["event"].max()) + 1

    # ---- events table ----
    E = {"event": np.arange(n_ev, dtype=np.int64)}
    first = np.full(n_ev, np.inf)
    last = np.full(n_ev, -np.inf)
    np.minimum.at(first, H["event"], H["t_ns"])
    np.maximum.at(last, H["event"], H["t_ns"])
    E["t_ns"] = first
    E["duration_ns"] = last - first
    E["n_hits"] = np.bincount(H["event"], minlength=n_ev).astype(np.int32)
    ab = H["kind"] == mc.HIT_AB
    layers_ab = np.zeros(n_ev, dtype=np.int64)
    paddles_ab = np.zeros(n_ev, dtype=np.int64)
    for layer in setup.layers:
        sel = H["layer"] == layer.index
        E[f"n_ab_{layer.name}"] = np.bincount(H["event"][sel & ab], minlength=n_ev).astype(np.int32)
        E[f"n_one_{layer.name}"] = np.bincount(H["event"][sel & ~ab], minlength=n_ev).astype(np.int32)
        layers_ab |= np.where(E[f"n_ab_{layer.name}"] > 0, 1 << layer.index, 0)
    if len(setup.paddles) <= 62:
        np.bitwise_or.at(paddles_ab, H["event"][ab], (1 << H["paddle"][ab].astype(np.int64)))
    else:
        paddles_ab[:] = -1
    E["layers_ab"] = layers_ab.astype(np.int32)
    E["paddles_ab"] = paddles_ab
    E["n_layers_ab"] = sum((layers_ab >> layer.index) & 1 for layer in setup.layers).astype(np.int32)

    events_path = os.path.join(folder, mc.EVENTS_FILE)
    mc.write_tree(events_path, "hits", H)
    mc.write_tree(events_path, "events", E, update=True)

    # ---- summary ----
    report()
    report("Hits per paddle (AB = both ends, A / B = one end only; both-ends share = AB / all hits)")
    head = f"{'paddle':<10}{'layer':<9}{'pulses A':>9}{'pulses B':>9}{'AB':>7}{'A only':>8}{'B only':>8}{'AB rate Hz':>12}{'both %':>8}"
    report(head)
    report("-" * len(head))
    for p in setup.paddles:
        c = per_paddle[p.name]
        allh = c["AB"] + c["A_only"] + c["B_only"]
        c["AB_rate_hz"] = c["AB"] / span_s if span_s > 0 else float("nan")
        c["both_percent"] = 100.0 * c["AB"] / allh if allh else float("nan")
        report(f"{p.name:<10}{p.layer.name:<9}{c['pulses_A']:>9}{c['pulses_B']:>9}{c['AB']:>7}{c['A_only']:>8}"
               f"{c['B_only']:>8}{c['AB_rate_hz']:>12.3f}{c['both_percent']:>8.1f}")

    report()
    report(f"{n_ev} events; by the layers with a two-end hit:")
    categories = {}
    for mask in sorted(set(int(m) for m in layers_ab)):
        names = [layer.name for layer in setup.layers if mask >> layer.index & 1]
        label = " + ".join(names) if names else "one-end hits only"
        n = int((layers_ab == mask).sum())
        categories[label] = n
        report(f"  {label:<30}{n:>8}{n / span_s:>12.4f} Hz")

    report()
    report(f"Coincidences between layers (two-end hits in both, within {window:g} ns), with the number "
           f"expected by chance:")
    coincidences = {}
    t_layer = {layer.index: np.sort(H["t_ns"][ab & (H["layer"] == layer.index)]) for layer in setup.layers}
    for i, li in enumerate(setup.layers):
        for lj in setup.layers[i + 1:]:
            both = int(((layers_ab >> li.index & 1) & (layers_ab >> lj.index & 1)).sum())
            t1, t2 = t_layer[li.index], t_layer[lj.index]
            shifted = [count_near(t1, t2, window, s * 1e3) for s in SHIFTS_US] if len(t1) and len(t2) else [0]
            chance_shift = float(np.mean(shifted))
            chance_calc = len(t1) * len(t2) * 2 * window * 1e-9 / span_s if span_s > 0 else float("nan")
            key = f"{li.name}+{lj.name}"
            coincidences[key] = {"events": both, "rate_hz": both / span_s if span_s > 0 else float("nan"),
                                 "chance_shifted": chance_shift, "chance_shifted_total": int(np.sum(shifted)),
                                 "chance_calculated": chance_calc}
            report(f"  {key:<22}{both:>7} events ({both / span_s:.4f} +- {np.sqrt(both) / span_s:.4f} Hz); "
                   f"by chance: {chance_shift:.3g} (shifted, {np.sum(shifted)} in {len(shifted)} shifts), "
                   f"{chance_calc:.3g} (N1 N2 2w / T)")

    pair_counts = {}
    multi = E["n_layers_ab"] >= 2
    if multi.any():
        sel = ab & multi[H["event"]]
        ev, pad, lay = H["event"][sel], H["paddle"][sel], H["layer"][sel]
        for e in np.unique(ev):
            m = ev == e
            ps = sorted(set(zip(lay[m].tolist(), pad[m].tolist())))
            for a in range(len(ps)):
                for b in range(a + 1, len(ps)):
                    if ps[a][0] != ps[b][0]:
                        key = f"{setup.paddles[ps[a][1]].name}+{setup.paddles[ps[b][1]].name}"
                        pair_counts[key] = pair_counts.get(key, 0) + 1
    if pair_counts:
        report()
        report("Events with two-end hits in two paddles of different layers, per paddle pair:")
        for key in sorted(pair_counts):
            report(f"  {key:<22}{pair_counts[key]:>7}")

    info = mc.stamp("build_events")
    info.update({"run": pinfo["run"], "span_ns": pinfo["span_ns"], "setup": {"name": setup.name, "path": setup.path},
                 "event_settings": setup.events, "v_eff_nominal": v_nom, "paddles": per_paddle,
                 "n_events": n_ev, "n_hits": int(len(H["t_ns"])), "categories": categories,
                 "coincidences": coincidences, "pair_counts": pair_counts})
    mc.write_info(events_path, info)
    with open(os.path.join(folder, "events_summary.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Saved {events_path} and events_summary.txt")
    return events_path


if __name__ == "__main__":
    main()
