#!/usr/bin/env python3
"""
make_fake_run.py  -  a simulated DAQ_YY run with known truth, to test the analysis

Simulates cosmic muons (I0 cos^2) and a low-energy background through the
paddles of a setup file, turns every crossing into pulses at both ends of the
paddle (light speed, attenuation, gains, time jitter, cable offsets) and the
pulses into self-trigger records of the AD9694 card (the measured pulse shape,
noise, threshold, pulse_pre / pulse_post). Writes, in the output folder:
    fake_run_YYYYMMDD_HHMMSS.root   in the raw2root_new format (trees ch<id>)
    fake_run_YYYYMMDD_HHMMSS.conf   the settings, as XIAOLAB DAQ writes them
    fake_run_YYYYMMDD_HHMMSS_truth.json and _truth.root   the true constants and muons
The file name carries the setup's valid_from date, so the analysis picks that setup.

Usage:
    python3 tests/make_fake_run.py /path/to/folder                 # 600 s, setup of 2026-11-12
    python3 tests/make_fake_run.py /path/to/folder --duration 1800 --seed 2
    python3 tests/make_fake_run.py /path/to/folder --swap-top-ends # wire TOP A and B the wrong way round
then
    python3 run_analysis.py /path/to/folder/fake_run_....root --out /path/to/folder/fake_muon
    python3 tests/check_fake_run.py /path/to/folder/fake_muon

The model is simple on purpose (no electronics response beyond one pulse shape,
no photon statistics beyond a smeared energy, no pile-up of a record with the
next pulse): it checks that the analysis recovers what was put in, signs and
frames included, not that the detector is understood.
"""

import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse
import json
import math
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import muon_common as mc  # noqa: E402

BASELINE = 8191.0
LATENCY = 6            # the trigger sample (2 x pulse_pre) comes 6 samples after the threshold crossing
TAU_RISE, TAU_DECAY = 4.0, 840.0     # ns: the measured average pulse (run 56/57)


def shape(t):
    """Pulse of height ~1 starting at t = 0 (ns)."""
    t = np.asarray(t, dtype=np.float64)
    s = np.where(t > 0, (1 - np.exp(-np.maximum(t, 0) / TAU_RISE)) * np.exp(-np.maximum(t, 0) / TAU_DECAY), 0.0)
    return s


PEAK = float(shape(np.linspace(0, 200, 20001)).max())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", help="Folder for the simulated run")
    ap.add_argument("--setup", default=None, help="Setup file (default: setups/two_layer_2026-11-12.toml)")
    ap.add_argument("--duration", type=float, default=600.0, help="Seconds (default 600)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--v", type=float, default=16.0, help="Signal speed, cm/ns (default 16)")
    ap.add_argument("--jitter", type=float, default=0.45, help="Time jitter of one end, ns (default 0.45)")
    ap.add_argument("--attenuation", type=float, default=38.0, help="Attenuation length, cm (default 38)")
    ap.add_argument("--mip", type=float, default=850.0,
                    help="sqrt(hA hB) of a vertical muon, ADC (default 850)")
    ap.add_argument("--background", type=float, default=8.0,
                    help="Low-energy hits per second and paddle (default 8)")
    ap.add_argument("--noise", type=float, default=20.0, help="Baseline noise, ADC rms (default 20)")
    ap.add_argument("--threshold", type=int, default=400, help="Threshold of every channel, ADC (default 400)")
    ap.add_argument("--pulse-pre", type=int, default=50)
    ap.add_argument("--pulse-post", type=int, default=150)
    ap.add_argument("--swap-top-ends", action="store_true",
                    help="Swap the A and B cables of the top paddles (to see the map check catch it)")
    args = ap.parse_args()
    setup = mc.load_setup(args.setup or os.path.join(mc.SETUPS, "two_layer_2026-11-12.toml"))
    rng = np.random.default_rng(args.seed)
    os.makedirs(args.out, exist_ok=True)
    day = (setup.valid_from or __import__("datetime").date(2026, 11, 12)).strftime("%Y%m%d")
    name = f"fake_run_{day}_120000"
    base = os.path.join(os.path.abspath(args.out), name)
    T = args.duration * 1e9                            # ns

    # ---- true constants ----
    truth = {"v_eff_cm_per_ns": args.v, "jitter_ns": args.jitter, "attenuation_cm": args.attenuation,
             "mip_adc": args.mip, "noise_adc": args.noise, "threshold": args.threshold,
             "duration_s": args.duration, "background_hz_per_paddle": args.background,
             "flux": setup.flux, "swap_top_ends": args.swap_top_ends, "paddles": {}, "channels": {}}
    for p in setup.paddles:
        gA, gB = math.exp(rng.normal(0, 0.12)), math.exp(rng.normal(0, 0.12))
        oA, oB = rng.normal(0, 0.4), rng.normal(0, 0.4)
        truth["paddles"][p.name] = {"gain_A": gA, "gain_B": gB, "offset_A_ns": oA, "offset_B_ns": oB,
                                    "tdiff_offset_ns": oA - oB, "gain_ratio_A_over_B": gA / gB}

    # ---- muons on a plane above everything, wide enough for every paddle ----
    z_top = max(p.z for p in setup.paddles)
    z_bot = min(p.z for p in setup.paddles)
    margin = 2.5 * (z_top - z_bot) + 10
    x0 = min(p.x - 0.5 * p.size_x for p in setup.paddles) - margin
    x1 = max(p.x + 0.5 * p.size_x for p in setup.paddles) + margin
    y0 = min(p.y - 0.5 * p.size_y for p in setup.paddles) - margin
    y1 = max(p.y + 0.5 * p.size_y for p in setup.paddles) + margin
    power = float(setup.flux["cos_power"])
    rate = float(setup.flux["vertical_intensity"]) * 1e-4 * (x1 - x0) * (y1 - y0) * 2 * math.pi / (power + 2)
    n = rng.poisson(rate * args.duration)
    tm = np.sort(rng.uniform(0, T, n))
    xs, ys = rng.uniform(x0, x1, n), rng.uniform(y0, y1, n)
    cos_t = rng.random(n) ** (1 / (power + 2))
    tan_t = np.sqrt(1 - cos_t ** 2) / cos_t
    ph = rng.uniform(0, 2 * math.pi, n)
    bx, by = tan_t * np.cos(ph), tan_t * np.sin(ph)
    crossings = []           # (time at the paddle, paddle index, u, energy factor)
    muon_rows = {k: [] for k in ("t_ns", "theta_deg", "phi_deg")}
    for layer in setup.layers:
        muon_rows[f"paddle_{layer.name}"] = []
        muon_rows[f"x_{layer.name}"] = []
        muon_rows[f"y_{layer.name}"] = []
    hit_any = np.zeros(n, dtype=bool)
    per_paddle = {}
    for p in setup.paddles:
        px, py = xs + bx * (p.z - z_top), ys + by * (p.z - z_top)
        inside = p.contains(px, py)
        per_paddle[p.index] = (inside, px, py)
        hit_any |= inside
    idx = np.where(hit_any)[0]
    for i in idx:
        for layer in setup.layers:
            muon_rows[f"paddle_{layer.name}"].append(-1)
            muon_rows[f"x_{layer.name}"].append(np.nan)
            muon_rows[f"y_{layer.name}"].append(np.nan)
        theta = math.degrees(math.atan(tan_t[i]))
        muon_rows["t_ns"].append(tm[i])
        muon_rows["theta_deg"].append(theta)
        muon_rows["phi_deg"].append(math.degrees(math.atan2(by[i], bx[i])) % 360)
        for p in setup.paddles:
            inside, px, py = per_paddle[p.index]
            if not inside[i]:
                continue
            path_to = (z_top - p.z) / cos_t[i]
            u = float(p.u_of(px[i], py[i]))
            crossings.append((tm[i] + path_to / mc.C_CM_PER_NS, p.index, u, 1.0 / cos_t[i]))
            muon_rows[f"paddle_{p.layer.name}"][-1] = p.index
            muon_rows[f"x_{p.layer.name}"][-1] = float(px[i])
            muon_rows[f"y_{p.layer.name}"][-1] = float(py[i])
    n_muon_crossings = len(crossings)
    # Landau-like energy loss (Moyal, inverse CDF), longer path when inclined
    from scipy.special import erfcinv
    z = -2.0 * np.log(math.sqrt(2.0) * erfcinv(rng.uniform(1e-9, 1 - 1e-9, n_muon_crossings)))
    energy = 1.0 + 0.09 * z
    crossings = [(t, k, u, path * e) for (t, k, u, path), e in zip(crossings, energy)]
    # low-energy background: evenly along each paddle, exponential spectrum
    for p in setup.paddles:
        nb = rng.poisson(args.background * args.duration)
        for t, u, e in zip(rng.uniform(0, T, nb), rng.uniform(-0.5 * p.length, 0.5 * p.length, nb),
                           rng.exponential(0.45, nb)):
            crossings.append((t, p.index, u, e))

    # ---- pulses at the ends ----
    pulses = {chid: [] for chid in setup.channels}
    lam, v = args.attenuation, args.v
    for t, k, u, e in crossings:
        p = setup.paddles[k]
        tp = truth["paddles"][p.name]
        for end in ("A", "B"):
            if end not in p.channels:
                continue
            d = 0.5 * p.length - u if end == "A" else 0.5 * p.length + u
            g = tp[f"gain_{end}"]
            h = args.mip * math.exp(0.5 * p.length / lam) * g * e * math.exp(-d / lam)
            te = t + d / v + tp[f"offset_{end}_ns"] + rng.normal(0, args.jitter)
            chid = p.channels[end]
            if args.swap_top_ends and p.layer.index == 0:
                chid = p.channels["B" if end == "A" else "A"]
            pulses[chid].append((te, h))

    # ---- records ----
    pre2, post2 = 2 * args.pulse_pre, 2 * args.pulse_post
    out = mc.ROOT.TFile(base + ".root", "RECREATE")
    MAX = 20000
    counts = {}
    bar = mc.Progress("Writing records", len(pulses))
    for chid in sorted(pulses):
        tree = mc.ROOT.TTree(f"ch{chid}", f"ch{chid}")
        buf = {"ttt": np.zeros(1, np.uint64), "tbusy": np.zeros(1, np.uint64), "nsamps": np.zeros(1, np.uint32),
               "x_numb": np.arange(MAX, dtype=np.uint32), "chdata": np.zeros(MAX, np.uint16),
               "t2prev": np.zeros(1, np.float64), "xmin": np.zeros(1, np.uint32), "xmax": np.zeros(1, np.uint32),
               "rms": np.zeros(1, np.float32), "baseline": np.zeros(1, np.float32), "area": np.zeros(1, np.float64)}
        for key, leaf in (("ttt", "ttt/l"), ("tbusy", "tbusy/l"), ("nsamps", "nsamps/i"),
                          ("x_numb", "x_numb[nsamps]/i"), ("chdata", "chdata[nsamps]/s"), ("t2prev", "t2prev/D"),
                          ("xmin", "xmin/i"), ("xmax", "xmax/i"), ("rms", "rms/F"), ("baseline", "baseline/F"),
                          ("area", "area/D")):
            tree.Branch(key, buf[key], leaf)
        plist = sorted(pulses[chid])
        tt = np.array([q[0] for q in plist])
        hh = np.array([q[1] for q in plist])
        busy_until, last_ttt, n_rec = -1.0, None, 0
        for te, h in plist:
            if h * PEAK <= args.threshold or te <= busy_until:
                continue
            # first 2 ns sample where the pulse crosses the threshold (noise-free), then the trigger sample
            grid = np.arange(math.floor(te / 2) - 2, math.floor(te / 2) + 40) * 2.0
            above = np.where(h * shape(grid - te) > args.threshold)[0]
            if not len(above):
                continue
            k_cross = int(round(grid[above[0]] / 2))
            k_trig = k_cross + LATENCY
            k_trig += k_trig % 2                       # ttt steps by 2 ticks
            t_end = te + TAU_DECAY * math.log(max(h * PEAK / args.threshold, 1.0001))
            k_below = int(t_end / 2) + 1
            k0 = k_trig - pre2
            k1 = max(k_below, k_trig) + post2 + 2
            ns = min(k1 - k0, MAX)
            ts = (k0 + np.arange(ns)) * 2.0
            w = BASELINE + h * shape(ts - te) + rng.normal(0, args.noise, ns)
            # later pulses of the same channel inside this record
            j0, j1 = np.searchsorted(tt, te, "right"), np.searchsorted(tt, ts[-1], "left")
            for te2, h2 in zip(tt[j0:j1], hh[j0:j1]):
                w += h2 * shape(ts - te2)
            w = np.clip(np.round(w), 0, mc.ADC_MAX).astype(np.uint16)
            buf["chdata"][:ns] = w
            buf["nsamps"][0] = ns
            buf["ttt"][0] = k_trig
            buf["t2prev"][0] = 0.0 if last_ttt is None else (k_trig - last_ttt) * 2e-9
            buf["baseline"][0] = float(w[:20].mean())
            buf["rms"][0] = float(w[:20].std(ddof=1))
            buf["xmin"][0], buf["xmax"][0] = 0, ns - 1
            buf["area"][0] = float((buf["baseline"][0] - w[35:101]).sum()) if ns > 101 else 0.0
            tree.Fill()
            last_ttt = k_trig
            busy_until = ts[-1]
            n_rec += 1
        tree.Write()
        counts[chid] = n_rec
        bar.update(1, f"ch{chid}")
    bar.close()
    out.Close()

    with open(base + ".conf", "w", encoding="utf-8") as f:
        f.write("# testDAQ configuration, written by tests/make_fake_run.py (simulated run)\n")
        f.write("# ---- note ----\n# simulated run for testing the analysis\n# ---- end of note ----\n")
        f.write("global_config = {\n")
        f.write(f"    pulse_pre={args.pulse_pre};\n    pulse_post={args.pulse_post};\n    ctrlmode=3;\n")
        f.write(f"    trigmode=0;\n    trig_pre=100;\n    trig_post=100;\n    daqruntime={int(args.duration)};\n}};\n")
        f.write("board_config = (\n")
        boards = sorted({mc.channel_board(c) for c in setup.channels})
        for b in boards:
            f.write(f"{{   brdID={b};\n    channelenable=0xFF;\n    pulsePolarity=0xFF;\n")
            for ch in range(8):
                f.write(f"    ch{ch}_threshold={args.threshold};\n")
            f.write("}" + ("," if b != boards[-1] else "") + "\n")
        f.write(");\n")
    truth["records"] = {str(k): v for k, v in counts.items()}
    truth["muons_crossing_a_paddle"] = len(muon_rows["t_ns"])
    truth["muon_crossings"] = n_muon_crossings
    mc.write_json(base + "_truth.json", truth)
    cols = {k: np.asarray(v, dtype=np.int32 if k.startswith("paddle_") else np.float64) for k, v in muon_rows.items()}
    mc.write_tree(base + "_truth.root", "muons", cols)
    print(f"{len(muon_rows['t_ns'])} muons crossed a paddle ({n_muon_crossings} crossings), "
          f"{sum(counts.values())} records in {len(counts)} channels")
    print(f"Saved {base}.root, .conf, _truth.json, _truth.root")


if __name__ == "__main__":
    main()
