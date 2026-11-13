#!/usr/bin/env python3
"""
reconstruct.py  -  step 4 of the paddle muon analysis

Reconstructs a straight muon track for every event with two-end hits in two
layers or more: the hit point in each layer, the track's zenith and azimuth
angles, its path length, time of flight and speed. Also predicts, with a toy
cosmic-muon simulation of the setup's geometry, the rates and angle
distributions expected for fully efficient paddles. Reads events.root (step 2)
and calibration.json (step 3); writes tracks.root and tracks_summary.txt.

Usage:
    python3 reconstruct.py                     # asks for the analysis folder
    python3 reconstruct.py /path/to/run00056_..._muon
    python3 reconstruct.py FOLDER --calibration /other/run_muon/calibration.json

Method (cm, ns, degrees; u = position along a paddle from its centre, positive
towards its A end):
    hit point    measured u = -(v / 2) x (tA - tB - offset) with the paddle's v and
                 offset from the calibration; its resolution is (v / 2) x
                 sigma(tA - tB), so near an end the measured u often lies beyond
                 the paddle (about 8 % of the hits in runs 56-58, as the resolution
                 predicts). The muon did cross inside, so the point used is the
                 mean of the measurement's Gaussian cut to the paddle, with that
                 cut Gaussian's spread as its error (always inside; the measured u
                 is kept as u_meas). Across the paddle the point is its centre line
                 (resolution width / sqrt 12), in height its centre plane.
    crossing     t = (tA + tB) / 2 - L / (2 v) - the paddle's time offset
    one layer    several two-end hits in neighbouring paddles (within 5 ns) are
                 merged into one point weighted by sqrt(hA hB) (flag ADJACENT);
                 hits in paddles that are not neighbours: the biggest group is
                 used (flag AMBIGUOUS)
    track        straight lines x(z) and y(z) through the points, weighted by
                 1 / sigma^2 (two layers: the line through both points)
    angles       zenith theta from the vertical; azimuth phi of the direction the
                 muon came from, counted from +x towards +y (0 = from +x, 90 =
                 from +y: the wall side in the first setup)
    speed        time of flight between the highest and the lowest layer, path
                 s = height difference x sqrt(1 + slope_x^2 + slope_y^2),
                 beta = s / (c x tof), as measured: no layer time offsets are
                 applied (the paddles' cable delays are unknown; setting them
                 from beta = 1 would force the result). With a 30 cm lever arm
                 the time resolution makes many tracks look faster than light.
    expected     toy muons, I(theta) = I0 cos^n(theta) (setup [flux]), through the
                 paddles; rates for 100 % efficient paddles, and their angles
                 reconstructed like the data (see muon_common.toy_cosmics)

Outputs (in the analysis folder):
    tracks.root          tree 'tracks' (one entry per reconstructed muon), tree
                         'points' (every two-end hit of the tracks' layers, each at
                         its own paddle; used = 1 if it went into the track's point),
                         tree 'toy' (the toy muons that cross two layers, with
                         weights in Hz), and 'info'
    tracks_summary.txt   tracks, angles, speed, measured against expected rates
"""

import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse
import math
import sys

import numpy as np

import muon_common as mc

MERGE_DT_NS = 5.0     # hits of neighbouring paddles within this time are one point


def hit_table(H, setup, calib):
    """Point, resolutions and crossing time of every two-end hit."""
    n = len(H["t_ns"])
    out = {k: np.full(n, np.nan) for k in ("x", "y", "z", "u", "u_meas", "t", "sx", "sy", "uamp", "hgeo", "hmip")}
    out["outside"] = np.zeros(n, dtype=bool)
    cal = calib["paddles"]
    for p in setup.paddles:
        sel = (H["paddle"] == p.index) & (H["kind"] == mc.HIT_AB)
        if not sel.any():
            continue
        c = cal.get(p.name, {})
        v = c.get("v_eff_cm_per_ns", setup.calibration["v_eff_cm_per_ns"])
        off = c.get("tdiff_offset_ns", 0.0)
        sig_u = mc.position_sigma(calib, p, setup)
        u_meas = -0.5 * v * (H["tdiff_ns"][sel] - off)
        # the muon crossed inside the paddle: best estimate and its error from the measurement
        u, s_u = mc.position_in_paddle(u_meas, 0.5 * p.length, sig_u)
        out["u"][sel], out["u_meas"][sel] = u, u_meas
        out["x"][sel], out["y"][sel], out["z"][sel] = p.point(u)
        out["t"][sel] = H["t_ns"][sel] - 0.5 * p.length / v - c.get("time_offset_ns", 0.0)
        along_x = p.axis == "x"
        out["sx"][sel] = s_u if along_x else p.width / math.sqrt(12)
        out["sy"][sel] = p.width / math.sqrt(12) if along_x else s_u
        out["outside"][sel] = np.abs(u_meas) > 0.5 * p.length + 3 * sig_u
        hgeo = np.sqrt(H["hA"][sel] * H["hB"][sel])
        out["hgeo"][sel] = hgeo
        if c.get("mip_adc"):
            out["hmip"][sel] = hgeo / c["mip_adc"]
        if c.get("attenuation_cm"):
            ratio = np.log(H["hA"][sel] / H["hB"][sel]) - c.get("log_ratio_offset", 0.0)
            out["uamp"][sel] = 0.5 * c["attenuation_cm"] * ratio
    return out


def layer_point(idx, H, T, setup):
    """One point for the two-end hits `idx` of one layer in one event, and its flags."""
    if len(idx) == 1:
        i = idx[0]
        return {k: T[k][i] for k in ("x", "y", "z", "sx", "sy", "u", "u_meas", "t", "uamp", "hgeo", "hmip")} | \
            {"paddle": int(H["paddle"][i]), "hit": int(H["hit"][i]), "outside": bool(T["outside"][i])}, 0, [i]
    # group neighbouring paddles that fired together
    groups = []
    for i in sorted(idx, key=lambda i: -T["hgeo"][i]):
        p = setup.paddles[int(H["paddle"][i])]
        for g in groups:
            if any(setup.paddles[int(H["paddle"][j])].neighbours(p) and
                   abs(T["t"][j] - T["t"][i]) < MERGE_DT_NS for j in g):
                g.append(i)
                break
        else:
            groups.append([i])
    flags = mc.TRACK_ADJACENT if any(len(g) > 1 for g in groups) else 0
    if len(groups) > 1:
        flags |= mc.TRACK_AMBIGUOUS
    best = max(groups, key=lambda g: sum(T["hgeo"][j] for j in g))
    w = np.array([T["hgeo"][j] for j in best])
    w = w / w.sum()
    pt = {k: float(np.sum(w * np.array([T[k][j] for j in best]))) for k in ("x", "y", "z", "u", "u_meas", "t")}
    for k in ("uamp", "hmip"):
        vals = np.array([T[k][j] for j in best])
        pt[k] = float(np.nansum(vals * w)) if np.isfinite(vals).any() else float("nan")
    pt["hgeo"] = float(np.sum([T["hgeo"][j] for j in best]))
    if len(best) > 1:   # the point lies between the paddles' centre lines: spread as resolution
        for k, s in (("x", "sx"), ("y", "sy")):
            spread = np.sqrt(np.sum(w * (np.array([T[k][j] for j in best]) - pt[k]) ** 2))
            pt[s] = float(max(np.max([T[s][j] for j in best]), spread))
    else:
        pt["sx"], pt["sy"] = float(T["sx"][best[0]]), float(T["sy"][best[0]])
    i0 = best[0]
    pt.update({"paddle": int(H["paddle"][i0]), "hit": int(H["hit"][i0]),
               "outside": bool(any(T["outside"][j] for j in best))})
    return pt, flags, list(best)


def weighted_mean(v, w=None):
    """Mean and its error (weights = toy weights; None = data, plain mean)."""
    v = np.asarray(v, dtype=np.float64)
    if w is None:
        return (float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v)))) if len(v) > 1 else (float("nan"),) * 2
    m = float(np.average(v, weights=w))
    return m, 0.0


def map_check(cols, toy, setup, calib):
    """Is the channel map (which paddle is where, which end is A) consistent with the tracks?

    Where a muon crosses a layer along its paddles depends on the paddle it crosses
    in another layer (it came down at an angle, from near above that paddle). For
    each layer L and each other layer M with two paddles or more, the mean position
    along L's paddles is taken for the tracks through each paddle of M, in the data
    and in the toy, and the data means are fitted as a + b x the toy means. With a
    right map b is about +1; b about -1 says that the direction along L's paddles
    (their A ends) or the order of M's paddles is the other way round. The position
    along L's paddles is taken from the timing and, independently, from the light
    sharing ln(hA / hB). a is the shift of the whole layer against the toy."""
    rows = []
    if toy is None:
        return rows
    ok = (cols["flags"] & mc.TRACK_ADJACENT) == 0
    for L in setup.layers:
        axes = {p.axis for p in L.paddles}
        if len(axes) != 1:
            continue
        axis = axes.pop()
        coord = f"{axis}_{L.name}"
        for M in setup.layers:
            if M is L or len(M.paddles) < 2:
                continue
            for method in ("timing", "light"):
                if method == "timing":
                    data = cols[coord]
                else:   # position from ln(hA / hB): centre + a_sign x u_amp along the axis
                    data = np.full(len(cols[coord]), np.nan)
                    for p in L.paddles:
                        sel = cols[f"paddle_{L.name}"] == p.index
                        centre = p.x if axis == "x" else p.y
                        data[sel] = centre + p.a_sign * cols[f"uamp_{L.name}"][sel]
                d_m, d_e, t_m, names = [], [], [], []
                for q in M.paddles:
                    sel = ok & (cols[f"paddle_{M.name}"] == q.index) & np.isfinite(data)
                    tsel = toy[f"paddle_{M.name}"] == q.index
                    if sel.sum() < 5 or not tsel.any():
                        continue
                    m, e = weighted_mean(data[sel])
                    tm, _ = weighted_mean(toy[coord][tsel], toy["weight_hz"][tsel])
                    d_m.append(m), d_e.append(e), t_m.append(tm), names.append((q.name, int(sel.sum())))
                if len(d_m) < 2:
                    continue
                d_m, d_e, t_m = np.array(d_m), np.array(d_e), np.array(t_m)
                w = 1.0 / np.maximum(d_e, 1e-6) ** 2
                tbar = np.sum(w * t_m) / w.sum()
                dbar = np.sum(w * d_m) / w.sum()
                sxx = np.sum(w * (t_m - tbar) ** 2)
                if sxx <= 0:
                    continue
                b = float(np.sum(w * (t_m - tbar) * (d_m - dbar)) / sxx)
                b_err = float(1.0 / math.sqrt(sxx))
                shift = float(dbar - tbar)
                shift_err = float(1.0 / math.sqrt(w.sum()))
                plus, zero, minus = (abs(b - s) > 2 * b_err for s in (1.0, 0.0, -1.0))   # excluded at 2 sigma
                check = (f"check the A/B cables of the {L.name} paddles, the order of the {M.name} paddles and "
                         f"their positions in the setup file")
                if not plus:
                    verdict = "as expected" if minus else "not decided yet (more tracks needed)"
                elif not minus and zero:
                    verdict = f"REVERSED (slope -1; +1 excluded at {abs(b - 1) / b_err:.1f} sigma): {check}"
                elif not zero:
                    verdict = (f"+1 excluded at {abs(b - 1) / b_err:.1f} sigma, slope near 0"
                               f"{' or -1' if not minus else ''}: the {L.name} positions do not follow the "
                               f"{M.name} paddles; {check}")
                else:
                    verdict = f"neither +1 nor -1 (slope {b:.2f}): check the geometry in the setup file"
                rows.append({"layer": L.name, "axis": axis, "by_layer": M.name, "method": method,
                             "paddles": [{"name": n, "tracks": k, "data_cm": float(dm), "data_err_cm": float(de),
                                          "toy_cm": float(tm)} for (n, k), dm, de, tm in zip(names, d_m, d_e, t_m)],
                             "slope": b, "slope_err": b_err, "shift_cm": shift, "shift_err_cm": shift_err,
                             "verdict": verdict})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", default=None,
                    help="Analysis folder with events.root and calibration.json (asked for if left out)")
    ap.add_argument("--setup", default=None, help="Setup file (default: the one step 1 used)")
    ap.add_argument("--calibration", default=None,
                    help="calibration.json to use (default: the one in the analysis folder)")
    ap.add_argument("--toy-muons", type=int, default=None,
                    help="Toy muons per layer for the expectation (default: the setup's [flux] toy_muons)")
    args = ap.parse_args()
    folder = mc.ask_analysis_folder(args.folder, mc.EVENTS_FILE, "build_events.py")
    run(folder, args.setup, args.calibration, args.toy_muons)


def run(folder, setup_path=None, calib_path=None, toy_muons=None):
    """Step 4 for one analysis folder; returns the path of tracks.root."""
    events_path = mc.need_file(folder, mc.EVENTS_FILE, "build_events.py")
    pinfo = mc.read_info(mc.need_file(folder, mc.PULSES_FILE, "extract_pulses.py"))
    einfo = mc.read_info(events_path)
    setup = mc.setup_of(folder, setup_path)
    calib_path = os.path.abspath(calib_path) if calib_path else mc.need_file(folder, mc.CALIB_FILE, "calibrate.py")
    if not os.path.isfile(calib_path):
        sys.exit(f"ERROR: no such calibration file: {calib_path}")
    calib = mc.read_json(calib_path)
    H = mc.read_tree(events_path, "hits")
    E = mc.read_tree(events_path, "events")
    span_s = pinfo["span_ns"] * 1e-9
    nl = len(setup.layers)

    lines = []

    def report(text=""):
        print(text)
        lines.append(text)

    report(f"Analysis folder : {folder}")
    report(f"Run             : {pinfo['run']} ({span_s:.3f} s)")
    report(f"Setup           : {setup.name}")
    report(f"Calibration     : {calib_path}" + ("" if os.path.dirname(calib_path) == os.path.abspath(folder)
                                                 else f" (from run {calib.get('info', {}).get('run', '?')})"))

    T = hit_table(H, setup, calib)
    cand = np.where(E["n_layers_ab"] >= 2)[0]
    ab_idx = np.where(H["kind"] == mc.HIT_AB)[0]
    by_event = {}
    for i in ab_idx:
        e = int(H["event"][i])
        if E["n_layers_ab"][e] >= 2:
            by_event.setdefault(e, []).append(i)
    n = len(cand)
    cols = {"event": cand.astype(np.int64), "t_ns": E["t_ns"][cand].astype(np.float64),
            "flags": np.zeros(n, dtype=np.int32), "n_layers": np.zeros(n, dtype=np.int32)}
    per_layer = {}
    for layer in setup.layers:
        for k in ("paddle", "hit"):
            per_layer[f"{k}_{layer.name}"] = np.full(n, -1, dtype=np.int64 if k == "hit" else np.int32)
        for k in ("x", "y", "z", "u", "u_meas", "t", "sx", "sy", "uamp", "hgeo", "hmip"):
            per_layer[f"{k}_{layer.name}"] = np.full(n, np.nan)
    pts = {k: [] for k in ("event", "track", "layer", "paddle", "hit", "x", "y", "z", "u", "u_meas", "hgeo", "used")}
    bar = mc.Progress("Tracks", max(n, 1))
    for row, e in enumerate(cand):
        idx = by_event.get(int(e), [])
        for layer in setup.layers:
            li = [i for i in idx if H["layer"][i] == layer.index]
            if not li:
                continue
            pt, fl, used = layer_point(li, H, T, setup)
            for i in li:      # every two-end hit of this layer, at its own paddle
                for k, v in (("event", int(e)), ("track", row), ("layer", layer.index), ("paddle", int(H["paddle"][i])),
                             ("hit", int(H["hit"][i])), ("x", T["x"][i]), ("y", T["y"][i]), ("z", T["z"][i]),
                             ("u", T["u"][i]), ("u_meas", T["u_meas"][i]), ("hgeo", T["hgeo"][i]),
                             ("used", 1 if i in used else 0)):
                    pts[k].append(v)
            cols["flags"][row] |= fl | (mc.TRACK_OUTSIDE if pt["outside"] else 0)
            cols["n_layers"][row] += 1
            for k in ("paddle", "hit", "x", "y", "z", "u", "u_meas", "t", "sx", "sy", "uamp", "hgeo", "hmip"):
                per_layer[f"{k}_{layer.name}"][row] = pt[k]
        bar.update(1)
    bar.close()
    cols.update(per_layer)

    names = [layer.name for layer in setup.layers]
    X = np.stack([cols[f"x_{nm}"] for nm in names], axis=1).reshape(n, nl)
    Y = np.stack([cols[f"y_{nm}"] for nm in names], axis=1).reshape(n, nl)
    Z = np.stack([cols[f"z_{nm}"] for nm in names], axis=1).reshape(n, nl)
    SX = np.stack([cols[f"sx_{nm}"] for nm in names], axis=1).reshape(n, nl)
    SY = np.stack([cols[f"sy_{nm}"] for nm in names], axis=1).reshape(n, nl)
    Tt = np.stack([cols[f"t_{nm}"] for nm in names], axis=1).reshape(n, nl)
    fit = mc.fit_lines(X, Y, Z, np.where(np.isfinite(SX), SX, 1.0), np.where(np.isfinite(SY), SY, 1.0))
    theta, phi, sig_theta = mc.direction(fit["bx"], fit["by"], fit["var_bx"], fit["var_by"])
    with np.errstate(invalid="ignore"):
        hi = np.argmax(np.where(np.isfinite(Z), Z, -np.inf), axis=1) if n else np.zeros(0, dtype=int)
        lo = np.argmin(np.where(np.isfinite(Z), Z, np.inf), axis=1) if n else np.zeros(0, dtype=int)
        r = np.arange(n)
        dz = Z[r, hi] - Z[r, lo]
        s = dz * np.sqrt(1 + fit["bx"] ** 2 + fit["by"] ** 2)
        tof = Tt[r, lo] - Tt[r, hi]
        beta = s / (mc.C_CM_PER_NS * tof)
        inv_beta = mc.C_CM_PER_NS * tof / s
    cols.update({"ax": fit["ax"], "ay": fit["ay"], "z0": fit["z0"], "bx": fit["bx"], "by": fit["by"],
                 "theta_deg": theta, "phi_deg": phi,
                 "sigma_theta_deg": sig_theta, "chi2": fit["chi2"], "ndf": fit["ndf"].astype(np.int32),
                 "s_cm": s, "tof_ns": tof, "tof_expected_ns": s / mc.C_CM_PER_NS, "beta": beta,
                 "inv_beta": inv_beta})
    cols = {k: (np.asarray(v, dtype=np.float64) if np.asarray(v).dtype.kind == "f" else np.asarray(v))
            for k, v in cols.items()}

    tracks_path = os.path.join(folder, mc.TRACKS_FILE)
    mc.write_tree(tracks_path, "tracks", cols)
    ints = ("event", "track", "hit")
    mc.write_tree(tracks_path, "points", {k: np.asarray(v, dtype=np.int64 if k in ints else
                                                       (np.int32 if k in ("layer", "paddle", "used") else np.float64))
                                          for k, v in pts.items()}, update=True)

    # ---- toy expectation ----
    print("Toy cosmic muons for the expectation ...")
    toy = mc.toy_cosmics(setup, calib, toy_muons)
    tcols = {}
    if toy is not None:
        tcols = {"weight_hz": toy["weight"], "theta_deg": toy["theta"], "phi_deg": toy["phi"],
                 "theta_true_deg": toy["theta_true"], "phi_true_deg": toy["phi_true"], "s_cm": toy["s"],
                 "tof_true_ns": toy["tof_true"], "multi": toy["multi"].astype(np.int32)}
        for layer in setup.layers:
            tcols[f"paddle_{layer.name}"] = toy["paddles"][:, layer.index].astype(np.int32)
            tcols[f"x_{layer.name}"] = toy["x"][:, layer.index]
            tcols[f"y_{layer.name}"] = toy["y"][:, layer.index]
            tcols[f"u_true_{layer.name}"] = toy["u_true"][:, layer.index]   # true position along its paddle
        mc.write_tree(tracks_path, "toy", tcols, update=True)

    # ---- summary ----
    ok = np.isfinite(theta)
    report()
    report(f"{n} tracks (events with two-end hits in {'two' if nl == 2 else 'two or more'} layers), "
           f"{n / span_s:.4f} +- {math.sqrt(n) / span_s:.4f} Hz")
    for bit, name, text in mc.TRACK_FLAGS:
        k = int(((cols["flags"] & bit) != 0).sum())
        if k:
            report(f"  {name:<10}{k:>6}  {text}")
    if ok.any():
        report(f"Zenith angle    : mean {theta[ok].mean():.1f} deg, median {np.median(theta[ok]):.1f} deg, "
               f"typical error {np.median(sig_theta[ok]):.1f} deg")
        good_tof = np.isfinite(inv_beta)
        if good_tof.any():
            report(f"Time of flight  : median {np.median(tof[good_tof]):.3f} ns over a median path of "
                   f"{np.median(s[good_tof]):.1f} cm (s/c = {np.median(s[good_tof]) / mc.C_CM_PER_NS:.3f} ns)")
            ib = inv_beta[good_tof]
            report(f"1/beta          : median {np.median(ib):.3f}, spread (robust sigma) "
                   f"{1.4826 * np.median(np.abs(ib - np.median(ib))):.3f}; tof - s/c robust sigma "
                   f"{1.4826 * np.median(np.abs((tof - s / mc.C_CM_PER_NS)[good_tof] - np.median((tof - s / mc.C_CM_PER_NS)[good_tof]))):.3f} ns")

    expected = {}
    if toy is not None:
        w = toy["weight"]
        report()
        report(f"Expected for 100 % efficient paddles (toy: I0 = {setup.flux['vertical_intensity']:g} /(m^2 s sr), "
               f"cos^{setup.flux['cos_power']:g}(theta), {len(w)} toy muons crossing two layers):")
        head = f"  {'paddles':<22}{'measured':>9}{'expected':>10}{'ratio':>8}   (counts in {span_s:.1f} s)"
        report(head)
        crossing = {}
        for i, pi in enumerate(setup.paddles):
            for pj in setup.paddles[i + 1:]:
                if pi.layer is pj.layer:
                    continue
                a, b = (pi, pj) if pi.layer.index < pj.layer.index else (pj, pi)
                sel = (toy["paddles"][:, a.layer.index] == a.index) & (toy["paddles"][:, b.layer.index] == b.index)
                rate = float(w[sel].sum())
                meas = int(((cols[f"paddle_{a.layer.name}"] == a.index) &
                            (cols[f"paddle_{b.layer.name}"] == b.index)).sum())
                crossing[f"{a.name}+{b.name}"] = {"expected_hz": rate, "measured": meas,
                                                  "expected_count": rate * span_s}
                ratio = meas / (rate * span_s) if rate > 0 else float("nan")
                report(f"  {a.name + ' + ' + b.name:<22}{meas:>9}{rate * span_s:>10.2f}{ratio:>8.2f}")
        total_rate = float(w.sum())
        ratio = n / (total_rate * span_s) if total_rate > 0 else float("nan")
        err = math.sqrt(max(n, 1)) / (total_rate * span_s) if total_rate > 0 else float("nan")
        report(f"  {'all':<22}{n:>9}{total_rate * span_s:>10.2f}{ratio:>8.2f} +- {err:.2f}")
        flux = setup.flux["vertical_intensity"] * ratio
        report(f"  -> measured / expected = {ratio:.2f} +- {err:.2f}: the product of the paddles' efficiencies and "
               f"the true flux over the assumed one")
        report(f"     (vertical intensity {flux:.0f} +- {setup.flux['vertical_intensity'] * err:.0f} /(m^2 s sr) "
               f"if the paddles were fully efficient)")
        th = toy["theta"]
        report(f"  toy zenith angle (reconstructed like the data): mean {np.average(th, weights=w):.1f} deg; "
               f"true mean {np.average(toy['theta_true'], weights=w):.1f} deg")
        report()
        report("Single paddles: two-end hit rate against muons expected through the paddle (any direction)")
        singles = {}
        for p in setup.paddles:
            exp_r = mc.singles_rate(p, setup)
            meas = einfo.get("paddles", {}).get(p.name, {}).get("AB_rate_hz", float("nan"))
            singles[p.name] = {"expected_hz": exp_r, "measured_hz": meas}
            report(f"  {p.name:<10}{meas:>8.3f} Hz measured, {exp_r:.3f} Hz muons expected, ratio {meas / exp_r:.2f}")
        report("  (a ratio above 1: other particles (electrons, gammas, radioactivity) or noise also give two-end hits)")
        expected = {"pairs": crossing, "total_hz": total_rate, "ratio": ratio, "ratio_err": err,
                    "flux_if_efficient": flux, "singles": singles,
                    "toy_mean_theta_deg": float(np.average(th, weights=w)),
                    "toy_true_mean_theta_deg": float(np.average(toy["theta_true"], weights=w))}

    checks = map_check(cols, tcols if toy is not None else None, setup, calib)
    if checks:
        report()
        report("Channel-map check: the mean position along a layer's paddles for the tracks through each paddle of")
        report("another layer, data against the toy; fitted as data = shift + slope x toy (slope +1 = map as the")
        report("setup file says, -1 = one direction or order the other way round):")
        for r in checks:
            parts = ", ".join(f"{q['name']} {q['data_cm']:+.1f} ± {q['data_err_cm']:.1f} (toy {q['toy_cm']:+.1f})"
                              for q in r["paddles"])
            report(f"  {r['layer']} layer, {r['axis']} from {r['method']:<6} by {r['by_layer']} paddle: {parts} cm")
            report(f"      slope {r['slope']:+.2f} ± {r['slope_err']:.2f}, shift {r['shift_cm']:+.1f} ± "
                   f"{r['shift_err_cm']:.1f} cm -> {r['verdict']}")

    info = mc.stamp("reconstruct")
    info.update({"run": pinfo["run"], "span_ns": pinfo["span_ns"], "setup": {"name": setup.name, "path": setup.path},
                 "calibration": calib_path, "n_tracks": int(n), "expected": expected, "flux": setup.flux,
                 "layers": names, "map_check": checks})
    mc.write_info(tracks_path, info)
    with open(os.path.join(folder, "tracks_summary.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Saved {tracks_path} and tracks_summary.txt")
    return tracks_path


if __name__ == "__main__":
    main()
