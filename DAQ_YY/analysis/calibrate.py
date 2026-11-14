#!/usr/bin/env python3
"""
calibrate.py  -  step 3 of the paddle muon analysis

Measures the constants the reconstruction needs from the hits of the run itself
(cosmic muons calibrate the paddles): for every paddle the offset and the spread
of tA - tB, the signal speed along the paddle, the light attenuation and the
height of a muon's signal; for every channel its muon peak and its threshold in
those units; and the time offsets between the layers. Reads events.root (step 2)
and writes calibration.json and calibration_summary.txt.

Usage:
    python3 calibrate.py                       # asks for the analysis folder
    python3 calibrate.py /path/to/run00056_..._muon
    python3 calibrate.py FOLDER --setup setups/two_layer_2026-11-12.toml

Another run's calibration.json can be used instead of this one in the next step
(reconstruct.py --calibration), e.g. a long run's constants for a short run.

Method (two-end hits only; u = position along the paddle from its centre, cm,
positive towards the A end; L = the paddle's length; v = signal speed). Two
passes: the first uses every two-end hit and is only good enough to find the
muons and their signal (the muon peak, below). The second repeats the speed,
offset, resolution and attenuation fits with the hits whose sqrt(hA x hB) is
above speed_min_mip (default 0.8) x the muon peak: single paddles also see
low-energy hits (electrons, gammas, radioactivity), and those pass both
thresholds mostly near the centre, which makes the spread of tA - tB look
narrower (in a simulated run v came out 7 % high without the cut). sqrt(hA x hB)
does not depend on where the paddle was hit, so the cut keeps the muons even
along the paddle.
    tA - tB      = -2 u / v + offset. Muons cross a paddle evenly along its length,
                   so tA - tB is spread evenly over 2 L / v, smeared by the time
                   resolution: a box from lo to hi convolved with a Gaussian
                   (sigma), on a small flat background. An unbinned fit gives
                   v = 2 L / (hi - lo), offset = (lo + hi) / 2, and sigma; the
                   position resolution is v / 2 x sigma.
                   With fewer than min_hits hits (setup [calibration]) the
                   defaults stay: v_eff_cm_per_ns, tdiff_sigma_ns, offset = median.
    attenuation  light falls as exp(-d / lambda) over a distance d to an end, so
                   ln(hA / hB) = 2 u / lambda + ln(gA / gB): a straight line in the
                   true u. The measured u carries the position resolution, which
                   would flatten a plain straight-line fit (lambda 20 % too long in
                   a simulated run), so the line is fitted against the mean true u
                   given the measured one (muons even along the paddle, Gaussian
                   resolution: "regression calibration"; points beyond 3 rms
                   dropped). It gives lambda and the gain ratio of the two ends.
    muon peak    sqrt(hA x hB) does not depend on u. Single paddles also see
                   electrons, gammas and noise, so the peak is taken from the hits
                   of the tracks (events with one two-end hit in each of two
                   layers), as measured: no correction for the track's angle (which
                   would assume that the signal grows with the path length; the
                   plot signal_vs_path tests that on the data). Its most probable
                   value (MPV) comes from a fit of a Moyal function (close to the
                   Landau shape of energy loss) from 0.75 x the peak up. Per
                   channel: the signal of that end for a track through the paddle's
                   centre, MIP x sqrt(gA / gB) for A and MIP / sqrt(gA / gB) for B
                   (a fit of the hits within L / 4 of the centre if the attenuation
                   was not measured). With fewer than min_hits track hits, all
                   two-end hits are used (and the summary says so). The threshold
                   of each channel is then also given in units of its muon peak.
    layer times  t = (tA + tB) / 2 - L / (2 v) is when the muon crossed the paddle.
                   The paddles' cable and electronics delays are NOT taken out:
                   that would need an independent measurement (a pulser into every
                   channel), and setting them so that the tracks have beta = 1
                   would force the result. Only as text, the summary quotes the
                   offsets that beta = 1 would need (the median of t(lower) -
                   t(upper) - s / c per layer, or per paddle by least squares when
                   every paddle has min_tof_events events); they are not applied.
    timing       with pulse times from the rising-edge fit (step 1, timing = "fit") the
                   spread of tA - tB of the same hits is also fitted with the crossings
                   between two samples, to show what the fit gains (tdiff_sigma_cfd_ns).

Outputs (in the analysis folder):
    calibration.json          the constants, with how each was found
    calibration_summary.txt   the same as a table
"""

import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse
import math

import numpy as np

import muon_common as mc


def fit_range(p, setup):
    max_dt = p.length / float(setup.calibration["v_eff_cm_per_ns"]) + float(setup.events["pair_margin_ns"])
    return -max_dt, max_dt


def calibrate_paddle(p, hits, setup, shared=None):
    """Constants of one paddle from its two-end hits. shared = (v, v_err, n_paddles):
    the signal speed fitted for all paddles together; then only the offset, the
    resolution and the background are fitted here."""
    cs = setup.calibration
    v_nom = float(cs["v_eff_cm_per_ns"])
    r0, r1 = fit_range(p, setup)
    td = hits["tdiff_ns"]
    out = {"n_hits": int(len(td)), "length_cm": p.length}
    free = mc.box_fit(td, r0, r1) if len(td) >= int(cs["min_hits"]) else None
    if free is not None:
        width = free["hi"] - free["lo"]
        cov = free["cov"]
        var_w = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
        out["v_free_cm_per_ns"] = 2 * p.length / width
        out["v_free_err_cm_per_ns"] = out["v_free_cm_per_ns"] * math.sqrt(var_w) / width if var_w > 0 else float("nan")
        out["tdiff_sigma_free_ns"] = free["sigma"]
        out["tdiff_sigma_free_err_ns"] = float(free["err"][2])
    if shared is not None:
        fit = mc.box_fit(td, r0, r1, half_width=p.length / shared[0]) if len(td) >= 10 else None
        how = f"fitted, signal speed shared by {shared[2]} paddles"
    else:
        fit, how = free, "fitted"
    ok = False
    if fit is not None:
        width = fit["hi"] - fit["lo"]
        cov = fit["cov"]
        if shared is not None:
            v, v_err = shared[0], shared[1]
        else:
            var_w = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
            v = 2 * p.length / width
            v_err = v * math.sqrt(var_w) / width if var_w > 0 else float("nan")
        off_err = (math.sqrt(cov[0, 0]) if shared is not None else
                   0.5 * math.sqrt(cov[0, 0] + cov[1, 1] + 2 * cov[0, 1])) if np.isfinite(cov).all() else float("nan")
        ok = 5.0 <= v <= 30.0 and 0.02 <= fit["sigma"] < 0.5 * width
        if ok:
            out.update({
                "tdiff_status": how,
                "tdiff_offset_ns": 0.5 * (fit["lo"] + fit["hi"]), "tdiff_offset_err_ns": off_err,
                "v_eff_cm_per_ns": v, "v_eff_err_cm_per_ns": v_err,
                "tdiff_sigma_ns": fit["sigma"], "tdiff_sigma_err_ns": float(fit["err"][2]),
                "tdiff_lo_ns": fit["lo"], "tdiff_hi_ns": fit["hi"], "tdiff_background": fit["bkg"],
                "tdiff_fit_range_ns": [r0, r1], "tdiff_fit_n": fit["n"],
            })
    if not ok:
        max_dt = r1
        why = (f"only {len(td)} two-end hits (min_hits {cs['min_hits']})" if len(td) < int(cs["min_hits"])
               else "fit failed or gave an implausible speed")
        out.update({
            "tdiff_status": f"default: {why}",
            "tdiff_offset_ns": float(np.median(td)) if len(td) else 0.0, "tdiff_offset_err_ns": float("nan"),
            "v_eff_cm_per_ns": shared[0] if shared is not None else v_nom,
            "v_eff_err_cm_per_ns": shared[1] if shared is not None else float("nan"),
            "tdiff_sigma_ns": float(cs["tdiff_sigma_ns"]), "tdiff_sigma_err_ns": float("nan"),
            "tdiff_fit_range_ns": [-max_dt, max_dt],
        })
    v = out["v_eff_cm_per_ns"]
    out["position_sigma_cm"] = 0.5 * v * out["tdiff_sigma_ns"]
    out["mean_time_sigma_ns"] = 0.5 * out["tdiff_sigma_ns"]
    u = -0.5 * v * (td - out["tdiff_offset_ns"])

    # attenuation: ln(hA / hB) against the expected true position given the measured one
    # ("regression calibration": a plain fit against the measured u would flatten the slope,
    # because u carries the position resolution)
    hA, hB = hits["hA"], hits["hB"]
    clean = ((hits["flagsA"] | hits["flagsB"]) & (mc.FLAG_SATURATED | mc.FLAG_PILEUP)) == 0
    sel = clean & (hA > 0) & (hB > 0)
    u_exp = mc.expected_position(u, 0.5 * p.length, out["position_sigma_cm"])
    line = mc.robust_line(u_exp[sel], np.log(hA[sel] / hB[sel])) if sel.sum() >= int(cs["min_hits"]) else None
    if line is not None and line["b"] > 0:
        lam = 2.0 / line["b"]
        out.update({"attenuation_status": "fitted", "attenuation_cm": lam,
                    "attenuation_err_cm": lam * line["b_err"] / line["b"],
                    "log_ratio_slope_per_cm": line["b"], "log_ratio_slope_err": line["b_err"],
                    "log_ratio_offset": line["a"], "gain_ratio_A_over_B": math.exp(line["a"]),
                    "log_ratio_rms": line["rms"], "attenuation_n": line["n"]})
    else:
        why = (f"only {sel.sum()} usable hits" if line is None else
               f"ln(hA/hB) does not rise with u (slope {line['b']:.4f} +- {line['b_err']:.4f} per cm)")
        out.update({"attenuation_status": f"not measured: {why}", "attenuation_cm": None,
                    "gain_ratio_A_over_B": math.exp(line["a"]) if line is not None else None})
        if line is not None:
            out.update({"log_ratio_slope_per_cm": line["b"], "log_ratio_slope_err": line["b_err"],
                        "log_ratio_offset": line["a"], "log_ratio_rms": line["rms"]})

    return out


def muon_peaks(H, E, setup, cal, thresholds):
    """Muon peak (MPV) of every paddle and channel, from the hits of muons that cross
    two layers or more, each times cos(zenith) of its track (the energy loss per
    vertical path). A paddle with fewer than min_hits such hits uses all its
    two-end hits instead (single paddles also see electrons, gammas and noise).
    Each end's height is moved to the paddle's centre with the attenuation length,
    so the channel's MPV is the signal of a vertical muon through the centre.
    Fills the mip_* entries of cal; returns the per-channel table."""
    min_hits = int(setup.calibration["min_hits"])
    events, x, y, z, u, t = muon_sample(H, E, setup, cal)
    cos_of = {}
    for e, idx in events.items():
        i, j = idx[0], idx[-1]
        d = math.dist((x[i], y[i], z[i]), (x[j], y[j], z[j]))
        for k in idx:
            cos_of[k] = abs(z[i] - z[j]) / d if d > 0 else 1.0
    channels = {}
    bad = mc.FLAG_SATURATED | mc.FLAG_PILEUP
    for p in setup.paddles:
        if not p.both_ends:
            continue
        c = cal[p.name]
        lam = c.get("attenuation_cm")
        muon = np.array(sorted(k for k in cos_of if H["paddle"][k] == p.index
                               and ((H["flagsA"][k] | H["flagsB"][k]) & bad) == 0), dtype=np.int64)
        if len(muon) >= min_hits:
            idx, cos = muon, np.ones(len(muon))      # as measured: no path-length correction
            sample = f"{len(muon)} hits of tracks (two layers), as measured"
            c["mip_sample"] = "muons"
        else:
            idx = np.where((H["paddle"] == p.index) & (H["kind"] == mc.HIT_AB) &
                           (((H["flagsA"] | H["flagsB"]) & bad) == 0))[0]
            cos = np.ones(len(idx))
            sample = (f"all {len(idx)} two-end hits (only {len(muon)} hits of muons crossing two layers, "
                      f"min_hits {min_hits}): not only muons")
            c["mip_sample"] = "all"
        hA, hB = H["hA"][idx] * cos, H["hB"][idx] * cos
        uu = u[idx]
        geo = np.sqrt(hA * hB)
        fit = mc.moyal_fit(geo)
        if fit is not None:
            c.update({"mip_status": "fitted", "mip_adc": fit["mpv"], "mip_err_adc": fit["mpv_err"],
                      "mip_width_adc": fit["width"], "mip_fit_range_adc": list(fit["range"]), "mip_n": fit["n"]})
        else:
            c.update({"mip_status": "median (too few hits for a fit)",
                      "mip_adc": float(np.median(geo)) if len(geo) else None})
        c["mip_from"] = sample
        gain = c.get("gain_ratio_A_over_B")
        for end, h in (("A", hA), ("B", hB)):
            chid = p.channels[end]
            thr = thresholds.get(chid)
            ch = {"paddle": p.name, "end": end, "threshold_adc": thr, "n_hits": int(len(h))}
            if lam and gain and c.get("mip_adc"):
                # sqrt(hA hB) = sqrt(gA gB) E exp(-L / 2 lambda), and at the centre hA = gA E exp(-L / 2 lambda):
                # so the end's signal for a muon through the centre is MIP x sqrt(gA / gB) (A) or / sqrt (B)
                f = math.sqrt(gain) if end == "A" else 1.0 / math.sqrt(gain)
                rel = (c.get("mip_err_adc") or 0.0) / c["mip_adc"]
                ch.update({"mip_status": "from the paddle's muon peak and the gain ratio of its ends",
                           "mip_adc": c["mip_adc"] * f, "mip_err_adc": c["mip_adc"] * f * rel,
                           "mip_from": f"{sample}; MIP x sqrt(gA/gB)^{'+' if end == 'A' else '-'}1"})
            else:
                centre = h[np.abs(uu) < 0.25 * p.length]
                fe = mc.moyal_fit(centre)
                ch["mip_from"] = f"{sample}; hits within L/4 of the paddle centre"
                if fe is not None:
                    ch.update({"mip_status": "fitted", "mip_adc": fe["mpv"], "mip_err_adc": fe["mpv_err"],
                               "mip_width_adc": fe["width"]})
                else:
                    ch.update({"mip_status": "median (too few hits for a fit)",
                               "mip_adc": float(np.median(centre)) if len(centre) else None})
            if thr is not None and ch.get("mip_adc"):
                ch["threshold_mip"] = thr / ch["mip_adc"]
            channels[str(chid)] = ch
    return channels


def muon_sample(H, E, setup, cal):
    """Events with one two-end hit in each of two layers or more (a muon crossing
    the setup): {event: hit indices from the highest to the lowest}, and the
    points (x, y, z), u and crossing times of all hits with this calibration."""
    x, y, z, u, t = hit_points(H, setup, cal)
    single = np.ones(len(E["event"]), dtype=bool)
    present = np.zeros(len(E["event"]), dtype=np.int32)
    for layer in setup.layers:
        n_ab = E[f"n_ab_{layer.name}"]
        single &= n_ab <= 1
        present += (n_ab == 1)
    good = single & (present >= 2)
    events = {}
    for i in np.where(H["kind"] == mc.HIT_AB)[0]:
        e = int(H["event"][i])
        if good[e]:
            events.setdefault(e, []).append(int(i))
    for e in events:
        events[e].sort(key=lambda i: -z[i])
    return events, x, y, z, u, t


def hit_points(H, setup, cal):
    """Position (x, y, z), u and crossing time of every two-end hit, with this calibration."""
    n = len(H["t_ns"])
    x, y, z, u, t = (np.full(n, np.nan) for _ in range(5))
    for p in setup.paddles:
        sel = H["paddle"] == p.index
        if not sel.any() or p.name not in cal:
            continue
        c = cal[p.name]
        v = c["v_eff_cm_per_ns"]
        u_meas = -0.5 * v * (H["tdiff_ns"][sel] - c["tdiff_offset_ns"])
        # the same track point as reconstruct.py (setup [tracks] points)
        u[sel], _ = mc.track_point(u_meas, p, c.get("position_sigma_cm", 0.0), setup)
        x[sel], y[sel], z[sel] = p.point(u[sel])
        t[sel] = H["t_ns"][sel] - 0.5 * p.length / v
    return x, y, z, u, t


def tof_offsets(H, E, setup, cal):
    """Time offsets of layers (or paddles) from events with one two-end hit per layer."""
    cs = setup.calibration
    by_event, x, y, z, u, t = muon_sample(H, E, setup, cal)
    layer_names = [layer.name for layer in setup.layers]
    rows = []        # (upper paddle, lower paddle, residual)
    for e, idx in by_event.items():
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                i, j = idx[a], idx[b]
                s = math.dist((x[i], y[i], z[i]), (x[j], y[j], z[j]))
                rows.append((int(H["paddle"][i]), int(H["paddle"][j]), t[j] - t[i] - s / mc.C_CM_PER_NS))
    out = {"n_events": len(by_event), "n_pairs": len(rows), "assumption": "beta = 1 (quoted only, not applied)",
           "layer_offsets_ns": {name: 0.0 for name in layer_names},
           "paddle_offsets_ns": {p.name: 0.0 for p in setup.paddles}}
    if not rows:
        out["status"] = "default: no events with one two-end hit in each of two layers"
        return out
    up = np.array([r[0] for r in rows])
    low = np.array([r[1] for r in rows])
    res = np.array([r[2] for r in rows])
    top = setup.layers[0]
    out["residual_raw_median_ns"] = float(np.median(res))
    out["residual_raw_sigma_ns"] = float(1.4826 * np.median(np.abs(res - np.median(res))))
    # layer offsets against the top layer
    lay_up = np.array([setup.paddles[k].layer.index for k in up])
    lay_low = np.array([setup.paddles[k].layer.index for k in low])
    n_layer = {}
    for layer in setup.layers[1:]:
        sel = (lay_up == top.index) & (lay_low == layer.index)
        n_layer[layer.name] = int(sel.sum())
        if sel.any():
            out["layer_offsets_ns"][layer.name] = float(np.median(res[sel]))
    out["layer_pairs"] = n_layer
    for p in setup.paddles:
        out["paddle_offsets_ns"][p.name] = out["layer_offsets_ns"][p.layer.name]
    counts = {p.index: int(((up == p.index) | (low == p.index)).sum()) for p in setup.paddles}
    involved = [p for p in setup.paddles if counts[p.index] > 0]
    if involved and all(counts[p.index] >= int(cs["min_tof_events"]) for p in involved):
        # least squares: res = o(lower) - o(upper); the mean over the top layer's paddles is 0
        index = {p.index: k for k, p in enumerate(involved)}
        A = np.zeros((len(res) + 1, len(involved)))
        for r, (i, j) in enumerate(zip(up, low)):
            A[r, index[j]] += 1.0
            A[r, index[i]] -= 1.0
        tops = [index[p.index] for p in involved if p.layer is top]
        A[-1, tops] = 1.0
        b = np.append(res, 0.0)
        keep = np.ones(len(b), dtype=bool)
        for _ in range(3):
            sol, *_ = np.linalg.lstsq(A[keep], b[keep], rcond=None)
            r = b - A @ sol
            rms = math.sqrt(np.mean(r[:-1][keep[:-1]] ** 2)) if keep[:-1].any() else 0.0
            new = np.append(np.abs(r[:-1]) <= 3 * max(rms, 1e-6), True)
            if np.array_equal(new, keep):
                break
            keep = new
        for p in involved:
            out["paddle_offsets_ns"][p.name] = float(sol[index[p.index]])
        out["status"] = f"per paddle (least squares, {int(keep[:-1].sum())} of {len(res)} pairs)"
    else:
        need = int(cs["min_tof_events"])
        out["status"] = (f"per layer (median against the top layer); per-paddle offsets need {need} events "
                         f"per paddle, fewest here {min(counts[p.index] for p in involved) if involved else 0}")
    corrected = res - np.array([out["paddle_offsets_ns"][setup.paddles[j].name] -
                                out["paddle_offsets_ns"][setup.paddles[i].name] for i, j in zip(up, low)])
    out["residual_sigma_ns"] = float(np.std(corrected)) if len(corrected) > 1 else float("nan")
    out["residual_robust_sigma_ns"] = float(1.4826 * np.median(np.abs(corrected - np.median(corrected))))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", default=None,
                    help="Analysis folder with events.root (asked for if left out)")
    ap.add_argument("--setup", default=None, help="Setup file (default: the one step 1 used)")
    args = ap.parse_args()
    folder = mc.ask_analysis_folder(args.folder, mc.EVENTS_FILE, "build_events.py")
    run(folder, args.setup)


def run(folder, setup_path=None):
    """Step 3 for one analysis folder; returns the path of calibration.json."""
    events_path = mc.need_file(folder, mc.EVENTS_FILE, "build_events.py")
    pinfo = mc.read_info(mc.need_file(folder, mc.PULSES_FILE, "extract_pulses.py"))
    setup = mc.setup_of(folder, setup_path)
    H = mc.read_tree(events_path, "hits")
    E = mc.read_tree(events_path, "events")
    thresholds = {int(k): v.get("threshold") for k, v in pinfo.get("channels", {}).items()}
    fit_timing = pinfo.get("pulse_settings", {}).get("timing") == "fit"
    if fit_timing:   # the same hits timed with the crossing instead of the edge fit, for the comparison
        P = mc.read_tree(os.path.join(folder, mc.PULSES_FILE), "pulses", ["edge", "edge_cfd"])
        shift = (P["edge_cfd"] - P["edge"]) * mc.NS_PER_SAMPLE
        ia, ib = H["pulseA"], H["pulseB"]
        H["tdiff_cfd_ns"] = np.where((ia >= 0) & (ib >= 0), H["tdiff_ns"] + shift[np.maximum(ia, 0)]
                                     - shift[np.maximum(ib, 0)], np.nan)

    lines = []

    def report(text=""):
        print(text)
        lines.append(text)

    report(f"Analysis folder : {folder}")
    report(f"Run             : {pinfo['run']}")
    report(f"Setup           : {setup.name}")
    mode = str(setup.calibration["v_eff_mode"]).lower()
    if mode not in ("shared", "per_paddle"):
        print(f"WARNING: v_eff_mode '{mode}' unknown (shared or per_paddle); per_paddle used")
        mode = "per_paddle"
    cut = float(setup.calibration["speed_min_mip"])
    all_hits = {}
    for p in setup.paddles:
        sel = (H["paddle"] == p.index) & (H["kind"] == mc.HIT_AB)
        all_hits[p.name] = {k: H[k][sel] for k in ("tdiff_ns", "hA", "hB", "flagsA", "flagsB", "tdiff_cfd_ns")
                            if k in H}

    def one_pass(hits, label):
        """Speed (shared or not), offsets, resolution and attenuation from these hits."""
        shared, shared_fit = None, None
        if mode == "shared":
            samples = [(hits[p.name]["tdiff_ns"], p.length, *fit_range(p, setup)) for p in setup.paddles
                       if p.both_ends and len(hits[p.name]["tdiff_ns"]) >= int(setup.calibration["min_hits"])]
            if samples:
                print(f"{label}: one signal speed for {len(samples)} paddles (profile likelihood) ...")
                shared_fit = mc.shared_speed_fit(samples)
                if shared_fit is not None and 5.0 <= shared_fit["v"] <= 30.0:
                    shared = (shared_fit["v"], shared_fit["err"], len(samples))
                else:
                    print("WARNING: the shared signal-speed fit failed; every paddle is fitted on its own")
        out = {}
        for p in setup.paddles:
            if p.both_ends:
                out[p.name] = calibrate_paddle(p, hits[p.name], setup, shared)
            else:
                out[p.name] = {"n_hits": 0, "tdiff_status": "default: one end only", "tdiff_offset_ns": 0.0,
                               "v_eff_cm_per_ns": float(setup.calibration["v_eff_cm_per_ns"]),
                               "tdiff_sigma_ns": float(setup.calibration["tdiff_sigma_ns"]),
                               "position_sigma_cm": p.length / math.sqrt(12), "mean_time_sigma_ns": float("nan"),
                               "attenuation_status": "not measured: one end only", "attenuation_cm": None,
                               "mip_status": "not measured", "mip_adc": None}
        return out, shared, shared_fit

    bar = mc.Progress("Calibrating", 4)
    # pass 1: every two-end hit -> a first speed, enough to find the muons and their peak
    first, _, _ = one_pass(all_hits, "Pass 1")
    bar.update(1, "first pass")
    muon_peaks(H, E, setup, first, thresholds)
    bar.update(1, "muon peaks")
    # pass 2: only hits with a muon-like signal. sqrt(hA hB) does not depend on where the paddle
    # was hit, so this cut removes low-energy background without making the muons uneven along it
    hits, cut_note = {}, {}
    for p in setup.paddles:
        h = all_hits[p.name]
        mip = first[p.name].get("mip_adc")
        if cut > 0 and mip and first[p.name].get("mip_sample") == "muons":
            keep = np.sqrt(h["hA"] * h["hB"]) > cut * mip
            hits[p.name] = {k: v[keep] for k, v in h.items()}
            cut_note[p.name] = f"sqrt(hA hB) > {cut:g} x {mip:.0f} ADC: {int(keep.sum())} of {len(keep)} hits"
            first[p.name]["speed_cut_adc"] = cut * mip
        else:
            hits[p.name] = h
            cut_note[p.name] = ("all two-end hits (no muon peak from tracks for the speed_min_mip cut)"
                                if cut > 0 else "all two-end hits (speed_min_mip = 0)")
    paddles, shared, shared_fit = one_pass(hits, "Pass 2")
    for p in setup.paddles:
        paddles[p.name]["speed_hits"] = cut_note[p.name]
        paddles[p.name]["speed_cut_adc"] = first[p.name].get("speed_cut_adc")
        td_cfd = hits[p.name].get("tdiff_cfd_ns")
        if (td_cfd is not None and "tdiff_sigma_free_ns" in paddles[p.name]
                and len(td_cfd) >= int(setup.calibration["min_hits"])):
            fc = mc.box_fit(td_cfd, *fit_range(p, setup))
            if fc is not None:
                paddles[p.name]["tdiff_sigma_cfd_ns"] = fc["sigma"]
                paddles[p.name]["tdiff_sigma_cfd_err_ns"] = float(fc["err"][2])
    bar.update(1, "second pass")
    tof = tof_offsets(H, E, setup, paddles)
    channels = muon_peaks(H, E, setup, paddles, thresholds)
    bar.update(1, "layer times, muon peaks")
    bar.close()

    report()
    calib_warning = None
    if shared is not None:
        report(f"Signal speed    : one speed for {shared[2]} paddles (setup v_eff_mode = shared): "
               f"v = {shared[0]:.2f} +{shared_fit['err_up']:.2f} -{shared_fit['err_down']:.2f} cm/ns")
        free = [(p.name, paddles[p.name].get("v_free_cm_per_ns"), paddles[p.name].get("v_free_err_cm_per_ns"))
                for p in setup.paddles if paddles[p.name].get("v_free_cm_per_ns")]
        if free:
            report("  each paddle on its own: " + ", ".join(f"{n} {v:.1f} +- {e:.1f}" for n, v, e in free))
            chi2 = sum(((v - shared[0]) / e) ** 2 for _, v, e in free if e and np.isfinite(e))
            from scipy.stats import chi2 as chi2_dist
            p_value = float(chi2_dist.sf(chi2, max(len(free) - 1, 1)))
            report(f"  agreement with the shared speed: chi2 = {chi2:.1f} for {len(free)} paddles "
                   f"(probability {p_value:.2g})")
            if p_value < 0.01:
                thr = {}
                for p in setup.paddles:
                    fr = [channels[str(c)].get("threshold_mip") for c in p.channels.values() if str(c) in channels]
                    fr = [f for f in fr if f]
                    if fr:
                        thr[p.name] = max(fr)
                report("  WARNING: the paddles disagree. A likely cause: thresholds near the muon signal. A muon near")
                report("  one end gives the far end only exp(-L/2/lambda) of its centre signal, so near the ends the")
                report("  far end often stays below threshold; the spread of tA - tB then looks narrower and v comes")
                report("  out too high, the more the higher the threshold. Highest threshold / muon peak per paddle: " +
                       ", ".join(f"{n} {t:.2f}" for n, t in thr.items()))
                report("  Positions near the paddle ends can then be off by about (dv / v) x L / 2. Lower thresholds")
                report("  remove the effect; v_eff_mode = \"per_paddle\" in the setup file shows the paddles' own values.")
                calib_warning = (f"signal speeds of the paddles disagree (chi2 {chi2:.0f} for {len(free)}): "
                                 "see calibration_summary.txt")
    else:
        report(f"Signal speed    : fitted for each paddle on its own (setup v_eff_mode = {mode})")
    report()
    report("Per paddle: tA - tB fit (offset, spread sigma, signal speed), position resolution, attenuation, muon peak")
    head = (f"{'paddle':<9}{'hits':>6}{'offset ns':>11}{'sigma ns':>10}{'v cm/ns':>14}{'pos cm':>8}"
            f"{'lambda cm':>11}{'gA/gB':>7}{'MIP ADC':>13}  status")
    report(head)
    report("-" * (len(head) + 20))
    for p in setup.paddles:
        c = paddles[p.name]
        v_txt = f"{c['v_eff_cm_per_ns']:.1f}" + (f"+-{c['v_eff_err_cm_per_ns']:.1f}"
                                                 if np.isfinite(c.get("v_eff_err_cm_per_ns", float("nan"))) else "")
        lam = c.get("attenuation_cm")
        mip = c.get("mip_adc")
        mip_txt = (f"{mip:.0f}" + (f"+-{c['mip_err_adc']:.0f}" if np.isfinite(c.get("mip_err_adc", float("nan")))
                                   else "")) if mip else "-"
        gain = c.get("gain_ratio_A_over_B")
        status = "; ".join(s for s in (c["tdiff_status"] if c["tdiff_status"] != "fitted" else "",
                                       c["attenuation_status"] if c["attenuation_status"] != "fitted" else "",
                                       c["mip_status"] if c["mip_status"] != "fitted" else "") if s) or "all fitted"
        report(f"{p.name:<9}{c['n_hits']:>6}{c['tdiff_offset_ns']:>11.3f}{c['tdiff_sigma_ns']:>10.3f}{v_txt:>14}"
               f"{c['position_sigma_cm']:>8.2f}{(f'{lam:.0f}' if lam else '-'):>11}"
               f"{(f'{gain:.2f}' if gain else '-'):>7}{mip_txt:>13}  {status}")
    report()
    report("Speed, offset, resolution and attenuation fitted on:")
    for p in setup.paddles:
        report(f"  {p.name:<8}{paddles[p.name].get('speed_hits', '-')}")
    if any("tdiff_sigma_cfd_ns" in paddles[p.name] for p in setup.paddles):
        report()
        report("Pulse times from the rising-edge fit (step 1). The same hits with the crossing between two samples")
        report("instead (free fits of both; the two share every hit, so their difference is better known than the")
        report("errors suggest):")
        report(f"  {'paddle':<9}{'sigma(tA-tB) fit':>18}{'crossing':>12}{'difference':>12}")
        for p in setup.paddles:
            c = paddles[p.name]
            if "tdiff_sigma_cfd_ns" in c:
                a, b = c["tdiff_sigma_free_ns"], c["tdiff_sigma_cfd_ns"]
                report(f"  {p.name:<9}{a:>15.3f} ns{b:>9.3f} ns{100 * (a - b) / b:>+11.1f} %")
    report()
    report("Muon peak (MIP) taken from:")
    for p in setup.paddles:
        if paddles[p.name].get("mip_from"):
            report(f"  {p.name:<8}{paddles[p.name]['mip_from']}")
    report()
    report("Per channel: muon peak of each end's height for a vertical muon through the paddle's centre, and the "
           "DAQ threshold in those units")
    head = f"{'channel':<10}{'end':<9}{'hits':>6}{'MIP ADC':>12}{'threshold':>11}{'thr/MIP':>9}"
    report(head)
    report("-" * len(head))
    for chid in sorted(channels, key=int):
        c = channels[chid]
        mip = c.get("mip_adc")
        mip_txt = (f"{mip:.0f}" + (f"+-{c['mip_err_adc']:.0f}" if np.isfinite(c.get("mip_err_adc", float("nan")))
                                   else "")) if mip else "-"
        thr = c.get("threshold_adc")
        thr_txt = "-" if thr is None else str(thr)
        frac_txt = f"{c['threshold_mip']:.2f}" if "threshold_mip" in c else "-"
        end_txt = f"{c['paddle']}-{c['end']}"
        report(f"ch{chid:<8}{end_txt:<9}{c['n_hits']:>6}{mip_txt:>12}{thr_txt:>11}{frac_txt:>9}")
    report()
    report("Layer times: no time offsets applied (the paddles' cable delays would need a pulser measurement).")
    report(f"  {tof['n_events']} events with one two-end hit in each of two layers or more ({tof['n_pairs']} layer pairs)")
    if tof["n_pairs"]:
        report(f"  measured t(lower) - t(upper) - s/c: median {tof['residual_raw_median_ns']:+.3f} ns, "
               f"robust sigma {tof['residual_raw_sigma_ns']:.3f} ns")
        report(f"  only as text, not applied: offsets that would give the tracks beta = 1 ({tof['status']}):")
        for name, value in tof["layer_offsets_ns"].items():
            report(f"    layer {name:<10} {value:+.3f} ns" + ("  (reference)" if name == setup.layers[0].name else ""))
        if tof["status"].startswith("per paddle"):
            for name, value in tof["paddle_offsets_ns"].items():
                report(f"    paddle {name:<9} {value:+.3f} ns")
    tof["applied"] = False
    for p in setup.paddles:
        paddles[p.name]["time_offset_ns"] = 0.0
        paddles[p.name]["beta1_offset_ns"] = tof["paddle_offsets_ns"][p.name]

    report()
    report("offset = centre of tA - tB (cable and electronics difference of the two ends); sigma = resolution of")
    report("tA - tB (one end: about sigma / sqrt 2; (tA + tB) / 2: sigma / 2); pos = position resolution along")
    report("the paddle, v / 2 x sigma; lambda = attenuation length; gA/gB = gain ratio of the ends; MIP = most")
    report("probable sqrt(hA hB) of a muon. thr/MIP = threshold as a fraction of the channel's muon peak: muons")
    report("crossing far from an end give less light there, so a threshold near 1 MIP loses many of them.")

    calib = {"info": {**mc.stamp("calibrate"), "run": pinfo["run"], "analysis_folder": folder,
                      "setup": {"name": setup.name, "path": setup.path}},
             "speed": ({"mode": "shared", "v_eff_cm_per_ns": shared[0], "err": shared[1],
                        "err_up": shared_fit["err_up"], "err_down": shared_fit["err_down"], "paddles": shared[2]}
                       if shared is not None else {"mode": mode}),
             "paddles": paddles, "channels": channels, "tof": tof,
             "warnings": [calib_warning] if calib_warning else []}
    path = os.path.join(folder, mc.CALIB_FILE)
    mc.write_json(path, calib)
    with open(os.path.join(folder, "calibration_summary.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Saved {path} and calibration_summary.txt")
    return path


if __name__ == "__main__":
    main()
