#!/usr/bin/env python3
"""
check_fake_run.py  -  compare the analysis of a simulated run with its truth

Reads the analysis folder of a run made by make_fake_run.py and the run's
_truth.json / _truth.root (found next to the raw file recorded in pulses.root),
and prints, for every constant and every track quantity, the true value, the
analysis result and their difference in units of the stated error.

Usage:
    python3 tests/check_fake_run.py /path/to/fake_muon
"""

import os
import sys
import math

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import muon_common as mc  # noqa: E402
from make_fake_run import PEAK  # noqa: E402  (the pulse shape's maximum: the measured height is h x PEAK)


def line(name, true, got, err=None, unit=""):
    pull = (got - true) / err if err and np.isfinite(err) and err > 0 else float("nan")
    e = f" +- {err:.3g}" if err is not None and np.isfinite(err) else ""
    p = f"{pull:+6.1f} sigma" if np.isfinite(pull) else ""
    print(f"  {name:<34}{true:>10.3f}{got:>10.3f}{e:<12} {unit:<6}{p}")


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    folder = os.path.abspath(sys.argv[1])
    pinfo = mc.read_info(mc.need_file(folder, mc.PULSES_FILE, "run_analysis.py"))
    base = os.path.splitext(pinfo["raw_file"])[0]
    truth = mc.read_json(base + "_truth.json")
    mu = mc.read_tree(base + "_truth.root", "muons")
    cal = mc.read_json(mc.need_file(folder, mc.CALIB_FILE, "calibrate.py"))
    T = mc.read_tree(mc.need_file(folder, mc.TRACKS_FILE, "reconstruct.py"), "tracks")
    tinfo = mc.read_info(os.path.join(folder, mc.TRACKS_FILE))
    setup = mc.setup_of(folder)
    print(f"Simulated run: {pinfo['raw_file']}")
    print(f"  {'quantity':<34}{'truth':>10}{'analysis':>10}")
    sp = cal.get("speed", {})
    if sp.get("mode") == "shared":
        line("signal speed (shared)", truth["v_eff_cm_per_ns"], sp["v_eff_cm_per_ns"], sp["err"], "cm/ns")
    for p in setup.paddles:
        c, t = cal["paddles"][p.name], truth["paddles"][p.name]
        print(f" {p.name}:")
        if "v_free_cm_per_ns" in c:
            line("signal speed (this paddle alone)", truth["v_eff_cm_per_ns"], c["v_free_cm_per_ns"],
                 c.get("v_free_err_cm_per_ns"), "cm/ns")
        line("offset of tA - tB", t["tdiff_offset_ns"], c["tdiff_offset_ns"], c.get("tdiff_offset_err_ns"), "ns")
        line("sigma(tA - tB) >= sqrt2 x jitter", math.sqrt(2) * truth["jitter_ns"], c["tdiff_sigma_ns"],
             c.get("tdiff_sigma_err_ns"), "ns")
        if c.get("attenuation_cm"):
            line("attenuation length", truth["attenuation_cm"], c["attenuation_cm"], c.get("attenuation_err_cm"), "cm")
        if c.get("gain_ratio_A_over_B"):
            line("gain ratio A / B", t["gain_ratio_A_over_B"], c["gain_ratio_A_over_B"], None)
        if c.get("mip_adc"):
            line(f"muon peak ({c.get('mip_sample', '?')})", truth["mip_adc"] * PEAK * math.sqrt(t["gain_A"] * t["gain_B"]),
                 c["mip_adc"], c.get("mip_err_adc"), "ADC")

    # ---- tracks matched to true muons by time ----
    origin = pinfo["ttt_origin"] * mc.NS_PER_TICK
    tmu = mu["t_ns"] - origin
    names = [layer.name for layer in setup.layers]
    two = np.ones(len(tmu), dtype=bool)
    for nm in names:
        two &= mu[f"paddle_{nm}"] >= 0
    order = np.argsort(tmu)
    ts = tmu[order]
    j = np.clip(np.searchsorted(ts, T["t_ns"]), 1, len(ts) - 1)
    near = np.where(np.abs(ts[j - 1] - T["t_ns"]) < np.abs(ts[j] - T["t_ns"]), j - 1, j)
    dt = T["t_ns"] - ts[near]
    m = order[near]
    matched = np.abs(dt) < 50
    print()
    print(f"Tracks: {len(T['t_ns'])} reconstructed, {matched.sum()} matched to a simulated muon (|dt| < 50 ns); "
          f"{two.sum()} simulated muons crossed {' and '.join(names)} -> {matched.sum() / max(two.sum(), 1):.1%} "
          f"reconstructed (losses: thresholds)")
    exp = tinfo.get("expected", {})
    if exp:
        print(f"  toy expectation for the run: {exp['total_hz'] * truth['duration_s']:.1f} muons through two layers "
              f"(simulated: {two.sum()})")
    ok = matched & ((T["flags"] & mc.TRACK_ADJACENT) == 0)
    print(f"  {'residual (analysis - truth)':<34}{'mean':>10}{'rms':>10}")
    for nm in names:
        for axis in ("x", "y"):
            r = T[f"{axis}_{nm}"][ok] - mu[f"{axis}_{nm}"][m[ok]]
            r = r[np.isfinite(r)]
            if len(r):
                print(f"  {axis + ' at the ' + nm + ' layer (cm)':<34}{r.mean():>10.2f}{r.std():>10.2f}")
    r = T["theta_deg"][ok] - mu["theta_deg"][m[ok]]
    print(f"  {'zenith angle (deg)':<34}{r.mean():>10.2f}{r.std():>10.2f}")
    dphi = (T["phi_deg"][ok] - mu["phi_deg"][m[ok]] + 180) % 360 - 180
    big = T["theta_deg"][ok] > 15
    print(f"  {'azimuth, zenith > 15 deg (deg)':<34}{dphi[big].mean():>10.2f}{dphi[big].std():>10.2f}")
    ib = T["inv_beta"][ok]
    ib = ib[np.isfinite(ib)]
    print(f"  {'1/beta (truth 1)':<34}{ib.mean() - 1:>10.3f}{ib.std():>10.3f}")
    print()
    print("Pass if: the constants agree within about 2 sigma (sigma(tA - tB) is larger than sqrt2 x jitter: the ADC")
    print("noise adds to the timing), the position residuals have a mean near 0 and an rms")
    print("near the position resolution (across a paddle: width / sqrt 12 = 1.2 cm for 4 cm), the angles no bias.")


if __name__ == "__main__":
    main()
