#!/usr/bin/env python3
"""
make_plots.py  -  step 5 of the paddle muon analysis

Draws the plots of the analysis (one plot per file, in plots/) and writes
analysis_summary.txt, which starts with the main numbers of the run and then
repeats the summaries of steps 1-4. Reads everything the earlier steps wrote.

Usage:
    python3 make_plots.py                      # asks for the analysis folder
    python3 make_plots.py /path/to/run00056_..._muon

Plots (plots/<name>.png):
    pulse_rates            used and rejected pulses per second, per channel
    baseline_noise         baseline noise per channel (median, 68 % and 95 % ranges)
    pulse_shape_<paddle>   average pulse of both ends, normalised to its height
    heights_<paddle>       pulse heights of both ends, with the DAQ thresholds
    tdiff_<paddle>         tA - tB with the box fit (gives the signal speed)
    position_<paddle>      measured position along the paddle; why some lie beyond its ends
    attenuation_<paddle>   ln(hA / hB) against the position along the paddle
    tdiff_vs_ratio_<paddle> ln(hA / hB) against tA - tB, both raw: no calibration at all
    mip_<paddle>           sqrt(hA x hB) with the muon-peak fit
    layer_delay_<a>_<b>    time between two-end hits of two layers, +-1000 ns and a zoom of +-10 ns
    event_types            events by the layers with a two-end hit
    event_quality          selection steps: all events, a paddle at both ends, a track, one
                           paddle per layer or several paddles of one layer
    coincidences_<a>_<b>   events per paddle pair, measured and expected
    rates_vs_time          hit rate per layer and coincidence rate over the run
    hitmap_<layer>         where the tracks cross each layer (top view)
    zenith, azimuth        track angles against the toy expectation
    tof_residual           measured time of flight minus s / c (no time offsets applied)
    beta                   beta of the tracks over its whole range, shaded where no particle going
                           down can be, against the simulation of beta = 1 particles
    energy_deposit         the tracks' signal per paddle as measured (ADC)
    energy_deposit_sim     simulation: the energy muons lose in the paddles (MeV, known values)
    signal_vs_path         signal against the path through the paddle (1 / cos zenith): data test
"""

import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import muon_common as mc
from muon_common import INK, INK2, MUTED, GRID, AXIS, SURFACE, CATEGORICAL, END_COLORS, FIT_COLOR

DATA_COLOR = INK


class Context:
    """Everything the earlier steps wrote, loaded once."""

    def __init__(self, folder, setup_path=None):
        self.folder = folder
        self.plots = os.path.join(folder, mc.PLOTS_DIR)
        os.makedirs(self.plots, exist_ok=True)
        pp = mc.need_file(folder, mc.PULSES_FILE, "extract_pulses.py")
        ep = mc.need_file(folder, mc.EVENTS_FILE, "build_events.py")
        self.setup = mc.setup_of(folder, setup_path)
        self.pinfo = mc.read_info(pp)
        self.P = mc.read_tree(pp, "pulses")
        self.shapes = mc.read_tree(pp, "shapes")
        self.einfo = mc.read_info(ep)
        self.H = mc.read_tree(ep, "hits")
        self.E = mc.read_tree(ep, "events")
        tp = os.path.join(folder, mc.TRACKS_FILE)
        self.has_tracks = os.path.isfile(tp)
        self.Pts = None
        # the calibration the reconstruction used (maybe another run's), else this folder's
        cal_path = os.path.join(folder, mc.CALIB_FILE)
        if self.has_tracks:
            self.tinfo = mc.read_info(tp)
            self.Tr = mc.read_tree(tp, "tracks")
            self.toy = mc.read_tree(tp, "toy") if self.tinfo.get("expected") else None
            f = mc.open_root(tp)
            has_points = bool(f.Get("points"))
            f.Close()
            self.Pts = mc.read_tree(tp, "points") if has_points else None
            if self.tinfo.get("calibration") and os.path.isfile(self.tinfo["calibration"]):
                cal_path = self.tinfo["calibration"]
        if not os.path.isfile(cal_path):
            mc.need_file(folder, mc.CALIB_FILE, "calibrate.py")
        self.calib = mc.read_json(cal_path)
        self.calib_path = cal_path
        self.run = self.pinfo["run"]
        self.short = mc.short_run(self.run)
        self.span_s = self.pinfo["span_ns"] * 1e-9
        self.footer = f"{self.run}  ·  setup {self.setup.name}  ·  {self.span_s:.1f} s"
        self.made = []

    def save(self, fig, name):
        mc.footer(fig, self.footer)
        path = os.path.join(self.plots, f"{name}.png")
        mc.save_figure(fig, path)
        self.made.append(path)


# ---- drawing helpers ----

def note(ax, text, loc="upper right"):
    x, ha = (0.98, "right") if "right" in loc else (0.02, "left")
    y, va = (0.97, "top") if "upper" in loc else (0.03, "bottom")
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=7, color=INK2, linespacing=1.4,
            bbox=dict(boxstyle="round,pad=0.45", fc=SURFACE, ec=GRID, lw=0.6), zorder=5)


def stairs(ax, values, bins, color, label=None, weights=None, wash=True, lw=1.1):
    counts, edges = np.histogram(values, bins=bins, weights=weights)
    if wash:
        ax.stairs(counts, edges, fill=True, color=color, alpha=0.10, lw=0)
    ax.stairs(counts, edges, color=color, lw=lw, label=label)
    return counts, edges


def points(ax, values, bins, label="data"):
    """Counts as points with approximate Poisson 68 % error bars; empty bins are not drawn."""
    counts, edges = np.histogram(values, bins=bins)
    centres = 0.5 * (edges[1:] + edges[:-1])
    down, up = mc.poisson_band(counts)
    m = counts > 0
    ax.errorbar(centres[m], counts[m], yerr=[down[m], up[m]], fmt="o", ms=4, color=DATA_COLOR, mec=SURFACE,
                mew=1.0, elinewidth=0.9, capsize=0, label=label, zorder=4)
    return counts, edges


def figure(w=6.4, h=4.2):
    return plt.subplots(figsize=(w, h))


def fmt_err(value, err, digits=2):
    if err is None or not np.isfinite(err):
        return f"{value:.{digits}f}"
    return f"{value:.{digits}f} ± {err:.{digits}f}"


# ---- pulse plots ----

def plot_pulse_rates(c):
    chans = sorted(int(k) for k in c.pinfo["channels"])
    info = c.pinfo["channels"]
    used = np.array([info[str(ch)]["used"] for ch in chans]) / c.span_s
    rej = np.array([info[str(ch)]["records"] - info[str(ch)]["used"] for ch in chans]) / c.span_s
    labels = [f"ch{ch}  {c.setup.channel_label(ch)}" for ch in chans]
    fig, ax = figure(6.4, 0.32 * len(chans) + 1.3)
    y = np.arange(len(chans))[::-1]
    ax.barh(y, used, height=0.6, color=CATEGORICAL[0], edgecolor=SURFACE, lw=1.0, label="used")
    ax.barh(y, rej, left=used, height=0.6, color=AXIS, edgecolor=SURFACE, lw=1.0,
            label="rejected (noise, too low, too short, no edge)")
    for yi, u, r in zip(y, used, rej):
        ax.text(u + r, yi, f"  {u:.2f}" + (f" + {r:.2f}" if r > 0 else ""), va="center", fontsize=7, color=INK2)
    ax.set_yticks(y, labels)
    ax.set_xlabel("pulses per second")
    ax.set_xlim(0, max(used + rej) * 1.3 if len(used) else 1)
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right")
    mc.title(ax, f"Pulses per channel: {c.short}")
    c.save(fig, "pulse_rates")


def plot_baseline_noise(c):
    chans = sorted(int(k) for k in c.pinfo["channels"])
    fig, ax = figure(6.4, 0.32 * len(chans) + 1.3)
    y = np.arange(len(chans))[::-1]
    for yi, ch in zip(y, chans):
        r = c.P["baseline_rms"][c.P["channel"] == ch]
        if not len(r):
            continue
        q = np.percentile(r, [2.5, 16, 50, 84, 97.5])
        ax.plot([q[0], q[4]], [yi, yi], color=CATEGORICAL[0], lw=1.0, solid_capstyle="round")
        ax.plot([q[1], q[3]], [yi, yi], color=CATEGORICAL[0], lw=4.0, solid_capstyle="round", alpha=0.45)
        ax.plot(q[2], yi, "o", ms=5, color=CATEGORICAL[0], mec=SURFACE, mew=1.2)
        ax.text(q[4], yi, f"  {q[2]:.1f}", va="center", fontsize=7, color=INK2)
    limit = float(c.pinfo["pulse_settings"]["max_baseline_rms"])
    ax.axvline(limit, color=MUTED, lw=0.9, ls="--")
    ax.set_yticks(y, [f"ch{ch}  {c.setup.channel_label(ch)}" for ch in chans])
    ax.set_xlabel("baseline noise (rms of the baseline samples, ADC)")
    ax.set_xlim(0, None)
    ax.grid(axis="y", visible=False)
    mc.title(ax, f"Baseline noise per waveform: {c.short}",
             f"dot = median, thick = 68 % of the waveforms, thin = 95 %; dashed: NOISY flag ({limit:g} ADC); "
             f"the card alone gives ~2.2 ADC")
    c.save(fig, "baseline_noise")


def plot_pulse_shape(c, p):
    fig, ax = figure()
    drawn = 0
    for end, chid in p.channels.items():
        sel = c.shapes["channel"] == chid
        if not sel.any():
            continue
        t, m = c.shapes["t_ns"][sel], c.shapes["mean"][sel]
        ax.fill_between(t, c.shapes["p16"][sel], c.shapes["p84"][sel], color=END_COLORS[end], alpha=0.12, lw=0)
        ax.plot(t, m, color=END_COLORS[end], lw=1.2,
                label=f"{end} end, ch{chid} (mean of {int(c.shapes['n'][sel][0])} pulses)")
        drawn += 1
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.set_xlabel(f"time from the pulse time ({c.pinfo['pulse_settings']['cfd_fraction']:.0%} of the height), ns")
    ax.set_ylabel("signal / pulse height")
    mc.headroom(ax, 1.15)
    if drawn > 1:
        ax.legend(loc="upper right")
    mc.title(ax, f"Average pulse shape, {p.name}: {c.short}", "line = mean, band = 16-84 % of the pulses")
    c.save(fig, f"pulse_shape_{p.name}")


def plot_heights(c, p):
    fig, ax = figure()
    used = (c.P["flags"] & mc.REJECT_MASK) == 0
    allh = c.P["height"][used & np.isin(c.P["channel"], list(p.channels.values()))]
    if not len(allh):
        plt.close(fig)
        return
    top = np.percentile(allh, 99.5) * 1.15
    bins = np.linspace(0, top, 70)
    for end, chid in p.channels.items():
        h = c.P["height"][used & (c.P["channel"] == chid)]
        stairs(ax, np.clip(h, 0, top * 0.999), bins, END_COLORS[end], label=f"{end} end, ch{chid}: {len(h)} pulses")
        thr = c.pinfo["channels"].get(str(chid), {}).get("threshold")
        if thr is not None:
            ax.axvline(thr, color=END_COLORS[end], lw=0.9, ls="--")
    ax.set_yscale("log")
    ax.set_ylim(0.7, None)
    mc.headroom(ax, log=True)
    ax.set_xlabel("pulse height (ADC above the baseline; the last bin holds the higher ones)")
    ax.set_ylabel("pulses per bin")
    ax.legend(loc="upper right")
    mc.title(ax, f"Pulse heights, {p.name}: {c.short}", "dashed: DAQ threshold of each channel")
    c.save(fig, f"heights_{p.name}")


def fitted_hits(c, p, cal):
    """The paddle's two-end hits, and which of them the calibration's speed fit used."""
    sel = (c.H["paddle"] == p.index) & (c.H["kind"] == mc.HIT_AB)
    cut = cal.get("speed_cut_adc")
    geo = np.sqrt(c.H["hA"][sel] * c.H["hB"][sel])
    used = geo > cut if cut else np.ones(sel.sum(), dtype=bool)
    return sel, used, cut


def plot_tdiff(c, p):
    cal = c.calib["paddles"].get(p.name, {})
    sel, used, cut = fitted_hits(c, p, cal)
    td = c.H["tdiff_ns"][sel]
    r0, r1 = cal.get("tdiff_fit_range_ns", [-8, 8])
    fig, ax = figure()
    nb = int(np.clip(used.sum() / 4, 16, 60))
    bins = np.linspace(r0, r1, nb + 1)
    if cut:
        ax.stairs(np.histogram(td, bins=bins)[0], bins, color=AXIS, lw=1.1, label=f"all {len(td)} two-end hits")
    stairs(ax, td[used], bins, CATEGORICAL[0],
           label=(f"{used.sum()} with sqrt(hA hB) > {cut:.0f} ADC (fitted)" if cut else f"{len(td)} two-end hits"))
    if str(cal.get("tdiff_status", "")).startswith("fitted"):
        t = np.linspace(r0, r1, 600)
        width = bins[1] - bins[0]
        curve = cal["tdiff_fit_n"] * width * mc.box_pdf(t, cal["tdiff_lo_ns"], cal["tdiff_hi_ns"],
                                                        cal["tdiff_sigma_ns"], cal["tdiff_background"], r0, r1)
        ax.plot(t, curve, color=FIT_COLOR, lw=1.3, label="fit: even spread along the paddle ⊗ resolution")
        note(ax, f"signal speed v = {fmt_err(cal['v_eff_cm_per_ns'], cal.get('v_eff_err_cm_per_ns'), 1)} cm/ns\n"
                 f"offset = {fmt_err(cal['tdiff_offset_ns'], cal.get('tdiff_offset_err_ns'), 3)} ns\n"
                 f"resolution σ = {fmt_err(cal['tdiff_sigma_ns'], cal.get('tdiff_sigma_err_ns'), 3)} ns\n"
                 f"position σ = {cal['position_sigma_cm']:.1f} cm", loc="upper left")
        v, off = cal["v_eff_cm_per_ns"], cal["tdiff_offset_ns"]
        sec = ax.secondary_xaxis("top", functions=(lambda t, v=v, off=off: -0.5 * v * (t - off),
                                                   lambda u, v=v, off=off: off - 2 * u / v))
        sec.set_xlabel(f"position along the paddle u (cm, + towards the A end at {p.a_end})", color=INK2)
    else:
        note(ax, f"no fit: {cal.get('tdiff_status', '?')}", loc="upper left")
    ax.set_xlabel("tA - tB (ns)")
    ax.set_ylabel("hits per bin")
    ax.set_ylim(0, None)
    mc.headroom(ax)
    ax.legend(loc="upper right")
    mc.title(ax, f"Time difference of the two ends, {p.name}: {c.short}")
    c.save(fig, f"tdiff_{p.name}")


def plot_attenuation(c, p):
    """ln(hA / hB) against the position from the timing: light loses exp(-d / lambda) on its way to
    each end, so the log of the two ends' height ratio is a straight line in u of slope 2 / lambda."""
    cal = c.calib["paddles"].get(p.name, {})
    sel, used, cut = fitted_hits(c, p, cal)
    pos = (c.H["hA"][sel] > 0) & (c.H["hB"][sel] > 0)
    v, off = cal.get("v_eff_cm_per_ns", 15.0), cal.get("tdiff_offset_ns", 0.0)
    u = -0.5 * v * (c.H["tdiff_ns"][sel] - off)
    lr = np.full(len(u), np.nan)
    lr[pos] = np.log(c.H["hA"][sel][pos] / c.H["hB"][sel][pos])
    fig, ax = figure()
    many = len(u) > 3000
    ms = 2.2 if many else 3.5
    if cut:
        ax.plot(u[~used], lr[~used], "o", ms=ms, color=AXIS, mec=SURFACE, mew=0.4,
                label=f"sqrt(hA hB) < {cut:.0f} ADC (not fitted)")
    ax.plot(u[used], lr[used], "o", ms=ms, color=CATEGORICAL[0], mec=SURFACE, mew=0.4,
            label=f"{used.sum()} hits fitted")
    if "log_ratio_slope_per_cm" in cal:
        uu = np.linspace(-0.75 * p.length, 0.75 * p.length, 200)
        ue = mc.expected_position(uu, 0.5 * p.length, cal.get("position_sigma_cm", 0))
        ax.plot(uu, cal["log_ratio_offset"] + cal["log_ratio_slope_per_cm"] * ue, color=FIT_COLOR, lw=1.3,
                label="fitted line, as it looks at a measured u")
        lam = cal.get("attenuation_cm")
        text = (f"attenuation length λ = {fmt_err(lam, cal.get('attenuation_err_cm'), 1)} cm\n"
                f"(light halves every λ ln 2 = {lam * math.log(2):.0f} cm)\n" if lam else "attenuation not measured\n") + \
            f"gain ratio A/B = {math.exp(cal['log_ratio_offset']):.2f} (at u = 0)\nscatter (rms) = {cal['log_ratio_rms']:.2f}"
        note(ax, text, loc="upper left")
    for edge in (-0.5 * p.length, 0.5 * p.length):
        ax.axvline(edge, color=MUTED, lw=0.8, ls="--")
    ax.set_xlabel(f"position along the paddle u from the timing (tA - tB), cm, + towards the A end at {p.a_end}")
    ax.set_ylabel("ln(hA / hB): log of the two ends' height ratio")
    mc.headroom(ax, 1.5)
    sec = ax.secondary_yaxis("right", functions=(np.exp, lambda r: np.log(np.maximum(r, 1e-9))))
    lo, hi = ax.get_ylim()
    ticks = [t for t in (0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50) if lo <= math.log(t) <= hi]
    sec.set_yticks(ticks, [f"{t:g}" for t in ticks])
    sec.set_ylabel("hA / hB", color=INK2)
    ax.legend(loc="upper right")
    mc.title(ax, f"Light attenuation, {p.name}: {c.short}",
             "hA ∝ exp(-(L/2 - u)/λ), hB ∝ exp(-(L/2 + u)/λ), so ln(hA/hB) = 2u/λ + ln(gA/gB):\n"
             "the slope gives λ (the tA - tB edges give only the speed v); dashed: the paddle's ends")
    c.save(fig, f"attenuation_{p.name}")


def plot_tdiff_vs_ratio(c, p):
    """Light sharing against timing with no calibration at all: ln(hA / hB) against tA - tB of
    every two-end hit. If the light fades along the paddle, the end nearer the hit sees more
    light and sees it earlier, so the two raw measurements must follow each other."""
    from matplotlib.colors import LogNorm
    from scipy.stats import spearmanr
    sel = (c.H["paddle"] == p.index) & (c.H["kind"] == mc.HIT_AB) & (c.H["hA"] > 0) & (c.H["hB"] > 0)
    td = c.H["tdiff_ns"][sel]
    lr = np.log(c.H["hA"][sel] / c.H["hB"][sel])
    if len(td) < 20:
        return
    fig, ax = figure(6.8, 4.8)
    xr = np.percentile(td, [0.3, 99.7]) + np.array([-0.5, 0.5])
    yr = np.percentile(lr, [0.3, 99.7]) + np.array([-0.3, 0.3])
    nb = int(np.clip(math.sqrt(len(td)) / 2, 20, 70))
    H2, xe, ye = np.histogram2d(td, lr, bins=[np.linspace(*xr, nb + 1), np.linspace(*yr, nb + 1)])
    im = ax.pcolormesh(xe, ye, np.where(H2 > 0, H2, np.nan).T, cmap=mc.SEQUENTIAL,
                       norm=LogNorm(vmin=1, vmax=max(H2.max(), 2)), shading="flat", rasterized=True)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.12)
    cb.set_label("hits per cell", color=INK2)
    cb.outline.set_visible(False)
    # profile: median ln(hA/hB) in slices of tA - tB
    pe = np.linspace(*xr, 25)
    px, py, pyerr = [], [], []
    for a, b in zip(pe[:-1], pe[1:]):
        s = (td >= a) & (td < b)
        if s.sum() >= 20:
            med = np.median(lr[s])
            px.append(0.5 * (a + b))
            py.append(med)
            pyerr.append(1.2533 * 1.4826 * np.median(np.abs(lr[s] - med)) / math.sqrt(s.sum()))
    ax.errorbar(px, py, yerr=pyerr, fmt="o", ms=4, color=INK, mec=SURFACE, mew=1.0, elinewidth=0.9, capsize=0,
                label="median in slices of tA - tB", zorder=4)
    rho = spearmanr(td, lr)[0]
    r = np.corrcoef(td, lr)[0, 1]
    note(ax, f"{len(td)} two-end hits\nrank correlation ρ = {rho:+.3f}\nlinear correlation r = {r:+.3f}\n"
             f"(no calibration used)", loc="upper right")
    ax.set_xlim(*xr)
    ax.set_ylim(*yr)
    ax.grid(False)
    sec = ax.secondary_yaxis("right", functions=(np.exp, lambda q: np.log(np.maximum(q, 1e-9))))
    ticks = [t for t in (0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20) if yr[0] <= math.log(t) <= yr[1]]
    sec.set_yticks(ticks, [f"{t:g}" for t in ticks])
    sec.set_ylabel("hA / hB", color=INK2)
    ax.set_xlabel("tA - tB (ns), as measured")
    ax.set_ylabel("ln(hA / hB), as measured")
    ax.legend(loc="lower left")
    mc.title(ax, f"Light sharing against timing, {p.name}: {c.short}",
             "a hit nearer end A gives A more light (up) and gives it earlier (left): both raw measurements see the position")
    c.save(fig, f"tdiff_vs_ratio_{p.name}")


def plot_mip(c, p):
    """sqrt(hA hB) of all two-end hits and of the tracks' hits (as measured), with the
    fit the calibration made, every histogram as a fraction of its hits."""
    sel = (c.H["paddle"] == p.index) & (c.H["kind"] == mc.HIT_AB)
    geo_all = np.sqrt(c.H["hA"][sel] * c.H["hB"][sel])
    cal = c.calib["paddles"].get(p.name, {})
    if not len(geo_all):
        return
    geo_mu = np.zeros(0)
    if c.has_tracks:
        L = p.layer.name
        tsel = (c.Tr[f"paddle_{L}"] == p.index) & ((c.Tr["flags"] & (mc.TRACK_ADJACENT | mc.TRACK_AMBIGUOUS)) == 0)
        geo_mu = c.Tr[f"hgeo_{L}"][tsel]
        geo_mu = geo_mu[np.isfinite(geo_mu)]
    fig, ax = figure()
    top = np.percentile(geo_all, 99) * 1.2
    bins = np.linspace(0, top, 50 if len(geo_mu) >= 200 or not len(geo_mu) else 30)
    width = bins[1] - bins[0]
    ax.stairs(np.histogram(np.clip(geo_all, 0, top * 0.999), bins=bins)[0] / len(geo_all), bins,
              color=AXIS, lw=1.2, label=f"all two-end hits ({len(geo_all)}): also electrons, gammas, noise")
    fitted_on = cal.get("mip_sample", "all")
    if len(geo_mu):
        stairs(ax, np.clip(geo_mu, 0, top * 0.999), bins, CATEGORICAL[0], weights=np.full(len(geo_mu), 1 / len(geo_mu)),
               label=f"hits of tracks (two layers), as measured ({len(geo_mu)})")
    if cal.get("mip_status") == "fitted":
        lo, hi = cal["mip_fit_range_adc"]
        x = np.linspace(lo, hi, 400)
        norm = mc.moyal_cdf(hi, cal["mip_adc"], cal["mip_width_adc"]) - mc.moyal_cdf(lo, cal["mip_adc"], cal["mip_width_adc"])
        n_ref = len(geo_mu) if fitted_on == "muons" and len(geo_mu) else len(geo_all)
        ax.plot(x, cal["mip_n"] / n_ref * width * mc.moyal_pdf(x, cal["mip_adc"], cal["mip_width_adc"]) / norm,
                color=FIT_COLOR, lw=1.3,
                label=f"Moyal (Landau-like) fit to the {'tracks' if fitted_on == 'muons' else 'two-end hits'}, "
                      f"{lo:.0f}-{hi:.0f} ADC")
        note(ax, f"most probable signal (MPV) = {fmt_err(cal['mip_adc'], cal.get('mip_err_adc'), 0)} ADC\n"
                 f"width = {cal['mip_width_adc']:.0f} ADC\nfrom: {'the tracks' if fitted_on == 'muons' else 'all two-end hits'}",
             loc="upper left")
    ax.set_xlabel("sqrt(hA x hB), ADC (does not depend on where the paddle was crossed)")
    ax.set_ylabel("fraction of the hits per bin")
    ax.set_ylim(0, None)
    mc.headroom(ax, 1.45)
    ax.legend(loc="upper right")
    mc.title(ax, f"Signal of the paddle, {p.name}: {c.short}")
    c.save(fig, f"mip_{p.name}")


def plot_position(c, p):
    """Measured position along the paddle (from tA - tB): why some lie beyond the ends."""
    cal = c.calib["paddles"].get(p.name, {})
    if not str(cal.get("tdiff_status", "")).startswith("fitted"):
        return
    sel, used, cut = fitted_hits(c, p, cal)
    v, off, sig = cal["v_eff_cm_per_ns"], cal["tdiff_offset_ns"], cal["tdiff_sigma_ns"]
    u = -0.5 * v * (c.H["tdiff_ns"][sel][used] - off)
    s_u = 0.5 * v * sig
    half = 0.5 * p.length
    r0, r1 = cal["tdiff_fit_range_ns"]
    umax = 0.5 * v * max(abs(r0), abs(r1))
    fig, ax = figure()
    bins = np.linspace(-umax, umax, int(np.clip(len(u) / 40, 30, 90)) + 1)
    width = bins[1] - bins[0]
    stairs(ax, u, bins, CATEGORICAL[0], label=f"{len(u)} two-end hits" + (" (fitted)" if cut else ""))
    # the fitted model, an even spread over the paddle smeared by the resolution, in cm
    uu = np.linspace(-umax, umax, 800)
    tt = off - 2 * uu / v
    pdf_t = mc.box_pdf(tt, cal["tdiff_lo_ns"], cal["tdiff_hi_ns"], sig, cal["tdiff_background"], r0, r1)
    ax.plot(uu, len(u) * width * pdf_t * 2 / v, color=FIT_COLOR, lw=1.3,
            label=f"even along it, smeared by σ = {s_u:.1f} cm")
    for sgn in (-1, 1):
        ax.axvspan(sgn * half, sgn * umax, color=GRID, alpha=0.5, lw=0)
        ax.axvline(sgn * half, color=MUTED, lw=0.9, ls="--")
    # fractions beyond the ends: data, and the smearing alone (the fit without its flat part)
    grid = np.linspace(r0, r1, 20001)
    dens = mc.box_pdf(grid, cal["tdiff_lo_ns"], cal["tdiff_hi_ns"], sig, 0.0, r0, r1)
    ug = -0.5 * v * (grid - off)
    step = grid[1] - grid[0]
    exp_out = float(np.sum(dens[np.abs(ug) > half]) * step)
    exp_far = float(np.sum(dens[np.abs(ug) > half + 3 * s_u]) * step)
    out = float(np.mean(np.abs(u) > half))
    far = float(np.mean(np.abs(u) > half + 3 * s_u))
    note(ax, f"beyond the ends: {out:.1%}\n  smearing alone: {exp_out:.1%}\n"
             f"beyond the ends + 3σ: {far:.2%}\n  smearing alone: {exp_far:.3%}", loc="upper left")
    ax.set_xlim(-umax, umax)
    ax.set_ylim(0, None)
    mc.headroom(ax, 1.6)
    ax.set_xlabel(f"measured position along the paddle u (cm, + towards the A end at {p.a_end})")
    ax.set_ylabel("hits per bin")
    ax.legend(loc="upper right")
    mc.title(ax, f"Measured position along {p.name}: {c.short}",
             "grey: beyond the paddle's ends, reached only through the time resolution; tracks use the "
             "best estimate inside the paddle")
    c.save(fig, f"position_{p.name}")


def event_quality(c):
    """Selection steps of the events, each with its count: [(label, count, kind)]; kind
    'step' for the chain, 'clean' / 'several' for the split of the tracks."""
    E, s = c.E, c.setup
    names = [layer.name for layer in s.layers]
    n_ab = sum(E[f"n_ab_{nm}"].astype(np.int64) for nm in names)
    rows = [("all events", len(E["event"]), "step"),
            ("a paddle seen at both ends", int((n_ab >= 1).sum()), "step")]
    multi = E["n_layers_ab"] >= 2
    both = " and ".join(names) if len(names) == 2 else "two layers or more"
    rows.append((f"both ends in the {both} layers (a track)" if len(names) == 2 else
                 f"both ends in {both} (a track)", int(multi.sum()), "step"))
    if len(names) > 2:
        rows.append(("... in every layer", int((E["n_layers_ab"] == len(names)).sum()), "step"))
    one = multi.copy()
    for nm in names:
        one &= E[f"n_ab_{nm}"] <= 1
    flags = np.zeros(len(E["event"]), dtype=np.int64)
    if c.has_tracks and len(c.Tr["event"]):
        flags[c.Tr["event"]] = c.Tr["flags"]
    inside = one & ((flags & mc.TRACK_OUTSIDE) == 0)
    adjacent = multi & ~one & ((flags & mc.TRACK_AMBIGUOUS) == 0)
    apart = multi & ~one & ((flags & mc.TRACK_AMBIGUOUS) != 0)
    rows += [("one paddle per layer", int(one.sum()), "clean"),
             ("one paddle per layer, every hit < 3σ beyond its paddle", int(inside.sum()), "clean"),
             ("several paddles of one layer, side by side", int(adjacent.sum()), "several"),
             ("several paddles of one layer, not side by side", int(apart.sum()), "several")]
    return rows


def plot_event_quality(c):
    rows = event_quality(c)
    n_all = max(rows[0][1], 1)
    fig, ax = figure(7.2, 0.42 * len(rows) + 1.5)
    y = np.arange(len(rows))[::-1]
    colors = {"step": CATEGORICAL[0], "clean": CATEGORICAL[2], "several": CATEGORICAL[1]}
    for yi, (label, n, kind) in zip(y, rows):
        ax.barh(yi, max(n, 0.5), height=0.6, color=colors[kind], edgecolor=SURFACE, lw=1.0)
        ax.text(max(n, 0.5), yi, f"  {n}  ({100 * n / n_all:.3g} %)", va="center", fontsize=7, color=INK2)
    ax.set_yticks(y, [("    " if kind != "step" else "") + label for label, n, kind in rows])
    ax.set_xscale("log")
    ax.set_xlim(0.4, n_all * 12)
    ax.set_xlabel("events (log scale); in brackets: % of all events")
    ax.grid(axis="y", visible=False)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=colors["step"], label="selection steps, each inside the one above"),
                       Patch(color=colors["clean"], label="tracks: one paddle per layer"),
                       Patch(color=colors["several"], label="tracks: several paddles of one layer")],
              loc="lower right")
    mc.title(ax, f"Event quality: {c.short}",
             f"all {n_all} events = hits within {c.setup.events['window_ns']:g} ns; most have one end of one paddle only")
    c.save(fig, "event_quality")


# ---- event plots ----

def plot_layer_delay(c, la, lb):
    ab = c.H["kind"] == mc.HIT_AB
    t1 = np.sort(c.H["t_ns"][ab & (c.H["layer"] == la.index)])
    t2 = np.sort(c.H["t_ns"][ab & (c.H["layer"] == lb.index)])
    if not len(t1) or not len(t2):
        return
    span = 1000.0
    lo = np.searchsorted(t2, t1 - span, "left")
    hi = np.searchsorted(t2, t1 + span, "right")
    dts = np.concatenate([t2[l:h] - t for t, l, h in zip(t1, lo, hi)]) if len(t1) else np.zeros(0)
    window = float(c.setup.events["window_ns"])
    fig = plt.figure(figsize=(11.6, 4.4))
    from matplotlib.gridspec import GridSpec
    gs = GridSpec(1, 2, figure=fig, wspace=0.22, top=0.84)
    ax, az = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    width = 4.0
    bins = np.arange(-span, span + width, width)
    stairs(ax, dts, bins, CATEGORICAL[0], label=f"{len(dts)} pairs of two-end hits")
    chance = len(t1) * len(t2) * width * 1e-9 / c.span_s
    if chance >= 0.05:
        ax.axhline(chance, color=FIT_COLOR, lw=1.2, label=f"chance level: N1 N2 x bin / T = {chance:.2g} per bin")
    else:
        note(ax, f"chance level N1 N2 x bin / T = {chance:.1e} per bin:\nbelow the axis, no chance pairs expected",
             loc="upper left")
    for w in (-window, window):
        ax.axvline(w, color=MUTED, lw=0.8, ls="--")
    ax.set_yscale("log")
    ax.set_ylim(min(chance, 0.5) * 0.3 if chance >= 0.05 else 0.5, None)
    mc.headroom(ax, log=True)
    ax.set_xlabel(f"t({lb.name}) - t({la.name}) of two-end hits (ns)")
    ax.set_ylabel(f"pairs per {width:g} ns")
    ax.legend(loc="upper right")
    ax.set_title(f"±{span:g} ns; dashed: the event window ±{window:g} ns", loc="left", fontsize=8, color=INK2)
    # zoom on the peak
    near = dts[np.abs(dts) < 10]
    zw = 0.2 if len(near) >= 300 else 0.5
    zb = np.arange(-10, 10 + zw, zw)
    stairs(az, near, zb, CATEGORICAL[0], label=f"{len(near)} pairs within ±10 ns")
    if len(near):
        med = np.median(near)
        rob = 1.4826 * np.median(np.abs(near - med))
        dz = abs(np.mean([p.z for p in la.paddles]) - np.mean([p.z for p in lb.paddles]))
        az.axvline(dz / mc.C_CM_PER_NS, color=FIT_COLOR, lw=1.1, ls="--",
                   label=f"light over the {dz:g} cm between the layers: {dz / mc.C_CM_PER_NS:.2f} ns")
        note(az, f"median {med:+.2f} ns\nrobust σ {rob:.2f} ns\n(mean times (tA + tB)/2,\nno offsets, no path)",
             loc="upper left")
    az.set_xlim(-10, 10)
    az.set_ylim(0, None)
    mc.headroom(az, 1.35)
    az.set_xlabel(f"t({lb.name}) - t({la.name}) of two-end hits (ns)")
    az.set_ylabel(f"pairs per {zw:g} ns")
    az.legend(loc="upper right")
    az.set_title("zoom: the peak of particles crossing both layers", loc="left", fontsize=8, color=INK2)
    fig.text(0.125, 0.975, f"Time between the {la.name} and {lb.name} layers: {c.short}", ha="left", va="top",
             fontsize=9.5, color=INK)
    fig.text(0.125, 0.935, "every pair of two-end hits in the two layers; particles crossing both make the peak",
             ha="left", va="top", fontsize=8, color=INK2)
    c.save(fig, f"layer_delay_{la.name}_{lb.name}")


def plot_event_types(c):
    cats = c.einfo["categories"]
    names = list(cats)
    counts = np.array([cats[k] for k in names])
    fig, ax = figure(6.4, 0.38 * len(names) + 1.3)
    y = np.arange(len(names))[::-1]
    ax.barh(y, counts, height=0.6, color=CATEGORICAL[0])
    for yi, k in zip(y, counts):
        ax.text(k, yi, f"  {k}  ({k / c.span_s:.3g} Hz)", va="center", fontsize=7, color=INK2)
    ax.set_yticks(y, names)
    ax.set_xscale("log")
    ax.set_xlim(0.8, counts.max() * 8 if len(counts) else 10)
    ax.set_xlabel("events")
    ax.grid(axis="y", visible=False)
    mc.title(ax, f"Events by the layers with a two-end hit: {c.short}",
             f"{c.einfo['n_events']} events; one-end hits only = no paddle fired at both ends")
    c.save(fig, "event_types")


def plot_coincidences(c, la, lb):
    exp = (c.tinfo.get("expected") or {}).get("pairs", {}) if c.has_tracks else {}
    rows, cols = la.paddles, lb.paddles
    M = np.zeros((len(rows), len(cols)))
    X = np.full((len(rows), len(cols)), np.nan)
    for i, pa in enumerate(rows):
        for j, pb in enumerate(cols):
            if c.has_tracks:
                M[i, j] = ((c.Tr[f"paddle_{la.name}"] == pa.index) & (c.Tr[f"paddle_{lb.name}"] == pb.index)).sum()
            else:
                M[i, j] = c.einfo.get("pair_counts", {}).get(f"{pa.name}+{pb.name}", 0)
            key = f"{pa.name}+{pb.name}"
            if key in exp:
                X[i, j] = exp[key]["expected_count"]
    fig, ax = figure(1.4 * len(cols) + 2.8, 0.95 * len(rows) + 1.9)
    im = ax.imshow(M, cmap=mc.SEQUENTIAL, vmin=0, vmax=max(M.max(), 1), aspect="auto")
    for i in range(len(rows)):
        for j in range(len(cols)):
            dark = M[i, j] > 0.55 * max(M.max(), 1)
            text = f"{M[i, j]:.0f}"
            if np.isfinite(X[i, j]):
                text += f"\nsimulation {X[i, j]:.1f}" + (f"\nratio {M[i, j] / X[i, j]:.2f}" if X[i, j] > 0 else "")
            ax.text(j, i, text, ha="center", va="center", fontsize=7.5, color="#ffffff" if dark else INK)
    ax.set_xticks(range(len(cols)), [p.name for p in cols])
    ax.set_yticks(range(len(rows)), [p.name for p in rows])
    ax.set_xlabel(f"{lb.name} layer")
    ax.set_ylabel(f"{la.name} layer")
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label("tracks measured", color=INK2)
    cb.outline.set_visible(False)
    flux = c.setup.flux
    mc.title(ax, f"Tracks per paddle pair: {c.short}",
             (f"simulation (known values): muons with I = I0 cos{'²' if flux['cos_power'] == 2 else '^' + format(flux['cos_power'], 'g')}θ, I0 = "
              f"{flux['vertical_intensity']:g} /(m² s sr) (sea level, PDG),\nstraight through the setup's paddles, every "
              f"paddle fully efficient, x the run time; ratio = measured / simulation") if exp else None)
    c.save(fig, f"coincidences_{la.name}_{lb.name}")


def plot_rates_vs_time(c):
    ab = c.H["kind"] == mc.HIT_AB
    nb = int(np.clip(c.span_s / 2, 4, 100))
    edges = np.linspace(0, c.span_s, nb + 1)
    width = edges[1] - edges[0]
    fig, ax = figure()
    for layer in c.setup.layers:
        t = c.H["t_ns"][ab & (c.H["layer"] == layer.index)] * 1e-9
        k, _ = np.histogram(t, bins=edges)
        ax.stairs(k / width, edges, color=mc.layer_color(layer.index), lw=1.2,
                  label=f"{layer.name} layer, two-end hits ({len(t) / c.span_s:.2f} Hz)")
    multi = c.E["n_layers_ab"] >= 2
    tc = c.E["t_ns"][multi] * 1e-9
    if len(tc):
        nc = int(np.clip(len(tc) // 5, 1, nb))
        ec = np.linspace(0, c.span_s, nc + 1)
        k, _ = np.histogram(tc, bins=ec)
        w = ec[1] - ec[0]
        down, up = mc.poisson_band(k)
        ax.errorbar(0.5 * (ec[1:] + ec[:-1]), k / w, xerr=0.5 * w, yerr=[down / w, up / w], fmt="o", ms=4,
                    color=DATA_COLOR, mec=SURFACE, mew=1.0, elinewidth=0.9, capsize=0,
                    label=f"coincidences of two layers or more ({len(tc) / c.span_s:.3f} Hz)")
    ax.set_yscale("log")
    ax.set_xlabel("time in the run (s)")
    ax.set_ylabel("rate (Hz)")
    ax.set_xlim(0, c.span_s)
    mc.headroom(ax, log=True)
    ax.legend(loc="center right")
    mc.title(ax, f"Rates over the run: {c.short}")
    c.save(fig, "rates_vs_time")


# ---- track plots ----

def outline(ax, p, color, fill=False, lw=0.8):
    rect = Rectangle((p.x - 0.5 * p.size_x, p.y - 0.5 * p.size_y), p.size_x, p.size_y, fill=False,
                     edgecolor=color, lw=lw)
    ax.add_patch(rect)
    if fill:
        ax.add_patch(Rectangle((p.x - 0.5 * p.size_x, p.y - 0.5 * p.size_y), p.size_x, p.size_y,
                               facecolor=color, alpha=0.10, lw=0))


def plot_hitmap(c, layer):
    """Where the tracks cross the layer. A track through two neighbouring paddles of the layer
    shows each paddle's own hit (never a point between them)."""
    fig, ax = figure(6.0, 5.6)
    for p in c.setup.paddles:
        if p.layer is layer:
            outline(ax, p, mc.layer_color(layer.index), fill=True, lw=1.0)
            if p.axis == "x":   # label just beyond the paddle's +x end
                ax.text(p.x + 0.5 * p.size_x + 0.6, p.y, p.name, fontsize=6.5, color=INK2, ha="left",
                        va="center", zorder=2)
            else:               # label just beyond its +y end, along its length
                ax.text(p.x, p.y + 0.5 * p.size_y + 0.6, p.name, fontsize=6.5, color=INK2, ha="center",
                        va="bottom", rotation=90, zorder=2)
        else:
            outline(ax, p, AXIS, lw=0.7)
    if c.Pts is not None:
        sel = c.Pts["layer"] == layer.index
        track = c.Pts["track"][sel]
        ids, counts = np.unique(track, return_counts=True)
        multi = np.isin(track, ids[counts > 1])
        x, y = c.Pts["x"][sel], c.Pts["y"][sel]
        ax.plot(x[~multi], y[~multi], "o", ms=4.5, color=DATA_COLOR, mec=SURFACE, mew=1.0,
                label=f"{(~multi).sum()} tracks, one paddle of this layer", zorder=4)
        if multi.any():
            ax.plot(x[multi], y[multi], "o", ms=4.5, color=CATEGORICAL[1], mec=SURFACE, mew=1.0, zorder=5,
                    label=f"{int((counts > 1).sum())} tracks through two neighbouring paddles: each paddle's hit")
    else:   # older tracks.root without the points tree
        s = c.Tr[f"paddle_{layer.name}"] >= 0
        ax.plot(c.Tr[f"x_{layer.name}"][s], c.Tr[f"y_{layer.name}"][s], "o", ms=4.5, color=DATA_COLOR, mec=SURFACE,
                mew=1.0, label=f"{s.sum()} tracks", zorder=4)
    ax.set_aspect("equal")
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.margins(0.08)
    ax.legend(loc="lower right")
    mc.title(ax, f"Where the tracks cross the {layer.name} layer: {c.short}",
             f"filled: {layer.name} paddles; grey: the other layers; across a paddle: its centre line, "
             f"along it: the best estimate inside the paddle")
    c.save(fig, f"hitmap_{layer.name}")


def plot_angle(c, key, name, unit_range, step_small, step_large):
    ok = np.isfinite(c.Tr[key])
    vals = c.Tr[key][ok]
    if not len(vals):
        return
    step = step_small if len(vals) >= 200 else step_large
    bins = np.arange(unit_range[0], unit_range[1] + step, step)
    fig, ax = figure()
    if c.toy is not None and len(c.toy["weight_hz"]):
        w = c.toy["weight_hz"]
        k, e = np.histogram(c.toy[key], bins=bins, weights=w)
        scale = len(vals) / k.sum() if k.sum() > 0 else 0
        ax.stairs(k * scale, e, color=mc.EXPECTED_COLOR, lw=1.3, label="toy cosmic muons, scaled to the data")
        ax.stairs(k * scale, e, fill=True, color=mc.EXPECTED_COLOR, alpha=0.10, lw=0)
    points(ax, vals, bins, label=f"{len(vals)} tracks")
    ax.set_xlim(*unit_range)
    ax.set_ylim(0, None)
    mc.headroom(ax)
    ax.set_xlabel(name)
    ax.set_ylabel(f"tracks per {step:g}°")
    ax.legend(loc="upper right")
    return fig, ax


def plot_zenith(c):
    out = plot_angle(c, "theta_deg", "zenith angle θ (degrees from the vertical)", (0, 90), 5, 10)
    if out:
        fig, ax = out
        mc.title(ax, f"Zenith angle of the tracks: {c.short}",
                 "the toy is reconstructed like the data (positions smeared, across a paddle at its centre)")
        c.save(fig, "zenith")


def plot_azimuth(c):
    out = plot_angle(c, "phi_deg", "azimuth φ of the arrival direction (0 = from +x, 90 = from +y)", (0, 360), 5, 5)
    if out:
        fig, ax = out
        ax.set_xticks(range(0, 361, 45))
        mc.title(ax, f"Azimuth of the tracks: {c.short}",
                 "5° bins; the setup's square acceptance favours the diagonals a little (see the simulation)")
        c.save(fig, "azimuth")


def resolution_tof(c):
    """Expected sigma of the time of flight from the paddles' mean-time resolutions."""
    sig = [c.calib["paddles"].get(p.name, {}).get("mean_time_sigma_ns") for p in c.setup.paddles]
    sig = [s for s in sig if s and np.isfinite(s)]
    return math.sqrt(2) * float(np.median(sig)) if sig else float("nan")


def plot_tof(c):
    ok = np.isfinite(c.Tr["tof_ns"]) & np.isfinite(c.Tr["s_cm"])
    d = (c.Tr["tof_ns"] - c.Tr["tof_expected_ns"])[ok]
    if not len(d):
        return
    fig, ax = figure()
    lim = max(4.0, np.percentile(np.abs(d), 98) * 1.3)
    bw = 0.25 if len(d) >= 100 else 0.5
    bins = np.arange(-lim, lim + bw, bw)
    points(ax, np.clip(d, -lim, lim - 1e-6), bins, label=f"{len(d)} tracks (measured)")
    s_exp = resolution_tof(c)
    if np.isfinite(s_exp):
        x = np.linspace(-lim, lim, 400)
        ax.plot(x, len(d) * bw * np.exp(-0.5 * (x / s_exp) ** 2) / (s_exp * math.sqrt(2 * math.pi)),
                color=mc.EXPECTED_COLOR, lw=1.3,
                label=f"simulation: β = 1, smeared by the measured resolution σ = {s_exp:.2f} ns")
    robust = 1.4826 * np.median(np.abs(d - np.median(d)))
    note(ax, f"median {np.median(d):+.3f} ns\nrobust σ {robust:.3f} ns\nmedian s/c {np.median(c.Tr['tof_expected_ns'][ok]):.3f} ns",
         loc="upper left")
    ax.set_xlabel("measured time of flight - s / c (ns); no layer time offsets applied")
    ax.set_ylabel(f"tracks per {bw:g} ns")
    ax.set_ylim(0, None)
    mc.headroom(ax)
    ax.legend(loc="upper right")
    mc.title(ax, f"Time of flight between the layers: {c.short}",
             "s = path between the two hit points; the paddles' cable delays are not known, so the median may be off 0")
    c.save(fig, "tof_residual")


def plot_beta(c):
    """beta = path / (c x time of flight) as measured, over its whole range: shaded where no
    particle going down can be (beta > 1: faster than light; beta < 0: the lower layer first)."""
    ok = np.isfinite(c.Tr["beta"]) & np.isfinite(c.Tr["s_cm"])
    b = c.Tr["beta"][ok]
    if not len(b):
        return
    s = c.Tr["s_cm"][ok]
    lo, hi = -1.0, 3.0
    width = 0.05 if len(b) >= 200 else 0.1
    edges = np.arange(lo, hi + width / 2, width)

    def counts(values, weights=None):
        w = np.ones(len(values)) if weights is None else weights
        inside = (values >= lo) & (values < hi)
        k, _ = np.histogram(values[inside], bins=edges, weights=w[inside])
        return k, float(w[values < lo].sum()), float(w[values >= hi].sum())

    k, under, over = counts(b)
    fig, ax = figure(7.4, 4.4)
    ax.axvspan(lo - 4 * width, 0, color=GRID, alpha=0.55, lw=0)
    ax.axvspan(1, hi + 4 * width, color=AXIS, alpha=0.35, lw=0)
    x_under, x_over = lo - 2.5 * width, hi + 2.5 * width
    centres = 0.5 * (edges[1:] + edges[:-1])
    allc = np.concatenate([[under], k, [over]])
    allx = np.concatenate([[x_under], centres, [x_over]])
    down, up = mc.poisson_band(allc)
    m = allc > 0
    ax.errorbar(allx[m], allc[m], yerr=[down[m], up[m]], fmt="o", ms=3.6, color=DATA_COLOR, mec=SURFACE, mew=0.8,
                elinewidth=0.8, capsize=0, label=f"{len(b)} tracks (measured, no layer time offsets)", zorder=4)
    s_exp = resolution_tof(c)
    if np.isfinite(s_exp):   # simulation: beta = 1 particles on the same paths, measured with this resolution
        rng = np.random.default_rng(3)
        reps = 400
        ss = np.repeat(s, reps)
        tof = ss / mc.C_CM_PER_NS + rng.normal(0, s_exp, len(ss))
        si, su, so = counts(ss / (mc.C_CM_PER_NS * tof), np.full(len(ss), 1.0 / reps))
        ax.stairs(si, edges, color=mc.EXPECTED_COLOR, lw=1.3,
                  label=f"simulation: β = 1, measured with σ(tof) = {s_exp:.2f} ns")
        for x0, val in ((x_under, su), (x_over, so)):
            ax.plot([x0 - width / 2, x0 + width / 2], [val, val], color=mc.EXPECTED_COLOR, lw=1.3)
    ax.axvline(1, color=MUTED, lw=0.9)
    ax.axvline(0, color=MUTED, lw=0.9)
    ticks = [-0.5, 0, 0.5, 1, 1.5, 2, 2.5]
    ax.set_xticks(ticks + [x_under, x_over], [f"{t:g}" for t in ticks] + [f"< {lo:g}", f"> {hi:g}"])
    ax.set_xlim(lo - 4 * width, hi + 4 * width)
    ax.set_ylim(0, None)
    mc.headroom(ax, 1.45)
    ytop = ax.get_ylim()[1]
    ax.text(1 + 0.04, ytop * 0.97, "β > 1: faster than light,\nonly the time resolution", fontsize=7, color=INK2,
            va="top")
    ax.text(lo - 3.6 * width, ytop * 0.97, "β < 0: lower layer first\n(upward, or resolution)", fontsize=7,
            color=INK2, va="top")
    ax.set_xlabel(f"β = path / (c x time of flight); the outer bins hold the tracks below {lo:g} and above {hi:g}")
    ax.set_ylabel(f"tracks per {width:g}")
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.86))
    n_fast, n_neg = int((b > 1).sum()), int((b < 0).sum())
    mc.title(ax, f"Speed of the tracks: {c.short}",
             f"{n_fast} of {len(b)} tracks ({n_fast / len(b):.0%}) measured at β > 1, {n_neg} below 0: the flight "
             f"({np.median(s) / mc.C_CM_PER_NS:.2f} ns) is only ~{(np.median(s) / mc.C_CM_PER_NS) / s_exp:.0f} x the "
             f"time resolution")
    c.save(fig, "beta")


def plot_energy_deposit(c):
    """The tracks' signal in each paddle, as measured (ADC): energy deposit x light yield x gain."""
    fig, ax = figure()
    vals = {}
    for p in c.setup.paddles[:len(CATEGORICAL)]:
        L = p.layer.name
        sel = (c.Tr[f"paddle_{L}"] == p.index) & ((c.Tr["flags"] & (mc.TRACK_ADJACENT | mc.TRACK_AMBIGUOUS)) == 0)
        v = c.Tr[f"hgeo_{L}"][sel]
        vals[p.name] = v[np.isfinite(v)]
    allv = np.concatenate(list(vals.values())) if vals else np.zeros(0)
    if not len(allv):
        plt.close(fig)
        return
    top = np.percentile(allv, 98) * 1.25
    bins = np.linspace(0, top, 41 if len(allv) >= 400 else 21)
    lines = []
    for k, p in enumerate(c.setup.paddles[:len(CATEGORICAL)]):
        v = vals[p.name]
        if not len(v):
            continue
        stairs(ax, np.clip(v, 0, top * 0.999), bins, CATEGORICAL[k], label=f"{p.name} ({len(v)} tracks)", wash=False)
        mpv = c.calib["paddles"].get(p.name, {}).get("mip_adc")
        if mpv:
            lines.append(f"{p.name}: most probable {mpv:.0f} ADC")
    sim = c.tinfo.get("energy_sim_mpv_mev") if c.has_tracks else None
    if sim:
        lines.append("simulation (known muon values): muons here lose")
        lines.append("most probably " + ", ".join(f"{k} {v:.2f}" for k, v in sim.items()) + " MeV")
        lines.append("(energy_deposit_sim.png); ADC to MeV needs an")
        lines.append("energy calibration, e.g. a gamma source")
    note(ax, "\n".join(lines), loc="upper right")
    ax.set_xlabel("sqrt(hA x hB) of the track's hit, ADC (as measured, no angle correction)")
    ax.set_ylabel("tracks per bin")
    ax.set_ylim(0, None)
    mc.headroom(ax, 1.5)
    ax.legend(loc="upper left")
    mc.title(ax, f"Energy deposit of the tracks, measured signal: {c.short}",
             "signal = deposited energy x light yield x gain: each paddle has its own scale")
    c.save(fig, "energy_deposit")


def energy_simulation(c, seed=4):
    """Energy the toy muons lose in each layer's paddle (MeV): PDG Landau most probable loss and width for
    the path thickness / cos(theta) through plastic, a Moyal-shaped spread. Weighted like the toy."""
    if c.toy is None or not len(c.toy["weight_hz"]):
        return None
    rng = np.random.default_rng(seed)
    out = {}
    cos_t = np.cos(np.radians(c.toy["theta_true_deg"]))
    for layer in c.setup.layers:
        sel = c.toy[f"paddle_{layer.name}"] >= 0
        if not sel.any():
            continue
        thick = np.array([c.setup.paddles[k].thickness for k in c.toy[f"paddle_{layer.name}"][sel]])
        path = thick / cos_t[sel]
        e = mc.sample_moyal(rng, mc.landau_mpv_mev(path), mc.landau_xi_mev(path))
        out[layer.name] = (e, c.toy["weight_hz"][sel])
    return out


def plot_energy_sim(c):
    sim = energy_simulation(c)
    if not sim:
        return
    fig, ax = figure()
    bins = np.linspace(0, 3.0, 61)
    mpvs = {}
    for layer in c.setup.layers:
        if layer.name not in sim:
            continue
        e, w = sim[layer.name]
        k, _ = np.histogram(np.clip(e, 0, 2.999), bins=bins, weights=w)
        k = k / max(k.sum(), 1e-30)
        ax.stairs(k, bins, color=mc.layer_color(layer.index), lw=1.3, label=f"{layer.name} layer")
        j = int(np.argmax(k))
        mpvs[layer.name] = 0.5 * (bins[j] + bins[j + 1])
    t0 = c.setup.paddles[0].thickness
    note(ax, f"SIMULATION with known muon values:\nminimum-ionising muons (βγ = 20), PDG Landau formula;\n"
             f"a vertical muon in {t0:g} cm of plastic: most probable {mc.landau_mpv_mev(t0):.2f} MeV;\n"
             f"paths = thickness / cos θ of the toy's tracks;\nno light yield, photo-statistics or noise",
         loc="upper right")
    ax.set_xlabel("energy lost in the paddle (MeV)")
    ax.set_ylabel("fraction of the tracks per 0.05 MeV")
    ax.set_ylim(0, None)
    mc.headroom(ax, 1.5)
    ax.legend(loc="upper left")
    mc.title(ax, f"Expected energy deposit of muon tracks (simulation): {c.short}",
             "the toy's tracks through the setup; compare the shape, not the scale, with energy_deposit.png")
    c.save(fig, "energy_deposit_sim")


def plot_signal_vs_path(c):
    """Does the signal grow with the path through the paddle? Each hit's signal over its paddle's most
    probable signal, against 1 / cos(zenith) of its track: a data test of the angle correction."""
    rows_f, rows_r = [], []
    cos_t = np.cos(np.radians(c.Tr["theta_deg"]))
    clean = (c.Tr["flags"] & (mc.TRACK_ADJACENT | mc.TRACK_AMBIGUOUS)) == 0
    for p in c.setup.paddles:
        mpv = c.calib["paddles"].get(p.name, {}).get("mip_adc")
        if not mpv:
            continue
        L = p.layer.name
        sel = clean & (c.Tr[f"paddle_{L}"] == p.index) & np.isfinite(c.Tr[f"hgeo_{L}"]) & (cos_t > 0)
        rows_f.append(1 / cos_t[sel])
        rows_r.append(c.Tr[f"hgeo_{L}"][sel] / mpv)
    if not rows_f:
        return
    f, r = np.concatenate(rows_f), np.concatenate(rows_r)
    if len(f) < 30:
        return
    edges = np.unique(np.percentile(f, np.linspace(0, 100, 9 if len(f) >= 400 else 5)))
    xm, ym, ye = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        s = (f >= a) & (f <= b)
        if s.sum() < 10:
            continue
        med = np.median(r[s])
        rs = 1.4826 * np.median(np.abs(r[s] - med))
        xm.append(np.median(f[s]))
        ym.append(med)
        ye.append(1.2533 * rs / math.sqrt(s.sum()))
    xm, ym, ye = map(np.array, (xm, ym, ye))
    fig, ax = figure()
    ax.errorbar(xm, ym, yerr=ye, fmt="o", ms=5, color=DATA_COLOR, mec=SURFACE, mew=1.0, elinewidth=0.9, capsize=0,
                label=f"median of {len(f)} hits per bin (tracks, every paddle)", zorder=4)
    w = 1 / np.maximum(ye, 1e-9) ** 2
    a_prop = float(np.sum(w * ym * xm) / np.sum(w * xm * xm))       # ym = a x f
    a_flat = float(np.sum(w * ym) / np.sum(w))                       # ym = a
    chi_prop = float(np.sum(w * (ym - a_prop * xm) ** 2))
    chi_flat = float(np.sum(w * (ym - a_flat) ** 2))
    xx = np.linspace(1, max(xm.max(), 1.05) * 1.02, 50)
    ax.plot(xx, a_prop * xx, color=FIT_COLOR, lw=1.3,
            label=f"signal ∝ path length (χ² = {chi_prop:.1f} for {len(xm) - 1})")
    ax.plot(xx, np.full_like(xx, a_flat), color=AXIS, lw=1.3,
            label=f"no dependence on the path (χ² = {chi_flat:.1f} for {len(xm) - 1})")
    ax.set_xlabel("1 / cos(zenith) of the track = path through the paddle / its thickness")
    ax.set_ylabel("signal / the paddle's most probable signal")
    mc.headroom(ax, 1.25)
    ax.legend(loc="upper left")
    mc.title(ax, f"Signal against the path through the paddle: {c.short}",
             "a particle crossing at θ goes through thickness / cos θ; the reconstructed θ (error ~7°) flattens the trend")
    c.save(fig, "signal_vs_path")


def sim_mpv(c):
    """Most probable simulated energy loss per layer (MeV), quoted as text on the data plot."""
    sim = energy_simulation(c)
    if not sim:
        return
    out = {}
    bins = np.linspace(0, 3.0, 121)
    for name, (e, w) in sim.items():
        k, _ = np.histogram(e, bins=bins, weights=w)
        j = int(np.argmax(k))
        out[name] = 0.5 * (bins[j] + bins[j + 1])
    c.tinfo["energy_sim_mpv_mev"] = out


# ---- summary ----

def headline(c):
    s = c.setup
    lines = [f"Paddle muon analysis of {c.run}", "=" * (29 + len(c.run)), ""]
    conf = c.pinfo.get("conf") or {}
    lines.append(f"Data      : {c.pinfo['raw_file']}")
    lines.append(f"Setup     : {s.name} ({len(s.paddles)} paddles, {len(s.channels)} channels, "
                 f"{len(s.layers)} layers)")
    if s.geometry_source:
        lines.append(f"Geometry  : {s.geometry_source}")
    lines.append(f"Run       : {c.span_s:.2f} s of data (daqruntime {conf.get('daqruntime', '?')} s), trigmode "
                 f"{conf.get('trigmode', '?')}, started {conf.get('started', '?')}")
    if conf.get("note"):
        lines.append(f"Note      : {conf['note']}")
    lines.append("")
    rej = sum(v["records"] - v["used"] for v in c.pinfo["channels"].values())
    tot = sum(v["records"] for v in c.pinfo["channels"].values())
    lines.append(f"Pulses    : {tot} waveforms, {tot - rej} used, {rej} rejected (mostly noise bursts)")
    for p in s.paddles:
        e = c.einfo["paddles"].get(p.name, {})
        cal = c.calib["paddles"].get(p.name, {})
        thr = []
        for end, chid in p.channels.items():
            cc = c.calib.get("channels", {}).get(str(chid), {})
            if "threshold_mip" in cc:
                thr.append(f"{end} {cc['threshold_mip']:.2f}")
        lines.append(f"  {p.name:<6} two-end hits {e.get('AB', 0):>6} ({e.get('AB_rate_hz', float('nan')):.2f} Hz, "
                     f"{e.get('both_percent', float('nan')):.0f} % of its hits), v = {cal.get('v_eff_cm_per_ns', float('nan')):.1f} cm/ns, "
                     f"σ(tA-tB) = {cal.get('tdiff_sigma_ns', float('nan')):.2f} ns, position σ = "
                     f"{cal.get('position_sigma_cm', float('nan')):.1f} cm, muon peak {cal.get('mip_adc') or float('nan'):.0f} ADC"
                     + (f", threshold/peak {', '.join(thr)}" if thr else ""))
    sp = c.calib.get("speed", {})
    if sp.get("mode") != "shared":
        lines.append("Speed     : each paddle its own signal speed: " + ", ".join(
            f"{p.name} {c.calib['paddles'][p.name]['v_eff_cm_per_ns']:.1f} ± "
            f"{c.calib['paddles'][p.name].get('v_eff_err_cm_per_ns', float('nan')):.1f}"
            for p in s.paddles if p.name in c.calib["paddles"]) + " cm/ns")
    if sp.get("mode") == "shared":
        lines.append(f"Speed     : one signal speed for {sp['paddles']} paddles, {sp['v_eff_cm_per_ns']:.2f} ± "
                     f"{sp['err']:.2f} cm/ns; each alone: " +
                     ", ".join(f"{p.name} {c.calib['paddles'][p.name]['v_free_cm_per_ns']:.1f}"
                               for p in s.paddles if c.calib['paddles'].get(p.name, {}).get('v_free_cm_per_ns')))
    rhos = []
    from scipy.stats import spearmanr
    for p in s.paddles:
        sel = (c.H["paddle"] == p.index) & (c.H["kind"] == mc.HIT_AB) & (c.H["hA"] > 0) & (c.H["hB"] > 0)
        if sel.sum() >= 20:
            rhos.append(f"{p.name} {spearmanr(c.H['tdiff_ns'][sel], np.log(c.H['hA'][sel] / c.H['hB'][sel]))[0]:+.2f}")
    if rhos:
        lines.append("Light/time: rank correlation of ln(hA/hB) with tA - tB, raw (no calibration): " + ", ".join(rhos))
    for w in c.calib.get("warnings", []):
        lines.append(f"WARNING   : {w}")
    lines.append("")
    rows = event_quality(c)
    n_all = max(rows[0][1], 1)
    lines.append("Event quality (each step inside the one above; the indented lines split the tracks):")
    for label, n, kind in rows:
        lines.append(f"  {('  ' if kind != 'step' else '') + label:<62}{n:>8}  {100 * n / n_all:7.3f} %")
    lines.append("")
    for key, co in c.einfo.get("coincidences", {}).items():
        lines.append(f"Coincid.  : {key}: {co['events']} events = {co['rate_hz']:.3f} ± "
                     f"{math.sqrt(co['events']) / c.span_s:.3f} Hz (by chance {co['chance_calculated']:.2g})")
    if c.has_tracks:
        tr = c.Tr
        ok = np.isfinite(tr["theta_deg"])
        exp = c.tinfo.get("expected") or {}
        lines.append(f"Tracks    : {c.tinfo['n_tracks']}; zenith mean {np.mean(tr['theta_deg'][ok]) if ok.any() else float('nan'):.1f}°"
                     + (f" (toy {exp['toy_mean_theta_deg']:.1f}°)" if exp else ""))
        if exp:
            lines.append(f"Expected  : {exp['total_hz'] * c.span_s:.1f} tracks for fully efficient paddles; "
                         f"measured / expected = {exp['ratio']:.2f} ± {exp['ratio_err']:.2f}")
        good = np.isfinite(tr["tof_ns"])
        if good.any():
            d = (tr["tof_ns"] - tr["tof_expected_ns"])[good]
            lines.append(f"TOF       : measured tof - s/c median {np.median(d):+.2f} ns, robust σ "
                         f"{1.4826 * np.median(np.abs(d - np.median(d))):.2f} ns (simulation of β = 1 with the "
                         f"measured resolution: σ {resolution_tof(c):.2f} ns) over a median path of "
                         f"{np.median(tr['s_cm'][good]):.1f} cm; no layer time offsets applied")
            lines.append(f"β         : median {np.nanmedian(tr['beta'][good]):.2f}; "
                         f"{int((tr['beta'][good] >= 1).sum())} of {int(good.sum())} tracks measured at β ≥ 1, "
                         f"{int((tr['beta'][good] < 0).sum())} below 0 (time resolution)")
        for r in c.tinfo.get("map_check", []):
            lines.append(f"Map check : {r['layer']} layer {r['axis']} ({r['method']}) by {r['by_layer']} paddle: "
                         f"slope {r['slope']:+.2f} ± {r['slope_err']:.2f} (+1 expected), shift "
                         f"{r['shift_cm']:+.1f} ± {r['shift_err_cm']:.1f} cm -> {r['verdict']}")
    lines.append("")
    lines.append(f"Plots     : {len(c.made)} files in {c.plots}")
    return lines


def write_summary(c):
    parts = headline(c)
    for name in ("pulses_summary.txt", "events_summary.txt", "calibration_summary.txt", "tracks_summary.txt"):
        path = os.path.join(c.folder, name)
        if os.path.isfile(path):
            with open(path, encoding="utf-8") as f:
                parts += ["", "", f"---- {name} " + "-" * max(0, 60 - len(name)), "", f.read().rstrip()]
    path = os.path.join(c.folder, "analysis_summary.txt")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(parts) + "\n")
    print("\n".join(headline(c)))
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", default=None,
                    help="Analysis folder with the outputs of steps 1-4 (asked for if left out)")
    ap.add_argument("--setup", default=None, help="Setup file (default: the one step 1 used)")
    args = ap.parse_args()
    folder = mc.ask_analysis_folder(args.folder, mc.EVENTS_FILE, "build_events.py")
    run(folder, args.setup)


def run(folder, setup_path=None):
    """Step 5 for one analysis folder; returns the list of plot files."""
    c = Context(folder, setup_path)
    s = c.setup
    jobs = [("pulse rates", lambda: plot_pulse_rates(c)), ("baseline noise", lambda: plot_baseline_noise(c))]
    for p in s.paddles:
        jobs += [(f"{p.name} shape", lambda p=p: plot_pulse_shape(c, p)),
                 (f"{p.name} heights", lambda p=p: plot_heights(c, p))]
        if p.both_ends:
            jobs += [(f"{p.name} tA-tB", lambda p=p: plot_tdiff(c, p)),
                     (f"{p.name} position", lambda p=p: plot_position(c, p)),
                     (f"{p.name} attenuation", lambda p=p: plot_attenuation(c, p)),
                     (f"{p.name} light vs timing", lambda p=p: plot_tdiff_vs_ratio(c, p)),
                     (f"{p.name} muon peak", lambda p=p: plot_mip(c, p))]
    for i, la in enumerate(s.layers):
        for lb in s.layers[i + 1:]:
            jobs.append((f"{la.name}-{lb.name} delay", lambda la=la, lb=lb: plot_layer_delay(c, la, lb)))
            jobs.append((f"{la.name}-{lb.name} pairs", lambda la=la, lb=lb: plot_coincidences(c, la, lb)))
    jobs += [("event types", lambda: plot_event_types(c)), ("event quality", lambda: plot_event_quality(c)),
             ("rates", lambda: plot_rates_vs_time(c))]
    if c.has_tracks and len(c.Tr["event"]):
        for layer in s.layers:
            jobs.append((f"{layer.name} hit map", lambda layer=layer: plot_hitmap(c, layer)))
        jobs += [("zenith", lambda: plot_zenith(c)), ("azimuth", lambda: plot_azimuth(c)),
                 ("time of flight", lambda: plot_tof(c)), ("beta", lambda: plot_beta(c)),
                 ("energy (simulation)", lambda: sim_mpv(c)), ("energy deposit", lambda: plot_energy_deposit(c)),
                 ("energy simulation", lambda: plot_energy_sim(c)), ("signal vs path", lambda: plot_signal_vs_path(c))]
    bar = mc.Progress("Plotting", len(jobs))
    for name, job in jobs:
        job()
        bar.update(1, name)
    bar.close()
    path = write_summary(c)
    print(f"Saved {len(c.made)} plots in {c.plots} and {path}")
    return c.made


if __name__ == "__main__":
    main()
