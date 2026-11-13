#!/usr/bin/env python3
"""
event_display.py  -  pictures of single events of the paddle muon analysis

For every chosen track it draws
    event_<n>.png            the setup from above (x-y) and from two sides (x-z,
                             y-z): the paddles that fired, the hit points with
                             their resolution, and the fitted straight line
    event_<n>_waveforms.png  the event's waveforms on one time axis (left) and their
                             rising edges divided by the height (right), each with
                             its pulse time marked (needs the raw ROOT file)
into event_displays/ in the analysis folder. n = the event number of events.root.

Usage:
    python3 event_display.py                       # asks for the analysis folder; the first 10 tracks
    python3 event_display.py FOLDER --first 20
    python3 event_display.py FOLDER --events 12,57,301
    python3 event_display.py FOLDER --all
    python3 event_display.py FOLDER --several 5     # also 5 tracks with several paddles of one layer
    python3 event_display.py FOLDER --no-waveforms
"""

import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse
import math
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

import muon_common as mc
from muon_common import INK, INK2, MUTED, AXIS, SURFACE, CATEGORICAL


def parse_events(text):
    out = []
    for part in str(text).replace(" ", "").split(","):
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


class Data:
    def __init__(self, folder, setup_path=None):
        self.folder = folder
        self.setup = mc.setup_of(folder, setup_path)
        pp = mc.need_file(folder, mc.PULSES_FILE, "extract_pulses.py")
        ep = mc.need_file(folder, mc.EVENTS_FILE, "build_events.py")
        tp = mc.need_file(folder, mc.TRACKS_FILE, "reconstruct.py")
        self.pinfo = mc.read_info(pp)
        self.tinfo = mc.read_info(tp)
        self.P = mc.read_tree(pp, "pulses", ["channel", "entry", "ttt", "t_ns", "edge", "baseline", "height", "flags"])
        self.H = mc.read_tree(ep, "hits")
        self.Tr = mc.read_tree(tp, "tracks")
        cal = self.tinfo.get("calibration") or os.path.join(folder, mc.CALIB_FILE)
        self.calib = mc.read_json(cal) if os.path.isfile(cal) else {"paddles": {}}
        self.run = self.pinfo["run"]
        self.short = mc.short_run(self.run)
        self.raw = None
        raw = self.pinfo.get("raw_file")
        if raw and os.path.isfile(raw):
            self.raw = mc.open_root(raw)

    def waveform(self, chid, entry):
        if self.raw is None:
            return None
        tree = self.raw.Get(f"ch{chid}")
        if not tree:
            return None
        tree.GetEntry(int(entry))
        return np.frombuffer(tree.chdata, dtype=np.uint16, count=tree.nsamps).astype(np.float64)


def draw_paddle(ax, p, view, color, fired, single):
    """One paddle in a view: 'xy' (from above), 'xz' or 'yz' (from the side)."""
    if view == "xy":
        x0, w, y0, h = p.x - 0.5 * p.size_x, p.size_x, p.y - 0.5 * p.size_y, p.size_y
    elif view == "xz":
        x0, w, y0, h = p.x - 0.5 * p.size_x, p.size_x, p.z - 0.5 * p.thickness, p.thickness
    else:
        x0, w, y0, h = p.y - 0.5 * p.size_y, p.size_y, p.z - 0.5 * p.thickness, p.thickness
    if fired:
        ax.add_patch(Rectangle((x0, y0), w, h, facecolor=color, alpha=0.22 if not single else 0.08, lw=0))
    ax.add_patch(Rectangle((x0, y0), w, h, fill=False, edgecolor=color if fired else AXIS,
                           lw=1.0 if fired else 0.7, ls="-" if not single else (0, (3, 2))))


def geometry_figure(d, row, out_dir):
    s = d.setup
    T = d.Tr
    e = int(T["event"][row])
    ev_hits = np.where(d.H["event"] == e)[0]
    fired = {}
    for i in ev_hits:
        k = int(d.H["paddle"][i])
        ab = d.H["kind"][i] == mc.HIT_AB
        fired[k] = fired.get(k, False) or ab
    fig = plt.figure(figsize=(10.2, 5.8))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.1, 1.0], hspace=0.42, wspace=0.22, top=0.86)
    ax_xy = fig.add_subplot(gs[:, 0])
    ax_xz = fig.add_subplot(gs[0, 1])
    ax_yz = fig.add_subplot(gs[1, 1])
    zs = [p.z for p in s.paddles]
    z_lo, z_hi = min(zs) - 8, max(zs) + 8
    zz = np.linspace(z_lo, z_hi, 50)
    # the fitted line: x(z) = ax + bx (z - z0), y(z) = ay + by (z - z0)
    lx = T["ax"][row] + T["bx"][row] * (zz - T["z0"][row])
    ly = T["ay"][row] + T["by"][row] * (zz - T["z0"][row])
    pts = [(T[f"x_{L.name}"][row], T[f"y_{L.name}"][row], T[f"z_{L.name}"][row], T[f"sx_{L.name}"][row],
            T[f"sy_{L.name}"][row], L) for L in s.layers if T[f"paddle_{L.name}"][row] >= 0]
    for ax, view in ((ax_xy, "xy"), (ax_xz, "xz"), (ax_yz, "yz")):
        for p in s.paddles:
            draw_paddle(ax, p, view, mc.layer_color(p.layer.index), p.index in fired,
                        p.index in fired and not fired[p.index])
    ax_xy.plot(lx, ly, color=INK, lw=1.0, zorder=3)
    ax_xz.plot(lx, zz, color=INK, lw=1.0, zorder=3)
    ax_yz.plot(ly, zz, color=INK, lw=1.0, zorder=3)
    for x, y, z, sx, sy, L in pts:
        col = mc.layer_color(L.index)
        kw = dict(fmt="o", ms=5.5, color=col, mec=SURFACE, mew=1.2, elinewidth=1.4, capsize=0, zorder=4)
        ax_xy.errorbar(x, y, xerr=sx, yerr=sy, label=f"{L.name} layer", **kw)
        ax_xz.errorbar(x, z, xerr=sx, **kw)
        ax_yz.errorbar(y, z, xerr=sy, **kw)
    ext_x = max(abs(p.x) + 0.5 * p.size_x for p in s.paddles) + 4
    ext_y = max(abs(p.y) + 0.5 * p.size_y for p in s.paddles) + 4
    ext = max(ext_x, ext_y)
    ax_xy.set_xlim(-ext, ext)
    ax_xy.set_ylim(-ext, ext)
    ax_xy.set_aspect("equal")
    ax_xy.set_xlabel("x (cm)")
    ax_xy.set_ylabel("y (cm)")
    ax_xy.legend(loc="lower left")
    for ax, lab in ((ax_xz, "x (cm)"), (ax_yz, "y (cm)")):
        ax.set_xlim(-ext, ext)
        ax.set_ylim(z_lo, z_hi)
        ax.set_xlabel(lab)
        ax.set_ylabel("z (cm)")
    flags = mc.flag_text(int(T["flags"][row]), mc.TRACK_FLAGS)
    lines = []
    for L in s.layers:
        k = int(T[f"paddle_{L.name}"][row])
        if k < 0:
            continue
        p = s.paddles[k]
        i = int(T[f"hit_{L.name}"][row])
        sig = mc.position_sigma(d.calib, p, s)
        um = T[f"u_meas_{L.name}"][row] if f"u_meas_{L.name}" in T else T[f"u_{L.name}"][row]
        lines.append(f"{L.name}: {p.name}  u = {T[f'u_{L.name}'][row]:+.1f} cm (measured {um:+.1f} ± {sig:.1f})  "
                     f"hA {d.H['hA'][i]:.0f}  hB {d.H['hB'][i]:.0f} ADC  tA-tB {d.H['tdiff_ns'][i]:+.2f} ns")
    ax_xy.text(0.0, -0.13, "\n".join(lines), transform=ax_xy.transAxes, fontsize=7, color=INK2, va="top",
               linespacing=1.5)
    t_s = T["t_ns"][row] * 1e-9
    sub = (f"θ = {T['theta_deg'][row]:.1f} ± {T['sigma_theta_deg'][row]:.1f}°, φ = {T['phi_deg'][row]:.0f}°, "
           f"path {T['s_cm'][row]:.1f} cm, tof {T['tof_ns'][row]:.2f} ns (s/c {T['tof_expected_ns'][row]:.2f}), "
           f"1/β = {T['inv_beta'][row]:.2f}; flags: {flags}")
    fig.text(0.125, 0.975, f"Event {e} of {d.short}, {t_s:.6f} s into the run", ha="left", va="top",
             fontsize=9.5, color=INK)
    fig.text(0.125, 0.94, sub, ha="left", va="top", fontsize=8, color=INK2)
    ax_xy.set_title("top view", loc="left", fontsize=8, color=INK2)
    ax_xz.set_title("side view from -y", loc="left", fontsize=8, color=INK2)
    ax_yz.set_title("side view from +x", loc="left", fontsize=8, color=INK2)
    mc.footer(fig, f"{d.run}  ·  setup {s.name}  ·  filled: two-end hit, dashed: one end only; "
                   f"bars: position resolution")
    path = os.path.join(out_dir, f"event_{e}.png")
    mc.save_figure(fig, path)
    return path


def waveform_figure(d, row, out_dir):
    """The event's waveforms: all of them (left) and their rising edges, each divided by its
    height (right), where the times of the 20 % crossings show the differences directly."""
    s = d.setup
    T = d.Tr
    e = int(T["event"][row])
    t_event = T["t_ns"][row]
    ps = d.pinfo["pulse_settings"]
    ref = int(d.pinfo["trigger_sample"])
    origin = int(d.pinfo["ttt_origin"])
    frac = float(ps["cfd_fraction"])
    mhz = float(ps.get("lowpass_mhz", 0) or 0)
    rows = []
    for i in np.where(d.H["event"] == e)[0]:
        for end in ("A", "B"):
            k = int(d.H[f"pulse{end}"][i])
            if k >= 0:
                rows.append((int(d.H["paddle"][i]), end, k, int(i)))
    rows.sort()
    if not rows or d.raw is None:
        return None
    from extract_pulses import lowpass
    fig = plt.figure(figsize=(11.6, 4.6))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.0, 1.0], wspace=0.22, top=0.84)
    ax, az = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    times, labels = [], []
    for n, (pi, end, k, hit) in enumerate(rows[:8]):
        chid = int(d.P["channel"][k])
        w = d.waveform(chid, d.P["entry"][k])
        if w is None:
            continue
        positive = d.pinfo["channels"].get(str(chid), {}).get("positive", True)
        sign = 1 if positive else -1
        y_raw = (w - d.P["baseline"][k]) * sign
        y = (lowpass(w[None, :], mhz, ps.get("lowpass_kind", "gauss"))[0] - d.P["baseline"][k]) * sign if mhz > 0 else y_raw
        t = ((int(d.P["ttt"][k]) - origin) * mc.NS_PER_TICK + (np.arange(len(w)) - ref) * mc.NS_PER_SAMPLE
             - s.time_offset(chid) - t_event)
        col = CATEGORICAL[n % len(CATEGORICAL)]
        p = s.paddles[pi]
        h = float(d.P["height"][k])
        tp = d.P["t_ns"][k] - t_event
        ax.plot(t, y, color=col, lw=0.9, label=f"{p.name}-{end}  ch{chid}  height {h:.0f}")
        ax.plot(tp, frac * h, "o", ms=5, color=col, mec=SURFACE, mew=1.2, zorder=4)
        # zoom: samples as dots, the (filtered) waveform as a line, both divided by the height
        az.plot(t, y_raw / h, "o", ms=2.6, color=col, alpha=0.45 if mhz > 0 else 0.9, mew=0)
        az.plot(t, y / h, color=col, lw=1.1)
        az.plot(tp, frac, "o", ms=6, color=col, mec=SURFACE, mew=1.2, zorder=5)
        az.axvline(tp, color=col, lw=0.7, alpha=0.6)
        times.append((p, end, tp))
    if not times:
        plt.close(fig)
        return None
    t0 = min(x[2] for x in times)
    t1 = max(x[2] for x in times)
    ax.set_xlim(t0 - 60, t0 + 340)
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.set_xlabel("time from the event's first hit (ns)")
    ax.set_ylabel("ADC above the baseline")
    mc.headroom(ax, 1.2)
    ax.legend(loc="upper right")
    ax.set_title("whole records (dots: pulse times)", loc="left", fontsize=8, color=INK2)
    az.set_xlim(t0 - 14, t1 + 10)
    az.set_ylim(-0.15, 1.15)
    az.axhline(frac, color=MUTED, lw=0.8, ls="--")
    az.axhline(0, color=AXIS, lw=0.8)
    az.set_xlabel("time from the event's first hit (ns)")
    az.set_ylabel("signal / pulse height")
    az.set_title(f"rising edges: dots = samples (2 ns), lines = {'the filtered waveform' if mhz > 0 else 'linear between samples'}; "
                 f"dashed: {frac:.0%}", loc="left", fontsize=8, color=INK2)
    # the differences the analysis uses
    lines = []
    by_paddle = {}
    for p, end, tp in times:
        by_paddle.setdefault(p.name, {})[end] = tp
    for name, ends in by_paddle.items():
        if "A" in ends and "B" in ends:
            lines.append(f"{name}: tA - tB = {ends['A'] - ends['B']:+.2f} ns")
    if np.isfinite(T["tof_ns"][row]):
        lines.append(f"time of flight {T['tof_ns'][row]:.2f} ns (s/c {T['tof_expected_ns'][row]:.2f})")
    az.text(0.98, 0.04, "\n".join(lines), transform=az.transAxes, ha="right", va="bottom", fontsize=7.5,
            color=INK2, linespacing=1.5, bbox=dict(boxstyle="round,pad=0.45", fc=SURFACE, ec=mc.GRID, lw=0.6))
    filt = f"{ps.get('lowpass_kind')} low-pass {mhz:g} MHz" if mhz > 0 else "no filter"
    fig.text(0.125, 0.975, f"Waveforms of event {e}, {d.short}", ha="left", va="top", fontsize=9.5, color=INK)
    fig.text(0.125, 0.935, f"each pulse's time = where its rising edge crosses {frac:.0%} of its height ({filt}); "
                           f"the ends of a paddle differ by -2u/v, the layers by the time of flight",
             ha="left", va="top", fontsize=8, color=INK2)
    mc.footer(fig, f"{d.run}  ·  setup {s.name}")
    path = os.path.join(out_dir, f"event_{e}_waveforms.png")
    mc.save_figure(fig, path)
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", default=None,
                    help="Analysis folder with tracks.root (asked for if left out)")
    ap.add_argument("--events", default=None, help="Event numbers, e.g. 12,57,300-310")
    ap.add_argument("--first", type=int, default=10, help="The first N tracks (default 10)")
    ap.add_argument("--all", action="store_true", help="Every track")
    ap.add_argument("--several", type=int, default=0,
                    help="Also the first N tracks with several paddles of one layer (default 0)")
    ap.add_argument("--no-waveforms", action="store_true", help="Leave out the waveform pictures")
    ap.add_argument("--setup", default=None, help="Setup file (default: the one step 1 used)")
    args = ap.parse_args()
    folder = mc.ask_analysis_folder(args.folder, mc.TRACKS_FILE, "reconstruct.py")
    events = parse_events(args.events) if args.events else None
    run(folder, events, None if args.all else args.first, not args.no_waveforms, args.setup, args.several)


def run(folder, events=None, first=10, waveforms=True, setup_path=None, several=0):
    """Event pictures for one analysis folder; returns the list of files."""
    d = Data(folder, setup_path)
    if events is not None:
        rows = [int(np.where(d.Tr["event"] == e)[0][0]) for e in events if (d.Tr["event"] == e).any()]
        missing = [e for e in events if not (d.Tr["event"] == e).any()]
        if missing:
            print(f"Note: no track in event(s) {', '.join(map(str, missing))}")
    else:
        rows = list(range(len(d.Tr["event"]) if first is None else min(first, len(d.Tr["event"]))))
    if several:
        odd = np.where((d.Tr["flags"] & (mc.TRACK_ADJACENT | mc.TRACK_AMBIGUOUS)) != 0)[0]
        rows += [int(r) for r in odd if int(r) not in rows][:several]
    out_dir = os.path.join(folder, mc.DISPLAY_DIR)
    os.makedirs(out_dir, exist_ok=True)
    if waveforms and d.raw is None:
        print(f"Note: raw file {d.pinfo.get('raw_file')} not found: no waveform pictures")
    made = []
    bar = mc.Progress("Event pictures", max(len(rows), 1))
    for row in rows:
        made.append(geometry_figure(d, row, out_dir))
        if waveforms:
            w = waveform_figure(d, row, out_dir)
            if w:
                made.append(w)
        bar.update(1, f"event {int(d.Tr['event'][row])}")
    bar.close()
    print(f"Saved {len(made)} pictures in {out_dir}")
    return made


if __name__ == "__main__":
    main()
