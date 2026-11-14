#!/usr/bin/env python3
"""
event_display.py  -  pictures of single events of the paddle muon analysis

For every chosen track it draws one page, event_<n>.png in event_displays/ of the
analysis folder (n = the event number of events.root):
    top view and two side views   the paddles (the ones that fired filled), each paddle
                                  end's readout as a small box coloured by its pulse
                                  height in this event (the light sharing at a glance),
                                  every hit at its own paddle with its ±1 sigma error as
                                  the track fit used it (open marker: where tA - tB alone
                                  puts the hit, when that is beyond the paddle's end), the
                                  fitted track as an arrow (the way the particle went),
                                  and a side of the room marked if the setup's [display]
                                  names one (the wall of the first setup)
    info table                    the track's numbers and every hit's
    waveforms                     the event's records on one time axis, and their rising
                                  edges divided by the height, where the time differences
                                  show, with the samples each edge fit used and the fitted
                                  template (needs the raw ROOT file)

Usage:
    python3 event_display.py                       # asks for the analysis folder; the first 10 tracks
    python3 event_display.py FOLDER --first 20
    python3 event_display.py FOLDER --events 12,57,301
    python3 event_display.py FOLDER --all            # every track (also --all-events)
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
        f = mc.open_root(pp)
        has_templates = bool(f.Get("templates"))
        f.Close()
        fit = has_templates and self.pinfo["pulse_settings"].get("timing") == "fit"
        self.P = mc.read_tree(pp, "pulses", ["channel", "entry", "ttt", "t_ns", "edge", "baseline", "height", "flags"]
                              + (["edge_cfd", "fit_amp"] if fit else []))
        self.templates = {}
        if fit:
            tm = mc.read_tree(pp, "templates")
            for chid in np.unique(tm["channel"]):
                sel = tm["channel"] == chid
                self.templates[int(chid)] = (tm["x_ns"][sel], tm["y"][sel])
        self.H = mc.read_tree(ep, "hits")
        self.Tr = mc.read_tree(tp, "tracks")
        f = mc.open_root(tp)
        has_points = bool(f.Get("points"))
        f.Close()
        self.Pts = mc.read_tree(tp, "points") if has_points else None
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


READOUT_LEN = 2.0      # cm: the readout box drawn beyond each paddle end
VIEW_MARGIN = 0.8      # cm: a hit measured beyond the view is drawn this far inside its edge


def end_heights(d, e):
    """{(paddle index, end): pulse height} of every pulse of event e (two-end and one-end hits)."""
    out = {}
    for i in np.where(d.H["event"] == e)[0]:
        for end in ("A", "B"):
            if int(d.H[f"pulse{end}"][i]) >= 0:
                out[(int(d.H["paddle"][i]), end)] = float(d.H[f"h{end}"][i])
    return out


def paddle_box(p, view):
    """(x0, y0, width, height) of the paddle in a view: 'xy' from above, 'xz' or 'yz' from the side."""
    if view == "xy":
        return p.x - 0.5 * p.size_x, p.y - 0.5 * p.size_y, p.size_x, p.size_y
    if view == "xz":
        return p.x - 0.5 * p.size_x, p.z - 0.5 * p.thickness, p.size_x, p.thickness
    return p.y - 0.5 * p.size_y, p.z - 0.5 * p.thickness, p.size_y, p.thickness


def readout_box(p, end, view, aside=0):
    """(x0, y0, w, h) of an end's readout box beyond the paddle, and its label's anchor (x, y, alignment;
    top view only); None when the view looks along the paddle (its ends point at the viewer).
    aside = +1 or -1: the label sits next to the box on that side across the paddle, which leaves the
    space beyond the box to a hit drawn there."""
    sign = p.a_sign if end == "A" else -p.a_sign          # the end sits at + or - along the axis
    half, size = 0.5 * p.length, READOUT_LEN
    if view == "xy":   # label just beyond the box, growing away from the end (or next to the box)
        if p.axis == "x":
            x0 = p.x + sign * half + (0 if sign > 0 else -size)
            box = (x0, p.y - 0.5 * p.width, size, p.width)
            if aside:
                return box, (x0 + 0.5 * size, p.y + aside * (0.5 * p.width + 0.3),
                             dict(ha="center", va="bottom" if aside > 0 else "top"))
            return box, (p.x + sign * (half + size + 0.3), p.y, dict(ha="left" if sign > 0 else "right", va="center"))
        y0 = p.y + sign * half + (0 if sign > 0 else -size)
        box = (p.x - 0.5 * p.width, y0, p.width, size)
        if aside:
            return box, (p.x + aside * (0.5 * p.width + 0.3), y0 + 0.5 * size,
                         dict(ha="left" if aside > 0 else "right", va="center", rotation=90))
        return box, (p.x, p.y + sign * (half + size + 0.3),
                     dict(ha="center", va="bottom" if sign > 0 else "top", rotation=90))
    horizontal = "x" if view == "xz" else "y"
    if p.axis != horizontal:
        return None
    c = p.x if horizontal == "x" else p.y
    x0 = c + sign * half + (0 if sign > 0 else -size)
    return (x0, p.z - 0.5 * size, size, size), None


def outer_side(p):
    """+1 or -1: the side across the paddle away from the other paddles of its layer."""
    others = [q.y if p.axis == "x" else q.x for q in p.layer.paddles if q is not p]
    own = p.y if p.axis == "x" else p.x
    return 1 if not others or own >= float(np.mean(others)) else -1


def hit_points(d, row):
    """Every two-end hit of the track's layers at its own paddle: (x, y, z, sx, sy, paddle, hit, u, u_meas),
    with the errors the track fit used."""
    s = d.setup
    T = d.Tr
    pts = []
    if d.Pts is not None:
        for j in np.where(d.Pts["track"] == row)[0]:
            p = s.paddles[int(d.Pts["paddle"][j])]
            # along the paddle (reconstruct.py): the measurement's Gaussian cut to the paddle, smaller
            # than the resolution near an end; across: the paddle's width / sqrt 12
            _, s_u = mc.position_in_paddle(d.Pts["u_meas"][j], 0.5 * p.length, mc.position_sigma(d.calib, p, s))
            across = p.width / math.sqrt(12)
            sx, sy = (float(s_u), across) if p.axis == "x" else (across, float(s_u))
            pts.append((d.Pts["x"][j], d.Pts["y"][j], d.Pts["z"][j], sx, sy, p, int(d.Pts["hit"][j]),
                        d.Pts["u"][j], d.Pts["u_meas"][j]))
    else:
        for L in s.layers:
            k = int(T[f"paddle_{L.name}"][row])
            if k >= 0:
                p = s.paddles[k]
                pts.append((T[f"x_{L.name}"][row], T[f"y_{L.name}"][row], T[f"z_{L.name}"][row],
                            T[f"sx_{L.name}"][row], T[f"sy_{L.name}"][row], p, int(T[f"hit_{L.name}"][row]),
                            T[f"u_{L.name}"][row], T[f"u_meas_{L.name}"][row]))
    return pts


def draw_views(d, row, axes, heights, norm):
    s = d.setup
    T = d.Tr
    e = int(T["event"][row])
    pts = hit_points(d, row)
    # a hit that tA - tB alone puts beyond its paddle's end is drawn there too; that end's label moves aside
    beyond = {(p.index, "A" if um > 0 else "B") for x, y, z, sx, sy, p, hit, u, um in pts
              if np.isfinite(um) and abs(um) > 0.5 * p.length}
    ab = {}
    for i in np.where(d.H["event"] == e)[0]:
        k = int(d.H["paddle"][i])
        ab[k] = ab.get(k, False) or d.H["kind"][i] == mc.HIT_AB
    ext = max(max(abs(p.x) + 0.5 * p.size_x, abs(p.y) + 0.5 * p.size_y) for p in s.paddles) + READOUT_LEN + 7.5
    zs = [p.z for p in s.paddles]
    z_lo, z_hi = min(zs) - 9, max(zs) + 9
    zz = np.array([z_hi, z_lo])
    lx = T["ax"][row] + T["bx"][row] * (zz - T["z0"][row])
    ly = T["ay"][row] + T["by"][row] * (zz - T["z0"][row])
    for ax, view in zip(axes, ("xy", "xz", "yz")):
        for p in s.paddles:
            col = mc.layer_color(p.layer.index)
            x0, y0, w, h = paddle_box(p, view)
            fired = p.index in ab
            ax.add_patch(Rectangle((x0, y0), w, h, facecolor=col if fired else mc.GRID,
                                   alpha=(0.28 if ab.get(p.index) else 0.12) if fired else 0.55, lw=0, zorder=1))
            ax.add_patch(Rectangle((x0, y0), w, h, fill=False, edgecolor=col if fired else AXIS, zorder=2,
                                   lw=1.1 if fired else 0.7, ls="-" if (not fired or ab.get(p.index)) else (0, (3, 2))))
            for end in ("A", "B"):
                box = readout_box(p, end, view, outer_side(p) if (p.index, end) in beyond else 0)
                if box is None:
                    continue
                (bx0, by0, bw, bh), label = box
                hgt = heights.get((p.index, end))
                face = mc.SEQUENTIAL(norm(hgt)) if hgt is not None else SURFACE
                ax.add_patch(Rectangle((bx0, by0), bw, bh, facecolor=face, edgecolor=INK2 if hgt is not None else AXIS,
                                       lw=0.7, zorder=3))
                if view == "xy":
                    tx, ty, align = label
                    ax.text(tx, ty, f"{end} {hgt:.0f}" if hgt is not None else end, fontsize=6,
                            color=INK2 if hgt is not None else MUTED, zorder=3, **align)
        # the track: a line through the detector with an arrow the way the particle went (down)
        if view == "xy":
            start, stop = (lx[0], ly[0]), (lx[1], ly[1])
        elif view == "xz":
            start, stop = (lx[0], zz[0]), (lx[1], zz[1])
        else:
            start, stop = (ly[0], zz[0]), (ly[1], zz[1])
        ax.annotate("", xy=stop, xytext=start, zorder=6,
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.3, mutation_scale=11, shrinkA=0, shrinkB=0))
    ax_xy, ax_xz, ax_yz = axes
    seen = set()
    shown = {"o": False, "edge": False}
    for x, y, z, sx, sy, p, hit, u, um in pts:
        col = mc.layer_color(p.layer.index)
        kw = dict(fmt="o", ms=6.5, color=col, mec=SURFACE, mew=1.4, elinewidth=1.6, capsize=0, zorder=7)
        ax_xy.errorbar(x, y, xerr=sx, yerr=sy, label=None if p.layer.name in seen else f"{p.layer.name} layer hit", **kw)
        seen.add(p.layer.name)
        ax_xz.errorbar(x, z, xerr=sx, **kw)
        ax_yz.errorbar(y, z, xerr=sy, **kw)
        if np.isfinite(um) and abs(um) > 0.5 * p.length:     # the timing alone puts it beyond the paddle's end
            c = p.x if p.axis == "x" else p.y
            lim = ext - VIEW_MARGIN
            pos = c + p.a_sign * um
            pulled = abs(pos) > lim                          # further out than the view: drawn at its edge
            pos = min(max(pos, -lim), lim)
            out = 1 if pos > c else -1
            shown["edge" if pulled else "o"] = True
            mk = dict(ms=5.5, mfc="none", mec=col, mew=1.2, ls="none", zorder=6)
            xm, ym = (pos, p.y) if p.axis == "x" else (p.x, pos)
            ax_xy.plot([x, xm], [y, ym], color=col, lw=0.8, ls=":", zorder=6)
            ax_xy.plot(xm, ym, marker=("o" if not pulled else (">" if out > 0 else "<") if p.axis == "x"
                                       else ("^" if out > 0 else "v")), **mk)
            side, a = (ax_xz, x) if p.axis == "x" else (ax_yz, y)   # the side view along the paddle
            side.plot([a, pos], [z, z], color=col, lw=0.8, ls=":", zorder=6)
            side.plot(pos, z, marker=("o" if not pulled else ">" if out > 0 else "<"), **mk)
    # the room side (the wall of the first setup)
    side, label = s.display.get("side", ""), s.display.get("side_label", "")
    if side in ("+x", "-x", "+y", "-y") and label:
        sgn = 1 if side[0] == "+" else -1
        band = 1.6
        if side[1] == "y":
            ax_xy.add_patch(Rectangle((-ext, sgn * ext - (band if sgn > 0 else 0)), 2 * ext, band, color=AXIS, lw=0, zorder=1))
            ax_xy.text(0, sgn * (ext - 0.5 * band), label, ha="center", va="center", fontsize=7.5, color=INK, zorder=2)
            ax_yz.axvspan(sgn * ext - (band if sgn > 0 else 0), sgn * ext + (0 if sgn > 0 else band), color=AXIS, lw=0)
            ax_yz.text(sgn * (ext - 0.5 * band), 0.5 * (z_lo + z_hi), label, rotation=90, ha="center", va="center",
                       fontsize=7.5, color=INK)
        else:
            ax_xy.add_patch(Rectangle((sgn * ext - (band if sgn > 0 else 0), -ext), band, 2 * ext, color=AXIS, lw=0, zorder=1))
            ax_xy.text(sgn * (ext - 0.5 * band), 0, label, rotation=90, ha="center", va="center", fontsize=7.5, color=INK)
            ax_xz.axvspan(sgn * ext - (band if sgn > 0 else 0), sgn * ext + (0 if sgn > 0 else band), color=AXIS, lw=0)
            ax_xz.text(sgn * (ext - 0.5 * band), 0.5 * (z_lo + z_hi), label, rotation=90, ha="center", va="center",
                       fontsize=7.5, color=INK)
    ax_xy.set_xlim(-ext, ext)
    ax_xy.set_ylim(-ext, ext)
    ax_xy.set_aspect("equal")
    ax_xy.set_xlabel("x (cm)")
    ax_xy.set_ylabel("y (cm)")
    ax_xy.grid(False)
    handles, labels = ax_xy.get_legend_handles_labels()
    from matplotlib.lines import Line2D
    for show, marker, text in ((shown["o"] or shown["edge"], "o", "from tA - tB alone (beyond the end)"),
                               (shown["edge"], ">", "the same, beyond the view")):
        if show:
            handles.append(Line2D([], [], marker=marker, ms=5.5, mfc="none", mec=INK2, mew=1.2, ls="none"))
            labels.append(text)
    ax_xy.legend(handles, labels, loc="lower left", fontsize=7)
    ax_xy.set_title("seen from above", loc="left", fontsize=8.5, color=INK2)
    for ax, lab, name in ((ax_xz, "x (cm)", "seen from the side (from -y)"), (ax_yz, "y (cm)", "seen from the side (from +x)")):
        ax.set_xlim(-ext, ext)
        ax.set_ylim(z_lo, z_hi)
        ax.set_xlabel(lab)
        ax.set_ylabel("z (cm)")
        ax.grid(False)
        ax.set_title(name, loc="left", fontsize=8.5, color=INK2)
        for layer in s.layers:
            ax.text(ext - 0.6, layer.z + 1.2, layer.name, ha="right", va="bottom", fontsize=7, color=INK2)
    return pts


def info_table(d, row, ax, pts, norm):
    s = d.setup
    T = d.Tr
    ax.axis("off")
    fl = mc.flag_text(int(T["flags"][row]), mc.TRACK_FLAGS)
    rows = [("zenith θ", f"{T['theta_deg'][row]:.1f} ± {T['sigma_theta_deg'][row]:.1f}°"),
            ("azimuth φ (arrival)", f"{T['phi_deg'][row]:.0f}°"),
            ("path between layers", f"{T['s_cm'][row]:.1f} cm"),
            ("time of flight", f"{T['tof_ns'][row]:.2f} ns"),
            ("light would need", f"{T['tof_expected_ns'][row]:.2f} ns"),
            ("β (no time offsets)", f"{T['beta'][row]:.2f}"),
            ("track flags", fl)]
    lines = [f"{k:<21}{v}" for k, v in rows]
    lines += ["", f"{'hit':<6}{'u ± σ cm':>11}{'meas.':>7}{'hA':>6}{'hB':>6}{'tA-tB':>7}"]
    for x, y, z, sx, sy, p, i, u, um in pts:
        su = sx if p.axis == "x" else sy
        lines.append(f"{p.name:<6}{f'{u:+.1f} ± {su:.1f}':>11}{um:>+7.1f}{d.H['hA'][i]:>6.0f}{d.H['hB'][i]:>6.0f}"
                     f"{d.H['tdiff_ns'][i]:>+7.2f}")
    lines += ["", "u ± σ: along the paddle, + towards end A:", "   best estimate inside the paddle and",
              "   its error, as the track fit used them", "meas.: from tA - tB alone",
              "hA, hB: pulse heights (ADC)", "tA-tB in ns"]
    ax.text(0.0, 1.0, "\n".join(lines), transform=ax.transAxes, ha="left", va="top", fontsize=7.4, color=INK,
            family="monospace", linespacing=1.45)
    cax = ax.inset_axes([0.02, 0.05, 0.9, 0.035])
    from matplotlib.colorbar import ColorbarBase
    cb = ColorbarBase(cax, cmap=mc.SEQUENTIAL, norm=norm, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.set_label("readout boxes: pulse height in this event (ADC)", color=INK2, fontsize=7)
    cax.tick_params(labelsize=6.5)


def draw_waveforms(d, row, ax, az):
    """The event's records (left) and their rising edges divided by the height (right)."""
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
                rows.append((int(d.H["paddle"][i]), end, k))
    rows.sort()
    from extract_pulses import lowpass
    from matplotlib.lines import Line2D
    times, handles = [], []
    for n, (pi, end, k) in enumerate(rows[:8]):
        chid = int(d.P["channel"][k])
        w = d.waveform(chid, d.P["entry"][k])
        if w is None:
            continue
        sign = 1 if d.pinfo["channels"].get(str(chid), {}).get("positive", True) else -1
        y_raw = (w - d.P["baseline"][k]) * sign
        y = ((lowpass(w[None, :], mhz, ps.get("lowpass_kind", "gauss"))[0] - d.P["baseline"][k]) * sign
             if mhz > 0 else y_raw)
        t = ((int(d.P["ttt"][k]) - origin) * mc.NS_PER_TICK + (np.arange(len(w)) - ref) * mc.NS_PER_SAMPLE
             - s.time_offset(chid) - t_event)
        col = CATEGORICAL[n % len(CATEGORICAL)]
        p = s.paddles[pi]
        h = float(d.P["height"][k])
        tp = d.P["t_ns"][k] - t_event
        amp = float(d.P["fit_amp"][k]) if "fit_amp" in d.P else float("nan")
        fitted = np.isfinite(amp) and chid in d.templates
        level = frac * (amp if fitted else h)          # the pulse time is where the (fitted) edge reaches this
        mark = "D" if fitted else "o"                  # the same marker for the pulse time in both panels
        ax.plot(t, y, color=col, lw=0.9)
        ax.plot(tp, level, mark, ms=5, color=col, mec=SURFACE, mew=1.2, zorder=4)
        handles.append(Line2D([], [], color=col, lw=1.2, marker=mark, ms=5, mec=SURFACE, mew=1.0,
                              label=f"{p.name}-{end}  ch{chid}  {h:.0f} ADC"))
        az.plot(t, y_raw / h, "o", ms=2.6, color=col, alpha=0.45 if mhz > 0 else 0.9, mew=0)
        if fitted:
            j = int(math.floor(d.P["edge_cfd"][k])) + np.arange(int(ps["fit_window"][0]), int(ps["fit_window"][1]) + 1)
            j = j[(j >= 0) & (j < len(t))]
            tx, ty = d.templates[chid]
            tt = np.linspace(t[j[0]], t[j[-1]], 200)
            az.plot(t, y / h, color=col, lw=0.8, alpha=0.75)
            az.plot(tt, amp * np.interp(tt - tp, tx, ty, left=0.0, right=ty[-1]) / h, color=col, lw=2.2,
                    alpha=0.5, solid_capstyle="round", zorder=3)
            az.plot(t[j], y[j] / h, "o", ms=3.8, color=col, mec=SURFACE, mew=0.6, zorder=4)
        else:
            az.plot(t, y / h, color=col, lw=1.2)
        az.plot(tp, level / h, mark, ms=6, color=col, mec=SURFACE, mew=1.2, zorder=5)
        az.axvline(tp, color=col, lw=0.7, alpha=0.6)
        times.append((p, end, tp))
    if not times:
        for a in (ax, az):
            a.text(0.5, 0.5, "raw waveforms not available", transform=a.transAxes, ha="center", color=MUTED)
        return
    t0, t1 = min(x[2] for x in times), max(x[2] for x in times)
    ax.set_xlim(t0 - 60, t0 + 340)
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.set_xlabel("time from the event's first hit (ns)")
    ax.set_ylabel("ADC above the baseline")
    mc.headroom(ax, 1.2)
    legend = dict(handles=handles, fontsize=6.5, ncol=2 if len(handles) > 4 else 1)   # the same in both panels
    ax.legend(loc="upper right", **legend)
    az.legend(loc="upper left", **legend)
    filt = f"{ps.get('lowpass_kind')} low-pass {mhz:g} MHz" if mhz > 0 else "no filter"
    ax.set_title(f"the records ({filt}); {'◆' if d.templates else '●'} = pulse time of each channel",
                 loc="left", fontsize=8.5, color=INK2)
    az.set_xlim(t0 - 14, t1 + 10)
    az.set_ylim(-0.15, 1.15)
    az.axhline(frac, color=MUTED, lw=0.8, ls="--")
    az.axhline(0, color=AXIS, lw=0.8)
    az.set_xlabel("time from the event's first hit (ns)")
    az.set_ylabel("signal / pulse height")
    if d.templates:
        az.set_title(f"rising edges / height: thin = {'filtered' if mhz > 0 else 'samples'}, big dots = the samples "
                     f"each fit uses, thick = fitted template;\n◆ = pulse time ({frac:.0%} of the fitted height); "
                     f"faint dots = raw samples (2 ns)", loc="left", fontsize=8, color=INK2)
    else:
        az.set_title(f"rising edges: dots = samples (2 ns), lines = {'filtered' if mhz > 0 else 'between samples'}; "
                     f"dashed: {frac:.0%}, where each pulse's time is taken", loc="left", fontsize=8.5, color=INK2)
    by_paddle = {}
    for p, end, tp in times:
        by_paddle.setdefault(p.name, {})[end] = tp
    lines = [f"{name}: tA - tB = {ends['A'] - ends['B']:+.2f} ns" for name, ends in by_paddle.items()
             if "A" in ends and "B" in ends]
    if np.isfinite(T["tof_ns"][row]):
        lines.append(f"time of flight {T['tof_ns'][row]:.2f} ns (light {T['tof_expected_ns'][row]:.2f})")
    az.text(0.98, 0.04, "\n".join(lines), transform=az.transAxes, ha="right", va="bottom", fontsize=7.5,
            color=INK2, linespacing=1.5, bbox=dict(boxstyle="round,pad=0.45", fc=SURFACE, ec=mc.GRID, lw=0.6))


def event_page(d, row, out_dir, waveforms=True):
    """One picture per event: views, info table, waveforms."""
    from matplotlib.colors import Normalize
    s = d.setup
    T = d.Tr
    e = int(T["event"][row])
    heights = end_heights(d, e)
    norm = Normalize(vmin=0, vmax=max(max(heights.values()) if heights else 1, 1))
    fig = plt.figure(figsize=(14.0, 9.6))
    gs = GridSpec(3, 3, figure=fig, width_ratios=[1.3, 1.0, 0.78], height_ratios=[1, 1, 1.08],
                  hspace=0.5, wspace=0.24, top=0.9, bottom=0.07, left=0.05, right=0.985)
    ax_xy = fig.add_subplot(gs[0:2, 0])
    ax_xz = fig.add_subplot(gs[0, 1])
    ax_yz = fig.add_subplot(gs[1, 1])
    ax_info = fig.add_subplot(gs[0:2, 2])
    ax_w = fig.add_subplot(gs[2, 0])
    ax_z = fig.add_subplot(gs[2, 1:3])
    pts = draw_views(d, row, (ax_xy, ax_xz, ax_yz), heights, norm)
    info_table(d, row, ax_info, pts, norm)
    if waveforms and d.raw is not None:
        draw_waveforms(d, row, ax_w, ax_z)
    else:
        for a in (ax_w, ax_z):
            a.axis("off")
            a.text(0.5, 0.5, "waveforms not drawn" + (" (no raw file)" if d.raw is None else ""),
                   transform=a.transAxes, ha="center", color=MUTED)
    t_s = T["t_ns"][row] * 1e-9
    fig.text(0.05, 0.975, f"Event {e}  ·  {d.short}  ·  {t_s:.6f} s into the run", ha="left", va="top",
             fontsize=12, color=INK)
    fig.text(0.05, 0.945, f"zenith {T['theta_deg'][row]:.1f} ± {T['sigma_theta_deg'][row]:.1f}°, "
                          f"azimuth {T['phi_deg'][row]:.0f}°, path "
                          f"{T['s_cm'][row]:.1f} cm, time of flight {T['tof_ns'][row]:.2f} ns, β {T['beta'][row]:.2f}"
                          f"  ·  arrow: the fitted track, the way the particle went",
             ha="left", va="top", fontsize=8.5, color=INK2)
    mc.footer(fig, f"{d.run}  ·  setup {s.name}  ·  filled: paddle seen at both ends (dashed outline: one end)\n"
                   f"bars: ±1σ of each hit as the track fit used it (along the paddle: from the timing, cut to the "
                   f"paddle, so smaller near an end; across: width / √12)  ·  open marker: where tA - tB alone puts "
                   f"a hit that lies beyond its paddle's end")
    path = os.path.join(out_dir, f"event_{e}.png")
    mc.save_figure(fig, path)
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", default=None,
                    help="Analysis folder with tracks.root (asked for if left out)")
    ap.add_argument("--events", default=None, help="Event numbers, e.g. 12,57,300-310")
    ap.add_argument("--first", type=int, default=10, help="The first N tracks (default 10)")
    ap.add_argument("--all", "--all-events", dest="all", action="store_true",
                    help="Every track (a 30-minute run: ~640 pages, ~13 minutes, ~200 MB)")
    ap.add_argument("--several", type=int, default=0,
                    help="Also the first N tracks with several paddles of one layer (default 0)")
    ap.add_argument("--no-waveforms", action="store_true", help="Leave out the waveforms")
    ap.add_argument("--setup", default=None, help="Setup file (default: the one step 1 used)")
    args = ap.parse_args()
    folder = mc.ask_analysis_folder(args.folder, mc.TRACKS_FILE, "reconstruct.py")
    events = parse_events(args.events) if args.events else None
    run(folder, events, None if args.all else args.first, not args.no_waveforms, args.setup, args.several)


def run(folder, events=None, first=10, waveforms=True, setup_path=None, several=0):
    """Event pages for one analysis folder; returns the list of files."""
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
    if len(rows) > 50:
        print(f"Drawing {len(rows)} event pages: about {len(rows) * 1.2 / 60:.0f} minutes and "
              f"{len(rows) * 0.32:.0f} MB")
    if waveforms and d.raw is None:
        print(f"Note: raw file {d.pinfo.get('raw_file')} not found: pages without waveforms")
    made = []
    bar = mc.Progress("Event pages", max(len(rows), 1))
    for row in rows:
        made.append(event_page(d, row, out_dir, waveforms))
        bar.update(1, f"event {int(d.Tr['event'][row])}")
    bar.close()
    print(f"Saved {len(made)} pages in {out_dir}")
    return made


if __name__ == "__main__":
    main()
