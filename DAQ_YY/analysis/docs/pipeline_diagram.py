#!/usr/bin/env python3
"""
pipeline_diagram.py  -  the paddle muon analysis on one slide

Draws the six steps of the analysis (extract_pulses.py ... event_display.py)
side by side. Every box carries the ID it has in the detailed block diagram of
../README.md ("Block diagram"), so either picture can be used to name a box
("in 3a, ..."). Writes pipeline.svg, pipeline.pdf and pipeline.png next to this
script: 16:9, 13.33 x 7.5 inches, the size of a widescreen slide.

    python3 docs/pipeline_diagram.py               # -> docs/pipeline.svg, .pdf, .png
    python3 docs/pipeline_diagram.py --out DIR --dpi 300

When a step changes, change its boxes in STEPS below and in the README diagram,
then run this again. Only matplotlib is needed (no ROOT). A label too long for
its box is reported (exit code 1): shorten it in STEPS.
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Rectangle  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# Colours of ../../scripts/daqyy_common.py (the palette of the analysis plots), so the slide matches them.
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
ACCENT = "#2a78d6"   # categorical slot 1: step numbers and box IDs
CARD = "#f2f1ec"     # background of a step
BOX = "#ffffff"

SANS = dict(fontfamily="DejaVu Sans", fontstretch="condensed")  # narrow, so that the labels fit
MONO = dict(fontfamily="DejaVu Sans Mono")

# ---- what the picture shows ----

INPUTS = {  # cylinders at the top
    "raw": ("raw ROOT file", "waveforms, 2 ns samples"),
    "conf": ("run .conf", "pulse polarity, thresholds"),
    "setup": ("setup file  setups/*.toml", "geometry, channel map, settings"),
}

# One entry per step: number, name, script, what it writes, what it takes from the setup file and the run's
# .conf (the dotted arrows), and its boxes: (ID, what, detail), with True added for a box that is simulation
# or only quoted, never applied to the data (drawn dashed). IDs as in the README diagram.
STEPS = [
    ("1", "Pulses", "extract_pulses.py", "pulses.root", "[pulse], .conf polarity", [
        ("1a", "Low-pass filter 80 MHz", "Butterworth, zero-phase"),
        ("1b", "Baseline and noise", "samples up to trigger − 20"),
        ("1c", "Height, 20 % crossing", "rise time, charge"),
        ("1d", "Shape tests → flags", "noise, no edge, low: left out"),
        ("1e", "Edge template", "per channel, 5000 clean pulses"),
        ("1f", "Pulse time", "edge fit, 20 % of fitted height"),
    ]),
    ("2", "Hits and events", "build_events.py", "events.root", "[[paddles]], [events]", [
        ("2a", "Pair A and B ends", "|tA − tB| ≤ L/v + 5 ns"),
        ("2b", "Hits", "tA − tB, mean time, hA, hB"),
        ("2c", "Events", "hits within 30 ns of the first"),
        ("2d", "Chance coincidences", "one layer shifted in time"),
    ]),
    ("3", "Calibration", "calibrate.py", "calibration.json", "[calibration], .conf thresholds", [
        ("3a", "Box fit of tA − tB", "v, offset, σ of each paddle"),
        ("3b", "Attenuation", "ln(hA/hB) against u: λ, gA/gB"),
        ("3c", "Muon peak", "√(hA·hB) of track hits, Moyal"),
        ("3d", "Second pass of 3a, 3b", "hits above 0.8 × muon peak"),
        ("3e", "β = 1 layer offsets", "quoted only, never applied", True),
        ("3f", "Neighbour times", "board-offset check, quoted only", True),
    ]),
    ("4", "Tracks", "reconstruct.py", "tracks.root", "geometry, [tracks], [flux]", [
        ("4a", "Point along paddle", "u from tA − tB, as measured"),
        ("4b", "Point across paddle", "centre line, error width/√12"),
        ("4c", "One point per layer", "neighbours within 5 ns merged"),
        ("4d", "Straight-line fit", "zenith ± error, azimuth"),
        ("4e", "Path, TOF, β", "no time offsets applied"),
        ("4f", "Toy cosmic muons", "I₀ cos²θ: expected rates, angles", True),
        ("4g", "Channel-map check", "data against the toy", True),
    ]),
    ("5", "Plots", "make_plots.py", "plots/, analysis_summary.txt", "geometry", [
        ("5a", "Plots, one per file", "3D tracks also as a GIF"),
        ("5b", "Run summary", "main numbers, steps 1–4"),
    ]),
    ("6", "Event pages", "event_display.py", "event_displays/", "[display]", [
        ("6a", "One page per track", "views, table, raw waveforms"),
    ]),
]

# ---- layout, in inches (1 data unit = 1 inch) ----

W, H = 13.333, 7.5          # a 16:9 slide
MARGIN, MARGIN_R = 0.2, 0.3   # the right one also holds the dotted line to step 6
COLS = 5                    # steps 5 and 6 share the last column: both start from tracks.root
GAP = 0.3                   # between the columns: the arrows
CW = (W - MARGIN - MARGIN_R - (COLS - 1) * GAP) / COLS
CYL_H, CYL_RY = 0.62, 0.065
BUS_Y = H - MARGIN - CYL_H  # the dotted line the settings run along
CARD_TOP = BUS_Y - 0.3
NAME_Y, SCRIPT_Y, OUT_Y, USES_Y, HEAD = 0.27, 0.53, 0.71, 0.89, 1.02   # below a card's top
BOX_H, BOX_GAP, PAD, INSET = 0.56, 0.08, 0.09, 0.08
STACK_GAP = 0.34            # between steps 5 and 6
LEGEND_Y = 0.36

FITS = []                   # (text, left, right, label): checked once everything is drawn


def text(ax, x, y, s, size, color=INK, weight="normal", font=SANS, ha="left", fit=None):
    """Text centred on y; fit = (left, right) edges it must stay within."""
    t = ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha, va="center", **font)
    if fit:
        FITS.append((t, fit[0], fit[1], s))
    return t


def arrow(ax, xy_from, xy_to, color=INK2, lw=1.4, scale=13):
    ax.annotate("", xy=xy_to, xytext=xy_from, arrowprops=dict(
        arrowstyle="-|>", color=color, lw=lw, mutation_scale=scale, shrinkA=0, shrinkB=0))


def dotted(ax, xs, ys):
    ax.plot(xs, ys, color=MUTED, lw=1.2, ls=(0, (0.1, 2.2)), dash_capstyle="round", solid_capstyle="round")


def dotted_arrow(ax, xy_from, xy_to):
    """Dotted line with a small solid head (a dotted head looks frayed)."""
    (x0, y0), (x1, y1) = xy_from, xy_to
    d = ((x1 - x0) ** 2 + (y1 - y0) ** 2) ** 0.5
    xm, ym = x1 - 0.09 * (x1 - x0) / d, y1 - 0.09 * (y1 - y0) / d
    dotted(ax, [x0, xm], [y0, ym])
    arrow(ax, (xm, ym), xy_to, color=MUTED, lw=1.2, scale=9)


def cylinder(ax, xc, w, title, detail):
    """An input file, standing on the line y = BUS_Y."""
    x0, y0, h, ry = xc - w / 2, BUS_Y, CYL_H, CYL_RY
    ax.add_patch(Ellipse((xc, y0 + ry), w, 2 * ry, fc=BOX, ec=AXIS, lw=0.9, zorder=2))
    ax.add_patch(Rectangle((x0, y0 + ry), w, h - 2 * ry, fc=BOX, ec="none", zorder=2.1))
    for x in (x0, x0 + w):
        ax.plot([x, x], [y0 + ry, y0 + h - ry], color=AXIS, lw=0.9, zorder=2.2)
    ax.add_patch(Ellipse((xc, y0 + h - ry), w, 2 * ry, fc=GRID, ec=AXIS, lw=0.9, zorder=2.3))
    yc = y0 + 0.5 * (h - ry)
    edges = (x0 + 0.06, x0 + w - 0.06)
    text(ax, xc, yc + 0.09, title, 10.5, INK, "bold", ha="center", fit=edges).set_zorder(3)
    text(ax, xc, yc - 0.1, detail, 8.8, INK2, ha="center", fit=edges).set_zorder(3)


def box(ax, x0, y0, w, entry):
    bid, what, detail = entry[:3]
    aside = len(entry) > 3 and entry[3]
    ax.add_patch(FancyBboxPatch((x0, y0), w, BOX_H, boxstyle="round,pad=0,rounding_size=0.06",
                                fc=SURFACE if aside else BOX, ec=MUTED if aside else GRID, lw=0.9,
                                ls=(0, (3.5, 2.2)) if aside else "-"))
    right = x0 + w - 0.07
    text(ax, x0 + 0.09, y0 + 0.66 * BOX_H, bid, 9.5, ACCENT, "bold", MONO)
    text(ax, x0 + 0.33, y0 + 0.66 * BOX_H, what, 10.5, INK2 if aside else INK, "bold", fit=(x0, right))
    text(ax, x0 + 0.33, y0 + 0.29 * BOX_H, detail, 8.8, INK2, fit=(x0, right))


def card(ax, x0, top, step):
    """One step: header, then its boxes top to bottom. Returns the card's bottom."""
    num, name, script, out, uses, boxes = step
    n = len(boxes)
    h = HEAD + n * BOX_H + (n - 1) * BOX_GAP + PAD
    ax.add_patch(FancyBboxPatch((x0, top - h), CW, h, boxstyle="round,pad=0,rounding_size=0.1", fc=CARD, ec="none"))
    right = x0 + CW - 0.08
    ax.add_patch(Circle((x0 + 0.27, top - NAME_Y), 0.15, fc=ACCENT, ec="none"))
    text(ax, x0 + 0.27, top - NAME_Y, num, 11, "white", "bold", ha="center")
    text(ax, x0 + 0.5, top - NAME_Y, name, 13, INK, "bold", fit=(x0, right))
    text(ax, x0 + 0.1, top - SCRIPT_Y, script, 8.5, INK2, font=MONO, fit=(x0, right))
    text(ax, x0 + 0.1, top - OUT_Y, "→ " + out, 8.5, INK2, font=MONO, fit=(x0, right))
    text(ax, x0 + 0.1, top - USES_Y, uses, 8.5, MUTED, fit=(x0, right))
    y = top - HEAD
    for entry in boxes:
        box(ax, x0 + INSET, y - BOX_H, CW - 2 * INSET, entry)
        y -= BOX_H + BOX_GAP
    return top - h


def legend(ax, renderer):
    """One line at the bottom: what the three kinds of lines and boxes mean."""
    x, y = MARGIN + 0.05, LEGEND_Y

    def label(x, s):
        t = text(ax, x, y, s, 8.8, INK2)
        return ax.transData.inverted().transform(t.get_window_extent(renderer))[1, 0] + 0.45

    arrow(ax, (x, y), (x + 0.42, y))
    x = label(x + 0.52, "data files, from step to step")
    dotted_arrow(ax, (x, y), (x + 0.42, y))
    x = label(x + 0.52, "settings: sections of the setup file [ ], and the run's .conf")
    ax.add_patch(FancyBboxPatch((x, y - 0.1), 0.42, 0.2, boxstyle="round,pad=0,rounding_size=0.04",
                                fc=SURFACE, ec=MUTED, lw=0.9, ls=(0, (3.5, 2.2))))
    label(x + 0.52, "simulation or quoted only: never applied to the data")


def draw():
    fig = plt.figure(figsize=(W, H), facecolor=SURFACE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")
    x0 = [MARGIN + k * (CW + GAP) for k in range(COLS)]

    # the steps; 5 and 6 one above the other in the last column
    for k, step in enumerate(STEPS[:5]):
        bottom = card(ax, x0[k], CARD_TOP, step)
    y6 = bottom - STACK_GAP - NAME_Y   # bottom: of step 5, the last card drawn; y6: step 6's name line
    card(ax, x0[4], y6 + NAME_Y, STEPS[5])

    # data flow: 1 -> 2 -> 3 -> 4 -> 5, and 4 -> 6
    for k in range(4):
        arrow(ax, (x0[k] + CW, CARD_TOP - NAME_Y), (x0[k + 1], CARD_TOP - NAME_Y))
    arrow(ax, (x0[3] + CW, y6), (x0[4], y6))

    # inputs: the raw file into step 1; the .conf and the setup file on the dotted line into every step
    raw_x = x0[0] + 0.38 * CW
    drops = [x0[0] + 0.82 * CW] + [x0[k] + 0.5 * CW for k in range(1, 5)]
    channel = x0[4] + CW + 0.12      # the line to step 6 runs down beside step 5
    cylinder(ax, raw_x, 1.75, *INPUTS["raw"])
    arrow(ax, (raw_x, BUS_Y), (raw_x, CARD_TOP))
    dotted(ax, [drops[0], channel], [BUS_Y, BUS_Y])
    cylinder(ax, 0.5 * (drops[0] + drops[1]), 1.72, *INPUTS["conf"])
    cylinder(ax, 0.5 * (drops[2] + drops[3]), 2.35, *INPUTS["setup"])
    for xd in drops:
        dotted_arrow(ax, (xd, BUS_Y), (xd, CARD_TOP))
    dotted(ax, [channel, channel], [BUS_Y, y6])
    dotted_arrow(ax, (channel, y6), (x0[4] + CW, y6))

    renderer = fig.canvas.get_renderer()
    legend(ax, renderer)
    return fig, ax


def too_long(fig, ax):
    """Labels that run past their box, as lines of text."""
    fig.canvas.draw()
    renderer, inv = fig.canvas.get_renderer(), ax.transData.inverted()
    out = []
    for t, left, right, s in FITS:
        (x0, _), (x1, _) = inv.transform(t.get_window_extent(renderer))
        over = max(x1 - right, left - x0)
        if over > 0.005:
            out.append(f"  '{s}': {over:.2f} in too wide")
    return out


def main():
    ap = argparse.ArgumentParser(description="Draw the analysis pipeline on one 16:9 slide "
                                             "(pipeline.svg, .pdf, .png).")
    ap.add_argument("--out", default=HERE, help="folder for the three files (default: next to this script)")
    ap.add_argument("--dpi", type=int, default=300, help="resolution of the PNG (default 300: 4000 x 2250 pixels)")
    args = ap.parse_args()

    plt.rcParams.update({"svg.hashsalt": "pipeline", "pdf.fonttype": 42})  # same files from the same picture
    fig, ax = draw()
    bad = too_long(fig, ax)
    os.makedirs(args.out, exist_ok=True)
    for ext, meta in (("svg", {"Date": None}), ("pdf", {"CreationDate": None}), ("png", {})):
        path = os.path.join(args.out, f"pipeline.{ext}")
        fig.savefig(path, dpi=args.dpi, facecolor=SURFACE, metadata=meta)
        print("wrote", path)
    if bad:
        print("Labels too long for their boxes (shorten them in STEPS):", *bad, sep="\n")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
