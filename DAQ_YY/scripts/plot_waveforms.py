#!/usr/bin/env python3
"""
plot_waveforms.py

Raw waveforms of one channel of the AD9694 boards, from a ROOT file written
by raw2root_new: one row per trigger, with the whole waveform and zooms on
its first and last 100 samples. Made to check short waveforms (the SHORT_WF
flag of analyze_baseline.py): where they stop, and whether the samples they
hold look like normal pedestal.

Usage:
    python3 plot_waveforms.py                  # asks for everything
    python3 plot_waveforms.py run.root --board 1 --channel 0
    python3 plot_waveforms.py run.root --board 1 --channel 0 --triggers all --out DIR

Leave out the ROOT file, --board or --channel and the script asks for it.
With no ROOT file on the command line it also asks for the output folder and
the triggers, and offers another channel after each plot. Enter keeps the
value shown in [brackets]; Tab completes paths.

--triggers: "short" (the default: every short, long or missing waveform, plus
the first full-length one to compare with), "all", or trigger numbers such
as "0,2,5-7". Trigger numbers count the triggers of the board, as in
analyze_baseline.py. The full length is the most common waveform length
among the channels of the board, unless --nsamps gives it.

Before plotting, the script prints one line per trigger of the channel:
entry, ttt, nsamps, status, and tbusy with its rise since the channel's
previous waveform (a channel that was busy shows it there).

Output (written into --out, default: <run name>_waveforms in the folder of the input file):
    waveforms_<tree>_<triggers>.png   10 triggers per page; further pages end in _p2, _p3, ...

Needs daqyy_common.py in the same folder.
"""

import os

# One thread is plenty here; keep numpy's BLAS from starting one per core.
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse
import math
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from daqyy_common import (
    BLUE, CRITICAL, GRID, open_root, channel_trees, read_branches, waveform, run_name,
    default_output, make_output_folder, enable_path_completion, ask, existing_file, output_folder,
    Progress,
)

ROWS_PER_PAGE = 10
ZOOM = 100                   # samples in the first/last zooms
DAQ_BASELINE_SAMPLES = 20    # raw2root_new: baseline and rms from the first 20 samples


def board_events(root_file, board_chans, chn, full_length=None):
    """One dict per trigger of the board, for channel chn: entry, nsamps,
    status (ok/short/long/missing) and tbusy. All channels of a board share
    the trigger time tag, so the board's triggers are the distinct ttt values
    of its channels. Returns (events, full length, tree of chn)."""
    trees = {c: root_file.Get(name) for c, name in board_chans.items()}
    bar = Progress("Scanning board", len(trees))
    index = {}
    for c, t in trees.items():
        index[c] = read_branches(t, ("ttt", "nsamps", "tbusy") if c == chn else ("ttt", "nsamps"))
        bar.update(note=board_chans[c])
    bar.close()
    triggers = np.unique(np.concatenate([ix["ttt"] for ix in index.values()]))
    if full_length is None:
        lengths, counts = np.unique(np.concatenate([ix["nsamps"] for ix in index.values()]),
                                    return_counts=True)
        full_length = int(lengths[np.argmax(counts)])
    mine = index[chn]
    entry_of = {}
    for i, t in enumerate(mine["ttt"]):
        entry_of.setdefault(int(t), i)
    events, last_tbusy = [], 0
    for k, t in enumerate(triggers):
        i = entry_of.get(int(t))
        if i is None:
            events.append({"trigger": k, "ttt": int(t), "status": "missing"})
            continue
        n, tbusy = int(mine["nsamps"][i]), int(mine["tbusy"][i])
        if n == full_length:
            status = "ok"
        else:
            status = "short" if n < full_length else "long"
        events.append({"trigger": k, "ttt": int(t), "entry": i, "nsamps": n, "status": status,
                       "tbusy": tbusy, "tbusy_step": tbusy - last_tbusy})
        last_tbusy = tbusy
    return events, full_length, trees[chn]


def print_events(name, board, chn, events, full_length):
    n_odd = sum(e["status"] != "ok" for e in events)
    print(f"\n{name} (board {board}, channel {chn}): {len(events)} triggers on the board, "
          f"{n_odd} without a full-length ({full_length}-sample) waveform")
    print(f"{'trigger':>7} {'entry':>5} {'ttt':>13} {'nsamps':>6}  {'status':<7} "
          f"{'tbusy':>11} {'tbusy rise':>11}")
    for e in events:
        if e["status"] == "missing":
            print(f"{e['trigger']:>7} {'-':>5} {e['ttt']:>13} {'-':>6}  {'MISSING':<7}")
            continue
        status = e["status"] if e["status"] == "ok" else e["status"].upper()
        print(f"{e['trigger']:>7} {e['entry']:>5} {e['ttt']:>13} {e['nsamps']:>6}  {status:<7} "
              f"{e['tbusy']:>11} {e['tbusy_step']:>11}")
    print()


def select_triggers(text, events):
    """Events to plot: "short", "all", or trigger numbers such as "0,2,5-7"."""
    text = text.strip().lower()
    if text == "all":
        return list(events)
    if text == "short":
        odd = [e for e in events if e["status"] != "ok"]
        reference = next((e for e in events if e["status"] == "ok"), None)
        return sorted(odd + ([reference] if reference else []), key=lambda e: e["trigger"])
    wanted = set()
    for part in text.split(","):
        first, _, last = part.strip().partition("-")
        try:
            wanted.update(range(int(first), int(last or first) + 1))
        except ValueError:
            raise ValueError(f"not a trigger list: {text!r} (use short, all, or e.g. 0,2,5-7)")
    rows = [e for e in events if e["trigger"] in wanted]
    if not rows:
        raise ValueError(f"no trigger {text!r}; this board has triggers 0-{len(events) - 1}")
    if len(rows) < len(wanted):
        print(f"  note: this board has triggers 0-{len(events) - 1}; plotting the "
              f"{len(rows)} of yours that exist")
    return rows


def one_of(options, what, allow_quit=False):
    """Converter for ask(): a board or channel number from options."""
    def convert(text):
        if allow_quit and text.lower() in ("", "q", "quit"):
            return None
        try:
            value = int(text)
        except ValueError:
            raise ValueError(f"not a {what} number: {text!r}")
        if value not in options:
            raise ValueError(f"no {what} {value}; choose from {', '.join(map(str, options))}")
        return value
    return convert


def draw_row(e, w, full_length, a_first, a_whole, a_last):
    """One trigger: zoom on the first samples, the whole waveform, zoom on the last."""
    k = e["trigger"]
    a_whole.set_xlim(-0.5, full_length - 0.5)
    if e["status"] == "missing":
        a_whole.axvspan(-0.5, full_length - 0.5, color=GRID, alpha=0.6, lw=0, zorder=0)
        a_whole.set_title(f"trigger {k}: no waveform recorded", loc="left")
        a_whole.set_title("⚠ missing", loc="right", color=CRITICAL)
        a_first.axis("off")
        a_last.axis("off")
        return
    n = len(w)
    x = np.arange(n)
    if n < full_length:
        a_whole.axvspan(n - 0.5, full_length - 0.5, color=GRID, alpha=0.6, lw=0, zorder=0)
    a_whole.plot(x, w, color=BLUE, lw=0.4)
    a_whole.set_xlim(-0.5, max(n, full_length) - 0.5)
    stats = f"mean {w.mean():.2f}, σ {w.std():.2f}" if n else "no samples"
    a_whole.set_title(f"trigger {k} (entry {e['entry']}): {n} of {full_length} samples, {stats}",
                      loc="left")
    if n < full_length:
        a_whole.set_title(f"⚠ short: samples {n}–{full_length - 1} not recorded (shaded)",
                          loc="right", color=CRITICAL)
    elif n > full_length:
        a_whole.set_title(f"⚠ long: {n - full_length} samples beyond the full length",
                          loc="right", color=CRITICAL)
    m = min(ZOOM, n)
    a_first.axvspan(-0.5, DAQ_BASELINE_SAMPLES - 0.5, color=BLUE, alpha=0.1, lw=0)
    a_first.plot(x[:m], w[:m], color=BLUE, lw=0.6, marker=".", ms=2.5)
    a_first.set_xlim(-1, ZOOM)
    a_last.plot(x[n - m:], w[n - m:], color=BLUE, lw=0.6, marker=".", ms=2.5)
    a_last.set_xlim(n - ZOOM - 1, n)
    for a in (a_first, a_last):
        a.xaxis.set_major_locator(MaxNLocator(3, integer=True))


def plot_channel(tree, name, board, chn, rows, full_length, run, path_base):
    """Pages of ROWS_PER_PAGE triggers, one row each. Returns the files written."""
    files = []
    pages = math.ceil(len(rows) / ROWS_PER_PAGE)
    bar = Progress("Plotting waveforms", len(rows) + pages)  # every row, and saving every page
    for page in range(pages):
        chunk = rows[page * ROWS_PER_PAGE:(page + 1) * ROWS_PER_PAGE]
        waves = {e["trigger"]: waveform(tree, e["entry"]) for e in chunk if e["status"] != "missing"}
        filled = [w for w in waves.values() if len(w)]
        lo = min(int(w.min()) for w in filled) - 2 if filled else 0
        hi = max(int(w.max()) for w in filled) + 2 if filled else 1
        height = 1.8 * len(chunk) + 1.0
        fig, axes = plt.subplots(len(chunk), 3, figsize=(14, height), squeeze=False, sharey=True,
                                 gridspec_kw={"width_ratios": [1, 5, 1], "wspace": 0.05, "hspace": 0.8})
        fig.subplots_adjust(top=1 - 0.8 / height, bottom=0.55 / height)  # margins in inches
        for e, (a_first, a_whole, a_last) in zip(chunk, axes):
            draw_row(e, waves.get(e["trigger"]), full_length, a_first, a_whole, a_last)
            bar.update(note=f"trigger {e['trigger']}")
        axes[0][0].set_ylim(lo, hi)  # shared by every panel, so pages compare row to row
        axes[0][0].set_title(f"first {ZOOM} samples (shaded: DAQ\nbaseline, samples "
                             f"0–{DAQ_BASELINE_SAMPLES - 1})", loc="left")
        axes[0][2].set_title(f"last {ZOOM} recorded\nsamples", loc="left")
        for a in axes[-1]:
            a.set_xlabel("sample")
        for row in axes:
            row[0].set_ylabel("ADC")
        page_text = f", page {page + 1} of {pages}" if pages > 1 else ""
        fig.suptitle(f"Raw waveforms of {name} (board {board}, channel {chn}){page_text}: {run}",
                     x=0.125, y=1 - 0.12 / height, va="top", ha="left", fontsize=10)
        path = f"{path_base}.png" if page == 0 else f"{path_base}_p{page + 1}.png"
        bar.update(0, f"saving page {page + 1} of {pages}")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        bar.update()
        files.append(path)
    bar.close()
    return files


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", nargs="?", default=None,
                    help="ROOT file written by raw2root_new (asked for if left out)")
    ap.add_argument("--board", type=int, default=None, help="Board number (asked for if left out)")
    ap.add_argument("--channel", type=int, default=None,
                    help="Channel of the board, 0-7 (asked for if left out)")
    ap.add_argument("--triggers", default=None,
                    help='"short" (default), "all", or trigger numbers such as 0,2,5-7')
    ap.add_argument("--out", default=None,
                    help="Output directory (default: <run name>_waveforms in the folder of the "
                         "input file)")
    ap.add_argument("--nsamps", type=int, default=None,
                    help="Full waveform length in samples (default: the most common length "
                         "on the board)")
    args = ap.parse_args()

    interactive = args.input is None
    if interactive:
        enable_path_completion()
        print("DAQ_YY raw waveforms. Enter keeps the value in [brackets]; Tab completes paths.")
        args.input = ask("ROOT file", convert=existing_file)
    root_file = open_root(args.input)
    if root_file is None:
        sys.exit(f"ERROR: cannot open {args.input}")
    layout = channel_trees(root_file)
    if not layout:
        sys.exit("ERROR: no ch<id> trees in the file.")
    run = run_name(args.input)
    out_dir = args.out or default_output(args.input, "waveforms")
    if interactive:
        out_dir = ask("Output folder", default=out_dir, convert=output_folder)
    make_output_folder(out_dir)
    print("Boards in the file (channels): "
          + ", ".join(f"{b} ({','.join(map(str, chs))})" for b, chs in layout.items()))

    board, chn = args.board, args.channel
    while True:
        if board is None:
            board = ask("Board", default=next(iter(layout)), convert=one_of(layout, "board"))
        elif board not in layout:
            sys.exit(f"ERROR: no board {board}; the file has boards {', '.join(map(str, layout))}")
        if chn is None:
            chn = ask("Channel", default=next(iter(layout[board])),
                      convert=one_of(layout[board], "channel"))
        elif chn not in layout[board]:
            sys.exit(f"ERROR: board {board} has no channel {chn}; it has "
                     f"{', '.join(map(str, layout[board]))}")
        name = layout[board][chn]
        events, full_length, tree = board_events(root_file, layout[board], chn, args.nsamps)
        print_events(name, board, chn, events, full_length)

        if interactive:
            text, rows = ask("Triggers to plot: short, all, or numbers such as 0,2,5-7",
                             default=args.triggers or "short",
                             convert=lambda t: (t, select_triggers(t, events)))
        else:
            text = args.triggers or "short"
            try:
                rows = select_triggers(text, events)
            except ValueError as err:
                sys.exit(f"ERROR: {err}")
        if text.strip().lower() == "short" and all(e["status"] == "ok" for e in rows):
            print("No short, long or missing waveforms on this channel; "
                  "plotting the first full-length one.")
        tag = text.strip().lower().replace(" ", "").replace(",", "_")
        if tag not in ("short", "all"):
            tag = "t" + tag
        for path in plot_channel(tree, name, board, chn, rows, full_length, run,
                                 os.path.join(out_dir, f"waveforms_{name}_{tag}")):
            print(f"Saved {path}")

        if not interactive:
            break
        board = ask("\nNext board (Enter to quit)", convert=one_of(layout, "board", allow_quit=True))
        if board is None:
            break
        chn = None


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("\nCancelled.")
