#!/usr/bin/env python3
"""
run_analysis.py  -  the whole paddle muon analysis of one run

Runs the steps one after the other on one DAQ_YY run:
    1 extract_pulses.py   waveforms -> pulses.root
    2 build_events.py     pulses -> hits (both ends of a paddle) -> events.root
    3 calibrate.py        signal speed, offsets, resolution, attenuation, muon peak -> calibration.json
    4 reconstruct.py      tracks, angles, time of flight, toy expectation -> tracks.root
    5 make_plots.py       plots/ and analysis_summary.txt
    6 event_display.py    event_displays/ (the first 10 tracks and 5 with several paddles of one layer)

Usage:
    python3 run_analysis.py                    # asks for the ROOT file and the output folder
    python3 run_analysis.py run.root           # output: <run name>_muon next to run.root
    python3 run_analysis.py run.root --out /path/to/folder
    python3 run_analysis.py run.root --calibration other_muon/calibration.json   # use these constants, skip step 3
    python3 run_analysis.py run.root --from 3  # redo steps 3-6 (after editing the setup's [calibration])

Without a ROOT file on the command line the script asks for it and for the
output folder; Enter keeps the folder shown in [brackets] (a folder named after
the run next to the data file), Tab completes paths. Every option is also
explained by --help of the step it belongs to. Needs muon_common.py, the
setups/ folder and ../scripts/daqyy_common.py.
"""

import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse
import sys
import time

import muon_common as mc
import extract_pulses
import build_events
import calibrate
import reconstruct
import make_plots
import event_display


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file", nargs="?", default=None, help="ROOT file of the run (asked for if left out)")
    ap.add_argument("--out", default=None, help="Analysis folder (default: <run name>_muon next to the ROOT file)")
    ap.add_argument("--setup", default=None, help="Setup file (default: the one in setups/ valid for the run's date)")
    ap.add_argument("--calibration", default=None,
                    help="Use this calibration.json (e.g. from a longer run) instead of calibrating this run")
    ap.add_argument("--from", dest="start", type=int, default=1, choices=range(1, 7),
                    help="First step to run (default 1); the earlier steps' files must be in the folder")
    ap.add_argument("--events", type=int, default=10, help="Event pictures of the first N tracks (default 10, 0: none)")
    ap.add_argument("--max-waveforms", type=int, default=None,
                    help="Step 1 measures only the first N waveforms of every channel (quick tests)")
    ap.add_argument("--pulse-pre", type=int, default=None, help="pulse_pre if the run has no .conf")
    ap.add_argument("--polarity", choices=("positive", "negative"), default=None,
                    help="Pulse polarity of every channel (default: from the .conf)")
    ap.add_argument("--toy-muons", type=int, default=None, help="Toy muons per layer (default: setup [flux])")
    ap.add_argument("--lowpass", type=float, default=None,
                    help="Low-pass filter (MHz, -3 dB) before measuring; 0 = off (default: setup [pulse])")
    ap.add_argument("--lowpass-kind", choices=("gauss", "butter"), default=None, help="Filter shape")
    ap.add_argument("--raw-timing", action="store_true", help="Pulse time from the raw samples even with the filter")
    args = ap.parse_args()

    asked = args.file is None
    path = mc.ask_data_file(args.file)
    default = mc.default_out(path)
    out = mc.ask_out_folder(args.out, default) if asked else mc.ask_out_folder(args.out or default, None)
    start = time.monotonic()

    def step(n, text):
        print(f"\n==== step {n}: {text} " + "=" * max(0, 60 - len(text)))

    if args.start <= 1:
        step(1, "pulses (extract_pulses.py)")
        extract_pulses.run(path, out, args.setup, args.pulse_pre, args.polarity, args.max_waveforms,
                           args.lowpass, args.lowpass_kind, args.raw_timing)
    if args.start <= 2:
        step(2, "hits and events (build_events.py)")
        build_events.run(out, args.setup)
    if args.start <= 3:
        if args.calibration:
            step(3, f"calibration skipped: using {args.calibration}")
        else:
            step(3, "calibration (calibrate.py)")
            calibrate.run(out, args.setup)
    if args.start <= 4:
        step(4, "tracks (reconstruct.py)")
        reconstruct.run(out, args.setup, args.calibration, args.toy_muons)
    if args.start <= 5:
        step(5, "plots and summary (make_plots.py)")
        make_plots.run(out, args.setup)
    if args.start <= 6 and args.events > 0:
        step(6, "event pictures (event_display.py)")
        event_display.run(out, first=args.events, setup_path=args.setup, several=5)
    print(f"\nDone in {mc.clock(time.monotonic() - start)}. Results in {out}")
    print(f"Start with {os.path.join(out, 'analysis_summary.txt')}")


if __name__ == "__main__":
    main()
