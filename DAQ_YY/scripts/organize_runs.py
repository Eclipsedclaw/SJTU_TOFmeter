#!/usr/bin/env python3
"""
organize_runs.py  -  move the runs of a flat raw data folder into a folder per run

Before: every run in one folder (DAQ_yy_data/raw):
    run00073_file0000_20261126_160206.dat  .conf  .root  _events.root
    run00073_file0001_20261126_170256.dat  .conf  .root
    run00073_file0000_20261126_160206_muon/  _muon.log   (outputs: _muon, _baseline, _waveforms, ...)
After: the layout the DAQ web app writes from v1.0.28 on (its setting "a folder for every run"):
    DAQ_yy_data/run00073/dat/    the .dat and .conf files
    DAQ_yy_data/run00073/root/   the .root and _events.root files
    DAQ_yy_data/run00073/muon/   the paddle analysis of the run, muon.log next to it
                                 (an output made from another file alone: muon_file0001/)
A run is file 0 of a run number and the files after it, up to the next file 0. A number used for two runs
(run00068 in November 2026) gives each run folder its start as well: run00068_20261119_175732/.

Nothing is copied: every file is renamed on the same disk, never over an existing one. A run stays where it
is while it may still change: a .dat without its .root (being converted), a conversion's .part file, or a
file changed in the last 10 minutes (a run or an analysis going). Every move is listed in
organize_runs_YYYYMMDD_HHMMSS.txt in the target folder; --undo reads that list and moves everything back.

ORDER: move the runs only after the DAQ app writes run folders (v1.0.28, setting on, raw data folder = the
target folder, e.g. DAQ_yy_data). An older app still gives testDAQ the flat folder: with no .dat left there,
testDAQ would number the next run 0.

Usage:
    python3 organize_runs.py /data/data/jiancheng/DAQ_yy_data/raw           # dry run: prints every move
    python3 organize_runs.py /data/data/jiancheng/DAQ_yy_data/raw --do      # moves (into the folder above)
    python3 organize_runs.py /data/data/jiancheng/DAQ_yy_data/raw --to /data/data/jiancheng/DAQ_yy_data --do
    python3 organize_runs.py --undo /data/data/jiancheng/DAQ_yy_data/organize_runs_20261126_190000.txt
"""

import argparse
import os
import re
import sys
import time

STEM = re.compile(r"^(\.?)(run(\d+)_file(\d+)_(\d{8}_\d{6}))(.*)$")
QUIET_S = 600                   # a run changed within this time is left where it is
CONVERT_S = 86400               # a .dat without .root younger than this may still be converted: left too
OUTPUT = re.compile(r"_([A-Za-z]\w*)")


def kind_of(rest, is_dir):
    """Where an entry of a run goes, from the text after its stem: (kind, name), or None (not known: left)."""
    if is_dir:
        m = OUTPUT.fullmatch(rest)
        return ("output", m.group(1)) if m else None
    if rest in (".dat", ".conf"):
        return ("dat", "")
    if rest in (".root", "_events.root"):
        return ("root", "")
    m = re.fullmatch(r"_([A-Za-z]\w*)\.log", rest)
    return ("log", m.group(1)) if m else None


def newest_change(path):
    """The newest modification time of a file, or of a folder and everything in it."""
    t = os.lstat(path).st_mtime
    if os.path.isdir(path) and not os.path.islink(path):
        for here, dirs, files in os.walk(path):
            for name in dirs + files:
                try:
                    t = max(t, os.lstat(os.path.join(here, name)).st_mtime)
                except OSError:
                    pass
    return t


def open_paths():
    """Paths that running processes (those this account may look at) have open or work in: an analysis
    reading a run's ROOT files for many minutes writes nothing meanwhile, so the time alone cannot tell."""
    out = set()
    for pid in filter(str.isdigit, os.listdir("/proc")):
        for link in [f"/proc/{pid}/cwd"] + [f"/proc/{pid}/fd/{fd}" for fd in _list(f"/proc/{pid}/fd")]:
            try:
                out.add(os.readlink(link))
            except OSError:
                pass
    return out


def _list(folder):
    try:
        return os.listdir(folder)
    except OSError:
        return []


def plan(folder, target):
    """The runs of the flat folder: [{"name", "dest", "moves": [(src, dst)], "why_not": [...]}], and the run files
    whose kind is not known (left)."""
    entries = {}
    unknown = []
    for name in sorted(os.listdir(folder)):
        m = STEM.match(name)
        if not m:
            continue
        hidden, stem, number, fileno, start, rest = m.groups()
        path = os.path.join(folder, name)
        if hidden:                                  # .NAME.root.xxxx.part: a conversion going
            entries.setdefault(stem, []).append((name, ("busy", "")))
            continue
        kind = kind_of(rest, os.path.isdir(path))
        if kind is None:
            unknown.append(name)
            continue
        entries.setdefault(stem, []).append((name, kind))
    key = {s: (int(STEM.match(s).group(3)), STEM.match(s).group(5), int(STEM.match(s).group(4))) for s in entries}
    runs = []
    for stem in sorted(entries, key=lambda s: key[s]):
        number, start, fileno = key[stem]
        if not runs or runs[-1]["number"] != number or fileno == 0:
            runs.append({"number": number, "start": start, "stems": []})
        runs[-1]["stems"].append(stem)
    per_number = {}
    for r in runs:
        per_number[r["number"]] = per_number.get(r["number"], 0) + 1
    now = time.time()
    in_use = open_paths()
    out = []
    for r in runs:
        name = f"run{r['number']:05d}" + (f"_{r['start']}" if per_number[r["number"]] > 1 else "")
        dest = os.path.join(target, name)
        first = r["stems"][0]
        moves, why = [], []
        for stem in r["stems"]:
            fileno = key[stem][2]
            names = {n for n, _ in entries[stem]}
            if (stem + ".dat" in names and stem + ".root" not in names
                    and now - os.lstat(os.path.join(folder, stem + ".dat")).st_mtime < CONVERT_S):
                why.append(f"{stem}.dat has no .root yet")       # (an older one was never converted: it moves)
            for n, (kind, sub) in entries[stem]:
                src = os.path.join(folder, n)
                if kind == "busy":
                    why.append(f"{n}: a conversion is going")
                    continue
                if now - newest_change(src) < QUIET_S:
                    why.append(f"{n} changed in the last {QUIET_S // 60} minutes")
                if any(p == src or p.startswith(src + os.sep) for p in in_use):
                    why.append(f"{n} is open in a running program")
                if kind in ("dat", "root"):
                    dst = os.path.join(dest, kind, n)
                else:
                    own = sub if stem == first else f"{sub}_file{fileno:04d}"
                    dst = os.path.join(dest, own + (".log" if kind == "log" else ""))
                if os.path.lexists(dst):
                    why.append(f"{dst} exists already")
                moves.append((src, dst))
        out.append({"name": name, "dest": dest, "moves": moves, "why_not": why})
    return out, unknown


def undo(listing):
    lines = [ln.rstrip("\n").split("\t") for ln in open(listing, encoding="utf-8") if "\t" in ln]
    back, made = 0, set()
    for src, dst in reversed(lines):
        if os.path.lexists(dst) and not os.path.lexists(src):
            os.rename(dst, src)
            back += 1
            made.update({os.path.dirname(dst), os.path.dirname(os.path.dirname(dst))})
    for d in sorted(made, key=len, reverse=True):       # the run folders, if empty now
        try:
            os.rmdir(d)
        except OSError:
            pass
    print(f"moved back {back} of {len(lines)} entries listed in {listing}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", help="The flat raw data folder (e.g. DAQ_yy_data/raw)")
    ap.add_argument("--to", default=None, help="Where the run folders go (default: the folder above)")
    ap.add_argument("--do", action="store_true", help="Move (without it: only print what would move)")
    ap.add_argument("--undo", default=None, help="Move back everything listed in this organize_runs_*.txt")
    args = ap.parse_args()
    if args.undo:
        undo(args.undo)
        return
    if not args.folder:
        ap.error("give the flat raw data folder (or --undo LIST)")
    folder = os.path.realpath(os.path.expanduser(args.folder))
    target = os.path.realpath(os.path.expanduser(args.to)) if args.to else os.path.dirname(folder)
    if not os.path.isdir(folder) or not os.path.isdir(target):
        sys.exit(f"ERROR: {folder} or {target} is not a folder")
    if os.stat(folder).st_dev != os.stat(target).st_dev:
        sys.exit(f"ERROR: {target} is on another disk than {folder}: this script only renames, it never copies")
    runs, unknown = plan(folder, target)
    ready = [r for r in runs if not r["why_not"]]
    print(f"{folder}: {len(runs)} runs, {sum(len(r['moves']) for r in runs)} entries; into {target}")
    for r in runs:
        kinds = {}
        for src, dst in r["moves"]:
            part = os.path.relpath(dst, r["dest"]).split(os.sep)[0]
            kinds[part] = kinds.get(part, 0) + 1
        what = ", ".join(f"{k} {n}" for k, n in sorted(kinds.items()))
        state = "moves" if not r["why_not"] else "STAYS: " + "; ".join(sorted(set(r["why_not"]))[:3])
        print(f"  {r['name']:<26} {what:<44} {state}")
    for n in unknown:
        print(f"  (left where it is, kind not known: {n})")
    print(f"{len(ready)} runs can move, {len(runs) - len(ready)} stay.")
    if not args.do:
        print("Dry run: nothing moved. Add --do to move (after the DAQ app writes run folders, see --help).")
        return
    listing = os.path.join(target, f"organize_runs_{time.strftime('%Y%m%d_%H%M%S')}.txt")
    moved = 0
    with open(listing, "w", encoding="utf-8") as log:
        log.write(f"# organize_runs.py {time.strftime('%Y-%m-%d %H:%M:%S')}: from {folder} to {target}; "
                  "every line: old<TAB>new (--undo moves them back)\n")
        for r in ready:
            for src, dst in r["moves"]:
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                if os.path.lexists(dst):            # appeared since the plan: never written over
                    print(f"  not moved, {dst} exists: {src}")
                    continue
                os.rename(src, dst)
                log.write(f"{src}\t{dst}\n")
                log.flush()
                moved += 1
    print(f"Moved {moved} entries of {len(ready)} runs; the list (for --undo): {listing}")


if __name__ == "__main__":
    main()
