"""
daqyy_common.py

Shared pieces of the DAQ_YY scripts (analyze_baseline.py, plot_waveforms.py):
reading the ROOT files written by raw2root_new, asking for input on the
terminal, and the plot style. Keep this file next to the scripts.
"""

import glob
import os
import re
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import ROOT

ROOT.gROOT.SetBatch(True)

ADC_CODES = 1 << 14          # AD9694 samples are 14 bit: codes 0..16383
CHANNELS_PER_BOARD = 8

# ---- plot style: colors from the dataviz reference palette, light surface ----
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
BLUE, ORANGE = "#2a78d6", "#eb6834"                         # categorical slots 1, 2
GOOD, WARNING, CRITICAL = "#0ca30c", "#fab219", "#d03b3b"   # status colors
SEQUENTIAL = LinearSegmentedColormap.from_list(
    "blues", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]).with_extremes(bad=SURFACE)
DIVERGING = LinearSegmentedColormap.from_list(
    "blue_red", [BLUE, "#f0efec", "#e34948"]).with_extremes(bad=SURFACE)

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "xtick.color": AXIS, "ytick.color": AXIS, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "text.color": INK, "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "legend.frameon": False, "legend.fontsize": 7,
})


# ---- reading raw2root_new files ----

def run_name(path):
    """Run name of a ROOT file: its file name without the extension."""
    return os.path.splitext(os.path.basename(path))[0]


def default_output(path, kind):
    """Default output folder: <folder of the input file>/<run name>_<kind>,
    so every run gets its own folder next to its data."""
    return os.path.join(os.path.dirname(os.path.abspath(path)), f"{run_name(path)}_{kind}")


def open_root(path):
    """Open a ROOT file; None if it cannot be read."""
    try:
        root_file = ROOT.TFile.Open(path)
    except OSError:  # PyROOT raises instead of returning a null pointer
        return None
    return root_file if root_file and not root_file.IsZombie() else None


def channel_trees(root_file):
    """{board: {channel: tree name}} of the ch<id> trees, id = 10000 + 10*board + channel."""
    layout = {}
    for key in root_file.GetListOfKeys():
        m = re.fullmatch(r"ch(\d+)", key.GetName())
        if key.GetClassName() != "TTree" or not m:
            continue
        tid = int(m.group(1))
        board = (tid - 10000) // 10          # same decoding as raw2root_new
        layout.setdefault(board, {})[tid - 10000 - board * 10] = key.GetName()
    return {b: dict(sorted(chs.items())) for b, chs in sorted(layout.items())}


def read_branches(tree, names):
    """Integer branches (such as ttt, nsamps, tbusy) of every entry, as numpy
    arrays; only those branches are read, so this is fast."""
    tree.SetBranchStatus("*", 0)
    for name in names:
        tree.SetBranchStatus(name, 1)
    values = {name: [] for name in names}
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        for name in names:
            values[name].append(int(getattr(tree, name)))
    tree.SetBranchStatus("*", 1)
    if tree.GetBranch("x_numb"):
        tree.SetBranchStatus("x_numb", 0)  # only holds 0..nsamps-1; skip reading it
    return {name: np.array(v, dtype=np.int64) for name, v in values.items()}


def waveform(tree, entry):
    """ADC samples of one entry. Afterwards the tree holds that entry, so its
    other branches (tbusy, baseline, ...) can be read from it."""
    tree.GetEntry(entry)
    return np.frombuffer(tree.chdata, dtype=np.uint16, count=tree.nsamps).astype(np.int64)


# ---- progress bar ----

def clock(seconds):
    """Duration as m:ss, or h:mm:ss from one hour."""
    m, s = divmod(int(round(seconds)), 60)
    return f"{m // 60}:{m % 60:02d}:{s:02d}" if m >= 60 else f"{m}:{s:02d}"


class Progress:
    """Progress bar on stderr, so printed tables stay clean:
        Reading triggers  [██████░░░░░░░░░░░░]  12/40  30%  0:05 left  board 3
    On a terminal it is redrawn in place, at most 10 times a second. When
    stderr goes to a file or pipe only the finished line is written."""

    WIDTH = 24        # characters of the bar itself
    LABEL_WIDTH = 18  # labels are padded so bars of successive steps line up

    def __init__(self, label, total):
        self.label, self.total, self.done, self.note = label, max(int(total), 0), 0, ""
        self.start = self._drawn = time.monotonic()
        self.live = sys.stderr.isatty()
        self._length = 0
        sys.stdout.flush()  # keep earlier printed lines ahead of the bar in logs
        if self.live:
            self._draw("")

    def update(self, n=1, note=None):
        """Count n more steps done; note (optional) says what is happening now."""
        self.done += n
        changed = note is not None and note != self.note
        if note is not None:
            self.note = note
        now = time.monotonic()
        if self.live and (changed or now - self._drawn >= 0.1 or self.done >= self.total):
            self._drawn = now
            self._draw(self._eta(now))

    def close(self):
        """Finish the bar on its own line, with the time it took."""
        self.note = ""
        self._draw(f"done in {clock(time.monotonic() - self.start)}", end="\n")

    def _eta(self, now):
        if not self.done or self.done >= self.total:
            return ""
        return f"{clock((now - self.start) / self.done * (self.total - self.done))} left"

    def _draw(self, timing, end=""):
        frac = self.done / self.total if self.total else 1.0
        filled = int(round(frac * self.WIDTH))
        parts = [f"{self.label:<{self.LABEL_WIDTH}}", f"[{'█' * filled}{'░' * (self.WIDTH - filled)}]",
                 f"{self.done:>{len(str(self.total))}}/{self.total} {frac:>4.0%}", timing, self.note]
        line = "  ".join(p for p in parts if p)
        pad = " " * max(self._length - len(line), 0)   # wipe a longer previous line
        self._length = len(line)
        sys.stderr.write(("\r" if self.live else "") + line + pad + end)
        sys.stderr.flush()


# ---- asking on the terminal ----

def enable_path_completion():
    """Let Tab complete file and folder names at the prompts."""
    try:
        import readline
    except ImportError:
        return

    def complete(text, state):
        matches = sorted(glob.glob(os.path.expanduser(text) + "*"))
        matches = [m + "/" if os.path.isdir(m) else m for m in matches]
        return matches[state] if state < len(matches) else None

    readline.set_completer_delims(" \t\n;")
    readline.set_completer(complete)
    readline.parse_and_bind("tab: complete")


def ask(question, default=None, convert=str):
    """Ask until convert(answer) succeeds. Enter takes the default shown in
    [brackets]. convert raises ValueError with a message for a bad answer."""
    shown = f" [{default}]" if default is not None else ""
    while True:
        try:
            answer = input(f"{question}{shown}: ").strip()
        except EOFError:
            sys.exit("\nERROR: input ended before every question was answered.")
        if not answer and default is not None:
            answer = str(default)
        try:
            return convert(answer)
        except ValueError as err:
            print(f"  {err}")


def existing_file(text):
    path = os.path.abspath(os.path.expanduser(text))
    if not os.path.isfile(path):
        raise ValueError(f"no such file: {path}")
    return path


def output_folder(text):
    """Folder path; created right away, so an unusable path is asked again."""
    path = os.path.abspath(os.path.expanduser(text))
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as err:
        raise ValueError(f"cannot create {path}: {err.strerror}")
    return path


def make_output_folder(path):
    """Create the output folder, or stop with a hint to choose another one."""
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as err:
        sys.exit(f"ERROR: cannot create {path}: {err.strerror}; choose another folder with --out")
