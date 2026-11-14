"""
muon_common.py

Code shared by the steps of the paddle muon analysis (DAQ_YY/analysis/):

  * the setup file (setups/*.toml): which channel reads which end of which
    paddle, and where the paddles are;
  * the run's .conf, which XIAOLAB DAQ writes next to every data file;
  * the analysis files the steps hand on to each other (pulses.root,
    events.root, calibration.json, tracks.root) and their 'info' record;
  * fit helpers, the straight-line track fit, and a toy cosmic-muon
    simulation that predicts rates and angles for the setup's geometry;
  * prompts, progress bars and the plot style, from ../scripts/daqyy_common.py.

Keep this file next to the step scripts. Units: times in ns, lengths in cm,
angles in degrees, pulse heights in ADC counts, charges in ADC x ns.
"""

import os

# One thread is plenty; keep numpy's BLAS from starting one per core.
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import datetime
import glob
import json
import math
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.normpath(os.path.join(HERE, os.pardir, "scripts"))
SETUPS = os.path.join(HERE, "setups")
sys.path.insert(0, SCRIPTS)

import daqyy_common as dc  # noqa: E402  (puts ROOT in batch mode and sets the plot style)
from daqyy_common import (  # noqa: E402,F401  (re-exported for the step scripts)
    ROOT, Progress, ask, existing_file, output_folder, enable_path_completion,
    open_root, channel_trees, run_name, make_output_folder, clock,
    SURFACE, INK, INK2, MUTED, GRID, AXIS, SEQUENTIAL,
)

try:
    import tomllib
except ImportError:  # Python < 3.11
    try:
        import tomli as tomllib
    except ImportError:
        tomllib = None

VERSION = "0.1.0"

# ---- units and constants ----
NS_PER_SAMPLE = 2.0        # AD9694: 500 MS/s
NS_PER_TICK = 2.0          # ttt ticks
NS_PER_BUSY_TICK = 4.0     # tbusy ticks (see the AD9694 notes in ../README.md)
C_CM_PER_NS = 29.9792458   # speed of light
ADC_MAX = 16383            # 14-bit samples

# ---- files the steps hand on ----
PULSES_FILE = "pulses.root"
EVENTS_FILE = "events.root"
CALIB_FILE = "calibration.json"
TRACKS_FILE = "tracks.root"
PLOTS_DIR = "plots"
DISPLAY_DIR = "event_displays"
KIND = "muon"              # default output folder: <run name>_muon next to the data file

# ---- pulse flags (bits of 'flags' in pulses.root) ----
FLAG_SHORT = 1             # record too short for the baseline or the peak search
FLAG_NO_EDGE = 2           # no rising edge found
FLAG_NOISE = 4             # noise burst, not a pulse (flat-top or undershoot test)
FLAG_LOW = 8               # height below min_height
FLAG_SATURATED = 16        # a sample at the ADC limit: the height is a lower limit
FLAG_PILEUP = 32           # a second pulse rises on the tail
FLAG_SLOPE = 64            # baseline window not flat (tail of an earlier pulse)
FLAG_NOISY = 128           # baseline noise above max_baseline_rms
FLAG_NO_FIT = 256          # timing = "fit", but no edge fit: the time is the 20 % crossing
REJECT_MASK = FLAG_SHORT | FLAG_NO_EDGE | FLAG_NOISE | FLAG_LOW
PULSE_FLAGS = [
    (FLAG_SHORT, "SHORT", "record too short for the baseline window or the peak search (rejected)"),
    (FLAG_NO_EDGE, "NO_EDGE", "no rising edge through the constant fraction of the height (rejected)"),
    (FLAG_NOISE, "NOISE", "noise burst: no flat top after the peak, or a deep swing below the baseline (rejected)"),
    (FLAG_LOW, "LOW", "height below min_height (rejected)"),
    (FLAG_SATURATED, "SATURATED", "a sample at the ADC limits or beyond max_raw_adc / min_raw_adc (the amplifier "
                                  "clips): the height is only a lower limit (kept, left out of the muon-peak fits)"),
    (FLAG_PILEUP, "PILEUP", "a second pulse rises on the tail of the first (kept)"),
    (FLAG_SLOPE, "SLOPE", "the baseline window is not flat, e.g. the tail of an earlier pulse (kept)"),
    (FLAG_NOISY, "NOISY", "baseline noise above max_baseline_rms (kept)"),
    (FLAG_NO_FIT, "NO_FIT", "no rising-edge fit (window beyond the record, fit more than fit_max_shift from the "
                            "crossing, or no template): the time is the 20 % crossing (kept)"),
]

# ---- hit kinds and track flags ----
HIT_AB, HIT_A_ONLY, HIT_B_ONLY = 0, 1, 2
TRACK_ADJACENT = 1         # a layer had hits in neighbouring paddles, merged into one point
TRACK_AMBIGUOUS = 2        # a layer had hits in paddles that are not neighbours: the biggest was used
TRACK_OUTSIDE = 4          # a hit's position along its paddle lies beyond the paddle's end (> 3 sigma)
TRACK_FLAGS = [
    (TRACK_ADJACENT, "ADJACENT", "hits in neighbouring paddles of one layer, merged into one point"),
    (TRACK_AMBIGUOUS, "AMBIGUOUS", "hits in paddles of one layer that are not neighbours; the biggest one was used"),
    (TRACK_OUTSIDE, "OUTSIDE", "a position along a paddle lies more than 3 sigma beyond the paddle's end"),
]

# ---- defaults of the setup file's settings sections ----
PULSE_DEFAULTS = {
    "polarity": "conf",          # "conf" (pulsePolarity of the run's .conf), "positive" or "negative"
    "baseline_end": -20,         # the baseline is the mean of samples 0 .. trigger sample + baseline_end
    "search": [-30, 50],         # the peak is searched in trigger sample + search[0] .. + search[1]
    "smooth": 5,                 # samples averaged for the height (running mean)
    "cfd_fraction": 0.2,         # the pulse time is where the rising edge crosses this fraction of the height
    "timing": "fit",             # "fit": the rising edge fitted with the channel's own template, the time where the
                                 # fit reaches cfd_fraction of its height; "cfd": the crossing between two samples
    "fit_window": [-4, 4],       # samples fitted, around the one at or before the crossing: the edge, not the top
    "fit_max_shift": 1.5,        # samples; a fit further than this from the crossing has failed (NO_FIT)
    "template_pulses": 5000,     # clean pulses (no flag) per channel averaged into its edge template
    "template_min_height": 300.0,  # ADC; smaller pulses are left out of the template
    "template_min_pulses": 50,   # fewer clean pulses: no template, the channel keeps the crossing (NO_FIT)
    "charge_window": [-5, 130],  # samples around the edge summed for the charge
    "min_height": 100.0,         # ADC; lower peaks are not pulses
    "plateau_window": [5, 45],   # samples after the peak where a real pulse stays high (flat-top test)
    "min_plateau": 0.5,          # flat-top test: mean there / height; 0 switches the test off
    "max_undershoot": 0.5,       # noise test: deepest swing below the baseline / height
    "pileup_fraction": 0.25,     # a later rise of 4 samples above this fraction of the height ...
    "pileup_min": 200.0,         # ... and above this many ADC counts is a second pulse
    "max_baseline_slope": 100.0, # ADC between the first and last 30 baseline samples
    "max_baseline_rms": 60.0,    # ADC
    "max_raw_adc": 16383,        # a raw sample at or above this: SATURATED (the amplifier clips lower)
    "min_raw_adc": 0,            # a raw sample at or below this: SATURATED
    "lowpass_mhz": 0.0,          # low-pass filter (-3 dB frequency, MHz) before measuring; 0 = off
    "lowpass_kind": "gauss",     # "gauss" (no overshoot) or "butter" (4th order); both zero-phase
    "lowpass_timing": True,      # the pulse time from the filtered waveform too (false: raw samples)
}
EVENT_DEFAULTS = {
    "pair_margin_ns": 5.0,       # A and B pulses pair up if |tA - tB| <= length / v_eff + this
    "window_ns": 30.0,           # hits within this time of an event's first hit belong to the event
}
CALIBRATION_DEFAULTS = {
    "v_eff_mode": "per_paddle",  # "per_paddle": each paddle its own signal speed; "shared": one for all
    "v_eff_cm_per_ns": 15.0,     # signal speed along a paddle until calibrate.py measures it
    "speed_min_mip": 0.8,        # speed and attenuation fits use hits with sqrt(hA hB) > this x muon peak
    "tdiff_sigma_ns": 0.6,       # resolution of tA - tB until calibrate.py measures it
    "min_hits": 50,              # fewer A&B hits in a paddle: keep the defaults
    "min_tof_events": 20,        # events per paddle needed for per-paddle time offsets
    "mip_mev": 0.0,              # most probable energy loss of a vertical muon in a paddle, MeV; 0 = from its
                                 # thickness (plastic scintillator, Landau formula: 0.95 MeV for 0.6 cm)
}
DISPLAY_DEFAULTS = {
    "side": "",                  # a side of the room to mark in the event pictures: "+x", "-x", "+y", "-y"
    "side_label": "",            # its name, e.g. "wall"
}
TRACK_DEFAULTS = {
    "points": "measured",        # a hit's position along its paddle for the track fit: "measured" = where tA - tB
                                 # puts it (error: the resolution; near an end it can lie beyond the paddle, as in
                                 # the hit maps), "inside" = the best estimate inside the paddle (the measurement's
                                 # Gaussian cut to the paddle, error: its spread); tracks.root keeps both
}
FLUX_DEFAULTS = {
    "vertical_intensity": 70.0,  # muons / (m^2 s sr) at sea level above ~1 GeV (PDG review)
    "cos_power": 2.0,            # I(theta) = I0 cos^n(theta)
    "toy_muons": 1000000,        # toy muons generated per layer
    "seed": 1,
}

# ---- plot colours: the dataviz reference palette, categorical slots in their fixed order ----
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
END_COLORS = {"A": CATEGORICAL[0], "B": CATEGORICAL[1]}
FIT_COLOR = CATEGORICAL[1]           # fit curves over a blue (slot 1) histogram
EXPECTED_COLOR = CATEGORICAL[0]      # toy expectation (line) under the data (ink points)


def layer_color(index):
    return CATEGORICAL[index % len(CATEGORICAL)]


class SetupError(ValueError):
    pass


# =====================================================================
# setup file
# =====================================================================

def _date(value, where):
    if value is None or isinstance(value, datetime.date):
        return value
    try:
        return datetime.date.fromisoformat(str(value))
    except ValueError:
        raise SetupError(f"{where}: '{value}' is not a date (write 2026-11-12)")


def parse_channel(value, where):
    """Channel id from 10010 or "ch10010" (DAQ_YY: 10000 + 10 x board + channel)."""
    text = str(value).strip().lower()
    if text.startswith("ch"):
        text = text[2:]
    if not text.isdigit() or int(text) < 10010:
        raise SetupError(f"{where}: channel '{value}' is not a DAQ_YY channel id such as 10010 "
                         f"(10000 + 10 x board + channel)")
    return int(text)


def channel_board(chid):
    return (chid - 10000) // 10


def channel_number(chid):
    return (chid - 10000) % 10


def _merge(defaults, given, where):
    out = dict(defaults)
    for key, value in (given or {}).items():
        if key not in defaults:
            print(f"WARNING: {where}: unknown setting '{key}' ignored (known: {', '.join(defaults)})")
            continue
        out[key] = value
    return out


class Layer:
    def __init__(self, name, z):
        self.name, self.z = name, z
        self.index = -1
        self.paddles = []


class Paddle:
    """One scintillator bar. Its position along its length is u (cm from the
    centre, positive towards the end that channel A reads)."""

    def __init__(self, index, spec, layers, where):
        def need(key):
            if key not in spec:
                raise SetupError(f"{where}: missing '{key}'")
            return spec[key]

        self.index = index
        self.name = str(need("name"))
        where = f"{where} ({self.name})"
        self.label = str(spec.get("label", self.name))
        layer = str(need("layer"))
        if layer not in layers:
            raise SetupError(f"{where}: layer '{layer}' is not one of the [[layers]] ({', '.join(layers)})")
        self.layer = layers[layer]
        self.axis = str(need("axis")).strip().lower()
        if self.axis not in ("x", "y"):
            raise SetupError(f"{where}: axis must be \"x\" or \"y\" (the paddle's long direction)")
        self.x = float(need("x_cm"))
        self.y = float(need("y_cm"))
        self.z = float(spec.get("z_cm", self.layer.z))
        self.length = float(need("length_cm"))
        self.width = float(need("width_cm"))
        self.thickness = float(spec.get("thickness_cm", 1.0))
        if min(self.length, self.width, self.thickness) <= 0:
            raise SetupError(f"{where}: sizes must be positive")
        a_end = str(need("A_end")).strip().lower()
        if a_end not in (f"+{self.axis}", f"-{self.axis}"):
            raise SetupError(f"{where}: A_end must be \"+{self.axis}\" or \"-{self.axis}\" for a paddle "
                             f"along {self.axis} (the end that channel A reads)")
        self.a_end = a_end
        self.a_sign = 1.0 if a_end[0] == "+" else -1.0
        channels = need("channels")
        self.channels = {}
        for end in ("A", "B"):
            if end in channels:
                self.channels[end] = parse_channel(channels[end], where)
        unknown = set(channels) - {"A", "B"}
        if unknown or not self.channels:
            raise SetupError(f"{where}: channels must be {{ A = <id>, B = <id> }}")

    @property
    def size_x(self):
        return self.length if self.axis == "x" else self.width

    @property
    def size_y(self):
        return self.width if self.axis == "x" else self.length

    @property
    def both_ends(self):
        return "A" in self.channels and "B" in self.channels

    def point(self, u):
        """(x, y, z) arrays of the points u cm from the centre along the paddle."""
        u = np.asarray(u, dtype=np.float64)
        if self.axis == "x":
            return self.x + self.a_sign * u, np.full_like(u, self.y), np.full_like(u, self.z)
        return np.full_like(u, self.x), self.y + self.a_sign * u, np.full_like(u, self.z)

    def u_of(self, x, y):
        """Position along the paddle of the point (x, y)."""
        return self.a_sign * ((np.asarray(x) - self.x) if self.axis == "x" else (np.asarray(y) - self.y))

    def contains(self, x, y):
        return (np.abs(np.asarray(x) - self.x) <= 0.5 * self.size_x) & \
               (np.abs(np.asarray(y) - self.y) <= 0.5 * self.size_y)

    def neighbours(self, other, tolerance=0.5):
        """Side by side in the same layer, along the same axis, touching or with a small gap."""
        if other.layer is not self.layer or other.axis != self.axis or other is self:
            return False
        if self.axis == "x":
            gap = abs(other.y - self.y) - 0.5 * (self.width + other.width)
            overlap = abs(other.x - self.x) < 0.5 * (self.length + other.length)
        else:
            gap = abs(other.x - self.x) - 0.5 * (self.width + other.width)
            overlap = abs(other.y - self.y) < 0.5 * (self.length + other.length)
        return overlap and gap <= tolerance

    def describe(self):
        chans = ", ".join(f"ch{c} ({e})" for e, c in self.channels.items())
        return (f"{self.name}: {self.layer.name} layer, along {self.axis}, centre ({self.x:g}, {self.y:g}, "
                f"{self.z:g}) cm, {self.length:g} x {self.width:g} x {self.thickness:g} cm, A end at "
                f"{self.a_end}; {chans}")


class Setup:
    """The setup file: layers, paddles, channel map and analysis settings."""

    def __init__(self, raw, path):
        self.path = os.path.abspath(path)
        where = os.path.basename(path)
        self.name = str(raw.get("name") or os.path.splitext(where)[0])
        self.description = str(raw.get("description", ""))
        self.geometry_source = str(raw.get("geometry_source", ""))
        self.valid_from = _date(raw.get("valid_from"), where)
        self.valid_until = _date(raw.get("valid_until"), where)
        layers = {}
        for i, spec in enumerate(raw.get("layers") or []):
            if "name" not in spec or "z_cm" not in spec:
                raise SetupError(f"{where}: [[layers]] number {i + 1} needs name and z_cm")
            name = str(spec["name"])
            if name in layers:
                raise SetupError(f"{where}: layer '{name}' is listed twice")
            layers[name] = Layer(name, float(spec["z_cm"]))
        if not layers:
            raise SetupError(f"{where}: no [[layers]]")
        self.layers = sorted(layers.values(), key=lambda layer: -layer.z)   # top first
        for i, layer in enumerate(self.layers):
            layer.index = i
        self.paddles = []
        for i, spec in enumerate(raw.get("paddles") or []):
            p = Paddle(i, spec, layers, f"{where}: [[paddles]] number {i + 1}")
            if any(q.name == p.name for q in self.paddles):
                raise SetupError(f"{where}: paddle name '{p.name}' is used twice")
            self.paddles.append(p)
            p.layer.paddles.append(p)
        if not self.paddles:
            raise SetupError(f"{where}: no [[paddles]]")
        self.channels = {}
        for p in self.paddles:
            for end, chid in p.channels.items():
                if chid in self.channels:
                    q, e = self.channels[chid]
                    raise SetupError(f"{where}: ch{chid} is both {q.name}-{e} and {p.name}-{end}")
                self.channels[chid] = (p, end)
        self.pulse = _merge(PULSE_DEFAULTS, raw.get("pulse"), f"{where} [pulse]")
        self.events = _merge(EVENT_DEFAULTS, raw.get("events"), f"{where} [events]")
        self.calibration = _merge(CALIBRATION_DEFAULTS, raw.get("calibration"), f"{where} [calibration]")
        self.flux = _merge(FLUX_DEFAULTS, raw.get("flux"), f"{where} [flux]")
        self.display = _merge(DISPLAY_DEFAULTS, raw.get("display"), f"{where} [display]")
        self.tracks = _merge(TRACK_DEFAULTS, raw.get("tracks"), f"{where} [tracks]")
        if self.tracks["points"] not in ("measured", "inside"):
            raise SetupError(f"{where} [tracks]: points must be \"measured\" or \"inside\", "
                             f"not {self.tracks['points']!r}")
        self.board_offsets = {int(k): float(v) for k, v in (raw.get("board_offsets_ns") or {}).items()}
        self.channel_offsets = {parse_channel(k, where): float(v)
                                for k, v in (raw.get("channel_offsets_ns") or {}).items()}
        for key in ("search", "charge_window", "plateau_window"):
            v = self.pulse[key]
            if not (isinstance(v, (list, tuple)) and len(v) == 2 and int(v[0]) < int(v[1])):
                raise SetupError(f"{where} [pulse]: {key} must be two numbers, low first")
        if not 0.02 <= float(self.pulse["cfd_fraction"]) <= 0.95:
            raise SetupError(f"{where} [pulse]: cfd_fraction must lie in 0.02 .. 0.95")
        if self.pulse["timing"] not in ("fit", "cfd"):
            raise SetupError(f"{where} [pulse]: timing must be \"fit\" or \"cfd\", not {self.pulse['timing']!r}")
        v = self.pulse["fit_window"]
        if not (isinstance(v, (list, tuple)) and len(v) == 2 and int(v[0]) <= 0 < int(v[1])):
            raise SetupError(f"{where} [pulse]: fit_window must be two numbers of samples, the first <= 0 < the second")
        if not 0 < float(self.pulse["fit_max_shift"]) <= 5:
            raise SetupError(f"{where} [pulse]: fit_max_shift must lie in 0 .. 5 samples")

    def paddle(self, name):
        for p in self.paddles:
            if p.name == name:
                return p
        raise KeyError(name)

    def channel_label(self, chid):
        if chid in self.channels:
            p, end = self.channels[chid]
            return f"{p.name}-{end}"
        return "not in setup"

    def time_offset(self, chid):
        """ns subtracted from every pulse time of the channel (board + channel offsets)."""
        return self.board_offsets.get(channel_board(chid), 0.0) + self.channel_offsets.get(chid, 0.0)

    def describe(self):
        lines = [f"setup {self.name} ({self.path})"]
        if self.description:
            lines.append(f"  {self.description}")
        if self.geometry_source:
            lines.append(f"  geometry: {self.geometry_source}")
        for layer in self.layers:
            lines.append(f"  layer {layer.name}: z = {layer.z:g} cm, {len(layer.paddles)} paddles")
        for p in self.paddles:
            lines.append(f"    {p.describe()}")
        return "\n".join(lines)


def load_setup(path):
    if tomllib is None:
        sys.exit("ERROR: reading the setup file needs Python 3.11 or newer (or: pip install tomli)")
    try:
        with open(path, "rb") as f:
            raw = tomllib.load(f)
    except OSError as err:
        sys.exit(f"ERROR: cannot read setup file {path}: {err.strerror}")
    except tomllib.TOMLDecodeError as err:
        sys.exit(f"ERROR: {path} is not valid TOML: {err}")
    try:
        return Setup(raw, path)
    except SetupError as err:
        sys.exit(f"ERROR: {err}")


def run_datetime(path):
    """Start of the run from its file name (..._20261112_221007), or None."""
    m = re.search(r"_(\d{8})_(\d{6})", os.path.basename(path))
    if not m:
        return None
    try:
        return datetime.datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H%M%S")
    except ValueError:
        return None


def choose_setup(explicit, data_path):
    """The given setup file, else the file in setups/ that is valid for the run's date
    (the latest valid_from on or before it)."""
    if explicit:
        return load_setup(os.path.abspath(os.path.expanduser(explicit)))
    when = run_datetime(data_path)
    files = sorted(glob.glob(os.path.join(SETUPS, "*.toml")))
    if not files:
        sys.exit(f"ERROR: no setup files in {SETUPS}; give one with --setup")
    valid = []
    for path in files:
        s = load_setup(path)
        if when is not None:
            day = when.date()
            if (s.valid_from and day < s.valid_from) or (s.valid_until and day > s.valid_until):
                continue
        valid.append(s)
    if not valid:
        sys.exit(f"ERROR: no setup file in {SETUPS} is valid for a run of {when:%Y-%m-%d}; "
                 f"give one with --setup")
    valid.sort(key=lambda s: (s.valid_from or datetime.date.min, s.path))
    return valid[-1]


# =====================================================================
# the run's .conf
# =====================================================================

class RunConf:
    """Settings of a run from <run>.conf, which XIAOLAB DAQ writes next to the data."""

    def __init__(self, path):
        self.path = path
        with open(path, encoding="utf-8", errors="replace") as f:
            self.text = f.read()
        body = "\n".join(re.sub(r"#.*", "", line) for line in self.text.splitlines())

        def number(name, text=body):
            m = re.search(rf"\b{name}\s*=\s*(-?0[xX][0-9a-fA-F]+|-?\d+)", text)
            return int(m.group(1), 0) if m else None

        self.pulse_pre_typed = number("pulse_pre")
        self.pulse_post = number("pulse_post")
        self.trigmode = number("trigmode")
        self.ctrlmode = number("ctrlmode")
        self.daqruntime = number("daqruntime")
        self.trig_pre = number("trig_pre")
        self.trig_post = number("trig_post")
        # the board library uses 20 for pulse_pre outside 5..1020 (as analyze_sync.py)
        p = self.pulse_pre_typed
        self.pulse_pre = p if p is not None and 5 <= p <= 1020 else (20 if p is not None else None)
        self.boards = {}
        m = re.search(r"board_config\s*=\s*\((.*)\)", body, re.S)
        for block in re.findall(r"\{(.*?)\}", m.group(1), re.S) if m else []:
            brd = number("brdID", block)
            if brd is None:
                continue
            thresholds = {int(c): int(v, 0) for c, v in
                          re.findall(r"\bch(\d+)_threshold\s*=\s*(-?0[xX][0-9a-fA-F]+|-?\d+)", block)}
            self.boards[brd] = {
                "enable": number("channelenable", block),
                "polarity": number("pulsePolarity", block),
                "thresholds": thresholds,
            }
        m = re.search(r"^# ---- note ----\n(.*?)^# ---- end of note ----", self.text, re.M | re.S)
        self.note = " ".join(line.lstrip("#").strip() for line in m.group(1).splitlines()
                             if line.strip("# ")) if m else ""
        m = re.search(r"^#\s*started by (\S+) at (.+)$", self.text, re.M)
        self.started_by, self.started = (m.group(1), m.group(2).strip()) if m else ("", "")

    def threshold(self, chid):
        b = self.boards.get(channel_board(chid))
        return None if b is None else b["thresholds"].get(channel_number(chid))

    def positive(self, chid):
        """True / False from pulsePolarity (bit = 1: positive), None if unknown."""
        b = self.boards.get(channel_board(chid))
        if b is None or b["polarity"] is None:
            return None
        return bool((b["polarity"] >> channel_number(chid)) & 1)

    def summary(self):
        return {"file": self.path, "pulse_pre": self.pulse_pre, "pulse_pre_typed": self.pulse_pre_typed,
                "pulse_post": self.pulse_post, "trigmode": self.trigmode, "ctrlmode": self.ctrlmode,
                "daqruntime": self.daqruntime, "trig_pre": self.trig_pre, "trig_post": self.trig_post,
                "note": self.note, "started": self.started, "started_by": self.started_by,
                "boards": {str(k): {"enable": v["enable"], "polarity": v["polarity"],
                                    "thresholds": {str(c): t for c, t in v["thresholds"].items()}}
                           for k, v in self.boards.items()}}


def read_conf(data_path):
    """RunConf of <run>.conf next to the data file, or None."""
    path = os.path.splitext(data_path)[0] + ".conf"
    return RunConf(path) if os.path.isfile(path) else None


# =====================================================================
# analysis files
# =====================================================================

_LEAF = {np.dtype("float64"): "D", np.dtype("float32"): "F", np.dtype("int64"): "L",
         np.dtype("int32"): "I", np.dtype("uint64"): "l", np.dtype("uint32"): "i"}


def write_tree(path, name, columns, update=False):
    """Write {column: numpy array} as TTree `name` (RECREATE, or UPDATE to add a tree)."""
    columns = {k: np.ascontiguousarray(v) for k, v in columns.items()}
    for k, v in columns.items():
        if v.dtype not in _LEAF:
            raise TypeError(f"column {k}: dtype {v.dtype} cannot be written (use float64/32, int64/32)")
    n = len(next(iter(columns.values())))
    if n:
        opts = ROOT.RDF.RSnapshotOptions()
        opts.fMode = "UPDATE" if update else "RECREATE"
        ROOT.RDF.FromNumpy(columns).Snapshot(name, path, list(columns), opts)
        return
    # RDF cannot take empty arrays: write the empty tree by hand
    f = ROOT.TFile(path, "UPDATE" if update else "RECREATE")
    tree = ROOT.TTree(name, name)
    buffers = []
    for k, v in columns.items():
        buf = np.zeros(1, dtype=v.dtype)
        buffers.append(buf)
        tree.Branch(k, buf, f"{k}/{_LEAF[v.dtype]}")
    tree.Write()
    f.Close()


def read_tree(path, name, columns=None):
    """{column: numpy array} of TTree `name` (all columns, or the listed ones)."""
    df = ROOT.RDataFrame(name, path)
    data = df.AsNumpy(columns) if columns else df.AsNumpy()
    return {k: np.asarray(v) for k, v in data.items()}


def write_info(path, info):
    """Store a JSON record as the TNamed 'info' in a ROOT file (UPDATE)."""
    f = ROOT.TFile(path, "UPDATE")
    ROOT.TNamed("info", json.dumps(info, default=_json_default)).Write("info", ROOT.TObject.kOverwrite)
    f.Close()


def read_info(path):
    f = open_root(path)
    if f is None:
        sys.exit(f"ERROR: cannot open {path}")
    obj = f.Get("info")
    info = json.loads(obj.GetTitle()) if obj else {}
    f.Close()
    return info


def _json_default(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (datetime.date, datetime.datetime)):
        return value.isoformat()
    raise TypeError(f"not JSON serialisable: {type(value)}")


def write_json(path, data):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=_json_default)
        f.write("\n")


def read_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def need_file(folder, name, made_by):
    """Path of an analysis file in `folder`; stop with a hint if it is missing."""
    path = os.path.join(folder, name)
    if not os.path.isfile(path):
        sys.exit(f"ERROR: no {name} in {folder}; run {made_by} first")
    return path


def now_text():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def stamp(step):
    return {"step": step, "version": VERSION, "created": now_text()}


def default_out(data_path):
    """Default output folder: <run name>_muon next to the data file (its real folder)."""
    return os.path.join(os.path.dirname(os.path.realpath(data_path)), f"{run_name(data_path)}_{KIND}")


def short_run(name):
    """run00056 for run00056_file0000_20261112_221007."""
    m = re.match(r"run\d+", name)
    return m.group(0) if m else name


def ask_data_file(given):
    """The raw ROOT file from the command line, or asked on the terminal."""
    if given:
        path = os.path.realpath(os.path.expanduser(given))
        if not os.path.isfile(path):
            sys.exit(f"ERROR: no such file: {path}")
        return path
    enable_path_completion()
    return os.path.realpath(ask("ROOT file (raw data of the run)", convert=existing_file))


def ask_out_folder(given, default):
    """Output folder from --out, or asked (Enter keeps the default next to the data)."""
    if given:
        path = os.path.abspath(os.path.expanduser(given))
        make_output_folder(path)
        return path
    enable_path_completion()
    return ask("Output folder", default=default, convert=output_folder)


def ask_analysis_folder(given, needed, made_by):
    """Folder with the outputs of the earlier steps, from the command line or asked."""
    def check(text):
        path = os.path.abspath(os.path.expanduser(text))
        if not os.path.isfile(os.path.join(path, needed)):
            raise ValueError(f"no {needed} in {path} (made by {made_by})")
        return path
    if given:
        try:
            return check(given)
        except ValueError as err:
            sys.exit(f"ERROR: {err}")
    enable_path_completion()
    return ask(f"Analysis folder (with {needed})", convert=check)


def setup_of(folder, explicit=None):
    """The setup the earlier steps used (recorded in pulses.root), or the given one."""
    if explicit:
        return load_setup(os.path.abspath(os.path.expanduser(explicit)))
    info = read_info(need_file(folder, PULSES_FILE, "extract_pulses.py"))
    path = info.get("setup", {}).get("path")
    if not path or not os.path.isfile(path):
        sys.exit(f"ERROR: the setup file recorded in {PULSES_FILE} ({path}) is gone; give one with --setup")
    return load_setup(path)


def flag_text(flags, table=PULSE_FLAGS):
    names = [name for bit, name, _ in table if flags & bit]
    return "+".join(names) if names else "ok"


# =====================================================================
# fits
# =====================================================================

def numerical_errors(nll, p, scale):
    """Errors and covariance from the numerical Hessian of a negative log-likelihood."""
    p = np.asarray(p, dtype=np.float64)
    n = len(p)
    h = np.asarray(scale, dtype=np.float64) * 1e-3 + 1e-9
    hess = np.zeros((n, n))
    f0 = nll(p)
    for i in range(n):
        for j in range(i, n):
            if i == j:
                pp, pm = p.copy(), p.copy()
                pp[i] += h[i]
                pm[i] -= h[i]
                hess[i, i] = (nll(pp) - 2 * f0 + nll(pm)) / h[i] ** 2
            else:
                ppp, ppm, pmp, pmm = p.copy(), p.copy(), p.copy(), p.copy()
                ppp[i] += h[i]; ppp[j] += h[j]
                ppm[i] += h[i]; ppm[j] -= h[j]
                pmp[i] -= h[i]; pmp[j] += h[j]
                pmm[i] -= h[i]; pmm[j] -= h[j]
                hess[i, j] = hess[j, i] = (nll(ppp) - nll(ppm) - nll(pmp) + nll(pmm)) / (4 * h[i] * h[j])
    try:
        cov = np.linalg.inv(hess)
        err = np.sqrt(np.where(np.diag(cov) > 0, np.diag(cov), np.nan))
    except np.linalg.LinAlgError:
        cov = np.full((n, n), np.nan)
        err = np.full(n, np.nan)
    return err, cov


def _phi(z):
    return np.exp(-0.5 * z * z) / math.sqrt(2 * math.pi)


def _Phi(z):
    from scipy.special import ndtr
    return ndtr(z)


def box_pdf(t, lo, hi, sigma, bkg, r0, r1):
    """Uniform distribution lo..hi smeared by a Gaussian of sigma, plus a flat fraction
    bkg, normalised on r0..r1 (the distribution of tA - tB along a paddle)."""
    def G(x, m):
        z = (x - m) / sigma
        return (x - m) * _Phi(z) + sigma * _phi(z)
    box = (_Phi((t - lo) / sigma) - _Phi((t - hi) / sigma)) / (hi - lo)
    norm = (G(r1, lo) - G(r0, lo) - G(r1, hi) + G(r0, hi)) / (hi - lo)
    return (1 - bkg) * box / max(norm, 1e-12) + bkg / (r1 - r0)


def box_fit(t, r0, r1, half_width=None):
    """Unbinned fit of box_pdf to t in r0..r1. Free: the box's ends lo, hi, the
    resolution sigma and the background fraction; with half_width given the box
    is centre +- half_width and only centre, sigma and background are free.
    Returns dict (lo, hi, sigma, bkg, err, cov, n, nll) or None."""
    from scipy.optimize import minimize
    t = np.asarray(t, dtype=np.float64)
    t = t[(t >= r0) & (t <= r1)]
    if len(t) < 10:
        return None

    def nll4(p):
        lo, hi, sigma, bkg = p
        if not (r0 <= lo < hi <= r1) or not 0.02 <= sigma <= 0.5 * (r1 - r0) or not 0 <= bkg < 0.9:
            return 1e30
        return -np.sum(np.log(np.maximum(box_pdf(t, lo, hi, sigma, bkg, r0, r1), 1e-300)))

    if half_width is None:
        nll = nll4
        p0 = [np.percentile(t, 4), np.percentile(t, 96), 0.4, 0.02]
        starts = (p0, [p0[0] + 0.3, p0[1] - 0.3, 0.6, 0.05], [p0[0] - 0.3, p0[1] + 0.3, 0.25, 0.01])
        scale = [1, 1, 0.3, 0.05]
    else:
        def nll(p):
            return nll4([p[0] - half_width, p[0] + half_width, p[1], p[2]])
        m = float(np.median(t))
        starts = ([m, 0.4, 0.02], [m, 0.8, 0.05], [m, 0.2, 0.01])
        scale = [0.3, 0.3, 0.05]
    best = None
    for start in starts:
        res = minimize(nll, start, method="Nelder-Mead",
                       options={"xatol": 1e-5, "fatol": 1e-7, "maxiter": 20000, "maxfev": 20000})
        if best is None or res.fun < best.fun:
            best = res
    if best.fun >= 1e29:
        return None
    err, cov = numerical_errors(nll, best.x, scale)
    if half_width is None:
        lo, hi, sigma, bkg = best.x
    else:
        centre, sigma, bkg = best.x
        lo, hi = centre - half_width, centre + half_width
        # in the lo, hi, sigma, bkg order of the free fit (lo and hi move together)
        c = cov
        cov = np.array([[c[0, 0], c[0, 0], c[0, 1], c[0, 2]], [c[0, 0], c[0, 0], c[0, 1], c[0, 2]],
                        [c[1, 0], c[1, 0], c[1, 1], c[1, 2]], [c[2, 0], c[2, 0], c[2, 1], c[2, 2]]])
        err = np.array([err[0], err[0], err[1], err[2]])
    return {"lo": lo, "hi": hi, "sigma": sigma, "bkg": bkg, "err": err, "cov": cov, "n": len(t),
            "nll": best.fun, "range": (r0, r1)}


def shared_speed_fit(samples, v_range=(6.0, 30.0)):
    """One signal speed for several paddles: samples = [(tdiff array, length, r0, r1)].
    Each paddle's box is offset +- length / v; for every v the paddles' offsets,
    resolutions and backgrounds are fitted, and v minimises the summed -log L
    (a profile likelihood). Returns dict (v, err, nll, per-paddle fits) or None."""
    from scipy.optimize import minimize_scalar

    def total(v):
        s = 0.0
        for t, length, r0, r1 in samples:
            f = box_fit(t, r0, r1, half_width=length / v)
            if f is None:
                return 1e30
            s += f["nll"]
        return s

    vs = np.linspace(v_range[0], v_range[1], 25)
    scan = np.array([total(v) for v in vs])
    k = int(np.argmin(scan))
    lo, hi = vs[max(k - 1, 0)], vs[min(k + 1, len(vs) - 1)]
    res = minimize_scalar(total, bounds=(lo, hi), method="bounded", options={"xatol": 1e-3})
    v = float(res.x)
    f0 = total(v)

    def crossing(direction):
        from scipy.optimize import brentq
        target = f0 + 0.5
        step, x = 0.05 * v, v
        for _ in range(60):
            x2 = x + direction * step
            if not v_range[0] * 0.5 <= x2 <= v_range[1] * 1.5:
                return float("nan")
            if total(x2) >= target:
                return brentq(lambda y: total(y) - target, min(x, x2), max(x, x2), xtol=1e-3)
            x = x2
        return float("nan")

    up, down = crossing(+1), crossing(-1)
    err = 0.5 * (up - down) if np.isfinite(up) and np.isfinite(down) else float("nan")
    fits = [box_fit(t, r0, r1, half_width=length / v) for t, length, r0, r1 in samples]
    return {"v": v, "err": err, "err_up": up - v, "err_down": v - down, "nll": f0, "fits": fits,
            "scan_v": vs, "scan_nll": scan}


def moyal_pdf(x, mpv, width):
    z = (np.asarray(x, dtype=np.float64) - mpv) / width
    return np.exp(-0.5 * (z + np.exp(-z))) / (width * math.sqrt(2 * math.pi))


def moyal_cdf(x, mpv, width):
    from scipy.special import erfc
    z = (np.asarray(x, dtype=np.float64) - mpv) / width
    return erfc(np.exp(-0.5 * z) / math.sqrt(2))


def mode_estimate(x):
    """Most frequent value: peak of a smoothed histogram in log(x)."""
    x = np.asarray(x, dtype=np.float64)
    x = x[x > 0]
    if len(x) < 5:
        return float(np.median(x)) if len(x) else float("nan")
    lx = np.log(x)
    bins = max(10, min(60, int(len(x) / 6)))
    counts, edges = np.histogram(lx, bins=bins, range=(np.percentile(lx, 1), np.percentile(lx, 99)))
    smooth = np.convolve(counts, [1, 2, 3, 2, 1], mode="same")
    k = int(np.argmax(smooth))
    return float(np.exp(0.5 * (edges[k] + edges[k + 1])))


def moyal_fit(x, lo=None, hi=None):
    """Unbinned fit of a Moyal distribution (a close approximation of the Landau
    shape of energy loss) to x in lo..hi, normalised on that range. MPV = the peak."""
    from scipy.optimize import minimize
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x) & (x > 0)]
    if len(x) < 20:
        return None
    mode = mode_estimate(x)
    lo = 0.75 * mode if lo is None else lo
    hi = np.percentile(x, 98) if hi is None else hi
    sel = x[(x >= lo) & (x <= hi)]
    if len(sel) < 15 or hi <= lo:
        return None

    def nll(p):
        mpv, width = p
        if width <= 0.005 * mode or not lo - 0.5 * mode < mpv < hi:
            return 1e30
        norm = moyal_cdf(hi, mpv, width) - moyal_cdf(lo, mpv, width)
        if norm <= 1e-12:
            return 1e30
        return -np.sum(np.log(np.maximum(moyal_pdf(sel, mpv, width) / norm, 1e-300)))

    res = minimize(nll, [mode, 0.12 * mode], method="Nelder-Mead",
                   options={"xatol": 1e-4, "fatol": 1e-7, "maxiter": 10000})
    if res.fun >= 1e29:
        return None
    err, _ = numerical_errors(nll, res.x, [0.05 * mode, 0.02 * mode])
    return {"mpv": res.x[0], "width": res.x[1], "mpv_err": err[0], "width_err": err[1],
            "n": len(sel), "n_all": len(x), "range": (lo, hi), "mode": mode}


def expected_position(u_measured, half_length, sigma):
    """Mean true position along a paddle given the measured one, for muons spread
    evenly over -half_length..+half_length and a Gaussian resolution sigma: the mean
    of N(u_measured, sigma) cut to the paddle. Fitting against it instead of the
    measured u keeps a straight-line slope unbiased ("regression calibration")."""
    from scipy.special import ndtr
    um = np.asarray(u_measured, dtype=np.float64)
    if not sigma or sigma <= 0:
        return np.clip(um, -half_length, half_length)
    a = (-half_length - um) / sigma
    b = (half_length - um) / sigma
    return um + sigma * (_phi(a) - _phi(b)) / np.maximum(ndtr(b) - ndtr(a), 1e-12)


def position_in_paddle(u_measured, half_length, sigma):
    """Best estimate of where along a paddle a muon crossed, knowing it crossed inside.

    The measured position (from tA - tB) scatters around the true one with the
    resolution sigma, so near an end it often lies beyond the paddle (about
    2 sigma / (sqrt(2 pi) L) of all hits for muons spread evenly). Combining the
    measurement with "the muon was inside -L/2..L/2" gives a Gaussian cut to the
    paddle: its mean is the estimate (always inside) and its spread the error
    (smaller near the ends). Returns (mean, sigma) arrays."""
    from scipy.stats import truncnorm
    um = np.asarray(u_measured, dtype=np.float64)
    if not sigma or sigma <= 0:
        return np.clip(um, -half_length, half_length), np.zeros_like(um)
    a = (-half_length - um) / sigma
    b = (half_length - um) / sigma
    mean, var = truncnorm.stats(a, b, loc=um, scale=sigma, moments="mv")
    mean = np.clip(np.nan_to_num(np.asarray(mean, dtype=np.float64), nan=0.0), -half_length, half_length)
    std = np.sqrt(np.clip(np.nan_to_num(np.asarray(var, dtype=np.float64), nan=0.0), 1e-6, None))
    return mean, std


def track_point(u_measured, paddle, sigma, setup):
    """Position along a paddle that a track goes through, and its error, as the setup's [tracks] points
    says: "measured" = the measured position itself with the resolution sigma (default), "inside" = the
    best estimate inside the paddle (position_in_paddle). Returns (u, error) arrays."""
    um = np.asarray(u_measured, dtype=np.float64)
    if setup.tracks["points"] == "inside":
        return position_in_paddle(um, 0.5 * paddle.length, sigma)
    return um, np.full_like(um, float(sigma or 0.0))


def robust_line(x, y, rounds=5, clip=3.0):
    """Straight line y = a + b x by least squares, dropping points beyond clip x rms
    in each round. Returns dict or None."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 5:
        return None
    for _ in range(rounds):
        xs, ys = x[keep], y[keep]
        A = np.vstack([np.ones_like(xs), xs]).T
        coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
        resid = y - (coef[0] + coef[1] * x)
        rms = np.sqrt(np.mean(resid[keep] ** 2))
        new = np.isfinite(x) & np.isfinite(y) & (np.abs(resid) <= clip * max(rms, 1e-12))
        if new.sum() < 5 or np.array_equal(new, keep):
            break
        keep = new
    xs, ys = x[keep], y[keep]
    n = len(xs)
    resid = ys - (coef[0] + coef[1] * xs)
    s2 = np.sum(resid ** 2) / max(n - 2, 1)
    sxx = np.sum((xs - xs.mean()) ** 2)
    b_err = math.sqrt(s2 / sxx) if sxx > 0 else float("nan")
    a_err = math.sqrt(s2 * (1 / n + xs.mean() ** 2 / sxx)) if sxx > 0 else float("nan")
    return {"a": coef[0], "b": coef[1], "a_err": a_err, "b_err": b_err, "rms": math.sqrt(s2),
            "n": n, "n_all": int(np.isfinite(x).sum()), "keep": keep}


# =====================================================================
# tracks
# =====================================================================

def fit_lines(x, y, z, sx, sy):
    """Straight lines x = ax + bx (z - z0), y = ay + by (z - z0), one per row, through
    the points of the row weighted by 1/sigma^2. NaN entries are skipped. Every row
    needs points at two different z at least. Returns a dict of arrays."""
    x, y, z = (np.asarray(a, dtype=np.float64) for a in (x, y, z))
    use = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    wz = np.where(use, 1.0, 0.0)
    z0 = np.sum(np.where(use, z, 0), axis=1) / np.maximum(wz.sum(axis=1), 1)
    dz = np.where(use, z - z0[:, None], 0.0)
    out = {"z0": z0, "npoints": use.sum(axis=1)}
    chi2 = np.zeros(len(z0))
    for name, v, s in (("x", x, sx), ("y", y, sy)):
        w = np.where(use, 1.0 / np.asarray(s, dtype=np.float64) ** 2, 0.0)
        v0 = np.where(use, v, 0.0)
        S, Sz, Szz = w.sum(1), (w * dz).sum(1), (w * dz * dz).sum(1)
        Sv, Szv = (w * v0).sum(1), (w * dz * v0).sum(1)
        D = S * Szz - Sz * Sz
        with np.errstate(invalid="ignore", divide="ignore"):
            b = (S * Szv - Sz * Sv) / D
            a = (Szz * Sv - Sz * Szv) / D
            var_b = S / D
        chi2 += np.sum(w * (v0 - a[:, None] - b[:, None] * dz) ** 2 * use, axis=1)
        out[f"a{name}"], out[f"b{name}"], out[f"var_b{name}"] = a, b, var_b
    out["chi2"] = chi2
    out["ndf"] = 2 * out["npoints"] - 4
    return out


def direction(bx, by, var_bx=None, var_by=None):
    """Zenith and azimuth (degrees) of the direction the muon came from, for lines
    x = ax + bx z, y = ay + by z (z up). Azimuth 0 = from +x, 90 = from +y.
    With variances: the zenith's error too."""
    bx, by = np.asarray(bx, dtype=np.float64), np.asarray(by, dtype=np.float64)
    r = np.hypot(bx, by)
    theta = np.degrees(np.arctan(r))
    phi = np.degrees(np.arctan2(by, bx)) % 360.0
    if var_bx is None:
        return theta, phi
    with np.errstate(invalid="ignore", divide="ignore"):
        s2 = np.where(r > 1e-9, (bx * bx * var_bx + by * by * var_by) / np.maximum(r * r, 1e-18),
                      0.5 * (var_bx + var_by))
        sigma_theta = np.degrees(np.sqrt(s2) / (1 + r * r))
    return theta, phi, sigma_theta


# plastic scintillator (polyvinyltoluene based), PDG atomic and nuclear properties
PVT = {"density": 1.032, "z_over_a": 0.54141, "i_ev": 64.7,
       "x0": 0.1464, "x1": 2.4855, "c_bar": 3.1997, "a": 0.16101, "k": 3.2393}


def _landau_terms(thickness_cm, beta_gamma, m):
    bg2 = beta_gamma ** 2
    b2 = bg2 / (1 + bg2)
    x = np.asarray(thickness_cm, dtype=np.float64) * m["density"]   # g / cm^2
    xi = 0.5 * 0.307075 * m["z_over_a"] * x / b2                     # MeV
    lx = math.log10(beta_gamma)
    if lx >= m["x1"]:
        delta = 2 * math.log(10) * lx - m["c_bar"]
    elif lx >= m["x0"]:
        delta = 2 * math.log(10) * lx - m["c_bar"] + m["a"] * (m["x1"] - lx) ** m["k"]
    else:
        delta = 0.0
    return xi, b2, bg2, delta


def landau_mpv_mev(thickness_cm, beta_gamma=20.0, m=PVT):
    """Most probable energy loss (MeV) of a minimum-ionising particle crossing a plastic
    slab of this thickness (cm; a number or an array) straight through (PDG review,
    'Passage of particles through matter': Delta_p = xi [ln(2 m c^2 b^2 g^2 / I) +
    ln(xi / I) + j - b^2 - delta], xi = K/2 Z/A x / b^2, j = 0.2). It is lower than the
    mean loss (about 2 MeV/cm), which the long Landau tail pulls up. Cosmic muons (beta
    gamma 10-100) give 0.945-0.964 MeV for 0.6 cm (0.953 at the beta gamma 20 used).
    Used only by the simulation (known values), never on the data."""
    xi, b2, bg2, delta = _landau_terms(thickness_cm, beta_gamma, m)
    i_mev = m["i_ev"] * 1e-6
    out = xi * (math.log(2 * 0.51099895 * bg2 / i_mev) + np.log(xi / i_mev) + 0.2 - b2 - delta)
    return float(out) if np.ndim(out) == 0 else out


def landau_xi_mev(thickness_cm, beta_gamma=20.0, m=PVT):
    """Width parameter xi (MeV) of the Landau energy-loss distribution (FWHM about 4 xi)."""
    xi = _landau_terms(thickness_cm, beta_gamma, m)[0]
    return float(xi) if np.ndim(xi) == 0 else xi


def sample_moyal(rng, loc, scale):
    """Random numbers from a Moyal distribution (close to Landau) with peak loc and width scale."""
    from scipy.special import erfcinv
    loc = np.asarray(loc, dtype=np.float64)
    u = rng.uniform(1e-9, 1 - 1e-9, loc.shape)
    return loc + np.asarray(scale) * (-2.0 * np.log(math.sqrt(2.0) * erfcinv(u)))


def mip_mev(paddle, setup):
    """A vertical muon's most probable energy loss in the paddle (MeV), for the simulation and
    for text on plots only: the data processing never converts ADC to MeV with it."""
    given = float(setup.calibration.get("mip_mev") or 0.0)
    return given if given > 0 else landau_mpv_mev(paddle.thickness)


def position_sigma(calib, paddle, setup):
    """Resolution (cm) of the position along a paddle: (v_eff / 2) x sigma(tA - tB)."""
    pc = (calib or {}).get("paddles", {}).get(paddle.name, {})
    v = pc.get("v_eff_cm_per_ns", setup.calibration["v_eff_cm_per_ns"])
    s = pc.get("tdiff_sigma_ns", setup.calibration["tdiff_sigma_ns"])
    return 0.5 * v * s


# =====================================================================
# toy cosmic-muon simulation
# =====================================================================

def toy_cosmics(setup, calib=None, n_per_layer=None, seed=None):
    """Toy cosmic muons through the setup's paddles, to predict rates and angles.

    For each layer as the highest layer crossed, muons are started on a rectangle
    around that layer's paddles, with I(theta) = I0 cos^n(theta) (through a
    horizontal plane: dN ~ cos^(n+1)(theta) dcos(theta) dphi). Each straight line is
    followed to every paddle's centre plane; a paddle is crossed if the point lies
    inside its outline. Muons that cross a higher layer are left to that layer's
    sample, so none is counted twice. Detection efficiency, the trigger threshold,
    the soft (electron, gamma) component, multiple scattering and the paddles'
    thickness are not simulated.

    Each muon crossing two layers or more is reconstructed as the data are: the
    position along each crossed paddle smeared by its resolution (and then taken
    inside the paddle if the setup's [tracks] points = "inside"), the position
    across it set to the paddle's centre, then the same line fit.

    Returns a dict of arrays (one row per such muon) with 'weight' in Hz per muon.
    """
    flux = setup.flux
    n = int(n_per_layer or flux["toy_muons"])
    rng = np.random.default_rng(flux["seed"] if seed is None else seed)
    power = float(flux["cos_power"])
    i0 = float(flux["vertical_intensity"]) * 1e-4          # per cm^2 s sr
    nl = len(setup.layers)
    rows = []
    for layer in setup.layers:
        if not layer.paddles:
            continue
        x0 = min(p.x - 0.5 * p.size_x for p in layer.paddles)
        x1 = max(p.x + 0.5 * p.size_x for p in layer.paddles)
        y0 = min(p.y - 0.5 * p.size_y for p in layer.paddles)
        y1 = max(p.y + 0.5 * p.size_y for p in layer.paddles)
        area = (x1 - x0) * (y1 - y0)
        weight = i0 * area * 2 * math.pi / (power + 2) / n
        xs = rng.uniform(x0, x1, n)
        ys = rng.uniform(y0, y1, n)
        cos_t = rng.random(n) ** (1.0 / (power + 2))
        tan_t = np.sqrt(1 - cos_t ** 2) / cos_t
        ph = rng.uniform(0, 2 * math.pi, n)
        bx, by = tan_t * np.cos(ph), tan_t * np.sin(ph)
        crossed = np.full((n, nl), -1, dtype=np.int64)
        multi = np.zeros(n, dtype=bool)
        for other in setup.layers:
            for p in other.paddles:
                px = xs + bx * (p.z - layer.z)
                py = ys + by * (p.z - layer.z)
                inside = p.contains(px, py)
                multi |= inside & (crossed[:, other.index] >= 0)
                crossed[:, other.index] = np.where(inside & (crossed[:, other.index] < 0),
                                                   p.index, crossed[:, other.index])
        keep = crossed[:, layer.index] >= 0
        if layer.index > 0:
            keep &= (crossed[:, :layer.index] < 0).all(axis=1)
        keep &= (crossed >= 0).sum(axis=1) >= 2
        idx = np.where(keep)[0]
        if not len(idx):
            continue
        rows.append({"weight": np.full(len(idx), weight), "paddles": crossed[idx], "multi": multi[idx],
                     "x0": xs[idx], "y0": ys[idx], "z0": np.full(len(idx), layer.z),
                     "bx_true": bx[idx], "by_true": by[idx]})
    if not rows:
        return None
    toy = {k: np.concatenate([r[k] for r in rows]) for k in rows[0]}
    m = len(toy["weight"])
    xm = np.full((m, nl), np.nan)
    ym = np.full((m, nl), np.nan)
    zm = np.full((m, nl), np.nan)
    sxm = np.full((m, nl), np.nan)
    sym = np.full((m, nl), np.nan)
    utm = np.full((m, nl), np.nan)                                     # true position along the crossed paddle
    for p in setup.paddles:
        rowsel = toy["paddles"][:, p.layer.index] == p.index
        if not rowsel.any():
            continue
        dz = p.z - toy["z0"][rowsel]
        xt = toy["x0"][rowsel] + toy["bx_true"][rowsel] * dz
        yt = toy["y0"][rowsel] + toy["by_true"][rowsel] * dz
        sig = position_sigma(calib, p, setup)
        j = p.layer.index
        utm[rowsel, j] = p.u_of(xt, yt)
        u_meas = utm[rowsel, j] + rng.normal(0, sig, rowsel.sum())
        u, s_u = track_point(u_meas, p, sig, setup)                    # as the data's reconstruction
        px, py, pz = p.point(u)
        xm[rowsel, j], ym[rowsel, j], zm[rowsel, j] = px, py, pz
        along_x = p.axis == "x"
        sxm[rowsel, j] = s_u if along_x else p.width / math.sqrt(12)
        sym[rowsel, j] = p.width / math.sqrt(12) if along_x else s_u
    fit = fit_lines(xm, ym, zm, np.where(np.isfinite(sxm), sxm, 1.0), np.where(np.isfinite(sym), sym, 1.0))
    theta, phi = direction(fit["bx"], fit["by"])
    theta_true, phi_true = direction(toy["bx_true"], toy["by_true"])
    zmax = np.nanmax(zm, axis=1)
    zmin = np.nanmin(zm, axis=1)
    s_meas = (zmax - zmin) * np.sqrt(1 + fit["bx"] ** 2 + fit["by"] ** 2)
    s_true = (zmax - zmin) * np.sqrt(1 + toy["bx_true"] ** 2 + toy["by_true"] ** 2)
    toy.update({"x": xm, "y": ym, "z": zm, "u_true": utm, "theta": theta, "phi": phi, "theta_true": theta_true,
                "phi_true": phi_true, "s": s_meas, "s_true": s_true, "tof_true": s_true / C_CM_PER_NS})
    return toy


def singles_rate(paddle, setup):
    """Expected rate (Hz) of muons crossing a horizontal paddle, any direction."""
    power = float(setup.flux["cos_power"])
    i0 = float(setup.flux["vertical_intensity"]) * 1e-4
    return i0 * paddle.length * paddle.width * 2 * math.pi / (power + 2)


# =====================================================================
# plots
# =====================================================================

def save_figure(fig, path):
    import matplotlib.pyplot as plt
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def footer(fig, text):
    """Small grey line (run, setup) under everything the axes draw, x labels included."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = [ax.get_tightbbox(renderer) for ax in fig.axes]
    boxes = [b for b in boxes if b is not None]
    bottom = min(b.y0 for b in boxes) if boxes else 0.0
    left = min(b.x0 for b in boxes) if boxes else 0.0
    gap = 4 * fig.dpi / 72
    fig.text(left / fig.bbox.width, (bottom - gap) / fig.bbox.height, text, fontsize=6.5, color=MUTED,
             ha="left", va="top")


def headroom(ax, factor=1.38, log=False):
    """Room above the data for the legend and the notes in the upper corners."""
    lo, hi = ax.get_ylim()
    if log:
        ax.set_ylim(lo, hi * factor ** 4)
    else:
        ax.set_ylim(lo, lo + (hi - lo) * factor)


def title(ax, text, sub=None):
    """Left-aligned title, with an optional second line in secondary ink."""
    ax.set_title(text if sub is None else f"{text}\n{sub}", loc="left", fontsize=9)


def poisson_band(counts):
    """Approximate 68 % error bars for counts (Gehrels-like: sqrt(n + 0.75) + 1 above, sqrt(n - 0.25) below)."""
    c = np.asarray(counts, dtype=np.float64)
    up = np.sqrt(c + 0.75) + 1.0
    down = np.where(c > 0, np.sqrt(np.maximum(c - 0.25, 0)), 0.0)
    return down, up
