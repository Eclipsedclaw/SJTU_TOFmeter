# Paddle muon analysis

From the raw waveforms of a DAQ_YY run to reconstructed cosmic muons: pulses, hits at both ends of every
scintillator paddle, events, calibration from the run itself, straight tracks with angles and time of
flight, a toy-simulation expectation, plots and a summary. Made for self-trigger runs of paddles read out
at both ends, starting with the two-layer test of 2026-11-12 (two top and two middle paddles on board 1,
channels 0–7). Every paddle, channel, layer and board comes from a setup file, so adding channels means
editing that file, not the code.

**No known muon values in the data processing.** Every number the analysis applies to the data is fitted
from the data themselves (speeds, offsets, resolutions, attenuation, signal peaks) or comes from the setup
file (geometry) or from signal processing (filter, timing). What muons are known to do (the cosmic flux and
its cos²θ, β = 1, the energy loss of a minimum-ionising particle) enters only the toy simulation, and every
plot that shows it says so ("simulation"). The setup is a test and may be wrong; an analysis that forced it
towards textbook values would hide that.

## Quick start

```bash
cd SJTU_TOFmeter/DAQ_YY/analysis
python3 run_analysis.py
```

```text
ROOT file (raw data of the run): /data/data/jiancheng/DAQ_yy_data/raw/run00058_file0000_20261113_132310.root
Output folder [/data/data/jiancheng/DAQ_yy_data/raw/run00058_file0000_20261113_132310_muon]:
```

Tab completes paths; Enter keeps the folder in brackets, a folder named after the run next to the data
file. To skip the questions:

```bash
python3 run_analysis.py /path/to/run.root                 # output next to the data
python3 run_analysis.py /path/to/run.root --out /path/to/folder
```

The results of runs 56–58 are in `/data/data/jiancheng/DAQ_yy_data/raw/<run>_muon/`, next to their data.
A 30-minute run (125 000 waveforms) takes about 6 minutes on one CPU core. **Passes if** the last lines
read `Done in ...` and `Start with .../analysis_summary.txt`.

Before the first run, check that Python finds ROOT (as for the scripts in [`../scripts`](../scripts)):
`python3 -c "import ROOT, numpy, scipy, matplotlib; print('OK')"`. If ROOT is missing,
`source /opt/root/root_install/bin/thisroot.sh`. The setup file needs Python 3.11 or newer (`tomllib`).

## What the analysis does

`run_analysis.py` runs six scripts one after the other; each can also be run on its own on the output
folder (they ask for it, or take it as the first argument), for example after changing a setting.

| Step | Script | Reads → writes | What it does |
|---|---|---|---|
| 1 | `extract_pulses.py` | raw ROOT file → `pulses.root` | Measures every waveform: baseline, height, time, charge, shape tests (flags) |
| 2 | `build_events.py` | → `events.root` | Pairs the pulses of a paddle's two ends into hits; groups hits within 30 ns into events |
| 3 | `calibrate.py` | → `calibration.json` | Signal speed, offsets and resolution of tA − tB, attenuation, muon peak, layer time offsets |
| 4 | `reconstruct.py` | → `tracks.root` | Hit points, straight tracks, angles, time of flight; toy expectation; channel-map check |
| 5 | `make_plots.py` | → `plots/`, `analysis_summary.txt` | 30-odd plots, one per file, and the summary |
| 6 | `event_display.py` | → `event_displays/` | One page per event (views, readout heights, info table, waveforms): the first 10 tracks and 5 with several paddles of one layer (`--events 12,57`, `--several N`) |

The ideas behind each step (the scripts' `--help` has the details):

1. **Pulse.** The front end integrates the charge: a pulse rises in ~10 ns, stays flat for ~20 ns and
   decays over ~1 µs. The waveform first goes through a low-pass filter (Butterworth, 80 MHz, run forwards
   and backwards so that no time moves) that removes the ringing at 110–190 MHz; see *The low-pass filter*.
   The *height* (running mean of 5 samples at the peak) measures the light. The *time* is where the rising
   edge crosses 20 % of the height ("constant fraction"), interpolated between samples: the same point on
   every pulse, small or big, so no time walk. Noise bursts (an oscillation that crossed the threshold) have
   no flat top and swing below the baseline; they are flagged and left out.
2. **Hit.** A muon crossing a paddle at u (cm from its centre, + towards end A) gives
   tA − tB = −2u/v + offset (v = signal speed along the paddle) and a mean time (tA + tB)/2 that does
   not depend on u. A pulse without a partner at the other end is a one-end hit.
3. **Calibration from the run itself.** Particles cross a paddle evenly along its length, so tA − tB
   spreads evenly over 2L/v: a box, smeared by the time resolution. Fitting its edges gives v (each paddle
   its own), the offset and the resolution. Light falls as exp(−d/λ) towards each end, so
   ln(hA/hB) = 2u/λ + ln(gA/gB) gives the attenuation length λ and the ends' gain ratio. The paddle's
   signal peak (most probable √(hA·hB), which does not depend on u) is taken from the tracks' hits as
   measured. No time offsets are applied between the layers: the paddles' cable delays would need a
   pulser measurement; the offsets that β = 1 would need are only quoted as text.
4. **Track.** The position along each paddle comes from the timing; across it, the paddle's centre line
   (resolution width/√12). The measured position scatters by the resolution (about 4 cm) and so lands
   beyond a paddle's end for muons near it; the point used is the best estimate *inside* the paddle (see
   *Positions beyond the paddle ends*). A straight line through the points of the layers gives the zenith θ and the
   azimuth φ of the arrival direction; the time between the layers gives the time of flight and β.
   A toy simulation (muons with I₀ cos²θ through the same geometry, reconstructed like the data) gives the
   rates and angle distributions expected for fully efficient paddles.

## The setup file

[`setups/two_layer_2026-11-12.toml`](setups/two_layer_2026-11-12.toml) describes the first test: which
channel reads which end of which paddle, where the paddles are, and the analysis settings. The analysis
uses the setup file whose `valid_from` is the latest on or before the run's date (from the file name).
**Whenever paddles, cables or positions change, copy the file, edit it and give it a new `valid_from`**;
keep the old one so old runs are still analysed with their own setup. `--setup file.toml` overrides the
choice.

Coordinates of the first setup, seen from above (slide 17, the wall at the top):

```text
                 wall  (+y)
            TOP1A(ch0)  TOP2A(ch2)
               |  |       |  |
 MID1A(ch4) ===|==|=======|==|=== MID1B(ch5)      x to the right (MID A ends on the left)
 MID2A(ch6) ===|==|=======|==|=== MID2B(ch7)      z up; the top layer 30 cm above the middle one
               |  |       |  |                    the layers cross at x = y = 0
            TOP1B(ch1)  TOP2B(ch3)
```

One `[[paddles]]` table per paddle:

```toml
[[paddles]]
name = "TOP1"
layer = "top"                        # one of the [[layers]] (each with name and z_cm)
axis = "y"                           # long direction: "x" or "y"
x_cm = -2.0                          # centre (z_cm too, if not the layer's)
y_cm = 0.0
length_cm = 40.0                     # scintillator length: tA - tB spans 2 x length / v
width_cm = 4.0
thickness_cm = 0.6
A_end = "+y"                         # where the end read by channel A sits
channels = { A = 10010, B = 10011 }  # DAQ_YY id = 10000 + 10 x board + channel
```

- **More paddles:** add `[[paddles]]` tables; more layers: add `[[layers]]` (sorted by `z_cm`, top first).
  Paddles of one layer side by side are recognised as neighbours (a muon crossing both is merged into one
  point).
- **More boards:** the boards' clocks start together in ctrlmode 3, but each board has a fixed time offset
  (up to ±5 ns, measured with [`../scripts/analyze_sync.py`](../scripts/analyze_sync.py)). Enter them under
  `[board_offsets_ns]`, e.g. `2 = 4.766`; single cable delays under `[channel_offsets_ns]`.
- **Settings** (`[pulse]`, `[events]`, `[calibration]`, `[flux]`, `[display]` for the side of the room the
  event pages mark) are explained line by line in the file;
  every key is optional. The geometry of the first setup is the simulation design
  (`simulation/src/DetectorConstruction.cc`): replace it with the measured values.

## Outputs

In the output folder:

| File | Contents |
|---|---|
| `analysis_summary.txt` | **Start here.** The main numbers, then the summaries of steps 1–4 |
| `pulses.root` | Tree `pulses`: one entry per waveform (channel, time, height, charge, flags, ...); tree `shapes`: average pulse per channel |
| `events.root` | Tree `hits` (one per hit, with its event number); tree `events` (one per event) |
| `calibration.json` | Every constant, with how it was found; reusable for other runs (`run_analysis.py --calibration`) |
| `tracks.root` | Tree `tracks` (one per reconstructed muon); tree `toy` (the toy muons, with weights in Hz) |
| `*_summary.txt` | The tables each step printed |
| `plots/*.png` | One plot per file, listed below |
| `event_displays/` | `event_<n>.png`, one page per event: the setup from above and from two sides, each paddle end's readout box coloured by its pulse height (the light sharing at a glance), every hit at its own paddle, the fitted track as an arrow, the room's wall; an info table; the records and their rising edges divided by the height, where the time differences show |

The ROOT trees can be browsed directly, e.g. `root -l tracks.root` then `tracks->Draw("theta_deg")`.

| Plot | Look for |
|---|---|
| `pulse_rates`, `baseline_noise` | Channels with many rejected (noise) pulses or a noisy baseline |
| `pulse_shape_<paddle>` | Both ends with the same shape |
| `position_<paddle>` | Measured position along the paddle: how many lie beyond the ends, against the smearing alone |
| `heights_<paddle>` | Where the thresholds cut the spectrum |
| `tdiff_<paddle>` | A box with smeared edges; the fit gives v, offset, resolution |
| `attenuation_<paddle>` | ln(hA/hB) rising with u (see *Reading the plots*); the right axis is hA/hB itself |
| `tdiff_vs_ratio_<paddle>` | The same check with no calibration: raw ln(hA/hB) against raw tA − tB, 2D, with their correlation |
| `mip_<paddle>` | The signal peak of the tracks (blue) against all two-end hits (grey, mostly not muons) |
| `layer_delay_<a>_<b>` | One peak above a flat chance level (left, ±1000 ns) and its shape (right, ±10 ns) |
| `event_quality` | The selection steps: all events → a paddle at both ends → a track → one paddle per layer, or several paddles of one layer |
| `event_types`, `coincidences_<a>_<b>`, `rates_vs_time` | Counts and rates; tracks per paddle pair against the toy |
| `hitmap_<layer>`, `zenith`, `azimuth` | Where the tracks cross, and their angles against the toy |
| `tof_residual` | Measured time of flight minus s/c, against the simulation of β = 1 with the measured resolution |
| `beta` | β over its whole measured range; shaded where no particle going down can be (β > 1, β < 0) |
| `energy_deposit` | The tracks' signal per paddle as measured (ADC); each paddle has its own scale |
| `energy_deposit_sim` | Simulation: the energy muons lose in the paddles (MeV, PDG Landau formula) |
| `signal_vs_path` | Signal against the path through the paddle (1/cos θ): does it grow with the path? |

## Read the results

`analysis_summary.txt` starts with a block like this (run 58):

```text
Coincid.  : top+middle: 640 events = 0.356 ± 0.014 Hz (by chance 0.013)
Tracks    : 640; zenith mean 23.5° (toy 22.8°)
Expected  : 730.8 tracks for fully efficient paddles; measured / expected = 0.88 ± 0.03
TOF       : measured tof - s/c median +0.01 ns, robust σ 0.34 ns (simulation of β = 1 with the measured
            resolution: σ 0.35 ns) over a median path of 33.0 cm; no layer time offsets applied
β         : median 0.98; 309 of 640 tracks measured at β ≥ 1, 6 below 0 (time resolution)
Map check : top layer y (timing) by middle paddle: slope -0.53 ± 0.44 (+1 expected), ...
```

- **Event quality**: each selection step with its count and share of all events. In run 58: 80 503
  events (hits within 30 ns), 47.8 % with a paddle seen at both ends, 640 tracks (0.8 %), of which 616 have
  one paddle per layer and 24 two neighbouring paddles of one layer (a muon crossing at their joint, or a
  shower); none had two paddles of a layer that are not neighbours.
- **threshold/peak** (per channel, in the paddle lines): the threshold as a fraction of that end's signal
  for a muon through the paddle's centre. A muon near the far end gives only exp(−L/2λ) of it (≈ 0.5 here),
  so thresholds above ~0.3 lose muons near the paddle ends.
- **both-ends share** (`events_summary.txt`): the fraction of a paddle's hits seen at both ends.
- **measured / expected**: the product of the paddles' efficiencies and of the true flux over the assumed
  one (70 /(m² s sr), cos²θ, sea level).
- **Speeds**: every paddle keeps its own signal speed (`v_eff_mode = "per_paddle"`; run 58 rejects one
  common speed). With `"shared"` the summary also shows each paddle's own fit and warns if they disagree.
- **Map check**: where a muon crosses the top layer along its paddles depends on which middle paddle it
  then crosses (and the other way round). The data's mean positions are fitted against the toy's as
  shift + slope × toy: slope +1 = the channel map is as the setup file says; −1 = one direction (the A/B
  cables) or one order (which paddle is 1 or 2) is the other way round. It uses the timing and,
  independently, the light sharing. A simulated run with the top cables swapped gives −1.4 ± 0.7.

## Reading the plots

- **Attenuation** (`attenuation_<paddle>`). A hit at u (cm from the centre, + towards end A) is L/2 − u
  from end A and L/2 + u from end B. The light reaching an end falls as exp(−distance/λ), so
  hA ∝ gA·exp(−(L/2 − u)/λ) and hB ∝ gB·exp(−(L/2 + u)/λ), and the energy deposit cancels in their ratio:
  ln(hA/hB) = 2u/λ + ln(gA/gB). The y axis is that log ratio (right axis: the ratio hA/hB itself), the x
  axis the position from the timing. The slope of the straight line gives λ (λ ≈ 30 cm here: the light
  halves every λ ln 2 ≈ 21 cm), its value at u = 0 the gain ratio of the ends. λ does **not** come from the
  edges of the tA − tB box: those give only the signal speed v, which sets the position scale. The fitted
  line looks S-shaped because measured positions beyond the ends come from hits near the ends.
  `tdiff_vs_ratio_<paddle>` makes the same comparison with nothing calibrated: the raw height ratio
  against the raw time difference. A hit nearer end A gives A more light and gives it earlier, so the two
  must follow each other if the light fades along the paddle. Runs 56–58: rank correlation −0.89 to −0.95
  in every paddle (−0.93 in run 58), and the medians lie on a straight line that flattens only beyond the
  ends. The scatter of ln(hA/hB) around the attenuation line (0.24–0.27) is about what the timing's own
  position error alone gives (2/λ × 4 cm ≈ 0.27), so the light sharing places a hit at least as well as the
  timing; telling the two apart needs a third measurement (a third layer's track, or a source).
- **Coincidences** (`coincidences_<a>_<b>`). Each cell: tracks measured through that paddle pair, the
  simulation's number, and their ratio. The simulation sends cosmic muons with the sea-level flux
  I = I₀ cos²θ, I₀ = 70 /(m² s sr) (PDG, muons above ~1 GeV) through the setup's paddles in straight lines,
  every paddle fully efficient, for the run's duration. The ratio is the pair's efficiency times the true
  over the assumed flux (the lab's floors above, the altitude and the soft component change the flux).
- **Two paddles of one layer**. The bars touch; a track crossing a layer near the joint of two bars at an
  angle goes through both inside the 0.6 cm thickness. The simulation predicts this for 4.2 % of the tracks
  (27 in run 58); 24 were measured. The hit map and the event pictures show each bar's own hit; the
  track's point for the fit is their mean, weighted by the signal.
- **β** (`beta`). β = path / (c × time of flight). The flight between the layers lasts ~1.1 ns and is
  measured with σ ≈ 0.35 ns, so half the tracks come out at β > 1: not faster than light, but measured too
  fast. The plot shows every track from −1 to 3 (the outer bins hold the few beyond), shaded where no
  particle going down can be: β > 1, and β < 0, where the lower layer fired first (a particle going up, or
  the resolution). The simulation of β = 1 particles with the same resolution lies over the data. The full
  TOF-meter (top to bottom 70 cm) and better timing will separate slow particles.
- **Energy and the track angle** (`energy_deposit`, `signal_vs_path`). A minimum-ionising particle loses
  energy at a nearly constant rate per cm of path (dE/dx ≈ 2 MeV/cm in plastic), and that rate does not
  depend on the angle. The energy deposited in a slab does: a track at zenith θ crosses thickness / cos θ
  of it, so it deposits 1/cos θ times more (1.41 times at 45°). `signal_vs_path` shows this on the data:
  the signal rises with 1/cos θ (run 58: from 1.01 to 1.30 times the peak; "∝ path" χ² 10 for 7,
  "no dependence" χ² 105 for 7). Without the angle the peak of the tracks is ~10 % higher and 15–25 %
  wider than for vertical tracks; the analysis quotes it as measured and leaves the angle correction to
  plots that say so.

## The low-pass filter

Most channels ring at 110–190 MHz (2–6 times more noise power there than at 2–40 MHz; ch6 also at
50 MHz), far above the 35–50 MHz a 10 ns edge needs. Run 58 analysed with each filter setting (σ(tA − tB)
errors ±0.015 ns):

| Filter | σ(tA − tB), ns: TOP1 / TOP2 / MID1 / MID2 | TOF σ | Positions > 3σ beyond the ends |
|---|---|---|---|
| none | 0.509 / 0.596 / 0.566 / 0.613 | 0.373 ns | 0.23 % |
| Gaussian 120 MHz | 0.499 / 0.546 / 0.546 / 0.563 | 0.356 ns | 0.14 % |
| Gaussian 80 MHz | 0.495 / 0.522 / 0.537 / 0.519 | 0.360 ns | 0.12 % |
| Gaussian 40 MHz | 0.528 / 0.528 / 0.550 / 0.544 | 0.366 ns | 0.08 % |
| **Butterworth 80 MHz** (used) | **0.504 / 0.502 / 0.534 / 0.505** | **0.338 ns** | 0.14 % |
| Gaussian 60 MHz, time from the raw samples | 0.510 / 0.592 / 0.566 / 0.607 | 0.374 ns | 0.21 % |

The filter improves the timing, most for the noisy channels (TOP2 and MID2 by 18 %), and halves the far
position outliers. Below ~60 MHz it starts to slow the edge (TOP1 gets worse). Both filters are run forwards
and backwards (zero phase), so they move no pulse in time; the baseline noise in `pulses_summary.txt` is
always measured on the raw samples. Settings: `[pulse]` `lowpass_mhz` (0 = off), `lowpass_kind`,
`lowpass_timing`; or `--lowpass 0` to compare without it. A simulated run (white noise, no ringing) gives
the same constants with and without the filter.

## Positions beyond the paddle ends

A muon cannot cross outside a paddle, but its *measured* position can lie beyond the end: the position comes
from tA − tB with a resolution of about 4 cm, so a muon 1 cm from the end is measured beyond it almost half
the time. For muons spread evenly along a 40 cm paddle that puts 2σ/(√(2π) L) ≈ 8 % of the hits beyond the
ends, and that is what the data show (`position_<paddle>.png`, run 58: 7.5–8.5 % measured, 7.5–8.0 % from
the resolution). Far beyond (more than 3σ) Gaussian smearing gives 0.01 %; there were 0.2 % without the
filter, small pulses whose edge time the ringing had pulled, and 0.06–0.13 % with it.

The track points therefore use the best estimate given that the muon crossed inside: the measurement's
Gaussian cut to the paddle, its mean as the point and its spread as the error. It always lies inside; in a
simulated run it and the filter together brought the position error along the paddle from 5.4 to 4.6 cm and
the zenith error from 8.3° to 7.3°. The measured position is kept as `u_meas_<layer>` in `tracks.root`, and the flag
OUTSIDE marks a measurement more than 3σ beyond its paddle.

## Check the analysis with a simulated run

[`tests/make_fake_run.py`](tests/make_fake_run.py) simulates a run with known constants (cosmic muons and a
low-energy background through the setup's paddles, the measured pulse shape, noise, self-trigger records)
and writes it in the raw2root_new format; [`tests/check_fake_run.py`](tests/check_fake_run.py) compares the
analysis with the truth:

```bash
python3 tests/make_fake_run.py /path/to/test_folder                  # 600 s; --swap-top-ends to test the map check
python3 run_analysis.py /path/to/test_folder/fake_run_20261112_120000.root --out /path/to/test_folder/fake_muon
python3 tests/check_fake_run.py /path/to/test_folder/fake_muon
```

**Passes if** the constants agree with the truth within about 2σ and the position residuals have a mean
near 0. On 2026-11-13 (filter on, no assumptions): each paddle's signal speed 15.6–16.1 (truth 16.0) cm/ns,
λ 37.4–38.4 (38.0) cm, offsets and signal peaks within 1.6σ (the peak as measured equals the truth times the
tracks' path factor, ×1.09–1.10), positions unbiased (rms 1.1 cm across a paddle, 4.4–4.8 cm along it),
zenith −0.9 ± 7.3°. With no time offsets applied, the measured time of flight carries the simulated cables'
delays (+0.11 ns); after them +0.04 ns remain.

## Reference: runs 56–58 (2026-11-12/13, board 1, thresholds 400 ADC; run 56: 600; filter on)

| | run 56 (40 s) | run 57 (10 min) | run 58 (30 min, longer window) |
|---|---|---|---|
| Events / with a paddle at both ends | 1 374 / 39 % | 26 968 / 49 % | 80 503 / 48 % |
| Tracks (top + middle) | 11 | 205 | 640 |
| Rate | 0.28 ± 0.08 Hz | 0.342 ± 0.024 Hz | 0.356 ± 0.014 Hz |
| Measured / expected | 0.68 ± 0.21 | 0.84 ± 0.06 | 0.88 ± 0.03 |
| Signal speed (each paddle) | 16.0–21.0 cm/ns (too few hits) | 14.6–16.6 cm/ns | 14.4–16.5 cm/ns |
| σ(tA − tB) | 0.43–1.12 ns | 0.46–0.62 ns | 0.42–0.61 ns |
| Position resolution along a paddle | 3.5–12 cm | 3.4–5.1 cm | 3.0–5.0 cm |
| Signal peak of the tracks √(hA·hB) | 796–982 ADC (all hits) | 671–964 ADC | 694–991 ADC |
| TOF − s/c, median and robust σ (no offsets) | −0.16, 0.52 ns (11 tracks) | +0.07, 0.36 ns | +0.01, 0.34 ns |
| Tracks measured at β ≥ 1 | 6 of 11 | 86 of 205 | 309 of 640 |
| Zenith mean (toy) | 23.8° (21.6°) | 24.5° (22.7°) | 23.5° (22.8°) |

- Two-end hits per paddle are 2–4 times the expected muon rate: single paddles also see electrons, gammas
  and radioactivity, mostly below the muon peak (`mip_<paddle>.png`).
- ch6 (MID2-A) flags up to 31 % of its triggers as noise bursts, ch2 (TOP2-A) 14 % (run 57).
- Thresholds were 0.43–0.69 of the muon signal at the paddle centre (run 58); TOP2-A is the highest.
- The amplifiers clip at raw 14 600–15 100 ADC (6300–6900 above the baseline): a few dozen pulses per
  channel above ~7 muon peaks pile up there. `max_raw_adc = 14500` flags them SATURATED.
- Map check (run 58): the middle layer follows TOP1/TOP2 as expected (slope +1.45 ± 0.45); the top layer's
  position along its paddles does **not** follow MID1/MID2 (slopes −0.52 ± 0.43 timing, −0.10 ± 0.45
  light; +1 excluded at 3.5 and 2.4σ). To check: are TOP1A/TOP2A (ch0, ch2) at the wall, and is BOTTOM1
  (ch4/ch5) the middle paddle nearer the wall? Also the top hits sit 1.1 ± 0.4 cm towards −y.

## Known limitations

- **Geometry not measured**: sizes and the 30 cm between the layers are the simulation design. The
  expected rates, angles, path lengths and β depend on them.
- **End losses can bias the signal speed.** With thresholds near the muon signal the far end misses
  particles near the paddle ends, the tA − tB box looks narrower and v comes out high, more for higher
  thresholds (run 58: TOP2 16.5 with the highest threshold, MID1 14.4 cm/ns with the lowest). Positions
  near the ends can then be off by about (Δv/v) × L/2 ≈ 1–2 cm. Lower thresholds remove it.
- **No time offsets between the paddles.** Their cable delays are unknown; the time of flight and β are
  as measured (the delays happen to be small here: TOF − s/c median +0.01 ns in run 58). A pulser into all
  channels would measure them without assuming anything about the particles.
- The box fit assumes muons evenly along the paddle. Low-energy hits (mostly near the centre) are removed
  by a cut on √(hA·hB) > 0.8 × the muon peak; without it a simulated run gave v 7 % too high.
- With two layers the track goes through both points, so nothing checks a single hit; a third layer
  would give a χ² per track and a direct calibration of the middle layer's positions.
- The filter was chosen on run 58's ringing; another front end (faster pulses) needs its own scan.
- Checked on self-trigger runs (trigmode 0) only. Other trigger modes may place the trigger elsewhere in
  the record.
- Reads DAQ_YY files only (raw2root_new format). Another DAQ needs its own reader for step 1.
- The toy has no efficiencies, no soft component and no building above the lab.

## Units and conventions

- Times in ns, lengths in cm, angles in degrees, speed as β = v/c; pulse heights in ADC counts above the
  baseline, charges in ADC × ns. Samples and ttt ticks are 2 ns.
- u: position along a paddle from its centre, positive towards the end that channel A reads.
- Zenith θ from the vertical; azimuth φ of the direction the muon came from, 0 = from +x, 90° = from +y.
- Pulse flags (`pulses.root`, column `flags`; explained in `pulses_summary.txt`): SHORT, NO_EDGE, NOISE,
  LOW are rejected; SATURATED, PILEUP, SLOPE, NOISY are kept.
- Track flags (`tracks.root`): ADJACENT (hits in neighbouring paddles merged), AMBIGUOUS (hits in
  separate paddles of one layer: the biggest used), OUTSIDE (a position > 3σ beyond a paddle's end).
- Python + PyROOT, one thread; progress bars on stderr; the plot style of `../scripts/daqyy_common.py`,
  which every script imports (keep the folders together).
