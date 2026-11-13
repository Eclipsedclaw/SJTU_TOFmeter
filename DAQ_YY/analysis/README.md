# Paddle muon analysis

From the raw waveforms of a DAQ_YY run to reconstructed cosmic muons: pulses, hits at both ends of every
scintillator paddle, events, calibration from the run itself, straight tracks with angles and time of
flight, a toy-simulation expectation, plots and a summary. Made for self-trigger runs of paddles read out
at both ends, starting with the two-layer test of 2026-11-12 (two top and two middle paddles on board 1,
channels 0–7). Every paddle, channel, layer and board comes from a setup file, so adding channels means
editing that file, not the code.

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

A 30-minute run (125 000 waveforms) takes about 4 minutes on one CPU core. **Passes if** the last lines
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
| 6 | `event_display.py` | → `event_displays/` | Pictures of single events (first 10 tracks; `--events 12,57` for others) |

The ideas behind each step (the scripts' `--help` has the details):

1. **Pulse.** The front end integrates the charge: a pulse rises in ~10 ns, stays flat for ~20 ns and
   decays over ~1 µs. Its *height* (running mean of 5 samples at the peak) measures the light. Its
   *time* is where the rising edge crosses 20 % of the height ("constant fraction"), interpolated between
   samples: the same point on every pulse, small or big, so no time walk. Noise bursts (an oscillation
   that crossed the threshold) have no flat top and swing below the baseline; they are flagged and left out.
2. **Hit.** A muon crossing a paddle at u (cm from its centre, + towards end A) gives
   tA − tB = −2u/v + offset (v = signal speed along the paddle) and a mean time (tA + tB)/2 that does
   not depend on u. A pulse without a partner at the other end is a one-end hit.
3. **Calibration from the run itself.** Muons cross a paddle evenly along its length, so tA − tB spreads
   evenly over 2L/v: a box, smeared by the time resolution. Fitting it gives v, the offset and the
   resolution. Light falls as exp(−d/λ) towards each end, so ln(hA/hB) = 2u/λ + ln(gA/gB) gives the
   attenuation length λ and the ends' gain ratio. The muon peak (most probable √(hA·hB), which does not
   depend on u) comes from the muons that cross both layers, × cos(zenith). Layer time offsets come from
   requiring β = 1 on average (cosmic muons are relativistic).
4. **Track.** The position along each paddle comes from the timing; across it, the paddle's centre line
   (resolution width/√12). A straight line through the points of the layers gives the zenith θ and the
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
- **Settings** (`[pulse]`, `[events]`, `[calibration]`, `[flux]`) are explained line by line in the file;
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
| `event_displays/` | `event_<n>.png` (top and side views) and `event_<n>_waveforms.png` |

The ROOT trees can be browsed directly, e.g. `root -l tracks.root` then `tracks->Draw("theta_deg")`.

| Plot | Look for |
|---|---|
| `pulse_rates`, `baseline_noise` | Channels with many rejected (noise) pulses or a noisy baseline |
| `pulse_shape_<paddle>` | Both ends with the same shape |
| `heights_<paddle>` | Where the thresholds cut the spectrum |
| `tdiff_<paddle>` | A box with smeared edges; the fit gives v, offset, resolution |
| `attenuation_<paddle>` | ln(hA/hB) rising with u; the S-shaped curve is what the fitted line predicts at a measured u |
| `mip_<paddle>` | The muon peak of the tracks (blue) against all two-end hits (grey, mostly not muons) |
| `layer_delay_<a>_<b>` | One peak at 0 ns above a flat chance level |
| `event_types`, `coincidences_<a>_<b>`, `rates_vs_time` | Counts and rates; tracks per paddle pair against the toy |
| `hitmap_<layer>`, `zenith`, `azimuth` | Where the tracks cross, and their angles against the toy |
| `tof_residual`, `inv_beta` | Time of flight minus s/c against the expected resolution; 1/β around 1 |
| `energy_deposit` | Signal per vertical path, peaking at 1 MIP in every layer |

## Read the results

`analysis_summary.txt` starts with a block like this (run 58):

```text
Coincid.  : top+middle: 640 events = 0.356 ± 0.014 Hz (by chance 0.013)
Tracks    : 640; zenith mean 25.4° (toy 24.3°)
Expected  : 730.8 tracks for fully efficient paddles; measured / expected = 0.88 ± 0.03
TOF       : tof - s/c robust σ 0.37 ns (expected 0.41 ns) over a median path of 33.4 cm
Map check : top layer y (timing) by middle paddle: slope -0.24 ± 0.45 (+1 expected), ...
```

- **threshold/peak** (per channel, in the paddle lines): the threshold as a fraction of that end's signal
  for a muon through the paddle's centre. A muon near the far end gives only exp(−L/2λ) of it (≈ 0.5 here),
  so thresholds above ~0.3 lose muons near the paddle ends.
- **both-ends share** (`events_summary.txt`): the fraction of a paddle's hits seen at both ends.
- **measured / expected**: the product of the paddles' efficiencies and of the true flux over the assumed
  one (70 /(m² s sr), cos²θ, sea level).
- **Speed warning**: with one signal speed for all paddles (`v_eff_mode = "shared"`), the summary also
  shows each paddle's own fit; if they disagree, see *Known limitations*.
- **Map check**: where a muon crosses the top layer along its paddles depends on which middle paddle it
  then crosses (and the other way round). The data's mean positions are fitted against the toy's as
  shift + slope × toy: slope +1 = the channel map is as the setup file says; −1 = one direction (the A/B
  cables) or one order (which paddle is 1 or 2) is the other way round. It uses the timing and,
  independently, the light sharing. A simulated run with the top cables swapped gives −1.4 ± 0.7.

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
near 0. On 2026-11-13: signal speed 15.98 ± 0.18 (truth 16.00) cm/ns, λ 37.8–38.3 (38.0) cm, offsets and
muon peaks within 1.5σ, positions unbiased (rms 1.1 cm across a paddle, 5.3–5.5 cm along it).

## Reference: runs 56–58 (2026-11-12/13, board 1, thresholds 400 ADC; run 56: 600)

| | run 56 (40 s) | run 57 (10 min) | run 58 (30 min, longer window) |
|---|---|---|---|
| Tracks (top + middle) | 11 | 205 | 640 |
| Rate | 0.28 ± 0.08 Hz | 0.342 ± 0.024 Hz | 0.356 ± 0.014 Hz |
| Measured / expected | 0.68 ± 0.21 | 0.84 ± 0.06 | 0.88 ± 0.03 |
| Signal speed (shared) | 16.8 ± 0.7 cm/ns | 15.3 ± 0.1 cm/ns | 15.24 ± 0.06 cm/ns |
| σ(tA − tB) | 0.50–0.95 ns | 0.55–0.62 ns | 0.51–0.61 ns |
| Position resolution along a paddle | 4–8 cm | 4.2–4.8 cm | 3.9–4.7 cm |
| Attenuation length | 27–34 cm | 27–31 cm | 27–32 cm |
| Muon peak √(hA·hB) | 775–981 ADC | 649–886 ADC | 628–893 ADC |
| TOF residual (robust σ) | 0.26 ns | 0.45 ns | 0.37 ns |

- Two-end hits per paddle are 2–4 times the expected muon rate: single paddles also see electrons, gammas
  and radioactivity, mostly below the muon peak (`mip_<paddle>.png`).
- ch6 (MID2-A) flags up to 31 % of its triggers as noise bursts, ch2 (TOP2-A) 14 % (run 57).
- Thresholds were 0.43–0.69 of the muon signal at the paddle centre (run 58); TOP2-A is the highest.
- Map check (runs 57, 58): the middle layer follows TOP1/TOP2 as expected (slope +1.3 ± 0.5); the top
  layer's position along its paddles does **not** follow MID1/MID2 (slopes −0.24 ± 0.45 timing,
  −0.08 ± 0.42 light; +1 excluded at 2.6–2.8σ). To check: are TOP1A/TOP2A (ch0, ch2) at the wall, and is
  BOTTOM1 (ch4/ch5) the middle paddle nearer the wall? Also the top hits sit 1.2 ± 0.4 cm towards −y.

## Known limitations

- **Geometry not measured**: sizes and the 30 cm between the layers are the simulation design. The
  expected rates, angles, path lengths and β depend on them.
- **End losses bias the signal speed.** With thresholds near the muon signal the far end misses muons near
  the paddle ends, the tA − tB box looks narrower and v comes out high, more for higher thresholds (run 58:
  TOP2 16.9, MID1 14.5 cm/ns on their own). The shared speed averages this out; positions near the ends
  can be off by about (Δv/v) × L/2 ≈ 1–2 cm. Lower thresholds remove it.
- The box fit assumes muons evenly along the paddle. Low-energy hits (mostly near the centre) are removed
  by a cut on √(hA·hB) > 0.8 × the muon peak; without it a simulated run gave v 7 % too high.
- With two layers the track goes through both points, so nothing checks a single hit; a third layer
  would give a χ² per track and a direct calibration of the middle layer's positions.
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
