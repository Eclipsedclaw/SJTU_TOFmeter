# TOF-meter Geant4 simulation

Cosmic-muon TOF-meter: plastic-bar tracker, degrader, target and the three-layer Compton camera
(Sun et al., arXiv:2608.21216). The geometry is `src/DetectorConstruction.cc`; the run framework
(`TOFmeter.cc`, `TOFmeter*` classes) derives from Tsuguo Aramaki's "Original" framework (ASCII step
output, `/update`, `/OutputFile` ...), extended with a cosmic-ray generator (EXPACS spectra of all
species, or a muon formula) and a ROOT event ntuple. Tested with Geant4 11.4.3 (sequential run manager).

## Build and run

```sh
cmake -S SJTU_TOFmeter/simulation -B TOFmeter-work \
      -DTOF_GEANT4_DATA_DIR=/Users/jiancheng/tool/geant/geant4-11.4.3-build/data
cmake --build TOFmeter-work -j8
cd TOFmeter-work
./TOFmeter macro/smoke.mac          # macros are copied into the build directory
```

`TOF_GEANT4_DATA_DIR` is optional: it is used only when `GEANT4_DATA_DIR` is not set in the
environment (the local Geant4 install has its datasets in the build tree, not the install prefix).
The physics list is QGSP_BERT; set `PHYSLIST` to choose another reference list, e.g.
`PHYSLIST=FTFP_BERT_EMZ ./TOFmeter ...` (EM option 4). Avoid `*_HP` and `Shielding` until the
G4NDL dataset is installed (its download is incomplete). Without arguments `./TOFmeter` opens the
Qt GUI (`/control/execute macro/vis_qt.mac` fills the viewer); in batch mode pictures go to files
(`macro/vis*.mac`).

| Macro | What it does |
|---|---|
| `smoke.mac` | 2000 vertical 1 GeV mu- over the bar stations; ASCII + ROOT |
| `cosmic.mac` | 10^6 cosmic muons through the setup as written (~14 min live time) |
| `cosmic_expacs.mac` | 10^6 cosmic rays of all EXPACS species (neutrons, photons, e+-, ...) |
| `camera_validation.mac` | camera alone looking up, 2.7x10^6 muons (~15.5 h); compare with `analysis/camera_validation.C` |
| `stopping_example.mac` | mu- below 0.5 GeV through the example layout; summarise with `analysis/stops_summary.C` |
| `vis_geometry.mac` | geometry check, no particles: position table, volume tree, 5 PNG views, VRML file |
| `vis_qt.mac` | in the Qt GUI: geometry and 10 cosmic muons; the view turns slowly until you click, scroll or press a key in it (`/tof/vis/spin`) |
| `acceptance.mac` | cosmic-muon rates and geometric acceptance A*Omega per layer and coincidence (~5 min) |
| `muon_spectrum_check.mac` | generated muon energy spectrum vs. the EXPACS table, the Guan formula and the scan's band fluxes; integral numbers vs. PDG (`analysis/muon_spectrum.C`, ~1 min) |
| `lead_energy_scan.mac` | top & middle & Compton versus muon energy (20 points, 0.1-100 GeV) and top-lead thickness (0-20 cm); in parallel with `scripts/lead_energy_scan.sh` |
| `lead_lowE_scan.mac` | the same for 100-500 MeV in 20 MeV steps; in parallel with `sh scripts/lead_energy_scan.sh macro/lead_lowE_scan_def.mac` |
| `vis.mac`, `vis_camera.mac`, `vis_example.mac` | PNG pictures (oblique and side view) with 30-40 muons |
| `geom_camera_only.mac`, `geom_tofmeter_example.mac` | geometry setups used by the macros above |
| `test.mac` | the example of the original framework (antiprotons, ASCII) |
| `README.mac` | command reference |

The ROOT macros in `analysis/` are copied to the build directory as well, e.g.
`root -l -b -q 'analysis/camera_validation.C("output/camera_validation.root")'`.
Files added to `macro/`, `analysis/` or `scripts/` need a re-run of `cmake` to be copied.

## Geometry

Default = `DetectorConstruction.cc`, z up (muons come from +z), cm, from the top:

| z range | Part | Switch |
|---|---|---|
| 21 to 31 | lead layer, Pb 90x90x10 out to the lead walls, so the top is fully covered (`AbsorberPhys`); `/tof/det/leadThickness` changes the 10 cm, the bottom stays at 21 | `/tof/det/absorber` |
| 20 to 21 | stainless-steel plate 80x80x1 | `/tof/det/frame` |
| -0.2 to 0.4 | top TOF: station 0, 9 bars along x | `/tof/det/stations` |
| -2.5 to -2.2 | plate Al1 under it | `/tof/det/frame` |
| -30.2 to -29.6 | middle TOF: station 1, bars along y | |
| -32.5 to -32.2 | plate Al2 | |
| -52.2 to -38.4 | Compton detector: Si/YSO stack, 20x20 cm, Si 6 mm x3 and YSO 3 mm x2 alternating | `/tof/det/stack` |
| -52.5 to -52.2 | plate Al3 under the stack | |
| -67.7 to -67.1 | bottom TOF: station 2, bars along x | |
| -70.0 to -69.7 | plate Al4 | |
| -70 to 21 | four stainless L beams (10 cm legs, 1 cm thick) on the outside of the steel-plate corners | `/tof/det/frame` |
| -70 to 31 | four 10 cm lead walls, inner faces at +-45 cm (5 cm from the steel-plate edges), tops level with the lead layer | `/tof/det/leadWalls` |

Bars are 40 x 4 x 0.6 cm plastic (G4_PLASTIC_SC_VINYLTOLUENE), 2 cm above their plate; the four
plates are 80 x 62 x 0.3 cm G4_Al between the corner L beams (long side = distance between the
beams' inner faces, short side = gap between the leg tips; Al2 turned by 90 deg: 62 x 80). Changes relative to the original file: z is mirrored (that file
had z pointing down, with the lead at negative z); the world is air and 4 x 4 x 3 m (was vacuum,
1 x 1 x 1.6 m) so that the cosmic source plane fits;
the 9 bars of a station sit in an air envelope, so every bar has its own ID; the world has copy
number -1. The Si-YSO face-to-face gap is 2.85 cm (the code's pitch adds only the YSO thickness).

Optional parts (off by default, set with macro commands then `/update`):
- **Compton camera** (`/tof/det/camera true`): CH2 and CH1 YSO 3 mm, CH0 LYSO 6 mm (7.2 g/cm3),
  all 47.04 mm square, 25 / 35 mm centre to centre. `/tof/camera/position` is the centre of CH2;
  `/tof/camera/rotateX` turns it (0 = looking up). Pixels, wrapping, electronics and the Al cover
  are not modelled.
- **Target** (`/tof/target/...`) and **degrader** (`/tof/degrader/...`): boxes with NIST or custom
  (YSO, LYSO, CdZnTe) materials.
- **Blocks** (`/tof/blocks/add name material x y z dx dy dz cm`): any passive box, e.g. lead bricks.

`macro/geom_tofmeter_example.mac` shows them in an illustrative (not optimised) layout: a
20x20x2 cm Fe target at z = -36.5 to -34.5 cm under the middle TOF, and the paper's camera in place
of the Si/YSO stack (CH2 at z = -40 cm). The lead layer on top acts as the degrader.
Check any setup with `macro/vis_geometry.mac` (run its geometry macro first).

## Outputs

`/OutputFormat 0` (default) writes the original ASCII step file `<dir>/<file>.dat`, format in
`README.txt`. Active volumes: bars, stack YSO, Pb absorber, camera layers. `copyNb1` is the mother
copy number: station 0-2 for bars (`copyNb` = bar 0-8), 10 for the camera (`copyNb` = CH index),
-1 for volumes placed in the world. `/TrackType` and `/TrackEdep` filter this file only.

`/OutputFormat 1` (or 2 for both) writes `<dir>/<file>.root` with two trees (MeV, ns, cm):
- `events`, one row per event (`/EventFilter 1` drops events that only crossed air, `/EventFilter 2`
  keeps only events with energy in a TOF bar or a Compton-detector layer):
  primary `pdg ekin costh phi x0 y0 z0`; where it ended `end_x end_y end_z end_t end_ke end_proc end_vol`
  (`end_vol` = OutOfWorld if it left; a mu- stop ends with `muMinusCaptureAtRest`);
  camera `e_ch0..2`, first-deposit times `t_ch0..2` (-1 if none), `cam_prim` (bit k: primary
  deposited in CHk) and every camera step `cam_layer cam_pdg cam_trk cam_e cam_t cam_x cam_y cam_z`;
  bars with energy `bar_st bar_id bar_e bar_t bar_x bar_y bar_z` (energy-weighted positions);
  `e_stack0 e_stack1 e_absorber`.
- `run`, one row per run: `n_events generator seed`, and for cosmic runs `flux` (/cm2/s through the
  plane), `area` (cm2), `live_time` (s), `e_min e_max theta_max plane_z model particles`, `e_mono`
  (MeV, 0 unless mono-energetic); `lead_cm`, the lead-layer thickness (0 if not built).

Energies are true deposits: apply the measured resolution (camera ~4% at 662 keV as sigma/E,
~5.5% for MIPs, ~2 mm position) in the analysis.

## Cosmic-ray generator

`/GeneralParticleSource 2`: particles start on a horizontal rectangle (`/tof/cosmic/planeCenter`,
`halfX`, `halfY`) and go down, with (E, theta) following the flux through that plane. Every event
has weight 1; the energy range, `thetaMax` and the particle selection change the normalisation,
not the weights. Live time = n_events / (flux x area), printed at the end of the run and stored in
`run`. The plane must be larger than the setup by about its height x tan(thetaMax) and lie inside
the world (checked at run start). Two spectrum models, `/tof/cosmic/model`:

- **`expacs` (default)**: energy spectra of every species in an EXPACS table, PARMA model
  (T. Sato, PLoS ONE 10(12): e0144679, 2015): neutron, proton, alpha, mu+, mu-, e-, e+, gamma
  and ions Li7 ... Ni58 (per nucleon). `data/expacs_spectrum.txt` comes from
  `claude_TOF_project/EXPACS303-eng.xlsx` (EXPACS 3.03): sea level (1047 g/cm2), latitude 35 deg,
  longitude 142 deg (cut-off rigidity 11.5 GV), 2016-07-20 (W = 46.5), ground with water fraction
  0.2. **This is not Shanghai**: set the conditions in Excel, save, and run
  `python3 scripts/expacs_to_table.py <workbook.xlsx>` (needs openpyxl) to regenerate the table;
  the conditions are in its header and printed at every run start. `/tof/cosmic/species` selects
  particles (default `mu+ mu-`, `all` for everything).
  EXPACS 3.x gives omnidirectional fluxes only, so the zenith law is a model
  (`/tof/cosmic/angular`, default `all 2 mu+ guan mu- guan`): muons follow the energy-dependent
  shape of Guan et al., other species I ~ cos^2(theta). The whole omnidirectional flux is taken as
  downward-going (no albedo from below). Spectra are interpolated as power laws between the table
  points. Sampled energy and zenith distributions agree with the table for every species
  (chi2/ndf ~ 1, 6x10^5 particles).
- **`guan`**: muons only, Guan et al. (arXiv:1509.06176, Gaisser's formula extended to low
  energies), mu+/mu- = `/tof/cosmic/chargeRatio` (1.27), `/tof/cosmic/charge both|mu-|mu+`.

**Mono-energetic mode** (`/tof/cosmic/monoEnergy E`, 0 = off): every particle gets the kinetic
energy E (per nucleon for ions), with the zenith distribution of that energy and the species
mixture of the band. The normalisation is still the spectrum integrated over `energyMin` to
`energyMax`, so set those to the energy band E stands for: the live time then turns counts into
the rate due to the particles of that band (`lead_energy_scan_e.mac` does this).

`muon_spectrum_check.mac` (10^6 muons, zenith < 75 deg): the sampled spectrum follows the
EXPACS table within statistics from 1 MeV to 1 TeV. EXPACS / Guan / PDG: vertical intensity above
1 GeV/c 61.6 / 60.4 / ~70 m^-2 s^-1 sr^-1, mean kinetic energy near the vertical 4.5 / 4.4 / ~4 GeV
(5.8 / 6.0 GeV over all zenith angles: inclined muons are harder), mu+/mu- above 1 GeV/c 1.17 / 1.27
(set) / ~1.27. The ~12% deficit against the PDG vertical intensity is shared by both models
(EXPACS site: cut-off 11.5 GV).

Both give 0.80 muons /cm2/min through a horizontal plane (EXPACS 0.800, Guan 0.805) but differ
below ~1 GeV: EXPACS has 0.34x the Guan flux under 100 MeV and 0.87x at 0.1-1 GeV, so it gives
fewer stopping muons (example layout: 34 instead of 56 mu- stops per hour in the target, 2800
instead of 5200 in the lead layer). All species at these conditions: 5.6 particles /cm2/min
through the plane, 77% photons (mostly 0.01-1 MeV), 14% muons, 5% neutrons, 3.4% e-/e+.
For neutrons below ~100 eV raise the physics list's 10 us neutron time cut
(`/physics_engine/neutron/timeLimit 10 ms`, as in `cosmic_expacs.mac`); low-energy electrons
started far above the setup lose energy in the air on the way, so keep the plane close.

## Cosmic-muon acceptance (`acceptance.mac`)

Runs on the current geometry (put a geometry macro first for other setups) and prints a table via
`analysis/acceptance.C`. Run 1: 3x10^6 realistic cosmic muons (~10 min live time) give the rate
of each selection, including muons lost in the lead layer. Run 2: 2x10^7 straight 10 GeV mu- with
an isotropic (cosine-law) flux up to 80 deg, scattering, radiative processes and secondaries off,
give the geometric acceptance A*Omega = N_sel / N_gen x A_plane x pi sin^2(80 deg). A TOF is hit
when one bar has > 0.2 MeV, the Compton detector when a stack YSO or camera layer has > 0.1 MeV
(arguments of `acceptance.C`). The macro enlarges the world to 13 x 13 x 4 m. Default geometry,
EXPACS muons (Guan muons give rates within 1-4%):

| Selection | Rate /min | A*Omega (cm2 sr) |
|---|---|---|
| top TOF | 1173 +- 11 | 4486 +- 31 |
| top & middle | 493 +- 7 | 1233 +- 16 |
| top & middle & bottom | 176 +- 4 | 376 +- 9 |
| top & middle & Compton | 120 +- 3 | 255 +- 8 |
| top & middle & Compton & bottom | 103 +- 3 | 219 +- 7 |
| top & middle & not bottom | 317 +- 6 | 857 +- 14 |

## Lead thickness x muon energy (`lead_energy_scan.mac`)

Top & middle & Compton (TMC) for lead layers of 0-20 cm (1 cm steps) and 20 mono-energetic muon
energies, 0.1-100 GeV on a log grid: 2x10^5 EXPACS mu+/mu- per point with the zenith law of their
energy, from a 2 x 2 m plane at z = 42 cm (zenith < 75 deg). Each point is normalised to the
EXPACS flux of its band, E / 1.2 to E x 1.2, so its live time gives the TMC rate due to the muons
of that band; the 20 bands tile 0.083-120 GeV. `sh scripts/lead_energy_scan.sh` runs it (8 jobs,
~1.7 h for the points plus ~6 CPU min per spectrum check run) and `analysis/lead_energy_scan.C`
writes the table (`output/leadscan/lead_energy_scan.csv`, every point: counts, ratio, A_eff,
rate) and maps (`lead_energy_scan.pdf`). Thresholds as in `acceptance.mac`.
A scan is defined by a small macro (`lead_energy_scan_def.mac`, `lead_lowE_scan_def.mac`: energy
list, band rule, output directory, seeds, check-run range); `lead_energy_scan_t.mac` runs one
thickness of any of them, `scripts/lead_energy_scan.sh [definition] [jobs] [muons] [check muons]`
runs all thicknesses in parallel, and `lead_energy_scan.C` reads the grid from the files of a
directory.

Flux calibration (page 3, `lead_energy_scan_expected.png`): each point gets the coefficient
w = J_band / J_max, its EXPACS band flux relative to the band with the most muons (w = 1 at
3.79 GeV; 0.06 at 0.1 GeV, 0.014 at 100 GeV). The band flux is the right weight because the grid
is logarithmic (a point stands for a band 0.36 E wide); dJ/dE per GeV, which peaks at 0.3 GeV,
would over-weight the low energies. ratio x w is the flux-weighted ratio, and ratio x w x N_ref,
N_ref = J_max x plane area x exposure (2.16x10^5 muons per hour), the expected number of TMC muons:
~7500 per hour without lead, 7000 with 20 cm. The 4th argument of `lead_energy_scan.C` sets the
exposure in hours (default 1); the CSV gives `flux_coeff`, `ratio_TMC_weighted` and
`N_TMC_per_hour`.

| lead (cm) | 0 | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| TMC /min, sum of the 20 points | 124.8 | 129.5 | 126.7 | 125.8 | 125.3 | 124.5 | 121.5 | 120.4 | 120.9 | 120.0 | 116.9 |
| TMC /min, spectrum check run | 123.2 | 126.5 | 128.4 | 126.6 | 128.9 | 124.2 | 122.7 | 124.0 | 121.0 | 121.2 | 118.1 |

Statistical errors are 1.3 and 2.0 /min. The sum over the points agrees with a run of the
spectrum itself over 0.083-120 GeV at every thickness (ratio 0.97-1.03, pulls within 1.5 sigma),
and the band fluxes add up to the spectrum flux exactly. The full spectrum (0-1000 GeV) gives
121.9 +- 2.0 /min at 10 cm, so the >120 GeV muons add only ~1%.
- The TMC rate hardly depends on the lead: flat at ~125 /min up to ~11 cm, -6% at 20 cm. 77% of
  it comes from 0.5-9.4 GeV muons (8 bands), which no thickness up to 20 cm stops.
- Per muon, the ratio is (4.0-4.5) x 10^-3 of the plane (A_eff ~ 170 cm2) for 0.2-2 GeV and
  ~3.0 x 10^-3 above 10 GeV without lead (flatter zenith distribution). A layer stops the muons
  below ~0.10 GeV (3-4 cm), 0.14 (6-7), 0.21 (10-12) and 0.30 GeV (17-19 cm); 0.43 GeV passes 20 cm.
- With lead, high-energy muons gain TMC from secondaries: at 33.6 GeV a third of the TMC events
  (half at 100 GeV, 20 cm) have a muon line that misses both YSO layers, against 17% (29%)
  without lead. These coincidences are real triggers but not clean muon tracks.
- TMC events whose muon stops in the Compton region (stack, Al3 and the space above it) come
  from 0.10 GeV at 3 cm, 0.14 GeV at 6 cm, 0.21 GeV at 10-11 cm and 0.30 GeV at 16-18 cm. The
  mono-energetic grid cannot give their rate (the stopping window, ~10 MeV for one zenith angle,
  is much narrower than a band); a spectrum run over ~0.05-1 GeV per thickness would.

### Fine low-energy scan (`lead_lowE_scan.mac`)

100-500 MeV in 20 MeV steps (each point normalised to its band E -+ 10 MeV), 0-20 cm of lead,
2x10^5 muons per point, check runs of 2x10^6 muons over 90-510 MeV (74 min on 8 jobs). Output in
`output/leadscan_lowE/`; `root -l -b -q 'analysis/lead_energy_scan.C("output/leadscan_lowE")'`.
- Range-out: the TMC ratio halves at 3.3 cm of lead for 100 MeV, then ~1.3 cm more per 20 MeV
  (9.8 cm at 200 MeV, 16.8 cm at 300 MeV, 19.1 cm at 340 MeV); 360 MeV and above pass 20 cm.
- TMC rate of the 90-510 MeV muons: 14.4 /min without lead, 10.0 at 10 cm, 4.8 at 20 cm; the sum
  of the points agrees with the check runs (ratio 0.95-1.03).
- TMC with the muon stopped in the Compton region (check runs): ~15 +- 1.2 /h for every thickness
  from 3 to 20 cm. More lead only moves the incident energy of the muons that stop there (the flux
  per MeV rises, the stopping power at the incident energy falls, and the two cancel). Below 3 cm
  the stopping muons have less than 90 MeV, outside this scan, so 0-2 cm are not covered. With 20
  MeV steps the sum of the points now follows the check runs (pulls mostly within 1.5 sigma).

## Validation: camera alone vs. the paper (`camera_validation.mac`, 15.5 h)

| | EXPACS muons | Guan muons | Paper (Sec. 4.1, Table 2) |
|---|---|---|---|
| CH1 x CH2 trigger rate (427.6 keV thresholds) | 10.62 /min | 10.56 /min | ~10 /min |
| three-layer muons (CH0 > 330 keV) per 15 h | 4409 | 4377 | ~4000 |
| MIP MPV, angle-corrected: CH0 / CH1 / CH2 (MeV) | 5.20 / 1.69 / 1.68 | 5.22 / 1.69 / 1.68 | 4.86 / 1.63 / 1.56 |

Rates agree to 6-9%. The simulated MPVs match the Bethe-Landau expectation for few-GeV muons in
these crystals; the measured ones are 4-8% lower, which detector effects not simulated here
(scintillator non-proportionality relative to the 662 keV calibration, SiPM saturation in CH0)
can account for.

## Known limitations

- Muonic X-rays come from Geant4's `G4EmCaptureCascade`: 2p-1s energies are 1.4-3.9% low
  (Al 334.7, Fe 1233.7, Cu 1490.8, Pb 5673.9 keV) and the Pb 2p doublet is missing. Fine for
  rates and geometry studies; inject tabulated lines for spectral identification.
- Sequential only (the framework keeps global state); run several jobs with different
  `/gun/seed` values in parallel instead.
- No intrinsic 176Lu activity, no optical photons, no digitisation yet.
