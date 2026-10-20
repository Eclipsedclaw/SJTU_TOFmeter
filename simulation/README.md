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
(`macro/vis*.mac`). The GUI shows the TOFmeter logo as its Dock and window icon (the marks in the
repository's `data/` folder, built into the executable; the reversed mark on a dark desktop),
when CMake finds the Qt SVG module.

| Macro | What it does |
|---|---|
| `smoke.mac` | 2000 vertical 1 GeV mu- over the bar stations; ASCII + ROOT |
| `cosmic.mac` | 10^6 cosmic muons through the setup as written (~14 min live time) |
| `cosmic_expacs.mac` | 10^6 cosmic rays of all EXPACS species (neutrons, photons, e+-, ...) |
| `camera_validation.mac` | camera alone looking up, 2.7x10^6 muons (~15.5 h); compare with `analysis/camera_validation.C` |
| `stopping_example.mac` | mu- below 0.5 GeV through the example layout; summarise with `analysis/stops_summary.C` |
| `vis_geometry.mac` | geometry check, no particles: position table, volume tree, 5 PNG views, VRML file |
| `geom_target_al.mac` | the default setup plus the muonic X-ray target: Al panel 20 x 20 x 1 cm centred above the Si/YSO stack at z = 35.5 to 36.5 cm, 3.9 cm above the stack and 1 cm below the Al2 plate (the free space between them is 5.9 cm); muons stopping in it end in `Target` (`end_vol`) |
| `vis_target.mac` | geometry check of `geom_target_al.mac`, no particles: position table and 4 PNG views (whole setup, target region, two edge-on views with labels) in `output/geometry_target/` |
| `vis_qt.mac` | in the Qt GUI: geometry and 10 cosmic muons; the view turns slowly until you click, scroll or press a key in it (`/tof/vis/spin`) |
| `acceptance.mac` | cosmic-muon rates and geometric acceptance A*Omega per layer and coincidence (~5 min) |
| `muon_spectrum_check.mac` | generated muon energy spectrum vs. the EXPACS table, the Guan formula and the scan's band fluxes; integral numbers vs. PDG (`analysis/muon_spectrum.C`, ~1 min) |
| `lead_energy_scan.mac` | top & middle & Compton versus muon energy (20 points, 0.1-100 GeV) and top-lead thickness (0-20 cm); in parallel with `scripts/lead_energy_scan.sh` |
| `trident_validation.mac` | muon tridents (mu -> mu mu+ mu-) in lead against Russell et al., PRL 26 (1971) 46: Geant4 cross section versus energy, biased 10.5 GeV/c beam on 114.5 g/cm2 of lead, `analysis/trident_validation.C`; run with `--preinit macro/physics_tridents.mac` (~10 s) |
| `vis_tridents_qt.mac`, `vis_tridents.mac` | only muon trident events, kept and drawn (Qt GUI: one more trident each time the Run button is pressed, the last one staying in the viewer, which turns slowly until the first click as in `vis_qt.mac`; batch: PNG files in `output/trident_views/`, each kept event drawn twice by `vis_trident_event.mac`); alias `tridentKeep` (1 every trident, 2 detectable), batch also `tridentMax` (20); with `--preinit macro/physics_tridents.mac` |
| `trident_cosmic.mac` | rate of muon tridents in the TOF-meter, biased x10^5: `sh scripts/trident_cosmic.sh [jobs] [muons per job] [lead cm]`, then `analysis/trident_rate.C` |
| `trident_scan_def.mac` | the same rate with F = 10^3 (10^5 gives rates ~5-6% low, see below) for lead 10-50 cm, 80 x 10^7 muons per thickness (519 h of exposure; 9-12 h per thickness on 10 cores, more when the computer is busy or asleep), source plane at z = 142 cm: `mkdir -p output/tridents_F1e3; nohup nice -n 15 caffeinate -is sh scripts/trident_scan.sh macro/trident_scan_def.mac 10 >> output/tridents_F1e3/scan.log 2>&1 &`; `trident_rate.C` runs on `output/tridents_F1e3/t<cm>/` when that thickness is done, `summary.txt` at the end; `touch output/tridents_F1e3/STOP` stops after the running processes, starting again (without STOP) goes on; `sh scripts/trident_remote.sh "40 50" [jobs]` in an unpacked copy of this folder builds and runs some thicknesses on another computer (same seeds), to be copied back with `scp -r`; `analysis/trident_thickness.C` then plots the detectable rate per year versus thickness (`trident_thickness.png/.pdf/.csv`) |
| `lead_lowE_scan.mac` | the same for 100-500 MeV in 20 MeV steps; in parallel with `sh scripts/lead_energy_scan.sh macro/lead_lowE_scan_def.mac` |
| `lead_thick_scan.mac` | the same for 60 energies (0.01-100 GeV) and lead 0-50 cm, 10^6 muons per point, source plane at z = 142 cm, with the muon energy after the lead (`pb_ke top_ke cmp_ke`): `caffeinate -is sh scripts/lead_energy_scan.sh macro/lead_thick_scan_def.mac 10` (about 2 days on 10 cores), then `analysis/lead_energy_scan.C` and `analysis/lead_after_spectrum.C` (run by the script) |
| `lead_filter.mac` | sea-level muon spectrum before and after a plain lead slab (8 x 8 m, in vacuum, nothing else), one thickness per process: `sh scripts/lead_filter.sh [jobs] [muons] [thicknesses in cm]` (default 6 jobs, 2x10^6 muons, 5-50 cm, ~25 min), then `analysis/lead_filter_spectrum.C`: spectra before and after, the probability to get through versus energy, and the flux of slow muons after the lead. `--below E` simulates only the muons that can leave the lead below E GeV (many more statistics there; `--below 2 --parts 4 10 15000000` gives 10x the default below 2 GeV in ~1 h, ~15 GB), `--parts n` splits each thickness into n processes |
| `vis.mac`, `vis_camera.mac`, `vis_example.mac` | PNG pictures (oblique and side view) with 30-40 muons |
| `lead_energy_scan_def.mac`, `lead_lowE_scan_def.mac`, `lead_thick_scan_def.mac` | definitions of the three lead x energy scans (output directory, energies, thicknesses, muons per point, source plane, event filter), read by `scripts/lead_energy_scan.sh` and executed by every process; the scans are built from `lead_energy_scan_t.mac` (one thickness), `lead_energy_scan_points_t.mac` (its energy points), `lead_energy_scan_e.mac` / `lead_energy_scan_e_lin.mac` (one point, log / equal bands), `lead_energy_scan_setup.mac` (geometry and source), `lead_spectrum_t.mac` (the spectrum check run) and `lead_energy_scan_check_t.mac` (only the check run, `--check-only`) |
| `geom_camera_only.mac`, `geom_tofmeter_example.mac` | geometry setups used by the macros above |
| `test.mac` | the example of the original framework (antiprotons, ASCII) |
| `README.mac` | command reference |

The ROOT macros in `analysis/` are copied to the build directory as well, e.g.
`root -l -b -q 'analysis/camera_validation.C("output/camera_validation.root")'`:

| Analysis macro | Input, what it gives |
|---|---|
| `acceptance.C` | `acceptance.mac` (runs it): rates and A*Omega per selection |
| `camera_validation.C` | `camera_validation.mac`: trigger rate, three-layer muons and MIP peaks against the paper |
| `stops_summary.C` | `stopping_example.mac`: where the primary stopped, per hour, and prompt deposits near the 2p-1s lines |
| `muon_spectrum.C` | `muon_spectrum_check.mac` (runs it): generated spectrum against EXPACS, Guan and PDG |
| `lead_energy_scan.C` | a lead x energy scan (the script runs it): coincidence rates per point and per thickness; `muonOnly` for the muon itself |
| `lead_after_spectrum.C` | the thick-lead scan (the script runs it): muon energy after the lead, before the top TOF and at the Compton detector, versus thickness |
| `lead_steps.C` | the thick-lead scan: the same muons step by step (sea level, top TOF, Compton detector) for 0-50 cm, as flux |
| `compton_arrival.C` | the check runs of the 0-20 cm log and fine scans: energy of the muons reaching the Compton detector |
| `lead_scan_status.C` | `scripts/lead_scan_status.sh` (runs it): chart of finished / running / pending points |
| `lead_filter_spectrum.C` | `scripts/lead_filter.sh` (runs it): spectra before and after a plain lead slab |
| `trident_validation.C` | `trident_validation.mac` (runs it): Geant4 tridents against Russell et al. |
| `trident_rate.C` | one folder of trident runs (the scripts run it): rates per day, where, and how many bars |
| `trident_thickness.C` | a trident scan (`output/tridents_F1e3`): detectable tridents per year versus lead thickness |

Files added to `macro/`, `analysis/`, `scripts/` or `data/` need a re-run of `cmake` to be copied.
`scripts/lead_energy_scan.sh` marks every finished process with `<dir>/log/<name>.done` and skips
those when started again, so an interrupted scan goes on where it stopped (delete the `.done`
files to redo it); a macro edited in `macro/` reaches the build directory with `make` or `cmake .`.
`sh scripts/lead_scan_stop.sh [definition]` stops a scan after the processes now running (none
is started any more; those running finish all their points), and
`sh scripts/lead_scan_resume.sh [definition] [jobs]` goes on from there in the background, first
waiting for processes still finishing (default `macro/lead_thick_scan_def.mac`, 10 jobs).
`sh scripts/lead_scan_status.sh [definition]` shows which points are finished (1), running (R) or
pending (0), as a chart per thickness and energy, `<dir>/scan_status.csv` and `<dir>/scan_status.png`.

`./TOFmeter --preinit <file.mac> [macro]` runs `file.mac` before the kernel is initialised, for the
commands Geant4 accepts only then. `macro/physics_tridents.mac` switches on muon tridents
(`muToMuonPairProd`, muons above 0.85 GeV) and photon conversion to mu+mu- (`GammaToMuPair`); both
are off otherwise. `/tof/physics/crossSection <process> <particle> <material> <E1 E2 ... GeV>` prints
a process's Geant4 cross section after a `/run/beamOn` (0 events is enough). Rare processes can be
made more frequent with `/process/em/setBiasingFactor <process> <F> true` and `/run/physicsModified`;
the affected tracks then carry weight 1/F (`tri_w`). Keep F times the probability per track well
below 1: a second biased interaction of a track already at 1/F gets 1/F^2, so rates come out low by
about half that probability (cosmic tridents with 10-20 cm of lead: 5.2-6.1% of the primary tridents
at F = 10^5, 20-30% above 100 GeV; ~0.05% at F = 10^3).

## Geometry

Default = `DetectorConstruction.cc`, z up (muons come from +z), cm. The setup stands on the ground
at z = 0, the bottom of Al4, of the corner beams and of the lead walls; x = y = 0 on its vertical
axis. From the top:

| z range | Part | Switch |
|---|---|---|
| 91 to 101 | lead layer, Pb 90x90x10 out to the lead walls, so the top is fully covered (`AbsorberPhys`); `/tof/det/leadThickness` changes the 10 cm, the bottom stays at 91 | `/tof/det/absorber` |
| 90 to 91 | stainless-steel plate 80x80x1 | `/tof/det/frame` |
| 69.8 to 70.4 | top TOF: station 0, 9 bars along x | `/tof/det/stations` |
| 67.5 to 67.8 | plate Al1 under it | `/tof/det/frame` |
| 39.8 to 40.4 | middle TOF: station 1, bars along y | |
| 37.5 to 37.8 | plate Al2 | |
| 17.8 to 31.6 | Compton detector: Si/YSO stack, 20x20 cm, Si 6 mm x3 and YSO 3 mm x2 alternating | `/tof/det/stack` |
| 17.5 to 17.8 | plate Al3 under the stack | |
| 2.3 to 2.9 | bottom TOF: station 2, bars along x | |
| 0 to 0.3 | plate Al4 | |
| 0 to 91 | four stainless L beams (10 cm legs, 1 cm thick) on the outside of the steel-plate corners | `/tof/det/frame` |
| 0 to 101 | four 10 cm lead walls, inner faces at +-45 cm (5 cm from the steel-plate edges), tops level with the lead layer | `/tof/det/leadWalls` |

Until 2026-10-20 z = 0 was at the top TOF, so every z was 70 cm lower (lead layer 21 to 31, Al4 at
-70 cm). Files made before have no run column `setup_bottom_z`; the analysis macros that use
positions (`lead_energy_scan.C`, `trident_rate.C`, `trident_thickness.C`) then take -70 cm, so old
and new files can be analysed alike. Macros of stand-alone setups (camera alone, plain lead slab,
the trident validation) keep their own origin.

Bars are 40 x 4 x 0.6 cm plastic (G4_PLASTIC_SC_VINYLTOLUENE), 2 cm above their plate; the four
plates are 80 x 62 x 0.3 cm G4_Al between the corner L beams (long side = distance between the
beams' inner faces, short side = gap between the leg tips; Al2 turned by 90 deg: 62 x 80). Changes relative to the original file: z is mirrored (that file
had z pointing down, with the lead at negative z) and starts at the ground; the world is air and 4 x 4 x 3 m (was vacuum,
1 x 1 x 1.6 m) so that the cosmic source plane fits; it is centred at z = 0, so it also holds 1.5 m of
air under the setup (no floor is modelled), and the source plane can be at most at z = 150 cm;
the 9 bars of a station sit in an air envelope, so every bar has its own ID; the world has copy
number -1. The Si-YSO face-to-face gap is 2.85 cm (the code's pitch adds only the YSO thickness).

Optional parts (off by default, set with macro commands then `/update`):
- **Compton camera** (`/tof/det/camera true`): CH2 and CH1 YSO 3 mm, CH0 LYSO 6 mm (7.2 g/cm3),
  all 47.04 mm square, 25 / 35 mm centre to centre. `/tof/camera/position` is the centre of CH2;
  `/tof/camera/rotateX` turns it (0 = looking up). Pixels, wrapping, electronics and the Al cover
  are not modelled.
- **Target** (`/tof/target/...`, position = centre) and **degrader** (`/tof/degrader/...`): boxes
  with NIST or custom (YSO, LYSO, CdZnTe) materials. `macro/geom_target_al.mac` puts the muonic X-ray
  target, 20 x 20 x 1 cm Al, at z = 35.5 to 36.5 cm above the Compton detector (the target volume is
  `Target`, so `end_vol == "Target"` marks muons stopped in it).
- **Blocks** (`/tof/blocks/add name material x y z dx dy dz cm`): any passive box, e.g. lead bricks.

`macro/geom_tofmeter_example.mac` shows them in an illustrative (not optimised) layout: a
20x20x2 cm Fe target at z = 33.5 to 35.5 cm under the middle TOF, and the paper's camera in place
of the Si/YSO stack (CH2 at z = 30 cm). The lead layer on top acts as the degrader.
Check any setup with `macro/vis_geometry.mac` (run its geometry macro first).

## Outputs

`/OutputFormat 0` (default) writes the original ASCII step file `<dir>/<file>.dat`: one line per step
in an active volume (bars, stack YSO, Pb absorber, camera layers), tab-separated: event, track ID,
parent ID, PDG code, `copyNb1`, `copyNb`, material, time (ns), kinetic energy at the start of the
step (MeV), energy deposit (MeV), creator process (`N/A` for the primary), process of the step, momentum
px py pz after the step (MeV), step length (cm), position x y z after the step (cm); `/OutputType 1`
writes 0 for both processes. `copyNb1` is the mother
copy number: station 0-2 for bars (`copyNb` = bar 0-8), 10 for the camera (`copyNb` = CH index),
-1 for volumes placed in the world. `/TrackType` and `/TrackEdep` filter this file only.

`/OutputFormat 1` (or 2 for both) writes `<dir>/<file>.root` with two trees (MeV, ns, cm):
- `events`, one row per event (`/EventFilter 1` drops events that only crossed air, `/EventFilter 2`
  keeps only events with energy in a TOF bar or a Compton-detector layer, `/EventFilter 4` only
  those with energy in both a top-TOF and a middle-TOF bar):
  primary `pdg ekin costh phi x0 y0 z0`; where it ended `end_x end_y end_z end_t end_ke end_proc end_vol`
  (`end_vol` = OutOfWorld if it left; a mu- stop ends with `muMinusCaptureAtRest`);
  camera `e_ch0..2`, first-deposit times `t_ch0..2` (-1 if none), `cam_prim` (bit k: primary
  deposited in CHk) and every camera step `cam_layer cam_pdg cam_trk cam_e cam_t cam_x cam_y cam_z`;
  bars with energy `bar_st bar_id bar_e bar_t bar_x bar_y bar_z` (energy-weighted positions);
  `e_stack0 e_stack1 e_absorber`; `cmp_ke`, the kinetic energy of the primary where it first enters a
  layer of the Compton detector (stack Si or YSO, camera), -1 if it never does; `top_ke` the same
  where it first enters a bar of the top TOF, and `pb_ke` where it first leaves the lead layer (-1
  if it stopped in it or went round it). Muon tridents, one
  entry per interaction: `tri_x tri_y tri_z` (vertex), `tri_ein tri_eout` (muon kinetic energy
  before and after), `tri_epair` (total energy of the pair), `tri_mll` (invariant mass of the two
  like-sign muons), `tri_w` (weight), `tri_pdg`; and per muon from a trident (`tmu_proc` 49) or
  from a photon conversion to mu+mu- (15): `tmu_pdg tmu_e0 tmu_e1` (kinetic energy at creation and
  at the end), `tmu_parent` (track ID of the muon that made it), `tri_parent` (the same for the
  trident), and the TOF stations crossed by each trident muon, `tmu_st` (bit 0 top, 1 middle, 2
  bottom), with its first top and middle bar `tmu_top tmu_mid`; `prim_st prim_top prim_mid` are the
  same for the primary. `/EventFilter 3` keeps only events with one of them.
  `/tof/tridents/keepForVis 1|2` keeps events with a trident (1) or a detectable one (2: the
  primary and its pair all cross the top and middle TOF) for `/vis/reviewKeptEvents`, and stops the
  run after `/tof/tridents/maxKept` (20) events, printing one line per kept event.
  `/tof/tridents/search true` turns `/run/beamOn N`, and so the Run button of the Qt toolbar
  (`/run/beamOn 1`), into a search: muons are run until N events have been kept as `keepForVis`
  says (1 if it is 0), at most `/tof/tridents/searchLimit` (10^6) muons; `false` gives normal runs
  again. With `/vis/drawOnlyToBeKeptEvents` the viewer shows the last trident until the next one.
  This is done by `TOFmeterRunManager`, a `G4RunManager` whose `BeamOn` is otherwise unchanged.
- `run`, one row per run: `n_events generator seed`, and for cosmic runs `flux` (/cm2/s through the
  plane), `area` (cm2), `live_time` (s), `e_min e_max theta_max plane_z model particles`, `e_mono`
  (MeV, 0 unless mono-energetic); `lead_cm`, the lead-layer thickness (0 if not built);
  `event_filter`, the `/EventFilter` of the run; `setup_bottom_z`, the z of the ground (0 cm; missing
  in files made with the old origin, where it was -70 cm).

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
  and ions Li7 ... Ni58 (per nucleon). `data/expacs_spectrum.txt` comes from EXPACS 3.03: sea level (1047 g/cm2), latitude 35 deg,
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


## Known limitations

- Muonic X-rays come from Geant4's `G4EmCaptureCascade`: 2p-1s energies are 1.4-3.9% low
  (Al 334.7, Fe 1233.7, Cu 1490.8, Pb 5673.9 keV) and the Pb 2p doublet is missing. Fine for
  rates and geometry studies; inject tabulated lines for spectral identification.
- Sequential only (the framework keeps global state); run several jobs with different
  `/gun/seed` values in parallel instead.
- No intrinsic 176Lu activity, no optical photons, no digitisation yet.
- No floor, walls or building: below z = 0 there is only air (the world reaches z = -150 cm).
- `/tof/vis/spin` exists only when TOFmeter runs the Qt GUI (started without a macro argument), so
  `vis_qt.mac` and `vis_tridents_qt.mac` work only there.
