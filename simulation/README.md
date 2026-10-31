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

On an HTCondor cluster with a shared filesystem, the scripts in `condor/` (copied to the build directory like
`scripts/`) run the CZT run as jobs: one core per process, at most 60 at a time, with a short test first; see
`condor/README.md`.

| Macro | What it does |
|---|---|
| `smoke.mac` | 2000 vertical 1 GeV mu- over the bar stations; ASCII + ROOT |
| `cosmic.mac` | 10^6 cosmic muons through the setup as written (~14 min live time) |
| `cosmic_expacs.mac` | 10^6 cosmic rays of all EXPACS species (neutrons, photons, e+-, ...) |
| `camera_validation.mac` | camera alone looking up, 2.7x10^6 muons (~15.5 h); compare with `analysis/camera_validation.C` |
| `stopping_example.mac` | mu- below 0.5 GeV through the example layout; summarise with `analysis/stops_summary.C` |
| `vis_geometry.mac` | geometry check, no particles: position table, volume tree, 5 PNG views, VRML file |
| `geom_target_al.mac` | the muonic X-ray target set explicitly (the default since 2026-10-23): Al panel 20 x 20 x 1 cm centred above the Si/YSO stack at z = 35.5 to 36.5 cm, 3.9 cm above the stack and 1 cm below the Al2 plate (the free space between them is 5.9 cm); copy it for another target; muons stopping in it end in `Target` (`end_vol`) |
| `muonic_xray.mac` | muonic X-rays of cosmic mu- stopped in the target, as an ideal detector sees them (every photon of the capture, energy at creation and when leaving the volume): `sh scripts/muonic_xray.sh [--half cm] [--band "Emin Emax"] [--geometry macro] [--dir path] [processes] [muons per process] [at a time]` (default `geom_target_al.mac`, muons of 0.1-2 GeV from a 4 x 4 m plane at z = 112 cm, 30 x 5 x 10^6 muons = 50 h of exposure, ~1.5 h on 10 cores; a 2 x 2 m plane misses ~10% of the target's stops and 15-30% of the plates'; stop with `touch <dir>/STOP`, start again with more processes to add statistics), then `analysis/muonic_xray.C` for the target and for the Al plates |
| `muonic_xray_stack.mac` | muonic X-rays as the Si/YSO stack sees them: sea-level cosmic mu+ and mu- of all energies (0-1000 GeV, 4 x 4 m plane at z = 112 cm) through the default setup with its Al target, only events with energy in the stack (`/EventFilter 6`, hits `stk_*`): `sh scripts/muonic_xray_stack.sh [--dir path] [processes] [muons per process] [at a time]` (default 385 x 10^7 muons = 3 weeks of exposure, at most 90% of the physical cores at a time (10 here), about 4 days and 12 GB; resumable like `muonic_xray.sh`, any finished part can be analysed; progress bar and time left: `sh scripts/muonic_xray_stack_status.sh [--watch 60]`) |
| `muonic_xray_target.mac` | the same with a thicker Al target (20 x 20 x d cm, top at 36.5 cm growing down, at least 1 mm above the stack: 5 cm at 31.7-36.7 cm) and only the muons below 1 GeV, which hold nearly every stop in the target (the faster muons' background comes from the 3-week run): `sh scripts/muonic_xray_target.sh [--dir path] [--emax GeV] ["thicknesses"] [processes per thickness] [muons] [at a time]` (default "2 3 4 5", 29 x 10^7 muons = 7.0 days of exposure per thickness, the thicknesses taking turns, at most 90% of the physical cores, about 12 h; one folder per thickness, `output/muonic_xray_target/t<d>cm`, in the layout of the 3-week run, so `scripts/muonic_xray_stack_status.sh --dir`, the skim and the cuts work on each) |
| `muonic_xray_czt.mac` | the 4 cm Al target (32.5-36.5 cm) with the Si/YSO stack taken out and the CZT plane (`/tof/czt/...`, 4 x 4 pixels of 25 x 25 x 10 mm) right under it, only events with energy in the CZT (`/EventFilter 7`), in two bands of cosmic muons, 0-1 GeV (the target's stops) and 1-1000 GeV (background; the old runs had no CZT): `sh scripts/muonic_xray_czt.sh [--dir path] [slow processes] [fast processes] [muons] [at a time]` (default 8 slow + 28 fast processes of 5 x 10^6 muons, about 1 day of exposure each, a ~2.7 h test on 10 cores; `output/muonic_xray_czt/slow` and `fast`, each in the layout of the stack runs); as HTCondor jobs on a cluster: `sh condor/muonic_xray_czt_submit.sh`, the same processes and outputs (`condor/README.md`) |
| `vis_target.mac` | geometry check of `geom_target_al.mac`, no particles: position table and 4 PNG views (whole setup, target region, two edge-on views with labels) in `output/geometry_target/` |
| `vis_czt.mac` | pictures of the CZT setup (the geometry of `muonic_xray_czt.mac`) for slides, no particles: the whole setup cut open at y = 1.25 cm (the CZT plane shows in the cut under the target) and the target with the CZT plane from below, every other pixel lighter, plus the position table; `sh scripts/geometry_czt.sh` runs it, finishes the pictures (`scripts/picture_finish.py`: shrunk for smooth edges, white border cut, colour legend) and draws the setup to scale from the table (`analysis/geometry_schematic.C`), all in `output/geometry_czt/` (the pictures as Geant4 drew them in `raw/`) |
| `vis_qt.mac` | in the Qt GUI: geometry and 10 cosmic muons, the Al plates as outlines so that the target and the stack show; the view turns slowly until you click, scroll or press a key in it (`/tof/vis/spin`) |
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
| `lead_target_stops.C` | the lead-slab runs (`output/lead_filter_below2GeV`): mu- stopping in the 20x20 cm Al target versus the top lead thickness (each muon leaving the slab carried through the steel plate, TOF stations and Al plates to the target with Bethe-Bloch CSDA ranges, checked against Geant4 in the lead; tracks through the top layer only, the share that would come through a wall printed), and for thicker targets the X-rays reaching the stack (Al self-absorption) |
| `lead_target_stops_check.C` | the same lead-slab runs: why those stops hardly change with the lead. The sea-level mu- spectrum per MeV and per unit of ln E, the slice of it that stops in the target for each thickness (j x dE0, CSDA against Geant4), j, S and j x S (the range spectrum), Geant4 stopping depths in the lead, the directions of the slow muons against the top opening, all relative to 10 cm with the stops of `lead_target_stops.C` |
| `trident_validation.C` | `trident_validation.mac` (runs it): Geant4 tridents against Russell et al. |
| `trident_rate.C` | one folder of trident runs (the scripts run it): rates per day, where, and how many bars |
| `trident_thickness.C` | a trident scan (`output/tridents_F1e3`): detectable tridents per year versus lead thickness |
| `muonic_xray.C` | `scripts/muonic_xray.sh` (runs it): muonic X-ray and nuclear-capture gamma spectra of mu- stopped in a volume, emitted and leaving it, per day; the lines with their transitions and yields per stop; a range and bin width zoom in, e.g. `muonic_xray.C("output/muonic_xray", "Target", "mux_", 13, 27, 100, 0.25)` for 0-100 keV; a last argument `true` gives the yield (muonic X-rays / stopped mu-, per bin) instead of the rate per day (`_yield` files) |
| `muonic_xray_stack.C` | `scripts/muonic_xray_stack.sh` (any finished part of it): the add-back spectrum of the Si/YSO stack with a detector response, per day of exposure: every hit smeared with sigma/E = 4% at 1 MeV scaling as 1/sqrt(E) (Si and YSO alike), hits above 10 keV, the hits within 100 ns of an event added up; split by what made them (target X-rays, muonic X-rays of other volumes, everything else), with S, B and S/sqrt(B) in +-2 sigma around the Al lines; no TOF selection |

| `muonic_xray_stack_check.C` | the same run: health check of its finished processes (agreement between processes against Poisson, seeds, impossible hit values, the largest hit) with plots of the hit energies, times and positions and the rate per process |
| `muonic_xray_stack_skim.C` | the same run: the stack response and add-back applied once, the 20-1000 keV candidates of events with a middle-TOF hit or a stopped muon kept with their TOF and truth columns (`xray_skim.root`, with `etrue`, the entry's energy before the response); run in parts in parallel and merge with `analysis/skim_merge.C` (or `hadd`, which keeps only the first part's conditions). Detector `czt` skims a band of `muonic_xray_czt.sh` (its `czt_*` hits, a pixel's centre as its position), with another energy resolution (sigma at 1 MeV, e.g. 6.37 keV for a FWHM/E of 1.5%); each of the 16 pixels is read out alone: an entry per pixel hit, no sum over pixels (`pix`, and `npix`, `eOther`: the pixels hit within 100 ns and their energy). For both detectors also the TOF as measured, in columns of its own: bar energies with 10% FWHM, station times with 300 ps FWHM, both bar ends read out (the position along a bar from the time difference, sigma vEff x sigma_t = 1.9 cm with vEff 15 cm/ns; across the bars the bar's centre), the trigger `trig` (top and middle station within 20 ns) and from it the track (`xTgt`, `yTgt` where it crosses the target's centre plane, default z = 34.5 cm), the speed `beta` and the muon's kinetic energy `tkin`, and the bottom station's largest bar and point (`barBotM`, `xBot`, `yBot`); the true columns stay for the older cut studies. Every condition is in the file (`conditions`, `detector` and `TParameter`s) |
| `muonic_xray_stack_cuts.C` | the skim: prompt entries of stopped muons by stop volume (truth); cut variables of target X-rays against the rest; a cut flow (middle TOF in coincidence, no bottom TOF, middle bars 2-6, slow muon > 3 MeV in the middle TOF, highest hit in the top two layers, top-to-middle time > 1.3 ns) with S, B and S/sqrt(B) at the Al lines; a fit of the counts after the cuts (no truth). Assumed TOF response: 0.2 ns, 15% |
| `muonic_xray_stack_corner.C` | a skim (default the 3-week run's; a thicker target's with the 3-week run's faster muons as the rest of its background): corner plot of the five cut variables in the Al L_alpha window for the entries with a middle-TOF hit in coincidence, target X-rays (points) against the rest (log density), each variable's shapes with the cut and the shares of S and B it keeps alone and after the other cuts (N-1); e.g. `muonic_xray_stack_corner.C("output/muonic_xray_target/t4cm/xray_skim.root", 1.75, 1.1)` for the 4 cm target with its own cuts |
| `muonic_xray_target_compare.C` | the skims of the 1 cm target (3-week run) and of `muonic_xray_target.sh`: the same cuts per thickness, with the fast-muon background (primaries above 1 GeV) of the 3-week run for all; mu- stops with a prompt entry and the share through a wall (`pbKe < 0`), S before and after the cuts, the backgrounds, S/sqrt(B) for 21 days and the days to 5 sigma at the Al lines, the cut efficiency versus the depth of the stop |
| `muonic_xray_target_cuts.C` | the same skims: why the cuts of the 1 cm target lose the stops deep in a thick target, and what each thickness gives with its own slow-muon cuts. Each cut's efficiency for the target's L_alpha entries by the depth of the stop (thickest target); the middle-TOF energy and top-to-middle time against that depth, beside the fast-muon background; a scan of the two slow-muon cuts per thickness (best L_alpha S/sqrt(B) in 21 days, cross-validated on even and odd events); the background of the 0.5-1 GeV muons per thickness, a check of the fast-muon background taken from the 1 cm run |
| `muonic_xray_target_truth.C` | the stack runs of every thickness (1 cm: 3-week run; 2-5 cm: `muonic_xray_target.sh`), MC truth before any response or cut, per day: mu- stopped in the target, its muonic X-rays at the Al lines emitted, leaving the target and leaving it unscattered, stack entries with the line's full energy (true add-back energy), the stack's truth spectrum of target X-rays, the same with the skim's response (entries in the +-2 sigma windows, how many hold the line's full energy) and S after the cuts of the two cut studies; fits of the Al L lines with the response (smooth background and Gaussian against Landau lines, Poisson likelihood). One plot per file. Each thickness's whole truth spectrum and line table: `muonic_xray.C(<run folder>, "Target", "stk_")`. Another energy resolution for the spectra with the response (sigma at 1 MeV, e.g. 6.37 keV for a CZT's FWHM/E of 1.5% at 1 MeV), its files with a tag |
| `muonic_xray_target_select.C` | one target's skims (default 4 cm with the 3-week run's faster muons), split into the target's muonic X-rays (red) and the rest (blue): corner plots of ten candidate variables (the muon's five, the stack entry's five) before and after the cuts, the most efficient cuts (greedy on S/sqrt(B) in the Al L_alpha window, cross-validated), the cut flow, the spectra before and after them, and simulated measurements fitted without truth (free or known line energy and width; smooth, known or control-region background shape) with null tests. Detector `czt` for the CZT plane's run (`muonic_xray_czt.sh`), with its resolution. Skims with the measured TOF: the TOF trigger is a fixed first cut (what it costs shows in the cut flow), the track at the target and the TOF kinetic energy are the muon's variables 1 and 3; the track's zenith angle and the bottom-TOF bar are muon variables too; single-pixel CZT skims: the pixels hit in the event, their energy, the pixel's place and its number 1-16 are the entry's (the bottom-TOF bar and the pixel number are shown in the corner plots, not cut on). Every condition (the skims' and its own) goes to `muonic_xray_select_conditions.txt`, the CSV (`condition` rows) and the foot of every plot |
| `skim_merge.C` | skims of `muonic_xray_stack_skim.C` made in parts: merges them as `hadd` does (`cand` and `info` in the given order; the skim's `TParameter`s keep the first part's value, merge mode 'f', and lose `IsConst()` if the parts differ, which `muonic_xray_target_select.C` reports) and merges their `conditions`: the lines all parts share once, the others (processes, exposure, seed) part by part, with the totals (`hadd` keeps only the first part's text) |
| `geometry_schematic.C` | a log with the table of `/tof/det/print` (default that of `vis_czt.mac`, `output/geometry_czt/geometry.log`; `scripts/geometry_czt.sh` runs it): drawings to scale for slides, PDF and PNG, one per file: a section through the setup at y = 1.25 cm with labels, the same around the target with its dimensions, the setup seen from above |

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
are off otherwise.

The atomic cascade of a stopped mu- (its muonic X-rays and Auger electrons) comes from
`TOFmeterEmCaptureCascade` by default: an (n,l)-resolved cascade from n = 14 with calculated
radiative (E1) and K/L Auger rates and an initial l distribution (2l+1) exp(0.16 l), whose rates use
the tabulated 1s level and whose emitted energies use pure hydrogen-like levels with the reduced mass
(the user's model, from a G4EmCaptureCascade of the same name, renamed because Geant4's own class has
that name). `/tof/physics/muonCascade geant4` switches back to Geant4's `G4EmCaptureCascade`;
`/tof/physics/muonCascadeZmax Z` keeps the first-principles model for elements up to Z only. For Al
it gives K_alpha 354.9 keV with 0.80 per stop (measured 346.8 keV, ~0.8; Geant4's model 334.7 keV,
0.56); for heavy elements its point-nucleus 1s level is far too deep (Pb K_alpha 14.2 MeV instead of
~5.9 MeV, Geant4 5.67 MeV; Fe ~13% high), so use `muonCascadeZmax` (e.g. 20) when the lead's photons
matter. Its rate tables take a few seconds per isotope at the first capture in it.
`/tof/physics/crossSection <process> <particle> <material> <E1 E2 ... GeV>` prints
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
| 35.5 to 36.5 | muonic X-ray target `Target`: Al 20x20x1 cm on the axis, 3.9 cm above the stack and 1 cm below Al2 (not supported in the model) | `/tof/target/enable` |
| 17.8 to 31.6 | Compton detector: Si/YSO stack, 20x20 cm, Si 6 mm x3 and YSO 3 mm x2 alternating | `/tof/det/stack` |
| 17.5 to 17.8 | plate Al3 under the stack | |
| 2.3 to 2.9 | bottom TOF: station 2, bars along x | |
| 0 to 0.3 | plate Al4 | |
| 0 to 91 | four stainless L beams (10 cm legs, 1 cm thick) on the outside of the steel-plate corners | `/tof/det/frame` |
| 0 to 101 | four 10 cm lead walls, inner faces at +-45 cm (5 cm from the steel-plate edges), tops level with the lead layer | `/tof/det/leadWalls` |

Until 2026-10-23 the target was off by default: the lead and trident scans and earlier runs have
none (the muonic X-ray run set the same panel with `geom_target_al.mac`). Until 2026-10-20 z = 0 was at the top TOF, so every z was 70 cm lower (lead layer 21 to 31, Al4 at
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

Configurable parts (set with macro commands, then `/update`):
- **Compton camera** (`/tof/det/camera true`): CH2 and CH1 YSO 3 mm, CH0 LYSO 6 mm (7.2 g/cm3),
  all 47.04 mm square, 25 / 35 mm centre to centre. `/tof/camera/position` is the centre of CH2;
  `/tof/camera/rotateX` turns it (0 = looking up). Pixels, wrapping, electronics and the Al cover
  are not modelled.
- **Target** (`/tof/target/...`, position = centre), **on by default**: the muonic X-ray target above
  (volume `Target`, so `end_vol == "Target"` marks muons stopped in it); `/tof/target/enable false`
  removes it, `material`, `size` and `position` change it (`macro/geom_target_al.mac` sets the default
  values explicitly). The stand-alone setups (camera alone, lead slab, trident validation, the
  spectrum check) switch it off. **Degrader** (`/tof/degrader/...`), off by default. Both are boxes
  with NIST or custom (YSO, LYSO, CdZnTe, CZT) materials.
- **CZT plane** (`/tof/czt/...`), off by default: `pixels` x `pixels` boxes `CZTPixel` (copy number
  ix * pixels + iy) of `pixelSize` x `pixelSize` x `thickness` CZT, side by side, centred on `position`; default the
  CZT layer of `compton_satellite_sim/compton_camera_DSSD`, 4 x 4 pixels of 25 x 25 x 10 mm of Cd0.9Zn0.1Te (material
  `CZT`: Cd 9, Zn 1, Te 10 atoms, 5.78 g/cm3), at z = 31.5 to 32.5 cm, right under a 4 cm target at 32.5 to 36.5 cm.
  It overlaps the Si/YSO stack: take that out (`/tof/det/stack false`), as `macro/muonic_xray_czt.mac` does.
- **Blocks** (`/tof/blocks/add name material x y z dx dy dz cm`): any passive box, e.g. lead bricks.

`macro/geom_tofmeter_example.mac` shows them in an illustrative (not optimised) layout: a
20x20x2 cm Fe target at z = 33.5 to 35.5 cm under the middle TOF, and the paper's camera in place
of the Si/YSO stack (CH2 at z = 30 cm). The lead layer on top acts as the degrader.
Check any setup with `macro/vis_geometry.mac` (run its geometry macro first).
Pictures of the CZT setup for slides (3D views and drawings to scale): `sh scripts/geometry_czt.sh`, see `vis_czt.mac` above.

## Outputs

`/OutputFormat 0` (default) writes the original ASCII step file `<dir>/<file>.dat`: one line per step
in an active volume (bars, stack Si and YSO, Pb absorber, camera layers), tab-separated: event, track ID,
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
  `e_stack0 e_stack1` (the YSO layers of the Si/YSO stack, from the bottom), `e_si0 e_si1 e_si2` (its Si
  layers, from the bottom), `e_absorber`; the hits of the stack, one per layer and time (the deposits in a
  layer within `/tof/stack/hitWindow`, 100 ns, of the first one): `stk_layer` (0, 2, 4 Si and 1, 3 YSO, from
  the bottom), `stk_src` (what made it, as bits of the tracks that deposited, each with its secondaries:
  1 muonic X-ray of the target, 2 muonic X-ray of another volume, 4 nuclear-capture gamma, 8 anything else),
  `stk_e`, `stk_t` (first deposit) and the energy-weighted `stk_x stk_y stk_z`; `/EventFilter 6` keeps only
  events with energy in the stack. `cmp_ke`, the kinetic energy of the primary where it first enters a
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
  Photons of mu- captures at rest, one entry each: `mux_e` (energy at creation), `mux_exit` (energy
  when it first left the volume it was born in, -1 if it was absorbed there), `mux_t` (creation time),
  `mux_x mux_y mux_z`, `mux_src` (0 muonic X-ray of the atomic cascade, 1 gamma of the nuclear
  capture) and `mux_parent` (track ID of the captured muon, 1 = primary). `/EventFilter 5` keeps only
  events with one of them.
  The hits of the CZT plane (`/tof/czt/...`), as the stack's, one per pixel and time: `czt_pix` (ix * pixels + iy),
  `czt_src` (the stack's bits), `czt_e`, `czt_t` and the energy-weighted `czt_x czt_y czt_z`; `/EventFilter 7` keeps
  only events with energy in the CZT plane.
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


## Muonic X-rays of the Al target, ideal detector (`muonic_xray.mac`)

`sh scripts/muonic_xray.sh 30` (first-principles cascade, 50 h of exposure, 1.3 h on 10 cores), every
photon of the capture counted, no detector. Per day of sea-level exposure:

| | Al target 20 x 20 x 1 cm | Al plates Al1-Al4 |
|---|---|---|
| mu- stopped and captured | 324 +- 12 | 4509 +- 46 |
| K_alpha 354.9 keV (0.82 / 0.80 per stop) emitted / leaving unscattered | 264 / 193 | 3624 / 3160 |
| L_alpha 65.7 keV emitted / leaving unscattered | 209 / 122 | 2951 / 2241 |
| M_alpha 23.0 keV emitted / leaving unscattered | 116 / 12 | 1692 / 438 |
| K_beta + K_gamma + K_delta emitted | 47 | 677 |
| nuclear-capture gammas (0.1-8 MeV, ~0.8 us after the stop) | 408 | 5933 |

The plates catch 14 times more mu- than the target and give the same lines. Geant4's own cascade
(15 h sample) gives K_alpha at 334.7 keV with 0.57 per stop and L_alpha 0.40 per stop instead.
Muons stopping in the target have 186-396 MeV at sea level. Plots and line lists:
`output/muonic_xray/muonic_xray_Target.*`, `muonic_xray_FrameworkAlPhys.*`.

## Known limitations

- Muonic X-ray energies are model lines without fine structure: the first-principles cascade
  (default) puts Al K_alpha at 354.9 keV (+2.3%) and heavy elements far too high (above); Geant4's
  model is 1.4-3.9% low (Al 334.7, Fe 1233.7, Cu 1490.8, Pb 5673.9 keV) and lacks the Pb 2p doublet.
  Fine for rates and geometry studies; inject tabulated lines for spectral identification.
- Sequential only (the framework keeps global state); run several jobs with different
  `/gun/seed` values in parallel instead.
- No intrinsic 176Lu activity, no optical photons, no digitisation yet.
- No floor, walls or building: below z = 0 there is only air (the world reaches z = -150 cm).
- `/tof/vis/spin` exists only when TOFmeter runs the Qt GUI (started without a macro argument), so
  `vis_qt.mac` and `vis_tridents_qt.mac` work only there.
