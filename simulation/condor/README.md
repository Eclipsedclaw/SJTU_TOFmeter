# HTCondor jobs: the CZT run on a cluster

The CZT run (`macro/muonic_xray_czt.mac`, run on the shared server by `scripts/muonic_xray_czt.sh`) as HTCondor
jobs, for a cluster with a shared filesystem such as SJTU's INPAC cluster. Each job is one process on one core
(TOFmeter runs sequentially). The macros, seeds, outputs and `.done` markers are those of the shared server's run
(each job macro also names its target material, `/tof/target/material G4_Al` by default), so processes run in either
place fit together process by process.

| File | What it does |
|---|---|
| `muonic_xray_czt_submit.sh` | Writes the macros of the processes not yet done and one HTCondor submit file, then submits it. |
| `muonic_xray_czt_job.sh` | One job: `./TOFmeter <dir>/<band>/log/job<k>.mac` (with `--preinit <file>` when the macro's first line is `# preinit <file>`, as for `--detector lyso`). It marks the process done only if the run reached its end, and writes host, time and speed to `job<k>.out`. |
| `muonic_xray_czt_compare.C` | Compares two outputs of the same process (same seed) event by event, or their distributions if the events differ. |
| `muonic_xray_czt_check.sh` | The state of each job of a folder, checked (exit code 0, output not empty, log without error): queued, running, held, finished OK, failed, the progress; for failed and held jobs the HTCondor ID, the reason and the log paths, and the commands to run them again. Read-only. |
| `muonic_xray_czt_status.sh` | Progress of a folder's jobs: bars overall, per band and per running job; exposure, failed jobs, time left. |
| `beam_xray_submit.sh`, `beam_xray_job.sh`, `beam_xray_check.sh` | (2026-11-07) The same three for the beam study of Al in Si (`macro/beam_xray_lyso.mac`): a mu- beam into the target on the LYSO panel. See "The beam study" below. |
| `beam_steel_submit.sh`, `beam_steel_job.sh`, `beam_steel_check.sh` | (2026-11-10) The same three for the steel study with U, Th and K: a mu- beam into a 1 cm steel plate on the CZT plane (band `beam`, `macro/beam_steel_lyso.mac`) and the plate's natural radioactivity (band `decay`, `macro/beam_steel_decay.mac`). See "The steel study" below. |
| `room_background_submit.sh`, `room_background_check.sh` | (2026-11-20) The steel study's room background, in the beam run's setup on a 30 cm concrete floor: the floor's natural decays (band `floor`), radon in the lead castle's air (`radon`) and sea-level cosmic rays (`cosmic`); the jobs run `beam_steel_job.sh`. See "The steel study's room background" below. |
| `muonic_library_submit.sh`, `muonic_library_check.sh` | (2026-11-23) The muonic X-ray library: every element from H to Mo as a pure plate, and two test alloys, 10^8 mu- each in the steel beam run's setup, each material at its own beam energy (muons stop 8.2424 mm deep in all); jobs 11-100 of each material `--slim`; the jobs run `beam_steel_job.sh`. See "The muonic X-ray library" below. |

CMake copies this folder's scripts (`*.sh`, `*.C`) into the build directory (`<build>/condor/`); run `cmake` again
after changing a file here.

## Before submitting

- Build TOFmeter on the cluster's shared filesystem. The jobs start in the build directory, read `./TOFmeter`,
  `macro/` and `data/` there, and write their output into it. Nothing is transferred (`should_transfer_files = NO`).
- The jobs take the environment of the shell you submit from (`getenv = True`), so set up Geant4 there.
- `condor_submit` may not start in some shells, e.g. with `error while loading shared libraries:
  libcrypto.so.3`, when the shell lacks libraries HTCondor needs. Submit from a shell where `condor_q` works; if
  the jobs need more than that shell's environment, add `--setup <file>`. Each job runs under bash and sources that
  file before TOFmeter starts. The script checks that `condor_submit` starts before it writes anything.

## Test

```sh
cd <build>
sh condor/muonic_xray_czt_submit.sh --test      # 2 jobs: slow and fast process 1, 2 x 10^5 muons each
                                                # (add --setup <file> as above if needed)
condor_q                                        # batch tof_muonic_xray_czt_condor_test
cat output/muonic_xray_czt_condor_test/*/log/job1.out
```

It works if both processes have `log/job1.done`. Each `job1.out` names the host and the Geant4 version, and gives
the speed and the time 5 x 10^6 muons would take (an upper limit: the slow band's first mu- captures take about a
minute once). `--no-submit` writes everything without submitting; `rm -r output/muonic_xray_czt_condor_test` before
testing again.

The test uses the seeds of the shared server's process 1, so its muons are that process's first 2 x 10^5. Copy the
shared server's `output/muonic_xray_czt/slow/stk_1.root` (and `fast/`) to the cluster, e.g. into `server/slow/`, then:

```sh
root -l -b -q 'condor/muonic_xray_czt_compare.C("output/muonic_xray_czt_condor_test/slow/stk_1.root", "server/slow/stk_1.root")'
```

`the same events, bit for bit` means the same physics as on the shared server. With another Geant4 version or
compiler, the events differ. The macro then compares the two files as independent samples and says whether they are
`consistent within the statistics`:

- the events written per muon (Poisson);
- the muon's end processes (chi2);
- the muon's energy and the event's CZT energy (Kolmogorov-Smirnov).

On the shared server, a 5 x 10^4-muon run matched process 1's first events (54 of 54 slow, 53 of 53 fast).

On the INPAC cluster (2026-10-31) the setup that works is:

- the CentOS 7 workers (423 slots; the 97 EL9 slots refuse these jobs);
- `LCG_105a/x86_64-centos7-gcc11-opt` (Geant4 11.2.0), built on the login node `bl-0` with `-DWITH_GEANT4_UIVIS=OFF`;
- submitted with `--setup <file sourcing that view> --os-major 7`.

Its test agreed with the shared server's process 1 (Geant4 11.2.2) within the statistics in both bands:

| Band | Events per 10^6 muons, cluster | Shared server | p (end processes) | p (KS, muon energy) | p (KS, CZT energy) |
|---|---|---|---|---|---|
| slow | 760 +- 62 | 779 +- 12 | 0.13 | 0.68 | 0.92 |
| fast | 1120 +- 75 | 1203 +- 16 | 0.50 | 0.56 | 0.51 |

The live time per muon was the same as on the shared server.

## Production

```sh
sh condor/muonic_xray_czt_submit.sh slow 79 88 fast 47 168
```

This runs the processes the shared server had not finished on 2026-10-31. The arguments are `<band> <first> <last>`,
repeated as needed. Of each range, the processes without `log/job<k>.done` are submitted; finished ones are skipped,
so the same command can be run again later.

| Option | Default | Meaning |
|---|---|---|
| `--target name` | `G4_Al` | the target's NIST material (`/tof/target/material` in each job macro), e.g. `G4_Si`; see Another target material |
| `--seed-base b` | by target | slow seeds b + k, fast b + 100000 + k: G4_Al 9500000, G4_Si 9700000; any other target needs one |
| `--dir path` | `output/muonic_xray_czt` (G4_Al), `output/muonic_xray_czt_<name>` (others) | the run's folder, laid out as on the shared server |
| `--muons n` | 5 x 10^6 | muons per process: ~35 min slow, ~50 min fast on the shared server |
| `--max n` | 60 (the most) | jobs at a time (`max_materialize`; one core each) |
| `--setup file` | | a file each job sources before TOFmeter starts (your Geant4 setup), on top of the submitting shell's environment |
| `--os-major n` | any | run only on workers of that OS major version (`TARGET.OpSysMajorVer`, e.g. 9 for EL9): the one TOFmeter was built for; `condor_status -af Machine OpSysAndVer` lists the workers |
| `--no-submit` | | write the macros and the submit file only |

Bands: `slow` is 0-1 GeV, seed 9500000 + k; `fast` is 1-1000 GeV, seed 9600000 + k (with G4_Al).

Each job asks for 1 core and 1 GB of memory (TOFmeter uses ~120 MB). A job that fails, or ends without the run's
last line `Cosmic rays: <muons> events` (e.g. a macro command not found, after which TOFmeter still exits with 0), is
put on hold. So is a job still running after 23 h (the limit is 24 h).

The script refuses to submit while jobs of the same folder are in the queue. The 60-job limit holds per submission,
and this keeps a process from running twice.

## Another target material

The same setup with a silicon target (the same 20 x 20 x 4 cm box, the same CZT plane and TOF), as large as the
aluminium run (200 slow and 684 fast processes of 5 x 10^6 muons, 24.1 days of exposure in each band; about 560
core-hours on the cluster, the aluminium jobs' mean of 0.47 h per slow and 0.68 h per fast process):

```sh
cd <build>
sh condor/muonic_xray_czt_submit.sh --test --target G4_Si --setup <file> --os-major 7    # 2 short jobs first
sh condor/muonic_xray_czt_check.sh output/muonic_xray_czt_si_condor_test
sh condor/muonic_xray_czt_submit.sh --target G4_Si --setup <file> --os-major 7 slow 1 200 fast 1 684
sh condor/muonic_xray_czt_check.sh output/muonic_xray_czt_si
```

- Each job macro sets `/tof/target/material G4_Si` before `macro/muonic_xray_czt.mac` rebuilds the geometry
  (`/update`); TOFmeter's log then says `Target: G4_Si, size (20,20,4) cm at (0,0,34.5) cm` (the `G4_Al ... 1 cm`
  line before it is the default geometry at start-up).
- The run goes into `output/muonic_xray_czt_si` (batch `tof_muonic_xray_czt_si`) with its own seeds: slow 9700000 + k,
  fast 9800000 + k, used by no other run. A folder holds the jobs of one target: the script refuses to add G4_Si jobs to
  a folder of G4_Al jobs, or the other way round.
- Tested on the shared server (2026-11-02): a slow job of 2 x 10^5 muons ran in 203 s; its one mu- captured in the
  target emitted silicon's cascade at the model's point-nucleus energies (M 26.7, L_alpha 76.2, K_alpha 411.7 keV; the
  measured K_alpha is 400.2 keV), and its CZT hits were tagged as the target's X-rays, so the skim and the analysis
  take the run as they take the aluminium one. The analysis' line energies are aluminium's so far (L_alpha 65.73 keV).

## The LYSO panel (`--detector lyso`)

The same setup and bands with the CZT plane replaced by the LYSO panel and its own Lu-176 (2026-11-03,
`macro/muonic_xray_lyso.mac`: since 2026-11-04 20 x 20 pixels of 6.25 x 6.25 x 3 mm on a 6.33 mm pitch with ESR and an Al frame, its Al face towards the target, its pixels' centre at cztZ 32.332 cm; 12.6 kBq
of Lu-176 decays within +-1 us of every cosmic-ray event). Each job runs `./TOFmeter --preinit macro/preinit_lyso.mac
<job.mac>`: the job macro's first line, `# preinit macro/preinit_lyso.mac`, tells `muonic_xray_czt_job.sh`. Two months
for each target, as the CZT run's 60 days (499 slow and 1704 fast processes of 5 x 10^6 muons). Rebuild TOFmeter on
the cluster first (`cmake` again, then `make`: there is a new source file, `src/IntrinsicSource.cc`, and new macros):

```sh
cd <build>
sh condor/muonic_xray_czt_submit.sh --test --detector lyso --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work --setup <file> --os-major 7    # 2 short jobs first
sh condor/muonic_xray_czt_check.sh /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work/output/muonic_xray_lyso_condor_test
sh condor/muonic_xray_czt_submit.sh --detector lyso --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work --setup <file> --os-major 7 --max 30 slow 1 499 fast 1 1704
sh condor/muonic_xray_czt_check.sh /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work/output/muonic_xray_lyso
sh condor/muonic_xray_czt_submit.sh --detector lyso --target G4_Si --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work --setup <file> --os-major 7 --max 30 slow 1 499 fast 1 1704
sh condor/muonic_xray_czt_check.sh /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work/output/muonic_xray_lyso_si
```

`--storage <path>` (2026-11-04) puts the run folder under `<path>/output/` (any run, any option) in place of the build
directory's `output/`: the job macros, TOFmeter's logs, the `.done` markers, the outputs and HTCondor's files. TOFmeter
itself still runs in the build directory. The path must be absolute and an existing folder you can write to. The
check and status scripts take the run folder as it is printed, e.g. `/lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work/output/muonic_xray_lyso`. A LYSO job writes
about 10-11 MB (5 x 10^6 muons; the CZT's 1.6-2.4 MB), so 60 days take about 24 GB per target.

- Folders `output/muonic_xray_lyso` (G4_Al) and `output/muonic_xray_lyso_si` (G4_Si), batches `tof_muonic_xray_lyso`
  and `tof_muonic_xray_lyso_si`; seeds slow 9900000 + k and fast 10000000 + k (Al), slow 10100000 + k and fast
  10200000 + k (Si), used by no other run. A folder holds one detector and one target: the script refuses to mix them.
  The two targets are separate folders and batches, so both may be in the queue at once; the 60-job limit holds per
  submission, so then give each `--max 30` to stay within 60 cores together.
- The test's seeds are those of production process 1, so its muons are that process's first 2 x 10^5.
- Speed (shared server, 2026-11-03): a slow job of 2 x 10^5 muons in 185 s, as the CZT run's (187 s; with the Auger
  electrons that the radioactive decay turns on it took 273 s, so the macro turns them off). On the cluster, as the CZT
  jobs, about 0.5 h per slow and 0.7 h per fast process: ~1400 core-hours per target, about a day at 60 cores.
  TOFmeter uses ~135 MB.
- The check script names the detector in its first line (`detector lyso, target G4_Al`) and, at the end, the Geant4
  version of the finished jobs (from their `job<k>.out`), with a warning when a folder mixes versions.
- Geant4 11.2.1 or later (2026-11-04): 11.2.0 (LCG_105a) samples beta decays wrongly (fixed in 11.2.1, problem
  report #2588), which makes the LYSO's Lu-176 deposits too soft; the CZT runs have no radioactive decay and are not
  affected. After a rebuild, a LYSO job's `.out` must read `with geant4-11-02-patch-01` or later; jobs of an older
  build go into a folder of their own (the submit script skips every process that has a `.done` marker).

## The CZT run with fluorescence (`--fluo`)

The CZT runs so far had fluorescence off (QGSP_BERT's default): no Cd/Te K X-rays, so no escape peaks (23.2 and 27.5 keV
below each line; 9% of the Al L_alpha photons in the CZT, 7% of Si's) and no Pb K X-rays of the shield (72.8-87.3 keV,
Pb K_alpha 75.0 keV 1.3 keV below Si's L_alpha). `--fluo` (2026-11-04) runs the same CZT setup with
`./TOFmeter --preinit macro/preinit_fluo.mac` (`/process/em/fluo true` and `/process/em/deexcitationIgnoreCut true`:
without the second, the 1 mm production cut would drop these X-rays again), into `output/muonic_xray_czt_fluo` (G4_Al)
and `output/muonic_xray_czt_fluo_si` (G4_Si), seeds slow 10300000 + k and fast 10400000 + k (Al), 10500000 + k and
10600000 + k (Si):

```sh
cd <build>
sh condor/muonic_xray_czt_submit.sh --test --fluo --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work --setup <file> --os-major 7
sh condor/muonic_xray_czt_check.sh /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work/output/muonic_xray_czt_fluo_condor_test
sh condor/muonic_xray_czt_submit.sh --fluo --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work --setup <file> --os-major 7 --max 30 slow 1 499 fast 1 1704
sh condor/muonic_xray_czt_check.sh /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work/output/muonic_xray_czt_fluo
sh condor/muonic_xray_czt_submit.sh --fluo --target G4_Si --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work --setup <file> --os-major 7 --max 30 slow 1 499 fast 1 1704
sh condor/muonic_xray_czt_check.sh /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work/output/muonic_xray_czt_fluo_si
```

The 60-core limit holds for all your jobs together: two runs at a time with `--max 30` each, or four with `--max 15`.

## The beam study (`beam_xray_*.sh`, 2026-11-07)

A mono-energetic mu- beam goes straight down into the 20 x 20 x 4 cm target lying on the LYSO panel
(`macro/beam_xray_lyso.mac`, `./TOFmeter --preinit macro/preinit_lyso.mac`).
- Defaults: PSI muE1's spot (39 x 25 mm FWHM), 29.95 MeV (85 MeV/c: every muon stops at 2.5 cm depth in Si, rms 0.9 mm),
  a G4_Si target with 1% of its atoms Al.
- `analysis/beam_impurity_spectrum.C` reweights the Al-caught and Si-caught muons (column `cap_z`) to any Al
  concentration and scales the run to one day of beam.
- 5 x 10^5 muons per job: about 25-40 min, 140 MB of output, 150 MB of memory (local test 2026-11-07: 2.7 ms per muon
  with Geant4 11.2.2).
- Folder `output/beam_xray_lyso_si_al0.01`, seeds 10700000 + k.

The code is new (`/tof/target/impurity`, the `cap_z` column, the new macro), so run `cmake` and `make` in the build
directory first. Then:

```sh
cd <build>
sh condor/beam_xray_submit.sh --test --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work --setup <file> --os-major 7      # 1 job of 2 x 10^4 muons first
sh condor/beam_xray_check.sh /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work/output/beam_xray_lyso_si_al0.01_condor_test
sh condor/beam_xray_submit.sh --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work --setup <file> --os-major 7 1 200       # 10^8 muons
sh condor/beam_xray_check.sh /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work/output/beam_xray_lyso_si_al0.01
```

- The check script follows the rules of `muonic_xray_czt_check.sh`: exit code 0, a non-empty `beam/beam_<k>.root`,
  and the log's `Number of events processed : <muons>` without an error. It is read-only and prints the commands that
  would run failed jobs again.
- Other setups get folders of their own: `--impurity none` for pure Si, another `--fraction`, `--ekin`, `--sigma-x`
  or `--sigma-y`. Give each a different `--seed-base`.
- To analyse, sync the folder into the work folder's `output/lyso/`, then:
  `root -l -b -q 'analysis/beam_impurity_spectrum.C("output/lyso/beam_xray_lyso_si_al0.01")'`.

## The steel study with U, Th and K (`beam_steel_*.sh`, 2026-11-10)

A mono-energetic mu- beam goes straight down into a 20 x 20 x 1 cm stainless-steel plate lying on the CZT plane (4 x 4
pixels of 25 x 25 x 10 mm, their top against the plate; `--detector lyso` for the LYSO panel), with
`./TOFmeter --preinit macro/preinit_lyso.mac` (radioactive decay and fluorescence; nothing of the LYSO).
- Two bands in one folder: `beam` (processes 1-N, 5 x 10^5 muons each, seeds 10900000 + k, every event written) and
  `decay` (processes 1-4, 5 x 10^6 decays each, seeds 11000000 + k: 1 U-238, 2 Th-232, 3 K-40, 4 U-235, each nucleus
  with its whole chain in the event, secular equilibrium).
- Target: `--composition Fe:rest,Cr:16,Ni:10,Mo:2,Mn:2,Si:0.75,N:0.10,C:0.08,P:0.045,S:0.03,U:0.1,Th:0.1,K:0.1` (mass
  percent; U, Th, K raised to 0.1% so that enough muons catch on them) at `--density 8.0`. The analysis reweights every
  element to the real steel (`cap_z`).
- The muonic X-rays get the measured (or interpolated) fine-structure lines of all 13 elements
  (`/tof/physics/muonXrayEnergies measured`, set by the macro).
- Speed (local test 2026-11-10, Geant4 11.2.2): a beam job spends ~3 min building the cascade models of ~40 isotopes,
  then ~1.3 ms per muon: ~15 min for 5 x 10^5 muons, 150 MB of memory, ~230 MB of output (460 B per muon: ~46 GB for
  10^8). A decay job: ~0.5 ms per chain, ~45 min for 5 x 10^6, ~280 MB of output.
- Folder `output/beam_steel_czt_trace` (with `--storage`, under `<path>/output/`).

The code is new (`/tof/target/composition`, the measured lines of the new elements, and a fix in the panel's hit
building: hits at the times of long-lived decays, ~1e25 ns, made the old `time < first + window` test loop forever), so
bring the sources to the cluster and run `cmake` and `make` in the build directory first. Then:

```sh
cd <build>
sh condor/beam_steel_submit.sh --test --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work --setup ~/tof_setup_c7.sh --os-major 7    # beam 1 x 2e4 muons, decay 4 x 2e4
sh condor/beam_steel_check.sh /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work/output/beam_steel_czt_trace_condor_test
sh condor/beam_steel_submit.sh --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work --setup ~/tof_setup_c7.sh --os-major 7 1 200    # 10^8 muons + the 4 decay jobs
sh condor/beam_steel_check.sh /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work/output/beam_steel_czt_trace
```

- The check script follows the rules of the others, per band: exit code 0, a non-empty `<band>/<file>.root` (the job
  macro's `bmFile`), and the log's `Number of events processed : <events>` without an error. It is read-only and
  prints the commands that would run failed jobs again (`--no-beam` for failed decay jobs only).
- To analyse, sync the folder into the work folder's `output/czt/`, then:
  `root -l -b -q 'analysis/beam_steel_trace_spectrum.C("output/czt/beam_steel_czt_trace")'`.
  Then the peaks, named and fitted (from the ROOT file that macro writes, in seconds):
  `root -l -b -q 'analysis/beam_steel_peak_fit.C("output/czt/beam_steel_czt_trace")'`.

## The steel study's room background (`room_background_*.sh`, 2026-11-20)

For the TOF-trigger ON/OFF subtraction at 50 ppm U, Th and K (`analysis/beam_steel_onoff.C`): the room's background in
the setup the steel beam run had (the plate on the CZT plane inside TOFmeter's default setup: a 10 cm lead roof and four
10 cm lead walls, open at the bottom), now standing on a 30 cm concrete floor (`G4_CONCRETE`, -30..0 cm), with
`--preinit macro/preinit_lyso.mac`. Three bands, one folder, a fixed table of 127 processes (`<band>/log/plan.txt`):
- `floor` (98 processes, seeds 11100000 + k, `macro/room_decay.mac`): U-238, Th-232, K-40 and U-235, each nucleus with its
  whole chain, in five regions of the floor: `inner` (x and y within +-55 cm, under the castle) and the ring around it out
  to +-100 cm (`px`, `nx`, `py`, `ny`). At the analysis's default activities (K-40 400, U-238 40, Th-232 30 Bq/kg) that is
  ~0.1 day of the inner region and ~0.02 day of the ring.
- `radon` (4 processes of 8 x 10^5, seeds 11200000 + k): Rn-222 with its short-lived daughters in the castle's air, ~1 day
  at 50 Bq/m3.
- `cosmic` (25 processes of 3.25 x 10^6, seeds 11300000 + k, `macro/room_cosmic.mac`): every EXPACS species at sea level,
  0.25 day.
- Speed (local test 2026-11-20, Geant4 11.2.2): U-238 and Th-232 chains ~2040 /s, K-40 ~25000 /s, radon ~140 /s, cosmic
  rays ~450 /s: ~1.6-2 h per job, ~210 core-hours in all. Memory ~150-570 MB per job.
- Folder `output/beam_steel_czt_room` (with `--storage`, under `<path>/output/`).

```sh
cd <build>
sh condor/room_background_submit.sh --test --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work --setup ~/tof_setup_c7.sh --os-major 7    # 6 jobs of 2e4 events
sh condor/room_background_check.sh /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work/output/beam_steel_czt_room_condor_test
sh condor/room_background_submit.sh --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work --setup ~/tof_setup_c7.sh --os-major 7    # the 127 jobs
sh condor/room_background_check.sh /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work/output/beam_steel_czt_room
```

- New macros only (no new code), so after bringing the sources over, `cmake .` in the build directory is enough.
- The check script follows the rules of the others, per band (exit code 0, a non-empty `<band>/<band>_<k>.root`, the
  log's end without an error) and lists what each job simulates with `--all`. It is read-only; the submit command it
  prints submits every process without a `.done` marker.
- To analyse, sync the folder into the work folder's `output/czt/`, then (the walls' estimate runs locally with
  `macro/room_lines.mac` into `output/czt/beam_steel_czt_room_walls/lines`):
  `root -l -b -q 'analysis/room_background_spectrum.C("output/czt/beam_steel_czt_room", "", "output/czt/beam_steel_czt_room_walls/lines")'`
  and then
  `root -l -b -q 'analysis/beam_steel_onoff.C("output/czt/beam_steel_czt_trace/analysis_50ppm/beam_steel_trace_spectrum.root", "output/czt/beam_steel_czt_room/analysis/room_background_spectrum.root")'`.

## The muonic X-ray library (`muonic_library_*.sh`, 2026-11-23)

For unmixing a compound's spectrum into its elements (the user: pure-element spectra first, "super precise", then
solve for the components' ratios): one spectrum per material of `data/muonic_library_materials.txt`, every element
from H to Mo as a pure 20 x 20 x 1 cm plate and the test alloys Ni80Cr20 and Fe50Ni50 (mass), in the steel beam run's
setup (`macro/beam_steel_lyso.mac` on the CZT plane, `--preinit macro/preinit_lyso.mac`, every event written).
- The gases as liquids (NIST's lH2, lN2, lO2, lAr, lKr, lBr; He 0.125, Ne 1.207, F 1.505, Cl 1.5625 g/cm3).
- Every material at its own kinetic energy, tuned by `scripts/muonic_library_depth.sh` (2026-11-23): the muons stop at
  the steel run's mean depth, 8.2424 mm below the plate's top, in every material (all 44 within 0.025 mm; 27-89 MeV/c).
- The X-ray lines: `/tof/physics/muonXrayEnergies measured`, since 2026-11-23 for every element up to Mo
  (`scripts/muonic_lines.py`).
- Per material (folder `<id>_<name>`, id = Z; 101, 102 the alloys) 100 processes of 10^6 muons, seeds
  12000000 + 1000 id + k; processes 1-10 in the full format, 11-100 `TOFmeter --slim` (a job macro line `# slim`).
- 4400 processes. Local check runs (Geant4 11.2.2, 10^4 muons): 41-72 s per job of 10^4 including the start-up; so
  about 20-40 min per job of 10^6 on the cluster, ~1500-3000 core-hours (1-2 days on 60 cores).
- Size: full ~340-465, slim ~125 bytes per muon, so ~16 GB per material, ~0.7 TB in all.
- Folder `output/muonic_library_czt` (with `--storage`, under `<path>/output/`).

```sh
cd <build>
sh condor/muonic_library_submit.sh --test --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work --setup ~/tof_setup_c7.sh --os-major 7    # 88 jobs of 1e4 muons: process 1 full, 2 slim
sh condor/muonic_library_check.sh /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work/output/muonic_library_czt_condor_test
sh condor/muonic_library_submit.sh --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work --setup ~/tof_setup_c7.sh --os-major 7    # the 4400 jobs
sh condor/muonic_library_check.sh /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work/output/muonic_library_czt
```

- New code (`--slim`, the measured lines of Z 1-42): rebuild TOFmeter (`make`) after bringing the sources over; the
  submit script refuses a build without them.
- `--materials Fe,Cr,...` submits only some materials; `--jobs`, `--muons` and `--full-jobs` change the plan (keep them
  for a folder: the seeds and the full/slim split follow the process number).
- The check script follows the rules of the others, per material (exit code 0, a non-empty `<id>_<name>/beam_<k>.root`,
  the log's end without an error) and also checks that each job wrote the ntuple its macro asked for (slim or full).
  It is read-only; the submit command it prints submits every process without a `.done` marker.
- The test's depths under the cluster's Geant4: sync a few test folders back and run
  `root -l -b -q 'analysis/muonic_library_depth.C("output/czt/muonic_library_czt_condor_test/26_Fe/beam_*.root")'`.

## Following the jobs

```sh
sh condor/muonic_xray_czt_check.sh <dir | cluster id> [--all]  # each job's state, checked; failed and held jobs
sh condor/muonic_xray_czt_status.sh [--dir path] [--watch 60]   # progress bars, exposure, time left
condor_q                                      # running, idle, held
condor_q -hold                                # held jobs and why; see <dir>/<band>/log/job<k>.err and .log
condor_release <id>                           # run a held job again after a fix
condor_rm <id>                                # remove; then submit again to rerun what is not done
```

`muonic_xray_czt_check.sh` checks rather than assumes. A job finished OK when HTCondor's log of its submission says
it ended with exit code 0, its `stk_<k>.root` exists and is not empty, and TOFmeter's log reached `Cosmic rays:
<muons> events` without an error (Geant4's `EEEE` errors, crashes, aborts; its warnings do not count). It prints the
states per job (as ranges of process numbers; `--all` one line each) and in total, the progress (e.g. `37/50 finished
OK`), each failed or held job with its HTCondor ID, the reason and its logs, and the commands that would run them again
(`condor_release`, or the submit command with the failed processes), which it only prints. Run it where you submitted:
`condor_q` there tells queued, running and held.

`muonic_xray_czt_status.sh` shows:

- a bar over all the folder's processes, and one per band with the days of exposure done;
- one bar per running job (from the last `>>> Event N` of its log, which is written in blocks and lags a little);
- the failed jobs (`FAILED` in `.err`) and the stopped ones (log unchanged for 30 min without `.done`: held or
  removed);
- HTCondor's counts for the batch, and the time left (each band's job speed, shared by the jobs running now).

It only reads files, so it also works on a copy of the folder.

## Outputs

The layout is that of the shared server's run:

- `<dir>/<band>/stk_<k>.root`;
- in `<dir>/<band>/log/`, for each process:
  - `job<k>.mac`, the macro;
  - `job<k>.log`, TOFmeter's log;
  - `job<k>.out` and `job<k>.err`, the job's own output and errors;
  - `job<k>.done`, the marker of a finished process;
- in `<dir>/condor/`, each submission's `submit_<date>.sub` and HTCondor's log of its jobs, `submit_<date>.log`.

To complete the shared server's run, copy `stk_<k>.root` and `log/` of the finished processes into its
`output/muonic_xray_czt/<band>/`. The skim (`analysis/muonic_xray_stack_skim.C`) and the analysis take them as they
are.
