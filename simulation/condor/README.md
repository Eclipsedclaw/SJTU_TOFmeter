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
S=/lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work     # --storage: the run folders go to $S/output/
sh condor/muonic_xray_czt_submit.sh --test --detector lyso --storage $S --setup <file> --os-major 7    # 2 short jobs first
sh condor/muonic_xray_czt_check.sh $S/output/muonic_xray_lyso_condor_test
sh condor/muonic_xray_czt_submit.sh --detector lyso --storage $S --setup <file> --os-major 7 --max 30 slow 1 499 fast 1 1704
sh condor/muonic_xray_czt_check.sh $S/output/muonic_xray_lyso
sh condor/muonic_xray_czt_submit.sh --detector lyso --target G4_Si --storage $S --setup <file> --os-major 7 --max 30 slow 1 499 fast 1 1704
sh condor/muonic_xray_czt_check.sh $S/output/muonic_xray_lyso_si
```

`--storage <path>` (2026-11-04) puts the run folder under `<path>/output/` (any run, any option) in place of the build
directory's `output/`: the job macros, TOFmeter's logs, the `.done` markers, the outputs and HTCondor's files. TOFmeter
itself still runs in the build directory. The path must be absolute and an existing folder you can write to. The
check and status scripts take the run folder as it is printed, e.g. `$S/output/muonic_xray_lyso`. A LYSO job writes
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
S=/lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work
sh condor/muonic_xray_czt_submit.sh --test --fluo --storage $S --setup <file> --os-major 7
sh condor/muonic_xray_czt_check.sh $S/output/muonic_xray_czt_fluo_condor_test
sh condor/muonic_xray_czt_submit.sh --fluo --storage $S --setup <file> --os-major 7 --max 30 slow 1 499 fast 1 1704
sh condor/muonic_xray_czt_check.sh $S/output/muonic_xray_czt_fluo
sh condor/muonic_xray_czt_submit.sh --fluo --target G4_Si --storage $S --setup <file> --os-major 7 --max 30 slow 1 499 fast 1 1704
sh condor/muonic_xray_czt_check.sh $S/output/muonic_xray_czt_fluo_si
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
S=/lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work
sh condor/beam_xray_submit.sh --test --storage $S --setup <file> --os-major 7      # 1 job of 2 x 10^4 muons first
sh condor/beam_xray_check.sh $S/output/beam_xray_lyso_si_al0.01_condor_test
sh condor/beam_xray_submit.sh --storage $S --setup <file> --os-major 7 1 200       # 10^8 muons
sh condor/beam_xray_check.sh $S/output/beam_xray_lyso_si_al0.01
```

- The check script follows the rules of `muonic_xray_czt_check.sh`: exit code 0, a non-empty `beam/beam_<k>.root`,
  and the log's `Number of events processed : <muons>` without an error. It is read-only and prints the commands that
  would run failed jobs again.
- Other setups get folders of their own: `--impurity none` for pure Si, another `--fraction`, `--ekin`, `--sigma-x`
  or `--sigma-y`. Give each a different `--seed-base`.
- To analyse, sync the folder into the work folder's `output/lyso/`, then:
  `root -l -b -q 'analysis/beam_impurity_spectrum.C("output/lyso/beam_xray_lyso_si_al0.01")'`.

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
