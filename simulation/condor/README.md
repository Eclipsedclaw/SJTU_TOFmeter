# HTCondor jobs: the CZT run on a cluster

The CZT run (`macro/muonic_xray_czt.mac`, run on the shared server by `scripts/muonic_xray_czt.sh`) as HTCondor
jobs, for a cluster with a shared filesystem such as SJTU's INPAC cluster. Each job is one process on one core
(TOFmeter runs sequentially). The macros, seeds, outputs and `.done` markers are those of the shared server's run,
so processes run in either place fit together process by process.

| File | What it does |
|---|---|
| `muonic_xray_czt_submit.sh` | Writes the macros of the processes not yet done and one HTCondor submit file, then submits it. |
| `muonic_xray_czt_job.sh` | One job: `./TOFmeter <dir>/<band>/log/job<k>.mac`. It marks the process done only if the run reached its end, and writes host, time and speed to `job<k>.out`. |
| `muonic_xray_czt_compare.C` | Compares two outputs of the same process (same seed) event by event. |

CMake copies this folder into the build directory (`<build>/condor/`); run `cmake` again after changing a file here.

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
| `--dir path` | `output/muonic_xray_czt` | the run's folder, laid out as on the shared server |
| `--muons n` | 5 x 10^6 | muons per process: ~35 min slow, ~50 min fast on the shared server |
| `--max n` | 60 (the most) | jobs at a time (`max_materialize`; one core each) |
| `--setup file` | | a file each job sources before TOFmeter starts (your Geant4 setup), on top of the submitting shell's environment |
| `--os-major n` | any | run only on workers of that OS major version (`TARGET.OpSysMajorVer`, e.g. 9 for EL9): the one TOFmeter was built for; `condor_status -af Machine OpSysAndVer` lists the workers |
| `--no-submit` | | write the macros and the submit file only |

Bands: `slow` is 0-1 GeV, seed 9500000 + k; `fast` is 1-1000 GeV, seed 9600000 + k.

Each job asks for 1 core and 1 GB of memory (TOFmeter uses ~120 MB). A job that fails, or ends without the run's
last line `Cosmic rays: <muons> events` (e.g. a macro command not found, after which TOFmeter still exits with 0), is
put on hold. So is a job still running after 23 h (the limit is 24 h).

The script refuses to submit while jobs of the same folder are in the queue. The 60-job limit holds per submission,
and this keeps a process from running twice.

## Following the jobs

```sh
condor_q                                      # running, idle, held
condor_q -hold                                # held jobs and why; see <dir>/<band>/log/job<k>.err and .log
condor_release <id>                           # run a held job again after a fix
condor_rm <id>                                # remove; then submit again to rerun what is not done
ls output/muonic_xray_czt/*/log/*.done | wc -l   # processes done
```

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
