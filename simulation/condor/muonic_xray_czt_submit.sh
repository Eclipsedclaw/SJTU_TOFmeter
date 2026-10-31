#!/bin/sh
# The CZT run of scripts/muonic_xray_czt.sh as HTCondor jobs, on a cluster with a shared filesystem such as SJTU's INPAC
# cluster (condor/README.md): one job of one core per process (TOFmeter runs sequentially), with the macros, seeds
# (9500000 + k slow, 9600000 + k fast), outputs and .done markers of the shared server's run, so that processes run in
# either place fit together process by process. At most 60 jobs at a time (max_materialize; 60 cores may be used), and
# no submission while jobs of the same folder are still in the queue (the limit is per submission, and no process may
# run twice). The jobs run muonic_xray_czt_job.sh, the one beside this script, in this build directory, and write into
# it. Their environment is that of this shell (getenv: set up Geant4 first); with --setup <file> each job also sources
# that file before TOFmeter starts, for when condor_submit does not start in a shell set up for Geant4 (submit from one
# without it). In the build directory:
#   sh condor/muonic_xray_czt_submit.sh [--dir path] [--muons n] [--max n] [--setup file] [--os-major n] [--no-submit] <band> <first> <last> [<band> <first> <last> ...]
#   sh condor/muonic_xray_czt_submit.sh --test [--setup file] [--os-major n] [--no-submit]
# <band> is slow (0-1 GeV) or fast (1-1000 GeV); of its processes first-last, those without <dir>/<band>/log/job<k>.done
# are submitted, e.g. slow 79 88 fast 47 168 (what the shared server's run had not finished on 2026-10-31). Defaults:
# output/muonic_xray_czt, 5 x 10^6 muons per process (~35 min slow, ~50 min fast on the shared server; a job must end
# within 24 h, one still running after 23 h is put on hold), 60 jobs at a time. --test: process 1 of each band with
# 2 x 10^5 muons into output/muonic_xray_czt_condor_test, a few minutes (the seeds of the shared server's process 1, so these
# are its first 2 x 10^5 muons: condor/muonic_xray_czt_compare.C). --no-submit: write the macros and the submit file
# only. Writes the macros as muonic_xray_czt.sh does (<dir>/<band>/log/job<k>.mac; beside them TOFmeter's log
# job<k>.log and the job's output and errors job<k>.out, job<k>.err) and per submission <dir>/condor/submit_<date>.sub
# with HTCondor's log of its jobs (submit_<date>.log). Follow them with condor_q; failed ones: condor_q -hold, then
# condor_release after a fix. --os-major n: run only on workers of that OS major version (TARGET.OpSysMajorVer,
# e.g. 9 for EL9, 7 for CentOS 7), the one TOFmeter was built on; condor_status -af Machine OpSysAndVer lists them.

dir=output/muonic_xray_czt
muons=5000000
max=60
setup=
osmajor=
submit=yes
test=no
while [ $# -gt 0 ]; do
  case "$1" in
    --dir) dir=$2; shift 2 ;;
    --muons) muons=$2; shift 2 ;;
    --max) max=$2; shift 2 ;;
    --setup) setup=$2; shift 2 ;;
    --os-major) osmajor=$2; shift 2 ;;
    --no-submit) submit=no; shift ;;
    --test) test=yes; shift ;;
    *) break ;;
  esac
done
if [ "$test" = yes ]; then
  dir=output/muonic_xray_czt_condor_test
  muons=200000
  set -- slow 1 1 fast 1 1
fi
if [ $# -eq 0 ] || [ $(($# % 3)) -ne 0 ]; then
  echo "usage: sh condor/muonic_xray_czt_submit.sh [--dir path] [--muons n] [--max n] [--setup file] [--os-major n] [--no-submit] <band> <first> <last> ..." >&2
  echo "       sh condor/muonic_xray_czt_submit.sh --test [--setup file] [--os-major n] [--no-submit]" >&2
  exit 1
fi
for v in "$muons" "$max"; do
  case "$v" in
    ''|*[!0-9]*) echo "--muons and --max take whole numbers" >&2; exit 1 ;;
  esac
done
if [ "$muons" -lt 1 ]; then
  echo "--muons: at least 1" >&2
  exit 1
fi
if [ "$max" -lt 1 ] || [ "$max" -gt 60 ]; then
  echo "--max: 1-60 jobs at a time (one core each, 60 cores at most)" >&2
  exit 1
fi
case "$osmajor" in
  *[!0-9]*) echo "--os-major: a whole number, e.g. 9" >&2; exit 1 ;;
esac
case "$dir$setup" in
  *[[:space:]]*) echo "--dir and --setup: paths without spaces" >&2; exit 1 ;;
esac
if [ -n "$setup" ]; then
  if [ ! -f "$setup" ] || [ ! -r "$setup" ]; then
    echo "--setup: cannot read $setup" >&2
    exit 1
  fi
  setup=$(cd "$(dirname "$setup")" && pwd -P)/$(basename "$setup")
fi
job=$(cd "$(dirname "$0")" && pwd -P)/muonic_xray_czt_job.sh
for f in ./TOFmeter macro/muonic_xray_czt.mac data/expacs_spectrum.txt "$job"; do
  if [ ! -e "$f" ]; then
    echo "Run this in the build directory ($f is missing)" >&2
    exit 1
  fi
done
name=tof_$(basename "$dir")
if [ "$submit" = yes ]; then
  if ! command -v condor_submit > /dev/null; then
    echo "no condor_submit here: run this on the cluster's login node (or with --no-submit)" >&2
    exit 1
  fi
  # HTCondor's programs must start in this shell; a shell set up for Geant4 (another system's libraries) can keep them
  # from it
  out=$(condor_submit -version 2>&1)
  case "$out" in
    *"error while loading shared libraries"*|*"not found (required by"*)
      printf '%s\n' "$out" >&2
      echo "condor_submit cannot start in this shell's environment; nothing was written or submitted. Run this in a shell" >&2
      echo "where it starts (condor_q is a quick test); if the jobs need more than that shell's environment, add" >&2
      echo "--setup <file>: each job sources the file before TOFmeter starts" >&2
      exit 1 ;;
  esac
  if condor_q -af JobBatchName 2> /dev/null | grep -qx "$name"; then
    echo "jobs of $dir (batch $name) are still in the queue (condor_q): submit when they are over, or condor_rm them" >&2
    exit 1
  fi
fi

# the macros of the processes to run, as scripts/muonic_xray_czt.sh writes them
items=
ranges=
n=0
finished=0
while [ $# -gt 0 ]; do
  band=$1
  first=$2
  last=$3
  shift 3
  case "$band" in
    slow) emin=0; emax=1; base=9500000 ;;
    fast) emin=1; emax=1000; base=9600000 ;;
    *) echo "band $band: slow or fast" >&2; exit 1 ;;
  esac
  for v in "$first" "$last"; do
    case "$v" in
      ''|*[!0-9]*) echo "$band $first $last: processes are numbered 1, 2, ..." >&2; exit 1 ;;
    esac
  done
  if [ "$first" -lt 1 ] || [ "$last" -lt "$first" ]; then
    echo "$band $first $last: first at least 1, last at least first" >&2
    exit 1
  fi
  sub="$dir/$band"
  mkdir -p "$sub/log"
  [ -e "$sub/started" ] || date +%s > "$sub/started"
  for k in $(seq "$first" "$last"); do
    if [ -e "$sub/log/job$k.done" ]; then
      finished=$((finished + 1))
      continue
    fi
    {
      echo "/control/verbose 2"
      echo "/control/alias czSeed $((base + k))"
      echo "/control/alias czEvents $muons"
      echo "/control/alias czDir $sub"
      echo "/control/alias czFile stk_$k"
      echo "/control/alias czEmin $emin"
      echo "/control/alias czEmax $emax"
      echo "/control/alias tgtCm 4"
      echo "/control/alias tgtZ 34.5"
      echo "/control/alias cztZ 32"
      echo "/control/execute macro/muonic_xray_czt.mac"
    } > "$sub/log/job$k.mac"
    items="$items$band $k
"
    n=$((n + 1))
  done
  ranges="$ranges${ranges:+, }$band $first-$last"
done
if [ "$n" -eq 0 ]; then
  echo "nothing to submit: $ranges of $dir are done ($finished .done markers)"
  [ "$test" = yes ] && echo "to test again: rm -r $dir"
  exit 0
fi

# the submit file: the jobs start in this directory (initialdir) and see it, nothing is transferred; bash, as setup
# files are mostly written for it
mkdir -p "$dir/condor"
here=$(pwd -P)
stamp=$(date +%Y%m%d_%H%M%S)
file="$dir/condor/submit_$stamp.sub"
{
  echo "# HTCondor jobs of the CZT run (macro/muonic_xray_czt.mac), by condor/muonic_xray_czt_submit.sh on $(date '+%F %T')"
  echo "# in $here: $ranges, $n processes of $muons muons ($finished done before), at most $max at a time;"
  echo "# environment: the submitting shell's${setup:+, then $setup sourced by each job}"
  echo "universe              = vanilla"
  echo "executable            = /bin/bash"
  echo "arguments             = \"$job $dir \$(band) \$(k)${setup:+ $setup}\""
  echo "transfer_executable   = False"
  echo "should_transfer_files = NO"
  echo "initialdir            = $here"
  echo "getenv                = True"
  [ -n "$osmajor" ] && echo "requirements          = (TARGET.OpSysMajorVer == $osmajor)"
  echo "request_cpus          = 1"
  echo "request_memory        = 1024"
  echo "notification          = Never"
  echo "batch_name            = $name"
  echo "max_materialize       = $max"
  echo "on_exit_hold          = ExitCode =!= 0"
  echo "on_exit_hold_reason   = \"TOFmeter failed or stopped early: see $dir/\$(band)/log/job\$(k).err and .log\""
  echo "periodic_hold         = (JobStatus == 2) && (time() - EnteredCurrentStatus > 23 * 3600)"
  echo "periodic_hold_reason  = \"running for over 23 h (jobs must end within 24 h)\""
  echo "log                   = $dir/condor/submit_$stamp.log"
  echo "output                = $dir/\$(band)/log/job\$(k).out"
  echo "error                 = $dir/\$(band)/log/job\$(k).err"
  echo "queue band, k from ("
  printf '%s' "$items"
  echo ")"
} > "$file"
echo "$(date '+%F %T')  $n processes ($ranges; $finished done before) of $muons muons, at most $max at a time: $file"
if [ "$submit" = no ]; then
  echo "not submitted (--no-submit); to submit: condor_submit $file"
  exit 0
fi
out=$(condor_submit "$file" 2>&1)
status=$?
printf '%s\n' "$out"
if [ "$status" -ne 0 ]; then
  echo "condor_submit FAILED (its message above)" >&2
  case "$out" in
    *aterializ*) echo "this HTCondor does not take max_materialize: delete that line of $file and submit at most 60 processes at once" >&2 ;;
  esac
  exit 1
fi
echo "follow: condor_q (batch $name); condor_q -hold (failed jobs, why); ls $dir/*/log/*.done | wc -l (processes done);"
echo "        cat $dir/<band>/log/job<k>.out (a job's host, time and speed)"
