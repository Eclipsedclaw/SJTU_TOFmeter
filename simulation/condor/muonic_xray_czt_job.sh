#!/bin/sh
# One HTCondor job of the CZT run, as condor/muonic_xray_czt_submit.sh writes and submits them: process <k> of the
# band <band> (slow or fast) of the run in <dir>, i.e. ./TOFmeter <dir>/<band>/log/job<k>.mac as scripts/muonic_xray_czt.sh
# runs it on the shared server, in the build directory (the job's initialdir; the cluster's filesystem is shared, so
# nothing is transferred), with the environment of the shell that submitted it (getenv) and, if given, the setup file
# sourced on top (a file of yours that sets up Geant4). TOFmeter's output goes to <dir>/<band>/log/job<k>.log. A process
# that ran to its end (the run's last line "Cosmic rays: <muons> events" in that log, and its stk_<k>.root) leaves
# <dir>/<band>/log/job<k>.done and is not run again; otherwise the exit status is not 0 (TOFmeter's, or 1) and the job
# goes on hold (condor_q -hold).
#   bash condor/muonic_xray_czt_job.sh <dir> <band> <k> [setup file]
dir=$1
band=$2
k=$3
setup=$4
job=$dir/$band/log/job$k
environment="the submitting shell's"
if [ -n "$setup" ]; then
  if [ ! -r "$setup" ]; then
    echo "cannot read the setup file $setup on $(uname -n)" >&2
    exit 1
  fi
  here=$(pwd)
  set --   # so that the setup file does not take this script's arguments for its own
  . "$setup"
  cd "$here" || exit 1
  environment="$environment and $setup"
fi
if [ ! -x ./TOFmeter ] || [ ! -e "$job.mac" ]; then
  echo "no ./TOFmeter or $job.mac in $(pwd) on $(uname -n)" >&2
  exit 1
fi
if [ -e "$job.done" ]; then
  echo "$band process $k is done already ($job.done)"
  exit 0
fi
muons=$(sed -n 's/^\/control\/alias czEvents //p' "$job.mac")
start=$(date +%s)
echo "$(date '+%F %T')  $band process $k: $muons muons on $(uname -n) in $(pwd), environment $environment"
./TOFmeter "$job.mac" > "$job.log" 2>&1
status=$?
s=$(($(date +%s) - start))
geant4=$(sed -n 's/^ *Geant4 version Name: *\([^ ]*\).*/\1/p' "$job.log" | head -n 1)
if [ "$status" -eq 0 ] && grep -q "^Cosmic rays: $muons events" "$job.log" && [ -s "$dir/$band/stk_$k.root" ]; then
  touch "$job.done"
  # the speed, to size the jobs (24 h at most); with the start-up (in the slow band ~1 min: the first captures compute
  # each isotope's cascade), so for more muons than this job's an upper limit
  awk -v n="$muons" -v s="$s" -v g="$geant4" -v d="$(date '+%F %T')" 'BEGIN { if (s < 1) s = 1
    printf "%s  done in %d s with %s: %.0f muons/s with the start-up, 5 x 10^6 muons in %.2f h or less\n",
           d, s, g, n / s, 5e6 / n * s / 3600 }'
  exit 0
fi
echo "$(date '+%F %T')  $band process $k FAILED after $s s (TOFmeter's exit status $status); the end of $job.log:" >&2
tail -n 20 "$job.log" >&2
[ "$status" -ne 0 ] && exit "$status"
exit 1
