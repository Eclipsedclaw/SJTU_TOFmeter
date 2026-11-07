#!/bin/sh
# One HTCondor job of the beam study (2026-11-07), as condor/beam_xray_submit.sh writes and submits them: process <k> of
# the run in <dir>, i.e. ./TOFmeter --preinit macro/preinit_lyso.mac <dir>/beam/log/job<k>.mac, in the build directory
# (the job's initialdir; the cluster's filesystem is shared, so nothing is transferred). The environment is the
# submitting shell's (getenv), with the setup file sourced on top if given. TOFmeter's output goes to
# <dir>/beam/log/job<k>.log. A process that ran to its end leaves <dir>/beam/log/job<k>.done and is not run again. Its
# end means the log's "Number of events processed : <muons>" and a non-empty <dir>/beam/beam_<k>.root. Otherwise the
# exit status is not 0 (TOFmeter's, or 1) and the job goes on hold (condor_q -hold).
#   bash condor/beam_xray_job.sh <dir> <band> <k> [setup file]      (<band> is always beam)
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
muons=$(sed -n 's/^\/control\/alias bmEvents //p' "$job.mac")
preinit=$(sed -n 's/^# preinit //p' "$job.mac")
start=$(date +%s)
echo "$(date '+%F %T')  $band process $k: $muons muons on $(uname -n) in $(pwd), environment $environment${preinit:+, --preinit $preinit}"
./TOFmeter ${preinit:+--preinit "$preinit"} "$job.mac" > "$job.log" 2>&1
status=$?
s=$(($(date +%s) - start))
geant4=$(sed -n 's/^ *Geant4 version Name: *\([^ ]*\).*/\1/p' "$job.log" | head -n 1)
if [ "$status" -eq 0 ] && grep -q "Number of events processed : $muons\$" "$job.log" && [ -s "$dir/$band/beam_$k.root" ]; then
  touch "$job.done"
  # the speed, to size the jobs (24 h at most), with the start-up included
  awk -v n="$muons" -v s="$s" -v g="$geant4" -v d="$(date '+%F %T')" 'BEGIN { if (s < 1) s = 1
    printf "%s  done in %d s with %s: %.0f muons/s with the start-up, 10^6 muons in %.2f h or less\n",
           d, s, g, n / s, 1e6 / n * s / 3600 }'
  exit 0
fi
echo "$(date '+%F %T')  $band process $k FAILED after $s s (TOFmeter's exit status $status); the end of $job.log:" >&2
tail -n 20 "$job.log" >&2
[ "$status" -ne 0 ] && exit "$status"
exit 1
