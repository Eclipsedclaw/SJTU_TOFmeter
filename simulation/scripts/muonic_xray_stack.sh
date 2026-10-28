#!/bin/sh
# Muonic X-rays as the Si/YSO stack sees them: many TOFmeter processes of macro/muonic_xray_stack.mac, each a separate
# stretch of sea-level exposure to cosmic muons of all energies through the default setup with its Al target, at most
# 90% of the physical cores at a time (hyperthreads do not count: 10 on the shared server, the default); a run can
# be taken over by a driver with another number at a time (stop the old driver's xargs, start the script again).
# Run it in the build directory:
#   sh scripts/muonic_xray_stack.sh [--dir path] [--all-cores] [processes] [muons per process] [processes at a time]
# Defaults: output/muonic_xray_stack, 385 processes of 10^7 muons (1.31 h of exposure each, 3 weeks in all), 90% of
# the physical cores at a time. Writes <dir>/stk_<k>.root with a macro and log per process in <dir>/log/. Measured
# on this server: about 0.7 ms per muon alone, ~2.7 h per process with 12 running, so the 3 weeks take about 4 days
# at 10 at a time; about 3 kB per second of exposure on disk (12 GB for 3 weeks).
# A process that finishes leaves <dir>/log/job<k>.done; those are skipped when the script is started again (so more
# processes can be added later, or an interrupted run goes on). touch <dir>/STOP: no new process starts; delete it
# and start again to go on. Every finished process is a complete sample of its own: analyse what is done at any
# time. Long runs, in the background:
#   nohup nice -n 10 sh scripts/muonic_xray_stack.sh >> output/muonic_xray_stack.log 2>&1 &

dir=output/muonic_xray_stack; allcores=0
while [ $# -gt 0 ]; do
  case "$1" in
    --dir) dir=$2; shift 2 ;;
    --all-cores) allcores=1; shift ;;   # up to every physical core: an exception, e.g. on an otherwise idle server
    *) break ;;
  esac
done
processes=${1:-385}
events=${2:-10000000}
# at most 90% of the physical cores (distinct core/socket pairs; hyperthreads do not count: more processes than
# physical cores run slower in total), all of them with --all-cores
cores=$(lscpu -p=Core,Socket | grep -v '^#' | sort -u | wc -l)
limit=$((cores * 9 / 10))
[ "$allcores" -eq 1 ] && limit=$cores
[ "$limit" -ge 1 ] || limit=1
parallel=${3:-$limit}
if [ "$parallel" -gt "$limit" ]; then
  echo "at most $limit processes at a time (90% of the physical cores; asked for $parallel)" >&2
  parallel=$limit
fi
for f in ./TOFmeter macro/muonic_xray_stack.mac; do
  if [ ! -e "$f" ]; then
    echo "Run this in the build directory ($f is missing)" >&2
    exit 1
  fi
done
mkdir -p "$dir/log"
if ! mkdir "$dir/driver.lock" 2>/dev/null; then
  echo "$dir/driver.lock exists: a run in $dir is going on (or one was killed: then delete the lock)" >&2
  exit 1
fi
trap 'rmdir "$dir/driver.lock" 2>/dev/null' EXIT
trap 'exit 1' INT TERM
[ -e "$dir/started" ] || date +%s > "$dir/started"   # start of the run, for scripts/muonic_xray_stack_status.sh

for k in $(seq 1 "$processes"); do
  {
    echo "/control/verbose 2"
    echo "/control/alias stkSeed $((9200000 + k))"
    echo "/control/alias stkEvents $events"
    echo "/control/alias stkDir $dir"
    echo "/control/alias stkFile stk_$k"
    echo "/control/execute macro/muonic_xray_stack.mac"
  } > "$dir/log/job$k.mac"
done
todo=$(for k in $(seq 1 "$processes"); do [ -e "$dir/log/job$k.done" ] || echo "$k"; done)
echo "$(date '+%F %T')  $processes x $events cosmic muons of all energies -> $dir;" \
     "$(echo $todo | wc -w | tr -d ' ') processes to run, $parallel at a time"
export TOF_STK_DIR="$dir" TOF_STK_MAX="$parallel"
# Each process starts only while fewer than TOF_STK_MAX processes of this run are going, those of an earlier
# driver included (so a run can be taken over by a driver with another number at a time without waiting), and
# not at all if an earlier driver is still running it. start.lock makes the starts one at a time.
echo $todo | tr ' ' '\n' | grep . | xargs -P "$parallel" -n 1 sh -c '
  [ -e "$TOF_STK_DIR/STOP" ] && exit 0
  [ -e "$TOF_STK_DIR/log/job$1.done" ] && exit 0
  # ^./TOFmeter: the processes themselves, not a shell or ps whose command line holds the same words
  pgrep -f "^\./TOFmeter $TOF_STK_DIR/log/job$1\.mac\$" > /dev/null && exit 0
  until mkdir "$TOF_STK_DIR/start.lock" 2> /dev/null; do sleep 1; done
  while [ "$(pgrep -c -f "^\./TOFmeter $TOF_STK_DIR/log/job")" -ge "$TOF_STK_MAX" ]; do sleep 30; done
  if [ -e "$TOF_STK_DIR/STOP" ]; then rmdir "$TOF_STK_DIR/start.lock"; exit 0; fi
  ./TOFmeter "$TOF_STK_DIR/log/job$1.mac" > "$TOF_STK_DIR/log/job$1.log" 2>&1 &
  pid=$!
  sleep 2   # so that the next start counts this one
  rmdir "$TOF_STK_DIR/start.lock"
  if wait $pid; then
    touch "$TOF_STK_DIR/log/job$1.done"
    echo "$(date "+%F %T")  job $1 done ($(ls "$TOF_STK_DIR/log" | grep -c "\.done$") finished)"
  else
    echo "$(date "+%F %T")  job $1 FAILED, see $TOF_STK_DIR/log/job$1.log"
  fi' sh

echo "$(date '+%F %T')  finished; $(ls "$dir/log" | grep -c '\.done$') of $processes processes done"
