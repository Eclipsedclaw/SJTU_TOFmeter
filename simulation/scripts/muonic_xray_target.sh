#!/bin/sh
# Muonic X-rays in the Si/YSO stack for thicker Al targets: TOFmeter processes of macro/muonic_xray_target.mac, as the
# 3-week run of scripts/muonic_xray_stack.sh (1 cm target) but with only the cosmic muons below 1 GeV, which hold
# nearly every muon that stops in the target. Several thicknesses in one run: their processes take turns, so every
# thickness has about the same exposure at any time. At most 90% of the physical cores at a time (hyperthreads do not
# count: 10 on the shared server, the default). Run it in the build directory:
#   sh scripts/muonic_xray_target.sh [--dir path] [--emax GeV] [thicknesses in cm] [processes per thickness]
#                                    [muons per process] [processes at a time]
# Defaults: output/muonic_xray_target, 1 GeV, "2 3 4 5", 29 processes of 10^7 muons per thickness (5.8 h of exposure
# each, 7.0 days per thickness), 90% of the physical cores. A target of d cm is 20 x 20 x d cm of Al with its top at
# 36.5 cm, as the default 1 cm one, so it grows down towards the stack; one that would come closer than 1 mm to the
# stack (its top Si ends at 31.6 cm) is raised to keep the 1 mm (5 cm: 31.7 - 36.7 cm). Writes <dir>/t<d>cm/stk_<k>.root
# with a macro and log per process in <dir>/t<d>cm/log/: the layout of output/muonic_xray_stack, so the status script,
# the skim and the cut macros work on each thickness, e.g.
#   sh scripts/muonic_xray_stack_status.sh --dir output/muonic_xray_target/t2cm
# Measured in the default run (10 processes at a time): about 61 min per process of 10^7 muons (0.37 ms per muon, the
# start-up included); a muon below 1 GeV stands for 2.08 ms of exposure, so a day of exposure takes about 25 min of
# 10 cores (12.1 h for the defaults); ~0.1 GB per day of exposure on disk.
# A finished process leaves <dir>/t<d>cm/log/job<k>.done and is skipped when the script is started again (so more
# processes can be added later, or an interrupted run goes on); touch <dir>/STOP: no new process starts (delete it and
# start again to go on). In the background:
#   nohup nice -n 10 sh scripts/muonic_xray_target.sh >> output/muonic_xray_target.log 2>&1 &

dir=output/muonic_xray_target; emax=1
while [ $# -gt 0 ]; do
  case "$1" in
    --dir) dir=$2; shift 2 ;;
    --emax) emax=$2; shift 2 ;;
    *) break ;;
  esac
done
thicknesses=${1:-"2 3 4 5"}
processes=${2:-29}
events=${3:-10000000}
# at most 90% of the physical cores (distinct core/socket pairs; hyperthreads do not count)
cores=$(lscpu -p=Core,Socket | grep -v '^#' | sort -u | wc -l)
limit=$((cores * 9 / 10))
[ "$limit" -ge 1 ] || limit=1
parallel=${4:-$limit}
if [ "$parallel" -gt "$limit" ]; then
  echo "at most $limit processes at a time (90% of the physical cores; asked for $parallel)" >&2
  parallel=$limit
fi
for f in ./TOFmeter macro/muonic_xray_target.mac; do
  if [ ! -e "$f" ]; then
    echo "Run this in the build directory ($f is missing)" >&2
    exit 1
  fi
done
mkdir -p "$dir"
if ! mkdir "$dir/driver.lock" 2>/dev/null; then
  echo "$dir/driver.lock exists: a run in $dir is going on (or one was killed: then delete the lock)" >&2
  exit 1
fi
trap 'rmdir "$dir/driver.lock" 2>/dev/null' EXIT
trap 'exit 1' INT TERM

for d in $thicknesses; do
  sub="$dir/t${d}cm"
  mkdir -p "$sub/log"
  [ -e "$sub/started" ] || date +%s > "$sub/started"   # for scripts/muonic_xray_stack_status.sh
  # centre: top at 36.5 cm, or 1 mm above the stack; seeds 9300000 + 1000 x (thickness in mm) + k
  z=$(awk -v d="$d" 'BEGIN { top = 36.5; if (31.7 + d > top) top = 31.7 + d; printf "%.4g", top - 0.5 * d }')
  base=$(awk -v d="$d" 'BEGIN { printf "%d", 9300000 + 1000 * int(10 * d + 0.5) }')
  for k in $(seq 1 "$processes"); do
    {
      echo "/control/verbose 2"
      echo "/control/alias stkSeed $((base + k))"
      echo "/control/alias stkEvents $events"
      echo "/control/alias stkDir $sub"
      echo "/control/alias stkFile stk_$k"
      echo "/control/alias tgtCm $d"
      echo "/control/alias tgtZ $z"
      echo "/control/alias tgtEmax $emax"
      echo "/control/execute macro/muonic_xray_target.mac"
    } > "$sub/log/job$k.mac"
  done
  echo "$(date '+%F %T')  $d cm target, centre at z = $z cm: $processes x $events muons below $emax GeV -> $sub"
done
# the thicknesses take turns: job 1 of each, then job 2 of each, ...
todo=$(for k in $(seq 1 "$processes"); do
         for d in $thicknesses; do [ -e "$dir/t${d}cm/log/job$k.done" ] || echo "t${d}cm $k"; done
       done)
echo "$(date '+%F %T')  $(printf '%s\n' "$todo" | grep -c .) processes to run, $parallel at a time"
export TOF_TGT_DIR="$dir" TOF_TGT_MAX="$parallel"
# Each process starts only while fewer than TOF_TGT_MAX processes of this run are going (those of an earlier driver
# included), and not at all if one is still running it; start.lock makes the starts one at a time.
printf '%s\n' "$todo" | grep . | xargs -P "$parallel" -n 2 sh -c '
  job="$TOF_TGT_DIR/$1/log/job$2"
  [ -e "$TOF_TGT_DIR/STOP" ] && exit 0
  [ -e "$job.done" ] && exit 0
  # ^./TOFmeter: the processes themselves, not a shell or ps whose command line holds the same words
  pgrep -f "^\./TOFmeter $job\.mac\$" > /dev/null && exit 0
  until mkdir "$TOF_TGT_DIR/start.lock" 2> /dev/null; do sleep 1; done
  while [ "$(pgrep -c -f "^\./TOFmeter $TOF_TGT_DIR/t")" -ge "$TOF_TGT_MAX" ]; do sleep 30; done
  if [ -e "$TOF_TGT_DIR/STOP" ]; then rmdir "$TOF_TGT_DIR/start.lock"; exit 0; fi
  ./TOFmeter "$job.mac" > "$job.log" 2>&1 &
  pid=$!
  sleep 2   # so that the next start counts this one
  rmdir "$TOF_TGT_DIR/start.lock"
  if wait $pid; then
    touch "$job.done"
    echo "$(date "+%F %T")  $1 job $2 done ($(ls "$TOF_TGT_DIR"/t*/log | grep -c "\.done$") finished)"
  else
    echo "$(date "+%F %T")  $1 job $2 FAILED, see $job.log"
  fi' sh

echo "$(date '+%F %T')  finished; $(ls "$dir"/t*/log | grep -c '\.done$') processes done"
