#!/bin/sh
# Muonic X-rays in the CZT plane under the 4 cm Al target: TOFmeter processes of macro/muonic_xray_czt.mac (Si/YSO
# stack out, CZT plane of 4 x 4 pixels of 25 x 25 x 10 mm in its place, its top against the target's bottom at 32.5 cm),
# in two bands of cosmic muons: the slow ones (0 - 1 GeV), which hold nearly every stop in the target, and the faster
# ones (1 - 1000 GeV), background only; the old runs had no CZT, so both are needed. Fast processes first (they take
# longest), at most 90% of the physical cores at a time (10 on the shared server, the default). In the build directory:
#   sh scripts/muonic_xray_czt.sh [--dir path] [slow processes] [fast processes] [muons per process] [at a time]
# Defaults: output/muonic_xray_czt, 8 slow and 28 fast processes of 5 x 10^6 muons: about 1 day of exposure each
# (5 x 10^6 muons of 0-1 GeV stand for 2.9 h, of 1-1000 GeV for 0.85 h), a test run of about 2.7 h on 10 cores
# (estimated from the stack runs: ~31 and ~48 min per process). Writes <dir>/slow/stk_<k>.root and <dir>/fast/stk_<k>.root
# with a macro and log per process in <dir>/slow/log/ and <dir>/fast/log/. A finished process leaves log/job<k>.done
# and is skipped when the script is started again (so more processes can be added later); touch <dir>/STOP: no new
# process starts. In the background:
#   nohup nice -n 10 sh scripts/muonic_xray_czt.sh >> output/muonic_xray_czt.log 2>&1 &

dir=output/muonic_xray_czt
while [ $# -gt 0 ]; do
  case "$1" in
    --dir) dir=$2; shift 2 ;;
    *) break ;;
  esac
done
slow=${1:-8}
fast=${2:-28}
events=${3:-5000000}
# at most 90% of the physical cores (distinct core/socket pairs; hyperthreads do not count)
cores=$(lscpu -p=Core,Socket | grep -v '^#' | sort -u | wc -l)
limit=$((cores * 9 / 10))
[ "$limit" -ge 1 ] || limit=1
parallel=${4:-$limit}
if [ "$parallel" -gt "$limit" ]; then
  echo "at most $limit processes at a time (90% of the physical cores; asked for $parallel)" >&2
  parallel=$limit
fi
for f in ./TOFmeter macro/muonic_xray_czt.mac; do
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

# 4 cm target at 32.5-36.5 cm, the CZT plane's centre 5 mm under its bottom face; seeds 9500000 (slow) or
# 9600000 (fast) + k
for band in fast slow; do
  sub="$dir/$band"
  mkdir -p "$sub/log"
  [ -e "$sub/started" ] || date +%s > "$sub/started"
  if [ "$band" = slow ]; then n=$slow; emin=0; emax=1; base=9500000; else n=$fast; emin=1; emax=1000; base=9600000; fi
  for k in $(seq 1 "$n"); do
    {
      echo "/control/verbose 2"
      echo "/control/alias czSeed $((base + k))"
      echo "/control/alias czEvents $events"
      echo "/control/alias czDir $sub"
      echo "/control/alias czFile stk_$k"
      echo "/control/alias czEmin $emin"
      echo "/control/alias czEmax $emax"
      echo "/control/alias tgtCm 4"
      echo "/control/alias tgtZ 34.5"
      echo "/control/alias cztZ 32"
      echo "/control/execute macro/muonic_xray_czt.mac"
    } > "$sub/log/job$k.mac"
  done
  echo "$(date '+%F %T')  $band: $n x $events muons of $emin-$emax GeV -> $sub"
done
todo=$(for k in $(seq 1 "$fast"); do [ -e "$dir/fast/log/job$k.done" ] || echo "fast $k"; done
       for k in $(seq 1 "$slow"); do [ -e "$dir/slow/log/job$k.done" ] || echo "slow $k"; done)
echo "$(date '+%F %T')  $(printf '%s\n' "$todo" | grep -c .) processes to run, $parallel at a time"
export TOF_CZT_DIR="$dir" TOF_CZT_MAX="$parallel"
# Each process starts only while fewer than TOF_CZT_MAX processes of this run are going (those of an earlier driver
# included), and not at all if one is still running it; start.lock makes the starts one at a time.
printf '%s\n' "$todo" | grep . | xargs -P "$parallel" -n 2 sh -c '
  job="$TOF_CZT_DIR/$1/log/job$2"
  [ -e "$TOF_CZT_DIR/STOP" ] && exit 0
  [ -e "$job.done" ] && exit 0
  # ^./TOFmeter: the processes themselves, not a shell or ps whose command line holds the same words
  pgrep -f "^\./TOFmeter $job\.mac\$" > /dev/null && exit 0
  until mkdir "$TOF_CZT_DIR/start.lock" 2> /dev/null; do sleep 1; done
  while [ "$(pgrep -c -f "^\./TOFmeter $TOF_CZT_DIR/")" -ge "$TOF_CZT_MAX" ]; do sleep 30; done
  if [ -e "$TOF_CZT_DIR/STOP" ]; then rmdir "$TOF_CZT_DIR/start.lock"; exit 0; fi
  ./TOFmeter "$job.mac" > "$job.log" 2>&1 &
  pid=$!
  sleep 2   # so that the next start counts this one
  rmdir "$TOF_CZT_DIR/start.lock"
  if wait $pid; then
    touch "$job.done"
    echo "$(date "+%F %T")  $1 job $2 done ($(ls "$TOF_CZT_DIR"/*/log | grep -c "\.done$") finished)"
  else
    echo "$(date "+%F %T")  $1 job $2 FAILED, see $job.log"
  fi' sh

echo "$(date '+%F %T')  finished; $(ls "$dir"/*/log | grep -c '\.done$') processes done"
