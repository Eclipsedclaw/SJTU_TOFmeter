#!/bin/sh
# Muonic X-rays of cosmic muons stopped in the target: many TOFmeter processes of
# macro/muonic_xray.mac, several at a time, then analysis/muonic_xray.C. Run it in the build directory:
#   sh scripts/muonic_xray.sh [--half cm] [--dir path] [--geometry macro] [--band "Emin Emax"]
#                             [processes] [muons per process] [processes at a time]
# Defaults: plane half size 200 cm (4 x 4 m at z = 112 cm, world 4.2 x 4.2 m), output/muonic_xray,
# macro/geom_target_al.mac, muons of 0.1 - 2 GeV, 30 processes of 5 x 10^6 muons (1.68 h of
# sea-level exposure each, 50 h in all, about 1.5 h on 10 cores), 10 at a time. Writes
# <dir>/mux_<k>.root with a macro and log per process in
# <dir>/log/; at the end analysis/muonic_xray.C makes <dir>/muonic_xray_*.png/.pdf/.csv for the
# target and for the Al plates of the frame.
# A process that finishes leaves <dir>/log/job<k>.done; those are skipped when the script is started
# again (so more processes can be added later, or an interrupted run goes on). touch <dir>/STOP: no
# new process starts; delete it and start again to go on. Long runs, in the background:
#   nohup nice -n 15 caffeinate -is sh scripts/muonic_xray.sh 60 >> output/muonic_xray.log 2>&1 &

half=200; dir=output/muonic_xray; geometry=macro/geom_target_al.mac; band="0.1 2.0"
while [ $# -gt 0 ]; do
  case "$1" in
    --half) half=$2; shift 2 ;;
    --dir) dir=$2; shift 2 ;;
    --geometry) geometry=$2; shift 2 ;;
    --band) band=$2; shift 2 ;;
    *) break ;;
  esac
done
processes=${1:-30}
events=${2:-5000000}
parallel=${3:-10}
for f in ./TOFmeter macro/muonic_xray.mac "$geometry" analysis/muonic_xray.C; do
  if [ ! -e "$f" ]; then
    echo "Run this in the build directory ($f is missing)" >&2
    exit 1
  fi
done
emin=${band% *}; emax=${band#* }
world=$(awk -v h="$half" 'BEGIN { w = 2 * h + 20; print (w > 400 ? w : 400) }')
mkdir -p "$dir/log"
if ! mkdir "$dir/driver.lock" 2>/dev/null; then
  echo "$dir/driver.lock exists: a run in $dir is going on (or one was killed: then delete the lock)" >&2
  exit 1
fi
trap 'rmdir "$dir/driver.lock" 2>/dev/null' EXIT
trap 'exit 1' INT TERM

for k in $(seq 1 "$processes"); do
  {
    echo "/control/verbose 2"
    echo "/control/alias muxSeed $((8000000 + k))"
    echo "/control/alias muxEvents $events"
    echo "/control/alias muxDir $dir"
    echo "/control/alias muxFile mux_$k"
    echo "/control/alias muxGeometry $geometry"
    echo "/control/alias muxHalf $half"
    echo "/control/alias muxWorld $world"
    echo "/control/alias muxEmin $emin"
    echo "/control/alias muxEmax $emax"
    echo "/control/execute macro/muonic_xray.mac"
  } > "$dir/log/job$k.mac"
done
todo=$(for k in $(seq 1 "$processes"); do [ -e "$dir/log/job$k.done" ] || echo "$k"; done)
echo "$(date '+%F %T')  $processes x $events muons of $emin-$emax GeV, plane $((2 * half)) x $((2 * half)) cm," \
     "$geometry -> $dir; $(echo $todo | wc -w | tr -d ' ') processes to run, $parallel at a time"
export TOF_MUX_DIR="$dir"
echo $todo | tr ' ' '\n' | grep . | xargs -P "$parallel" -n 1 sh -c '
  [ -e "$TOF_MUX_DIR/STOP" ] && exit 0
  [ -e "$TOF_MUX_DIR/log/job$1.done" ] && exit 0
  if ./TOFmeter "$TOF_MUX_DIR/log/job$1.mac" > "$TOF_MUX_DIR/log/job$1.log" 2>&1; then
    touch "$TOF_MUX_DIR/log/job$1.done"
    echo "$(date "+%F %T")  job $1 done ($(ls "$TOF_MUX_DIR/log" | grep -c "\.done$") finished)"
  else
    echo "$(date "+%F %T")  job $1 FAILED, see $TOF_MUX_DIR/log/job$1.log"
  fi' sh

echo "$(date '+%F %T')  finished; $(ls "$dir/log" | grep -c '\.done$') of $processes processes done"
command -v root > /dev/null || exit 0
root -l -b -q "analysis/muonic_xray.C(\"$dir\", \"Target\")" > "$dir/muonic_xray_Target.txt" 2>&1
root -l -b -q "analysis/muonic_xray.C(\"$dir\", \"FrameworkAlPhys\")" > "$dir/muonic_xray_FrameworkAlPhys.txt" 2>&1
echo "tables: $dir/muonic_xray_Target.txt, $dir/muonic_xray_FrameworkAlPhys.txt"
