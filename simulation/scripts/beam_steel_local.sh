#!/bin/sh
# Local runs of the steel beam study (2026-11-09): mu- of 85 MeV/c into a 1 cm stainless-steel plate lying on the LYSO
# panel (macro/beam_steel_lyso.mac, TOFmeter --preinit macro/preinit_lyso.mac), several TOFmeter processes at the same
# time on this computer, one thread each, at most 90% of the physical cores.
#   job k: <muons> muons with the seed <seed-base> + k into <dir>/beam/beam_<k>.root; its macro and log in
#   <dir>/beam/log/job<k>.mac and .log, and job<k>.done once it ran to its end ("Number of events processed : <muons>"
#   in the log and a non-empty ROOT file), as the Condor jobs of the Si study. Jobs with a .done are skipped.
# Every event is written (/EventFilter 0; --filter 8 keeps only those with energy in the panel), so that
# analysis/beam_steel_spectrum.C counts the muons each element caught. Refuses to start at 85% memory use or more.
# In the build directory, with Geant4 and ROOT set up:
#   sh scripts/beam_steel_local.sh [--jobs 10] [--first 1] [--muons 1000000] [--seed-base 10800000]
#        [--material G4_STAINLESS-STEEL] [--ekin 29.95] [--filter 0] [--dir output/lyso/beam_steel_lyso]
#        [--detector lyso|czt]
# --detector czt: the CZT plane (macro/panel_czt.mac, 4 x 4 pixels of 25 x 25 x 10 mm, its top against the plate), into
# output/czt/beam_steel_czt by default.
# Then: root -l -b -q 'analysis/beam_steel_spectrum.C("output/lyso/beam_steel_lyso")'
jobs=10; first=1; muons=1000000; seedbase=10800000; material=G4_STAINLESS-STEEL; ekin=29.95; filter=0
dir=""; detector=lyso; sigmax=1.656; sigmay=1.062
while [ $# -gt 0 ]; do
  case "$1" in
    --jobs) jobs=$2; shift 2 ;;
    --first) first=$2; shift 2 ;;
    --muons) muons=$2; shift 2 ;;
    --seed-base) seedbase=$2; shift 2 ;;
    --material) material=$2; shift 2 ;;
    --ekin) ekin=$2; shift 2 ;;
    --filter) filter=$2; shift 2 ;;
    --dir) dir=$2; shift 2 ;;
    --detector) detector=$2; shift 2 ;;
    *) echo "usage: sh scripts/beam_steel_local.sh [--jobs n] [--first k] [--muons n] [--seed-base b] [--material name] [--ekin MeV] [--filter 0|8] [--dir path] [--detector lyso|czt]" >&2
       exit 1 ;;
  esac
done
for v in "$jobs" "$first" "$muons" "$seedbase" "$filter"; do
  case "$v" in
    ''|*[!0-9]*) echo "--jobs, --first, --muons, --seed-base and --filter take whole numbers" >&2; exit 1 ;;
  esac
done
case "$detector" in
  lyso) panelz=32.332 ;;   # the centre of the 3 mm LYSO layer under the 0.1 mm Al and 0.08 mm ESR
  czt) panelz=32.0 ;;      # the centre of the 10 mm CZT, its top against the plate at 32.5 cm
  *) echo "--detector is lyso or czt" >&2; exit 1 ;;
esac
[ -n "$dir" ] || dir=output/$detector/beam_steel_$detector
limit=$(( $(lscpu -p=Core,Socket | grep -v '^#' | sort -u | wc -l) * 9 / 10 ))
if [ "$jobs" -lt 1 ] || [ "$jobs" -gt "$limit" ]; then
  echo "--jobs must be 1 to $limit (90% of the physical cores)" >&2
  exit 1
fi
for f in ./TOFmeter macro/beam_steel_lyso.mac macro/preinit_lyso.mac macro/panel_$detector.mac; do
  if [ ! -e "$f" ]; then
    echo "Run this in the build directory ($f is missing)" >&2
    exit 1
  fi
done
used=$(awk '/^MemTotal/{t=$2} /^MemAvailable/{a=$2} END{printf "%d", 100*(t-a)/t}' /proc/meminfo)
if [ "$used" -ge 85 ]; then
  echo "memory: $used% in use, not starting (85% or more)" >&2
  exit 1
fi

log=$dir/beam/log
mkdir -p "$log"
start=$(date +%s)
last=$((first + jobs - 1))
pids=""
for k in $(seq "$first" "$last"); do
  job=$log/job$k
  if [ -e "$job.done" ]; then
    echo "job $k is done already ($job.done)"
    continue
  fi
  {
    echo "/control/verbose 2"
    echo "/control/alias bmSeed $((seedbase + k))"
    echo "/control/alias bmEvents $muons"
    echo "/control/alias bmDir $dir/beam"
    echo "/control/alias bmFile beam_$k"
    echo "/control/alias bmEkin $ekin"
    echo "/control/alias bmSigX $sigmax"
    echo "/control/alias bmSigY $sigmay"
    echo "/control/alias bmZ 33.55"
    echo "/control/alias bmFilter $filter"
    echo "/control/alias bmPanel $detector"
    echo "/control/alias bmPanelZ $panelz"
    echo "/tof/target/material $material"
    echo "/control/execute macro/beam_steel_lyso.mac"
  } > "$job.mac"
  (
    OMP_NUM_THREADS=1 ./TOFmeter --preinit macro/preinit_lyso.mac "$job.mac" > "$job.log" 2>&1
    status=$?
    if [ "$status" -eq 0 ] && grep -q "Number of events processed : $muons\$" "$job.log" && [ -s "$dir/beam/beam_$k.root" ]; then
      touch "$job.done"
      exit 0
    fi
    echo "job $k FAILED (TOFmeter's exit status $status): see $job.log" >&2
    exit 1
  ) &
  pids="$pids $!"
done
echo "$(date '+%F %T')  started jobs $first-$last: $muons muons of $ekin MeV each into 1 cm of $material on the $detector panel ($dir/beam, logs in $log)"
status=0
for p in $pids; do
  wait "$p" || status=1
done
echo "$(date '+%F %T')  finished after $(( $(date +%s) - start )) s: $(ls "$log"/job*.done 2>/dev/null | wc -l) jobs done in $log"
exit $status
