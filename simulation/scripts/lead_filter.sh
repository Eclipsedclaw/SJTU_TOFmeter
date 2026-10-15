#!/bin/sh
# Muon spectrum before and after lead slabs of several thicknesses (macro/lead_filter.mac), then
# analysis/lead_filter_spectrum.C. Run it in the build directory:
#   sh scripts/lead_filter.sh [--below E] [--parts n] [--dir path] [jobs] [muons per thickness] [thicknesses in cm]
# Defaults: 6 jobs, 2000000 muons, 5 10 20 30 40 50 cm, the whole spectrum 9.25 MeV - 108 GeV,
# output in output/lead_filter. Every event is written, about 0.17 GB per 10^6 muons.
# --below E (GeV): only the muons that can leave the lead below E, i.e. up to
#   E + 0.2 GeV + 0.07 GeV per cm of lead (the most a muon at 75 deg loses), for many more
#   statistics below E at the same cost; the analysis then shows the spectrum below E only.
#   Output in output/lead_filter_below<E>GeV.
# --parts n: each thickness as n processes with different seeds (the muons are shared), which the
#   analysis adds up; with more cores than thicknesses this makes the thick slabs finish sooner.
# Example, 10x the statistics of the default run below 2 GeV:
#   sh scripts/lead_filter.sh --below 2 --parts 4 10 14000000

below=""; parts=1; dir=""
while [ $# -gt 0 ]; do
  case "$1" in
    --below) below=$2; shift 2 ;;
    --parts) parts=$2; shift 2 ;;
    --dir) dir=$2; shift 2 ;;
    *) break ;;
  esac
done
jobs=${1:-6}; [ $# -gt 0 ] && shift
muons=${1:-2000000}; [ $# -gt 0 ] && shift
thicknesses=${*:-5 10 20 30 40 50}
if [ ! -x ./TOFmeter ] || [ ! -f macro/lead_filter.mac ]; then
  echo "Run this in the build directory (TOFmeter and macro/lead_filter.mac must be there)" >&2
  exit 1
fi
if [ -z "$dir" ]; then
  if [ -n "$below" ]; then dir=output/lead_filter_below${below}GeV; else dir=output/lead_filter; fi
fi
mkdir -p "$dir/log"
perPart=$(( (muons + parts - 1) / parts ))
list=""
for t in $thicknesses; do
  if [ -n "$below" ]; then emax=$(awk -v e="$below" -v t="$t" 'BEGIN { printf "%.4g", e + 0.2 + 0.07 * t }')
  else emax=108.118; fi
  p=0
  while [ $p -lt "$parts" ]; do
    name=t$t; [ "$parts" -gt 1 ] && name=t${t}_p$p
    {
      echo "/control/verbose 2"
      echo "/control/alias leadT $t"
      echo "/control/alias leadFilterEvents $perPart"
      echo "/control/alias leadFilterSeed $((30000 + 100 * p + t))"
      echo "/control/alias leadFilterEmax $emax"
      echo "/control/alias leadFilterDir $dir"
      echo "/control/alias leadFilterFile lead_filter_$name"
      echo "/control/execute macro/lead_filter.mac"
    } > "$dir/log/$name.mac"
    list="$name $list"   # thickest first
    p=$((p + 1))
  done
done
echo "$(date '+%F %T')  lead slabs of $thicknesses cm, $muons muons each${below:+ (only those that can leave below $below GeV)}," \
     "in $parts process(es) per thickness, $jobs at a time -> $dir"
export TOF_FILTER_LOG="$dir/log"
echo $list | tr ' ' '\n' | grep . | xargs -P "$jobs" -n 1 sh -c '
  if ./TOFmeter "$TOF_FILTER_LOG/$1.mac" > "$TOF_FILTER_LOG/$1.log" 2>&1; then
    echo "$(date "+%F %T")  $1 done"
  else
    echo "$(date "+%F %T")  $1 FAILED, see $TOF_FILTER_LOG/$1.log"
  fi' sh
root -l -b -q "analysis/lead_filter_spectrum.C(\"$dir\", ${below:-0})"
