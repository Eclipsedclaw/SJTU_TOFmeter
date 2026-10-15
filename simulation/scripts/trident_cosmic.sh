#!/bin/sh
# Rate of muon tridents in the TOF-meter: several TOFmeter jobs of macro/trident_cosmic.mac with
# different seeds, then analysis/trident_rate.C. Run it in the build directory:
#   sh scripts/trident_cosmic.sh [jobs, default 8] [muons per job, default 1000000] [lead cm, default 10]
# Writes output/tridents_t<lead>/trident_cosmic_<job>.root (earlier files of this script there are
# replaced) and one macro and log per job in output/tridents_t<lead>/log/.

jobs=${1:-8}
events=${2:-1000000}
lead=${3:-10}
if [ ! -x ./TOFmeter ] || [ ! -f macro/trident_cosmic.mac ] || [ ! -f macro/physics_tridents.mac ]; then
  echo "Run this in the build directory (TOFmeter and macro/ must be there)" >&2
  exit 1
fi
dir=output/tridents_t$lead
mkdir -p "$dir/log"
rm -f "$dir"/trident_cosmic_*.root

for k in $(seq 1 "$jobs"); do
  {
    echo "/control/verbose 2"
    echo "/control/alias tridentSeed $((7000 + 100 * lead + k))"
    echo "/control/alias tridentEvents $events"
    echo "/control/alias tridentDir $dir"
    echo "/control/alias tridentFile trident_cosmic_$k"
    echo "/control/alias tridentLead $lead"
    echo "/control/execute macro/trident_cosmic.mac"
  } > "$dir/log/job$k.mac"
done

echo "$(date '+%F %T')  $jobs jobs x $events muons, lead $lead cm -> $dir"
export TOF_TRIDENT_LOG="$dir/log/job"
seq 1 "$jobs" | xargs -P "$jobs" -n 1 sh -c '
  if ./TOFmeter --preinit macro/physics_tridents.mac "$TOF_TRIDENT_LOG$1.mac" > "$TOF_TRIDENT_LOG$1.log" 2>&1; then
    echo "$(date "+%F %T")  job $1 done"
  else
    echo "$(date "+%F %T")  job $1 FAILED, see $TOF_TRIDENT_LOG$1.log"
  fi' sh

root -l -b -q "analysis/trident_rate.C(\"$dir\")"
