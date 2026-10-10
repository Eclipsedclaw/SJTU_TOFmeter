#!/bin/sh
# Lead-thickness x muon-energy scan, one TOFmeter process per lead thickness (0-20 cm), several at
# a time. Run it in the build directory:
#   sh scripts/lead_energy_scan.sh [definition] [jobs] [muons per point] [muons per check run]
# definition: macro/lead_energy_scan_def.mac (default, 0.1-100 GeV log grid) or
# macro/lead_lowE_scan_def.mac (100-500 MeV in 20 MeV steps); jobs: default 8; the muon numbers
# default to those of the definition. Writes <dir>/leadscan_t<cm>_E<GeV>.root (one file per
# point), leadspec_t<cm>.root (check run) and one macro and log per thickness in <dir>/log/,
# then runs analysis/lead_energy_scan.C on <dir> (leadScanDir of the definition).

def=macro/lead_energy_scan_def.mac
case "$1" in *.mac) def=$1; shift ;; esac
jobs=${1:-8}
events=$2
spectrum=$3
if [ ! -x ./TOFmeter ] || [ ! -f macro/lead_energy_scan_t.mac ] || [ ! -f "$def" ]; then
  echo "Run this in the build directory (TOFmeter, macro/lead_energy_scan_t.mac and $def must be there)" >&2
  exit 1
fi
dir=$(sed -n 's|^/control/alias leadScanDir *||p' "$def")
if [ -z "$dir" ]; then
  echo "$def defines no leadScanDir" >&2
  exit 1
fi
mkdir -p "$dir/log"

for t in $(seq 0 20); do
  {
    echo "/control/verbose 2"
    echo "/control/execute $def"
    [ -n "$events" ] && echo "/control/alias leadScanEvents $events"
    [ -n "$spectrum" ] && echo "/control/alias leadSpectrumEvents $spectrum"
    echo "/control/alias leadT $t"
    echo "/control/execute macro/lead_energy_scan_t.mac"
  } > "$dir/log/t$t.mac"
done

# thickest lead first: those processes take longest
echo "$(date '+%F %T')  $def -> $dir, $jobs jobs${events:+, $events muons per point}${spectrum:+, $spectrum per check run}"
export TOF_SCAN_DIR="$dir"
seq 20 -1 0 | xargs -P "$jobs" -n 1 sh -c '
  if ./TOFmeter "$TOF_SCAN_DIR/log/t$1.mac" > "$TOF_SCAN_DIR/log/t$1.log" 2>&1; then
    echo "$(date "+%F %T")  lead $1 cm done"
  else
    echo "$(date "+%F %T")  lead $1 cm FAILED, see $TOF_SCAN_DIR/log/t$1.log"
  fi' sh

root -l -b -q "analysis/lead_energy_scan.C(\"$dir\")"
