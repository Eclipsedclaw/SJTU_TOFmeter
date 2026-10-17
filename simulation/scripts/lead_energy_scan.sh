#!/bin/sh
# Lead-thickness x muon-energy scan, several TOFmeter processes at a time. Run it in the build
# directory:
#   sh scripts/lead_energy_scan.sh [--check-only] [definition] [jobs] [muons per point] [muons per check run]
# definition: macro/lead_energy_scan_def.mac (default, 0.1-100 GeV log grid, lead 0-20 cm),
# macro/lead_lowE_scan_def.mac (100-500 MeV in 20 MeV steps) or macro/lead_thick_scan_def.mac
# (50 energies, lead 0-50 cm); jobs: default 8; the muon numbers default to those of the
# definition. The definition gives the thicknesses (0 to leadScanTmax cm) and leadScanGroups: each
# thickness runs as that many processes, each with a part of the energies (and its own seed when
# there are several), plus one process for its spectrum check run. Writes
# <dir>/leadscan_t<cm>_E<GeV>.root (one file per point), leadspec_t<cm>.root (check run) and a
# macro and log per process in <dir>/log/, then runs analysis/lead_energy_scan.C and
# analysis/lead_after_spectrum.C on <dir> (leadScanDir of the definition).
# A process that finishes leaves <dir>/log/<name>.done, and those are skipped when the script is
# started again: an interrupted scan goes on where it stopped (delete the .done files to redo).
# --check-only repeats only the spectrum check runs (same seeds, so the same events), e.g. after
# an output column was added.
# Long scans: keep the computer awake and on power, e.g. caffeinate -is sh scripts/lead_energy_scan.sh ...
# or start it with scripts/lead_scan_resume.sh, which does that in the background; stop it after
# the running processes with scripts/lead_scan_stop.sh, and go on later with lead_scan_resume.sh.

check_only=false
if [ "$1" = "--check-only" ]; then check_only=true; shift; fi
def=macro/lead_energy_scan_def.mac
case "$1" in *.mac) def=$1; shift ;; esac
jobs=${1:-8}
events=$2
spectrum=$3
for f in ./TOFmeter macro/lead_energy_scan_points_t.mac macro/lead_energy_scan_check_t.mac "$def"; do
  if [ ! -e "$f" ]; then
    echo "Run this in the build directory ($f is missing)" >&2
    exit 1
  fi
done
alias_value() { sed -n "s|^/control/alias $1 *||p" "$def" | tr -d '"'; }
dir=$(alias_value leadScanDir)
if [ -z "$dir" ]; then
  echo "$def defines no leadScanDir" >&2
  exit 1
fi
tmax=$(alias_value leadScanTmax); tmax=${tmax:-20}
groups=$(alias_value leadScanGroups); groups=${groups:-1}
seed0=$(alias_value leadScanSeed0)
energies=$(alias_value leadScanEnergies)
nE=$(echo $energies | wc -w | tr -d ' ')
size=$(( (nE + groups - 1) / groups ))
logdir="$dir/log"
mkdir -p "$logdir"

# one macro per process: <name>.mac, names t<cm> (or t<cm>_g<group>) and check_t<cm>
header() {
  echo "/control/verbose 2"
  echo "/control/execute $def"
  [ -n "$events" ] && echo "/control/alias leadScanEvents $events"
  [ -n "$spectrum" ] && echo "/control/alias leadSpectrumEvents $spectrum"
  echo "/control/alias leadT $1"
}
list="$logdir/jobs.txt"
: > "$list"
if ! $check_only; then
  # most expensive first: highest energies, thickest lead
  g=$((groups - 1))
  while [ $g -ge 0 ]; do
    for t in $(seq "$tmax" -1 0); do
      name=t$t
      [ "$groups" -gt 1 ] && name=t${t}_g$g
      {
        header "$t"
        if [ "$groups" -gt 1 ]; then
          part=$(echo $energies | tr ' ' '\n' | sed -n "$((g * size + 1)),$(((g + 1) * size))p" | tr '\n' ' ')
          echo "/control/alias leadScanEnergies \"${part% }\""
          echo "/control/alias leadScanSeed0 $((seed0 + 1000 * g))"
        fi
        echo "/control/execute macro/lead_energy_scan_points_t.mac"
      } > "$logdir/$name.mac"
      echo "$name" >> "$list"
    done
    g=$((g - 1))
  done
fi
for t in $(seq "$tmax" -1 0); do
  $check_only && rm -f "$logdir/check_t$t.done"
  { header "$t"; echo "/control/execute macro/lead_energy_scan_check_t.mac"; } > "$logdir/check_t$t.mac"
  echo "check_t$t" >> "$list"
done

todo=$(while read -r name; do [ -e "$logdir/$name.done" ] || echo "$name"; done < "$list")
total=$(wc -l < "$list" | tr -d ' ')
left=$(echo $todo | wc -w | tr -d ' ')
echo "$(date '+%F %T')  $def -> $dir: lead 0-$tmax cm, $nE energies in $groups group(s), $jobs jobs," \
     "$left of $total processes to run${events:+, $events muons per point}${spectrum:+, $spectrum per check run}"
export TOF_SCAN_LOG="$logdir"
echo $todo | tr ' ' '\n' | grep . | xargs -P "$jobs" -n 1 sh -c '
  if pgrep -qf "[.]/TOFmeter $TOF_SCAN_LOG/$1[.]mac"; then
    echo "$(date "+%F %T")  $1 is still running from an earlier start; not started again"
  elif ./TOFmeter "$TOF_SCAN_LOG/$1.mac" > "$TOF_SCAN_LOG/$1.log" 2>&1; then
    touch "$TOF_SCAN_LOG/$1.done"
    echo "$(date "+%F %T")  $1 done"
  else
    echo "$(date "+%F %T")  $1 FAILED, see $TOF_SCAN_LOG/$1.log"
  fi' sh

missing=$(while read -r name; do [ -e "$logdir/$name.done" ] || echo "$name"; done < "$list" | wc -l | tr -d ' ')
echo "$(date '+%F %T')  finished; $missing of $total processes not done"
root -l -b -q "analysis/lead_energy_scan.C(\"$dir\")"
root -l -b -q "analysis/lead_after_spectrum.C(\"$dir\")"
