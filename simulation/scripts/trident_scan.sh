#!/bin/sh
# Muon-trident rate versus lead thickness: many TOFmeter processes of macro/trident_cosmic.mac, several
# at a time, then analysis/trident_rate.C for each thickness. Run it in the build directory:
#   sh scripts/trident_scan.sh [definition] [jobs at a time]   (default macro/trident_scan_def.mac, 10)
# The definition gives the output directory, the thicknesses (run in their order), the processes per
# thickness, the muons per process, the bias factor F and the source-plane height (see
# macro/trident_scan_def.mac). Writes <dir>/t<cm>/trident_cosmic_<k>.root with a macro and log per
# process in <dir>/t<cm>/log/; when the last process of a thickness is done, analysis/trident_rate.C
# makes <dir>/t<cm>/trident_rate.png/.pdf/.txt; at the end <dir>/summary.txt has the rate tables.
# A process that finishes leaves <dir>/t<cm>/log/job<k>.done, and those are skipped when the script
# is started again, so an interrupted scan goes on where it stopped (delete the .done files to redo).
# touch <dir>/STOP: no new process starts (those running finish and leave their .done); delete
# STOP and start the script again to go on. Long scans: start it in the background, at low priority
# and keeping the computer awake (on power):
#   nohup nice -n 15 caffeinate -is sh scripts/trident_scan.sh macro/trident_scan_def.mac 10 >> <dir>/scan.log 2>&1 &

def=${1:-macro/trident_scan_def.mac}
jobs=${2:-10}
for f in ./TOFmeter macro/trident_cosmic.mac macro/physics_tridents.mac analysis/trident_rate.C "$def"; do
  if [ ! -e "$f" ]; then
    echo "Run this in the build directory ($f is missing)" >&2
    exit 1
  fi
done
alias_value() { sed -n "s|^/control/alias $1 *||p" "$def" | tr -d '"'; }
dir=$(alias_value tridentScanDir)
leads=$(alias_value tridentScanLead)
perLead=$(alias_value tridentScanJobs)
events=$(alias_value tridentEvents)
bias=$(alias_value tridentBias)
planeZ=$(alias_value tridentPlaneZ)
seed0=$(alias_value tridentScanSeed0)
for v in dir leads perLead events bias planeZ seed0; do
  if [ -z "$(eval echo \$$v)" ]; then
    echo "$def does not define all of tridentScanDir tridentScanLead tridentScanJobs tridentEvents tridentBias tridentPlaneZ tridentScanSeed0" >&2
    exit 1
  fi
done
# the muons start on the plane: it must be above the lead layer (bottom at 21 cm)
for t in $leads; do
  if ! awk -v t="$t" -v z="$planeZ" 'BEGIN { exit !(t + 21 < z) }'; then
    echo "lead $t cm reaches z = $(awk -v t="$t" 'BEGIN { print t + 21 }') cm, not below the source plane at $planeZ cm" >&2
    exit 1
  fi
done
mkdir -p "$dir"
if ! mkdir "$dir/driver.lock" 2>/dev/null; then
  echo "$dir/driver.lock exists: a scan of $dir is running (or one was killed: then delete the lock)" >&2
  exit 1
fi
trap 'rmdir "$dir/driver.lock" 2>/dev/null' EXIT
trap 'exit 1' INT TERM

# one macro per process, t<cm>/log/job<k>.mac
list="$dir/jobs.txt"
: > "$list"
for t in $leads; do
  mkdir -p "$dir/t$t/log"
  for k in $(seq 1 "$perLead"); do
    {
      echo "/control/verbose 2"
      echo "/control/execute $def"
      echo "/control/alias tridentSeed $((seed0 + 1000 * t + k))"
      echo "/control/alias tridentDir $dir/t$t"
      echo "/control/alias tridentFile trident_cosmic_$k"
      echo "/control/alias tridentLead $t"
      echo "/control/execute macro/trident_cosmic.mac"
    } > "$dir/t$t/log/job$k.mac"
    echo "$t:$k" >> "$list"
  done
done
total=$(wc -l < "$list" | tr -d ' ')
todo=$(while IFS=: read -r t k; do [ -e "$dir/t$t/log/job$k.done" ] || echo "$t:$k"; done < "$list")
left=$(echo $todo | wc -w | tr -d ' ')
echo "$(date '+%F %T')  $def -> $dir: lead $leads cm, $perLead x $events muons per thickness, F = $bias," \
     "plane at z = $planeZ cm, $jobs jobs at a time; $left of $total processes to run"
[ -e "$dir/STOP" ] && echo "$dir/STOP exists: nothing will start until it is deleted"

export TOF_SCAN_DIR="$dir" TOF_SCAN_PER_LEAD="$perLead" TOF_SCAN_TOTAL="$total"
echo $todo | tr ' ' '\n' | grep . | xargs -P "$jobs" -n 1 sh -c '
  t=${1%:*}; k=${1#*:}
  d=$TOF_SCAN_DIR/t$t
  mac=$d/log/job$k.mac
  [ -e "$TOF_SCAN_DIR/STOP" ] && exit 0
  [ -e "$d/log/job$k.done" ] && exit 0
  if pgrep -f "[.]/TOFmeter --preinit macro/physics_tridents.mac $mac" > /dev/null; then
    echo "$(date "+%F %T")  t$t job $k is still running from an earlier start; not started again"
    exit 0
  fi
  if ./TOFmeter --preinit macro/physics_tridents.mac "$mac" > "$d/log/job$k.log" 2>&1; then
    touch "$d/log/job$k.done"
    ndone=$(ls "$TOF_SCAN_DIR"/t*/log/ | grep -c "\.done$")
    echo "$(date "+%F %T")  t$t job $k done ($ndone of $TOF_SCAN_TOTAL)"
  else
    echo "$(date "+%F %T")  t$t job $k FAILED, see $d/log/job$k.log"
    exit 0
  fi
  # the last process of this thickness: its rate plot and table
  # (without ROOT, e.g. on another computer, analyse later where ROOT is)
  if [ "$(ls "$d/log/" | grep -c "\.done$")" -ge "$TOF_SCAN_PER_LEAD" ] && command -v root > /dev/null \
     && mkdir "$d/log/analysis.lock" 2>/dev/null; then
    root -l -b -q "analysis/trident_rate.C(\"$d\")" > "$d/trident_rate.txt" 2>&1
    rmdir "$d/log/analysis.lock"
    echo "$(date "+%F %T")  t$t complete: $d/trident_rate.png, $d/trident_rate.txt"
  fi' sh

missing=$(while IFS=: read -r t k; do [ -e "$dir/t$t/log/job$k.done" ] || echo x; done < "$list" | wc -l | tr -d ' ')
note=""
[ -e "$dir/STOP" ] && note=" (stopped by $dir/STOP)"
echo "$(date '+%F %T')  finished; $missing of $total processes not done$note"
for t in $leads; do
  [ -f "$dir/t$t/trident_rate.txt" ] || continue
  echo "=== lead $t cm"
  grep -E "^[0-9]+ cosmic muons|^bias factor|^all tridents|^  in the|^  elsewhere|^detectable|^  \.\.\.|^median" "$dir/t$t/trident_rate.txt"
done > "$dir/summary.txt"
echo "rate tables of the finished thicknesses: $dir/summary.txt"
