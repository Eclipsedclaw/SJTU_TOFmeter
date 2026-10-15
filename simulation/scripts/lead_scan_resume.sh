#!/bin/sh
# Go on with a lead x energy scan stopped by scripts/lead_scan_stop.sh, or interrupted in any other
# way: starts scripts/lead_energy_scan.sh in the background, at low priority and keeping the
# computer awake (on power). It skips every process that has its .done marker, so the scan goes on
# exactly where it stopped; processes still finishing from before are waited for first. Run it in
# the build directory:
#   sh scripts/lead_scan_resume.sh [definition] [jobs]   (default macro/lead_thick_scan_def.mac, 10)
# Progress is appended to <dir>/scan.log; stop again with scripts/lead_scan_stop.sh.

def=${1:-macro/lead_thick_scan_def.mac}
jobs=${2:-10}
if [ ! -f "$def" ] || [ ! -f scripts/lead_energy_scan.sh ]; then
  echo "no $def or scripts/lead_energy_scan.sh: run this in the build directory" >&2
  exit 1
fi
dir=$(sed -n 's|^/control/alias leadScanDir *||p' "$def" | tr -d '"')
if pgrep -qf "lead_energy_scan.sh $def"; then
  echo "a scan of $def is already running or waiting to start (stop it with scripts/lead_scan_stop.sh)" >&2
  exit 1
fi
mkdir -p "$dir"
# "[.]/TOFmeter" matches the processes but not the command line of this pattern itself
finishing=$(pgrep -f "[.]/TOFmeter $dir/log/" | wc -l | tr -d ' ')
nohup nice -n 10 caffeinate -is sh -c \
  "while pgrep -qf '[.]/TOFmeter $dir/log/'; do sleep 30; done; exec sh scripts/lead_energy_scan.sh $def $jobs" \
  >> "$dir/scan.log" 2>&1 &
if [ "$finishing" -gt 0 ]; then
  echo "$finishing processes of the earlier start are still finishing; the scan goes on when they are done"
else
  echo "the scan of $def goes on now, $jobs processes at a time"
fi
ndone=$(ls "$dir"/log/ 2>/dev/null | grep -c '\.done$')
total=$(cat "$dir/log/jobs.txt" 2>/dev/null | wc -l | tr -d ' ')
echo "$ndone of $total processes done; progress: tail -f $dir/scan.log (stop: sh scripts/lead_scan_stop.sh $def)"
