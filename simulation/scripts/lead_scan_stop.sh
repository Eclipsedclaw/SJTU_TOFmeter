#!/bin/sh
# Stop a lead x energy scan of scripts/lead_energy_scan.sh after the processes now running: no new
# process is started, and every running one finishes its macro (all its energy points) and leaves
# its .done marker. scripts/lead_scan_resume.sh later goes on exactly from there. Run it in the
# build directory:
#   sh scripts/lead_scan_stop.sh [definition]     (default macro/lead_thick_scan_def.mac)
# The computer is kept awake (on power) until the running processes have finished.

def=${1:-macro/lead_thick_scan_def.mac}
if [ ! -f "$def" ]; then
  echo "no $def: run this in the build directory" >&2
  exit 1
fi
dir=$(sed -n 's|^/control/alias leadScanDir *||p' "$def" | tr -d '"')

# the driver script (or a resume waiting to start it) and its xargs, which starts the processes;
# the TOFmeter processes themselves keep running
drivers=$(pgrep -f "lead_energy_scan.sh $def")
feeders=""
for d in $drivers; do feeders="$feeders $(pgrep -P "$d" xargs)"; done
if [ -z "$drivers" ]; then
  echo "no scan of $def is running or waiting to start"
else
  kill $drivers 2>/dev/null
  [ -n "$(echo $feeders)" ] && kill $feeders 2>/dev/null
  echo "$(date '+%F %T')  stopped: no new processes of $def will start" | tee -a "$dir/scan.log"
fi

# "[.]/TOFmeter" matches the processes but not the command line of this pattern itself
running=$(pgrep -f "[.]/TOFmeter $dir/log/")
if [ -n "$running" ]; then
  echo "still finishing, each with its last energy points (then its .done marker):"
  for p in $running; do
    j=$(ps -o command= -p "$p" | sed 's|.*/log/||; s|\.mac$||')
    echo "  $j: $(grep -c 'Run end' "$dir/log/$j.log") points finished, $(grep '^/tof/cosmic/monoEnergy' "$dir/log/$j.log" | tail -1 | sed 's|/tof/cosmic/monoEnergy|now at|')"
  done
  nohup caffeinate -is sh -c "while pgrep -qf '[.]/TOFmeter $dir/log/'; do sleep 60; done" > /dev/null 2>&1 &
fi
ndone=$(ls "$dir"/log/ 2>/dev/null | grep -c '\.done$')
total=$(cat "$dir/log/jobs.txt" 2>/dev/null | wc -l | tr -d ' ')
echo "$ndone of $total processes done (more as the running ones finish). Go on later with:"
echo "  sh scripts/lead_scan_resume.sh $def"
