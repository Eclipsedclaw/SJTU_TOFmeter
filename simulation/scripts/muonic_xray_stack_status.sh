#!/bin/sh
# Progress of scripts/muonic_xray_stack.sh: a bar over all its muons, one per running process, and the time left.
# It only reads files. Run it in the build directory:
#   sh scripts/muonic_xray_stack_status.sh [--dir path] [--watch seconds]
# Muons done: those of the finished processes (<dir>/log/job<k>.done) plus, for the running ones, the last
# ">>> Event N" line of their log (printed every 1000 events). The time left comes from the speed now: the muons of
# each running process over its running time (ps), added up, so a change in the number at a time shows within
# minutes (without running processes: the processes finished in the last 6 h, or the average since <dir>/started).
# The exposure comes from the live time a finished process printed.
# --watch 60 redraws every 60 s (Ctrl-C ends it).

dir=output/muonic_xray_stack; watch=0
while [ $# -gt 0 ]; do
  case "$1" in
    --dir) dir=$2; shift 2 ;;
    --watch) watch=$2; shift 2 ;;
    *) echo "unknown argument $1" >&2; exit 1 ;;
  esac
done
if [ ! -d "$dir/log" ]; then
  echo "no $dir/log: run this in the build directory" >&2
  exit 1
fi

bar() {   # fraction, width
  awk -v f="$1" -v w="$2" 'BEGIN { n = int(f * w + 0.5); s = ""; for (i = 0; i < w; ++i) s = s (i < n ? "█" : "░"); printf "%s", s }'
}

duration() {   # seconds -> "1 d 13 h", "5 h 10 min", "12 min"
  awk -v s="$1" 'BEGIN { d = int(s / 86400); h = int((s % 86400) / 3600); m = int((s % 3600) / 60);
                         if (d > 0) printf "%d d %d h", d, h; else if (h > 0) printf "%d h %d min", h, m; else printf "%d min", m }'
}

status() {
  jobs=$(ls "$dir/log" | grep -c '^job[0-9]*\.mac$')
  events=$(sed -n 's/^\/control\/alias stkEvents //p' "$dir/log/job1.mac")
  finished=$(ls "$dir/log" | grep -c '\.done$')
  done_events=$((finished * events))
  running=0; stopped=0; lines=""; speed=0   # muons per second of the running processes
  for log in "$dir"/log/job*.log; do
    [ -e "$log" ] || continue
    k=${log##*/job}; k=${k%.log}
    [ -e "$dir/log/job$k.done" ] && continue
    n=$(grep '^>>> Event' "$log" | tail -1 | awk '{ print $3 }')
    n=${n:-0}
    done_events=$((done_events + n))
    f=$(awk -v n="$n" -v e="$events" 'BEGIN { printf "%.4f", n / e }')
    pid=$(pgrep -f "^\./TOFmeter $dir/log/job$k\.mac\$" | head -1)
    if [ -n "$pid" ]; then
      running=$((running + 1))
      seconds=$(ps -o etimes= -p "$pid" | tr -d ' ')
      speed=$(awk -v s="$speed" -v n="$n" -v t="${seconds:-0}" 'BEGIN { print s + (t > 0 ? n / t : 0) }')
      lines="$lines$(printf '  job %-4s [%s] %3.0f%%' "$k" "$(bar "$f" 20)" "$(awk -v f="$f" 'BEGIN { print 100 * f }')")\n"
    else
      stopped=$((stopped + 1))
      lines="$lines$(printf '  job %-4s stopped at %s muons without finishing (see log/job%s.log)' "$k" "$n" "$k")\n"
    fi
  done
  total=$((jobs * events))
  fraction=$(awk -v d="$done_events" -v t="$total" 'BEGIN { printf "%.4f", d / t }')
  # start of the run: <dir>/started (written by the driver), else the oldest job macro
  if [ -e "$dir/started" ]; then start=$(cat "$dir/started"); else start=$(stat -c %Y "$(ls -tr "$dir"/log/job*.mac | head -1)"); fi
  now=$(date +%s)
  elapsed=$((now - start))
  # speed: the processes finished in the last 6 h if there are enough of them, else the average since the start
  recent=$(find "$dir/log" -name 'job*.done' -mmin -360 | wc -l)
  printf '%s  %s\n' "$(date '+%F %T')" "$dir"
  printf '[%s] %5.1f%%   %d of %d processes done, %d running' "$(bar "$fraction" 40)" \
         "$(awk -v f="$fraction" 'BEGIN { print 100 * f }')" "$finished" "$jobs" "$running"
  [ "$stopped" -gt 0 ] && printf ', %d stopped' "$stopped"
  [ -e "$dir/STOP" ] && printf ' (STOP file: no new process starts)'
  printf '\n'
  printf '%s of %s muons' "$(awk -v n="$done_events" 'BEGIN { printf "%.3g", n }')" "$(awk -v n="$total" 'BEGIN { printf "%.3g", n }')"
  # exposure: hours of live time per process, from a finished one
  # "... = 1.3077 h of equivalent live time": the number before "h"
  hours=$(grep -h 'h of equivalent live time' "$dir"/log/job*.log 2>/dev/null | head -1 | awk '{ print $(NF-5) }')
  if [ -n "$hours" ]; then
    awk -v h="$hours" -v d="$done_events" -v t="$total" -v e="$events" \
      'BEGIN { printf "; %.1f of %.1f days of sea-level exposure", h * d / e / 24, h * t / e / 24 }'
  fi
  printf '\n'
  if [ "$done_events" -gt 0 ] && [ "$elapsed" -gt 0 ]; then
    if awk -v s="$speed" 'BEGIN { exit !(s > 0) }'; then
      left=$(awk -v s="$speed" -v d="$done_events" -v t="$total" 'BEGIN { printf "%d", (t - d) / s }')
      now_speed=$(awk -v s="$speed" -v r="$running" 'BEGIN { printf "; now %.0f muons/s with %d running", s, r }')
    elif [ "$recent" -ge 12 ]; then
      left=$(awk -v r="$recent" -v e="$events" -v d="$done_events" -v t="$total" 'BEGIN { printf "%d", (t - d) / (r * e / 21600) }')
      now_speed=""
    else
      left=$(awk -v d="$done_events" -v t="$total" -v s="$elapsed" 'BEGIN { printf "%d", (t - d) * s / d }')
      now_speed=""
    fi
    printf 'running for %s%s; about %s left (done around %s)\n' "$(duration "$elapsed")" "$now_speed" \
           "$(duration "$left")" "$(date -d "@$((now + left))" '+%a %d %b %H:%M')"
  fi
  [ -n "$lines" ] && printf '%b' "$lines"   # %b: the \n become new lines, a % stays a %
}

if [ "$watch" -gt 0 ] 2>/dev/null; then
  while :; do
    printf '\033[H\033[2J'
    status
    sleep "$watch"
  done
else
  status
fi
