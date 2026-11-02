#!/bin/sh
# Progress of the CZT run's HTCondor jobs (condor/muonic_xray_czt_submit.sh) in one folder: a bar over all its
# processes (those with a macro, <dir>/<band>/log/job<k>.mac), one per band with the exposure done, one per running
# job, the failed and stopped ones, what HTCondor reports (where condor_q runs) and the time left. It only reads
# files, so it also works on a copy of the folder. In the build directory:
#   sh condor/muonic_xray_czt_status.sh [--dir path] [--watch seconds]
# Muons done: those of the finished processes (job<k>.done) and, for the running ones, the last ">>> Event N" line of
# their log (every 1000 muons; the log is written in blocks, so it lags a little). A process runs while its log
# changes (within 30 min); one whose log stopped changing without a .done stopped (held or removed: condor_q -hold),
# one whose job<k>.err says FAILED failed. The time left: the muons still to run in each band at the speed of that
# band's jobs (the finished ones' from job<k>.out, else the running ones'), shared by the jobs running now.
# --watch 60 redraws every 60 s (Ctrl-C ends it).

dir=output/muonic_xray_czt
watch=0
while [ $# -gt 0 ]; do
  case "$1" in
    --dir) dir=$2; shift 2 ;;
    --watch) watch=$2; shift 2 ;;
    *) echo "unknown argument $1" >&2; exit 1 ;;
  esac
done
if [ ! -d "$dir/slow/log" ] && [ ! -d "$dir/fast/log" ]; then
  echo "no $dir/slow/log or $dir/fast/log: run this in the build directory, or give --dir" >&2
  exit 1
fi
name=tof_$(basename "$dir")

# one line per process: band, k, state (done, running, stopped, failed or waiting), its muons, the muons done, its
# running time in s (or 0) and, for a finished one, its live time in h
processes() {
  for band in slow fast; do
    log=$dir/$band/log
    [ -d "$log" ] || continue
    for mac in $(ls "$log"/job*.mac 2> /dev/null | sort -V); do
      k=${mac##*/job}
      k=${k%.mac}
      job=$log/job$k
      m=$(sed -n 's/^\/control\/alias czEvents //p' "$mac")
      n=0
      s=0
      h=0
      if [ -e "$job.done" ]; then
        state=done
        n=$m
        s=$(sed -n 's/.* done in \([0-9]*\) s .*/\1/p' "$job.out" 2> /dev/null)
        # "... = 2.88857 h of equivalent live time": the number before "h"
        h=$(tail -n 20 "$job.log" 2> /dev/null | grep 'h of equivalent live time' | awk '{ print $(NF-5) }')
      elif [ -e "$job.log" ]; then
        n=$(tail -n 20 "$job.log" | grep '^>>> Event' | tail -n 1 | awk '{ print $3 }')
        if grep -q FAILED "$job.err" 2> /dev/null; then
          state=failed
        elif [ $((now - $(stat -c %Y "$job.log"))) -lt 1800 ]; then
          state=running
          [ -e "$job.out" ] && s=$((now - $(stat -c %Y "$job.out")))
        else
          state=stopped
        fi
      else
        state=waiting
      fi
      echo "$band $k $state $m ${n:-0} ${s:-0} ${h:-0}"
    done
  done
}

status() {
  now=$(date +%s)
  q=
  if command -v condor_q > /dev/null; then
    q=$(condor_q -constraint "JobBatchName == \"$name\"" -af JobStatus 2> /dev/null |
        awk '{ c[$1]++ } END { printf "%d running, %d idle, %d held", c[2], c[1], c[5] }')
  fi
  # the report; its first line "LEFT <s>" gives the time left (-1: unknown) for the end's date, put in by the shell
  out=$(processes | awk -v date="$(date '+%F %T')" -v dir="$dir" -v name="$name" -v q="$q" '
    function bar(f, w,   i, n, s) { n = int(f * w + 0.5); s = ""; for (i = 0; i < w; ++i) s = s (i < n ? "█" : "░"); return s }
    function dur(s,   d, h, m) { d = int(s / 86400); h = int((s % 86400) / 3600); m = int((s % 3600) / 60)
                                 if (d > 0) return sprintf("%d d %d h", d, h); if (h > 0) return sprintf("%d h %d min", h, m)
                                 return sprintf("%d min", m) }
    function num(x) { return sprintf("%.3g", x) }
    {
      b = $1; st = $3; m = $4; n = $5; s = $6; h = $7
      procs[b]++; total[b] += m; all += m; count[st]++; bandcount[b, st]++
      if (st == "done" || st == "running") done[b] += n
      if (st == "done" && s > 0) { v[b] += m / s; nv[b]++ }
      if (st == "done" && h > 0) hours[b] = h / m
      if (st == "running" && s > 60 && n > 0) { vr[b] += n / s; nvr[b]++ }
      if (st == "running")
        run[++nrun] = sprintf("  %s %-4s [%s] %3.0f%%  %s of %s muons%s", b, $2, bar(n / m, 20), 100 * n / m, num(n), num(m),
                              s > 0 ? ", running " dur(s) : "")
      if (st == "failed" || st == "stopped")
        bad[++nbad] = sprintf("  %s %-4s %s at %s of %s muons: %s/%s/log/job%s.err and .log", b, $2, st, num(n), num(m), dir, b, $2)
    }
    END {
      if (NR == 0) { print "LEFT -1"; printf "%s  %s: no processes (no job macros)\n", date, dir; exit }
      split("slow fast", bands, " ")
      sum = 0; work = 0; known = 1
      for (i = 1; i <= 2; ++i) {
        b = bands[i]
        if (!procs[b]) continue
        sum += done[b]
        speed = nv[b] ? v[b] / nv[b] : (nvr[b] ? vr[b] / nvr[b] : 0)   # muons per second of one job
        if (total[b] > done[b]) { if (speed > 0) work += (total[b] - done[b]) / speed; else known = 0 }
      }
      nrunning = count["running"] + 0
      print "LEFT " (known && nrunning > 0 && work > 0 ? int(work / nrunning) : -1)
      printf "%s  %s (HTCondor batch %s)\n", date, dir, name
      printf "[%s] %5.1f%%   %d of %d processes done, %d running, %d waiting", bar(sum / all, 40), 100 * sum / all,
             count["done"], NR, nrunning, count["waiting"]
      if (count["failed"]) printf ", %d failed", count["failed"]
      if (count["stopped"]) printf ", %d stopped", count["stopped"]
      printf "\n"
      if (q != "") print "HTCondor: " q
      for (i = 1; i <= 2; ++i) {
        b = bands[i]
        if (!procs[b]) continue
        printf "  %s [%s] %5.1f%%  %d of %d done, %d running", b, bar(done[b] / total[b], 20), 100 * done[b] / total[b],
               bandcount[b, "done"], procs[b], bandcount[b, "running"]
        if (hours[b] > 0) printf "; %.1f of %.1f days of exposure", hours[b] * done[b] / 24, hours[b] * total[b] / 24
        printf "\n"
      }
      if (sum >= all) print "all done"
      else if (known && nrunning > 0 && work > 0) printf "about %s left (done around @WHEN@) with %d running\n", dur(work / nrunning), nrunning
      else if (nrunning == 0) print "no job running now (waiting for free slots, or not submitted: condor_q)"
      if (nrun) { print "running:"; for (i = 1; i <= nrun; ++i) print run[i] }
      if (nbad) { print "failed or stopped (condor_q -hold; condor_release, or submit again):"; for (i = 1; i <= nbad; ++i) print bad[i] }
    }')
  left=$(printf '%s\n' "$out" | sed -n '1s/^LEFT //p')
  when=
  [ "${left:--1}" -ge 0 ] && when=$(date -d "@$((now + left))" '+%a %d %b %H:%M')
  printf '%s\n' "$out" | sed -e 1d -e "s/@WHEN@/$when/"
}

if [ "$watch" -gt 0 ] 2> /dev/null; then
  while :; do
    printf '\033[H\033[2J'
    status
    sleep "$watch"
  done
else
  status
fi
