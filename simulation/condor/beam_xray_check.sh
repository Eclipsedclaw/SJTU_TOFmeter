#!/bin/sh
# The state of every HTCondor job of a beam-study run (condor/beam_xray_submit.sh) in one folder, checked rather than
# assumed: per job and in total queued, running, held, finished OK and failed, and the progress (e.g. 37/200 finished OK).
# A job finished OK when all of these hold:
#   - HTCondor's log of its submission (<dir>/condor/submit_<date>.log) says it ended with exit code 0;
#   - its output <dir>/beam/beam_<k>.root exists and is not empty;
#   - TOFmeter's log <dir>/beam/log/job<k>.log reached the run's end ("Number of events processed : <muons>") without an
#     error (a Geant4 error, a crash, an abort; Geant4's warnings do not count).
# For each failed or held job it prints the HTCondor ID, the reason and the log paths. At the end it prints the commands
# that would run them again. It only reads: it never submits, releases, removes or edits anything. Run it where you
# submitted (condor_q there tells queued, running and held), in the build directory:
#   sh condor/beam_xray_check.sh <output folder | HTCondor cluster id> [--all]
# e.g. output/beam_xray_lyso_si_al0.01, or the cluster id condor_submit printed. The jobs of each state are listed as
# ranges of process numbers; --all lists every job on its own line. A job macro never run shows as not submitted.

all=no
arg=
for a in "$@"; do
  case "$a" in
    --all) all=yes ;;
    -*) echo "unknown option $a" >&2; exit 1 ;;
    *) arg=$a ;;
  esac
done
if [ -z "$arg" ]; then
  echo "usage: sh condor/beam_xray_check.sh <output folder | HTCondor cluster id> [--all]" >&2
  exit 1
fi
if [ -d "$arg" ]; then
  dir=${arg%/}
else
  case "$arg" in
    *[!0-9]*) echo "$arg: neither a folder nor an HTCondor cluster id (run this in the build directory)" >&2; exit 1 ;;
  esac
  log=$(grep -l -E "^[0-9]{3} \(0*$arg\.[0-9]+\.[0-9]+\)" output/*/condor/submit_*.log 2> /dev/null | head -n 1)
  if [ -z "$log" ]; then
    echo "no submission of cluster $arg in output/*/condor/submit_*.log: give its output folder" >&2
    exit 1
  fi
  dir=${log%/condor/*}
fi
if [ ! -d "$dir/beam/log" ]; then
  echo "$dir: no beam/log (run this in the build directory, with the run's output folder)" >&2
  exit 1
fi
name=tof_$(basename "$dir")
first=$(ls "$dir"/beam/log/job*.mac 2> /dev/null | head -n 1)
tgt=; imp=; frac=; ekin=
if [ -n "$first" ]; then
  tgt=$(sed -n 's/^\/tof\/target\/material //p' "$first")
  imp=$(sed -n 's/^\/tof\/target\/impurity //p' "$first")
  frac=$(sed -n 's/^\/tof\/target\/impurityFraction //p' "$first")
  ekin=$(sed -n 's/^\/control\/alias bmEkin //p' "$first")
fi
setupline="target ${tgt:-?}${imp:+ with $imp at atom fraction $frac}, mu- of ${ekin:-?} MeV, "

# HTCondor's view of the jobs still in the queue (and of submissions whose jobs are not all created yet)
qok=no
queue=
factories=
if command -v condor_q > /dev/null 2>&1; then
  if queue=$(condor_q -constraint "JobBatchName == \"$name\"" -af ClusterId ProcId JobStatus HoldReason 2> /dev/null); then
    qok=yes
    factories=$(condor_q -factory -constraint "JobBatchName == \"$name\"" -af ClusterId 2> /dev/null)
  fi
fi

# the command that submits again: the latest submission's options
latest=$(ls "$dir"/condor/submit_*.sub 2> /dev/null | sort | tail -n 1)
resub=
[ -n "$latest" ] && resub=$(sed -n 's/^# resubmit: //p' "$latest")

# TOFmeter's log: an error, not a warning (Geant4's errors print EEEE, its warnings WWWW)
errors='EEEE|FatalException|\*\*\* Break \*\*\*|[Ss]egmentation (violation|fault)|terminate called|core dumped|bad_alloc|Error in <|^ERROR|Aborted|Killed'
{
  # each job macro with the checks of its files: .done, output not empty, the run's end in the log, the first error
  jl=$dir/beam/log
  for mac in "$jl"/job*.mac; do
    [ -e "$mac" ] || continue
    k=${mac##*/job}
    k=${k%.mac}
    job=$jl/job$k
    m=$(sed -n 's/^\/control\/alias bmEvents //p' "$mac")
    d=0; [ -e "$job.done" ] && d=1
    s=0; [ -s "$dir/beam/beam_$k.root" ] && s=1
    e=0; [ -e "$job.log" ] && grep -q "Number of events processed : $m\$" "$job.log" && e=1
    err=
    [ -e "$job.log" ] && err=$(grep -m 1 -E "$errors" "$job.log")
    [ -z "$err" ] && [ -s "$job.err" ] && err=$(grep -m 1 -E 'FAILED|cannot read the setup file|no \./TOFmeter' "$job.err")
    printf 'M\tbeam\t%s\t%s\t%s\t%s\t%s\t%s\n' "$k" "$m" "$d" "$s" "$e" "$(printf '%s' "$err" | tr '\t' ' ' | cut -c1-160)"
  done
  # the submissions in time order: their jobs (ProcId = place in the queue list) and HTCondor's log of them
  for sub in $(ls "$dir"/condor/submit_*.sub 2> /dev/null | sort); do
    stamp=${sub##*/submit_}
    stamp=${stamp%.sub}
    awk -v st="$stamp" '/^queue band, k from \(/ { on = 1; next } on && /^\)/ { on = 0 }
                        on && NF == 2 { printf "S\t%s\t%d\t%s\t%s\n", st, i++, $1, $2 }' "$sub"
    [ -e "$dir/condor/submit_$stamp.log" ] &&
      awk -v st="$stamp" '{ gsub(/\t/, " "); printf "L\t%s\t%s\n", st, $0 }' "$dir/condor/submit_$stamp.log"
  done
  printf '%s\n' "$queue" | awk 'NF >= 3 { r = ""; for (i = 4; i <= NF; ++i) r = r (i > 4 ? " " : "") $i
                                          printf "Q\t%s\t%s\t%s\t%s\n", $1, $2, $3, r }'
  printf '%s\n' "$factories" | awk 'NF { printf "A\t%s\n", $1 }'
} | awk -F '\t' -v dir="$dir" -v name="$name" -v setupline="$setupline" -v qok="$qok" -v all="$all" -v resub="$resub" \
        -v date="$(date '+%F %T')" '
  function add(s, t) { return s == "" ? t : s "; " t }
  function bar(f, w,   i, n, s) { n = int(f * w + 0.5); s = ""; for (i = 0; i < w; ++i) s = s (i < n ? "#" : "."); return s }
  $1 == "M" { j = $3 + 0; job[j] = 1; muons[j] = $4; done[j] = $5; out[j] = $6; end[j] = $7; err[j] = $8
              if (j > maxk) maxk = j; next }
  $1 == "S" { j = $5 + 0; sstamp[j] = $2; sidx[j] = $3; next }   # a later submission replaces an earlier one
  $1 == "L" {
    line = $3
    if (pend != "") {   # the line after an end (its exit code) or a hold (its reason)
      if (pendcode == "005") {
        if (match(line, /return value [0-9]+/)) ev_exit[pend] = substr(line, RSTART + 13, RLENGTH - 13) + 0
        else if (match(line, /signal [0-9]+/)) ev_exit[pend] = substr(line, RSTART, RLENGTH)
      } else { r = line; sub(/^ +/, "", r); ev_hold[pend] = r }
      pend = ""
    }
    if (match(line, /^[0-9][0-9][0-9] \([0-9]+\.[0-9]+\.[0-9]+\)/)) {
      code = substr(line, 1, 3)
      split(substr(line, 6, RLENGTH - 6), p, ".")
      key = (p[1] + 0) "." (p[2] + 0)
      stampcl[$2] = p[1] + 0
      ev_last[key] = code
      if (code == "005" || code == "012") { pend = key; pendcode = code }
    }
    next
  }
  $1 == "Q" { key = ($2 + 0) "." ($3 + 0); q_st[key] = $4 + 0; q_why[key] = $5; alive[$2 + 0] = 1; next }
  $1 == "A" { alive[$2 + 0] = 1; next }
  END {
    split("ok running queued held failed unknown notsub", order, " ")
    label["ok"] = "finished OK"; label["running"] = "running"; label["queued"] = "queued"; label["held"] = "held"
    label["failed"] = "failed"; label["unknown"] = "unknown"; label["notsub"] = "not submitted"
    total = 0
    for (j in job) {
      ++total
      k = j
      files = ""
      if (!out[j]) files = add(files, "beam_" k ".root missing or empty")
      if (!end[j]) files = add(files, "TOFmeter'\''s log lacks the run'\''s end")
      if (err[j] != "") files = add(files, "error in the log: " err[j])
      st = ""; why = ""; hid = "-"
      if (j in sstamp) {
        if (!(sstamp[j] in stampcl)) {   # no HTCondor log of its submission (yet)
          if (qok == "yes") { st = "failed"; why = "no HTCondor log of its submission (" dir "/condor/submit_" sstamp[j] ".log)" }
          else st = "unknown"
        } else {
          key = stampcl[sstamp[j]] "." sidx[j]; hid = key
          if (key in q_st) {
            s = q_st[key]
            if (s == 1) st = "queued"
            else if (s == 2 || s == 6 || s == 7) st = "running"
            else if (s == 5) { st = "held"; why = q_why[key] }
            else if (s == 3) { st = "failed"; why = "being removed from the queue" }
          }
          if (st == "") {
            if (key in ev_last) {
              c = ev_last[key]
              if (c == "005") {
                if (ev_exit[key] "" == "0") { if (files == "") st = "ok"; else { st = "failed"; why = "exit code 0, but " files } }
                else { st = "failed"; why = "exit code " ev_exit[key] (err[j] != "" ? "; " err[j] : "") }
              } else if (c == "009") { st = "failed"; why = "removed from the queue (condor_rm)" }
              else if (c == "012") {
                if (qok == "yes") { st = "failed"; why = "held (" ev_hold[key] "), then removed from the queue" }
                else { st = "held"; why = ev_hold[key] }
              } else if (qok != "yes") st = "unknown"
              else if (files == "" && done[j]) { st = "ok"; why = "no end in HTCondor'\''s log; exit code 0 by its .done marker" }
              else { st = "failed"; why = "left the queue without an end in HTCondor'\''s log" (files != "" ? "; " files : "") }
            } else if (qok != "yes") st = "unknown"
            else if (alive[stampcl[sstamp[j]]]) st = "queued"   # its submission is still in the queue: not started yet
            else { st = "failed"; why = "its submission left the queue before it started" }
          }
        }
      } else if (done[j] && files == "") { st = "ok"; why = "run outside HTCondor: exit code 0 by its .done marker" }
      else if (done[j]) { st = "failed"; why = "run outside HTCondor with a .done marker, but " files }
      else st = "notsub"
      state[j] = st; reason[j] = why; htid[j] = hid
      count[st]++
      if (st == "ok") okmuons += muons[j]
    }
    if (total == 0) { printf "%s  %s: no job macros (%s/beam/log/job<k>.mac)\n", date, dir, dir; exit 1 }
    printf "%s  %s: %sHTCondor batch %s%s\n", date, dir, setupline, name, qok == "yes" ? "" : " (no condor_q here: queued, running and held unknown)"
    printf "progress [%s] %d/%d finished OK (%.1f%%), %.4g muons simulated\n", bar(count["ok"] / total, 40), count["ok"], total,
           100 * count["ok"] / total, okmuons
    line = ""
    for (i = 1; i <= 7; ++i) if (count[order[i]]) line = line (line == "" ? "" : ", ") count[order[i]] " " label[order[i]]
    printf "totals: %s\n", line
    # per job: the process numbers of each state, as ranges
    printf "per job (process numbers):\n"
    for (i = 1; i <= 7; ++i) {
      s = order[i]; o = ""; start = -1
      for (k = 1; k <= maxk + 1; ++k) {
        hit = (k in state) && state[k] == s
        if (hit && start < 0) start = k
        if (!hit && start >= 0) { o = o (o == "" ? "" : ", ") (start == k - 1 ? start : start "-" (k - 1)); start = -1 }
      }
      if (o != "") printf "  %-14s %s\n", label[s], o
    }
    if (all == "yes") {
      printf "every job:\n"
      for (k = 1; k <= maxk; ++k) {
        if (!(k in state)) continue
        printf "  %4d  %-14s job %-12s %s\n", k, label[state[k]], htid[k], reason[k]
      }
    }
    # the failed and held ones, and the commands that would run them again (printed only)
    nbad = 0; held = ""; again = ""; stale = ""
    start = -1
    for (k = 1; k <= maxk + 1; ++k) {
      if ((k in state) && (state[k] == "failed" || state[k] == "held")) {
        if (!nbad++) printf "failed or held:\n"
        printf "  %4d  job %-12s %s: %s\n        logs: %s/beam/log/job%d.err, .log, .out\n", k, htid[k], label[state[k]],
               reason[k] == "" ? "(no reason given)" : reason[k], dir, k
        if (state[k] == "held") held = held " " htid[k]
        if (state[k] == "failed" && done[k]) stale = stale " " dir "/beam/log/job" k ".done"
      }
      hit = (k in state) && state[k] == "failed"
      if (hit && start < 0) start = k
      if (!hit && start >= 0) { again = again "\n    " (resub != "" ? resub : "sh condor/beam_xray_submit.sh --dir " dir) " " start " " (k - 1); start = -1 }
    }
    if (nbad) {
      printf "to run them again (printed only; this script changed nothing):\n"
      if (held != "") {
        printf "  held, after fixing the cause:   condor_release%s\n", held
        printf "  or take them out of the queue:  condor_rm%s   (then submit them again as below)\n", held
      }
      if (stale != "")
        printf "  first remove the .done markers of these failed jobs (the submit script skips a process that has one):\n    rm%s\n", stale
      if (again != "")
        printf "  failed, once batch %s has left the queue (condor_q):%s\n", name, again
    } else if (count["ok"] == total) printf "all %d jobs finished OK\n", total
  }'

# the Geant4 version of each finished job (its job<k>.out: "done in N s with geant4-...:"): one folder should hold one
# version, as a run mixed from two builds differs where the versions differ
versions=$(cat "$dir"/beam/log/job*.out 2> /dev/null | sed -n 's/.* done in [0-9]* s with \([^:]*\):.*/\1/p' | sort | uniq -c)
if [ -n "$versions" ]; then
  printf 'Geant4 of the finished jobs: %s\n' "$(printf '%s\n' "$versions" | awk '{ printf "%s%s (%d jobs)", (NR > 1 ? ", " : ""), $2, $1 }')"
  [ "$(printf '%s\n' "$versions" | wc -l)" -gt 1 ] && echo "WARNING: this folder mixes Geant4 versions; run each in a folder of its own"
fi
exit 0
