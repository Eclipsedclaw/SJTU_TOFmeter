#!/bin/sh
# Status of every point of a lead x energy scan of scripts/lead_energy_scan.sh, per lead thickness
# and muon energy, and of the check runs: 1 finished, R running, 0 pending. Read-only; it prints a
# chart, writes <dir>/scan_status.csv and draws <dir>/scan_status.png (analysis/lead_scan_status.C).
# Run it in the build directory:
#   sh scripts/lead_scan_status.sh [definition]     (default macro/lead_thick_scan_def.mac)
# A point counts as finished when its run has ended ("Run end" in the log of its process). The
# finished points of a process that was stopped half-way (no .done marker) are simulated again,
# with the same events, when the scan is resumed; they show as 1 until then.

def=${1:-macro/lead_thick_scan_def.mac}
if [ ! -f "$def" ]; then
  echo "no $def: run this in the build directory" >&2
  exit 1
fi
alias_value() { sed -n "s|^/control/alias $1 *||p" "$def" | tr -d '"'; }
dir=$(alias_value leadScanDir)
tmax=$(alias_value leadScanTmax); tmax=${tmax:-20}
groups=$(alias_value leadScanGroups); groups=${groups:-1}
energies=$(alias_value leadScanEnergies)
log="$dir/log"
running=" $(ps -Ao command= | sed -n "s|^\./TOFmeter $log/\(.*\)\.mac$|\1|p" | tr '\n' ' ') "

# one line per process: thickness, group (-1 = check run), state, runs ended
list_states() {
  for t in $(seq "$tmax" -1 0); do
    g=0
    while [ $g -le "$groups" ]; do
      if [ $g -eq "$groups" ]; then name=check_t$t; gg=-1
      elif [ "$groups" -gt 1 ]; then name=t${t}_g$g; gg=$g
      else name=t$t; gg=$g; fi
      n=0
      [ -e "$log/$name.log" ] && n=$(grep -c 'Run end' "$log/$name.log")
      if [ -e "$log/$name.done" ]; then s=done
      else
        case "$running" in
          *" $name "*) s=running ;;
          *) if [ -e "$log/$name.log" ]; then s=stopped; else s=pending; fi ;;
        esac
      fi
      echo "$t $gg $s $n"
      g=$((g + 1))
    done
  done
}
states=$(list_states)

csv="$dir/scan_status.csv"
echo "Scan $def ($dir), $(date '+%F %T'): 1 finished, R running, 0 pending"
echo "$states" | awk -v energies="$energies" -v groups="$groups" -v csv="$csv" '
  BEGIN {
    nE = split(energies, E, " ")
    size = int((nE + groups - 1) / groups)
    print "lead_cm,E_GeV,status" > csv      # status 1 finished, 0.5 running, 0 pending; E_GeV "check" for the check runs
    printf "%8s", "E (GeV)"
    for (g = 0; g < groups; ++g) {
      i = g * size + 1; j = (g + 1) * size; if (j > nE) j = nE
      printf " %-" (j - i + 1) "s", sprintf("%.3g-%.3g", E[i], E[j])
    }
    print "  check"
  }
  {
    t = $1; g = $2; s = $3; n = $4
    if (!(t in seen)) { seen[t] = 1; order[++nt] = t }
    if (g < 0) {
      c = (s == "done") ? "1" : (s == "running" ? "R" : "0")
      v = (s == "done") ? 1 : (s == "running" ? 0.5 : 0)
      check[t] = c; print t ",check," v > csv
      if (s == "done") ++chk
      next
    }
    cells = ""
    for (k = 0; k < size && g * size + k < nE; ++k) {
      e = E[g * size + k + 1]
      if (s == "done" || (s != "pending" && k < n)) { c = "1"; v = 1; ++fin }
      else if (s == "running" && k == n) { c = "R"; v = 0.5; ++run }
      else { c = "0"; v = 0 }
      cells = cells c
      print t "," e "," v > csv
    }
    row[t] = row[t] " " cells
  }
  END {
    for (r = 1; r <= nt; ++r) printf "%5s cm%s  %s\n", order[r], row[order[r]], check[order[r]]
    printf "points: %d of %d finished, %d running; check runs: %d of %d finished\n", fin, nt * nE, run, chk, nt
  }'
echo "processes: $(ls "$log" 2>/dev/null | grep -c '\.done$') of $(cat "$log/jobs.txt" 2>/dev/null | wc -l | tr -d ' ') done; table: $csv"
if command -v root > /dev/null 2>&1 && [ -f analysis/lead_scan_status.C ]; then
  root -l -b -q "analysis/lead_scan_status.C(\"$dir\")" > /dev/null 2>&1 && echo "chart: $dir/scan_status.png"
fi
