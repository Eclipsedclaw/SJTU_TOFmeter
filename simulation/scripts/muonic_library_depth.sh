#!/bin/sh
# Tunes the beam energy of each material of the muonic X-ray library (2026-11-23, data/muonic_library_materials.txt):
# the kinetic energy at which the steel study's mu- beam (macro/beam_steel_lyso.mac, 20 x 20 x 1 cm plate on the CZT
# plane) stops at the steel run's mean depth, 8.2424 mm below the plate's top (the user: the same depth for every
# material, "a material property"). Local TOFmeter runs on this computer, several materials at the same time, one
# thread each, at most 90% of the physical cores. Each material in passes of --muons muons, without --preinit
# (radioactive decay and fluorescence do not change where a muon stops):
#   pass 1 at the table's ekin; pass 2 at ekin x (target / depth)^(1 / 1.75) (range ~ E^1.75); later passes on the
#   straight line through the last two passes in log(depth) against log(ekin); until the depth is within --tolerance
#   of the target, at most 6 passes (a pass whose muons mostly left the plate steps 15% down instead);
# then a check run at the tuned energy with the production setup (--preinit macro/preinit_lyso.mac, as the Condor
# jobs), with --check-muons muons: its depth, its speed (the Geant4 run timer: event loop only) and its output size.
# Files, per material: <dir>/<name>/pass<k>.mac, .log, .root, check.mac, .log, .root, passes.txt, result.txt; at the end
# <dir>/materials_tuned.txt (the table with the tuned ekin, p, depth and rms: copy it over the table after a look) and
# <dir>/summary.txt. Materials with a result.txt are not run again (delete it to redo one).
# Refuses to start at 85% memory use or more. In the build directory, with Geant4 and ROOT set up:
#   sh scripts/muonic_library_depth.sh [--jobs 10] [--muons 2000] [--check-muons 2000] [--tolerance 0.02]
#        [--target 8.2424] [--only name,name] [--table data/muonic_library_materials.txt]
#        [--dir output/muonic_library_czt/tuning]

# ---- one material (run by the driver below through xargs; the settings come in the LIB_* variables) ----
if [ "$1" = --worker ]; then
  name=$2
  table=$LIB_TABLE; dir=$LIB_DIR; muons=$LIB_MUONS; checkmuons=$LIB_CHECK; tol=$LIB_TOL; target=$LIB_TARGET
  row=$(awk -v n="$name" '!/^#/ && $1 == n { print; exit }' "$table")
  if [ -z "$row" ]; then
    echo "$name: not in $table" >&2
    exit 1
  fi
  set -- $row
  id=$2; comp=$3; dens=$4; e=$5
  mkdir -p "$dir/$name" || exit 1
  d=$(cd "$dir/$name" && pwd)
  [ -s "$d/result.txt" ] && exit 0
  seed=$((12900000 + 100 * id))
  # run <base> <ekin> <seed> <events> [preinit]: one TOFmeter run into <base>.root
  run() {
    {
      echo "/control/verbose 2"
      echo "/control/alias bmSeed $3"
      echo "/control/alias bmEvents $4"
      echo "/control/alias bmDir $(dirname "$1")"
      echo "/control/alias bmFile $(basename "$1")"
      echo "/control/alias bmEkin $2"
      echo "/control/alias bmSigX 1.656"
      echo "/control/alias bmSigY 1.062"
      echo "/control/alias bmZ 33.55"
      echo "/control/alias bmFilter 0"
      echo "/control/alias bmPanel czt"
      echo "/control/alias bmPanelZ 32.0"
      echo "/tof/target/composition $comp"
      echo "/tof/target/density $dens g/cm3"
      echo "/control/execute macro/beam_steel_lyso.mac"
    } > "$1.mac"
    ./TOFmeter ${5:+--preinit "$5"} "$1.mac" > "$1.log" 2>&1
    grep -q "Number of events processed : $4\$" "$1.log" && [ -s "$1.root" ]
  }
  depth() {   # the depth line of one output
    ROOT_HIST=0 root -l -b -q "analysis/muonic_library_depth.C(\"$1\")" 2>/dev/null | grep '^depth '
  }
  : > "$d/passes.txt"
  k=0; pe=; pd=; done=no
  while [ $k -lt 6 ]; do
    k=$((k + 1))
    if ! run "$d/pass$k" "$e" $((seed + k)) "$muons"; then
      echo "$name: pass $k failed, see $d/pass$k.log" >&2
      exit 1
    fi
    line=$(depth "$d/pass$k.root")
    if [ -z "$line" ]; then
      echo "$name: no depth from $d/pass$k.root (ROOT)" >&2
      exit 1
    fi
    echo "pass $k ekin $e $line" >> "$d/passes.txt"
    set -- $line
    dep=$2; stopped=${10}; all=${12}
    verdict=$(awk -v a="$dep" -v t="$target" -v tol="$tol" -v s="$stopped" -v n="$all" 'BEGIN {
      if (s < 0.9 * n || a != a + 0) print "deep"; else if (a - t < tol && t - a < tol) print "ok"; else print "step" }')
    if [ "$verdict" = ok ]; then
      done=yes
      break
    elif [ "$verdict" = deep ]; then
      new=$(awk -v e="$e" 'BEGIN { printf "%.4f", 0.85 * e }')
      pe=; pd=
    elif [ -z "$pe" ]; then
      new=$(awk -v e="$e" -v a="$dep" -v t="$target" 'BEGIN { printf "%.4f", e * exp(log(t / a) / 1.75) }')
      pe=$e; pd=$dep
    else
      new=$(awk -v e1="$pe" -v a1="$pd" -v e2="$e" -v a2="$dep" -v t="$target" 'BEGIN {
        n = (e2 != e1) ? (log(a2) - log(a1)) / (log(e2) - log(e1)) : 1.75
        if (n != n || n < 1.2 || n > 2.5) n = 1.75
        printf "%.4f", e2 * exp(log(t / a2) / n) }')
      pe=$e; pd=$dep
    fi
    e=$new
  done
  if [ "$done" != yes ]; then
    echo "$name: not within $tol mm of $target mm after $k passes, see $d/passes.txt" >&2
    exit 1
  fi
  if ! run "$d/check" "$e" $((seed + 99)) "$checkmuons" macro/preinit_lyso.mac; then
    echo "$name: check run failed, see $d/check.log" >&2
    exit 1
  fi
  line=$(depth "$d/check.root")
  real=$(sed -n 's/.*User=\([0-9.]*\)s Real=\([0-9.]*\)s.*/\2/p' "$d/check.log" | tail -n 1)
  bytes=$(wc -c < "$d/check.root")
  awk -v n="$name" -v id="$id" -v c="$comp" -v r="$dens" -v e="$e" -v k="$k" -v t="${real:-nan}" -v b="$bytes" \
      -v m="$checkmuons" -v l="$line" 'BEGIN {
    split(l, f, " ")
    printf "%s %s %s %s %.4f %.2f %s %s %s %s %s %s %d %.3g %.0f\n", n, id, c, r, e, sqrt(e * e + 2 * e * 105.6583755),
           f[2], f[4], f[6], f[8], f[10], f[12], k, t / m * 1000, b / m }' > "$d/result.txt"
  exit 0
fi

# ---- the driver ----
jobs=10; muons=2000; checkmuons=2000; tol=0.02; target=8.2424; only=
table=data/muonic_library_materials.txt; dir=output/muonic_library_czt/tuning
while [ $# -gt 0 ]; do
  case "$1" in
    --jobs) jobs=$2; shift 2 ;;
    --muons) muons=$2; shift 2 ;;
    --check-muons) checkmuons=$2; shift 2 ;;
    --tolerance) tol=$2; shift 2 ;;
    --target) target=$2; shift 2 ;;
    --only) only=$2; shift 2 ;;
    --table) table=$2; shift 2 ;;
    --dir) dir=${2%/}; shift 2 ;;
    *) echo "usage: sh scripts/muonic_library_depth.sh [--jobs n] [--muons n] [--check-muons n] [--tolerance mm] [--target mm] [--only name,name] [--table file] [--dir path]" >&2
       exit 1 ;;
  esac
done
for v in "$jobs" "$muons" "$checkmuons"; do
  case "$v" in
    ''|*[!0-9]*) echo "--jobs, --muons and --check-muons take whole numbers" >&2; exit 1 ;;
  esac
done
limit=$(( $(lscpu -p=Core,Socket | grep -v '^#' | sort -u | wc -l) * 9 / 10 ))
if [ "$jobs" -lt 1 ] || [ "$jobs" -gt "$limit" ]; then
  echo "--jobs must be 1 to $limit (90% of the physical cores)" >&2
  exit 1
fi
for f in ./TOFmeter macro/beam_steel_lyso.mac macro/preinit_lyso.mac macro/panel_czt.mac analysis/muonic_library_depth.C "$table"; do
  if [ ! -e "$f" ]; then
    echo "Run this in the build directory ($f is missing)" >&2
    exit 1
  fi
done
used=$(awk '/^MemTotal/{t=$2} /^MemAvailable/{a=$2} END{printf "%d", 100*(t-a)/t}' /proc/meminfo)
if [ "$used" -ge 85 ]; then
  echo "memory: $used% in use, not starting (85% or more)" >&2
  exit 1
fi
names=$(awk '!/^#/ && NF >= 5 { print $1 }' "$table")
if [ -n "$only" ]; then
  names=$(printf '%s\n' "$only" | tr ',' '\n')
fi
mkdir -p "$dir"
echo "$(date '+%F %T')  tuning $(printf '%s\n' $names | wc -l) materials to $target mm (+-$tol), $muons muons a pass, $jobs at a time, in $dir"
export LIB_TABLE="$table" LIB_DIR="$dir" LIB_MUONS="$muons" LIB_CHECK="$checkmuons" LIB_TOL="$tol" LIB_TARGET="$target"
printf '%s\n' $names | xargs -P "$jobs" -n 1 sh "$0" --worker
status=$?

# the tuned table and the summary, from the results there are
awk -v dir="$dir" '
  /^#/ { print; next }
  NF < 5 { print; next }
  { note = ""; for (i = 9; i <= NF; ++i) note = note (i > 9 ? " " : "") $i
    e = $5; p = $6; dep = $7; rms = $8
    f = dir "/" $1 "/result.txt"; line = ""
    if ((getline line < f) > 0) { split(line, r, " "); e = r[5]; p = r[6]; dep = r[7]; rms = r[9] }
    close(f)
    printf "%-9s %-4s %-13s %-8s %-8s %-8s %-8s %-8s %s\n", $1, $2, $3, $4, e, p, dep, rms, note }' "$table" > "$dir/materials_tuned.txt"
{
  echo "# name id composition density ekin(MeV) p(MeV/c) | check run: depth err rms median (mm) stopped events | passes ms/muon bytes/muon"
  for n in $names; do
    if [ -s "$dir/$n/result.txt" ]; then cat "$dir/$n/result.txt"; else echo "$n not done"; fi
  done
} > "$dir/summary.txt"
echo "$(date '+%F %T')  done (status $status): $dir/summary.txt, $dir/materials_tuned.txt"
exit $status
