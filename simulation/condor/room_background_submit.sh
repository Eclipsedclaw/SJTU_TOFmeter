#!/bin/sh
# The room background of the steel study (2026-11-20: can a TOF-trigger ON/OFF subtraction show U, Th and K at 50 ppm
# in steel, with the room's background included?) as HTCondor jobs, on a cluster with a shared filesystem such as
# SJTU's INPAC cluster (condor/README.md). The setup is the steel beam run's (condor/beam_steel_submit.sh: the
# 20 x 20 x 1 cm plate on the CZT plane inside TOFmeter's default setup, i.e. under a 10 cm lead roof and between four
# 10 cm lead walls, open at the bottom), plus a 30 cm concrete floor (G4_CONCRETE) under it, at -30..0 cm. One core per
# job, three bands:
#   floor:  the floor's natural radioactivity (macro/room_decay.mac), each event one nucleus with its whole chain
#           (secular equilibrium): U-238, Th-232, K-40 and U-235, in five regions of the floor: inner (x and y within
#           +-55 cm: under the castle's opening and walls; 0.1 day of decays at the analysis's default activities) and
#           the ring around it out to +-100 cm (px, nx: x beyond +-55 cm; py, ny: y beyond +-55 cm within |x| < 55;
#           0.02 day), where the paths to the panel run through more concrete;
#   radon:  Rn-222 with its short-lived daughters in the air inside the castle (macro/room_decay.mac; ~1 day at
#           50 Bq/m3);
#   cosmic: sea-level cosmic rays of every species from EXPACS, open sky (macro/room_cosmic.mac; 0.25 day).
# The analysis (analysis/room_background_spectrum.C) scales each band to its activity or live time; the number of
# decays per job only sets the statistics. The process table (band: process numbers -> source, region, events):
#   floor   1-36   U-238 chain, 1.2e7 each: inner 1-24, px 25-28, nx 29-32, py 33-34, ny 35-36
#          37-64   Th-232 chain, 1.2e7 each: inner 37-54, px 55-57, nx 58-60, py 61-62, ny 63-64
#          65-93   K-40, 1.5e8 each: inner 65-83, px 84-86, nx 87-89, py 90-91, ny 92-93
#          94-98   U-235 chain: inner 94 (1.2e7), px 95, nx 96, py 97, ny 98 (2e6 each)
#   radon   1-4    Rn-222, 8e5 each
#   cosmic  1-25   3.25e6 cosmic-ray particles each
# Speed (local test 2026-11-20, Geant4 11.2.2, event loop): U-238 and Th-232 chains in the floor ~2040 /s, K-40
# ~25000 /s, radon in the castle's air ~140 chains/s (betas and alphas tracked through metres of air), cosmic rays
# ~450 /s (each primary through the lead roof); so ~1.6-2 h per job, ~210 core-hours in all (127 jobs).
# Seeds: floor 11100000 + k, radon 11200000 + k, cosmic 11300000 + k (--seed-base b: b + k, b + 100000 + k,
# b + 200000 + k). The plate: --composition and --density as the beam run (its absorption only).
# The run's folder is output/beam_steel_<detector>_room (--dir to choose; --storage <path> puts it under
# <path>/output/, e.g. on the cluster's Lustre storage). At most 60 jobs at a time (max_materialize), and no submission
# while jobs of the same folder are still in the queue. The jobs run beam_steel_job.sh, the one beside this script (it
# runs one job macro of any band), in this build directory. Their environment is this shell's (getenv: set up Geant4
# first); with --setup <file> each job also sources that file before TOFmeter starts. In the build directory:
#   sh condor/room_background_submit.sh [--bands floor,radon,cosmic] [--composition list] [--density g/cm3]
#        [--detector czt|lyso] [--seed-base b] [--dir path] [--storage path] [--max n] [--setup file] [--os-major n]
#        [--no-submit]
#   sh condor/room_background_submit.sh --test [the same options]   (floor: the 4 sources once in the inner region,
#        radon once, cosmic once, 2 x 10^4 events each, into <the folder>_condor_test)
# e.g. on the INPAC cluster's CentOS 7 workers:
#   sh condor/room_background_submit.sh --setup ~/tof_setup_c7.sh --os-major 7 --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work
# Of the processes of the bands asked for, those without <dir>/<band>/log/job<k>.done are submitted.
# Files:
#   - <dir>/<band>/log/job<k>.mac, and beside them TOFmeter's log job<k>.log and the job's output and errors .out, .err;
#   - the outputs <dir>/<band>/<band>_<k>.root;
#   - per submission <dir>/condor/submit_<date>.sub, with HTCondor's log of its jobs (submit_<date>.log).
# Check them with sh condor/room_background_check.sh <dir> (each job's state, verified; read-only). --os-major n: run
# only on workers of that OS major version (TARGET.OpSysMajorVer, 7 for CentOS 7), the one TOFmeter was built on.

dir=
storage=
bandlist=floor,radon,cosmic
composition=Fe:rest,Cr:16,Ni:10,Mo:2,Mn:2,Si:0.75,N:0.10,C:0.08,P:0.045,S:0.03,U:0.1,Th:0.1,K:0.1
density=8.0
detector=czt
seedbase=11100000
max=60
setup=
osmajor=
submit=yes
test=no
while [ $# -gt 0 ]; do
  case "$1" in
    --bands) bandlist=$2; shift 2 ;;
    --dir) dir=${2%/}; shift 2 ;;
    --storage) storage=${2%/}; shift 2 ;;
    --composition) composition=$2; shift 2 ;;
    --density) density=$2; shift 2 ;;
    --detector) detector=$2; shift 2 ;;
    --seed-base) seedbase=$2; shift 2 ;;
    --max) max=$2; shift 2 ;;
    --setup) setup=$2; shift 2 ;;
    --os-major) osmajor=$2; shift 2 ;;
    --no-submit) submit=no; shift ;;
    --test) test=yes; shift ;;
    -*) echo "unknown option $1" >&2; exit 1 ;;
    *) echo "no arguments besides the options (the process table is fixed): $1" >&2; exit 1 ;;
  esac
done

if ! printf '%s' "$composition" | grep -Eq '^([A-Z][a-z]?:([0-9]+\.?[0-9]*|\.[0-9]+|rest),)*[A-Z][a-z]?:([0-9]+\.?[0-9]*|\.[0-9]+|rest)$'; then
  echo "--composition: element:mass-percent pairs separated by commas, rest for the balance of one, e.g. Fe:rest,Cr:16,Ni:10" >&2
  exit 1
fi
if ! printf '%s' "$density" | grep -Eq '^([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$' || ! awk -v v="$density" 'BEGIN { exit !(v > 0) }'; then
  echo "--density takes a positive number (g/cm3)" >&2
  exit 1
fi
case "$detector" in
  czt) panelz=32.0 ;;     # the centre of the 10 mm CZT, its top against the plate at 32.5 cm
  lyso) panelz=32.332 ;;  # the centre of the 3 mm LYSO layer under the 0.1 mm Al and 0.08 mm ESR
  *) echo "--detector: czt or lyso" >&2; exit 1 ;;
esac
for v in "$max" "$seedbase"; do
  case "$v" in
    ''|*[!0-9]*) echo "--max and --seed-base take whole numbers" >&2; exit 1 ;;
  esac
done
if [ "$max" -lt 1 ] || [ "$max" -gt 60 ]; then
  echo "--max: 1-60 jobs at a time (one core each, 60 cores at most)" >&2
  exit 1
fi
case "$osmajor" in
  *[!0-9]*) echo "--os-major: a whole number, e.g. 7" >&2; exit 1 ;;
esac
bands=
for b in $(printf '%s' "$bandlist" | tr ',' ' '); do
  case "$b" in
    floor|radon|cosmic) case " $bands " in *" $b "*) ;; *) bands="$bands $b" ;; esac ;;
    *) echo "--bands: a comma-separated list of floor, radon, cosmic (got $b)" >&2; exit 1 ;;
  esac
done
bands=${bands# }
[ -n "$bands" ] || { echo "--bands: no band given" >&2; exit 1; }
if [ -n "$storage" ]; then
  case "$storage" in
    /*) ;;
    *) echo "--storage: an absolute path, e.g. /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work" >&2; exit 1 ;;
  esac
  if [ ! -d "$storage" ] || [ ! -w "$storage" ]; then
    echo "--storage: $storage is not a folder you can write to (make it first)" >&2
    exit 1
  fi
fi

# the run's folder
home=${storage:+$storage/}output/beam_steel_${detector}_room
[ -n "$dir" ] || dir=$home
basedir=$dir   # the resubmit line names this one (with --test for a test folder)
retest=
if [ "$test" = yes ]; then dir=${dir}_condor_test; retest="--test "; fi
case "$dir$setup" in
  *[[:space:]]*) echo "--dir, --setup and --storage: paths without spaces" >&2; exit 1 ;;
esac
if [ -n "$setup" ]; then
  if [ ! -f "$setup" ] || [ ! -r "$setup" ]; then
    echo "--setup: cannot read $setup" >&2
    exit 1
  fi
  setup=$(cd "$(dirname "$setup")" && pwd -P)/$(basename "$setup")
fi
decaymac=macro/room_decay.mac
cosmicmac=macro/room_cosmic.mac
preinit=macro/preinit_lyso.mac
job=$(cd "$(dirname "$0")" && pwd -P)/beam_steel_job.sh
for f in ./TOFmeter "$decaymac" "$cosmicmac" "$preinit" "macro/panel_$detector.mac" data/expacs_spectrum.txt "$job"; do
  if [ ! -e "$f" ]; then
    echo "Run this in the build directory ($f is missing; after a new macro, run cmake again)" >&2
    exit 1
  fi
done
if ! grep -a -q "Material005" ./TOFmeter; then
  echo "./TOFmeter has no /tof/target/composition: build it again (make) from the current sources" >&2
  exit 1
fi
batch=tof_$(basename "$dir")
if [ "$submit" = yes ]; then
  if ! command -v condor_submit > /dev/null; then
    echo "no condor_submit here: run this on the cluster's login node (or with --no-submit)" >&2
    exit 1
  fi
  # HTCondor's programs must start in this shell; a shell set up for Geant4 (another system's libraries) can keep them
  # from it
  out=$(condor_submit -version 2>&1)
  case "$out" in
    *"error while loading shared libraries"*|*"not found (required by"*)
      printf '%s\n' "$out" >&2
      echo "condor_submit cannot start in this shell's environment; nothing was written or submitted. Run this in a shell" >&2
      echo "where it starts (condor_q is a quick test); add --setup <file> for the jobs' Geant4" >&2
      exit 1 ;;
  esac
  if condor_q -af JobBatchName 2> /dev/null | grep -qx "$batch"; then
    echo "jobs of $dir (batch $batch) are still in the queue (condor_q): submit when they are over, or condor_rm them" >&2
    exit 1
  fi
fi

# the process table: one line per process, "k events what", what = decay <Z> <A> <nuclide> <region> or cosmic
region() {   # <region>: x y z halfx halfy halfz (cm) of the box the decays are placed in, and its volume
  case "$1" in
    inner) echo "0 0 -15 55 55 15 Floor" ;;
    px) echo "77.5 0 -15 22.5 100 15 Floor" ;;
    nx) echo "-77.5 0 -15 22.5 100 15 Floor" ;;
    py) echo "0 77.5 -15 55 22.5 15 Floor" ;;
    ny) echo "0 -77.5 -15 55 22.5 15 Floor" ;;
    castle) echo "0 0 45.5 45 45 45.5 WorldPhys" ;;   # the castle's inside: 90 x 90 cm between the walls, 0-91 cm
  esac
}
plan() {   # <band>
  awk -v band="$1" -v test="$test" 'BEGIN {
    k = 0
    if (band == "floor") {
      # nuclide, events per job, jobs per region (inner px nx py ny)
      n = split("92 238 U-238 1.2e7 24 4 4 2 2|90 232 Th-232 1.2e7 18 3 3 2 2|19 40 K-40 1.5e8 19 3 3 2 2|92 235 U-235 1.2e7 1 0 0 0 0", src, "|")
      split("inner px nx py ny", reg, " ")
      for (i = 1; i <= n; ++i) {
        split(src[i], f, " ")
        for (r = 1; r <= 5; ++r) {
          jobs = f[4 + r]; ev = f[4]
          if (f[3] == "U-235" && r > 1) { jobs = 1; ev = 2e6 }
          if (test == "yes") { if (r > 1) continue; jobs = 1; ev = 20000 }
          for (j = 0; j < jobs; ++j) printf "%d %.0f decay %s %s %s %s\n", ++k, ev, f[1], f[2], f[3], reg[r]
        }
      }
    } else if (band == "radon") {
      jobs = test == "yes" ? 1 : 4
      for (j = 0; j < jobs; ++j) printf "%d %.0f decay 86 222 Rn-222 castle\n", ++k, test == "yes" ? 20000 : 8e5
    } else {
      jobs = test == "yes" ? 1 : 25
      for (j = 0; j < jobs; ++j) printf "%d %.0f cosmic\n", ++k, test == "yes" ? 20000 : 3.25e6
    }
  }'
}

# the job macros' setup lines: every macro of a band must have the same ones
setuplines() {   # <band>
  echo "/control/alias bmPanel $detector"
  echo "/control/alias bmPanelZ $panelz"
  echo "/tof/target/composition $composition"
  echo "/tof/target/density $density g/cm3"
  if [ "$1" = cosmic ]; then echo "/control/execute $cosmicmac"; else echo "/control/execute $decaymac"; fi
}
for band in $bands; do
  mkdir -p "$dir/$band/log"
  want=$(setuplines "$band")
  for mac in "$dir/$band"/log/job*.mac; do
    [ -e "$mac" ] || continue
    have=$(grep -E '^/control/alias bm(Panel|PanelZ) |^/tof/target/|^/control/execute ' "$mac")
    if [ "$have" != "$want" ]; then
      echo "$dir holds $band jobs of another setup ($mac):" >&2
      printf '%s\n' "$have" | sed 's/^/  /' >&2
      echo "give another --dir for this one (composition $composition, $density g/cm3, detector $detector)" >&2
      exit 1
    fi
  done
done

# the macros of the processes to run
items=
n=0
finished=0
summary=
for band in $bands; do
  sub=$dir/$band
  [ -e "$sub/started" ] || date +%s > "$sub/started"
  case "$band" in
    floor) base=$seedbase ;;
    radon) base=$((seedbase + 100000)) ;;
    cosmic) base=$((seedbase + 200000)) ;;
  esac
  nb=0
  plan "$band" > "$sub/log/plan.txt.new"
  while read -r k events kind z a nuclide reg; do
    if [ -e "$sub/log/job$k.done" ]; then
      finished=$((finished + 1))
      continue
    fi
    {
      echo "# preinit $preinit"   # read by beam_steel_job.sh: TOFmeter --preinit <file>
      if [ "$kind" = decay ]; then echo "# source $nuclide in $reg"; else echo "# source cosmic rays"; fi
      echo "/control/verbose 2"
      echo "/control/alias bmSeed $((base + k))"
      echo "/control/alias bmEvents $events"
      echo "/control/alias bmDir $sub"
      echo "/control/alias bmFile ${band}_$k"
      if [ "$kind" = decay ]; then
        set -- $(region "$reg")
        echo "/control/alias bmIonZ $z"
        echo "/control/alias bmIonA $a"
        echo "/control/alias rmRegion $reg"
        echo "/control/alias rmX $1"
        echo "/control/alias rmY $2"
        echo "/control/alias rmZ $3"
        echo "/control/alias rmHalfX $4"
        echo "/control/alias rmHalfY $5"
        echo "/control/alias rmHalfZ $6"
        echo "/control/alias rmVolume $7"
        if [ "$band" = radon ]; then echo "/control/alias rmThreshold 1"; else echo "/control/alias rmThreshold 1.0e+60"; fi
      fi
      setuplines "$band"
    } > "$sub/log/job$k.mac"
    items="${items}$band $k
"
    n=$((n + 1))
    nb=$((nb + 1))
  done < "$sub/log/plan.txt.new"
  mv "$sub/log/plan.txt.new" "$sub/log/plan.txt"
  summary="$summary $band $nb,"
done
summary=${summary%,}
if [ "$n" -eq 0 ]; then
  echo "nothing to submit: the processes of $dir asked for are done ($finished .done markers)"
  [ "$test" = yes ] && echo "to test again: rm -r $dir"
  exit 0
fi

# the submit file: the jobs start in this directory (initialdir) and see it, nothing is transferred; bash, as setup
# files are mostly written for it
mkdir -p "$dir/condor"
here=$(pwd -P)
stamp=$(date +%Y%m%d_%H%M%S)
file="$dir/condor/submit_$stamp.sub"
{
  echo "# HTCondor jobs of the steel study's room background ($decaymac, $cosmicmac, --preinit $preinit), by condor/room_background_submit.sh on $(date '+%F %T')"
  echo "# in $here: processes per band$summary ($finished done before), at most $max at a time;"
  echo "# plate $composition (mass percent) at $density g/cm3 on the $detector panel, 30 cm concrete floor; seeds $seedbase + k (floor), + 100000 + k (radon), + 200000 + k (cosmic);"
  echo "# environment: the submitting shell's${setup:+, then $setup sourced by each job}"
  echo "# resubmit: sh condor/room_background_submit.sh ${retest}--bands $(printf '%s' "$bands" | tr ' ' ',') --composition $composition --density $density --detector $detector --seed-base $seedbase --dir $basedir --max $max${setup:+ --setup $setup}${osmajor:+ --os-major $osmajor}"
  echo "universe              = vanilla"
  echo "executable            = /bin/bash"
  echo "arguments             = \"$job $dir \$(band) \$(k)${setup:+ $setup}\""
  echo "transfer_executable   = False"
  echo "should_transfer_files = NO"
  echo "initialdir            = $here"
  echo "getenv                = True"
  [ -n "$osmajor" ] && echo "requirements          = (TARGET.OpSysMajorVer == $osmajor)"
  echo "request_cpus          = 1"
  echo "request_memory        = 1024"
  echo "notification          = Never"
  echo "batch_name            = $batch"
  echo "max_materialize       = $max"
  echo "on_exit_hold          = ExitCode =!= 0"
  echo "on_exit_hold_reason   = \"TOFmeter failed or stopped early: see $dir/\$(band)/log/job\$(k).err and .log\""
  echo "periodic_hold         = (JobStatus == 2) && (time() - EnteredCurrentStatus > 23 * 3600)"
  echo "periodic_hold_reason  = \"running for over 23 h (jobs must end within 24 h)\""
  echo "log                   = $dir/condor/submit_$stamp.log"
  echo "output                = $dir/\$(band)/log/job\$(k).out"
  echo "error                 = $dir/\$(band)/log/job\$(k).err"
  echo "queue band, k from ("
  printf '%s' "$items"
  echo ")"
} > "$file"
echo "$(date '+%F %T')  steel study's room background on the $detector panel: processes per band$summary ($finished done before), at most $max at a time: $file"
if [ "$submit" = no ]; then
  echo "not submitted (--no-submit); to submit: condor_submit $file"
  exit 0
fi
out=$(condor_submit "$file" 2>&1)
status=$?
printf '%s\n' "$out"
if [ "$status" -ne 0 ]; then
  echo "condor_submit FAILED (its message above)" >&2
  case "$out" in
    *aterializ*) echo "this HTCondor does not take max_materialize: delete that line of $file and submit at most 60 processes at once" >&2 ;;
  esac
  exit 1
fi
echo "check:  sh condor/room_background_check.sh $dir   (each job: queued, running, held, finished OK, failed; read-only)"
