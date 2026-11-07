#!/bin/sh
# The beam study of 2026-11-07 as HTCondor jobs, on a cluster with a shared filesystem such as SJTU's INPAC cluster
# (condor/README.md). A mono-energetic mu- beam goes into the 20 x 20 x 4 cm target lying on the LYSO panel
# (macro/beam_xray_lyso.mac, TOFmeter --preinit macro/preinit_lyso.mac); one job of one core per process.
# Defaults:
#   - beam: PSI muE1's spot (39 x 25 mm FWHM, sigma 1.656 x 1.062 cm), 29.95 MeV kinetic energy (85 MeV/c; it stops
#     at 2.5 cm depth in Si, 1.5 cm above the panel);
#   - target: G4_Si with 1% of its atoms Al (--impurity Al --fraction 0.01; analysis/beam_impurity_spectrum.C reweights
#     the Al- and Si-captured muons to any concentration; --impurity none for a pure target);
#   - 5 x 10^5 muons per process (~25-40 min), seeds 10700000 + k (--seed-base b: b + k).
# The run's folder is output/beam_xray_lyso_<target>[_<impurity><fraction>], e.g. output/beam_xray_lyso_si_al0.01
# (--dir to choose; --storage <path> puts it under <path>/output/, e.g. on the cluster's Lustre storage). One folder holds
# one setup: target, impurity, beam energy and spot.
# At most 60 jobs at a time (max_materialize; 60 cores may be used), and no submission while jobs of the same folder are
# still in the queue. The jobs run beam_xray_job.sh, the one beside this script, in this build directory. Their
# environment is this shell's (getenv: set up Geant4 first); with --setup <file> each job also sources that file before
# TOFmeter starts. In the build directory:
#   sh condor/beam_xray_submit.sh [--target G4_Si] [--impurity Al|none] [--fraction 0.01] [--ekin MeV] [--sigma-x cm]
#        [--sigma-y cm] [--muons n] [--seed-base b] [--dir path] [--storage path] [--max n] [--setup file]
#        [--os-major n] [--no-submit] <first> <last>
#   sh condor/beam_xray_submit.sh --test [the same options]   (process 1, 2 x 10^4 muons, into <the folder>_condor_test)
# e.g., 10^8 muons in 200 processes on the INPAC cluster's CentOS 7 workers:
#   sh condor/beam_xray_submit.sh --setup ~/tof_setup_c7.sh --os-major 7 1 200
# Of the processes first-last, those without <dir>/beam/log/job<k>.done are submitted.
# Files:
#   - <dir>/beam/log/job<k>.mac, and beside them TOFmeter's log job<k>.log and the job's output and errors job<k>.out, .err;
#   - the outputs <dir>/beam/beam_<k>.root;
#   - per submission <dir>/condor/submit_<date>.sub, with HTCondor's log of its jobs (submit_<date>.log).
# Check them with sh condor/beam_xray_check.sh <dir> (each job's state, verified; read-only). --os-major n: run only on
# workers of that OS major version (TARGET.OpSysMajorVer, 7 for CentOS 7), the one TOFmeter was built on.

dir=
storage=
target=G4_Si
impurity=Al
fraction=0.01
ekin=29.95
sigmax=1.656
sigmay=1.062
seedbase=10700000
muons=500000
max=60
setup=
osmajor=
submit=yes
test=no
while [ $# -gt 0 ]; do
  case "$1" in
    --dir) dir=${2%/}; shift 2 ;;
    --storage) storage=${2%/}; shift 2 ;;
    --target) target=$2; shift 2 ;;
    --impurity) impurity=$2; shift 2 ;;
    --fraction) fraction=$2; shift 2 ;;
    --ekin) ekin=$2; shift 2 ;;
    --sigma-x) sigmax=$2; shift 2 ;;
    --sigma-y) sigmay=$2; shift 2 ;;
    --seed-base) seedbase=$2; shift 2 ;;
    --muons) muons=$2; shift 2 ;;
    --max) max=$2; shift 2 ;;
    --setup) setup=$2; shift 2 ;;
    --os-major) osmajor=$2; shift 2 ;;
    --no-submit) submit=no; shift ;;
    --test) test=yes; shift ;;
    -*) echo "unknown option $1" >&2; exit 1 ;;
    *) break ;;
  esac
done

# a positive decimal number (energy, spot, fraction)
positive() { printf '%s' "$1" | grep -Eq '^([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$' && awk -v v="$1" 'BEGIN { exit !(v > 0) }'; }
case "$target" in
  G4_[A-Za-z]*) ;;
  *) echo "--target: a NIST material name such as G4_Si" >&2; exit 1 ;;
esac
case "$target" in
  *[!A-Za-z0-9_]*) echo "--target: a NIST material name such as G4_Si" >&2; exit 1 ;;
esac
if [ "$impurity" != none ]; then
  printf '%s' "$impurity" | grep -Eq '^[A-Z][a-z]?$' || { echo "--impurity: an element symbol such as Al, or none" >&2; exit 1; }
  if ! positive "$fraction" || ! awk -v f="$fraction" 'BEGIN { exit !(f < 1) }'; then
    echo "--fraction: an atom fraction between 0 and 1, e.g. 0.01" >&2
    exit 1
  fi
fi
for v in "$ekin" "$sigmax" "$sigmay"; do
  positive "$v" || { echo "--ekin, --sigma-x and --sigma-y take positive numbers (MeV, cm)" >&2; exit 1; }
done
for v in "$muons" "$max" "$seedbase"; do
  case "$v" in
    ''|*[!0-9]*) echo "--muons, --max and --seed-base take whole numbers" >&2; exit 1 ;;
  esac
done
if [ "$muons" -lt 1 ]; then
  echo "--muons: at least 1" >&2
  exit 1
fi
if [ "$max" -lt 1 ] || [ "$max" -gt 60 ]; then
  echo "--max: 1-60 jobs at a time (one core each, 60 cores at most)" >&2
  exit 1
fi
case "$osmajor" in
  *[!0-9]*) echo "--os-major: a whole number, e.g. 7" >&2; exit 1 ;;
esac
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
name=$(printf '%s' "${target#G4_}" | tr 'A-Z' 'a-z')
[ "$impurity" != none ] && name=${name}_$(printf '%s' "$impurity" | tr 'A-Z' 'a-z')$fraction
home=${storage:+$storage/}output/beam_xray_lyso_$name
[ -n "$dir" ] || dir=$home
if [ "$test" = yes ]; then
  dir=${dir}_condor_test
  muons=20000
  set -- 1 1
fi
if [ $# -ne 2 ]; then
  echo "usage: sh condor/beam_xray_submit.sh [--target G4_Si] [--impurity Al|none] [--fraction 0.01] [--ekin MeV] [--sigma-x cm] [--sigma-y cm] [--muons n] [--seed-base b] [--dir path] [--storage path] [--max n] [--setup file] [--os-major n] [--no-submit] <first> <last>" >&2
  echo "       sh condor/beam_xray_submit.sh --test [the same options]" >&2
  exit 1
fi
first=$1
last=$2
for v in "$first" "$last"; do
  case "$v" in
    ''|*[!0-9]*) echo "$first $last: processes are numbered 1, 2, ..." >&2; exit 1 ;;
  esac
done
if [ "$first" -lt 1 ] || [ "$last" -lt "$first" ]; then
  echo "$first $last: first at least 1, last at least first" >&2
  exit 1
fi
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
mainmac=macro/beam_xray_lyso.mac
preinit=macro/preinit_lyso.mac
job=$(cd "$(dirname "$0")" && pwd -P)/beam_xray_job.sh
for f in ./TOFmeter "$mainmac" "$preinit" macro/panel_lyso.mac "$job"; do
  if [ ! -e "$f" ]; then
    echo "Run this in the build directory ($f is missing; after a new macro, run cmake again)" >&2
    exit 1
  fi
done
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

# the job macro's setup lines: every macro of a folder must have the same ones
setuplines() {
  echo "/control/alias bmEkin $ekin"
  echo "/control/alias bmSigX $sigmax"
  echo "/control/alias bmSigY $sigmay"
  echo "/control/alias bmZ 36.55"
  echo "/tof/target/material $target"
  if [ "$impurity" != none ]; then
    echo "/tof/target/impurity $impurity"
    echo "/tof/target/impurityFraction $fraction"
  fi
  echo "/control/execute $mainmac"
}
sub=$dir/beam
mkdir -p "$sub/log"
want=$(setuplines)
for mac in "$sub"/log/job*.mac; do
  [ -e "$mac" ] || continue
  have=$(grep -E '^/control/alias bm(Ekin|SigX|SigY|Z) |^/tof/target/|^/control/execute ' "$mac")
  if [ "$have" != "$want" ]; then
    echo "$dir holds jobs of another setup ($mac):" >&2
    printf '%s\n' "$have" | sed 's/^/  /' >&2
    echo "give another --dir for this one (target $target, impurity $impurity $fraction, $ekin MeV, spot $sigmax x $sigmay cm)" >&2
    exit 1
  fi
done

# the macros of the processes to run
items=
n=0
finished=0
[ -e "$sub/started" ] || date +%s > "$sub/started"
for k in $(seq "$first" "$last"); do
  if [ -e "$sub/log/job$k.done" ]; then
    finished=$((finished + 1))
    continue
  fi
  {
    echo "# preinit $preinit"   # read by beam_xray_job.sh: TOFmeter --preinit <file>
    echo "/control/verbose 2"
    echo "/control/alias bmSeed $((seedbase + k))"
    echo "/control/alias bmEvents $muons"
    echo "/control/alias bmDir $sub"
    echo "/control/alias bmFile beam_$k"
    setuplines
  } > "$sub/log/job$k.mac"
  items="${items}beam $k
"
  n=$((n + 1))
done
if [ "$n" -eq 0 ]; then
  echo "nothing to submit: processes $first-$last of $dir are done ($finished .done markers)"
  [ "$test" = yes ] && echo "to test again: rm -r $dir"
  exit 0
fi

# the submit file: the jobs start in this directory (initialdir) and see it, nothing is transferred; bash, as setup
# files are mostly written for it
mkdir -p "$dir/condor"
here=$(pwd -P)
stamp=$(date +%Y%m%d_%H%M%S)
file="$dir/condor/submit_$stamp.sub"
imptext="no impurity"
imptopt=" --impurity none"
if [ "$impurity" != none ]; then
  imptext="$impurity at atom fraction $fraction"
  imptopt=" --impurity $impurity --fraction $fraction"
fi
{
  echo "# HTCondor jobs of the beam study ($mainmac, --preinit $preinit), by condor/beam_xray_submit.sh on $(date '+%F %T')"
  echo "# in $here: processes $first-$last, $n processes of $muons muons ($finished done before), at most $max at a time;"
  echo "# target $target with $imptext; mu- of $ekin MeV, spot sigma $sigmax x $sigmay cm; seeds $seedbase + k;"
  echo "# environment: the submitting shell's${setup:+, then $setup sourced by each job}"
  echo "# resubmit: sh condor/beam_xray_submit.sh --target $target$imptopt --ekin $ekin --sigma-x $sigmax --sigma-y $sigmay --seed-base $seedbase --dir $dir --muons $muons --max $max${setup:+ --setup $setup}${osmajor:+ --os-major $osmajor}"
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
echo "$(date '+%F %T')  beam study, $target with $imptext: $n processes ($first-$last; $finished done before) of $muons muons, at most $max at a time: $file"
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
echo "check:  sh condor/beam_xray_check.sh $dir   (each job: queued, running, held, finished OK, failed; read-only)"
