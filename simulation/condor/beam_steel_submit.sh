#!/bin/sh
# The steel study with U, Th and K impurities (2026-11-10) as HTCondor jobs, on a cluster with a shared filesystem such
# as SJTU's INPAC cluster (condor/README.md): a 20 x 20 x 1 cm plate of stainless steel lying on the CZT plane
# (macro/panel_czt.mac: 4 x 4 pixels of 25 x 25 x 10 mm, their top against the plate), one core per job, in two bands:
#   beam:  mu- of a mono-energetic beam into the plate (macro/beam_steel_lyso.mac): every event written (/EventFilter 0),
#          so that analysis/beam_steel_trace_spectrum.C counts the muons each element caught (cap_z) and reweights them;
#   decay: the plate's natural radioactivity (macro/beam_steel_decay.mac), one source per process, each with its whole
#          chain in an event (secular equilibrium): 1 U-238, 2 Th-232, 3 K-40, 4 U-235; only events with energy in
#          the panel are written.
# Defaults:
#   - target: --composition Fe:rest,Cr:16,Ni:10,Mo:2,Mn:2,Si:0.75,N:0.10,C:0.08,P:0.045,S:0.03,U:0.1,Th:0.1,K:0.1 (mass
#     percent; the user's steel, with U, Th and K at 0.1% each so that enough muons catch on them; the analysis reweights
#     every element to its real share), --density 8.0 (g/cm3);
#   - beam: PSI muE1's spot (39 x 25 mm FWHM, sigma 1.656 x 1.062 cm), 29.95 MeV kinetic energy (85 MeV/c; it stops about
#     0.8 cm deep, 2 mm above the panel), starting 0.5 mm above the plate;
#   - 5 x 10^5 muons per beam process (~15 min plus ~3 min building the cascade models of ~40 isotopes), seeds
#     10900000 + k; 5 x 10^6 decays per decay process (~45 min), seeds 11000000 + k (--seed-base b: b + k and
#     b + 100000 + k);
#   - --detector czt (lyso: the LYSO panel instead).
# The run's folder is output/beam_steel_<detector>_trace (--dir to choose; --storage <path> puts it under <path>/output/,
# e.g. on the cluster's Lustre storage). One folder holds one setup: composition, density, detector, beam energy and spot.
# At most 60 jobs at a time (max_materialize; 60 cores may be used), and no submission while jobs of the same folder are
# still in the queue. The jobs run beam_steel_job.sh, the one beside this script, in this build directory. Their
# environment is this shell's (getenv: set up Geant4 first); with --setup <file> each job also sources that file before
# TOFmeter starts. In the build directory:
#   sh condor/beam_steel_submit.sh [--composition list] [--density g/cm3] [--detector czt|lyso] [--ekin MeV]
#        [--sigma-x cm] [--sigma-y cm] [--muons n] [--decays n] [--no-decays] [--no-beam] [--seed-base b] [--dir path]
#        [--storage path] [--max n] [--setup file] [--os-major n] [--no-submit] <first> <last>
#   sh condor/beam_steel_submit.sh --test [the same options]   (beam process 1 with 2 x 10^4 muons and the 4 decay
#                                                                processes with 2 x 10^4 decays, into <the folder>_condor_test)
# e.g., 10^8 muons in 200 processes and the 4 decay processes on the INPAC cluster's CentOS 7 workers:
#   sh condor/beam_steel_submit.sh --setup ~/tof_setup_c7.sh --os-major 7 --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work 1 200
# Of the beam processes first-last and the decay processes 1-4, those without <dir>/<band>/log/job<k>.done are
# submitted; --no-decays leaves out the decay band, --no-beam the beam band (then no <first> <last>).
# Files:
#   - <dir>/<band>/log/job<k>.mac, and beside them TOFmeter's log job<k>.log and the job's output and errors .out, .err;
#   - the outputs <dir>/beam/beam_<k>.root and <dir>/decay/decay_<k>.root;
#   - per submission <dir>/condor/submit_<date>.sub, with HTCondor's log of its jobs (submit_<date>.log).
# Check them with sh condor/beam_steel_check.sh <dir> (each job's state, verified; read-only). --os-major n: run only on
# workers of that OS major version (TARGET.OpSysMajorVer, 7 for CentOS 7), the one TOFmeter was built on.

dir=
storage=
composition=Fe:rest,Cr:16,Ni:10,Mo:2,Mn:2,Si:0.75,N:0.10,C:0.08,P:0.045,S:0.03,U:0.1,Th:0.1,K:0.1
density=8.0
detector=czt
ekin=29.95
sigmax=1.656
sigmay=1.062
seedbase=10900000
muons=500000
decays=5000000
dobeam=yes
dodecays=yes
max=60
setup=
osmajor=
submit=yes
test=no
while [ $# -gt 0 ]; do
  case "$1" in
    --dir) dir=${2%/}; shift 2 ;;
    --storage) storage=${2%/}; shift 2 ;;
    --composition) composition=$2; shift 2 ;;
    --density) density=$2; shift 2 ;;
    --detector) detector=$2; shift 2 ;;
    --ekin) ekin=$2; shift 2 ;;
    --sigma-x) sigmax=$2; shift 2 ;;
    --sigma-y) sigmay=$2; shift 2 ;;
    --seed-base) seedbase=$2; shift 2 ;;
    --muons) muons=$2; shift 2 ;;
    --decays) decays=$2; shift 2 ;;
    --no-decays) dodecays=no; shift ;;
    --no-beam) dobeam=no; shift ;;
    --max) max=$2; shift 2 ;;
    --setup) setup=$2; shift 2 ;;
    --os-major) osmajor=$2; shift 2 ;;
    --no-submit) submit=no; shift ;;
    --test) test=yes; shift ;;
    -*) echo "unknown option $1" >&2; exit 1 ;;
    *) break ;;
  esac
done

# a positive decimal number (energy, spot, density)
positive() { printf '%s' "$1" | grep -Eq '^([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$' && awk -v v="$1" 'BEGIN { exit !(v > 0) }'; }
if ! printf '%s' "$composition" | grep -Eq '^([A-Z][a-z]?:([0-9]+\.?[0-9]*|\.[0-9]+|rest),)*[A-Z][a-z]?:([0-9]+\.?[0-9]*|\.[0-9]+|rest)$'; then
  echo "--composition: element:mass-percent pairs separated by commas, rest for the balance of one, e.g. Fe:rest,Cr:16,Ni:10" >&2
  exit 1
fi
for v in "$density" "$ekin" "$sigmax" "$sigmay"; do
  positive "$v" || { echo "--density, --ekin, --sigma-x and --sigma-y take positive numbers (g/cm3, MeV, cm)" >&2; exit 1; }
done
case "$detector" in
  czt) panelz=32.0 ;;     # the centre of the 10 mm CZT, its top against the plate at 32.5 cm
  lyso) panelz=32.332 ;;  # the centre of the 3 mm LYSO layer under the 0.1 mm Al and 0.08 mm ESR
  *) echo "--detector: czt or lyso" >&2; exit 1 ;;
esac
for v in "$muons" "$decays" "$max" "$seedbase"; do
  case "$v" in
    ''|*[!0-9]*) echo "--muons, --decays, --max and --seed-base take whole numbers" >&2; exit 1 ;;
  esac
done
if [ "$muons" -lt 1 ] || [ "$decays" -lt 1 ]; then
  echo "--muons, --decays: at least 1" >&2
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
if [ "$dobeam" = no ] && [ "$dodecays" = no ]; then
  echo "--no-beam and --no-decays together: nothing to submit" >&2
  exit 1
fi

# the run's folder
home=${storage:+$storage/}output/beam_steel_${detector}_trace
[ -n "$dir" ] || dir=$home
if [ "$test" = yes ]; then
  dir=${dir}_condor_test
  muons=20000
  decays=20000
  [ "$dobeam" = yes ] && set -- 1 1
fi
if [ "$dobeam" = yes ]; then
  if [ $# -ne 2 ]; then
    echo "usage: sh condor/beam_steel_submit.sh [--composition list] [--density g/cm3] [--detector czt|lyso] [--ekin MeV] [--sigma-x cm] [--sigma-y cm] [--muons n] [--decays n] [--no-decays] [--no-beam] [--seed-base b] [--dir path] [--storage path] [--max n] [--setup file] [--os-major n] [--no-submit] <first> <last>" >&2
    echo "       sh condor/beam_steel_submit.sh --test [the same options]" >&2
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
elif [ $# -ne 0 ]; then
  echo "--no-beam: no beam processes, so no <first> <last>" >&2
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
beammac=macro/beam_steel_lyso.mac
decaymac=macro/beam_steel_decay.mac
preinit=macro/preinit_lyso.mac
job=$(cd "$(dirname "$0")" && pwd -P)/beam_steel_job.sh
for f in ./TOFmeter "$beammac" "$decaymac" "$preinit" "macro/panel_$detector.mac" "$job"; do
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

# the job macros' setup lines: every macro of a band must have the same ones
targetlines() {
  echo "/control/alias bmPanel $detector"
  echo "/control/alias bmPanelZ $panelz"
  echo "/tof/target/composition $composition"
  echo "/tof/target/density $density g/cm3"
}
setuplines() {   # <band>
  if [ "$1" = beam ]; then
    echo "/control/alias bmEkin $ekin"
    echo "/control/alias bmSigX $sigmax"
    echo "/control/alias bmSigY $sigmay"
    echo "/control/alias bmZ 33.55"
    echo "/control/alias bmFilter 0"
    targetlines
    echo "/control/execute $beammac"
  else
    targetlines
    echo "/control/execute $decaymac"
  fi
}
# the decay sources: process k -> Z A name
decaysource() {
  case "$1" in
    1) echo "92 238 U-238" ;;
    2) echo "90 232 Th-232" ;;
    3) echo "19 40 K-40" ;;
    4) echo "92 235 U-235" ;;
  esac
}
bands=
[ "$dobeam" = yes ] && bands=beam
[ "$dodecays" = yes ] && bands="$bands decay"
for band in $bands; do
  mkdir -p "$dir/$band/log"
  want=$(setuplines "$band")
  for mac in "$dir/$band"/log/job*.mac; do
    [ -e "$mac" ] || continue
    have=$(grep -E '^/control/alias bm(Ekin|SigX|SigY|Z|Filter|Panel|PanelZ) |^/tof/target/|^/control/execute ' "$mac")
    if [ "$have" != "$want" ]; then
      echo "$dir holds $band jobs of another setup ($mac):" >&2
      printf '%s\n' "$have" | sed 's/^/  /' >&2
      echo "give another --dir for this one (composition $composition, $density g/cm3, detector $detector, $ekin MeV, spot $sigmax x $sigmay cm)" >&2
      exit 1
    fi
  done
done

# the macros of the processes to run
items=
n=0
nbeam=0
ndecay=0
finished=0
for band in $bands; do
  sub=$dir/$band
  [ -e "$sub/started" ] || date +%s > "$sub/started"
  if [ "$band" = beam ]; then list=$(seq "$first" "$last"); else list="1 2 3 4"; fi
  for k in $list; do
    if [ -e "$sub/log/job$k.done" ]; then
      finished=$((finished + 1))
      continue
    fi
    {
      echo "# preinit $preinit"   # read by beam_steel_job.sh: TOFmeter --preinit <file>
      echo "/control/verbose 2"
      if [ "$band" = beam ]; then
        echo "/control/alias bmSeed $((seedbase + k))"
        echo "/control/alias bmEvents $muons"
        echo "/control/alias bmDir $sub"
        echo "/control/alias bmFile beam_$k"
      else
        set -- $(decaysource "$k")
        echo "# source $3"
        echo "/control/alias bmSeed $((seedbase + 100000 + k))"
        echo "/control/alias bmEvents $decays"
        echo "/control/alias bmDir $sub"
        echo "/control/alias bmFile decay_$k"
        echo "/control/alias bmIonZ $1"
        echo "/control/alias bmIonA $2"
      fi
      setuplines "$band"
    } > "$sub/log/job$k.mac"
    items="${items}$band $k
"
    n=$((n + 1))
    if [ "$band" = beam ]; then nbeam=$((nbeam + 1)); else ndecay=$((ndecay + 1)); fi
  done
done
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
  echo "# HTCondor jobs of the steel study ($beammac, $decaymac, --preinit $preinit), by condor/beam_steel_submit.sh on $(date '+%F %T')"
  echo "# in $here: $nbeam beam processes of $muons muons${first:+ (of $first-$last)}, $ndecay decay processes of $decays decays ($finished done before), at most $max at a time;"
  echo "# target $composition (mass percent) at $density g/cm3 on the $detector panel; mu- of $ekin MeV, spot sigma $sigmax x $sigmay cm; seeds $seedbase + k (beam), $seedbase + 100000 + k (decay);"
  echo "# environment: the submitting shell's${setup:+, then $setup sourced by each job}"
  echo "# resubmit: sh condor/beam_steel_submit.sh --composition $composition --density $density --detector $detector --ekin $ekin --sigma-x $sigmax --sigma-y $sigmay --seed-base $seedbase --dir $dir --muons $muons --decays $decays --max $max${setup:+ --setup $setup}${osmajor:+ --os-major $osmajor}"
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
echo "$(date '+%F %T')  steel study, $composition on the $detector panel: $nbeam beam processes of $muons muons, $ndecay decay processes of $decays decays ($finished done before), at most $max at a time: $file"
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
echo "check:  sh condor/beam_steel_check.sh $dir   (each job: queued, running, held, finished OK, failed; read-only)"
