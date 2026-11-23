#!/bin/sh
# The muonic X-ray library (2026-11-23: every element from H to Mo as a pure plate, and two test alloys, each the CZT
# spectrum of 10^8 mu- of the steel study's beam, so that a compound's spectrum can be unmixed into its elements) as
# HTCondor jobs, on a cluster with a shared filesystem such as SJTU's INPAC cluster (condor/README.md). The setup is
# the steel beam run's (condor/beam_steel_submit.sh, macro/beam_steel_lyso.mac: a 20 x 20 x 1 cm plate on the CZT plane
# of 4 x 4 pixels of 25 x 25 x 10 mm, PSI muE1's spot, measured line energies, --preinit macro/preinit_lyso.mac, every
# event written), with the plate made of each material of data/muonic_library_materials.txt (composition, density) and
# the beam's kinetic energy tuned there by scripts/muonic_library_depth.sh, so that the muons stop at the steel run's
# mean depth, 8.2424 mm, in every material. One core per job; per material, in the folder <id>_<name>:
#   --jobs processes (100) of --muons mu- (10^6): 10^8 muons. Processes 1 to --full-jobs (10) write TOFmeter's full
#   events ntuple, the others run TOFmeter --slim (eventID, ekin, end_x/y/z, end_ke, end_vol, cap_z, mux_e, mux_src,
#   czt_pix, czt_src, czt_e, czt_t: about a quarter of the size; the user's choice: slim for 90%, full for 10%).
# Seeds 12000000 + 1000 id + k, id the element's Z (101 and 102 the alloys), k < 1000.
# Cost (local check runs, Geant4 11.2.2, in the tuning's summary.txt; the cluster's cores were ~1.5x faster for the
# steel): see condor/README.md; 44 materials x 10^8 muons, about 0.65 TB.
# The run's folder is output/muonic_library_czt (--dir to choose; --storage <path> puts it under <path>/output/, e.g. on
# the cluster's Lustre storage). At most 60 jobs at a time (max_materialize), and no submission while jobs of the same
# folder are still in the queue. The jobs run beam_steel_job.sh, the one beside this script (it runs one job macro of
# any folder), in this build directory. Their environment is this shell's (getenv: set up Geant4 first); with
# --setup <file> each job also sources that file before TOFmeter starts. In the build directory:
#   sh condor/muonic_library_submit.sh [--materials name,name] [--jobs 100] [--muons 1000000] [--full-jobs 10]
#        [--table data/muonic_library_materials.txt] [--dir path] [--storage path] [--max n] [--setup file]
#        [--os-major n] [--no-submit]
#   sh condor/muonic_library_submit.sh --test [the same options]   (per material process 1 in the full format and
#        process 2 slim, 10^4 muons each, into <the folder>_condor_test)
# e.g. on the INPAC cluster's CentOS 7 workers:
#   sh condor/muonic_library_submit.sh --setup ~/tof_setup_c7.sh --os-major 7 --storage /lustre/neutrino/zengjiancheng/SJTU_TOFmeter-work
# Of the processes of the materials asked for (all by default), those without <dir>/<id>_<name>/log/job<k>.done are
# submitted.
# Files:
#   - <dir>/<id>_<name>/log/job<k>.mac, and beside them TOFmeter's log job<k>.log and the job's output and errors .out,
#     .err;
#   - the outputs <dir>/<id>_<name>/beam_<k>.root;
#   - per submission <dir>/condor/submit_<date>.sub, with HTCondor's log of its jobs (submit_<date>.log).
# Check them with sh condor/muonic_library_check.sh <dir> (each job's state, verified; read-only). --os-major n: run
# only on workers of that OS major version (TARGET.OpSysMajorVer, 7 for CentOS 7), the one TOFmeter was built on.

dir=
storage=
materials=
table=data/muonic_library_materials.txt
jobs=100
muons=1000000
fulljobs=10
max=60
setup=
osmajor=
submit=yes
test=no
while [ $# -gt 0 ]; do
  case "$1" in
    --materials) materials=$2; shift 2 ;;
    --table) table=$2; shift 2 ;;
    --jobs) jobs=$2; shift 2 ;;
    --muons) muons=$2; shift 2 ;;
    --full-jobs) fulljobs=$2; shift 2 ;;
    --dir) dir=${2%/}; shift 2 ;;
    --storage) storage=${2%/}; shift 2 ;;
    --max) max=$2; shift 2 ;;
    --setup) setup=$2; shift 2 ;;
    --os-major) osmajor=$2; shift 2 ;;
    --no-submit) submit=no; shift ;;
    --test) test=yes; shift ;;
    -*) echo "unknown option $1" >&2; exit 1 ;;
    *) echo "no arguments besides the options (the materials come from --table): $1" >&2; exit 1 ;;
  esac
done

for v in "$jobs" "$muons" "$fulljobs" "$max"; do
  case "$v" in
    ''|*[!0-9]*) echo "--jobs, --muons, --full-jobs and --max take whole numbers" >&2; exit 1 ;;
  esac
done
if [ "$jobs" -lt 1 ] || [ "$jobs" -gt 999 ] || [ "$muons" -lt 1 ]; then
  echo "--jobs: 1-999 processes per material (the seeds leave room for 999); --muons: at least 1" >&2
  exit 1
fi
if [ "$max" -lt 1 ] || [ "$max" -gt 60 ]; then
  echo "--max: 1-60 jobs at a time (one core each, 60 cores at most)" >&2
  exit 1
fi
case "$osmajor" in
  *[!0-9]*) echo "--os-major: a whole number, e.g. 7" >&2; exit 1 ;;
esac
if [ ! -r "$table" ]; then
  echo "--table: cannot read $table (run this in the build directory)" >&2
  exit 1
fi
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

# the materials: name id composition density ekin p depth rms, tuned (a number in the depth column)
rows=$(awk '!/^#/ && NF >= 8' "$table")
if [ -n "$materials" ]; then
  for m in $(printf '%s' "$materials" | tr ',' ' '); do
    printf '%s\n' "$rows" | awk -v m="$m" '$1 == m { f = 1 } END { exit !f }' ||
      { echo "--materials: $m is not in $table" >&2; exit 1; }
  done
  rows=$(printf '%s\n' "$rows" | awk -v list=",$materials," 'index(list, "," $1 ",") > 0')
fi
untuned=$(printf '%s\n' "$rows" | awk '$7 !~ /^[0-9.]+$/ { printf " %s", $1 }')
if [ -n "$untuned" ]; then
  echo "not tuned in $table (no depth):$untuned; run scripts/muonic_library_depth.sh and copy its materials_tuned.txt" >&2
  exit 1
fi
bad=$(printf '%s\n' "$rows" | awk '$2 !~ /^[0-9]+$/ || $2 > 999 || $3 !~ /^([A-Z][a-z]?:([0-9]+\.?[0-9]*|\.[0-9]+|rest),)*[A-Z][a-z]?:([0-9]+\.?[0-9]*|\.[0-9]+|rest)$/ ||
                                       $4 !~ /^[0-9.]+$/ || $5 !~ /^[0-9.]+$/ { printf " %s", $1 }')
if [ -n "$bad" ]; then
  echo "$table: rows with a bad id, composition, density or ekin:$bad" >&2
  exit 1
fi

# the run's folder
home=${storage:+$storage/}output/muonic_library_czt
[ -n "$dir" ] || dir=$home
basedir=$dir   # the resubmit line names this one (with --test for a test folder)
retest=
if [ "$test" = yes ]; then
  dir=${dir}_condor_test
  retest="--test "
  jobs=2
  fulljobs=1
  muons=10000
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
preinit=macro/preinit_lyso.mac
job=$(cd "$(dirname "$0")" && pwd -P)/beam_steel_job.sh
for f in ./TOFmeter "$beammac" "$preinit" macro/panel_czt.mac "$job"; do
  if [ ! -e "$f" ]; then
    echo "Run this in the build directory ($f is missing; after a new macro, run cmake again)" >&2
    exit 1
  fi
done
if ! grep -a -q "Material005" ./TOFmeter || ! grep -a -q "TOFmeter --slim: the events ntuple" ./TOFmeter ||
   ! grep -a -q "measured line energies (Z 1-42, Th, U)" ./TOFmeter; then
  echo "./TOFmeter lacks /tof/target/composition, --slim or the measured lines of Z 1-42 (2026-11-23): build it again (make) from the current sources" >&2
  exit 1
fi
mkdir -p "$dir"
dir=$(cd "$dir" && pwd -P)
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

# the job macros' setup lines: every macro of a material must have the same ones
setuplines() {   # <composition> <density> <ekin>
  echo "/control/alias bmEkin $3"
  echo "/control/alias bmSigX 1.656"
  echo "/control/alias bmSigY 1.062"
  echo "/control/alias bmZ 33.55"
  echo "/control/alias bmFilter 0"
  echo "/control/alias bmPanel czt"
  echo "/control/alias bmPanelZ 32.0"
  echo "/tof/target/composition $1"
  echo "/tof/target/density $2 g/cm3"
  echo "/control/execute $beammac"
}
printf '%s\n' "$rows" > "$dir/.rows.$$"
while read -r name id comp dens ekin p depth rms note; do
  band=${id}_$name
  [ "$id" -lt 10 ] && band=0${id}_$name
  mkdir -p "$dir/$band/log"
  want=$(setuplines "$comp" "$dens" "$ekin")
  for mac in "$dir/$band"/log/job*.mac; do
    [ -e "$mac" ] || continue
    have=$(grep -E '^/control/alias bm(Ekin|SigX|SigY|Z|Filter|Panel|PanelZ) |^/tof/target/|^/control/execute ' "$mac")
    if [ "$have" != "$want" ]; then
      echo "$dir/$band holds jobs of another setup ($mac):" >&2
      printf '%s\n' "$have" | sed 's/^/  /' >&2
      echo "give another --dir for this one ($name: $comp at $dens g/cm3, $ekin MeV)" >&2
      rm -f "$dir/.rows.$$"
      exit 1
    fi
  done
done < "$dir/.rows.$$"

# the macros of the processes to run
items=
n=0
nfull=0
finished=0
nmat=0
while read -r name id comp dens ekin p depth rms note; do
  band=${id}_$name
  [ "$id" -lt 10 ] && band=0${id}_$name
  sub=$dir/$band
  nmat=$((nmat + 1))
  [ -e "$sub/started" ] || date +%s > "$sub/started"
  k=0
  while [ $k -lt "$jobs" ]; do
    k=$((k + 1))
    if [ -e "$sub/log/job$k.done" ]; then
      finished=$((finished + 1))
      continue
    fi
    kind=slim
    [ $k -le "$fulljobs" ] && kind=full
    {
      echo "# preinit $preinit"   # read by beam_steel_job.sh: TOFmeter --preinit <file>
      [ "$kind" = slim ] && echo "# slim"   # read by beam_steel_job.sh: TOFmeter --slim
      echo "# material $name ($kind): $comp (mass %) at $dens g/cm3, mu- of $ekin MeV ($p MeV/c), tuned depth $depth mm"
      echo "/control/verbose 2"
      echo "/control/alias bmSeed $((12000000 + 1000 * id + k))"
      echo "/control/alias bmEvents $muons"
      echo "/control/alias bmDir $sub"
      echo "/control/alias bmFile beam_$k"
      setuplines "$comp" "$dens" "$ekin"
    } > "$sub/log/job$k.mac"
    items="${items}$band $k
"
    n=$((n + 1))
    [ "$kind" = full ] && nfull=$((nfull + 1))
  done
done < "$dir/.rows.$$"
rm -f "$dir/.rows.$$"
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
  echo "# HTCondor jobs of the muonic X-ray library ($beammac, --preinit $preinit), by condor/muonic_library_submit.sh on $(date '+%F %T')"
  echo "# in $here: $n processes ($nfull full, $((n - nfull)) slim) of $muons muons in $nmat materials of $table ($finished done before), at most $max at a time;"
  echo "# seeds 12000000 + 1000 id + k; environment: the submitting shell's${setup:+, then $setup sourced by each job}"
  echo "# resubmit: sh condor/muonic_library_submit.sh ${retest}${materials:+--materials $materials }--table $table --jobs $jobs --muons $muons --full-jobs $fulljobs --dir $basedir --max $max${setup:+ --setup $setup}${osmajor:+ --os-major $osmajor}"
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
echo "$(date '+%F %T')  muonic X-ray library: $n processes ($nfull full, $((n - nfull)) slim) of $muons muons in $nmat materials ($finished done before), at most $max at a time: $file"
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
echo "check:  sh condor/muonic_library_check.sh $dir   (each job: queued, running, held, finished OK, failed; read-only)"
