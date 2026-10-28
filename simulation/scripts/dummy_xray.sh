#!/bin/sh
# Dummy muonic X-ray spectra, to check the atomic cascade: mu- of 1 keV stopped directly in a target
# of each material (macro/dummy_xray.mac), one TOFmeter process per material and cascade model, all at
# the same time (one core each, at most 90% of the physical cores at a time), then analysis/dummy_xray_compare.C overlays the
# emitted X-ray yields of all of them. Run it in the build directory:
#   sh scripts/dummy_xray.sh [--dir path] [--materials "G4_Al G4_Si"] [--models "firstPrinciples"]
#                            [--range keV] [--ymax top] [muons per run]
# Defaults: output/dummy_xray, G4_Al and G4_Si, the first-principles cascade (add geant4 for Geant4's
# G4EmCaptureCascade), 0-120 keV, y axis 0-1, 100000 muons per run. Writes <dir>/<material>_<model>/mux.root,
# with the macro and log of each process in <dir>/log/, and the overlay
# <dir>/dummy_xray_compare.png/.pdf/.txt with _lines.csv and _bins.csv. After a change of the cascade
# code: make, then run this again (the files of the last run are replaced).

dir=output/dummy_xray; materials="G4_Al G4_Si"; models=firstPrinciples; range=120; ymax=1
while [ $# -gt 0 ]; do
  case "$1" in
    --dir) dir=$2; shift 2 ;;
    --materials) materials=$2; shift 2 ;;
    --models) models=$2; shift 2 ;;
    --range) range=$2; shift 2 ;;
    --ymax) ymax=$2; shift 2 ;;
    *) break ;;
  esac
done
events=${1:-100000}
for f in ./TOFmeter macro/dummy_xray.mac macro/geom_target_al.mac analysis/dummy_xray_compare.C; do
  if [ ! -e "$f" ]; then
    echo "Run this in the build directory ($f is missing)" >&2
    exit 1
  fi
done
# the ROOT of /opt/root, not the snap one, which writes to ~/snap (outside the project)
case "$(command -v root)" in
  "") echo "root not found: source /opt/root/root_install/bin/thisroot.sh first" >&2
      exit 1 ;;
  /snap/*) echo "root is the snap ($(command -v root)), which writes to ~/snap:" \
                "source /opt/root/root_install/bin/thisroot.sh first" >&2
           exit 1 ;;
esac
mkdir -p "$dir/log"
runs=""
seed=9100001
for material in $materials; do
  for model in $models; do
    run=${material}_$model
    mkdir -p "$dir/$run"
    {
      echo "/control/verbose 2"
      echo "/control/alias dummyMaterial $material"
      echo "/control/alias dummyModel $model"
      echo "/control/alias dummySeed $seed"
      echo "/control/alias dummyEvents $events"
      echo "/control/alias dummyDir $dir/$run"
      echo "/control/alias dummyFile mux"
      echo "/control/execute macro/dummy_xray.mac"
    } > "$dir/log/$run.mac"
    runs="$runs $run"
    seed=$((seed + 1))
  done
done
parallel=$(echo $runs | wc -w)
limit=$(( $(lscpu -p=Core,Socket | grep -v '^#' | sort -u | wc -l) * 9 / 10 ))   # 90% of the physical cores
[ "$limit" -ge 1 ] || limit=1
[ "$parallel" -gt "$limit" ] && parallel=$limit

echo "$(date '+%F %T')  $events mu- of 1 keV per run, $parallel at a time ->$dir:$runs"
export TOF_DUMMY_DIR="$dir"
if ! echo $runs | tr ' ' '\n' | xargs -P "$parallel" -n 1 sh -c '
  if ./TOFmeter "$TOF_DUMMY_DIR/log/$1.mac" > "$TOF_DUMMY_DIR/log/$1.log" 2>&1; then
    echo "$(date "+%F %T")  $1 done"
  else
    echo "$(date "+%F %T")  $1 FAILED, see $TOF_DUMMY_DIR/log/$1.log"
    exit 1
  fi' sh; then
  exit 1
fi

root -l -b -q "analysis/dummy_xray_compare.C(\"$dir\", $range, 1, $ymax)" > "$dir/dummy_xray_compare.txt" 2>&1
cat "$dir/dummy_xray_compare.txt"
