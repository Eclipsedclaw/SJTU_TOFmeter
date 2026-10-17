#!/bin/sh
# Part of the trident scan (macro/trident_scan_def.mac) on another computer: builds TOFmeter there and
# runs scripts/trident_scan.sh for some of the thicknesses, with the same settings and seeds as the
# full scan, so the files can join those of the first computer. Run it in the unpacked simulation
# folder (the one with CMakeLists.txt), with Geant4 set up (geant4-config in PATH, datasets found):
#   sh scripts/trident_remote.sh ["40 50"] [jobs at a time, default all cores]
# Needs Geant4 (11.4.3 as on the first computer, for the same physics), CMake and a C++17 compiler;
# ROOT is not needed (trident_rate.C runs later on the first computer). Builds in ./build-remote and
# runs there in the background: build-remote/output/tridents_F1e3/t<cm>/, progress in
# build-remote/output/tridents_F1e3/scan.log; stop after the running processes with
# touch build-remote/output/tridents_F1e3/STOP, go on by starting this script again.
# Afterwards, on the first computer (in its build directory):
#   scp -r <user>@<host>:<this folder>/build-remote/output/tridents_F1e3/t40 output/tridents_F1e3/

leads=${1:-40 50}
jobs=$2
if [ -z "$jobs" ]; then
  jobs=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
fi
if [ ! -f CMakeLists.txt ] || [ ! -f macro/trident_scan_def.mac ]; then
  echo "Run this in the simulation folder (CMakeLists.txt and macro/ must be there)" >&2
  exit 1
fi
if ! command -v geant4-config > /dev/null; then
  echo "geant4-config not found: set up Geant4 first (source <geant4-install>/bin/geant4.sh)" >&2
  exit 1
fi
version=$(geant4-config --version)
echo "Geant4 $version"
[ "$version" = "11.4.3" ] || echo "warning: the first computer used Geant4 11.4.3; another version can give slightly different rates"
if [ -z "$GEANT4_DATA_DIR" ] && geant4-config --check-datasets 2>&1 | grep -q NOTFOUND; then
  echo "Geant4 datasets not found: install them or set GEANT4_DATA_DIR" >&2
  exit 1
fi

cmake -S . -B build-remote -DCMAKE_BUILD_TYPE=Release > build-remote.log 2>&1 \
  && cmake --build build-remote -j "$jobs" >> build-remote.log 2>&1
if [ ! -x build-remote/TOFmeter ]; then
  echo "build failed, see build-remote.log" >&2
  exit 1
fi
cd build-remote || exit 1

# the scan definition with only these thicknesses; settings and seeds otherwise unchanged
def=macro/trident_scan_def_remote.mac
sed "s|^/control/alias tridentScanLead .*|/control/alias tridentScanLead \"$leads\"|" macro/trident_scan_def.mac > "$def"
dir=$(sed -n 's|^/control/alias tridentScanDir *||p' "$def")
mkdir -p "$dir"
awake=""
command -v caffeinate > /dev/null && awake="caffeinate -is"
nohup nice -n 15 $awake sh scripts/trident_scan.sh "$def" "$jobs" >> "$dir/scan.log" 2>&1 &
sleep 3
cat "$dir/scan.log"
echo "running in the background, $jobs at a time: tail -f $(pwd)/$dir/scan.log"
echo "when done, on the first computer, in its build directory:"
for t in $leads; do echo "  scp -r $(whoami)@$(hostname):$(pwd)/$dir/t$t output/tridents_F1e3/"; done
