#!/bin/sh
# Pictures of the CZT setup for slides, one per file, in output/geometry_czt/: two views drawn by Geant4
# (macro/vis_czt.mac; scripts/picture_finish.py shrinks them for smooth edges, cuts the white border and adds a legend)
# and three drawings to scale from the table of placed volumes the same run prints (analysis/geometry_schematic.C):
#   geometry_3d_setup.png        the whole setup cut open at y = 1.25 cm, the CZT plane in the cut under the target
#   geometry_3d_target.png       the target and the CZT plane from below, with the middle TOF
#   geometry_side.pdf/.png       section through the setup at y = 1.25 cm, labelled
#   geometry_side_zoom.pdf/.png  the same around the target, with its dimensions
#   geometry_top.pdf/.png        seen from above
# geometry.log: the Geant4 output with the table; raw/: the pictures as Geant4 drew them. One process at a time.
# In the build directory, with Geant4 and ROOT set up:
#   sh scripts/geometry_czt.sh
out=output/geometry_czt
for f in ./TOFmeter macro/vis_czt.mac analysis/geometry_schematic.C scripts/picture_finish.py; do
  if [ ! -e "$f" ]; then
    echo "Run this in the build directory ($f is missing)" >&2
    exit 1
  fi
done
mkdir -p "$out/raw"
if ! ./TOFmeter macro/vis_czt.mac > "$out/geometry.log" 2>&1; then
  echo "TOFmeter failed: see $out/geometry.log" >&2
  exit 1
fi
lead='Pb=#5a5f66' steel='steel=#b5bac2' bars='plastic TOF bars=#4db050' target='Al target=#f29c12' czt='CZT plane=#7a40a1'
python3 -I scripts/picture_finish.py "$out/raw/setup.png" "$out/geometry_3d_setup.png" \
  "$lead" "$steel" 'Al plates=#a8cce8' "$bars" "$target" "$czt" || exit 1
python3 -I scripts/picture_finish.py "$out/raw/target.png" "$out/geometry_3d_target.png" \
  'plastic bars of the middle TOF=#4db050' "$target" "$czt" || exit 1
root -l -b -q "analysis/geometry_schematic.C(\"$out/geometry.log\", \"$out\")"
