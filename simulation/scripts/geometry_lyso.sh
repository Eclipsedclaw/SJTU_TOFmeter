#!/bin/sh
# Pictures of the LYSO setup for slides, one per file, in output/geometry_lyso/ (as scripts/geometry_czt.sh for the CZT
# plane): three views drawn by Geant4 (macro/vis_lyso.mac; scripts/picture_finish.py shrinks them for smooth edges, cuts
# the white border and adds a legend) and three drawings to scale from the table of placed volumes the same run prints
# (analysis/geometry_schematic.C):
#   geometry_3d_setup.png        the whole setup cut open at y = 1.0 cm, the LYSO panel in the cut under the target
#   geometry_3d_target.png       the target and the LYSO panel from below, with the middle TOF
#   geometry_3d_panel.png        a corner of the LYSO panel close up from below, under the target: pixels and gaps
#   geometry_side.pdf/.png       section through the setup at y = 1.0 cm, labelled
#   geometry_side_zoom.pdf/.png  the same around the target, with its dimensions
#   geometry_top.pdf/.png        seen from above
# geometry.log: the Geant4 output with the table; raw/: the pictures as Geant4 drew them. One process at a time.
# In the build directory, with Geant4 and ROOT set up:
#   sh scripts/geometry_lyso.sh
out=output/geometry_lyso
for f in ./TOFmeter macro/vis_lyso.mac analysis/geometry_schematic.C scripts/picture_finish.py; do
  if [ ! -e "$f" ]; then
    echo "Run this in the build directory ($f is missing)" >&2
    exit 1
  fi
done
mkdir -p "$out/raw"
if ! ./TOFmeter macro/vis_lyso.mac > "$out/geometry.log" 2>&1; then
  echo "TOFmeter failed: see $out/geometry.log" >&2
  exit 1
fi
lead='Pb=#5a5f66' steel='steel=#b5bac2' bars='plastic TOF bars=#4db050' target='Al target=#f29c12' lyso='LYSO panel=#7a40a1'
python3 -I scripts/picture_finish.py "$out/raw/setup.png" "$out/geometry_3d_setup.png" \
  "$lead" "$steel" 'Al plates=#a8cce8' "$bars" "$target" "$lyso" || exit 1
python3 -I scripts/picture_finish.py "$out/raw/target.png" "$out/geometry_3d_target.png" \
  'plastic bars of the middle TOF=#4db050' "$target" "$lyso" || exit 1
python3 -I scripts/picture_finish.py "$out/raw/panel.png" "$out/geometry_3d_panel.png" \
  'LYSO pixels, 6.25 x 6.25 x 3 mm, 0.08 mm ESR gaps=#7a40a1' "$target" || exit 1
ROOT_HIST=0 root -l -b -q "analysis/geometry_schematic.C(\"$out/geometry.log\", \"$out\", 1.0)"
