#!/bin/sh
# Pictures of muon trident events in the LYSO setup for slides, one per file, in output/trident_views_lyso/:
#   trident_<n>_full.png   event n in the whole setup, cut open in front of the TOF stations
#   trident_<n>_zoom.png   the same event closer: from the vertex in the lead roof to the middle TOF
#   trident_views.txt      what each event is (energies, vertex, bars hit) and what the searches took
# Pictures 1-3 are the detectable tridents whose three muons are the most separated (the most different bars of the
# top and middle TOF), picture 4 a typical one (the three muons in one bar of each TOF). Only the muon tracks are
# drawn (mu- red, mu+ blue); the TOF bars with more than 0.5 MeV are yellow; the lead roof is see-through.
# How: <jobs> searches with different seeds (7100 + job) find <kept> detectable tridents each, without drawing
# (macro/vis_tridents_lyso.mac, tridents made 10^5 times more often); analysis/trident_views_select.C picks the events
# and writes what to draw for each; the searches with chosen events run again with the offscreen viewer
# (macro/vis_tridents_lyso_viewer.mac; the same seed gives the same events, which the script checks) and draw them
# (macro/vis_trident_event_lyso.mac); scripts/thicken_tracks.py and scripts/picture_finish.py finish the pictures.
# jobs/ keeps the job macros, logs and ROOT files, raw/ the pictures as Geant4 drew them, marks/ the drawing files.
# In the build directory, with Geant4, ROOT, numpy and Pillow set up:
#   sh scripts/trident_views_lyso.sh [jobs, default 8] [kept per job, default 8]
# The jobs run in parallel, one thread each, at most 90% of the physical cores; 8 x 8 takes about 3 min per pass.
out=output/trident_views_lyso
jobs=${1:-8}
kept=${2:-8}
seed0=7100
limit=$(( $(lscpu -p=Core,Socket | grep -v '^#' | sort -u | wc -l) * 9 / 10 ))
if [ "$jobs" -lt 1 ] || [ "$jobs" -gt "$limit" ]; then
  echo "jobs must be 1 to $limit (90% of the physical cores)" >&2
  exit 1
fi
for f in ./TOFmeter macro/physics_tridents.mac macro/panel_lyso.mac macro/vis_tridents_lyso.mac \
         macro/vis_tridents_lyso_viewer.mac macro/vis_trident_event_lyso.mac analysis/trident_views_select.C \
         scripts/thicken_tracks.py scripts/picture_finish.py; do
  if [ ! -e "$f" ]; then
    echo "Run this in the build directory ($f is missing)" >&2
    exit 1
  fi
done
mkdir -p "$out/jobs" "$out/raw" "$out/marks"
rm -f "$out/redraw.txt"

# job macro: the aliases of macro/vis_tridents_lyso.mac, and for a drawing run the review of the kept events
job_macro() {  # <job> <draw 0|1> <name>
  {
    echo "/control/alias tvSeed $((seed0 + $1))"
    echo "/control/alias tvKept $kept"
    echo "/control/alias tvDir $out"
    echo "/control/alias tvFile $3"
    echo "/control/alias tvJob $1"
    echo "/control/alias tvDraw $2"
    echo "/control/execute macro/vis_tridents_lyso.mac"
    if [ "$2" = 1 ]; then
      echo "/vis/enable"
      echo "/control/alias tvEvent 0"
      echo "/vis/tsg/offscreen/set/file $out/raw/scratch_$1.png"
      echo "/vis/reviewKeptEvents macro/vis_trident_event_lyso.mac"
    fi
  } > "$out/jobs/$3.mac"
}

run_jobs() {  # <search|draw> <job> ...
  kind=$1
  shift
  draw=0
  [ "$kind" = draw ] && draw=1
  pids=""
  for j in "$@"; do
    job_macro "$j" "$draw" "${kind}_$j"
    OMP_NUM_THREADS=1 ./TOFmeter --preinit macro/physics_tridents.mac "$out/jobs/${kind}_$j.mac" \
      > "$out/jobs/${kind}_$j.log" 2>&1 &
    pids="$pids $!"
  done
  status=0
  for p in $pids; do
    wait "$p" || status=1
  done
  for j in "$@"; do
    if ! grep -q "Kept trident event $kept/$kept:" "$out/jobs/${kind}_$j.log"; then
      echo "$kind $j did not keep $kept tridents: see $out/jobs/${kind}_$j.log" >&2
      status=1
    fi
  done
  return $status
}

start=$(date +%s)
echo "Searches: $jobs x $kept detectable tridents ($out/jobs/search_*.log)"
run_jobs search $(seq 1 "$jobs") || exit 1
echo "  done in $(( $(date +%s) - start )) s"
ROOT_HIST=0 root -l -b -q "analysis/trident_views_select.C(\"$out\", $jobs, $kept, $seed0)" || exit 1
if [ ! -s "$out/redraw.txt" ]; then
  echo "no events chosen: see the output of analysis/trident_views_select.C above" >&2
  exit 1
fi
redraw=$(awk '{print $2}' "$out/redraw.txt" | sort -un | tr '\n' ' ')
echo "Drawing: searches $redraw($out/jobs/draw_*.log)"
run_jobs draw $redraw || exit 1
for j in $redraw; do
  grep 'Kept trident event' "$out/jobs/search_$j.log" > "$out/jobs/kept_search_$j.txt"
  grep 'Kept trident event' "$out/jobs/draw_$j.log" > "$out/jobs/kept_draw_$j.txt"
  if ! cmp -s "$out/jobs/kept_search_$j.txt" "$out/jobs/kept_draw_$j.txt"; then
    echo "search $j: the drawing run kept other events than the search (compare $out/jobs/kept_*_$j.txt)" >&2
    exit 1
  fi
done
echo "  done in $(( $(date +%s) - start )) s"

mum='#e31a1c' mup='#1f59f2'
while read -r n j k e; do
  for view in full zoom; do
    raw="$out/raw/trident_${n}_${view}.png"
    for f in "$raw" "$out/raw/trident_${n}_${view}_tracks.png"; do
      if [ ! -s "$f" ]; then
        echo "$f was not drawn: see $out/jobs/draw_$j.log" >&2
        exit 1
      fi
    done
    python3 -I scripts/thicken_tracks.py "$raw" "$out/raw/trident_${n}_${view}_wide.png" \
      --mask-from "$out/raw/trident_${n}_${view}_tracks.png" "$mum" "$mup" > /dev/null || exit 1
    python3 -I scripts/picture_finish.py "$out/raw/trident_${n}_${view}_wide.png" "$out/trident_${n}_${view}.png" \
      "μ⁻ track=$mum" "μ⁺ track=$mup" 'TOF bar hit (> 0.5 MeV)=#ffd900' 'plastic TOF bars=#4db050' \
      'Pb (roof see-through)=#5a5f66' 'steel=#b5bac2' 'Al plates=#a8cce8' 'Al target=#f29c12' || exit 1
  done
done < "$out/redraw.txt"
echo "Pictures in $out (trident_views.txt says what they show), $(( $(date +%s) - start )) s in all"
