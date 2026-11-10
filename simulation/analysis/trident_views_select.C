// Picks the muon trident events to draw for scripts/trident_views_lyso.sh and writes what to draw for each.
//   root -l -b -q 'analysis/trident_views_select.C("output/trident_views_lyso", 8, 8)'
// Input: the searches of macro/vis_tridents_lyso.mac, <dir>/jobs/search_<j>.root and .log for j = 1..jobs, each with
// <kept> kept events, all detectable tridents (the primary after its trident and the produced mu+ and mu- each crossed
// a bar of the top and of the middle TOF). The kept events come from the lines "Kept trident event k/K: event E" of
// the log and are looked up in the tree by eventID; the macro checks with the tree that each is detectable.
// The choice: pictures 1 to nSeparated are the clean events (one trident, no other muon) whose three muons crossed
// the most different bars of the top and middle TOF (then the widest spread of bar numbers, then the earliest); the
// last picture is a typical event, the first clean one whose three muons all crossed one bar of the top and one of
// the middle TOF. The bars with more than `threshold` MeV are highlighted (a muon crossing a 0.6 cm bar leaves about
// 1.2 MeV).
// Writes:
//   <dir>/marks/job<j>_event<k>.mac  one per kept event, executed for it by macro/vis_trident_event_lyso.mac: for a
//                                    chosen one the bars to highlight and two pictures, raw/trident_<n>_full.png (the
//                                    whole setup) and raw/trident_<n>_zoom.png (from the vertex to the middle TOF),
//                                    each also drawn with only the lead walls and steel beams (_tracks.png), where
//                                    the tracks are found; for the others nothing
//   <dir>/redraw.txt                 one line per picture: "<n> <job> <k> <eventID>" (written last)
//   <dir>/trident_views.txt          what each chosen event is, and what the searches took

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {
struct TridentMuon {
  int pdg = 0;
  double ekin = 0.;  // MeV, just after the trident
  int top = -1, mid = -1;  // first bar crossed in the top and middle TOF, -1 if none
};
struct BarHit {
  int station = 0, bar = 0;
  double e = 0., x = 0., y = 0., z = 0.;  // MeV, cm
};
struct TridentEvent {
  int job = 0, k = 0, eventID = 0, pdg = 0;
  double ekin = 0., costh = 1., x0 = 0., y0 = 0., z0 = 0.;  // the primary where it started: MeV, cm
  double vx = 0., vy = 0., vz = 0.;                         // the primary's trident: cm
  double ein = 0., eout = 0., epair = 0., mll = 0., weight = 0.;
  int tridents = 0, otherMuons = 0;
  bool detectable = false;
  std::vector<TridentMuon> muons;  // the primary after its trident, then the produced mu+ and mu-
  std::vector<BarHit> bars;        // the bars above the threshold
  int distinct = 0, spread = 0;    // different bars of the three muons (top + middle), spread of their numbers
  bool Clean() const { return tridents == 1 && otherMuons == 0 && muons.size() == 3; }
};

const char* kStationName[3] = {"top", "middle", "bottom"};
const char* kStationTag[3] = {"Al1", "Al2", "Al4"};  // StationPhys_<tag>, copy number = station
const double kRoofTop = 101., kRoofBottom = 91., kRoofHalf = 45.;  // the 10 cm lead layer, cm

const char* MuonName(int pdg) { return pdg > 0 ? "mu-" : "mu+"; }

std::string Where(const TridentEvent& ev)
{
  const double ax = std::fabs(ev.vx), ay = std::fabs(ev.vy);
  if (ev.vz >= kRoofBottom && ev.vz <= kRoofTop && ax <= kRoofHalf && ay <= kRoofHalf)
    return Form("in the lead roof at (%.1f, %.1f, %.1f) cm, %.1f cm under its top", ev.vx, ev.vy, ev.vz, kRoofTop - ev.vz);
  if (ev.vz <= kRoofTop && std::max(ax, ay) > kRoofHalf && std::max(ax, ay) <= kRoofHalf + 10.)
    return Form("in a lead wall at (%.1f, %.1f, %.1f) cm", ev.vx, ev.vy, ev.vz);
  if (ev.vz >= 90. && ev.vz < kRoofBottom && ax <= 40. && ay <= 40.)
    return Form("in the steel plate at (%.1f, %.1f, %.1f) cm", ev.vx, ev.vy, ev.vz);
  return Form("at (%.1f, %.1f, %.1f) cm, outside the lead", ev.vx, ev.vy, ev.vz);
}

std::string BarList(const TridentEvent& ev, int station)
{
  std::string list;
  for (const auto& b : ev.bars)
    if (b.station == station) list += Form("%s%d (%.1f MeV)", list.empty() ? "" : ", ", b.bar, b.e);
  return list.empty() ? "none" : list;
}
}  // namespace

void trident_views_select(const char* dir = "output/trident_views_lyso", int jobs = 8, int kept = 8,
                          int seed0 = 7100, double threshold = 0.5, int nSeparated = 3)
{
  std::vector<TridentEvent> events;
  Long64_t muons = 0;
  double liveTime = 0.;
  for (int job = 1; job <= jobs; ++job) {
    // the kept events, in order, from the log
    const TString logPath = Form("%s/jobs/search_%d.log", dir, job);
    std::ifstream log(logPath.Data());
    if (!log) {
      printf("cannot read %s\n", logPath.Data());
      return;
    }
    std::map<int, int> keptIndex;  // eventID -> k
    std::string line;
    while (std::getline(log, line)) {
      const auto pos = line.find("Kept trident event ");
      int k = 0, total = 0, id = 0;
      if (pos != std::string::npos && sscanf(line.c_str() + pos, "Kept trident event %d/%d: event %d", &k, &total, &id) == 3)
        keptIndex[id] = k;
    }
    if (static_cast<int>(keptIndex.size()) != kept) {
      printf("%s: %zu kept events, expected %d\n", logPath.Data(), keptIndex.size(), kept);
      return;
    }

    const TString rootPath = Form("%s/jobs/search_%d.root", dir, job);
    TFile file(rootPath);
    auto tree = dynamic_cast<TTree*>(file.Get("events"));
    auto run = dynamic_cast<TTree*>(file.Get("run"));
    if (!tree || !run) {
      printf("%s: no events or run tree\n", rootPath.Data());
      return;
    }
    int nEvents = 0;
    double lt = 0.;
    run->SetBranchAddress("n_events", &nEvents);
    run->SetBranchAddress("live_time", &lt);
    for (Long64_t i = 0; i < run->GetEntries(); ++i) {
      run->GetEntry(i);
      muons += nEvents;
      liveTime += lt;
    }

    TTreeReader reader(tree);
    TTreeReaderValue<int> eventID(reader, "eventID"), pdg(reader, "pdg");
    TTreeReaderValue<int> primSt(reader, "prim_st"), primTop(reader, "prim_top"), primMid(reader, "prim_mid");
    TTreeReaderValue<double> ekin(reader, "ekin"), costh(reader, "costh");
    TTreeReaderValue<double> x0(reader, "x0"), y0(reader, "y0"), z0(reader, "z0");
    TTreeReaderArray<int> barSt(reader, "bar_st"), barId(reader, "bar_id");
    TTreeReaderArray<double> barE(reader, "bar_e"), barX(reader, "bar_x"), barY(reader, "bar_y"), barZ(reader, "bar_z");
    TTreeReaderArray<double> triX(reader, "tri_x"), triY(reader, "tri_y"), triZ(reader, "tri_z");
    TTreeReaderArray<double> triEin(reader, "tri_ein"), triEout(reader, "tri_eout"), triEpair(reader, "tri_epair");
    TTreeReaderArray<double> triMll(reader, "tri_mll"), triW(reader, "tri_w");
    TTreeReaderArray<int> triParent(reader, "tri_parent");
    TTreeReaderArray<int> tmuPdg(reader, "tmu_pdg"), tmuProc(reader, "tmu_proc"), tmuParent(reader, "tmu_parent");
    TTreeReaderArray<int> tmuSt(reader, "tmu_st"), tmuTop(reader, "tmu_top"), tmuMid(reader, "tmu_mid");
    TTreeReaderArray<double> tmuE0(reader, "tmu_e0");
    int found = 0;
    while (reader.Next()) {
      const auto it = keptIndex.find(*eventID);
      if (it == keptIndex.end()) continue;
      ++found;
      TridentEvent ev;
      ev.job = job;
      ev.k = it->second;
      ev.eventID = *eventID;
      ev.pdg = *pdg;
      ev.ekin = *ekin;
      ev.costh = *costh;
      ev.x0 = *x0;
      ev.y0 = *y0;
      ev.z0 = *z0;
      ev.tridents = triX.GetSize();
      int first = -1;  // the primary's (first) trident
      for (std::size_t i = 0; i < triParent.GetSize(); ++i)
        if (triParent[i] == 1) {
          first = static_cast<int>(i);
          break;
        }
      if (first >= 0) {
        ev.vx = triX[first];
        ev.vy = triY[first];
        ev.vz = triZ[first];
        ev.ein = triEin[first];
        ev.eout = triEout[first];
        ev.epair = triEpair[first];
        ev.mll = triMll[first];
        ev.weight = triW[first];
      }
      ev.muons.push_back({*pdg, ev.eout, *primTop, *primMid});
      bool pairCrossed = true;
      int produced = 0;
      for (std::size_t i = 0; i < tmuPdg.GetSize(); ++i) {
        if (tmuProc[i] == 49 && tmuParent[i] == 1) {  // made by a trident of the primary
          ++produced;
          if ((tmuSt[i] & 3) != 3) pairCrossed = false;
          ev.muons.push_back({tmuPdg[i], tmuE0[i], tmuTop[i], tmuMid[i]});
        } else {
          ++ev.otherMuons;
        }
      }
      std::sort(ev.muons.begin() + 1, ev.muons.end(),
                [](const TridentMuon& a, const TridentMuon& b) { return a.pdg < b.pdg; });  // mu+ (-13) first
      ev.detectable = (*primSt & 3) == 3 && first >= 0 && pairCrossed && produced >= 2;
      for (std::size_t i = 0; i < barSt.GetSize(); ++i)
        if (barE[i] > threshold) ev.bars.push_back({barSt[i], barId[i], barE[i], barX[i], barY[i], barZ[i]});
      std::set<int> top, mid;
      for (const auto& m : ev.muons) {
        top.insert(m.top);
        mid.insert(m.mid);
      }
      ev.distinct = static_cast<int>(top.size() + mid.size());
      ev.spread = (*top.rbegin() - *top.begin()) + (*mid.rbegin() - *mid.begin());
      events.push_back(ev);
    }
    if (found != kept) {
      printf("%s: found %d of the %d kept events\n", rootPath.Data(), found, kept);
      return;
    }
  }
  std::sort(events.begin(), events.end(),
            [](const TridentEvent& a, const TridentEvent& b) { return a.job != b.job ? a.job < b.job : a.k < b.k; });
  for (const auto& ev : events)
    if (!ev.detectable) printf("warning: search %d event %d was kept but is not detectable by its columns\n", ev.job, ev.eventID);

  // the choice
  std::vector<int> order(events.size());
  for (std::size_t i = 0; i < order.size(); ++i) order[i] = static_cast<int>(i);
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    const auto &A = events[a], &B = events[b];
    if (A.Clean() != B.Clean()) return A.Clean();
    if (A.distinct != B.distinct) return A.distinct > B.distinct;
    if (A.spread != B.spread) return A.spread > B.spread;
    return a < b;
  });
  std::vector<int> chosen(order.begin(), order.begin() + std::min<std::size_t>(nSeparated, order.size()));
  int typical = -1;
  for (std::size_t i = 0; i < events.size(); ++i)
    if (events[i].Clean() && events[i].distinct == 2 &&
        std::find(chosen.begin(), chosen.end(), static_cast<int>(i)) == chosen.end()) {
      typical = static_cast<int>(i);
      chosen.push_back(typical);
      break;
    }

  // camera of the close-up: the viewpoint of macro/vis_tridents_lyso_viewer.mac, theta 75 deg, phi -62 deg
  const double theta = 75. * TMath::DegToRad(), phi = -62. * TMath::DegToRad();
  const double toCamera[3] = {std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)};
  const double up[3] = {-std::cos(theta) * std::cos(phi), -std::cos(theta) * std::sin(phi), std::sin(theta)};
  const double right[3] = {-std::sin(phi), std::cos(phi), 0.};
  const double sceneRadius = 92.4;  // cm: half the height shown at zoom 1, measured on the pictures of the whole setup
  const double aspect = 4. / 3., margin = 6.;  // picture width / height; cm around the points
  auto dot = [](const double* a, const double* p) { return a[0] * p[0] + a[1] * p[1] + a[2] * p[2]; };

  gSystem->mkdir(Form("%s/marks", dir), true);
  for (std::size_t i = 0; i < events.size(); ++i) {
    const auto& ev = events[i];
    std::ofstream mac(Form("%s/marks/job%d_event%d.mac", dir, ev.job, ev.k));
    const auto pos = std::find(chosen.begin(), chosen.end(), static_cast<int>(i));
    if (pos == chosen.end()) {
      mac << "# search " << ev.job << ", kept event " << ev.k << ", event " << ev.eventID << ": not drawn\n";
      continue;
    }
    const int n = static_cast<int>(pos - chosen.begin()) + 1;
    // the close-up: the vertex, where the primary entered the roof, and the bars hit in the top and middle TOF
    std::vector<std::vector<double>> points = {{ev.vx, ev.vy, ev.vz}};
    if (ev.z0 > kRoofTop && ev.vz < kRoofTop) {
      const double t = (ev.z0 - kRoofTop) / (ev.z0 - ev.vz);
      points.push_back({ev.x0 + t * (ev.vx - ev.x0), ev.y0 + t * (ev.vy - ev.y0), kRoofTop});
    }
    for (const auto& b : ev.bars)
      if (b.station < 2) points.push_back({b.x, b.y, b.z});
    double u0 = 1e9, u1 = -1e9, r0 = 1e9, r1 = -1e9, depth = 0.;
    for (const auto& p : points) {
      u0 = std::min(u0, dot(up, p.data()));
      u1 = std::max(u1, dot(up, p.data()));
      r0 = std::min(r0, dot(right, p.data()));
      r1 = std::max(r1, dot(right, p.data()));
      depth += dot(toCamera, p.data()) / points.size();
    }
    const double half = std::max(0.5 * (u1 - u0), 0.5 * (r1 - r0) / aspect) + margin;
    const double zoom = std::min(6., std::max(1.6, sceneRadius / half));
    const double cu = 0.5 * (u0 + u1), cr = 0.5 * (r0 + r1);
    double target[3];
    for (int a = 0; a < 3; ++a) target[a] = cu * up[a] + cr * right[a] + depth * toCamera[a];

    const TString scratch = Form("%s/raw/scratch_%d.png", dir, ev.job);
    mac << "# picture " << n << ": search " << ev.job << ", kept event " << ev.k << ", event " << ev.eventID
        << " (written by analysis/trident_views_select.C)\n";
    mac << "/vis/tsg/offscreen/set/file " << scratch << "\n";
    mac << "/vis/viewer/clearVisAttributesModifiers\n";
    for (const auto& b : ev.bars) {
      mac << "/vis/set/touchable WorldPhys 0 StationPhys_" << kStationTag[b.station] << " " << b.station
          << " ScintillatorPhys_" << kStationTag[b.station] << "_" << b.bar + 1 << " " << b.bar << "\n";
      mac << "/vis/touchable/set/colour 1 0.85 0 1\n";
    }
    // each view twice: as it is, and with only the lead walls and the steel beams, from which scripts/thicken_tracks.py
    // takes the tracks, so that they show over the roof (whose see-through faces the renderer does not draw lines
    // behind), the plates and the bars (which hide the tracks below them for tens of cm, seen from 15 deg above)
    const char* hidden[] = {"AbsorberLogic", "FrameworkFeLogic", "FrameworkAlLogic", "ScintillatorLogicX",
                            "ScintillatorLogicY", "TargetLogic", "LYSOPixelLogic"};
    auto picture = [&](const char* view) {
      mac << "/vis/tsg/offscreen/set/file " << dir << "/raw/trident_" << n << "_" << view << ".png\n/vis/viewer/rebuild\n";
      mac << "/vis/tsg/offscreen/set/file " << scratch << "\n";
      for (const char* lv : hidden) mac << "/vis/geometry/set/visibility " << lv << " 0 false\n";
      mac << "/vis/tsg/offscreen/set/file " << dir << "/raw/trident_" << n << "_" << view << "_tracks.png\n/vis/viewer/rebuild\n";
      mac << "/vis/tsg/offscreen/set/file " << scratch << "\n";
      for (const char* lv : hidden) mac << "/vis/geometry/set/visibility " << lv << " 0 true\n";
    };
    mac << "/vis/viewer/set/targetPoint 0 0 50 cm\n/vis/viewer/zoomTo 1.3\n";
    picture("full");
    mac << Form("/vis/viewer/set/targetPoint %.2f %.2f %.2f cm\n/vis/viewer/zoomTo %.3f\n", target[0], target[1], target[2], zoom);
    picture("zoom");
    mac << "/vis/viewer/set/targetPoint 0 0 50 cm\n/vis/viewer/zoomTo 1.3\n/vis/viewer/clearVisAttributesModifiers\n";
  }

  // what the pictures show
  const double factor = events.empty() || events[0].weight <= 0. ? 0. : 1. / events[0].weight;
  int six = 0, two = 0;
  for (const auto& ev : events) {
    if (ev.distinct == 6) ++six;
    if (ev.distinct == 2) ++two;
  }
  const double perYear = liveTime > 0. && factor > 0. ? events.size() / factor / liveTime * 3.15576e7 : 0.;
  FILE* txt = fopen(Form("%s/trident_views.txt", dir), "w");
  fprintf(txt, "Muon trident events for slides: scripts/trident_views_lyso.sh (macro/vis_tridents_lyso.mac, "
               "analysis/trident_views_select.C)\n");
  fprintf(txt, "Setup: the LYSO setup of macro/vis_lyso.mac (Si/YSO stack out, 4 cm Al target, LYSO panel 12.69 x 12.69 cm "
               "under it), 10 cm lead roof and walls.\n");
  fprintf(txt, "Muons: sea-level cosmic mu+ and mu- above 0.85 GeV (EXPACS spectrum, Guan zenith law) from a 2 x 2 m plane at "
               "z = 112 cm, zenith < 75 deg.\n");
  fprintf(txt, "Searches: %d (seeds %d to %d), each until its %d-th detectable trident (the primary after its trident and the "
               "produced mu+ and mu- all cross a bar of the top and of the middle TOF), tridents made %.0g times more often.\n",
          jobs, seed0 + 1, seed0 + jobs, kept, factor);
  fprintf(txt, "  %lld muons = %.2f h of sea-level exposure for %zu detectable tridents: without the biasing about %.1f per "
               "year, one every %.0f days\n",
          muons, liveTime / 3600., events.size(), perYear, perYear > 0. ? 365.25 / perYear : 0.);
  fprintf(txt, "  (at this biasing the rates come out ~5-6%% low, see README; analysis/trident_rate.C gives them properly)\n");
  fprintf(txt, "  three different bars of both the top and the middle TOF: %d of %zu; all three in one bar of each: %d of %zu\n",
          six, events.size(), two, events.size());
  fprintf(txt, "Pictures: only the muon tracks are drawn, mu- red and mu+ blue, over the roof, plates and bars they pass "
               "(hidden only behind the lead walls); the TOF bars with more than %.1f MeV are yellow; the lead roof is drawn "
               "see-through; the setup is cut open at y = -21 cm, in front of the TOF bars.\n",
          threshold);
  FILE* redraw = nullptr;
  std::vector<std::string> lines;
  for (std::size_t c = 0; c < chosen.size(); ++c) {
    const auto& ev = events[chosen[c]];
    const int n = static_cast<int>(c) + 1;
    const char* kind = chosen[c] == typical ? "a typical detectable trident: the three muons in one bar of each TOF"
                                            : "one of the most separated detectable tridents";
    fprintf(txt, "\ntrident_%d_full.png, trident_%d_zoom.png: %s\n", n, n, kind);
    fprintf(txt, "  search %d (seed %d), event %d (kept event %d)\n", ev.job, seed0 + ev.job, ev.eventID, ev.k);
    fprintf(txt, "  incoming %s: %.3g GeV at the source plane, zenith %.1f deg; %.3g GeV at the trident\n", MuonName(ev.pdg),
            ev.ekin / 1000., std::acos(std::min(1., ev.costh)) * TMath::RadToDeg(), ev.ein / 1000.);
    fprintf(txt, "  trident %s\n", Where(ev).c_str());
    if (ev.muons.size() >= 3)
      fprintf(txt, "  after it: the %s %.3g GeV, the produced %s %.3g GeV and %s %.3g GeV (pair %.3g GeV in total), "
                   "like-sign mass %.0f MeV\n",
              MuonName(ev.muons[0].pdg), ev.muons[0].ekin / 1000., MuonName(ev.muons[1].pdg), ev.muons[1].ekin / 1000.,
              MuonName(ev.muons[2].pdg), ev.muons[2].ekin / 1000., ev.epair / 1000., ev.mll);
    std::string top, mid;
    for (const auto& m : ev.muons) {
      top += Form(" %d", m.top);
      mid += Form(" %d", m.mid);
    }
    fprintf(txt, "  first bar of each muon (primary, mu+, mu-): top TOF%s, middle TOF%s (%d different bars in all)\n",
            top.c_str(), mid.c_str(), ev.distinct);
    for (int s = 0; s < 3; ++s) fprintf(txt, "  %s TOF bars above %.1f MeV: %s\n", kStationName[s], threshold, BarList(ev, s).c_str());
    if (!ev.Clean()) fprintf(txt, "  (%d tridents, %d other muons in the event)\n", ev.tridents, ev.otherMuons);
    lines.push_back(Form("%d %d %d %d", n, ev.job, ev.k, ev.eventID));
  }
  fclose(txt);
  redraw = fopen(Form("%s/redraw.txt", dir), "w");
  for (const auto& l : lines) fprintf(redraw, "%s\n", l.c_str());
  fclose(redraw);
  printf("%zu events chosen of %zu (%d with three different bars in both TOFs): %s/trident_views.txt, %s/redraw.txt\n",
         chosen.size(), events.size(), six, dir, dir);
}
