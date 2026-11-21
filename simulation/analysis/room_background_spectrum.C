// The room background of the steel study (2026-11-20: can a TOF-trigger ON/OFF subtraction show U, Th and K at 50 ppm
// in steel, with the room's background included?), per day on the pixel panel, from the jobs of
// condor/room_background_submit.sh in one folder:
//   <dir>/floor   the 30 cm concrete floor's natural decays (U-238, Th-232, K-40, U-235, each with its chain) in five
//                 regions: inner (x and y within +-55 cm, under the lead castle) and the ring around it out to +-100 cm
//                 (px, nx, py, ny);
//   <dir>/radon   Rn-222 with its short-lived daughters in the air inside the castle;
//   <dir>/cosmic  sea-level cosmic rays, every EXPACS species, open sky;
// and, if linesDir is given, local runs of macro/room_lines.mac (lines_<keV>.root, the energy with two decimals): the
// walls' and ceiling's strongest gamma lines through the castle's 10 cm of lead, an estimate.
//
// Normalisation, per day (seconds):
//   floor:  each region's mass (its box x 2.3 g/cm3, G4_CONCRETE) x the nuclide's activity per kg (arguments: K-40 400,
//           U-238 40, Th-232 30 Bq/kg, EU RP-112's typical concrete; U-235 the natural share of U-238's, 568.5 / 12355,
//           as beam_steel_trace_spectrum.C), each chain in secular equilibrium (one event = one decay of every
//           member); an event of a region weighs A T / N, N that region's simulated decays of that nuclide;
//   radon:  radonBqM3 x the castle's air, 0.7147 m3 (the 90 x 90 x 91 cm inside, 0.7371 m3, minus the Al plates
//           5952 cm3, the steel frame plate 6400, the four L beams 6916, the three bar stations 2592, the plate 400 and
//           the panel 100 cm3 in it), with its short-lived daughters in equilibrium (equilibrium factor 1);
//   cosmic: T / the simulated live time (the run trees' live_time);
//   walls (estimate): a line of energy E and yield y per decay of a chain of activity a (Bq/g) fills a closed concrete
//           room with an isotropic uncollided flux phi = a y / (mu/rho)(E) photons /cm2/s (the flux inside a cavity in
//           a thick, uniform emitter; mu/rho of NIST's ordinary concrete); the sphere of macro/room_lines.mac (95 cm
//           around the castle, 86850 cm2 of it above the floor) emits phi / 4 per cm2 and second inward. The walls'
//           scattered photons are left out (lower energies, which hardly cross 10 cm of lead): a lower estimate.
// Entries: each pixel hit (czt_e) is smeared by sigma1MeV x sqrt(E / MeV) and kept from threshold (20 keV) up, as in
// beam_steel_trace_spectrum.C (with its own random seed).
// Outputs in outDir (default <dir>/analysis): room_background_spectrum.root (per day in 0.5 keV bins, the binning of
// beam_steel_trace_spectrum.root: room_floor_<nuclide>_<region>, room_floor_<nuclide>, room_floor, room_radon,
// room_cosmic, room_walls, room_total; nuclides U238, Th232, K40, U235), .txt and .csv (the numbers), and one plot,
// room_background_spectrum (the parts per day, 20-7000 keV). Only reads the runs; writes into outDir. In the work folder:
//   root -l -b -q 'analysis/room_background_spectrum.C("output/czt/beam_steel_czt_room", "", "output/czt/beam_steel_czt_room_lines")'
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
const double kFillKeV = 0.5, kMaxKeV = 8000.;   // beam_steel_trace_spectrum.root's bins
const int kBins = int(kMaxKeV / kFillKeV);
const double kConcrete = 2.3;                    // g/cm3, G4_CONCRETE
const double kCastleAirM3 = 0.7147;              // the castle's air (see the head)
const double kU235PerU238 = 568.5 / 12355.;      // activity ratio in natural U
const double kSphereArea = 2. * M_PI * 95. * (95. + 50.5);   // cm2 of macro/room_lines.mac's sphere above the floor

void Save(TCanvas* c, const std::string& base)
{
  c->Print((base + ".png").c_str());
  c->Print((base + ".pdf").c_str());
}

TH1D* NewHist(const std::string& name)
{
  auto h = new TH1D(name.c_str(), "", kBins, 0., kMaxKeV);
  h->Sumw2();
  h->SetDirectory(nullptr);
  return h;
}

TH1D* Copy(const TH1D* h, const std::string& name, double scale, int rebin = 1)
{
  auto c = dynamic_cast<TH1D*>(h->Clone(name.c_str()));
  c->SetDirectory(nullptr);
  c->Scale(scale);
  if (rebin > 1) c->Rebin(rebin);
  return c;
}

std::string Sci(double v, int digits = 2)
{
  if (v == 0.) return "0";
  const int e = int(std::floor(std::log10(std::fabs(v))));
  char s[64];
  snprintf(s, sizeof s, "%.*f#times10^{%d}", digits - 1, v / std::pow(10., e), e);
  return s;
}

double WindowSum(const TH1D* h, double lo, double hi)
{
  double sum = 0.;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    const double x = h->GetBinCenter(b);
    if (x >= lo && x <= hi) sum += h->GetBinContent(b);
  }
  return sum;
}

std::vector<int> FinishedJobs(const std::string& bandDir, const char* prefix)
{
  std::vector<int> ks;
  if (void* d = gSystem->OpenDirectory(bandDir.c_str())) {
    while (const char* entry = gSystem->GetDirEntry(d)) {
      int k = 0;
      char tail[8] = "";
      if (sscanf(entry, (std::string(prefix) + "_%d.%5s").c_str(), &k, tail) == 2 && std::string(tail) == "root" &&
          !gSystem->AccessPathName(Form("%s/log/job%d.done", bandDir.c_str(), k)))
        ks.push_back(k);
    }
    gSystem->FreeDirectory(d);
  }
  std::sort(ks.begin(), ks.end());
  return ks;
}

std::string MacroValue(const std::string& file, const std::string& key)   // a job macro's line "key value"
{
  std::ifstream in(file);
  std::string l;
  while (std::getline(in, l))
    if (l.rfind(key + " ", 0) == 0) return l.substr(key.size() + 1);
  return "";
}

double MuRhoConcrete(double keV)   // cm2/g, NIST XCOM "Concrete, Ordinary" (= G4_CONCRETE), log-log between its points
{
  static const double e[] = {500., 600., 800., 1000., 1250., 1500., 2000., 3000.};
  static const double m[] = {0.08712, 0.08058, 0.07076, 0.06361, 0.05690, 0.05183, 0.04480, 0.03668};
  int i = 0;
  while (i < 6 && keV > e[i + 1]) ++i;
  const double f = std::log(keV / e[i]) / std::log(e[i + 1] / e[i]);
  return m[i] * std::pow(m[i + 1] / m[i], f);
}

struct WallLine {   // the room's strongest lines from 0.58 MeV up (yields per decay of the chain's parent)
  double keV, yield;
  int chain;   // 0 U-238, 1 Th-232, 2 K-40
  const char* nuclide;
};
const std::vector<WallLine> kWallLines = {
    {583.19, 0.3055, 1, "Tl-208"},  {609.31, 0.4549, 0, "Bi-214"},  {911.20, 0.258, 1, "Ac-228"},
    {968.97, 0.158, 1, "Ac-228"},   {1120.29, 0.1491, 0, "Bi-214"}, {1238.11, 0.05831, 0, "Bi-214"},
    {1377.67, 0.03968, 0, "Bi-214"}, {1460.82, 0.1066, 2, "K-40"},  {1588.19, 0.0322, 1, "Ac-228"},
    {1729.60, 0.02878, 0, "Bi-214"}, {1764.49, 0.1531, 0, "Bi-214"}, {1847.42, 0.02025, 0, "Bi-214"},
    {2204.21, 0.04913, 0, "Bi-214"}, {2614.51, 0.3585, 1, "Tl-208"}};
}  // namespace

void room_background_spectrum(const char* dir, const char* outDir = "", const char* linesDir = "", double kBqKg = 400.,
                              double uBqKg = 40., double thBqKg = 30., double radonBqM3 = 50., double seconds = 86400.,
                              double sigma1MeV = 6.36991, double threshold = 20.)
{
  gStyle->SetOptStat(0);
  const std::string base = dir;
  const std::string out = std::string(outDir).empty() ? base + "/analysis" : std::string(outDir);
  gSystem->mkdir(out.c_str(), true);
  std::ostringstream txt, csv;
  char line[600];
  TRandom3 random(4359);
  // the hits of one file into h (smeared, from threshold up); returns the entries filled
  auto Fill = [&](TTree* ev, TH1D* h) {
    std::vector<double>* e = nullptr;
    ev->SetBranchStatus("*", 0);
    ev->SetBranchStatus("czt_e", 1);
    ev->SetBranchAddress("czt_e", &e);
    double n = 0.;
    for (Long64_t i = 0; i < ev->GetEntries(); ++i) {
      ev->GetEntry(i);
      for (double eMeV : *e) {
        if (eMeV <= 0.) continue;
        const double kev = 1000. * eMeV + random.Gaus(0., sigma1MeV * std::sqrt(eMeV));
        if (kev >= threshold) { h->Fill(kev); n += 1.; }
      }
    }
    return n;
  };
  auto RunSum = [](TTree* run, const char* column) {
    double s = 0.;
    run->Draw(column, "", "goff");
    for (Long64_t i = 0; i < run->GetSelectedRows(); ++i) s += run->GetV1()[i];
    return s;
  };

  // the floor: one group per nuclide and region
  struct Group {
    std::string nuclide, tag, region;
    int jobs = 0;
    double events = 0., grams = 0., bqPerGram = 0., hits = 0.;
    TH1D* h = nullptr;
  };
  const std::vector<std::string> nuclides = {"U-238", "Th-232", "K-40", "U-235"}, regions = {"inner", "px", "nx", "py", "ny"};
  std::map<std::string, Group> groups;   // "<nuclide> <region>"
  const std::string floorDir = base + "/floor";
  for (int k : FinishedJobs(floorDir, "floor")) {
    const std::string mac = Form("%s/log/job%d.mac", floorDir.c_str(), k);
    const int z = std::atoi(MacroValue(mac, "/control/alias bmIonZ").c_str()), a = std::atoi(MacroValue(mac, "/control/alias bmIonA").c_str());
    const std::string region = MacroValue(mac, "/control/alias rmRegion");
    const double hx = std::atof(MacroValue(mac, "/control/alias rmHalfX").c_str()), hy = std::atof(MacroValue(mac, "/control/alias rmHalfY").c_str()),
                 hz = std::atof(MacroValue(mac, "/control/alias rmHalfZ").c_str());
    std::string nuclide, tag;
    double bq = 0.;
    if (z == 92 && a == 238) { nuclide = "U-238"; tag = "U238"; bq = uBqKg / 1000.; }
    else if (z == 90 && a == 232) { nuclide = "Th-232"; tag = "Th232"; bq = thBqKg / 1000.; }
    else if (z == 19 && a == 40) { nuclide = "K-40"; tag = "K40"; bq = kBqKg / 1000.; }
    else if (z == 92 && a == 235) { nuclide = "U-235"; tag = "U235"; bq = uBqKg * kU235PerU238 / 1000.; }
    else { printf("floor_%d.root: a source Z %d A %d this analysis does not know; skipped\n", k, z, a); continue; }
    if (std::find(regions.begin(), regions.end(), region) == regions.end() || hx <= 0. || hy <= 0. || hz <= 0.) {
      printf("floor job %d: region \"%s\" (%g x %g x %g) unknown; skipped\n", k, region.c_str(), hx, hy, hz);
      continue;
    }
    TFile f(Form("%s/floor_%d.root", floorDir.c_str(), k));
    auto run = dynamic_cast<TTree*>(f.Get("run"));
    auto ev = dynamic_cast<TTree*>(f.Get("events"));
    if (!run || !ev) { printf("skipped floor_%d.root: no run or events tree\n", k); continue; }
    Group& g = groups[nuclide + " " + region];
    const double grams = 8. * hx * hy * hz * kConcrete;
    if (!g.h) {
      g.nuclide = nuclide; g.tag = tag; g.region = region; g.grams = grams; g.bqPerGram = bq;
      g.h = NewHist("room_floor_" + tag + "_" + region);
    } else if (std::fabs(g.grams - grams) > 1e-6 * grams) {
      printf("floor job %d: region %s with another size than its other jobs; skipped\n", k, region.c_str());
      continue;
    }
    g.jobs++;
    g.events += RunSum(run, "n_events");
    g.hits += Fill(ev, g.h);
  }
  std::map<std::string, TH1D*> perNuclide;
  for (size_t i = 0; i < nuclides.size(); ++i) {
    const char* tags[] = {"U238", "Th232", "K40", "U235"};
    perNuclide[nuclides[i]] = NewHist(std::string("room_floor_") + tags[i]);
  }
  TH1D* floorAll = NewHist("room_floor");
  for (auto& [key, g] : groups) {
    if (g.events > 0.) g.h->Scale(g.grams * g.bqPerGram * seconds / g.events);
    perNuclide[g.nuclide]->Add(g.h);
    floorAll->Add(g.h);
  }

  // radon in the castle's air
  TH1D* radon = NewHist("room_radon");
  int radonJobs = 0;
  double radonEvents = 0., radonHits = 0.;
  for (int k : FinishedJobs(base + "/radon", "radon")) {
    const std::string mac = Form("%s/radon/log/job%d.mac", base.c_str(), k);
    if (MacroValue(mac, "/control/alias bmIonZ") != "86" || MacroValue(mac, "/control/alias bmIonA") != "222") {
      printf("radon job %d is not Rn-222; skipped\n", k);
      continue;
    }
    TFile f(Form("%s/radon/radon_%d.root", base.c_str(), k));
    auto run = dynamic_cast<TTree*>(f.Get("run"));
    auto ev = dynamic_cast<TTree*>(f.Get("events"));
    if (!run || !ev) { printf("skipped radon_%d.root: no run or events tree\n", k); continue; }
    ++radonJobs;
    radonEvents += RunSum(run, "n_events");
    radonHits += Fill(ev, radon);
  }
  const double radonBq = radonBqM3 * kCastleAirM3;
  if (radonEvents > 0.) radon->Scale(radonBq * seconds / radonEvents);

  // cosmic rays
  TH1D* cosmic = NewHist("room_cosmic");
  int cosmicJobs = 0;
  double cosmicEvents = 0., cosmicLive = 0., cosmicHits = 0.;
  for (int k : FinishedJobs(base + "/cosmic", "cosmic")) {
    TFile f(Form("%s/cosmic/cosmic_%d.root", base.c_str(), k));
    auto run = dynamic_cast<TTree*>(f.Get("run"));
    auto ev = dynamic_cast<TTree*>(f.Get("events"));
    if (!run || !ev) { printf("skipped cosmic_%d.root: no run or events tree\n", k); continue; }
    ++cosmicJobs;
    cosmicEvents += RunSum(run, "n_events");
    cosmicLive += RunSum(run, "live_time");
    cosmicHits += Fill(ev, cosmic);
  }
  if (cosmicLive > 0.) cosmic->Scale(seconds / cosmicLive);

  // the walls and the ceiling through the lead (estimate): the strongest lines
  TH1D* walls = NewHist("room_walls");
  struct WallRun { WallLine l; double events = 0., hits = 0., flux = 0., perDay = 0., entries = 0.; };
  std::vector<WallRun> wallRuns;
  std::vector<WallLine> wallMissing;
  const bool withWalls = !std::string(linesDir).empty();
  if (withWalls)
    for (const WallLine& l : kWallLines) {
      const std::string file = Form("%s/lines_%.2f.root", linesDir, l.keV);
      if (gSystem->AccessPathName(file.c_str())) { wallMissing.push_back(l); continue; }
      TFile f(file.c_str());
      auto run = dynamic_cast<TTree*>(f.Get("run"));
      auto ev = dynamic_cast<TTree*>(f.Get("events"));
      if (!run || !ev) { printf("skipped %s: no run or events tree\n", file.c_str()); wallMissing.push_back(l); continue; }
      WallRun w;
      w.l = l;
      w.events = RunSum(run, "n_events");
      TH1D* h = NewHist(Form("walls_%.0f", l.keV));
      w.hits = Fill(ev, h);
      const double a = (l.chain == 0 ? uBqKg : l.chain == 1 ? thBqKg : kBqKg) / 1000.;
      w.flux = a * l.yield / MuRhoConcrete(l.keV);
      w.perDay = w.flux / 4. * kSphereArea * seconds;
      if (w.events > 0.) h->Scale(w.perDay / w.events);
      w.entries = h->Integral();
      walls->Add(h);
      delete h;
      wallRuns.push_back(w);
    }

  TH1D* total = NewHist("room_total");
  for (TH1D* h : {floorAll, radon, cosmic, walls}) total->Add(h);

  // the ROOT file: every histogram per day in 0.5 keV bins
  {
    TFile hf((out + "/room_background_spectrum.root").c_str(), "RECREATE");
    for (auto& [key, g] : groups) g.h->Write();
    for (auto& [n, h] : perNuclide) h->Write();
    for (TH1D* h : {floorAll, radon, cosmic, walls, total}) h->Write();
    TNamed("settings", Form("sigma1MeV %g threshold %g seconds %g K40 %g U238 %g Th232 %g Bq/kg radon %g Bq/m3 castle air %g m3 walls %s",
                            sigma1MeV, threshold, seconds, kBqKg, uBqKg, thBqKg, radonBqM3, kCastleAirM3,
                            withWalls ? linesDir : "none")).Write();
  }

  // numbers
  auto Entries = [](const TH1D* h) { return WindowSum(h, 0., kMaxKeV); };
  txt << "Room background of the steel study on the CZT 10 x 10 cm panel (sigma " << sigma1MeV << " keV at 1 MeV, from " << threshold
      << " keV), per day of " << seconds << " s: " << base << "\n"
      << "  setup: the steel beam run's (plate on the panel under a 10 cm lead roof, between four 10 cm lead walls, open at the bottom)"
         " on a 30 cm concrete floor\n";
  snprintf(line, sizeof line, "  concrete (G4_CONCRETE, %.1f g/cm3): K-40 %g, U-238 %g (chain in equilibrium), Th-232 %g Bq/kg, U-235 %.3g Bq/kg; "
           "radon %g Bq/m3 in the castle's %.4g m3 of air (%.3g Bq)\n", kConcrete, kBqKg, uBqKg, thBqKg, uBqKg * kU235PerU238, radonBqM3,
           kCastleAirM3, radonBq);
  txt << line;
  csv << "part,nuclide,region,jobs,events,activity_Bq,simulated_days,simulated_entries,entries_per_day,entries_per_s\n";
  txt << "\n  floor: nuclide region | jobs, simulated decays, activity (Bq), simulated days | panel entries: simulated, per day, per s\n";
  for (const std::string& n : nuclides)
    for (const std::string& r : regions) {
      auto it = groups.find(n + " " + r);
      if (it == groups.end()) { txt << Form("    %-6s %-5s | no finished job\n", n.c_str(), r.c_str()); continue; }
      const Group& g = it->second;
      const double act = g.grams * g.bqPerGram, simDays = act > 0. ? g.events / act / 86400. : 0., perDay = Entries(g.h);
      snprintf(line, sizeof line, "    %-6s %-5s | %3d %10.4g %10.4g %9.4g | %9.0f %10.4g %8.4g\n", n.c_str(), r.c_str(), g.jobs, g.events, act,
               simDays, g.hits, perDay, perDay / seconds);
      txt << line;
      csv << Form("floor,%s,%s,%d,%.6g,%.6g,%.6g,%.0f,%.6g,%.6g\n", n.c_str(), r.c_str(), g.jobs, g.events, act, simDays, g.hits, perDay,
                  perDay / seconds);
    }
  for (const std::string& n : nuclides) {
    double inner = 0., all = 0.;
    for (const std::string& r : regions) {
      auto it = groups.find(n + " " + r);
      if (it == groups.end()) continue;
      const double e = Entries(it->second.h);
      all += e;
      if (r == "inner") inner += e;
    }
    snprintf(line, sizeof line, "    %-6s all regions: %.4g entries per day (%.4g /s), %.1f%% from the inner region\n", n.c_str(), all, all / seconds,
             all > 0. ? 100. * inner / all : 0.);
    txt << line;
  }
  const double floorDay = Entries(floorAll), radonDay = Entries(radon), cosmicDay = Entries(cosmic), wallsDay = Entries(walls),
               totalDay = Entries(total);
  snprintf(line, sizeof line, "    the floor: %.4g entries per day (%.4g /s)\n", floorDay, floorDay / seconds);
  txt << line;
  snprintf(line, sizeof line, "\n  radon: %d jobs, %.4g simulated decays of Rn-222 (%.3g days) | %.0f simulated entries, %.4g per day (%.4g /s)\n",
           radonJobs, radonEvents, radonBq > 0. ? radonEvents / radonBq / 86400. : 0., radonHits, radonDay, radonDay / seconds);
  txt << line;
  csv << Form("radon,Rn-222,castle,%d,%.6g,%.6g,%.6g,%.0f,%.6g,%.6g\n", radonJobs, radonEvents, radonBq,
              radonBq > 0. ? radonEvents / radonBq / 86400. : 0., radonHits, radonDay, radonDay / seconds);
  snprintf(line, sizeof line, "  cosmic rays: %d jobs, %.4g simulated particles, %.4g s of live time (%.3g days) | %.0f simulated entries, %.4g per day (%.4g /s)\n",
           cosmicJobs, cosmicEvents, cosmicLive, cosmicLive / 86400., cosmicHits, cosmicDay, cosmicDay / seconds);
  txt << line;
  csv << Form("cosmic,all,sky,%d,%.6g,,%.6g,%.0f,%.6g,%.6g\n", cosmicJobs, cosmicEvents, cosmicLive / 86400., cosmicHits, cosmicDay, cosmicDay / seconds);
  if (withWalls) {
    txt << "\n  walls and ceiling through the lead (estimate: the uncollided lines of a closed concrete room; sphere of " << Form("%.0f", kSphereArea)
        << " cm2):\n    line (keV) nuclide yield | flux (/cm2/s) | photons per day on the sphere, simulated | entries: simulated, per day\n";
    for (const WallRun& w : wallRuns) {
      snprintf(line, sizeof line, "    %8.2f %-7s %6.4f | %8.4g | %10.4g %10.4g | %6.0f %9.4g\n", w.l.keV, w.l.nuclide, w.l.yield, w.flux, w.perDay,
               w.events, w.hits, w.entries);
      txt << line;
      csv << Form("walls,%s,%.2f,1,%.6g,%.6g,,%.0f,%.6g,%.6g\n", w.l.nuclide, w.l.keV, w.events, w.perDay, w.hits, w.entries, w.entries / seconds);
    }
    for (const WallLine& l : wallMissing) txt << Form("    %8.2f %-7s: no run (lines_%.2f.root) - left out\n", l.keV, l.nuclide, l.keV);
    snprintf(line, sizeof line, "    walls and ceiling: %.4g entries per day (%.4g /s), %.2f%% of the floor's\n", wallsDay, wallsDay / seconds,
             floorDay > 0. ? 100. * wallsDay / floorDay : 0.);
    txt << line;
  } else txt << "\n  walls and ceiling: not estimated (no linesDir)\n";
  snprintf(line, sizeof line, "\n  the room in all: %.4g entries per day (%.4g /s): floor %.1f%%, radon %.2f%%, cosmic rays %.2f%%, walls and ceiling %.2f%%\n",
           totalDay, totalDay / seconds, totalDay > 0. ? 100. * floorDay / totalDay : 0., totalDay > 0. ? 100. * radonDay / totalDay : 0.,
           totalDay > 0. ? 100. * cosmicDay / totalDay : 0., totalDay > 0. ? 100. * wallsDay / totalDay : 0.);
  txt << line;
  csv << Form("total,all,all,,,,,,%.6g,%.6g\n", totalDay, totalDay / seconds);
  // the room's natural lines (+-1.5 sigma): what a beam-off screening of the plate has under its own lines
  txt << "\n  the room's entries per day at the natural lines (+-1.5 sigma): line | floor, radon, cosmic, walls, total\n";
  const std::vector<std::pair<double, const char*>> natLines = {{351.93, "Pb-214"}, {609.31, "Bi-214"}, {1460.82, "K-40"}, {2614.51, "Tl-208"}};
  for (const auto& [e0, nuc] : natLines) {
    const double sg = sigma1MeV * std::sqrt(e0 / 1000.), lo = e0 - 1.5 * sg, hi = e0 + 1.5 * sg;
    snprintf(line, sizeof line, "    %8.2f %-7s | %9.4g %9.4g %9.4g %9.4g %9.4g\n", e0, nuc, WindowSum(floorAll, lo, hi), WindowSum(radon, lo, hi),
             WindowSum(cosmic, lo, hi), WindowSum(walls, lo, hi), WindowSum(total, lo, hi));
    txt << line;
    csv << Form("natural_line,%s,%.2f,,,,,,%.6g,\n", nuc, e0, WindowSum(total, lo, hi));
  }

  // the plot: the room's parts per day
  {
    const int rebin = 4;   // 2 keV
    auto c = new TCanvas("c_room", "", 1100, 750);
    c->SetLogy(); c->SetLogx();
    c->SetLeftMargin(0.11); c->SetRightMargin(0.03); c->SetTopMargin(0.08); c->SetBottomMargin(0.12);
    TH1D* t = Copy(total, "plot_total", 1., rebin);
    t->GetXaxis()->SetTitleOffset(1.2);
    t->SetTitle(Form("The room's background in the CZT 10 x 10 cm panel, one day: the steel run's lead castle on a 30 cm concrete floor;"
                     "energy in one pixel (keV);entries per %d keV in one day", int(rebin * kFillKeV)));
    t->GetXaxis()->SetRangeUser(20., 7000.);
    t->GetXaxis()->SetMoreLogLabels();
    t->SetLineColor(kBlack); t->SetLineWidth(2);
    const double top = t->GetMaximum();
    t->SetMaximum(top * 300.); t->SetMinimum(0.3);
    t->Draw("hist");
    struct Part { TH1D* h; std::string label; int color; int style; };   // labels as strings: Form's buffer is reused
    std::vector<Part> parts = {
        {perNuclide["U-238"], Form("floor, U-238 chain (%g Bq/kg)", uBqKg), kViolet, 1},
        {perNuclide["Th-232"], Form("floor, Th-232 chain (%g Bq/kg)", thBqKg), kAzure + 1, 1},
        {perNuclide["K-40"], Form("floor, K-40 (%g Bq/kg)", kBqKg), kOrange + 1, 1},
        {perNuclide["U-235"], Form("floor, U-235 chain (%.2g Bq/kg)", uBqKg * kU235PerU238), kMagenta, 1},
        {radon, Form("radon in the castle's air (%g Bq/m^{3})", radonBqM3), kGreen + 2, 1},
        {cosmic, "cosmic rays (sea level, open sky)", kGray + 2, 1}};
    if (withWalls) parts.push_back({walls, "walls and ceiling through the lead (estimate)", kRed, 2});
    auto leg = new TLegend(0.13, 0.73, 0.97, 0.91);
    leg->SetNColumns(2); leg->SetMargin(0.12); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.022);
    leg->AddEntry(t, Form("the room in all (%s entries, %.3g /s)", Sci(totalDay).c_str(), totalDay / seconds), "l");
    for (const Part& p : parts) {
      TH1D* h = Copy(p.h, std::string("plot_") + p.h->GetName(), 1., rebin);
      h->GetXaxis()->SetRangeUser(20., 7000.);
      h->SetLineColor(p.color); h->SetLineWidth(2); h->SetLineStyle(p.style);
      h->Draw("hist same");
      leg->AddEntry(h, (p.label + ": " + Sci(Entries(p.h))).c_str(), "l");
    }
    t->Draw("hist same");
    leg->Draw();
    Save(c, out + "/room_background_spectrum");
  }

  std::ofstream(out + "/room_background_spectrum.txt") << txt.str();
  std::ofstream(out + "/room_background_spectrum.csv") << csv.str();
  printf("%s", txt.str().c_str());
  printf("written into %s: room_background_spectrum.root, .txt, .csv, room_background_spectrum (.png, .pdf)\n", out.c_str());
}
