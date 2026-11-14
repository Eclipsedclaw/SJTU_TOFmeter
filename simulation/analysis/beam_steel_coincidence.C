// The cascade coincidences of the steel study (2026-11-12, the user: study 3 of the suggestions for an element
// spectrometer, "Also do 3"; all 13 elements, each with its best pair). One muon's cascade sends its K, L, M X-rays out
// together: two hits of one muon in two pixels of the CZT panel, both within promptNs (10 ns) of the muon, per day at the
// run's rate, as a map of the two energies. For each element its best pair of lines: its own X-ray pairs in the box of
// +-1.5 resolutions around both lines (S) against all other pairs (B), compared with its best single line in the
// tagged-beam spectrum (the windows and terms of beam_steel_trace_spectrum.C: S the element's own X-ray entries, B the rest).
//
// Input: the beam band (beam_<k>.root: cap_z, czt_pix, czt_e, czt_t, czt_src), read as beam_steel_trace_spectrum.C reads
// it: each element's muons reweighted to the real steel (Z law), each hit smeared by sigma1MeV x sqrt(E / MeV) with the
// same random numbers in the same order (the singles here are the trace's prompt view: checked), threshold 20 keV.
// From analysis/beam_steel_trace_spectrum.root: the emitted lines (truth) and the singles (prompt, beam on).
// Accidental pairs (a hit of another muon in the window): prompt singles x beam-on singles per muon x rate x window.
// The panel: 4 x 4 pixels of 25 x 25 x 10 mm; two photons in one pixel make one hit (their sum), not a pair.
// Outputs in outDir (default <dir>/analysis): beam_steel_coincidence.txt, .csv, .root (the maps per day), and the plots
//   coincidence_map        all pairs of one day, the higher energy against the lower one (log axes)
//   coincidence_k_spot     around K K_alpha x K L_alpha: all pairs and K's own X-ray pairs (two pads), the box
//   coincidence_k_gate     the K K_alpha region: the singles, and the pairs gated on K L_alpha (two pads)
//   coincidence_s_over_b   per element: S/B of its best single line and of its best pair
//   coincidence_limits     per element: the mass fraction seen at 3 sigma in one day, single line and pair, without and
//                          with a 1% uncertainty on the background
// Run in the work folder after beam_steel_trace_spectrum.C (reads the 42 GB, about 12 min, one thread):
//   root -l -b -q 'analysis/beam_steel_coincidence.C("output/czt/beam_steel_czt_trace")'
#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace coincidence {
struct Element {   // the real steel, mass percent (as beam_steel_trace_spectrum.C)
  int z;
  const char* name;
  double a;
  double percent;
};

std::vector<int> FinishedJobs(const std::string& bandDir, const char* prefix)   // as beam_steel_trace_spectrum.C
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

const int kSources = 14;   // the 13 elements of the steel, then the muons caught elsewhere

struct Line {   // an emitted line (truth; within half a resolution merged); its own singles and all singles in +-1.5 sigma
  double e = 0., perDay = 0., own = 0., total = 0.;
};
struct Box {   // two lines of one element: the windows; the simulated pairs inside, per source and kind
  int el = 0;
  double ea = 0., eb = 0.;   // the higher line, the lower one
  double hiLo = 0., hiHi = 0., loLo = 0., loHi = 0.;
  double n[kSources][2] = {};   // [source][both hits X-rays of the plate]
  double S = 0., B = 0., acc = 0., sSim = 0., bSim = 0.;
};
// the mass fraction whose signal S' gives S' / sqrt(S' + B + (u B)^2) = 3 (S at the fraction f, B stays)
double Fraction3(double f, double S, double B, double u)
{
  const double b = B + u * u * B * B;
  return S > 0. ? f * 0.5 * (9. + std::sqrt(81. + 36. * b)) / S : 0.;
}
double WindowSum(const TH1D* h, double lo, double hi)   // the bins with centres in [lo, hi] (as the trace)
{
  double s = 0.;
  for (int b = std::max(1, h->GetXaxis()->FindFixBin(lo) - 1); b <= std::min(h->GetNbinsX(), h->GetXaxis()->FindFixBin(hi) + 1); ++b) {
    const double x = h->GetBinCenter(b);
    if (x >= lo && x <= hi) s += h->GetBinContent(b);
  }
  return s;
}
}  // namespace coincidence

void beam_steel_coincidence(const char* dir, const char* outDir = "", double fracU = 1e-7, double fracTh = 1e-7,
                            double fracK = 5e-5, double rate = 1.32e5, double seconds = 86400., double promptNs = 10.,
                            double sigma1MeV = 6.36991, double threshold = 20.)
{
  using namespace coincidence;
  gStyle->SetOptStat(0);
  const std::string base = dir;
  const std::string out = std::string(outDir).empty() ? base + "/analysis" : std::string(outDir);
  gSystem->mkdir(out.c_str(), true);
  auto Sigma = [&](double kev) { return sigma1MeV * std::sqrt(std::max(kev, 1.) / 1000.); };

  // the real steel and each element's share of the muons (Z law), as beam_steel_trace_spectrum.C
  std::vector<Element> steel = {{26, "Fe", 55.845, 0.}, {24, "Cr", 51.996, 16.}, {28, "Ni", 58.693, 10.},
                                {42, "Mo", 95.95, 2.}, {25, "Mn", 54.938, 2.}, {14, "Si", 28.086, 0.75},
                                {7, "N", 14.007, 0.10}, {6, "C", 12.011, 0.08}, {15, "P", 30.974, 0.045},
                                {16, "S", 32.06, 0.03}, {19, "K", 39.098, 100. * fracK}, {90, "Th", 232.038, 100. * fracTh},
                                {92, "U", 238.029, 100. * fracU}};
  double rest = 100.;
  for (size_t i = 1; i < steel.size(); ++i) rest -= steel[i].percent;
  steel[0].percent = rest;
  double zSum = 0.;
  for (const Element& el : steel) zSum += el.percent / el.a * el.z;
  const int nEl = steel.size();
  std::map<int, int> index;   // Z -> source
  for (int i = 0; i < nEl; ++i) index[steel[i].z] = i;
  auto SourceName = [&](int s) { return std::string(s < nEl ? steel[s].name : "other"); };

  // the trace's singles (per day, 0.5 keV) and the emitted lines
  const std::string traceName = base + "/analysis/beam_steel_trace_spectrum.root";
  std::unique_ptr<TFile> trace(TFile::Open(traceName.c_str()));
  if (!trace || trace->IsZombie()) { printf("cannot open %s: run analysis/beam_steel_trace_spectrum.C first\n", traceName.c_str()); return; }
  auto Get = [&](const std::string& name) -> TH1D* {
    auto h = dynamic_cast<TH1D*>(trace->Get(name.c_str()));
    if (!h) { printf("no %s in %s\n", name.c_str(), traceName.c_str()); return nullptr; }
    h = dynamic_cast<TH1D*>(h->Clone((name + "_c").c_str()));
    h->SetDirectory(nullptr);
    return h;
  };
  TH1D* viewPrompt = Get("view_prompt");
  TH1D* viewOn = Get("view_on");
  TH1D* accPrompt = Get("accidental_prompt");
  if (!viewPrompt || !viewOn || !accPrompt) return;
  auto beamPrompt = dynamic_cast<TH1D*>(viewPrompt->Clone("beam_prompt"));   // the prompt view without the accidentals
  beamPrompt->SetDirectory(nullptr);
  beamPrompt->Add(accPrompt, -1.);
  std::vector<TH1D*> ownX(nEl);
  for (int i = 0; i < nEl; ++i)
    if (!(ownX[i] = Get(std::string("day_") + steel[i].name + "_promptx"))) return;
  if (auto s = dynamic_cast<TNamed*>(trace->Get("settings"))) {
    double s1 = 0., th = 0., pn = 0., ra = 0., se = 0.;
    if (sscanf(s->GetTitle(), "sigma1MeV %lf threshold %lf promptNs %lf rate %lf seconds %lf", &s1, &th, &pn, &ra, &se) == 5 &&
        (std::fabs(s1 / sigma1MeV - 1.) > 1e-6 || std::fabs(th / threshold - 1.) > 1e-6 || std::fabs(pn / promptNs - 1.) > 1e-6 ||
         std::fabs(ra / rate - 1.) > 1e-6 || std::fabs(se / seconds - 1.) > 1e-6))
      printf("warning: the trace ran with other settings (%s): its singles and these differ\n", s->GetTitle());
  }
  std::vector<std::vector<Line>> lines(nEl);
  {
    auto t = dynamic_cast<TTree*>(trace->Get("lines"));
    if (!t) { printf("no lines tree in %s\n", traceName.c_str()); return; }
    int z = 0;
    double e = 0., pd = 0.;
    t->SetBranchAddress("z", &z); t->SetBranchAddress("e_keV", &e); t->SetBranchAddress("per_day", &pd);
    std::vector<std::vector<std::pair<double, double>>> raw(nEl);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      const auto it = index.find(z);
      if (it != index.end()) raw[it->second].push_back({e, pd});
    }
    for (int i = 0; i < nEl; ++i) {
      std::sort(raw[i].begin(), raw[i].end());
      std::vector<Line> merged;
      for (const auto& [e0, p] : raw[i])
        if (!merged.empty() && e0 - merged.back().e < 0.5 * Sigma(e0)) {
          Line& m = merged.back();
          m.e = (m.e * m.perDay + e0 * p) / (m.perDay + p);
          m.perDay += p;
        } else {
          Line l;
          l.e = e0; l.perDay = p;
          merged.push_back(l);
        }
      // each line's own singles; the eight largest of each element are the candidates for its pairs
      std::vector<Line> kept;
      for (Line& l : merged) {
        const double sg = Sigma(l.e);
        if (l.e - 1.5 * sg < threshold) continue;
        l.own = WindowSum(ownX[i], l.e - 1.5 * sg, l.e + 1.5 * sg);
        l.total = WindowSum(viewPrompt, l.e - 1.5 * sg, l.e + 1.5 * sg);
        if (l.own > 0.) kept.push_back(l);
      }
      std::sort(kept.begin(), kept.end(), [](const Line& a, const Line& b) { return a.own > b.own; });
      if (kept.size() > 8) kept.resize(8);
      std::sort(kept.begin(), kept.end(), [](const Line& a, const Line& b) { return a.e > b.e; });
      lines[i] = kept;
    }
  }
  // the boxes: every two candidate lines of one element whose windows do not overlap (else the two hits could swap)
  std::vector<Box> boxes;
  for (int i = 0; i < nEl; ++i)
    for (size_t a = 0; a < lines[i].size(); ++a)
      for (size_t b = a + 1; b < lines[i].size(); ++b) {
        const double ea = lines[i][a].e, eb = lines[i][b].e, sa = Sigma(ea), sb = Sigma(eb);
        if (ea - 1.5 * sa <= eb + 1.5 * sb) continue;
        Box x;
        x.el = i; x.ea = ea; x.eb = eb;
        x.hiLo = ea - 1.5 * sa; x.hiHi = ea + 1.5 * sa; x.loLo = eb - 1.5 * sb; x.loHi = eb + 1.5 * sb;
        boxes.push_back(x);
      }
  const double kStep = 2.;   // the boxes by the higher energy in 2 keV steps, to find a pair's boxes quickly
  std::vector<std::vector<int>> byHi(int(8000. / kStep) + 2);
  for (size_t k = 0; k < boxes.size(); ++k)
    for (int st = int(boxes[k].hiLo / kStep); st <= int(boxes[k].hiHi / kStep) && st < int(byHi.size()); ++st) byHi[st].push_back(k);
  // K: K_alpha and L_alpha among its candidates
  const int iK = index[19];
  double eKa = 0., eLa = 0.;
  for (const Line& l : lines[iK]) {
    if (std::fabs(l.e - 712.7) < 5. && (eKa == 0. || std::fabs(l.e - 712.7) < std::fabs(eKa - 712.7))) eKa = l.e;
    if (std::fabs(l.e - 141.8) < 3. && (eLa == 0. || std::fabs(l.e - 141.8) < std::fabs(eLa - 141.8))) eLa = l.e;
  }
  if (eKa == 0. || eLa == 0.) { printf("K K_alpha (712.7 keV) or L_alpha (141.8 keV) is not among K's lines: stopping\n"); return; }
  const double gLo = eLa - 1.5 * Sigma(eLa), gHi = eLa + 1.5 * Sigma(eLa), kLo = eKa - 1.5 * Sigma(eKa), kHi = eKa + 1.5 * Sigma(eKa);

  // the histograms (simulated pairs, unweighted; per day below): the map in cells of half a resolution each way
  std::vector<double> edges = {threshold};
  while (edges.back() < 3000.) edges.push_back(edges.back() + std::max(0.5, 0.5 * std::round(Sigma(edges.back()))));
  const int nE = edges.size() - 1;
  std::vector<TH2F*> pairMap(kSources);
  std::vector<std::array<TH2F*, 2>> spot(kSources);
  std::vector<std::array<TH1D*, 2>> gate(kSources);
  std::vector<TH1D*> singles(kSources);
  for (int s = 0; s < kSources; ++s) {
    pairMap[s] = new TH2F(("map_" + SourceName(s)).c_str(), "", nE, edges.data(), nE, edges.data());
    pairMap[s]->SetDirectory(nullptr);
    for (int x = 0; x < 2; ++x) {
      spot[s][x] = new TH2F(Form("spot_%s_%d", SourceName(s).c_str(), x), "", 35, 125., 160., 28, 680., 750.);
      gate[s][x] = new TH1D(Form("gate_%s_%d", SourceName(s).c_str(), x), "", 75, 640., 790.);
      spot[s][x]->SetDirectory(nullptr);
      gate[s][x]->SetDirectory(nullptr);
    }
    singles[s] = new TH1D(("singles_" + SourceName(s)).c_str(), "", viewPrompt->GetNbinsX(), viewPrompt->GetXaxis()->GetXmin(),
                          viewPrompt->GetXaxis()->GetXmax());
    singles[s]->SetDirectory(nullptr);
  }

  // the beam band, event by event
  const std::string beamDir = base + "/beam";
  const std::vector<int> beamJobs = FinishedJobs(beamDir, "beam");
  if (beamJobs.empty()) { printf("no finished beam jobs (beam_<k>.root with log/job<k>.done) in %s\n", beamDir.c_str()); return; }
  std::vector<double> events(kSources, 0.);
  double nSim = 0.;
  long nPrompt = 0, nPairs = 0, nSamePixel = 0;
  TRandom3 random(4357);   // the trace's seed, the trace's order: the same smeared energies
  const double tMax = seconds * 1e9;
  struct Hit { double kev; int pix; bool x; };
  std::vector<Hit> hits;
  int done = 0;
  for (int k : beamJobs) {
    TFile f(Form("%s/beam_%d.root", beamDir.c_str(), k));
    auto run = dynamic_cast<TTree*>(f.Get("run"));
    auto ev = dynamic_cast<TTree*>(f.Get("events"));
    if (!run || !ev) { printf("skipped beam_%d.root: no run or events tree\n", k); continue; }
    run->Draw("n_events:event_filter", "", "goff");
    for (Long64_t i = 0; i < run->GetSelectedRows(); ++i) {
      nSim += run->GetV1()[i];
      if (run->GetV2()[i] != 0) { printf("beam_%d.root was run with /EventFilter other than 0: stopping\n", k); return; }
    }
    int capZ = 0;
    std::vector<double>*e = nullptr, *t = nullptr;
    std::vector<int>*src = nullptr, *pix = nullptr;
    ev->SetBranchStatus("*", 0);
    for (const char* b : {"cap_z", "czt_e", "czt_t", "czt_src", "czt_pix"}) ev->SetBranchStatus(b, 1);
    ev->SetBranchAddress("cap_z", &capZ);
    ev->SetBranchAddress("czt_e", &e);
    ev->SetBranchAddress("czt_t", &t);
    ev->SetBranchAddress("czt_src", &src);
    ev->SetBranchAddress("czt_pix", &pix);
    for (Long64_t i = 0; i < ev->GetEntries(); ++i) {
      ev->GetEntry(i);
      const auto it = index.find(capZ);
      const int s = it != index.end() ? it->second : nEl;
      events[s] += 1.;
      hits.clear();
      for (size_t h = 0; h < e->size(); ++h) {
        const double eMeV = (*e)[h], tNs = std::max(0., (*t)[h]);
        if (eMeV <= 0. || tNs >= 2. * tMax) continue;
        const double kev = 1000. * eMeV + random.Gaus(0., sigma1MeV * std::sqrt(eMeV));   // as the trace, in its order
        if (kev < threshold || tNs >= promptNs) continue;
        hits.push_back({kev, (*pix)[h], ((*src)[h] & 1) != 0});
        singles[s]->Fill(kev);
      }
      nPrompt += hits.size();
      for (size_t a = 0; a < hits.size(); ++a)
        for (size_t b = a + 1; b < hits.size(); ++b) {
          if (hits[a].pix == hits[b].pix) { ++nSamePixel; continue; }
          const Hit& hi = hits[a].kev >= hits[b].kev ? hits[a] : hits[b];
          const Hit& lo = hits[a].kev >= hits[b].kev ? hits[b] : hits[a];
          const int x = hi.x && lo.x;
          ++nPairs;
          pairMap[s]->Fill(lo.kev, hi.kev);
          if (lo.kev >= 125. && lo.kev < 160. && hi.kev >= 680. && hi.kev < 750.) spot[s][x]->Fill(lo.kev, hi.kev);
          if (lo.kev >= gLo && lo.kev <= gHi) gate[s][x]->Fill(hi.kev);
          const int st = int(hi.kev / kStep);
          if (st < int(byHi.size()))
            for (int bx : byHi[st]) {
              Box& bo = boxes[bx];
              if (hi.kev >= bo.hiLo && hi.kev <= bo.hiHi && lo.kev >= bo.loLo && lo.kev <= bo.loHi) bo.n[s][x] += 1.;
            }
        }
    }
    if (++done % 20 == 0) { printf("  %d of %zu beam files\n", done, beamJobs.size()); fflush(stdout); }
  }

  // per day: each element's muons to the real steel, as the trace
  const double nDay = rate * seconds;
  std::vector<double> w(kSources, 0.);
  for (int i = 0; i < nEl; ++i) {
    if (events[i] <= 0.) { printf("no muon caught by %s (Z %d): stopping\n", steel[i].name, steel[i].z); return; }
    w[i] = nDay * steel[i].percent / steel[i].a * steel[i].z / zSum / events[i];
  }
  w[nEl] = nDay / nSim;
  // the check: these prompt singles, weighted, against the trace's prompt view without the accidentals
  auto mine = dynamic_cast<TH1D*>(singles[0]->Clone("singles_all"));
  mine->SetDirectory(nullptr);
  mine->Reset();
  for (int s = 0; s < kSources; ++s) mine->Add(singles[s], w[s]);
  double maxDev = 0.;
  for (int b = 1; b <= mine->GetNbinsX(); ++b)
    if (beamPrompt->GetBinContent(b) > 0.) maxDev = std::max(maxDev, std::fabs(mine->GetBinContent(b) / beamPrompt->GetBinContent(b) - 1.));
  const double accScale = rate * promptNs * 1e-9 / nDay;   // a foreign hit per muon's window, per beam-on entry per day
  for (Box& bo : boxes) {
    double all = 0., allSim = 0.;
    for (int s = 0; s < kSources; ++s)
      for (int x = 0; x < 2; ++x) { all += bo.n[s][x] * w[s]; allSim += bo.n[s][x]; }
    bo.sSim = bo.n[bo.el][1];
    bo.S = bo.sSim * w[bo.el];
    bo.acc = accScale * (WindowSum(beamPrompt, bo.hiLo, bo.hiHi) * WindowSum(viewOn, bo.loLo, bo.loHi) +
                         WindowSum(beamPrompt, bo.loLo, bo.loHi) * WindowSum(viewOn, bo.hiLo, bo.hiHi));
    bo.B = all - bo.S + bo.acc;
    bo.bSim = allSim - bo.sSim;
  }
  // each element's best single line and best pair (by S / sqrt(S + B) in one day; a pair needs 10 simulated signal pairs)
  auto Z = [](double S, double B) { return S + B > 0. ? S / std::sqrt(S + B) : 0.; };
  std::vector<int> bestLine(nEl, -1), bestBox(nEl, -1);
  std::vector<bool> fewPairs(nEl, false);
  for (int i = 0; i < nEl; ++i) {
    for (size_t l = 0; l < lines[i].size(); ++l) {
      const Line& li = lines[i][l];
      if (bestLine[i] < 0 || Z(li.own, li.total - li.own) > Z(lines[i][bestLine[i]].own, lines[i][bestLine[i]].total - lines[i][bestLine[i]].own))
        bestLine[i] = l;
    }
    for (size_t k = 0; k < boxes.size(); ++k) {
      const Box& bo = boxes[k];
      if (bo.el != i || bo.sSim < 10.) continue;
      if (bestBox[i] < 0 || Z(bo.S, bo.B) > Z(boxes[bestBox[i]].S, boxes[bestBox[i]].B)) bestBox[i] = k;
    }
    if (bestBox[i] < 0) {   // under 10 simulated signal pairs everywhere: the box with the most
      fewPairs[i] = true;
      for (size_t k = 0; k < boxes.size(); ++k)
        if (boxes[k].el == i && (bestBox[i] < 0 || boxes[k].sSim > boxes[bestBox[i]].sSim)) bestBox[i] = k;
    }
  }
  int iKBox = -1;
  for (size_t k = 0; k < boxes.size(); ++k)
    if (boxes[k].el == iK && boxes[k].ea == eKa && boxes[k].eb == eLa) iKBox = k;

  // the numbers
  const double pairsDay = [&] { double v = 0.; for (int s = 0; s < kSources; ++s) v += pairMap[s]->GetEntries() * w[s]; return v; }();
  std::ostringstream txt, csv;
  txt << Form("Cascade coincidences of the steel study: two hits of one muon in two pixels, both within %g ns of the muon, one "
              "day at %.3g mu-/s: %s\n", promptNs, rate, base.c_str())
      << Form("  %zu beam files, %.6g simulated muons: %ld prompt hits, %ld pairs in two pixels (%ld in one pixel: none expected)\n",
              beamJobs.size(), nSim, nPrompt, nPairs, nSamePixel)
      << Form("  check: these prompt hits, weighted, against the trace's prompt view (without accidentals): largest deviation %.2g\n",
              maxDev)
      << Form("  per day: %.4g prompt hits, %.4g pairs in two pixels (%.3g per muon); resolution sigma = %g keV x sqrt(E / MeV)\n",
              mine->Integral(), pairsDay, pairsDay / nDay, sigma1MeV)
      << "  windows +-1.5 sigma; S the element's own X-ray entries (a pair: both hits), B everything else (pairs: with the "
         "accidental pairs);\n  f3: the mass fraction seen at 3 sigma in one day, S' / sqrt(S' + B + (u B)^2) = 3, u = 0 "
         "(statistics) and 1% (the background known to 1%)\n\n";
  txt << "element mass frac. | best single line: keV       S/day     B/day      S/B  S/sqrt(S+B)  f3 (stat)  f3 (1%) | best pair: "
         "keV x keV          S/day     B/day   [simulated S, B]      S/B  S/sqrt(S+B)  f3 (stat)  f3 (1%) | S/B x\n";
  csv << "element,mass_fraction,line_keV,line_S,line_B,line_S_over_B,line_z,line_f3_stat,line_f3_1pct,pair_hi_keV,pair_lo_keV,"
         "pair_S,pair_B,pair_S_sim,pair_B_sim,pair_accidental,pair_S_over_B,pair_z,pair_f3_stat,pair_f3_1pct,few_simulated_pairs\n";
  for (int i = 0; i < nEl; ++i) {
    const double f = steel[i].percent / 100.;
    const Line& li = lines[i][bestLine[i]];
    const double sS = li.own, sB = li.total - li.own;
    const Box& bo = boxes[bestBox[i]];
    txt << Form("  %-2s %11.3g | %18.2f %9.3g %9.3g %8.2g %12.3g %10.3g %8.3g | %8.2f x %7.2f %9.3g %9.3g [%7.0f %9.0f] %8.2g %12.3g "
                "%10.3g %8.3g | %6.3g%s\n", steel[i].name, f, li.e, sS, sB, sS / sB, Z(sS, sB), Fraction3(f, sS, sB, 0.),
                Fraction3(f, sS, sB, 0.01), bo.ea, bo.eb, bo.S, bo.B, bo.sSim, bo.bSim, bo.S / bo.B, Z(bo.S, bo.B),
                Fraction3(f, bo.S, bo.B, 0.), Fraction3(f, bo.S, bo.B, 0.01), (bo.S / bo.B) / (sS / sB),
                fewPairs[i] ? "  (under 10 simulated signal pairs)" : "");
    csv << Form("%s,%.6g,%.2f,%.6g,%.6g,%.4g,%.4g,%.4g,%.4g,%.2f,%.2f,%.6g,%.6g,%.0f,%.0f,%.4g,%.4g,%.4g,%.4g,%.4g,%d\n", steel[i].name,
                f, li.e, sS, sB, sS / sB, Z(sS, sB), Fraction3(f, sS, sB, 0.), Fraction3(f, sS, sB, 0.01), bo.ea, bo.eb, bo.S, bo.B,
                bo.sSim, bo.bSim, bo.acc, bo.S / bo.B, Z(bo.S, bo.B), Fraction3(f, bo.S, bo.B, 0.), Fraction3(f, bo.S, bo.B, 0.01),
                int(fewPairs[i]));
  }
  if (iKBox >= 0) {
    const Box& bo = boxes[iKBox];
    txt << Form("\nK K_alpha %.2f keV x K L_alpha %.2f keV (box %.1f-%.1f x %.1f-%.1f keV): S %.4g per day (%.0f simulated pairs), "
                "B %.4g per day (%.0f simulated pairs; accidentals %.3g)\n", eKa, eLa, bo.hiLo, bo.hiHi, bo.loLo, bo.loHi, bo.S, bo.sSim,
                bo.B, bo.bSim, bo.acc);
    txt << "  B by source (per day):";
    for (int s = 0; s < kSources; ++s) {
      const double v = (bo.n[s][0] + bo.n[s][1]) * w[s] - (s == iK ? bo.S : 0.);
      if (v > 0.) txt << Form(" %s %.3g", SourceName(s).c_str(), v);
    }
    txt << "\n";
  }
  txt << "\nall candidate pairs (each element's 8 strongest lines in the singles, two windows apart): element, keV x keV, S, B per "
         "day [simulated], S/B, S/sqrt(S+B)\n";
  for (const Box& bo : boxes)
    txt << Form("  %-2s %8.2f x %8.2f  %9.3g %9.3g [%7.0f %9.0f] %9.2g %9.3g\n", steel[bo.el].name, bo.ea, bo.eb, bo.S, bo.B, bo.sSim,
                bo.bSim, bo.S / bo.B, Z(bo.S, bo.B));
  std::ofstream(out + "/beam_steel_coincidence.txt") << txt.str();
  std::ofstream(out + "/beam_steel_coincidence.csv") << csv.str();
  printf("%s", txt.str().c_str());

  // the maps per day
  auto total = dynamic_cast<TH2F*>(pairMap[0]->Clone("map_total"));
  total->SetDirectory(nullptr);
  total->Reset();
  for (int s = 0; s < kSources; ++s) total->Add(pairMap[s], w[s]);
  auto spotAll = dynamic_cast<TH2F*>(spot[0][0]->Clone("spot_all"));
  spotAll->SetDirectory(nullptr);
  spotAll->Reset();
  auto gateAll = dynamic_cast<TH1D*>(gate[0][0]->Clone("gate_all"));
  gateAll->SetDirectory(nullptr);
  gateAll->Reset();
  for (int s = 0; s < kSources; ++s)
    for (int x = 0; x < 2; ++x) { spotAll->Add(spot[s][x], w[s]); gateAll->Add(gate[s][x], w[s]); }
  auto spotK = dynamic_cast<TH2F*>(spot[iK][1]->Clone("spot_K"));
  spotK->SetDirectory(nullptr);
  spotK->Scale(w[iK]);
  auto gateK = dynamic_cast<TH1D*>(gate[iK][1]->Clone("gate_K"));
  gateK->SetDirectory(nullptr);
  gateK->Scale(w[iK]);
  {
    TFile rf((out + "/beam_steel_coincidence.root").c_str(), "RECREATE");
    total->Write();
    for (int s = 0; s < kSources; ++s) {
      auto c = dynamic_cast<TH2F*>(pairMap[s]->Clone(("map_" + SourceName(s) + "_day").c_str()));
      c->Scale(w[s]);
      c->Write();
      delete c;
    }
    spotAll->Write(); spotK->Write(); gateAll->Write(); gateK->Write();
    TNamed("settings", Form("sigma1MeV %g threshold %g promptNs %g rate %g seconds %g (maps: pairs per day; x the lower energy)",
                            sigma1MeV, threshold, promptNs, rate, seconds)).Write();
  }

  // the plots
  const int ex = int(std::floor(std::log10(rate)));
  const std::string beamText = Form("%.2g#times10^{%d} #mu^{-}/s", rate / std::pow(10., ex), ex);
  {   // all pairs
    auto c = new TCanvas("c_map", "", 1300, 1150);
    c->SetLogx(); c->SetLogy(); c->SetLogz();
    c->SetLeftMargin(0.13); c->SetRightMargin(0.16); c->SetTopMargin(0.1); c->SetBottomMargin(0.1);
    total->SetTitle(Form("Two hits of one muon in two pixels within %g ns, one day at %s, 1 cm steel, CZT 4#times4 pixels;lower "
                         "energy (keV);higher energy (keV);pairs per day in a cell (half a resolution each way)", promptNs, beamText.c_str()));
    total->GetXaxis()->SetMoreLogLabels(); total->GetXaxis()->SetNoExponent();
    total->GetYaxis()->SetMoreLogLabels(); total->GetYaxis()->SetNoExponent();
    total->GetXaxis()->SetTitleOffset(1.2); total->GetYaxis()->SetTitleOffset(1.7); total->GetZaxis()->SetTitleOffset(1.4);
    total->GetXaxis()->SetRangeUser(threshold, 3000.); total->GetYaxis()->SetRangeUser(threshold, 3000.);
    total->SetMinimum(std::max(1., 1e-7 * total->GetMaximum()));
    total->Draw("colz");
    TBox kb(gLo, kLo, gHi, kHi);
    kb.SetFillStyle(0); kb.SetLineColor(kRed + 1); kb.SetLineWidth(2);
    kb.Draw();
    TLatex lab;
    lab.SetTextSize(0.028); lab.SetTextColor(kRed + 1); lab.SetTextAlign(32);
    lab.DrawLatex(0.9 * gLo, eKa, "K K_{#alpha}#timesK L_{#alpha}");
    // orientation: Fe's spots (its lines in both pixels) and the arc of Fe K_alpha shared by two pixels (E1 + E2 = E)
    lab.SetTextColor(kBlack); lab.SetTextAlign(12); lab.SetTextSize(0.024);
    lab.DrawLatex(285., 1330., "Fe K_{#alpha}#timesL_{#alpha}");
    lab.DrawLatex(100., 1330., "Fe K_{#alpha}#timesM_{#alpha}");
    lab.DrawLatex(100., 285., "Fe L_{#alpha}#timesM_{#alpha}");
    lab.SetTextAlign(11);
    lab.DrawLatex(25., 2350., "arcs E_{1}+E_{2} = a line's energy: one photon shared by two pixels");
    c->Print((out + "/coincidence_map.png").c_str());
    c->Print((out + "/coincidence_map.pdf").c_str());
  }
  {   // around K's spot: all pairs, K's own pairs
    auto c = new TCanvas("c_spot", "", 1800, 800);
    c->Divide(2, 1, 0.002, 0.002);
    TH2F* h[2] = {spotAll, spotK};
    const char* what[2] = {"All pairs, one day (box: K K_{#alpha}#timesK L_{#alpha}, #pm1.5#sigma)", "K's own X-ray pairs, one day"};
    for (int p = 0; p < 2; ++p) {
      c->cd(p + 1);
      gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.24); gPad->SetTopMargin(0.1); gPad->SetBottomMargin(0.12);
      h[p]->SetTitle(Form("%s;lower energy (keV);higher energy (keV);pairs per day per 1#times2.5 keV", what[p]));
      h[p]->GetZaxis()->SetTitleOffset(2.1); h[p]->GetYaxis()->SetTitleOffset(1.5);
      h[p]->Draw("colz");
      TBox kb(gLo, kLo, gHi, kHi);
      kb.SetFillStyle(0); kb.SetLineColor(kRed + 1); kb.SetLineWidth(2);
      kb.DrawClone();
    }
    c->cd();
    c->Print((out + "/coincidence_k_spot.png").c_str());
    c->Print((out + "/coincidence_k_spot.pdf").c_str());
  }
  {   // K K_alpha: singles and pairs gated on K L_alpha
    auto c = new TCanvas("c_gate", "", 1500, 1300);
    c->Divide(1, 2, 0.002, 0.002);
    auto sAll = dynamic_cast<TH1D*>(viewPrompt->Clone("gate_singles_all"));
    auto sK = dynamic_cast<TH1D*>(ownX[iK]->Clone("gate_singles_K"));
    for (TH1D* h : {sAll, sK}) { h->SetDirectory(nullptr); h->Rebin(4); h->GetXaxis()->SetRangeUser(640., 790.); }
    TH1D* all[2] = {sAll, gateAll};
    TH1D* own[2] = {sK, gateK};
    const char* title[2] = {Form("Singles (tagged beam, %g ns), one day;energy in one pixel (keV);entries per 2 keV per day", promptNs),
                            Form("Pairs whose other hit is in K L_{#alpha} (%.1f-%.1f keV), one day;energy of the higher hit (keV);"
                                 "pairs per 2 keV per day", gLo, gHi)};
    for (int p = 0; p < 2; ++p) {
      c->cd(p + 1);
      gPad->SetLeftMargin(0.1); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.09); gPad->SetBottomMargin(0.12);
      double topAll = 0., topOwn = 0.;
      for (int b = all[p]->FindBin(640.); b <= all[p]->FindBin(789.); ++b) {
        topAll = std::max(topAll, all[p]->GetBinContent(b));
        topOwn = std::max(topOwn, own[p]->GetBinContent(b));
      }
      // K's part scaled by a round factor to about a third of the total's height (it is far smaller)
      const double want = topOwn > 0. ? 0.35 * topAll / topOwn : 1.;
      const double factor = want > 1. ? std::pow(10., std::floor(std::log10(want))) : 1.;
      auto k = dynamic_cast<TH1D*>(own[p]->Clone(Form("own_scaled_%d", p)));
      k->SetDirectory(nullptr);
      k->Scale(factor);
      all[p]->SetTitle(title[p]);
      all[p]->SetLineColor(kBlack); all[p]->SetLineWidth(1);
      all[p]->SetMinimum(0.);
      all[p]->SetMaximum(1.65 * std::max(topAll, topOwn * factor));   // the data in the lower 60%, the legend above
      all[p]->GetXaxis()->SetTitleSize(0.05); all[p]->GetXaxis()->SetLabelSize(0.045);
      all[p]->GetYaxis()->SetTitleSize(0.05); all[p]->GetYaxis()->SetLabelSize(0.045); all[p]->GetYaxis()->SetTitleOffset(0.95);
      all[p]->Draw("hist");
      k->SetLineColor(kRed + 1); k->SetLineWidth(2);
      k->Draw("hist same");
      auto wl = new TLine();
      wl->SetLineStyle(2); wl->SetLineColor(kGray + 2);
      for (double x : {kLo, kHi}) wl->DrawLine(x, 0., x, 0.62 * all[p]->GetMaximum());
      auto leg = new TLegend(0.12, 0.68, 0.5, 0.89);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.045); leg->SetTextFont(42);
      leg->AddEntry(all[p], p == 0 ? "all entries" : "all pairs", "l");
      leg->AddEntry(k, Form("K's own X-rays%s", factor > 1. ? Form(" #times %g", factor) : ""), "l");
      leg->AddEntry(wl, Form("K K_{#alpha} window %.1f-%.1f keV", kLo, kHi), "l");
      leg->Draw();
      double wS = 0., wB = 0.;
      for (int b = 1; b <= all[p]->GetNbinsX(); ++b) {
        const double x = all[p]->GetBinCenter(b);
        if (x >= kLo && x <= kHi) { wS += own[p]->GetBinContent(b); wB += all[p]->GetBinContent(b); }
      }
      TLatex note;
      note.SetNDC(); note.SetTextSize(0.042); note.SetTextAlign(31); note.SetTextFont(42);
      note.DrawLatex(0.96, 0.84, "in the window, per day:");
      note.DrawLatex(0.96, 0.77, Form("K %.3g, the rest %.3g (S/B %.2g)", wS, wB - wS, wS / std::max(1e-30, wB - wS)));
    }
    c->cd();
    c->Print((out + "/coincidence_k_gate.png").c_str());
    c->Print((out + "/coincidence_k_gate.pdf").c_str());
  }
  // per element: S/B, and the mass fraction seen at 3 sigma
  auto Frame = [&](const char* name, const char* title, const char* yTitle, double lo, double hi) {   // drawn as axes: the title as text
    auto h = new TH1D(name, Form(";;%s", yTitle), nEl, 0., nEl);
    h->SetDirectory(nullptr);
    for (int i = 0; i < nEl; ++i) h->GetXaxis()->SetBinLabel(i + 1, steel[i].name);
    h->GetXaxis()->SetLabelSize(0.055);
    h->SetMinimum(lo); h->SetMaximum(hi);
    h->Draw("axis");
    TLatex t;
    t.SetNDC(); t.SetTextSize(0.036); t.SetTextAlign(22); t.SetTextFont(42);
    t.DrawLatex(0.53, 0.95, title);
    return h;
  };
  auto Range = [](const std::vector<double>& v, double& lo, double& hi) {
    lo = 1e300; hi = 0.;
    for (double x : v) if (x > 0.) { lo = std::min(lo, x); hi = std::max(hi, x); }
    lo = std::pow(10., std::floor(std::log10(lo)) - 0.3); hi = std::pow(10., std::ceil(std::log10(hi)) + 0.3);
  };
  {
    std::vector<double> sb1, sb2;
    for (int i = 0; i < nEl; ++i) {
      const Line& li = lines[i][bestLine[i]];
      sb1.push_back(li.own / (li.total - li.own));
      sb2.push_back(boxes[bestBox[i]].S / boxes[bestBox[i]].B);
    }
    std::vector<double> both = sb1;
    both.insert(both.end(), sb2.begin(), sb2.end());
    double lo, hi;
    Range(both, lo, hi);
    auto c = new TCanvas("c_sb", "", 1400, 850);
    c->SetLogy(); c->SetGridy();
    c->SetLeftMargin(0.1); c->SetRightMargin(0.03); c->SetTopMargin(0.1); c->SetBottomMargin(0.1);
    Frame("frame_sb", "S/B of each element in one day: its best single line and its best pair (windows of 1.5 resolutions)", "S/B", lo, hi);
    auto g1 = new TGraph(), g2 = new TGraph();
    for (int i = 0; i < nEl; ++i) {
      if (sb1[i] > 0.) g1->AddPoint(i + 0.5, sb1[i]);
      if (sb2[i] > 0. && !fewPairs[i]) g2->AddPoint(i + 0.5, sb2[i]);
    }
    g1->SetMarkerStyle(24); g1->SetMarkerSize(1.8); g1->SetMarkerColor(kBlack);
    g2->SetMarkerStyle(20); g2->SetMarkerSize(1.8); g2->SetMarkerColor(kRed + 1);
    g1->Draw("p same"); g2->Draw("p same");
    TLatex lab;
    lab.SetTextSize(0.022); lab.SetTextAlign(21); lab.SetTextColor(kRed + 1);
    for (int i = 0; i < nEl; ++i)   // under 10 simulated signal pairs (C, N: their L lines are below the threshold): none
      if (sb2[i] > 0. && !fewPairs[i]) lab.DrawLatex(i + 0.5, sb2[i] * 1.45, Form("%.0f#times%.0f", boxes[bestBox[i]].ea, boxes[bestBox[i]].eb));
      else lab.DrawLatex(i + 0.5, sb1[i] * 3., "none");
    auto leg = new TLegend(0.45, 0.77, 0.97, 0.89);   // framed: its markers stand over a column
    leg->SetBorderSize(1); leg->SetFillStyle(1001); leg->SetFillColor(kWhite); leg->SetTextSize(0.03); leg->SetTextFont(42);
    leg->AddEntry(g1, "best single line", "p");
    leg->AddEntry(g2, "best pair (keV#timeskeV; none: under 10 simulated pairs)", "p");
    leg->Draw();
    c->Print((out + "/coincidence_s_over_b.png").c_str());
    c->Print((out + "/coincidence_s_over_b.pdf").c_str());
  }
  {
    std::vector<double> f1, f1u, f2, f2u, nominal;
    for (int i = 0; i < nEl; ++i) {
      const double f = steel[i].percent / 100.;
      const Line& li = lines[i][bestLine[i]];
      const Box& bo = boxes[bestBox[i]];
      f1.push_back(Fraction3(f, li.own, li.total - li.own, 0.)); f1u.push_back(Fraction3(f, li.own, li.total - li.own, 0.01));
      f2.push_back(fewPairs[i] ? 0. : Fraction3(f, bo.S, bo.B, 0.)); f2u.push_back(fewPairs[i] ? 0. : Fraction3(f, bo.S, bo.B, 0.01));
      nominal.push_back(f);
    }
    std::vector<double> all = f1;
    for (const auto* v : {&f1u, &f2, &f2u, &nominal}) all.insert(all.end(), v->begin(), v->end());
    double lo, hi;
    Range(all, lo, hi);
    auto c = new TCanvas("c_limits", "", 1400, 850);
    c->SetLogy(); c->SetGridy();
    c->SetLeftMargin(0.1); c->SetRightMargin(0.03); c->SetTopMargin(0.1); c->SetBottomMargin(0.1);
    Frame("frame_limits", "The mass fraction each element needs for 3 sigma in one day: its best single line or its best pair",
          "mass fraction", lo, hi);
    struct Set { const std::vector<double>* v; int style, color; double dx; const char* label; };
    const std::vector<Set> sets = {{&f1, 24, kBlack, -0.24, "best single line, statistics"},
                                   {&f1u, 25, kBlack, -0.08, "best single line, background known to 1%"},
                                   {&f2, 20, kRed + 1, 0.08, "best pair, statistics"},
                                   {&f2u, 21, kRed + 1, 0.24, "best pair, background known to 1%"}};
    auto leg = new TLegend(0.12, 0.13, 0.56, 0.33);
    leg->SetBorderSize(0); leg->SetFillStyle(1001); leg->SetFillColor(kWhite); leg->SetTextSize(0.03); leg->SetTextFont(42);
    for (const Set& st : sets) {
      auto g = new TGraph();
      for (int i = 0; i < nEl; ++i) if ((*st.v)[i] > 0.) g->AddPoint(i + 0.5 + st.dx, (*st.v)[i]);
      g->SetMarkerStyle(st.style); g->SetMarkerSize(1.6); g->SetMarkerColor(st.color);
      g->Draw("p same");
      leg->AddEntry(g, st.label, "p");
    }
    TLine nl;
    nl.SetLineColor(kBlue + 1); nl.SetLineWidth(3);
    for (int i = 0; i < nEl; ++i) nl.DrawLine(i + 0.1, nominal[i], i + 0.9, nominal[i]);
    auto gn = new TGraph();
    gn->SetLineColor(kBlue + 1); gn->SetLineWidth(3);
    leg->AddEntry(gn, "the steel's own fraction (seen where a marker is below it)", "l");
    leg->Draw();
    TLatex none;
    none.SetTextSize(0.022); none.SetTextAlign(21); none.SetTextColor(kRed + 1);
    bool anyNone = false;
    for (int i = 0; i < nEl; ++i)
      if (f2[i] <= 0.) { none.DrawLatex(i + 0.5, f1[i] * 3., "none"); anyNone = true; }
    if (anyNone) {
      none.SetNDC(); none.SetTextAlign(33);
      none.DrawLatex(0.96, 0.87, "none: no usable pair (under 10 simulated pairs; C, N: their L lines are below 20 keV)");
    }
    c->Print((out + "/coincidence_limits.png").c_str());
    c->Print((out + "/coincidence_limits.pdf").c_str());
  }
  printf("written into %s: beam_steel_coincidence.txt, .csv, .root, coincidence_map, _k_spot, _k_gate, _s_over_b, _limits (.png, .pdf)\n",
         out.c_str());
}
