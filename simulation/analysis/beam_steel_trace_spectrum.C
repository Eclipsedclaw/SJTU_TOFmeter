// The pixel panel's spectra of one day of a mu- beam into a 1 cm plate of stainless steel with U, Th and K impurities
// (2026-11-10: can the CZT plane see the steel's elements, and the impurities U, Th and K?), from the jobs of
// condor/beam_steel_submit.sh in one folder: <dir>/beam (muons, macro/beam_steel_lyso.mac) and <dir>/decay (the plate's
// natural radioactivity, macro/beam_steel_decay.mac).
//
// Composition. The simulated plate holds every element of the steel, with U, Th and K at 0.1% by mass each so that
// enough muons catch on them. Each muon is reweighted by the element that caught it (cap_z) to the real steel: the
// muons an element catches follow the Z law (as in Geant4), share_i = x_i Z_i / sum_j x_j Z_j with x the atom fractions
// of the real mass percentages (Fe the balance; U, Th, K the arguments, as mass fractions). A muon caught by element i
// weighs w_i = N_day share_i / n_i, n_i the simulated muons it caught (every event is written, so they are counted).
// Natural radioactivity: each decay-band event is one decay of every member of a chain in secular equilibrium (U-238,
// Th-232, K-40, U-235); it weighs A T / N_events, A the activity of the impurity in the plate (plateGrams x the mass
// fraction x the specific activity: natural U 12355 Bq/g of U-238 and 568.5 of U-235, Th 4057, natural K 31.72 Bq/g).
//
// Entries: each pixel hit of the panel (czt_*) is smeared by sigma1MeV x sqrt(E / MeV) (CZT 6.36991 keV: FWHM/E 1.5% at
// 1 MeV, as the CZT studies) and kept from threshold (20 keV) up. Four views of the same simulation, per day:
//   on:     beam on, every hit (no timing); a hit at time t after its muon counts with 1 - t / T (T the measuring time:
//           the chance that a nucleus made at a random moment decays within it), 0 beyond T; plus the natural decays
//   prompt: beam on, a tagged beam (a muon counter in front: at 1.32e5 mu-/s the muons come 7.6 us apart on average),
//           only the hits within promptNs (10 ns) of their muon, plus the accidental hits of the other muons and of
//           the natural decays inside that window (the "on" rate x promptNs per muon)
//   after:  the day after the beam, beam off: the activation the day of beam left (a hit at t counts with t / T for
//           t < T, (2T - t) / T for T < t < 2T) and the natural decays
//   nobeam: no beam at all: the natural decays alone (the plate's own radioactivity, as a screening measurement sees
//           it; the room's background is not simulated)
// For each line (the muonic X-rays each element emits, from the simulation's truth, lines closer than one resolution
// merged; and the natural gamma lines) the window +-1.5 sigma: the entries of that line's own source (S), the others
// (B), the significance S/sqrt(S + B) (S/sqrt(B) where B dominates; finite where B is small), the contrast (the
// content at the line over the straight line between +-2.5 sigma) and, for U, Th and K, the mass fraction that one day
// would see at 3 sigma (S grows with the fraction, B does not). Counting statistics only: knowing B's shape to a
// fraction of S/B is the real limit wherever S/B is small.
// The muonic X-ray intensities are the cascade model's (Fe K_alpha ~0.93 per muon against a measured 0.745); no pile-up
// (at 1.32e5 /s a 25 mm pixel takes ~1.4e4 entries/s, under 1% within its 100 ns); the CZT's hole trapping (low tails
// under each line) is not simulated.
//
// Outputs in outDir (default <dir>/analysis): beam_steel_trace_spectrum.txt and .csv (the numbers), one plot per file:
//   trace_spectrum_on       the day's spectrum, beam on, no timing, 20-7000 keV, by what made the entries
//   trace_spectrum_prompt   the same with the prompt window (tagged beam)
//   trace_spectrum_after    the day after the beam: activation and the natural U, Th, K lines
//   trace_k_region          600-800 keV, prompt: K K_alpha (713 keV) under Mo's 3d -> 2p lines (694, 718 keV)
//   trace_u_th_region       5800-6700 keV, prompt: the K_alpha lines of U and Th
// Only reads the run; writes into outDir. Run in the work folder, after the jobs are back (rsync):
//   root -l -b -q 'analysis/beam_steel_trace_spectrum.C("output/czt/beam_steel_czt_trace")'
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
const double kBinKeV = 2., kMaxKeV = 8000.;
const int kBins = int(kMaxKeV / kBinKeV);

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

struct Element {   // the real steel, mass percent (U, Th, K from the arguments; Fe the balance)
  int z;
  const char* name;
  double a;
  double percent;
};

struct Sample {   // panel entries of the muons one element caught, per 2 keV
  TH1D *on = nullptr, *onX = nullptr, *prompt = nullptr, *promptX = nullptr, *after = nullptr, *afterX = nullptr;
  double events = 0.;
  std::map<long, double> lines;   // emitted target X-rays (truth): energy in 10 eV -> photons
  void Make(const std::string& n)
  {
    on = NewHist(n + "_on"); onX = NewHist(n + "_onx"); prompt = NewHist(n + "_prompt");
    promptX = NewHist(n + "_promptx"); after = NewHist(n + "_after"); afterX = NewHist(n + "_afterx");
  }
};

struct Source {   // a natural decay chain of the decay band
  int z, a;
  const char* name;
  TH1D* h = nullptr;
  double events = 0., activity = 0.;
};

struct NaturalLine {
  double e;
  const char* nuclide;
  int source;   // index into the sources
};

double WindowSum(const TH1D* h, double lo, double hi)
{
  double sum = 0.;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    const double x = h->GetBinCenter(b);
    if (x >= lo && x <= hi) sum += h->GetBinContent(b);
  }
  return sum;
}

double Contrast(const TH1D* h, double e0, double sigma)
{
  auto at = [h](double e) { return h->GetBinContent(h->GetXaxis()->FindFixBin(e)); };
  const double side = 0.5 * (at(e0 - 2.5 * sigma) + at(e0 + 2.5 * sigma));
  return side > 0. ? at(e0) / side : 0.;
}

std::string Sci(double v, int digits = 2)
{
  if (v == 0.) return "0";
  const int e = int(std::floor(std::log10(std::fabs(v))));
  char s[64];
  snprintf(s, sizeof s, "%.*f#times10^{%d}", digits - 1, v / std::pow(10., e), e);
  return s;
}

TH1D* Copy(const TH1D* h, const std::string& name, double scale, int rebin = 1)
{
  auto c = dynamic_cast<TH1D*>(h->Clone(name.c_str()));
  c->SetDirectory(nullptr);
  c->Scale(scale);
  if (rebin > 1) c->Rebin(rebin);
  return c;
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
}  // namespace

void beam_steel_trace_spectrum(const char* dir, const char* outDir = "", double fracU = 1e-7, double fracTh = 1e-7,
                               double fracK = 5e-5, double rate = 1.32e5, double seconds = 86400., double promptNs = 10.,
                               double sigma1MeV = 6.36991, double threshold = 20., double plateGrams = 3200.)
{
  gStyle->SetOptStat(0);
  const std::string base = dir;
  const std::string out = std::string(outDir).empty() ? base + "/analysis" : std::string(outDir);
  gSystem->mkdir(out.c_str(), true);
  std::ostringstream txt, csv;
  char line[700];
  auto Sigma = [&](double kev) { return sigma1MeV * std::sqrt(kev / 1000.); };

  // the real steel (mass percent): the user's 2026-11-10 composition, U, Th, K from the arguments, Fe the balance
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
  std::map<int, double> share;
  std::map<int, const Element*> byZ;
  for (const Element& el : steel) { share[el.z] = el.percent / el.a * el.z / zSum; byZ[el.z] = &el; }

  // the beam band
  const std::string beamDir = base + "/beam";
  const std::vector<int> beamJobs = FinishedJobs(beamDir, "beam");
  if (beamJobs.empty()) { printf("no finished beam jobs (beam_<k>.root with log/job<k>.done) in %s\n", beamDir.c_str()); return; }
  const std::string simComposition = MacroValue(Form("%s/log/job%d.mac", beamDir.c_str(), beamJobs[0]), "/tof/target/composition");
  const std::string panelName = MacroValue(Form("%s/log/job%d.mac", beamDir.c_str(), beamJobs[0]), "/control/alias bmPanel");
  const char* panel = panelName == "lyso" ? "LYSO 12.69 x 12.69 cm" : "CZT 10 x 10 cm";
  std::map<int, Sample> samples;
  for (const Element& el : steel) samples[el.z].Make(el.name);
  Sample other;
  other.Make("other");
  TRandom3 random(4357);
  const double tMax = seconds * 1e9;   // ns
  double nSim = 0.;
  long nHit = 0;
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
    std::vector<double>*e = nullptr, *t = nullptr, *muxE = nullptr;
    std::vector<int>*src = nullptr, *muxSrc = nullptr;
    ev->SetBranchStatus("*", 0);
    for (const char* b : {"cap_z", "czt_e", "czt_t", "czt_src", "mux_e", "mux_src"}) ev->SetBranchStatus(b, 1);
    ev->SetBranchAddress("cap_z", &capZ);
    ev->SetBranchAddress("czt_e", &e);
    ev->SetBranchAddress("czt_t", &t);
    ev->SetBranchAddress("czt_src", &src);
    ev->SetBranchAddress("mux_e", &muxE);
    ev->SetBranchAddress("mux_src", &muxSrc);
    for (Long64_t i = 0; i < ev->GetEntries(); ++i) {
      ev->GetEntry(i);
      auto it = samples.find(capZ);
      Sample& s = it != samples.end() ? it->second : other;
      s.events += 1.;
      for (size_t p = 0; p < muxE->size(); ++p)
        if ((*muxSrc)[p] == 0) s.lines[std::lround((*muxE)[p] * 1e5)] += 1.;   // mux_src: 0 cascade X-ray, 1 capture gamma
      for (size_t h = 0; h < e->size(); ++h) {
        const double eMeV = (*e)[h], tNs = std::max(0., (*t)[h]);
        if (eMeV <= 0. || tNs >= 2. * tMax) continue;
        const double kev = 1000. * eMeV + random.Gaus(0., sigma1MeV * std::sqrt(eMeV));
        if (kev < threshold) continue;
        const bool isX = (*src)[h] & 1;
        const double wOn = tNs < tMax ? 1. - tNs / tMax : 0.;
        const double wAfter = tNs < tMax ? tNs / tMax : (2. * tMax - tNs) / tMax;
        ++nHit;
        if (wOn > 0.) { s.on->Fill(kev, wOn); if (isX) s.onX->Fill(kev, wOn); }
        if (wAfter > 0.) { s.after->Fill(kev, wAfter); if (isX) s.afterX->Fill(kev, wAfter); }
        if (tNs < promptNs) { s.prompt->Fill(kev); if (isX) s.promptX->Fill(kev); }
      }
    }
  }
  const double nDay = rate * seconds;
  std::map<int, double> wDay;
  for (const Element& el : steel) {
    if (samples[el.z].events <= 0.) { printf("no muon caught by %s (Z %d) in %s: stopping\n", el.name, el.z, beamDir.c_str()); return; }
    wDay[el.z] = nDay * share[el.z] / samples[el.z].events;
  }
  const double wOther = nDay / nSim;

  // the decay band: one decay of every member of the chain per event
  std::vector<Source> sources = {{92, 238, "U-238 chain"}, {90, 232, "Th-232 chain"}, {19, 40, "K-40"}, {92, 235, "U-235 chain"}};
  sources[0].activity = plateGrams * fracU * 12355.;
  sources[1].activity = plateGrams * fracTh * 4057.;
  sources[2].activity = plateGrams * fracK * 31.72;
  sources[3].activity = plateGrams * fracU * 568.5;
  for (Source& so : sources) so.h = NewHist(Form("nat_%d_%d", so.z, so.a));
  const std::string decayDir = base + "/decay";
  for (int k : FinishedJobs(decayDir, "decay")) {
    const std::string mac = Form("%s/log/job%d.mac", decayDir.c_str(), k);
    const int z = std::atoi(MacroValue(mac, "/control/alias bmIonZ").c_str()), a = std::atoi(MacroValue(mac, "/control/alias bmIonA").c_str());
    Source* so = nullptr;
    for (Source& s : sources) if (s.z == z && s.a == a) so = &s;
    if (!so) { printf("decay_%d.root: a source Z %d A %d this analysis does not know; skipped\n", k, z, a); continue; }
    TFile f(Form("%s/decay_%d.root", decayDir.c_str(), k));
    auto run = dynamic_cast<TTree*>(f.Get("run"));
    auto ev = dynamic_cast<TTree*>(f.Get("events"));
    if (!run || !ev) { printf("skipped decay_%d.root: no run or events tree\n", k); continue; }
    run->Draw("n_events", "", "goff");
    for (Long64_t i = 0; i < run->GetSelectedRows(); ++i) so->events += run->GetV1()[i];
    std::vector<double>* e = nullptr;
    ev->SetBranchStatus("*", 0);
    ev->SetBranchStatus("czt_e", 1);
    ev->SetBranchAddress("czt_e", &e);
    for (Long64_t i = 0; i < ev->GetEntries(); ++i) {
      ev->GetEntry(i);
      for (double eMeV : *e) {
        if (eMeV <= 0.) continue;
        const double kev = 1000. * eMeV + random.Gaus(0., sigma1MeV * std::sqrt(eMeV));
        if (kev >= threshold) so->h->Fill(kev);
      }
    }
  }
  TH1D* natural = NewHist("natural");
  for (Source& so : sources) {
    if (so.events > 0.) so.h->Scale(so.activity * seconds / so.events);
    else printf("no finished decay job of the %s: the spectra leave it out\n", so.name);
    natural->Add(so.h);
  }

  // the views, per day
  auto Beam = [&](TH1D* Sample::*part, const std::string& name, int rebin) {
    TH1D* h = nullptr;
    for (const Element& el : steel) {
      TH1D* c = Copy(samples[el.z].*part, name + "_" + el.name, wDay[el.z], rebin);
      if (!h) h = c; else { h->Add(c); delete c; }
    }
    TH1D* o = Copy(other.*part, name + "_o", wOther, rebin);
    h->Add(o); delete o;
    return h;
  };
  TH1D* on = Beam(&Sample::on, "v_on", 1);
  on->Add(natural);
  TH1D* prompt = Beam(&Sample::prompt, "v_prompt", 1);
  TH1D* accidental = Copy(on, "v_accidental", rate * promptNs * 1e-9);
  prompt->Add(accidental);
  TH1D* after = Beam(&Sample::after, "v_after", 1);
  after->Add(natural);
  TH1D* nobeam = Copy(natural, "v_nobeam", 1.);
  struct View { const char* name; TH1D* total; TH1D* Sample::*own; double scaleNatural; };
  // the element's own X-rays in a view, and how much of the natural decays a view holds
  const std::vector<View> views = {{"on", on, &Sample::onX, 1.}, {"prompt", prompt, &Sample::promptX, rate * promptNs * 1e-9},
                                   {"after", after, &Sample::afterX, 1.}};

  // numbers
  double entriesOn = 0.;
  for (int b = 1; b <= on->GetNbinsX(); ++b) entriesOn += on->GetBinContent(b);
  txt << "Steel with U, Th, K on the " << panel << " panel (sigma " << sigma1MeV << " keV at 1 MeV, from " << threshold
      << " keV): " << base << "\n  beam: " << beamJobs.size() << " finished jobs, " << Form("%.6g", nSim)
      << " simulated mu- into 1 cm of " << simComposition << " (mass percent), " << nHit << " panel hits\n";
  snprintf(line, sizeof line, "  one day: %.3g mu-/s x %.0f s = %.4g muons (%.4g times the simulated); plate %.0f g\n", rate, seconds,
           nDay, nDay / nSim, plateGrams);
  txt << line << "  muons caught, simulated -> the real steel (Z law; weight per simulated muon), mass percent:\n";
  for (const Element& el : steel) {
    snprintf(line, sizeof line, "    %-2s (Z %2d) %9.4g%%: %9.0f simulated (%.3g%%) -> %.4g%% of the muons, %.4g per day (weight %.4g)\n",
             el.name, el.z, el.percent, samples[el.z].events, 100. * samples[el.z].events / nSim, 100. * share[el.z],
             share[el.z] * nDay, wDay[el.z]);
    txt << line;
  }
  snprintf(line, sizeof line, "    other (not caught in the plate's atoms, or no capture): %.0f simulated\n", other.events);
  txt << line;
  txt << "  natural radioactivity of the plate (secular equilibrium):\n";
  for (const Source& so : sources) {
    snprintf(line, sizeof line, "    %-13s %.4g Bq (%.4g decays per day), %.0f simulated chains\n", so.name, so.activity,
             so.activity * seconds, so.events);
    txt << line;
  }
  const double perPixel = entriesOn / seconds / 16.;
  snprintf(line, sizeof line, "  panel entries per day, beam on: %.4g (%.4g /s; %.3g per pixel and second: a second entry within "
           "+-100 ns %.2f%%, within +-1 us %.1f%%; no pile-up simulated)\n", entriesOn, entriesOn / seconds, perPixel,
           100. * (1. - std::exp(-2. * perPixel * 100e-9)), 100. * (1. - std::exp(-2. * perPixel * 1e-6)));
  txt << line;
  csv << "kind,element,line_keV,yield_per_capture,view,window_lo,window_hi,S,B,S_over_sqrt_SplusB,contrast,fraction,fraction_3sigma\n";

  // the X-ray lines of each element (truth), merged within one resolution, yield >= 1% per capture, above 50 keV
  txt << "\nmuonic X-ray lines (the cascade model's yields per captured muon; +-1.5 sigma windows; S its own X-ray entries "
         "per day, B everything else):\n";
  txt << "  element  line (keV) yield | on: S, B, S/sqrt(S+B), contrast | prompt: S, B, S/sqrt(S+B), contrast | U, Th, K: mass "
         "fraction seen at 3 sigma in a day (on, prompt)\n";
  for (const Element& el : steel) {
    Sample& s = samples[el.z];
    std::vector<std::pair<double, double>> raw;   // keV, photons per capture
    for (const auto& [k, n] : s.lines) raw.push_back({k / 100., n / s.events});
    std::vector<std::pair<double, double>> merged;
    for (const auto& [e0, y] : raw) {
      if (!merged.empty() && e0 - merged.back().first < Sigma(e0)) {
        auto& m = merged.back();
        m.first = (m.first * m.second + e0 * y) / (m.second + y);
        m.second += y;
      } else merged.push_back({e0, y});
    }
    std::vector<std::pair<double, double>> kept;
    for (const auto& m : merged) if (m.second >= 0.01 && m.first >= 50.) kept.push_back(m);
    std::sort(kept.begin(), kept.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
    if (kept.size() > 8) kept.resize(8);
    std::sort(kept.begin(), kept.end());
    const bool trace = el.z == 19 || el.z == 90 || el.z == 92;
    for (const auto& [e0, y] : kept) {
      const double sg = Sigma(e0), lo = e0 - 1.5 * sg, hi = e0 + 1.5 * sg;
      txt << Form("  %-2s %10.2f %6.3f", el.name, e0, y);
      double f3s[2] = {0., 0.};
      for (int v = 0; v < 2; ++v) {
        const View& vw = views[v];
        const double S = WindowSum(s.*(vw.own), lo, hi) * wDay[el.z];
        const double B = std::max(0., WindowSum(vw.total, lo, hi) - S);
        const double z = S + B > 0. ? S / std::sqrt(S + B) : 0.;
        const double c = Contrast(vw.total, e0, sg);
        txt << Form(" | %9.3g %9.3g %8.3g %5.2f", S, B, z, c);
        const double fr = el.percent / 100.;
        // the fraction whose S' gives S' / sqrt(S' + B) = 3 (B stays): S' = (9 + sqrt(81 + 36 B)) / 2
        const double f3 = trace && S > 0. ? fr * 0.5 * (9. + std::sqrt(81. + 36. * B)) / S : 0.;
        if (trace) f3s[v] = f3;
        csv << Form("xray,%s,%.2f,%.4f,%s,%.2f,%.2f,%.6g,%.6g,%.6g,%.4f,%.4g,%.4g\n", el.name, e0, y, vw.name, lo, hi, S, B, z, c, fr, f3);
      }
      if (trace) txt << Form(" | %.3g, %.3g", f3s[0], f3s[1]);
      txt << "\n";
    }
  }

  // the natural gamma lines
  const std::vector<NaturalLine> naturalLines = {
    {238.63, "Pb-212", 1}, {295.22, "Pb-214", 0}, {351.93, "Pb-214", 0}, {583.19, "Tl-208", 1}, {609.32, "Bi-214", 0},
    {911.20, "Ac-228", 1}, {968.97, "Ac-228", 1}, {1120.29, "Bi-214", 0}, {1460.82, "K-40", 2}, {1764.49, "Bi-214", 0},
    {2614.51, "Tl-208", 1}, {185.72, "U-235", 3}};
  txt << "\nnatural gamma lines of U, Th, K (S the line's chain in the window, B everything else; per day):\n";
  txt << "  line (keV) nuclide | on: S, B, S/sqrt(S+B) | after (beam off, next day): S, B, S/sqrt(S+B) | no beam: S, B, "
         "S/sqrt(S+B), contrast\n";
  for (const NaturalLine& nl : naturalLines) {
    const double sg = Sigma(nl.e), lo = nl.e - 1.5 * sg, hi = nl.e + 1.5 * sg;
    const double S = WindowSum(sources[nl.source].h, lo, hi);
    txt << Form("  %8.2f %-7s", nl.e, nl.nuclide);
    const std::vector<std::pair<const char*, TH1D*>> vv = {{"on", on}, {"after", after}, {"nobeam", nobeam}};
    for (const auto& [vn, h] : vv) {
      const double B = std::max(0., WindowSum(h, lo, hi) - S);
      const double z = S + B > 0. ? S / std::sqrt(S + B) : 0.;
      txt << Form(" | %9.3g %9.3g %8.3g", S, B, z);
      if (std::string(vn) == "nobeam") txt << Form(" %5.2f", Contrast(h, nl.e, sg));
      csv << Form("natural,%s,%.2f,,%s,%.2f,%.2f,%.6g,%.6g,%.6g,%.4f,,\n", nl.nuclide, nl.e, vn, lo, hi, S, B, z, Contrast(h, nl.e, sg));
    }
    txt << "\n";
  }

  // plots: the components of a view
  const std::string beamText = Form("%s #mu^{-}/s, 85 MeV/c", Sci(rate).c_str());
  struct Group { const char* label; std::vector<int> zs; int color; };
  const std::vector<Group> groups = {{"Fe", {26}, kRed}, {"Cr", {24}, kBlue}, {"Ni", {28}, kGreen + 2},
                                     {"Mo", {42}, kMagenta + 1}, {"Mn", {25}, kCyan + 2},
                                     {"C, N, Si, P, S", {6, 7, 14, 15, 16}, kGray + 2}, {"K", {19}, kOrange},
                                     {"Th", {90}, kAzure + 1}, {"U", {92}, kViolet}};
  auto GroupHist = [&](TH1D* Sample::*part, const Group& g, const std::string& tag, int rebin) {
    TH1D* h = nullptr;
    for (int z : g.zs) {
      TH1D* c = Copy(samples[z].*part, tag + "_" + std::to_string(z), wDay[z], rebin);
      if (!h) h = c; else { h->Add(c); delete c; }
    }
    h->SetLineColor(g.color); h->SetLineWidth(2);
    return h;
  };
  auto Spectrum = [&](const std::string& tag, TH1D* total, TH1D* Sample::*xpart, TH1D* Sample::*allpart, double natScale,
                      const std::string& title, double lo, double hi, int rebin, bool logY, bool withAccidentals,
                      double headroom) {
    TH1D* t = Copy(total, tag + "_total", 1., rebin);
    std::vector<TH1D*> parts;
    std::vector<std::string> labels;
    for (const Group& g : groups) { parts.push_back(GroupHist(xpart, g, tag + "_x" + g.label, rebin)); labels.push_back(std::string("muonic X-rays of ") + g.label); }
    TH1D* allBeam = Beam(allpart, tag + "_beam", rebin);
    TH1D* xBeam = Beam(xpart, tag + "_xbeam", rebin);
    allBeam->Add(xBeam, -1.);   // the beam's entries that are not muonic X-rays of the plate
    allBeam->SetLineColor(kOrange + 7); allBeam->SetLineWidth(2);
    parts.push_back(allBeam); labels.push_back("beam: decay e^{-}, capture #gamma, n, activation");
    if (withAccidentals) {   // the other muons' hits in the window (the natural decays' share drawn on its own)
      TH1D* acc = Copy(accidental, tag + "_acc", 1., rebin);
      TH1D* accNat = Copy(natural, tag + "_accnat", natScale, rebin);
      acc->Add(accNat, -1.); delete accNat;
      acc->SetLineColor(kGray + 1); acc->SetLineWidth(2); acc->SetLineStyle(2);
      parts.push_back(acc); labels.push_back(Form("accidentals: other muons' hits within %g ns", promptNs));
    }
    TH1D* nat = Copy(natural, tag + "_nat", natScale, rebin);
    nat->SetLineColor(kYellow + 3); nat->SetLineWidth(2);
    parts.push_back(nat); labels.push_back(withAccidentals ? "natural U, Th, K decays within the window" : "natural U, Th, K decays");
    auto c = new TCanvas(("c_" + tag).c_str(), "", 1100, 750);
    if (logY) c->SetLogy();
    c->SetLeftMargin(0.11); c->SetRightMargin(0.03); c->SetTopMargin(0.08);
    t->SetTitle((title + ";energy in one pixel (keV);" + Form("entries per %g keV in one day", rebin * kBinKeV)).c_str());
    t->GetXaxis()->SetRangeUser(lo, hi);
    t->SetLineColor(kBlack); t->SetLineWidth(2);
    double ymin = 1e300;
    for (int b = t->FindBin(lo); b <= t->FindBin(hi); ++b) if (t->GetBinContent(b) > 0.) ymin = std::min(ymin, t->GetBinContent(b));
    if (logY) { t->SetMaximum(t->GetMaximum() * headroom); t->SetMinimum(std::max(1e-3, ymin * 1e-4)); }
    else { t->SetMinimum(0.); t->SetMaximum(1.35 * t->GetMaximum()); }
    t->Draw("hist");
    for (TH1D* h : parts) { h->GetXaxis()->SetRangeUser(lo, hi); h->Draw("hist same"); }
    t->Draw("hist same");
    auto leg = new TLegend(0.56, 0.52, 0.97, 0.91);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.023);
    double sum = 0.;
    for (int b = t->FindBin(lo); b <= t->FindBin(hi); ++b) sum += t->GetBinContent(b);
    leg->AddEntry(t, Form("all entries (%s in %g-%g keV)", Sci(sum).c_str(), lo, hi), "l");
    for (size_t i = 0; i < parts.size(); ++i) leg->AddEntry(parts[i], labels[i].c_str(), "l");
    leg->Draw();
    return std::make_pair(c, t);
  };
  auto Mark = [&](TH1D* t, double e, const char* text, int color, double yFactor) {
    TLatex tex; tex.SetTextSize(0.026); tex.SetTextAlign(21); tex.SetTextColor(color);
    tex.DrawLatex(e, t->GetBinContent(t->FindBin(e)) * yFactor, text);
  };

  // 1. beam on, no timing
  {
    auto [c, t] = Spectrum("on", on, &Sample::onX, &Sample::on, 1., Form("One day at %s: 1 cm steel with U, Th, K, %s (no timing)",
                           beamText.c_str(), panel), 20., 7000., 4, true, false, 300.);
    Mark(t, 1257.2, "Fe K_{#alpha}", kRed, 3.);
    Mark(t, 2707.5, "Mo K_{#alpha}", kMagenta + 1, 3.);
    Save(c, out + "/trace_spectrum_on");
  }
  // 2. prompt window
  {
    auto [c, t] = Spectrum("prompt", prompt, &Sample::promptX, &Sample::prompt, rate * promptNs * 1e-9,
                           Form("One day at %s, tagged beam: hits within %g ns of their muon (1 cm steel with U, Th, K, %s)",
                                beamText.c_str(), promptNs, panel), 20., 7000., 4, true, true, 300.);
    Mark(t, 1257.2, "Fe K_{#alpha}", kRed, 3.);
    Mark(t, 2707.5, "Mo K_{#alpha}", kMagenta + 1, 3.);
    Save(c, out + "/trace_spectrum_prompt");
  }
  // 3. the day after: activation and natural lines
  {
    TH1D* t = Copy(after, "after_total", 1., 2);
    TH1D* act = Beam(&Sample::after, "after_act", 2);
    act->SetLineColor(kOrange + 7); act->SetLineWidth(2);
    std::vector<TH1D*> nat;
    const int colors[4] = {kViolet, kAzure + 1, kOrange, kMagenta};
    for (size_t i = 0; i < sources.size(); ++i) {
      TH1D* h = Copy(sources[i].h, Form("after_nat%zu", i), 1., 2);
      h->SetLineColor(colors[i]); h->SetLineWidth(2);
      nat.push_back(h);
    }
    auto c = new TCanvas("c_after", "", 1100, 750);
    c->SetLogy(); c->SetLeftMargin(0.11); c->SetRightMargin(0.03); c->SetTopMargin(0.08);
    t->SetTitle(Form("The day after one day at %s: activation of the steel and its natural U, Th, K (%s);energy in one pixel "
                     "(keV);entries per %g keV in one day", beamText.c_str(), panel, 2 * kBinKeV));
    t->GetXaxis()->SetRangeUser(20., 3000.);
    t->SetLineColor(kBlack); t->SetLineWidth(2);
    t->SetMaximum(t->GetMaximum() * 300.); t->SetMinimum(0.3);
    t->Draw("hist"); act->Draw("hist same");
    for (TH1D* h : nat) h->Draw("hist same");
    t->Draw("hist same");
    auto leg = new TLegend(0.50, 0.68, 0.97, 0.91);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.026);
    leg->AddEntry(t, "all entries (beam off)", "l");
    leg->AddEntry(act, "activation by the day of beam", "l");
    for (size_t i = 0; i < sources.size(); ++i)
      leg->AddEntry(nat[i], Form("%s (%.3g Bq)", sources[i].name, sources[i].activity), "l");
    leg->Draw();
    for (const NaturalLine& nl : naturalLines)
      if (nl.e == 609.32 || nl.e == 1460.82 || nl.e == 2614.51) Mark(t, nl.e, Form("%.0f %s", nl.e, nl.nuclide), colors[nl.source], 4.);
    Save(c, out + "/trace_spectrum_after");
  }
  // 4. and 5. the regions of the impurities' K lines, prompt window, linear
  for (int region = 0; region < 2; ++region) {
    const bool kRegion = region == 0;
    const double lo = kRegion ? 600. : 5800., hi = kRegion ? 800. : 6700.;
    const int rebin = kRegion ? 1 : 4;
    auto [c, t] = Spectrum(kRegion ? "kreg" : "ureg", prompt, &Sample::promptX, &Sample::prompt, rate * promptNs * 1e-9,
                           Form("One day at %s, tagged beam, hits within %g ns: %s (%s)", beamText.c_str(), promptNs,
                                kRegion ? "K K_{#alpha} under Mo 3d #rightarrow 2p" : "the K_{#alpha} lines of U and Th", panel),
                           lo, hi, rebin, true, true, 3e4);
    if (kRegion) {
      Mark(t, 712.7, "K K_{#alpha} 712.7", kOrange + 1, 12.);
      Mark(t, 694.4, "Mo 694.4", kMagenta + 1, 1.6);
      Mark(t, 718.4, "Mo 718.4", kMagenta + 1, 1.6);
    } else {
      Mark(t, 6133., "U K_{#alpha2}", kViolet, 3.);
      Mark(t, 6420., "U K_{#alpha1}", kViolet, 3.);
      Mark(t, 6103.6, "Th K_{#alpha2}", kAzure + 1, 8.);
      Mark(t, 6367.7, "Th K_{#alpha1}", kAzure + 1, 8.);
    }
    Save(c, out + (kRegion ? "/trace_k_region" : "/trace_u_th_region"));
  }

  std::ofstream(out + "/beam_steel_trace_spectrum.txt") << txt.str();
  std::ofstream(out + "/beam_steel_trace_spectrum.csv") << csv.str();
  printf("%s", txt.str().c_str());
  printf("written into %s: beam_steel_trace_spectrum.txt, .csv, trace_spectrum_on, trace_spectrum_prompt, trace_spectrum_after, "
         "trace_k_region, trace_u_th_region (.png, .pdf)\n", out.c_str());
}
