// The pixel panel's (LYSO, or the CZT plane) spectrum of one day of a mu- beam into a 1 cm plate of stainless steel 316 lying on it (2026-11-09,
// the beam study of steel: can the panel tell the Fe, Cr and Ni peaks apart), from the jobs of
// scripts/beam_steel_local.sh (macro/beam_steel_lyso.mac). Every simulated muon stands for N_day / N_sim real ones
// (N_day = beam rate x measuring time).
//
// Composition. The simulated plate is Geant4's G4_STAINLESS-STEEL (Fe, Cr and Ni only, like 304). Each muon is
// reweighted by the element that caught it (cap_z) to the shares of 316: the muons an element catches follow the Z law
// (as in Geant4), share_i = x_i Z_i / sum_j x_j Z_j with x the atom fractions from the weight percentages given (Fe the
// rest). The weight of a muon caught by element i is w_i = N_day share_i / n_i, with n_i the simulated muons it caught:
// every event is written (/EventFilter 0), so they are counted exactly. Mo, Mn, Si and C (about 5% of the muons) are
// not simulated: their share is left out of every spectrum.
//
// Entries, as in analysis/beam_impurity_spectrum.C: each pixel hit of the panel (czt_*) is smeared by
// sigma1MeV x sqrt(E / MeV) and kept from 20 keV up. A hit at time t after its muon counts with the weight 1 - t / T,
// and 0 beyond T (T = the measuring time): the chance that a nucleus made at a random moment decays within the run.
// The panel's own Lu-176 is added from a decay-only run (luActivity decays per second for T). No pile-up.
// Two views of the same simulation:
//   no timing: at 1.32e8 mu-/s a muon comes every 8 ns, so no window can pick out one muon's prompt hits: every hit
//              counts (the realistic spectrum of that beam; per day)
//   prompt:    only the hits within promptNs (10 ns) of the muon's arrival, as a muon counter in front and a beam slow
//              enough to see each muon (below ~1e6 /s) would allow (per muon of the beam)
// Line energies are the measured ones (/tof/physics/muonXrayEnergies measured; Cr interpolated). Their intensities are
// the cascade model's: Fe K_alpha 0.93 per muon (measured 0.745), K_beta 0.027 (0.075), the higher K lines about 0.02
// (0.18), so the K_alpha peaks come out about 25% too strong and the K lines above them too weak.
//
// Outputs in outDir (default <dir>/analysis): beam_steel_spectrum.txt (the numbers), and one plot per file:
//   beam_steel_spectrum_day   the day's spectrum (no timing), 20-3000 keV, by what made the entries
//   beam_steel_k_lines        the day's spectrum, 900-1800 keV: the K_alpha lines of Cr, Fe and Ni
//   beam_steel_l_lines        the day's spectrum, 150-450 keV: their L_alpha lines (3 -> 2)
//   beam_steel_prompt         entries per muon, every hit against the prompt hits only, 150-1800 keV
//   beam_steel_shapes         entries per caught muon for Cr, Fe and Ni, 150-1800 keV: where each element's lines are
// Only reads the run; writes into outDir. Run in the work folder:
//   root -l -b -q 'analysis/beam_steel_spectrum.C("output/lyso/beam_steel_lyso")'
// The panel comes from the jobs' macros (/control/alias bmPanel; none: lyso): with the CZT plane (--detector czt,
// 2026-11-09) the labels, its 16 pixels, sigma1MeV 6.36991 keV (FWHM/E 1.5% at 1 MeV, as the CZT studies) and no Lu-176.
// sigma1MeV < 0 and luActivity < 0 (the defaults) mean the panel's own: LYSO 26.85 keV and 12557 decays/s.
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
#include <tuple>
#include <vector>

namespace {
void Save(TCanvas* c, const std::string& base)
{
  c->Print((base + ".png").c_str());
  c->Print((base + ".pdf").c_str());
}

const double kBinKeV = 2., kMaxKeV = 4000.;

// the lines as simulated (keV): K_alpha and L_alpha (3 -> 2) of /tof/physics/muonXrayEnergies measured, K_beta = their sum
struct Element {
  int z;
  const char* name;
  double ka, la;
  int color;
};
const std::vector<Element> kElements = {{24, "Cr", 1090.0, 226.7, kBlue}, {26, "Fe", 1255.81, 266.62, kRed},
                                        {28, "Ni", 1429.4, 309.8, kGreen + 2}};

struct Sample {        // panel entries of the muons one element caught, per 2 keV
  TH1D* all = nullptr;      // every hit, with the time weight
  TH1D* xray = nullptr;     // hits with energy of a muonic X-ray of the target (source bit 1)
  TH1D* rest = nullptr;     // the others: decay electrons, capture gammas, neutrons, activation
  TH1D* late = nullptr;     // hits more than 1 ms after the muon (activation)
  TH1D* prompt = nullptr;   // hits within the prompt window (weight 1)
  TH1D* promptX = nullptr;  // of these, muonic X-rays
  double events = 0.;       // simulated muons this element caught
  double entries = 0.;      // time-weighted entries above the threshold
  double promptEntries = 0.;
  void Make(const std::string& name)
  {
    const int n = int(kMaxKeV / kBinKeV);
    all = new TH1D((name + "_all").c_str(), "", n, 0., kMaxKeV);
    xray = new TH1D((name + "_xray").c_str(), "", n, 0., kMaxKeV);
    rest = new TH1D((name + "_rest").c_str(), "", n, 0., kMaxKeV);
    late = new TH1D((name + "_late").c_str(), "", n, 0., kMaxKeV);
    prompt = new TH1D((name + "_prompt").c_str(), "", n, 0., kMaxKeV);
    promptX = new TH1D((name + "_promptx").c_str(), "", n, 0., kMaxKeV);
    for (TH1D* h : {all, xray, rest, late, prompt, promptX}) { h->Sumw2(); h->SetDirectory(nullptr); }
  }
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

// how much a line stands out: the content at the line against the straight line between the sidebands at +-2.5 sigma
// (1 = no bump at all; a shoulder on a falling slope gives little more than 1)
double Contrast(const TH1D* h, double e0, double sigma)
{
  auto at = [h](double e) { return h->GetBinContent(h->GetXaxis()->FindFixBin(e)); };
  return at(e0) / (0.5 * (at(e0 - 2.5 * sigma) + at(e0 + 2.5 * sigma)));
}
}  // namespace

void beam_steel_spectrum(const char* dir, const char* outDir = "", double wCr = 17., double wNi = 12., double wMo = 2.5,
                         double wMn = 2., double wSi = 0.75, double wC = 0.08, double rate = 1.32e8,
                         double seconds = 86400., double promptNs = 10., double sigma1MeV = -1.,
                         const char* decayFile = "output/lyso20_concept/decay/decay_1.root", double luActivity = -1.,
                         double threshold = 20.)
{
  gStyle->SetOptStat(0);
  const std::string base = dir;
  const std::string beamDir = base + "/beam";
  const std::string out = std::string(outDir).empty() ? base + "/analysis" : std::string(outDir);
  gSystem->mkdir(out.c_str(), true);
  std::ostringstream txt;
  char line[600];

  // the finished jobs: beam_<k>.root with log/job<k>.done
  std::map<int, std::string> files;
  if (void* d = gSystem->OpenDirectory(beamDir.c_str())) {
    while (const char* entry = gSystem->GetDirEntry(d)) {
      int k = 0;
      char tail[8] = "";
      if (sscanf(entry, "beam_%d.%5s", &k, tail) == 2 && std::string(tail) == "root" &&
          !gSystem->AccessPathName(Form("%s/log/job%d.done", beamDir.c_str(), k)))
        files[k] = beamDir + "/" + entry;
    }
    gSystem->FreeDirectory(d);
  }
  if (files.empty()) { printf("no finished jobs (beam_<k>.root with log/job<k>.done) in %s\n", beamDir.c_str()); return; }
  std::string detector;   // the panel of the jobs' macros, the same for all
  for (const auto& [k, name] : files) {
    std::ifstream mac(Form("%s/log/job%d.mac", beamDir.c_str(), k));
    std::string l, d = "lyso";
    while (std::getline(mac, l)) if (l.rfind("/control/alias bmPanel ", 0) == 0) d = l.substr(23);
    if (!detector.empty() && d != detector) { printf("jobs of two panels (%s, %s) in %s; stopping\n", detector.c_str(), d.c_str(), beamDir.c_str()); return; }
    detector = d;
  }
  const bool czt = detector == "czt";
  if (!czt && detector != "lyso") { printf("unknown panel %s in the jobs' macros; stopping\n", detector.c_str()); return; }
  const char* panel = czt ? "CZT 10 x 10 cm" : "LYSO 12.69 x 12.69 cm";
  const double nPixels = czt ? 16. : 400.;
  if (sigma1MeV < 0.) sigma1MeV = czt ? 6.36991 : 26.85;
  if (luActivity < 0.) luActivity = czt ? 0. : 12557.;

  // the 316 shares of the muons (Z law over the atom fractions)
  const double wFe = 100. - wCr - wNi - wMo - wMn - wSi - wC;
  const std::vector<std::tuple<int, const char*, double, double>> ss316 = {
    {26, "Fe", 55.845, wFe}, {24, "Cr", 51.996, wCr}, {28, "Ni", 58.693, wNi}, {42, "Mo", 95.95, wMo},
    {25, "Mn", 54.938, wMn}, {14, "Si", 28.086, wSi}, {6, "C", 12.011, wC}};
  double zSum = 0.;
  for (const auto& [z, name, a, w] : ss316) zSum += w / a * z;
  std::map<int, double> share;
  for (const auto& [z, name, a, w] : ss316) share[z] = w / a * z / zSum;

  std::map<int, Sample> samples;
  for (const Element& el : kElements) samples[el.z].Make(el.name);
  Sample other;
  other.Make("other");
  TRandom3 random(4357);
  const double tMax = seconds * 1e9;   // ns
  double nSim = 0.;
  long nHit = 0;
  for (const auto& [k, name] : files) {
    TFile f(name.c_str());
    auto run = dynamic_cast<TTree*>(f.Get("run"));
    auto ev = dynamic_cast<TTree*>(f.Get("events"));
    if (!run || !ev) { printf("skipped %s: no run or events tree\n", name.c_str()); continue; }
    run->Draw("n_events:event_filter", "", "goff");
    bool allWritten = true;
    for (Long64_t i = 0; i < run->GetSelectedRows(); ++i) {
      nSim += run->GetV1()[i];
      if (run->GetV2()[i] != 0) allWritten = false;
    }
    if (!allWritten) {
      printf("%s was run with /EventFilter other than 0: the captures of each element cannot be counted; stopping\n", name.c_str());
      return;
    }
    int capZ = 0;
    std::vector<double>*e = nullptr, *t = nullptr;
    std::vector<int>* src = nullptr;
    ev->SetBranchStatus("*", 0);
    for (const char* b : {"cap_z", "czt_e", "czt_t", "czt_src"}) ev->SetBranchStatus(b, 1);
    ev->SetBranchAddress("cap_z", &capZ);
    ev->SetBranchAddress("czt_e", &e);
    ev->SetBranchAddress("czt_t", &t);
    ev->SetBranchAddress("czt_src", &src);
    for (Long64_t i = 0; i < ev->GetEntries(); ++i) {
      ev->GetEntry(i);
      auto it = samples.find(capZ);
      Sample& s = it != samples.end() ? it->second : other;
      s.events += 1.;
      for (size_t h = 0; h < e->size(); ++h) {
        const double eMeV = (*e)[h], tNs = (*t)[h];
        if (eMeV <= 0. || tNs >= tMax) continue;
        const double kev = 1000. * eMeV + random.Gaus(0., sigma1MeV * std::sqrt(eMeV));
        if (kev < threshold) continue;
        const double w = 1. - std::max(0., tNs) / tMax;
        const bool isX = (*src)[h] & 1;
        ++nHit;
        s.all->Fill(kev, w);
        (isX ? s.xray : s.rest)->Fill(kev, w);
        if (tNs > 1e6) s.late->Fill(kev, w);
        s.entries += w;
        if (tNs < promptNs) {
          s.prompt->Fill(kev);
          if (isX) s.promptX->Fill(kev);
          s.promptEntries += 1.;
        }
      }
    }
  }
  for (const Element& el : kElements)
    if (samples[el.z].events <= 0.) { printf("no muon caught by %s (Z %d) in %s\n", el.name, el.z, beamDir.c_str()); return; }

  // the panel's own Lu-176 over the time T
  TH1D* lu = new TH1D("lu", "", int(kMaxKeV / kBinKeV), 0., kMaxKeV);
  lu->Sumw2(); lu->SetDirectory(nullptr);
  double nDecays = 0., luEntries = 0.;
  if (luActivity > 0.) {
    TFile f(decayFile);
    auto run = dynamic_cast<TTree*>(f.Get("run"));
    auto ev = dynamic_cast<TTree*>(f.Get("events"));
    if (run && ev) {
      run->Draw("n_events", "", "goff");
      for (Long64_t i = 0; i < run->GetSelectedRows(); ++i) nDecays += run->GetV1()[i];
      std::vector<double>* e = nullptr;
      std::vector<int>* src = nullptr;
      ev->SetBranchStatus("*", 0);
      ev->SetBranchStatus("czt_e", 1); ev->SetBranchStatus("czt_src", 1);
      ev->SetBranchAddress("czt_e", &e); ev->SetBranchAddress("czt_src", &src);
      for (Long64_t i = 0; i < ev->GetEntries(); ++i) {
        ev->GetEntry(i);
        for (size_t h = 0; h < e->size(); ++h) {
          if (!((*src)[h] & 16) || (*e)[h] <= 0.) continue;
          const double kev = 1000. * (*e)[h] + random.Gaus(0., sigma1MeV * std::sqrt((*e)[h]));
          if (kev < threshold) continue;
          lu->Fill(kev);
          luEntries += 1.;
        }
      }
    } else printf("no Lu-176 decay-only run in %s: the spectra leave it out\n", decayFile);
  }
  if (nDecays > 0.) lu->Scale(luActivity * seconds / nDecays);

  // weights: per day (no timing) and per muon of the beam (prompt)
  const double nDay = rate * seconds;
  std::map<int, double> wDay;
  for (const Element& el : kElements) wDay[el.z] = nDay * share[el.z] / samples[el.z].events;
  const double wOther = nDay / nSim;
  double simulatedShare = 0.;
  for (const Element& el : kElements) simulatedShare += share[el.z];
  auto Sum = [&](TH1D* Sample::*part, const std::string& name, double perDay, int rebin) {   // the elements' sum
    TH1D* h = nullptr;
    for (const Element& el : kElements) {
      TH1D* c = Copy(samples[el.z].*part, name + "_" + el.name, wDay[el.z] * perDay, rebin);
      if (!h) h = c; else { h->Add(c); delete c; }
    }
    TH1D* o = Copy(other.*part, name + "_o", wOther * perDay, rebin);
    h->Add(o); delete o;
    return h;
  };
  auto Total = [&](const std::string& name, int rebin) {   // the day's spectrum, no timing, with Lu-176
    TH1D* h = Sum(&Sample::all, name, 1., rebin);
    TH1D* l = Copy(lu, name + "_lu", 1., rebin);
    h->Add(l); delete l;
    return h;
  };

  // numbers
  txt << "Beam study of stainless steel 316 on the " << panel << " panel (sigma " << sigma1MeV << " keV at 1 MeV): " << base << " (" << files.size() << " finished jobs), "
      << Form("%.6g", nSim) << " simulated mu- into 1 cm of G4_STAINLESS-STEEL, " << nHit << " panel entries of " << threshold
      << " keV and up\n";
  txt << "  muons caught, simulated (share) -> 316 (Z law over the atom fractions; weight per simulated muon for one day):\n";
  for (const Element& el : kElements) {
    snprintf(line, sizeof line, "    %s (Z %d): %.0f (%.2f%%) -> %.2f%% (%.4g)\n", el.name, el.z, samples[el.z].events,
             100. * samples[el.z].events / nSim, 100. * share[el.z], wDay[el.z]);
    txt << line;
  }
  snprintf(line, sizeof line, "    other (%.0f muons), and Mo %.2f%%, Mn %.2f%%, Si %.2f%%, C %.2f%% of the 316 muons not simulated (left out: %.2f%%)\n",
           other.events, 100. * share[42], 100. * share[25], 100. * share[14], 100. * share[6], 100. * (1. - simulatedShare));
  txt << line;
  snprintf(line, sizeof line, "one day: %.3g mu-/s x %.0f s = %.4g muons, %.4g times the simulated number\n", rate, seconds, nDay, nDay / nSim);
  txt << line;
  double entriesDay = 0., promptPerMuon = 0.;
  for (const Element& el : kElements) {
    entriesDay += samples[el.z].entries * wDay[el.z];
    promptPerMuon += samples[el.z].promptEntries * wDay[el.z] / nDay;
  }
  const double luPerSecond = nDecays > 0. ? luEntries / nDecays * luActivity : 0.;
  const double perPixel = (entriesDay / seconds + luPerSecond) / nPixels;
  snprintf(line, sizeof line, "  panel entries per muon (20 keV and up): %.4f with every hit (time-weighted), %.4f within %g ns; %.4g entries/s, %.4g per pixel and second\n",
           entriesDay / nDay, promptPerMuon, promptNs, entriesDay / seconds + luPerSecond, perPixel);
  txt << line;
  snprintf(line, sizeof line, "  pile-up left out: a pixel's entry has another within +-100 ns with probability %.1f%%\n",
           100. * (1. - std::exp(-2. * perPixel * 100e-9)));
  txt << line;

  // the lines: what fills a window of +-1.5 sigma around each (no timing: per day; prompt: per muon)
  TH1D* tot = Total("total_day", 1);
  TH1D* totPrompt = Sum(&Sample::prompt, "total_prompt", 1. / nDay, 1);
  const double sFe = sigma1MeV * std::sqrt(1.25581), sFeL = sigma1MeV * std::sqrt(0.26662);
  snprintf(line, sizeof line, "the lines (+-1.5 sigma windows; the K_alpha lines are %.1f-%.1f sigma apart, the L_alpha lines %.1f-%.1f):\n",
           (1255.81 - 1090.) / sFe, (1429.4 - 1255.81) / sFe, (266.62 - 226.7) / sFeL, (309.8 - 266.62) / sFeL);
  txt << line;
  txt << "  line         window (keV) | no timing, one day: entries | that element's X-rays | the other two's X-rays | not X-rays (activation) | Lu-176 | S/sqrt(B) | prompt (< " << promptNs << " ns): entries per muon | that element's X-rays\n";
  for (int series = 0; series < 2; ++series) {
    for (const Element& el : kElements) {
      const double e0 = series == 0 ? el.ka : el.la;
      const double half = 1.5 * sigma1MeV * std::sqrt(e0 / 1000.);
      const double lo = e0 - half, hi = e0 + half;
      const double all = WindowSum(tot, lo, hi);
      double own = 0., others = 0., rest = 0., late = 0., ownP = 0.;
      for (const Element& o : kElements) {
        const double x = WindowSum(samples[o.z].xray, lo, hi) * wDay[o.z];
        (o.z == el.z ? own : others) += x;
        rest += WindowSum(samples[o.z].rest, lo, hi) * wDay[o.z];
        late += WindowSum(samples[o.z].late, lo, hi) * wDay[o.z];
        if (o.z == el.z) ownP = WindowSum(samples[o.z].promptX, lo, hi) * wDay[o.z] / nDay;
      }
      const double luW = WindowSum(lu, lo, hi);
      const double allP = WindowSum(totPrompt, lo, hi);
      snprintf(line, sizeof line, "  %s %-8s %6.1f-%6.1f | %10.4g | %5.1f%% | %5.1f%% | %5.1f%% (%4.1f%%) | %5.2f%% | %9.3g | %10.4g | %5.1f%%\n",
               el.name, series == 0 ? "K_alpha" : "L_alpha", lo, hi, all, 100. * own / all, 100. * others / all,
               100. * rest / all, 100. * late / all, 100. * luW / all, own / std::sqrt(std::max(1., all - own)), allP,
               100. * ownP / allP);
      txt << line;
    }
  }
  {
    const int cRebin = czt ? 1 : 4;   // bins of a third of sigma or less
    const bool lApart = 4.5 * sFeL < 40.;   // sidebands clear of the neighbouring L_alpha peaks
    TH1D* t8 = Total("total_day_8", cRebin);
    TH1D* p8 = Sum(&Sample::prompt, "total_prompt_8", 1. / nDay, cRebin);
    txt << "  how much each line stands out of the whole spectrum (" << cRebin * kBinKeV << " keV bins): its content / the straight line between\n"
        << (lApart ? "  +-2.5 sigma (K_alpha, then L_alpha)\n" : "  +-2.5 sigma (not for the L_alpha lines: their neighbours lie within 2.5 sigma)\n");
    for (int series = 0; series < (lApart ? 2 : 1); ++series) {
      txt << "   ";
      for (const Element& el : kElements) {
        const double e0 = series == 0 ? el.ka : el.la, s = sigma1MeV * std::sqrt(e0 / 1000.);
        snprintf(line, sizeof line, "  %s %.2f (prompt %.2f)", el.name, Contrast(t8, e0, s), Contrast(p8, e0, s));
        txt << line;
      }
      txt << "\n";
    }
    delete t8; delete p8;
  }

  const std::string beamText = Form("%s #mu^{-}/s, 85 MeV/c", Sci(rate).c_str());
  // the lines' names above the curve h; on the wide log plots the three L_alpha lines (40 keV apart) get one name
  auto Labels = [&](TH1D* h, double lo, double hi, bool kSeries, bool lSeries, double yFactor, bool logY) {
    TLatex tex; tex.SetTextSize(0.028); tex.SetTextAlign(21);
    TLine mark; mark.SetLineStyle(3); mark.SetLineColor(kGray + 2);
    for (const Element& el : kElements) {
      for (int series = 0; series < 2; ++series) {
        if ((series == 0 && !kSeries) || (series == 1 && !lSeries)) continue;
        const double x = series == 0 ? el.ka : el.la;
        if (x < lo || x > hi || (logY && series == 1 && el.z != 26)) continue;
        const double y = h->GetBinContent(h->FindBin(x)) * yFactor;
        tex.SetTextColor(logY && series == 1 ? kBlack : el.color);
        tex.DrawLatex(x, y, logY && series == 1 ? "L_{#alpha}: Cr, Fe, Ni" : Form("%s %s", el.name, series == 0 ? "K_{#alpha}" : "L_{#alpha}"));
        if (!logY) mark.DrawLine(x, 0., x, 0.97 * y);
      }
    }
  };
  auto Components = [&](const std::string& tag, int rebin, double perDay) {   // the elements' X-rays and the rest
    std::vector<TH1D*> v;
    for (const Element& el : kElements) {
      TH1D* x = Copy(samples[el.z].xray, tag + "_x_" + el.name, wDay[el.z] * perDay, rebin);
      x->SetLineColor(el.color); x->SetLineWidth(2);
      v.push_back(x);
    }
    TH1D* r = Sum(&Sample::rest, tag + "_rest", perDay, rebin);
    r->SetLineColor(kOrange + 7); r->SetLineWidth(2);
    v.push_back(r);
    return v;
  };

  // 1. the day's spectrum, no timing, by what made the entries
  {
    TH1D* total = Total("day_total", 2);
    auto parts = Components("day", 2, 1.);
    TH1D* l = Copy(lu, "day_lu", 1., 2);
    l->SetLineColor(kViolet); l->SetLineWidth(2);
    auto c = new TCanvas("c_day", "", 1000, 700);
    c->SetLogy(); c->SetLeftMargin(0.11); c->SetRightMargin(0.03); c->SetTopMargin(0.08);
    total->SetTitle(Form("One day at %s: 1 cm stainless steel 316, %s;energy in one pixel (keV);"
                         "entries per %g keV in one day", beamText.c_str(), panel, 2. * kBinKeV));
    total->GetXaxis()->SetRangeUser(20., 3000.);
    total->SetMinimum(1e6); total->SetMaximum(total->GetMaximum() * 3e2);
    total->SetLineColor(kBlack); total->SetLineWidth(2);
    total->Draw("hist");
    for (TH1D* h : parts) h->Draw("hist same");
    if (luActivity > 0.) l->Draw("hist same");
    total->Draw("hist same");
    {
      TLatex g; g.SetTextSize(0.026); g.SetTextAlign(21); g.SetTextColor(kOrange + 7);
      g.DrawLatex(846.8, total->GetBinContent(total->FindBin(846.8)) * 2.5, "^{56}Fe 847 keV #gamma");
    }
    auto leg = new TLegend(0.47, 0.66, 0.97, 0.91);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.028);
    leg->AddEntry(total, Form("all entries (%s per day)", Sci(entriesDay + luPerSecond * seconds).c_str()), "l");
    for (size_t i = 0; i < kElements.size(); ++i)
      leg->AddEntry(parts[i], Form("muonic X-rays of %s (%.0f%% of the muons)", kElements[i].name, 100. * share[kElements[i].z]), "l");
    leg->AddEntry(parts.back(), "decay e^{-}, capture #gamma, n, activation", "l");
    if (luActivity > 0.) leg->AddEntry(l, "the LYSO's own Lu-176", "l");
    leg->Draw();
    Labels(total, 20., 3000., true, true, 2.5, true);
    Save(c, out + "/beam_steel_spectrum_day");
  }

  // 2. and 3. the K_alpha and L_alpha regions of the day's spectrum, linear
  for (int series = 0; series < 2; ++series) {
    const bool k = series == 0;
    const double lo = k ? 900. : 150., hi = k ? 1800. : 450.;
    const int rebin = czt ? 1 : (k ? 4 : 2);
    TH1D* total = Total(k ? "k_total" : "l_total", rebin);
    auto parts = Components(k ? "k" : "l", rebin, 1.);
    auto c = new TCanvas(k ? "c_k" : "c_l", "", 1000, 700);
    c->SetLeftMargin(0.11); c->SetRightMargin(0.03); c->SetTopMargin(0.08);
    total->SetTitle(Form("One day at %s: 1 cm stainless steel 316, %s;energy in one pixel (keV);"
                         "entries per %g keV in one day", beamText.c_str(), panel, rebin * kBinKeV));
    total->GetXaxis()->SetRangeUser(lo, hi);
    total->SetMinimum(0.);
    total->SetMaximum(1.35 * total->GetMaximum());
    total->SetLineColor(kBlack); total->SetLineWidth(3);
    total->GetYaxis()->SetTitleOffset(1.3);
    total->Draw("hist");
    for (TH1D* h : parts) { h->GetXaxis()->SetRangeUser(lo, hi); h->Draw("hist same"); }
    total->Draw("hist same");
    auto leg = new TLegend(0.62, 0.62, 0.97, 0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.03);
    leg->AddEntry(total, "all entries (no timing)", "l");
    for (size_t i = 0; i < kElements.size(); ++i)
      leg->AddEntry(parts[i], Form("muonic X-rays of %s", kElements[i].name), "l");
    leg->AddEntry(parts.back(), "decay e^{-}, capture #gamma, n, activation", "l");
    leg->Draw();
    Labels(total, lo, hi, k, !k, 1.06, false);
    if (k) {   // the gamma line under Ni K_alpha: 52V decay (after capture on 52Cr) and 52Cr(n,n')
      TLatex g; g.SetTextSize(0.026); g.SetTextAlign(21); g.SetTextColor(kOrange + 7);
      const double y = czt ? 1.45 * total->GetBinContent(total->FindBin(1434.1))   // narrow peaks: above Ni K_alpha's name
                           : 0.55 * parts.back()->GetBinContent(parts.back()->FindBin(1434.1));
      g.DrawLatex(1434.1, y, "1434 keV #gamma (^{52}V, ^{52}Cr)");
    }
    Save(c, out + (k ? "/beam_steel_k_lines" : "/beam_steel_l_lines"));
  }

  // 4. entries per muon: every hit against the prompt hits only
  {
    TH1D* every = Sum(&Sample::all, "pm_all", 1. / nDay, 2);
    TH1D* prompt = Sum(&Sample::prompt, "pm_prompt", 1. / nDay, 2);
    TH1D* l = Copy(lu, "pm_lu", 1. / nDay, 2);
    every->Add(l);
    auto c = new TCanvas("c_prompt", "", 1000, 700);
    c->SetLogy(); c->SetLeftMargin(0.11); c->SetRightMargin(0.03); c->SetTopMargin(0.08);
    every->SetTitle(Form("Entries per muon: every hit against those within %g ns of the muon (1 cm steel 316, %s);"
                         "energy in one pixel (keV);entries per %g keV per muon", promptNs, panel, 2. * kBinKeV));
    every->GetXaxis()->SetRangeUser(150., 1800.);
    prompt->GetXaxis()->SetRangeUser(150., 1800.);
    every->SetMaximum(every->GetMaximum() * 8.);
    every->SetMinimum(std::max(1e-7, 0.3 * prompt->GetMinimum(1e-12)));
    every->SetLineColor(kBlack); every->SetLineWidth(2);
    prompt->SetLineColor(kRed + 1); prompt->SetLineWidth(2);
    every->Draw("hist"); prompt->Draw("hist same");
    auto leg = new TLegend(0.40, 0.76, 0.97, 0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
    leg->AddEntry(every, Form("every hit, as at %s (no timing)", beamText.c_str()), "l");
    leg->AddEntry(prompt, Form("hits within %g ns: a tagged beam below ~10^{6} #mu^{-}/s", promptNs), "l");
    leg->Draw();
    Labels(every, 150., 1800., true, true, 1.5, true);
    Save(c, out + "/beam_steel_prompt");
  }

  // 5. entries per caught muon for each element: where its lines are
  {
    std::vector<TH1D*> per;
    double ymax = 0.;
    for (const Element& el : kElements) {
      TH1D* h = Copy(samples[el.z].all, std::string("per_") + el.name, 1. / samples[el.z].events, 2);
      h->GetXaxis()->SetRangeUser(150., 1800.);
      h->SetLineColor(el.color); h->SetLineWidth(2);
      ymax = std::max(ymax, h->GetMaximum());
      per.push_back(h);
    }
    auto c = new TCanvas("c_shapes", "", 1000, 700);
    c->SetLogy(); c->SetLeftMargin(0.11); c->SetRightMargin(0.03); c->SetTopMargin(0.08);
    per[0]->SetTitle(Form("Entries per caught muon, every hit: 1 cm steel, beam 85 MeV/c, %s;"
                          "energy in one pixel (keV);entries per %g keV per caught muon", panel, 2. * kBinKeV));
    per[0]->SetMaximum(ymax * 8.);
    per[0]->SetMinimum(1e-5);
    per[0]->Draw("hist");
    for (TH1D* h : per) h->Draw("hist same");
    auto leg = new TLegend(0.74, 0.72, 0.97, 0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
    for (size_t i = 0; i < kElements.size(); ++i) leg->AddEntry(per[i], Form("caught by %s", kElements[i].name), "l");
    leg->Draw();
    TLatex tex; tex.SetTextSize(0.028); tex.SetTextAlign(21);
    for (size_t i = 0; i < kElements.size(); ++i) {
      for (double x : {kElements[i].ka, kElements[i].la, kElements[i].ka + kElements[i].la}) {
        tex.SetTextColor(kElements[i].color);
        const std::string label = x == kElements[i].ka ? "K_{#alpha}" : x == kElements[i].la ? "L_{#alpha}" : "K_{#beta}";
        if (x <= 1800.) tex.DrawLatex(x, per[i]->GetBinContent(per[i]->FindBin(x)) * 1.6, label.c_str());
      }
    }
    Save(c, out + "/beam_steel_shapes");
  }

  std::ofstream(out + "/beam_steel_spectrum.txt") << txt.str();
  printf("%s", txt.str().c_str());
  printf("written into %s: beam_steel_spectrum.txt, beam_steel_spectrum_day, beam_steel_k_lines, beam_steel_l_lines, "
         "beam_steel_prompt, beam_steel_shapes (.png, .pdf)\n", out.c_str());
}
