// The LYSO panel's spectrum of one day of a mu- beam into a Si target with a trace of Al (2026-11-07, the beam study of
// Al in Si), from the jobs of condor/beam_xray_submit.sh. Every simulated muon stands for N_day / N_sim real ones
// (N_day = beam rate x measuring time).
//
// Reweighting to a real concentration c. The simulated target holds the impurity at an atom fraction fSim (1%). For a
// dilute impurity, the muons it caught (cap_z = its Z) and those the host caught do not change with c, so each kind
// is reweighted:
//   w_imp(c) = N_day P(c) / (N_sim P(fSim)),   w_host(c) = N_day (1 - P(c)) / (N_sim (1 - P(fSim))),
//   P(f) = f A / (1 - f + f A): the chance that a muon is caught by the impurity.
// A is the capture strength per impurity atom against a host atom, from the written events:
// A = (n_imp / n_host) (1 - fSim) / fSim. Geant4 takes it from the Z law, so A = 13/14 for Al in Si.
//
// Entries. Each pixel hit of the panel (czt_*) is smeared by sigma1MeV x sqrt(E / MeV) and kept from 20 keV up.
// A hit at time t after its muon counts with the weight 1 - t / T, and 0 beyond T (T = the measuring time): the
// chance that a nucleus made at a random moment of the run decays within it. Prompt hits get 1; the minutes-long
// decays of 28Al and 27Mg about 1; 24Na (15 h) about 0.6; 26Al 0.
// At these rates every hit counts: a muon comes every few ns, so no coincidence window can pick out one. The panel's
// own Lu-176 counts too: its entries come from a decay-only run, at luActivity decays per second for the time T.
// Each muon is alone in its event (no pile-up). The text output estimates the pile-up from the entry rate.
//
// Signal and noise. In the Al K_alpha window (+-2 sigma) the signal is the excess of (Si with c Al) over Si alone,
// S = N_day P(c) (s_imp - s_host), with s the window's entries per caught muon. Its noise is sqrt(B) of the day's
// Si-alone entries B. Z = S / sqrt(B) for statistics alone. Z = S / sqrt(B + (delta B)^2) if the shape of the
// Si-alone spectrum is known only to a fraction delta.
//
// Line energies are those of the cascade model (TOFmeterEmCaptureCascade): Al K_alpha 354.9 and Si K_alpha 411.7 keV.
// The measured ones are 2-3% lower: 346.83 and 400.18 keV. With a Lu K_alpha X-ray (54.07 keV) lost from the pixel,
// Si K_alpha appears at 357.6 keV in the simulation, and at 346.1 keV in reality: on top of Al K_alpha either way.
//
// Outputs in outDir (default <dir>/analysis): beam_impurity_spectrum.txt (the numbers), and one plot per file:
//   beam_spectrum_day     the day's spectrum, 20-1000 keV, Si with 1 ppm Al, by what made the entries
//   beam_shapes           entries per caught muon, Si against Al, 250-500 keV: where the lines are
//   beam_residual         (Si with c Al) - (Si alone) in units of the day's noise, per 2 keV, 250-500 keV:
//                         a pseudo-measurement with 1 ppm, and the expected excess for 1, 10, 100 ppm
//   beam_significance     Z of the Al K_alpha window against the concentration
// Only reads the run; writes into outDir. Run in the work folder:
//   root -l -b -q 'analysis/beam_impurity_spectrum.C("output/lyso/beam_xray_lyso_si_al0.01")'
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
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

// the cascade model's lines (keV) and the Lu K_alpha1 X-ray that a pixel can lose
const double kAlKa = 354.9, kSiKa = 411.7, kAlLa = 65.73, kSiLa = 76.24, kLuKa1 = 54.07;
const double kBinKeV = 2., kMaxKeV = 2000.;

struct Sample {        // panel entries of one kind of muon, per 2 keV, with the time weight
  TH1D* all = nullptr;
  TH1D* xray = nullptr;   // entries with energy of a muonic X-ray of the target (source bit 1)
  TH1D* rest = nullptr;   // all others: decay electrons, capture gammas, neutrons, activation
  TH1D* late = nullptr;   // hits more than 1 ms after the muon (activation)
  double events = 0.;     // written events
  double entries = 0.;    // weighted entries above the threshold, all energies
  void Make(const std::string& name)
  {
    const int n = int(kMaxKeV / kBinKeV);
    all = new TH1D((name + "_all").c_str(), "", n, 0., kMaxKeV);
    xray = new TH1D((name + "_xray").c_str(), "", n, 0., kMaxKeV);
    rest = new TH1D((name + "_rest").c_str(), "", n, 0., kMaxKeV);
    late = new TH1D((name + "_late").c_str(), "", n, 0., kMaxKeV);
    for (TH1D* h : {all, xray, rest, late}) { h->Sumw2(); h->SetDirectory(nullptr); }
  }
};

double WindowSum(const TH1D* h, double lo, double hi, double* err = nullptr)
{
  double sum = 0., e2 = 0.;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    const double x = h->GetBinCenter(b);
    if (x < lo || x > hi) continue;
    sum += h->GetBinContent(b);
    e2 += h->GetBinError(b) * h->GetBinError(b);
  }
  if (err) *err = std::sqrt(e2);
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

std::string Num(double v)   // two significant digits, without an exponent from 0.001 to 10^7
{
  char s[64];
  if (v >= 100. && v < 1e7) snprintf(s, sizeof s, "%.0f", std::round(v / std::pow(10., std::floor(std::log10(v)) - 1.)) * std::pow(10., std::floor(std::log10(v)) - 1.));
  else if (v >= 1e-3 && v < 100.) snprintf(s, sizeof s, "%.2g", v);
  else snprintf(s, sizeof s, "%.1e", v);
  return s;
}
}  // namespace

void beam_impurity_spectrum(const char* dir, const char* outDir = "", double fSim = 0.01, int zImp = 13, int zHost = 14,
                            double rate = 1.32e8, double seconds = 86400., double sigma1MeV = 26.85,
                            const char* decayFile = "output/lyso20_concept/decay/decay_1.root", double luActivity = 12557.,
                            double threshold = 20.)
{
  gStyle->SetOptStat(0);
  const std::string base = dir;
  const std::string beamDir = base + "/beam";
  const std::string out = std::string(outDir).empty() ? base + "/analysis" : std::string(outDir);
  gSystem->mkdir(out.c_str(), true);
  std::ostringstream txt;

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

  Sample host, imp, other;
  host.Make("host"); imp.Make("imp"); other.Make("other");
  TRandom3 random(4357);
  const double tMax = seconds * 1e9;   // ns
  double nSim = 0.;
  long nHit = 0;
  for (const auto& [k, name] : files) {
    TFile f(name.c_str());
    auto run = dynamic_cast<TTree*>(f.Get("run"));
    auto ev = dynamic_cast<TTree*>(f.Get("events"));
    if (!run || !ev) { printf("skipped %s: no run or events tree\n", name.c_str()); continue; }
    run->Draw("n_events", "", "goff");
    for (Long64_t i = 0; i < run->GetSelectedRows(); ++i) nSim += run->GetV1()[i];
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
      Sample& s = capZ == zHost ? host : capZ == zImp ? imp : other;
      s.events += 1.;
      for (size_t h = 0; h < e->size(); ++h) {
        const double eMeV = (*e)[h], tNs = (*t)[h];
        if (eMeV <= 0. || tNs >= tMax) continue;
        const double w = 1. - std::max(0., tNs) / tMax;
        const double kev = 1000. * eMeV + random.Gaus(0., sigma1MeV * std::sqrt(eMeV));
        if (kev < threshold) continue;
        ++nHit;
        s.all->Fill(kev, w);
        ((*src)[h] & 1 ? s.xray : s.rest)->Fill(kev, w);
        if (tNs > 1e6) s.late->Fill(kev, w);
        s.entries += w;
      }
    }
  }
  if (host.events <= 0. || imp.events <= 0.) { printf("no events of Z %d or Z %d in %s\n", zHost, zImp, beamDir.c_str()); return; }

  // the panel's own Lu-176 over the time T
  TH1D* lu = new TH1D("lu", "", int(kMaxKeV / kBinKeV), 0., kMaxKeV);
  lu->Sumw2(); lu->SetDirectory(nullptr);
  double nDecays = 0., luEntries = 0.;
  {
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

  // weights
  const double nDay = rate * seconds;
  const double A = (imp.events / host.events) * (1. - fSim) / fSim;
  auto P = [A](double f) { return f * A / (1. - f + f * A); };
  const double pSim = P(fSim);
  auto wHost = [&](double c) { return nDay * (1. - P(c)) / (nSim * (1. - pSim)); };
  auto wImp = [&](double c) { return nDay * P(c) / (nSim * pSim); };
  const double wOther = nDay / nSim;
  auto Spectrum = [&](double c, const std::string& name) {
    auto h = dynamic_cast<TH1D*>(host.all->Clone(name.c_str()));
    h->SetDirectory(nullptr);
    h->Scale(wHost(c));
    h->Add(imp.all, wImp(c));
    h->Add(other.all, wOther);
    h->Add(lu);
    return h;
  };

  // the Al K_alpha window and the per-muon numbers in it
  const double half = 2. * sigma1MeV * std::sqrt(kAlKa / 1000.);
  const double lo = kAlKa - half, hi = kAlKa + half;
  double eHost = 0., eImp = 0.;
  const double perHost = WindowSum(host.all, lo, hi, &eHost) / (nSim * (1. - pSim));
  const double perImp = WindowSum(imp.all, lo, hi, &eImp) / (nSim * pSim);
  const double perDiff = perImp - perHost;
  const double perDiffErr = std::hypot(eImp / (nSim * pSim), eHost / (nSim * (1. - pSim)));
  TH1D* h0 = Spectrum(0., "si_alone");
  const double B = WindowSum(h0, lo, hi);
  const double bLu = WindowSum(lu, lo, hi);
  const double bLate = WindowSum(host.late, lo, hi) * wHost(0.);
  const double bOther = WindowSum(other.all, lo, hi) * wOther;
  auto S = [&](double c) { return nDay * P(c) * perDiff; };
  auto Zstat = [&](double c) { return S(c) / std::sqrt(B); };
  auto Zsys = [&](double c, double delta) { return S(c) / std::sqrt(B + delta * delta * B * B); };
  auto Limit = [&](double z, double delta) { return z * std::sqrt(B + delta * delta * B * B) / (nDay * A * perDiff); };

  // the pile-up the simulation leaves out: entries per second and pixel, all energies above the threshold
  const double entriesPerMuon = host.entries / (nSim * (1. - pSim));
  const double luPerSecond = nDecays > 0. ? luEntries / nDecays * luActivity : 0.;
  const double perPixel = (entriesPerMuon * rate + luPerSecond) / 400.;
  const double pileup = 1. - std::exp(-2. * perPixel * 100e-9);

  char line[512];
  txt << "Beam study of Al in Si: " << base << " (" << files.size() << " finished jobs), "
      << Form("%.6g", nSim) << " simulated mu-, " << nHit << " panel entries of " << threshold << " keV and up\n";
  txt << "  written events: Si-caught (Z " << zHost << ") " << host.events << ", Al-caught (Z " << zImp << ") " << imp.events
      << ", other " << other.events << "\n";
  snprintf(line, sizeof line, "  capture strength per Al atom against a Si atom A = %.4f (Z law 13/14 = %.4f); in the simulation %.3f%% of the muons were caught by Al (atom fraction %.3g)\n",
           A, 13. / 14., 100. * pSim, fSim);
  txt << line;
  snprintf(line, sizeof line, "one day: %.3g mu-/s x %.0f s = %.4g muons, %.4g times the simulated number\n", rate, seconds, nDay, nDay / nSim);
  txt << line;
  snprintf(line, sizeof line, "  LYSO entries per Si-caught muon (20 keV and up, time-weighted): %.4f; with Lu-176 %.4g entries/s in all, %.4g per pixel and second\n",
           entriesPerMuon, entriesPerMuon * rate + luPerSecond, perPixel);
  txt << line;
  snprintf(line, sizeof line, "  pile-up left out: a pixel's entry has another within +-100 ns with probability %.1f%%\n", 100. * pileup);
  txt << line;
  snprintf(line, sizeof line, "Al K_alpha window %.1f-%.1f keV (%.1f keV +-2 sigma, sigma %.2f keV; simulated line energies, measured Al K_alpha 346.83 keV)\n",
           lo, hi, kAlKa, half / 2.);
  txt << line;
  snprintf(line, sizeof line, "  entries in it per Si-caught muon %.5f, per Al-caught muon %.5f: an Al atom in place of a Si atom adds %.5f +- %.5f (MC)\n",
           perHost, perImp, perDiff, perDiffErr);
  txt << line;
  snprintf(line, sizeof line, "  Si alone, one day: B = %.4g entries (Lu-176 %.3g, %.2g%%; activation, later than 1 ms, %.3g, %.2g%%; other muons %.3g)\n",
           B, bLu, 100. * bLu / B, bLate, 100. * bLate / B, bOther);
  txt << line;
  {
    const double bX = WindowSum(host.xray, lo, hi) * wHost(0.), bR = WindowSum(host.rest, lo, hi) * wHost(0.);
    snprintf(line, sizeof line, "    of B: the Si-caught muons' muonic X-rays (partial deposits, escape peaks, K_alpha's tail) %.1f%%, their other entries (decay electrons, capture gammas, neutrons, activation) %.1f%%\n",
             100. * bX / B, 100. * bR / B);
    txt << line;
  }
  txt << "  concentration (atom fraction) | Al-caught muons/day | excess S | S/B | Z stat | Z with 0.1% | Z with 1% of B known\n";
  for (double c : {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2}) {
    snprintf(line, sizeof line, "  %8.0e (%6g ppm) | %10.3g | %10.3g | %9.2e | %8.3g | %8.3g | %8.3g\n", c, c * 1e6, nDay * P(c), S(c),
             S(c) / B, Zstat(c), Zsys(c, 1e-3), Zsys(c, 1e-2));
    txt << line;
  }
  snprintf(line, sizeof line, "  lowest concentration seen in one day: 3 sigma %.3g ppm, 5 sigma %.3g ppm (statistics alone); with the Si-alone shape known to 0.1%%: %.3g ppm, %.3g ppm; to 1%%: %.3g ppm, %.3g ppm\n",
           1e6 * Limit(3., 0.), 1e6 * Limit(5., 0.), 1e6 * Limit(3., 1e-3), 1e6 * Limit(5., 1e-3), 1e6 * Limit(3., 1e-2),
           1e6 * Limit(5., 1e-2));
  txt << line;
  snprintf(line, sizeof line, "  MC statistics of the excess per Al atom: %.1f%% (relative)\n", 100. * perDiffErr / std::fabs(perDiff));
  txt << line;

  // 1. the day's spectrum with 1 ppm Al, by what made the entries
  const double c1 = 1e-6;
  {
    TH1D* total = Spectrum(c1, "total_1ppm");
    auto hx = dynamic_cast<TH1D*>(host.xray->Clone("si_xray")); hx->Scale(wHost(c1));
    auto hr = dynamic_cast<TH1D*>(host.rest->Clone("si_rest")); hr->Scale(wHost(c1)); hr->Add(other.all, wOther);
    auto ha = dynamic_cast<TH1D*>(imp.all->Clone("al_1ppm")); ha->Scale(wImp(c1));
    auto c = new TCanvas("c_day", "", 1000, 700);
    c->SetLogy(); c->SetLeftMargin(0.11); c->SetRightMargin(0.03); c->SetTopMargin(0.08);
    total->SetTitle(Form("One day at %s #mu^{-}/s, 85 MeV/c: 4 cm Si with 1 ppm Al, LYSO 12.69 x 12.69 cm;"
                         "energy in one pixel (keV);entries per %g keV in one day", Sci(rate).c_str(), kBinKeV));
    total->GetXaxis()->SetRangeUser(20., 1000.);
    total->SetMinimum(1e2); total->SetMaximum(total->GetMaximum() * 3e3);   // room for the legend above the curves
    total->SetLineColor(kBlack); total->SetLineWidth(2);
    hx->SetLineColor(kBlue); hr->SetLineColor(kGreen + 2); ha->SetLineColor(kRed); lu->SetLineColor(kViolet);
    total->Draw("hist");
    for (TH1D* h : {hx, hr, ha, lu}) { h->SetLineWidth(2); h->Draw("hist same"); }
    auto leg = new TLegend(0.36, 0.71, 0.97, 0.91);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.03);
    leg->AddEntry(total, Form("all entries (%s in one day, 20-1000 keV)", Sci(total->Integral(11, 500)).c_str()), "l");
    leg->AddEntry(hx, "Si-caught muons: their muonic X-rays", "l");
    leg->AddEntry(hr, "Si-caught muons: decay e^{-}, #gamma, n, activation", "l");
    leg->AddEntry(ha, Form("Al-caught muons at 1 ppm (%s of them per day)", Sci(nDay * P(c1)).c_str()), "l");
    leg->AddEntry(lu, "the LYSO's own Lu-176", "l");
    leg->Draw();
    TLatex tex; tex.SetTextSize(0.026); tex.SetTextAngle(90); tex.SetTextAlign(12);
    for (auto [x, label] : std::vector<std::pair<double, std::string>>{{0.5 * (kAlLa + kSiLa), "L_{#alpha} (Al, Si)"},
                                                                       {kAlKa, "Al K_{#alpha}"}, {kSiKa, "Si K_{#alpha}"}}) {
      const double y = total->GetBinContent(total->FindBin(x)) * 2.;
      tex.DrawLatex(x, y, label.c_str());
    }
    Save(c, out + "/beam_spectrum_day");
  }

  // 2. entries per caught muon: Si against Al, where the lines are
  {
    auto hs = dynamic_cast<TH1D*>(host.all->Clone("per_si")); hs->Scale(1. / (nSim * (1. - pSim)));
    auto ha = dynamic_cast<TH1D*>(imp.all->Clone("per_al")); ha->Scale(1. / (nSim * pSim));
    auto c = new TCanvas("c_shapes", "", 1000, 700);
    c->SetLeftMargin(0.12); c->SetRightMargin(0.03); c->SetTopMargin(0.08);
    hs->SetTitle(Form("Entries per caught muon: 4 cm Si target, beam 85 MeV/c, LYSO 12.69 x 12.69 cm;energy in one pixel (keV);"
                      "entries per %g keV per caught muon", kBinKeV));
    hs->GetXaxis()->SetRangeUser(250., 500.);
    ha->GetXaxis()->SetRangeUser(250., 500.);   // the maxima within the drawn range
    hs->SetMaximum(1.3 * std::max(hs->GetMaximum(), ha->GetMaximum()));
    hs->SetMinimum(0.);
    hs->SetLineColor(kBlack); hs->SetLineWidth(2); ha->SetLineColor(kRed); ha->SetLineWidth(2);
    hs->GetYaxis()->SetTitleOffset(1.4);
    hs->Draw("hist"); ha->Draw("hist same");
    auto box = new TH1D("window", "", 1, lo, hi);
    box->SetBinContent(1, hs->GetMaximum()); box->SetFillColorAlpha(kOrange, 0.25); box->SetLineWidth(0);
    box->Draw("hist same"); hs->Draw("hist same"); ha->Draw("hist same");
    auto leg = new TLegend(0.66, 0.52, 0.97, 0.72);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
    leg->AddEntry(hs, "muon caught by Si", "l");
    leg->AddEntry(ha, "muon caught by Al", "l");
    leg->AddEntry(box, Form("Al K_{#alpha} window %.0f-%.0f keV", lo, hi), "f");
    leg->Draw();
    // the lines (the cascade model's energies), marked from the axis up to their labels
    const double ymax = hs->GetMaximum();
    TLatex tex; tex.SetTextSize(0.028); tex.SetTextAlign(21);
    TLine mark; mark.SetLineStyle(3); mark.SetLineColor(kGray + 2);
    const std::vector<std::tuple<double, double, std::string, int>> marks = {
      {kSiKa, 0.95, "Si K_{#alpha}", kBlack}, {kSiKa - kLuKa1, 0.88, "Si K_{#alpha} - Lu K_{#alpha} (escape)", kBlack},
      {kAlKa, 0.81, "Al K_{#alpha}", kRed}};
    for (const auto& [x, f, label, color] : marks) {
      mark.DrawLine(x, 0., x, (f - 0.01) * ymax);
      tex.SetTextColor(color);
      tex.DrawLatex(x, f * ymax, label.c_str());
    }
    Save(c, out + "/beam_shapes");
  }

  // 3. the excess over Si alone in units of the day's noise: a pseudo-measurement with 1 ppm, and the expected excess
  {
    TH1D* data = Spectrum(c1, "pseudo_1ppm");
    const int b1 = h0->FindBin(250.), b2 = h0->FindBin(500.);
    auto pulls = new TH1D("pulls", "", b2 - b1 + 1, h0->GetBinLowEdge(b1), h0->GetBinLowEdge(b2 + 1));
    std::vector<TH1D*> expect;
    const std::vector<double> cs = {1e-6, 1e-5, 1e-4};
    const std::vector<int> colors = {kRed, kOrange + 1, kBlue};
    for (size_t i = 0; i < cs.size(); ++i) {
      auto h = new TH1D(Form("expect_%zu", i), "", b2 - b1 + 1, h0->GetBinLowEdge(b1), h0->GetBinLowEdge(b2 + 1));
      TH1D* hc = Spectrum(cs[i], Form("spec_%zu", i));
      for (int b = b1; b <= b2; ++b) h->SetBinContent(b - b1 + 1, (hc->GetBinContent(b) - h0->GetBinContent(b)) / std::sqrt(h0->GetBinContent(b)));
      h->SetLineColor(colors[i]); h->SetLineWidth(2);
      expect.push_back(h);
    }
    TRandom3 noise(2026);
    for (int b = b1; b <= b2; ++b) {
      const double n = h0->GetBinContent(b);
      const double measured = data->GetBinContent(b) + noise.Gaus(0., std::sqrt(data->GetBinContent(b)));
      pulls->SetBinContent(b - b1 + 1, (measured - n) / std::sqrt(n));
      pulls->SetBinError(b - b1 + 1, 1.);
    }
    auto c = new TCanvas("c_residual", "", 1000, 700);
    c->SetLeftMargin(0.11); c->SetRightMargin(0.03); c->SetTopMargin(0.08);
    double ymax = 4.;
    for (TH1D* h : expect) ymax = std::max(ymax, 1.6 * h->GetMaximum());   // room for the legend
    double ymin = -4.;
    for (TH1D* h : expect) ymin = std::min(ymin, 1.15 * h->GetMinimum());
    pulls->SetTitle(Form("One day of %s #mu^{-}/s: (Si with Al) - (Si alone), in units of the day's noise;energy in one pixel (keV);"
                         "excess per %g keV / #sqrt{entries}", Sci(rate).c_str(), kBinKeV));
    pulls->SetMinimum(ymin); pulls->SetMaximum(ymax);
    pulls->SetMarkerStyle(20); pulls->SetMarkerSize(0.6); pulls->SetLineColor(kBlack);
    pulls->Draw("e1");
    TLine zero; zero.SetLineStyle(2); zero.DrawLine(250., 0., 500., 0.);
    for (TH1D* h : expect) h->Draw("hist same");
    pulls->Draw("e1 same");
    auto leg = new TLegend(0.14, 0.66, 0.60, 0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
    leg->AddEntry(pulls, "a measurement with 1 ppm Al (1#sigma error bars)", "pe");
    for (size_t i = 0; i < cs.size(); ++i)
      leg->AddEntry(expect[i], Form("expected excess, %g ppm (the Al K_{#alpha} window: %s#sigma)", cs[i] * 1e6, Num(Zstat(cs[i])).c_str()), "l");
    leg->Draw();
    Save(c, out + "/beam_residual");
  }

  // 4. the window's significance against the concentration
  {
    const int n = 61;
    auto gStat = new TGraph(n), gSys3 = new TGraph(n), gSys2 = new TGraph(n);
    for (int i = 0; i < n; ++i) {
      const double c = std::pow(10., -8. + 6. * i / (n - 1.));
      gStat->SetPoint(i, c * 1e6, Zstat(c));
      gSys3->SetPoint(i, c * 1e6, Zsys(c, 1e-3));
      gSys2->SetPoint(i, c * 1e6, Zsys(c, 1e-2));
    }
    auto c = new TCanvas("c_sig", "", 1000, 700);
    c->SetLogx(); c->SetLogy(); c->SetLeftMargin(0.11); c->SetRightMargin(0.03); c->SetTopMargin(0.08);
    auto frame = c->DrawFrame(1e-2, 1e-3, 1e4, 3e4);
    frame->SetTitle(Form("Al in the Si target, one day of %s #mu^{-}/s: the Al K_{#alpha} window (%.0f-%.0f keV);"
                         "Al concentration (ppm of the atoms);excess / noise (#sigma)", Sci(rate).c_str(), lo, hi));
    gStat->SetLineColor(kBlack); gSys3->SetLineColor(kBlue); gSys2->SetLineColor(kRed);
    for (TGraph* g : {gStat, gSys3, gSys2}) { g->SetLineWidth(3); g->Draw("l same"); }
    TLine l; l.SetLineStyle(2);
    l.DrawLine(1e-2, 3., 1e4, 3.); l.DrawLine(1e-2, 5., 1e4, 5.);
    l.SetLineColor(kGray + 2); l.DrawLine(1., 1e-3, 1., 3e4);
    TLatex tex; tex.SetTextSize(0.03);
    tex.DrawLatex(2e-2, 3.4, "3#sigma"); tex.DrawLatex(2e-2, 5.8, "5#sigma");
    tex.SetTextColor(kGray + 2); tex.DrawLatex(1.15, 3e-2, "1 ppm");
    auto leg = new TLegend(0.14, 0.70, 0.70, 0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
    leg->AddEntry(gStat, Form("statistics alone (3#sigma at %s ppm)", Num(1e6 * Limit(3., 0.)).c_str()), "l");
    leg->AddEntry(gSys3, Form("Si-alone spectrum known to 0.1%% (3#sigma at %s ppm)", Num(1e6 * Limit(3., 1e-3)).c_str()), "l");
    leg->AddEntry(gSys2, Form("Si-alone spectrum known to 1%% (3#sigma at %s ppm)", Num(1e6 * Limit(3., 1e-2)).c_str()), "l");
    leg->Draw();
    Save(c, out + "/beam_significance");
  }

  std::ofstream(out + "/beam_impurity_spectrum.txt") << txt.str();
  printf("%s", txt.str().c_str());
  printf("written into %s: beam_impurity_spectrum.txt, beam_spectrum_day, beam_shapes, beam_residual, beam_significance (.png, .pdf)\n",
         out.c_str());
}
