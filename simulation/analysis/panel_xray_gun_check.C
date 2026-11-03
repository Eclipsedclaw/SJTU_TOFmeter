// What the pixel panel itself does to an X-ray line (macro/panel_xray_gun.mac: photons of one energy into the panel from
// above, the target out), for the LYSO panel against the CZT plane: per photon shot, how often it leaves its full energy
// in one pixel (what the single-pixel analysis can use), shares it between pixels, leaves only part of it (a K X-ray of
// the crystal or a scattered photon escaped), or nothing (passed or scattered out without a deposit: not written). Then,
// with the panel's response (sigma = sigma1MeV x sqrt(E / MeV), single pixels, 10 keV threshold), the share of photons
// with an entry in the analysis window (line +- 2 sigma). Reads <dir>/gun_<config>_<tag>.root for the configurations
// lyso_fluo (LYSO, fluorescence on, as the LYSO runs), czt_fluo (CZT with it on) and czt_nofluo (CZT as the CZT runs,
// fluorescence off); those missing are skipped. The line's energy is the primary's (ekin). Writes
// <dir>/panel_xray_gun_<tag>.txt (also printed) and panel_xray_gun_<tag>.png/.pdf: the true single-pixel deposits per
// photon, the configurations overlaid (one comparison plot).
//   root -l -b -q 'analysis/panel_xray_gun_check.C("output/lyso_concept/gun", "al")'

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

void panel_xray_gun_check(const char* dir = "output/lyso_concept/gun", const char* tag = "al", double sigmaLyso = 26.85,
                          double sigmaCzt = 6.36991, long long shot = 100000)
{
  struct Config { std::string name, title; double sigma; int colour, style; };
  const std::vector<Config> configs = {{"lyso_fluo", "LYSO, fluorescence on (the LYSO runs)", sigmaLyso, kRed + 1, 1},
                                       {"czt_fluo", "CZT, fluorescence on", sigmaCzt, kBlue + 1, 1},
                                       {"czt_nofluo", "CZT, fluorescence off (the CZT runs)", sigmaCzt, kBlue + 1, 2}};
  gStyle->SetOptStat(0);
  std::ostringstream out;
  std::vector<TH1D*> spectra;
  std::vector<std::string> titles;
  double line = 0.;
  TRandom3 random(6573);
  for (const auto& c : configs) {
    const std::string path = Form("%s/gun_%s_%s.root", dir, c.name.c_str(), tag);
    if (gSystem->AccessPathName(path.c_str())) continue;
    TFile file(path.c_str());
    auto events = dynamic_cast<TTree*>(file.Get("events"));
    auto run = dynamic_cast<TTree*>(file.Get("run"));
    if (!events || !run) continue;
    int nEvents = 0;
    run->SetBranchAddress("n_events", &nEvents);
    run->GetEntry(0);
    double ekin = 0.;
    std::vector<double>* e = nullptr;
    events->SetBranchAddress("ekin", &ekin);
    events->SetBranchAddress("czt_e", &e);
    auto h = new TH1D(Form("h_%s", c.name.c_str()), "", 100, 0., 100.);
    h->SetDirectory(nullptr);
    long long one = 0, shared = 0, part = 0, inWindow = 0;
    for (Long64_t i = 0; i < events->GetEntries(); ++i) {
      events->GetEntry(i);
      line = ekin * 1000.;
      const double sigma = c.sigma * std::sqrt(line / 1000.);
      double sum = 0., largest = 0.;
      bool entry = false;
      for (double d : *e) {
        const double keV = d * 1000.;
        sum += keV;
        largest = std::max(largest, keV);
        h->Fill(keV);
        const double measured = keV + random.Gaus(0., c.sigma * std::sqrt(d));
        entry = entry || (measured >= 10. && std::abs(measured - line) <= 2. * sigma);
      }
      if (largest >= line - 0.5) ++one;
      else if (sum >= line - 0.5) ++shared;
      else ++part;
      if (entry) ++inWindow;
    }
    const double n = nEvents > 0 ? nEvents : shot;
    h->Scale(1. / n);
    h->SetLineColor(c.colour);
    h->SetLineStyle(c.style);
    h->SetLineWidth(2);
    spectra.push_back(h);
    titles.push_back(c.title);
    out << Form("%-40s %.2f keV, %.0f photons: full energy in one pixel %5.1f%%, shared by pixels %4.1f%%, part of it %5.1f%%, "
                "none %5.1f%%; an entry within +-2 sigma (%.2f keV) %5.1f%%\n",
                c.title.c_str(), line, n, 100. * one / n, 100. * shared / n, 100. * part / n,
                100. * (n - events->GetEntries()) / n, 2. * c.sigma * std::sqrt(line / 1000.), 100. * inWindow / n);
  }
  if (spectra.empty()) {
    printf("no gun files %s/gun_*_%s.root\n", dir, tag);
    return;
  }
  printf("%s", out.str().c_str());
  std::ofstream(Form("%s/panel_xray_gun_%s.txt", dir, tag)) << out.str();

  auto canvas = new TCanvas("gun", "gun", 1000, 700);
  canvas->SetLogy();
  canvas->SetLeftMargin(0.11);
  canvas->SetRightMargin(0.04);
  double top = 0.;
  for (auto h : spectra) top = std::max(top, h->GetMaximum());
  auto frame = spectra[0];
  frame->SetTitle(Form(";true energy in a pixel (keV);pixel deposits per photon shot, per keV"));
  frame->SetMaximum(top * 8.);
  frame->SetMinimum(1e-5);
  frame->Draw("hist");
  for (std::size_t i = 1; i < spectra.size(); ++i) spectra[i]->Draw("hist same");
  auto legend = new TLegend(0.13, 0.72, 0.62, 0.89);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  for (std::size_t i = 0; i < spectra.size(); ++i) legend->AddEntry(spectra[i], titles[i].c_str(), "l");
  legend->Draw();
  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextSize(0.033);
  text.DrawLatex(0.13, 0.66, Form("%.2f keV photons into the panel from above, no target; 1 keV bins", line));
  canvas->Print(Form("%s/panel_xray_gun_%s.png", dir, tag));
  canvas->Print(Form("%s/panel_xray_gun_%s.pdf", dir, tag));
}
