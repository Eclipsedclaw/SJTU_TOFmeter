// Thicker Al targets, the MC truth before any detector response or cut: how many muonic X-rays the mu- stopped in the
// target make, how many get out of it and how many reach the Si/YSO stack, per day, and what the response and the cuts
// leave of them. From the TOFmeter files of the 3-week run (1 cm target) and of scripts/muonic_xray_target.sh (2-5 cm):
//   root -l -b -q 'analysis/muonic_xray_target_truth.C'
// For every primary mu- captured at rest in the target (end_vol Target, end_proc muMinusCaptureAtRest):
// 1. its muonic X-rays (mux_src 0, mux_parent 1) at the Al lines: emitted (energy at creation), leaving the target
//    (mux_exit >= 0) and leaving it unscattered (the exit energy is the energy at creation);
// 2. what reaches the stack, true energies: the stack hits of the event grouped as analysis/muonic_xray_stack_skim.C
//    does (within 100 ns of the first hit), the entries with a target X-ray in them (stk_src bit 1); full energy: an
//    entry within 0.5 keV of a line, nothing else added to it;
// 3. the same entries with the stack response of the skim (sigma = 0.04 MeV x sqrt(E / MeV) per hit, hits below 10 keV
//    dropped; other random numbers than the skims): those in the +-2 sigma windows, i.e. S before the cuts, and how many
//    of them hold the line's full energy; after the cuts, S of muonic_xray_target_compare.C (cuts of the 1 cm target)
//    and of muonic_xray_target_cuts.C (each thickness's own cuts), from their CSVs if they are there.
// The runs wrote only the events with energy in the stack (/EventFilter 6): a stop that put nothing at all in the stack
// is missing from the stops and from 1 (1 cm target: 149 stops per day written, 159 per day from
// analysis/lead_target_stops.C); 2 and 3 are complete. analysis/muonic_xray.C with the prefix stk_ gives each
// thickness's whole truth spectrum and its line table.
// 4. Fits of the Al L lines in the spectrum with the response, 30-160 keV, per thickness: a smooth background exp(pol2)
//    and the L_alpha and L_beta lines as Gaussians, what the response makes of a line fully absorbed (L_alpha amount,
//    position and width free; L_beta amount free, 23 keV above L_alpha, width scaled as sqrt(E)); Poisson likelihood
//    on the counts, Baker-Cousins chi2. Landau lines (the energy loss of a charged particle crossing a thin layer) were
//    tried on 2026-10-30 with the same 7 free parameters and fitted worse at every thickness (1 cm: chi2 94/58 against
//    82/58), as the physics expects.
// Writes <dir>/muonic_xray_target_truth.csv and one plot per file, <dir>/muonic_xray_target_truth_<name>.png/.pdf:
// lines (made in the target, leaving it unscattered), stack (into the stack with the full energy), spectrum (the
// stack's truth spectrum, 0-200 keV), response (the same with the response, the Gaussian fits dashed), window (what the
// Al L_alpha window holds), funnel (Al L_alpha from emission to the window). The spectra's legends give each
// histogram's entries, as a ROOT statistics box does.
// Another energy resolution: sigma1MeV, the energy sigma at 1 MeV in keV (default 40, the skim's; scaling as sqrt(E)),
// for the smearing, the windows and the fits, with tag added to every file name; e.g. a CZT's FWHM/E of 1.5% at 1 MeV,
// sigma = 15 keV / 2.3548 = 6.37 keV:
//   root -l -b -q 'analysis/muonic_xray_target_truth.C("output/muonic_xray_stack", "output/muonic_xray_target", "2 3 4 5", 4357, 6.37, "_fwhm15")'
// The cut studies' S after the cuts are only shown for the default resolution, the one they were made with.

#include <algorithm>
#include <cstring>
#include <fstream>
#include <functional>
#include <sstream>

namespace {
// the stack response of analysis/muonic_xray_stack_skim.C and the windows of analysis/muonic_xray_stack_cuts.C: the
// energy sigma at 1 MeV (keV), set by the macro's argument
double gSigma1MeV = 40.;
double SigmaKeV(double eKeV) { return gSigma1MeV * std::sqrt(eKeV / 1000.); }

struct Window { const char* name; double eKeV; };
const std::vector<Window> kWindows = {{"Al L_alpha", 65.73}, {"Al L_beta", 88.73}, {"Al K_alpha", 354.9}};

bool InWindow(double eKeV, const Window& w) { return std::abs(eKeV - w.eKeV) <= 2. * SigmaKeV(w.eKeV); }

struct Sample {
  double cm = 1., days = 0.;
  TString dir;
  double stops = 0., xrays = 0.;   // mu- captured in the target; their muonic X-rays above 10 keV
  // per line: emitted; left the target; left it unscattered; stack entries of its full energy (truth); stack entries
  // with a target X-ray in its window with the response, and those of them with the line's full energy
  double emitted[3] = {}, left[3] = {}, unscattered[3] = {}, full[3] = {}, window[3] = {}, windowFull[3] = {};
  double entries = 0., responseEntries = 0.;   // stack entries with a target X-ray, 20-1000 keV: truth, with the response
  double before[3] = {-1., -1., -1.}, cut1cm[3] = {-1., -1., -1.}, cutOwn[3] = {-1., -1., -1.};   // per day, from the CSVs
  TH1D *truth = nullptr, *response = nullptr;
};
}  // namespace

void muonic_xray_target_truth(const char* ref = "output/muonic_xray_stack", const char* dir = "output/muonic_xray_target",
                              const char* thicknesses = "2 3 4 5", unsigned seed = 4357, double sigma1MeV = 40.,
                              const char* tag = "")
{
  TH1::AddDirectory(false);
  gSigma1MeV = sigma1MeV;
  std::vector<Sample> samples(1);
  samples[0].dir = ref;
  std::istringstream list(thicknesses);
  for (double d; list >> d;) {
    Sample s;
    s.cm = d;
    s.dir = Form("%s/t%gcm", dir, d);
    samples.push_back(s);
  }

  TRandom3 random(seed);
  struct Hit { double t, e, measured; int src; };   // ns, keV (true, with the response)
  std::vector<Hit> hits;
  for (auto& s : samples) {
    s.truth = new TH1D(Form("truth%g", s.cm), "", 500, 0., 500.);         // 1 keV
    s.response = new TH1D(Form("response%g", s.cm), "", 250, 0., 500.);   // 2 keV
    for (int k = 1; !gSystem->AccessPathName(Form("%s/stk_%d.root", s.dir.Data(), k)); ++k) {
      if (gSystem->AccessPathName(Form("%s/log/job%d.done", s.dir.Data(), k))) continue;   // not finished
      TFile file(Form("%s/stk_%d.root", s.dir.Data(), k));
      auto run = dynamic_cast<TTree*>(file.Get("run"));
      auto events = dynamic_cast<TTree*>(file.Get("events"));
      if (!run || !events || run->GetEntries() != 1) continue;
      double lt = 0.;
      run->SetBranchAddress("live_time", &lt);
      run->GetEntry(0);
      s.days += lt / 86400.;
      int pdg = 0;
      char endProc[128] = "", endVol[128] = "";
      std::vector<double>*muxE = nullptr, *muxExit = nullptr, *stkE = nullptr, *stkT = nullptr;
      std::vector<int>*muxSrc = nullptr, *muxParent = nullptr, *stkSrc = nullptr;
      events->SetBranchAddress("pdg", &pdg);
      events->SetBranchAddress("end_proc", endProc);
      events->SetBranchAddress("end_vol", endVol);
      events->SetBranchAddress("mux_e", &muxE); events->SetBranchAddress("mux_exit", &muxExit);
      events->SetBranchAddress("mux_src", &muxSrc); events->SetBranchAddress("mux_parent", &muxParent);
      events->SetBranchAddress("stk_e", &stkE); events->SetBranchAddress("stk_t", &stkT);
      events->SetBranchAddress("stk_src", &stkSrc);
      // read the rest of an event only for a mu- captured in the target
      TBranch *bPdg = events->GetBranch("pdg"), *bProc = events->GetBranch("end_proc"), *bVol = events->GetBranch("end_vol");
      std::vector<TBranch*> rest;
      for (const char* b : {"mux_e", "mux_exit", "mux_src", "mux_parent", "stk_e", "stk_t", "stk_src"}) rest.push_back(events->GetBranch(b));
      for (Long64_t i = 0; i < events->GetEntries(); ++i) {
        bPdg->GetEntry(i);
        if (pdg != 13) continue;
        bProc->GetEntry(i);
        if (std::strcmp(endProc, "muMinusCaptureAtRest") != 0) continue;
        bVol->GetEntry(i);
        if (std::strncmp(endVol, "Target", 6) != 0) continue;
        for (auto b : rest) b->GetEntry(i);
        s.stops += 1.;
        // 1. its muonic X-rays
        for (std::size_t p = 0; p < muxE->size(); ++p) {
          if ((*muxParent)[p] != 1 || (*muxSrc)[p] != 0) continue;
          const double eKeV = 1000. * (*muxE)[p], exitKeV = 1000. * (*muxExit)[p];
          if (eKeV > 10.) s.xrays += 1.;
          for (int w = 0; w < 3; ++w) {
            if (std::abs(eKeV - kWindows[w].eKeV) > 0.5) continue;
            s.emitted[w] += 1.;
            if (exitKeV < 0.) continue;   // absorbed in the target
            s.left[w] += 1.;
            if (std::abs(exitKeV - eKeV) < 1e-6) s.unscattered[w] += 1.;
          }
        }
        // 2. and 3. its stack entries, true and with the response
        hits.clear();
        for (std::size_t h = 0; h < stkE->size(); ++h) {
          const double e = (*stkE)[h];   // MeV
          hits.push_back({(*stkT)[h], 1000. * e, 1000. * (e + random.Gaus(0., gSigma1MeV / 1000. * std::sqrt(e))), (*stkSrc)[h]});
        }
        std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) { return a.t < b.t; });
        for (std::size_t h = 0; h < hits.size();) {
          const double t0 = hits[h].t;
          double eTrue = 0., eMeasured = 0.;
          int srcTrue = 0, srcMeasured = 0;
          for (; h < hits.size() && hits[h].t < t0 + 100.; ++h) {
            eTrue += hits[h].e;
            srcTrue |= hits[h].src;
            if (hits[h].measured < 10.) continue;   // the threshold of the skim
            eMeasured += hits[h].measured;
            srcMeasured |= hits[h].src;
          }
          if (srcTrue & 1) {
            s.truth->Fill(eTrue);
            if (eTrue >= 20. && eTrue <= 1000.) s.entries += 1.;
            for (int w = 0; w < 3; ++w) if (std::abs(eTrue - kWindows[w].eKeV) < 0.5) s.full[w] += 1.;
          }
          if (!(srcMeasured & 1)) continue;
          s.response->Fill(eMeasured);                           // every entry, down to the 10 keV threshold
          if (eMeasured < 20. || eMeasured > 1000.) continue;   // the candidates of the skim
          s.responseEntries += 1.;
          for (int w = 0; w < 3; ++w) {
            if (!InWindow(eMeasured, kWindows[w])) continue;
            s.window[w] += 1.;
            if (std::abs(eTrue - kWindows[w].eKeV) < 0.5) s.windowFull[w] += 1.;
          }
        }
      }
    }
    printf("%g cm: %s, %.2f days, %.0f mu- captured in the target\n", s.cm, s.dir.Data(), s.days, s.stops);
  }
  samples.erase(std::remove_if(samples.begin(), samples.end(), [](const Sample& s) { return s.days <= 0.; }), samples.end());
  if (samples.empty()) {
    printf("no finished stk_<k>.root in %s or %s/t<d>cm\n", ref, dir);
    return;
  }
  const int nS = samples.size();

  // S per day before and after the cuts from the CSVs of the cut studies: compare (before, 1 cm cuts), cuts (own cuts)
  auto readCsv = [&](const char* path, int column, int windowColumn, int cutsColumn) {
    std::ifstream in(path);
    std::string line;
    if (!std::getline(in, line)) return false;   // header
    while (std::getline(in, line)) {
      std::vector<std::string> f;
      std::stringstream fields(line);
      for (std::string x; std::getline(fields, x, ',');) f.push_back(x);
      if (int(f.size()) <= std::max({column, windowColumn, cutsColumn, 6})) continue;
      if (cutsColumn >= 0 && f[cutsColumn] != "own") continue;
      for (auto& s : samples) {
        if (std::abs(std::stod(f[0]) - s.cm) > 1e-6) continue;
        for (int w = 0; w < 3; ++w) {
          if (f[windowColumn] != kWindows[w].name) continue;
          if (cutsColumn >= 0) s.cutOwn[w] = std::stod(f[column]);
          else if (column == 5) s.before[w] = std::stod(f[column]);
          else s.cut1cm[w] = std::stod(f[column]);
        }
      }
    }
    return true;
  };
  const bool skimResponse = std::abs(sigma1MeV - 40.) < 1e-9;   // the cut studies were made with the skim's response
  const bool haveCompare = skimResponse && readCsv(Form("%s/muonic_xray_target_compare.csv", dir), 5, 4, -1) &&
                           readCsv(Form("%s/muonic_xray_target_compare.csv", dir), 6, 4, -1);
  const bool haveCuts = skimResponse && readCsv(Form("%s/muonic_xray_target_cuts.csv", dir), 6, 5, 2);

  // ---- tables
  printf("\n=== Thicker Al targets, MC truth per day (no detector response, no cut) ===\n");
  printf("response: energy sigma %.1f keV at 1 MeV, scaling as sqrt(E) (L_alpha: sigma %.2f keV)\n", gSigma1MeV, SigmaKeV(kWindows[0].eKeV));
  printf("mu- captured in the target (in the written events) and their muonic X-rays above 10 keV:\n");
  printf(" target    days     mu- stops per day    X-rays per day   per stop\n");
  for (const auto& s : samples)
    printf("  %4.1f   %6.2f    %7.1f +- %5.1f      %7.1f         %4.2f\n", s.cm, s.days, s.stops / s.days,
           std::sqrt(s.stops) / s.days, s.xrays / s.days, s.xrays / std::max(1., s.stops));
  std::ofstream csv(Form("%s/muonic_xray_target_truth%s.csv", dir, tag));
  csv << "target_cm,days,mu_stops_day,line,keV,emitted_day,left_day,left_unscattered_day,full_energy_stack_day,"
         "window_response_day,window_full_energy_day,window_skim_day,after_1cm_cuts_day,after_own_cuts_day\n";
  for (int w = 0; w < 3; ++w) {
    printf("\n%s (%.2f keV), per day; in brackets the share of the line's emitted X-rays:\n", kWindows[w].name, kWindows[w].eKeV);
    printf(" target   emitted   left the target   left unscattered   full energy     window (+-%.0f keV)     after the cuts\n",
           2. * SigmaKeV(kWindows[w].eKeV));
    printf("   cm                                                     in the stack    with the response      1 cm / own\n");
    for (const auto& s : samples) {
      const double d = s.days, n = std::max(1., s.emitted[w]);
      printf("  %4.1f    %6.1f    %6.1f (%3.0f%%)    %6.1f (%3.0f%%)     %5.1f (%4.1f%%)   %5.1f (full %4.1f)   %5.1f / %5.1f\n", s.cm,
             s.emitted[w] / d, s.left[w] / d, 100. * s.left[w] / n, s.unscattered[w] / d, 100. * s.unscattered[w] / n,
             s.full[w] / d, 100. * s.full[w] / n, s.window[w] / d, s.windowFull[w] / d, s.cut1cm[w], s.cutOwn[w]);
      csv << s.cm << "," << s.days << "," << s.stops / d << "," << kWindows[w].name << "," << kWindows[w].eKeV << ","
          << s.emitted[w] / d << "," << s.left[w] / d << "," << s.unscattered[w] / d << "," << s.full[w] / d << ","
          << s.window[w] / d << "," << s.windowFull[w] / d << "," << s.before[w] << "," << s.cut1cm[w] << "," << s.cutOwn[w]
          << "\n";
    }
  }
  printf("\nstack entries with a target X-ray, 20-1000 keV, per day: truth / with the response; window entries of the\n"
         "skims (S before the cuts of muonic_xray_target_compare.C) for the same response with other random numbers:\n");
  for (const auto& s : samples) {
    printf("  %4.1f cm: %6.1f / %6.1f;", s.cm, s.entries / s.days, s.responseEntries / s.days);
    for (int w = 0; w < 3; ++w) printf("   %s %5.1f (skim %5.1f)", kWindows[w].name, s.window[w] / s.days, s.before[w]);
    printf("\n");
  }
  if (!skimResponse) printf("(another response than the cut studies': no S after their cuts)\n");
  else if (!haveCompare || !haveCuts) printf("(no %s of the cut studies: no S after the cuts)\n", haveCompare ? "cuts CSV" : "CSVs");

  // ---- plots, one per file (CLAUDE.md rule 8)
  gStyle->SetOptStat(0);
  const int kFitStyle = 11;
  gStyle->SetLineStyleString(kFitStyle, "40 20");   // the legend's sample of the dashed fits: in the PNG it looks as drawDashed's
  const int colour[7] = {kGray + 2, kAzure + 2, kGreen + 2, kOrange + 1, kRed + 1, kMagenta + 1, kCyan + 2};   // as the compare
  const int lineColour[3] = {kRed + 1, kOrange + 1, kAzure + 2};
  const char* lineNames[3] = {"Al L_{#alpha}", "Al L_{#beta}", "Al K_{#alpha}"};
  auto canvas = [](const char* name, bool logy) {
    auto c = new TCanvas(name, name, 900, 700);
    c->SetLogy(logy);
    c->SetLeftMargin(0.12); c->SetRightMargin(0.03); c->SetTopMargin(0.09); c->SetBottomMargin(0.12);
    return c;
  };
  auto save = [&](TCanvas* c, const char* name, const char* title) {
    c->cd();
    TLatex heading;
    heading.SetNDC(); heading.SetTextSize(0.045);
    heading.DrawLatex(0.12, 0.935, title);
    c->Print(Form("%s/muonic_xray_target_truth%s_%s.pdf", dir, tag, name));
    c->SaveAs(Form("%s/muonic_xray_target_truth%s_%s.png", dir, tag, name));
    delete c;
  };
  auto frame = [](double x1, double y1, double x2, double y2, const char* titles) {
    auto f = gPad->DrawFrame(x1, y1, x2, y2, titles);
    f->GetXaxis()->SetTitleSize(0.045); f->GetYaxis()->SetTitleSize(0.045);
    f->GetXaxis()->SetTitleOffset(1.1); f->GetYaxis()->SetTitleOffset(1.25);
    return f;
  };
  auto legendAt = [](double x1, double y1, double x2, double y2, double size) {
    auto l = new TLegend(x1, y1, x2, y2);
    l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(size);
    return l;
  };
  // a quantity per day against the thickness
  auto perDay = [&](std::function<double(const Sample&)> count, int col, int marker, int style) {
    auto g = new TGraphErrors();
    for (const auto& s : samples) {
      const double n = count(s);
      const int p = g->GetN();
      g->SetPoint(p, s.cm, n / s.days);
      g->SetPointError(p, 0., std::sqrt(std::max(0., n)) / s.days);
    }
    g->SetMarkerStyle(marker); g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineStyle(style); g->SetLineWidth(2);
    g->SetMarkerSize(1.2);
    g->Draw("PL SAME");
    return g;
  };
  const double xMax = samples.back().cm + 0.5;
  double yMax = 0., eMax = 0.;
  for (const auto& s : samples) {
    yMax = std::max(yMax, s.xrays / s.days);
    eMax = std::max(eMax, s.entries / s.days);
  }

  // the Al lines made in the target and leaving it unscattered
  {
    auto c = canvas("lines", true);
    frame(0.5, 5., xMax, 15. * yMax, ";Al target thickness (cm);per day");
    auto legend = legendAt(0.14, 0.69, 0.97, 0.89, 0.034);
    legend->SetNColumns(2);
    legend->AddEntry(perDay([](const Sample& s) { return s.xrays; }, kGray + 2, 21, 1), "all muonic X-rays > 10 keV", "lp");
    legend->AddEntry(perDay([](const Sample& s) { return s.stops; }, kBlack, 24, 1), "#mu^{-} stopped in the target", "lp");
    for (int w : {2, 0, 1}) {
      legend->AddEntry(perDay([w](const Sample& s) { return s.emitted[w]; }, lineColour[w], 24, 2), Form("%s emitted", lineNames[w]), "lp");
      legend->AddEntry(perDay([w](const Sample& s) { return s.unscattered[w]; }, lineColour[w], 20, 1),
                       Form("%s leaving unscattered", lineNames[w]), "lp");
    }
    legend->Draw();
    save(c, "lines", "Made in the target, and leaving it unscattered");
  }

  // into the stack with the full energy
  {
    auto c = canvas("stack", true);
    frame(0.5, 0.5, xMax, 30. * eMax, ";Al target thickness (cm);stack entries per day (true energies)");
    auto legend = legendAt(0.14, 0.72, 0.97, 0.89, 0.034);
    legend->SetNColumns(2);
    legend->AddEntry(perDay([](const Sample& s) { return s.entries; }, kGray + 2, 21, 1), "with a target X-ray, 20-1000 keV", "lp");
    for (int w : {0, 1, 2})
      legend->AddEntry(perDay([w](const Sample& s) { return s.full[w]; }, lineColour[w], 20, 1), Form("%s, its full energy", lineNames[w]), "lp");
    legend->Draw();
    save(c, "stack", "Target X-rays into the stack (truth)");
  }

  // 4. fits of the Al L lines with the response, on the counts (before the plots scale them per day)
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");   // classic Minuit, as muonic_xray_stack_cuts.C
  const double fitLow = 30., fitHigh = 160., shiftB = kWindows[1].eKeV - kWindows[0].eKeV;
  const double widthB = std::sqrt(kWindows[1].eKeV / kWindows[0].eKeV);   // L_beta width / L_alpha width
  auto gaussian = [](double e, double pos, double width) {   // unit area
    return std::exp(-0.5 * std::pow((e - pos) / width, 2)) / (std::sqrt(2. * M_PI) * width);
  };
  std::vector<TF1*> fits(nS);
  std::vector<double> chi2(nS, 0.);
  std::vector<int> ndf(nS, 0);
  for (int k = 0; k < nS; ++k) {
    TH1D* h = samples[k].response;
    const double bin = h->GetBinWidth(1);
    auto f = new TF1(Form("fit%d", k), [=](double* x, double* p) {
      const double e = x[0];
      return std::exp(p[0] + p[1] * (e - 80.) + p[2] * (e - 80.) * (e - 80.)) +
             bin * (p[3] * gaussian(e, p[4], p[5]) + p[6] * gaussian(e, p[4] + shiftB, p[5] * widthB));
    }, fitLow, fitHigh, 7);
    const double nAlpha = samples[k].windowFull[0];
    f->SetParameters(std::log(std::max(1., 0.5 * h->GetBinContent(h->FindBin(80.)))), -0.01, 0., nAlpha, kWindows[0].eKeV,
                     SigmaKeV(kWindows[0].eKeV), 0.15 * nAlpha);
    f->SetParLimits(3, 0., 1e7); f->SetParLimits(4, 50., 80.); f->SetParLimits(5, 1., 30.); f->SetParLimits(6, 0., 1e7);
    h->Fit(f, "LQNR");
    // Baker-Cousins chi2 of a Poisson likelihood fit, 2 sum (mu - n + n ln(n / mu)), at the bin centres as the fit
    int bins = 0;
    for (int b = h->FindBin(fitLow + 1e-6); b <= h->FindBin(fitHigh - 1e-6); ++b, ++bins) {
      const double n = h->GetBinContent(b), mu = f->Eval(h->GetBinCenter(b));
      chi2[k] += 2. * (mu - n + (n > 0. ? n * std::log(n / mu) : 0.));
    }
    fits[k] = f;
    ndf[k] = bins - f->GetNumberFreeParameters();
  }
  printf("\nGaussian fits of the Al L lines with the response, %.0f-%.0f keV (Poisson likelihood; Baker-Cousins chi2):\n",
         fitLow, fitHigh);
  printf(" target   L_alpha per day   position (keV)   width (keV; response %4.1f)   chi2/ndf   truth: L_alpha full energy "
         "per day\n", SigmaKeV(kWindows[0].eKeV));
  for (int k = 0; k < nS; ++k) {
    const TF1* g = fits[k];
    const double d = samples[k].days;
    printf("  %4.1f      %6.1f +- %4.1f     %5.1f +- %3.1f          %5.1f +- %3.1f              %5.1f/%d          %5.1f\n",
           samples[k].cm, g->GetParameter(3) / d, g->GetParError(3) / d, g->GetParameter(4), g->GetParError(4),
           g->GetParameter(5), g->GetParError(5), chi2[k], ndf[k], samples[k].full[0] / d);
  }
  // a fit per day
  auto perDayCurve = [&](const TF1* f, double days, int col, const char* name) {
    auto curve = new TF1(name, [=](double* x, double*) { return f->Eval(x[0]) / days; }, fitLow, fitHigh, 0);
    curve->SetLineColor(col); curve->SetLineWidth(3);
    return curve;
  };
  // a curve drawn as separate dashes, dashPx long with gaps of gapPx, measured along the curve as it shows on the pad:
  // the PNG output starts a line style's pattern anew at every point of a curve, so a dashed TF1 comes out solid (or
  // as dense dots)
  auto drawDashed = [](const TF1* curve, double dashPx, double gapPx) {
    const double xScale = gPad->GetWw() * gPad->GetAbsWNDC() / (gPad->GetX2() - gPad->GetX1());   // pixels per unit
    const double yScale = gPad->GetWh() * gPad->GetAbsHNDC() / (gPad->GetY2() - gPad->GetY1());
    auto yPad = [](double y) { return gPad->GetLogy() ? std::log10(std::max(y, 1e-30)) : y; };
    TGraph* dash = nullptr;
    double length = 0., uLast = curve->GetXmin(), yLast = curve->Eval(uLast);
    for (double u = curve->GetXmin(); u <= curve->GetXmax() + 1e-9; u += 0.1) {
      const double y = curve->Eval(u);
      length += std::hypot((u - uLast) * xScale, (yPad(y) - yPad(yLast)) * yScale);
      uLast = u; yLast = y;
      if (std::fmod(length, dashPx + gapPx) < dashPx) {
        if (!dash) {
          dash = new TGraph();
          dash->SetLineColor(curve->GetLineColor()); dash->SetLineWidth(curve->GetLineWidth());
        }
        dash->SetPoint(dash->GetN(), u, y);
      } else if (dash) {
        dash->Draw("L SAME");
        dash = nullptr;
      }
    }
    if (dash) dash->Draw("L SAME");
  };

  // the stack spectra of the target's X-rays, true and with the response
  // the range xLow-xHigh; the legend from legendX on, each histogram with its entries (under- and overflow included,
  // as TH1::GetEntries and a statistics box count them)
  auto spectra = [&](const char* name, bool logy, TH1D* (*get)(Sample&), double xLow, double xHigh, double legendX,
                     const char* titles, const char* title, const std::vector<TF1*>& curves) {
    auto c = canvas(name, logy);
    double hMax = 0., hMin = 1e30;
    for (auto& s : samples) {
      TH1D* h = get(s);
      h->Scale(1. / s.days);
      for (int b = h->FindBin(xLow + 1e-6); b <= h->FindBin(xHigh - 1e-6); ++b) {
        hMax = std::max(hMax, h->GetBinContent(b));
        if (h->GetBinContent(b) > 0.) hMin = std::min(hMin, h->GetBinContent(b));
      }
    }
    const double yLow = logy ? 0.5 * hMin : 0., yHigh = logy ? 40. * hMax : 1.3 * hMax;
    frame(xLow, yLow, xHigh, yHigh, titles);
    TLatex mark;
    mark.SetTextSize(0.03); mark.SetTextColor(kGray + 2); mark.SetTextAngle(90); mark.SetTextAlign(12);
    for (int w = 0; w < 3; ++w) {   // where the lines are
      if (kWindows[w].eKeV > xHigh) continue;
      auto l = new TLine(kWindows[w].eKeV, yLow, kWindows[w].eKeV, logy ? 1.5 * hMax : 1.05 * hMax);
      l->SetLineColor(kGray + 1); l->SetLineStyle(3);
      l->Draw();
      mark.DrawLatex(kWindows[w].eKeV + 3., logy ? 2. * hMax : 1.08 * hMax, lineNames[w]);
    }
    auto legend = legendAt(legendX, 0.6, legendX + 0.41, 0.89, 0.034);
    for (int k = nS - 1; k >= 0; --k) {   // the thickest first, the 1 cm target on top
      TH1D* h = get(samples[k]);
      h->SetLineColor(colour[k % 7]); h->SetLineWidth(2);
      h->Draw("HIST SAME");
    }
    for (auto curve : curves) drawDashed(curve, 12., 6.);
    for (int k = 0; k < nS; ++k)
      legend->AddEntry(get(samples[k]), Form("%g cm target: %.0f entries", samples[k].cm, get(samples[k])->GetEntries()), "l");
    if (!curves.empty()) {
      auto dashed = new TLine();   // the legend's sample of the fits, in the colours of their spectra
      dashed->SetLineStyle(kFitStyle); dashed->SetLineWidth(3);
      legend->AddEntry(dashed, "fits with Gaussian lines", "l");
    }
    legend->Draw();
    save(c, name, title);
  };
  spectra("spectrum", true, [](Sample& s) { return s.truth; }, 0., 200., 0.56,   // the legend over 130-200 keV
          ";true add-back energy (keV);entries with a target X-ray, per day per keV",
          "The stack's MC truth spectrum of the target's X-rays", {});
  std::vector<TF1*> gaussCurves;
  for (int k = 0; k < nS; ++k)
    gaussCurves.push_back(perDayCurve(fits[k], samples[k].days, colour[k % 7], Form("gaussDay%d", k)));
  spectra("response", false, [](Sample& s) { return s.response; }, 0., 200., 0.56,   // the legend over 130-200 keV
          ";add-back energy with the response (keV);entries with a target X-ray, per day per 2 keV",
          Form("Stack response, FWHM/E %.1f%% at 1 MeV, before any cut", 2.3548 * gSigma1MeV / 10.), gaussCurves);

  // what the Al L_alpha window holds with the response: the line's full energy and the other target X-rays
  {
    std::vector<double> edges;
    for (const auto& s : samples) { edges.push_back(s.cm - 0.3); edges.push_back(s.cm + 0.3); }
    auto hFull = new TH1D("hFull", "", int(edges.size()) - 1, edges.data());
    auto hAll = new TH1D("hAll", "", int(edges.size()) - 1, edges.data());
    auto gRef = new TGraph(), gOwn = new TGraph();
    double top = 0.;
    for (const auto& s : samples) {
      const int b = hFull->FindBin(s.cm);
      hFull->SetBinContent(b, s.windowFull[0] / s.days);
      hAll->SetBinContent(b, s.window[0] / s.days);
      if (s.cut1cm[0] > 0.) gRef->SetPoint(gRef->GetN(), s.cm - 0.1, s.cut1cm[0]);
      if (s.cutOwn[0] > 0.) gOwn->SetPoint(gOwn->GetN(), s.cm + 0.1, s.cutOwn[0]);
      top = std::max({top, s.window[0] / s.days, s.cut1cm[0], s.cutOwn[0]});
    }
    auto c = canvas("window", false);
    frame(0.5, 0., xMax, 1.6 * top, ";Al target thickness (cm);Al L-alpha window entries per day");
    hAll->SetFillColor(kRed - 9); hAll->SetLineColor(kRed - 9);
    hFull->SetFillColor(kRed + 1); hFull->SetLineColor(kRed + 1);
    hAll->Draw("HIST SAME");
    hFull->Draw("HIST SAME");
    gRef->SetMarkerStyle(25); gRef->SetMarkerSize(1.6); gRef->SetMarkerColor(kBlack);
    gOwn->SetMarkerStyle(20); gOwn->SetMarkerSize(1.6); gOwn->SetMarkerColor(kBlack);
    auto legend = legendAt(0.14, 0.68, 0.97, 0.89, 0.034);
    legend->AddEntry(hFull, "the line's full energy (65.73 keV, truth)", "f");
    legend->AddEntry(hAll, "other target X-rays: scattered, part of their energy, summed", "f");
    if (gRef->GetN() > 0) {   // the cut studies' S, with the response they were made with
      gRef->Draw("P SAME");
      gOwn->Draw("P SAME");
      legend->AddEntry(gRef, "S after the cuts of the 1 cm target (every target X-ray entry)", "p");
      legend->AddEntry(gOwn, "S after each thickness's own cuts (every target X-ray entry)", "p");
    }
    legend->Draw();
    gPad->RedrawAxis();
    save(c, "window", "What the Al L-alpha window holds, before any cut");
  }

  // Al L_alpha from emission to the window
  {
    auto c = canvas("funnel", true);
    c->SetBottomMargin(0.17);
    double fMax = 0., fMin = 1e30;
    for (const auto& s : samples) {
      fMax = std::max(fMax, s.emitted[0] / s.days);
      fMin = std::min(fMin, s.windowFull[0] / s.days);
    }
    auto f = frame(0.5, 0.4 * fMin, 4.5, 4. * fMax, ";;Al L-alpha per day");
    f->GetXaxis()->SetLabelSize(0.);   // the stages are written under the axis
    f->GetXaxis()->SetTickLength(0.);
    const char* stages[4] = {"emitted", "#splitline{left the target}{unscattered}", "#splitline{full energy}{in the stack}",
                             "#splitline{full energy in the}{window (response)}"};
    TLatex label;
    label.SetNDC(); label.SetTextAlign(23); label.SetTextSize(0.036); label.SetTextFont(42);
    for (int b = 0; b < 4; ++b) label.DrawLatex(0.12 + (b + 0.5) / 4. * 0.85, 0.155, stages[b]);
    auto legend = legendAt(0.7, 0.62, 0.97, 0.89, 0.036);
    for (int k = 0; k < nS; ++k) {
      const auto& s = samples[k];
      const double values[4] = {s.emitted[0], s.unscattered[0], s.full[0], s.windowFull[0]};
      auto g = new TGraphErrors();
      for (int b = 0; b < 4; ++b) {
        g->SetPoint(b, b + 1. + 0.04 * (k - 0.5 * (nS - 1)), values[b] / s.days);
        g->SetPointError(b, 0., std::sqrt(values[b]) / s.days);
      }
      g->SetMarkerStyle(20); g->SetMarkerSize(1.2); g->SetMarkerColor(colour[k % 7]); g->SetLineColor(colour[k % 7]);
      g->SetLineWidth(2);
      g->Draw("PL SAME");
      legend->AddEntry(g, Form("%g cm target", s.cm), "lp");
    }
    legend->Draw();
    save(c, "funnel", "Al L-alpha from emission to the stack (truth)");
  }
  printf("\nplots: %s/muonic_xray_target_truth%s_{lines,stack,spectrum,response,window,funnel}.png (and .pdf); table: "
         "%s/muonic_xray_target_truth%s.csv\n", dir, tag, dir, tag);
}
