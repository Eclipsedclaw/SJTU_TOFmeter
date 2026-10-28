// Selection of the muonic X-rays of the target in the Si/YSO stack, from the skim of analysis/muonic_xray_stack_skim.C:
//   root -l -b -q 'analysis/muonic_xray_stack_cuts.C("output/muonic_xray_stack/xray_skim.root")'
// 1. Stopped muons (truth): the add-back entries within -5..+50 ns of the capture of the primary mu-, by the volume it
//    stopped in: where muonic X-rays come from at all, and how often.
// 2. Cut variables: shapes of signal (entries with an X-ray of the target, stk_src bit 1) and background, 40-420 keV,
//    for entries within -5..+30 ns of a middle-TOF hit.
// 3. Cuts, one after the other, with S, B and S/sqrt(B) in +-2 sigma around Al L_alpha, L_beta and K_alpha:
//    middle-TOF hit in coincidence; no bottom-TOF hit within -5..+50 ns of it; the muon in a middle bar over the target
//    (2-6); a slow muon: more than eMidCut in the middle-TOF bar; the highest stack hit in the top two layers; if the
//    top TOF was hit, a top-to-middle time above tofCut. Assumed TOF response (the simulation has true values): time
//    sigma sigmaT per station, energy sigma/E sigmaE, the same smearing for every entry of an event.
// Arguments: skim, TOF time (0.2 ns) and energy resolution (15%), the slow-muon (3 MeV: best L_alpha S/sqrt(B) of a scan
// of 1-3.5 MeV) and time-of-flight (1.3 ns) cuts. Finally a fit of the counts after all cuts (no truth): smooth
// background and the Al L_alpha and L_beta lines. Writes, next to the skim, muonic_xray_stack_stopped.png/.pdf,
// _cutvars.png/.pdf, _cuts.png/.pdf, _fit.png/.pdf and _cutflow.csv.

#include <set>

namespace {
// the stack response as in the skim: sigma (keV) of an energy (keV)
double SigmaKeV(double eKeV) { return 40. * std::sqrt(eKeV / 1000.); }

struct Window { const char* name; double eKeV; };
const std::vector<Window> kWindows = {{"Al L_alpha", 65.73}, {"Al L_beta", 88.73}, {"Al K_alpha", 354.9}};

bool InWindow(double eKeV, const Window& w) { return std::abs(eKeV - w.eKeV) <= 2. * SigmaKeV(w.eKeV); }
}  // namespace

void muonic_xray_stack_cuts(const char* skim = "output/muonic_xray_stack/xray_skim.root", double sigmaT = 0.2,
                            double sigmaE = 0.15, double eMidCut = 3.0, double tofCut = 1.3)
{
  TFile file(skim);
  auto cand = dynamic_cast<TTree*>(file.Get("cand"));
  auto info = dynamic_cast<TTree*>(file.Get("info"));
  if (!cand || !info) {
    printf("no skim in %s (run analysis/muonic_xray_stack_skim.C)\n", skim);
    return;
  }
  double lt = 0., liveTime = 0.;
  info->SetBranchAddress("live_time", &lt);
  for (Long64_t i = 0; i < info->GetEntries(); ++i) { info->GetEntry(i); liveTime += lt; }
  const double days = liveTime / 86400., perDay = 1. / days;
  int proc, event, n, layers, topLayer, src, pdg, endVol, endProc, tofTop, tofMid, barTop, barMid;
  float e, emax, rmax, t, tTop, tMid, eTop, eMid, botFirst, ekin, endT, endZ;
  cand->SetBranchAddress("proc", &proc); cand->SetBranchAddress("event", &event);
  cand->SetBranchAddress("e", &e); cand->SetBranchAddress("n", &n); cand->SetBranchAddress("topLayer", &topLayer);
  cand->SetBranchAddress("rmax", &rmax); cand->SetBranchAddress("t", &t); cand->SetBranchAddress("src", &src);
  cand->SetBranchAddress("pdg", &pdg); cand->SetBranchAddress("endVol", &endVol); cand->SetBranchAddress("endProc", &endProc);
  cand->SetBranchAddress("endT", &endT);
  cand->SetBranchAddress("tofTop", &tofTop); cand->SetBranchAddress("tofMid", &tofMid);
  cand->SetBranchAddress("tTop", &tTop); cand->SetBranchAddress("tMid", &tMid);
  cand->SetBranchAddress("barMid", &barMid); cand->SetBranchAddress("eMid", &eMid); cand->SetBranchAddress("botFirst", &botFirst);

  // ---- histograms
  const int nBins = 250;   // 2 keV, 0-500 keV
  // 1. stopped mu- by volume: 0 target, 1 Al plates, 2 stack, 3 TOF bars, 4 lead, 5 steel and the rest
  const char* stopNames[] = {"in the Al target", "in the Al plates", "in the stack (Si, YSO)", "in the TOF bars",
                             "in the lead", "elsewhere (steel, air)"};
  auto stopGroup = [](int volume) { return volume == 0 ? 0 : volume == 1 ? 1 : (volume == 2 || volume == 3) ? 2 : volume == 4 ? 3
                                                   : volume == 5 ? 4 : 5; };
  TH1D* hStop[6];
  long stopCount[6] = {0, 0, 0, 0, 0, 0};
  std::set<long long> stopEvents[6];
  for (int g = 0; g < 6; ++g) hStop[g] = new TH1D(Form("hStop%d", g), "", nBins, 0., 500.);
  // 2. cut variables (40-420 keV, prompt with a middle-TOF hit), signal and background
  auto vEMid = new TH1D("vEMidS", "", 60, 0., 12.), vEMidB = new TH1D("vEMidB", "", 60, 0., 12.);
  auto vBar = new TH1D("vBarS", "", 9, -0.5, 8.5), vBarB = new TH1D("vBarB", "", 9, -0.5, 8.5);
  auto vBot = new TH1D("vBotS", "", 2, -0.5, 1.5), vBotB = new TH1D("vBotB", "", 2, -0.5, 1.5);
  auto vTof = new TH1D("vTofS", "", 60, 0., 6.), vTofB = new TH1D("vTofB", "", 60, 0., 6.);
  auto vTop = new TH1D("vTopS", "", 5, -0.5, 4.5), vTopB = new TH1D("vTopB", "", 5, -0.5, 4.5);
  auto vR = new TH1D("vRS", "", 20, 0., 10.), vRB = new TH1D("vRB", "", 20, 0., 10.);
  // 3. stages: 0 middle TOF in coincidence, 1 + bottom veto, 2 + bar over the target, 3 + slow muon,
  //    4 + top layers, 5 + time of flight where measured
  const int nStages = 6;
  // std::string: Form() reuses a few buffers, so its pointers would be overwritten by the next calls
  const std::vector<std::string> stageNames = {"middle-TOF hit, entry within -5..+30 ns", "+ no bottom-TOF hit within 50 ns",
                                               "+ muon in middle bars 2-6 (over the target)",
                                               Form("+ slow muon: middle TOF > %.1f MeV", eMidCut),
                                               "+ highest hit in the top two layers",
                                               Form("+ top-to-middle time > %.1f ns (if top hit)", tofCut)};
  TH1D *hAll[nStages], *hSig[nStages];
  double sW[nStages][3] = {}, bW[nStages][3] = {};
  for (int s = 0; s < nStages; ++s) {
    hAll[s] = new TH1D(Form("hAll%d", s), "", nBins, 0., 500.);
    hSig[s] = new TH1D(Form("hSig%d", s), "", nBins, 0., 500.);
  }
  // after all cuts: counts in 30-160 keV for the fit (no truth), and what the background of the L_alpha window is
  auto hFit = new TH1D("hFit", "", 65, 30., 160.);
  hFit->Sumw2();   // keep sqrt(N) of the counts, also after the background is subtracted from a copy
  std::map<std::string, double> composition;
  // slow-muon threshold scan after stages 0-2: S and B in the L_alpha window
  const std::vector<double> scan = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5};
  std::vector<double> scanS(scan.size(), 0.), scanB(scan.size(), 0.);

  TRandom3 random;
  for (Long64_t i = 0; i < cand->GetEntries(); ++i) {
    cand->GetEntry(i);
    const bool signal = src & 1;
    // 1. stopped mu- (truth), prompt entries
    if (endProc == 1 && t - endT > -5. && t - endT < 50.) {
      const int g = stopGroup(endVol);
      hStop[g]->Fill(e);
      stopEvents[g].insert(1000000000LL * proc + event);
    }
    if (!tofMid || t - tMid < -5. || t - tMid > 30.) continue;
    // the assumed TOF response, the same for every entry of an event
    random.SetSeed(1000003ULL * proc + event + 1);
    const double eMidMeasured = eMid * (1. + sigmaE * random.Gaus());
    const double dtTof = tofTop ? (tMid + sigmaT * random.Gaus()) - (tTop + sigmaT * random.Gaus()) : -1.;
    const bool topBefore = tofTop && tTop < tMid + 5.;   // the top hit belongs to the muon (not a later particle)
    const bool bottom = botFirst >= 0. && botFirst - tMid > -5. && botFirst - tMid < 50.;
    const bool overTarget = barMid >= 2 && barMid <= 6;
    const bool slow = eMidMeasured > eMidCut;
    const bool upper = topLayer >= 3;
    const bool tofSlow = !topBefore || dtTof > tofCut;
    // 2. cut variables
    if (e >= 40. && e <= 420.) {
      (signal ? vEMid : vEMidB)->Fill(eMidMeasured);
      (signal ? vBar : vBarB)->Fill(barMid);
      (signal ? vBot : vBotB)->Fill(bottom ? 1 : 0);
      if (topBefore) (signal ? vTof : vTofB)->Fill(dtTof);
      (signal ? vTop : vTopB)->Fill(topLayer);
      (signal ? vR : vRB)->Fill(rmax);
    }
    // 3. stages
    const bool pass[nStages] = {true, !bottom, overTarget, slow, upper, tofSlow};
    if (!bottom && overTarget)
      for (std::size_t k = 0; k < scan.size(); ++k)
        if (eMidMeasured > scan[k] && InWindow(e, kWindows[0])) (signal ? scanS : scanB)[k] += 1.;
    bool all = true;
    for (int s = 0; s < nStages; ++s) {
      if (!pass[s]) { all = false; break; }
      hAll[s]->Fill(e);
      if (signal) hSig[s]->Fill(e);
      for (int w = 0; w < 3; ++w)
        if (InWindow(e, kWindows[w])) (signal ? sW : bW)[s][w] += 1.;
    }
    if (all) {
      hFit->Fill(e);
      if (!signal && InWindow(e, kWindows[0])) {
        // muonic X-rays of another volume (bit 2), by where the primary stopped; else what the primary did
        const char* volumes[] = {"target", "Al plates", "stack Si", "stack YSO", "TOF bars", "lead", "steel", "air"};
        if (src & 2) composition[std::string("muonic X-rays, muon stopped in the ") + volumes[endVol]] += 1.;
        else if (src & 4) composition["nuclear-capture gammas"] += 1.;
        else if (endProc == 1) composition[std::string("other energy, mu- stopped in the ") + volumes[endVol]] += 1.;
        else if (endProc == 2) composition[std::string("other energy, mu+ stopped in the ") + volumes[endVol]] += 1.;
        else composition["other energy, muon passed through"] += 1.;
      }
    }
  }
  for (int g = 0; g < 6; ++g) stopCount[g] = stopEvents[g].size();

  // ---- tables
  TString base = skim;
  base.ReplaceAll("xray_skim.root", "muonic_xray_stack");
  printf("\n=== X-ray selection in the Si/YSO stack (%s): %.2f days of exposure ===\n", skim, days);
  printf("assumed TOF response: time sigma %.2f ns per station, energy sigma/E %.0f%%\n", sigmaT, 100. * sigmaE);
  printf("\n1. mu- stopped (truth), with a prompt add-back entry of 20-1000 keV in the stack:\n");
  for (int g = 0; g < 6; ++g) {
    const int b1 = hStop[g]->FindBin(kWindows[0].eKeV - 2. * SigmaKeV(kWindows[0].eKeV));
    const int b2 = hStop[g]->FindBin(kWindows[0].eKeV + 2. * SigmaKeV(kWindows[0].eKeV));
    printf("   %-24s %7ld stops (%7.1f per day), %8.0f entries, in the L_alpha window %6.0f\n", stopNames[g], stopCount[g],
           stopCount[g] * perDay, hStop[g]->GetEntries(), hStop[g]->Integral(b1, b2));
  }
  printf("\n3. cut flow (counts of the %.1f days; S = entries with an X-ray of the target):\n", days);
  printf("%-50s", "");
  for (const auto& w : kWindows) printf(" %24s", Form("%s %.0f keV", w.name, w.eKeV));
  printf("\n%-50s", "");
  for (int w = 0; w < 3; ++w) printf(" %8s %8s %6s", "S", "B", "S/vB");
  printf("\n");
  std::ofstream csv((base + "_cutflow.csv").Data());
  csv << "stage,cut";
  for (const auto& w : kWindows) csv << "," << w.name << " S," << w.name << " B," << w.name << " S/sqrt(B)";
  csv << "\n";
  for (int s = 0; s < nStages; ++s) {
    printf("%d %-48s", s + 1, stageNames[s].c_str());
    csv << s + 1 << ",\"" << stageNames[s] << "\"";
    for (int w = 0; w < 3; ++w) {
      const double z = bW[s][w] > 0. ? sW[s][w] / std::sqrt(bW[s][w]) : 0.;
      printf(" %8.0f %8.0f %6.2f", sW[s][w], bW[s][w], z);
      csv << "," << sW[s][w] << "," << bW[s][w] << "," << z;
    }
    printf("\n");
    csv << "\n";
  }
  printf("\nbackground of the L_alpha window after all cuts (%.0f entries):\n", bW[nStages - 1][0]);
  std::vector<std::pair<double, std::string>> parts;
  for (const auto& [name, count] : composition) parts.push_back({count, name});
  std::sort(parts.rbegin(), parts.rend());
  for (const auto& [count, name] : parts)
    printf("   %-50s %6.0f (%4.1f%%)\n", name.c_str(), count, 100. * count / bW[nStages - 1][0]);

  // ---- the lines found without truth: smooth background exp(pol2) and two Gaussians of the known energy and width,
  // log-likelihood fit to the counts of 30-160 keV after all cuts
  const double binWidth = hFit->GetBinWidth(1);
  auto lines = [&](double* x, double* p) {
    const double gA = std::exp(-0.5 * std::pow((x[0] - p[4]) / SigmaKeV(kWindows[0].eKeV), 2)) /
                      (std::sqrt(2. * M_PI) * SigmaKeV(kWindows[0].eKeV));
    const double gB = std::exp(-0.5 * std::pow((x[0] - kWindows[1].eKeV) / SigmaKeV(kWindows[1].eKeV), 2)) /
                      (std::sqrt(2. * M_PI) * SigmaKeV(kWindows[1].eKeV));
    return std::exp(p[0] + p[1] * (x[0] - 80.) + p[2] * std::pow(x[0] - 80., 2)) + binWidth * (p[3] * gA + p[5] * gB);
  };
  auto model = new TF1("model", lines, 30., 160., 6);
  model->SetParameters(std::log(std::max(1., hFit->GetBinContent(hFit->FindBin(80.)))), -0.005, 0., 100., kWindows[0].eKeV, 50.);
  model->SetParLimits(4, kWindows[0].eKeV - 10., kWindows[0].eKeV + 10.);
  model->SetNpx(400);
  // classic Minuit: Minuit2 of the ROOT on the shared server starts MPI (its hardware probe and PMIx then complain)
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");
  hFit->Fit(model, "LQ0");
  const double nAlpha = model->GetParameter(3), dAlpha = model->GetParError(3);
  const double nBeta = model->GetParameter(5), dBeta = model->GetParError(5);
  printf("\nfit of the counts after all cuts (30-160 keV, no truth): L_alpha %.0f +- %.0f at %.1f keV (%.1f sigma), L_beta %.0f +- %.0f "
         "(%.1f sigma)\n", nAlpha, dAlpha, model->GetParameter(4), dAlpha > 0. ? nAlpha / dAlpha : 0., nBeta, dBeta,
         dBeta > 0. ? nBeta / dBeta : 0.);
  // the same as a likelihood ratio, the usual discovery significance: the fit again without the L_alpha line,
  // Z = sqrt(2 x the log-likelihood it loses); N: the counts keep the first fit
  TF1 noAlpha("noAlpha", lines, 30., 160., 6);
  noAlpha.SetParameters(model->GetParameters());
  noAlpha.FixParameter(3, 0.);
  noAlpha.FixParameter(4, model->GetParameter(4));
  hFit->Fit(&noAlpha, "LQ0N");
  auto nll = [&](const TF1& f) {   // Poisson -log L of the bins, at the bin centres as the fit takes it
    double sum = 0.;
    for (int b = 1; b <= hFit->GetNbinsX(); ++b) {
      const double mu = f.Eval(hFit->GetBinCenter(b));
      sum += mu - hFit->GetBinContent(b) * std::log(mu);
    }
    return sum;
  };
  printf("likelihood ratio against the fit without L_alpha: %.1f sigma\n",
         std::sqrt(std::max(0., 2. * (nll(noAlpha) - nll(*model)))));

  printf("\nslow-muon threshold after stages 1-3, L_alpha window:\n");
  for (std::size_t k = 0; k < scan.size(); ++k)
    printf("   middle TOF > %.1f MeV: S %5.0f  B %7.0f  S/sqrt(B) %.2f\n", scan[k], scanS[k], scanB[k],
           scanB[k] > 0. ? scanS[k] / std::sqrt(scanB[k]) : 0.);

  // ---- plots
  const int blue = TColor::GetColor("#2a78d6"), orange = TColor::GetColor("#eb6834"), aqua = TColor::GetColor("#1baf7a");
  const int yellow = TColor::GetColor("#eda100"), magenta = TColor::GetColor("#e87ba4"), green = TColor::GetColor("#008300");
  const int grey = TColor::GetColor("#8c8a83"), inkSecondary = TColor::GetColor("#52514e");
  const int colours[6] = {orange, blue, aqua, yellow, magenta, grey};
  gStyle->SetOptStat(0);
  auto axes = [&](TH1* h, const char* titles, double size) {
    h->SetTitle(titles);
    for (auto axis : {h->GetXaxis(), h->GetYaxis()}) {
      axis->SetTitleSize(size);
      axis->SetLabelSize(0.9 * size);
      axis->SetTitleColor(inkSecondary);
    }
  };
  auto lineMarkers = [&](double yLow, double yHigh, double textY, double textSize, bool withEnergy = true) {
    auto label = new TLatex();
    label->SetTextFont(42);
    label->SetTextSize(textSize);
    label->SetTextColor(inkSecondary);
    label->SetTextAngle(90);
    label->SetTextAlign(12);
    const char* latex[] = {"Al L_{#alpha}", "Al L_{#beta}", "Al K_{#alpha}"};
    for (int w = 0; w < 3; ++w) {
      auto line = new TLine(kWindows[w].eKeV, yLow, kWindows[w].eKeV, yHigh);
      line->SetLineStyle(2);
      line->SetLineColor(inkSecondary);
      line->Draw();
      label->DrawLatex(kWindows[w].eKeV + 3., textY, withEnergy ? Form("%s %.1f", latex[w], kWindows[w].eKeV) : latex[w]);
    }
  };
  auto legendAt = [&](double x1, double y1, double x2, double y2, double size) {
    auto legend = new TLegend(x1, y1, x2, y2);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(size);
    return legend;
  };
  auto title = [&](const char* text) {
    auto latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(62);
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.08, 0.935, text);
  };

  // 1. stopped muons
  {
    TCanvas c("c1", "stopped muons", 1400, 800);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.08); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.1); gPad->SetBottomMargin(0.11);
    auto all = static_cast<TH1D*>(hStop[0]->Clone("allStops"));
    for (int g = 1; g < 6; ++g) all->Add(hStop[g]);
    for (auto h : hStop) h->Scale(perDay);
    all->Scale(perDay);
    axes(all, ";add-back energy (keV);entries per day per 2 keV", 0.04);
    double minimum = 1e30;
    for (int b = 1; b <= nBins; ++b)
      for (auto h : hStop) if (h->GetBinContent(b) > 0.) minimum = std::min(minimum, h->GetBinContent(b));
    const double maximum = all->GetMaximum();
    all->SetMinimum(0.5 * minimum);
    all->SetMaximum(3e4 * maximum);   // the legend above the line labels
    all->SetLineColor(inkSecondary);
    all->SetLineWidth(3);
    all->Draw("HIST");
    for (int g = 0; g < 6; ++g) {
      hStop[g]->SetLineColor(colours[g]);
      hStop[g]->SetLineWidth(2);
      hStop[g]->Draw("HIST SAME");
    }
    lineMarkers(0.5 * minimum, 2. * maximum, 2.5 * maximum, 0.032, false);
    auto legend = legendAt(0.3, 0.66, 0.97, 0.88, 0.03);
    legend->SetNColumns(2);
    // the Greek letter first: ROOT 6.32 leaves a gap before one that follows other text
    legend->AddEntry(all, Form("#mu^{-} stopped anywhere (%.0f/day)", (stopCount[0] + stopCount[1] + stopCount[2] + stopCount[3] + stopCount[4] + stopCount[5]) * perDay), "l");
    for (int g = 0; g < 6; ++g) legend->AddEntry(hStop[g], Form("%s (%.0f/day)", stopNames[g], stopCount[g] * perDay), "l");
    legend->Draw();
    title(Form("Prompt stack entries of stopped muons (truth), %.1f days at sea level", days));
    c.Print(base + "_stopped.pdf");
    c.SaveAs(base + "_stopped.png");
  }

  // 2. cut variables, unit area
  {
    TCanvas c("c2", "cut variables", 1400, 850);
    auto pads = new TPad("pads", "", 0., 0., 1., 0.93);
    pads->Draw();
    pads->Divide(3, 2, 0.01, 0.01);
    struct Panel { TH1D *s, *b; const char* titles; double cut; bool logy; };
    const Panel panels[] = {{vEMid, vEMidB, ";middle-TOF bar energy (MeV);fraction", eMidCut, true},
                            {vTof, vTofB, ";top-to-middle time (ns);fraction", tofCut, true},
                            {vBar, vBarB, ";middle-TOF bar;fraction", -1., false},
                            {vBot, vBotB, ";bottom-TOF hit within 50 ns (0 no, 1 yes);fraction", -1., false},
                            {vTop, vTopB, ";highest stack layer hit (4 = top Si);fraction", -1., false},
                            {vR, vRB, ";largest |x| or |y| of the hits (cm);fraction", -1., false}};
    for (int p = 0; p < 6; ++p) {
      pads->cd(p + 1);
      gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.04); gPad->SetRightMargin(0.03);
      if (panels[p].logy) gPad->SetLogy();
      for (auto h : {panels[p].s, panels[p].b}) if (h->Integral() > 0.) h->Scale(1. / h->Integral());
      axes(panels[p].b, panels[p].titles, 0.055);
      panels[p].b->SetLineColor(grey); panels[p].b->SetLineWidth(2);
      panels[p].s->SetLineColor(orange); panels[p].s->SetLineWidth(2);
      panels[p].b->SetMaximum((panels[p].logy ? 5. : 1.25) * std::max(panels[p].s->GetMaximum(), panels[p].b->GetMaximum()));
      if (panels[p].logy) panels[p].b->SetMinimum(1e-4);
      else panels[p].b->SetMinimum(0.);
      panels[p].b->Draw("HIST");
      panels[p].s->Draw("HIST SAME");
      if (panels[p].cut > 0.) {
        auto line = new TLine(panels[p].cut, panels[p].logy ? 1e-4 : 0., panels[p].cut, panels[p].b->GetMaximum());
        line->SetLineStyle(2); line->SetLineColor(inkSecondary); line->Draw();
      }
      if (p == 0) {
        auto legend = legendAt(0.45, 0.78, 0.97, 0.94, 0.05);
        legend->AddEntry(panels[p].s, "X-ray of the target", "l");
        legend->AddEntry(panels[p].b, "everything else", "l");
        legend->Draw();
      }
    }
    c.cd();
    title(Form("Cut variables, 40-420 keV entries within -5..+30 ns of a middle-TOF hit (TOF: %.1f ns, %.0f%%)", sigmaT, 100. * sigmaE));
    c.Print(base + "_cutvars.pdf");
    c.SaveAs(base + "_cutvars.png");
  }

  // 3. spectra after the first and the last stage
  {
    TCanvas c("c3", "after the cuts", 1400, 850);
    auto pPlot = new TPad("pPlot", "", 0., 0.16, 1., 1.);
    auto pText = new TPad("pText", "", 0., 0., 1., 0.16);
    pPlot->Draw(); pText->Draw();
    pPlot->cd();
    gPad->SetLogy();
    gPad->SetLeftMargin(0.08); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.1); gPad->SetBottomMargin(0.11);
    const int last = nStages - 1;
    for (auto h : {hAll[0], hSig[0], hAll[last], hSig[last]}) h->Scale(perDay);
    axes(hAll[0], ";add-back energy (keV);entries per day per 2 keV", 0.042);
    double minimum = 1e30;
    for (int b = 1; b <= nBins; ++b)
      for (auto h : {hSig[0], hSig[last], hAll[last]}) if (h->GetBinContent(b) > 0.) minimum = std::min(minimum, h->GetBinContent(b));
    const double maximum = hAll[0]->GetMaximum();
    hAll[0]->SetMinimum(0.5 * minimum);
    hAll[0]->SetMaximum(300. * maximum);
    hAll[0]->SetLineColor(grey); hAll[0]->SetLineWidth(2); hAll[0]->SetLineStyle(2);
    hSig[0]->SetLineColor(orange); hSig[0]->SetLineWidth(1); hSig[0]->SetLineStyle(2);
    hAll[last]->SetLineColor(blue); hAll[last]->SetLineWidth(3);
    hSig[last]->SetLineColor(orange); hSig[last]->SetLineWidth(2);
    hAll[0]->Draw("HIST");
    hSig[0]->Draw("HIST SAME");
    hAll[last]->Draw("HIST SAME");
    hSig[last]->Draw("HIST SAME");
    lineMarkers(0.5 * minimum, 2. * maximum, 2.5 * maximum, 0.035);
    auto legend = legendAt(0.3, 0.68, 0.7, 0.88, 0.032);
    legend->AddEntry(hAll[0], "middle-TOF coincidence only: all entries", "l");
    legend->AddEntry(hSig[0], "   ... with an X-ray of the target", "l");
    legend->AddEntry(hAll[last], "after all cuts: all entries", "l");
    legend->AddEntry(hSig[last], "   ... with an X-ray of the target", "l");
    legend->Draw();
    title(Form("Stack add-back spectrum after the X-ray selection, %.1f days at sea level", days));
    pText->cd();
    auto box = new TPaveText(0.01, 0.05, 0.99, 0.95, "NDC");
    box->SetFillColor(kWhite); box->SetBorderSize(1); box->SetTextAlign(12); box->SetTextFont(42);
    box->SetTextSize(0.14); box->SetTextColor(inkSecondary);
    box->AddText(Form("Cuts: middle TOF in coincidence, no bottom TOF, middle bars 2-6, middle TOF > %.1f MeV, highest hit in the top two "
                      "layers, top-to-middle > %.1f ns if measured.", eMidCut, tofCut));
    std::string summary;
    for (int w = 0; w < 3; ++w)
      summary += Form("%s%s: S %.0f, B %.0f, S/#sqrt{B} %.1f", w ? ";  " : "", kWindows[w].name, sW[last][w], bW[last][w],
                      bW[last][w] > 0. ? sW[last][w] / std::sqrt(bW[last][w]) : 0.);
    box->AddText(("After all cuts, +-2 sigma, counts of the exposure: " + summary).c_str());
    box->AddText(Form("Stack: sigma/E 4%% at 1 MeV (1/sqrt(E)), 2 mm, 10 keV threshold. Assumed TOF: %.1f ns per station, %.0f%% energy.",
                      sigmaT, 100. * sigmaE));
    box->Draw();
    c.Print(base + "_cuts.pdf");
    c.SaveAs(base + "_cuts.png");
  }
  // 4. the fit after all cuts, linear: counts, fit, its smooth background, and the counts minus that background
  {
    TCanvas c("c4", "fit", 1400, 950);
    auto pTop = new TPad("pTop", "", 0., 0.38, 1., 1.);
    auto pBottom = new TPad("pBottom", "", 0., 0., 1., 0.38);
    pTop->Draw(); pBottom->Draw();
    auto background = new TF1("background", [&](double* x, double* p) {
      return std::exp(p[0] + p[1] * (x[0] - 80.) + p[2] * std::pow(x[0] - 80., 2)); }, 30., 160., 3);
    for (int p = 0; p < 3; ++p) background->SetParameter(p, model->GetParameter(p));
    background->SetNpx(400);
    pTop->cd();
    gPad->SetLeftMargin(0.1); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.12); gPad->SetBottomMargin(0.02);
    axes(hFit, Form(";;entries per %.0f keV in %.1f days", binWidth, days), 0.05);
    hFit->GetXaxis()->SetLabelSize(0.);
    hFit->SetMarkerStyle(20); hFit->SetMarkerSize(0.8); hFit->SetMarkerColor(inkSecondary); hFit->SetLineColor(inkSecondary);
    hFit->SetMinimum(0.);
    hFit->SetMaximum(1.3 * hFit->GetMaximum());
    hFit->Draw("E1");
    model->SetLineColor(blue); model->SetLineWidth(3);
    background->SetLineColor(grey); background->SetLineWidth(2); background->SetLineStyle(2);
    background->Draw("SAME");
    model->Draw("SAME");
    auto legend = legendAt(0.55, 0.6, 0.97, 0.86, 0.045);
    legend->AddEntry(hFit, "counts after all cuts", "pe");
    legend->AddEntry(model, "fit: background + Al L_{#alpha} + Al L_{#beta}", "l");
    legend->AddEntry(background, "its smooth background", "l");
    legend->Draw();
    auto heading = new TLatex();
    heading->SetNDC(); heading->SetTextFont(62); heading->SetTextSize(0.055);
    heading->DrawLatex(0.1, 0.92, Form("Al L lines found by a fit, no truth: L_{#alpha} %.0f +/- %.0f (%.1f sigma), L_{#beta} %.0f +/- %.0f",
                                        nAlpha, dAlpha, dAlpha > 0. ? nAlpha / dAlpha : 0., nBeta, dBeta));
    pBottom->cd();
    gPad->SetLeftMargin(0.1); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.22);
    auto residual = static_cast<TH1D*>(hFit->Clone("residual"));
    for (int b = 1; b <= residual->GetNbinsX(); ++b)
      residual->SetBinContent(b, hFit->GetBinContent(b) - background->Eval(hFit->GetBinCenter(b)));
    auto truthSignal = new TH1D("truthSignal", "", 65, 30., 160.);   // entries with a target X-ray, after all cuts
    for (int b = 1; b <= truthSignal->GetNbinsX(); ++b)
      truthSignal->SetBinContent(b, hSig[nStages - 1]->GetBinContent(hSig[nStages - 1]->FindBin(truthSignal->GetBinCenter(b))) * days);
    axes(residual, ";add-back energy (keV);minus background", 0.08);
    residual->GetYaxis()->SetTitleOffset(0.45);
    residual->GetYaxis()->SetNdivisions(505);
    const double span = 1.4 * std::max(std::abs(residual->GetMaximum()), std::abs(residual->GetMinimum()));
    residual->SetMinimum(-0.5 * span); residual->SetMaximum(span);
    residual->Draw("E1");
    auto zero = new TLine(30., 0., 160., 0.);
    zero->SetLineColor(grey); zero->Draw();
    auto signalOnly = new TF1("signalOnly", [&](double* x, double*) { return model->Eval(x[0]) - background->Eval(x[0]); }, 30., 160., 0);
    signalOnly->SetNpx(400); signalOnly->SetLineColor(blue); signalOnly->SetLineWidth(3);
    signalOnly->Draw("SAME");
    truthSignal->SetLineColor(orange); truthSignal->SetLineWidth(2);
    truthSignal->Draw("HIST SAME");
    auto legend2 = legendAt(0.62, 0.62, 0.97, 0.95, 0.07);
    legend2->AddEntry(signalOnly, "fitted lines", "l");
    legend2->AddEntry(truthSignal, "truth: entries with a target X-ray", "l");
    legend2->Draw();
    c.Print(base + "_fit.pdf");
    c.SaveAs(base + "_fit.png");
  }
  printf("\nplots: %s_stopped.png, _cutvars.png, _cuts.png, _fit.png (and .pdf); cut flow: %s_cutflow.csv\n", base.Data(), base.Data());
}
