// Energy spectrum of the muons when they reach the Compton detector, for the top & middle &
// Compton coincidences (TMC, the selection of the lead x energy scans), versus lead thickness.
// The arrival energy is cmp_ke, the kinetic energy of the primary muon where it first enters a
// layer of the Compton detector (the Si/YSO stack, or the camera when built).
// Input: the spectrum check runs of the two lead x energy scans (macro/lead_energy_scan.mac and
// lead_lowE_scan.mac), repeated with the cmp_ke column by
//   sh scripts/lead_energy_scan.sh --check-only macro/lead_energy_scan_def.mac
//   sh scripts/lead_energy_scan.sh --check-only macro/lead_lowE_scan_def.mac
// - fine scan, <lowDir>/leadspec_t<cm>.root: sea-level muons of 90-510 MeV at the source plane
// - log scan, <logDir>/leadspec_t<cm>.root: muons of 83 MeV-120 GeV, used only outside the
//   fine scan's range, so that every incident energy is counted once
// Each event weighs 1 / live time of its run, so the sums are rates of the EXPACS sea-level
// muons of 83 MeV-120 GeV. TMC events whose muon never entered the Compton detector (the
// coincidence came from secondaries) have no arrival energy; the table gives their share.
//   root -l -b -q analysis/compton_arrival.C
//   root -l -b -q 'analysis/compton_arrival.C("output/leadscan", "output/leadscan_lowE", 0.2, 0.1)'
// Writes output/compton_arrival.csv (spectrum per thickness), output/compton_arrival.pdf and the
// same page as output/compton_arrival_spectra.png.

namespace {
const int kNT = 21;           // lead 0-20 cm
const int kNLog = 51;         // arrival energy 1 MeV - 126 GeV, 10 bins per decade
const int kShown[] = {0, 5, 10, 15, 20};
const int kColours[] = {kBlack, kBlue + 1, kGreen + 2, kOrange + 7, kRed + 1};

struct Sample {
  double eKin = 0., eArr = 0., w = 0.;   // MeV, MeV, per hour
};

// TMC events of one check run with weight 3600 / live time; keep(ekin) selects the incident range
bool ReadRun(const TString& path, const TString& tmc, std::vector<Sample>& reached, double& tmcRate,
             double& eMin, double& eMax, std::function<bool(double)> keep)
{
  TFile file(path);
  auto run = dynamic_cast<TTree*>(file.Get("run"));
  auto events = dynamic_cast<TTree*>(file.Get("events"));
  if (!run || !events || run->GetEntries() != 1) {
    printf("missing or incomplete %s\n", path.Data());
    return false;
  }
  if (!events->GetBranch("cmp_ke")) {
    printf("%s has no cmp_ke column: repeat the check runs with --check-only\n", path.Data());
    return false;
  }
  double liveTime = 0.;
  run->SetBranchAddress("live_time", &liveTime);
  run->SetBranchAddress("e_min", &eMin);
  run->SetBranchAddress("e_max", &eMax);
  run->GetEntry(0);
  const double w = 3600. / liveTime;
  events->SetEstimate(events->GetEntries() + 1);
  const Long64_t n = events->Draw("ekin:cmp_ke", tmc, "goff");
  for (Long64_t i = 0; i < n; ++i) {
    const double eKin = events->GetV1()[i], eArr = events->GetV2()[i];
    if (!keep(eKin)) continue;
    tmcRate += w;
    if (eArr >= 0.) reached.push_back({eKin, eArr, w});
  }
  return true;
}
}  // namespace

void compton_arrival(const char* logDir = "output/leadscan", const char* lowDir = "output/leadscan_lowE",
                     double barThr = 0.2, double comptonThr = 0.1, const char* outBase = "output/compton_arrival")
{
  const TString tmc = Form("Sum$(bar_st==0 && bar_e>%g)>0 && Sum$(bar_st==1 && bar_e>%g)>0 && "
                           "(e_stack0>%g || e_stack1>%g || e_ch0>%g || e_ch1>%g || e_ch2>%g)",
                           barThr, barThr, comptonThr, comptonThr, comptonThr, comptonThr, comptonThr);

  double logEdges[kNLog + 1];
  for (int k = 0; k <= kNLog; ++k) logEdges[k] = std::pow(10., k / 10.);  // MeV
  std::vector<TH1D*> hLog(kNT), hLin(kNT);
  std::vector<std::vector<Sample>> reached(kNT);
  std::vector<double> tmcRate(kNT, 0.);
  std::vector<bool> ok(kNT, false);
  double fineMin = 0., fineMax = 0., logMin = 0., logMax = 0.;
  for (int t = 0; t < kNT; ++t) {
    hLog[t] = new TH1D(Form("hLog%d", t), "", kNLog, logEdges);
    hLin[t] = new TH1D(Form("hLin%d", t), "", 50, 0., 1000.);
    hLog[t]->Sumw2();
    hLin[t]->Sumw2();
    double lo = 0., hi = 0.;
    // fine scan first: its range decides which part of the log scan is used
    if (!ReadRun(Form("%s/leadspec_t%d.root", lowDir, t), tmc, reached[t], tmcRate[t], lo, hi,
                 [](double) { return true; }))
      continue;
    fineMin = lo;
    fineMax = hi;
    if (!ReadRun(Form("%s/leadspec_t%d.root", logDir, t), tmc, reached[t], tmcRate[t], logMin, logMax,
                 [&](double e) { return e < fineMin || e >= fineMax; }))
      continue;
    ok[t] = true;
    for (const auto& s : reached[t]) {
      hLog[t]->Fill(s.eArr, s.w);
      hLin[t]->Fill(s.eArr, s.w);
    }
  }
  int nOk = 0;
  for (bool b : ok) nOk += b;
  if (nOk == 0) return;

  // ---- table
  printf("\n=== Muon energy at the Compton detector, top & middle & Compton coincidences ===\n");
  printf("sea-level muons %.4g MeV - %.4g GeV at the plane: %.4g-%.4g MeV from the fine scan's check runs,\n"
         "the rest from the log scan's; TOF bar > %.2f MeV, Compton layer > %.2f MeV\n",
         logMin, logMax / 1000., fineMin, fineMax, barThr, comptonThr);
  printf("%5s %10s %14s %10s %11s %11s %11s %10s %10s %10s %13s\n", "lead", "TMC /min", "muon reached", "share",
         "mean E_arr", "median", "mean loss", "<100 MeV", "<200 MeV", "<500 MeV", "dR/dE 0-50");
  printf("%5s %10s %14s %10s %11s %11s %11s %10s %10s %10s %13s\n", "cm", "", "/min", "", "GeV", "GeV", "MeV",
         "/h", "/h", "/h", "/h/MeV");
  std::vector<double> below100(kNT), below200(kNT), below500(kNT), err100(kNT), err200(kNT), err500(kNT);
  for (int t = 0; t < kNT; ++t) {
    if (!ok[t]) continue;
    auto& v = reached[t];
    double sw = 0., swE = 0., swLoss = 0., s50 = 0.;
    double b[3] = {}, b2[3] = {};
    const double cut[3] = {100., 200., 500.};
    for (const auto& s : v) {
      sw += s.w;
      swE += s.w * s.eArr;
      swLoss += s.w * (s.eKin - s.eArr);
      if (s.eArr < 50.) s50 += s.w;
      for (int k = 0; k < 3; ++k)
        if (s.eArr < cut[k]) {
          b[k] += s.w;
          b2[k] += s.w * s.w;
        }
    }
    std::sort(v.begin(), v.end(), [](const Sample& a, const Sample& c) { return a.eArr < c.eArr; });
    double run = 0., median = 0.;
    for (const auto& s : v) {
      run += s.w;
      if (run >= 0.5 * sw) {
        median = s.eArr;
        break;
      }
    }
    below100[t] = b[0];
    below200[t] = b[1];
    below500[t] = b[2];
    err100[t] = std::sqrt(b2[0]);
    err200[t] = std::sqrt(b2[1]);
    err500[t] = std::sqrt(b2[2]);
    printf("%5d %10.2f %14.2f %10.3f %11.3f %11.3f %11.1f %10.2f %10.2f %10.1f %13.4f\n", t, tmcRate[t] / 60.,
           sw / 60., tmcRate[t] > 0. ? sw / tmcRate[t] : 0., swE / sw / 1000., median / 1000., swLoss / sw, b[0], b[1],
           b[2], s50 / 50.);
  }
  printf("share: TMC events whose muon itself entered the Compton detector (the others were fired by\n"
         "secondaries). Mean loss: incident minus arrival energy. Statistical errors of the <100/<200 MeV\n"
         "rates are typically 5-15%%; see the plots.\n");

  // ---- CSV: spectrum per thickness
  if (FILE* csv = fopen(TString(outBase) + ".csv", "w")) {
    fprintf(csv, "lead_cm,E_arr_low_MeV,E_arr_high_MeV,rate_per_hour,rate_err_per_hour,dRdE_per_min_per_GeV\n");
    for (int t = 0; t < kNT; ++t) {
      if (!ok[t]) continue;
      for (int i = 1; i <= kNLog; ++i) {
        const double width = hLog[t]->GetBinWidth(i) / 1000.;
        fprintf(csv, "%d,%.6g,%.6g,%.6g,%.6g,%.6g\n", t, hLog[t]->GetBinLowEdge(i), hLog[t]->GetBinLowEdge(i + 1),
                hLog[t]->GetBinContent(i), hLog[t]->GetBinError(i), hLog[t]->GetBinContent(i) / 60. / width);
      }
    }
    fclose(csv);
    printf("spectra: %s.csv\n", outBase);
  }

  // ---- plots
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  gStyle->SetNumberContours(64);
  TCanvas c1("c1", "arrival spectra", 1500, 1100);
  c1.Divide(2, 2);

  // (a) dR/dE, log-log
  c1.cd(1);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  auto legendA = new TLegend(0.15, 0.15, 0.42, 0.42);
  legendA->SetHeader("lead");
  double yMax = 0., yMin = 1e30;
  std::vector<TH1D*> perGeV;
  for (int k = 0; k < 5; ++k) {
    const int t = kShown[k];
    if (!ok[t]) continue;
    auto h = static_cast<TH1D*>(hLog[t]->Clone(Form("hPerGeV%d", t)));
    for (int i = 1; i <= kNLog; ++i) {
      const double width = h->GetBinWidth(i) / 1000.;
      h->SetBinContent(i, hLog[t]->GetBinContent(i) / 60. / width);
      h->SetBinError(i, hLog[t]->GetBinError(i) / 60. / width);
      if (h->GetBinContent(i) > 0.) {
        yMax = std::max(yMax, h->GetBinContent(i));
        yMin = std::min(yMin, h->GetBinContent(i));
      }
    }
    h->SetLineColor(kColours[k]);
    h->SetMarkerColor(kColours[k]);
    h->SetMarkerStyle(20 + k);
    h->SetMarkerSize(0.8);
    h->SetLineWidth(2);
    perGeV.push_back(h);
    legendA->AddEntry(h, Form("%d cm", t), "pl");
  }
  for (std::size_t k = 0; k < perGeV.size(); ++k) {
    auto h = perGeV[k];
    if (k == 0) {
      h->SetTitle("muons at the Compton detector (TMC);arrival kinetic energy (MeV);dR/dE (min^{-1} GeV^{-1})");
      h->GetYaxis()->SetRangeUser(std::max(0.5 * yMin, 1e-5 * yMax), 3. * yMax);
      h->Draw("E1 HIST");
    } else {
      h->Draw("E1 HIST SAME");
    }
  }
  legendA->Draw();

  // (b) slow muons, linear 0-1000 MeV
  c1.cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  auto legendB = new TLegend(0.62, 0.62, 0.89, 0.89);
  legendB->SetHeader("lead");
  double yMaxB = 0.;
  std::vector<TH1D*> perMeV;
  for (int k = 0; k < 5; ++k) {
    const int t = kShown[k];
    if (!ok[t]) continue;
    auto h = static_cast<TH1D*>(hLin[t]->Clone(Form("hPerMeV%d", t)));
    h->Scale(1. / h->GetBinWidth(1));  // per hour per MeV
    yMaxB = std::max(yMaxB, h->GetMaximum());
    h->SetLineColor(kColours[k]);
    h->SetMarkerColor(kColours[k]);
    h->SetMarkerStyle(20 + k);
    h->SetMarkerSize(0.8);
    h->SetLineWidth(2);
    perMeV.push_back(h);
    legendB->AddEntry(h, Form("%d cm", t), "pl");
  }
  for (std::size_t k = 0; k < perMeV.size(); ++k) {
    auto h = perMeV[k];
    if (k == 0) {
      h->SetTitle("slow muons at the Compton detector (TMC);arrival kinetic energy (MeV);dR/dE (h^{-1} MeV^{-1})");
      h->GetYaxis()->SetRangeUser(0., 1.15 * yMaxB);
      h->Draw("E1 HIST");
    } else {
      h->Draw("E1 HIST SAME");
    }
  }
  legendB->Draw();

  // (c) map: arrival energy versus lead thickness, per hour per decade
  c1.cd(3);
  gPad->SetLogx();
  gPad->SetRightMargin(0.16);
  gPad->SetBottomMargin(0.12);
  auto hMap = new TH2D("hMap", "muons at the Compton detector (TMC) per hour per decade;arrival kinetic energy (MeV);lead thickness (cm)",
                       kNLog, logEdges, kNT, -0.5, kNT - 0.5);
  for (int t = 0; t < kNT; ++t)
    if (ok[t])
      for (int i = 1; i <= kNLog; ++i) hMap->SetBinContent(i, t + 1, hLog[t]->GetBinContent(i) * 10.);
  hMap->GetXaxis()->SetTitleOffset(1.2);
  hMap->GetXaxis()->SetRangeUser(10., logEdges[kNLog]);
  hMap->Draw("COLZ");

  // (d) slow-muon rates versus thickness
  c1.cd(4);
  gPad->SetLogy();
  gPad->SetGridy();
  gPad->SetBottomMargin(0.12);
  auto mg = new TMultiGraph("mg", "muons arriving slower than ...;lead thickness (cm);rate (h^{-1})");
  auto legendD = new TLegend(0.6, 0.15, 0.89, 0.33);
  const std::vector<double>* rates[3] = {&below100, &below200, &below500};
  const std::vector<double>* errs[3] = {&err100, &err200, &err500};
  const char* names[3] = {"< 100 MeV", "< 200 MeV", "< 500 MeV"};
  const int coloursD[3] = {kRed + 1, kOrange + 7, kBlue + 1};
  for (int k = 0; k < 3; ++k) {
    auto g = new TGraphErrors();
    for (int t = 0; t < kNT; ++t)
      if (ok[t] && (*rates[k])[t] > 0.) {
        const int n = g->GetN();
        g->SetPoint(n, t, (*rates[k])[t]);
        g->SetPointError(n, 0., (*errs[k])[t]);
      }
    g->SetMarkerStyle(20 + k);
    g->SetMarkerColor(coloursD[k]);
    g->SetLineColor(coloursD[k]);
    mg->Add(g, "PL");
    legendD->AddEntry(g, names[k], "pl");
  }
  mg->Draw("A");
  legendD->Draw();

  c1.Print(TString(outBase) + ".pdf");
  c1.SaveAs(TString(outBase) + "_spectra.png");
  printf("plots: %s.pdf, %s_spectra.png\n", outBase, outBase);
}
