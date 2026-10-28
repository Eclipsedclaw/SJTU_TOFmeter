// Health check of a run of scripts/muonic_xray_stack.sh, on its finished processes (log/job<k>.done): do the processes
// agree with each other, and do the stack hits look right?
//   root -l -b -q 'analysis/muonic_xray_stack_check.C("output/muonic_xray_stack")'
// Per process: events written and mu- stopped in the target (with energy in the stack) per hour of exposure, against
// a Poisson spread: chi2/ndf near 1 when the processes agree; all seeds different. Over all hits: impossible values
// (energy <= 0, time < 0 or not a number, position outside its layer or the 20 x 20 cm of the stack), the largest
// energy and the most hits of an event. Plots: energy of the Si and of the YSO hits (log, 1 keV - 100 MeV: crossing
// muons make the MIP peaks), time of the hits after the first hit of their event (prompt, then mu- decays and
// captures), x-y of the hits, events written per hour for every process.
// Writes <dir>/muonic_xray_stack_check.png/.pdf.

#include <set>

namespace {
// Poisson chi2 / ndf of counts that should all have the same mean
double PoissonChi2(const std::vector<double>& counts)
{
  if (counts.size() < 2) return 0.;
  double mean = 0.;
  for (double c : counts) mean += c;
  mean /= counts.size();
  double chi2 = 0.;
  for (double c : counts) chi2 += (c - mean) * (c - mean) / mean;
  return chi2 / (counts.size() - 1);
}

// most probable value above an energy: the middle of the fullest bin (log bins: the geometric middle)
double PeakAbove(const TH1D* h, double low)
{
  int best = 0;
  for (int b = h->GetXaxis()->FindFixBin(low); b <= h->GetNbinsX(); ++b)
    if (best == 0 || h->GetBinContent(b) > h->GetBinContent(best)) best = b;
  return best > 0 ? std::sqrt(h->GetBinLowEdge(best) * h->GetBinLowEdge(best + 1)) : 0.;
}
}  // namespace

void muonic_xray_stack_check(const char* dir = "output/muonic_xray_stack")
{
  // the finished processes
  std::vector<std::pair<int, std::string>> files;
  if (void* dirp = gSystem->OpenDirectory(dir)) {
    while (const char* entry = gSystem->GetDirEntry(dirp)) {
      const std::string name = entry;
      if (name.rfind("stk_", 0) != 0 || name.size() <= 9 || name.substr(name.size() - 5) != ".root") continue;
      const std::string k = name.substr(4, name.size() - 9);
      if (gSystem->AccessPathName(Form("%s/log/job%s.done", dir, k.c_str()))) continue;   // true: not finished
      files.push_back({std::stoi(k), std::string(dir) + "/" + name});
    }
    gSystem->FreeDirectory(dirp);
  }
  std::sort(files.begin(), files.end());
  if (files.empty()) {
    printf("no finished process in %s\n", dir);
    return;
  }

  // the five layers, from the bottom: 0, 2, 4 Si and 1, 3 YSO (z ranges in cm)
  const double zLow[5] = {17.8, 21.25, 24.4, 27.85, 31.0}, zHigh[5] = {18.4, 21.55, 25.0, 28.15, 31.6};
  const int nLog = 100;
  std::vector<double> edges(nLog + 1);
  for (int i = 0; i <= nLog; ++i) edges[i] = std::pow(10., -3. + 5. * i / nLog);   // 1 keV - 100 MeV
  auto hSi = new TH1D("hSi", "", nLog, edges.data());
  auto hYSO = new TH1D("hYSO", "", nLog, edges.data());
  auto hTime = new TH1D("hTime", "", 200, 0., 10.);   // us after the first hit of the event
  auto hXY = new TH2D("hXY", "", 50, -12.5, 12.5, 50, -12.5, 12.5);
  auto rates = new TGraphErrors();
  std::vector<double> written, stops;
  std::set<int> seeds;
  long long hits = 0, badEnergy = 0, badTime = 0, outside = 0, notANumber = 0, duplicatedSeeds = 0;
  double liveTime = 0., maxEnergy = 0., maxEkin = 0.;
  std::size_t maxHits = 0;
  long long above100 = 0;   // hits above 100 MeV
  int maxProcess = 0, maxEvent = 0, maxPdg = 0, maxEventHits = 0;   // the event of the largest hit
  std::string maxVolume, maxProcessName;
  for (const auto& [k, path] : files) {
    TFile file(path.c_str());
    auto run = dynamic_cast<TTree*>(file.Get("run"));
    auto events = dynamic_cast<TTree*>(file.Get("events"));
    if (!run || !events || run->GetEntries() != 1) {
      printf("process %d: incomplete file\n", k);
      continue;
    }
    double lt = 0.;
    int seed = 0;   // an int column of the run tree
    run->SetBranchAddress("live_time", &lt);
    run->SetBranchAddress("seed", &seed);
    run->GetEntry(0);
    liveTime += lt;
    if (!seeds.insert(seed).second) ++duplicatedSeeds;
    int pdg = 0, eventID = 0;
    double ekin = 0.;
    char endProc[128] = "", endVol[128] = "";
    std::vector<int>* layer = nullptr;
    std::vector<double>*e = nullptr, *t = nullptr, *x = nullptr, *y = nullptr, *z = nullptr;
    events->SetBranchStatus("*", 0);
    for (const char* column : {"eventID", "pdg", "ekin", "end_proc", "end_vol", "stk_layer", "stk_e", "stk_t", "stk_x", "stk_y",
                               "stk_z"})
      events->SetBranchStatus(column, 1);
    events->SetBranchAddress("eventID", &eventID);
    events->SetBranchAddress("ekin", &ekin);
    events->SetBranchAddress("pdg", &pdg);
    events->SetBranchAddress("end_proc", endProc);
    events->SetBranchAddress("end_vol", endVol);
    events->SetBranchAddress("stk_layer", &layer);
    events->SetBranchAddress("stk_e", &e);
    events->SetBranchAddress("stk_t", &t);
    events->SetBranchAddress("stk_x", &x);
    events->SetBranchAddress("stk_y", &y);
    events->SetBranchAddress("stk_z", &z);
    long stopped = 0;
    for (Long64_t i = 0; i < events->GetEntries(); ++i) {
      events->GetEntry(i);
      if (pdg == 13 && std::string(endProc) == "muMinusCaptureAtRest" && std::string(endVol).rfind("Target", 0) == 0) ++stopped;
      maxHits = std::max(maxHits, e->size());
      double first = 1e30;
      for (double ti : *t) first = std::min(first, ti);
      for (std::size_t j = 0; j < e->size(); ++j) {
        ++hits;
        const int l = (*layer)[j];
        const double ej = (*e)[j], tj = (*t)[j], xj = (*x)[j], yj = (*y)[j], zj = (*z)[j];
        if (!std::isfinite(ej) || !std::isfinite(tj) || !std::isfinite(xj) || !std::isfinite(yj) || !std::isfinite(zj)) {
          ++notANumber;
          continue;
        }
        if (ej <= 0.) ++badEnergy;
        if (tj < 0.) ++badTime;
        if (l < 0 || l > 4 || zj < zLow[l] - 1e-6 || zj > zHigh[l] + 1e-6 || std::abs(xj) > 10. + 1e-6 || std::abs(yj) > 10. + 1e-6)
          ++outside;
        if (ej > 100.) ++above100;
        if (ej > maxEnergy) {
          maxEnergy = ej;
          maxProcess = k;
          maxEvent = eventID;
          maxPdg = pdg;
          maxEkin = ekin;
          maxVolume = endVol;
          maxProcessName = endProc;
          maxEventHits = int(e->size());
        }
        (l % 2 ? hYSO : hSi)->Fill(ej);
        hTime->Fill((tj - first) / 1000.);
        hXY->Fill(xj, yj);
      }
    }
    const double hours = lt / 3600.;
    written.push_back(events->GetEntries());
    stops.push_back(stopped);
    const int p = rates->GetN();
    rates->SetPoint(p, k, events->GetEntries() / hours);
    rates->SetPointError(p, 0., std::sqrt(double(events->GetEntries())) / hours);
  }

  // ---- summary
  const double hoursPerProcess = liveTime / 3600. / written.size();
  double sumWritten = 0., sumStops = 0.;
  for (std::size_t i = 0; i < written.size(); ++i) { sumWritten += written[i]; sumStops += stops[i]; }
  printf("\n=== Health check of %s: %zu finished processes, %.2f days of exposure ===\n", dir, written.size(), liveTime / 86400.);
  printf("seeds: %zu different, %lld used twice\n", seeds.size(), duplicatedSeeds);
  printf("events written: %.0f per process (%.0f per hour), spread chi2/ndf = %.2f (1 = Poisson)\n", sumWritten / written.size(),
         sumWritten / written.size() / hoursPerProcess, PoissonChi2(written));
  printf("mu- stopped in the target with energy in the stack: %.0f in all, %.1f per day, spread chi2/ndf = %.2f\n", sumStops,
         sumStops / (liveTime / 86400.), PoissonChi2(stops));
  printf("hits: %lld; not a number %lld, energy <= 0: %lld, time < 0: %lld, outside their layer: %lld\n", hits, notANumber,
         badEnergy, badTime, outside);
  printf("largest hit energy %.1f MeV; most hits in one event %zu\n", maxEnergy, maxHits);
  printf("most probable hit energy above 0.5 MeV (crossing muons): Si %.2f MeV, YSO %.2f MeV\n", PeakAbove(hSi, 0.5),
         PeakAbove(hYSO, 0.5));
  printf("hits above 100 MeV: %lld; the largest in process %d, event %d: a %s of %.1f GeV that ended in %s (%s), %d hits\n",
         above100, maxProcess, maxEvent, maxPdg == 13 ? "mu-" : maxPdg == -13 ? "mu+" : "particle", maxEkin / 1000.,
         maxVolume.c_str(), maxProcessName.c_str(), maxEventHits);

  // ---- plots
  const int blue = TColor::GetColor("#2a78d6"), orange = TColor::GetColor("#eb6834"), inkSecondary = TColor::GetColor("#52514e");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  TCanvas c("c", "stack run check", 1400, 1000);
  auto title = new TLatex();
  title->SetNDC();
  title->SetTextFont(62);
  title->SetTextSize(0.025);
  title->DrawLatex(0.05, 0.965, Form("Health check of the Si/YSO stack run: %zu processes, %.1f days of exposure, %lld hits",
                                     written.size(), liveTime / 86400., hits));
  auto pads = new TPad("pads", "", 0., 0., 1., 0.94);
  pads->Draw();
  pads->Divide(2, 2, 0.01, 0.01);
  auto style = [&](TH1* h, const char* titles) {
    h->SetTitle(titles);
    for (auto axis : {h->GetXaxis(), h->GetYaxis()}) {
      axis->SetTitleSize(0.05);
      axis->SetLabelSize(0.045);
      axis->SetTitleColor(inkSecondary);
    }
  };
  pads->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.13);
  style(hSi, ";hit energy (MeV);hits per bin");
  hSi->SetLineColor(blue);
  hYSO->SetLineColor(orange);
  hSi->SetLineWidth(2);
  hYSO->SetLineWidth(2);
  hSi->SetMaximum(3. * std::max(hSi->GetMaximum(), hYSO->GetMaximum()));
  hSi->Draw("HIST");
  hYSO->Draw("HIST SAME");
  auto legend = new TLegend(0.16, 0.78, 0.5, 0.9);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.045);
  legend->AddEntry(hSi, "Si layers (6 mm)", "l");
  legend->AddEntry(hYSO, "YSO layers (3 mm)", "l");
  legend->Draw();
  pads->cd(2);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.13);
  style(hTime, ";time after the first hit of the event (#mus);hits per 50 ns");
  hTime->SetLineColor(blue);
  hTime->SetLineWidth(2);
  hTime->Draw("HIST");
  pads->cd(3);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.14);
  gPad->SetBottomMargin(0.13);
  style(hXY, ";hit x (cm);hit y (cm)");
  hXY->Draw("COLZ");
  pads->cd(4);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.13);
  auto frame = gPad->DrawFrame(0., 0., files.back().first + 1., 1.);
  frame->SetTitle(";process;events written per hour of exposure");
  for (auto axis : {frame->GetXaxis(), frame->GetYaxis()}) {
    axis->SetTitleSize(0.05);
    axis->SetLabelSize(0.045);
    axis->SetTitleColor(inkSecondary);
  }
  const double mean = sumWritten / written.size() / hoursPerProcess;
  frame->SetMinimum(0.9 * mean);
  frame->SetMaximum(1.1 * mean);
  frame->Draw("AXIS");
  rates->SetMarkerStyle(20);
  rates->SetMarkerSize(0.5);
  rates->SetMarkerColor(blue);
  rates->SetLineColor(blue);
  rates->Draw("P");
  auto meanLine = new TLine(0., mean, files.back().first + 1., mean);
  meanLine->SetLineColor(inkSecondary);
  meanLine->SetLineStyle(2);
  meanLine->Draw();
  auto note = new TLatex();
  note->SetNDC();
  note->SetTextFont(42);
  note->SetTextSize(0.045);
  note->SetTextColor(inkSecondary);
  note->DrawLatex(0.18, 0.85, Form("spread chi2/ndf = %.2f (1 = Poisson)", PoissonChi2(written)));
  TString base = Form("%s/muonic_xray_stack_check", dir);
  c.Print(base + ".pdf");
  c.SaveAs(base + ".png");
  printf("plots: %s.png, .pdf\n", base.Data());
}
