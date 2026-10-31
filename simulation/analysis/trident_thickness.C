// Detectable muon tridents per year versus lead thickness, from a scan of scripts/trident_scan.sh
// (macro/trident_scan_def.mac): every <dir>/t<cm>/ with trident_cosmic_*.root is one point.
//   root -l -b -q 'analysis/trident_thickness.C("output/tridents_F1e3")'
// Detectable as in analysis/trident_rate.C: a trident of the primary muon whose three muons (the
// primary after it and the produced mu+ and mu-) each cross a bar of the top TOF and a bar of the
// middle TOF; also the parts with at least 2 and with 3 different bars hit in each TOF. Rates are
// sums of weights (1/F) over the live time, errors statistical. The Geant4 version of each thickness
// is read from its first job log (log/job1.log); points of another version get open markers.
// Writes <dir>/trident_thickness.png, .pdf and .csv.

namespace {
struct ThicknessPoint {
  double lead = 0., liveTime = 0., bias = 0.;
  Long64_t generated = 0;
  double det = 0., det2 = 0., two = 0., two2 = 0., three = 0., three2 = 0., leadLayer = 0.;
  long nDet = 0, nTwo = 0, nThree = 0;
  std::string geant4 = "?";
};

// "geant4-11-04-patch-03" in the run log -> "11.4.3"
std::string Geant4Version(const std::string& log)
{
  std::ifstream in(log);
  std::string line;
  while (std::getline(in, line)) {
    const auto pos = line.find("Geant4 version Name: geant4-");
    if (pos == std::string::npos) continue;
    int major = 0, minor = 0, patch = 0;
    if (std::sscanf(line.c_str() + pos + 28, "%d-%d-patch-%d", &major, &minor, &patch) == 3)
      return Form("%d.%d.%d", major, minor, patch);
    return line.substr(pos + 21, line.find(' ', pos + 21) - pos - 21);
  }
  return "?";
}

// 8e8 -> "8 #times 10^{8}" (TLatex), numbers below 10^5 as they are: 1000 -> "1000"
std::string PowerOfTen(double x)
{
  if (x < 1e5) return Form("%g", x);
  const int exponent = int(std::floor(std::log10(x) + 1e-9));
  const double mantissa = x / std::pow(10., exponent);
  if (std::abs(mantissa - 1.) < 1e-6) return Form("10^{%d}", exponent);
  return Form("%.3g #times 10^{%d}", mantissa, exponent);
}
}  // namespace

void trident_thickness(const char* dir = "output/tridents_F1e3")
{
  // the thickness folders t<cm>
  std::vector<double> leads;
  if (void* dirp = gSystem->OpenDirectory(dir)) {
    while (const char* entry = gSystem->GetDirEntry(dirp)) {
      const std::string name = entry;
      if (name.size() > 1 && name[0] == 't' && name.find_first_not_of("0123456789.", 1) == std::string::npos)
        leads.push_back(std::stod(name.substr(1)));
    }
    gSystem->FreeDirectory(dirp);
  }
  std::sort(leads.begin(), leads.end());

  const double year = 365. * 86400.;
  std::vector<ThicknessPoint> points;
  for (double t : leads) {
    const std::string sub = Form("%s/t%g", dir, t);
    std::vector<std::string> files;
    if (void* dirp = gSystem->OpenDirectory(sub.c_str())) {
      while (const char* entry = gSystem->GetDirEntry(dirp)) {
        const std::string name = entry;
        if (name.rfind("trident_cosmic_", 0) == 0 && name.size() > 5 && name.substr(name.size() - 5) == ".root")
          files.push_back(sub + "/" + name);
      }
      gSystem->FreeDirectory(dirp);
    }
    if (files.empty()) continue;
    ThicknessPoint p;
    p.lead = t;
    p.geant4 = Geant4Version(sub + "/log/job1.log");
    double maxWeight = 0.;
    for (const auto& path : files) {
      TFile file(path.c_str());
      auto run = dynamic_cast<TTree*>(file.Get("run"));
      auto events = dynamic_cast<TTree*>(file.Get("events"));
      if (!run || !events || run->GetEntries() != 1 || !events->GetBranch("tmu_st")) {
        printf("skipping %s: incomplete or without TOF information per muon\n", path.c_str());
        continue;
      }
      int n = 0;
      double lt = 0., leadCm = 0.;
      run->SetBranchAddress("n_events", &n);
      run->SetBranchAddress("live_time", &lt);
      run->SetBranchAddress("lead_cm", &leadCm);
      double ground = -70.;   // z of the ground (setup bottom); -70 cm in files made before it was put at 0
      if (run->GetBranch("setup_bottom_z")) run->SetBranchAddress("setup_bottom_z", &ground);
      run->GetEntry(0);
      const double leadBottom = ground + 91.;   // the lead layer starts 91 cm above the ground
      if (std::abs(leadCm - t) > 1e-6) printf("warning: %s has %g cm of lead, not %g\n", path.c_str(), leadCm, t);
      p.generated += n;
      p.liveTime += lt;

      std::vector<double>*w = nullptr, *x = nullptr, *y = nullptr, *z = nullptr;
      std::vector<int>*parent = nullptr, *proc = nullptr, *muParent = nullptr, *st = nullptr, *top = nullptr, *mid = nullptr;
      int primSt = 0, primTop = -1, primMid = -1;
      events->SetBranchAddress("tri_w", &w);
      events->SetBranchAddress("tri_x", &x);
      events->SetBranchAddress("tri_y", &y);
      events->SetBranchAddress("tri_z", &z);
      events->SetBranchAddress("tri_parent", &parent);
      events->SetBranchAddress("tmu_proc", &proc);
      events->SetBranchAddress("tmu_parent", &muParent);
      events->SetBranchAddress("tmu_st", &st);
      events->SetBranchAddress("tmu_top", &top);
      events->SetBranchAddress("tmu_mid", &mid);
      events->SetBranchAddress("prim_st", &primSt);
      events->SetBranchAddress("prim_top", &primTop);
      events->SetBranchAddress("prim_mid", &primMid);
      for (Long64_t i = 0; i < events->GetEntries(); ++i) {
        events->GetEntry(i);
        int primaryTrident = -1;
        for (std::size_t k = 0; k < w->size(); ++k) {
          maxWeight = std::max(maxWeight, (*w)[k]);
          // every trident in the lead layer (as in trident_rate.C): checks the cross section per cm
          if (std::max(std::abs((*x)[k]), std::abs((*y)[k])) < 45. && (*z)[k] > leadBottom && (*z)[k] < leadBottom + t)
            p.leadLayer += (*w)[k];
          if ((*parent)[k] == 1 && primaryTrident < 0) primaryTrident = int(k);
        }
        // the primary and its produced pair all cross the top (bit 0) and middle (bit 1) TOF
        if (primaryTrident < 0 || (primSt & 3) != 3) continue;
        std::set<int> topBars = {primTop}, midBars = {primMid};
        int produced = 0;
        bool ok = true;
        for (std::size_t k = 0; k < proc->size(); ++k) {
          if ((*proc)[k] != 49 || (*muParent)[k] != 1) continue;
          ++produced;
          if (((*st)[k] & 3) != 3) ok = false;
          topBars.insert((*top)[k]);
          midBars.insert((*mid)[k]);
        }
        if (!ok || produced < 2) continue;
        const double wk = (*w)[primaryTrident];
        p.det += wk;
        p.det2 += wk * wk;
        ++p.nDet;
        const int separate = int(std::min(topBars.size(), midBars.size()));
        if (separate >= 2) { p.two += wk; p.two2 += wk * wk; ++p.nTwo; }
        if (separate >= 3) { p.three += wk; p.three2 += wk * wk; ++p.nThree; }
      }
    }
    if (p.liveTime <= 0.) continue;
    p.bias = maxWeight > 0. ? 1. / maxWeight : 1.;
    points.push_back(p);
  }
  if (points.empty()) {
    printf("no t<cm>/trident_cosmic_*.root in %s (run scripts/trident_scan.sh)\n", dir);
    return;
  }

  // table and CSV, rates per year
  const std::string base = Form("%s/trident_thickness", dir);
  std::ofstream csv(base + ".csv");
  csv << "lead_cm,geant4,muons,live_time_h,bias_F,detectable_per_year,err,two_bars_per_year,err,three_bars_per_year,err,"
         "detectable_simulated\n";
  printf("\n=== Detectable muon tridents per year versus lead thickness (%s) ===\n", dir);
  printf("%8s %8s %10s %10s %18s %18s %18s %10s\n", "lead cm", "Geant4", "muons", "exposure h", "detectable /yr",
         ">=2 bars /yr", "3 bars /yr", "simulated");
  for (const auto& p : points) {
    const double s = year / p.liveTime;
    printf("%8g %8s %10.3g %10.1f %10.2f +- %5.2f %10.2f +- %5.2f %10.3f +- %5.3f %10ld\n", p.lead, p.geant4.c_str(),
           double(p.generated), p.liveTime / 3600., p.det * s, std::sqrt(p.det2) * s, p.two * s, std::sqrt(p.two2) * s,
           p.three * s, std::sqrt(p.three2) * s, p.nDet);
    csv << p.lead << "," << p.geant4 << "," << p.generated << "," << p.liveTime / 3600. << "," << p.bias << ","
        << p.det * s << "," << std::sqrt(p.det2) * s << "," << p.two * s << "," << std::sqrt(p.two2) * s << ","
        << p.three * s << "," << std::sqrt(p.three2) * s << "," << p.nDet << "\n";
  }

  // the plot: one series, markers filled for the most common Geant4 version, open for others
  std::map<std::string, int> versionCount;
  for (const auto& p : points) ++versionCount[p.geant4];
  std::string mainVersion = points.front().geant4;
  for (const auto& [v, n] : versionCount)
    if (n > versionCount[mainVersion]) mainVersion = v;
  const int blue = TColor::GetColor("#2a78d6");      // series colour
  const int inkSecondary = TColor::GetColor("#52514e");
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetGridColor(TColor::GetColor("#e3e2dc"));
  gStyle->SetGridStyle(1);
  TCanvas c("c", "detectable tridents versus lead thickness", 1200, 900);
  auto padPlot = new TPad("padPlot", "", 0., 0.3, 1., 1.);
  auto padText = new TPad("padText", "", 0., 0., 1., 0.3);
  for (auto pad : {padPlot, padText}) pad->Draw();

  padPlot->cd();
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gPad->SetGridy();
  double yMax = 0., xMax = 0.;
  for (const auto& p : points) {
    yMax = std::max(yMax, (p.det + std::sqrt(p.det2)) * year / p.liveTime);
    xMax = std::max(xMax, p.lead);
  }
  auto frame = gPad->DrawFrame(0., 0., 1.1 * xMax, 1.3 * yMax);
  frame->SetTitle(";lead layer thickness (cm);detectable tridents per year");
  for (auto axis : {frame->GetXaxis(), frame->GetYaxis()}) {
    axis->SetTitleSize(0.05);
    axis->SetLabelSize(0.045);
    axis->SetTitleColor(inkSecondary);
  }
  frame->GetXaxis()->SetTitleOffset(1.05);
  frame->GetYaxis()->SetTitleOffset(0.9);
  auto legend = new TLegend(0.13, 0.66, 0.6, 0.86);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.042);
  legend->SetTextFont(42);
  auto label = new TLatex();
  label->SetTextSize(0.04);
  label->SetTextFont(42);
  label->SetTextColor(inkSecondary);
  label->SetTextAlign(21);
  for (bool main : {true, false}) {
    auto graph = new TGraphErrors();
    std::string versions;
    for (const auto& p : points) {
      if ((p.geant4 == mainVersion) != main) continue;
      const double s = year / p.liveTime;
      const int i = graph->GetN();
      graph->SetPoint(i, p.lead, p.det * s);
      graph->SetPointError(i, 0., std::sqrt(p.det2) * s);
      if (versions.find(p.geant4) == std::string::npos) versions += (versions.empty() ? "" : ", ") + p.geant4;
    }
    if (graph->GetN() == 0) continue;
    graph->SetMarkerStyle(main ? 20 : 24);
    graph->SetMarkerSize(2.);
    graph->SetMarkerColor(blue);
    graph->SetLineColor(blue);
    graph->SetLineWidth(2);
    graph->Draw("P SAME");
    std::string thicknesses;
    for (int i = 0; i < graph->GetN(); ++i) thicknesses += Form("%s%g", i ? ", " : "", graph->GetX()[i]);
    legend->AddEntry(graph, Form("Geant4 %s (%s cm)", versions.c_str(), thicknesses.c_str()), "pe");
  }
  // values at the two ends of the scan
  for (const auto* p : {&points.front(), &points.back()}) {
    const double s = year / p->liveTime;
    label->DrawLatex(p->lead, (p->det + std::sqrt(p->det2)) * s + 0.04 * yMax, Form("%.1f / year", p->det * s));
  }
  if (versionCount.size() > 1) legend->Draw();
  auto title = new TLatex();
  title->SetNDC();
  title->SetTextSize(0.055);
  title->DrawLatex(0.1, 0.92, "Detectable muon tridents versus lead thickness");

  padText->cd();
  auto box = new TPaveText(0.01, 0.04, 0.99, 0.96, "NDC");
  box->SetFillColor(kWhite);
  box->SetBorderSize(1);
  box->SetTextAlign(12);
  box->SetTextSize(0.072);
  box->SetTextFont(42);
  box->SetTextColor(inkSecondary);
  std::string rows;
  for (const auto& p : points) {
    const double s = year / p.liveTime;
    rows += Form("%s%g cm: %.1f #pm %.1f", rows.empty() ? "" : ";   ", p.lead, p.det * s, std::sqrt(p.det2) * s);
  }
  box->AddText(Form("Detectable per year: %s", rows.c_str()));
  rows.clear();
  for (const auto& p : points) {
    const double s = year / p.liveTime;
    rows += Form("%s%g cm: %.1f", rows.empty() ? "" : ";   ", p.lead, p.two * s);
  }
  box->AddText(Form("... of them with #geq 2 different bars hit in each TOF: %s", rows.c_str()));
  rows.clear();
  for (const auto& p : points) rows += Form("%s%g cm: %.4f", rows.empty() ? "" : ";   ", p.lead, p.leadLayer * 86400. / p.liveTime / p.lead);
  box->AddText(Form("All tridents in the lead layer, per cm of lead and per day: %s", rows.c_str()));
  box->AddText("Detectable: the muon after the trident and its #mu^{+}#mu^{-} pair each cross a bar of the top and of the middle TOF.");
  box->AddText(Form("Each point: %s cosmic muons above 0.85 GeV = %.0f h at sea level (EXPACS), 2 x 2 m plane, trident bias F = %s;",
                    PowerOfTen(double(points.front().generated)).c_str(), points.front().liveTime / 3600.,
                    PowerOfTen(points.front().bias).c_str()));
  box->AddText("statistical errors (smaller than the markers). Geant4 trident model: total cross section checked (Russell et al. 1971), energy sharing not.");
  box->Draw();
  c.Print((base + ".pdf").c_str());
  c.SaveAs((base + ".png").c_str());
  printf("plots: %s.png, %s.pdf; table: %s.csv\n", base.c_str(), base.c_str(), base.c_str());
}
