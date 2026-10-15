// Rate of muon tridents in the TOF-meter, from scripts/trident_cosmic.sh (macro/trident_cosmic.mac):
// sea-level cosmic muons above 0.85 GeV, tridents biased (weight tri_w, 1e-5), only events with a
// trident written. Sums of weights over the live time are the true rates.
//   root -l -b -q 'analysis/trident_rate.C("output/tridents_t10")'
// Detectable: a trident of the primary muon whose three muons, the primary after it and the produced
// mu+ and mu-, each deposit energy in a bar of the top TOF and in a bar of the middle TOF. Also
// counted: how many different (4 cm) bars the three muons hit in each of the two TOFs; three
// different bars in both means three separate tracks. Where the tridents happen: lead layer, lead
// walls, steel plate, elsewhere (TOF bars, Al plates, corner beams, stack, air).
// The bias factor F is read from the weights (F = 1 / largest trident weight); every rate printed
// or plotted is corrected for it, while the simulated counts are F times the real numbers.
// Writes <dir>/trident_rate.pdf and .png.

namespace {
struct WeightSum {
  double w = 0., w2 = 0.;
  long n = 0;
  void Add(double x) { w += x; w2 += x * x; ++n; }
};
}  // namespace

void trident_rate(const char* dir = "output/tridents_t10")
{
  std::vector<TString> files;
  if (void* dirp = gSystem->OpenDirectory(dir)) {
    while (const char* entry = gSystem->GetDirEntry(dirp)) {
      const TString name = entry;
      if (name.BeginsWith("trident_cosmic_") && name.EndsWith(".root")) files.push_back(Form("%s/%s", dir, entry));
    }
    gSystem->FreeDirectory(dirp);
  }
  if (files.empty()) {
    printf("no trident_cosmic_*.root in %s (run scripts/trident_cosmic.sh)\n", dir);
    return;
  }

  Long64_t generated = 0;
  double liveTime = 0., lead = -1.;
  WeightSum all, inLead, inWalls, inSteel, elsewhere, detectable, twoBars, threeBars;
  double maxWeight = 0.;  // weight of a first trident = 1 / F
  double logEdges[31];
  for (int i = 0; i <= 30; ++i) logEdges[i] = 0.85 * std::pow(1000. / 0.85, i / 30.);  // GeV
  auto hAll = new TH1D("hAll", "", 30, logEdges);
  auto hDetectable = new TH1D("hDetectable", "", 30, logEdges);
  auto hBars = new TH1D("hBars", "", 3, 0.5, 3.5);
  for (auto h : {hAll, hDetectable, hBars}) h->Sumw2();
  std::vector<std::pair<double, double>> detectableEnergies;  // incoming energy, weight

  for (const auto& path : files) {
    TFile file(path);
    auto run = dynamic_cast<TTree*>(file.Get("run"));
    auto events = dynamic_cast<TTree*>(file.Get("events"));
    if (!run || !events || !events->GetBranch("tmu_st")) {
      printf("skipping %s: incomplete or without TOF information per muon\n", path.Data());
      continue;
    }
    int n = 0;
    double lt = 0., leadCm = 0.;
    run->SetBranchAddress("n_events", &n);
    run->SetBranchAddress("live_time", &lt);
    run->SetBranchAddress("lead_cm", &leadCm);
    run->GetEntry(0);
    generated += n;
    liveTime += lt;
    if (lead >= 0. && std::abs(leadCm - lead) > 1e-6) printf("warning: %s has %g cm of lead, not %g\n", path.Data(), leadCm, lead);
    lead = leadCm;

    std::vector<double>*x = nullptr, *y = nullptr, *z = nullptr, *w = nullptr, *ein = nullptr;
    std::vector<int>*parent = nullptr, *proc = nullptr, *muParent = nullptr, *st = nullptr, *top = nullptr, *mid = nullptr;
    int primSt = 0, primTop = -1, primMid = -1;
    events->SetBranchAddress("tri_x", &x);
    events->SetBranchAddress("tri_y", &y);
    events->SetBranchAddress("tri_z", &z);
    events->SetBranchAddress("tri_w", &w);
    events->SetBranchAddress("tri_ein", &ein);
    events->SetBranchAddress("tri_parent", &parent);
    events->SetBranchAddress("tmu_proc", &proc);
    events->SetBranchAddress("tmu_parent", &muParent);
    events->SetBranchAddress("tmu_st", &st);
    events->SetBranchAddress("tmu_top", &top);
    events->SetBranchAddress("tmu_mid", &mid);
    events->SetBranchAddress("prim_st", &primSt);
    events->SetBranchAddress("prim_top", &primTop);
    events->SetBranchAddress("prim_mid", &primMid);
    const double leadTop = 21. + lead;
    for (Long64_t i = 0; i < events->GetEntries(); ++i) {
      events->GetEntry(i);
      int primaryTrident = -1;
      for (std::size_t k = 0; k < w->size(); ++k) {
        const double wk = (*w)[k];
        maxWeight = std::max(maxWeight, wk);
        all.Add(wk);
        hAll->Fill((*ein)[k] / 1000., wk);
        const double r = std::max(std::abs((*x)[k]), std::abs((*y)[k]));
        const double zk = (*z)[k];
        if (r < 45. && zk > 21. && zk < leadTop) inLead.Add(wk);
        else if (r >= 45. && r <= 55. && zk > -70. && zk < leadTop) inWalls.Add(wk);
        else if (r < 40. && zk > 20. && zk < 21.) inSteel.Add(wk);
        else elsewhere.Add(wk);
        if ((*parent)[k] == 1 && primaryTrident < 0) primaryTrident = int(k);
      }
      // detectable: the primary and its produced pair all cross the top (bit 0) and middle (bit 1) TOF
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
      detectable.Add(wk);
      hDetectable->Fill((*ein)[primaryTrident] / 1000., wk);
      detectableEnergies.push_back({(*ein)[primaryTrident] / 1000., wk});
      const int separate = int(std::min(topBars.size(), midBars.size()));
      hBars->Fill(std::min(separate, 3), wk);
      if (separate >= 2) twoBars.Add(wk);
      if (separate >= 3) threeBars.Add(wk);
    }
  }
  if (liveTime <= 0. || all.n == 0) {
    printf("no tridents in %s\n", dir);
    return;
  }

  const double perDay = 86400. / liveTime;
  const double biasFactor = maxWeight > 0. ? 1. / maxWeight : 1.;
  printf("\n=== Muon tridents in the TOF-meter, lead layer %g cm (Geant4 muToMuonPairProd, biased) ===\n", lead);
  printf("%lld cosmic muons above 0.85 GeV through the 2 x 2 m plane = %.2f h of sea-level exposure (EXPACS),\n"
         "%ld biased tridents in %zu files\n", generated, liveTime / 3600., all.n, files.size());
  printf("bias factor F = %.3g: tridents were simulated F times more often, each with weight 1/F. The rates\n"
         "below are real rates, corrected for F; the simulated counts in brackets are F times the real numbers.\n",
         biasFactor);
  printf("%-62s %18s %12s %10s\n", "", "per day", "per 30 days", "per year");
  auto row = [&](const char* label, const WeightSum& s) {
    const double r = s.w * perDay, e = std::sqrt(s.w2) * perDay;
    printf("%-62s %9.4f +- %6.4f %12.3f %10.1f   (%ld simulated)\n", label, r, e, 30. * r, 365. * r, s.n);
  };
  row("all tridents in the setup", all);
  row("  in the lead layer", inLead);
  row("  in the lead walls", inWalls);
  row("  in the steel plate", inSteel);
  row("  elsewhere (TOFs, Al plates, beams, stack, air)", elsewhere);
  row("detectable: all three muons cross the top and middle TOF", detectable);
  row("  ... in at least 2 different bars of each TOF", twoBars);
  row("  ... in 3 different bars of each TOF (three separate tracks)", threeBars);
  if (!detectableEnergies.empty()) {
    std::sort(detectableEnergies.begin(), detectableEnergies.end());
    double cum = 0., median = 0.;
    for (const auto& [e, wk] : detectableEnergies) {
      cum += wk;
      if (cum >= 0.5 * detectable.w) {
        median = e;
        break;
      }
    }
    printf("median energy of the muons with a detectable trident: %.1f GeV\n", median);
  }
  printf("Rates use Geant4's trident model: total rate checked against Russell et al. at 10.5 GeV\n"
         "(analysis/trident_validation.C), but pairs too hard and opening angles too small, so the bar\n"
         "separation is pessimistic and the detectable fraction uncertain.\n");

  // plots: rates per day per energy bin
  gStyle->SetOptStat(0);
  for (auto h : {hAll, hDetectable, hBars}) h->Scale(perDay);
  TCanvas c("c", "trident rate", 1400, 720);
  auto padSpectrum = new TPad("padSpectrum", "", 0., 0.2, 0.5, 1.);
  auto padBars = new TPad("padBars", "", 0.5, 0.2, 1., 1.);
  auto padText = new TPad("padText", "", 0., 0., 1., 0.2);
  for (auto pad : {padSpectrum, padBars, padText}) pad->Draw();
  padSpectrum->cd();
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLeftMargin(0.13);
  hAll->SetTitle(Form("tridents in the setup, %g cm lead;incoming muon energy (GeV);tridents per day (corrected for F)", lead));
  hAll->SetLineColor(kBlack);
  hAll->SetLineWidth(2);
  hAll->SetMaximum(6. * hAll->GetMaximum());  // room for the legend
  hAll->Draw("HIST E1");
  hDetectable->SetLineColor(kRed + 1);
  hDetectable->SetLineWidth(2);
  hDetectable->Draw("HIST E1 SAME");
  auto legend = new TLegend(0.5, 0.78, 0.89, 0.89);
  legend->AddEntry(hAll, "all tridents", "l");
  legend->AddEntry(hDetectable, "detectable (3 muons, top & middle TOF)", "l");
  legend->Draw();
  padBars->cd();
  gPad->SetLeftMargin(0.14);
  hBars->SetTitle("detectable tridents: different bars hit (fewer of the two TOFs);bars;tridents per day (corrected for F)");
  hBars->GetXaxis()->SetNdivisions(3);
  hBars->SetFillColor(kAzure - 9);
  hBars->Draw("HIST");
  // totals and the bias factor, so that simulated counts are not taken for real ones
  padText->cd();
  auto box = new TPaveText(0.01, 0.06, 0.99, 0.94, "NDC");
  box->SetFillColor(kWhite);
  box->SetBorderSize(1);
  box->SetTextAlign(12);
  box->SetTextSize(0.16);
  box->AddText(Form("%lld muons simulated (above 0.85 GeV) = %.2f h of sea-level exposure; %ld tridents simulated, %ld of them detectable",
                    generated, liveTime / 3600., all.n, detectable.n));
  box->AddText(Form("Bias factor F = %.3g: tridents were made F times more often and each carries weight 1/F. Every rate shown is corrected for F.",
                    biasFactor));
  box->AddText(Form("Real rates: %.3g tridents per day in the setup, %.3g per day detectable (%.1f per year).",
                    all.w * perDay, detectable.w * perDay, 365. * detectable.w * perDay));
  box->Draw();
  const TString base = Form("%s/trident_rate", dir);
  c.Print(base + ".pdf");
  c.SaveAs(base + ".png");
  printf("plots: %s.pdf, %s.png\n", base.Data(), base.Data());
}
