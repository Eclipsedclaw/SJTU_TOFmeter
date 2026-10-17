// Muon tridents in lead from macro/trident_validation.mac against Russell et al., Phys. Rev. Lett.
// 26 (1971) 46: 10.5 GeV/c muons on 114.5 g/cm2 of lead, 4.428e9 incident muons.
//   root -l -b -q analysis/trident_validation.C
//   root -l -b -q 'analysis/trident_validation.C("output/trident_validation.root")'
// Tridents are biased (tri_w = 1 / biasing factor), so sums of weights are true numbers.
// - effective total cross section: sum of weights / (muons x nuclei per cm2 of the target)
// - the paper's energy cut: all three final-state muons (the incoming one after the interaction and
//   the produced pair) with total energy above 1.5 GeV when they leave the target (vacuum around
//   it, so the end-of-track energy is the exit energy); the paper quotes 16% of the total for it
// - expected events in the paper's run: 4.428e9 x tridents per muon x fraction above 1.5 GeV x 82%
//   (their geometrical efficiency for these events, footnote 15), against their 75.2 +- 10.5
//   observed, 88.9 +- 2.6 predicted with Fermi statistics and 120.9 +- 2.6 without the exchange
//   term (Born approximation 101.6 and 138.4); Geant4's model has no exchange term
// - invariant mass of the two like-sign muons (the paper's Fig. 2), in the paper's run
// Writes <input>.pdf and .png.

void trident_validation(const char* fileName = "output/trident_validation.root")
{
  const double muonMass = 105.658;               // MeV
  const double targetGrams = 114.5;              // g/cm2 of lead
  const double nucleiPerCm2 = targetGrams * 6.02214e23 / 207.2;
  const double paperMuons = 4.428e9;
  const double paperGeometric = 0.82;
  const double energyCut = 1500.;                // MeV, total energy of each final-state muon

  TFile file(fileName);
  auto run = dynamic_cast<TTree*>(file.Get("run"));
  auto events = dynamic_cast<TTree*>(file.Get("events"));
  if (!run || !events || !events->GetBranch("tri_w")) {
    printf("%s: no trident data (run macro/trident_validation.mac with --preinit macro/physics_tridents.mac)\n",
           fileName);
    return;
  }
  Long64_t generated = 0;
  {
    int n = 0;
    run->SetBranchAddress("n_events", &n);
    for (Long64_t i = 0; i < run->GetEntries(); ++i) {
      run->GetEntry(i);
      generated += n;
    }
  }

  std::vector<double>*triW = nullptr, *triEin = nullptr, *triEout = nullptr, *triEpair = nullptr, *triMll = nullptr,
                       *triZ = nullptr, *tmuE0 = nullptr, *tmuE1 = nullptr;
  std::vector<int>*triPdg = nullptr, *tmuPdg = nullptr, *tmuProc = nullptr;
  int pdg = 0;
  double endKe = 0.;
  events->SetBranchAddress("tri_w", &triW);
  events->SetBranchAddress("tri_ein", &triEin);
  events->SetBranchAddress("tri_eout", &triEout);
  events->SetBranchAddress("tri_epair", &triEpair);
  events->SetBranchAddress("tri_mll", &triMll);
  events->SetBranchAddress("tri_z", &triZ);
  events->SetBranchAddress("tri_pdg", &triPdg);
  events->SetBranchAddress("tmu_pdg", &tmuPdg);
  events->SetBranchAddress("tmu_proc", &tmuProc);
  events->SetBranchAddress("tmu_e0", &tmuE0);
  events->SetBranchAddress("tmu_e1", &tmuE1);
  events->SetBranchAddress("pdg", &pdg);
  events->SetBranchAddress("end_ke", &endKe);

  const double toPaper = paperMuons / generated;  // scale a weight sum of this run to the paper's run
  auto hMll = new TH1D("hMll", "", 8, 211., 611.);
  auto hFrac = new TH1D("hFrac", "", 40, -3., 0.);
  auto hDepth = new TH1D("hDepth", "", 20, -5.044, 5.044);
  auto hEnergy = new TH1D("hEnergy", "", 50, 0., 10.5);
  for (auto h : {hMll, hFrac, hDepth, hEnergy}) h->Sumw2();

  long nTridents = 0, nSelected = 0, nPseudo = 0, nMultiple = 0;
  double sumW = 0., sumW2 = 0., selW = 0., selW2 = 0.;
  for (Long64_t i = 0; i < events->GetEntries(); ++i) {
    events->GetEntry(i);
    for (std::size_t k = 0; k < tmuProc->size(); ++k)
      if ((*tmuProc)[k] == 15) ++nPseudo;
    const std::size_t n = triW->size();
    if (n == 0) continue;
    if (n > 1) ++nMultiple;
    for (std::size_t k = 0; k < n; ++k) {
      const double w = (*triW)[k];
      ++nTridents;
      sumW += w;
      sumW2 += w * w;
      hFrac->Fill(std::log10((*triEpair)[k] / (*triEin)[k]), w * toPaper);
      hDepth->Fill((*triZ)[k], w * toPaper);
    }
    // paper's cut, for events with one trident of the beam muon: the beam muon and the two
    // produced muons all leave the target with more than 1.5 GeV
    if (n != 1 || (*triPdg)[0] != pdg) continue;
    std::vector<double> exitEnergies = {endKe + muonMass};
    for (std::size_t k = 0; k < tmuProc->size(); ++k)
      if ((*tmuProc)[k] == 49) exitEnergies.push_back((*tmuE1)[k] + muonMass);
    if (exitEnergies.size() != 3) continue;
    for (double e : exitEnergies) hEnergy->Fill(e / 1000., (*triW)[0] * toPaper);
    if (*std::min_element(exitEnergies.begin(), exitEnergies.end()) < energyCut) continue;
    ++nSelected;
    selW += (*triW)[0];
    selW2 += (*triW)[0] * (*triW)[0];
    hMll->Fill((*triMll)[0], (*triW)[0] * toPaper * paperGeometric);
  }
  if (nTridents == 0) {
    printf("no tridents in %s\n", fileName);
    return;
  }

  const double perMuon = sumW / generated, perMuonErr = std::sqrt(sumW2) / generated;
  const double sigma = perMuon / nucleiPerCm2 / 1e-30, sigmaErr = perMuonErr / nucleiPerCm2 / 1e-30;  // ub
  const double fraction = selW / sumW;
  const double fractionErr = std::sqrt(fraction * (1. - fraction) / nTridents);
  const double inAcceptance = paperMuons * perMuon * fraction * paperGeometric;
  const double inAcceptanceErr = paperMuons * std::sqrt(selW2) / generated * paperGeometric;
  printf("\n=== Muon tridents in lead: Geant4 (muToMuonPairProd) vs Russell et al., PRL 26 (1971) 46 ===\n");
  printf("%lld beam mu- of 10.5 GeV/c on %.1f g/cm2 of lead: %ld biased tridents (weight %.3g), %ld events with more\n"
         "than one, %ld muons from photon conversions to mu+mu- (not biased)\n",
         generated, targetGrams, nTridents, sumW / nTridents, nMultiple, nPseudo);
  printf("%-58s %18s   %s\n", "", "Geant4", "paper");
  printf("%-58s %10.3g +- %.2g   QED ~0.46 (60 nb / 0.13); measured ~0.39 (51 nb / 0.13)\n",
         "total trident cross section per lead nucleus (ub)", sigma, sigmaErr);
  printf("%-58s %10.3g +- %.2g   %.3g x 4.428e9 = %.0f in total\n", "tridents per incident muon (true)", perMuon,
         perMuonErr, perMuon, paperMuons * perMuon);
  printf("%-58s %10.3f +- %.3f   0.16\n", "fraction with all three muons above 1.5 GeV", fraction, fractionErr);
  printf("%-58s %10.1f +- %.1f   observed 75.2 +- 10.5; Fermi 88.9 +- 2.6; no exchange 120.9 +- 2.6\n",
         "events in the paper's acceptance (x 82% geometric)", inAcceptance, inAcceptanceErr);
  printf("(no exchange term in Geant4: compare with the paper's no-exchange numbers, Born 138.4 before their\n"
         " -13%% higher-order correction; the paper also removed 6.7 pseudotrident events)\n");
  printf("like-sign mu mu invariant mass, events in the paper's acceptance per 50 MeV bin:\n  ");
  for (int b = 1; b <= hMll->GetNbinsX(); ++b) printf(" %.0f-%.0f: %.1f", hMll->GetBinLowEdge(b), hMll->GetBinLowEdge(b + 1), hMll->GetBinContent(b));
  printf("   (overflow %.1f)\n", hMll->GetBinContent(hMll->GetNbinsX() + 1));

  gStyle->SetOptStat(0);
  TCanvas c("c", "trident validation", 1300, 1000);
  c.Divide(2, 2);
  c.cd(1);
  hMll->SetTitle("like-sign #mu#mu mass, paper's acceptance and run;invariant mass of the two like-sign muons (MeV);events / 50 MeV");
  hMll->SetLineWidth(2);
  hMll->SetLineColor(kBlue + 1);
  hMll->Draw("HIST E1");
  auto note = new TLatex();
  note->SetNDC();
  note->SetTextSize(0.04);
  note->DrawLatex(0.4, 0.82, Form("Geant4: %.1f #pm %.1f events", inAcceptance, inAcceptanceErr));
  note->DrawLatex(0.4, 0.76, "paper: 75.2 #pm 10.5 observed");
  note->DrawLatex(0.4, 0.70, "QED 88.9 (Fermi), 120.9 (no exchange)");
  c.cd(2);
  gPad->SetLogy();
  hFrac->SetTitle("pair energy / incoming muon energy (all tridents);log_{10}(E_{pair} / E_{#mu});tridents in the paper's run");
  hFrac->Draw("HIST");
  c.cd(3);
  hEnergy->SetTitle("final-state muon energies (one-trident events);total energy at the target exit (GeV);muons in the paper's run");
  hEnergy->Draw("HIST");
  auto cut = new TLine(1.5, 0., 1.5, hEnergy->GetMaximum());
  cut->SetLineColor(kRed + 1);
  cut->SetLineStyle(2);
  cut->Draw();
  c.cd(4);
  hDepth->SetTitle("trident vertex in the target (all tridents);z (cm);tridents in the paper's run");
  hDepth->SetMinimum(0.);
  hDepth->Draw("HIST E1");
  TString base = fileName;
  base.ReplaceAll(".root", "");
  c.Print(base + ".pdf");
  c.SaveAs(base + ".png");
  printf("plots: %s.pdf, %s.png\n", base.Data(), base.Data());
}
