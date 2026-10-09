// Cosmic-muon acceptance of the setup, from the two runs of macro/acceptance.mac:
//  - cosmic run (sea-level spectrum, both charges, full physics): rate of each selection;
//    muons that stop in the lead layer or scatter away are lost as in reality
//  - straight-line run (10 GeV mu-, isotropic cosine-law flux, no scattering or radiative
//    processes): geometric acceptance A*Omega = N_sel / N_gen * A_plane * pi sin^2(theta_max)
// A TOF counts as hit when one of its bars has more than barThr (MeV); the Compton detector
// when one of its active layers (stack YSO or camera CH0-2) has more than comptonThr.
//   root -l -b -q 'analysis/acceptance.C(600, 80)'   (straight-line plane half size cm, theta max deg)

struct RunInfo {
  Long64_t generated = 0;
  double liveTime = 0.;  // s, cosmic runs only
};

RunInfo ReadRun(TFile& file)
{
  RunInfo info;
  auto run = static_cast<TTree*>(file.Get("run"));
  int n = 0;
  double lt = 0.;
  run->SetBranchAddress("n_events", &n);
  run->SetBranchAddress("live_time", &lt);
  for (Long64_t i = 0; i < run->GetEntries(); ++i) {
    run->GetEntry(i);
    info.generated += n;
    info.liveTime += lt;
  }
  return info;
}

void acceptance(double geoHalfCm = 600., double geoThetaMaxDeg = 80.,
                const char* cosmicFile = "output/acceptance_cosmic.root",
                const char* geoFile = "output/acceptance_geometric.root",
                double barThr = 0.2, double comptonThr = 0.1)
{
  TFile fc(cosmicFile), fg(geoFile);
  auto cosmic = static_cast<TTree*>(fc.Get("events"));
  auto straight = static_cast<TTree*>(fg.Get("events"));
  if (!cosmic || !straight || !fc.Get("run") || !fg.Get("run")) {
    printf("missing trees in %s or %s\n", cosmicFile, geoFile);
    return;
  }
  const RunInfo rc = ReadRun(fc), rg = ReadRun(fg);
  const double sin2 = std::pow(std::sin(geoThetaMaxDeg * TMath::DegToRad()), 2);
  const double perEvent = 4. * geoHalfCm * geoHalfCm * TMath::Pi() * sin2 / rg.generated;  // cm2 sr

  const TString top = Form("Sum$(bar_st==0 && bar_e>%g)>0", barThr);
  const TString middle = Form("Sum$(bar_st==1 && bar_e>%g)>0", barThr);
  const TString bottom = Form("Sum$(bar_st==2 && bar_e>%g)>0", barThr);
  const TString compton = Form("(e_stack0>%g || e_stack1>%g || e_ch0>%g || e_ch1>%g || e_ch2>%g)",
                               comptonThr, comptonThr, comptonThr, comptonThr, comptonThr);
  const TString tm = top + " && " + middle;
  struct Selection { const char* name; TString cut; };
  const std::vector<Selection> selections = {
    {"lead layer", "e_absorber > 0"},
    {"top TOF", top},
    {"middle TOF", middle},
    {"bottom TOF", bottom},
    {"Compton detector", compton},
    {"top & middle", tm},
    {"top & middle & bottom", tm + " && " + bottom},
    {"top & middle & Compton", tm + " && " + compton},
    {"top & middle & Compton & bottom", tm + " && " + compton + " && " + bottom},
    {"top & middle & not bottom", tm + " && !(" + bottom + ")"},
  };

  printf("\n=== Cosmic-muon acceptance ===\n");
  printf("cosmic run:   %lld muons = %.0f s of live time (%s)\n", rc.generated, rc.liveTime, cosmicFile);
  printf("straight run: %lld muons, plane %.0f x %.0f cm, zenith < %.0f deg (%s)\n", rg.generated,
         2 * geoHalfCm, 2 * geoHalfCm, geoThetaMaxDeg, geoFile);
  printf("hit thresholds: TOF bar > %.2f MeV, Compton layer > %.2f MeV\n\n", barThr, comptonThr);
  printf("%-34s %9s %18s %10s %20s\n", "selection", "muons", "rate /min", "per day", "A*Omega (cm2 sr)");
  for (const auto& sel : selections) {
    const Long64_t nc = cosmic->GetEntries(sel.cut);
    const Long64_t ng = straight->GetEntries(sel.cut);
    const double rate = rc.liveTime > 0 ? nc / rc.liveTime * 60. : 0.;
    const double rateErr = rc.liveTime > 0 ? std::sqrt(double(nc)) / rc.liveTime * 60. : 0.;
    printf("%-34s %9lld %9.2f +- %6.2f %10.0f %11.1f +- %6.1f\n", sel.name, nc, rate, rateErr,
           rate * 60. * 24., ng * perEvent, std::sqrt(double(ng)) * perEvent);
  }
  printf("\nRates: realistic sea-level muons (statistical errors only). A*Omega: straight lines;\n"
         "\"not bottom\" counts muons that stop or scatter out before the bottom TOF in the rates,\n"
         "and tracks that miss the bottom TOF geometrically in A*Omega.\n");
}
