// Camera-alone cosmic-muon run vs. Sun et al., arXiv:2608.21216, Sec. 4.1 and Table 2.
//   root -l -b -q 'camera_validation.C("output/camera_validation.root")'
// Energies are true deposits (no resolution smearing); "corrected" multiplies by
// cos(zenith), i.e. normal incidence, like the paper's incident-angle correction.

double FitMPV(TH1D& h, double& width)
{
  const double peak = h.GetBinCenter(h.GetMaximumBin());
  TF1 f("f", "landau", 0.7 * peak, 2.0 * peak);
  f.SetParameters(h.GetMaximum(), peak, 0.07 * peak);
  h.Fit(&f, "QNR");
  width = f.GetParameter(2);
  return f.GetParameter(1);  // ROOT's Landau location parameter, as a "standard Landau fit" reports
}

void camera_validation(const char* fname = "output/camera_validation.root")
{
  TFile file(fname);
  auto events = static_cast<TTree*>(file.Get("events"));
  auto run = static_cast<TTree*>(file.Get("run"));
  if (!events || !run) { printf("no events/run trees in %s\n", fname); return; }

  // Normalisation from the run tree (the events tree may be filtered, /EventFilter 1)
  double liveTime = 0., lt = 0.;
  int nRun = 0;
  Long64_t nGenerated = 0;
  run->SetBranchAddress("live_time", &lt);
  run->SetBranchAddress("n_events", &nRun);
  for (Long64_t i = 0; i < run->GetEntries(); ++i) { run->GetEntry(i); liveTime += lt; nGenerated += nRun; }

  // Working mode: CH1 and CH2 above 1000 DAC (427.6 keV); CH0 at 500 DAC (330 keV)
  const TString trigger = "e_ch1 > 0.4276 && e_ch2 > 0.4276";
  const TString threeLayer = trigger + " && e_ch0 > 0.330";
  const Long64_t nTrigger = events->GetEntries(trigger);
  const Long64_t nThree = events->GetEntries(threeLayer);
  const Long64_t nThreePrimary = events->GetEntries(threeLayer + " && cam_prim == 7");

  printf("\n=== Camera alone, cosmic muons (%s) ===\n", fname);
  printf("generated muons %lld (%lld stored), live time %.0f s = %.2f h\n", nGenerated,
         events->GetEntries(), liveTime, liveTime / 3600.);
  printf("CH1*CH2 triggers      %8lld  -> %6.2f /min        (paper: ~10 /min)\n",
         nTrigger, nTrigger / liveTime * 60.);
  printf("three-layer muons     %8lld  -> %6.0f per 15 h     (paper: ~4000 in ~15 h)\n",
         nThree, nThree / liveTime * 15. * 3600.);
  printf("  of which the primary muon crossed all three layers: %lld\n", nThreePrimary);

  const char* name[3] = {"CH0", "CH1", "CH2"};
  const double paperMPV[3] = {4.862, 1.631, 1.563};  // Table 2, track corrected, MeV
  printf("\nMIP peak of three-layer muons (MeV)   raw MPV   corrected MPV  width/MPV   paper\n");
  for (int ch = 0; ch < 3; ++ch) {
    const double hi = (ch == 0) ? 12. : 5.;
    TH1D raw(Form("raw%d", ch), "", 240, 0., hi), cor(Form("cor%d", ch), "", 240, 0., hi);
    events->Draw(Form("e_ch%d>>raw%d", ch, ch), threeLayer, "goff");
    events->Draw(Form("e_ch%d*costh>>cor%d", ch, ch), threeLayer, "goff");
    double wRaw = 0., wCor = 0.;
    const double mRaw = FitMPV(raw, wRaw), mCor = FitMPV(cor, wCor);
    printf("  %s                                 %7.3f   %7.3f        %5.1f%%      %5.3f\n",
           name[ch], mRaw, mCor, 100. * wCor / mCor, paperMPV[ch]);
  }
  printf("(paper widths sigma/MPV ~5.4-5.6%% include detector resolution; the simulation has none)\n");
}
