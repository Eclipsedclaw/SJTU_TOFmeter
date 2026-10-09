// Where the primary muons stop, and the prompt camera signal for stops in the target.
//   root -l -b -q 'stops_summary.C("output/stopping_example.root")'
// Rates use the live time in the run tree (cosmic source, /GeneralParticleSource 2).
// "Prompt" = camera deposits within +-window ns of the stop time by particles other
// than the primary: muonic X-rays are emitted at the stop, capture products later
// (mu- lifetime ~200 ns in Fe, ~75 ns in Pb). Energies are true deposits.

void stops_summary(const char* fname = "output/stopping_example.root",
                   const char* target = "Target", double window = 10.)
{
  TFile file(fname);
  auto events = static_cast<TTree*>(file.Get("events"));
  auto run = static_cast<TTree*>(file.Get("run"));
  if (!events || !run) { printf("no events/run trees in %s\n", fname); return; }

  double liveTime = 0., lt = 0.;
  int nRun = 0;
  Long64_t nGenerated = 0;
  run->SetBranchAddress("live_time", &lt);
  run->SetBranchAddress("n_events", &nRun);
  for (Long64_t i = 0; i < run->GetEntries(); ++i) { run->GetEntry(i); liveTime += lt; nGenerated += nRun; }
  const double hours = liveTime / 3600.;
  printf("\n=== %s: %lld muons generated, live time %.2f h ===\n", fname, nGenerated, hours);

  int pdg = 0;
  double endT = 0., endKE = 0.;
  char endProc[64], endVol[64];
  std::vector<int>* camTrk = nullptr;
  std::vector<double>* camE = nullptr;
  std::vector<double>* camT = nullptr;
  std::vector<int>* camLayer = nullptr;
  events->SetBranchAddress("pdg", &pdg);
  events->SetBranchAddress("end_t", &endT);
  events->SetBranchAddress("end_ke", &endKE);
  events->SetBranchAddress("end_proc", endProc);
  events->SetBranchAddress("end_vol", endVol);
  events->SetBranchAddress("cam_trk", &camTrk);
  events->SetBranchAddress("cam_e", &camE);
  events->SetBranchAddress("cam_t", &camT);
  events->SetBranchAddress("cam_layer", &camLayer);

  std::map<std::string, long> stops;
  TH1D prompt("prompt", "prompt camera energy, stops in target;E (MeV);events", 400, 0., 8.);
  long nTargetStops = 0, nWithPrompt = 0;
  for (Long64_t i = 0; i < events->GetEntries(); ++i) {
    events->GetEntry(i);
    if (endKE > 0. || std::string(endVol) == "OutOfWorld") continue;  // did not stop
    stops[std::string(endVol) + "  (" + endProc + ")"]++;
    if (std::string(endVol) != target) continue;
    ++nTargetStops;
    double e = 0.;
    for (size_t k = 0; k < camE->size(); ++k)
      if (camTrk->at(k) != 1 && std::fabs(camT->at(k) - endT) < window) e += camE->at(k);
    if (e > 0.) { ++nWithPrompt; prompt.Fill(e); }
  }

  printf("\nprimary stopped in               events   per hour\n");
  for (const auto& [where, n] : stops) printf("  %-36s %7ld   %8.2f\n", where.c_str(), n, n / hours);

  printf("\nstops in %s: %ld (%.2f /h, %.1f /day)\n", target, nTargetStops, nTargetStops / hours, nTargetStops / hours * 24.);
  printf("  with prompt camera energy (|t - t_stop| < %.0f ns): %ld (%.2f /day)\n", window, nWithPrompt,
         nWithPrompt / hours * 24.);
  // Muonic K-alpha as produced by Geant4 11.4 (G4EmCaptureCascade) and as measured
  struct Line { const char* name; double g4, measured; };
  const Line lines[] = {{"Al", 0.3347, 0.3468}, {"Fe", 1.2337, 1.2547}, {"Cu", 1.4908, 1.5122}, {"Pb", 5.6739, 5.902}};
  printf("  prompt deposits near the 2p-1s lines (+-30 keV around the Geant4 energy):\n");
  for (const auto& l : lines) {
    const int b1 = prompt.FindBin(l.g4 - 0.03), b2 = prompt.FindBin(l.g4 + 0.03);
    printf("    %-2s  Geant4 %.4f MeV (measured ~%.3f): %4.0f events\n", l.name, l.g4, l.measured,
           prompt.Integral(b1, b2));
  }
  TString out(fname);
  out.ReplaceAll(".root", "_prompt.root");
  prompt.SaveAs(out);
}
