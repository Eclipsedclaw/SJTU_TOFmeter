// Muonic X-rays of cosmic mu- stopped in one part of the setup, as an ideal detector would see them:
// every photon of the capture, without any detector (scripts/muonic_xray.sh, macro/muonic_xray.mac).
//   root -l -b -q 'analysis/muonic_xray.C("output/muonic_xray")'                    the Al target
//   root -l -b -q 'analysis/muonic_xray.C("output/muonic_xray", "FrameworkAlPhys")'  the four Al plates
// Arguments: directory, volume (every end_vol that starts with it), file name prefix (default mux_),
// Z and A of the element (for the names of the lines; Al 13, 27), and the X-ray range and bin width
// (default 0-500 keV in 1 keV bins). A smaller range is a zoom: one panel with the X-rays and the
// nuclear-capture gammas leaving the volume in that range, files <base>_0-<range>keV.*, e.g.
//   root -l -b -q 'analysis/muonic_xray.C("output/muonic_xray", "Target", "mux_", 13, 27, 100, 0.25)'
// Last argument yield = true: the spectra per stopped mu- instead of per day, i.e. the yield
// (muonic X-rays / stopped mu-) in each bin, so a line's height is its yield; files <base>_yield.*
// Selection: primary mu- that stopped and was captured at rest (end_proc muMinusCaptureAtRest) in the
// volume; its photons (mux_parent 1): muonic X-rays of the atomic cascade (mux_src 0) and gammas of
// the nuclear capture (mux_src 1). Two spectra each: emitted (energy at creation) and leaving the
// volume (energy when the photon first left it; absorbed photons are missing, scattered ones lower).
// Both cascade models give exact line energies: Geant4's G4EmCaptureCascade hydrogen-like levels
// (n >= 2) with a tabulated 1s level, the first-principles TOFmeterEmCaptureCascade (default,
// /tof/physics/muonCascade) pure hydrogen-like levels with the reduced mass; the macro names the
// lines n -> n' from either level scheme. Rates are per day of sea-level exposure (the live time).
// Writes <dir>/muonic_xray_<volume>[_0-<range>keV][_yield].png/.pdf and _lines.csv.

namespace {
// Binding energies of the muonic atom (MeV), index n: as in G4EmCaptureCascade (tabulated 1s level,
// hydrogen-like above) or, with hydrogenic = true, pure hydrogen-like with the reduced mass as the
// emitted energies of TOFmeterEmCaptureCascade
std::vector<double> CascadeLevels(int Z, int A, bool hydrogenic = false)
{
  const int nK = 28;
  const int listK[nK] = {1, 2, 4, 6, 8, 11, 14, 17, 18, 21, 24, 26, 29, 32, 38, 40, 41, 44, 49, 53, 55,
                         60, 65, 70, 75, 81, 85, 92};
  const double listE[nK] = {0.00275, 0.011, 0.043, 0.098, 0.173, 0.326, 0.524, 0.765, 0.853, 1.146, 1.472,
                            1.708, 2.081, 2.475, 3.323, 3.627, 3.779, 4.237, 5.016, 5.647, 5.966,
                            6.793, 7.602, 8.421, 9.249, 10.222, 10.923, 11.984};
  double eK = 0.;
  for (int i = 0; i < nK; ++i) {
    if (listK[i] == Z) { eK = listE[i]; break; }
    if (i > 0 && listK[i - 1] < Z && Z < listK[i]) {
      const double y1 = listE[i - 1] / (listK[i - 1] * listK[i - 1]), y2 = listE[i] / (listK[i] * listK[i]);
      eK = (y1 + (y2 - y1) * (Z - listK[i - 1]) / double(listK[i] - listK[i - 1])) * Z * Z;
      break;
    }
  }
  const double muMass = 105.6583755, eMass = 0.51099895;
  const double nucleus = A * 931.494 - Z * eMass;   // close enough for the reduced mass
  const double e = 13.6e-6 * Z * Z * (muMass * nucleus / (muMass + nucleus)) / eMass;
  std::vector<double> levels(15, 0.);
  levels[1] = hydrogenic ? e : eK;
  for (int n = 2; n <= 14; ++n) levels[n] = e / (n * n);
  return levels;
}

std::string LineName(int n, int n2)
{
  const char* series = "KLMNOPQ";
  const char* greek[] = {"#alpha", "#beta", "#gamma", "#delta", "#epsilon"};
  if (n2 <= 7 && n - n2 <= 5) return Form("%c_{%s}", series[n2 - 1], greek[n - n2 - 1]);
  return Form("%d#rightarrow%d", n, n2);
}

struct Line { double emitted = 0., leftFull = 0., left = 0.; };
}  // namespace

void muonic_xray(const char* dir = "output/muonic_xray", const char* volume = "Target",
                 const char* prefix = "mux_", int Z = 13, int A = 27, double xMaxKeV = 500., double binKeV = 1.,
                 bool yield = false)
{
  std::vector<std::string> files;
  if (void* dirp = gSystem->OpenDirectory(dir)) {
    while (const char* entry = gSystem->GetDirEntry(dirp)) {
      const std::string name = entry;
      if (name.rfind(prefix, 0) == 0 && name.size() > 5 && name.substr(name.size() - 5) == ".root")
        files.push_back(std::string(dir) + "/" + name);
    }
    gSystem->FreeDirectory(dirp);
  }
  if (files.empty()) {
    printf("no %s*.root in %s (run scripts/muonic_xray.sh)\n", prefix, dir);
    return;
  }
  const std::string vol = volume;
  const std::string display = vol == "Target" ? "the Al target" : vol == "FrameworkAlPhys" ? "the Al plates Al1-Al4" : vol;

  const bool zoom = xMaxKeV < 499.;
  const int nX = int(std::lround(xMaxKeV / binKeV));
  auto hXEmit = new TH1D("hXEmit", "", nX, 0., xMaxKeV);   // keV
  auto hXLeft = new TH1D("hXLeft", "", nX, 0., xMaxKeV);
  auto hGZoom = new TH1D("hGZoom", "", nX, 0., xMaxKeV);   // nuclear-capture gammas leaving, X-ray binning
  auto hGEmit = new TH1D("hGEmit", "", 400, 0., 8000.);
  auto hGLeft = new TH1D("hGLeft", "", 400, 0., 8000.);
  auto hDelay = new TH1D("hDelay", "", 100, 0., 5000.);    // ns after the stop, nuclear-capture gammas
  std::map<long, Line> lines;                              // key: energy in 0.01 keV
  std::map<std::string, long> stopsByVolume;
  double liveTime = 0.;
  Long64_t generated = 0;
  long stops = 0, xrays = 0, xraysAbove10 = 0, gammas = 0;
  for (const auto& path : files) {
    TFile file(path.c_str());
    auto run = dynamic_cast<TTree*>(file.Get("run"));
    auto events = dynamic_cast<TTree*>(file.Get("events"));
    if (!run || !events || run->GetEntries() != 1 || !events->GetBranch("mux_e")) {
      printf("skipping %s: incomplete or without the mux_* columns\n", path.c_str());
      continue;
    }
    int n = 0;
    double lt = 0.;
    run->SetBranchAddress("n_events", &n);
    run->SetBranchAddress("live_time", &lt);
    run->GetEntry(0);
    generated += n;
    liveTime += lt;
    int pdg = 0;
    double endT = 0.;
    char endProc[128] = "", endVol[128] = "";
    std::vector<double>*e = nullptr, *exitE = nullptr, *t = nullptr;
    std::vector<int>*src = nullptr, *parent = nullptr;
    events->SetBranchAddress("pdg", &pdg);
    events->SetBranchAddress("end_t", &endT);
    events->SetBranchAddress("end_proc", endProc);
    events->SetBranchAddress("end_vol", endVol);
    events->SetBranchAddress("mux_e", &e);
    events->SetBranchAddress("mux_exit", &exitE);
    events->SetBranchAddress("mux_t", &t);
    events->SetBranchAddress("mux_src", &src);
    events->SetBranchAddress("mux_parent", &parent);
    for (Long64_t i = 0; i < events->GetEntries(); ++i) {
      events->GetEntry(i);
      if (pdg != 13 || std::string(endProc) != "muMinusCaptureAtRest" || std::string(endVol).rfind(vol, 0) != 0) continue;
      ++stops;
      ++stopsByVolume[endVol];
      for (std::size_t k = 0; k < e->size(); ++k) {
        if ((*parent)[k] != 1) continue;   // photons of another muon (rare)
        const double eKeV = (*e)[k] * 1000., exitKeV = (*exitE)[k] * 1000.;
        if ((*src)[k] == 0) {
          ++xrays;
          if (eKeV > 10.) ++xraysAbove10;
          hXEmit->Fill(eKeV);
          Line& line = lines[std::lround(eKeV * 100.)];
          line.emitted += 1.;
          if (exitKeV >= 0.) {
            hXLeft->Fill(exitKeV);
            line.left += 1.;
            if (std::abs(exitKeV - eKeV) < 1e-6) line.leftFull += 1.;
          }
        } else {
          ++gammas;
          hGEmit->Fill(eKeV);
          if (exitKeV >= 0.) {
            hGLeft->Fill(exitKeV);
            hGZoom->Fill(exitKeV);
          }
          hDelay->Fill((*t)[k] - endT);
        }
      }
    }
  }
  if (stops == 0) {
    printf("no mu- stopped in %s* in %s (exposure %.1f h)\n", volume, dir, liveTime / 3600.);
    return;
  }
  // cosmic runs: per day of exposure; with yield = true, or for runs without live time (a beam), per
  // stopped mu-: the yield (photons / stopped mu-)
  const bool cosmic = liveTime > 0.;
  const bool perStop = yield || !cosmic;
  const double perDay = perStop ? 1. / stops : 86400. / liveTime;   // the normalisation of every spectrum
  const double dayScale = cosmic ? 86400. / liveTime : 0.;
  const char* unit = perStop ? "stop" : "day";

  // names of the lines: n -> n' of the same level formula
  const auto levelsG4 = CascadeLevels(Z, A), levelsH = CascadeLevels(Z, A, true);
  auto nameOf = [&](double eKeV) -> std::string {
    for (const auto* levels : {&levelsG4, &levelsH})
      for (int n = 2; n <= 14; ++n)
        for (int n2 = 1; n2 < n; ++n2)
          if (std::abs(((*levels)[n2] - (*levels)[n]) * 1000. - eKeV) < 0.06) return LineName(n, n2);
    return "?";
  };

  printf("\n=== Muonic X-rays of mu- stopped in %s* (%s, Z = %d), no detector ===\n", volume, dir, Z);
  if (cosmic)
    printf("%lld muons, %.1f h of sea-level exposure; mu- stopped and captured: %ld = %.2f per hour = %.1f per day\n",
           generated, liveTime / 3600., stops, stops / (liveTime / 3600.), stops * dayScale);
  else
    printf("%lld particles (no live time: rates per stopped mu-); mu- stopped and captured: %ld\n", generated, stops);
  if (stopsByVolume.size() > 1)
    for (const auto& [v, k] : stopsByVolume) {
      if (cosmic) printf("   %-20s %6ld = %.1f per day\n", v.c_str(), k, k * dayScale);
      else printf("   %-20s %6ld\n", v.c_str(), k);
    }
  printf("photons per stopped mu-: %.3f muonic X-rays (%.3f above 10 keV), %.3f nuclear-capture gammas\n",
         double(xrays) / stops, double(xraysAbove10) / stops, double(gammas) / stops);
  const double lineMin = zoom ? 1. : 10.;              // lowest line energy in the table (keV)
  const double printMin = zoom ? 0.01 : 0.02;          // printed from this yield per stop
  printf("lines of %.0f-%.0f keV with at least %.0f%% per stop (all from 0.5%% in the CSV):\n", lineMin, xMaxKeV,
         100. * printMin);
  if (perStop)
    printf("%10s %8s %16s %26s %20s\n", "line (keV)", "n -> n'", "yield (emitted)", "yield leaving unscattered",
           "yield leaving (all)");
  else
    printf("%10s %8s %10s %12s %16s %12s\n", "line (keV)", "n -> n'", "per stop", "emitted/day", "left unscattered",
           "left/day");
  TString base = zoom ? Form("%s/muonic_xray_%s_0-%gkeV", dir, volume, xMaxKeV) : Form("%s/muonic_xray_%s", dir, volume);
  if (yield && cosmic) base += "_yield";
  std::ofstream csv((base + "_lines.csv").Data());
  csv << Form("energy_keV,transition,per_stop,emitted_per_%s,left_unscattered_per_%s,left_per_%s\n", unit, unit, unit);
  std::vector<std::pair<double, double>> labelled;   // energy, emitted per day, for the plot labels
  for (const auto& [key, line] : lines) {
    const double eKeV = key / 100.;
    if (line.emitted / stops < 0.005 || eKeV < lineMin || eKeV > xMaxKeV) continue;   // at least 0.5% per stop
    const bool print = line.emitted / stops >= printMin;
    std::string name = nameOf(eKeV);
    TString plain = name.c_str();
    plain.ReplaceAll("#alpha", "a"); plain.ReplaceAll("#beta", "b"); plain.ReplaceAll("#gamma", "g");
    plain.ReplaceAll("#delta", "d"); plain.ReplaceAll("#epsilon", "e"); plain.ReplaceAll("#rightarrow", "->");
    plain.ReplaceAll("_{", ""); plain.ReplaceAll("}", "");
    if (print && perStop)
      printf("%10.1f %8s %16.3f %18.3f (%3.0f%%) %20.3f\n", eKeV, plain.Data(), line.emitted / stops,
             line.leftFull / stops, 100. * line.leftFull / line.emitted, line.left / stops);
    else if (print)
      printf("%10.1f %8s %10.3f %12.1f %9.1f (%3.0f%%) %12.1f\n", eKeV, plain.Data(), line.emitted / stops,
             line.emitted * perDay, line.leftFull * perDay, 100. * line.leftFull / line.emitted, line.left * perDay);
    csv << eKeV << "," << plain << "," << line.emitted / stops << "," << line.emitted * perDay << ","
        << line.leftFull * perDay << "," << line.left * perDay << "\n";
    if (eKeV >= (zoom ? 5. : 20.) && line.emitted / stops >= (zoom ? 0.02 : 0.03)) labelled.push_back({eKeV, line.emitted * perDay});
  }
  if (Z == 13) printf("measured Al K_alpha (2p -> 1s): 346.8 keV, about 0.8 per stopped mu-; models: Geant4 334.7 keV, "
                      "first-principles 354.9 keV\n");
  printf("nuclear-capture gammas: %.2f per %s emitted, %.2f per %s leaving the volume\n", hGEmit->GetEntries() * perDay, unit,
         hGLeft->GetEntries() * perDay, unit);

  // ---- plots: rates per day per keV (X-rays) and per day per 20 keV (gammas)
  for (auto h : {hXEmit, hXLeft, hGEmit, hGLeft, hGZoom}) h->Scale(perDay);
  const int blue = TColor::GetColor("#2a78d6"), orange = TColor::GetColor("#eb6834");
  const int inkSecondary = TColor::GetColor("#52514e");
  gStyle->SetOptStat(0);
  TCanvas c("c", "muonic X-rays", 1400, 1150);
  auto pX = new TPad("pX", "", 0., zoom ? 0.13 : 0.47, 1., 1.);   // a zoom has no gamma panel
  auto pG = new TPad("pG", "", 0., 0.13, 1., 0.47);
  auto pT = new TPad("pT", "", 0., 0., 1., 0.13);
  pX->Draw();
  if (!zoom) pG->Draw();
  pT->Draw();
  auto styleHist = [&](TH1D* h, int colour, int width) {
    h->SetLineColor(colour);
    h->SetLineWidth(width);
  };
  pX->cd();
  const double k = zoom ? 0.6 : 1.;   // ROOT text scales with the pad height: the zoom pad is taller
  gPad->SetLogy(perStop ? 0 : 1);     // a yield: linear, 0 to 1
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.03);
  gPad->SetTopMargin(0.1);
  styleHist(hXEmit, blue, 2);
  styleHist(hXLeft, orange, 2);
  const TString perBin = std::abs(binKeV - 1.) < 1e-9 ? TString("keV") : TString(Form("%g keV", binKeV));
  if (perStop)
    hXEmit->SetTitle(Form(";photon energy (keV);yield = muonic X-rays / stopped #mu^{-} (per %s)", perBin.Data()));
  else
    hXEmit->SetTitle(Form(";photon energy (keV);photons per day per %s", perBin.Data()));
  if (perStop) {   // a yield is at most 1 per stopped mu-
    hXEmit->SetMinimum(0.);
    hXEmit->SetMaximum(1.);
  } else {
    hXEmit->SetMinimum(0.5 * perDay);
    hXEmit->SetMaximum(hXEmit->GetMaximum() * (zoom ? 40. : 8.));   // room for the legend
  }
  for (auto axis : {hXEmit->GetXaxis(), hXEmit->GetYaxis()}) {
    axis->SetTitleSize(0.045 * k);
    axis->SetLabelSize(0.04 * k);
    axis->SetTitleColor(inkSecondary);
  }
  hXEmit->GetYaxis()->SetTitleOffset(zoom ? 1.1 : 0.85);
  const int grey = TColor::GetColor("#8c8a83");
  hGZoom->SetLineColor(grey);
  hGZoom->SetLineWidth(2);
  hXEmit->Draw("HIST");
  if (zoom && !perStop) hGZoom->Draw("HIST SAME");   // too small to see on a linear yield axis
  hXLeft->Draw("HIST SAME");
  auto label = new TLatex();
  label->SetTextFont(42);
  label->SetTextSize(0.04 * k);
  label->SetTextColor(inkSecondary);
  label->SetTextAlign(21);
  std::sort(labelled.begin(), labelled.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
  std::vector<double> placed;
  for (const auto& entry : labelled) {
    const double eKeV = entry.first, rate = entry.second;
    const double spacing = 0.07 * xMaxKeV;   // labels at least this far apart
    if (placed.size() >= (zoom ? 8u : 6u)) break;
    if (std::any_of(placed.begin(), placed.end(), [eKeV, spacing](double x) { return std::abs(x - eKeV) < spacing; })) continue;
    placed.push_back(eKeV);
    // centred on the line, but kept inside the frame near its ends
    label->SetTextAlign(eKeV > 0.93 * xMaxKeV ? 31 : eKeV < 0.07 * xMaxKeV ? 11 : 21);
    // above the bin (several lines can share one)
    const double top = std::max(rate, hXEmit->GetBinContent(hXEmit->FindBin(eKeV)));
    label->DrawLatex(eKeV, perStop ? top + 0.02 : top * 1.6, Form("%s %.1f", nameOf(eKeV).c_str(), eKeV));
  }
  auto title = new TLatex();
  title->SetNDC();
  title->SetTextFont(62);
  title->SetTextSize(0.05 * k);
  title->DrawLatex(0.08, 0.93, Form("Muonic X-ray %s of #mu^{-} stopped in %s, ideal detector%s", perStop ? "yield" : "spectrum",
                                    display.c_str(), zoom ? Form(", 0-%g keV", xMaxKeV) : ""));
  // full range: between the L and the K lines; zoom: across the top, above the labels
  auto legend = zoom ? new TLegend(0.38, 0.79, 0.96, 0.9) : new TLegend(0.3, 0.74, 0.66, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.038 * k);
  legend->AddEntry(hXEmit, "emitted (energy at creation)", "l");
  legend->AddEntry(hXLeft, "leaving the volume (energy at its surface)", "l");
  if (zoom && !perStop) legend->AddEntry(hGZoom, "nuclear-capture gammas leaving the volume (~800 ns later)", "l");
  legend->Draw();

  if (!zoom) {   // the nuclear-capture gammas, 0-8 MeV
    pG->cd();
    gPad->SetLogy();
    gPad->SetLeftMargin(0.08);
    gPad->SetRightMargin(0.03);
    gPad->SetBottomMargin(0.17);
    styleHist(hGEmit, blue, 2);
    styleHist(hGLeft, orange, 2);
    hGEmit->SetTitle(Form(";photon energy (keV);gammas per %s per 20 keV", perStop ? "stopped #mu^{-}" : "day"));
    hGEmit->SetMinimum(0.5 * perDay);
    hGEmit->SetMaximum(std::max(10. * perDay, hGEmit->GetMaximum() * 5.));
    for (auto axis : {hGEmit->GetXaxis(), hGEmit->GetYaxis()}) {
      axis->SetTitleSize(0.07);
      axis->SetLabelSize(0.06);
      axis->SetTitleColor(inkSecondary);
  }
  hGEmit->GetYaxis()->SetTitleOffset(0.55);
  hGEmit->Draw("HIST");
  hGLeft->Draw("HIST SAME");
  auto gTitle = new TLatex();
  gTitle->SetNDC();
  gTitle->SetTextFont(42);
  gTitle->SetTextSize(0.07);
  gTitle->SetTextColor(inkSecondary);
  gTitle->DrawLatex(0.42, 0.82, Form("nuclear-capture gammas: %.3g per %s, delayed (mean %.0f ns)",
                                     hGEmit->GetEntries() > 0 ? hGEmit->Integral() : 0., unit, hDelay->GetMean()));
  }

  pT->cd();
  auto box = new TPaveText(0.01, 0.05, 0.99, 0.95, "NDC");
  box->SetFillColor(kWhite);
  box->SetBorderSize(1);
  box->SetTextAlign(12);
  box->SetTextFont(42);
  box->SetTextSize(zoom ? (perStop ? 0.13 : 0.16) : 0.19);
  box->SetTextColor(inkSecondary);
  if (perStop && cosmic)
    box->AddText(Form("%.0f h of sea-level exposure: %ld #mu^{-} stopped (%.0f per day); %.2f X-rays per stop above 10 keV. "
                      "A line's height is its yield.", liveTime / 3600., stops, stops * dayScale, double(xraysAbove10) / stops));
  else
    box->AddText(cosmic ? Form("%.0f h of sea-level exposure: %ld #mu^{-} stopped = %.1f per day; %.2f X-rays per stop above 10 keV",
                               liveTime / 3600., stops, stops * dayScale, double(xraysAbove10) / stops)
                        : Form("%ld #mu^{-} stopped; %.2f X-rays per stop above 10 keV", stops, double(xraysAbove10) / stops));
  if (zoom)
    box->AddText(Form("%g keV bins. Lines below 10 keV (0.4, 0.5, 5.8, 10.7 keV, ...) are absorbed where they are made.", binKeV));
  if (zoom && perStop)
    box->AddText(Form("Not shown: nuclear-capture gammas below %g keV leaving the volume, %.3f per stop (~800 ns later).",
                      xMaxKeV, hGZoom->Integral()));
  if (zoom)
    box->AddText("Exact model lines (hydrogen-like L and M energies), no detector, no resolution");
  else
    box->AddText(Z == 13 ? "Exact model lines, no detector. Al K_{#alpha}: 346.8 keV measured; 354.9 (first-principles), 334.7 (Geant4 model)"
                         : "Exact model lines (hydrogen-like levels, a few % off or more), no detector");
  box->Draw();
  c.Print(base + ".pdf");
  c.SaveAs(base + ".png");
  printf("plots: %s.png, .pdf; lines: %s_lines.csv\n", base.Data(), base.Data());
}
