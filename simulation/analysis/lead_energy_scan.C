// Top & middle & Compton coincidences versus muon energy and top-lead thickness, from a scan of
// macro/lead_energy_scan_t.mac in one directory: leadscan_t<cm>_E<GeV>.root (one per point) and
// the spectrum check runs leadspec_t<cm>.root. The energies, their bands and the thicknesses are
// read from the files, so the same macro serves the log scan (lead_energy_scan.mac, 0.1-100 GeV)
// and the fine low-energy scan (lead_lowE_scan.mac, 100-500 MeV in 20 MeV steps).
//   root -l -b -q analysis/lead_energy_scan.C
//   root -l -b -q 'analysis/lead_energy_scan.C("output/leadscan_lowE")'
//   root -l -b -q 'analysis/lead_energy_scan.C("output/leadscan", 0.2, 0.1, 24)'   (24 h exposure)
//
// For every point (lead thickness t, muon kinetic energy E):
//   ratio  N(TMC) / N(generated), the fraction of the muons of energy E (sea-level zenith
//          distribution, spread over the 2 x 2 m source plane) that fire top & middle & Compton
//   A_eff  ratio x plane area (cm2); the rate is A_eff times the muon flux per cm2
//   rate   N(TMC) / live time: the rate due to the sea-level muons of the energy band (EXPACS)
//          that E stands for. Summed over the points it is the rate of the muons of the whole
//          band range, which the spectrum check run of the same thickness measures directly.
// A TOF counts as hit when one of its bars has more than barThr (MeV), the Compton detector
// when a stack YSO or camera layer has more than comptonThr.
// TMC+stop: TMC events whose muon stopped in the Compton region (-52.5 < z < -32.5 cm inside the
// frame, |x|, |y| < 40 cm: the stack, its plate Al3 and the space above it up to Al2). For one
// energy the stopping window is narrow (~10 MeV for one zenith angle), so the sum over the points
// gives its rate only if the grid is fine enough; the check run gives it in any case.
// Flux calibration: every point is weighted with w = J_band / J_max, its band flux relative to
// the band with the most muons (w = 1 there). On a log grid a point stands for a band 0.36 E wide,
// which is why the band flux, not dJ/dE, is the weight. ratio x w is the flux-weighted ratio,
// and ratio x w x N_ref, with N_ref = J_max x plane area x exposure the muons of that band
// crossing the plane, the expected number of TMC muons in the exposure (= rate x exposure).
// Writes <dir>/lead_energy_scan.csv (every point) and <dir>/lead_energy_scan.pdf (maps, totals
// per thickness, ratio versus thickness per energy, flux calibration) with a .png per page.

namespace {
enum { kTop, kTM, kTMC, kTMCnoB, kTMCstop, kNSel };
const char* kSelName[kNSel] = {"top", "top_middle", "TMC", "TMC_noBottom", "TMC_stopped"};

struct Counts {
  bool filled = false;
  double flux = 0., area = 0., liveTime = 0., planeZ = 0.;  // /cm2/s through the plane, cm2, s, cm
  double eMono = 0., eMin = 0., eMax = 0., lead = 0.;       // MeV, cm
  Long64_t generated = 0;
  Long64_t n[kNSel] = {};
  double Ratio(int s) const { return generated > 0 ? double(n[s]) / generated : 0.; }
  double RatioErr(int s) const { return generated > 0 ? std::sqrt(double(n[s])) / generated : 0.; }
  double Rate(int s) const { return liveTime > 0 ? n[s] / liveTime * 60. : 0.; }  // per minute
  double RateErr(int s) const { return liveTime > 0 ? std::sqrt(double(n[s])) / liveTime * 60. : 0.; }
};

// Run summary and selection counts of one file; false if the file is incomplete
bool ReadCounts(const TString& path, const TString* cuts, Counts& c)
{
  TFile file(path);
  auto run = dynamic_cast<TTree*>(file.Get("run"));
  auto events = dynamic_cast<TTree*>(file.Get("events"));
  if (!run || !events || !run->GetBranch("e_mono") || !run->GetBranch("lead_cm") || run->GetEntries() != 1)
    return false;
  int nEvents = 0;
  run->SetBranchAddress("n_events", &nEvents);
  run->SetBranchAddress("flux", &c.flux);
  run->SetBranchAddress("area", &c.area);
  run->SetBranchAddress("live_time", &c.liveTime);
  run->SetBranchAddress("plane_z", &c.planeZ);
  run->SetBranchAddress("e_min", &c.eMin);
  run->SetBranchAddress("e_max", &c.eMax);
  run->SetBranchAddress("e_mono", &c.eMono);
  run->SetBranchAddress("lead_cm", &c.lead);
  run->GetEntry(0);
  c.generated = nEvents;
  for (int s = 0; s < kNSel; ++s) c.n[s] = events->GetEntries(cuts[s]);
  c.filled = true;
  return true;
}

// Index of x in the sorted list v (relative tolerance), -1 if absent
int FindValue(const std::vector<double>& v, double x)
{
  for (std::size_t i = 0; i < v.size(); ++i)
    if (std::abs(v[i] - x) <= 1e-6 * std::max(1., std::abs(x))) return int(i);
  return -1;
}

// Bin edges around sorted centres: midpoints, half a step beyond the ends
std::vector<double> CentreEdges(const std::vector<double>& c)
{
  std::vector<double> e(c.size() + 1);
  const double first = c.size() > 1 ? c[1] - c[0] : 1., last = c.size() > 1 ? c[c.size() - 1] - c[c.size() - 2] : 1.;
  e[0] = c[0] - 0.5 * first;
  for (std::size_t i = 1; i < c.size(); ++i) e[i] = 0.5 * (c[i - 1] + c[i]);
  e[c.size()] = c.back() + 0.5 * last;
  return e;
}
}  // namespace

void lead_energy_scan(const char* dir = "output/leadscan", double barThr = 0.2, double comptonThr = 0.1,
                      double exposureHours = 1.)
{
  const TString top = Form("Sum$(bar_st==0 && bar_e>%g)>0", barThr);
  const TString middle = Form("Sum$(bar_st==1 && bar_e>%g)>0", barThr);
  const TString bottom = Form("Sum$(bar_st==2 && bar_e>%g)>0", barThr);
  const TString compton = Form("(e_stack0>%g || e_stack1>%g || e_ch0>%g || e_ch1>%g || e_ch2>%g)",
                               comptonThr, comptonThr, comptonThr, comptonThr, comptonThr);
  const TString stopped = "(end_ke < 0.001 && end_z > -52.5 && end_z < -32.5 && abs(end_x) < 40 && abs(end_y) < 40)";
  const TString tmc = top + " && " + middle + " && " + compton;
  const TString cuts[kNSel] = {top, top + " && " + middle, tmc, tmc + " && !(" + bottom + ")",
                               tmc + " && " + stopped};

  // ---- read the scan points and the spectrum check runs
  std::vector<Counts> points, checks;
  void* dirp = gSystem->OpenDirectory(dir);
  if (!dirp) {
    printf("cannot open %s\n", dir);
    return;
  }
  while (const char* entry = gSystem->GetDirEntry(dirp)) {
    const TString name = entry;
    const bool isPoint = name.BeginsWith("leadscan_t"), isSpectrum = name.BeginsWith("leadspec_t");
    if (!(isPoint || isSpectrum) || !name.EndsWith(".root")) continue;
    Counts c;
    if (!ReadCounts(Form("%s/%s", dir, entry), cuts, c)) {
      printf("skipping %s: no complete run summary\n", entry);
      continue;
    }
    if (isSpectrum != (c.eMono == 0.)) {
      printf("skipping %s: %s\n", entry, isSpectrum ? "not a spectrum run" : "not mono-energetic");
      continue;
    }
    (isSpectrum ? checks : points).push_back(c);
  }
  gSystem->FreeDirectory(dirp);
  if (points.empty()) {
    printf("no scan points in %s\n", dir);
    return;
  }

  // ---- the grid: thicknesses (cm) and energies (MeV) found in the files
  std::vector<double> leads, energies;
  for (const auto& p : points) {
    if (FindValue(leads, p.lead) < 0) leads.push_back(p.lead);
    if (FindValue(energies, p.eMono) < 0) energies.push_back(p.eMono);
  }
  std::sort(leads.begin(), leads.end());
  std::sort(energies.begin(), energies.end());
  const int nT = int(leads.size()), nE = int(energies.size());
  std::vector<std::vector<Counts>> grid(nT, std::vector<Counts>(nE));
  for (const auto& p : points) {
    Counts& slot = grid[FindValue(leads, p.lead)][FindValue(energies, p.eMono)];
    if (slot.filled) printf("warning: a second file for lead %g cm, E = %g MeV replaces the first\n", p.lead, p.eMono);
    slot = p;
  }
  std::vector<Counts> spectrum(nT);
  for (const auto& c : checks) {
    const int it = FindValue(leads, c.lead);
    if (it >= 0) spectrum[it] = c;
    else printf("spectrum check run for %g cm has no scan points; ignored\n", c.lead);
  }
  int missing = 0;
  for (auto& row : grid) for (auto& p : row) if (!p.filled) ++missing;

  // bands of each energy (MeV), from any thickness
  std::vector<double> bandLow(nE, 0.), bandHigh(nE, 0.), bandFlux(nE, 0.);
  for (int ie = 0; ie < nE; ++ie)
    for (int it = 0; it < nT; ++it)
      if (grid[it][ie].filled) {
        bandLow[ie] = grid[it][ie].eMin;
        bandHigh[ie] = grid[it][ie].eMax;
        bandFlux[ie] = grid[it][ie].flux;
        break;
      }
  bool equalBands = true, contiguous = true;
  for (int ie = 1; ie < nE; ++ie) {
    const double w0 = bandHigh[0] - bandLow[0], w = bandHigh[ie] - bandLow[ie];
    if (std::abs(w - w0) > 1e-3 * w0) equalBands = false;
    if (std::abs(bandLow[ie] - bandHigh[ie - 1]) > 1e-3 * bandLow[ie]) contiguous = false;
  }
  // display: MeV and a linear axis for a grid below 1 GeV of equal bands, else GeV and log axis
  const bool mev = energies.back() < 1000.;
  const bool logAxis = !equalBands;
  const char* unit = mev ? "MeV" : "GeV";
  auto show = [&](double eMeV) { return mev ? eMeV : eMeV / 1000.; };
  const TString rangeText = Form("%.3g-%.3g %s", show(bandLow[0]), show(bandHigh[nE - 1]), unit);
  auto energyText = [&](int ie) { return mev ? TString(Form("%.0f", energies[ie])) : TString(Form("%.3g", energies[ie] / 1000.)); };
  auto leadText = [&](int it) { return TString(Form("%g", leads[it])); };

  const Counts* any = nullptr;
  for (auto& row : grid) for (auto& p : row) if (p.filled) any = &p;
  printf("\n=== Lead thickness x muon energy scan (%s): %d energies x %d thicknesses, %d of %d points "
         "missing, %d spectrum check runs ===\n", dir, nE, nT, missing, nE * nT, int(checks.size()));
  printf("energies %s to %s %s in bands of %s width covering %s%s; source plane %.0f cm2 at z = %g cm, "
         "%lld muons per point\nTOF bar > %.2f MeV, Compton layer > %.2f MeV\n", energyText(0).Data(),
         energyText(nE - 1).Data(), unit, equalBands ? "equal" : "log", rangeText.Data(),
         contiguous ? "" : " (with gaps)", any->area, any->planeZ, any->generated, barThr, comptonThr);
  for (const auto& c : checks)
    if (std::abs(c.eMin - bandLow[0]) > 1e-3 * bandLow[0] || std::abs(c.eMax - bandHigh[nE - 1]) > 1e-3 * bandHigh[nE - 1])
      printf("warning: the check run for %g cm covers %g-%g MeV, not the bands of the points\n", c.lead, c.eMin, c.eMax);

  // ---- tables: rows = lead thickness, columns = energy
  auto printMatrix = [&](const char* title, auto value, const char* format) {
    printf("\n%s\nlead\\E %s", title, unit);
    for (int ie = 0; ie < nE; ++ie) printf(" %6s", energyText(ie).Data());
    printf("\n");
    for (int it = 0; it < nT; ++it) {
      printf("%5s cm  ", leadText(it).Data());
      for (int ie = 0; ie < nE; ++ie) {
        if (grid[it][ie].filled) printf(format, value(grid[it][ie]));
        else printf(" %6s", "-");
      }
      printf("\n");
    }
  };
  printMatrix("Ratio: top & middle & Compton / muons generated (x 1e-3)",
              [](const Counts& p) { return 1e3 * p.Ratio(kTMC); }, " %6.3f");
  printMatrix("Rate: top & middle & Compton due to the muons of each energy band (per minute)",
              [](const Counts& p) { return p.Rate(kTMC); }, " %6.3f");
  printMatrix("Top & middle & Compton with the muon stopped in the Compton region / muons generated (x 1e-3)",
              [](const Counts& p) { return 1e3 * p.Ratio(kTMCstop); }, " %6.3f");

  // ---- flux calibration: w = J_band / J_max and the expected number of TMC muons
  const int iMax = int(std::max_element(bandFlux.begin(), bandFlux.end()) - bandFlux.begin());
  std::vector<double> coeff(nE, 0.), perGeV(nE, 0.);
  double perGeVMax = 0.;
  for (int ie = 0; ie < nE; ++ie) {
    coeff[ie] = bandFlux[iMax] > 0. ? bandFlux[ie] / bandFlux[iMax] : 0.;
    perGeV[ie] = bandHigh[ie] > bandLow[ie] ? bandFlux[ie] / ((bandHigh[ie] - bandLow[ie]) / 1000.) : 0.;
    perGeVMax = std::max(perGeVMax, perGeV[ie]);
  }
  const double nRef = bandFlux[iMax] * any->area * exposureHours * 3600.;  // muons of the max band
  auto energyIndex = [&](const Counts& p) { return FindValue(energies, p.eMono); };
  auto weighted = [&](const Counts& p) { return p.Ratio(kTMC) * coeff[energyIndex(p)]; };
  auto expected = [&](const Counts& p) { return weighted(p) * nRef; };
  auto expectedErr = [&](const Counts& p) { return p.RatioErr(kTMC) * coeff[energyIndex(p)] * nRef; };

  printf("\nFlux coefficients w = J_band / J_max (EXPACS mu+ + mu- of each band through the plane)\n");
  printf("%9s %21s %15s %8s %14s\n", Form("E (%s)", unit), Form("band (%s)", unit), "J_band /cm2/s", "w",
         "(dJ/dE)/max");
  for (int ie = 0; ie < nE; ++ie)
    printf("%9s %9.4g - %-9.4g %15.4e %8.4f %14.4f\n", energyText(ie).Data(), show(bandLow[ie]), show(bandHigh[ie]),
           bandFlux[ie], coeff[ie], perGeVMax > 0. ? perGeV[ie] / perGeVMax : 0.);
  if (equalBands)
    printf("w = 1 at %s %s. The bands are equally wide, so w is also the spectrum dJ/dE relative to its maximum.\n",
           energyText(iMax).Data(), unit);
  else
    printf("w = 1 at %s %s. The last column, the spectrum per GeV, is for comparison only: on this log\n"
           "grid each point stands for a band 0.36 E wide, so the muons it represents scale with J_band.\n",
           energyText(iMax).Data(), unit);
  printMatrix(Form("Flux-weighted ratio: ratio x w (x 1e-3), equal to the ratio at %s %s", energyText(iMax).Data(), unit),
              [&](const Counts& p) { return 1e3 * weighted(p); }, " %6.3f");
  double maxN = 0.;
  for (auto& row : grid) for (auto& p : row) if (p.filled) maxN = std::max(maxN, expected(p));
  printMatrix(Form("Expected top & middle & Compton muons in %g h: N = ratio x w x N_ref, N_ref = %.4g muons of "
                   "the %s %s band through the plane", exposureHours, nRef, energyText(iMax).Data(), unit),
              expected, maxN < 9999.95 ? " %6.1f" : (maxN < 999999.5 ? " %6.0f" : " %6.2g"));

  // ---- totals over the energy bands, per thickness, and the spectrum checks
  std::vector<double> total(nT, 0.), totalErr(nT, 0.), totalStop(nT, 0.), totalStopErr(nT, 0.);
  for (int it = 0; it < nT; ++it) {
    for (int ie = 0; ie < nE; ++ie) {
      total[it] += grid[it][ie].Rate(kTMC);
      totalErr[it] += std::pow(grid[it][ie].RateErr(kTMC), 2);
      totalStop[it] += grid[it][ie].Rate(kTMCstop) * 60.;
      totalStopErr[it] += std::pow(grid[it][ie].RateErr(kTMCstop) * 60., 2);
    }
    totalErr[it] = std::sqrt(totalErr[it]);
    totalStopErr[it] = std::sqrt(totalStopErr[it]);
  }
  printf("\nTop & middle & Compton rate of the sea-level muons %s, per minute\n", rangeText.Data());
  printf("%7s %22s %8s %22s %14s %8s %22s\n", "lead", Form("sum of the %d bands", nE), "vs first", "spectrum check run",
         "sum/spectrum", "pull", Form("expected N in %g h", exposureHours));
  int best = 0;
  for (int it = 0; it < nT; ++it) {
    printf("%4s cm %13.2f +- %5.2f %8.3f", leadText(it).Data(), total[it], totalErr[it], total[0] > 0 ? total[it] / total[0] : 0.);
    const Counts& s = spectrum[it];
    if (s.filled && s.n[kTMC] > 0) {
      const double r = s.Rate(kTMC), e = s.RateErr(kTMC);
      printf(" %13.2f +- %5.2f %14.3f %8.2f", r, e, total[it] / r, (total[it] - r) / std::hypot(totalErr[it], e));
    } else {
      printf(" %22s %14s %8s", "", "", "");
    }
    printf(" %13.0f +- %6.0f\n", total[it] * 60. * exposureHours, totalErr[it] * 60. * exposureHours);
    if (total[it] > total[best]) best = it;
  }
  printf("highest top & middle & Compton rate: %s cm of lead\n", leadText(best).Data());

  printf("\nTop & middle & Compton with the muon stopped in the Compton region, per hour\n");
  printf("%7s %22s %22s %14s %8s\n", "lead", Form("sum of the %d bands", nE), "spectrum check run", "sum/spectrum", "pull");
  int bestStop = 0, bestStopCheck = -1;
  for (int it = 0; it < nT; ++it) {
    printf("%4s cm %13.2f +- %5.2f", leadText(it).Data(), totalStop[it], totalStopErr[it]);
    const Counts& s = spectrum[it];
    if (s.filled) {
      const double r = s.Rate(kTMCstop) * 60., e = s.RateErr(kTMCstop) * 60.;
      printf(" %13.2f +- %5.2f", r, e);
      if (s.n[kTMCstop] > 0)
        printf(" %14.3f %8.2f", totalStop[it] / r, (totalStop[it] - r) / std::hypot(totalStopErr[it], e));
      if (bestStopCheck < 0 || r > spectrum[bestStopCheck].Rate(kTMCstop) * 60.) bestStopCheck = it;
    }
    printf("\n");
    if (totalStop[it] > totalStop[bestStop]) bestStop = it;
  }
  printf("highest rate of stopped muons: %s cm (sum of the points)", leadText(bestStop).Data());
  if (bestStopCheck >= 0) printf(", %s cm (spectrum check runs)", leadText(bestStopCheck).Data());
  printf("\n(the sum resolves the stopping window only if the energy grid is fine enough; where they\n"
         "disagree, the spectrum check run is the measurement)\n");

  // ---- every point to CSV
  const TString csvName = Form("%s/lead_energy_scan.csv", dir);
  if (FILE* csv = fopen(csvName, "w")) {
    fprintf(csv, "lead_cm,E_GeV,E_low_GeV,E_high_GeV,generated,live_time_s,band_flux_cm2s,plane_cm2");
    for (int s = 0; s < kNSel; ++s) fprintf(csv, ",n_%s", kSelName[s]);
    fprintf(csv, ",ratio_TMC,ratio_TMC_err,Aeff_TMC_cm2,rate_TMC_min,rate_TMC_err_min,rate_TMC_day,"
                 "flux_coeff,ratio_TMC_weighted,N_TMC_per_hour,N_TMC_per_hour_err\n");
    for (int it = 0; it < nT; ++it)
      for (int ie = 0; ie < nE; ++ie) {
        const Counts& p = grid[it][ie];
        if (!p.filled) continue;
        fprintf(csv, "%g,%.6g,%.6g,%.6g,%lld,%.6g,%.6g,%.6g", leads[it], p.eMono / 1000., p.eMin / 1000.,
                p.eMax / 1000., p.generated, p.liveTime, p.flux, p.area);
        for (int s = 0; s < kNSel; ++s) fprintf(csv, ",%lld", p.n[s]);
        fprintf(csv, ",%.6g,%.6g,%.6g,%.6g,%.6g,%.6g", p.Ratio(kTMC), p.RatioErr(kTMC),
                p.Ratio(kTMC) * p.area, p.Rate(kTMC), p.RateErr(kTMC), p.Rate(kTMC) * 60. * 24.);
        fprintf(csv, ",%.6g,%.6g,%.6g,%.6g\n", coeff[ie], weighted(p), expected(p) / exposureHours,
                expectedErr(p) / exposureHours);
      }
    fclose(csv);
    printf("table of every point: %s\n", csvName.Data());
  }

  // ---- plots
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  gStyle->SetNumberContours(64);
  std::vector<double> eEdges;
  if (contiguous) {
    for (int ie = 0; ie < nE; ++ie) eEdges.push_back(show(bandLow[ie]));
    eEdges.push_back(show(bandHigh[nE - 1]));
  } else {
    std::vector<double> centres;
    for (double e : energies) centres.push_back(show(e));
    eEdges = CentreEdges(centres);
  }
  const std::vector<double> tEdges = CentreEdges(leads);
  const TString energyTitle = Form("muon kinetic energy (%s)", unit);
  auto makeMap = [&](const char* name, const char* title, auto value) {
    auto h = new TH2D(name, title, nE, eEdges.data(), nT, tEdges.data());
    for (int it = 0; it < nT; ++it)
      for (int ie = 0; ie < nE; ++ie)
        if (grid[it][ie].filled) h->SetBinContent(ie + 1, it + 1, value(grid[it][ie]));
    h->GetXaxis()->SetTitle(energyTitle);
    h->GetYaxis()->SetTitle("lead thickness (cm)");
    if (logAxis) h->GetXaxis()->SetMoreLogLabels();
    h->GetXaxis()->SetTitleOffset(1.2);
    return h;
  };
  auto drawMap = [&](TH2D* h) {
    gPad->SetLogx(logAxis);
    gPad->SetRightMargin(0.16);
    gPad->SetBottomMargin(0.12);
    h->Draw("COLZ");
  };
  TH2D* hRatio = makeMap("hRatio", "ratio: top & middle & Compton / muons generated (10^{-3})",
                         [](const Counts& p) { return 1e3 * p.Ratio(kTMC); });
  TH2D* hRate = makeMap("hRate", "rate: top & middle & Compton per energy band (min^{-1})",
                        [](const Counts& p) { return p.Rate(kTMC); });
  TH2D* hStop = makeMap("hStop", "same, muon stopped in the Compton region / generated (10^{-3})",
                        [](const Counts& p) { return 1e3 * p.Ratio(kTMCstop); });

  TCanvas c1("c1", "lead energy scan", 1500, 1100);
  c1.Divide(2, 2);
  int pad = 0;
  for (TH2D* h : {hRatio, hRate, hStop}) {
    c1.cd(++pad);
    drawMap(h);
  }
  c1.cd(4);
  gPad->SetGridy();
  gPad->SetBottomMargin(0.12);
  auto mg = new TMultiGraph("mg", Form("top & middle & Compton, muons %s;lead thickness (cm);rate (min^{-1})", rangeText.Data()));
  auto legend = new TLegend(0.4, 0.15, 0.89, 0.29);
  auto gSum = new TGraphErrors();
  auto gSpec = new TGraphErrors();
  for (int it = 0; it < nT; ++it) {
    gSum->SetPoint(gSum->GetN(), leads[it], total[it]);
    gSum->SetPointError(gSum->GetN() - 1, 0., totalErr[it]);
    if (spectrum[it].filled) {
      gSpec->SetPoint(gSpec->GetN(), leads[it], spectrum[it].Rate(kTMC));
      gSpec->SetPointError(gSpec->GetN() - 1, 0., spectrum[it].RateErr(kTMC));
    }
  }
  gSum->SetMarkerStyle(20);
  mg->Add(gSum, "PL");
  legend->AddEntry(gSum, Form("sum of the %d mono-energetic points", nE), "pl");
  if (gSpec->GetN() > 0) {
    gSpec->SetMarkerStyle(24);
    gSpec->SetMarkerColor(kRed + 1);
    gSpec->SetLineColor(kRed + 1);
    mg->Add(gSpec, "P");
    legend->AddEntry(gSpec, "spectrum check run", "pl");
  }
  mg->Draw("A");
  mg->GetYaxis()->SetRangeUser(0., 1.15 * std::max(total[best], gSpec->GetN() > 0 ? TMath::MaxElement(gSpec->GetN(), gSpec->GetY()) : 0.));
  legend->Draw();

  // ratio versus thickness, one curve per energy. The energies below 500 MeV, which range out
  // within the scanned thicknesses, get thick lines in distinct colours and a label at their drop;
  // the others overlap on the plateau and get thin lines in blue-grey shades with open markers.
  TCanvas c2("c2", "ratio versus thickness", 1200, 820);
  gPad->SetRightMargin(0.17);
  gPad->SetLeftMargin(0.09);
  gPad->SetGridy();
  auto mg2 = new TMultiGraph("mg2", "top & middle & Compton / muons generated;lead thickness (cm);ratio (10^{-3})");
  auto legend2 = new TLegend(0.84, 0.1, 0.995, 0.9);
  legend2->SetHeader(Form("E (%s)", unit));
  legend2->SetTextSize(nE > 20 ? 0.024 : 0.026);
  const int boldColours[] = {kRed + 1, kOrange + 7, kGreen + 2, kBlue + 1, kMagenta + 1, kCyan + 2,
                             kOrange - 3, kViolet + 2, kSpring - 6, kAzure + 7, kPink + 9, kGray + 3};
  const int boldMarkers[] = {20, 21, 22, 23, 33, 34, 29, 47, 43, 45, 39, 41};
  const int openMarkers[] = {24, 25, 26, 32, 27, 28, 30, 46, 42, 44};
  int nBold = 0;
  for (double e : energies) if (e < 500.) ++nBold;
  const int nPlain = nE - nBold;
  std::vector<TLatex*> labels;
  std::vector<TGraph*> boldGraphs;
  double maxRatio = 0.;
  for (int ie = 0; ie < nE; ++ie) {
    auto g = new TGraphErrors();
    for (int it = 0; it < nT; ++it) {
      const Counts& p = grid[it][ie];
      if (!p.filled) continue;
      const int n = g->GetN();
      g->SetPoint(n, leads[it], 1e3 * p.Ratio(kTMC));
      g->SetPointError(n, 0., 1e3 * p.RatioErr(kTMC));
      maxRatio = std::max(maxRatio, 1e3 * (p.Ratio(kTMC) + p.RatioErr(kTMC)));
    }
    if (ie < nBold) {
      const int colour = boldColours[ie % 12];
      g->SetLineColor(colour);
      g->SetMarkerColor(colour);
      g->SetLineWidth(nBold > 8 ? 3 : 4);
      g->SetMarkerStyle(boldMarkers[ie % 12]);
      g->SetMarkerSize(nBold > 8 ? 1.1 : 1.4);
      boldGraphs.push_back(g);
      // label just right of where the curve falls below half of its value at the thinnest lead,
      // in staggered rows so that neighbouring labels do not overlap
      for (int k = 1; k < g->GetN(); ++k) {
        if (g->GetPointY(k) < 0.5 * g->GetPointY(0)) {
          const int row = int(labels.size());
          const double y = nBold > 5 ? 0.6 + 0.55 * (row % 6) : 1.0 + 0.45 * row;
          const double x = std::min(g->GetPointX(k) + 0.3, tEdges.back() - (mev ? 1.0 : 2.5));  // inside the frame
          auto label = new TLatex(x, y, mev ? energyText(ie) : energyText(ie) + " GeV");
          label->SetTextColor(colour);
          label->SetTextSize(nBold > 8 ? 0.03 : 0.035);
          label->SetTextFont(62);
          labels.push_back(label);
          break;
        }
      }
    }
    else {
      // light to dark blue-grey with energy
      const double f = nPlain > 1 ? double(ie - nBold) / (nPlain - 1) : 0.;
      const int colour = TColor::GetColor(Float_t(0.70 - 0.55 * f), Float_t(0.75 - 0.55 * f), Float_t(0.85 - 0.35 * f));
      g->SetLineColor(colour);
      g->SetMarkerColor(colour);
      g->SetLineWidth(1);
      g->SetMarkerStyle(openMarkers[(ie - nBold) % 10]);
      g->SetMarkerSize(0.9);
    }
    mg2->Add(g, "PL");
    legend2->AddEntry(g, energyText(ie), "pl");
  }
  mg2->Draw("A");
  mg2->GetXaxis()->SetLimits(tEdges.front(), tEdges.back());
  mg2->GetYaxis()->SetRangeUser(0., 1.08 * maxRatio);
  // draw the low-energy curves again on top of the plateau curves
  for (auto it = boldGraphs.rbegin(); it != boldGraphs.rend(); ++it) (*it)->Draw("PL SAME");
  for (auto label : labels) label->Draw();
  legend2->Draw();

  // flux calibration: coefficients, weighted ratio, expected numbers
  TCanvas c3("c3", "flux calibration", 1500, 1100);
  c3.Divide(2, 2);
  c3.cd(1);
  gPad->SetLogx(logAxis);
  gPad->SetGridx();
  gPad->SetGridy();
  auto frame3 = gPad->DrawFrame(logAxis ? 0.7 * eEdges.front() : eEdges.front(), 0.,
                                logAxis ? 1.25 * eEdges.back() : eEdges.back(), 1.12);
  frame3->SetTitle(Form("flux coefficient w = J_{band} / J_{max};%s;w", energyTitle.Data()));
  auto gCoeff = new TGraph(), gPerGeV = new TGraph();
  for (int ie = 0; ie < nE; ++ie) {
    gCoeff->SetPoint(ie, show(energies[ie]), coeff[ie]);
    if (perGeVMax > 0.) gPerGeV->SetPoint(ie, show(energies[ie]), perGeV[ie] / perGeVMax);
  }
  gCoeff->SetMarkerStyle(20);
  gCoeff->SetMarkerSize(1.2);
  gCoeff->SetLineWidth(3);
  gCoeff->SetLineColor(kBlue + 1);
  gCoeff->SetMarkerColor(kBlue + 1);
  gCoeff->Draw("PL");
  auto legend3 = new TLegend(0.56, 0.76, 0.89, 0.89);
  legend3->AddEntry(gCoeff, equalBands ? "w: band flux / max (= dJ/dE / max)" : "w: band flux / max (used)", "pl");
  if (!equalBands) {
    gPerGeV->SetLineStyle(2);
    gPerGeV->SetLineColor(kGray + 2);
    gPerGeV->SetMarkerStyle(24);
    gPerGeV->SetMarkerColor(kGray + 2);
    gPerGeV->Draw("PL");
    legend3->AddEntry(gPerGeV, "dJ/dE / max (per GeV, for comparison)", "pl");
  }
  legend3->Draw();
  auto hWeighted = makeMap("hWeighted", Form("flux-weighted ratio: ratio #times w (10^{-3}), w = 1 at %s %s",
                                             energyText(iMax).Data(), unit),
                           [&](const Counts& p) { return 1e3 * weighted(p); });
  auto hExpected = makeMap("hExpected", Form("expected top & middle & Compton muons in %g h", exposureHours), expected);
  c3.cd(2);
  drawMap(hWeighted);
  c3.cd(3);
  drawMap(hExpected);
  c3.cd(4);
  gPad->SetGridy();
  gPad->SetBottomMargin(0.12);
  auto mg3 = new TMultiGraph("mg3", Form("expected top & middle & Compton muons in %g h;lead thickness (cm);N", exposureHours));
  auto gAll = new TGraphErrors(), gLow = new TGraphErrors();
  for (int it = 0; it < nT; ++it) {
    double nLow = 0., nLowErr = 0.;
    for (int ie = 0; ie < nE; ++ie)
      if (grid[it][ie].filled && energies[ie] < 500.) {
        nLow += expected(grid[it][ie]);
        nLowErr += std::pow(expectedErr(grid[it][ie]), 2);
      }
    gAll->SetPoint(it, leads[it], total[it] * 60. * exposureHours);
    gAll->SetPointError(it, 0., totalErr[it] * 60. * exposureHours);
    gLow->SetPoint(it, leads[it], nLow);
    gLow->SetPointError(it, 0., std::sqrt(nLowErr));
  }
  gAll->SetMarkerStyle(20);
  mg3->Add(gAll, "PL");
  auto legend4 = new TLegend(0.4, 0.4, 0.89, 0.55);
  legend4->AddEntry(gAll, Form("all %d energies (%s)", nE, rangeText.Data()), "pl");
  if (nPlain > 0 && nBold > 0) {  // only when the low energies are a subset
    gLow->SetMarkerStyle(21);
    gLow->SetMarkerColor(kRed + 1);
    gLow->SetLineColor(kRed + 1);
    mg3->Add(gLow, "PL");
    legend4->AddEntry(gLow, Form("the %d energies below 500 MeV", nBold), "pl");
  }
  mg3->Draw("A");
  mg3->GetYaxis()->SetRangeUser(0., 1.15 * total[best] * 60. * exposureHours);
  legend4->Draw();

  const TString base = Form("%s/lead_energy_scan", dir);
  c1.Print(base + ".pdf(");
  c2.Print(base + ".pdf");
  c3.Print(base + ".pdf)");
  c1.SaveAs(base + "_maps.png");
  c2.SaveAs(base + "_ratio_vs_lead.png");
  c3.SaveAs(base + "_expected.png");
  printf("plots: %s.pdf, %s_maps.png, %s_ratio_vs_lead.png, %s_expected.png\n", base.Data(), base.Data(),
         base.Data(), base.Data());
}
