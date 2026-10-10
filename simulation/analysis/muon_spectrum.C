// Initial-energy spectrum of the generated muons against the atmospheric muon flux.
//   root -l -b -q analysis/muon_spectrum.C
//   root -l -b -q 'analysis/muon_spectrum.C("output/muon_spectrum.root", "output/leadscan")'
// - spectrumFile (macro/muon_spectrum_check.mac, every event kept): kinetic energies of the
//   generated mu+ and mu-, as the flux dJ/dE through the horizontal source plane (zenith < 75 deg)
// - curves: the EXPACS table (data/expacs_spectrum.txt, mu+ + mu-, omnidirectional flux times the
//   fraction crossing the plane for the Guan zenith shape, as the generator does) and, independent
//   of EXPACS, the Guan et al. sea-level formula (arXiv:1509.06176) for the same plane and angles
// - scanDir (macro/lead_energy_scan.mac): the band flux every mono-energetic point stands for, and
//   a check that each point's muons all have its energy; plus the energy spectrum of the
//   top & middle & Compton coincidences for 0, 10 and 20 cm of lead
// Prints integral numbers next to standard sea-level values (PDG Review of Particle Physics,
// "Cosmic rays": vertical intensity above 1 GeV/c ~ 70 /m2/s/sr, ~1 /cm2/min through a horizontal
// surface, mean energy ~4 GeV, mu+/mu- ~ 1.27). Writes <spectrum file>.pdf and .png.

namespace {
const double kMuMass = 105.658;  // MeV

// Guan et al. dI/dE in 1/(cm2 s sr GeV), total energy in GeV (as in CosmicRaySource.cc)
double CosThetaStar(double c)
{
  const double p1 = 0.102573, p2 = -0.068287, p3 = 0.958633, p4 = 0.0407253, p5 = 0.817285;
  return std::sqrt((c * c + p1 * p1 + p2 * std::pow(c, p3) + p4 * std::pow(c, p5)) / (1. + p1 * p1 + p2 + p4));
}
double GuanIntensity(double e, double c)
{
  const double cs = CosThetaStar(c);
  return 0.14 * std::pow(e * (1. + 3.64 / (e * std::pow(cs, 1.29))), -2.7)
         * (1. / (1. + 1.1 * e * cs / 115.) + 0.054 / (1. + 1.1 * e * cs / 850.));
}
// Integral over cos(theta) in [a, 1] of I(e, c) c^k
double GuanCosIntegral(double e, double a, int k)
{
  const int n = 400;
  double sum = 0.;
  for (int i = 0; i <= n; ++i) {
    const double c = a + (1. - a) * i / n;
    const double w = (i == 0 || i == n) ? 1. : ((i % 2) ? 4. : 2.);
    sum += w * GuanIntensity(e, c) * std::pow(c, k);
  }
  return sum * (1. - a) / n / 3.;
}

struct RunSummary {
  double liveTime = 0., area = 0., flux = 0., eMin = 0., eMax = 0., eMono = 0., lead = 0.;
  int n = 0;
};
bool ReadRun(TFile& file, RunSummary& r)
{
  auto run = dynamic_cast<TTree*>(file.Get("run"));
  if (!run || run->GetEntries() != 1) return false;
  run->SetBranchAddress("n_events", &r.n);
  run->SetBranchAddress("live_time", &r.liveTime);
  run->SetBranchAddress("area", &r.area);
  run->SetBranchAddress("flux", &r.flux);
  run->SetBranchAddress("e_min", &r.eMin);
  run->SetBranchAddress("e_max", &r.eMax);
  if (run->GetBranch("e_mono")) run->SetBranchAddress("e_mono", &r.eMono);
  if (run->GetBranch("lead_cm")) run->SetBranchAddress("lead_cm", &r.lead);
  run->GetEntry(0);
  return true;
}
}  // namespace

void muon_spectrum(const char* spectrumFile = "output/muon_spectrum.root", const char* scanDir = "output/leadscan",
                   const char* expacsFile = "data/expacs_spectrum.txt")
{
  // ---- generated sample
  TFile file(spectrumFile);
  auto events = dynamic_cast<TTree*>(file.Get("events"));
  RunSummary run;
  if (!events || !ReadRun(file, run) || run.eMono != 0.) {
    printf("%s: no spectrum run (run macro/muon_spectrum_check.mac first)\n", spectrumFile);
    return;
  }
  if (events->GetEntries() != run.n) printf("warning: %lld events stored of %d generated (event filter?)\n",
                                            events->GetEntries(), run.n);
  const double cosMin = std::cos(75. * TMath::DegToRad());
  const int nBins = 60;
  double edges[nBins + 1];
  for (int i = 0; i <= nBins; ++i) edges[i] = std::pow(10., -3. + 6. * i / nBins);  // GeV
  auto hGen = new TH1D("hGen", "", nBins, edges);
  events->Draw("ekin/1000.>>hGen", "", "goff");
  const double norm = 1. / (run.liveTime * run.area);  // per cm2 per s
  for (int i = 1; i <= nBins; ++i) {
    const double w = hGen->GetBinWidth(i);
    hGen->SetBinContent(i, hGen->GetBinContent(i) * norm / w);
    hGen->SetBinError(i, hGen->GetBinError(i) * norm / w);
  }

  // ---- curves: EXPACS table and the Guan formula, flux through the plane in 1/(cm2 s GeV)
  auto gExpacs = new TGraph();
  {
    std::ifstream in(expacsFile);
    std::string line;
    int iPlus = -1, iMinus = -1;
    while (std::getline(in, line)) {
      if (line.rfind("# columns:", 0) == 0) {
        std::istringstream names(line.substr(10));
        std::string name;
        for (int k = 0; names >> name; ++k) {
          if (name == "mu+") iPlus = k;
          if (name == "mu-") iMinus = k;
        }
        continue;
      }
      if (line.empty() || line[0] == '#' || iPlus < 0) continue;
      std::istringstream values(line);
      std::vector<double> v;
      for (double x; values >> x;) v.push_back(x);
      if (int(v.size()) <= std::max(iPlus, iMinus)) continue;
      const double eGeV = v[0] / 1000., phi = (v[iPlus] + v[iMinus]) * 1000.;  // 1/(cm2 s GeV), omnidirectional
      if (phi <= 0. || eGeV < 1e-3 || eGeV > 1e3) continue;
      const double etot = eGeV + kMuMass / 1000.;
      const double kappa = GuanCosIntegral(etot, cosMin, 1) / GuanCosIntegral(etot, 0., 0);
      gExpacs->SetPoint(gExpacs->GetN(), eGeV, phi * kappa);
    }
  }
  auto gGuan = new TGraph();
  for (int i = 0; i <= 120; ++i) {
    const double eGeV = std::pow(10., -3. + 6. * i / 120.);
    gGuan->SetPoint(i, eGeV, TMath::TwoPi() * GuanCosIntegral(eGeV + kMuMass / 1000., cosMin, 1));
  }
  if (gExpacs->GetN() < 2) {
    printf("cannot read the mu+/mu- columns of %s\n", expacsFile);
    return;
  }
  auto gLogExpacs = new TGraph();  // log J versus log E, for interpolation
  for (int i = 0; i < gExpacs->GetN(); ++i)
    gLogExpacs->SetPoint(i, std::log(gExpacs->GetPointX(i)), std::log(gExpacs->GetPointY(i)));
  auto expacsAt = [&](double e) { return std::exp(gLogExpacs->Eval(std::log(e))); };

  // ---- integral numbers: simulation, the Guan formula (same plane and angles), reference
  const TString p = Form("sqrt(ekin*ekin + 2*ekin*%g)", kMuMass);
  const double vertCone = TMath::TwoPi() * 0.05 * 0.975;  // cos(theta) > 0.95: solid angle x mean cos
  const double nVert = events->GetEntries("costh > 0.95 && " + p + " > 1000");
  const double iVert = nVert / (run.liveTime * run.area * vertCone) * 1e4;  // /m2/s/sr
  const double nPlus = events->GetEntries("pdg == -13 && " + p + " > 1000");
  const double nMinus = events->GetEntries("pdg == 13 && " + p + " > 1000");
  auto meanOf = [&](const char* cut) {
    const Long64_t n = events->Draw("ekin/1000.", cut, "goff");
    double sum = 0.;
    for (Long64_t i = 0; i < n; ++i) sum += events->GetV1()[i];
    return n > 0 ? sum / n : 0.;
  };
  const double meanAll = meanOf(""), meanVert = meanOf("costh > 0.95");
  // Guan: flux, mean energies and vertical intensity by integration over E (log steps)
  double gFlux = 0., gSumE = 0., gVertN = 0., gVertE = 0., gVertI = 0.;
  for (int i = 0; i < 3000; ++i) {
    const double e = std::pow(10., -3. + 6. * (i + 0.5) / 3000.), de = e * std::log(10.) * 6. / 3000.;
    const double etot = e + kMuMass / 1000.;
    const double j = TMath::TwoPi() * GuanCosIntegral(etot, cosMin, 1) * de;
    const double jv = TMath::TwoPi() * GuanCosIntegral(etot, 0.95, 1) * de;
    gFlux += j;
    gSumE += e * j;
    gVertN += jv;
    gVertE += e * jv;
    if (std::sqrt(e * e + 2. * e * kMuMass / 1000.) > 1.) gVertI += GuanIntensity(etot, 1.) * de;
  }
  printf("\n=== Generated muon spectrum: %d EXPACS muons, %.0f s live time, plane %.0f cm2, zenith < 75 deg ===\n",
         run.n, run.liveTime, run.area);
  printf("%-48s %12s %10s   %s\n", "", "simulation", "Guan", "sea-level reference (PDG)");
  printf("%-48s %12.3f %10.3f   ~1 for a horizontal detector\n", "flux through the horizontal plane (/cm2/min)",
         run.flux * 60., gFlux * 60.);
  printf("%-48s %7.1f+-%3.1f %10.1f   ~70\n", "vertical intensity, p > 1 GeV/c (/m2/s/sr)", iVert,
         iVert / std::sqrt(std::max(nVert, 1.)), gVertI * 1e4);
  printf("%-48s %12.2f %10.2f   ~4 (mean energy at ground)\n", "mean kinetic energy, near vertical (GeV)", meanVert,
         gVertE / gVertN);
  printf("%-48s %12.2f %10.2f\n", "mean kinetic energy, zenith < 75 deg (GeV)", meanAll, gSumE / gFlux);
  printf("%-48s %12.3f %10s   ~1.27 (1-100 GeV/c)\n", "mu+ / mu-, p > 1 GeV/c", nMinus > 0 ? nPlus / nMinus : 0.,
         "1.27 (set)");
  printf("near vertical = cos(theta) > 0.95. EXPACS here is for 35N 142E (cut-off 11.5 GV), not Shanghai.\n");

  // ---- scan points: band flux each stands for, and their energies
  auto gBand = new TGraphErrors();
  std::map<double, TGraph*> tmcSpectra;  // lead cm -> rate per GeV of top & middle & Compton
  int nFiles = 0, nBad = 0;
  double maxDev = 0.;
  if (void* dirp = gSystem->OpenDirectory(scanDir)) {
    const TString tmc = "Sum$(bar_st==0 && bar_e>0.2)>0 && Sum$(bar_st==1 && bar_e>0.2)>0 && "
                        "(e_stack0>0.1 || e_stack1>0.1 || e_ch0>0.1 || e_ch1>0.1 || e_ch2>0.1)";
    std::map<double, std::vector<std::pair<double, double>>> tmcPoints;
    while (const char* entry = gSystem->GetDirEntry(dirp)) {
      const TString name = entry;
      if (!name.BeginsWith("leadscan_t") || !name.EndsWith(".root")) continue;
      TFile f(Form("%s/%s", scanDir, entry));
      RunSummary r;
      auto t = dynamic_cast<TTree*>(f.Get("events"));
      if (!t || !ReadRun(f, r) || r.eMono <= 0.) continue;
      ++nFiles;
      const double lo = t->GetMinimum("ekin"), hi = t->GetMaximum("ekin");
      const double dev = std::max(std::abs(lo / r.eMono - 1.), std::abs(hi / r.eMono - 1.));
      maxDev = std::max(maxDev, dev);
      if (dev > 1e-6) {
        ++nBad;
        printf("  %s: ekin %g - %g MeV, expected %g\n", entry, lo, hi, r.eMono);
      }
      const double width = (r.eMax - r.eMin) / 1000.;  // GeV
      if (std::abs(r.lead) < 1e-6) gBand->SetPoint(gBand->GetN(), r.eMono / 1000., r.flux / width);
      if (std::abs(r.lead) < 1e-6 || std::abs(r.lead - 10.) < 1e-6 || std::abs(r.lead - 20.) < 1e-6)
        tmcPoints[r.lead].push_back({r.eMono / 1000., t->GetEntries(tmc) / r.liveTime * 60. / width});
    }
    gSystem->FreeDirectory(dirp);
    for (auto& [lead, pts] : tmcPoints) {
      std::sort(pts.begin(), pts.end());
      auto g = new TGraph();
      for (auto& [e, rate] : pts) if (rate > 0.) g->SetPoint(g->GetN(), e, rate);
      tmcSpectra[lead] = g;
    }
  }
  if (nFiles > 0)
    printf("\nscan points in %s: %d files, %d with a muon energy other than their nominal one "
           "(largest relative deviation %.1e)\n", scanDir, nFiles, nBad, maxDev);
  if (gBand->GetN() > 0) {
    gBand->Sort();
    printf("band flux of the 20 points / EXPACS curve at the same energy:");
    for (int i = 0; i < gBand->GetN(); ++i) printf(" %.3f", gBand->GetPointY(i) / expacsAt(gBand->GetPointX(i)));
    printf("\n  (the band flux is an average over E/1.2 - Ex1.2, so it sits slightly off a curved spectrum)\n");
  }

  // ---- plots
  gStyle->SetOptStat(0);
  TCanvas c1("c1", "muon spectrum", 1100, 1000);
  auto top = new TPad("top", "", 0., 0.32, 1., 1.);
  auto bottom = new TPad("bottom", "", 0., 0., 1., 0.32);
  top->SetBottomMargin(0.02);
  bottom->SetTopMargin(0.03);
  bottom->SetBottomMargin(0.3);
  for (auto pad : {top, bottom}) {
    pad->SetLogx();
    pad->SetGridx();
    pad->SetGridy();
    pad->SetLeftMargin(0.12);
    pad->Draw();
  }
  top->cd();
  top->SetLogy();
  auto frame = top->DrawFrame(1e-3, 1e-9, 1e3, 2e-2);
  frame->SetTitle("sea-level muons through the source plane (zenith < 75#circ);;dJ/dE (cm^{-2} s^{-1} GeV^{-1})");
  frame->GetXaxis()->SetLabelSize(0);
  gExpacs->SetLineColor(kBlue + 1);
  gExpacs->SetLineWidth(3);
  gExpacs->Draw("L");
  gGuan->SetLineColor(kRed + 1);
  gGuan->SetLineWidth(2);
  gGuan->SetLineStyle(2);
  gGuan->Draw("L");
  hGen->SetMarkerStyle(20);
  hGen->SetMarkerSize(0.8);
  hGen->SetMarkerColor(kBlack);
  hGen->SetLineColor(kBlack);
  hGen->Draw("E1 SAME");
  if (gBand->GetN() > 0) {
    gBand->SetMarkerStyle(21);
    gBand->SetMarkerSize(1.3);
    gBand->SetMarkerColor(kGreen + 2);
    gBand->SetLineColor(kGreen + 2);
    gBand->Draw("P");
  }
  auto legend = new TLegend(0.15, 0.06, 0.62, 0.32);
  legend->AddEntry(hGen, Form("generated, %d muons", run.n), "pe");
  legend->AddEntry(gExpacs, "EXPACS table (mu^{+} + mu^{-})", "l");
  legend->AddEntry(gGuan, "Guan et al. formula", "l");
  if (gBand->GetN() > 0) legend->AddEntry(gBand, "band flux of the 20 scan points", "p");
  legend->Draw();

  bottom->cd();
  auto hRatio = static_cast<TH1D*>(hGen->Clone("hRatio"));
  auto gGuanRatio = new TGraph();
  for (int i = 1; i <= nBins; ++i) {
    const double e = hRatio->GetBinCenter(i);
    // average of the curve over the bin (integral in log steps / width), as the histogram is
    const int sub = 20;
    const double dlog = std::log(edges[i] / edges[i - 1]) / sub;
    double integral = 0.;
    for (int k = 0; k < sub; ++k) {
      const double ek = edges[i - 1] * std::exp((k + 0.5) * dlog);
      integral += expacsAt(ek) * ek * dlog;
    }
    const double avg = integral / hRatio->GetBinWidth(i);
    if (avg > 0. && hGen->GetBinContent(i) > 0.) {
      hRatio->SetBinContent(i, hGen->GetBinContent(i) / avg);
      hRatio->SetBinError(i, hGen->GetBinError(i) / avg);
    } else {
      hRatio->SetBinContent(i, 0.);
      hRatio->SetBinError(i, 0.);
    }
    gGuanRatio->SetPoint(gGuanRatio->GetN(), e, gGuan->Eval(e) / expacsAt(e));
  }
  hRatio->GetYaxis()->SetRangeUser(0., 2.);
  hRatio->GetYaxis()->SetTitle("/ EXPACS");
  hRatio->GetXaxis()->SetTitle("muon kinetic energy (GeV)");
  for (auto axis : {hRatio->GetXaxis(), hRatio->GetYaxis()}) {
    axis->SetLabelSize(0.09);
    axis->SetTitleSize(0.1);
  }
  hRatio->GetYaxis()->SetTitleOffset(0.5);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->Draw("E1");
  gGuanRatio->SetLineColor(kRed + 1);
  gGuanRatio->SetLineWidth(2);
  gGuanRatio->SetLineStyle(2);
  gGuanRatio->Draw("L");
  TLine one(1e-3, 1., 1e3, 1.);
  one.SetLineColor(kBlue + 1);
  one.Draw();

  // energy spectrum of the coincidences
  TCanvas c2("c2", "coincidence spectrum", 1100, 750);
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();
  auto frame2 = gPad->DrawFrame(0.07, 0., 150., 1.);
  frame2->SetTitle("muons giving top & middle & Compton;muon kinetic energy (GeV);dR/dE (min^{-1} GeV^{-1})");
  auto legend2 = new TLegend(0.6, 0.7, 0.89, 0.89);
  legend2->SetHeader("lead");
  const int colours[] = {kBlack, kBlue + 1, kRed + 1};
  int k = 0;
  double ymax = 0.;
  for (auto& [lead, g] : tmcSpectra) {
    for (int i = 0; i < g->GetN(); ++i) ymax = std::max(ymax, g->GetPointY(i));
    g->SetLineColor(colours[k % 3]);
    g->SetMarkerColor(colours[k % 3]);
    g->SetMarkerStyle(20 + k);
    g->SetLineWidth(2);
    g->Draw("PL");
    legend2->AddEntry(g, Form("%.0f cm", lead), "pl");
    ++k;
  }
  frame2->GetYaxis()->SetRangeUser(0., 1.1 * ymax);
  legend2->Draw();

  TString base = spectrumFile;
  base.ReplaceAll(".root", "");
  c1.Print(base + ".pdf(");
  c2.Print(base + ".pdf)");
  c1.SaveAs(base + ".png");
  c2.SaveAs(base + "_coincidences.png");
  printf("plots: %s.pdf, %s.png, %s_coincidences.png\n", base.Data(), base.Data(), base.Data());
}
