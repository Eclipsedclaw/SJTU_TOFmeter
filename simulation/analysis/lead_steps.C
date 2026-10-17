// The same muons followed step by step: top & middle & Compton coincidences made by the muon itself
// (it crossed the top and the middle TOF and entered the Compton detector), from the thick-lead scan
// (macro/lead_thick_scan.mac), for several lead thicknesses. Kinetic energy of each muon
//   before the lead      at sea level (ekin; only air between the source and the lead)
//   before the top TOF   entering a top-TOF bar (top_ke), after the lead and the steel plate
//   final state          entering the Compton detector (cmp_ke)
// (the energy entering the middle TOF is not recorded, so it is not shown)
//   root -l -b -q 'analysis/lead_steps.C("output/leadscan_thick")'
// Rates per hour as in analysis/lead_after_spectrum.C: each event weighs 3600 / live time of its
// mono-energetic point, and is spread over the point's band (above 50 MeV, the part below folded
// back; below 50 MeV kept at its own energy). Unlike there, not evenly but with the slope of the
// sea-level spectrum in the band (from the fluxes of the bands on either side), and stretched by
// dx/dE between neighbouring points (from the medians of the energy x, away from the threshold),
// which removes the steps and ripples the grid of energies otherwise leaves in fine bins.
// Also the energy when they leave the lead layer (pb_ke), compared over the thicknesses with the
// spectrum of the same kind of muons without lead (0 cm, at sea level), as the plain lead slab of
// analysis/lead_filter_spectrum.C but for these coincidence muons only.
// The spectra are drawn as flux per m2 s sr GeV: each rate divided by the coincidence's geometric
// factor G(E) = rate / vertical flux at sea level, fitted to the 0 cm points (see below).
// Writes <dir>/lead_steps.pdf with <dir>/lead_steps_flux.png (spectra of the three steps, one panel
// per thickness, over the sea-level spectrum of all muons in grey; _flux_zoom.png: 10-1000 MeV), <dir>/lead_steps_energy.png and _energy_loglog.png (energy after the lead, all
// thicknesses together, linear below 2 GeV and log-log), <dir>/lead_steps_quantiles.png (median and
// spread of the energy at each step), <dir>/lead_steps.csv (quantiles, means), lead_steps_after_lead.csv
// (the linear spectra, rates and fluxes) and lead_steps_geometric_factor.csv (G of each 0 cm point).

namespace {
const int kSteps = 3;
const char* kStepNames[kSteps] = {"before the lead (sea level)", "before the top TOF", "entering the Compton detector (final)"};
const char* kStepShort[kSteps] = {"before lead", "before top TOF", "Compton entry"};

// An event of a mono-energetic point (energy e0, band [lo, hi] at sea level, where the spectrum goes
// as E^-g) found with the energy x (MeV): a muon of the band with energy E would be found at
// x + j (E - e0), with j = jLow below e0 and jHigh above (both 1 at sea level, where x = e0)
struct Spread {
  double x, e0, lo, hi, g, jLow, jHigh;
  double A() const { return x - jLow * (e0 - lo); }   // the energies the event is spread over
  double B() const { return x + jHigh * (hi - e0); }
  double Sea(double y) const { return e0 + (y - x) / (y < x ? jLow : jHigh); }   // energy in the band found at y
  double Integral(double e) const { return std::abs(g - 1.) < 1e-6 ? std::log(e) : std::pow(e, 1. - g) / (1. - g); }
  double Part(double y1, double y2) const   // the fraction of the event found in [y1, y2]
  {
    return (Integral(Sea(std::min(y2, B()))) - Integral(Sea(std::max(y1, A())))) / (Integral(hi) - Integral(lo));
  }
};

// the part [c, d] of the spread, at y or (mirror) at 2 kSpreadMin - y, each bin with its own part
const double kSpreadMin = 50.;   // MeV
void FillPart(TH1D* h, const Spread& sp, double c, double d, double w, bool mirror = false)
{
  const double m = 2. * kSpreadMin, za = mirror ? m - d : c, zb = mirror ? m - c : d;
  auto part = [&](double z1, double z2) { return mirror ? sp.Part(m - z2, m - z1) : sp.Part(z1, z2); };
  const TAxis* axis = h->GetXaxis();
  const double xmin = axis->GetXmin(), xmax = axis->GetXmax();
  if (za < xmin) h->Fill(xmin - 1., w * part(za, std::min(zb, xmin)));
  if (zb > xmax) h->Fill(2. * xmax, w * part(std::max(za, xmax), zb));
  const int first = std::max(1, axis->FindFixBin(za)), last = std::min(axis->GetNbins(), axis->FindFixBin(zb));
  for (int k = first; k <= last; ++k) {
    const double lo = std::max(za, axis->GetBinLowEdge(k)), hi = std::min(zb, axis->GetBinUpEdge(k));
    if (hi > lo) h->Fill(axis->GetBinCenter(k), w * part(lo, hi));
  }
}

// an event after the lead: below 50 MeV kept at x; above, spread, the part that falls below 50 MeV
// folded back above 50 MeV (as analysis/lead_after_spectrum.C, which spreads evenly with j = 1)
void FillBandShift(TH1D* h, const Spread& sp, double w)
{
  const double a = sp.A(), b = sp.B();
  if (sp.x < kSpreadMin || b <= a) {
    h->Fill(sp.x, w);
    return;
  }
  FillPart(h, sp, std::max(a, kSpreadMin), b, w);
  if (a < kSpreadMin) FillPart(h, sp, a, kSpreadMin, w, true);
}

// the selected muons of one mono-energetic point: kinetic energy leaving the lead (pb_ke), before
// the top TOF (top_ke) and entering the Compton detector (cmp_ke), negative if not recorded
struct ScanPoint {
  int it = -1;
  std::size_t band = 0;
  double eMono = 0., eMin = 0., eMax = 0., w = 0.;
  double accept = 0.;     // cm2 sr per selected muon: area x 2 pi PlaneFactor / muons simulated
  double vertical = 0.;   // the generator's vertical flux in the band, per m2 s sr GeV
  std::vector<std::array<double, 3>> muons;
};

// Zenith shape of the generator's muons (as src/CosmicRaySource.cc, Guan et al., arXiv:1509.06176):
// intensity at total energy e (GeV) and cos(theta) c, up to a constant
double GuanIntensity(double e, double c)
{
  const double p1 = 0.102573, p2 = -0.068287, p3 = 0.958633, p4 = 0.0407253, p5 = 0.817285;
  const double cs = std::sqrt((c * c + p1 * p1 + p2 * std::pow(c, p3) + p4 * std::pow(c, p5)) / (1. + p1 * p1 + p2 + p4));
  return std::pow(e * (1. + 3.64 / (e * std::pow(cs, 1.29))), -2.7)
         * (1. / (1. + 1.1 * e * cs / 115.) + 0.054 / (1. + 1.1 * e * cs / 850.));
}

// Int_cosMin^1 I(c) / I(1) c dc at kinetic energy eKin (MeV): the flux through the source plane (per
// cm2 s) is the vertical flux (per cm2 s sr) times 2 pi this
double PlaneFactor(double eKin, double cosMin)
{
  const double e = eKin / 1000. + 0.1056583755;
  const int n = 400;   // Simpson
  double sum = 0.;
  for (int k = 0; k <= n; ++k) {
    const double c = cosMin + (1. - cosMin) * k / n;
    sum += (k == 0 || k == n ? 1. : (k % 2 ? 4. : 2.)) * GuanIntensity(e, c) * c;
  }
  return sum * (1. - cosMin) / (3. * n) / GuanIntensity(e, 1.);
}

double Quantile(const TH1D* h, double q)   // MeV, from the bin contents (log-interpolated in a bin)
{
  const double total = h->Integral();
  double sum = 0.;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    const double c = h->GetBinContent(b);
    if (c > 0. && sum + c >= q * total) {
      const double lo = h->GetBinLowEdge(b), hi = h->GetBinLowEdge(b + 1);
      return lo * std::pow(hi / lo, (q * total - sum) / c);
    }
    sum += c;
  }
  return 0.;
}

TH1D* PerGeV(const TH1D* h, const char* name)
{
  auto d = static_cast<TH1D*>(h->Clone(name));
  for (int b = 1; b <= d->GetNbinsX(); ++b) {
    const double width = h->GetBinWidth(b) / 1000.;
    d->SetBinContent(b, h->GetBinContent(b) / width);
    d->SetBinError(b, h->GetBinError(b) / width);
  }
  return d;
}
}  // namespace

void lead_steps(const char* dir = "output/leadscan_thick", double barThr = 0.2, double comptonThr = 0.1)
{
  TH1::AddDirectory(false);
  const double wanted[] = {0., 5., 10., 20., 30., 40., 50.};
  const int nShow = 7;
  // energy bins: the bands of the scan's points (each "before the lead" point then fills one bin),
  // continued down to 1 MeV with the same width in log E
  std::vector<std::array<double, 3>> bands;   // e_min, e_max, flux
  if (void* d0 = gSystem->OpenDirectory(dir)) {
    while (const char* entry = gSystem->GetDirEntry(d0)) {
      const TString name = entry;
      if (!name.BeginsWith("leadscan_t0_E") || !name.EndsWith(".root")) continue;
      TFile file(Form("%s/%s", dir, entry));
      auto run = dynamic_cast<TTree*>(file.Get("run"));
      if (!run || run->GetEntries() != 1) continue;
      double lo = 0., hi = 0., flux = 0.;
      run->SetBranchAddress("e_min", &lo);
      run->SetBranchAddress("e_max", &hi);
      run->SetBranchAddress("flux", &flux);
      run->GetEntry(0);
      bands.push_back({lo, hi, flux});
    }
    gSystem->FreeDirectory(d0);
  }
  std::sort(bands.begin(), bands.end());
  if (bands.size() < 3) {
    printf("no 0 cm points in %s\n", dir);
    return;
  }
  // the sea-level spectrum in each band as E^-g, g from the flux per MeV of the bands on either
  // side: each event is spread over its band with this shape, so the grid of energies leaves no steps
  std::vector<double> lows, slopes;
  const double top = bands.back()[1];
  for (std::size_t j = 0; j < bands.size(); ++j) {
    const auto& below = bands[j == 0 ? j : j - 1];
    const auto& above = bands[j + 1 == bands.size() ? j : j + 1];
    lows.push_back(bands[j][0]);
    slopes.push_back(-std::log(above[2] / (above[1] - above[0]) / (below[2] / (below[1] - below[0])))
                     / std::log(std::sqrt(above[0] * above[1] / (below[0] * below[1]))));
  }
  const double ratio = lows[1] / lows[0];
  std::vector<double> edges;
  for (double e = lows[0] / ratio; e > 1.; e /= ratio) edges.insert(edges.begin(), e);
  edges.insert(edges.begin(), edges.empty() ? lows[0] / ratio : edges.front() / ratio);
  for (double e : lows) edges.push_back(e);
  edges.push_back(top);
  // rates per hour (quantiles, table, CSV) and, with I in the name, flux per m2 s sr (plots)
  TH1D *steps[nShow][kSteps], *stepsI[nShow][kSteps];
  TH1D *afterLin[nShow], *seaLin[nShow];                       // leaving the lead; sea level, linear
  TH1D *afterLogI[nShow], *afterLinI[nShow], *seaLinI[nShow];
  double rate[nShow] = {}, meanE[nShow][kSteps] = {};
  for (int i = 0; i < nShow; ++i)
    for (int s = 0; s < kSteps; ++s) {
      steps[i][s] = new TH1D(Form("step_%d_%d", i, s), "", int(edges.size()) - 1, edges.data());
      stepsI[i][s] = new TH1D(Form("stepI_%d_%d", i, s), "", int(edges.size()) - 1, edges.data());
      steps[i][s]->Sumw2();
      stepsI[i][s]->Sumw2();
    }
  for (int i = 0; i < nShow; ++i) {
    afterLogI[i] = new TH1D(Form("afterLogI_%d", i), "", int(edges.size()) - 1, edges.data());
    afterLin[i] = new TH1D(Form("afterLin_%d", i), "", 100, 0., 2000.);   // 20 MeV bins
    seaLin[i] = new TH1D(Form("seaLin_%d", i), "", 100, 0., 2000.);
    afterLinI[i] = new TH1D(Form("afterLinI_%d", i), "", 100, 0., 2000.);
    seaLinI[i] = new TH1D(Form("seaLinI_%d", i), "", 100, 0., 2000.);
    TH1D* hs[] = {afterLogI[i], afterLin[i], seaLin[i], afterLinI[i], seaLinI[i]};
    for (TH1D* h : hs) h->Sumw2();
  }
  long nEvents = 0;
  int nFiles = 0;
  std::vector<ScanPoint> points;
  std::vector<double> monos(lows.size(), 0.);   // e_mono of each band

  void* d = gSystem->OpenDirectory(dir);
  if (!d) {
    printf("cannot open %s\n", dir);
    return;
  }
  while (const char* entry = gSystem->GetDirEntry(d)) {
    const TString name = entry;
    if (!name.BeginsWith("leadscan_t") || !name.EndsWith(".root")) continue;
    TFile file(Form("%s/%s", dir, entry));
    auto run = dynamic_cast<TTree*>(file.Get("run"));
    auto events = dynamic_cast<TTree*>(file.Get("events"));
    if (!run || !events || run->GetEntries() != 1 || !events->GetBranch("top_ke") || !events->GetBranch("prim_st")) continue;
    double lead = 0., eMono = 0., eMin = 0., eMax = 0., liveTime = 0., flux = 0., area = 0., thetaMax = 0.;
    int simulated = 0;
    run->SetBranchAddress("lead_cm", &lead);
    run->SetBranchAddress("e_mono", &eMono);
    run->SetBranchAddress("e_min", &eMin);
    run->SetBranchAddress("e_max", &eMax);
    run->SetBranchAddress("live_time", &liveTime);
    run->SetBranchAddress("flux", &flux);
    run->SetBranchAddress("area", &area);
    run->SetBranchAddress("theta_max", &thetaMax);
    run->SetBranchAddress("n_events", &simulated);
    run->GetEntry(0);
    int it = -1;
    for (int i = 0; i < nShow; ++i)
      if (std::abs(lead - wanted[i]) < 1e-6) it = i;
    if (it < 0 || eMono <= 0. || liveTime <= 0.) continue;
    const auto band = std::lower_bound(lows.begin(), lows.end(), eMin * (1. - 1e-9));
    if (band == lows.end() || std::abs(*band - eMin) > 1e-6 * eMin) continue;   // not a band of the 0 cm points
    ++nFiles;
    points.emplace_back();
    ScanPoint& point = points.back();
    point.it = it;
    point.band = band - lows.begin();
    point.eMono = eMono;
    point.eMin = eMin;
    point.eMax = eMax;
    point.w = 3600. / liveTime;
    const double solid = TMath::TwoPi() * PlaneFactor(eMono, std::cos(thetaMax * TMath::DegToRad()));   // sr
    point.accept = simulated > 0 ? area * solid / simulated : 0.;
    point.vertical = flux / solid / ((eMax - eMin) / 1000.) * 1e4;
    monos[point.band] = eMono;
    std::vector<int>* barSt = nullptr;
    std::vector<double>* barE = nullptr;
    double top = 0., cmp = 0., pb = 0., s0 = 0., s1 = 0., c0 = 0., c1 = 0., c2 = 0.;
    int primSt = 0;
    events->SetBranchStatus("*", false);
    const char* used[] = {"bar_st", "bar_e", "top_ke", "cmp_ke", "pb_ke", "e_stack0", "e_stack1", "e_ch0", "e_ch1", "e_ch2", "prim_st"};
    for (const char* b : used) events->SetBranchStatus(b, true);
    events->SetBranchAddress("bar_st", &barSt);
    events->SetBranchAddress("bar_e", &barE);
    events->SetBranchAddress("top_ke", &top);
    events->SetBranchAddress("pb_ke", &pb);
    events->SetBranchAddress("cmp_ke", &cmp);
    events->SetBranchAddress("e_stack0", &s0);
    events->SetBranchAddress("e_stack1", &s1);
    events->SetBranchAddress("e_ch0", &c0);
    events->SetBranchAddress("e_ch1", &c1);
    events->SetBranchAddress("e_ch2", &c2);
    events->SetBranchAddress("prim_st", &primSt);
    for (Long64_t i = 0; i < events->GetEntries(); ++i) {
      events->GetEntry(i);
      bool hitTop = false, hitMiddle = false;
      for (std::size_t k = 0; k < barSt->size(); ++k) {
        if ((*barE)[k] <= barThr) continue;
        hitTop = hitTop || (*barSt)[k] == 0;
        hitMiddle = hitMiddle || (*barSt)[k] == 1;
      }
      const bool compton = s0 > comptonThr || s1 > comptonThr || c0 > comptonThr || c1 > comptonThr || c2 > comptonThr;
      if (!hitTop || !hitMiddle || !compton) continue;
      if ((primSt & 3) != 3 || cmp < 0. || top < 0.) continue;   // the muon itself through all three
      point.muons.push_back({pb, top, cmp});
    }
  }
  gSystem->FreeDirectory(d);

  // An event of a point is spread over the point's band moved to its energy x and stretched by
  // dx/dE, the slope between the medians of x of this point and the next point down (below e_mono)
  // or up (above): the moved bands of neighbouring points then meet. Without the stretch, they
  // overlap by about 3% of their width after thick lead (the faster muons lose a little more), which
  // makes bumps of about 5% at the band edges. Only between points whose muons all came out well
  // above the ~40 MeV needed to reach the Compton detector (5% quantile of x above 100 MeV): nearer
  // that, the medians are pushed up by the muons that did not make it, and stretching by them makes
  // the spectra jump by ~10% from one thickness to the next (with the check runs of all thicknesses).
  const double kStretchMin = 100.;   // MeV
  const std::size_t nBands = lows.size();
  std::vector<std::vector<std::array<double, 3>>> medians(nShow, std::vector<std::array<double, 3>>(nBands, {-1., -1., -1.}));
  auto lowest = medians;   // 5% quantiles
  for (const ScanPoint& p : points)
    for (int q = 0; q < 3; ++q) {
      std::vector<double> x;
      for (const auto& m : p.muons)
        if (m[q] >= 0.) x.push_back(m[q]);
      if (x.size() < 20) continue;
      std::nth_element(x.begin(), x.begin() + x.size() / 2, x.end());
      medians[p.it][p.band][q] = x[x.size() / 2];
      std::nth_element(x.begin(), x.begin() + x.size() / 20, x.end());
      lowest[p.it][p.band][q] = x[x.size() / 20];
    }
  auto stretch = [&](int it, std::size_t j, std::size_t k, int q) -> double {   // between the points j < k
    if (k >= nBands || lowest[it][j][q] < kStretchMin || lowest[it][k][q] < kStretchMin || monos[k] <= monos[j]) return 1.;
    return std::min(std::max((medians[it][k][q] - medians[it][j][q]) / (monos[k] - monos[j]), 0.5), 2.);
  };
  // Flux per sr: the geometric factor of the coincidence, G(E) = rate / vertical flux of the
  // muons at sea level (cm2 sr), from the 0 cm points. It grows with E, ~200 cm2 sr at 0.06 GeV to
  // ~260 at 100 GeV, as the zenith distribution widens; a cubic in ln E fits it from 55 MeV, where
  // the coincidence is fully efficient. Each muon's rate is divided by G at its sea-level energy
  // (which sets its zenith angle), so the 0 cm curves are the generator's vertical flux; below
  // 55 MeV G is taken at 55 MeV, leaving in the curves the muons that do not reach the Compton detector.
  const double kGeoMin = 55.;   // MeV
  TGraphErrors geo;
  for (const ScanPoint& p : points)
    if (p.it == 0 && p.eMono >= kGeoMin && p.muons.size() >= 100) {
      const int k = geo.GetN();
      geo.SetPoint(k, std::log(p.eMono), p.muons.size() * p.accept);
      geo.SetPointError(k, 0., std::sqrt(double(p.muons.size())) * p.accept);
    }
  if (geo.GetN() < 8) {
    printf("too few 0 cm points above %g MeV for the geometric factor in %s\n", kGeoMin, dir);
    return;
  }
  TF1 geoFit("geoFit", "pol3", std::log(kGeoMin), std::log(2e5));
  geo.Fit(&geoFit, "Q0");
  auto geometric = [&](double e) { return geoFit.Eval(std::log(std::max(e, kGeoMin))); };   // cm2 sr

  for (const ScanPoint& p : points) {
    const int it = p.it;
    const double w = p.w, wI = w / 3600. / (1e-4 * geometric(p.eMono)), g = slopes[p.band];
    double jLow[3], jHigh[3];
    for (int q = 0; q < 3; ++q) {
      jLow[q] = p.band > 0 ? stretch(it, p.band - 1, p.band, q) : 1.;
      jHigh[q] = stretch(it, p.band, p.band + 1, q);
    }
    const Spread sea{p.eMono, p.eMono, p.eMin, p.eMax, g, 1., 1.};
    for (const auto& m : p.muons) {
      ++nEvents;
      rate[it] += w;
      FillPart(steps[it][0], sea, p.eMin, p.eMax, w);
      FillPart(stepsI[it][0], sea, p.eMin, p.eMax, wI);
      FillPart(seaLin[it], sea, p.eMin, p.eMax, w);
      FillPart(seaLinI[it], sea, p.eMin, p.eMax, wI);
      if (m[0] >= 0.) {
        const Spread after{m[0], p.eMono, p.eMin, p.eMax, g, jLow[0], jHigh[0]};
        FillBandShift(afterLogI[it], after, wI);
        FillBandShift(afterLin[it], after, w);
        FillBandShift(afterLinI[it], after, wI);
      }
      for (int s = 1; s < kSteps; ++s) {
        const Spread step{m[s], p.eMono, p.eMin, p.eMax, g, jLow[s], jHigh[s]};
        FillBandShift(steps[it][s], step, w);
        FillBandShift(stepsI[it][s], step, wI);
      }
      meanE[it][0] += w * p.eMono;
      meanE[it][1] += w * m[1];
      meanE[it][2] += w * m[2];
    }
  }
  if (nEvents == 0) {
    printf("no muon-only TMC events with top_ke / prim_st in %s\n", dir);
    return;
  }

  // ---- table
  printf("\n=== The same muons step by step: top & middle & Compton made by the muon itself (%s) ===\n", dir);
  printf("%d files, %ld events; TOF bar > %.2f MeV, Compton layer > %.2f MeV. Energies in MeV (median [10%%, 90%%], mean)\n",
         nFiles, nEvents, barThr, comptonThr);
  printf("%5s %9s", "lead", "rate/h");
  for (int s = 0; s < kSteps; ++s) printf(" %28s", kStepShort[s]);
  printf("\n");
  const TString csvName = Form("%s/lead_steps.csv", dir);
  FILE* csv = fopen(csvName, "w");
  if (csv) fprintf(csv, "lead_cm,step,rate_h,mean_MeV,q10_MeV,q25_MeV,median_MeV,q75_MeV,q90_MeV\n");
  for (int i = 0; i < nShow; ++i) {
    if (rate[i] <= 0.) continue;
    printf("%5g %9.0f", wanted[i], rate[i]);
    for (int s = 0; s < kSteps; ++s) {
      const double m = meanE[i][s] / rate[i];
      printf("   %6.0f [%5.0f, %6.0f] %6.0f", Quantile(steps[i][s], 0.5), Quantile(steps[i][s], 0.1), Quantile(steps[i][s], 0.9), m);
      if (csv)
        fprintf(csv, "%g,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", wanted[i], kStepShort[s], rate[i], m, Quantile(steps[i][s], 0.1),
                Quantile(steps[i][s], 0.25), Quantile(steps[i][s], 0.5), Quantile(steps[i][s], 0.75), Quantile(steps[i][s], 0.9));
    }
    printf("\n");
  }
  if (csv) fclose(csv);
  printf("mean energy lost: in the lead and the steel plate (before lead - before top TOF), and in the two TOF bars\n"
         "and two Al plates down to the Compton detector (before top TOF - Compton entry)\n");
  for (int i = 0; i < nShow; ++i) {
    if (rate[i] <= 0.) continue;
    printf("  %2g cm: lead + steel %6.1f MeV, TOF bars + Al plates %5.1f MeV\n", wanted[i],
           (meanE[i][0] - meanE[i][1]) / rate[i], (meanE[i][1] - meanE[i][2]) / rate[i]);
  }
  printf("geometric factor G = rate / vertical flux at sea level (cm2 sr), cubic in ln E (chi2 %.0f, %d points):\n ",
         geoFit.GetChisquare(), geo.GetN());
  for (double e : {60., 300., 1000., 3000., 1e4, 3e4, 1e5}) printf("  %g GeV %.0f", e / 1000., geometric(e));
  printf("\n");
  std::vector<const ScanPoint*> zero;
  for (const ScanPoint& p : points)
    if (p.it == 0) zero.push_back(&p);
  std::sort(zero.begin(), zero.end(), [](const ScanPoint* a, const ScanPoint* b) { return a->eMono < b->eMono; });
  if (FILE* out = fopen(Form("%s/lead_steps_geometric_factor.csv", dir), "w")) {
    fprintf(out, "E_MeV,G_cm2sr,G_err_cm2sr,G_used_cm2sr,generator_vertical_flux_m2ssrGeV\n");
    for (const ScanPoint* p : zero)
      fprintf(out, "%g,%.6g,%.6g,%.6g,%.6g\n", p->eMono, p->muons.size() * p->accept, std::sqrt(double(p->muons.size())) * p->accept,
              geometric(p->eMono), p->vertical);
    fclose(out);
  }

  // ---- page 1: flux spectra of the three steps, one panel per thickness
  gStyle->SetOptStat(0);
  const char* fluxTitle = "flux [m^{2} s sr GeV/n]^{-1}";
  const int stepColours[kSteps] = {kBlack, kBlue + 1, kRed + 1};
  const int stepStyles[kSteps] = {1, 1, 2};   // the Compton entry is only ~6 MeV below the top TOF: dashed, on top
  TCanvas cf("cf", "steps: spectra", 2000, 1050);
  cf.Divide(4, 2);
  double yMax = 0., yMin = 1e30;
  std::vector<std::vector<TH1D*>> curves(nShow);
  for (int i = 0; i < nShow; ++i)
    for (int s = 0; s < kSteps; ++s) {
      TH1D* h = PerGeV(stepsI[i][s], Form("flux_%d_%d", i, s));
      h->SetLineColor(stepColours[s]);
      h->SetLineStyle(stepStyles[s]);
      h->SetLineWidth(2);
      for (int b = 1; b <= h->GetNbinsX(); ++b)
        if (h->GetBinContent(b) > 0.) {
          yMax = std::max(yMax, h->GetBinContent(b));
          yMin = std::min(yMin, h->GetBinContent(b));
        }
      curves[i].push_back(h);
    }
  // the reference in every panel: the sea-level spectrum of all the muons, the generator's vertical
  // flux in each band (the black line is only the muons that make the coincidence at that thickness)
  auto seaAll = static_cast<TH1D*>(stepsI[0][0]->Clone("seaAll"));
  seaAll->Reset();
  for (const ScanPoint& p : points) seaAll->SetBinContent(seaAll->FindFixBin(std::sqrt(p.eMin * p.eMax)), p.vertical);
  seaAll->SetLineColor(kGray);
  seaAll->SetLineWidth(6);   // a light band under the other lines
  for (int i = 0; i < nShow; ++i) {
    cf.cd(i + 1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetBottomMargin(0.12);
    gPad->SetLeftMargin(0.14);
    auto frame = gPad->DrawFrame(1., std::max(0.5 * yMin, 1e-4 * yMax), 1.2e5, 3. * yMax);
    frame->SetTitle(Form("%g cm lead: %.0f muons per hour;kinetic energy (MeV);%s", wanted[i], rate[i], fluxTitle));
    frame->GetXaxis()->SetTitleOffset(1.2);
    seaAll->Draw("HIST SAME ][");
    for (TH1D* h : curves[i]) h->Draw("HIST SAME ][");
  }
  cf.cd(8);
  auto legend = new TLegend(0.05, 0.25, 0.95, 0.85);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.045);
  legend->SetHeader("the same muons (muon only)", "L");
  for (int s = 0; s < kSteps; ++s) legend->AddEntry(curves[0][s], kStepNames[s], "l");
  legend->AddEntry(seaAll, "sea level, all muons (EXPACS)", "l");
  legend->Draw();

  // ---- the same page zoomed to 10-1000 MeV, the y range from the curves in that window
  const double zoomLo = 10., zoomHi = 1000.;
  double zMax = 0., zMin = 1e30;
  for (int i = 0; i < nShow; ++i)
    for (TH1D* h : curves[i])
      for (int b = h->FindFixBin(zoomLo); b <= h->FindFixBin(zoomHi); ++b)
        if (h->GetBinContent(b) > 0.) {
          zMax = std::max(zMax, h->GetBinContent(b));
          zMin = std::min(zMin, h->GetBinContent(b));
        }
  TCanvas cz("cz", "steps: spectra, 10-1000 MeV", 2000, 1050);
  cz.Divide(4, 2);
  for (int i = 0; i < nShow; ++i) {
    cz.cd(i + 1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetBottomMargin(0.12);
    gPad->SetLeftMargin(0.14);
    auto frame = gPad->DrawFrame(zoomLo, std::max(0.5 * zMin, 3e-3 * zMax), zoomHi, 2. * zMax);
    frame->SetTitle(Form("%g cm lead: %.0f muons per hour;kinetic energy (MeV);%s", wanted[i], rate[i], fluxTitle));
    frame->GetXaxis()->SetTitleOffset(1.2);
    frame->GetXaxis()->SetMoreLogLabels();
    frame->GetXaxis()->SetNoExponent();
    seaAll->Draw("HIST SAME ][");
    for (TH1D* h : curves[i]) h->Draw("HIST SAME ][");
  }
  cz.cd(8);
  legend->Draw();

  // ---- the energy after the lead, all thicknesses together (as analysis/lead_filter_spectrum.C)
  const int thickColours[] = {kBlack, kBlue + 1, kAzure + 7, kGreen + 2, kOrange + 7, kRed + 1, kMagenta + 2};
  TCanvas cc("cc", "after the lead, linear", 1100, 800);
  gPad->SetGridy();
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.11);
  TH1D* ref = PerGeV(seaLinI[0], "refLin");
  ref->SetTitle(Form("TMC muons (muon only): kinetic energy leaving the lead, below 2 GeV;kinetic energy (MeV);%s", fluxTitle));
  ref->SetLineColor(kBlack);
  ref->SetLineWidth(3);
  ref->SetMinimum(0.);
  double linTop = ref->GetMaximum();
  std::vector<TH1D*> linCurves;
  for (int i = 1; i < nShow; ++i) {
    if (rate[i] <= 0.) continue;
    TH1D* h = PerGeV(afterLinI[i], Form("afterLinGeV_%d", i));
    h->SetLineColor(thickColours[i]);
    h->SetLineWidth(2);
    h->SetTitle(Form("after %g cm of lead", wanted[i]));
    linTop = std::max(linTop, h->GetMaximum());
    linCurves.push_back(h);
  }
  ref->SetMaximum(1.3 * linTop);
  ref->GetXaxis()->SetTitleOffset(1.2);
  ref->Draw("HIST");
  for (TH1D* h : linCurves) h->Draw("HIST SAME");
  ref->Draw("HIST SAME");
  auto legendC = new TLegend(0.55, 0.62, 0.89, 0.89);
  legendC->SetTextSize(0.03);
  legendC->AddEntry(ref, "no lead (0 cm), at sea level", "l");
  for (TH1D* h : linCurves) legendC->AddEntry(h, h->GetTitle(), "l");
  legendC->Draw();

  TCanvas ccl("ccl", "after the lead, log-log", 1100, 800);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.11);
  TH1D* refLog = PerGeV(stepsI[0][0], "refLog");
  double logMax = refLog->GetMaximum(), logMin = 1e30;
  std::vector<TH1D*> logCurves;
  for (int i = 1; i < nShow; ++i) {
    if (rate[i] <= 0.) continue;
    TH1D* h = PerGeV(afterLogI[i], Form("afterLogGeV_%d", i));
    h->SetLineColor(thickColours[i]);
    h->SetLineWidth(2);
    logMax = std::max(logMax, h->GetMaximum());
    for (int b = 1; b <= h->GetNbinsX(); ++b)
      if (h->GetBinContent(b) > 0.) logMin = std::min(logMin, h->GetBinContent(b));
    logCurves.push_back(h);
  }
  auto frameL = gPad->DrawFrame(10., std::max(0.5 * logMin, 1e-4 * logMax), 1.2e5, 3. * logMax);
  frameL->SetTitle(Form("TMC muons (muon only): kinetic energy leaving the lead;kinetic energy (MeV);%s", fluxTitle));
  frameL->GetXaxis()->SetTitleOffset(1.2);
  refLog->SetLineColor(kBlack);
  refLog->SetLineWidth(3);
  for (TH1D* h : logCurves) h->Draw("HIST SAME ][");
  refLog->Draw("HIST SAME ][");
  auto legendL = new TLegend(0.14, 0.15, 0.45, 0.42);
  legendL->SetTextSize(0.03);
  legendL->AddEntry(refLog, "no lead (0 cm), at sea level", "l");
  for (std::size_t k = 0; k < logCurves.size(); ++k) legendL->AddEntry(logCurves[k], linCurves[k]->GetTitle(), "l");
  legendL->Draw();
  if (FILE* out = fopen(Form("%s/lead_steps_after_lead.csv", dir), "w")) {   // the linear spectra: rates, fluxes
    fprintf(out, "E_low_MeV,E_high_MeV,no_lead_sea_level_hGeV");
    for (int i = 1; i < nShow; ++i)
      if (rate[i] > 0.) fprintf(out, ",after_%gcm_hGeV", wanted[i]);
    fprintf(out, ",no_lead_sea_level_m2ssrGeV");
    for (int i = 1; i < nShow; ++i)
      if (rate[i] > 0.) fprintf(out, ",after_%gcm_m2ssrGeV", wanted[i]);
    fprintf(out, "\n");
    for (int b = 1; b <= ref->GetNbinsX(); ++b) {
      const double width = ref->GetBinWidth(b) / 1000.;
      fprintf(out, "%g,%g,%.6g", ref->GetBinLowEdge(b), ref->GetBinLowEdge(b + 1), seaLin[0]->GetBinContent(b) / width);
      for (int i = 1; i < nShow; ++i)
        if (rate[i] > 0.) fprintf(out, ",%.6g", afterLin[i]->GetBinContent(b) / width);
      fprintf(out, ",%.6g", ref->GetBinContent(b));
      for (TH1D* h : linCurves) fprintf(out, ",%.6g", h->GetBinContent(b));
      fprintf(out, "\n");
    }
    fclose(out);
  }

  // ---- the energy at each step, median with 25-75% (thick) and 10-90% (thin)
  TCanvas ce("ce", "steps: energies", 2000, 1050);
  ce.Divide(4, 2);
  double eMinAll = 1e30, eMaxAll = 0.;
  for (int i = 0; i < nShow; ++i)
    for (int s = 0; s < kSteps; ++s) {
      if (rate[i] <= 0.) continue;
      eMinAll = std::min(eMinAll, Quantile(steps[i][s], 0.1));
      eMaxAll = std::max(eMaxAll, Quantile(steps[i][s], 0.9));
    }
  for (int i = 0; i < nShow; ++i) {
    ce.cd(i + 1);
    gPad->SetLogy();
    gPad->SetGridy();
    gPad->SetBottomMargin(0.16);
    gPad->SetLeftMargin(0.14);
    auto frame = new TH1F(Form("eframe%d", i), Form("%g cm lead;;kinetic energy (MeV)", wanted[i]), kSteps, 0.5, kSteps + 0.5);
    for (int s = 0; s < kSteps; ++s) frame->GetXaxis()->SetBinLabel(s + 1, kStepShort[s]);
    frame->GetXaxis()->SetLabelSize(0.055);
    frame->SetMinimum(0.5 * eMinAll);
    frame->SetMaximum(8. * eMaxAll);   // room for the text above the whiskers
    frame->Draw();   // empty: only the frame and its title on the log scale
    if (rate[i] > 0.) {
      auto losses = new TLatex();
      losses->SetNDC();
      losses->SetTextSize(0.045);
      losses->DrawLatex(0.17, 0.85, Form("mean loss in lead + steel: %.0f MeV", (meanE[i][0] - meanE[i][1]) / rate[i]));
      losses->DrawLatex(0.17, 0.79, Form("in the TOFs and Al plates: %.1f MeV", (meanE[i][1] - meanE[i][2]) / rate[i]));
    }
    auto wide = new TGraphAsymmErrors();
    auto box = new TGraphAsymmErrors();
    auto mean = new TGraph();
    for (int s = 0; s < kSteps; ++s) {
      if (rate[i] <= 0.) continue;
      const double q10 = Quantile(steps[i][s], 0.1), q25 = Quantile(steps[i][s], 0.25), q50 = Quantile(steps[i][s], 0.5),
                   q75 = Quantile(steps[i][s], 0.75), q90 = Quantile(steps[i][s], 0.9);
      wide->SetPoint(s, s + 1, q50);
      wide->SetPointError(s, 0., 0., q50 - q10, q90 - q50);
      box->SetPoint(s, s + 1, q50);
      box->SetPointError(s, 0.12, 0.12, q50 - q25, q75 - q50);
      mean->SetPoint(s, s + 1, meanE[i][s] / rate[i]);
    }
    wide->SetLineWidth(1);
    wide->SetLineColor(kGray + 2);
    wide->Draw("Z SAME");
    box->SetFillColor(kAzure - 9);
    box->SetLineColor(kAzure + 2);
    box->Draw("2 SAME");
    box->SetMarkerStyle(20);
    box->SetMarkerColor(kBlack);
    box->SetLineWidth(2);
    box->Draw("PL SAME");
    mean->SetMarkerStyle(24);
    mean->SetMarkerColor(kRed + 1);
    mean->SetMarkerSize(1.2);
    mean->Draw("P SAME");
  }
  ce.cd(8);
  auto legend2 = new TLegend(0.05, 0.35, 0.95, 0.85);
  legend2->SetBorderSize(0);
  legend2->SetTextSize(0.055);
  legend2->SetHeader("energy of the same muons (muon only)", "L");
  auto sampleBox = new TGraphAsymmErrors();
  sampleBox->SetFillColor(kAzure - 9);
  sampleBox->SetLineColor(kAzure + 2);
  sampleBox->SetMarkerStyle(20);
  auto sampleWide = new TGraph();
  sampleWide->SetLineColor(kGray + 2);
  auto sampleMean = new TGraph();
  sampleMean->SetMarkerStyle(24);
  sampleMean->SetMarkerColor(kRed + 1);
  legend2->AddEntry(sampleBox, "median, middle 50% (box)", "pf");
  legend2->AddEntry(sampleWide, "middle 80% (10-90%)", "l");
  legend2->AddEntry(sampleMean, "mean", "p");
  legend2->Draw();

  const TString base = Form("%s/lead_steps", dir);
  cf.Print(base + ".pdf(");
  cz.Print(base + ".pdf");
  cc.Print(base + ".pdf");
  ccl.Print(base + ".pdf");
  ce.Print(base + ".pdf)");
  cf.SaveAs(base + "_flux.png");
  cz.SaveAs(base + "_flux_zoom.png");
  cc.SaveAs(base + "_energy.png");
  ccl.SaveAs(base + "_energy_loglog.png");
  ce.SaveAs(base + "_quantiles.png");
  printf("plots: %s.pdf, %s_flux.png, %s_flux_zoom.png, %s_energy.png, %s_energy_loglog.png, %s_quantiles.png; tables: %s, %s_after_lead.csv\n",
         base.Data(), base.Data(), base.Data(), base.Data(), base.Data(), base.Data(), csvName.Data(), base.Data());
}
