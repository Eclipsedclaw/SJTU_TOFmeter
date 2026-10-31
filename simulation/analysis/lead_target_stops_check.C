// Why the mu- stops in the target hardly change with the top lead (analysis/lead_target_stops.C: within 4% from 0 to
// 30 cm), checked on the lead-slab runs of scripts/lead_filter.sh (EXPACS sea-level muons through a plain lead slab in
// vacuum, Geant4). A muon stops where its range ends, so a thin target takes one narrow slice of the sea-level spectrum,
// and the lead above only chooses which slice:
//     stops per day = j(E0) x dE0,   dE0 = (the range the target spans) x S(E0),
// j the flux per MeV at the energy E0 whose range ends in the target, S the stopping power there. A slice at higher E0
// has more muons per MeV but is narrower in MeV, as S falls towards its minimum; from ~0.1 to ~1 GeV the product j x S
// (the range spectrum: muons per g/cm2 of range) changes little. All for mu- (the ones with muonic X-rays):
//  1. the sea-level spectrum per MeV and per unit of ln E (equal steps on a log axis, wider in MeV at higher energy);
//  2. the slice: energy before the lead of the near-vertical (cos > 0.9) mu- that leave it with the energies a vertical
//     mu- needs to stop in the 1 cm target (CSDA through the steel plate, TOF stations, Al plates and air);
//  3. j, S and j x S of the near-vertical mu- against E0, with the Geant4 number in each slice at its E0;
//  4. Geant4 stopping depths, all angles: the mu- per day that end their range in each cm of lead, against straight CSDA
//     ranges of the same muons;
//  5. directions of the slow mu- leaving the lead (all angles): scattered in it, more come at angles the top opening of
//     the walls (seen from the target) does not cover; analysis/lead_target_stops.C leaves those tracks out;
//  6. all relative to 10 cm, with the stops of analysis/lead_target_stops.C (read from <out>/lead_target_stops.csv);
//     printed also the slow mu- through the top opening onto the target under 10 cm, as a check of its absolute scale.
// Geant4 has fewer slow mu- in a slice than the straight CSDA ranges, the more so the thicker the lead: scattered, a
// slow muon zigzags, so near the end of its range its energy falls faster with depth than along a straight line.
//   root -l -b -q 'analysis/lead_target_stops_check.C("output/lead_filter_below2GeV")'
// Writes <out>/lead_target_stops_check.png/.pdf and prints the numbers.

#include <fstream>
#include <functional>
#include <sstream>

namespace {
constexpr double kMuon = 105.6584, kElectron = 0.51099895;   // MeV

double Interpolate(const std::vector<double>& x, const std::vector<double>& y, double v)
{
  std::size_t j = std::lower_bound(x.begin(), x.end(), v) - x.begin();
  j = std::clamp<std::size_t>(j, 1, x.size() - 1);
  return y[j - 1] + (y[j] - y[j - 1]) * (v - x[j - 1]) / (x[j] - x[j - 1]);
}

// CSDA ranges as in analysis/lead_target_stops.C: Bethe-Bloch with the density effect (Sternheimer parameters, PDG)
struct Medium {
  double za, iEv, rho, cBar, x0, x1, a, k, delta0;
  std::vector<double> logT, logR;
  double DEdx(double T) const   // MeV cm2/g
  {
    const double g = 1. + T / kMuon, bg2 = g * g - 1., b2 = bg2 / (g * g);
    const double wMax = 2. * kElectron * bg2 / (1. + 2. * g * kElectron / kMuon + std::pow(kElectron / kMuon, 2));
    const double x = 0.5 * std::log10(bg2), ln10 = std::log(10.);
    const double delta = x >= x1 ? 2. * ln10 * x - cBar
                         : x >= x0 ? 2. * ln10 * x - cBar + a * std::pow(x1 - x, k)
                                   : delta0 * std::pow(10., 2. * (x - x0));
    const double i = iEv * 1e-6;
    return 0.307075 * za / b2 * (0.5 * std::log(2. * kElectron * bg2 * wMax / (i * i)) - b2 - 0.5 * delta);
  }
  void Build()   // range table 1 MeV - 20 GeV; below 1 MeV (< 0.03 g/cm2) the range is taken as T / (2 dE/dx)
  {
    const int n = 4000;
    const double step = std::pow(2e4, 1. / n);
    double t = 1., r = 0.5 / DEdx(1.);
    for (int j = 0; j <= n; ++j, t *= step) {
      logT.push_back(std::log(t));
      logR.push_back(std::log(r));
      r += 0.5 * (1. / DEdx(t) + 1. / DEdx(t * step)) * t * (step - 1.);
    }
  }
  double Range(double T) const { return std::exp(Interpolate(logT, logR, std::log(T))); }    // g/cm2
  double Energy(double r) const { return std::exp(Interpolate(logR, logT, std::log(r))); }   // MeV, of range r g/cm2
  double After(double T, double x) const   // kinetic energy after x g/cm2, -1 if it stops in them
  {
    const double r = Range(T) - x;
    return r > 0. ? Energy(r) : -1.;
  }
  double Before(double T, double x) const { return Energy(Range(T) + x); }   // kinetic energy before x g/cm2, to leave with T
};
}  // namespace

void lead_target_stops_check(const char* dir = "output/lead_filter_below2GeV", const char* out = "output/lead_target_stops",
                             unsigned seed = 4357)
{
  TH1::AddDirectory(false);
  Medium steel{0.4658, 286.0, 8.00, 4.2911, -0.0012, 3.1531, 0.14680, 2.9632, 0.12};   // G4_STAINLESS-STEEL (Fe parameters)
  Medium plastic{0.54141, 64.7, 1.032, 3.1997, 0.1464, 2.4855, 0.16101, 3.2393, 0.};   // G4_PLASTIC_SC_VINYLTOLUENE
  Medium aluminium{0.48181, 166.0, 2.699, 4.2395, 0.1708, 3.0127, 0.08024, 3.6345, 0.12};
  Medium air{0.49919, 85.7, 1.205e-3, 10.5961, 1.7418, 4.2759, 0.10914, 3.3994, 0.};
  Medium lead{0.39575, 823.0, 11.35, 6.2018, 0.3776, 3.8073, 0.09359, 3.1608, 0.14};
  for (Medium* m : {&steel, &plastic, &aluminium, &air, &lead}) m->Build();

  // the energies leaving the lead (MeV) with which a vertical mu- stops in the 1 cm target: it gets through the steel
  // plate, top TOF, Al1, middle TOF, Al2 and the air (cm, as analysis/lead_target_stops.C), but not through the target
  const std::vector<std::pair<const Medium*, double>> above = {{&steel, 1.},   {&plastic, 0.6}, {&aluminium, 0.3},
                                                               {&plastic, 0.6}, {&aluminium, 0.3}, {&air, 51.7}};
  double aboveGcm2 = 0.;
  for (const auto& [medium, cm] : above) aboveGcm2 += cm * medium->rho;
  auto atTarget = [&](double T) {   // kinetic energy at the target top, -1 if it stops before
    for (const auto& [medium, cm] : above)
      if ((T = medium->After(T, cm * medium->rho)) <= 0.) return -1.;
    return T;
  };
  auto lowest = [](const std::function<bool(double)>& passes) {   // the lowest energy that passes (bisection)
    double a = 0.1, b = 2000.;
    for (int i = 0; i < 60; ++i) {
      const double m = 0.5 * (a + b);
      if (passes(m)) b = m; else a = m;
    }
    return b;
  };
  const double tLo = lowest([&](double T) { return atTarget(T) > 0.; });
  const double tHi = lowest([&](double T) {
    const double e = atTarget(T);
    return e > 0. && aluminium.After(e, aluminium.rho) > 0.;
  });
  const double targetTop = 36.5, targetHalf = 10., leadBottom = 91., wallInner = 45., deg = 180. / M_PI;
  // the top opening seen from the target centre: to the middle of an edge and to a corner
  const double edgeAngle = std::atan(wallInner / (leadBottom - targetTop)) * deg;
  const double cornerAngle = std::atan(std::sqrt(2.) * wallInner / (leadBottom - targetTop)) * deg;

  // the files of each thickness, with their exposure and the top of their simulated spectrum
  const std::vector<int> ts = {0, 5, 10, 20, 30, 40, 50};   // 0: no lead (the muons before the slabs)
  const int nT = ts.size();
  auto indexOf = [&](int t) { return int(std::find(ts.begin(), ts.end(), t) - ts.begin()); };
  const int iRef = indexOf(10), iThick = nT - 1;
  struct File { int k; TString path; double live, area, eMax; };
  std::vector<File> files;
  double liveAll = 0.;
  std::vector<double> liveT(nT, 0.);
  if (void* d = gSystem->OpenDirectory(dir)) {
    while (const char* entry = gSystem->GetDirEntry(d)) {
      TString name = entry;
      if (!name.BeginsWith("lead_filter_t") || !name.EndsWith(".root")) continue;
      TString rest = name(13, name.Length() - 18);
      if (rest.Index("_p") != kNPOS) rest = rest(0, rest.Index("_p"));
      const int k = indexOf(rest.Atoi());
      if (k == 0 || k >= nT) continue;
      File f{k, Form("%s/%s", dir, entry), 0., 0., 1e9};
      TFile file(f.path);
      auto run = dynamic_cast<TTree*>(file.Get("run"));
      if (!run || run->GetEntries() != 1) continue;
      run->GetEntry(0);
      auto value = [&](const char* leaf, double fallback) {
        auto l = run->GetLeaf(leaf);
        return l ? l->GetValue() : fallback;
      };
      f.live = value("live_time", 0.);
      f.area = value("area", 0.);
      f.eMax = value("e_max", 1e9);
      if (f.eMax < 100.) f.eMax *= 1000.;   // given in GeV
      files.push_back(f);
      liveAll += f.live;
      liveT[k] += f.live;
    }
    gSystem->FreeDirectory(d);
  }
  if (files.empty()) { printf("no lead_filter_t<cm>_p<k>.root in %s\n", dir); return; }

  // the sea-level mu- spectrum, log bins up to the top of the thickest slab's runs; every file has the same spectrum
  // below its e_max, so it counts only in the bins below that, and adds its live time to them
  const int nE = 120;
  const double eLow = 9.25, eHigh = 5700.;
  std::vector<double> edges(nE + 1);
  for (int b = 0; b <= nE; ++b) edges[b] = eLow * std::pow(eHigh / eLow, double(b) / nE);
  TH1D countAll("countAll", "", nE, edges.data()), countVert("countVert", "", nE, edges.data());
  std::vector<double> exposure(nE + 2, 0.);
  for (const auto& f : files)
    for (int b = 1; b <= nE; ++b)
      if (edges[b] <= f.eMax * (1. + 1e-9)) exposure[b] += f.live;
  std::vector<TH1D*> slice(nT), angle(nT), depth(nT), depthCsda(nT);
  for (int k = 0; k < nT; ++k) {
    slice[k] = new TH1D(Form("slice%d", k), "", 275, 0., 1100.);   // energy before the lead, 4 MeV
    angle[k] = new TH1D(Form("angle%d", k), "", 45, 0., 90.);     // zenith angle leaving the lead, 2 deg
    depth[k] = new TH1D(Form("depth%d", k), "", 50, 0., 50.);     // stopping depth in the lead, 1 cm
    depthCsda[k] = new TH1D(Form("depthCsda%d", k), "", 50, 0., 50.);
  }
  // per m2 per day: the slow mu- (leaving the lead with tLo - tHi) near vertical, their sum of w^2, the same from CSDA;
  // the slow mu- at all angles, those of them through the top opening; the mu- stopped in the lead
  std::vector<double> slowV(nT, 0.), slowV2(nT, 0.), csdaV(nT, 0.), slowAll(nT, 0.), slowTop(nT, 0.), stopped(nT, 0.);
  TRandom3 random(seed);
  // through the top opening of the walls: back along its direction from a random point of the target top face to the
  // bottom of the lead (as analysis/lead_target_stops.C)
  auto fromTop = [&](double ux, double uy, double uz) {
    const double xe = random.Uniform(-targetHalf, targetHalf), ye = random.Uniform(-targetHalf, targetHalf);
    const double back = (leadBottom - targetTop) / uz;   // uz < 0
    return std::abs(xe + ux * back) < wallInner && std::abs(ye + uy * back) < wallInner;
  };

  for (const auto& f : files) {
    TFile file(f.path);
    auto events = dynamic_cast<TTree*>(file.Get("events"));
    if (!events) continue;
    const int k = f.k, t = ts[k];
    int pdg = 0;
    double ekin = 0., costh = 0., phi = 0., x0 = 0., y0 = 0., z0 = 0., endX = 0., endY = 0., endZ = 0., endKe = 0.;
    events->SetBranchStatus("*", false);
    for (const char* b : {"pdg", "ekin", "costh", "phi", "x0", "y0", "z0", "end_x", "end_y", "end_z", "end_ke"})
      events->SetBranchStatus(b, true);
    events->SetBranchAddress("pdg", &pdg); events->SetBranchAddress("ekin", &ekin);
    events->SetBranchAddress("costh", &costh); events->SetBranchAddress("phi", &phi);
    events->SetBranchAddress("x0", &x0); events->SetBranchAddress("y0", &y0); events->SetBranchAddress("z0", &z0);
    events->SetBranchAddress("end_x", &endX); events->SetBranchAddress("end_y", &endY);
    events->SetBranchAddress("end_z", &endZ); events->SetBranchAddress("end_ke", &endKe);
    const double wBefore = 1e4 * 86400. / (f.area * liveAll);   // per m2 per day: before the lead, all files together
    const double wAfter = 1e4 * 86400. / (f.area * liveT[k]);   // after the slab: the files of its thickness
    const double half = 0.5 * t;
    const Long64_t n = events->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
      events->GetEntry(i);
      if (pdg != 13) continue;   // mu- only
      // before the lead
      const int b = countAll.FindFixBin(ekin);
      if (b >= 1 && b <= nE && edges[b] <= f.eMax * (1. + 1e-9)) {
        countAll.Fill(ekin);
        if (costh > 0.9) countVert.Fill(ekin);
      }
      if (ekin >= tLo && ekin < tHi) {   // with no lead it would stop in the target
        const double s = std::sqrt(std::max(0., 1. - costh * costh));
        slowAll[0] += wBefore;
        if (fromTop(s * std::cos(phi), s * std::sin(phi), -costh)) slowTop[0] += wBefore;
        angle[0]->Fill(std::acos(costh) * deg, wBefore);
        if (costh > 0.9) {
          slowV[0] += wBefore; slowV2[0] += wBefore * wBefore; csdaV[0] += wBefore;
          slice[0]->Fill(ekin, wBefore);
        }
      }
      // through the slab: CSDA for the same muon (straight), then what Geant4 did
      const double csdaDepth = lead.Range(ekin) / lead.rho * costh;   // cm
      if (csdaDepth < t) depthCsda[k]->Fill(csdaDepth, wAfter);
      if (costh > 0.9) {
        const double e = lead.After(ekin, t * lead.rho / costh);
        if (e >= tLo && e < tHi) csdaV[k] += wAfter;
      }
      if (endZ >= -half - 1e-3) {   // stopped in the slab
        if (endZ <= half + 1e-3) { depth[k]->Fill(half - endZ, wAfter); stopped[k] += wAfter; }
        continue;
      }
      if (endKe < tLo || endKe >= tHi) continue;   // left the lead too fast to stop in the target (or too slow)
      const double dx = endX - x0, dy = endY - y0, dz = endZ - z0, norm = std::sqrt(dx * dx + dy * dy + dz * dz);
      slowAll[k] += wAfter;
      if (fromTop(dx / norm, dy / norm, dz / norm)) slowTop[k] += wAfter;
      angle[k]->Fill(std::acos(-dz / norm) * deg, wAfter);
      if (costh > 0.9) { slowV[k] += wAfter; slowV2[k] += wAfter * wAfter; slice[k]->Fill(ekin, wAfter); }
    }
    printf("%s: %lld muons\n", f.path.Data(), n);
  }

  // the spectrum per MeV per m2 per day, and per unit of ln E (x E)
  const double area = files.front().area;
  auto perMeV = [&](const TH1D& counts, const char* name) {
    auto h = new TH1D(name, "", nE, edges.data());
    for (int b = 1; b <= nE; ++b) {
      if (exposure[b] <= 0.) continue;
      const double scale = 1e4 * 86400. / (area * exposure[b] * h->GetBinWidth(b));
      h->SetBinContent(b, counts.GetBinContent(b) * scale);
      h->SetBinError(b, counts.GetBinError(b) * scale);
    }
    return h;
  };
  auto perLn = [&](const TH1D* h, const char* name) {
    auto p = static_cast<TH1D*>(h->Clone(name));
    for (int b = 1; b <= nE; ++b) {
      p->SetBinContent(b, h->GetBinContent(b) * h->GetBinCenter(b));
      p->SetBinError(b, h->GetBinError(b) * h->GetBinCenter(b));
    }
    return p;
  };
  TH1D* jAll = perMeV(countAll, "jAll");
  TH1D* jVert = perMeV(countVert, "jVert");
  TH1D* lnAll = perLn(jAll, "lnAll");
  TH1D* lnVert = perLn(jVert, "lnVert");
  auto peak = [](const TH1D* h) { return h->GetBinCenter(h->GetMaximumBin()); };
  auto within = [](const TH1D* h, double fraction, double& lo, double& hi) {   // where h is within a fraction of its top
    const double top = h->GetMaximum();
    lo = hi = 0.;
    for (int b = 1; b <= h->GetNbinsX(); ++b)
      if (h->GetBinContent(b) >= fraction * top) {
        if (lo == 0.) lo = h->GetBinCenter(b);
        hi = h->GetBinCenter(b);
      }
  };
  // the slice of each thickness: vertical CSDA, the energies before the lead that leave it with tLo - tHi
  std::vector<double> e0Lo(nT), e0Hi(nT), e0Mid(nT);
  for (int k = 0; k < nT; ++k) {
    e0Lo[k] = ts[k] == 0 ? tLo : lead.Before(tLo, ts[k] * lead.rho);
    e0Hi[k] = ts[k] == 0 ? tHi : lead.Before(tHi, ts[k] * lead.rho);
    e0Mid[k] = 0.5 * (e0Lo[k] + e0Hi[k]);
  }
  auto depthAt = [&](int t) { return depth[iThick]->GetBinContent(std::min(t + 2, 50)); };   // the cm 1 - 2 cm below t
  std::map<int, std::pair<double, double>> model;   // analysis/lead_target_stops.C: mu- stops per day and error
  {
    std::ifstream csv(Form("%s/lead_target_stops.csv", out));
    std::string line;
    std::getline(csv, line);   // header
    while (std::getline(csv, line) && !line.empty()) {
      std::istringstream row(line);
      std::vector<double> v;
      for (std::string cell; std::getline(row, cell, ',');) v.push_back(std::atof(cell.c_str()));
      if (v.size() >= 3) model[int(v[0] + 0.5)] = {v[1], v[2]};
    }
  }

  // ---- numbers
  printf("\n=== Why the mu- stops in the target hardly change with the top lead (%s) ===\n", dir);
  printf("mu- only (they give the muonic X-rays), per m2 per day through a horizontal surface; before the lead: all %zu "
         "files (%.2f days); after a slab: its own files (", files.size(), liveAll / 86400.);
  for (int k = 1; k < nT; ++k) printf("%s%d cm %.2f d", k > 1 ? ", " : "", ts[k], liveT[k] / 86400.);
  printf(")\n");
  double lo, hi, loV, hiV;
  within(jAll, 0.95, lo, hi);
  within(jVert, 0.95, loV, hiV);
  const double topAll = jAll->GetMaximum();
  printf("\n1. the sea-level spectrum, all zenith angles (near vertical, cos > 0.9, in brackets)\n");
  printf("   per MeV:          highest at %.0f MeV (%.0f MeV), within 5%% of that from %.0f to %.0f MeV (%.0f to %.0f MeV);\n",
         peak(jAll), peak(jVert), lo, hi, loV, hiV);
  printf("                     at 1 GeV %.0f%%, at 2 GeV %.0f%%, at 5 GeV %.0f%% of the highest\n",
         100. * jAll->Interpolate(1000.) / topAll, 100. * jAll->Interpolate(2000.) / topAll,
         100. * jAll->Interpolate(5000.) / topAll);
  printf("   per unit of ln E: highest at %.2f GeV (%.2f GeV)%s\n", peak(lnAll) / 1000., peak(lnVert) / 1000.,
         lnAll->GetMaximumBin() == nE || lnVert->GetMaximumBin() == nE ? "; the top bin, so still rising at 5.7 GeV" : "");
  printf("\n2. leaving the lead with %.1f - %.1f MeV a vertical mu- stops in the 1 cm target (it gets through the %.2f g/cm2\n"
         "   of steel plate, TOF stations, Al plates and air, not through the 2.70 g/cm2 of Al); that slice of the spectrum:\n",
         tLo, tHi, aboveGcm2);
  printf(" lead   E0 of the slice    dE0     j(E0), cos > 0.9    j x dE0       mu- in the slice, cos > 0.9      Geant4\n");
  printf("  cm    (vertical, MeV)    (MeV)   per MeV/m2/day      per m2/day    CSDA (each muon)   Geant4        / CSDA\n");
  for (int k = 0; k < nT; ++k) {
    const double jv = jVert->Interpolate(e0Mid[k]), dE0 = e0Hi[k] - e0Lo[k];
    printf("  %3d   %6.1f - %6.1f   %5.2f   %8.1f            %8.0f      %8.0f           %6.0f +- %3.0f   %.3f\n", ts[k], e0Lo[k],
           e0Hi[k], dE0, jv, jv * dE0, csdaV[k], slowV[k], std::sqrt(slowV2[k]), slowV[k] / csdaV[k]);
  }
  const double jvRef = jVert->Interpolate(e0Mid[iRef]), dRef = e0Hi[iRef] - e0Lo[iRef];
  printf("\n3. relative to 10 cm of lead\n");
  printf(" lead   j(E0)   dE0    j x dE0   Geant4: mu- in   Geant4: mu- ending   slow mu- through    their share    lead_target_stops.C:\n");
  printf("  cm                             the slice        per cm of lead,      the top opening     through the    mu- stops in the\n");
  printf("                                 (cos > 0.9)      1-2 cm deeper        (all angles)        top opening    1 cm target\n");
  for (int k = 0; k < nT; ++k) {
    const double jv = jVert->Interpolate(e0Mid[k]), dE0 = e0Hi[k] - e0Lo[k];
    const bool hasModel = model.count(ts[k]) && model.count(10);
    printf("  %3d   %5.3f   %5.3f  %5.3f     %5.3f            %5.3f                %5.3f               %4.1f%%          %s\n", ts[k],
           jv / jvRef, dE0 / dRef, jv * dE0 / (jvRef * dRef), slowV[k] / slowV[iRef], depthAt(ts[k]) / depthAt(10),
           slowTop[k] / slowTop[iRef], 100. * slowTop[k] / slowAll[k],
           hasModel ? Form("%.3f", model[ts[k]].first / model[10].first) : "-");
  }
  printf("\n4. where the slowed muons go: mu- stopped in the lead itself, per m2 per day (Geant4):");
  for (int k = 1; k < nT; ++k) printf(" %d cm %.0f%s", ts[k], stopped[k], k + 1 < nT ? "," : "\n");
  printf("   per cm of lead in the 50 cm slab, Geant4 (straight CSDA ranges):");
  for (int d : {0, 2, 5, 10, 20, 30, 40, 49})
    printf(" %d-%d cm %.0f (%.0f)%s", d, d + 1, depth[iThick]->GetBinContent(d + 1), depthCsda[iThick]->GetBinContent(d + 1),
           d < 49 ? "," : "\n");
  printf("   slow mu- leaving the lead at all angles, per m2 per day:");
  for (int k = 0; k < nT; ++k) printf(" %d cm %.0f%s", ts[k], slowAll[k], k + 1 < nT ? "," : "\n");
  printf("\n5. absolute check of lead_target_stops.C under 10 cm of lead: the slow mu- that come through the top opening\n"
         "   onto the 0.04 m2 of the target are %.0f per day; lead_target_stops.C, which also follows the angle through the\n"
         "   layers below the lead and in the target, has %s mu- stops per day\n", 0.04 * slowTop[iRef],
         model.count(10) ? Form("%.1f", model[10].first) : "-");

  // ---- plots
  gStyle->SetOptStat(0);
  TCanvas c("c", "why the stops do not grow with the lead", 1800, 1000);
  c.Divide(3, 2);
  const int colour[7] = {kBlack, kBlue + 1, kAzure + 1, kGreen + 2, kOrange + 1, kRed + 1, kMagenta + 2};
  TLatex heading, text;
  heading.SetNDC(); heading.SetTextSize(0.05);
  text.SetTextSize(0.035);
  auto pad = [&](int p, bool logx) {
    c.cd(p);
    gPad->SetLogx(logx);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.1); gPad->SetBottomMargin(0.12);
  };
  auto frame = [](double x1, double y1, double x2, double y2, const char* titles) {
    auto f = gPad->DrawFrame(x1, y1, x2, y2, titles);
    f->GetXaxis()->SetTitleSize(0.045); f->GetYaxis()->SetTitleSize(0.045);
    f->GetXaxis()->SetTitleOffset(1.15); f->GetYaxis()->SetTitleOffset(1.45);
    return f;
  };
  auto styled = [](TH1* h, int col, int width, int style) {
    h->SetLineColor(col); h->SetLineWidth(width); h->SetLineStyle(style);
    return h;
  };
  auto legendAt = [](double x1, double y1, double x2, double y2, double size) {
    auto l = new TLegend(x1, y1, x2, y2);
    l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(size);
    return l;
  };

  // 1. the spectrum, per MeV and per unit of ln E
  pad(1, true);
  frame(eLow, 0., eHigh, 1.3, ";#mu^{-} kinetic energy at sea level E (MeV);relative to its highest value");
  TH1D* curves[4] = {jAll, lnAll, jVert, lnVert};
  const char* names1[4] = {"per MeV: dN/dE", "per unit of ln E: E dN/dE", "per MeV, near vertical",
                           "per unit of ln E, near vertical"};
  auto legend1 = legendAt(0.16, 0.7, 0.7, 0.89, 0.034);
  for (int q = 0; q < 4; ++q) {
    auto h = styled(static_cast<TH1D*>(curves[q]->Clone(Form("norm%d", q))), q % 2 ? kRed + 1 : kAzure + 2, q < 2 ? 3 : 2,
                    q < 2 ? 1 : 2);
    h->Scale(1. / h->GetMaximum());
    h->Draw("HIST SAME");
    legend1->AddEntry(h, names1[q], "l");
  }
  legend1->Draw();
  for (int k : {0, iRef, iThick}) {   // the slices that stop in the target under 0, 10 and 50 cm of lead
    auto line = new TLine(e0Mid[k], 0., e0Mid[k], 0.6);
    line->SetLineColor(kGray + 2); line->SetLineStyle(7); line->SetLineWidth(2);
    line->Draw();
    text.SetTextAlign(21);
    text.DrawLatex(e0Mid[k], 0.62, Form("%d cm", ts[k]));
  }
  text.SetTextAlign(11);
  text.DrawLatex(11., 0.05, "dashed: E_{0} that stops in the target under that much lead");
  heading.DrawLatex(0.14, 0.93, "#mu^{-} at sea level: the peak depends on the binning");

  // 2. the slices
  pad(2, false);
  double sMax = 0.;
  for (int k = 0; k < nT; ++k) {
    slice[k]->Scale(1. / slice[k]->GetBinWidth(1));
    sMax = std::max(sMax, slice[k]->GetMaximum());
  }
  frame(0., 0., 1100., 1.35 * sMax, ";#mu^{-} kinetic energy before the lead E_{0} (MeV);#mu^{-} per MeV (m^{-2} day^{-1})");
  auto legend2 = legendAt(0.4, 0.5, 0.97, 0.89, 0.036);
  legend2->SetHeader("near-vertical muons in it (m^{-2} day^{-1})");
  for (int k = 0; k < nT; ++k) {
    styled(slice[k], colour[k], 2, 1)->Draw("HIST SAME");
    legend2->AddEntry(slice[k], Form("%2d cm of lead: %.0f +/- %.0f", ts[k], slowV[k], std::sqrt(slowV2[k])), "l");
  }
  legend2->Draw();
  heading.DrawLatex(0.14, 0.93, "#mu^{-} stopping in the 1 cm target: one slice each");

  // 3. j, S and j x S
  pad(3, true);
  const double eRef = e0Mid[iRef], jRef = jVert->Interpolate(eRef), sRef = lead.DEdx(eRef);
  frame(40., 0., 3000., 2., ";#mu^{-} kinetic energy before the lead E_{0} (MeV);relative to the E_{0} of 10 cm of lead");
  auto gJ = new TGraph(), gS = new TGraph(), gJS = new TGraph();
  for (int b = 1; b <= nE; ++b) {
    const double e = jVert->GetBinCenter(b);
    if (e < 40. || e > 3000. || jVert->GetBinContent(b) <= 0.) continue;
    gJ->SetPoint(gJ->GetN(), e, jVert->GetBinContent(b) / jRef);
    gS->SetPoint(gS->GetN(), e, lead.DEdx(e) / sRef);
    gJS->SetPoint(gJS->GetN(), e, jVert->GetBinContent(b) * lead.DEdx(e) / (jRef * sRef));
  }
  gJ->SetLineColor(kAzure + 2); gJ->SetLineWidth(2);
  gS->SetLineColor(kGray + 2); gS->SetLineWidth(2); gS->SetLineStyle(2);
  gJS->SetLineColor(kBlack); gJS->SetLineWidth(3);
  for (auto g : {gJ, gS, gJS}) g->Draw("L SAME");
  auto gSlice = new TGraphErrors();
  text.SetTextAlign(21);
  for (int k = 0; k < nT; ++k) {
    gSlice->SetPoint(k, e0Mid[k], slowV[k] / slowV[iRef]);
    gSlice->SetPointError(k, 0., std::sqrt(slowV2[k]) / slowV[iRef]);
    text.DrawLatex(e0Mid[k], slowV[k] / slowV[iRef] + 0.07, Form("%d", ts[k]));
  }
  gSlice->SetMarkerStyle(20); gSlice->SetMarkerColor(kRed + 1); gSlice->SetLineColor(kRed + 1); gSlice->SetMarkerSize(1.2);
  gSlice->Draw("P SAME");
  auto legend3 = legendAt(0.16, 0.66, 0.9, 0.89, 0.034);
  legend3->AddEntry(gJ, "j(E_{0}): near-vertical muons per MeV", "l");
  legend3->AddEntry(gS, "S(E_{0}): energy loss in lead per g/cm^{2}", "l");
  legend3->AddEntry(gJS, "j S: muons per g/cm^{2} of range", "l");
  legend3->AddEntry(gSlice, "Geant4: muons in the slice (labels: cm of lead)", "p");
  legend3->Draw();
  heading.DrawLatex(0.14, 0.93, "More muons per MeV, but a narrower slice");

  // 4. stopping depths in the lead
  pad(4, false);
  double dMax = 0.;
  for (int k = 1; k < nT; ++k)
    for (int b = 1; b <= ts[k]; ++b) dMax = std::max({dMax, depth[k]->GetBinContent(b), depthCsda[k]->GetBinContent(b)});
  frame(0., 0., 50., 1.35 * dMax, ";depth in the lead (cm);#mu^{-} stopping per cm of lead (m^{-2} day^{-1})");
  for (int k = 1; k < nT; ++k) {   // a line through the bins, without the drop to zero at the bottom of the slab
    depth[k]->GetXaxis()->SetRange(1, ts[k]);
    styled(depth[k], colour[k], k == iThick ? 3 : 1, 1)->Draw("L SAME");
  }
  styled(depthCsda[iThick], kBlack, 2, 2)->Draw("L SAME");
  auto legend4 = legendAt(0.16, 0.7, 0.9, 0.89, 0.034);
  legend4->AddEntry(depth[iThick], "Geant4, 50 cm slab (all angles)", "l");
  legend4->AddEntry(depth[1], "Geant4, the thinner slabs (colours of panel 2)", "l");
  legend4->AddEntry(depthCsda[iThick], "straight CSDA ranges of the same muons", "l");
  legend4->Draw();
  heading.DrawLatex(0.14, 0.93, "Geant4: stops per cm of lead hardly change");

  // 5. directions of the slow muons leaving the lead
  pad(5, false);
  const int shown[4] = {0, iRef, indexOf(30), iThick};
  double aMax = 0.;
  for (int k : shown) {
    if (angle[k]->Integral() > 0.) angle[k]->Scale(1. / angle[k]->Integral());
    aMax = std::max(aMax, angle[k]->GetMaximum());
  }
  frame(0., 0., 90., 1.45 * aMax, ";zenith angle of the slow muons leaving the lead (deg);share per 2 deg");
  auto band = new TBox(edgeAngle, 0., cornerAngle, 1.45 * aMax);
  band->SetFillColor(kGray); band->SetFillStyle(1001);
  band->Draw();
  for (int k : shown) styled(angle[k], colour[k], 2, 1)->Draw("HIST SAME");
  auto legend5 = legendAt(0.5, 0.58, 0.97, 0.89, 0.034);
  legend5->SetHeader("share through the top opening");
  for (int k : shown) {
    const TString label = k == 0 ? TString("no lead (sky directions)") : TString::Format("%d cm of lead", ts[k]);
    legend5->AddEntry(angle[k], TString::Format("%s: %.0f%%", label.Data(), 100. * slowTop[k] / slowAll[k]), "l");
  }
  legend5->AddEntry(band, "edge to corner of the opening", "f");
  legend5->Draw();
  gPad->RedrawAxis();
  heading.DrawLatex(0.14, 0.93, "Slow muons leave the lead spread out in angle");

  // 6. all relative to 10 cm
  pad(6, false);
  frame(-2., 0.2, 52., 1.3, ";top lead thickness (cm);relative to 10 cm of lead");
  auto relative = [](int marker, int col, int line) {
    auto g = new TGraphErrors();
    g->SetMarkerStyle(marker); g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineStyle(line); g->SetLineWidth(2);
    g->SetMarkerSize(1.3);
    return g;
  };
  auto g6Slice = relative(20, kRed + 1, 1), g6Csda = relative(24, kGray + 2, 2), g6Depth = relative(21, kGreen + 2, 1),
       g6Top = relative(22, kOrange + 1, 1), g6Model = relative(33, kAzure + 2, 1);
  for (int k = 0; k < nT; ++k) {
    g6Slice->SetPoint(k, ts[k], slowV[k] / slowV[iRef]);
    g6Slice->SetPointError(k, 0., std::sqrt(slowV2[k]) / slowV[iRef]);
    g6Csda->SetPoint(k, ts[k], csdaV[k] / csdaV[iRef]);
    g6Depth->SetPoint(k, ts[k], depthAt(ts[k]) / depthAt(10));
    g6Top->SetPoint(k, ts[k], slowTop[k] / slowTop[iRef]);
    if (model.count(ts[k]) && model.count(10)) {
      const int p = g6Model->GetN();
      g6Model->SetPoint(p, ts[k], model[ts[k]].first / model[10].first);
      g6Model->SetPointError(p, 0., model[ts[k]].second / model[10].first);
    }
  }
  for (auto g : {g6Csda, g6Slice, g6Depth, g6Top, g6Model})
    if (g->GetN() > 0) g->Draw("PL SAME");
  auto legend6 = legendAt(0.16, 0.14, 0.9, 0.4, 0.034);
  legend6->AddEntry(g6Slice, "Geant4: muons in the slice (near vertical)", "lp");
  legend6->AddEntry(g6Csda, "CSDA: j dE_{0} of the same muons", "lp");
  legend6->AddEntry(g6Depth, "Geant4: stops per cm of lead, 1-2 cm deeper", "lp");
  legend6->AddEntry(g6Top, "Geant4: slow muons through the top opening", "lp");
  if (g6Model->GetN() > 0) legend6->AddEntry(g6Model, "lead_target_stops.C: stops in the target", "lp");
  legend6->Draw();
  heading.DrawLatex(0.14, 0.93, "All relative to 10 cm of lead");

  gSystem->mkdir(out, true);
  c.Print(Form("%s/lead_target_stops_check.pdf", out));
  c.SaveAs(Form("%s/lead_target_stops_check.png", out));
  printf("\nplots: %s/lead_target_stops_check.png (and .pdf)\n", out);
}
