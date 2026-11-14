// The steel study's tagged-beam spectrum fitted with element templates (2026-11-12, the user: study 2 of the suggestions
// for an element spectrometer: "Element-template fit prototype plus fake datasets: detection limit and false-alarm rate
// per element"; cases chosen: ideal templates, the simulation's own noise, a calibration mismatch).
//
// One day of data is a sum of templates, each with a free factor c (1 = the run's steel): each element's own muonic X-ray
// entries (its lines with their escapes, sums and Compton continuum, as the panel records them: day_<El>_promptx of
// analysis/beam_steel_trace_spectrum.root), the beam's other entries (capture gammas, neutrons, electrons), the X-rays of
// muons caught outside the plate, and the accidentals. The fit is linear least squares with the expected one-day Poisson
// variance (the counts are large) in bins of about half a resolution, eMin-eMax (20-7000 keV).
//   detection limit  the mass fraction one day sees at 3 sigma: 3 sigma(c) x the steel's fraction (the Fisher matrix with
//                    the element absent and every other factor free); beside it the best single line's, with the formula
//                    of beam_steel_trace_spectrum.C (S' / sqrt(S' + B) = 3)
//   false alarm      fake days (Poisson) with the element left out of the truth, fitted with every template: the share
//                    claiming it at >= 3 sigma (0.135% by chance alone). Also a steel of Fe, Cr, Ni, Mo, Mn alone: the
//                    share of days claiming any of the 8 others
// Cases:
//   a   ideal: the templates are the truth
//   b   the simulation's own noise (1e8 muons; one day is 1.14e10): each fake day's truth is the templates plus a fresh draw
//       of that noise (as another simulation of the same size would give). b1: the fit ignores it; b2: the fit adds it to
//       the variance. And sigma(c) with 1 to 10000 times the simulation: how much simulation a one-day fit needs
//   c   calibration: the fake days' energy scale (1 + gain) and resolution (1 + widen) times the templates' (0.1%, 3%).
//       c1: the fit ignores it; c2: a gain and a width in the fit as two more templates (the model's derivatives), one
//       linear step from the templates' calibration; c3: the gain and the width fitted to convergence (the templates
//       transformed anew at each step), the days fitted linearly around that point
// Outputs in outDir (default <dir>/analysis): template_fit.txt and .csv, and the plots (.png, .pdf)
//   template_fit_limits          per element: the 3 sigma one-day mass fraction, single line and template (a, b, c2)
//   template_fit_false_alarms    per element: the false-alarm share in each case
//   template_fit_mc_statistics   per element: sigma(c) with k times the simulation over the ideal one
//   template_fit_correlations    the fit's correlations at the run's steel (ideal templates)
// Run in the work folder after beam_steel_trace_spectrum.C (a minute):
//   root -l -b -q 'analysis/beam_steel_template_fit.C("output/czt/beam_steel_czt_trace")'
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMatrixDSym.h>
#include <TNamed.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace templatefit {
using Matrix = std::vector<std::vector<double>>;   // [bin][template]

struct Element {   // the real steel, mass percent (as beam_steel_trace_spectrum.C)
  int z;
  const char* name;
  double a;
  double percent;
};

// the covariance of a linear least-squares fit (the columns of A, weight W per bin), equilibrated before the inversion
bool Covariance(const Matrix& A, const std::vector<double>& W, Matrix& cov)
{
  const int nb = A.size(), nt = A.empty() ? 0 : A[0].size();
  Matrix F(nt, std::vector<double>(nt, 0.));
  for (int j = 0; j < nb; ++j)
    for (int t = 0; t < nt; ++t) {
      const double at = A[j][t] * W[j];
      if (at != 0.)
        for (int u = 0; u <= t; ++u) F[t][u] += at * A[j][u];
    }
  std::vector<double> d(nt, 0.);
  for (int t = 0; t < nt; ++t) d[t] = F[t][t] > 0. ? 1. / std::sqrt(F[t][t]) : 0.;
  TMatrixDSym m(nt);
  for (int t = 0; t < nt; ++t)
    for (int u = 0; u <= t; ++u) m(t, u) = m(u, t) = F[t][u] * d[t] * d[u];
  double det = 0.;
  m.Invert(&det);
  if (det == 0. || !std::isfinite(det)) return false;
  cov.assign(nt, std::vector<double>(nt, 0.));
  for (int t = 0; t < nt; ++t)
    for (int u = 0; u < nt; ++u) cov[t][u] = m(t, u) * d[t] * d[u];
  return true;
}
// one row of the fit as weights on the data: c_row = sum_j g_j data_j
std::vector<double> EstimatorRow(const Matrix& A, const std::vector<double>& W, const Matrix& cov, int row)
{
  std::vector<double> g(A.size(), 0.);
  for (size_t j = 0; j < A.size(); ++j) {
    double s = 0.;
    for (size_t u = 0; u < cov.size(); ++u) s += cov[row][u] * A[j][u];
    g[j] = s * W[j];
  }
  return g;
}
double GaussIn(double lo, double hi, double mu, double s)
{
  return 0.5 * (std::erf((hi - mu) / (M_SQRT2 * s)) - std::erf((lo - mu) / (M_SQRT2 * s)));
}
// a spectrum in fine bins (width w from 0) with its energy scale (1 + gain) times and its resolution (1 + widen) times the
// original's: the entries moved (linear within a bin), then spread by the extra Gaussian (sqrt((1 + widen)^2 - 1) sigma)
std::vector<double> Transformed(const std::vector<double>& f, double w, double gain, double widen, double s1MeV)
{
  const int n = f.size();
  std::vector<double> cum(n + 1, 0.), g(n, 0.);
  for (int k = 0; k < n; ++k) cum[k + 1] = cum[k] + f[k];
  auto C = [&](double x) {
    const double u = x / w;
    if (u <= 0.) return 0.;
    if (u >= n) return cum[n];
    const int k = int(u);
    return cum[k] + (u - k) * f[k];
  };
  for (int k = 0; k < n; ++k) g[k] = C((k + 1) * w / (1. + gain)) - C(k * w / (1. + gain));
  if (widen <= 0.) return g;
  std::vector<double> out(n, 0.);
  const double extra = std::sqrt((1. + widen) * (1. + widen) - 1.);
  for (int k = 0; k < n; ++k) {
    if (g[k] == 0.) continue;
    const double x = (k + 0.5) * w, s = extra * s1MeV * std::sqrt(std::max(x, 1.) / 1000.);
    if (s < 0.05 * w) { out[k] += g[k]; continue; }
    const int r = int(std::ceil(5. * s / w));
    const int m0 = std::max(0, k - r), m1 = std::min(n - 1, k + r);
    double norm = 0.;
    for (int m = m0; m <= m1; ++m) norm += GaussIn(m * w, (m + 1) * w, x, s);
    for (int m = m0; m <= m1; ++m) out[m] += g[k] * GaussIn(m * w, (m + 1) * w, x, s) / norm;
  }
  return out;
}
}  // namespace templatefit

void beam_steel_template_fit(const char* dir, const char* outDir = "", int toys = 10000, double eMin = 20., double eMax = 7000.,
                             double gain = 0.001, double widen = 0.03, double fracU = 1e-7, double fracTh = 1e-7, double fracK = 5e-5)
{
  using namespace templatefit;
  gStyle->SetOptStat(0);
  const std::string base = dir;
  const std::string out = std::string(outDir).empty() ? base + "/analysis" : std::string(outDir);
  gSystem->mkdir(out.c_str(), true);
  const std::string inName = base + "/analysis/beam_steel_trace_spectrum.root";
  std::unique_ptr<TFile> in(TFile::Open(inName.c_str()));
  if (!in || in->IsZombie()) { printf("cannot open %s: run analysis/beam_steel_trace_spectrum.C first\n", inName.c_str()); return; }
  double s1MeV = 6.36991, rate = 1.32e5, seconds = 86400., promptNs = 10.;
  if (auto s = dynamic_cast<TNamed*>(in->Get("settings"))) {
    double th = 0.;
    sscanf(s->GetTitle(), "sigma1MeV %lf threshold %lf promptNs %lf rate %lf seconds %lf", &s1MeV, &th, &promptNs, &rate, &seconds);
  }
  auto Res = [&](double e) { return s1MeV * std::sqrt(std::max(e, 1.) / 1000.); };

  // the real steel (as beam_steel_trace_spectrum.C)
  std::vector<Element> steel = {{26, "Fe", 55.845, 0.}, {24, "Cr", 51.996, 16.}, {28, "Ni", 58.693, 10.},
                                {42, "Mo", 95.95, 2.}, {25, "Mn", 54.938, 2.}, {14, "Si", 28.086, 0.75},
                                {7, "N", 14.007, 0.10}, {6, "C", 12.011, 0.08}, {15, "P", 30.974, 0.045},
                                {16, "S", 32.06, 0.03}, {19, "K", 39.098, 100. * fracK}, {90, "Th", 232.038, 100. * fracTh},
                                {92, "U", 238.029, 100. * fracU}};
  double rest = 100.;
  for (size_t i = 1; i < steel.size(); ++i) rest -= steel[i].percent;
  steel[0].percent = rest;
  const int nEl = steel.size();

  // the templates in the input's fine bins: the 13 elements' X-rays, the beam's other entries, other muons' X-rays, accidentals
  auto Get = [&](const std::string& name) -> TH1D* {
    auto h = dynamic_cast<TH1D*>(in->Get(name.c_str()));
    if (!h) printf("no %s in %s\n", name.c_str(), inName.c_str());
    return h;
  };
  TH1D* total = Get("view_prompt");
  if (!total) return;
  const int nFine = total->GetNbinsX();
  const double w = total->GetBinWidth(1);
  const int iBeam = nEl, iOther = nEl + 1, iAcc = nEl + 2;
  int nT = nEl + 3;
  std::vector<std::string> names, shortNames;
  for (const Element& el : steel) { names.push_back(el.name); shortNames.push_back(el.name); }
  names.push_back("the beam's other entries"); shortNames.push_back("beam");
  names.push_back("other muons' X-rays"); shortNames.push_back("other #mu^{-}");
  names.push_back("the accidentals"); shortNames.push_back("accid.");
  Matrix F(nT, std::vector<double>(nFine, 0.)), E2(nT, std::vector<double>(nFine, 0.));   // [template][fine bin]
  auto Take = [&](const std::string& xName, const std::string& allName, int t) {   // the X-rays to t, the rest to the beam
    TH1D* x = Get(xName);
    TH1D* all = Get(allName);
    if (!x || !all) return false;
    for (int k = 0; k < nFine; ++k) {
      const double cx = x->GetBinContent(k + 1), ex = std::pow(x->GetBinError(k + 1), 2);
      F[t][k] = cx; E2[t][k] = ex;
      F[iBeam][k] += all->GetBinContent(k + 1) - cx;
      E2[iBeam][k] += std::max(0., std::pow(all->GetBinError(k + 1), 2) - ex);   // a subset: the variances subtract
    }
    return true;
  };
  for (int i = 0; i < nEl; ++i)
    if (!Take(std::string("day_") + steel[i].name + "_promptx", std::string("day_") + steel[i].name + "_prompt", i)) return;
  if (!Take("day_other_promptx", "day_other_prompt", iOther)) return;
  TH1D* acc = Get("accidental_prompt");
  if (!acc) return;
  for (int k = 0; k < nFine; ++k) { F[iAcc][k] = acc->GetBinContent(k + 1); E2[iAcc][k] = std::pow(acc->GetBinError(k + 1), 2); }

  // the fit's bins: about half a resolution, whole fine bins, eMin-eMax
  std::vector<int> fitBin(nFine, -1);
  std::vector<double> lowEdge;
  for (double e = eMin; e < eMax - 1e-9;) {
    const double width = w * std::max(1., std::round(Res(e) / (2. * w)));
    for (int k = int(std::lround(e / w)); k < int(std::lround((e + width) / w)) && k < nFine; ++k) fitBin[k] = lowEdge.size();
    lowEdge.push_back(e);
    e += width;
  }
  const int nb = lowEdge.size();
  // a template without entries in the fit's range is left out (the other muons' X-rays are not the plate's: not flagged as
  // X-rays, they are in the beam's other entries)
  for (int t = nT - 1; t >= nEl; --t) {
    double s = 0.;
    for (int k = 0; k < nFine; ++k) if (fitBin[k] >= 0) s += F[t][k];
    if (s > 0.) continue;
    printf("template \"%s\" has no entries in %g-%g keV: left out\n", names[t].c_str(), eMin, eMax);
    F.erase(F.begin() + t); E2.erase(E2.begin() + t); names.erase(names.begin() + t); shortNames.erase(shortNames.begin() + t);
    --nT;
  }
  auto Merge = [&](const std::vector<double>& fine) {
    std::vector<double> m(nb, 0.);
    for (int k = 0; k < nFine; ++k) if (fitBin[k] >= 0) m[fitBin[k]] += fine[k];
    return m;
  };
  Matrix A(nb, std::vector<double>(nT, 0.)), Vt(nb, std::vector<double>(nT, 0.));   // contents, simulation variances
  for (int t = 0; t < nT; ++t) {
    const std::vector<double> m = Merge(F[t]), v = Merge(E2[t]);
    for (int j = 0; j < nb; ++j) { A[j][t] = m[j]; Vt[j][t] = v[j]; }
  }
  double maxDev = 0.;
  {
    std::vector<double> tot(nFine);
    for (int k = 0; k < nFine; ++k) tot[k] = total->GetBinContent(k + 1);
    const std::vector<double> mt = Merge(tot);
    for (int j = 0; j < nb; ++j) {
      double s = 0.;
      for (int t = 0; t < nT; ++t) s += A[j][t];
      if (mt[j] > 0.) maxDev = std::max(maxDev, std::fabs(s / mt[j] - 1.));
    }
  }
  // the calibration mismatch: every template's fake-day version; the model's derivatives for the gain and the width
  Matrix Ac(nb, std::vector<double>(nT, 0.));
  for (int t = 0; t < nT; ++t) {
    const std::vector<double> m = Merge(Transformed(F[t], w, gain, widen, s1MeV));
    for (int j = 0; j < nb; ++j) Ac[j][t] = m[j];
  }
  std::vector<double> fineTotal(nFine, 0.);
  for (int t = 0; t < nT; ++t)
    for (int k = 0; k < nFine; ++k) fineTotal[k] += F[t][k];
  const double dGain = 0.0005, dWiden = 0.01;   // the derivatives' steps (smaller than the mismatch: the fit is linearized)
  const std::vector<double> mUp = Merge(Transformed(fineTotal, w, dGain, 0., s1MeV)), mDown = Merge(Transformed(fineTotal, w, -dGain, 0., s1MeV));
  const std::vector<double> mWide = Merge(Transformed(fineTotal, w, 0., dWiden, s1MeV)), m0 = Merge(fineTotal);
  Matrix A2 = A;   // with the gain and the width, linearized at the templates' calibration (one step)
  for (int j = 0; j < nb; ++j) {
    A2[j].push_back((mUp[j] - mDown[j]) / (2. * dGain));
    A2[j].push_back((mWide[j] - m0[j]) / dWiden);
  }
  // c3: the gain and the width fitted to convergence on the fake days' expectation (each step transforms the templates
  // and the derivatives anew); the fake days are then fitted linearly around that point
  auto ModelAt = [&](double g, double wd, Matrix& M) {
    M.assign(nb, std::vector<double>(nT + 2, 0.));
    for (int t = 0; t < nT; ++t) {
      const std::vector<double> m = Merge(Transformed(F[t], w, g, wd, s1MeV));
      for (int j = 0; j < nb; ++j) M[j][t] = m[j];
    }
    const std::vector<double> up = Merge(Transformed(fineTotal, w, g + dGain, wd, s1MeV)), down = Merge(Transformed(fineTotal, w, g - dGain, wd, s1MeV));
    const std::vector<double> wide = Merge(Transformed(fineTotal, w, g, wd + dWiden, s1MeV)), here = Merge(Transformed(fineTotal, w, g, wd, s1MeV));
    for (int j = 0; j < nb; ++j) { M[j][nT] = (up[j] - down[j]) / (2. * dGain); M[j][nT + 1] = (wide[j] - here[j]) / dWiden; }
  };
  double gHat = 0., wHat = 0.;
  int steps = 0;
  Matrix A3;
  {
    std::vector<double> expect(nb, 0.), mu(nb, 0.);
    for (int j = 0; j < nb; ++j)
      for (int t = 0; t < nT; ++t) { expect[j] += Ac[j][t]; mu[j] += A[j][t]; }
    std::vector<double> W(nb);
    for (int j = 0; j < nb; ++j) W[j] = 1. / std::max(1., mu[j]);
    for (steps = 1; steps <= 8; ++steps) {
      ModelAt(gHat, wHat, A3);
      Matrix cov;
      if (!Covariance(A3, W, cov)) break;
      double dg = 0., dw = 0.;
      const std::vector<double> rg = EstimatorRow(A3, W, cov, nT), rw = EstimatorRow(A3, W, cov, nT + 1);
      for (int j = 0; j < nb; ++j) { dg += rg[j] * expect[j]; dw += rw[j] * expect[j]; }
      gHat += dg;
      wHat = std::max(0., wHat + dw);
      if (std::fabs(dg) < 1e-8 && std::fabs(dw) < 1e-6) break;
    }
    ModelAt(gHat, wHat, A3);
  }

  // the models: the run's steel, and without element X; their variances
  auto Without = [&](const Matrix& M, const std::vector<int>& gone) {   // the sum of the templates not in gone
    std::vector<double> s(nb, 0.);
    for (int j = 0; j < nb; ++j)
      for (int t = 0; t < nT; ++t)
        if (std::find(gone.begin(), gone.end(), t) == gone.end()) s[j] += M[j][t];
    return s;
  };
  auto Weights = [&](const std::vector<double>& mu, const std::vector<double>* vSim, double k) {
    std::vector<double> W(nb);
    for (int j = 0; j < nb; ++j) W[j] = 1. / std::max(1., mu[j] + (vSim ? (*vSim)[j] / k : 0.));
    return W;
  };
  const std::vector<double> mu0 = Without(A, {}), vSim0 = Without(Vt, {});
  Matrix cov0;
  if (!Covariance(A, Weights(mu0, nullptr, 1.), cov0)) { printf("the fit's matrix cannot be inverted\n"); return; }
  const std::vector<double> ks = {1., 3., 10., 30., 100., 300., 1000., 3000., 10000.};
  std::vector<double> sigA(nEl), sigC3(nEl), precision(nEl);
  Matrix sigB(nEl, std::vector<double>(ks.size()));
  for (int x = 0; x < nEl; ++x) {
    precision[x] = std::sqrt(cov0[x][x]);
    const std::vector<double> mu = Without(A, {x}), vs = Without(Vt, {x});
    Matrix cov;
    Covariance(A, Weights(mu, nullptr, 1.), cov);
    sigA[x] = std::sqrt(cov[x][x]);
    for (size_t q = 0; q < ks.size(); ++q) {
      Covariance(A, Weights(mu, &vs, ks[q]), cov);
      sigB[x][q] = std::sqrt(cov[x][x]);
    }
    Covariance(A3, Weights(mu, nullptr, 1.), cov);
    sigC3[x] = std::sqrt(cov[x][x]);
  }
  // the best single line (the trace's windows and formula): each element's emitted lines merged within half a resolution
  std::vector<double> lineE(nEl, 0.), lineF3(nEl, 0.);
  if (auto t = dynamic_cast<TTree*>(in->Get("lines"))) {
    int z = 0;
    double e = 0., pd = 0.;
    t->SetBranchAddress("z", &z); t->SetBranchAddress("e_keV", &e); t->SetBranchAddress("per_day", &pd);
    std::map<int, std::vector<std::pair<double, double>>> raw;
    for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); raw[z].push_back({e, pd}); }
    auto WindowSum = [&](const std::vector<double>& fine, double lo, double hi) {
      double s = 0.;
      for (int k = std::max(0, int(lo / w) - 1); k <= std::min(nFine - 1, int(hi / w) + 1); ++k) {
        const double c = (k + 0.5) * w;
        if (c >= lo && c <= hi) s += fine[k];
      }
      return s;
    };
    for (int i = 0; i < nEl; ++i) {
      auto& v = raw[steel[i].z];
      std::sort(v.begin(), v.end());
      std::vector<std::pair<double, double>> merged;
      for (const auto& [e0, p] : v)
        if (!merged.empty() && e0 - merged.back().first < 0.5 * Res(e0)) {
          auto& m = merged.back();
          m.first = (m.first * m.second + e0 * p) / (m.second + p);
          m.second += p;
        } else merged.push_back({e0, p});
      double best = -1.;
      for (const auto& [e0, p] : merged) {
        const double s = Res(e0), lo = e0 - 1.5 * s, hi = e0 + 1.5 * s;
        if (lo < eMin) continue;
        const double S = WindowSum(F[i], lo, hi), B = WindowSum(fineTotal, lo, hi) - S;
        if (S <= 0.) continue;
        const double zz = S / std::sqrt(S + B);
        if (zz > best) { best = zz; lineE[i] = e0; lineF3[i] = steel[i].percent / 100. * 0.5 * (9. + std::sqrt(81. + 36. * B)) / S; }
      }
    }
  }

  // the fake days
  TRandom3 rng(20261112);
  auto Poisson = [&](double m) { return m > 1000. ? rng.Gaus(m, std::sqrt(m)) : double(rng.Poisson(std::max(0., m))); };
  auto Dot = [](const std::vector<double>& g, const std::vector<double>& d) { double s = 0.; for (size_t j = 0; j < g.size(); ++j) s += g[j] * d[j]; return s; };
  enum { kA, kB1, kB2, kC1, kC2, kC3, kCases };
  const char* caseName[kCases] = {"a ideal", "b1 sim. noise, ignored", "b2 sim. noise, in the variance", "c1 calibration off, ignored",
                                  "c2 calibration, one linear step", "c3 calibration, fitted (iterated)"};
  Matrix claims(nEl, std::vector<double>(kCases, 0.)), zMean(nEl, std::vector<double>(kCases, 0.));
  double gainFit = 0., widenFit = 0.;
  std::vector<double> data(nb), truth(nb);
  // the fits' rows for the elements in `test`, with the truth `gone` removed
  struct Fit { std::vector<std::vector<double>> g; std::vector<double> s; };
  auto MakeFit = [&](const Matrix& M, const std::vector<double>& W, const std::vector<int>& test) {
    Fit f;
    Matrix cov;
    Covariance(M, W, cov);
    for (int x : test) { f.g.push_back(EstimatorRow(M, W, cov, x)); f.s.push_back(std::sqrt(cov[x][x])); }
    return f;
  };
  auto RunToys = [&](const std::vector<int>& gone, const std::vector<int>& test, std::vector<std::vector<double>>& count,
                     std::vector<std::vector<double>>& zsum, std::vector<double>& any, bool nuisances) {
    const std::vector<double> mu = Without(A, gone), vs = Without(Vt, gone), muC = Without(Ac, gone);
    const std::vector<double> W = Weights(mu, nullptr, 1.), Wb = Weights(mu, &vs, 1.);
    const Fit fa = MakeFit(A, W, test), fb = MakeFit(A, Wb, test), fc = MakeFit(A2, W, test), fd = MakeFit(A3, W, test);
    std::vector<std::vector<double>> gNuis;
    if (nuisances) {
      Matrix cov;
      Covariance(A2, W, cov);
      gNuis = {EstimatorRow(A2, W, cov, nT), EstimatorRow(A2, W, cov, nT + 1)};
    }
    count.assign(test.size(), std::vector<double>(kCases, 0.));
    zsum.assign(test.size(), std::vector<double>(kCases, 0.));
    any.assign(kCases, 0.);
    for (int toy = 0; toy < toys; ++toy) {
      bool hit[kCases] = {};
      auto Score = [&](int c, const Fit& f, size_t q) {
        const double z = Dot(f.g[q], data) / f.s[q];
        zsum[q][c] += z;
        if (z >= 3.) { count[q][c] += 1.; hit[c] = true; }
      };
      for (int j = 0; j < nb; ++j) data[j] = Poisson(mu[j]);   // a
      for (size_t q = 0; q < test.size(); ++q) Score(kA, fa, q);
      for (int j = 0; j < nb; ++j) data[j] = Poisson(std::max(0., mu[j] + std::sqrt(vs[j]) * rng.Gaus()));   // b
      for (size_t q = 0; q < test.size(); ++q) { Score(kB1, fa, q); Score(kB2, fb, q); }
      for (int j = 0; j < nb; ++j) data[j] = Poisson(muC[j]);   // c
      for (size_t q = 0; q < test.size(); ++q) { Score(kC1, fa, q); Score(kC2, fc, q); Score(kC3, fd, q); }
      if (nuisances) { gainFit += Dot(gNuis[0], data) / toys; widenFit += Dot(gNuis[1], data) / toys; }
      for (int c = 0; c < kCases; ++c) if (hit[c]) any[c] += 1.;
    }
  };
  for (int x = 0; x < nEl; ++x) {
    std::vector<std::vector<double>> cnt, zs;
    std::vector<double> any;
    RunToys({x}, {x}, cnt, zs, any, x == nEl - 3);   // the fitted gain and width once (K left out)
    for (int c = 0; c < kCases; ++c) { claims[x][c] = cnt[0][c] / toys; zMean[x][c] = zs[0][c] / toys; }
    printf("  fake days without %s: done\n", steel[x].name);
    fflush(stdout);
  }
  std::vector<int> minor;
  for (int i = 5; i < nEl; ++i) minor.push_back(i);   // Si, N, C, P, S, K, Th, U
  std::vector<std::vector<double>> mCnt, mZ;
  std::vector<double> mAny;
  RunToys(minor, minor, mCnt, mZ, mAny, false);

  // the numbers
  std::ostringstream txt, csv;
  txt << Form("Template fit of the tagged-beam spectrum (hits within %g ns of their muon), one day at %.3g mu-/s: %s\n", promptNs, rate, base.c_str())
      << Form("  %d templates: the 13 elements' own X-ray entries", nT);
  for (int t = nEl; t < nT; ++t) txt << ", " << names[t];
  txt << Form(" (their sum is the prompt view within %.1g)\n", maxDev)
      << Form("  %d bins of about half a resolution (sigma %g keV x sqrt(E / MeV)), %g-%g keV; least squares with the expected "
              "Poisson variance of one day\n", nb, s1MeV, eMin, eMax)
      << Form("  %d fake days per case; calibration mismatch: energy scale +%g%%, resolution +%g%% (derivative steps %g%%, %g%%);\n"
              "  fitted: c2 (one linear step) gain %+.4f%%, width %+.2f%% on average; c3 (iterated, %d steps) gain %+.5f%%, width %+.3f%%\n\n",
              toys, 100. * gain, 100. * widen, 100. * dGain, 100. * dWiden, 100. * gainFit, 100. * widenFit, steps, 100. * gHat, 100. * wHat);
  txt << "the mass fraction one day sees at 3 sigma (and the precision of c at the run's steel, ideal templates):\n"
         "  element  steel     sigma(c)  | single line (keV)        | template: a ideal  b sim. noise x1    x10   x100  x1000  | c3 calibration fitted\n";
  for (int x = 0; x < nEl; ++x) {
    const double f = steel[x].percent / 100.;
    txt << Form("  %-2s  %10.3g  %9.3g  | %9.3g (%8.2f)  | %9.3g  %9.3g %9.3g %9.3g %9.3g  | %9.3g\n", steel[x].name, f, precision[x],
                lineF3[x], lineE[x], 3. * sigA[x] * f, 3. * sigB[x][0] * f, 3. * sigB[x][2] * f, 3. * sigB[x][4] * f, 3. * sigB[x][6] * f,
                3. * sigC3[x] * f);
  }
  txt << Form("\nfalse alarms: fake days without the element (majors: a steel without it), claimed at >= 3 sigma (%.3g%% by "
              "chance), in %d days; the mean z beside\n  element", 100. * 0.5 * std::erfc(3. / M_SQRT2), toys);
  for (int c = 0; c < kCases; ++c) txt << Form(" | %-30s", caseName[c]);
  txt << "\n";
  for (int x = 0; x < nEl; ++x) {
    txt << Form("  %-7s", steel[x].name);
    for (int c = 0; c < kCases; ++c) txt << Form(" | %8.3g%% (z %+8.2f)%10s", 100. * claims[x][c], zMean[x][c], "");
    txt << "\n";
  }
  txt << "\na steel of Fe, Cr, Ni, Mo, Mn alone: the days claiming any of Si, N, C, P, S, K, Th, U at >= 3 sigma\n ";
  for (int c = 0; c < kCases; ++c) txt << Form(" | %s: %.3g%%", caseName[c], 100. * mAny[c] / toys);
  txt << "\n  per element (case c1, the calibration ignored):";
  for (size_t q = 0; q < minor.size(); ++q) txt << Form(" %s %.3g%% (z %+.1f)", steel[minor[q]].name, 100. * mCnt[q][kC1] / toys, mZ[q][kC1] / toys);
  txt << "\n";
  csv << "element,mass_fraction,sigma_c_run,single_line_keV,f3_single_line,f3_template_ideal,f3_sim_x1,f3_sim_x10,f3_sim_x100,"
         "f3_sim_x1000,f3_calibration_fitted,false_a,false_b1,false_b2,false_c1,false_c2,false_c3,z_a,z_b1,z_b2,z_c1,z_c2,z_c3\n";
  for (int x = 0; x < nEl; ++x) {
    const double f = steel[x].percent / 100.;
    csv << Form("%s,%.6g,%.5g,%.2f,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g", steel[x].name, f, precision[x], lineE[x], lineF3[x], 3. * sigA[x] * f,
                3. * sigB[x][0] * f, 3. * sigB[x][2] * f, 3. * sigB[x][4] * f, 3. * sigB[x][6] * f, 3. * sigC3[x] * f);
    for (int c = 0; c < kCases; ++c) csv << Form(",%.5g", claims[x][c]);
    for (int c = 0; c < kCases; ++c) csv << Form(",%.4g", zMean[x][c]);
    csv << "\n";
  }
  std::ofstream(out + "/template_fit.txt") << txt.str();
  std::ofstream(out + "/template_fit.csv") << csv.str();
  printf("%s", txt.str().c_str());

  // the plots
  auto Frame = [&](const char* name, const char* title, const char* yTitle, double lo, double hi) {   // axes; the title as text
    auto h = new TH1D(name, Form(";;%s", yTitle), nEl, 0., nEl);
    h->SetDirectory(nullptr);
    for (int i = 0; i < nEl; ++i) h->GetXaxis()->SetBinLabel(i + 1, steel[i].name);
    h->GetXaxis()->SetLabelSize(0.055);
    h->SetMinimum(lo); h->SetMaximum(hi);
    h->Draw("axis");
    TLatex t;
    t.SetNDC(); t.SetTextSize(0.036); t.SetTextAlign(22); t.SetTextFont(42);
    t.DrawLatex(0.53, 0.95, title);
    return h;
  };
  struct Set { std::vector<double> v; int style, color; double dx; std::string label; };
  auto Markers = [&](const std::vector<Set>& sets, TLegend* leg, double floor) {
    for (const Set& st : sets) {
      auto g = new TGraph();
      for (int i = 0; i < nEl; ++i) if (st.v[i] > 0. || floor > 0.) g->AddPoint(i + 0.5 + st.dx, std::max(st.v[i], floor));
      g->SetMarkerStyle(st.style); g->SetMarkerSize(1.6); g->SetMarkerColor(st.color);
      g->Draw("p same");
      leg->AddEntry(g, st.label.c_str(), "p");
    }
  };
  {
    std::vector<Set> sets(4);
    sets[0] = {{}, 24, kBlack, -0.27, "best single line"};
    sets[1] = {{}, 20, kRed + 1, -0.09, "template fit, ideal templates"};
    sets[2] = {{}, 22, kMagenta + 2, 0.09, "template fit, this simulation's noise in the variance"};
    sets[3] = {{}, 21, kAzure + 2, 0.27, "template fit, gain and width fitted (calibration off)"};
    double lo = 1e300, hi = 0.;
    for (int x = 0; x < nEl; ++x) {
      const double f = steel[x].percent / 100.;
      const double v[4] = {lineF3[x], 3. * sigA[x] * f, 3. * sigB[x][0] * f, 3. * sigC3[x] * f};
      for (int s = 0; s < 4; ++s) {
        sets[s].v.push_back(v[s]);
        if (v[s] > 0.) { lo = std::min(lo, v[s]); hi = std::max(hi, v[s]); }
      }
      lo = std::min(lo, f); hi = std::max(hi, f);
    }
    auto c = new TCanvas("c_limits", "", 1400, 850);
    c->SetLogy(); c->SetGridy();
    c->SetLeftMargin(0.1); c->SetRightMargin(0.03); c->SetTopMargin(0.1); c->SetBottomMargin(0.1);
    Frame("frame_limits", "The mass fraction one day sees at 3 sigma: single line against the template fit",
          "mass fraction", std::pow(10., std::floor(std::log10(lo)) - 1.3), std::pow(10., std::ceil(std::log10(hi)) + 0.3));
    auto leg = new TLegend(0.12, 0.13, 0.62, 0.35);
    leg->SetBorderSize(0); leg->SetFillStyle(1001); leg->SetFillColor(kWhite); leg->SetTextSize(0.03); leg->SetTextFont(42);
    Markers(sets, leg, 0.);
    TLine nl;
    nl.SetLineColor(kGreen + 3); nl.SetLineWidth(3);
    for (int x = 0; x < nEl; ++x) nl.DrawLine(x + 0.1, steel[x].percent / 100., x + 0.9, steel[x].percent / 100.);
    auto gn = new TGraph();
    gn->SetLineColor(kGreen + 3); gn->SetLineWidth(3);
    leg->AddEntry(gn, "the steel's own fraction (seen where a marker is below it)", "l");
    leg->Draw();
    c->Print((out + "/template_fit_limits.png").c_str());
    c->Print((out + "/template_fit_limits.pdf").c_str());
  }
  {
    const double floor = 0.3 / toys, nominal = 0.5 * std::erfc(3. / M_SQRT2);
    std::vector<Set> sets(kCases);
    const int style[kCases] = {20, 26, 22, 25, 28, 21};
    const int color[kCases] = {kBlack, kMagenta + 2, kMagenta + 2, kAzure + 2, kAzure + 2, kAzure + 2};
    for (int c = 0; c < kCases; ++c) {
      sets[c] = {{}, style[c], color[c], -0.3 + 0.12 * c, caseName[c]};
      for (int x = 0; x < nEl; ++x) sets[c].v.push_back(claims[x][c]);
    }
    auto c = new TCanvas("c_false", "", 1400, 850);
    c->SetLogy(); c->SetGridy();
    c->SetLeftMargin(0.1); c->SetRightMargin(0.03); c->SetTopMargin(0.1); c->SetBottomMargin(0.1);
    Frame("frame_false", Form("False alarms: the share of %d fake days without the element that claim it at 3 sigma", toys),
          "share of the days", 0.1 / toys, 3.);
    auto leg = new TLegend(0.12, 0.5, 0.95, 0.67);   // between the claims and the chance level
    leg->SetNColumns(2);
    leg->SetBorderSize(0); leg->SetFillStyle(1001); leg->SetFillColor(kWhite); leg->SetTextSize(0.028); leg->SetTextFont(42);
    Markers(sets, leg, floor);
    auto ln = new TLine(0., nominal, nEl, nominal);
    ln->SetLineStyle(2); ln->SetLineColor(kGray + 2);
    ln->Draw();
    leg->AddEntry(ln, "by chance alone (0.135%)", "l");
    TLatex t;
    t.SetTextSize(0.026); t.SetTextColor(kGray + 2); t.SetTextAlign(13);
    t.DrawLatex(nEl - 2.9, floor * 0.75, Form("bottom row: none in %d", toys));
    leg->Draw();
    c->Print((out + "/template_fit_false_alarms.png").c_str());
    c->Print((out + "/template_fit_false_alarms.pdf").c_str());
  }
  {
    auto c = new TCanvas("c_mc", "", 1300, 900);
    c->SetLogx(); c->SetLogy(); c->SetGridy();
    c->SetLeftMargin(0.11); c->SetRightMargin(0.03); c->SetTopMargin(0.1); c->SetBottomMargin(0.12);
    double hi = 1.;
    for (int x = 0; x < nEl; ++x) for (size_t q = 0; q < ks.size(); ++q) hi = std::max(hi, sigB[x][q] / sigA[x]);
    auto fr = new TH1D("frame_mc", ";simulated muons / the run's (1e8);#sigma(c) over that of ideal templates", 100, 0.7, 1.4e4);
    fr->SetDirectory(nullptr);
    fr->SetMinimum(0.9); fr->SetMaximum(1.5 * hi);
    fr->GetXaxis()->SetNoExponent(); fr->GetXaxis()->SetMoreLogLabels();
    fr->GetXaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleSize(0.045);
    fr->Draw("axis");
    TLatex t;
    t.SetNDC(); t.SetTextSize(0.036); t.SetTextAlign(22); t.SetTextFont(42);
    t.DrawLatex(0.55, 0.95, "How much simulation a one-day template fit needs: its noise against the day's own");
    // the 13 elements lie on one curve (every bin holds the simulation's weight, ~114 per muon): K drawn, the spread
    // given, and sqrt(1 + w / k) with w from the run's own simulation (k = 1)
    const int iK = nEl - 3;
    double spread = 0.;
    for (int x = 0; x < nEl; ++x)
      for (size_t q = 0; q < ks.size(); ++q) spread = std::max(spread, std::fabs(sigB[x][q] / sigA[x] / (sigB[iK][q] / sigA[iK]) - 1.));
    const double wEff = std::pow(sigB[iK][0] / sigA[iK], 2) - 1.;
    auto g = new TGraph(), ga = new TGraph();
    for (size_t q = 0; q < ks.size(); ++q) g->AddPoint(ks[q], sigB[iK][q] / sigA[iK]);
    for (double k = 0.8; k <= 1.2e4; k *= 1.1) ga->AddPoint(k, std::sqrt(1. + wEff / k));
    ga->SetLineStyle(2); ga->SetLineColor(kGray + 2); ga->SetLineWidth(2);
    ga->Draw("l same");
    g->SetLineColor(kRed + 1); g->SetMarkerColor(kRed + 1); g->SetMarkerStyle(20); g->SetMarkerSize(1.4); g->SetLineWidth(2);
    g->Draw("lp same");
    auto leg = new TLegend(0.3, 0.72, 0.97, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(1001); leg->SetFillColor(kWhite); leg->SetTextSize(0.032); leg->SetTextFont(42);
    leg->SetMargin(0.12);
    leg->AddEntry(g, Form("K (the other 12 elements within %.1f%%)", 100. * spread), "lp");
    leg->AddEntry(ga, Form("#sqrt{1 + %.0f / k}: the simulation's variance is %.0f times a day's", wEff, wEff), "l");
    leg->Draw();
    c->Print((out + "/template_fit_mc_statistics.png").c_str());
    c->Print((out + "/template_fit_mc_statistics.pdf").c_str());
  }
  {
    auto h = new TH2D("corr", "The template fit's correlations at the run's steel (ideal templates, one day);;", nT, 0., nT, nT, 0., nT);
    h->SetDirectory(nullptr);
    for (int t = 0; t < nT; ++t)
      for (int u = 0; u < nT; ++u) h->SetBinContent(t + 1, u + 1, cov0[t][u] / std::sqrt(cov0[t][t] * cov0[u][u]));
    for (int t = 0; t < nT; ++t) {
      h->GetXaxis()->SetBinLabel(t + 1, shortNames[t].c_str());
      h->GetYaxis()->SetBinLabel(t + 1, shortNames[t].c_str());
    }
    auto c = new TCanvas("c_corr", "", 1300, 1150);
    c->SetLeftMargin(0.12); c->SetRightMargin(0.13); c->SetTopMargin(0.08); c->SetBottomMargin(0.1);
    gStyle->SetPaintTextFormat(".2f");
    h->SetMinimum(-1.); h->SetMaximum(1.);
    h->GetXaxis()->SetLabelSize(0.035); h->GetYaxis()->SetLabelSize(0.035);
    h->SetMarkerSize(0.9);
    h->Draw("colz text");
    c->Print((out + "/template_fit_correlations.png").c_str());
    c->Print((out + "/template_fit_correlations.pdf").c_str());
  }
  printf("written into %s: template_fit.txt, .csv, template_fit_limits, _false_alarms, _mc_statistics, _correlations (.png, .pdf)\n",
         out.c_str());
}
