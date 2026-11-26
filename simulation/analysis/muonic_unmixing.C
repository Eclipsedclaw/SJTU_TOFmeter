// Unmixing a sample's muonic X-ray spectrum into the library's pure spectra (2026-11-25, the user's idea from
// fluorescence imaging: the pure components' spectra, a "phasor" picture of them, and the shares that rebuild the
// total; their choices: the prompt window, the phasor as a picture and a cross-check, the test samples Ni80Cr20,
// Fe50Ni50 and SS316).
//
// Templates: the pure plates' prompt spectra per muon (analysis/muonic_library_templates.C: 1e8 muons each, every
// plate tuned to the same mean stopping depth, 8.24 mm). A sample's spectrum of N muons is modelled as
//     mu_b = N x sum_i a_i T_ib       (b a bin, i an element, T_ib the pure plate's entries in b per muon)
// with a_i the share of the sample's muons that element i catches. The fit is Poisson maximum likelihood (iterated
// weighted least squares: weights 1 / (mu_b + the templates' own variance)), the shares free (negative allowed, so no
// pull toward zero), in bins of about half a resolution, eMin-eMax (20-7000 keV). N is known (a tagged beam counts the
// muons), so sum_i a_i = 1 is a check, not a constraint.
//   known fit    only the elements the sample holds
//   blind fit    every element whose pure plate gives at least minPerMuon (0.2) prompt entries per muon: not H, He
//                and Li (0.01-0.12: their K lines, 1.9-18.7 keV, are below the 20 keV threshold, and what they give is
//                mostly the decay electrons of the first 10 ns, alike for every element)
// The truth: the muons each element caught (cap_z), less the setup's own captures (the plastic, aluminium and lead
// around the plate catch ~0.014% of the muons; every pure template carries them, so they are not a sample's).
// Mass fractions from the shares: the capture law of Geant4 (Fermi-Teller: share ~ atoms x Z; O 0.56 Z, halogens
// 0.66 Z) turned around, mass_i ~ a_i A_i / Z_i.
// Checks:
//   precision      fits of 100 measurements of 1e5 muons, 100 of 1e6, 10 of 1e7 and the whole sample (~1e8)
//   honest errors  the pulls (fit - truth) / error of those measurements: their width ~1 means the errors are right,
//                  their centre the bias
//   false claims   the blind fit: the elements a sample does not hold that come out at 3 sigma or more
//   matrix effect  each element's X-ray entries inside the sample (its own muons) over those of its pure plate, line
//                  by line: the X-rays leave through a different material (density, absorption), so a pure template
//                  is somewhat off inside a mixture; at high statistics that bias beats the statistical error
//   self-test      each pure plate's even jobs fitted (blind) with every plate's odd jobs: its own share ~1, the others
//                  ~0 (the method and the templates' noise, no matrix effect)
//   phasor         a spectrum S(E) as one point (G, S): G = sum S(E) cos(phi) / sum S(E), S the same with sin, phi =
//                  2 pi n (E - eMin) / (phasorMax - eMin), harmonic n. A mixture's point is the average of its
//                  components' points weighted by their entries: on the line between two, inside the triangle of three.
//                  The entry fractions read from it (two components; three for SS316: Fe, Cr, Ni), turned into muon
//                  shares, against the full fit
// Outputs in outDir (default <dir>/analysis): muonic_unmixing.txt and .csv, and one plot per file (.png, .pdf):
//   unmixing_spectrum_<sample>   the sample's spectrum, the fitted sum and its components (known fit)
//   unmixing_fit_vs_truth        fitted share over the true share, every element of the three samples (whole samples)
//   unmixing_pulls               the pulls of the 1e5 and 1e6 measurements against a unit Gaussian
//   unmixing_precision           SS316: each element's relative error against the number of muons
//   unmixing_matrix              the line intensities inside the samples over the pure plates'
//   unmixing_blind               the blind fit of the whole samples: the absent elements' shares over their errors
//   unmixing_phasor              the phasor plane: every pure plate and the samples
//   unmixing_phasor_zoom         around the samples: their components' line and triangle
// Run in the work folder after the templates (compiled; a few minutes):
//   root -l -b -q 'analysis/muonic_unmixing.C+("output/spectrum_decode/muonic_library_czt")'
#include <TCanvas.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMatrixDSym.h>
#include <TNamed.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace unmix {
using Matrix = std::vector<std::vector<double>>;   // [bin][template]
const int kFine = 14000;                           // the templates' 0.5 keV bins
const double kFineKeV = 0.5;
const char* kSym[43] = {"n", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
                        "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
                        "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo"};
const double kMass[43] = {0., 1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.180, 22.990,
                          24.305, 26.982, 28.086, 30.974, 32.06, 35.45, 39.948, 39.098, 40.078, 44.956, 47.867, 50.942,
                          51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.630, 74.922, 78.971, 79.904,
                          83.798, 85.468, 87.62, 88.906, 91.224, 92.906, 95.95};   // standard atomic weights

double CaptureWeight(int z) { return (z == 8 ? 0.56 : (z == 9 || z == 17 || z == 35) ? 0.66 : 1.) * z; }
int ZOf(const std::string& s)
{
  for (int z = 1; z <= 42; ++z)
    if (s == kSym[z]) return z;
  if (s == "Pb") return 82;
  if (s == "Th") return 90;
  if (s == "U") return 92;
  return 0;
}
std::string Sym(int z) { return z >= 1 && z <= 42 ? std::string(kSym[z]) : std::string(Form("Z%d", z)); }

// a composition "Fe:rest,Cr:16,Ni:10" (mass percent) as Z -> percent, "rest" the balance
std::map<int, double> ParseComposition(const std::string& c)
{
  std::map<int, double> m;
  int restZ = 0;
  double sum = 0.;
  std::stringstream s(c);
  std::string part;
  while (std::getline(s, part, ',')) {
    const size_t colon = part.find(':');
    if (colon == std::string::npos) continue;
    const int z = ZOf(part.substr(0, colon));
    const std::string v = part.substr(colon + 1);
    if (z == 0) continue;
    if (v == "rest") { restZ = z; continue; }
    m[z] = std::atof(v.c_str());
    sum += m[z];
  }
  if (restZ) m[restZ] = 100. - sum;
  return m;
}

std::string Setting(const std::string& settings, const std::string& key)   // "key value" in a TNamed's title
{
  std::istringstream s(settings);
  std::string k, v;
  while (s >> k >> v)
    if (k == key) return v;
  return "";
}

std::vector<double> Contents(TFile& f, const std::string& name)
{
  std::vector<double> v(kFine, 0.);
  auto h = dynamic_cast<TH1*>(f.Get(name.c_str()));
  if (!h) { printf("no %s in %s\n", name.c_str(), f.GetName()); return v; }
  for (int k = 0; k < kFine && k < h->GetNbinsX(); ++k) v[k] = h->GetBinContent(k + 1);
  return v;
}
std::map<int, double> Captures(TFile& f, const std::string& name)   // cap_z -> muons
{
  std::map<int, double> m;
  if (auto h = dynamic_cast<TH1*>(f.Get(name.c_str())))
    for (int b = 1; b <= h->GetNbinsX(); ++b)
      if (h->GetBinContent(b) > 0.) m[b - 1] = h->GetBinContent(b);
  return m;
}
double Total(const std::map<int, double>& m)
{
  double s = 0.;
  for (const auto& [z, n] : m) s += n;
  return s;
}

// the covariance of a weighted least-squares fit, (A^T W A)^-1, equilibrated before the inversion (as
// analysis/beam_steel_template_fit.C)
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

struct Fit {
  std::vector<double> a, err;
  Matrix cov;
  double chi2 = 0.;
  int ndf = 0, iterations = 0;
  bool ok = false;
};

// mu_b = muons x sum_t a_t T_bt fitted to the counts d_b: iterated weighted least squares with weights
// 1 / (mu_b + sum_t (muons a_t)^2 VT_bt), the variance of the counts plus that of the templates (VT: their variance per
// muon^2). At convergence the shares solve the Poisson likelihood's equations (the templates' noise added to the
// variance); the covariance is (A^T W A)^-1 at the solution. floor keeps a weight finite where the model is ~0.
Fit FitShares(const Matrix& T, const Matrix& VT, const std::vector<double>& d, double muons, double floor = 0.01)
{
  Fit r;
  const int nb = d.size(), nt = nb > 0 ? T[0].size() : 0;
  if (nt == 0 || muons <= 0.) return r;
  Matrix A(nb, std::vector<double>(nt, 0.));
  for (int b = 0; b < nb; ++b)
    for (int t = 0; t < nt; ++t) A[b][t] = muons * T[b][t];
  std::vector<double> W(nb), a(nt, 0.);
  for (int b = 0; b < nb; ++b) W[b] = 1. / std::max(d[b], 1.);   // the start: the data's own variance
  Matrix cov;
  for (int it = 0; it < 60; ++it) {
    if (!Covariance(A, W, cov)) return r;
    std::vector<double> y(nt, 0.);
    for (int b = 0; b < nb; ++b) {
      const double wd = W[b] * d[b];
      if (wd != 0.)
        for (int t = 0; t < nt; ++t) y[t] += A[b][t] * wd;
    }
    double change = 0.;
    for (int t = 0; t < nt; ++t) {
      double s = 0.;
      for (int u = 0; u < nt; ++u) s += cov[t][u] * y[u];
      change = std::max(change, std::fabs(s - a[t]) / std::sqrt(cov[t][t]));
      a[t] = s;
    }
    for (int b = 0; b < nb; ++b) {
      double mu = 0., vt = 0.;
      for (int t = 0; t < nt; ++t) {
        mu += A[b][t] * a[t];
        vt += muons * a[t] * muons * a[t] * VT[b][t];
      }
      W[b] = 1. / (std::max(mu, floor) + vt);
    }
    r.iterations = it + 1;
    if (it > 0 && change < 1e-4) break;
  }
  if (!Covariance(A, W, cov)) return r;
  r.a = a;
  r.cov = cov;
  r.err.resize(nt);
  for (int t = 0; t < nt; ++t) r.err[t] = std::sqrt(cov[t][t]);
  int used = 0;
  for (int b = 0; b < nb; ++b) {
    double mu = 0.;
    for (int t = 0; t < nt; ++t) mu += A[b][t] * a[t];
    if (mu <= 0. && d[b] <= 0.) continue;
    r.chi2 += (d[b] - mu) * (d[b] - mu) * W[b];
    ++used;
  }
  r.ndf = used - nt;
  r.ok = true;
  return r;
}

struct Phase {   // a phasor, and the entries in its range
  double g = 0., s = 0., entries = 0.;
};
}  // namespace unmix

void muonic_unmixing(const char* dir = "output/spectrum_decode/muonic_library_czt",
                     const char* table = "data/muonic_library_materials.txt", const char* outDir = "", double eMin = 20.,
                     double eMax = 7000., double phasorMax = 3500., int harmonic = 1, double minPerMuon = 0.2)
{
  using namespace unmix;
  gStyle->SetOptStat(0);
  const std::string base = dir, tdir = base + "/analysis/templates";
  const std::string out = std::string(outDir).empty() ? base + "/analysis" : std::string(outDir);
  gSystem->mkdir(out.c_str(), true);

  // the library's materials (data/muonic_library_materials.txt: name, id, composition, ...)
  struct Row { std::string name, composition, band; int id; };
  std::vector<Row> rows;
  {
    std::ifstream in(table);
    if (!in) { printf("cannot read %s\n", table); return; }
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream s(line);
      Row r;
      if (!(s >> r.name >> r.id >> r.composition)) continue;
      r.band = (r.id < 10 ? "0" : "") + std::to_string(r.id) + "_" + r.name;
      rows.push_back(r);
    }
  }

  // the pure plates: their spectra (all, two halves), X-ray entries of their own muons, muons by cap_z
  struct Pure {
    int z = 0;
    std::string band;
    double muons = 0., halfMuons[2] = {0., 0.};
    std::vector<double> all, half[2], ownX;
    std::map<int, double> capture;
  };
  std::map<int, Pure> pure;
  std::string settings0;
  double s1MeV = 6.36991, threshold = 20., promptNs = 10.;
  auto CheckSettings = [&](const std::string& settings, const std::string& what) {   // every file smeared alike
    const std::string mine = Setting(settings, "sigma1MeV") + " " + Setting(settings, "threshold") + " " + Setting(settings, "promptNs");
    if (settings0.empty()) {
      settings0 = mine;
      s1MeV = std::atof(Setting(settings, "sigma1MeV").c_str());
      threshold = std::atof(Setting(settings, "threshold").c_str());
      promptNs = std::atof(Setting(settings, "promptNs").c_str());
      return true;
    }
    if (mine == settings0) return true;
    printf("%s: smearing, threshold, window \"%s\" differ from the others' \"%s\": stopping\n", what.c_str(), mine.c_str(), settings0.c_str());
    return false;
  };
  for (const Row& r : rows) {
    const auto comp = ParseComposition(r.composition);
    if (comp.size() != 1) continue;
    const int z = comp.begin()->first;
    TFile f((tdir + "/" + r.band + ".root").c_str());
    if (f.IsZombie()) { printf("no templates of %s in %s: run analysis/muonic_library_templates.C first\n", r.band.c_str(), tdir.c_str()); return; }
    auto st = dynamic_cast<TNamed*>(f.Get("settings"));
    if (!st || !CheckSettings(st->GetTitle(), r.band)) return;
    Pure p;
    p.z = z;
    p.band = r.band;
    p.all = Contents(f, "prompt");
    p.half[0] = Contents(f, "prompt_half0");
    p.half[1] = Contents(f, "prompt_half1");
    p.ownX = Contents(f, Form("prompt_x_z%d", z));
    p.capture = Captures(f, "capture_z");
    p.muons = Total(p.capture);
    p.halfMuons[0] = Total(Captures(f, "capture_z_half0"));
    p.halfMuons[1] = Total(Captures(f, "capture_z_half1"));
    pure[z] = p;
  }
  printf("%zu pure plates read from %s (sigma %g keV at 1 MeV, threshold %g keV, prompt %g ns)\n", pure.size(), tdir.c_str(), s1MeV, threshold, promptNs);

  // the fit's bins: about half a resolution, whole fine bins, eMin-eMax (as analysis/beam_steel_template_fit.C)
  auto Res = [&](double e) { return s1MeV * std::sqrt(std::max(e, 1.) / 1000.); };
  std::vector<int> fitBin(kFine, -1);
  std::vector<double> edges;
  for (double e = eMin; e < eMax - 1e-9;) {
    const double width = kFineKeV * std::max(1., std::round(Res(e) / (2. * kFineKeV)));
    for (int k = int(std::lround(e / kFineKeV)); k < int(std::lround((e + width) / kFineKeV)) && k < kFine; ++k) fitBin[k] = edges.size();
    edges.push_back(e);
    e += width;
  }
  const int nb = edges.size();
  edges.push_back(std::min(eMax, edges.back() + kFineKeV * std::max(1., std::round(Res(edges.back()) / (2. * kFineKeV)))));
  auto Merge = [&](const std::vector<double>& fine) {
    std::vector<double> m(nb, 0.);
    for (int k = 0; k < kFine; ++k)
      if (fitBin[k] >= 0) m[fitBin[k]] += fine[k];
    return m;
  };
  auto InFit = [&](const std::vector<double>& fine) {
    double s = 0.;
    for (int k = 0; k < kFine; ++k)
      if (fitBin[k] >= 0) s += fine[k];
    return s;
  };

  // the elements a fit can see: at least minPerMuon prompt entries per muon in eMin-eMax
  std::vector<int> visible, invisible;
  std::map<int, double> perMuon;
  for (auto& [z, p] : pure) {
    perMuon[z] = InFit(p.all) / p.muons;
    (perMuon[z] >= minPerMuon ? visible : invisible).push_back(z);
  }
  // the setup's own captures: the share of each pure plate's muons caught by another element (or none, cap_z 0)
  std::map<int, std::map<int, double>> setup;
  for (auto& [z, p] : pure)
    for (const auto& [cz, n] : p.capture)
      if (cz != z) setup[z][cz] = n / p.muons;

  // templates per muon for the elements zs: half -1 the whole plates, 0 or 1 one half
  auto Templates = [&](const std::vector<int>& zs, int half, Matrix& T, Matrix& VT) {
    T.assign(nb, std::vector<double>(zs.size(), 0.));
    VT.assign(nb, std::vector<double>(zs.size(), 0.));
    for (size_t t = 0; t < zs.size(); ++t) {
      const Pure& p = pure.at(zs[t]);
      const double n = half < 0 ? p.muons : p.halfMuons[half];
      const std::vector<double> m = Merge(half < 0 ? p.all : p.half[half]);
      for (int b = 0; b < nb; ++b) { T[b][t] = m[b] / n; VT[b][t] = m[b] / (n * n); }
    }
  };
  auto PhasorOf = [&](const std::vector<double>& fine) {
    Phase ph;
    for (int k = 0; k < kFine; ++k) {
      const double e = (k + 0.5) * kFineKeV;
      if (e < eMin || e >= phasorMax || fine[k] == 0.) continue;
      const double phi = 2. * M_PI * harmonic * (e - eMin) / (phasorMax - eMin);
      ph.g += fine[k] * std::cos(phi);
      ph.s += fine[k] * std::sin(phi);
      ph.entries += fine[k];
    }
    if (ph.entries > 0.) { ph.g /= ph.entries; ph.s /= ph.entries; }
    return ph;
  };
  std::map<int, Phase> purePhase;
  for (int z : visible) purePhase[z] = PhasorOf(pure[z].all);

  // the samples: the two test alloys of the library and SS316 (the steel trace run without its U, Th, K muons)
  struct Measurement {   // one measurement: its counts in the fit's bins, its muons by cap_z, its phasor
    std::vector<double> d;
    std::map<int, double> capture;
    double muons = 0.;
    Phase phase;
  };
  struct Sample {
    std::string file, label;
    std::vector<int> zs;              // the elements it holds
    std::map<int, double> mass;       // their mass percent
    Measurement whole;
    std::map<int, std::vector<double>> ownX;   // each element's X-ray entries (fine bins)
    std::vector<Measurement> sets[3];          // measurements of 1e5, 1e6, 1e7 muons
  };
  std::vector<Sample> samples;
  for (const auto& [file, label] : std::vector<std::pair<std::string, std::string>>{{"101_Ni80Cr20", "Ni80Cr20"}, {"102_Fe50Ni50", "Fe50Ni50"}, {"SS316", "SS316"}}) {
    TFile f((tdir + "/" + file + ".root").c_str());
    if (f.IsZombie()) { printf("no templates of the sample %s in %s: skipped\n", file.c_str(), tdir.c_str()); continue; }
    auto st = dynamic_cast<TNamed*>(f.Get("settings"));
    if (!st || !CheckSettings(st->GetTitle(), file)) return;
    Sample s;
    s.file = file;
    s.label = label;
    s.mass = ParseComposition(Setting(st->GetTitle(), "composition"));
    {
      std::stringstream ds(Setting(st->GetTitle(), "drop"));
      std::string z;
      while (std::getline(ds, z, ',')) s.mass.erase(std::atoi(z.c_str()));
      double sum = 0.;
      for (const auto& [z2, m] : s.mass) sum += m;
      for (auto& [z2, m] : s.mass) m *= 100. / sum;
    }
    for (const auto& [z, m] : s.mass) {
      if (!pure.count(z) || perMuon[z] < minPerMuon) { printf("%s: %s has no visible pure plate: left out of its fit\n", label.c_str(), Sym(z).c_str()); continue; }
      s.zs.push_back(z);
    }
    std::sort(s.zs.begin(), s.zs.end(), [&](int x, int y) { return s.mass[x] > s.mass[y]; });
    const std::vector<double> all = Contents(f, "prompt");
    s.whole.d = Merge(all);
    s.whole.capture = Captures(f, "capture_z");
    s.whole.muons = Total(s.whole.capture);
    s.whole.phase = PhasorOf(all);
    for (int z : s.zs) s.ownX[z] = Contents(f, Form("prompt_x_z%d", z));
    for (int set = 0; set < 2; ++set) {
      const char* prefix = set == 0 ? "slice5" : "slice6";
      const int n = std::atoi(Setting(st->GetTitle(), set == 0 ? "slices5" : "slices6").c_str());
      auto capz = dynamic_cast<TH2*>(f.Get(Form("%s_capz", prefix)));
      if (!capz) continue;
      for (int j = 0; j < n; ++j) {
        Measurement m;
        const std::vector<double> fine = Contents(f, Form("%s_%d", prefix, j));
        m.d = Merge(fine);
        m.phase = PhasorOf(fine);
        for (int zb = 1; zb <= capz->GetNbinsY(); ++zb)
          if (capz->GetBinContent(j + 1, zb) > 0.) m.capture[zb - 1] = capz->GetBinContent(j + 1, zb);
        m.muons = Total(m.capture);
        if (set == 1 && j / 10 < 10) {   // ten 1e6 measurements make one of 1e7
          if (j % 10 == 0) s.sets[2].push_back(Measurement());
          Measurement& big = s.sets[2].back();
          if (big.d.empty()) big.d.assign(nb, 0.);
          for (int b = 0; b < nb; ++b) big.d[b] += m.d[b];
          for (const auto& [z, c] : m.capture) big.capture[z] += c;
          big.muons += m.muons;
          big.phase.g += m.phase.g * m.phase.entries;   // entry-weighted: divided below
          big.phase.s += m.phase.s * m.phase.entries;
          big.phase.entries += m.phase.entries;
        }
        s.sets[set].push_back(m);
      }
    }
    if (s.sets[1].size() % 10 != 0 && !s.sets[2].empty()) s.sets[2].pop_back();   // an incomplete 1e7 measurement
    for (Measurement& m : s.sets[2]) { m.phase.g /= m.phase.entries; m.phase.s /= m.phase.entries; }
    printf("%s: %.0f muons, elements", label.c_str(), s.whole.muons);
    for (int z : s.zs) printf(" %s %.4g%%", Sym(z).c_str(), s.mass[z]);
    printf("; measurements: %zu of 1e5, %zu of 1e6, %zu of 1e7\n", s.sets[0].size(), s.sets[1].size(), s.sets[2].size());
    samples.push_back(s);
  }
  if (samples.empty()) { printf("no sample to unmix\n"); return; }
  const char* setName[3] = {"1e5", "1e6", "1e7"};

  // the truth of a measurement: the share of its muons each of the sample's elements caught, less the setup's
  // captures its templates carry (per element of the sample, weighted by its muons)
  auto Truth = [&](const Sample& s, const Measurement& m) {
    std::map<int, double> f;
    double fTotal = 0.;
    for (int z : s.zs) {
      const double w = m.capture.count(z) ? m.capture.at(z) / m.muons : 0.;
      if (setup.count(z))
        for (const auto& [cz, share] : setup.at(z)) f[cz] += w * share;
    }
    for (const auto& [cz, v] : f) fTotal += v;
    std::vector<double> t;
    for (int z : s.zs) {
      const double n = m.capture.count(z) ? m.capture.at(z) : 0.;
      t.push_back((n - m.muons * (f.count(z) ? f[z] : 0.)) / (m.muons * (1. - fTotal)));
    }
    return t;
  };
  // mass fractions from shares (the capture law turned around) and their covariance
  auto MassFractions = [&](const std::vector<int>& zs, const Fit& fit, std::vector<double>& w, std::vector<double>& werr) {
    const int n = zs.size();
    std::vector<double> g(n);
    double U = 0.;
    for (int i = 0; i < n; ++i) { g[i] = kMass[zs[i]] / CaptureWeight(zs[i]); U += fit.a[i] * g[i]; }
    w.assign(n, 0.);
    werr.assign(n, 0.);
    for (int i = 0; i < n; ++i) w[i] = fit.a[i] * g[i] / U;
    for (int i = 0; i < n; ++i) {
      double v = 0.;
      for (int k = 0; k < n; ++k)
        for (int l = 0; l < n; ++l)
          v += g[k] * ((i == k) - w[i]) / U * fit.cov[k][l] * g[l] * ((i == l) - w[i]) / U;
      werr[i] = std::sqrt(std::max(v, 0.));
    }
  };

  // the self-test: every visible plate's even jobs, fitted with every visible plate's odd jobs
  struct SelfTest { int z; double own, ownErr, worst; int worstZ, claims; bool ok; };
  std::vector<SelfTest> self;
  {
    Matrix T0, VT0;
    Templates(visible, 0, T0, VT0);
    for (size_t i = 0; i < visible.size(); ++i) {
      const int z = visible[i];
      const Fit fit = FitShares(T0, VT0, Merge(pure[z].half[1]), pure[z].halfMuons[1]);
      SelfTest st{z, 0., 0., 0., 0, 0, fit.ok};
      if (fit.ok) {
        st.own = fit.a[i];
        st.ownErr = fit.err[i];
        for (size_t q = 0; q < visible.size(); ++q) {
          if (q == i) continue;
          const double sig = fit.a[q] / fit.err[q];
          if (std::fabs(sig) > std::fabs(st.worst)) { st.worst = sig; st.worstZ = visible[q]; }
          if (sig >= 3.) ++st.claims;
        }
      }
      self.push_back(st);
    }
  }
  printf("self-test done (%zu plates)\n", self.size());

  // the fits of the samples
  struct Result { double muons; std::vector<double> truth, a, err; Phase phase; };   // one measurement, known fit
  struct SampleFits {
    Fit known, blind;
    std::vector<double> truth, mass, massErr;
    std::vector<Result> sets[3];
    int blindClaims1e6 = 0, blindSlices1e6 = 0;   // the 1e6 measurements with an absent element at >= 3 sigma
  };
  std::vector<SampleFits> fits(samples.size());
  Matrix TB, VTB;
  Templates(visible, -1, TB, VTB);
  for (size_t si = 0; si < samples.size(); ++si) {
    const Sample& s = samples[si];
    SampleFits& sf = fits[si];
    Matrix T, VT;
    Templates(s.zs, -1, T, VT);
    sf.known = FitShares(T, VT, s.whole.d, s.whole.muons);
    sf.blind = FitShares(TB, VTB, s.whole.d, s.whole.muons);
    sf.truth = Truth(s, s.whole);
    if (!sf.known.ok) { printf("%s: the known fit failed\n", s.label.c_str()); continue; }
    MassFractions(s.zs, sf.known, sf.mass, sf.massErr);
    for (int set = 0; set < 3; ++set)
      for (const Measurement& m : s.sets[set]) {
        const Fit fit = FitShares(T, VT, m.d, m.muons);
        if (!fit.ok) continue;
        sf.sets[set].push_back({m.muons, Truth(s, m), fit.a, fit.err, m.phase});
      }
    std::set<int> held(s.zs.begin(), s.zs.end());
    for (const Measurement& m : s.sets[1]) {
      const Fit fit = FitShares(TB, VTB, m.d, m.muons);
      if (!fit.ok) continue;
      ++sf.blindSlices1e6;
      for (size_t q = 0; q < visible.size(); ++q)
        if (!held.count(visible[q]) && fit.a[q] / fit.err[q] >= 3.) { ++sf.blindClaims1e6; break; }
    }
    printf("%s fitted (%d iterations; %zu + %zu + %zu measurements)\n", s.label.c_str(), sf.known.iterations, sf.sets[0].size(), sf.sets[1].size(), sf.sets[2].size());
  }

  // the matrix effect, line by line: an element's strongest lines in its pure plate's X-ray entries (local maxima of
  // the fit's bins, at least 3% of the highest, 3 resolutions apart, up to 6, below phasorMax), each in +-1.5 sigma
  struct LineRatio { std::string sample; int z; double e, ratio, err; };
  std::vector<LineRatio> ratios;
  auto Lines = [&](int z) {
    const std::vector<double> m = Merge(pure[z].ownX);
    double top = 0.;
    for (int b = 0; b < nb; ++b) top = std::max(top, m[b]);
    std::vector<std::pair<double, double>> peaks;   // content, energy
    for (int b = 1; b + 1 < nb; ++b)
      if (m[b] >= 0.03 * top && m[b] > m[b - 1] && m[b] >= m[b + 1] && edges[b] < phasorMax)
        peaks.push_back({m[b], 0.5 * (edges[b] + edges[b + 1])});
    std::sort(peaks.rbegin(), peaks.rend());
    std::vector<double> lines;
    for (const auto& [c, e] : peaks) {
      bool near = false;
      for (double l : lines) near = near || std::fabs(l - e) < 3. * Res(e);
      if (!near) lines.push_back(e);
      if (lines.size() == 6) break;
    }
    std::sort(lines.begin(), lines.end());
    return lines;
  };
  for (size_t si = 0; si < samples.size(); ++si) {
    const Sample& s = samples[si];
    for (int z : s.zs) {
      const double nS = s.whole.capture.count(z) ? s.whole.capture.at(z) : 0., nP = pure[z].capture.count(z) ? pure[z].capture.at(z) : 0.;
      if (nS <= 0. || nP <= 0.) continue;
      for (double e : Lines(z)) {
        const int k0 = std::max(0, int((e - 1.5 * Res(e)) / kFineKeV)), k1 = std::min(kFine - 1, int((e + 1.5 * Res(e)) / kFineKeV));
        double cs = 0., cp = 0.;
        for (int k = k0; k <= k1; ++k) { cs += s.ownX.at(z)[k]; cp += pure[z].ownX[k]; }
        if (cs < 100. || cp < 100.) continue;
        const double r = (cs / nS) / (cp / nP);
        ratios.push_back({s.label, z, e, r, r * std::sqrt(1. / cs + 1. / cp)});
      }
    }
  }

  // the phasor cross-check: two components (the alloys) and three (SS316: Fe, Cr, Ni; its other elements ignored)
  auto PerMuonInPhasor = [&](int z) { return purePhase[z].entries / pure[z].muons; };
  auto BinaryShare = [&](const Phase& m, int z1, int z2, double& f1, double& off) {   // the muon share of z1
    const Phase &p1 = purePhase[z1], &p2 = purePhase[z2];
    const double dx = p1.g - p2.g, dy = p1.s - p2.s, l2 = dx * dx + dy * dy;
    f1 = ((m.g - p2.g) * dx + (m.s - p2.s) * dy) / l2;
    off = ((m.g - p2.g) * dy - (m.s - p2.s) * dx) / std::sqrt(l2);
    const double u1 = f1 / PerMuonInPhasor(z1), u2 = (1. - f1) / PerMuonInPhasor(z2);
    return u1 / (u1 + u2);
  };
  auto TernaryShares = [&](const Phase& m, const int z[3], double f[3]) {   // entry fractions, then muon shares
    const Phase &a = purePhase[z[0]], &b = purePhase[z[1]], &c = purePhase[z[2]];
    const double det = (b.g - a.g) * (c.s - a.s) - (c.g - a.g) * (b.s - a.s);
    f[1] = ((m.g - a.g) * (c.s - a.s) - (c.g - a.g) * (m.s - a.s)) / det;
    f[2] = ((b.g - a.g) * (m.s - a.s) - (m.g - a.g) * (b.s - a.s)) / det;
    f[0] = 1. - f[1] - f[2];
    std::vector<double> u(3);
    double su = 0.;
    for (int i = 0; i < 3; ++i) { u[i] = f[i] / PerMuonInPhasor(z[i]); su += u[i]; }
    return std::vector<double>{u[0] / su, u[1] / su, u[2] / su};
  };
  struct PhasorCheck { std::string sample; std::vector<int> zs; std::vector<double> share, rms, truth, fitRms; double off = 0.; int n = 0; };
  std::vector<PhasorCheck> phasorChecks;
  for (size_t si = 0; si < samples.size(); ++si) {
    const Sample& s = samples[si];
    std::vector<int> zs = s.zs.size() == 2 ? s.zs : std::vector<int>{26, 24, 28};
    bool have = true;
    for (int z : zs) have = have && purePhase.count(z) && std::find(s.zs.begin(), s.zs.end(), z) != s.zs.end();
    if (!have || (s.zs.size() != 2 && s.label != "SS316")) continue;
    PhasorCheck pc;
    pc.sample = s.label;
    pc.zs = zs;
    auto Shares = [&](const Phase& ph, double& off) {
      if (zs.size() == 2) {
        double f1 = 0.;
        const double a1 = BinaryShare(ph, zs[0], zs[1], f1, off);
        return std::vector<double>{a1, 1. - a1};
      }
      double f[3];
      const int z3[3] = {zs[0], zs[1], zs[2]};
      off = 0.;
      return TernaryShares(ph, z3, f);
    };
    pc.share = Shares(s.whole.phase, pc.off);
    // the truth among these components, and the per-measurement spread (1e6) of the phasor's and the fit's shares
    const std::vector<double> t = fits[si].truth;
    double sum = 0.;
    for (int z : zs) sum += t[std::find(s.zs.begin(), s.zs.end(), z) - s.zs.begin()];
    for (int z : zs) pc.truth.push_back(t[std::find(s.zs.begin(), s.zs.end(), z) - s.zs.begin()] / sum);
    std::vector<double> m1(zs.size(), 0.), m2(zs.size(), 0.), f1(zs.size(), 0.), f2(zs.size(), 0.);
    for (const Result& r : fits[si].sets[1]) {
      double off = 0.;
      const std::vector<double> sh = Shares(r.phase, off);
      double fsum = 0.;
      for (int z : zs) fsum += r.a[std::find(s.zs.begin(), s.zs.end(), z) - s.zs.begin()];
      for (size_t i = 0; i < zs.size(); ++i) {
        const double fa = r.a[std::find(s.zs.begin(), s.zs.end(), zs[i]) - s.zs.begin()] / fsum;
        m1[i] += sh[i]; m2[i] += sh[i] * sh[i]; f1[i] += fa; f2[i] += fa * fa;
      }
      ++pc.n;
    }
    for (size_t i = 0; i < zs.size(); ++i) {
      pc.rms.push_back(pc.n > 1 ? std::sqrt(std::max(0., m2[i] / pc.n - std::pow(m1[i] / pc.n, 2))) : 0.);
      pc.fitRms.push_back(pc.n > 1 ? std::sqrt(std::max(0., f2[i] / pc.n - std::pow(f1[i] / pc.n, 2))) : 0.);
    }
    phasorChecks.push_back(pc);
  }

  // ---------------------------------------------------------------------------------------------- the numbers
  std::ofstream txt(out + "/muonic_unmixing.txt"), csv(out + "/muonic_unmixing.csv");
  txt << "Unmixing muonic X-ray spectra with the library's pure plates: " << base << "\n";
  txt << "  prompt entries (< " << promptNs << " ns), sigma " << s1MeV << " keV at 1 MeV, from " << threshold << " keV; fit "
      << eMin << "-" << eMax << " keV in " << nb << " bins of about half a resolution; Poisson likelihood, shares free\n";
  txt << "  share a_i: the part of a sample's muons element i catches (truth: cap_z, less the setup's own captures); "
         "mass % from the shares by the capture law (share ~ atoms x Z)\n\n";
  txt << "pure plates: prompt entries per muon in the fit's range (blind fit: those with >= " << minPerMuon << ")\n ";
  for (const auto& [z, v] : perMuon) txt << " " << Sym(z) << " " << Form("%.3g", v) << (v < minPerMuon ? " (not seen)" : "");
  txt << "\n\nself-test: each plate's even jobs (5e7 muons) fitted with all " << visible.size()
      << " visible plates' odd jobs: own share, the largest other in sigma, others at >= 3 sigma\n";
  int selfClaims = 0;
  double selfWorst = 0.;
  for (const SelfTest& st : self) {
    txt << Form("  %-3s %s %.5f +- %.5f   largest other %+.2f sigma (%s)   claims %d\n", Sym(st.z).c_str(), st.ok ? "" : "FAILED",
                st.own, st.ownErr, st.worst, Sym(st.worstZ).c_str(), st.claims);
    selfClaims += st.claims;
    selfWorst = std::max(selfWorst, std::fabs(st.worst));
  }
  txt << Form("  in all: %d claims at >= 3 sigma among %zu x %zu others (chance alone: %.1f); largest |other| %.2f sigma\n",
              selfClaims, self.size(), visible.size() - 1, 0.00135 * self.size() * (visible.size() - 1), selfWorst);
  csv << "sample,measurement,index,muons,element,truth,share,error\n";
  for (size_t si = 0; si < samples.size(); ++si) {
    const Sample& s = samples[si];
    const SampleFits& sf = fits[si];
    if (!sf.known.ok) continue;
    txt << "\n" << s.label << Form(": %.0f muons; known fit chi2 / ndf = %.0f / %d = %.2f", s.whole.muons, sf.known.chi2, sf.known.ndf, sf.known.chi2 / sf.known.ndf) << "\n";
    txt << "  element   mass %: true   fitted            share: true      fitted                (fit - true) / error   fit / true\n";
    double sumA = 0., sumV = 0.;
    for (size_t i = 0; i < s.zs.size(); ++i) {
      const int z = s.zs[i];
      sumA += sf.known.a[i];
      for (size_t k = 0; k < s.zs.size(); ++k) sumV += sf.known.cov[i][k];
      txt << Form("  %-4s  %10.4f  %9.4f +- %-8.4f  %10.6f  %10.6f +- %-9.6f  %+8.2f   %8.4f +- %.4f\n", Sym(z).c_str(), s.mass.at(z),
                  100. * sf.mass[i], 100. * sf.massErr[i], sf.truth[i], sf.known.a[i], sf.known.err[i],
                  (sf.known.a[i] - sf.truth[i]) / sf.known.err[i], sf.known.a[i] / sf.truth[i], sf.known.err[i] / sf.truth[i]);
      csv << s.label << ",whole,0," << s.whole.muons << "," << Sym(z) << "," << sf.truth[i] << "," << sf.known.a[i] << "," << sf.known.err[i] << "\n";
    }
    txt << Form("  sum of the shares %.5f +- %.5f (1 if the elements explain every muon)\n", sumA, std::sqrt(sumV));
    // the blind fit
    txt << Form("  blind fit (%zu elements free): chi2 / ndf = %.0f / %d;", visible.size(), sf.blind.chi2, sf.blind.ndf);
    std::set<int> held(s.zs.begin(), s.zs.end());
    int nClaims = 0;
    std::string claims;
    for (size_t q = 0; q < visible.size() && sf.blind.ok; ++q) {
      const int z = visible[q];
      if (held.count(z)) continue;
      const double sig = sf.blind.a[q] / sf.blind.err[q];
      if (std::fabs(sig) >= 3.) claims += Form(" %s %+.1f", Sym(z).c_str(), sig);
      if (sig >= 3.) ++nClaims;
    }
    txt << " absent elements at |share| >= 3 sigma:" << (claims.empty() ? " none" : claims) << "\n";
    txt << "    held:";
    for (size_t q = 0; q < visible.size() && sf.blind.ok; ++q)
      if (held.count(visible[q])) {
        const size_t i = std::find(s.zs.begin(), s.zs.end(), visible[q]) - s.zs.begin();
        txt << Form(" %s %.6f +- %.6f (known %.6f)", Sym(visible[q]).c_str(), sf.blind.a[q], sf.blind.err[q], sf.known.a[i]);
      }
    txt << "\n";
    txt << Form("  blind fits of the %d measurements of 1e6 muons: %d claim an absent element at >= 3 sigma\n", sf.blindSlices1e6, sf.blindClaims1e6);
    // the measurements: precision and pulls
    txt << "  measurements (known fit): mean error / true share, and the pulls (fit - true) / error: mean +- its error, rms\n";
    for (int set = 0; set < 3; ++set) {
      const std::vector<Result>& rs = sf.sets[set];
      if (rs.empty()) continue;
      txt << Form("    %zu x %s:", rs.size(), setName[set]);
      for (size_t i = 0; i < s.zs.size(); ++i) {
        double e = 0., t = 0., p1 = 0., p2 = 0.;
        for (const Result& r : rs) {
          e += r.err[i]; t += r.truth[i];
          const double p = (r.a[i] - r.truth[i]) / r.err[i];
          p1 += p; p2 += p * p;
        }
        const double n = rs.size(), mean = p1 / n, rms = std::sqrt(std::max(0., p2 / n - mean * mean));
        txt << Form("  %s %.3g%% [%+.2f+-%.2f, %.2f]", Sym(s.zs[i]).c_str(), 100. * e / t, mean, rms / std::sqrt(n), rms);
      }
      txt << "\n";
      int index = 0;
      for (const Result& r : rs) {
        for (size_t i = 0; i < s.zs.size(); ++i)
          csv << s.label << "," << setName[set] << "," << index << "," << r.muons << "," << Sym(s.zs[i]) << "," << r.truth[i] << "," << r.a[i] << "," << r.err[i] << "\n";
        ++index;
      }
    }
  }
  txt << "\nmatrix effect: X-ray entries per own muon in a peak's +-1.5 sigma (a line, or the sum of two X-rays in one pixel), "
         "inside the sample over the pure plate\n";
  for (const LineRatio& lr : ratios) txt << Form("  %-9s %-3s %8.2f keV  %.4f +- %.4f\n", lr.sample.c_str(), Sym(lr.z).c_str(), lr.e, lr.ratio, lr.err);
  txt << "\nphasor (harmonic " << harmonic << ", " << eMin << "-" << phasorMax << " keV): G, S of the pure plates\n ";
  for (const auto& [z, ph] : purePhase) txt << Form(" %s (%.4f, %.4f)", Sym(z).c_str(), ph.g, ph.s);
  txt << "\n";
  for (const Sample& s : samples) txt << Form("  %s (%.5f, %.5f)\n", s.label.c_str(), s.whole.phase.g, s.whole.phase.s);
  txt << "  shares read from the phasor (entry fractions over each plate's entries per muon), among the components named; "
         "spread over the 1e6 measurements, the full fit's beside it\n";
  for (const PhasorCheck& pc : phasorChecks) {
    txt << "  " << pc.sample << (pc.zs.size() == 2 ? Form(" (off the line by %.2g)", pc.off) : " (Fe, Cr, Ni only)") << ":";
    for (size_t i = 0; i < pc.zs.size(); ++i)
      txt << Form("  %s phasor %.5f (true %.5f; spread %.5f, fit's %.5f)", Sym(pc.zs[i]).c_str(), pc.share[i], pc.truth[i], pc.rms[i], pc.fitRms[i]);
    txt << "\n";
  }
  txt.close();
  csv.close();

  // ---------------------------------------------------------------------------------------------- the plots
  auto Title = [](const char* text) {
    TLatex t;
    t.SetNDC(); t.SetTextSize(0.036); t.SetTextAlign(22); t.SetTextFont(42);
    t.DrawLatex(0.53, 0.955, text);
  };
  auto Save = [&](TCanvas* c, const std::string& name) {
    c->Print((out + "/" + name + ".png").c_str());
    c->Print((out + "/" + name + ".pdf").c_str());
  };
  const int sampleMarker[3] = {20, 21, 33}, sampleColor[3] = {kRed + 1, kAzure + 2, kGreen + 3};
  const int palette[10] = {kRed + 1, kAzure + 2, kGreen + 3, kOrange + 7, kMagenta + 2, kCyan + 2, kViolet + 1, kGray + 2, kYellow + 3, kPink + 6};
  const int compColor[10] = {kAzure + 2, kGreen + 3, kOrange + 7, kMagenta + 2, kCyan + 2, kViolet + 1, kYellow + 3, kPink + 6, kSpring + 4, kGray + 2};
  auto ElementColor = [&](int z) {   // one colour per element in every plot
    static const std::map<int, int> c = {{26, kRed + 1}, {24, kAzure + 2}, {28, kGreen + 3}, {25, kOrange + 7}, {42, kMagenta + 2},
                                         {14, kCyan + 2}, {7, kViolet + 1}, {6, kGray + 2}, {15, kYellow + 3}, {16, kPink + 6}};
    return c.count(z) ? c.at(z) : kBlack;
  };

  // the spectra and their fits
  for (size_t si = 0; si < samples.size(); ++si) {
    const Sample& s = samples[si];
    const SampleFits& sf = fits[si];
    if (!sf.known.ok) continue;
    Matrix T, VT;
    Templates(s.zs, -1, T, VT);
    auto Hist = [&](const std::string& name, const std::vector<double>& m) {
      auto h = new TH1D(name.c_str(), ";energy (keV);entries per keV", nb, edges.data());
      h->SetDirectory(nullptr);
      for (int b = 0; b < nb; ++b) h->SetBinContent(b + 1, m[b] / (edges[b + 1] - edges[b]));
      return h;
    };
    auto c = new TCanvas(Form("c_spec_%zu", si), "", 1400, 850);
    c->SetLogy();
    c->SetLeftMargin(0.09); c->SetRightMargin(0.03); c->SetTopMargin(0.09); c->SetBottomMargin(0.1);
    TH1D* data = Hist("data_" + s.label, s.whole.d);
    data->GetXaxis()->SetRangeUser(eMin, phasorMax);
    double top = 0.;
    for (int b = 0; b < nb; ++b) if (edges[b] < phasorMax) top = std::max(top, data->GetBinContent(b + 1));
    data->SetMaximum(3. * top);
    data->SetMinimum(1e-4 * top);
    data->SetLineColor(kBlack);
    data->SetLineWidth(1);
    data->Draw("hist");
    std::vector<double> sum(nb, 0.);
    auto leg = new TLegend(0.6, 0.6, 0.96, 0.89);
    leg->SetBorderSize(0); leg->SetFillStyle(1001); leg->SetFillColor(kWhite); leg->SetTextSize(0.03); leg->SetTextFont(42);
    leg->AddEntry(data, Form("%s, %.3g muons (simulated)", s.label.c_str(), s.whole.muons), "l");
    const size_t shown = std::min<size_t>(s.zs.size(), 5);
    std::vector<double> rest(nb, 0.);
    for (size_t i = 0; i < s.zs.size(); ++i) {
      std::vector<double> m(nb);
      for (int b = 0; b < nb; ++b) { m[b] = s.whole.muons * sf.known.a[i] * T[b][i]; sum[b] += m[b]; }
      if (i < shown) {
        TH1D* h = Hist(Form("comp_%s_%d", s.label.c_str(), s.zs[i]), m);
        h->SetLineColor(compColor[i]);
        h->SetLineWidth(1);
        h->Draw("hist same");
        leg->AddEntry(h, Form("%s: share %.4g%%", Sym(s.zs[i]).c_str(), 100. * sf.known.a[i]), "l");
      } else
        for (int b = 0; b < nb; ++b) rest[b] += m[b];
    }
    if (shown < s.zs.size()) {
      TH1D* h = Hist("rest_" + s.label, rest);
      h->SetLineColor(kGray + 1);
      h->Draw("hist same");
      std::string names;
      for (size_t i = shown; i < s.zs.size(); ++i) names += (names.empty() ? "" : ", ") + Sym(s.zs[i]);
      leg->AddEntry(h, ("the others: " + names).c_str(), "l");
    }
    TH1D* fitted = Hist("fit_" + s.label, sum);
    fitted->SetLineColor(kRed);
    fitted->SetLineStyle(2);
    fitted->SetLineWidth(2);
    fitted->Draw("hist same");
    data->Draw("hist same");   // on top: where the fit agrees, the black hides the red dashes
    leg->AddEntry(fitted, Form("the fitted sum (#chi^{2}/ndf %.2f)", sf.known.chi2 / sf.known.ndf), "l");
    leg->Draw();
    Title(Form("%s: its prompt spectrum as a sum of the pure plates' spectra (fit %g-%g keV)", s.label.c_str(), eMin, eMax));
    Save(c, "unmixing_spectrum_" + s.label);
  }

  // fitted over true shares
  {
    auto c = new TCanvas("c_fvt", "", 1300, 900);
    c->SetLogx(); c->SetGridy();
    c->SetLeftMargin(0.1); c->SetRightMargin(0.03); c->SetTopMargin(0.09); c->SetBottomMargin(0.11);
    double lo = 1., hi = 1., xlo = 1.;
    std::vector<TGraphErrors*> gs;
    for (size_t si = 0; si < samples.size(); ++si) {
      if (!fits[si].known.ok) continue;
      auto g = new TGraphErrors();
      for (size_t i = 0; i < samples[si].zs.size(); ++i) {
        const double t = fits[si].truth[i], r = fits[si].known.a[i] / t, e = fits[si].known.err[i] / t;
        g->AddPoint(t, r);
        g->SetPointError(g->GetN() - 1, 0., e);
        lo = std::min(lo, r - e); hi = std::max(hi, r + e); xlo = std::min(xlo, t);
      }
      g->SetMarkerStyle(sampleMarker[si]); g->SetMarkerSize(1.6); g->SetMarkerColor(sampleColor[si]); g->SetLineColor(sampleColor[si]);
      gs.push_back(g);
    }
    const double pad = 0.1 * (hi - lo);
    auto fr = new TH1D("frame_fvt", ";true share of the muons;fitted share / true share", 100, 0.3 * xlo, 2.);
    fr->SetDirectory(nullptr);
    fr->SetMinimum(lo - pad); fr->SetMaximum(hi + 2. * pad);
    fr->GetXaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleSize(0.045);
    fr->Draw("axis");
    TLine one;
    one.SetLineStyle(2); one.SetLineColor(kGray + 2);
    one.DrawLine(0.3 * xlo, 1., 2., 1.);
    auto leg = new TLegend(0.6, 0.15, 0.96, 0.29);
    leg->SetBorderSize(0); leg->SetFillStyle(1001); leg->SetFillColor(kWhite); leg->SetTextSize(0.032); leg->SetTextFont(42);
    TLatex lab;   // the labels: left of the points (Ni80Cr20), above (Fe50Ni50), right (SS316)
    lab.SetTextSize(0.026); lab.SetTextFont(42); lab.SetTextAlign(12);
    size_t gi = 0;
    for (size_t si = 0; si < samples.size(); ++si) {
      if (!fits[si].known.ok) continue;
      gs[gi]->Draw("p same");
      leg->AddEntry(gs[gi], Form("%s (%.2g muons)", samples[si].label.c_str(), samples[si].whole.muons), "p");
      for (size_t i = 0; i < samples[si].zs.size(); ++i) {
        const double x = fits[si].truth[i], y = fits[si].known.a[i] / fits[si].truth[i];
        lab.SetTextColor(sampleColor[si]);
        lab.SetTextAlign(si % 3 == 0 ? 32 : si % 3 == 1 ? 21 : 12);
        lab.DrawLatex(si % 3 == 0 ? x / 1.15 : si % 3 == 1 ? x : x * 1.15, si % 3 == 1 ? y + 0.035 * (hi - lo + 3. * pad) : y, Sym(samples[si].zs[i]).c_str());
      }
      ++gi;
    }
    leg->Draw();
    Title(Form("Each element's fitted share over its true share (whole samples, fit %g-%g keV)", eMin, eMax));
    Save(c, "unmixing_fit_vs_truth");
  }

  // the pulls of the measurements
  {
    auto c = new TCanvas("c_pulls", "", 1300, 900);
    c->SetLeftMargin(0.1); c->SetRightMargin(0.03); c->SetTopMargin(0.09); c->SetBottomMargin(0.11);
    TH1D* h[2];
    double top = 0.;
    for (int set = 0; set < 2; ++set) {
      h[set] = new TH1D(Form("pulls_%d", set), ";(fitted share - the mean of its measurements) / its error;measurements x elements", 48, -6., 6.);
      h[set]->SetDirectory(nullptr);
      for (size_t si = 0; si < samples.size(); ++si) {
        const std::vector<Result>& rs = fits[si].sets[set];
        if (rs.size() < 2) continue;
        for (size_t i = 0; i < samples[si].zs.size(); ++i) {
          double mean = 0.;
          for (const Result& r : rs) mean += r.a[i] / rs.size();
          for (const Result& r : rs) h[set]->Fill((r.a[i] - mean) / r.err[i] / std::sqrt(1. - 1. / rs.size()));
        }
      }
      top = std::max({top, h[set]->GetMaximum(), h[set]->GetEntries() * h[set]->GetBinWidth(1) / std::sqrt(2. * M_PI)});
    }
    h[0]->SetLineColor(kAzure + 2); h[0]->SetLineWidth(2);
    h[1]->SetLineColor(kRed + 1); h[1]->SetLineWidth(2);
    h[0]->SetMaximum(1.3 * top);
    h[0]->GetXaxis()->SetTitleSize(0.045); h[0]->GetYaxis()->SetTitleSize(0.045);
    h[0]->Draw("hist");
    h[1]->Draw("hist same");
    auto leg = new TLegend(0.13, 0.7, 0.6, 0.89);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.03); leg->SetTextFont(42);
    for (int set = 0; set < 2; ++set) {
      auto g = new TGraph();
      for (double x = -6.; x <= 6.001; x += 0.05) g->AddPoint(x, h[set]->GetEntries() * h[set]->GetBinWidth(1) * std::exp(-0.5 * x * x) / std::sqrt(2. * M_PI));
      g->SetLineColor(h[set]->GetLineColor()); g->SetLineStyle(2);
      g->Draw("l same");
      leg->AddEntry(h[set], Form("%.0f shares from measurements of %s muons: rms %.2f", h[set]->GetEntries(), setName[set], h[set]->GetRMS()), "l");
    }
    auto gl = new TGraph();
    gl->SetLineStyle(2); gl->SetLineColor(kGray + 2);
    leg->AddEntry(gl, "a unit Gaussian: the error bars match the scatter", "l");
    leg->Draw();
    Title("Repeated measurements: each share's scatter around its mean, in units of its error");
    Save(c, "unmixing_pulls");
  }

  // the precision against the number of muons (SS316)
  for (size_t si = 0; si < samples.size(); ++si) {
    const Sample& s = samples[si];
    if (s.label != "SS316" || !fits[si].known.ok) continue;
    auto c = new TCanvas("c_prec", "", 1300, 900);
    c->SetLogx(); c->SetLogy(); c->SetGridy();
    c->SetLeftMargin(0.11); c->SetRightMargin(0.16); c->SetTopMargin(0.09); c->SetBottomMargin(0.11);
    double lo = 1e9, hi = 0.;
    std::vector<TGraph*> gs;
    for (size_t i = 0; i < s.zs.size(); ++i) {
      auto g = new TGraph();
      for (int set = 0; set < 3; ++set) {
        const std::vector<Result>& rs = fits[si].sets[set];
        if (rs.empty()) continue;
        double e = 0., t = 0., n = 0.;
        for (const Result& r : rs) { e += r.err[i]; t += r.truth[i]; n += r.muons; }
        g->AddPoint(n / rs.size(), e / t);
      }
      g->AddPoint(s.whole.muons, fits[si].known.err[i] / fits[si].truth[i]);
      for (int p = 0; p < g->GetN(); ++p) { lo = std::min(lo, g->GetPointY(p)); hi = std::max(hi, g->GetPointY(p)); }
      g->SetLineColor(palette[i]); g->SetMarkerColor(palette[i]); g->SetMarkerStyle(20); g->SetMarkerSize(1.3); g->SetLineWidth(2);
      gs.push_back(g);
    }
    auto fr = new TH1D("frame_prec", ";muons in the measurement;error of the share / the share", 100, 5e4, 2e8);
    fr->SetDirectory(nullptr);
    fr->SetMinimum(0.5 * lo); fr->SetMaximum(2. * hi);
    fr->GetXaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleSize(0.045);
    fr->Draw("axis");
    auto leg = new TLegend(0.85, 0.15, 0.99, 0.89);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032); leg->SetTextFont(42);
    for (size_t i = 0; i < s.zs.size(); ++i) {
      gs[i]->Draw("lp same");
      leg->AddEntry(gs[i], Form("%s %.3g%%", Sym(s.zs[i]).c_str(), s.mass.at(s.zs[i])), "lp");
    }
    leg->Draw();
    TLatex t;
    t.SetNDC(); t.SetTextSize(0.03); t.SetTextFont(42); t.SetTextColor(kGray + 2);
    t.DrawLatex(0.14, 0.15, "the last point: the whole sample, with the templates' own noise (1e8 muons each)");
    Title("SS316: how the error of each element's share falls with the muons measured");
    Save(c, "unmixing_precision");
  }

  // the matrix effect
  {
    auto c = new TCanvas("c_matrix", "", 1300, 900);
    c->SetLogx(); c->SetGridy();
    c->SetLeftMargin(0.1); c->SetRightMargin(0.17); c->SetTopMargin(0.09); c->SetBottomMargin(0.11);
    double lo = 1., hi = 1.;
    for (const LineRatio& lr : ratios) { lo = std::min(lo, lr.ratio - lr.err); hi = std::max(hi, lr.ratio + lr.err); }
    const double pad = 0.1 * (hi - lo);
    auto fr = new TH1D("frame_matrix", ";line energy (keV);entries per muon: inside the sample / pure plate", 100, eMin, phasorMax);
    fr->SetDirectory(nullptr);
    fr->SetMinimum(lo - pad); fr->SetMaximum(hi + pad);
    fr->GetXaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleSize(0.045);
    fr->GetXaxis()->SetMoreLogLabels(); fr->GetXaxis()->SetNoExponent();
    fr->Draw("axis");
    TLine one;
    one.SetLineStyle(2); one.SetLineColor(kGray + 2);
    one.DrawLine(eMin, 1., phasorMax, 1.);
    auto leg = new TLegend(0.84, 0.12, 0.99, 0.9);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.026); leg->SetTextFont(42);
    std::map<std::pair<std::string, int>, TGraphErrors*> graphs;
    std::vector<std::pair<std::string, int>> order;
    for (const LineRatio& lr : ratios) {
      auto key = std::make_pair(lr.sample, lr.z);
      if (!graphs.count(key)) { graphs[key] = new TGraphErrors(); order.push_back(key); }
      auto g = graphs[key];
      g->AddPoint(lr.e, lr.ratio);
      g->SetPointError(g->GetN() - 1, 0., lr.err);
    }
    for (const auto& key : order) {
      size_t si = 0;
      while (si < samples.size() && samples[si].label != key.first) ++si;
      auto g = graphs[key];
      g->SetMarkerStyle(sampleMarker[si % 3]); g->SetMarkerSize(sampleMarker[si % 3] == 33 ? 2.0 : 1.4);
      g->SetMarkerColor(ElementColor(key.second)); g->SetLineColor(ElementColor(key.second));
      g->Draw("p same");
      leg->AddEntry(g, Form("%s in %s", Sym(key.second).c_str(), key.first.c_str()), "p");
    }
    leg->Draw();
    Title("Matrix effect: each element's peaks (lines, sums) inside the sample over its pure plate's");
    Save(c, "unmixing_matrix");
  }

  // the blind fit: the absent elements in units of their errors
  {
    auto c = new TCanvas("c_blind", "", 1400, 850);
    c->SetGridy();
    c->SetLeftMargin(0.08); c->SetRightMargin(0.03); c->SetTopMargin(0.09); c->SetBottomMargin(0.1);
    const int nv = visible.size();
    double lo = -4., hi = 4.;
    std::vector<TGraph*> gs;
    for (size_t si = 0; si < samples.size(); ++si) {
      auto g = new TGraph();
      std::set<int> held(samples[si].zs.begin(), samples[si].zs.end());
      if (fits[si].blind.ok)
        for (int q = 0; q < nv; ++q) {
          if (held.count(visible[q])) continue;
          const double sig = fits[si].blind.a[q] / fits[si].blind.err[q];
          g->AddPoint(q + 0.5 + 0.2 * (int(si) - 1), sig);
          lo = std::min(lo, sig); hi = std::max(hi, sig);
        }
      g->SetMarkerStyle(sampleMarker[si]); g->SetMarkerSize(1.4); g->SetMarkerColor(sampleColor[si]);
      gs.push_back(g);
    }
    auto fr = new TH1D("frame_blind", ";;fitted share / its error", nv, 0., nv);
    fr->SetDirectory(nullptr);
    for (int q = 0; q < nv; ++q) fr->GetXaxis()->SetBinLabel(q + 1, Sym(visible[q]).c_str());
    fr->GetXaxis()->SetLabelSize(0.045);
    fr->GetYaxis()->SetTitleSize(0.04);
    fr->SetMinimum(lo - 1.); fr->SetMaximum(hi + 0.3 * (hi - lo));
    fr->Draw("axis");
    TLine l3;
    l3.SetLineStyle(2); l3.SetLineColor(kGray + 2);
    l3.DrawLine(0., 3., nv, 3.);
    l3.DrawLine(0., -3., nv, -3.);
    auto leg = new TLegend(0.1, 0.78, 0.55, 0.89);
    leg->SetBorderSize(0); leg->SetFillStyle(1001); leg->SetFillColor(kWhite); leg->SetTextSize(0.03); leg->SetTextFont(42);
    leg->SetNColumns(3);
    for (size_t si = 0; si < samples.size(); ++si) {
      gs[si]->Draw("p same");
      leg->AddEntry(gs[si], samples[si].label.c_str(), "p");
    }
    leg->Draw();
    Title(Form("Blind fit (%d elements free, %g-%g keV): the elements the samples do not hold", nv, eMin, eMax));
    Save(c, "unmixing_blind");
  }

  // the phasor plane
  auto PhasorPlot = [&](const char* name, bool zoom) {
    auto c = new TCanvas(name, "", 1100, 1000);
    c->SetLeftMargin(0.12); c->SetRightMargin(0.04); c->SetTopMargin(0.09); c->SetBottomMargin(0.1);
    double x0 = 1e9, x1 = -1e9, y0 = 1e9, y1 = -1e9;
    auto Extend = [&](double x, double y) { x0 = std::min(x0, x); x1 = std::max(x1, x); y0 = std::min(y0, y); y1 = std::max(y1, y); };
    std::set<int> near;
    for (const Sample& s : samples) for (int z : s.zs) if (!zoom || s.zs.size() <= 3 || z == 24 || z == 26 || z == 28) near.insert(z);
    for (const auto& [z, ph] : purePhase) if (!zoom || near.count(z)) Extend(ph.g, ph.s);
    for (const Sample& s : samples) Extend(s.whole.phase.g, s.whole.phase.s);
    const double span = std::max(x1 - x0, y1 - y0), mx = 0.5 * (x0 + x1), my = 0.5 * (y0 + y1);
    const double half = 0.5 * span * (zoom ? 1.35 : 1.12);
    auto fr = new TH2D(Form("frame_%s", name), ";G (cosine part);S (sine part)", 10, mx - half, mx + half, 10, my - half, my + half);
    fr->SetDirectory(nullptr);
    fr->GetXaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleSize(0.045);
    fr->Draw("axis");
    if (!zoom) {   // the pure plates in Z order
      auto line = new TGraph();
      for (const auto& [z, ph] : purePhase) line->AddPoint(ph.g, ph.s);
      line->SetLineColor(kGray); line->SetLineWidth(1);
      line->Draw("l same");
    }
    TLatex lab;
    lab.SetTextFont(42); lab.SetTextSize(zoom ? 0.034 : 0.024); lab.SetTextAlign(11); lab.SetTextColor(kGray + 3);
    auto pts = new TGraph();
    for (const auto& [z, ph] : purePhase) {
      if (ph.g < mx - half || ph.g > mx + half || ph.s < my - half || ph.s > my + half) continue;
      pts->AddPoint(ph.g, ph.s);
      lab.DrawLatex(ph.g + 0.01 * half, ph.s + 0.01 * half, Sym(z).c_str());
    }
    pts->SetMarkerStyle(20); pts->SetMarkerSize(zoom ? 1.3 : 0.9); pts->SetMarkerColor(kGray + 3);
    pts->Draw("p same");
    auto leg = new TLegend(0.15, 0.13, 0.62, 0.3);
    leg->SetBorderSize(0); leg->SetFillStyle(1001); leg->SetFillColor(kWhite); leg->SetTextSize(0.03); leg->SetTextFont(42);
    leg->AddEntry(pts, "pure plates", "p");
    std::vector<size_t> drawOrder;   // the triangle first, under the two-plate lines
    for (size_t si = 0; si < samples.size(); ++si) if (samples[si].zs.size() != 2) drawOrder.push_back(si);
    for (size_t si = 0; si < samples.size(); ++si) if (samples[si].zs.size() == 2) drawOrder.push_back(si);
    for (size_t si : drawOrder) {
      const Sample& s = samples[si];
      std::vector<int> comps = s.zs.size() == 2 ? s.zs : std::vector<int>{26, 24, 28};
      auto poly = new TGraph();
      for (int z : comps) if (purePhase.count(z)) poly->AddPoint(purePhase[z].g, purePhase[z].s);
      if (comps.size() == 3 && poly->GetN() == 3) poly->AddPoint(poly->GetPointX(0), poly->GetPointY(0));
      poly->SetLineColor(sampleColor[si]); poly->SetLineStyle(comps.size() == 3 ? 3 : 2); poly->SetLineWidth(2);
      poly->Draw("l same");
      auto star = new TGraph();
      star->AddPoint(s.whole.phase.g, s.whole.phase.s);
      star->SetMarkerStyle(29); star->SetMarkerSize(zoom ? 3.2 : 2.4); star->SetMarkerColor(sampleColor[si]);
      star->Draw("p same");
      leg->AddEntry(star, Form("%s (%s)", s.label.c_str(), comps.size() == 2 ? "dashed: its two plates" : "dotted: the Fe, Cr, Ni triangle"), "p");
    }
    leg->Draw();
    Title(Form("Phasor picture %s (harmonic %d, %g-%g keV)", zoom ? "around the samples" : "of the prompt spectra", harmonic, eMin, phasorMax));
    Save(c, name);
  };
  PhasorPlot("unmixing_phasor", false);
  PhasorPlot("unmixing_phasor_zoom", true);

  printf("written into %s: muonic_unmixing.txt, .csv and the unmixing_* plots (.png, .pdf)\n", out.c_str());
}
