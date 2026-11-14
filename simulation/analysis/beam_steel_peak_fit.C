// The peaks of the steel study's tagged-beam spectrum, found blind, fitted and named (2026-11-10, the user: "For the
// simulated spectrum, try to identify all the peaks with fit and their estimated information, show it in a more clear
// plot. If necessary separate panels to show detailed peaks").
//
// Input: <dir>/analysis/beam_steel_trace_spectrum.root of analysis/beam_steel_trace_spectrum.C: the tagged-beam view
// (the hits within 10 ns of their muon), one day at the run's rate, in 0.5 keV bins with the simulation's own errors; its
// parts (each element's muonic X-rays, the beam's other entries, the accidentals); the emitted lines (truth).
// Bins (2026-11-11, the user: finer at low energy): the fits use about half a resolution, 0.5 keV up to 100 keV, 1 keV up
// to 400 keV, 2 keV above (BinWidth); each plot one width, that of its lowest energy.
//  1. Search, blind: at every bin a Gaussian of the CZT's resolution (sigma1MeV x sqrt(E / MeV)) on a straight line over
//     +-4 sigma, a weighted linear fit; each local maximum of its significance (area / error, in the simulation's
//     statistics) above zMin is a candidate.
//  2. Fit, in sections of 30 resolutions (fit windows 6 resolutions wider on each side; a peak is reported by the section
//     that holds it): the candidates as Gaussians (integrated over each bin) with widths at the resolution, free areas and
//     positions, on a cubic background (Legendre), chi2 with the simulation's errors. The residuals are searched the same
//     way and a peak above zMin added (up to 8 times per section); peaks under 3 sigma are dropped. A fit with free widths
//     is printed beside: sigma_free / sigma_res is 1 for one line, more for a blend.
//  3. Uncertainties for one day of data: the simulation's errors are sqrt(w) times a day's Poisson errors (w the mean
//     weight of an entry in the window, about 114 = one day of muons / the simulated ones): the fit's errors / sqrt(w).
//  4. Names: the same Gaussians (positions and widths kept) fitted to each part of the spectrum give each part's share of
//     a peak (the shares add up to the peak). The largest names it. An element's X-rays: the strongest emitted line within
//     2 sigma, named by its transition (K, L, M, ... by the lower level; alpha, beta, ... by Delta n) from the cascade's
//     levels (measured K_alpha and L_alpha, M_alpha of Th and U; Dirac x 1.0047 from n = 3, as TOFmeterEmCaptureCascade),
//     else a Cd or Te K X-ray escape (23.0-31.0 keV below a line) or two lines summed in one pixel. The beam's other
//     entries: a capture gamma (truth) or a known line (annihilation, (n,n')). The accidentals: an activation line.
//  5. Lines that should be seen (2026-11-11): each element's emitted lines (truth) x a full-energy efficiency (Eff x a factor
//     from the fitted K_alpha, L_alpha, M_alpha, N_alpha lines), kept where the search would find them alone (>= zMin); at the
//     end of peaks_prompt.txt with the peak each lies under, grey in the linear plots.
// Outputs in outDir (default <dir>/analysis): peaks_prompt.txt and .csv (every peak), and the plots (.png, .pdf)
//   peaks_prompt_overview   the whole spectrum, every fitted peak numbered (log energy axis)
//   peaks_prompt_panels     six energy regions, each peak named with its fitted energy
//   peaks_prompt_linear_<lo>_<hi>   (2026-11-11) one plot per energy range, linear: the simulated spectrum minus the baseline
//                           (the fit's background and structures), the lines' fit over it, each peak named; grey dashed
//                           the lines that should be seen (5.), named in grey where no red label has their name
// Run in the work folder after beam_steel_trace_spectrum.C:
//   root -l -b -q 'analysis/beam_steel_peak_fit.C("output/czt/beam_steel_czt_trace")'
#include <Math/MinimizerOptions.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMatrixDSym.h>
#include <TNamed.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVectorD.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace peakfit {
double gS1 = 6.36991;   // the resolution at 1 MeV (keV), read from the ROOT file's settings
double Res(double e) { return gS1 * std::sqrt(std::max(e, 1.) / 1000.); }
double GaussIn(double lo, double hi, double mu, double s)   // the share of a Gaussian (mu, s) in [lo, hi]
{
  return 0.5 * (std::erf((hi - mu) / (M_SQRT2 * s)) - std::erf((lo - mu) / (M_SQRT2 * s)));
}
// the fit's bins: 0.5 keV up to 100 keV, 1 keV up to 400 keV, 2 keV above (half a resolution or less; the input's bins
// merged), never finer than the input's bins (gFineKeV)
double gFineKeV = 0.5;
const double kRefKeV = 2.;   // the backgrounds' unit: entries per 2 keV (the former bins: the same numbers in 2 keV bins)
double BinWidth(double e)
{
  const double w = e < 100. ? 0.5 : e < 400. ? 1. : 2.;
  return std::max(gFineKeV, std::ceil(w / gFineKeV - 1e-6) * gFineKeV);
}
double Legendre(int k, double u)
{
  return k == 0 ? 1. : k == 1 ? u : k == 2 ? 0.5 * (3. * u * u - 1.) : 0.5 * (5. * u * u * u - 3. * u);
}

// the measured lines of TOFmeterEmCaptureCascade's table (keV): K_alpha 2p3/2, 2p1/2 -> 1s; L_alpha 3d5/2 -> 2p3/2,
// 3d3/2 -> 2p1/2, 3d3/2 -> 2p3/2; M_alpha (Th, U) 4f7/2 -> 3d5/2, 4f5/2 -> 3d3/2, 4f5/2 -> 3d5/2
struct ElementInfo {
  int z;
  const char* name;
  double a;
  double ka[2], la[3], ma[3];
};
const std::vector<ElementInfo> kElements = {
  {26, "Fe", 55.845, {1257.19, 1253.06}, {265.30, 268.90, 264.87}, {0., 0., 0.}},
  {24, "Cr", 51.996, {1094.46, 1091.45}, {226.23, 228.93, 225.92}, {0., 0., 0.}},
  {28, "Ni", 58.693, {1431.83, 1426.31}, {308.23, 313.04, 307.66}, {0., 0., 0.}},
  {42, "Mo", 95.95, {2720.05, 2694.97}, {694.36, 718.37, 691.44}, {0., 0., 0.}},
  {25, "Mn", 54.938, {1175.60, 1172.06}, {245.57, 248.74, 245.21}, {0., 0., 0.}},
  {14, "Si", 28.086, {400.29, 399.94}, {76.61, 76.93, 76.57}, {0., 0., 0.}},
  {7, "N", 14.007, {102.41, 102.39}, {19.03, 19.05, 19.03}, {0., 0., 0.}},
  {6, "C", 12.011, {75.26, 75.25}, {13.96, 13.97, 13.96}, {0., 0., 0.}},
  {15, "P", 30.974, {457.49, 457.03}, {88.00, 88.43, 87.95}, {0., 0., 0.}},
  {16, "S", 32.06, {517.63, 517.04}, {100.17, 100.72, 100.11}, {0., 0., 0.}},
  {19, "K", 39.098, {713.07, 711.92}, {141.48, 142.57, 141.36}, {0., 0., 0.}},
  {90, "Th", 232.038, {6367.65, 6103.60}, {2949.06, 3148.50, 2883.54}, {1133.82, 1186.38, 1120.51}},
  {92, "U", 238.029, {6420.00, 6133.00}, {3060.37, 3275.56, 2988.56}, {1185.82, 1243.45, 1171.27}}};

double WeightedMean(const double* v, const double* w, int n)
{
  double s = 0., sw = 0.;
  for (int i = 0; i < n; ++i) { s += v[i] * w[i]; sw += w[i]; }
  return s / sw;
}
// the binding energy of the circular level n (keV), as the cascade builds it
double Level(const ElementInfo& el, int n)
{
  static const double wk[2] = {2., 1.}, wl[3] = {9., 5., 1.};
  if (n == 2) return Level(el, 3) + WeightedMean(el.la, wl, 3);
  if (n == 1) return Level(el, 2) + WeightedMean(el.ka, wk, 2);
  const double mmu = 105658.3755, nucleus = el.a * 931494.10242, m = mmu * nucleus / (mmu + nucleus), za = el.z / 137.035999084;
  double s = 0., sw = 0.;
  for (double j : {n - 0.5, n - 1.5}) {
    const double k = j + 0.5, d = k - std::sqrt(k * k - za * za), w = 2. * j + 1.;
    s += w * m * (1. - 1. / std::sqrt(1. + std::pow(za / (n - d), 2)));
    sw += w;
  }
  return 1.0047 * s / sw;
}
// the fine-structure splitting of the level n, l (keV; Dirac, point nucleus): j = l - 1/2 below j = l + 1/2; 0 for s
double Splitting(const ElementInfo& el, int n, int l)
{
  if (l < 1) return 0.;
  const double mmu = 105658.3755, nucleus = el.a * 931494.10242, m = mmu * nucleus / (mmu + nucleus), za = el.z / 137.035999084;
  auto Binding = [&](double k) { const double d = k - std::sqrt(k * k - za * za); return m * (1. - 1. / std::sqrt(1. + std::pow(za / (n - d), 2))); };
  return Binding(l) - Binding(l + 1.);   // k = j + 1/2
}
// the distance of a transition's two main components (keV): the lower level's splitting less the upper one's (the lower
// level circular, l = lo - 1; the upper one l = lo, a p state for the K series)
double DoubletGap(const ElementInfo& el, int up, int lo)
{
  return std::fabs(Splitting(el, lo, lo - 1) - Splitting(el, up, lo == 1 ? 1 : lo));
}
// the transition nUp -> nLow whose photon is nearest to e: the measured components for 2 -> 1, 3 -> 2 (and 4 -> 3 of Th
// and U), the levels' difference otherwise
void Transition(const ElementInfo& el, double e, int& nUp, int& nLow)
{
  double best = 1e300;
  for (int lo = 1; lo <= 14; ++lo)
    for (int up = lo + 1; up <= lo + 7; ++up) {
      std::vector<double> es;
      if (up == 2 && lo == 1) es.assign(el.ka, el.ka + 2);
      else if (up == 3 && lo == 2) es.assign(el.la, el.la + 3);
      else if (up == 4 && lo == 3 && el.ma[0] > 0.) es.assign(el.ma, el.ma + 3);
      else es.push_back(Level(el, lo) - Level(el, up));
      for (double v : es)
        if (std::fabs(v - e) < best) { best = std::fabs(v - e); nUp = up; nLow = lo; }
    }
}
std::string SeriesName(int nUp, int nLow, bool latex)
{
  static const char* series = "KLMNOPQRSTUVWX";
  static const char* greekLatex[] = {"#alpha", "#beta", "#gamma", "#delta", "#varepsilon", "#zeta", "#eta"};
  static const char* greekText[] = {"alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta"};
  const int dn = std::min(7, nUp - nLow);
  return latex ? Form("%c_{%s}", series[nLow - 1], greekLatex[dn - 1]) : Form("%c_%s (%d->%d)", series[nLow - 1], greekText[dn - 1], nUp, nLow);
}
// a line's name; the measured fine-structure components numbered where the CZT can tell them apart (> 0.75 sigma):
// K_alpha1 2p3/2 -> 1s, K_alpha2 2p1/2 -> 1s; L_alpha1 3d5/2 -> 2p3/2, L_alpha2 3d3/2 -> 2p1/2, L_alpha3 3d3/2 -> 2p3/2
std::string LineName(const ElementInfo& el, double e, bool latex)
{
  int up = 0, lo = 0;
  Transition(el, e, up, lo);
  const double* comp = nullptr;
  int n = 0;
  if (up == 2 && lo == 1) { comp = el.ka; n = 2; }
  else if (up == 3 && lo == 2) { comp = el.la; n = 3; }
  else if (up == 4 && lo == 3 && el.ma[0] > 0.) { comp = el.ma; n = 3; }
  if (!comp || std::fabs(comp[0] - comp[1]) < 0.75 * Res(e)) return SeriesName(up, lo, latex);
  int best = 0;
  for (int i = 1; i < n; ++i) if (std::fabs(comp[i] - e) < std::fabs(comp[best] - e)) best = i;
  static const char* states[3][3] = {{"2p3/2->1s", "2p1/2->1s", ""}, {"3d5/2->2p3/2", "3d3/2->2p1/2", "3d3/2->2p3/2"},
                                     {"4f7/2->3d5/2", "4f5/2->3d3/2", "4f5/2->3d5/2"}};
  return latex ? std::string(Form("%c_{#alpha%d}", "KLM"[lo - 1], best + 1))
               : std::string(Form("%c_alpha%d (%s)", "KLM"[lo - 1], best + 1, states[lo - 1][best]));
}

// a rough full-energy efficiency, only to rank the possible names (it enters no number): the steel in front of the panel
// (0.27 cm of iron: the muons stop 1.8 mm from that face, seen at an angle) x the CZT's absorption (10 mm), x (300 keV /
// E)^1.3 above 300 keV for the full-energy share. Attenuation (cm2/g, NIST XCOM, rounded).
double LogLog(const std::vector<std::pair<double, double>>& t, double e)
{
  if (e <= t.front().first) return t.front().second;
  for (size_t i = 1; i < t.size(); ++i)
    if (e <= t[i].first) {
      const double f = std::log(e / t[i - 1].first) / std::log(t[i].first / t[i - 1].first);
      return std::exp(std::log(t[i - 1].second) + f * std::log(t[i].second / t[i - 1].second));
    }
  return t.back().second;
}
const std::vector<std::pair<double, double>> kMuFe = {{10., 170.6}, {15., 57.08}, {20., 25.68}, {30., 8.176}, {40., 3.629},
  {50., 1.958}, {60., 1.205}, {80., 0.5952}, {100., 0.3717}, {150., 0.1964}, {200., 0.146}, {300., 0.1099}, {500., 0.084},
  {1000., 0.0599}, {2000., 0.0426}, {3000., 0.0362}, {5000., 0.0314}, {8000., 0.0299}};
const std::vector<std::pair<double, double>> kMuCdTe = {{10., 150.}, {20., 22.}, {26.6, 11.}, {26.8, 34.}, {31.7, 22.}, {31.9, 38.},
  {40., 20.}, {50., 11.}, {60., 6.5}, {80., 3.1}, {100., 1.85}, {150., 0.68}, {200., 0.35}, {300., 0.16}, {500., 0.092},
  {1000., 0.057}, {2000., 0.041}, {3000., 0.037}, {5000., 0.037}, {8000., 0.039}};
double Eff(double e)
{
  return std::exp(-7.87 * 0.27 * LogLog(kMuFe, e)) * (1. - std::exp(-5.85 * LogLog(kMuCdTe, e))) * std::min(1., std::pow(300. / e, 1.3));
}
// the chance that a Cd or Te K X-ray escapes the pixel: the absorption within 0.25 mm of the entrance face, x 0.25
double EscapeProb(double e)
{
  const double m = 5.85 * LogLog(kMuCdTe, e);
  return 0.25 * (1. - std::exp(-0.025 * m)) / (1. - std::exp(-m));
}
struct KXray {   // the CZT's K X-rays (keV), their share of the escapes, the K edge they need
  double e, weight, edge;
  const char* latex;
  const char* text;
};
const std::vector<KXray> kCdTe = {{23.174, 0.45 * 0.56, 26.71, "Cd K_{#alpha}", "Cd K_alpha1"},
                                  {22.984, 0.45 * 0.29, 26.71, "Cd K_{#alpha}", "Cd K_alpha2"},
                                  {26.095, 0.45 * 0.15, 26.71, "Cd K_{#beta}", "Cd K_beta"},
                                  {27.472, 0.55 * 0.55, 31.81, "Te K_{#alpha}", "Te K_alpha1"},
                                  {27.202, 0.55 * 0.29, 31.81, "Te K_{#alpha}", "Te K_alpha2"},
                                  {30.996, 0.55 * 0.16, 31.81, "Te K_{#beta}", "Te K_beta"}};
// the kinds of names
enum Kind { kindNone = 0, kindLine, kindEscape, kindPairEscape, kindSum, kindFluorescence, kindGamma, kindStructure };
const char* KindName(int k)
{
  static const char* names[] = {"?", "line", "Cd/Te escape", "pair escape", "sum", "Cd/Te fluorescence", "gamma", "structure"};
  return names[k];
}

struct KnownLine {
  double e;
  const char* latex;
  const char* text;
};
const std::vector<KnownLine> kKnown = {
  {511.0, "e^{+}e^{-}", "e+e- annihilation"},         {846.8, "^{56}Fe(n,n')", "56Fe(n,n') 846.8 (and 56Mn decay)"},
  {1238.3, "^{56}Fe(n,n')", "56Fe(n,n') 1238.3"},      {1408.2, "^{54}Fe(n,n')", "54Fe(n,n') 1408.2"},
  {1434.1, "^{52}Cr(n,n')", "52Cr(n,n') 1434.1 (and 52V decay)"}, {1332.5, "^{60}Ni(n,n')", "60Ni(n,n') 1332.5"},
  {1454.3, "^{58}Ni(n,n')", "58Ni(n,n') 1454.3"},      {1810.7, "^{56}Mn", "56Mn decay 1810.7"},
  {2113.1, "^{56}Mn", "56Mn decay 2113.1"},            {558.5, "^{114}Cd", "113Cd(n,g)114Cd 558.5"},
  {1460.8, "^{40}K", "40K decay 1460.8"},              {609.3, "^{214}Bi", "214Bi decay 609.3"},
  {2614.5, "^{208}Tl", "208Tl decay 2614.5"}};

// weighted linear least squares: y ~ sum_j c_j cols[j]; the weights are scaled to their mean (the covariance back)
bool LinearFit(const std::vector<std::vector<double>>& cols, const std::vector<double>& y, const std::vector<double>& w,
               std::vector<double>& c, TMatrixDSym* cov = nullptr)
{
  const int m = cols.size(), n = y.size();
  double wMean = 0.;
  for (double v : w) wMean += v / n;
  TMatrixDSym nm(m);
  TVectorD r(m);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j) {
      const double wi = w[i] / wMean;
      r[j] += wi * cols[j][i] * y[i];
      for (int k = 0; k <= j; ++k) nm(j, k) += wi * cols[j][i] * cols[k][i];
    }
  for (int j = 0; j < m; ++j)
    for (int k = 0; k < j; ++k) nm(k, j) = nm(j, k);
  double det = 0.;
  nm.Invert(&det);
  if (det == 0. || !std::isfinite(det)) return false;
  c.assign(m, 0.);
  for (int j = 0; j < m; ++j)
    for (int k = 0; k < m; ++k) c[j] += nm(j, k) * r[k];
  if (cov) { cov->ResizeTo(m, m); *cov = nm; *cov *= 1. / wMean; }
  return true;
}

struct Candidate {   // a peak of the search (or of a residual search), updated by the fits
  double mu = 0., area = 0., z = 0.;
  bool alive = true, reported = false;
};
struct FitPeak {   // a peak in a section's fit; external peaks (outside the window) stay fixed
  double mu = 0., area = 0., sigma = 0.;
  double dmu = 0., darea = 0.;
  double sFree = 0., dsFree = 0., areaFree = 0., dareaFree = 0.;
  bool free = true;
  int candidate = -1;
  double mu0 = 0.;   // where it entered the section: its position may move +-2.5 resolutions from here (no drift)
};
struct Section {   // a fitted section: its window, its background, its peaks
  double c0 = 0., c1 = 0., w0 = 0., w1 = 0.;
  std::vector<double> bg;
  std::vector<FitPeak> peaks;
  double chi2 = 0.;
  int ndf = 0, status = 0;
  double Model(double x, double bw, bool withPeaks = true) const   // the entries of a bin bw keV wide at x
  {
    const double u = (x - 0.5 * (w0 + w1)) / (0.5 * (w1 - w0));
    double v = 0.;
    for (size_t k = 0; k < bg.size(); ++k) v += bg[k] * Legendre(k, u);
    v *= bw / kRefKeV;   // the background per 2 keV
    if (withPeaks)
      for (const FitPeak& p : peaks) v += p.area * GaussIn(x - 0.5 * bw, x + 0.5 * bw, p.mu, p.sigma);
    return v;
  }
};
struct Peak {   // a reported peak
  int number = 0;
  FitPeak fit;
  double dmuDay = 0., dareaDay = 0., zDay = 0., zSim = 0., sb = 0.;
  int section = 0;
  std::vector<std::pair<int, double>> shares;   // part index -> share of the area
  std::string latex, text, source;
  double eLine = 0., emitted = 0., own = 0.;
  int kind = kindNone;
  bool broad = false;   // free width over 1.5 resolutions, by more than 2 errors
  int z = 0, up = 0, lo = 0;   // a line's element and transition
};
}  // namespace peakfit

void beam_steel_peak_fit(const char* dir, const char* outDir = "", double zMin = 5., double eMin = 20., double eMax = 7000.)
{
  using namespace peakfit;
  gStyle->SetOptStat(0);
  TF1::DefaultAddToGlobalList(false);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(200000);
  const std::string base = dir;
  const std::string out = std::string(outDir).empty() ? base + "/analysis" : std::string(outDir);
  const std::string inName = base + "/analysis/beam_steel_trace_spectrum.root";
  std::unique_ptr<TFile> in(TFile::Open(inName.c_str()));
  if (!in || in->IsZombie()) { printf("cannot open %s: run analysis/beam_steel_trace_spectrum.C first\n", inName.c_str()); return; }
  gSystem->mkdir(out.c_str(), true);
  auto Get = [&](const std::string& name) -> TH1D* {
    auto h = dynamic_cast<TH1D*>(in->Get(name.c_str()));
    if (!h) { printf("no %s in %s\n", name.c_str(), inName.c_str()); return nullptr; }
    h = dynamic_cast<TH1D*>(h->Clone((name + "_pf").c_str()));
    h->SetDirectory(nullptr);
    return h;
  };
  std::string settings;
  if (auto s = dynamic_cast<TNamed*>(in->Get("settings"))) settings = s->GetTitle();
  auto Setting = [&](const char* key, double def) {
    const auto p = settings.find(std::string(key) + " ");
    double v = def;
    if (p != std::string::npos) sscanf(settings.c_str() + p + std::strlen(key) + 1, "%lf", &v);
    return v;
  };
  gS1 = Setting("sigma1MeV", gS1);
  const double rate = Setting("rate", 1.32e5), promptNs = Setting("promptNs", 10.), seconds = Setting("seconds", 86400.);
  TH1D* fine = Get("view_prompt");   // in the input's bins, for the plots
  if (!fine) return;
  gFineKeV = fine->GetBinWidth(1);
  std::vector<double> edges = {fine->GetXaxis()->GetXmin()};   // the fit's bins (BinWidth): whole input bins merged
  while (edges.back() < fine->GetXaxis()->GetXmax() - 1e-6) edges.push_back(edges.back() + BinWidth(edges.back() + 1e-6));
  auto Merge = [&](TH1D* h) {
    auto m = dynamic_cast<TH1D*>(h->Rebin(edges.size() - 1, (std::string(h->GetName()) + "_v").c_str(), edges.data()));
    m->SetDirectory(nullptr);
    return m;
  };
  auto GetMerged = [&](const std::string& name) -> TH1D* {
    TH1D* h = Get(name);
    if (!h) return nullptr;
    TH1D* m = Merge(h);
    delete h;
    return m;
  };
  TH1D* total = Merge(fine);

  // the parts of the spectrum: each element's X-rays, the other muonic X-rays, the beam's other entries, the accidentals
  struct Part { std::string name; TH1D* h; int z; };
  std::vector<Part> parts;
  TH1D* beamOther = nullptr;
  for (const ElementInfo& el : kElements) {
    TH1D* x = GetMerged(std::string("day_") + el.name + "_promptx");
    TH1D* all = GetMerged(std::string("day_") + el.name + "_prompt");
    if (!x || !all) return;
    parts.push_back({el.name, x, el.z});
    all->Add(x, -1.);
    if (!beamOther) beamOther = all;
    else { beamOther->Add(all); delete all; }
  }
  TH1D* otherX = GetMerged("day_other_promptx");
  TH1D* otherAll = GetMerged("day_other_prompt");
  TH1D* acc = GetMerged("accidental_prompt");
  if (!otherX || !otherAll || !acc) return;
  otherAll->Add(otherX, -1.);
  beamOther->Add(otherAll);
  const int iOtherX = parts.size();
  parts.push_back({"other mu-", otherX, 0});
  const int iBeam = parts.size();
  parts.push_back({"beam", beamOther, -1});
  const int iAcc = parts.size();
  parts.push_back({"accidentals", acc, -2});
  double maxDev = 0.;
  for (int b = total->FindBin(eMin); b <= total->FindBin(eMax); ++b) {
    double s = 0.;
    for (const Part& p : parts) s += p.h->GetBinContent(b);
    if (total->GetBinContent(b) > 0.) maxDev = std::max(maxDev, std::fabs(s / total->GetBinContent(b) - 1.));
  }
  // the emitted lines (truth) and capture gammas
  std::map<int, std::vector<std::pair<double, double>>> lines;   // z -> (keV, photons per day)
  std::map<int, double> capturesPerDay;
  if (auto t = dynamic_cast<TTree*>(in->Get("lines"))) {
    int z = 0;
    double e = 0., perCapture = 0., perDay = 0.;
    t->SetBranchAddress("z", &z); t->SetBranchAddress("e_keV", &e);
    t->SetBranchAddress("per_capture", &perCapture); t->SetBranchAddress("per_day", &perDay);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      lines[z].push_back({e, perDay});
      if (perCapture > 0.) capturesPerDay[z] = perDay / perCapture;
    }
  }
  std::map<int, TH1D*> gammas;
  for (const ElementInfo& el : kElements) gammas[el.z] = Get(std::string("day_") + el.name + "_gammas");

  // the bins
  const int nb = total->GetNbinsX();
  std::vector<double> X(nb + 2, 0.), Y(nb + 2, 0.), E2(nb + 2, 0.), L(nb + 2, 0.), W(nb + 2, 0.);   // center, content, error^2, low edge, width
  double e2Floor = 1e300;
  for (int b = 1; b <= nb; ++b) {
    X[b] = total->GetBinCenter(b); Y[b] = total->GetBinContent(b); E2[b] = std::pow(total->GetBinError(b), 2);
    L[b] = total->GetBinLowEdge(b); W[b] = total->GetBinWidth(b);
    if (E2[b] > 0.) e2Floor = std::min(e2Floor, E2[b]);
  }
  for (int b = 1; b <= nb; ++b) if (E2[b] <= 0.) E2[b] = e2Floor;
  auto First = [&](double e) { const int b = total->FindBin(e); return X[b] < e ? b + 1 : b; };   // the bins with centers in [e, f]
  auto Last = [&](double f) { const int b = total->FindBin(f); return X[b] > f ? b - 1 : b; };
  const int bFirst = First(eMin), bLast = Last(eMax);

  // 1. the search: candidates (keV, area, significance) where a Gaussian of the resolution stands out of a straight line
  auto Scan = [&](const std::vector<double>& y, double from, double to, double zCut, std::vector<std::array<double, 3>>& found) {
    std::vector<double> z(nb + 2, 0.), a(nb + 2, 0.);
    const int b0 = std::max(bFirst, total->FindBin(from)), b1 = std::min(bLast, total->FindBin(to));
    for (int b = b0; b <= b1; ++b) {
      const double mu = X[b], s = Res(mu);
      const double reach = std::max(3., std::ceil(4. * s / W[b])) * W[b] + 1e-6;   // +-4 sigma in whole bins, at least 3
      std::vector<std::vector<double>> cols(3);
      std::vector<double> yy, ww;
      for (int i = total->FindBin(mu - reach); i <= total->FindBin(mu + reach); ++i) {
        if (i < bFirst || i > bLast || std::fabs(X[i] - mu) > reach) continue;
        cols[0].push_back(GaussIn(L[i], L[i] + W[i], mu, s));
        cols[1].push_back(W[i] / kRefKeV);   // the straight line per 2 keV
        cols[2].push_back(W[i] / kRefKeV * (X[i] - mu) / s);
        yy.push_back(y[i]);
        ww.push_back(1. / E2[i]);
      }
      std::vector<double> c;
      TMatrixDSym cov(3);
      if (yy.size() < 5 || !LinearFit(cols, yy, ww, c, &cov) || cov(0, 0) <= 0.) continue;
      a[b] = c[0];
      z[b] = c[0] / std::sqrt(cov(0, 0));
    }
    for (int b = b0; b <= b1; ++b) {
      if (z[b] < zCut || a[b] <= 0.) continue;
      const double r = std::max(1., double(std::lround(0.75 * Res(X[b]) / W[b]))) * W[b] + 1e-6;   // the highest within r
      bool top = true;
      for (int i = std::max(b0, total->FindBin(X[b] - r)); i <= std::min(b1, total->FindBin(X[b] + r)) && top; ++i)
        if (i != b && (z[i] > z[b] || (z[i] == z[b] && i < b))) top = false;
      if (top) found.push_back({X[b], a[b], z[b]});
    }
  };
  std::vector<Candidate> candidates;
  {
    std::vector<std::array<double, 3>> found;
    Scan(Y, eMin, eMax, zMin, found);
    for (const auto& f : found) { Candidate c; c.mu = f[0]; c.area = f[1]; c.z = f[2]; candidates.push_back(c); }
  }
  const size_t nSearch = candidates.size();

  // 2. the fits, section by section
  std::vector<std::pair<double, double>> cores;
  for (double c0 = eMin; c0 < eMax - 1.;) {
    double c1 = c0 + 30. * Res(c0 + 15. * Res(c0));
    if (eMax - c1 < 15. * Res(eMax)) c1 = eMax;
    else {   // the boundary moved into the widest gap between the candidates within +-5 resolutions (not on a peak)
      const double s = Res(c1);
      double best = c1, bestGap = -1.;
      for (double x = c1 - 5. * s; x <= c1 + 5. * s; x += 0.25 * s) {
        double gap = 1e300;
        for (const Candidate& c : candidates) gap = std::min(gap, std::fabs(c.mu - x) / Res(c.mu));
        if (gap > bestGap) { bestGap = gap; best = x; }
      }
      c1 = best;
    }
    cores.push_back({c0, c1});
    c0 = c1;
  }
  const int bgDeg = 3, nbg = bgDeg + 1;
  // fit the peaks of a window: the free ones (area, position, width at the resolution or free), the fixed ones as given
  auto FitWindow = [&](Section& sec, bool freeWidth) -> bool {
    std::vector<int> idx;
    for (size_t k = 0; k < sec.peaks.size(); ++k) if (sec.peaks[k].free) idx.push_back(k);
    const int per = freeWidth ? 3 : 2, npar = nbg + per * idx.size();
    const double wc = 0.5 * (sec.w0 + sec.w1), wh = 0.5 * (sec.w1 - sec.w0);
    std::vector<FitPeak>& pk = sec.peaks;
    auto model = [&pk, wc, wh, per, freeWidth, nbg](double* xx, double* p) {
      const double x = xx[0], u = (x - wc) / wh, bw = BinWidth(x);   // x: a bin's center
      double v = 0.;
      for (int k = 0; k < nbg; ++k) v += p[k] * Legendre(k, u);
      v *= bw / kRefKeV;   // the background per 2 keV
      int q = nbg;
      for (const FitPeak& f : pk) {
        double a = f.area, m = f.mu, s = f.sigma;
        if (f.free) { a = p[q]; m = p[q + 1]; s = freeWidth ? p[q + 2] : Res(m); q += per; }
        v += a * GaussIn(x - 0.5 * bw, x + 0.5 * bw, m, s);
      }
      return v;
    };
    // start: a linear fit of the areas and the background at the given positions and widths
    const int b0 = First(sec.w0), b1 = Last(sec.w1);
    std::vector<std::vector<double>> cols(nbg + idx.size());
    std::vector<double> yy, ww;
    for (int b = b0; b <= b1; ++b) {
      double fixedPart = 0.;
      for (const FitPeak& f : pk) if (!f.free) fixedPart += f.area * GaussIn(L[b], L[b] + W[b], f.mu, f.sigma);
      for (int k = 0; k < nbg; ++k) cols[k].push_back(W[b] / kRefKeV * Legendre(k, (X[b] - wc) / wh));
      for (size_t j = 0; j < idx.size(); ++j) {
        const FitPeak& f = pk[idx[j]];
        cols[nbg + j].push_back(GaussIn(L[b], L[b] + W[b], f.mu, freeWidth ? f.sigma : Res(f.mu)));
      }
      yy.push_back(Y[b] - fixedPart);
      ww.push_back(1. / E2[b]);
    }
    std::vector<double> start;
    if (!LinearFit(cols, yy, ww, start)) return false;
    TF1 f("peakfit_model", model, sec.w0, sec.w1, npar);
    for (int k = 0; k < nbg; ++k) f.SetParameter(k, start[k]);
    for (size_t j = 0; j < idx.size(); ++j) {
      const FitPeak& p = pk[idx[j]];
      const int q = nbg + per * j;
      const double m0 = p.mu0 > 0. ? p.mu0 : p.mu, s = Res(m0);
      f.SetParameter(q, std::max(start[nbg + j], 1e-3 * std::fabs(start[0]) * BinWidth(wc) / kRefKeV + 1.));   // 1e-3 of a bin's background
      f.SetParameter(q + 1, std::clamp(p.mu, m0 - 2.5 * s, m0 + 2.5 * s));
      f.SetParLimits(q + 1, m0 - 2.5 * s, m0 + 2.5 * s);
      if (freeWidth) { f.SetParameter(q + 2, p.sigma > 0. ? p.sigma : s); f.SetParLimits(q + 2, 0.4 * s, 4. * s); }
    }
    TFitResultPtr r = total->Fit(&f, "QRSN0");
    if (!r.Get()) return false;
    for (int k = 0; k < npar; ++k) if (!std::isfinite(f.GetParameter(k))) return false;
    if (!freeWidth) sec.status = r->Status();   // Minuit2: 0 fine, 1 covariance forced positive, 2-5 a problem (printed)
    sec.bg.assign(nbg, 0.);
    for (int k = 0; k < nbg; ++k) sec.bg[k] = f.GetParameter(k);
    for (size_t j = 0; j < idx.size(); ++j) {
      FitPeak& p = pk[idx[j]];
      const int q = nbg + per * j;
      if (freeWidth) {
        p.areaFree = f.GetParameter(q); p.dareaFree = f.GetParError(q);
        p.sFree = f.GetParameter(q + 2); p.dsFree = f.GetParError(q + 2);
      } else {
        p.area = f.GetParameter(q); p.darea = f.GetParError(q);
        p.mu = f.GetParameter(q + 1); p.dmu = f.GetParError(q + 1);
        p.sigma = Res(p.mu);
      }
    }
    if (!freeWidth) { sec.chi2 = r->Chi2(); sec.ndf = r->Ndf(); }
    return true;
  };

  std::vector<Section> sections;
  std::vector<Peak> peaks;
  int nAdded = 0, nDropped = 0, nMerged = 0;
  // two free peaks that converged within half a resolution: one, their sum (the weaker's candidate dropped)
  auto MergeClose = [&](Section& sec) {
    for (size_t a = 0; a < sec.peaks.size(); ++a)
      for (size_t b = a + 1; b < sec.peaks.size(); ++b) {
        FitPeak& p = sec.peaks[a];
        const FitPeak& q = sec.peaks[b];
        if (!p.free || !q.free || std::fabs(p.mu - q.mu) > 0.5 * Res(p.mu)) continue;
        if (p.area + q.area > 0. && p.area > 0. && q.area > 0.) p.mu = (p.mu * p.area + q.mu * q.area) / (p.area + q.area);
        p.area += q.area;
        candidates[q.candidate].alive = false;
        sec.peaks.erase(sec.peaks.begin() + b);
        ++nMerged;
        return true;
      }
    return false;
  };
  for (const auto& [c0, c1] : cores) {
    Section sec;
    sec.c0 = c0; sec.c1 = c1;
    sec.w0 = std::max(eMin, c0 - 6. * Res(c0));
    sec.w1 = std::min(eMax, c1 + 6. * Res(c1));
    for (size_t i = 0; i < candidates.size(); ++i) {
      const Candidate& c = candidates[i];
      if (!c.alive) continue;
      const double s = Res(c.mu);
      FitPeak p;
      p.mu = c.mu; p.area = c.area; p.sigma = s; p.candidate = i; p.mu0 = c.mu;
      if (c.mu >= sec.w0 && c.mu <= sec.w1) sec.peaks.push_back(p);
      else if (c.mu > sec.w0 - 4. * s && c.mu < sec.w1 + 4. * s) { p.free = false; sec.peaks.push_back(p); }
    }
    bool ok = FitWindow(sec, false);
    // residual peaks
    for (int iter = 0; ok && iter < 8; ++iter) {
      std::vector<double> res(nb + 2, 0.);
      for (int b = First(sec.w0); b <= Last(sec.w1); ++b) res[b] = Y[b] - sec.Model(X[b], W[b]);
      std::vector<std::array<double, 3>> found;
      Scan(res, sec.w0 + 2. * Res(sec.w0), sec.w1 - 2. * Res(sec.w1), zMin, found);
      std::sort(found.begin(), found.end(), [](const auto& a, const auto& b) { return a[2] > b[2]; });
      bool added = false;
      for (const auto& c : found) {
        bool near = false;
        for (const FitPeak& p : sec.peaks) if (std::fabs(p.mu - c[0]) < 0.8 * Res(c[0])) near = true;
        if (near) continue;
        Candidate nc;
        nc.mu = c[0]; nc.area = c[1]; nc.z = c[2];
        candidates.push_back(nc);
        FitPeak p;
        p.mu = c[0]; p.area = c[1]; p.sigma = Res(c[0]); p.candidate = candidates.size() - 1; p.mu0 = c[0];
        sec.peaks.push_back(p);
        ++nAdded;
        added = true;
        break;
      }
      if (!added) break;
      ok = FitWindow(sec, false);
      while (ok && MergeClose(sec)) ok = FitWindow(sec, false);
    }
    while (ok && MergeClose(sec)) ok = FitWindow(sec, false);
    // drop the free peaks under 3 sigma (the weakest first)
    while (ok) {
      int weakest = -1;
      double zWeak = 3.;
      for (size_t k = 0; k < sec.peaks.size(); ++k) {
        const FitPeak& p = sec.peaks[k];
        const double z = p.darea > 0. ? p.area / p.darea : 0.;
        if (p.free && z < zWeak) { zWeak = z; weakest = k; }
      }
      if (weakest < 0) break;
      candidates[sec.peaks[weakest].candidate].alive = false;
      sec.peaks.erase(sec.peaks.begin() + weakest);
      ++nDropped;
      ok = FitWindow(sec, false);
      while (ok && MergeClose(sec)) ok = FitWindow(sec, false);
    }
    if (!ok) { printf("section %.1f-%.1f keV: the fit failed, left out\n", c0, c1); continue; }
    // the next sections see the fitted positions and areas
    for (const FitPeak& p : sec.peaks)
      if (p.free) { candidates[p.candidate].mu = p.mu; candidates[p.candidate].area = p.area; }
    Section free = sec;
    if (FitWindow(free, true))
      for (size_t k = 0; k < sec.peaks.size(); ++k) {
        sec.peaks[k].sFree = free.peaks[k].sFree; sec.peaks[k].dsFree = free.peaks[k].dsFree;
        sec.peaks[k].areaFree = free.peaks[k].areaFree; sec.peaks[k].dareaFree = free.peaks[k].dareaFree;
      }
    // each part's share of each peak: the same Gaussians and background, fitted (linearly) to the part
    const int b0 = First(sec.w0), b1 = Last(sec.w1);
    const double wc = 0.5 * (sec.w0 + sec.w1), wh = 0.5 * (sec.w1 - sec.w0);
    std::vector<std::vector<double>> cols(nbg + sec.peaks.size());
    std::vector<double> ww;
    double sumE2 = 0., sumY = 0.;
    for (int b = b0; b <= b1; ++b) {
      for (int k = 0; k < nbg; ++k) cols[k].push_back(W[b] / kRefKeV * Legendre(k, (X[b] - wc) / wh));
      for (size_t j = 0; j < sec.peaks.size(); ++j)
        cols[nbg + j].push_back(GaussIn(L[b], L[b] + W[b], sec.peaks[j].mu, sec.peaks[j].sigma));
      ww.push_back(1. / E2[b]);
      if (Y[b] > 0.) { sumE2 += E2[b]; sumY += Y[b]; }
    }
    std::vector<std::vector<double>> amp(parts.size());
    for (size_t ip = 0; ip < parts.size(); ++ip) {
      std::vector<double> yy;
      for (int b = b0; b <= b1; ++b) yy.push_back(parts[ip].h->GetBinContent(b));
      LinearFit(cols, yy, ww, amp[ip]);
    }
    const double weight = sumY > 0. ? sumE2 / sumY : 1.;   // the mean weight of an entry: a day's errors are sqrt(w) smaller
    for (size_t j = 0; j < sec.peaks.size(); ++j) {
      const FitPeak& p = sec.peaks[j];
      // reported once, by the section whose core holds it, and not where another section already reported a peak within
      // half a resolution (the same peak, found twice)
      if (!p.free || p.mu < c0 || p.mu >= c1 || candidates[p.candidate].reported) continue;
      candidates[p.candidate].reported = true;
      bool twice = false;
      for (const Peak& q : peaks) if (std::fabs(q.fit.mu - p.mu) < 0.5 * p.sigma) twice = true;
      if (twice) continue;
      Peak pk;
      pk.fit = p;
      pk.section = sections.size();
      pk.dmuDay = p.dmu / std::sqrt(weight);
      pk.dareaDay = p.darea / std::sqrt(weight);
      pk.zSim = p.darea > 0. ? p.area / p.darea : 0.;
      pk.zDay = pk.dareaDay > 0. ? p.area / pk.dareaDay : 0.;
      double in = 0.;
      for (int b = total->FindBin(p.mu - 1.5 * p.sigma); b <= total->FindBin(p.mu + 1.5 * p.sigma); ++b) in += Y[b];
      const double sIn = p.area * GaussIn(p.mu - 1.5 * p.sigma, p.mu + 1.5 * p.sigma, p.mu, p.sigma);
      pk.sb = in - sIn > 0. ? sIn / (in - sIn) : 0.;
      double sum = 0.;
      for (size_t ip = 0; ip < parts.size(); ++ip) if (!amp[ip].empty()) sum += amp[ip][nbg + j];
      for (size_t ip = 0; ip < parts.size(); ++ip)
        if (!amp[ip].empty() && sum != 0.) pk.shares.push_back({int(ip), amp[ip][nbg + j] / sum});
      std::sort(pk.shares.begin(), pk.shares.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
      peaks.push_back(pk);
    }
    sections.push_back(sec);
  }
  std::sort(peaks.begin(), peaks.end(), [](const Peak& a, const Peak& b) { return a.fit.mu < b.fit.mu; });

  // 4. the names
  struct Hyp {
    double score = 0.;
    int kind = kindNone;
    std::string latex, text;
    double eLine = 0.;   // the energy the peak should have under this name
    int z = 0, up = 0, lo = 0;
  };
  auto Find = [&](int z) -> const ElementInfo* {
    for (const ElementInfo& e : kElements) if (e.z == z) return &e;
    return nullptr;
  };
  // the best explanation of a peak at mu (resolution s, free width / s) by one element's X-rays, each candidate scored by
  // its emitted photons per day x the rough efficiency x its chance (escape, pair escape, sum, structure) x its closeness
  auto Explain = [&](const ElementInfo& el, double mu, double s, double width, bool structuresOnly = false) {
    Hyp best;
    const auto& L = lines[el.z];
    const double caps = capturesPerDay.count(el.z) ? capturesPerDay[el.z] : 1.;
    auto Prox = [&](double e, double w) {
      const double d = (e - mu) / (w * s);
      return std::fabs(e - mu) < 2. * std::max(1., w) * s ? std::exp(-0.5 * d * d) : 0.;
    };
    auto Consider = [&](double score, int kind, const std::string& latex, const std::string& text, double eLine, int up = 0, int lo = 0) {
      if (structuresOnly && kind != kindStructure) return;
      if (score > best.score) {
        best.score = score; best.kind = kind; best.latex = latex; best.text = text; best.eLine = eLine;
        best.z = el.z; best.up = up; best.lo = lo;
      }
    };
    const std::string n = el.name;
    static const double wk[2] = {2., 1.}, wl[3] = {9., 5., 1.}, wm[3] = {20., 14., 1.};
    // the energy a line's peak should have: an unresolved measured doublet (components under 0.75 sigma apart) its mean
    auto Reference = [&](double e, int up, int lo) {
      const double* comp = up == 2 && lo == 1 ? el.ka : up == 3 && lo == 2 ? el.la : up == 4 && lo == 3 && el.ma[0] > 0. ? el.ma : nullptr;
      if (!comp || std::fabs(comp[0] - comp[1]) >= 0.75 * Res(e)) return e;
      return up == 2 ? WeightedMean(el.ka, wk, 2) : up == 3 ? WeightedMean(el.la, wl, 3) : WeightedMean(el.ma, wm, 3);
    };
    double wMax = 0., escapeTotal = 0.;
    for (const auto& [e, pd] : L) wMax = std::max(wMax, pd * Eff(e));
    struct Strong { double e, w; int up, lo; };
    std::vector<Strong> strong;   // the strong lines: keV, photons x efficiency, the transition
    for (const auto& [e, pd] : L) {
      const double w = pd * Eff(e);
      int up = 0, lo = 0;
      Transition(el, e, up, lo);
      if (w >= 0.03 * wMax) strong.push_back({e, w, up, lo});
      Consider(w * Prox(e, 0.5), kindLine, n + " " + LineName(el, e, true), n + " " + LineName(el, e, false), Reference(e, up, lo), up, lo);
      for (const KXray& k : kCdTe) {
        if (e <= k.edge) continue;
        escapeTotal += w * EscapeProb(e) * k.weight;
        Consider(w * EscapeProb(e) * k.weight * Prox(e - k.e, 0.5), kindEscape, n + " " + LineName(el, e, true) + " esc.",
                 Form("%s escape of %s %s (%.2f keV)", k.text, el.name, LineName(el, e, false).c_str(), e), e - k.e);
      }
      if (e > 1100.) {
        Consider(w * 0.02 * Prox(e - 511., 0.5), kindPairEscape, n + " " + LineName(el, e, true) + " SE",
                 Form("single escape (-511 keV) of %s %s (%.2f keV)", el.name, LineName(el, e, false).c_str(), e), e - 511.);
        Consider(w * 0.005 * Prox(e - 1022., 0.5), kindPairEscape, n + " " + LineName(el, e, true) + " DE",
                 Form("double escape (-1022 keV) of %s %s (%.2f keV)", el.name, LineName(el, e, false).c_str(), e), e - 1022.);
      }
    }
    // two X-rays of one muon in one pixel: transitions one after the other in the cascade (0.3: the same pixel, both whole);
    // where the sum is a transition itself (K_alpha + L_alpha = K_beta: the levels add up), both are named
    for (size_t a = 0; a < strong.size(); ++a)
      for (size_t b = a + 1; b < strong.size(); ++b) {
        const Strong& p = strong[a].up >= strong[b].up ? strong[a] : strong[b];   // the upper transition
        const Strong& q = &p == &strong[a] ? strong[b] : strong[a];
        if (p.lo < q.up) continue;   // not one after the other in one cascade
        const double e = p.e + q.e, score = p.w * q.w / caps * 0.3 * Prox(e, 0.5);
        if (score <= 0.) continue;
        double wDirect = 0., eDirect = 0.;
        for (const auto& [e2, pd2] : L)
          if (std::fabs(e2 - e) < 0.5 * s && pd2 * Eff(e2) > wDirect) { wDirect = pd2 * Eff(e2); eDirect = e2; }
        const double direct = eDirect > 0. ? wDirect * Prox(eDirect, 0.5) : 0., share = score / (score + direct);
        const std::string sumLatex = LineName(el, q.e, true) + "+" + LineName(el, p.e, true);
        const std::string sumText = LineName(el, q.e, false) + " + " + LineName(el, p.e, false);
        if (eDirect > 0. && share >= 0.2) {
          int up = 0, lo = 0;
          Transition(el, eDirect, up, lo);
          Consider(score + direct, kindSum, n + " " + LineName(el, eDirect, true) + "/" + sumLatex,
                   Form("%s %s and the sum %s in one pixel (the same energy; the sum about %.0f%%)", el.name,
                        LineName(el, eDirect, false).c_str(), sumText.c_str(), 100. * share), e, up, lo);
        } else if (eDirect <= 0.)
          Consider(score, kindSum, n + " " + sumLatex, Form("sum in one pixel: %s %s", el.name, sumText.c_str()), e);
      }
    if (mu < 35.)   // the escaping K X-rays absorbed in a neighbouring pixel (about half of them)
      for (const KXray& k : kCdTe)
        Consider(escapeTotal * k.weight * 0.5 * Prox(k.e, 0.5), kindFluorescence, k.latex,
                 Form("%s from a neighbouring pixel (escapes of %s X-rays)", k.text, el.name), k.e);
    for (const Strong& st : strong) {   // the continuum's structures: Compton edge, backscatter, the Compton shoulder
      // x the Compton share of the interactions, 1 / (1 + (150 keV / E)^3): below ~150 keV photo-absorption wins
      const double e = st.e, w = st.w / (1. + std::pow(150. / e, 3));
      const double eC = e * (1. - 1. / (1. + 2. * e / 511.)), eB = e / (1. + 2. * e / 511.);
      Consider(w * 0.01 * Prox(eC, 1.), kindStructure, n + " " + LineName(el, e, true) + " Compton edge",
               Form("Compton edge (%.0f keV) of %s %s", eC, el.name, LineName(el, e, false).c_str()), eC);
      Consider(w * 0.01 * Prox(eB, std::max(1., width)), kindStructure, n + " " + LineName(el, e, true) + " backscatter",
               Form("backscatter (%.0f keV) of %s %s", eB, el.name, LineName(el, e, false).c_str()), eB);
      if (mu > eC - 150. && mu < eC)
        Consider(w * 0.0005, kindStructure, n + " " + LineName(el, e, true) + " Compton",
                 Form("Compton continuum below the edge (%.0f keV) of %s %s", eC, el.name, LineName(el, e, false).c_str()), 0.);
    }
    return best;
  };
  for (size_t i = 0; i < peaks.size(); ++i) {
    Peak& pk = peaks[i];
    pk.number = i + 1;
    const double mu = pk.fit.mu, s = pk.fit.sigma;
    const double width = pk.fit.sFree > 0. ? pk.fit.sFree / s : 1.;
    pk.broad = pk.fit.sFree > 1.5 * s && pk.fit.sFree - 1.5 * s > 2. * pk.fit.dsFree;
    if (pk.shares.empty()) { pk.latex = pk.text = "?"; pk.kind = kindStructure; continue; }
    const int top = pk.shares[0].first;
    const Part& part = parts[top];
    pk.own = pk.fit.area * pk.shares[0].second;
    if (part.z > 0) {   // an element's muonic X-rays
      const ElementInfo* el = Find(part.z);
      pk.source = std::string(el->name) + " X-rays";
      Hyp h = Explain(*el, mu, s, width);
      if (h.kind == kindNone) {
        h.kind = kindStructure;
        h.latex = std::string(el->name) + " X-rays";
        h.text = Form("%s X-rays: no line, escape, sum or edge within 2 sigma", el->name);
      }
      // the plot's short label: the name before a "/" (a line that is also a sum: marked *)
      auto Short = [](const std::string& latex) { const auto p = latex.find('/'); return p == std::string::npos ? latex : latex.substr(0, p) + "*"; };
      std::string shortLatex = Short(h.latex);
      // a blend: a second element with 15% of the peak or more and a line (escape, sum) there, appended; with 30% or
      // more it is named first when the largest part offers only a structure
      if (pk.shares.size() > 1 && pk.shares[1].second >= 0.15 && parts[pk.shares[1].first].z > 0) {
        const ElementInfo* el2 = Find(parts[pk.shares[1].first].z);
        Hyp h2 = Explain(*el2, mu, s, width);
        const bool line1 = h.kind >= kindLine && h.kind <= kindSum, line2 = h2.kind >= kindLine && h2.kind <= kindSum;
        if (!line1 && line2 && pk.shares[1].second >= 0.3) {
          shortLatex = Short(h2.latex) + " +" + el->name;
          h2.latex += " + " + h.latex; h2.text += " + " + h.text; h = h2;
        } else if (line1 && line2) {
          shortLatex += std::string(" +") + el2->name;
          h.latex += " + " + h2.latex; h.text += " + " + h2.text;
        }
      }
      // wider than two resolutions (by more than 2 errors, or at the free width's limit of 4): a structure, whatever line
      // lies there
      if (h.kind != kindStructure && (pk.fit.sFree - 2. * s > 2. * pk.fit.dsFree || pk.fit.sFree > 3.9 * s)) {
        h.kind = kindStructure;
        h.text = Form("broad structure (free width %.1f keV, %.1f resolutions) at %s", pk.fit.sFree, width, h.text.c_str());
        h.latex = "broad: " + h.latex;
        shortLatex = "broad: " + shortLatex;
      }
      pk.kind = h.kind; pk.latex = shortLatex; pk.text = h.text; pk.eLine = h.eLine;
      pk.z = h.z; pk.up = h.up; pk.lo = h.lo;
      if (h.kind == kindLine)
        for (const auto& [e, pd] : lines[part.z]) if (std::fabs(e - mu) < 1.5 * s) pk.emitted += pd;
    } else if (top == iBeam || top == iAcc) {   // the beam's other entries or the accidentals: known lines, capture gammas
      pk.source = top == iBeam ? "beam, not X-rays" : "accidentals";
      double dBest = 2.;
      const KnownLine* known = nullptr;
      for (const KnownLine& kl : kKnown) if (std::fabs(kl.e - mu) / s < dBest) { dBest = std::fabs(kl.e - mu) / s; known = &kl; }
      double gBest = 0., gE = 0.;
      const char* gEl = nullptr;
      if (top == iBeam)
        for (const ElementInfo& el : kElements) {
          TH1D* g = gammas[el.z];
          if (!g) continue;
          auto Sum = [&](double a, double b) { return g->Integral(g->FindBin(a), g->FindBin(b)); };
          const double inW = Sum(mu - 1.5 * s, mu + 1.5 * s);
          const double side = Sum(mu - 4. * s, mu - 2.5 * s) + Sum(mu + 2.5 * s, mu + 4. * s);   // 3 sigma of sidebands
          if (inW - side > gBest) {
            gBest = inW - side;
            gEl = el.name;
            double sw = 0., swx = 0.;
            for (int b = g->FindBin(mu - 1.5 * s); b <= g->FindBin(mu + 1.5 * s); ++b) { sw += g->GetBinContent(b); swx += g->GetBinContent(b) * g->GetBinCenter(b); }
            gE = sw > 0. ? swx / sw : mu;
          }
        }
      // a capture gamma when its emitted line (truth) holds more than the peak (detection takes far less), else a known line
      if (gEl && gBest > pk.own) {
        pk.kind = kindGamma;
        pk.latex = Form("#mu^{-}%s #gamma", gEl);
        pk.text = Form("gamma after mu- capture on %s (emitted line at %.1f keV, %.3g per day)", gEl, gE, gBest);
        pk.eLine = gE;
      } else if (known) {
        pk.kind = kindGamma;
        pk.latex = known->latex;
        pk.text = known->text;
        pk.eLine = known->e;
      } else {
        pk.kind = kindStructure;
        pk.latex = top == iBeam ? "beam" : "accidental";
        pk.text = top == iBeam ? "the beam's other entries, no known line" : "accidentals, no known line";
      }
    } else {
      pk.kind = kindLine;
      pk.source = "other muonic X-rays";
      pk.latex = "other #mu^{-}";
      pk.text = "muonic X-rays of muons outside the plate's atoms";
    }
    if (pk.broad && pk.kind != kindStructure) pk.text += " (broad: blend or structure)";
  }
  // an unresolved doublet that the fit split in two (two peaks under one resolution apart, both one measured transition of
  // one element): the lower peak named by the lower component, the upper by the upper
  for (size_t i = 0; i + 1 < peaks.size(); ++i) {
    Peak& a = peaks[i];
    Peak& b = peaks[i + 1];
    if (a.kind != kindLine || b.kind != kindLine || a.z == 0 || a.z != b.z || a.up != b.up || a.lo != b.lo) continue;
    if (b.fit.mu - a.fit.mu > a.fit.sigma) continue;
    const ElementInfo* el = Find(a.z);
    const double* comp = a.up == 2 && a.lo == 1 ? el->ka : a.up == 3 && a.lo == 2 ? el->la : a.up == 4 && a.lo == 3 && el->ma[0] > 0. ? el->ma : nullptr;
    if (!comp) continue;
    static const char* states[3][2] = {{"2p3/2->1s", "2p1/2->1s"}, {"3d5/2->2p3/2", "3d3/2->2p1/2"}, {"4f7/2->3d5/2", "4f5/2->3d3/2"}};
    const int iLo = comp[0] < comp[1] ? 0 : 1, iHi = 1 - iLo;
    const char series = "KLM"[a.lo - 1];
    a.latex = Form("%s %c_{#alpha%d}", el->name, series, iLo + 1);
    a.text = Form("%s %c_alpha%d (%s): the lower half of an unresolved doublet (#%d the upper)", el->name, series, iLo + 1, states[a.lo - 1][iLo], b.number);
    a.eLine = comp[iLo];
    b.latex = Form("%s %c_{#alpha%d}", el->name, series, iHi + 1);
    b.text = Form("%s %c_alpha%d (%s): the upper half of an unresolved doublet (#%d the lower)", el->name, series, iHi + 1, states[a.lo - 1][iHi], a.number);
    b.eLine = comp[iHi];
  }
  // two neighbours under one resolution apart named the same transition of one element, too far apart for its two
  // components (DoubletGap under half their distance): the line and the shape around it (a Compton edge or a backscatter
  // under it; seen in fine bins). The one nearer the line keeps the name, the other takes its best structure.
  for (size_t i = 0; i + 1 < peaks.size(); ++i) {
    Peak& a = peaks[i];
    Peak& b = peaks[i + 1];
    if (a.kind != b.kind || (a.kind != kindLine && a.kind != kindSum) || a.z == 0 || a.z != b.z || a.up == 0 || a.up != b.up || a.lo != b.lo) continue;
    const double gap = b.fit.mu - a.fit.mu;
    const ElementInfo* el = Find(a.z);
    if (!el || gap > a.fit.sigma || DoubletGap(*el, a.up, a.lo) >= 0.5 * gap) continue;
    const bool aFar = std::fabs(a.fit.mu - a.eLine) > std::fabs(b.fit.mu - b.eLine);
    Peak& far = aFar ? a : b;
    const Peak& near = aFar ? b : a;
    const Hyp h = Explain(*el, far.fit.mu, far.fit.sigma, far.fit.sFree > 0. ? far.fit.sFree / far.fit.sigma : 1., true);
    const bool named = h.kind == kindStructure;
    far.kind = kindStructure;
    far.emitted = 0.;
    far.eLine = named ? h.eLine : 0.;
    far.latex = named ? h.latex : std::string(el->name) + " shoulder";
    far.text = Form("%s, beside #%d %s: a second Gaussian for the shape around that line", named ? h.text.c_str() : "shoulder",
                    near.number, near.text.c_str());
  }

  // 5. the lines that should be seen (2026-11-11, the user: "show lines that are supposed to be seen at location with grey
  // color"; grey in the linear plots): each element's emitted lines (truth; within half a resolution merged) x the
  // full-energy efficiency, seen when the search would find it on its own (expected significance >= zMin, as in 1.). The
  // efficiency: the rough one (Eff) x a factor from the fitted main lines (Delta n = 1: K_alpha, L_alpha, M_alpha, N_alpha;
  // not broad): their detected per emitted photon, interpolated in log energy.
  std::vector<std::pair<double, double>> calibration;   // keV -> detected per emitted / Eff
  {
    std::map<std::array<int, 3>, std::pair<double, double>> own;   // (z, up, lo) -> the peaks' own entries, error^2
    for (const Peak& pk : peaks)
      if (pk.kind == kindLine && pk.z > 0 && pk.up == pk.lo + 1 && !pk.broad && pk.fit.area > 0.) {
        auto& o = own[{pk.z, pk.up, pk.lo}];
        o.first += pk.own;
        o.second += std::pow(pk.fit.darea * pk.own / pk.fit.area, 2);
      }
    for (const auto& [key, o] : own) {
      const ElementInfo* el = Find(key[0]);
      double sum = 0., sumE = 0.;
      for (const auto& [e, pd] : lines[key[0]]) {
        int up = 0, lo = 0;
        Transition(*el, e, up, lo);
        if (up == key[1] && lo == key[2]) { sum += pd; sumE += pd * e; }
      }
      if (sum > 0. && o.first > 20. * std::sqrt(o.second)) calibration.push_back({sumE / sum, o.first / sum / Eff(sumE / sum)});
    }
    std::sort(calibration.begin(), calibration.end());
  }
  auto Efficiency = [&](double e) {   // full energy, per emitted photon
    if (calibration.empty()) return 0.;
    double k = e <= calibration.front().first ? calibration.front().second : calibration.back().second;
    for (size_t i = 1; i < calibration.size(); ++i)
      if (e > calibration[i - 1].first && e <= calibration[i].first)
        k = calibration[i - 1].second * std::pow(calibration[i].second / calibration[i - 1].second,
                                                 std::log(e / calibration[i - 1].first) / std::log(calibration[i].first / calibration[i - 1].first));
    return k * Eff(e);
  };
  auto AmplitudeError = [&](double mu) {   // the search's error on a Gaussian's area at mu (its straight line over +-4 sigma)
    const double s = Res(mu), w = BinWidth(mu), reach = std::max(3., std::ceil(4. * s / w)) * w + 1e-6;
    std::vector<std::vector<double>> cols(3);
    std::vector<double> yy, ww, c;
    for (int i = total->FindBin(mu - reach); i <= total->FindBin(mu + reach); ++i) {
      if (i < bFirst || i > bLast || std::fabs(X[i] - mu) > reach) continue;
      cols[0].push_back(GaussIn(L[i], L[i] + W[i], mu, s));
      cols[1].push_back(W[i] / kRefKeV);
      cols[2].push_back(W[i] / kRefKeV * (X[i] - mu) / s);
      yy.push_back(0.);
      ww.push_back(1. / E2[i]);
    }
    TMatrixDSym cov(3);
    return yy.size() >= 5 && LinearFit(cols, yy, ww, c, &cov) && cov(0, 0) > 0. ? std::sqrt(cov(0, 0)) : 0.;
  };
  struct Expected { int z = 0, up = 0, lo = 0, peak = 0; double e = 0., perDay = 0., detected = 0., zSim = 0.; bool carried = false; };
  std::vector<Expected> expected;
  for (const ElementInfo& el : kElements) {
    auto ls = lines[el.z];
    std::sort(ls.begin(), ls.end());
    std::vector<std::pair<double, double>> merged;   // within half a resolution: the photons' mean energy, their number
    for (const auto& [e, pd] : ls)
      if (!merged.empty() && e - merged.back().first < 0.5 * Res(e)) {
        auto& m = merged.back();
        m.first = (m.first * m.second + e * pd) / (m.second + pd);
        m.second += pd;
      } else merged.push_back({e, pd});
    for (const auto& [e, pd] : merged) {
      if (e < eMin || e > eMax) continue;
      Expected x;
      x.z = el.z; x.e = e; x.perDay = pd; x.detected = pd * Efficiency(e);
      const double de = AmplitudeError(e);
      x.zSim = de > 0. ? x.detected / de : 0.;
      if (x.zSim < zMin) continue;
      Transition(el, e, x.up, x.lo);
      // the reported peak it lies under (the nearest within 1.5 resolutions); carried: a peak there has its name
      double best = 1.5 * Res(e);
      for (const Peak& pk : peaks) {
        const double d = std::fabs(pk.fit.mu - e);
        if (d < best) { best = d; x.peak = pk.number; }
        if (d < 1.5 * Res(e) && pk.kind != kindStructure && pk.z == x.z && pk.up == x.up && pk.lo == x.lo) x.carried = true;
      }
      expected.push_back(x);
    }
  }
  std::sort(expected.begin(), expected.end(), [](const Expected& a, const Expected& b) { return a.e < b.e; });

  // the numbers
  std::ostringstream txt, csv;
  txt << "Peaks of the tagged-beam spectrum (hits within " << promptNs << " ns of their muon), one day at " << Form("%.3g", rate)
      << " mu-/s: " << base << "\n"
      << Form("  CZT resolution sigma = %.5g keV x sqrt(E / MeV); search and fits %g-%g keV in bins of %g keV up to 100 keV, %g keV "
              "up to 400 keV, %g keV above; the parts add up to the total within %.1g\n", gS1, eMin, eMax, BinWidth(50.), BinWidth(200.),
              BinWidth(1e4), maxDev)
      << Form("  %zu candidates from the search (significance >= %g in the simulation's statistics), %d added from the residuals, "
              "%d dropped (under 3 sigma), %d merged (two within half a resolution): %zu peaks in %zu sections\n", nSearch, zMin, nAdded,
              nDropped, nMerged, peaks.size(), sections.size())
      << "  columns: fitted energy +- one day's error [the simulation's]; area (entries per day) +- one day [simulation]; the free-"
         "width fit's area; significance in one day; S/B in +-1.5 sigma;\n  the resolution and the free fit's width; the kind and "
         "the name (from the part with the largest share); the line (truth) and the fit minus it; shares; detected per emitted\n"
         "  photon of the element's lines within +-1.5 sigma (lines only). Where the fixed- and free-width areas differ, the fit's "
         "model (not the statistics) limits the area. Kinds: line; Cd/Te escape\n  (a line minus a Cd or Te K X-ray that left the "
         "pixel); pair escape (-511 or -1022 keV); sum (two X-rays of one muon in one pixel); Cd/Te fluorescence (in a neighbouring "
         "pixel);\n  gamma (annihilation, (n,n'), capture); structure (Compton edge, backscatter, Compton continuum: no line)\n\n";
  csv << "number,e_keV,de_day,de_sim,area,darea_day,darea_sim,area_free,z_day,s_over_b,sigma_res,sigma_free,dsigma_free,kind,name,"
         "source,line_keV,e_minus_line,share1,part1,share2,part2,detected_per_emitted,section_chi2,section_ndf\n";
  for (const Peak& pk : peaks) {
    const FitPeak& p = pk.fit;
    std::string sh;
    for (size_t k = 0; k < pk.shares.size() && k < 3; ++k)
      if (k == 0 || pk.shares[k].second >= 0.1) sh += Form("%s%s %.0f%%", k ? ", " : "", parts[pk.shares[k].first].name.c_str(), 100. * pk.shares[k].second);
    const double eff = pk.emitted > 0. ? pk.own / pk.emitted : 0.;
    const Section& sec = sections[pk.section];
    txt << Form("%3d %8.2f +-%6.3f [%5.2f]  %9.3g +-%8.2g [%8.2g] %9.3g  %8.3g  %8.3g  %5.2f %5.2f+-%4.2f  %-12s %-26s", pk.number,
                p.mu, pk.dmuDay, p.dmu, p.area, pk.dareaDay, p.darea, p.areaFree, pk.zDay, pk.sb, p.sigma, p.sFree, p.dsFree,
                KindName(pk.kind), pk.text.c_str());
    if (pk.eLine > 0.) txt << Form("  line %8.2f (%+6.2f)", pk.eLine, p.mu - pk.eLine);
    txt << "  | " << sh;
    if (eff > 0.) txt << Form("  | %.3g", eff);
    txt << "\n";
    const auto Share = [&](size_t k, bool name) { return k < pk.shares.size() ? (name ? parts[pk.shares[k].first].name : Form("%.4f", pk.shares[k].second)) : std::string(); };
    csv << Form("%d,%.3f,%.3f,%.3f,%.6g,%.4g,%.4g,%.6g,%.4g,%.4g,%.3f,%.3f,%.3f,%s,\"%s\",\"%s\",%.2f,%.3f,%s,%s,%s,%s,%.4g,%.1f,%d\n",
                pk.number, p.mu, pk.dmuDay, p.dmu, p.area, pk.dareaDay, p.darea, p.areaFree, pk.zDay, pk.sb, p.sigma, p.sFree, p.dsFree,
                KindName(pk.kind), pk.text.c_str(), pk.source.c_str(), pk.eLine, pk.eLine > 0. ? p.mu - pk.eLine : 0.,
                Share(0, false).c_str(), Share(0, true).c_str(), Share(1, false).c_str(), Share(1, true).c_str(), eff, sec.chi2, sec.ndf);
  }
  txt << "\nsections (fit windows; chi2 with the simulation's errors):\n";
  for (const Section& sec : sections) {
    int n = 0;
    for (const FitPeak& p : sec.peaks) if (p.free && p.mu >= sec.c0 && p.mu < sec.c1) ++n;
    txt << Form("  %7.1f-%7.1f keV (window %7.1f-%7.1f): %2d peaks, chi2/ndf %.1f/%d%s\n", sec.c0, sec.c1, sec.w0, sec.w1, n, sec.chi2,
                sec.ndf, sec.status > 1 ? Form(" (Minuit2 status %d)", sec.status) : "");
  }
  txt << Form("\nlines that should be seen (grey in the linear plots): emitted (truth; within half a resolution merged) x the "
              "full-energy efficiency,\n  expected significance >= %g in the simulation's statistics (as the search). The efficiency: "
              "Eff x a factor from the fitted main lines (keV: factor)\n ", zMin);
  for (const auto& [e, k] : calibration) txt << Form(" %.0f: %.2f", e, k);
  txt << "\n";
  for (const Expected& x : expected) {
    const ElementInfo* el = Find(x.z);
    txt << Form("  %8.2f keV  %-2s %-28s emitted %9.3g/day  detected %9.3g/day  z %7.1f  %s\n", x.e, el->name,
                LineName(*el, x.e, false).c_str(), x.perDay, x.detected, x.zSim,
                x.peak ? Form("peak #%d%s", x.peak, x.carried ? "" : " (named otherwise)") : "no peak within 1.5 sigma");
  }
  std::ofstream(out + "/peaks_prompt.txt") << txt.str();
  std::ofstream(out + "/peaks_prompt.csv") << csv.str();
  printf("%s", txt.str().c_str());

  // the plots
  // the sections' weights at x: 1 inside a core, cross-faded over +-1 resolution at each boundary between two sections
  auto Weights = [&](double x) {
    std::vector<std::pair<const Section*, double>> ws;
    double sum = 0.;
    for (size_t k = 0; k < sections.size(); ++k) {
      const Section& sec = sections[k];
      const double d0 = Res(sec.c0), d1 = Res(sec.c1);
      const double rise = k == 0 ? 1. : std::clamp((x - (sec.c0 - d0)) / (2. * d0), 0., 1.);
      const double fall = k + 1 == sections.size() ? (x < sec.c1 ? 1. : 0.) : std::clamp((sec.c1 + d1 - x) / (2. * d1), 0., 1.);
      if (k == 0 && x < sec.c0) continue;
      if (rise * fall > 0.) { ws.push_back({&sec, rise * fall}); sum += rise * fall; }
    }
    if (sum > 0.)
      for (auto& w : ws) w.second /= sum;
    return ws;
  };
  auto ModelAt = [&](double x, double bw, bool withPeaks) {   // the sections' fit at x: the entries of a bin bw keV wide
    double v = 0.;
    for (const auto& [sec, w] : Weights(x)) v += w * sec->Model(x, bw, withPeaks);
    return v;
  };
  auto Rebinned = [&](double bw, const std::string& name) {   // the spectrum in bins bw keV wide (the input's merged)
    auto h = dynamic_cast<TH1D*>(fine->Clone(name.c_str()));
    h->SetDirectory(nullptr);
    const int k = std::max(1, int(std::lround(bw / gFineKeV)));
    if (k > 1) h->Rebin(k);
    return h;
  };
  // a section's lines: the line peaks it reported itself; everything else it holds (its structures, the neighbours'
  // peaks in its margins, peaks found twice) is its baseline
  std::map<int, std::pair<int, int>> reportOf;   // candidate -> (kind, the reporting section)
  for (const Peak& pk : peaks) reportOf[pk.fit.candidate] = {pk.kind, pk.section};
  auto IsLine = [&](const Section* sec, const FitPeak& p) {
    const auto it = reportOf.find(p.candidate);
    return p.free && it != reportOf.end() && it->second.second == int(sec - sections.data()) && it->second.first != kindStructure;
  };
  auto Part = [&](double x, double bw, bool lines) {   // the sections' (cross-faded) lines, or their baseline: a bin bw keV wide at x
    double v = 0.;
    for (const auto& [sec, w] : Weights(x)) {
      double u = lines ? 0. : sec->Model(x, bw, false);
      for (const FitPeak& p : sec->peaks)
        if (IsLine(sec, p) == lines) u += p.area * GaussIn(x - 0.5 * bw, x + 0.5 * bw, p.mu, p.sigma);
      v += w * u;
    }
    return v;
  };
  auto Baseline = [&](double x, double bw) { return Part(x, bw, false); };
  auto Lines = [&](double x, double bw) { return Part(x, bw, true); };
  // the fit as curves (the sections cross-faded), leaving out points <= 0 (log scales): with its peaks, or its baseline;
  // the entries of a bin bw keV wide, or (bw 0) per keV in the fit's bins
  auto DrawCurves = [&](double lo, double hi, double bw, bool withPeaks, int style, int width) {
    TGraph* first = nullptr;
    TGraph* g = nullptr;
    for (double x = lo; x <= hi; x += 0.25 * std::min(bw > 0. ? bw : BinWidth(x), Res(x))) {
      const double w = bw > 0. ? bw : BinWidth(x);
      const double v = (withPeaks ? ModelAt(x, w, true) : Baseline(x, w)) / (bw > 0. ? 1. : w);   // without the peaks: the baseline
      if (v <= 0.) { if (g) g->Draw("l same"); g = nullptr; continue; }
      if (!g) {
        g = new TGraph();
        g->SetLineColor(kRed + 1); g->SetLineStyle(style); g->SetLineWidth(width);
        if (!first) first = g;
      }
      g->AddPoint(x, v);
    }
    if (g) g->Draw("l same");
    return first;
  };
  const int ex = int(std::floor(std::log10(rate)));
  const std::string beamText = Form("%.2g#times10^{%d} #mu^{-}/s", rate / std::pow(10., ex), ex);
  double xMax = eMin;
  for (const Peak& pk : peaks) xMax = std::max(xMax, pk.fit.mu);
  xMax = std::min(eMax, 1.15 * xMax);
  {   // the whole spectrum, every peak numbered
    auto c = new TCanvas("c_overview", "", 1800, 800);
    c->SetLogy(); c->SetLogx();
    c->SetLeftMargin(0.07); c->SetRightMargin(0.02); c->SetTopMargin(0.08); c->SetBottomMargin(0.12);
    auto t = dynamic_cast<TH1D*>(total->Clone("overview_total"));
    t->SetDirectory(nullptr);
    t->Scale(1., "width");   // per keV: the fit's bins differ in width
    t->GetXaxis()->SetRangeUser(eMin, xMax);
    t->SetLineColor(kBlack); t->SetLineWidth(1);
    t->SetTitle(Form("Tagged-beam spectrum (hits within %g ns of their muon), one day at %s, 1 cm steel, CZT 10 x 10 cm: the %zu "
                     "fitted peaks;energy in one pixel (keV);entries per keV in one day", promptNs, beamText.c_str(), peaks.size()));
    t->GetXaxis()->SetNoExponent();
    t->GetXaxis()->SetTitleOffset(1.2); t->GetYaxis()->SetTitleOffset(0.85);
    double yLow = 1e300;
    for (int b = t->FindBin(eMin); b <= t->FindBin(xMax); ++b) if (t->GetBinContent(b) > 0.) yLow = std::min(yLow, t->GetBinContent(b));
    t->SetMaximum(t->GetMaximum() * 12.);
    t->SetMinimum(0.3 * yLow);
    t->Draw("hist");
    DrawCurves(eMin, xMax, 0., true, 1, 2);
    t->Draw("hist same");
    // numbers: lines red, structures grey; close neighbours (under 12 pixels apart) one, two or three steps higher
    TLatex tex;
    tex.SetTextSize(0.016); tex.SetTextAngle(90); tex.SetTextAlign(12);
    const double frameW = 1800. * (1. - 0.07 - 0.02);
    double pxLast = -1e9;
    int level = 0;
    for (const Peak& pk : peaks) {
      const double px = std::log(pk.fit.mu / eMin) / std::log(xMax / eMin) * frameW;
      level = px - pxLast < 12. ? (level + 1) % 3 : 0;
      pxLast = px;
      tex.SetTextColor(pk.kind == kindStructure ? kGray + 1 : kRed + 1);
      const double w = BinWidth(pk.fit.mu);
      tex.DrawLatex(pk.fit.mu, ModelAt(pk.fit.mu, w, true) / w * 1.3 * std::pow(2.2, level), Form("%d", pk.number));
    }
    // the energy labels, drawn here (ROOT's own crowd 80, 90, 100)
    t->GetXaxis()->SetLabelSize(0.);
    TLatex xl;
    xl.SetTextSize(0.04); xl.SetTextAlign(23); xl.SetTextFont(42);
    for (double e : {20., 30., 50., 100., 200., 300., 500., 1000., 2000., 3000., 5000.})
      if (e >= eMin && e <= xMax)
        xl.DrawLatexNDC(0.07 + std::log(e / eMin) / std::log(xMax / eMin) * (1. - 0.07 - 0.02), 0.11, Form("%g", e));
    TLatex note;
    note.SetTextSize(0.03); note.SetTextAlign(33);
    note.DrawLatexNDC(0.975, 0.9, "black: simulated   red: fit; peak numbers (grey: a structure, no line) as in peaks_prompt.txt and _panels");
    note.SetTextAlign(11);
    note.DrawLatexNDC(0.08, 0.15, Form("bins: %g keV up to 100 keV, %g keV up to 400 keV, %g keV above", BinWidth(50.), BinWidth(200.), BinWidth(1e4)));
    c->Print((out + "/peaks_prompt_overview.png").c_str());
    c->Print((out + "/peaks_prompt_overview.pdf").c_str());
  }
  {   // the regions in panels, each peak named
    const std::vector<std::pair<double, double>> regions = {{20., 100.}, {100., 200.}, {200., 450.}, {450., 1000.}, {1000., 1800.}, {1800., 3600.}};
    auto c = new TCanvas("c_panels", "", 2400, 2500);
    auto head = new TPad("head", "", 0., 0.952, 1., 1.);
    auto body = new TPad("body", "", 0., 0., 1., 0.952);
    head->Draw(); body->Draw();
    head->cd();
    TLatex title;
    title.SetTextSize(0.28); title.SetTextAlign(12); title.SetTextFont(42);
    title.DrawLatexNDC(0.01, 0.7, Form("Tagged beam (hits within %g ns of their muon), one day at %s, 1 cm steel, CZT 10 x 10 cm: "
                                       "%zu fitted peaks (number, name, keV)", promptNs, beamText.c_str(), peaks.size()));
    TLatex key;
    key.SetTextSize(0.22); key.SetTextAlign(12); key.SetTextFont(42);
    key.DrawLatexNDC(0.01, 0.25, "red: a line (esc.: Cd/Te escape; SE/DE: single/double escape; +X: with element X; *: also two X-rays in one pixel); grey: a structure; dashed: background + structures");
    auto lData = new TGraph(), lFit = new TGraph(), lBg = new TGraph();
    lData->SetLineColor(kBlack);
    lFit->SetLineColor(kRed + 1); lFit->SetLineWidth(2);
    lBg->SetLineColor(kRed + 1); lBg->SetLineStyle(2);
    auto leg = new TLegend(0.82, 0.04, 0.995, 0.96);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.22); leg->SetTextFont(42);
    leg->AddEntry(lData, "simulated", "l");
    leg->AddEntry(lFit, "fit (widths: resolution)", "l");
    leg->AddEntry(lBg, "its baseline", "l");
    leg->Draw();
    body->cd();
    body->Divide(2, 3, 0.001, 0.001);
    for (size_t r = 0; r < regions.size(); ++r) {
      const double lo = regions[r].first, hi = std::min(regions[r].second, eMax);
      body->cd(r + 1);
      gPad->SetLogy();
      const double lm = 0.09, rm = 0.02, tm = 0.07, bm = 0.11;
      gPad->SetLeftMargin(lm); gPad->SetRightMargin(rm); gPad->SetTopMargin(tm); gPad->SetBottomMargin(bm);
      const double pw = BinWidth(lo);   // the panel's bins
      auto t = Rebinned(pw, Form("panel_total_%zu", r));
      t->GetXaxis()->SetRangeUser(lo, hi);
      t->SetLineColor(kBlack); t->SetLineWidth(1);
      // the data in the lower 52% of the frame (log), the labels in a row above them
      double yMin = 1e300, yData = 0.;
      for (int b = t->FindBin(lo); b <= t->FindBin(hi); ++b) {
        const double y = t->GetBinContent(b);
        if (y > 0.) { yMin = std::min(yMin, y); yData = std::max(yData, y); }
      }
      yMin *= 0.5;
      const double range = std::log10(2. * yData / yMin) / 0.52;
      t->SetMinimum(yMin);
      t->SetMaximum(yMin * std::pow(10., range));
      t->SetTitle(Form("%g-%g keV;energy in one pixel (keV);entries per %g keV in one day", lo, hi, pw));
      t->GetXaxis()->SetTitleSize(0.045); t->GetXaxis()->SetLabelSize(0.04);
      t->GetYaxis()->SetTitleSize(0.045); t->GetYaxis()->SetLabelSize(0.04); t->GetYaxis()->SetTitleOffset(1.0);
      t->Draw("hist");
      DrawCurves(lo, hi, pw, false, 2, 1);
      DrawCurves(lo, hi, pw, true, 1, 2);
      t->Draw("hist same");
      // labels: rotated, at least 1.2 text heights apart, a red line down to each peak
      const double textSize = 0.028;
      const double frameW = gPad->GetWw() * gPad->GetAbsWNDC() * (1. - lm - rm), padH = gPad->GetWh() * gPad->GetAbsHNDC();
      const double dxMin = 1.2 * textSize * padH / frameW * (hi - lo);
      const double yLab = yMin * std::pow(10., 0.545 * range);
      double xLast = lo - 0.2 * dxMin;   // the first label most of a text height inside the frame
      TLatex tex;
      tex.SetTextSize(textSize); tex.SetTextAngle(90); tex.SetTextAlign(12);
      TLine l;
      for (const Peak& pk : peaks) {
        if (pk.fit.mu < lo || pk.fit.mu > hi) continue;
        const int color = pk.kind == kindStructure ? kGray + 1 : kRed + 1;
        const double yTop = ModelAt(pk.fit.mu, pw, true), xLab = std::max(pk.fit.mu, xLast + dxMin);
        xLast = xLab;
        l.SetLineColor(color);
        tex.SetTextColor(color);
        if (yTop > 0.) l.DrawLine(pk.fit.mu, yTop * 1.12, xLab, yLab * 0.9);
        tex.DrawLatex(xLab, yLab, Form("%d %s %.1f", pk.number, pk.latex.c_str(), pk.fit.mu));
      }
    }
    c->cd();
    c->Print((out + "/peaks_prompt_panels.png").c_str());
    c->Print((out + "/peaks_prompt_panels.pdf").c_str());
  }
  // the peaks over the baseline, linear, one energy range per plot (the user, 2026-11-11). The baseline is the fit's
  // continuum: each section's background and its structures (Compton edges, Compton continuum, backscatter; grey above).
  // What is left are the lines (with the escapes, sums and the 511 keV line), drawn with the fit: all of them and each one.
  // ranges of peaks of comparable heights (Fe K_alpha on its own: 30 times the peaks beside it)
  const std::vector<std::pair<double, double>> ranges = {{20., 75.},     {75., 125.},    {125., 200.},    {200., 290.},
                                                         {290., 460.},   {460., 800.},   {800., 1220.},   {1220., 1290.},
                                                         {1290., 1480.}, {1480., 1900.}, {1900., 3000.}};
  std::string linearFiles;
  for (const auto& [lo, hi] : ranges) {
    const std::string tag = Form("peaks_prompt_linear_%g_%g", lo, hi);
    linearFiles += (linearFiles.empty() ? "" : ", ") + tag;
    auto c = new TCanvas(("c_" + tag).c_str(), "", 1600, 900);
    const double lm = 0.09, rm = 0.02, tm = 0.17, bm = 0.11;
    c->SetLeftMargin(lm); c->SetRightMargin(rm); c->SetTopMargin(tm); c->SetBottomMargin(bm);
    const double bw = BinWidth(lo);   // the plot's bins
    auto r = Rebinned(bw, "r_" + tag);   // minus the baseline, the errors kept
    for (int b = r->FindBin(lo); b <= r->FindBin(hi); ++b) r->SetBinContent(b, r->GetBinContent(b) - Baseline(r->GetBinCenter(b), bw));
    r->GetXaxis()->SetRangeUser(lo, hi);
    double top = 0., bottom = 0.;
    for (int b = r->FindBin(lo); b <= r->FindBin(hi); ++b) {
      top = std::max(top, r->GetBinContent(b));
      bottom = std::min(bottom, r->GetBinContent(b) - r->GetBinError(b));
    }
    for (double x = lo; x <= hi; x += 0.25 * std::min(bw, Res(x))) top = std::max(top, Lines(x, bw));   // the curve's tops too
    const double yMax = top / 0.52;   // the data in the lower half, the labels above
    r->SetMinimum(std::min(-0.06 * yMax, 1.15 * bottom));
    r->SetMaximum(yMax);
    r->SetLineColor(kBlack); r->SetLineWidth(1); r->SetMarkerSize(0);
    r->SetTitle(Form("Tagged beam, one day at %s, 1 cm steel, CZT 10 x 10 cm: the peaks over the baseline, %g-%g keV;energy in one "
                     "pixel (keV);entries per %g keV per day, minus the baseline", beamText.c_str(), lo, hi, bw));
    r->GetXaxis()->SetTitleSize(0.045); r->GetXaxis()->SetLabelSize(0.04);
    r->GetYaxis()->SetTitleSize(0.042); r->GetYaxis()->SetLabelSize(0.04); r->GetYaxis()->SetTitleOffset(1.0);
    r->Draw("hist");
    r->Draw("e same");
    TLine zero;
    zero.SetLineColor(kGray + 1); zero.SetLineStyle(3);
    zero.DrawLine(lo, 0., hi, 0.);
    TLine seen;   // the lines that should be seen (truth): grey, dashed, under the fit
    seen.SetLineColor(kGray + 1); seen.SetLineStyle(2);
    for (const Expected& x : expected) if (x.e >= lo && x.e <= hi) seen.DrawLine(x.e, r->GetMinimum(), x.e, 0.53 * yMax);
    TGraph* gEach = nullptr;
    for (const Peak& pk : peaks) {
      const FitPeak& p = pk.fit;
      if (pk.kind == kindStructure || p.mu < lo - 3. * p.sigma || p.mu > hi + 3. * p.sigma) continue;
      auto g = new TGraph();
      for (double x = std::max(lo, p.mu - 4. * p.sigma); x <= std::min(hi, p.mu + 4. * p.sigma); x += 0.25 * std::min(bw, p.sigma))
        g->AddPoint(x, p.area * GaussIn(x - 0.5 * bw, x + 0.5 * bw, p.mu, p.sigma));
      g->SetLineColor(kRed - 7); g->SetLineWidth(1);
      if (g->GetN() > 1) g->Draw("l same");
      if (!gEach) gEach = g;
    }
    auto gAll = new TGraph();
    for (double x = lo; x <= hi; x += 0.25 * std::min(bw, Res(x))) gAll->AddPoint(x, Lines(x, bw));
    gAll->SetLineColor(kRed + 1); gAll->SetLineWidth(2);
    gAll->Draw("l same");
    r->Draw("hist same");
    // labels: rotated, in one row above the data, at least 1.2 text heights apart, a line down to each: red the fit's lines,
    // grey a line that should be seen where no red label has its name
    struct Label { double x, y; std::string text; int color; };
    std::vector<Label> labels;
    for (const Peak& pk : peaks)
      if (pk.kind != kindStructure && pk.fit.mu >= lo && pk.fit.mu <= hi)
        labels.push_back({pk.fit.mu, std::max(Lines(pk.fit.mu, bw), r->GetBinContent(r->FindBin(pk.fit.mu))) + 0.02 * yMax,
                          Form("%d %s %.1f", pk.number, pk.latex.c_str(), pk.fit.mu), kRed + 1});
    for (const Expected& x : expected)
      if (!x.carried && x.e >= lo && x.e <= hi)
        labels.push_back({x.e, 0.53 * yMax, Form("%s %s %.1f", Find(x.z)->name, LineName(*Find(x.z), x.e, true).c_str(), x.e), kGray + 1});
    std::stable_sort(labels.begin(), labels.end(), [](const Label& a, const Label& b) { return a.x < b.x; });
    const double textSize = 0.028, frameW = 1600. * (1. - lm - rm), padH = 900.;
    const double dxMin = 1.2 * textSize * padH / frameW * (hi - lo), yLab = 0.56 * yMax;
    double xLast = lo - 0.2 * dxMin;
    TLatex tex;
    tex.SetTextSize(textSize); tex.SetTextAngle(90); tex.SetTextAlign(12);
    TLine l;
    for (const Label& lb : labels) {
      const double xLab = std::max(lb.x, xLast + dxMin);
      xLast = xLab;
      l.SetLineColor(lb.color);
      tex.SetTextColor(lb.color);
      l.DrawLine(lb.x, lb.y, xLab, yLab - 0.015 * yMax);
      tex.DrawLatex(xLab, yLab, lb.text.c_str());
    }
    TLatex note;
    note.SetTextSize(0.026); note.SetTextAlign(12); note.SetTextFont(42);
    note.DrawLatexNDC(lm, 0.915, "baseline: the fit's background and its structures (Compton edges, Compton continuum, backscatter)");
    auto leg = new TLegend(0.40, 0.85, 0.98, 0.89);
    leg->SetNColumns(4); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.026); leg->SetTextFont(42);
    leg->AddEntry(r, "simulated", "le");
    leg->AddEntry(gAll, "fit: the lines", "l");
    if (gEach) leg->AddEntry(gEach, "each line", "l");
    auto gSeen = new TGraph();
    gSeen->SetLineColor(kGray + 1); gSeen->SetLineStyle(2);
    leg->AddEntry(gSeen, "should be seen (truth)", "l");
    leg->Draw();
    c->Print((out + "/" + tag + ".png").c_str());
    c->Print((out + "/" + tag + ".pdf").c_str());
  }
  printf("written into %s: peaks_prompt.txt, .csv, peaks_prompt_overview, peaks_prompt_panels, %s (.png, .pdf)\n", out.c_str(),
         linearFiles.c_str());
}
