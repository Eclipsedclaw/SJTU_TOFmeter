// One thicker Al target (default 4 cm), its simulated data split into the target's muonic X-rays (signal, red) and
// everything else (background, blue): which measured variables tell them apart, the most efficient cuts on them, and
// whether a measurement could still find the Al lines after the cuts. From the skims (analysis/muonic_xray_stack_skim.C)
// of the target's run of scripts/muonic_xray_target.sh and of the 3-week run:
//   root -l -b -q 'analysis/muonic_xray_target_select.C'
// The data: every stack entry of 20-1000 keV with a middle-TOF hit in coincidence (-5..+30 ns), from the target's run
// (cosmic muons up to eSplit) and from the 3-week run's faster muons, each scaled to the measurement shown, measureDays
// (by default since 2026-11-02 the target's run's own days: its entries count as simulated, the faster muons' are
// scaled to those days; the 3-week run's slower muons belong to its 1 cm target). Signal: an entry with a muonic X-ray of the target (stk_src bit 1); the
// background holds the muonic X-rays of the other volumes (Al plates, carbon of the TOF bars), which make the same or
// nearby lines. The TOF response of muonic_xray_stack_cuts.C (sigmaT per station, sigmaE on the bar energies, its
// seeds; the top-TOF energy smeared alike).
// 1. Corner plots of the candidate variables (for the CZT run of single pixels with the measured TOF one plot since
//    2026-11-02, corner_after: the muon's and the entry's variables together, each range a cut rejects shaded grey),
//    the muon's (bottom-TOF hit, middle-TOF bar, middle-TOF energy,
//    top-to-middle time, top-TOF energy) and the stack entry's (highest and lowest layer hit, hits in the entry, largest
//    hit / entry, outermost hit), signal as red points, background as a blue density. Since 2026-11-02 on the captured
//    data at every energy, as a measurement records it (the TOF trigger, an entry in time; one point per event,
//    signal an event with a muonic X-ray of the target), with the cuts chosen in the L_alpha window drawn. Every
//    event inside the axes (asked for on 2026-11-02): with the measured TOF the bar energies and the kinetic energy on
//    log axes, the kinetic energy's last bin the events with beta >= 1 (no kinetic energy); with single pixels the
//    other pixels' energy on a log axis, its first bin the entries without any; since a later request the same day finer:
//    the other pixels' energy zoomed in on 0.01-30 MeV at 50 bins a decade, its last bin the entries of 30 MeV or more, its
//    diagonal on a log scale, and the entry's energy on a log axis at 2% wide bins (the CZT resolves the lines there).
// 2. The most efficient cuts: greedy, one cut at a time the variable and threshold that raise S/sqrt(B) in the L_alpha
//    window the most (by at least 1%), then each threshold again with the others fixed; chosen on the even events and
//    scored on the odd ones and back, to see the optimisation's own gain. The CZT run of single pixels with the measured
//    TOF uses its own cuts instead (2026-11-02: chosen by hand, then optimised for the L_alpha peak with the work folder's
//    output/muonic_xray_czt/pixel_tof/cut_optimization.C), the most efficient ones printed and fitted beside them. The
//    cut flow.
// 3. The spectrum before and after the cuts by what made it (20-500 keV since 2026-11-05, Al's 20-450 keV before;
//    entries per 2 keV in measureDays days, drawn
//    as ROOT histograms: the simulated total with its statistical error bars, each class a line), and fits without truth
//    to the simulated spectrum after the cuts itself, what a measurement of measureDays days sees (the MC truth;
//    since 2026-11-02 in place of simulated measurements drawn from it): a smooth background, exp(pol3) (exp(pol2) with
//    backgroundDegree 2: for a panel as coarse as the LYSO, where the cubic term takes in a peak of its width), and
//    Gaussian peaks: two at free positions, from the bins alone, their widths the detector's resolution at their
//    positions (measurement_fit; since 2026-11-05, the resolution known from calibration; the same with free widths is
//    still fitted beside it), or L_alpha and L_beta at the known line energies (only the amounts free); each with the
//    likelihood-ratio significance of its L_alpha peak. A null test: the same fits to the expected spectrum without any
//    muonic X-ray. Since 2026-11-05 (asked for: the fitting all the way to 500 keV, to identify the K lines) also a
//    K-line fit at 160-500 keV by the same method, with three peaks, drawn in spectrum_after beside the L-line fit's
//    30-160 keV. The cut sets compared, each with the known-energy fit (with the stack's
//    skims also those of the 1 cm target and the earlier own cuts of muonic_xray_target_cuts.C).
// Writes <out>/muonic_xray_select_<plot>.png/.pdf (one plot per file) and muonic_xray_select.csv.
// The CZT plane under the 4 cm target (scripts/muonic_xray_czt.sh, its bands skimmed with detector "czt"): detector
// "czt" and the CZT's resolution, the slow band's skim as the target's run, the fast band's as the faster muons; the
// entry's variables are then the pixel hits, the largest hit / entry and whether it lies in the 4 central pixels:
//   root -l -b -q 'analysis/muonic_xray_target_select.C("output/muonic_xray_czt/pixel_tof/full_run/xray_skim_slow.root", "output/muonic_xray_czt/pixel_tof/full_run/xray_skim_fast.root", "output/muonic_xray_czt/pixel_tof/full_run", "4 cm target, CZT", "czt", 6.36991)'
// (as the work folder's output/muonic_xray_czt/pixel_tof/run_analysis.sh runs it, after the skims of both bands).
// Skims with the measured TOF (branch trig, analysis/muonic_xray_stack_skim.C since 2026-10-30) are used as measured,
// sigmaT and sigmaE are then not used: the dataset is what the TOF trigger takes (top and middle station hit), a fixed
// step after the middle-TOF coincidence in the cut flow, so that it shows what the trigger costs; the muon's variables 1
// and 3 are then the reconstructed track where it crosses the target's centre plane (max(|x|,|y|)) and the kinetic
// energy from the time of flight; the middle-TOF coincidence and the bottom-TOF veto use the measured times. CZT skims
// of single pixels (branch npix): an entry is one pixel, its variables the pixels hit in the event within 100 ns (5),
// their energy without it (6) and the pixel's place (9). Every condition, of the skims (their TNamed "conditions") and
// of this study, is printed first, written to muonic_xray_select_conditions.txt and as "condition" rows of the CSV, and
// summed up in two lines at the foot of every plot. Also with the measured TOF: the zenith angle of the TOF track (11,
// cut as "within x deg of the vertical") and, for skims with the bottom station's bar (barBotM), the bar of a bottom-TOF
// hit (10, 0 if none); with single pixels the pixel's number 1-16 (12). 10 and 12 are not cut on: a set of bars or
// pixels is not a range, and their useful parts are variables 0 (no bottom hit) and 9 (central pixels). The muon's corner
// plots leave out the bottom-TOF hit and bar (0, 10; asked for on 2026-11-02, as their densities tell little): 0 stays a
// cut the optimisation can choose, and a chosen one is named under the plot as not drawn (in the 4 cm CZT run the veto,
// after the other fixed cuts, still took B from 3.3 to 2.3 per day at no signal cost).

#include <TError.h>
#include <TLatex.h>
#include <TVirtualPad.h>
#include <TParameter.h>

#include "CutProfile.h"   // the analysis cuts of a detector configuration (2026-11-05)

#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>

namespace {
// the detector's energy sigma at 1 MeV (keV; the stack's 40, scaling as sqrt(E)) and which detector, set by the macro
double gSigma1MeV = 40.;
bool gCZT = false;
bool gTOF = false;     // skims with the measured TOF: the trigger, variables 1 and 3 the track and kinetic energy
bool gPixel = false;   // CZT skims of single pixels: variables 5 and 6 the pixels hit and their energy
bool gBottomBar = false;   // skims with the bottom-TOF bar as measured (barBotM): variable 10
bool gOtherPanel = false;  // single-pixel skims of a panel other than the 16-pixel CZT plane (2026-11-04: the LYSO
                           // panel), known by the skim's pixelPitch
double SigmaKeV(double eKeV) { return gSigma1MeV * std::sqrt(eKeV / 1000.); }
// a title at (x, y) in the current pad's NDC, its text size lowered (by 4% steps, to half at most) until it ends left of
// xMax (2026-11-05: the LYSO panel's label, "LYSO 12.69 x 12.69 cm", ran past the corner plot's edge)
void DrawFitted(TLatex& t, double x, double y, double xMax, const char* text)
{
  t.SetNDC();
  t.SetTitle(text);
  const double room = (xMax - x) * gPad->GetWw() * gPad->GetAbsWNDC();   // pixels
  for (int i = 0; i < 17; ++i) {
    UInt_t w = 0, h = 0;
    t.GetBoundingBox(w, h);
    if (w <= room) break;
    t.SetTextSize(t.GetTextSize() * 0.96);
  }
  t.DrawLatex(x, y, text);
}
// the target's muonic X-ray lines as simulated (the cascade model's energies, keV): L_alpha, L_beta, K_alpha, K_beta;
// Al's unless the target argument names Si (2026-11-03)
double kAlpha = 65.73, kBeta = 88.73, kKAlpha = 354.9, kKBeta = 420.6;
const char* gTarget = "Al";
bool InAlpha(double eKeV) { return std::abs(eKeV - kAlpha) <= 2. * SigmaKeV(kAlpha); }

// the candidate cut variables, all measured (no truth); the CZT has no layers (variables 5 and 6), and a single pixel
// no hits to add up (7 and 8); 10 and 11, the bottom-TOF bar and the track's zenith angle, need the measured TOF, 12,
// the pixel, single CZT pixels; 13, the entry's energy with the response, is shown in the corner plots (since
// 2026-11-02), not cut on (the L_alpha window defines the signal). 0-4, 10 and 11 are the muon's, the others the entry's
const int kNV = 14;
// a log axis (log) from lo to hi; split > 0: a bin of its own at the axis' end for the values at or above split (or,
// splitLow, below it), such as a placeholder, drawn at the bin's centre and named splitLabel in the plots; over > 0
// (with no split or a low one): a last bin of its own for the values at or above over, named overLabel, so that the
// axis can zoom in on the range below it; logShape: the corner plot's diagonal draws the shapes on a log scale
struct Variable {
  const char* name; int bins; double lo, hi; bool integer;
  bool log = false; double split = 0.; bool splitLow = false; const char* splitLabel = "";
  double over = 0.; const char* overLabel = ""; bool logShape = false;
};
bool Used(int i)
{
  if (i == 10) return gBottomBar;
  if (i == 11) return gTOF;
  if (i == 12) return gCZT && gPixel;
  if (i == 13) return true;
  return !gCZT || (gPixel ? i != 7 && i != 8 : i != 5 && i != 6);
}
bool IsMuon(int i) { return i < 5 || i == 10 || i == 11; }
Variable kVars[kNV] = {{"bottom-TOF hit", 2, -0.5, 1.5, true},
                             {"middle-TOF bar", 9, -0.5, 8.5, true},
                             {"middle-TOF energy (MeV)", 40, 0., 10., false},
                             {"top-to-middle time (ns)", 40, 0., 4., false},
                             {"top-TOF energy (MeV)", 40, 0., 10., false},
                             {"highest stack layer", 5, -0.5, 4.5, true},
                             {"lowest stack layer", 5, -0.5, 4.5, true},
                             {"hits in the entry", 5, 0.5, 5.5, true},
                             {"largest hit / entry", 26, 0., 1.04, false},
                             {"outermost hit (cm)", 22, 0., 11., false},
                             {"bottom-TOF bar (0 = none)", 10, -0.5, 9.5, true},
                             {"track zenith angle (deg)", 30, 0., 60., false},
                             {"CZT pixel (1-16)", 16, 0.5, 16.5, true},
                             // on a log axis, 2% wide bins (1.3 keV at 66 keV, near the CZT's sigma of 1.6 keV there;
                             // asked for on 2026-11-02), up to 1300 keV so that the label "1000" fits in the pad
                             {"entry energy (keV)", 210, 20., 1300., false, true}};

// a variable's bin edges: even, or even in log on a log axis, the split bin at its end
std::vector<double> Edges(const Variable& v)
{
  std::vector<double> e;
  if (!v.log) {
    for (int k = 0; k <= v.bins; ++k) e.push_back(v.lo + (v.hi - v.lo) * k / v.bins);
    return e;
  }
  const bool split = v.split > 0., over = v.over > 0.;
  const double a = split && v.splitLow ? v.split : v.lo;   // the part even in log: a to b
  const double b = over ? v.over : split && !v.splitLow ? v.split : v.hi;
  const int n = v.bins - (split ? 1 : 0) - (over ? 1 : 0);
  if (split && v.splitLow) e.push_back(v.lo);
  for (int k = 0; k <= n; ++k) e.push_back(a * std::pow(b / a, double(k) / n));
  if ((split && !v.splitLow) || over) e.push_back(v.hi);
  return e;
}
// where a value is drawn: beyond the split (or at or above over) at that bin's centre, else within the axis' range
double ShownAt(double x, const Variable& v)
{
  const bool split = v.split > 0.;
  if (split && (v.splitLow ? x < v.split : x >= v.split)) return v.splitLow ? std::sqrt(v.lo * v.split) : std::sqrt(v.split * v.hi);
  if (v.over > 0. && x >= v.over) return v.log ? std::sqrt(v.over * v.hi) : 0.5 * (v.over + v.hi);
  const double lo = split && v.splitLow ? v.split : v.lo, hi = v.over > 0. ? v.over : split && !v.splitLow ? v.split : v.hi;
  return v.log ? std::min(std::max(x, lo * (1. + 1e-9)), hi * (1. - 1e-9)) : std::min(std::max(x, lo + 1e-6), hi - 1e-6);
}

struct Entry {
  double w;        // the events of the measurement shown that it stands for: measureDays / its sample's days
  int kind, half;  // 0 muonic X-ray of the target (signal), 1 muonic X-ray of another volume, 2 the rest; event parity
  int sample;      // 0 the target's run (muons up to eSplit), 1 the 3-week run's faster muons
  bool trig;       // the TOF trigger fired (always true for skims without the measured TOF)
  bool timed;      // the top TOF timed the muon: variables 3 and 4 exist (and with the measured TOF the track, 1)
  double e;        // add-back energy with the response (keV)
  double v[kNV];
  long long id;    // the event: (sample x 1000 + process) x 10^8 + event number
  int passed = -1; // one point per event (corner plots): bit i, one of its entries passes cut i by itself; -1 an entry
  bool lu = false;  // energy of the pixel panel's own radioactivity (src bit 16, /tof/intrinsic/: the LYSO's Lu-176), no
                    // muonic X-ray: a decay alone or piled up on a muon's hit; part of the rest
};

// a cut keeps lo <= v <= hi; an inactive cut keeps everything, and so does a cut on a variable the entry lacks
struct Cut { bool active = false; double lo = -1e9, hi = 1e9; };
using Cuts = std::array<Cut, kNV>;

bool Defined(const Entry& x, int i) { return (i != 3 && i != 4 && i != 11 && (i != 1 || !gTOF)) || x.timed; }
bool Passes(const Entry& x, const Cuts& cuts, int i)
{
  const Cut& c = cuts[i];
  return !c.active || !Defined(x, i) || (x.v[i] >= c.lo && x.v[i] <= c.hi);
}
bool PassesAll(const Entry& x, const Cuts& cuts)
{
  for (int i = 0; i < kNV; ++i) if (!Passes(x, cuts, i)) return false;
  return true;
}

// the thresholds the optimisation tries for variable i
std::vector<Cut> Options(int i)
{
  std::vector<Cut> o;
  auto range = [&](double lo, double hi) { Cut c; c.active = true; c.lo = lo; c.hi = hi; o.push_back(c); };
  switch (i) {
    case 0: range(-1e9, 0.5); break;   // no bottom-TOF hit
    case 1:
      if (gTOF) for (double x : {6., 8., 10., 12., 14., 17., 20.}) range(-1e9, x);   // track at the target within x cm
      else                                                                           // middle-TOF bars a-b
        for (double a : {-0.5, 0.5, 1.5, 2.5})
          for (double b : {5.5, 6.5, 7.5, 8.5}) if (a > 0. || b < 8.) range(a, b);
      break;
    case 2: for (int k = 4; k <= 16; ++k) range(k / 4., 1e9); break;     // middle-TOF energy above 1-4 MeV
    case 3:
      if (gTOF)   // kinetic energy from the TOF from lo to hi MeV
        for (double lo : {0., 15., 25.}) for (double hi : {40., 50., 60., 70., 85., 100., 150.}) range(lo, hi);
      else for (int k = 8; k <= 18; ++k) range(k / 10., 1e9);            // top-to-middle time above 0.8-1.8 ns
      break;
    case 4: for (int k = 2; k <= 12; ++k) range(k / 4., 1e9); break;     // top-TOF energy above 0.5-3 MeV
    case 5:
      if (gCZT) for (double k : {1.5, 2.5, 3.5}) range(-1e9, k);         // at most 1-3 pixels hit in the event
      else for (double k : {0.5, 1.5, 2.5, 3.5}) range(k, 1e9);          // layer at least 1-4
      break;
    case 6:
      if (gCZT) for (double x : {0.05, 0.1, 0.2, 0.5, 1., 2.}) range(-1e9, x);   // the other pixels below x MeV
      else for (double k : {0.5, 1.5, 2.5, 3.5}) range(k, 1e9);
      break;
    case 7: for (double k : {1.5, 2.5, 3.5}) range(-1e9, k); break;      // at most 1-3 hits
    case 8: for (double x : {0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99}) range(x, 1e9); break;
    case 9:   // outermost hit within x cm; the CZT's pixel centres are at 1.25 and 3.75 cm: the 4 central pixels only
      if (gCZT && !gOtherPanel) range(-1e9, 2.5);
      else for (double x : {9.5, 9., 8., 7., 6., 5.}) range(-1e9, x);
      break;
    case 11: for (double a : {15., 20., 25., 30., 35., 40.}) range(-1e9, a); break;   // track within a deg of the vertical
    default: break;   // 10, 12 and 13, the bottom-TOF bar, the CZT pixel and the energy: shown, not cut
  }
  return o;
}

std::string Describe(const Cuts& cuts, int i)
{
  const Cut& c = cuts[i];
  if (!c.active) return "none";
  switch (i) {
    case 0: return "no bottom-TOF hit";
    case 1:
      if (gTOF) return Form("track at the target within %.0f cm", c.hi);
      return Form("middle-TOF bar %d-%d", int(std::ceil(c.lo)), int(std::floor(c.hi)));
    case 2: return c.hi < 1e8 ? Form("middle-TOF energy %.1f-%.1f MeV", c.lo, c.hi) : Form("middle-TOF energy > %.2f MeV", c.lo);
    case 3:
      if (gTOF) return c.lo > 0. ? Form("TOF kinetic energy %.0f-%.0f MeV", c.lo, c.hi) : Form("TOF kinetic energy < %.0f MeV", c.hi);
      return Form("top-to-middle time > %.1f ns", c.lo);
    case 4: return c.hi < 1e8 ? Form("top-TOF energy %.1f-%.1f MeV", c.lo, c.hi) : Form("top-TOF energy > %.2f MeV", c.lo);
    case 5:
      if (gCZT) return Form("at most %d pixels hit", int(std::floor(c.hi)));
      return Form("highest layer >= %d", int(std::ceil(c.lo)));
    case 6:
      if (gCZT) return Form("other pixels < %.2f MeV", c.hi);
      return Form("lowest layer >= %d", int(std::ceil(c.lo)));
    case 7: return Form("at most %d hits", int(std::floor(c.hi)));
    case 8: return Form("largest hit > %.2f of the entry", c.lo);
    case 9: return gCZT && !gOtherPanel ? std::string("one of the 4 central pixels")
                                        : std::string(Form("outermost hit < %.1f cm", c.hi));
    case 11: return Form("track within %.0f deg of the vertical", c.hi);
    case 13: return c.lo > 0. ? Form("entry energy %.0f-%.0f keV", c.lo, c.hi) : Form("entry energy < %.0f keV", c.hi);
    default: return "none";
  }
}
}  // namespace

void muonic_xray_target_select(const char* slowSkim = "output/muonic_xray_target/t4cm/xray_skim.root",
                               const char* ref = "output/muonic_xray_stack/xray_skim_pbke.root",
                               const char* outDir = "output/muonic_xray_target/t4cm", const char* label = "4 cm target",
                               const char* detector = "stack", double sigma1MeV = 40., double eSplit = 1000.,
                               double measureDays = 0., int toys = 300, double sigmaT = 0.2, double sigmaE = 0.15,
                               unsigned seed = 4357, const char* target = "Al", const char* cutProfile = "czt",
                               int backgroundDegree = 3)
// target: the target's element, Al or Si (2026-11-03: the same geometry, its own lines and names).
// backgroundDegree (2026-11-05): the fits' background exp(pol3), or exp(pol2) with 2 (its cubic term held at 0).
// cutProfile (2026-11-05): the analysis cuts of the pixel skims with the measured TOF, a profile of analysis/cuts/ by name
// (czt: the CZT plane's, as before; lyso: the LYSO panel's, 20 x 20 pixels, 12.69 x 12.69 cm) or a path
// (analysis/CutProfile.h).
// measureDays: the measurement every number and plot shows (days); 0, the default since 2026-11-02 (asked for: the
// actual simulated counts, no 21-day equivalents), the target's run's own days, the faster muons' entries scaled to them.
// toys and seed: not used since 2026-11-02 (the fits take the simulated spectrum, no simulated measurements); kept so
// that the arguments stay in their places
{
  TH1::AddDirectory(false);
  const TString out = outDir;
  gSigma1MeV = sigma1MeV;
  gTarget = target;
  if (std::string(target) == "Si") { kAlpha = 76.24; kBeta = 102.93; kKAlpha = 411.7; kKBeta = 487.9; }
  gCZT = std::string(detector) == "czt";
  const char* entryName = gCZT ? "CZT entry" : "stack entry";
  const char* energyTitle = "add-back energy with the response (keV)";   // single pixels: the pixel energy (below)
  // the skims: their kind (the measured TOF, single CZT pixels), the same for both, and their conditions
  std::vector<std::string> conditions;   // everything the study depends on, for its outputs
  std::string stamp[2];                  // the same in two lines at the foot of the plots
  double tofFwhmT = -1., tofFwhmE = -1., vEff = -1.;
  for (int s = 0; s < 2; ++s) {
    const char* path = s ? ref : slowSkim;
    TFile file(path);
    auto cand = dynamic_cast<TTree*>(file.Get("cand"));
    if (!cand) {
      printf("no skim in %s (run analysis/muonic_xray_stack_skim.C)\n", path);
      return;
    }
    const bool tof = cand->GetBranch("trig") != nullptr, single = gCZT && cand->GetBranch("npix") != nullptr;
    const bool bottomBar = cand->GetBranch("barBotM") != nullptr;
    if (s == 0) { gTOF = tof; gPixel = single; gBottomBar = bottomBar; }
    else if (tof != gTOF || single != gPixel || bottomBar != gBottomBar) {
      printf("%s and %s are skims of different kinds (measured TOF, single pixels): skim both again\n", slowSkim, ref);
      return;
    }
    conditions.push_back(Form("skim %s:", path));
    if (auto text = dynamic_cast<TNamed*>(file.Get("conditions"))) {
      std::istringstream lines(text->GetTitle());
      for (std::string line; std::getline(lines, line);) conditions.push_back("  " + line);
    } else {
      conditions.push_back("  no conditions recorded (a skim made before 2026-10-30)");
    }
    auto parameter = [&](const char* name) {
      auto p = dynamic_cast<TParameter<double>*>(file.Get(name));
      if (p && !p->IsConst()) {   // a skim merged from parts made with different settings
        conditions.push_back(Form("  WARNING: the parts merged into this skim differ in %s (the first part's %g used)", name,
                                  p->GetVal()));
        printf("WARNING: %s: the merged parts differ in %s\n", path, name);
      }
      return p ? p->GetVal() : -1.;
    };
    if (s == 0) { tofFwhmT = parameter("tofFwhmT"); tofFwhmE = parameter("tofFwhmE"); vEff = parameter("vEff"); }
    const double skimSigma = parameter("sigma1MeV");
    if (skimSigma > 0. && std::abs(skimSigma - sigma1MeV) > 1e-3) {
      conditions.push_back(Form("  WARNING: this skim's energy sigma is %g keV at 1 MeV, the study assumes %g", skimSigma, sigma1MeV));
      printf("WARNING: %s has an energy sigma of %g keV at 1 MeV, this study assumes %g\n", path, skimSigma, sigma1MeV);
    }
  }
  // the pixel panel (2026-11-04): the CZT plane (16 pixels, 2.5 cm pitch) or another, e.g. the LYSO panel's 30 x 30 on
  // 0.667 cm: its pitch from the skim (pixelPitch; none in the skims before it: the CZT), its pixels from the highest one
  double pitch = 2.5;
  int nPixels = 16;
  if (gCZT && gPixel) {
    TFile file(slowSkim);
    if (auto p = dynamic_cast<TParameter<double>*>(file.Get("pixelPitch"))) pitch = p->GetVal();
    if (auto n = dynamic_cast<TParameter<double>*>(file.Get("pixelsPerSide")))   // the skim's, from the run's log
      nPixels = int(n->GetVal() * n->GetVal());
    else if (auto cand = dynamic_cast<TTree*>(file.Get("cand")))   // older skims: from the highest pixel hit
      nPixels = std::max(16, int(cand->GetMaximum("pix")) + 1);
    gOtherPanel = std::abs(pitch - 2.5) > 1e-3;
    if (gOtherPanel) entryName = "pixel entry";
  }
  if (gTOF) {   // the reconstructed track and kinetic energy in place of the middle-TOF bar and the time of flight
    kVars[1] = {"track at the target (cm)", 30, 0., 30., false};   // max(|x|,|y|) where it crosses the centre plane
    // every event inside the axes: the bar energies (a hit above 0.2 MeV) on log axes to 10^2.5 MeV, 12 bins a decade; the
    // kinetic energy to 100 GeV, 8 a decade, its last bin the skim's placeholder 1e5 MeV of the events with beta >= 1
    kVars[2] = {"middle-TOF energy (MeV)", 42, 0.1, std::pow(10., 2.5), false, true};
    kVars[3] = {"TOF kinetic energy (MeV)", 49, 0.1, std::pow(10., 5.5), false, true, 1e5, false, "#beta#geq1"};
    kVars[4] = {"top-TOF energy (MeV)", 42, 0.1, std::pow(10., 2.5), false, true};
  }
  if (gCZT && gPixel) {   // a single pixel: the pixels hit in the event and their energy, where the pixel is
    energyTitle = "pixel energy with the response (keV)";
    kVars[5] = {"pixels hit in the event", 16, 0.5, 16.5, true};
    // on a log axis zoomed in on 0.01-30 MeV (asked for on 2026-11-02), 50 bins a decade (4.7% wide, 3 keV at 66 keV: the
    // muonic X-rays that reach other pixels, below 0.5 MeV, apart from a muon crossing them, ~7 MeV), its shapes on a log
    // scale; its first bin the entries without energy in the other pixels (a pixel counts from 10 keV; below that only the
    // rounding of a sum, |value| < 1e-6 MeV), its last bin those of 30 MeV or more (the fast band's muons reach 2 GeV)
    kVars[6] = {"other pixels' energy (MeV)", 177, std::pow(10., -2.5), 30. * std::sqrt(10.), false, true, 0.01, true, "0",
                30., "#geq30", true};
    kVars[9] = {"pixel centre (cm)", 5, 0., 5., false};   // max(|x|,|y|): 1.25 inner, 3.75 outer
    if (gOtherPanel) {   // more and smaller pixels: their centres out to the panel's half width, their numbers
      static std::string pixelAxis;
      pixelAxis = Form("pixel (1-%d)", nPixels);
      kVars[5] = {"pixels hit in the event", 20, 0.5, 20.5, true};
      kVars[9] = {"pixel centre (cm)", 22, 0., 0.5 * std::sqrt(double(nPixels)) * pitch + 0.5, false};
      kVars[12] = {pixelAxis.c_str(), 60, 0.5, nPixels + 0.5, false};
    }
  } else if (gCZT) {      // the CZT entry: pixel hits, the largest pixel hit, the pixel centres (1.25 or 3.75 cm out)
    kVars[7] = {"pixel hits in the entry", 5, 0.5, 5.5, true};
    kVars[8] = {"largest hit / entry", 26, 0., 1.04, false};
    kVars[9] = {"outermost pixel centre (cm)", 5, 0., 5., false};
  }
  if (gOtherPanel)
    stamp[0] = Form("pixel panel of %d pixels on a %.3g cm pitch, each read out alone; sigma/E %.1f%% at 662 keV, 10 keV threshold",
                    nPixels, pitch, 100. * sigma1MeV * std::sqrt(0.662) / 662.);
  else
  stamp[0] = gCZT ? Form("CZT plane, %s; FWHM/E %.1f%% at 1 MeV, 10 keV threshold", gPixel ? "each of its 16 pixels read out alone"
                                                                                           : "its pixels summed within 100 ns",
                         2.3548 * sigma1MeV / 10.)
                  : Form("Si/YSO stack, add-back within 100 ns; FWHM/E %.1f%% at 1 MeV", 2.3548 * sigma1MeV / 10.);
  stamp[1] = gTOF ? Form("TOF: %.0f ps FWHM per station, bar energy %.0f%% FWHM, both bar ends read out (%.1f cm along a bar); "
                         "trigger: top and middle", 1000. * tofFwhmT, 100. * tofFwhmE, vEff * tofFwhmT / 2.3548)
                  : Form("TOF response of this study: sigma %.2f ns per station, sigma %.0f%% on the bar energies; middle-TOF "
                         "coincidence, no trigger", sigmaT, 100. * sigmaE);
  std::vector<Entry> entries;
  std::vector<Entry> sideLu;   // the panel's decays alone in the trigger's sidebands (2026-11-05; none without /tof/intrinsic/)
  const double kSideScale = 35. / 1930.;   // the window's 35 ns over the sidebands' 980 + 950 ns
  // the entries of a skim: all of them, or (fastOnly) the background of its primaries above eSplit; returns its days
  auto load = [&](const char* path, bool fastOnly) {
    TFile file(path);
    auto cand = dynamic_cast<TTree*>(file.Get("cand"));
    auto info = dynamic_cast<TTree*>(file.Get("info"));
    if (!cand || !info) {
      printf("no skim in %s (run analysis/muonic_xray_stack_skim.C)\n", path);
      return 0.;
    }
    double lt = 0., skimDays = 0.;
    info->SetBranchAddress("live_time", &lt);
    for (Long64_t i = 0; i < info->GetEntries(); ++i) { info->GetEntry(i); skimDays += lt / 86400.; }
    int proc, event, n, layers, topLayer, src, tofTop, tofMid, barMid, trig = 0, npix = 0;
    float e, emax, rmax, t, tTop, tMid, eTop, eMid, botFirst, ekin;
    float tMidM = -1.f, eTopM = 0.f, eMidM = 0.f, botM = -1.f, tkin = 1e5f, xTgt = 1e3f, yTgt = 1e3f, eOther = 0.f;
    float slopeX = 0.f, slopeY = 0.f;
    int barBotM = -1, pix = -1;
    cand->SetBranchAddress("proc", &proc); cand->SetBranchAddress("event", &event);
    cand->SetBranchAddress("e", &e); cand->SetBranchAddress("n", &n); cand->SetBranchAddress("layers", &layers);
    cand->SetBranchAddress("topLayer", &topLayer); cand->SetBranchAddress("emax", &emax); cand->SetBranchAddress("rmax", &rmax);
    cand->SetBranchAddress("t", &t); cand->SetBranchAddress("src", &src); cand->SetBranchAddress("ekin", &ekin);
    cand->SetBranchAddress("tofTop", &tofTop); cand->SetBranchAddress("tofMid", &tofMid);
    cand->SetBranchAddress("tTop", &tTop); cand->SetBranchAddress("tMid", &tMid);
    cand->SetBranchAddress("eTop", &eTop); cand->SetBranchAddress("eMid", &eMid);
    cand->SetBranchAddress("barMid", &barMid); cand->SetBranchAddress("botFirst", &botFirst);
    if (gTOF) {
      cand->SetBranchAddress("trig", &trig); cand->SetBranchAddress("tMidM", &tMidM); cand->SetBranchAddress("botM", &botM);
      cand->SetBranchAddress("eTopM", &eTopM); cand->SetBranchAddress("eMidM", &eMidM); cand->SetBranchAddress("tkin", &tkin);
      cand->SetBranchAddress("xTgt", &xTgt); cand->SetBranchAddress("yTgt", &yTgt);
      cand->SetBranchAddress("slopeX", &slopeX); cand->SetBranchAddress("slopeY", &slopeY);
    }
    if (gBottomBar) cand->SetBranchAddress("barBotM", &barBotM);
    if (gPixel) {
      cand->SetBranchAddress("npix", &npix); cand->SetBranchAddress("eOther", &eOther); cand->SetBranchAddress("pix", &pix);
    }
    TRandom3 random;
    for (Long64_t i = 0; i < cand->GetEntries(); ++i) {
      cand->GetEntry(i);
      const int kind = (src & 1) ? 0 : (src & 2) ? 1 : 2;
      const long long id = (static_cast<long long>(fastOnly ? 1 : 0) * 1000 + proc) * 100000000LL + event;
      int lowest = 0;
      while (lowest < 4 && !(layers >> lowest & 1)) ++lowest;
      if (gTOF) {   // as measured: a middle-TOF hit in coincidence; the trigger and the reconstruction from the skim
        if (tMidM < 0.) continue;
        if (fastOnly && (kind == 0 || ekin <= eSplit)) continue;
        // in time with the middle TOF; or (2026-11-05) a decay of the panel's own radioactivity alone in the sidebands
        // -990..-10 and +40..+990 ns, as panel_intrinsic_check.C: the same decays at 55 times the statistics, weighted by
        // 35 / 1930 (they come at random times, unrelated to the muon)
        const double dt = t - tMidM;
        const bool inTime = dt >= -5. && dt <= 30.;
        const bool side = src == 16 && ((dt >= -990. && dt <= -10.) || (dt >= 40. && dt <= 990.));
        if (!inTime && !side) continue;
        const bool bottom = botM >= 0. && botM - tMidM > -5. && botM - tMidM < 50.;
        (inTime ? entries : sideLu)
            .push_back({(inTime ? 1. : kSideScale) / skimDays, kind, event & 1, fastOnly ? 1 : 0, trig != 0, trig != 0, e,
                        {bottom ? 1. : 0., std::max(std::abs(xTgt), std::abs(yTgt)), eMidM, tkin, trig ? eTopM : -1.,
                         gPixel ? double(npix) : double(topLayer), gPixel ? eOther / 1000. : double(lowest), double(n),
                         e > 0. ? emax / e : 0., rmax, bottom && barBotM >= 0 ? barBotM + 1. : 0.,
                         trig ? std::atan(std::hypot(slopeX, slopeY)) * 180. / M_PI : -1., pix + 1., e}, id});
        (inTime ? entries : sideLu).back().lu = (src & 16) && kind == 2;
        continue;
      }
      if (!tofMid || t - tMid < -5. || t - tMid > 30.) continue;   // a middle-TOF hit in coincidence
      if (fastOnly && (kind == 0 || ekin <= eSplit)) continue;
      // the TOF response of muonic_xray_stack_cuts.C with its seeds, then the top-TOF energy
      random.SetSeed(1000003ULL * proc + event + 1);
      const double eMidMeasured = eMid * (1. + sigmaE * random.Gaus());
      const double dtTof = tofTop ? (tMid + sigmaT * random.Gaus()) - (tTop + sigmaT * random.Gaus()) : -1.;
      const double eTopMeasured = tofTop ? eTop * (1. + sigmaE * random.Gaus()) : -1.;
      const bool timed = tofTop && tTop < tMid + 5.;
      const bool bottom = botFirst >= 0. && botFirst - tMid > -5. && botFirst - tMid < 50.;
      entries.push_back({1. / skimDays, kind, event & 1, fastOnly ? 1 : 0, true, timed, e,
                         {bottom ? 1. : 0., double(barMid), eMidMeasured, dtTof, eTopMeasured, double(topLayer),
                          double(lowest), double(n), e > 0. ? emax / e : 0., rmax, 0., 0., 0., e}, id});
    }
    return skimDays;
  };
  const double thickDays = load(slowSkim, false), refDays = load(ref, true);
  if (thickDays <= 0. || refDays <= 0.) return;
  // the measurement shown: by default the target's run's own days; every entry (loaded per day of its sample) then
  // counts the events of it that it stands for: 1 in the target's run, thickDays / refDays for the faster muons (until
  // their run has as many days)
  if (measureDays <= 0.) measureDays = thickDays;
  for (auto& x : entries) x.w *= measureDays;
  for (auto& x : sideLu) x.w *= measureDays;
  const std::string inDays = Form("in %.2f days", measureDays);
  // the TOF trigger takes the dataset out of the middle-TOF coincidences: what it keeps in the L_alpha window (in
  // measureDays days), then only what it keeps
  double coinS = 0., coinB = 0., trigS = 0., trigB = 0.;
  for (const auto& x : entries) {
    if (!InAlpha(x.e)) continue;
    (x.kind == 0 ? coinS : coinB) += x.w;
    if (x.trig) (x.kind == 0 ? trigS : trigB) += x.w;
  }
  if (gTOF) entries.erase(std::remove_if(entries.begin(), entries.end(), [](const Entry& x) { return !x.trig; }), entries.end());
  if (gTOF) sideLu.erase(std::remove_if(sideLu.begin(), sideLu.end(), [](const Entry& x) { return !x.trig; }), sideLu.end());
  const bool hasLu = !sideLu.empty() || std::any_of(entries.begin(), entries.end(), [](const Entry& x) { return x.lu; });
  std::vector<Entry> alpha;   // the target's L_alpha window
  for (const auto& x : entries) if (InAlpha(x.e)) alpha.push_back(x);
  const double fitLow = 30., fitHigh = 160.;   // the range of the L-line fits to simulated measurements (keV)
  const double kFitLow = fitHigh, kFitHigh = 500.;   // the K-line fit's (2026-11-05), from the L-line fit's end
  // the spectra's drawn range (keV): 20-500 keV for both targets since 2026-11-05, so that the K-line fit shows whole
  // (Al ended at 450 keV from 2026-11-02; its higher K lines lie at 443.6-473.2 keV, Si's K_beta at 487.9 keV)
  const double plotLow = 20., plotHigh = 500.;

  // the conditions of this study, after its skims'; printed first, and kept with the outputs
  conditions.push_back(Form("study: analysis/muonic_xray_target_select.C, %s, detector %s, written %s; outputs in %s", label,
                            detector, TDatime().AsSQLString(), outDir));
  conditions.push_back(Form("  samples: %s, all entries (%.2f days, muons up to %g MeV); %s, the background of primaries above "
                            "%g MeV only (%.2f days)", slowSkim, thickDays, eSplit, ref, eSplit, refDays));
  conditions.push_back(Form("  shown: a measurement of %.2f days; a simulated entry stands for %.4g of its events (muons up to %g "
                            "MeV) or %.4g (faster muons)", measureDays, measureDays / thickDays, eSplit, measureDays / refDays));
  conditions.push_back(Form("  entries: %s entries of 20-1000 keV within -5..+30 ns of the %s middle-TOF time", entryName,
                            gTOF ? "measured" : "true"));
  if (gTOF)
    conditions.push_back("  dataset: the TOF trigger (top and middle station hit as measured, within 20 ns), a fixed first cut");
  else
    conditions.push_back(Form("  TOF response applied in this study: sigma %g ns per station, sigma %g on the bar energies, seeds "
                              "1000003 x process + event + 1; no trigger", sigmaT, sigmaE));
  conditions.push_back(Form("  signal: entries with a muonic X-ray of the target (src bit 1); L_alpha window %.2f +- %.2f keV (2 sigma "
                            "of %g keV x sqrt(E / MeV))", kAlpha, 2. * SigmaKeV(kAlpha), sigma1MeV));
  {
    std::string muon = "  muon variables:", entry = Form("  %s variables:", entryName);
    for (int i = 0; i < kNV; ++i) {
      if (!Used(i)) continue;
      std::string& list = IsMuon(i) ? muon : entry;
      list += Form("%s %d %s", list.back() == ':' ? "" : ",", i, kVars[i].name);
    }
    conditions.push_back(muon + " (0: within -5..+50 ns of the middle-TOF time" +
                         (gTOF ? "; 1: max(|x|,|y|) of the track where it crosses the target's centre plane; 3: from the "
                                 "top-to-middle time of flight, for a muon" : "") +
                         (gBottomBar ? "; 10: the largest bar of a bottom-TOF hit in that window, 0 if none, shown, not cut" : "") +
                         (gTOF ? "; 11: the zenith angle of the TOF track)" : ")"));
    conditions.push_back(entry + (gCZT && gPixel ? " (12: ix x 4 + iy + 1, ix and iy counting along x and y from the -x, -y "
                                                   "corner; shown, not cut; " : " (") +
                         "13: the entry's energy with the response, shown, not cut)");
  }
  conditions.push_back(Form("  cuts: greedy on S/sqrt(B) %s in the L_alpha window, each step at least 1%% better and keeping "
                            "at least 30 signal entries (60 on all events), chosen on even and odd events alone to cross-validate",
                            inDays.c_str()));
  conditions.push_back(Form("  fits: to the simulated spectrum after the cuts %s (the MC truth itself), %g-%g keV (the L lines) and "
                            "%g-%g keV (the K lines)", inDays.c_str(), fitLow, fitHigh, kFitLow, kFitHigh));
  // the analysis cuts of the pixel skims with the measured TOF: a cut profile (2026-11-05), named here with where it
  // comes from (its notes up to the format's)
  CutProfile profile;
  if (gTOF && gCZT && gPixel) {
    if (!profile.Load(cutProfile)) return;
    conditions.push_back(Form("  analysis cuts: the cut profile %s (%s): %s; the most efficient cuts printed beside them",
                              profile.name.c_str(), profile.path.c_str(), profile.Describe().c_str()));
    for (const auto& note : profile.notes) {
      if (note.rfind("A cut keeps", 0) == 0) break;
      conditions.push_back("    " + note);
    }
  }
  printf("\n=== conditions ===\n");
  for (const auto& line : conditions) printf("%s\n", line.c_str());
  {
    std::ofstream file((out + "/muonic_xray_select_conditions.txt").Data());
    for (const auto& line : conditions) file << line << "\n";
  }

  const char* kindNames[3] = {"muonic X-rays of the target", "muonic X-rays of other volumes", "the rest"};
  printf("\n=== %s: simulated %s entries with a middle-TOF hit in coincidence%s ===\n", label,
         gOtherPanel ? "pixel" : gCZT ? "CZT" : "stack",
         gTOF ? " and the TOF trigger" : "");
  if (gTOF)
    printf("the trigger (top and middle) keeps in the L_alpha window %s: target X-rays %.1f of %.1f (%.0f%%), everything else "
           "%.1f of %.1f (%.0f%%)\n", inDays.c_str(), trigS, coinS, 100. * trigS / std::max(1e-12, coinS), trigB, coinB,
           100. * trigB / std::max(1e-12, coinB));
  printf("%s (%.2f days, muons up to %.0f MeV) and the faster muons of %s (%.2f days, scaled to %.2f); energy sigma %.2f keV at "
         "1 MeV (L_alpha window +-%.1f keV); %s (simulated entries):\n", slowSkim, thickDays, eSplit, ref, refDays, measureDays,
         gSigma1MeV, 2. * SigmaKeV(kAlpha), inDays.c_str());
  for (int k = 0; k < 3; ++k) {
    double all = 0., inAlpha = 0.;
    long nAll = 0, nAlpha = 0;
    for (const auto& x : entries) if (x.kind == k) { all += x.w; ++nAll; if (InAlpha(x.e)) { inAlpha += x.w; ++nAlpha; } }
    printf("   %-32s %10.1f (%7ld entries); in the L_alpha window %9.1f (%6ld entries)\n", kindNames[k], all, nAll, inAlpha,
           nAlpha);
  }

  // ---- S/sqrt(B) in measureDays days in the L_alpha window, from half h (0, 1; each with twice its weight) or both (2);
  // S and B: their entries in measureDays days
  auto score = [&](const Cuts& cuts, int h, double* sOut = nullptr, double* bOut = nullptr, long* nSignal = nullptr) {
    double s = 0., b = 0.;
    long n = 0;
    for (const auto& x : alpha) {
      if (h < 2 && x.half != h) continue;
      if (!PassesAll(x, cuts)) continue;
      const double w = h < 2 ? 2. * x.w : x.w;
      if (x.kind == 0) { s += w; ++n; }
      else b += w;
    }
    if (sOut) *sOut = s;
    if (bOut) *bOut = b;
    if (nSignal) *nSignal = n;
    return b > 0. ? s / std::sqrt(b) : 0.;
  };
  // the greedy optimisation on half h (or both): its cuts and the order they were chosen in
  auto optimise = [&](int h, std::vector<int>& order) {
    Cuts cuts{};
    order.clear();
    const long minSignal = h < 2 ? 30 : 60;   // signal entries a cut must keep
    double fom = score(cuts, h);
    for (;;) {
      int bestVar = -1;
      Cut bestCut;
      double bestF = fom;
      for (int i = 0; i < kNV; ++i) {
        if (cuts[i].active || !Used(i)) continue;
        for (const auto& option : Options(i)) {
          Cuts trial = cuts;
          trial[i] = option;
          long n = 0;
          const double f = score(trial, h, nullptr, nullptr, &n);
          if (n >= minSignal && f > bestF) { bestF = f; bestVar = i; bestCut = option; }
        }
      }
      if (bestVar < 0 || bestF < 1.01 * fom) break;   // less than 1% better: stop
      cuts[bestVar] = bestCut;
      order.push_back(bestVar);
      fom = bestF;
    }
    // each chosen threshold again with the others fixed (dropping the cut if that is better), until nothing changes
    for (int pass = 0; pass < 5; ++pass) {
      bool changed = false;
      for (int i : order) {
        Cut best = cuts[i];
        double bestF = score(cuts, h);
        std::vector<Cut> options = Options(i);
        options.push_back(Cut{});
        for (const auto& option : options) {
          Cuts trial = cuts;
          trial[i] = option;
          long n = 0;
          const double f = score(trial, h, nullptr, nullptr, &n);
          if (n >= minSignal && f > bestF * (1. + 1e-9)) { bestF = f; best = option; }
        }
        if (best.active != cuts[i].active || best.lo != cuts[i].lo || best.hi != cuts[i].hi) { cuts[i] = best; changed = true; }
      }
      if (!changed) break;
    }
    order.erase(std::remove_if(order.begin(), order.end(), [&](int i) { return !cuts[i].active; }), order.end());
    return cuts;
  };
  std::vector<int> order0, order1, bestOrder;
  const Cuts cuts0 = optimise(0, order0), cuts1 = optimise(1, order1), best = optimise(2, bestOrder);
  const double zCv = 0.5 * (score(cuts0, 1) + score(cuts1, 0));
  // the cuts of the analysis: for the CZT run of single pixels with the measured TOF, chosen by hand on 2026-11-02 and
  // optimised for the L_alpha peak the same day (the work folder's output/muonic_xray_czt/pixel_tof/cut_optimization.C: the
  // counting significance of the true counts under the peak; the features both halves of the events chose alone, added
  // to the cuts chosen by hand before), in the order of the corner plot; else the most efficient cuts
  const bool handChosen = gTOF && gCZT && gPixel;
  Cuts cuts = best;
  std::vector<int> order = bestOrder;
  if (handChosen) {   // since 2026-11-05 from the cut profile, loaded above (czt: the cuts above, as before)
    cuts = Cuts{};
    order = {1, 11, 4, 2, 3, 13, 6, 5};   // the vector of the corner plot
    auto keep = [&](int i, const char* key, double halfStep) {
      if (!profile.Has(key)) return;
      cuts[i].active = true; cuts[i].lo = profile.Lo(key); cuts[i].hi = profile.Hi(key) + halfStep;
    };
    keep(1, "track", 0.);    // the track at the target (cm)
    keep(11, "zenith", 0.);  // its zenith angle (deg)
    keep(4, "topE", 0.);     // top-TOF energy (MeV)
    keep(2, "midE", 0.);     // middle-TOF energy (MeV)
    keep(3, "tkin", 0.);     // TOF kinetic energy (MeV)
    keep(13, "entryE", 0.);  // the pixel entry (keV)
    keep(6, "otherE", 0.);   // the other pixels (MeV)
    keep(5, "npix", 0.5);    // pixels hit in the event: at most hi (between the bins, as the optimisation's)
  }
  const std::string cutsLabel = handChosen ? "the cuts of this analysis (profile " + profile.name + ")" : "";
  const char* cutsName = handChosen ? cutsLabel.c_str() : "the most efficient cuts (chosen in the L-alpha window)";
  // the cut sets of the earlier studies: the 1 cm target's, and muonic_xray_target_cuts.C's own for this target
  auto earlier = [](double eMidCut, double tofCut) {
    Cuts c{};
    c[0].active = true; c[0].hi = 0.5;
    c[1].active = true; c[1].lo = 1.5; c[1].hi = 6.5;
    c[2].active = true; c[2].lo = eMidCut;
    c[3].active = true; c[3].lo = tofCut;
    if (Used(5)) { c[5].active = true; c[5].lo = 2.5; }   // the highest stack layer: none in the CZT
    return c;
  };
  const Cuts cuts1cm = earlier(3.0, 1.3), cutsOwn = earlier(1.75, 1.1);

  printf("\nthe most efficient cuts (greedy, S/sqrt(B) %s in the L_alpha window; even and odd events alone):\n", inDays.c_str());
  auto list = [&](const char* label, const Cuts& c, const std::vector<int>& o) {
    printf("   %-12s", label);
    for (std::size_t k = 0; k < o.size(); ++k) printf("%s%s", k ? "; " : "", Describe(c, o[k]).c_str());
    printf("\n");
  };
  list("all events:", best, bestOrder);
  list("even events:", cuts0, order0);
  list("odd events:", cuts1, order1);
  printf("   S/sqrt(B) on the events they were chosen on %.2f; cross-validated (even's cuts on the odd events and back) %.2f\n",
         score(best, 2), zCv);
  if (handChosen) {
    printf("the cuts of this analysis (optimised for the Al L_alpha peak on 2026-11-02, cut_optimization.C%s):\n",
           std::string(target) == "Al" ? "" : "; used unchanged here");
    list("", cuts, order);
  }

  // the cut flow of the chosen cuts, in the order they were chosen; with the measured TOF the trigger first
  std::vector<std::string> flowNames = {"middle-TOF hit in coincidence"};
  std::vector<double> flowS, flowB, flowZ;
  {
    Cuts partial{};
    double s, b;
    if (gTOF) {
      flowS.push_back(coinS); flowB.push_back(coinB);
      flowZ.push_back(coinB > 0. ? coinS / std::sqrt(coinB) : 0.);
      flowNames.push_back("TOF trigger: top and middle");
    }
    flowZ.push_back(score(partial, 2, &s, &b)); flowS.push_back(s); flowB.push_back(b);
    for (int i : order) {
      partial[i] = cuts[i];
      flowZ.push_back(score(partial, 2, &s, &b)); flowS.push_back(s); flowB.push_back(b);
      flowNames.push_back(Describe(cuts, i));
    }
  }
  printf("\ncut flow in the L_alpha window, %s:\n   %-40s %9s %9s %16s\n", inDays.c_str(), "", "S", "B", "S/sqrt(B)");
  std::ofstream csv((out + "/muonic_xray_select.csv").Data());
  // S and B in measureDays days (the conditions give them); the fits' rows: their L_alpha significance, amount and error
  csv << "table,step,cut,S,B,S_sqrtB,Z_fit,amount,error\n";
  for (std::size_t k = 0; k < conditions.size(); ++k) {   // the conditions first, one row per line
    std::string line = conditions[k];
    std::replace(line.begin(), line.end(), '"', '\'');
    csv << "condition," << k << ",\"" << line << "\",,,,,,\n";
  }
  for (std::size_t k = 0; k < flowS.size(); ++k) {
    printf("   %zu %-38s %9.1f %9.1f %16.2f\n", k, flowNames[k].c_str(), flowS[k], flowB[k], flowZ[k]);
    csv << "cutflow," << k << ",\"" << flowNames[k] << "\"," << flowS[k] << "," << flowB[k] << "," << flowZ[k] << ",,,\n";
  }

  // ---- the spectra by what made them, in measureDays days, 0-200 keV in 2 keV
  auto spectrum = [&](const Cuts& c, const char* tag, int bins = 100, double hi = 200.) {
    std::array<TH1D*, 3> h;
    for (int k = 0; k < 3; ++k) {
      h[k] = new TH1D(Form("spectrum_%s_%d", tag, k), "", bins, 0., hi);
      h[k]->Sumw2();
    }
    for (const auto& x : entries) if (PassesAll(x, c)) h[x.kind]->Fill(x.e, x.w);
    return h;
  };
  const auto after = spectrum(cuts, "after");   // the simulated measurements' (the plots': at 0-1000 keV, below)

  // what the L_alpha window holds after the cuts, in measureDays days: by kind, and the rest by sample
  {
    double kinds[3] = {}, rest[2] = {};
    for (const auto& x : alpha) {
      if (!PassesAll(x, cuts)) continue;
      kinds[x.kind] += x.w;
      if (x.kind == 2) rest[x.sample] += x.w;
    }
    printf("\nthe L_alpha window after the cuts, %s: %s %.1f, %s %.1f, %s %.1f (from muons up to %.0f MeV %.1f, from faster "
           "ones %.1f)\n", inDays.c_str(), kindNames[0], kinds[0], kindNames[1], kinds[1], kindNames[2], kinds[2], eSplit, rest[0],
           rest[1]);
  }
  // the pixel panel's own radioactivity (2026-11-05, the LYSO's Lu-176 of /tof/intrinsic/; part of the rest): its
  // entries in the coincidence window as counted, and as the decays alone in the sidebands estimate them, before the
  // cuts, after them and in the L_alpha window after them, beside the rest
  if (hasLu) {
    struct Sum { double w = 0.; long raw = 0; };
    auto sum = [&](const std::vector<Entry>& v, bool luOnly, const Cuts& c, bool window) {
      Sum r;
      for (const auto& x : v)
        if ((!luOnly || x.lu) && (luOnly || x.kind == 2) && PassesAll(x, c) && (!window || InAlpha(x.e))) { r.w += x.w; ++r.raw; }
      return r;
    };
    printf("\nthe panel's own radioactivity (/tof/intrinsic/, truth bit 16: a decay alone or piled up on a muon's hit; part of the "
           "rest), %s (simulated entries in brackets):\n   %-42s %22s %22s %22s\n", inDays.c_str(), "", "counted in -5..+30 ns",
           "from the sidebands", "the rest (all of it)");
    const char* stepNames[3] = {"before the cuts (20-1000 keV)", "after the cuts", "after the cuts, in the L_alpha window"};
    for (int k = 0; k < 3; ++k) {
      const Cuts& c = k ? cuts : Cuts{};
      const Sum a = sum(entries, true, c, k == 2), b = sum(sideLu, false, c, k == 2), r = sum(entries, false, c, k == 2);
      printf("   %-42s %12.1f (%7ld) %12.1f (%7ld) %12.1f (%7ld)\n", stepNames[k], a.w, a.raw, b.w, b.raw, r.w, r.raw);
      csv << "intrinsic," << k << ",\"" << stepNames[k] << "\"," << a.w << "," << b.w << "," << r.w << ",,,\n";
    }
    printf("   (the sidebands: the decays alone at -990..-10 and +40..+990 ns of the middle TOF, times %.4f; after the cuts they may\n"
           "    overestimate a little: their other-pixel cuts see the hits around them, not the muon's prompt ones)\n", kSideScale);
  }

  // ---- fits without truth to the simulated spectrum itself: what the simulation gives in measureDays days, bin by bin
  // (the MC truth after the cuts; asked for on 2026-11-02 in place of simulated measurements drawn from it; by default
  // the target's run as simulated, the faster muons scaled to its days), in three ways:
  //  bins:  background exp(pol3) and two Gaussian peaks, each with its own amount, position and width, all from the
  //         simulated bins alone (asked for on 2026-11-03 in place of the free fit with L_beta tied to L_alpha: no line
  //         energy or spacing known in advance; the fit of measurement_fit): its start values from the bins (the
  //         background fitted alone; peak 1 at the largest excess over it, fitted with the background; peak 2 where a
  //         scan of its position over the bins, at peak 1's width, finds the best likelihood), its limits the fit range
  //         for the positions and, since 2026-11-05, for the widths the detector's resolution at each peak's own
  //         position (kWidthMin-kWidthMax times it; 1-1: fixed at it; the resolution known from calibration, not
  //         from these data: LYSO's is 4 times CZT's); peak 1 the larger; each peak's significance against the same
  //         fit without it;
  //  free:  the same with free widths, half a bin to 20 keV (the bins fit before 2026-11-05), for comparison;
  //  known: background exp(pol3), the Gaussian lines at their known energies with the stack's known width, only their
  //         amounts free (as a calibrated measurement fits them; an amount may come out negative).
  // Left out on 2026-11-02: a background shape known in advance (the simulated rest after the cuts: a cosmic background
  // cannot be predicted that exactly, and its template missed the target's own continuum), and the background shape of
  // the measurement's own control region (middle-TOF energy below its cut), which holds more of the target's X-rays
  // than background: it found L_alpha in the null test (median 4.2 sigma) and 23 sigma in a simulated measurement.
  // The L_alpha significance: the likelihood ratio against the same fit without L_alpha, with the sign of its amount.
  // A null test: the expected spectrum without any muonic X-ray, where an honest fit finds no L_alpha. Fit range:
  // fitLow-fitHigh. The background exp(pol2) instead with backgroundDegree 2. The fits take the simulated spectrum as exact: the simulation's own statistics (the plots' error bars) are
  // not in them.
  auto gaussian = [](double e, double pos, double width) {
    return std::exp(-0.5 * std::pow((e - pos) / width, 2)) / (std::sqrt(2. * M_PI) * width);
  };
  // the bins fit's widths (2026-11-05, asked for: the detector's energy resolution is known in advance, from its
  // calibration sources, not from these data): res, each peak's width the resolution at its own position,
  // SigmaKeV(position), times a parameter within kWidthMin-kWidthMax (equal: fixed). Fixed at 1, as a muonic X-ray line
  // is far narrower than any detector's resolution: a single-bin spike can no longer pass for a peak, nor can a broad
  // shoulder (with 1-4 times the resolution the LYSO's peak took in the target X-rays' own Compton shoulder below the
  // line and the Lu K X-rays of the next pixels: 1054 entries at 69.7 keV for 322 full-energy Si L_alpha entries).
  // free: half a bin to 20 keV (the bins fit before 2026-11-05)
  const double kWidthMin = 1., kWidthMax = 1.;   // equal: the widths fixed at the resolution
  auto binsModel = [&](const char* name, bool res) {   // p4-6, p7-9: the peaks' amounts (entries), positions and widths
    return new TF1(name, [=](double* x, double* p) {
      const double u = (x[0] - 80.) / 50.;
      const double w1 = res ? p[6] * SigmaKeV(p[5]) : p[6], w2 = res ? p[9] * SigmaKeV(p[8]) : p[9];
      return std::exp(p[0] + p[1] * u + p[2] * u * u + p[3] * u * u * u) +
             2. * (p[4] * gaussian(x[0], p[5], w1) + p[7] * gaussian(x[0], p[8], w2));
    }, fitLow, fitHigh, 10);
  };
  auto knownModel = [&](const char* name) {
    return new TF1(name, [=](double* x, double* p) {
      const double u = (x[0] - 80.) / 50.;
      return std::exp(p[0] + p[1] * u + p[2] * u * u + p[3] * u * u * u) +
             2. * (p[4] * gaussian(x[0], kAlpha, SigmaKeV(kAlpha)) + p[5] * gaussian(x[0], kBeta, SigmaKeV(kBeta)));
    }, fitLow, fitHigh, 6);
  };
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");   // classic Minuit, as muonic_xray_stack_cuts.C
  const int backDegree = backgroundDegree == 2 ? 2 : 3;
  const char* backName = backDegree == 2 ? "exp(pol2)" : "exp(pol3)";
  auto doFit = [&](TH1D* h, TF1* f) {   // every fit: the Poisson likelihood in the fit range; exp(pol2): p3 held at 0
    if (backDegree == 2) f->FixParameter(3, 0.);
    h->Fit(f, "LQNR");
  };
  TF1 *fBins = binsModel("bins", true), *fBinsNull = binsModel("binsNull", true), *fBinsBack = binsModel("binsBackground", true);
  TF1 *fFree = binsModel("free", false), *fFreeNull = binsModel("freeNull", false), *fFreeBack = binsModel("freeBackground", false);
  TF1 *fKnown = knownModel("known"), *fKnownNull = knownModel("knownNull");
  fKnownNull->FixParameter(4, 0.);
  auto nll = [&](TH1D* h, TF1* f) {   // Poisson -log L of the fit range (f's), at the bin centres as the fit
    double sum = 0.;
    for (int b = h->FindBin(f->GetXmin() + 1e-6); b <= h->FindBin(f->GetXmax() - 1e-6); ++b) {
      const double mu = std::max(1e-12, f->Eval(h->GetBinCenter(b)));
      sum += mu - h->GetBinContent(b) * std::log(mu);
    }
    return sum;
  };
  struct LineFit {
    double z, amount, error, position, width;                    // L_alpha (known), or the bins fit's peak 1
    double positionError = 0., widthError = 0.;
    double amount2 = 0., error2 = 0., position2 = 0., positionError2 = 0., width2 = 0., widthError2 = 0.;   // its peak 2
    double z2 = 0.;                                                                                      // and its significance
  };
  // the fit from the bins alone with the functions F (the fit), FN (its variants) and FB (the background alone); res: the
  // widths bounded by the resolution, else free
  auto fitBins = [&](TH1D* h, TF1* F, TF1* FN, TF1* FB, bool res) {
    const int bFirst = h->FindBin(fitLow + 1e-6), bLast = h->FindBin(fitHigh - 1e-6);
    const double binW = h->GetBinWidth(1), total = std::max(1., h->Integral(bFirst, bLast));
    const double wLo = res ? kWidthMin : 0.5 * binW, wHi = res ? kWidthMax : 20.;   // the width parameters' limits
    auto limitWidth = [&](TF1* f, int p) {   // a width parameter within its limits, or fixed when they are one value
      if (wHi > wLo) f->SetParLimits(p, wLo, wHi);
      else f->FixParameter(p, wLo);
    };
    auto widthPar = [&](double keV, double pos) { return std::min(wHi, std::max(wLo, res ? keV / SigmaKeV(pos) : keV)); };
    auto widthKeV = [&](double par, double pos) { return res ? par * SigmaKeV(pos) : par; };
    // the background alone first, the peaks' parameters fixed at no amount
    FB->SetParameters(std::log(std::max(1e-3, total / (bLast - bFirst + 1))), 0., 0., 0., 0., fitLow, widthPar(binW, fitLow), 0.,
                      fitHigh, widthPar(binW, fitHigh));
    for (int p = 4; p < 10; ++p) FB->FixParameter(p, FB->GetParameter(p));
    doFit(h, FB);
    // peak 1's start: the bin of the largest excess over that background, its width from the bins around it above
    // half that excess, its amount their excess (76% of a Gaussian lies within its half maximum), its position their
    // excess-weighted centre
    auto excess = [&](int b) { return h->GetBinContent(b) - FB->Eval(h->GetBinCenter(b)); };
    int top = bFirst;
    for (int b = bFirst; b <= bLast; ++b) if (excess(b) > excess(top)) top = b;
    int lo = top, hi = top;
    while (lo > bFirst && excess(lo - 1) > 0.5 * excess(top)) --lo;
    while (hi < bLast && excess(hi + 1) > 0.5 * excess(top)) ++hi;
    double sum = 0., centre = 0.;
    for (int b = lo; b <= hi; ++b) { const double x = std::max(0., excess(b)); sum += x; centre += x * h->GetBinCenter(b); }
    // the background and peak 1 first, peak 2 at no amount
    for (int p = 0; p < 10; ++p) F->ReleaseParameter(p);
    for (int p = 0; p < 4; ++p) F->SetParameter(p, FB->GetParameter(p));
    const double pos1 = sum > 0. ? centre / sum : h->GetBinCenter(top);
    F->SetParameter(4, std::max(1e-3, sum / 0.76)); F->SetParameter(5, pos1);
    F->SetParameter(6, widthPar((hi - lo + 1) * binW / 2.3548, pos1));
    F->SetParLimits(4, 0., total); F->SetParLimits(5, fitLow, fitHigh); limitWidth(F, 6);
    F->FixParameter(7, 0.); F->FixParameter(8, fitHigh); F->FixParameter(9, F->GetParameter(6));
    doFit(h, F);
    // peak 2's start: its position scanned over every bin more than 3 of peak 1's widths and a bin from it, with peak
    // 1's width (parameter), its amount and the background refitted at each (peak 1 kept); the best likelihood starts
    // the full fit (a scan, not the largest excess: a single-bin spike must not win over a peak)
    const double mu1 = F->GetParameter(5), w1 = F->GetParameter(6);
    double amount2 = 1e-3, mu2 = fitHigh, bestNll = 1e300;
    for (int b = bFirst; b <= bLast; ++b) {
      const double c = h->GetBinCenter(b);
      if (std::abs(c - mu1) < 3. * widthKeV(w1, mu1) + binW) continue;
      for (int p = 0; p < 10; ++p) { FN->ReleaseParameter(p); FN->SetParameter(p, F->GetParameter(p)); }
      for (int p = 4; p < 7; ++p) FN->FixParameter(p, F->GetParameter(p));
      FN->SetParameter(7, std::max(1., h->GetBinContent(b) - F->Eval(c))); FN->SetParLimits(7, 0., total);
      FN->FixParameter(8, c); FN->FixParameter(9, w1);
      doFit(h, FN);
      const double value = nll(h, FN);
      if (value < bestNll) { bestNll = value; amount2 = std::max(1e-3, FN->GetParameter(7)); mu2 = c; }
    }
    // the full fit, every parameter free within its limits
    for (int p = 7; p < 10; ++p) F->ReleaseParameter(p);
    F->SetParameter(7, amount2); F->SetParameter(8, mu2); F->SetParameter(9, w1);
    for (int p : {4, 7}) F->SetParLimits(p, 0., total);
    for (int p : {5, 8}) F->SetParLimits(p, fitLow, fitHigh);
    for (int p : {6, 9}) limitWidth(F, p);
    doFit(h, F);
    const int a = F->GetParameter(7) > F->GetParameter(4) ? 7 : 4, o = a == 4 ? 7 : 4;   // peak 1, the larger
    // each peak's significance: the same fit without it, the other peak held at its position and width (free, it would
    // move over to take the removed peak's place)
    auto without = [&](int k) {
      const int other = k == 4 ? 7 : 4;
      for (int p = 0; p < 10; ++p) { FN->ReleaseParameter(p); FN->SetParameter(p, F->GetParameter(p)); }
      for (int j = 0; j < 3; ++j) FN->FixParameter(k + j, j ? F->GetParameter(k + j) : 0.);
      FN->SetParLimits(other, 0., total);
      FN->FixParameter(other + 1, F->GetParameter(other + 1)); FN->FixParameter(other + 2, F->GetParameter(other + 2));
      doFit(h, FN);
      return std::sqrt(std::max(0., 2. * (nll(h, FN) - nll(h, F))));
    };
    LineFit r{without(a), F->GetParameter(a), F->GetParError(a), F->GetParameter(a + 1),
              widthKeV(F->GetParameter(a + 2), F->GetParameter(a + 1))};
    r.z2 = without(o);
    r.positionError = F->GetParError(a + 1); r.widthError = widthKeV(F->GetParError(a + 2), F->GetParameter(a + 1));
    r.amount2 = F->GetParameter(o); r.error2 = F->GetParError(o);
    r.position2 = F->GetParameter(o + 1); r.positionError2 = F->GetParError(o + 1);
    r.width2 = widthKeV(F->GetParameter(o + 2), F->GetParameter(o + 1));
    r.widthError2 = widthKeV(F->GetParError(o + 2), F->GetParameter(o + 1));
    return r;
  };
  auto fitMeasurement = [&](TH1D* h, int how) {   // 0 the bins alone, widths from the resolution; 2 free widths; 1 known
    if (how == 0) return fitBins(h, fBins, fBinsNull, fBinsBack, true);
    if (how == 2) return fitBins(h, fFree, fFreeNull, fFreeBack, false);
    double window = 0.;
    for (int b = h->FindBin(kAlpha - 2. * SigmaKeV(kAlpha)); b <= h->FindBin(kAlpha + 2. * SigmaKeV(kAlpha)); ++b) window += h->GetBinContent(b);
    const double back = std::log(std::max(1., h->GetBinContent(h->FindBin(80.))));
    fKnown->SetParameters(back, -0.3, 0., 0., 0.1 * window, 0.02 * window);
    fKnown->SetParLimits(4, -0.3 * window, 2. * window); fKnown->SetParLimits(5, -0.3 * window, 2. * window);
    doFit(h, fKnown);
    for (int p : {0, 1, 2, 3, 5}) fKnownNull->SetParameter(p, fKnown->GetParameter(p));
    fKnownNull->SetParLimits(5, -0.3 * window, 2. * window);
    doFit(h, fKnownNull);
    const double z = std::sqrt(std::max(0., 2. * (nll(h, fKnownNull) - nll(h, fKnown))));
    return LineFit{fKnown->GetParameter(4) < 0. ? -z : z, fKnown->GetParameter(4), fKnown->GetParError(4), kAlpha, SigmaKeV(kAlpha)};
  };
  // the simulated spectrum in measureDays days (entries per 2 keV), with or without its muonic X-rays
  auto expected = [&](const std::array<TH1D*, 3>& h, bool withXrays, const char* name) {
    auto m = new TH1D(name, "", h[0]->GetNbinsX(), h[0]->GetXaxis()->GetXmin(), h[0]->GetXaxis()->GetXmax());
    for (int b = 1; b <= m->GetNbinsX(); ++b)
      m->SetBinContent(b, h[2]->GetBinContent(b) + (withXrays ? h[0]->GetBinContent(b) + h[1]->GetBinContent(b) : 0.));
    return m;
  };
  double truthAlpha[3] = {}, truthFit[3] = {};   // muonic X-ray entries after the cuts in measureDays days: window, fit range
  for (const auto& x : entries) {
    if (!PassesAll(x, cuts)) continue;
    if (InAlpha(x.e)) truthAlpha[x.kind] += x.w;
    if (x.e >= fitLow && x.e < fitHigh) truthFit[x.kind] += x.w;
  }
  printf("\nfits without truth to the simulated spectrum after the cuts %s (the MC truth):\n", inDays.c_str());
  printf("   truth after the cuts: L_alpha window: target X-rays %.0f, other muonic X-rays %.0f, the rest %.0f; %.0f-%.0f keV: "
         "target X-rays %.0f, other muonic X-rays %.0f, the rest %.0f\n", truthAlpha[0], truthAlpha[1], truthAlpha[2], fitLow, fitHigh,
         truthFit[0], truthFit[1], truthFit[2]);
  const std::string widthWords = kWidthMax > kWidthMin ? Form("widths %g-%g x resolution", kWidthMin, kWidthMax)
                                                       : std::string("widths = resolution");
  const std::string modelNames[3] = {Form("bins alone, %s, %s", widthWords.c_str(), backName),
                                     Form("known energies and width, background %s", backName),
                                     Form("bins alone, free widths, background %s", backName)};
  TH1D *expect = expected(after, true, "expected"), *expectNull = expected(after, false, "expectedNull");
  // how well a fit f (as last fitted) describes a simulated spectrum, against the statistics of measureDays days
  // (Baker-Cousins: twice the log-likelihood ratio to a perfect description; a fluctuation of the simulation larger than
  // those of measureDays days counts too)
  const int nBinsFit = expect->FindBin(fitHigh - 1e-6) - expect->FindBin(fitLow + 1e-6) + 1;
  // the free parameters: the known-energy fit 6, the bins fits 10; fewer with exp(pol2) (p3 held) and the widths fixed
  const int backHeld = backDegree == 2 ? 1 : 0, widthsHeld = kWidthMax > kWidthMin ? 0 : 2;
  const int ndfFit = nBinsFit - 6 + backHeld, ndfFree = nBinsFit - 10 + backHeld, ndfBins = ndfFree + widthsHeld;
  auto chi2Of = [&](TH1D* h, TF1* f) {   // over f's range
    double chi2 = 0.;
    for (int b = h->FindBin(f->GetXmin() + 1e-6); b <= h->FindBin(f->GetXmax() - 1e-6); ++b) {
      const double nB = h->GetBinContent(b), fB = std::max(1e-12, f->Eval(h->GetBinCenter(b)));
      chi2 += 2. * (fB - nB + (nB > 0. ? nB * std::log(nB / fB) : 0.));
    }
    return chi2;
  };
  const int noteLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kWarning;   // not Minuit's notes on a start value outside its limits (the null fits' fixed L_alpha)
  // the cut sets compared, each with the known-energy fit to its own expected spectrum
  printf("\nthe cut sets compared: S, B and S/sqrt(B) in the L_alpha window, and the fit with the known energies and width "
         "(background %s) to the simulated spectrum, %s:\n", backName, inDays.c_str());
  printf("   %-46s %8s %8s %10s %18s %14s %13s\n", "", "S", "B", "S/sqrt(B)", "L_alpha amount", "significance", "chi2/ndf");
  struct Set { const char* name; const Cuts* cuts; };
  std::vector<Set> sets;
  if (!gTOF) {   // the earlier sets cut on the middle-TOF bar and the time of flight, which the measured TOF replaces
    sets.push_back({"cuts of the 1 cm target", &cuts1cm});
    sets.push_back({"own cuts of muonic_xray_target_cuts.C", &cutsOwn});
  }
  sets.push_back({"the most efficient cuts (this optimisation)", &best});
  if (handChosen) sets.push_back({cutsName, &cuts});
  for (std::size_t k = 0; k < sets.size(); ++k) {
    double s, b;
    const double z = score(*sets[k].cuts, 2, &s, &b);
    TH1D* e = expected(spectrum(*sets[k].cuts, Form("set%zu", k)), true, Form("expectedSet%zu", k));
    const LineFit r = fitMeasurement(e, 1);
    const double c2 = chi2Of(e, fKnown);
    printf("   %-46s %8.1f %8.1f %10.2f %10.1f +- %4.1f %8.1f sigma %7.1f/%d\n", sets[k].name, s, b, z, r.amount, r.error, r.z, c2,
           ndfFit);
    csv << "compare," << k << ",\"" << sets[k].name << "\"," << s << "," << b << "," << z << "," << r.z << "," << r.amount << ","
        << r.error << "\n";
    delete e;
  }
  // the null test, then the fits of this analysis (each model its own functions, so fBins and fKnown keep their fits to
  // the simulated spectrum for the plots)
  LineFit fit[3], fitNull[3];
  for (int how : {0, 2, 1}) fitNull[how] = fitMeasurement(expectNull, how);
  for (int how : {0, 2, 1}) fit[how] = fitMeasurement(expect, how);
  gErrorIgnoreLevel = noteLevel;
  printf("\nthe fits to the simulated spectrum of this analysis (%.0f-%.0f keV, %s; the bins fit's peak 1):\n   %-50s %-22s %24s %18s\n",
         fitLow, fitHigh, inDays.c_str(), "", "", "L_alpha amount", "significance");
  for (int how : {0, 2, 1})
    for (int n = 0; n < 2; ++n) {
      const LineFit& r = n ? fitNull[how] : fit[how];
      printf("   %-50s %-22s %14.1f +- %6.1f %12.1f sigma\n", modelNames[how].c_str(), n ? "without muonic X-rays" : "as simulated",
             r.amount, r.error, r.z);
      csv << (n ? "null" : "expected") << "," << how << ",\"" << modelNames[how] << "\",,,," << r.z << "," << r.amount << ","
          << r.error << "\n";
    }
  const LineFit fitShown = fit[0];   // measurement_fit shows the fit from the bins alone (2026-11-03), widths from the resolution
  const double chi2Fit = chi2Of(expect, fKnown), chi2Bins = chi2Of(expect, fBins), chi2Free = chi2Of(expect, fFree);
  for (int how : {0, 2}) {
    const LineFit& r = fit[how];
    printf("   %s%s: peak 1 at %.2f +- %.2f keV, sigma %.2f +- %.2f keV (%.2f x resolution), %.1f +- %.1f entries, %.1f sigma; peak 2 "
           "at %.2f +- %.2f keV, sigma %.2f +- %.2f keV (%.2f x resolution), %.1f +- %.1f entries, %.1f sigma; chi2/ndf %.1f/%d\n",
           how ? "the fit with free widths" : "the fit from the bins alone, widths from the resolution",
           how ? "" : " (measurement_fit)", r.position, r.positionError, r.width, r.widthError, r.width / SigmaKeV(r.position),
           r.amount, r.error, r.z, r.position2, r.positionError2, r.width2, r.widthError2, r.width2 / SigmaKeV(r.position2), r.amount2,
           r.error2, r.z2, how ? chi2Free : chi2Bins, how ? ndfFree : ndfBins);
  }
  // what the simulation's own fluctuations add to it: their variance over that of measureDays days, bin by bin (a
  // simulated entry counts measureDays / its band's days entries of a measurement, so a band of fewer simulated days
  // than measureDays is noisier than the measurement; at its own days a band is exactly as noisy)
  double chi2Mc = 0.;
  for (int b = expect->FindBin(fitLow + 1e-6); b <= expect->FindBin(fitHigh - 1e-6); ++b) {
    double var = 0.;
    for (int k = 0; k < 3; ++k) var += std::pow(after[k]->GetBinError(b), 2);
    if (expect->GetBinContent(b) > 0.) chi2Mc += var / expect->GetBinContent(b);
  }
  printf("   the known-energy fit: L_beta %.1f +- %.1f; chi2/ndf %.1f/%d (Baker-Cousins, against the statistics of %.2f days); the "
         "simulation's own fluctuations alone make about %.0f of it\n", fKnown->GetParameter(5), fKnown->GetParError(5), chi2Fit,
         ndfFit, measureDays, chi2Mc);

  // ---- the K-line fit (2026-11-05, asked for: the fitting all the way to 500 keV, as the K lines may be identified
  // there): the simulated spectrum after the cuts at kFitLow-kFitHigh by the bins fit's method (the Poisson likelihood
  // at the bin centres, Gaussian peaks of the resolution's width at free positions, background exp(pol3), or exp(pol2)
  // with backgroundDegree 2, its u = (E - centre) / half the range), with kPeaks peaks: peak 1 at the largest excess
  // over the background alone, each further one where a scan of its position over the bins (more than 3 resolutions and
  // a bin from the earlier ones) finds the best likelihood, the earlier ones held, then all refitted; each peak's
  // significance against the same fit without it, the others' positions held. Three peaks: this region holds more than
  // the L region's two, K_alpha and K_beta, for Al the higher K lines (443.6-473.2 keV: as many entries as K_beta) and,
  // in the LYSO, K_alpha's escape peak (a Lu K X-ray leaving the pixel: about 50 keV below K_alpha, found at 4.6 sigma in
  // the tests of 2026-11-05; fitted apart, it no longer pulls K_alpha 2 keV low). A peak is named after the target's line
  // within 2 resolutions of it. Its own range and background: the continuum under the K lines (their Compton scatters,
  // the 511 keV line's tail) is not the L region's, and one smooth function would have to hold over a factor 17 in energy.
  // Beside it: the null test (the rest alone) and the fit at the known K_alpha and K_beta energies (only the amounts free).
  const int kPeaks = 3;
  const double kUC = 0.5 * (kFitLow + kFitHigh), kUS = 0.5 * (kFitHigh - kFitLow);
  auto kBackgroundOf = [=](const TF1* f, double e) {   // the K-line fit's smooth background
    const double u = (e - kUC) / kUS;
    return std::exp(f->GetParameter(0) + f->GetParameter(1) * u + f->GetParameter(2) * u * u + f->GetParameter(3) * u * u * u);
  };
  auto kModel = [&](const char* name, int nPeaks) {   // p0-3 the background; p(4 + 2k), p(5 + 2k): peak k's amount and position
    return new TF1(name, [=](double* x, double* p) {
      const double u = (x[0] - kUC) / kUS;
      double v = std::exp(p[0] + p[1] * u + p[2] * u * u + p[3] * u * u * u);
      for (int k = 0; k < nPeaks; ++k) v += 2. * p[4 + 2 * k] * gaussian(x[0], p[5 + 2 * k], SigmaKeV(p[5 + 2 * k]));
      return v;
    }, kFitLow, kFitHigh, 4 + 2 * nPeaks);
  };
  struct KPeak { double mu, muError, amount, error, z; };
  // the fit F (FN: its variants, FB: the background alone) of nPeaks peaks to h; the peaks by energy
  auto fitK = [&](TH1D* h, TF1* F, TF1* FN, TF1* FB, int nPeaks) {
    const int nPar = 4 + 2 * nPeaks, bFirst = h->FindBin(kFitLow + 1e-6), bLast = h->FindBin(kFitHigh - 1e-6);
    const double binW = h->GetBinWidth(1), total = std::max(1., h->Integral(bFirst, bLast));
    auto from = [&](TF1* f, const TF1* g) { for (int p = 0; p < nPar; ++p) { f->ReleaseParameter(p); f->SetParameter(p, g->GetParameter(p)); } };
    for (int p = 0; p < nPar; ++p) FB->ReleaseParameter(p);
    FB->SetParameter(0, std::log(std::max(1e-3, total / (bLast - bFirst + 1))));
    for (int p = 1; p < 4; ++p) FB->SetParameter(p, 0.);
    for (int k = 0; k < nPeaks; ++k) { FB->FixParameter(4 + 2 * k, 0.); FB->FixParameter(5 + 2 * k, kFitLow); }
    doFit(h, FB);
    // peak 1: its start at the bins of the largest excess over that background (as the bins fit's)
    auto excess = [&](int b) { return h->GetBinContent(b) - FB->Eval(h->GetBinCenter(b)); };
    int top = bFirst;
    for (int b = bFirst; b <= bLast; ++b) if (excess(b) > excess(top)) top = b;
    int lo = top, hi = top;
    while (lo > bFirst && excess(lo - 1) > 0.5 * excess(top)) --lo;
    while (hi < bLast && excess(hi + 1) > 0.5 * excess(top)) ++hi;
    double sum = 0., centre = 0.;
    for (int b = lo; b <= hi; ++b) { const double x = std::max(0., excess(b)); sum += x; centre += x * h->GetBinCenter(b); }
    from(F, FB);
    F->SetParameter(4, std::max(1e-3, sum / 0.76)); F->SetParameter(5, sum > 0. ? centre / sum : h->GetBinCenter(top));
    F->SetParLimits(4, 0., total); F->SetParLimits(5, kFitLow, kFitHigh);
    for (int k = 1; k < nPeaks; ++k) { F->FixParameter(4 + 2 * k, 0.); F->FixParameter(5 + 2 * k, kFitHigh); }
    doFit(h, F);
    // peaks 2 and on: each at the best position of a scan, the earlier ones held, then all of them refitted
    for (int k = 1; k < nPeaks; ++k) {
      double bestNll = 1e300, bestMu = kFitHigh, bestAmount = 1e-3;
      for (int b = bFirst; b <= bLast; ++b) {
        const double c = h->GetBinCenter(b);
        bool near = false;
        for (int j = 0; j < k; ++j) near = near || std::abs(c - F->GetParameter(5 + 2 * j)) < 3. * SigmaKeV(F->GetParameter(5 + 2 * j)) + binW;
        if (near) continue;
        from(FN, F);
        for (int j = 0; j < k; ++j) { FN->FixParameter(4 + 2 * j, F->GetParameter(4 + 2 * j)); FN->FixParameter(5 + 2 * j, F->GetParameter(5 + 2 * j)); }
        FN->SetParameter(4 + 2 * k, std::max(1., h->GetBinContent(b) - F->Eval(c))); FN->SetParLimits(4 + 2 * k, 0., total);
        FN->FixParameter(5 + 2 * k, c);
        for (int j = k + 1; j < nPeaks; ++j) { FN->FixParameter(4 + 2 * j, 0.); FN->FixParameter(5 + 2 * j, kFitHigh); }
        doFit(h, FN);
        const double value = nll(h, FN);
        if (value < bestNll) { bestNll = value; bestMu = c; bestAmount = std::max(1e-3, FN->GetParameter(4 + 2 * k)); }
      }
      F->ReleaseParameter(4 + 2 * k); F->ReleaseParameter(5 + 2 * k);
      F->SetParameter(4 + 2 * k, bestAmount); F->SetParameter(5 + 2 * k, bestMu);
      for (int j = 0; j <= k; ++j) { F->SetParLimits(4 + 2 * j, 0., total); F->SetParLimits(5 + 2 * j, kFitLow, kFitHigh); }
      doFit(h, F);
    }
    const double best = nll(h, F);
    std::vector<KPeak> peaks;
    for (int k = 0; k < nPeaks; ++k) {   // without peak k, the others' positions held
      from(FN, F);
      FN->FixParameter(4 + 2 * k, 0.); FN->FixParameter(5 + 2 * k, F->GetParameter(5 + 2 * k));
      for (int j = 0; j < nPeaks; ++j)
        if (j != k) { FN->SetParLimits(4 + 2 * j, 0., total); FN->FixParameter(5 + 2 * j, F->GetParameter(5 + 2 * j)); }
      doFit(h, FN);
      peaks.push_back({F->GetParameter(5 + 2 * k), F->GetParError(5 + 2 * k), F->GetParameter(4 + 2 * k), F->GetParError(4 + 2 * k),
                       std::sqrt(std::max(0., 2. * (nll(h, FN) - best)))});
    }
    std::sort(peaks.begin(), peaks.end(), [](const KPeak& a, const KPeak& b) { return a.mu < b.mu; });
    return peaks;
  };
  // the target's line within 2 resolutions of a peak: 1 K_alpha, 2 K_beta, 3 the higher K lines (the cascade model's
  // levels go as 1/n^2: K_gamma at 1.25 and the series' limit at 4/3 times K_alpha), 0 none
  auto kLineOf = [&](double mu) {
    const double tol = 2. * SigmaKeV(mu);
    if (std::abs(mu - kKAlpha) <= tol) return 1;
    if (std::abs(mu - kKBeta) <= tol) return 2;
    if (mu >= 1.25 * kKAlpha - tol && mu <= 4. / 3. * kKAlpha + tol) return 3;
    return 0;
  };
  const char* kLineText[4] = {"no target line", "K_alpha", "K_beta", "K_gamma and up"};
  const char* kLineLatex[4] = {"no target line", "K_{#alpha}", "K_{#beta}", "K_{#gamma} and up"};
  const auto afterFull = spectrum(cuts, "afterFull", 500, 1000.);   // at 0-1000 keV (also spectrum_after's)
  TH1D *expectK = expected(afterFull, true, "expectedK"), *expectKNull = expected(afterFull, false, "expectedKNull");
  TF1 *fK = kModel("kLines", kPeaks), *fKNull = kModel("kLinesNull", kPeaks), *fKVariant = kModel("kLinesVariant", kPeaks);
  TF1 *fKBack = kModel("kLinesBackground", kPeaks), *fKKnown = kModel("kKnown", 2), *fKKnownVariant = kModel("kKnownVariant", 2);
  gErrorIgnoreLevel = kWarning;
  const std::vector<KPeak> kNullPeaks = fitK(expectKNull, fKNull, fKVariant, fKBack, kPeaks);
  const std::vector<KPeak> kFound = fitK(expectK, fK, fKVariant, fKBack, kPeaks);
  double kKnownZ[2] = {};
  {   // K_alpha and K_beta at the known energies, only their amounts (which may come out negative) and the background free
    const double total = std::max(1., expectK->Integral(expectK->FindBin(kFitLow + 1e-6), expectK->FindBin(kFitHigh - 1e-6)));
    fKKnown->SetParameters(std::log(std::max(1., expectK->GetBinContent(expectK->FindBin(kUC)))), 0., 0., 0., 0.05 * total, kKAlpha,
                           0.01 * total, kKBeta);
    fKKnown->FixParameter(5, kKAlpha); fKKnown->FixParameter(7, kKBeta);
    fKKnown->SetParLimits(4, -total, total); fKKnown->SetParLimits(6, -total, total);
    doFit(expectK, fKKnown);
    const double best = nll(expectK, fKKnown);
    for (int k = 0; k < 2; ++k) {
      for (int p = 0; p < 8; ++p) { fKKnownVariant->ReleaseParameter(p); fKKnownVariant->SetParameter(p, fKKnown->GetParameter(p)); }
      fKKnownVariant->FixParameter(5, kKAlpha); fKKnownVariant->FixParameter(7, kKBeta);
      fKKnownVariant->FixParameter(4 + 2 * k, 0.); fKKnownVariant->SetParLimits(6 - 2 * k, -total, total);
      doFit(expectK, fKKnownVariant);
      kKnownZ[k] = (fKKnown->GetParameter(4 + 2 * k) < 0. ? -1. : 1.) * std::sqrt(std::max(0., 2. * (nll(expectK, fKKnownVariant) - best)));
    }
  }
  gErrorIgnoreLevel = noteLevel;
  {
    // the truth: the target's X-ray entries within 2 sigma of each K line, and in the higher K lines' range
    const double gammaLow = 1.25 * kKAlpha, limit = 4. / 3. * kKAlpha;
    double tAlpha = 0., tBeta = 0., tHigher = 0.;
    for (const auto& x : entries) {
      if (x.kind != 0 || !PassesAll(x, cuts)) continue;
      if (std::abs(x.e - kKAlpha) <= 2. * SigmaKeV(kKAlpha)) tAlpha += x.w;
      if (std::abs(x.e - kKBeta) <= 2. * SigmaKeV(kKBeta)) tBeta += x.w;
      if (x.e >= gammaLow - 2. * SigmaKeV(gammaLow) && x.e <= std::min(kFitHigh, limit + 2. * SigmaKeV(limit))) tHigher += x.w;
    }
    const int nBinsK = expectK->FindBin(kFitHigh - 1e-6) - expectK->FindBin(kFitLow + 1e-6) + 1;
    const int ndfK = nBinsK - (backDegree + 1) - 2 * kPeaks, ndfKnownK = nBinsK - (backDegree + 1) - 2;
    double chi2McK = 0.;
    for (int b = expectK->FindBin(kFitLow + 1e-6); b <= expectK->FindBin(kFitHigh - 1e-6); ++b) {
      double var = 0.;
      for (int k = 0; k < 3; ++k) var += std::pow(afterFull[k]->GetBinError(b), 2);
      if (expectK->GetBinContent(b) > 0.) chi2McK += var / expectK->GetBinContent(b);
    }
    printf("\nthe K-line fit to the simulated spectrum after the cuts (%.0f-%.0f keV, %s; %d peaks of the resolution's width at free "
           "positions, background %s):\n", kFitLow, kFitHigh, inDays.c_str(), kPeaks, backName);
    printf("   truth after the cuts, the target's X-ray entries: within 2 sigma of K_alpha (%.1f keV) %.1f, of K_beta (%.1f keV) %.1f",
           kKAlpha, tAlpha, kKBeta, tBeta);
    if (gammaLow < kFitHigh) printf(", at %.1f-%.1f keV (the higher K lines, %.1f-%.1f keV, +-2 sigma) %.1f", gammaLow - 2. * SigmaKeV(gammaLow),
                                    std::min(kFitHigh, limit + 2. * SigmaKeV(limit)), gammaLow, limit, tHigher);
    printf("\n");
    for (std::size_t k = 0; k < kFound.size(); ++k) {
      const KPeak& p = kFound[k];
      const std::string what = kLineOf(p.mu) ? std::string(gTarget) + " " + kLineText[kLineOf(p.mu)] : std::string(kLineText[0]);
      double around[3] = {};   // the simulated entries within 1 resolution of it, by class (MC truth): what makes the peak
      for (int b = 1; b <= afterFull[0]->GetNbinsX(); ++b)
        if (std::abs(afterFull[0]->GetBinCenter(b) - p.mu) <= SigmaKeV(p.mu))
          for (int c = 0; c < 3; ++c) around[c] += afterFull[c]->GetBinContent(b);
      printf("   peak at %.2f +- %.2f keV (%s): %.1f +- %.1f entries, %.1f sigma; within 1 resolution of it as simulated: target X-rays "
             "%.0f, other muonic X-rays %.0f, the rest %.0f\n", p.mu, p.muError, what.c_str(), p.amount, p.error, p.z, around[0], around[1],
             around[2]);
      csv << "kfit," << k << ",\"peak at " << Form("%.2f", p.mu) << " keV (" << what << ")\",,,," << p.z << "," << p.amount << "," << p.error
          << "\n";
    }
    printf("   chi2/ndf %.1f/%d (Baker-Cousins, against the statistics of %.2f days); the simulation's own fluctuations alone make about "
           "%.0f of it\n", chi2Of(expectK, fK), ndfK, measureDays, chi2McK);
    printf("   null test, the same fit to the spectrum without any muonic X-ray: peaks at");
    for (std::size_t k = 0; k < kNullPeaks.size(); ++k) {
      const KPeak& p = kNullPeaks[k];
      printf("%s %.1f keV (%.1f +- %.1f entries, %.1f sigma)", k ? "," : "", p.mu, p.amount, p.error, p.z);
      csv << "kfit_null," << k << ",\"peak at " << Form("%.2f", p.mu) << " keV\",,,," << p.z << "," << p.amount << "," << p.error << "\n";
    }
    printf("\n   the known energies (only the amounts free): K_alpha %.1f +- %.1f entries (%.1f sigma), K_beta %.1f +- %.1f (%.1f sigma); "
           "chi2/ndf %.1f/%d\n", fKKnown->GetParameter(4), fKKnown->GetParError(4), kKnownZ[0], fKKnown->GetParameter(6),
           fKKnown->GetParError(6), kKnownZ[1], chi2Of(expectK, fKKnown), ndfKnownK);
    for (int k = 0; k < 2; ++k)
      csv << "kfit_known," << k << ",\"" << (k ? "K_beta" : "K_alpha") << " at " << (k ? kKBeta : kKAlpha) << " keV\",,,," << kKnownZ[k] << ","
          << fKKnown->GetParameter(4 + 2 * k) << "," << fKKnown->GetParError(4 + 2 * k) << "\n";
  }

  // ---- plots, one per file (CLAUDE.md rule 8)
  gStyle->SetOptStat(0);
  const int red = kRed + 1, blue = TColor::GetColor("#2171b5"), lightBlue = TColor::GetColor("#9ecae1");
  const int kindColour[3] = {red, lightBlue, blue};
  // the conditions in two lines at the foot of a plot (size: of the canvas height)
  auto footer = [&](TVirtualPad* c, double size, double gap) {
    c->cd();
    TLatex text;
    text.SetNDC(); text.SetTextFont(42); text.SetTextSize(size); text.SetTextColor(kGray + 2); text.SetTextAlign(11);
    text.DrawLatex(0.01, 0.006 + gap, stamp[0].c_str());
    text.DrawLatex(0.01, 0.006, stamp[1].c_str());
  };
  auto save = [&](TCanvas* c, const char* name, double footerSize = 0.02) {
    footer(c, footerSize, 1.3 * footerSize);
    c->Print(out + Form("/muonic_xray_select_%s.pdf", name));
    c->SaveAs(out + Form("/muonic_xray_select_%s.png", name));
    delete c;
  };
  auto canvas = [](const char* name, bool logy, int width = 900, int height = 700) {
    auto c = new TCanvas(name, name, width, height);
    c->SetLogy(logy);
    c->SetLeftMargin(0.12); c->SetRightMargin(0.03); c->SetTopMargin(0.09); c->SetBottomMargin(0.17);   // the footer under
    return c;
  };
  // the titles' and texts' sizes: of a 900 x 700 canvas, times textScale on a larger one (the spectra since 2026-11-05)
  auto heading = [](const char* text, double textScale = 1., double x = 0.12, double y = 0.935) {
    TLatex latex;
    latex.SetNDC(); latex.SetTextSize(0.045 * textScale);
    DrawFitted(latex, x, y, 0.98, text);
  };
  auto frame = [](double x1, double y1, double x2, double y2, const char* titles, double textScale = 1.) {
    auto f = gPad->DrawFrame(x1, y1, x2, y2, titles);
    f->GetXaxis()->SetTitleSize(0.045 * textScale); f->GetYaxis()->SetTitleSize(0.045 * textScale);
    if (textScale != 1.) {
      f->GetXaxis()->SetLabelSize(0.035 * textScale); f->GetYaxis()->SetLabelSize(0.035 * textScale);
      f->GetXaxis()->SetTickLength(0.02); f->GetYaxis()->SetTickLength(0.015);
    }
    f->GetXaxis()->SetTitleOffset(1.1); f->GetYaxis()->SetTitleOffset(textScale != 1. ? 1.0 : 1.25);
    return f;
  };
  auto legendAt = [](double x1, double y1, double x2, double y2, double size) {
    auto l = new TLegend(x1, y1, x2, y2);
    l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(size);
    return l;
  };
  // a curve drawn as separate dashes, dashPx long with gaps of gapPx along the curve as it shows on the pad (the PNG
  // output starts a line style's pattern anew at every point of a curve, so a dashed TF1 comes out solid)
  auto drawDashed = [](const TF1* curve, double dashPx, double gapPx, int col, int width) {
    const double xScale = gPad->GetWw() * gPad->GetAbsWNDC() / (gPad->GetX2() - gPad->GetX1());   // pixels per unit
    const double yScale = gPad->GetWh() * gPad->GetAbsHNDC() / (gPad->GetY2() - gPad->GetY1());
    // within the pad's range: a curve far outside it (a fit run away) would spoil the lengths and the drawing
    auto yPad = [](double y) {
      return std::min(std::max(gPad->GetLogy() ? std::log10(std::max(y, 1e-30)) : y, gPad->GetY1()), gPad->GetY2());
    };
    auto clip = [&](double y) { return gPad->GetLogy() ? std::pow(10., yPad(y)) : yPad(y); };
    TGraph* dash = nullptr;
    double length = 0., uLast = curve->GetXmin(), yLast = clip(curve->Eval(uLast));
    for (double u = curve->GetXmin(); u <= curve->GetXmax() + 1e-9; u += 0.01) {   // fine enough for dots on a steep line
      const double y = clip(curve->Eval(u));
      length += std::hypot((u - uLast) * xScale, (yPad(y) - yPad(yLast)) * yScale);
      uLast = u; yLast = y;
      if (std::fmod(length, dashPx + gapPx) < dashPx) {
        if (!dash) { dash = new TGraph(); dash->SetLineColor(col); dash->SetLineWidth(width); }
        dash->SetPoint(dash->GetN(), u, y);
      } else if (dash) {
        dash->Draw("L SAME");
        dash = nullptr;
      }
    }
    if (dash) dash->Draw("L SAME");
  };

  // 1. and the corner plots after the cuts: the variables vars for the L_alpha window entries that pass `keep`, with
  // the cuts drawn; signal as red points, background as a blue density (log), each variable's shapes on the diagonal
  {
    double stops[2] = {0., 1.}, r[2] = {0.9, 0.03}, g[2] = {0.95, 0.19}, b[2] = {1., 0.55};   // light to dark blue
    TColor::CreateGradientColorTable(2, stops, r, g, b, 64);
  }
  auto thousands = [](double v) {   // 318227 -> 318,227
    const std::string digits = Form("%.0f", v);
    std::string grouped;
    for (std::size_t i = 0; i < digits.size(); ++i) {
      if (i > 0 && (digits.size() - i) % 3 == 0) grouped += ',';
      grouped += digits[i];
    }
    return grouped;
  };
  // data: the records shown, CZT or stack entries or (byEvent) events; hidden: variables of the same group that are not
  // drawn, whose cuts (if chosen) are named under the plot
  auto corner = [&](const std::vector<Entry>& data, bool byEvent, const std::vector<int>& vars, std::function<bool(const Entry&)> keep,
                    const char* name, const char* title, const char* subtitle, const std::vector<int>& hidden) {
    const int nV = vars.size();
    std::vector<TH1D*> hS(nV), hB(nV), hL(nV);
    std::vector<std::vector<TH2D*>> hB2(nV, std::vector<TH2D*>(nV, nullptr));
    std::vector<std::vector<TGraph*>> gS(nV, std::vector<TGraph*>(nV, nullptr)), gL = gS;
    for (int a = 0; a < nV; ++a) {
      const std::vector<double> ea = Edges(kVars[vars[a]]);
      hS[a] = new TH1D(Form("%s_s%d", name, a), "", ea.size() - 1, ea.data());
      hB[a] = new TH1D(Form("%s_b%d", name, a), "", ea.size() - 1, ea.data());
      hL[a] = new TH1D(Form("%s_l%d", name, a), "", ea.size() - 1, ea.data());
      for (int c = 0; c < a; ++c) {
        const std::vector<double> ec = Edges(kVars[vars[c]]);
        hB2[a][c] = new TH2D(Form("%s_b%d%d", name, a, c), "", ec.size() - 1, ec.data(), ea.size() - 1, ea.data());
        gS[a][c] = new TGraph();
        gL[a][c] = new TGraph();
      }
    }
    // the pixel panel's own radioactivity (2026-11-05, asked for: where the LYSO's Lu-176 lies, to see how to cut it):
    // its entries in the coincidence window (one point per entry, its muon variables those of the event's muon; part of
    // B), violet under the signal's points, its shape (dashed) on the diagonals and its share in "cut keeps"
    double lAll = 0.;
    std::vector<double> lKept(nV, 0.);
    long nLu = 0;
    TRandom3 spreadLu(11);
    if (hasLu)
      for (const auto& x : entries) {
        if (!x.lu || !keep(x)) continue;
        lAll += x.w;
        ++nLu;
        for (int a = 0; a < nV; ++a) {
          const int i = vars[a];
          if (Passes(x, cuts, i)) lKept[a] += x.w;
          if (!Defined(x, i)) continue;
          const double va = ShownAt(x.v[i], kVars[i]);
          hL[a]->Fill(va, x.w);
          for (int c = 0; c < a; ++c) {
            const int j = vars[c];
            if (!Defined(x, j)) continue;
            gL[a][c]->SetPoint(gL[a][c]->GetN(), ShownAt(x.v[j], kVars[j]) + (kVars[j].integer ? spreadLu.Uniform(-0.3, 0.3) : 0.),
                               va + (kVars[i].integer ? spreadLu.Uniform(-0.3, 0.3) : 0.));
          }
        }
      }
    const int violet = TColor::GetColor("#9467bd");
    double sAll = 0., bAll = 0.;
    std::vector<double> sKept(nV, 0.), bKept(nV, 0.);
    long nSig = 0, nBack = 0;
    TRandom3 spread(7);
    for (const auto& x : data) {
      if (!keep(x)) continue;
      const bool sig = x.kind == 0;
      (sig ? sAll : bAll) += x.w;
      ++(sig ? nSig : nBack);
      for (int a = 0; a < nV; ++a) {
        const int i = vars[a];
        if (x.passed >= 0 ? (x.passed >> i & 1) : Passes(x, cuts, i)) (sig ? sKept : bKept)[a] += x.w;
        if (!Defined(x, i)) continue;
        const double va = ShownAt(x.v[i], kVars[i]);
        (sig ? hS : hB)[a]->Fill(va, x.w);
        for (int c = 0; c < a; ++c) {
          const int j = vars[c];
          if (!Defined(x, j)) continue;
          const double vc = ShownAt(x.v[j], kVars[j]);
          if (!sig) hB2[a][c]->Fill(vc, va, x.w);
          else gS[a][c]->SetPoint(gS[a][c]->GetN(), vc + (kVars[j].integer ? spread.Uniform(-0.3, 0.3) : 0.),
                                  va + (kVars[i].integer ? spread.Uniform(-0.3, 0.3) : 0.));
        }
      }
    }
    // where the cuts on a variable are, and what they reject
    auto bounds = [&](int i) {
      std::vector<double> v;
      if (cuts[i].active && cuts[i].lo > kVars[i].lo) v.push_back(cuts[i].lo);
      if (cuts[i].active && cuts[i].hi < kVars[i].hi) v.push_back(cuts[i].hi);
      return v;
    };
    const int width = std::max(1600, 250 * nV);   // cells of at least 250 px
    TCanvas c("corner", "corner plot", width, int(width * 1660. / 1600.));
    const double foot = 0.036, cell = 1. / nV, cellY = (1. - foot) / nV;   // square cells over the footer
    auto axisStyle = [](TAxis* axis, bool shown, const char* t) {
      axis->SetLabelSize(shown ? 0.075 : 0.); axis->SetTitleSize(shown ? 0.08 : 0.); axis->SetTitle(shown ? t : "");
      axis->SetTitleOffset(1.15); axis->SetNdivisions(505);
    };
    auto line = [](double x1, double y1, double x2, double y2) {
      auto l = new TLine(x1, y1, x2, y2);
      l->SetLineColor(kBlack); l->SetLineStyle(2); l->SetLineWidth(2);
      l->Draw();
      return l;
    };
    // a variable's split bin (e.g. beta >= 1) and its bin above over (e.g. >= 30 MeV): a grey line where each starts and
    // its name in it, the variable along x or y
    auto splitMark = [](int v, bool alongX, double top = 0.97) {   // top: where a name along x hangs from
      const Variable& var = kVars[v];
      if (var.split <= 0. && var.over <= 0.) return;
      gPad->Update();
      auto at = [](double f, bool x) {   // a fraction f of the frame along x or y, in the pad's coordinates
        const double u = x ? gPad->GetUxmin() + f * (gPad->GetUxmax() - gPad->GetUxmin())
                           : gPad->GetUymin() + f * (gPad->GetUymax() - gPad->GetUymin());
        return (x ? gPad->GetLogx() : gPad->GetLogy()) ? std::pow(10., u) : u;
      };
      auto mark = [&](double edge, double centre, const char* label) {
        auto l = alongX ? new TLine(edge, at(0., false), edge, at(1., false)) : new TLine(at(0., true), edge, at(1., true), edge);
        l->SetLineColor(kGray + 2); l->SetLineStyle(2);
        l->Draw();
        TLatex text;
        text.SetTextFont(42); text.SetTextSize(0.06); text.SetTextColor(kGray + 2); text.SetTextAlign(32);
        text.SetTextAngle(alongX ? 90. : 0.);   // along x: hanging from the frame's top
        if (alongX) text.DrawLatex(centre, at(top, false), label);
        else text.DrawLatex(at(0.97, true), centre, label);
      };
      if (var.split > 0.) mark(var.split, ShownAt(var.splitLow ? var.lo : var.hi, var), var.splitLabel);
      if (var.over > 0.) mark(var.over, ShownAt(var.hi, var), var.overLabel);
    };
    TLine* cutLine = nullptr;
    for (int row = 0; row < nV; ++row)
      for (int col = 0; col <= row; ++col) {
        c.cd();
        auto pad = new TPad(Form("pad%d%d", row, col), "", col * cell, 1. - (row + 1) * cellY, (col + 1) * cell, 1. - row * cellY);
        pad->SetLeftMargin(0.2); pad->SetBottomMargin(0.2); pad->SetRightMargin(0.03); pad->SetTopMargin(0.03);
        pad->Draw();
        pad->cd();
        const int i = vars[row], j = vars[col];
        pad->SetLogx(kVars[j].log);
        if (row != col) pad->SetLogy(kVars[i].log);
        const bool bottomRow = row == nV - 1;
        if (row == col) {   // the shapes, unit area
          TH1D *hb = hB[row], *hs = hS[row], *hl = hL[row];
          for (auto h : {hb, hs, hl}) if (h->Integral() > 0.) h->Scale(1. / h->Integral());
          const double top = std::max({hb->GetMaximum(), hs->GetMaximum(), hl->GetMaximum()});
          // on a log scale (logShape) four decades under the highest bin, else from 0
          const double floor = kVars[i].logShape ? 1e-4 * top : 0.;
          pad->SetLogy(kVars[i].logShape);
          hb->SetMinimum(floor); hb->SetMaximum(kVars[i].logShape ? 50. * top : 1.45 * top);
          hb->SetFillColor(lightBlue); hb->SetLineColor(blue);
          hs->SetLineColor(red); hs->SetLineWidth(2);
          axisStyle(hb->GetXaxis(), bottomRow, kVars[i].name);
          axisStyle(hb->GetYaxis(), false, "");
          if (i == 0 && bottomRow) { hb->GetXaxis()->SetBinLabel(1, "none"); hb->GetXaxis()->SetBinLabel(2, "hit"); }
          if (i == 13 && kVars[i].log) { hb->GetXaxis()->SetMoreLogLabels(); hb->GetXaxis()->SetNoExponent(); }
          hb->Draw("HIST");
          if (cuts[i].active) {
            for (auto [x1, x2] : {std::make_pair(kVars[i].lo, std::max(kVars[i].lo, cuts[i].lo)),
                                  std::make_pair(std::min(kVars[i].hi, cuts[i].hi), kVars[i].hi)}) {
              if (x2 <= x1) continue;
              auto box = new TBox(x1, floor, x2, top);   // shaded, the shapes still seen through it
              box->SetFillStyle(1001); box->SetFillColorAlpha(kGray + 2, 0.4); box->SetLineWidth(0);
              box->Draw();
            }
          }
          hs->Draw("HIST SAME");
          if (nLu > 0) { hl->SetLineColor(violet); hl->SetLineWidth(2); hl->SetLineStyle(2); hl->Draw("HIST SAME"); }
          for (double x : bounds(i)) cutLine = line(x, floor, x, top);
          splitMark(i, true, 0.82);   // under the "cut keeps" line
          TLatex text;
          text.SetNDC(); text.SetTextSize(0.068);
          const std::string luKept = nLu > 0 ? Form(", Lu %.0f%%", 100. * lKept[row] / std::max(1e-12, lAll)) : "";
          DrawFitted(text, 0.23, 0.89, 0.98, cuts[i].active ? Form("cut keeps S %.0f%%, B %.0f%%%s", 100. * sKept[row] / std::max(1e-12, sAll),
                                                                   100. * bKept[row] / std::max(1e-12, bAll), luKept.c_str())
                                                            : "no cut");
          gPad->RedrawAxis();
        } else {            // background density and signal points: x the column's variable, y the row's
          TH2D* h = hB2[row][col];
          pad->SetLogz();
          axisStyle(h->GetXaxis(), bottomRow, kVars[j].name);
          axisStyle(h->GetYaxis(), col == 0, kVars[i].name);
          if (j == 0 && bottomRow) { h->GetXaxis()->SetBinLabel(1, "none"); h->GetXaxis()->SetBinLabel(2, "hit"); }
          if (j == 13 && kVars[j].log) { h->GetXaxis()->SetMoreLogLabels(); h->GetXaxis()->SetNoExponent(); }
          h->SetMinimum(1e-3 * h->GetMaximum());   // three decades of background density
          h->Draw("COL");
          // what a cut rejects shaded as on the diagonal (asked for on 2026-11-02): beyond a cut on the column's
          // variable over the whole height, beyond one on the row's over the whole width (darker where both reject),
          // semi-transparent so that the density stays visible; the signal drawn over it
          auto shade = [](double x1, double y1, double x2, double y2) {
            auto box = new TBox(x1, y1, x2, y2);
            box->SetFillStyle(1001); box->SetFillColorAlpha(kGray + 2, 0.4); box->SetLineWidth(0);
            box->Draw();
          };
          const double xLo = kVars[j].lo, xHi = kVars[j].hi, yLo = kVars[i].lo, yHi = kVars[i].hi;
          if (cuts[j].active && cuts[j].lo > xLo) shade(xLo, yLo, std::min(cuts[j].lo, xHi), yHi);
          if (cuts[j].active && cuts[j].hi < xHi) shade(std::max(cuts[j].hi, xLo), yLo, xHi, yHi);
          if (cuts[i].active && cuts[i].lo > yLo) shade(xLo, yLo, xHi, std::min(cuts[i].lo, yHi));
          if (cuts[i].active && cuts[i].hi < yHi) shade(xLo, std::max(cuts[i].hi, yLo), xHi, yHi);
          // the panel's own decays under the signal, as medium dots (style 7: ROOT's PNG painter draws them as pixels, 0.3 s
          // for 40 000 where filled circles take 170 s)
          gL[row][col]->SetMarkerStyle(7); gL[row][col]->SetMarkerColor(violet);
          if (gL[row][col]->GetN() > 0) gL[row][col]->Draw("P SAME");
          gS[row][col]->SetMarkerStyle(20); gS[row][col]->SetMarkerSize(0.45); gS[row][col]->SetMarkerColor(red);
          if (gS[row][col]->GetN() > 0) gS[row][col]->Draw("P SAME");
          for (double x : bounds(j)) line(x, kVars[i].lo, x, kVars[i].hi);
          for (double y : bounds(i)) line(kVars[j].lo, y, kVars[j].hi, y);
          splitMark(j, true);
          splitMark(i, false);
          gPad->RedrawAxis();
        }
      }
    // the free upper right: what is shown, the cuts, the numbers
    c.cd();
    const double infoCells = nV >= 7 ? 3. : nV >= 5 ? 2. : 1.;   // the free cells of the upper right (text scales with its height)
    auto info = new TPad("info", "", infoCells * cell + 0.03, 1. - infoCells * cellY + 0.03, 0.99, 0.99);
    info->Draw();
    info->cd();
    TLatex text;
    text.SetNDC(); text.SetTextFont(62); text.SetTextSize(0.075);
    DrawFitted(text, 0.02, 0.9, 0.99, title);
    text.SetTextFont(42); text.SetTextSize(0.052);
    text.DrawLatex(0.02, 0.81, subtitle);
    std::vector<std::string> drawn = {"cuts drawn:"};   // wrapped at the commas
    for (int i : vars) {
      if (!cuts[i].active) continue;
      const std::string next = Describe(cuts, i);
      if (drawn.back().size() + next.size() > 58) drawn.push_back(next);
      else drawn.back() += (drawn.back().back() == ':' ? " " : ", ") + next;
    }
    if (drawn.size() == 1 && drawn[0] == "cuts drawn:") drawn[0] += " none on these variables";
    for (int i : hidden)
      if (cuts[i].active) drawn.push_back("not drawn: " + Describe(cuts, i));
    for (int i : vars)   // a diagonal on a log scale, under the cuts
      if (kVars[i].logShape && drawn.size() < 6) drawn.push_back(std::string(kVars[i].name) + ": its diagonal on a log scale");
    const bool many = drawn.size() > 4;   // up to 6 lines under the numbers: the legend and the numbers higher up
    auto marker = static_cast<TGraph*>(gS[1][0]->Clone("marker"));   // a legend marker big enough to see
    marker->SetMarkerSize(1.4);
    auto area = static_cast<TH1D*>(hB[0]->Clone("area"));
    auto legend = new TLegend(0.02, many ? 0.56 : 0.46, 0.9, many ? 0.79 : 0.77);
    legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(0.052);
    const char* unit = byEvent ? "events" : "entries";
    legend->AddEntry(marker, Form("%s with a muonic X-ray of the target (S)", unit), "p");
    legend->AddEntry(area, Form("other %s (B)", unit), "f");
    if (nLu > 0) {   // a legend marker big enough to see, not transparent
      auto luMarker = new TGraph(); luMarker->SetMarkerStyle(20); luMarker->SetMarkerSize(1.4); luMarker->SetMarkerColor(violet);
      luMarker->SetLineColor(violet); luMarker->SetLineStyle(2); luMarker->SetLineWidth(2);
      legend->AddEntry(luMarker, "entries of the panel's Lu-176 in time (part of B)", "pl");
    }
    if (cutLine) legend->AddEntry(cutLine, cutsName, "l");
    legend->Draw();
    // the simulated events themselves, of all the simulated days (asked for on 2026-11-02: no 21-day equivalents here;
    // event_flow.txt has them in brackets)
    const std::string luSimulated = nLu > 0 ? Form("; Lu-176 %s entries", thousands(nLu).c_str()) : "";
    text.DrawLatex(0.02, many ? 0.49 : 0.38, Form("simulated: S %s, B %s %s%s", thousands(nSig).c_str(), thousands(nBack).c_str(), unit,
                                                  luSimulated.c_str()));
    text.SetTextSize(0.042);
    text.DrawLatex(0.02, many ? 0.435 : 0.325, Form("all simulated days: %.2f of muons below %.0f GeV, %.2f above", thickDays,
                                                    eSplit / 1000., refDays));
    text.SetTextSize(many ? 0.042 : 0.045);
    for (std::size_t k = 0; k < drawn.size() && k < 6; ++k)
      text.DrawLatex(0.02, (many ? 0.37 : 0.27) - (many ? 0.05 : 0.055) * k, drawn[k].c_str());
    text.DrawLatex(0.02, 0.05, "shaded grey = rejected by a cut; diagonal: unit area; B log density below it");
    const TString file = out + Form("/muonic_xray_select_%s", name);
    footer(&c, 0.0105, 0.0135);
    c.Print(file + ".pdf");
    c.SaveAs(file + ".png");
  };
  // the corner plots' variables: those this kind of skim has, related ones next to each other
  auto usedOf = [](std::vector<int> v) {
    v.erase(std::remove_if(v.begin(), v.end(), [](int i) { return !Used(i); }), v.end());
    return v;
  };
  // the entry's variables (for its representative entry), and the one corner plot's (since 2026-11-02, asked for that day:
  // the muon's and the entry's together, in the order asked for; without the bottom-TOF hit and bar, 0 and 10, the CZT
  // pixel and the pixel centre, 12 and 9, whose cuts, if chosen, are named under the plot as not drawn)
  const std::vector<int> stackVars = usedOf(gCZT ? std::vector<int>{13, 5, 6, 9, 12, 7, 8} : std::vector<int>{13, 5, 6, 7, 8, 9});
  const std::vector<int> cornerVars = usedOf(gCZT && gPixel ? std::vector<int>{1, 11, 4, 2, 3, 13, 6, 5}
                                                             : std::vector<int>{1, 11, 2, 3, 4, 13, 5, 6, 7, 8, 9});
  std::vector<int> cornerHidden;
  for (int i : usedOf({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13}))
    if (std::find(cornerVars.begin(), cornerVars.end(), i) == cornerVars.end()) cornerHidden.push_back(i);
  auto passesGroup = [&](const Entry& x, const std::vector<int>& group) {
    for (int i : group) if (!Passes(x, cuts, i)) return false;
    return true;
  };
  // the corner plots show the captured data, as event_flow.C counts it (asked for on 2026-11-02): the TOF-triggered
  // events with a CZT (stack) entry in time, at every energy (step 2a of event_flow.C), the cuts as chosen in the L_alpha
  // window drawn, each cut's share by itself (event_flow_cuts.png); one point per event, as event_flow.C counts them:
  // signal an event with a muonic X-ray of the target among its entries; the entry's variables those of its
  // representative entry, a muonic X-ray of the target first, then one that passes the entry's cuts (so that the event
  // passes them if any entry does), then one in the L_alpha window; else the first; each cut's "cut keeps" counts an
  // event that any of its entries passes, as event_flow.C. Drawn in one plot, corner_after
  // (since 2026-11-02: no plots after the other group's cuts, which no step of the flow shows)
  auto eventsOf = [&](std::function<bool(const Entry&)> keep) {
    std::vector<Entry> events;
    std::vector<int> passed;   // per event, bit i: one of its entries passes cut i by itself
    std::map<long long, std::size_t> where;
    auto better = [&](const Entry& a, const Entry& b) {
      if (a.kind != b.kind) return a.kind < b.kind;
      const bool passA = passesGroup(a, stackVars), passB = passesGroup(b, stackVars);
      if (passA != passB) return passA;
      return InAlpha(a.e) && !InAlpha(b.e);
    };
    for (const auto& x : entries) {
      if (!keep(x)) continue;
      auto [it, isNew] = where.emplace(x.id, events.size());
      if (isNew) { events.push_back(x); passed.push_back(0); }
      else if (better(x, events[it->second])) events[it->second] = x;
      for (int i = 0; i < kNV; ++i) if (Passes(x, cuts, i)) passed[it->second] |= 1 << i;
    }
    for (std::size_t k = 0; k < events.size(); ++k) events[k].passed = passed[k];
    return events;
  };
  const std::vector<Entry> captured = eventsOf([](const Entry&) { return true; });
  const char* anEntry = gOtherPanel ? "a pixel entry" : gCZT ? "a CZT entry" : "a stack entry";
  const char* eventsName = gTOF ? "triggered events" : "events with a middle-TOF hit";
  const std::string capturedAll = Form("%s with %s, all energies (step 2a; one point per event)", eventsName, anEntry);
  corner(captured, true, cornerVars, [](const Entry&) { return true; }, "corner_after",
         Form("The muon and the %s, %s", entryName, label), capturedAll.c_str(), cornerHidden);

  // 2. the cut flow
  {
    auto c = canvas("cutflow", true);
    c->SetBottomMargin(0.42);
    const int nSteps = flowS.size();
    auto axis = new TH1D("steps", Form(";;entries in the %s L-alpha window, %.2f days", gTarget, measureDays), nSteps, -0.5, nSteps - 0.5);
    for (int k = 0; k < nSteps; ++k) axis->GetXaxis()->SetBinLabel(k + 1, (k ? "+ " + flowNames[k] : flowNames[k]).c_str());
    axis->GetXaxis()->LabelsOption("v");
    axis->GetXaxis()->SetLabelSize(0.032);
    axis->GetYaxis()->SetTitleSize(0.045); axis->GetYaxis()->SetTitleOffset(1.25);
    // the range of both curves (the background can end below the signal), its smallest value above zero
    double low = 1e30, high = 0.;
    for (const auto* flow : {&flowS, &flowB})
      for (double v : *flow) { if (v > 0.) low = std::min(low, v); high = std::max(high, v); }
    axis->SetMinimum(0.3 * low);
    axis->SetMaximum(40. * high);
    axis->Draw("AXIS");
    auto gSig = new TGraph(), gBack = new TGraph();
    for (int k = 0; k < nSteps; ++k) { gSig->SetPoint(k, k, flowS[k]); gBack->SetPoint(k, k, flowB[k]); }
    gSig->SetMarkerStyle(20); gSig->SetMarkerSize(1.3); gSig->SetMarkerColor(red); gSig->SetLineColor(red); gSig->SetLineWidth(2);
    gBack->SetMarkerStyle(21); gBack->SetMarkerSize(1.3); gBack->SetMarkerColor(blue); gBack->SetLineColor(blue); gBack->SetLineWidth(2);
    gBack->Draw("PL SAME");
    gSig->Draw("PL SAME");
    TLatex numbers;   // not `label`: the study's label goes into the heading below
    numbers.SetTextSize(0.032); numbers.SetTextAlign(21);
    for (int k = 0; k < nSteps; ++k) numbers.DrawLatex(k, 1.8 * std::max(flowS[k], flowB[k]), Form("%.1f", flowZ[k]));
    auto legend = legendAt(0.14, 0.83, 0.97, 0.9, 0.034);
    legend->SetNColumns(2);
    legend->AddEntry(gSig, "muonic X-rays of the target (S)", "lp");
    legend->AddEntry(gBack, "everything else (B)", "lp");
    legend->Draw();
    TLatex note;
    note.SetNDC(); note.SetTextSize(0.032); note.SetTextAlign(31);
    note.DrawLatex(0.96, 0.8, Form("numbers: S/sqrt(B) %s after the step", inDays.c_str()));
    heading(Form("Cut flow of %s, %s", handChosen ? "the analysis cuts" : "the most efficient cuts", label));
    save(c, "cutflow");
  }

  // 3. the spectra by what made them: the rest (blue) at the bottom, other muonic X-rays (light blue), the target's red
  const char* shortNames[3] = {"target muonic X-rays", "other muonic X-rays", "the rest"};   // for the legends right of the lines
  auto drawStack = [&](const std::array<TH1D*, 3>& h, double scale, TLegend* legend, const char* tag) {
    auto all = static_cast<TH1D*>(h[2]->Clone(Form("%s_all", tag)));
    all->Add(h[1]); all->Add(h[0]);
    auto lower = static_cast<TH1D*>(h[2]->Clone(Form("%s_lower", tag)));
    lower->Add(h[1]);
    auto rest = static_cast<TH1D*>(h[2]->Clone(Form("%s_rest", tag)));
    TH1D* layer[3] = {all, lower, rest};   // drawn from the top down, each in the colour of its top component
    for (int k = 0; k < 3; ++k) {
      layer[k]->Scale(scale);
      layer[k]->SetFillColor(kindColour[k]); layer[k]->SetLineColor(kindColour[k]); layer[k]->SetFillStyle(1001);
      layer[k]->Draw("HIST SAME");
      if (legend) legend->AddEntry(layer[k], Form("%s: %.0f entries", shortNames[k], h[k]->GetEntries()), "f");
    }
    return all;
  };
  // the target's lines as simulated (the cascade model's energies; the measured Al K_alpha is 346.8 keV), named dx right
  // of them: L_alpha and L_beta, with kLines the K lines too (L_alpha and L_beta then named together)
  auto lineMarks = [&](double yLow, double yHigh, bool logy, bool kLines, double dx, double textScale = 1.) {
    TLatex mark;
    mark.SetTextSize(0.03 * textScale); mark.SetTextColor(kGray + 2); mark.SetTextAngle(90); mark.SetTextAlign(12);
    const std::string el = gTarget;
    std::vector<std::pair<double, std::string>> marks = {{kAlpha, el + (kLines ? " L_{#alpha}, L_{#beta}" : " L_{#alpha}")},
                                                         {kBeta, kLines ? "" : el + " L_{#beta}"}};
    if (kLines) { marks.push_back({kKAlpha, el + " K_{#alpha}"}); marks.push_back({kKBeta, el + " K_{#beta}"}); }
    for (const auto& [eKeV, nameLine] : marks) {
      auto l = new TLine(eKeV, yLow, eKeV, yHigh);
      l->SetLineColor(kGray + 2); l->SetLineStyle(3);
      l->Draw();
      if (!nameLine.empty()) mark.DrawLatex(eKeV + dx, logy ? 1.3 * yHigh : yHigh, nameLine.c_str());
    }
  };
  // what a step holds, as event_flow.C counts it: each class's entries at eMin-eMax keV (the bins whose centre lies
  // there), in measureDays days and simulated, and the step's events at every energy (signal: a muonic X-ray of the
  // target among the entries that pass)
  struct StepCount { double total[3]; long raw[3]; double sEvents, bEvents; long sRaw, bRaw; };
  auto stepCount = [&](const std::array<TH1D*, 3>& h, const Cuts& cut, double eMin = 0., double eMax = 1e9) {
    StepCount n{};
    for (int k = 0; k < 3; ++k) {
      double sum = 0.;
      for (int b = 1; b <= h[k]->GetNbinsX(); ++b)
        if (h[k]->GetBinCenter(b) >= eMin && h[k]->GetBinCenter(b) < eMax) sum += h[k]->GetBinContent(b);
      n.total[k] = sum;
    }
    std::map<long long, std::pair<double, bool>> events;   // the event's weight, signal
    for (const auto& x : entries) {
      if (!PassesAll(x, cut)) continue;
      if (x.e >= eMin && x.e < eMax) ++n.raw[x.kind];
      auto& event = events[x.id];
      event.first = x.w;
      event.second = event.second || x.kind == 0;
    }
    for (const auto& [id, event] : events) {
      if (event.second) { n.sEvents += event.first; ++n.sRaw; }
      else { n.bEvents += event.first; ++n.bRaw; }
    }
    return n;
  };
  auto count = [&](double v) { return v < 100. ? std::string(Form("%.1f", v)) : thousands(v); };
  // the spectrum of the entries that pass the cut set c, per 2 keV in measureDays days at plotLow-plotHigh, as ROOT draws a histogram
  // (asked for on 2026-11-02): the simulated total as a line with its statistical error bars (the simulation's own), each
  // class as a thin line in its colour; with the step's events (at every energy) and each class's entries in the range
  // drawn, in measureDays days as event_flow.C counts them; mark: the L-line fit's and (2026-11-05) the K-line fit's
  // ranges, curves and results; perDay: the
  // histograms drawn per day (asked for on 2026-11-03 for spectrum_after; the legend's numbers stay those of measureDays days)
  auto spectrumPlot = [&](const std::array<TH1D*, 3>& h, const Cuts& cut, bool logy, bool mark, const char* name, const char* title,
                          bool perDay = false) {
    // 1600 x 1150 since 2026-11-05 (asked for: the plots got crowded), its texts as large as on the 900 x 700 canvas before
    const double ts = 0.62;
    auto c = canvas(name, logy, 1600, 1150);
    c->SetLeftMargin(0.075); c->SetRightMargin(0.02); c->SetTopMargin(0.065); c->SetBottomMargin(0.105);
    const double scale = perDay ? 1. / measureDays : 1.;
    auto total = static_cast<TH1D*>(h[2]->Clone(Form("%s_total", name)));
    total->Add(h[1]); total->Add(h[0]);
    total->Scale(scale);
    double top = 0., low = 1e30;
    for (int b = 1; b <= total->GetNbinsX(); ++b) {
      if (total->GetBinCenter(b) < plotLow || total->GetBinCenter(b) >= plotHigh) continue;
      top = std::max(top, total->GetBinContent(b) + total->GetBinError(b));
      for (int k = 0; k < 3; ++k) if (h[k]->GetBinContent(b) > 0.) low = std::min(low, scale * h[k]->GetBinContent(b));
    }
    const StepCount n = stepCount(h, cut, plotLow, plotHigh);
    const double headroom = logy ? 1000. : 1.6;   // the frame's top over the highest bar: room for the legend
    frame(plotLow, logy ? 0.5 * low : 0., plotHigh, headroom * top,
          perDay ? Form(";%s;entries per %.0f keV per day", energyTitle, total->GetBinWidth(1))
                 : Form(";%s;entries per %.0f keV in %.2f days", energyTitle, total->GetBinWidth(1), measureDays), ts);
    TBox *band = nullptr, *bandK = nullptr;
    if (mark) {   // the fits' ranges, behind everything: the L-line fit's, the K-line fit's lighter (2026-11-05)
      band = new TBox(fitLow, 0., fitHigh, headroom * top);
      band->SetFillColor(TColor::GetColor("#ececec")); band->SetLineColor(TColor::GetColor("#ececec"));
      band->Draw();
      bandK = new TBox(kFitLow, 0., kFitHigh, headroom * top);
      bandK->SetFillColor(TColor::GetColor("#f6f6f6")); bandK->SetLineColor(TColor::GetColor("#f6f6f6"));
      bandK->Draw();
    }
    TH1D* line[3];
    for (int k = 2; k >= 0; --k) {
      line[k] = static_cast<TH1D*>(h[k]->Clone(Form("%s_line%d", name, k)));
      line[k]->Scale(scale);
      line[k]->SetLineColor(kindColour[k]); line[k]->SetLineWidth(1);
      line[k]->Draw("HIST SAME");
    }
    // the pixel panel's own radioactivity (2026-11-05, asked for: how the LYSO's Lu-176 compares): its entries in the
    // coincidence window (part of the rest) in violet, and the sidebands' estimate of them dashed
    TH1D *lu = nullptr, *luSide = nullptr;
    double luSum[2] = {}, luRaw[2] = {};   // counted, from the sidebands: in measureDays days at plotLow-plotHigh, simulated
    if (hasLu) {
      const int nb = h[0]->GetNbinsX();
      const double x1 = h[0]->GetXaxis()->GetXmin(), x2 = h[0]->GetXaxis()->GetXmax();
      lu = new TH1D(Form("%s_lu", name), "", nb, x1, x2);
      luSide = new TH1D(Form("%s_luSide", name), "", nb, x1, x2);
      for (int k = 0; k < 2; ++k)
        for (const auto& x : k ? sideLu : entries) {
          if ((!k && !x.lu) || !PassesAll(x, cut)) continue;
          (k ? luSide : lu)->Fill(x.e, x.w);
          const int b = h[0]->FindBin(x.e);
          if (h[0]->GetBinCenter(b) >= plotLow && h[0]->GetBinCenter(b) < plotHigh) { luSum[k] += x.w; ++luRaw[k]; }
        }
      const int violet = TColor::GetColor("#9467bd");
      for (TH1D* l : {lu, luSide}) { l->Scale(scale); l->SetLineColor(violet); l->SetLineWidth(2); }
      luSide->SetLineStyle(2);
      lu->Draw("HIST SAME");
      luSide->Draw("HIST SAME");
    }
    total->SetLineColor(kBlack); total->SetLineWidth(1); total->SetMarkerSize(0);   // thin: the classes show under it
    total->Draw("HIST SAME");
    total->Draw("E1 SAME");
    // the L-line fit over its range (2026-11-05, asked for: with a coarse resolution its background shows better here
    // than in measurement_fit): the fit (fBins as last fitted: to this spectrum after the cuts) and its smooth background
    TF1 *fitAll = nullptr, *fitBack = nullptr;
    if (mark) {
      fitAll = new TF1(Form("%s_fit", name), [=](double* x, double*) { return scale * fBins->Eval(x[0]); }, fitLow, fitHigh, 0);
      fitBack = new TF1(Form("%s_fitBackground", name), [=](double* x, double*) {
        const double u = (x[0] - 80.) / 50.;
        return scale * std::exp(fBins->GetParameter(0) + fBins->GetParameter(1) * u + fBins->GetParameter(2) * u * u +
                                fBins->GetParameter(3) * u * u * u);
      }, fitLow, fitHigh, 0);
      drawDashed(fitBack, 9., 6., kRed, 2);
      fitAll->SetLineColor(kRed); fitAll->SetLineWidth(2); fitAll->SetNpx(1000);
      fitAll->Draw("SAME");
      // the K-line fit over its range, alike (2026-11-05, asked for: the fitting all the way to 500 keV)
      auto kAll = new TF1(Form("%s_kfit", name), [=](double* x, double*) { return scale * fK->Eval(x[0]); }, kFitLow, kFitHigh, 0);
      auto kBack = new TF1(Form("%s_kfitBackground", name), [=](double* x, double*) { return scale * kBackgroundOf(fK, x[0]); },
                           kFitLow, kFitHigh, 0);
      drawDashed(kBack, 9., 6., kRed, 2);
      kAll->SetLineColor(kRed); kAll->SetLineWidth(2); kAll->SetNpx(1000);
      kAll->Draw("SAME");
    }
    lineMarks(logy ? 0.5 * low : 0., logy ? 3. * top : 1.15 * top, logy, true, 0.006 * (plotHigh - plotLow), ts);
    // between the L lines and the K lines, as tall as its rows
    const int rows = 7 + (hasLu ? 2 : 0) + (mark ? 5 + int(kFound.size()) : 0);
    auto legend = legendAt(0.25, 0.92 - rows * 0.026 * ts * 1.55, 0.74, 0.92, 0.026 * ts);
    legend->SetHeader(Form("%s:", inDays.c_str()));
    legend->AddEntry((TObject*)nullptr, Form("events: S %s, B %s (simulated %s and %s)", count(n.sEvents).c_str(), count(n.bEvents).c_str(),
                                             thousands(n.sRaw).c_str(), thousands(n.bRaw).c_str()), "");
    legend->AddEntry((TObject*)nullptr, Form("entries at %.0f-%.0f keV (simulated):", plotLow, plotHigh), "");
    auto numbers = [&](double v, long r) { return count(v) + " (" + thousands(r) + ")"; };
    legend->AddEntry(total, Form("total: %s", numbers(n.total[0] + n.total[1] + n.total[2], n.raw[0] + n.raw[1] + n.raw[2]).c_str()), "le");
    for (int k = 0; k < 3; ++k) legend->AddEntry(line[k], Form("%s: %s", shortNames[k], numbers(n.total[k], n.raw[k]).c_str()), "l");
    if (hasLu) {
      legend->AddEntry(lu, Form("   of it the panel's Lu-176: %s", numbers(luSum[0], long(luRaw[0])).c_str()), "l");
      legend->AddEntry(luSide, Form("   Lu-176 from the sidebands x %.4f: %s", kSideScale, numbers(luSum[1], long(luRaw[1])).c_str()), "l");
    }
    if (mark) {
      legend->AddEntry(band, Form("L-line fit from the bins alone, %.0f-%.0f keV:", fitLow, fitHigh), "f");
      legend->AddEntry((TObject*)nullptr, Form("   peak at %.2f keV: %.0f +/- %.0f entries, %.1f sigma", fitShown.position, fitShown.amount,
                                               fitShown.error, fitShown.z), "");
      legend->AddEntry(bandK, Form("K-line fit from the bins alone, %.0f-%.0f keV:", kFitLow, kFitHigh), "f");
      for (const auto& p : kFound) {   // each peak named after the target's line within 2 resolutions of it
        const int line = kLineOf(p.mu);
        const std::string what = line ? std::string(gTarget) + " " + kLineLatex[line] : std::string(kLineLatex[0]);
        legend->AddEntry((TObject*)nullptr, Form("   peak at %.1f keV (%s): %.0f +/- %.0f entries, %.1f sigma", p.mu, what.c_str(),
                                                 p.amount, p.error, p.z), "");
      }
      legend->AddEntry(fitAll, kWidthMax > kWidthMin ? Form("   the fits: peaks, the L fit's widths #geq resolution (L 2, K %d)", kPeaks)
                                                     : Form("   the fits: peaks of the resolution's width (L 2, K %d)", kPeaks), "l");
      auto dash = new TLine(); dash->SetLineColor(kRed); dash->SetLineStyle(2); dash->SetLineWidth(2);
      legend->AddEntry(dash, Form("   their smooth backgrounds, %s each", backName), "l");
    }
    legend->Draw();
    gPad->RedrawAxis();
    heading(title, ts, 0.075, 0.952);
    save(c, name, 0.02 * ts);
  };
  const auto beforeFull = spectrum(Cuts{}, "beforeFull", 500, 1000.);   // afterFull: made for the K-line fit
  spectrumPlot(beforeFull, Cuts{}, true, false, "spectrum_before",
               Form("%s, before the cuts: %s", gTOF ? "Captured" : "Middle-TOF hit", label));
  spectrumPlot(afterFull, cuts, false, true, "spectrum_after", Form("After the corner-plot cuts: %s", label), true);
  {
    const StepCount all = stepCount(afterFull, cuts), drawn = stepCount(afterFull, cuts, plotLow, plotHigh);
    printf("\nafter the cuts, %.2f days: events S %.1f, B %.1f (simulated %ld, %ld)\n"
           "   entries at every energy: target muonic X-rays %.1f, other muonic X-rays %.1f, the rest %.1f (simulated %ld, %ld, %ld)\n"
           "   of these at %.0f-%.0f keV (drawn): %.1f, %.1f, %.1f (simulated %ld, %ld, %ld)\n", measureDays, all.sEvents, all.bEvents,
           all.sRaw, all.bRaw, all.total[0], all.total[1], all.total[2], all.raw[0], all.raw[1], all.raw[2], plotLow, plotHigh,
           drawn.total[0], drawn.total[1], drawn.total[2], drawn.raw[0], drawn.raw[1], drawn.raw[2]);
  }

  // the L-line fit of the simulated spectrum after the cuts (step 3 of event_flow.C) over the fit's range, as ROOT draws
  // a fitted histogram (2026-11-02): the simulated total per day and 2 keV (asked for on 2026-11-03: this plot alone per
  // day; the fit and its numbers stay those of measureDays days) as a line with its statistical
  // error bars (the simulation's), the classes as thin lines, the fit as a red curve with its smooth background dashed, a
  // box with the fit's numbers (bigger than the other plots); in a panel below, the simulated spectrum minus the fitted
  // smooth background, with the fitted lines: the excess that the significance measures
  {
    const StepCount n = stepCount(after, cuts, fitLow, fitHigh);
    const int bLow = after[0]->FindBin(fitLow + 1e-6), bHigh = after[0]->FindBin(fitHigh - 1e-6);
    auto total = static_cast<TH1D*>(after[2]->Clone("fitTotal"));   // per day and 2 keV, with the simulation's errors
    total->Add(after[1]); total->Add(after[0]);
    total->Scale(1. / measureDays);
    auto c = new TCanvas("measurement_fit", "measurement_fit", 1600, 1600);
    auto upper = new TPad("upper", "", 0., 0.34, 1., 1.);
    upper->SetLeftMargin(0.12); upper->SetRightMargin(0.03); upper->SetTopMargin(0.12); upper->SetBottomMargin(0.015);
    auto lower = new TPad("lower", "", 0., 0., 1., 0.34);
    lower->SetLeftMargin(0.12); lower->SetRightMargin(0.03); lower->SetTopMargin(0.03); lower->SetBottomMargin(0.4);
    upper->Draw();
    lower->Draw();
    // the upper panel: the simulated spectrum, its classes, the fit
    upper->cd();
    double top = 0.;
    for (int b = bLow; b <= bHigh; ++b) top = std::max(top, total->GetBinContent(b) + total->GetBinError(b));
    auto background = new TF1("background", [=](double* x, double*) {
      const double u = (x[0] - 80.) / 50.;
      return std::exp(fBins->GetParameter(0) + fBins->GetParameter(1) * u + fBins->GetParameter(2) * u * u +
                      fBins->GetParameter(3) * u * u * u);
    }, fitLow, fitHigh, 0);
    auto drawn = [&](TF1* f, const char* name) {   // a curve of the fit (entries per 2 keV in measureDays days) per day
      return new TF1(name, [=](double* x, double*) { return f->Eval(x[0]) / measureDays; }, fitLow, fitHigh, 0);
    };
    // the fit from the bins alone (fBins as last fitted: to this spectrum)
    TF1 *fitCurve = drawn(fBins, "fitDrawn"), *backgroundCurve = drawn(background, "backgroundDrawn");
    const double curves = fitCurve->GetMaximum(fitLow, fitHigh);
    auto upperFrame = gPad->DrawFrame(fitLow, 0., fitHigh, 2. * std::max(top, curves), ";;entries per day per 2 keV");
    upperFrame->GetXaxis()->SetLabelSize(0.);   // shared with the panel below
    upperFrame->GetYaxis()->SetTitleSize(0.052); upperFrame->GetYaxis()->SetLabelSize(0.045); upperFrame->GetYaxis()->SetTitleOffset(1.05);
    TH1D* line[3];
    for (int k = 2; k >= 0; --k) {
      line[k] = static_cast<TH1D*>(after[k]->Clone(Form("fitLine%d", k)));
      line[k]->Scale(1. / measureDays);
      line[k]->SetLineColor(kindColour[k]); line[k]->SetLineWidth(2);
      line[k]->Draw("HIST SAME");
    }
    total->SetLineColor(kBlack); total->SetLineWidth(3); total->SetMarkerSize(0);
    total->Draw("HIST SAME");
    total->Draw("E1 SAME");
    drawDashed(backgroundCurve, 14., 10., kRed, 4);
    fitCurve->SetLineColor(kRed); fitCurve->SetLineWidth(4); fitCurve->SetNpx(2000);
    fitCurve->Draw("SAME");
    {   // the simulated lines, dotted, each named just above the fit there (not at one height: the legend takes the right)
      TLatex mark;
      mark.SetTextSize(0.03); mark.SetTextColor(kGray + 2); mark.SetTextAngle(90); mark.SetTextAlign(12);
      const std::string el = gTarget;
      for (const auto& [eKeV, name] : {std::pair<double, std::string>{kAlpha, el + " L_{#alpha}"}, {kBeta, el + " L_{#beta}"}}) {
        double yTop = fitCurve->Eval(eKeV);   // above the fit and the data's error bars within 4 keV
        for (int b = total->FindBin(eKeV - 4.); b <= total->FindBin(eKeV + 4.); ++b)
          yTop = std::max(yTop, total->GetBinContent(b) + total->GetBinError(b));
        yTop += 0.05 * std::max(top, curves);
        auto l = new TLine(eKeV, 0., eKeV, yTop);
        l->SetLineColor(kGray + 2); l->SetLineStyle(3);
        l->Draw();
        mark.DrawLatex(eKeV + 2., yTop, name.c_str());
      }
    }
    // the fitted peaks' positions mu (asked for on 2026-11-03), in both panels: a dashed red line up to the fit and its
    // value, not in the legend. Upper panel: the values high above the fit with arrows down to the peaks; lower panel:
    // beside each peak, on the side away from the other
    const double mu[2] = {fitShown.position, fitShown.position2}, muWidth[2] = {fitShown.width, fitShown.width2};
    const double muError[2] = {fitShown.positionError, fitShown.positionError2};
    const int side[2] = {mu[0] > mu[1] ? 1 : -1, mu[1] > mu[0] ? 1 : -1};
    auto muGap = [&](int k) { return k ? 1. : 2. * muWidth[k] + 0.5; };   // from the dashed line to its value
    auto muLine = [&](double x, double yLow, double yHigh) {
      auto l = new TLine(x, yLow, x, yHigh);
      l->SetLineColor(kRed); l->SetLineStyle(2); l->SetLineWidth(3);
      l->Draw();
    };
    {
      // the values high above the fit, an arrow from each down to its peak (asked for on 2026-11-03): peak 1's centred
      // over its peak, peak 2's a row lower, ending 3 keV right of its peak; inside the frame, left of the fit's box
      const double frameTop = 2. * std::max(top, curves);
      const double boxLeft = fitLow + (0.545 - 0.12) / 0.85 * (fitHigh - fitLow);   // the box's left edge (keV)
      for (int k = 0; k < 2; ++k) {
        const double height = fitCurve->Eval(mu[k]);
        muLine(mu[k], 0., height);
        const std::string text = Form("#mu_{%d} = %.2f #pm %.2f keV", k + 1, mu[k], muError[k]);
        TLatex label(0., 0., text.c_str());
        label.SetTextColor(kRed); label.SetTextSize(0.03); label.SetTextAlign(31);
        const double width = label.GetXsize(), y = (k ? 0.62 : 0.72) * frameTop;
        const double xEnd = std::min(std::max(k ? mu[k] + 3. : mu[k] + 0.5 * width, fitLow + 1. + width), boxLeft - 1.);
        label.DrawLatex(xEnd, y, text.c_str());
        auto arrow = new TArrow(mu[k], y - 0.015 * frameTop, mu[k], height + 0.02 * frameTop, 0.012, "|>");
        arrow->SetLineColor(kRed); arrow->SetFillColor(kRed); arrow->SetLineWidth(2);
        arrow->Draw();
      }
    }
    // the fit's numbers, in a box as ROOT shows a fit
    auto box = new TPaveText(0.545, 0.59, 0.965, 0.87, "NDC");
    box->SetBorderSize(1); box->SetFillColor(kWhite); box->SetTextAlign(12); box->SetTextFont(42); box->SetTextSize(0.026);
    box->AddText(Form("fit to the simulated bins alone, %.2f days", measureDays));
    if (kWidthMax > kWidthMin)
      box->AddText(Form("widths %g-%g #times resolution (%.2f keV at #mu_{1}), %s", kWidthMin, kWidthMax, SigmaKeV(fitShown.position),
                        backName));
    else
      box->AddText(Form("widths fixed: #sigma = %.4g keV #sqrt{E/MeV}; %s", gSigma1MeV, backName));
    auto widthText = [&](double w, double e) {   // a string of its own: Form's buffer is reused by the Form it goes into
      return std::string(kWidthMax > kWidthMin ? Form("width %.2f #pm %.2f keV", w, e) : Form("width %.2f keV", w));
    };
    box->AddText(Form("peak 1: #mu = %.2f #pm %.2f keV, %s", fitShown.position, fitShown.positionError,
                      widthText(fitShown.width, fitShown.widthError).c_str()));
    box->AddText(Form("   %.1f #pm %.1f entries, significance %.1f #sigma", fitShown.amount, fitShown.error, fitShown.z));
    box->AddText(Form("peak 2: #mu = %.2f #pm %.2f keV, %s", fitShown.position2, fitShown.positionError2,
                      widthText(fitShown.width2, fitShown.widthError2).c_str()));
    box->AddText(Form("   %.1f #pm %.1f entries, significance %.1f #sigma", fitShown.amount2, fitShown.error2, fitShown.z2));
    box->AddText(Form("#chi^{2} / ndf   %.1f / %d", chi2Bins, ndfBins));
    box->Draw();
    auto dashed = [](int col, int width) { auto l = new TLine(); l->SetLineStyle(2); l->SetLineWidth(width); l->SetLineColor(col); return l; };
    auto legend = legendAt(0.57, 0.25, 0.96, 0.57, 0.03);
    legend->AddEntry(total, Form("simulated total: %s entries", count(n.total[0] + n.total[1] + n.total[2]).c_str()), "le");
    for (int k = 0; k < 3; ++k) legend->AddEntry(line[k], Form("%s: %s entries", shortNames[k], count(n.total[k]).c_str()), "l");
    legend->AddEntry(fitCurve, kWidthMax > kWidthMin ? "fit: two peaks, widths #geq resolution" : "fit: two peaks of the resolution's width",
                     "l");
    legend->AddEntry(dashed(kRed, 4), "its smooth background", "l");
    auto dotted = new TLine(); dotted->SetLineStyle(3); dotted->SetLineColor(kGray + 2);
    legend->AddEntry(dotted, "dotted: simulated lines, not fitted", "l");
    legend->Draw();
    TLatex text;
    text.SetNDC(); text.SetTextSize(0.032);
    text.DrawLatex(0.15, 0.84, Form("after the corner-plot cuts, %.2f days:", measureDays));
    text.DrawLatex(0.15, 0.795, Form("   %s target events, %s others", count(n.sEvents).c_str(), count(n.bEvents).c_str()));
    text.DrawLatex(0.15, 0.75, Form("   (their entries at %.0f-%.0f keV: the legend)", fitLow, fitHigh));
    gPad->RedrawAxis();
    // the lower panel: the simulated spectrum minus the fitted smooth background (the simulation's errors), and the fitted lines
    lower->cd();
    auto residual = new TH1D("fitResidual", "", bHigh - bLow + 1, after[0]->GetBinLowEdge(bLow), after[0]->GetBinLowEdge(bHigh + 1));
    residual->Sumw2();
    double rLow = 0., rHigh = 0.;
    for (int b = bLow; b <= bHigh; ++b) {
      const double r = total->GetBinContent(b) - backgroundCurve->Eval(total->GetBinCenter(b)), error = total->GetBinError(b);
      residual->SetBinContent(b - bLow + 1, r);
      residual->SetBinError(b - bLow + 1, error);
      rLow = std::min(rLow, r - error); rHigh = std::max(rHigh, r + error);
    }
    const double yLow = 1.25 * rLow, yHigh = 1.6 * rHigh;   // room above for the caption
    auto lowerFrame = gPad->DrawFrame(fitLow, yLow, fitHigh, yHigh, Form(";%s;residual", energyTitle));
    lowerFrame->GetXaxis()->SetLabelSize(0.087); lowerFrame->GetXaxis()->SetTitleSize(0.1); lowerFrame->GetXaxis()->SetTitleOffset(1.05);
    lowerFrame->GetYaxis()->SetLabelSize(0.087); lowerFrame->GetYaxis()->SetTitleSize(0.1); lowerFrame->GetYaxis()->SetTitleOffset(0.55);
    lowerFrame->GetYaxis()->SetNdivisions(505);
    auto zero = new TLine(fitLow, 0., fitHigh, 0.);
    zero->SetLineColor(kGray + 2); zero->SetLineStyle(2); zero->SetLineWidth(2);
    zero->Draw();
    for (double e : {kAlpha, kBeta}) {
      auto mark = new TLine(e, yLow, e, yLow + 0.75 * (yHigh - yLow));   // below the caption
      mark->SetLineColor(kGray + 2); mark->SetLineStyle(3);
      mark->Draw();
    }
    residual->SetLineColor(kBlack); residual->SetLineWidth(3); residual->SetMarkerSize(0);
    residual->Draw("HIST SAME");
    residual->Draw("E1 SAME");
    auto lines = new TF1("lines", [=](double* x, double*) { return fitCurve->Eval(x[0]) - backgroundCurve->Eval(x[0]); }, fitLow, fitHigh, 0);
    lines->SetLineColor(kRed); lines->SetLineWidth(4); lines->SetNpx(2000);
    lines->Draw("SAME");
    {
      TLatex label;
      label.SetTextColor(kRed); label.SetTextSize(0.06);
      for (int k = 0; k < 2; ++k) {
        const double height = lines->Eval(mu[k]);
        muLine(mu[k], yLow, height);
        label.SetTextAlign(side[k] > 0 ? 12 : 32);
        label.DrawLatex(mu[k] + side[k] * muGap(k), k ? height + 0.15 * (yHigh - yLow) : 0.55 * height,
                        Form("#mu_{%d} = %.2f keV", k + 1, mu[k]));
      }
    }
    TLatex panel;
    panel.SetNDC(); panel.SetTextSize(0.075);
    panel.DrawLatex(0.15, 0.85, "the simulated spectrum minus the fitted background; red: the fitted peaks");
    gPad->RedrawAxis();
    c->cd();
    TLatex title;
    title.SetNDC(); title.SetTextSize(0.036);
    DrawFitted(title, 0.12, 0.945, 0.98, Form("The %s L-line fit, %s", gTarget, label));
    save(c, "measurement_fit", 0.0155);
  }
  printf("\nplots: %s/muonic_xray_select_{corner_after,cutflow,spectrum_before,spectrum_after,"
         "measurement_fit}.png (and .pdf); table: %s/muonic_xray_select.csv\n", out.Data(), out.Data());
}
