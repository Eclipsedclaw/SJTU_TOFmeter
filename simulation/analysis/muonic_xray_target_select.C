// One thicker Al target (default 4 cm), its simulated data split into the target's muonic X-rays (signal, red) and
// everything else (background, blue): which measured variables tell them apart, the most efficient cuts on them, and
// whether a measurement could still find the Al lines after the cuts. From the skims (analysis/muonic_xray_stack_skim.C)
// of the target's run of scripts/muonic_xray_target.sh and of the 3-week run:
//   root -l -b -q 'analysis/muonic_xray_target_select.C'
// The data: every stack entry of 20-1000 keV with a middle-TOF hit in coincidence (-5..+30 ns), from the target's run
// (cosmic muons up to eSplit) and from the 3-week run's faster muons, each per day of its own exposure (the 3-week run's
// slower muons belong to its 1 cm target). Signal: an entry with a muonic X-ray of the target (stk_src bit 1); the
// background holds the muonic X-rays of the other volumes (Al plates, carbon of the TOF bars), which make the same or
// nearby lines. The TOF response of muonic_xray_stack_cuts.C (sigmaT per station, sigmaE on the bar energies, its
// seeds; the top-TOF energy smeared alike).
// 1. Corner plots of the candidate variables, the muon's (bottom-TOF hit, middle-TOF bar, middle-TOF energy,
//    top-to-middle time, top-TOF energy) and the stack entry's (highest and lowest layer hit, hits in the entry, largest
//    hit / entry, outermost hit), signal as red points, background as a blue density. Since 2026-11-02 on the captured
//    data at every energy, as a measurement records it (the TOF trigger, an entry in time; one point per event,
//    signal an event with a muonic X-ray of the target), with the cuts chosen in the L_alpha window drawn. Every
//    event inside the axes (asked for on 2026-11-02): with the measured TOF the bar energies and the kinetic energy on
//    log axes, the kinetic energy's last bin the events with beta >= 1 (no kinetic energy); with single pixels the
//    other pixels' energy on a log axis, its first bin the entries without any.
// 2. The most efficient cuts: greedy, one cut at a time the variable and threshold that raise S/sqrt(B) in the L_alpha
//    window the most (by at least 1%), then each threshold again with the others fixed; chosen on the even events and
//    scored on the odd ones and back, to see the optimisation's own gain. The CZT run of single pixels with the measured
//    TOF uses cuts chosen by hand instead (2026-11-02: keeping more of the signal), the most efficient ones printed and
//    fitted beside them. The cut flow.
// 3. The spectrum before and after the cuts by what made it (20-450 keV, per day and 2 keV, drawn as ROOT histograms:
//    the simulated total with its statistical error bars, each class a line), and fits without truth to the expected
//    spectrum after the cuts itself, what a measurement of measureDays days expects (the MC truth, the Asimov data set;
//    since 2026-11-02 in place of simulated measurements drawn from it): smooth background and Gaussian Al L_alpha and
//    L_beta, either with free position and width (background exp(pol2), as muonic_xray_target_truth.C fits the target's
//    X-rays alone) or at the known line energies with the stack's known width (background exp(pol3), only the amounts
//    free), with the likelihood-ratio significance of L_alpha (measurement_fit). A null test: the same fits to the
//    expected spectrum without any muonic X-ray. The cut sets compared, each with the known-energy fit (with the stack's
//    skims also those of the 1 cm target and the earlier own cuts of muonic_xray_target_cuts.C).
// Writes <out>/muonic_xray_select_<plot>.png/.pdf (one plot per file) and muonic_xray_select.csv.
// The CZT plane under the 4 cm target (scripts/muonic_xray_czt.sh, its bands skimmed with detector "czt"): detector
// "czt" and the CZT's resolution, the slow band's skim as the target's run, the fast band's as the faster muons; the
// entry's variables are then the pixel hits, the largest hit / entry and whether it lies in the 4 central pixels:
//   root -l -b -q 'analysis/muonic_xray_target_select.C("output/muonic_xray_czt/slow/xray_skim.root", "output/muonic_xray_czt/fast/xray_skim.root", "output/muonic_xray_czt", "4 cm target, CZT", "czt", 6.37)'
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
#include <TParameter.h>

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
double SigmaKeV(double eKeV) { return gSigma1MeV * std::sqrt(eKeV / 1000.); }
const double kAlpha = 65.73, kBeta = 88.73;                              // Al L_alpha and L_beta (keV)
bool InAlpha(double eKeV) { return std::abs(eKeV - kAlpha) <= 2. * SigmaKeV(kAlpha); }

// the candidate cut variables, all measured (no truth); the CZT has no layers (variables 5 and 6), and a single pixel
// no hits to add up (7 and 8); 10 and 11, the bottom-TOF bar and the track's zenith angle, need the measured TOF, 12,
// the pixel, single CZT pixels; 13, the entry's energy with the response, is shown in the corner plots (since
// 2026-11-02), not cut on (the L_alpha window defines the signal). 0-4, 10 and 11 are the muon's, the others the entry's
const int kNV = 14;
// a log axis (log) from lo to hi; split > 0: a bin of its own at the axis' end for the values at or above split (or,
// splitLow, below it), such as a placeholder, drawn at the bin's centre and named splitLabel in the plots
struct Variable {
  const char* name; int bins; double lo, hi; bool integer;
  bool log = false; double split = 0.; bool splitLow = false; const char* splitLabel = "";
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
                             {"entry energy (keV)", 104, 0., 1040., false}};   // room for the axis label "1000"

// a variable's bin edges: even, or even in log on a log axis, the split bin at its end
std::vector<double> Edges(const Variable& v)
{
  std::vector<double> e;
  if (!v.log) {
    for (int k = 0; k <= v.bins; ++k) e.push_back(v.lo + (v.hi - v.lo) * k / v.bins);
    return e;
  }
  const bool split = v.split > 0.;
  const double a = split && v.splitLow ? v.split : v.lo, b = split && !v.splitLow ? v.split : v.hi;   // the part even in log
  const int n = split ? v.bins - 1 : v.bins;
  if (split && v.splitLow) e.push_back(v.lo);
  for (int k = 0; k <= n; ++k) e.push_back(a * std::pow(b / a, double(k) / n));
  if (split && !v.splitLow) e.push_back(v.hi);
  return e;
}
// where a value is drawn: beyond the split at the split bin's centre, else within the axis' range
double ShownAt(double x, const Variable& v)
{
  const bool split = v.split > 0.;
  if (split && (v.splitLow ? x < v.split : x >= v.split)) return v.splitLow ? std::sqrt(v.lo * v.split) : std::sqrt(v.split * v.hi);
  const double lo = split && v.splitLow ? v.split : v.lo, hi = split && !v.splitLow ? v.split : v.hi;
  return v.log ? std::min(std::max(x, lo * (1. + 1e-9)), hi * (1. - 1e-9)) : std::min(std::max(x, lo + 1e-6), hi - 1e-6);
}

struct Entry {
  double w;        // per day of its sample's exposure
  int kind, half;  // 0 muonic X-ray of the target (signal), 1 muonic X-ray of another volume, 2 the rest; event parity
  int sample;      // 0 the target's run (muons up to eSplit), 1 the 3-week run's faster muons
  bool trig;       // the TOF trigger fired (always true for skims without the measured TOF)
  bool timed;      // the top TOF timed the muon: variables 3 and 4 exist (and with the measured TOF the track, 1)
  double e;        // add-back energy with the response (keV)
  double v[kNV];
  long long id;    // the event: (sample x 1000 + process) x 10^8 + event number
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
      if (gCZT) range(-1e9, 2.5);
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
    case 2: return Form("middle-TOF energy > %.2f MeV", c.lo);
    case 3:
      if (gTOF) return c.lo > 0. ? Form("TOF kinetic energy %.0f-%.0f MeV", c.lo, c.hi) : Form("TOF kinetic energy < %.0f MeV", c.hi);
      return Form("top-to-middle time > %.1f ns", c.lo);
    case 4: return Form("top-TOF energy > %.2f MeV", c.lo);
    case 5:
      if (gCZT) return Form("at most %d pixels hit", int(std::floor(c.hi)));
      return Form("highest layer >= %d", int(std::ceil(c.lo)));
    case 6:
      if (gCZT) return Form("other pixels < %.2f MeV", c.hi);
      return Form("lowest layer >= %d", int(std::ceil(c.lo)));
    case 7: return Form("at most %d hits", int(std::floor(c.hi)));
    case 8: return Form("largest hit > %.2f of the entry", c.lo);
    case 9: return gCZT ? std::string("one of the 4 central pixels") : std::string(Form("outermost hit < %.1f cm", c.hi));
    case 11: return Form("track within %.0f deg of the vertical", c.hi);
    default: return "none";
  }
}
}  // namespace

void muonic_xray_target_select(const char* slowSkim = "output/muonic_xray_target/t4cm/xray_skim.root",
                               const char* ref = "output/muonic_xray_stack/xray_skim_pbke.root",
                               const char* outDir = "output/muonic_xray_target/t4cm", const char* label = "4 cm target",
                               const char* detector = "stack", double sigma1MeV = 40., double eSplit = 1000.,
                               double measureDays = 21., int toys = 300, double sigmaT = 0.2, double sigmaE = 0.15,
                               unsigned seed = 4357)
// toys and seed: not used since 2026-11-02 (the fits take the expected spectrum, no simulated measurements); kept so that
// the arguments stay in their places
{
  TH1::AddDirectory(false);
  const TString out = outDir;
  gSigma1MeV = sigma1MeV;
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
    // on a log axis to 10^3.5 MeV, 8 bins a decade; its first bin the entries without energy in the other pixels (a pixel
    // counts from 10 keV; below that only the rounding of a sum, |value| < 1e-6 MeV)
    kVars[6] = {"other pixels' energy (MeV)", 45, std::pow(10., -2.5), std::pow(10., 3.5), false, true, 0.01, true, "0"};
    kVars[9] = {"pixel centre (cm)", 5, 0., 5., false};   // max(|x|,|y|): 1.25 inner, 3.75 outer
  } else if (gCZT) {      // the CZT entry: pixel hits, the largest pixel hit, the pixel centres (1.25 or 3.75 cm out)
    kVars[7] = {"pixel hits in the entry", 5, 0.5, 5.5, true};
    kVars[8] = {"largest hit / entry", 26, 0., 1.04, false};
    kVars[9] = {"outermost pixel centre (cm)", 5, 0., 5., false};
  }
  stamp[0] = gCZT ? Form("CZT plane, %s; FWHM/E %.1f%% at 1 MeV, 10 keV threshold", gPixel ? "each of its 16 pixels read out alone"
                                                                                           : "its pixels summed within 100 ns",
                         2.3548 * sigma1MeV / 10.)
                  : Form("Si/YSO stack, add-back within 100 ns; FWHM/E %.1f%% at 1 MeV", 2.3548 * sigma1MeV / 10.);
  stamp[1] = gTOF ? Form("TOF: %.0f ps FWHM per station, bar energy %.0f%% FWHM, both bar ends read out (%.1f cm along a bar); "
                         "trigger: top and middle", 1000. * tofFwhmT, 100. * tofFwhmE, vEff * tofFwhmT / 2.3548)
                  : Form("TOF response of this study: sigma %.2f ns per station, sigma %.0f%% on the bar energies; middle-TOF "
                         "coincidence, no trigger", sigmaT, 100. * sigmaE);
  std::vector<Entry> entries;
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
        if (tMidM < 0. || t - tMidM < -5. || t - tMidM > 30.) continue;
        if (fastOnly && (kind == 0 || ekin <= eSplit)) continue;
        const bool bottom = botM >= 0. && botM - tMidM > -5. && botM - tMidM < 50.;
        entries.push_back({1. / skimDays, kind, event & 1, fastOnly ? 1 : 0, trig != 0, trig != 0, e,
                           {bottom ? 1. : 0., std::max(std::abs(xTgt), std::abs(yTgt)), eMidM, tkin, trig ? eTopM : -1.,
                            gPixel ? double(npix) : double(topLayer), gPixel ? eOther / 1000. : double(lowest), double(n),
                            e > 0. ? emax / e : 0., rmax, bottom && barBotM >= 0 ? barBotM + 1. : 0.,
                            trig ? std::atan(std::hypot(slopeX, slopeY)) * 180. / M_PI : -1., pix + 1., e}, id});
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
  // the TOF trigger takes the dataset out of the middle-TOF coincidences: what it keeps in the L_alpha window (per
  // day), then only what it keeps
  double coinS = 0., coinB = 0., trigS = 0., trigB = 0.;
  for (const auto& x : entries) {
    if (!InAlpha(x.e)) continue;
    (x.kind == 0 ? coinS : coinB) += x.w;
    if (x.trig) (x.kind == 0 ? trigS : trigB) += x.w;
  }
  if (gTOF) entries.erase(std::remove_if(entries.begin(), entries.end(), [](const Entry& x) { return !x.trig; }), entries.end());
  std::vector<Entry> alpha;   // the Al L_alpha window
  for (const auto& x : entries) if (InAlpha(x.e)) alpha.push_back(x);
  const double fitLow = 30., fitHigh = 160.;   // the range of the fits to simulated measurements (keV)
  const double plotLow = 20., plotHigh = 450.;   // the spectra's drawn range (keV; little lies above 450 keV, 2026-11-02)

  // the conditions of this study, after its skims'; printed first, and kept with the outputs
  conditions.push_back(Form("study: analysis/muonic_xray_target_select.C, %s, detector %s, written %s; outputs in %s", label,
                            detector, TDatime().AsSQLString(), outDir));
  conditions.push_back(Form("  samples: %s, all entries (%.2f days, muons up to %g MeV); %s, the background of primaries above "
                            "%g MeV only (%.2f days); each per day of its own exposure", slowSkim, thickDays, eSplit, ref, eSplit,
                            refDays));
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
  conditions.push_back(Form("  cuts: greedy on S/sqrt(B) in %g days in the L_alpha window, each step at least 1%% better and keeping "
                            "at least 30 signal entries (60 on all events), chosen on even and odd events alone to cross-validate",
                            measureDays));
  conditions.push_back(Form("  fits: to the expected spectrum of %g days (the MC truth after the cuts, the Asimov data set), %g-%g keV",
                            measureDays, fitLow, fitHigh));
  if (gTOF && gCZT && gPixel)
    conditions.push_back("  analysis cuts: chosen by hand on 2026-11-02 (track at the target within 10 cm, within 40 deg of the "
                         "vertical, middle-TOF energy > 1.6 MeV, TOF kinetic energy 10-80 MeV, top-TOF energy > 1.6 MeV, other "
                         "pixels < 0.5 MeV); the most efficient cuts printed beside them");
  printf("\n=== conditions ===\n");
  for (const auto& line : conditions) printf("%s\n", line.c_str());
  {
    std::ofstream file((out + "/muonic_xray_select_conditions.txt").Data());
    for (const auto& line : conditions) file << line << "\n";
  }

  const char* kindNames[3] = {"muonic X-rays of the target", "muonic X-rays of other volumes", "the rest"};
  printf("\n=== %s: simulated %s entries with a middle-TOF hit in coincidence%s ===\n", label, gCZT ? "CZT" : "stack",
         gTOF ? " and the TOF trigger" : "");
  if (gTOF)
    printf("the trigger (top and middle) keeps in the L_alpha window: target X-rays %.2f of %.2f per day (%.0f%%), everything else "
           "%.1f of %.1f (%.0f%%)\n", trigS, coinS, 100. * trigS / std::max(1e-12, coinS), trigB, coinB,
           100. * trigB / std::max(1e-12, coinB));
  printf("%s (%.2f days, muons up to %.0f MeV) and the faster muons of %s (%.2f days); energy sigma %.2f keV at 1 MeV "
         "(L_alpha window +-%.1f keV)\n", slowSkim, thickDays, eSplit, ref, refDays, gSigma1MeV, 2. * SigmaKeV(kAlpha));
  for (int k = 0; k < 3; ++k) {
    double all = 0., inAlpha = 0.;
    long nAll = 0, nAlpha = 0;
    for (const auto& x : entries) if (x.kind == k) { all += x.w; ++nAll; if (InAlpha(x.e)) { inAlpha += x.w; ++nAlpha; } }
    printf("   %-32s %8.1f per day (%7ld entries); in the L_alpha window %7.1f per day (%6ld entries)\n", kindNames[k], all,
           nAll, inAlpha, nAlpha);
  }

  // ---- S/sqrt(B) in measureDays days in the L_alpha window, from half h (0, 1; each with twice its weight) or both (2)
  auto score = [&](const Cuts& cuts, int h, double* sDay = nullptr, double* bDay = nullptr, long* nSignal = nullptr) {
    double s = 0., b = 0.;
    long n = 0;
    for (const auto& x : alpha) {
      if (h < 2 && x.half != h) continue;
      if (!PassesAll(x, cuts)) continue;
      const double w = h < 2 ? 2. * x.w : x.w;
      if (x.kind == 0) { s += w; ++n; }
      else b += w;
    }
    if (sDay) *sDay = s;
    if (bDay) *bDay = b;
    if (nSignal) *nSignal = n;
    return b > 0. ? std::sqrt(measureDays) * s / std::sqrt(b) : 0.;
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
  // the cuts of the analysis: for the CZT run of single pixels with the measured TOF, chosen by hand on 2026-11-02 to keep
  // more of the signal than the most efficient cuts (which give up signal to reject more background), in this order;
  // else the most efficient cuts
  const bool handChosen = gTOF && gCZT && gPixel;
  Cuts cuts = best;
  std::vector<int> order = bestOrder;
  if (handChosen) {
    cuts = Cuts{};
    order = {1, 11, 2, 3, 4, 6};
    auto keep = [&](int i, double lo, double hi) { cuts[i].active = true; cuts[i].lo = lo; cuts[i].hi = hi; };
    keep(1, -1e9, 10.);   // the track at the target within 10 cm
    keep(11, -1e9, 40.);  // within 40 deg of the vertical
    keep(2, 1.6, 1e9);    // middle-TOF energy above 1.6 MeV
    keep(3, 10., 80.);    // TOF kinetic energy 10-80 MeV
    keep(4, 1.6, 1e9);    // top-TOF energy above 1.6 MeV
    keep(6, -1e9, 0.5);   // the other pixels below 0.5 MeV (the most efficient cuts' CZT cut)
  }
  const char* cutsName = handChosen ? "the cuts of this analysis (chosen by hand)" : "the most efficient cuts (chosen in the L-alpha window)";
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

  printf("\nthe most efficient cuts (greedy, S/sqrt(B) in %.0f days in the L_alpha window; even and odd events alone):\n",
         measureDays);
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
    printf("the cuts of this analysis (chosen by hand on 2026-11-02, to keep more of the signal):\n");
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
      flowZ.push_back(coinB > 0. ? std::sqrt(measureDays) * coinS / std::sqrt(coinB) : 0.);
      flowNames.push_back("TOF trigger: top and middle");
    }
    flowZ.push_back(score(partial, 2, &s, &b)); flowS.push_back(s); flowB.push_back(b);
    for (int i : order) {
      partial[i] = cuts[i];
      flowZ.push_back(score(partial, 2, &s, &b)); flowS.push_back(s); flowB.push_back(b);
      flowNames.push_back(Describe(cuts, i));
    }
  }
  printf("\ncut flow in the L_alpha window, per day:\n   %-40s %9s %9s %16s\n", "", "S", "B", Form("S/sqrt(B) %.0f d", measureDays));
  std::ofstream csv((out + "/muonic_xray_select.csv").Data());
  csv << "table,step,cut,S_day,B_day,S_sqrtB,median_Z_fit,Z16,Z84\n";
  for (std::size_t k = 0; k < conditions.size(); ++k) {   // the conditions first, one row per line
    std::string line = conditions[k];
    std::replace(line.begin(), line.end(), '"', '\'');
    csv << "condition," << k << ",\"" << line << "\",,,,,,\n";
  }
  for (std::size_t k = 0; k < flowS.size(); ++k) {
    printf("   %zu %-38s %9.2f %9.1f %16.2f\n", k, flowNames[k].c_str(), flowS[k], flowB[k], flowZ[k]);
    csv << "cutflow," << k << ",\"" << flowNames[k] << "\"," << flowS[k] << "," << flowB[k] << "," << flowZ[k] << ",,,\n";
  }

  // ---- the spectra by what made them, per day, 0-200 keV in 2 keV
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

  // what the L_alpha window holds after the cuts, per day: by kind, and the rest by sample
  {
    double kinds[3] = {}, rest[2] = {};
    for (const auto& x : alpha) {
      if (!PassesAll(x, cuts)) continue;
      kinds[x.kind] += x.w;
      if (x.kind == 2) rest[x.sample] += x.w;
    }
    printf("\nthe L_alpha window after the cuts, per day: %s %.1f, %s %.1f, %s %.1f (from muons up to %.0f MeV %.1f, from faster "
           "ones %.1f)\n", kindNames[0], kinds[0], kindNames[1], kinds[1], kindNames[2], kinds[2], eSplit, rest[0], rest[1]);
  }

  // ---- fits without truth to the expected spectrum itself: what the simulation expects in measureDays days, bin by bin
  // (the MC truth after the cuts, the Asimov data set: a fit to it gives what a measurement expects, its median result;
  // asked for on 2026-11-02 in place of simulated measurements drawn from it), in two ways:
  //  free:  smooth background exp(pol2), Gaussian L_alpha and L_beta with free position and width (as
  //         muonic_xray_target_truth.C fits the target's X-rays alone);
  //  known: background exp(pol3), the Gaussian lines at their known energies with the stack's known width, only their
  //         amounts free (as a calibrated measurement fits them; an amount may come out negative).
  // Left out on 2026-11-02: a background shape known in advance (the simulated rest after the cuts: a cosmic background
  // cannot be predicted that exactly, and its template missed the target's own continuum), and the background shape of
  // the measurement's own control region (middle-TOF energy below its cut), which holds more of the target's X-rays
  // than background: it found L_alpha in the null test (median 4.2 sigma) and 23 sigma in a simulated measurement.
  // The L_alpha significance: the likelihood ratio against the same fit without L_alpha, with the sign of its amount.
  // A null test: the expected spectrum without any muonic X-ray, where an honest fit finds no L_alpha. Fit range:
  // fitLow-fitHigh. The fits take the expectation as exact: the simulation's own statistics (the plots' error bars) are
  // not in them.
  auto gaussian = [](double e, double pos, double width) {
    return std::exp(-0.5 * std::pow((e - pos) / width, 2)) / (std::sqrt(2. * M_PI) * width);
  };
  auto freeModel = [&](const char* name) {
    return new TF1(name, [=](double* x, double* p) {
      const double e = x[0];
      return std::exp(p[0] + p[1] * (e - 80.) + p[2] * (e - 80.) * (e - 80.)) +
             2. * (p[3] * gaussian(e, p[4], p[5]) + p[6] * gaussian(e, p[4] + kBeta - kAlpha, p[5] * std::sqrt(kBeta / kAlpha)));
    }, fitLow, fitHigh, 7);
  };
  auto knownModel = [&](const char* name) {
    return new TF1(name, [=](double* x, double* p) {
      const double u = (x[0] - 80.) / 50.;
      return std::exp(p[0] + p[1] * u + p[2] * u * u + p[3] * u * u * u) +
             2. * (p[4] * gaussian(x[0], kAlpha, SigmaKeV(kAlpha)) + p[5] * gaussian(x[0], kBeta, SigmaKeV(kBeta)));
    }, fitLow, fitHigh, 6);
  };
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");   // classic Minuit, as muonic_xray_stack_cuts.C
  TF1 *fFree = freeModel("free"), *fFreeNull = freeModel("freeNull"), *fKnown = knownModel("known"), *fKnownNull = knownModel("knownNull");
  fFreeNull->FixParameter(3, 0.); fFreeNull->FixParameter(4, kAlpha); fFreeNull->FixParameter(5, SigmaKeV(kAlpha));
  fKnownNull->FixParameter(4, 0.);
  auto nll = [&](TH1D* h, TF1* f) {   // Poisson -log L of the fit range, at the bin centres as the fit
    double sum = 0.;
    for (int b = h->FindBin(fitLow + 1e-6); b <= h->FindBin(fitHigh - 1e-6); ++b) {
      const double mu = std::max(1e-12, f->Eval(h->GetBinCenter(b)));
      sum += mu - h->GetBinContent(b) * std::log(mu);
    }
    return sum;
  };
  struct LineFit { double z, amount, error, position, width; };
  auto fitMeasurement = [&](TH1D* h, int how) {   // 0 free, 1 known
    double window = 0.;
    for (int b = h->FindBin(kAlpha - 2. * SigmaKeV(kAlpha)); b <= h->FindBin(kAlpha + 2. * SigmaKeV(kAlpha)); ++b) window += h->GetBinContent(b);
    const double back = std::log(std::max(1., h->GetBinContent(h->FindBin(80.))));
    if (how == 1) {
      fKnown->SetParameters(back, -0.3, 0., 0., 0.1 * window, 0.02 * window);
      fKnown->SetParLimits(4, -0.3 * window, 2. * window); fKnown->SetParLimits(5, -0.3 * window, 2. * window);
      h->Fit(fKnown, "LQNR");
      for (int p : {0, 1, 2, 3, 5}) fKnownNull->SetParameter(p, fKnown->GetParameter(p));
      fKnownNull->SetParLimits(5, -0.3 * window, 2. * window);
      h->Fit(fKnownNull, "LQNR");
      const double z = std::sqrt(std::max(0., 2. * (nll(h, fKnownNull) - nll(h, fKnown))));
      return LineFit{fKnown->GetParameter(4) < 0. ? -z : z, fKnown->GetParameter(4), fKnown->GetParError(4), kAlpha, SigmaKeV(kAlpha)};
    }
    fFree->SetParameters(back, -0.01, 0., 0.2 * window, kAlpha, SigmaKeV(kAlpha), 0.05 * window);
    fFree->SetParLimits(3, 0., 1e8); fFree->SetParLimits(4, 50., 80.); fFree->SetParLimits(5, 3., 30.); fFree->SetParLimits(6, 0., 1e8);
    h->Fit(fFree, "LQNR");
    for (int p : {0, 1, 2, 6}) fFreeNull->SetParameter(p, fFree->GetParameter(p));
    fFreeNull->SetParLimits(6, 0., 1e8);
    h->Fit(fFreeNull, "LQNR");
    return LineFit{std::sqrt(std::max(0., 2. * (nll(h, fFreeNull) - nll(h, fFree)))), fFree->GetParameter(3), fFree->GetParError(3),
                   fFree->GetParameter(4), fFree->GetParameter(5)};
  };
  // the expected spectrum of measureDays days (entries per 2 keV), with or without its muonic X-rays
  auto expected = [&](const std::array<TH1D*, 3>& h, bool withXrays, const char* name) {
    auto m = new TH1D(name, "", h[0]->GetNbinsX(), h[0]->GetXaxis()->GetXmin(), h[0]->GetXaxis()->GetXmax());
    for (int b = 1; b <= m->GetNbinsX(); ++b)
      m->SetBinContent(b, measureDays * (h[2]->GetBinContent(b) + (withXrays ? h[0]->GetBinContent(b) + h[1]->GetBinContent(b) : 0.)));
    return m;
  };
  double truthAlpha[3] = {}, truthFit[3] = {};   // muonic X-ray entries after the cuts in measureDays days: window, fit range
  for (const auto& x : entries) {
    if (!PassesAll(x, cuts)) continue;
    if (InAlpha(x.e)) truthAlpha[x.kind] += measureDays * x.w;
    if (x.e >= fitLow && x.e < fitHigh) truthFit[x.kind] += measureDays * x.w;
  }
  printf("\nfits without truth to the expected spectrum of %.0f days after the cuts (the MC truth, the Asimov data set):\n", measureDays);
  printf("   truth after the cuts: L_alpha window: target X-rays %.0f, other muonic X-rays %.0f, the rest %.0f; %.0f-%.0f keV: "
         "target X-rays %.0f, other muonic X-rays %.0f, the rest %.0f\n", truthAlpha[0], truthAlpha[1], truthAlpha[2], fitLow, fitHigh,
         truthFit[0], truthFit[1], truthFit[2]);
  const char* modelNames[2] = {"free position and width, background exp(pol2)", "known energies and width, background exp(pol3)"};
  TH1D *expect = expected(after, true, "expected"), *expectNull = expected(after, false, "expectedNull");
  // how well the known-energy fit (fKnown as last fitted) describes an expected spectrum, against the statistics of
  // measureDays days (Baker-Cousins: twice the log-likelihood ratio to a perfect description; a fluctuation of the
  // simulation larger than those of measureDays days counts too)
  const int ndfFit = expect->FindBin(fitHigh - 1e-6) - expect->FindBin(fitLow + 1e-6) + 1 - 6;
  auto chi2Of = [&](TH1D* h) {
    double chi2 = 0.;
    for (int b = h->FindBin(fitLow + 1e-6); b <= h->FindBin(fitHigh - 1e-6); ++b) {
      const double nB = h->GetBinContent(b), fB = std::max(1e-12, fKnown->Eval(h->GetBinCenter(b)));
      chi2 += 2. * (fB - nB + (nB > 0. ? nB * std::log(nB / fB) : 0.));
    }
    return chi2;
  };
  const int noteLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kWarning;   // not Minuit's notes on a start value outside its limits (the null fits' fixed L_alpha)
  // the cut sets compared, each with the known-energy fit to its own expected spectrum
  printf("\nthe cut sets compared: S/sqrt(B) in the L_alpha window, and the fit with the known energies and width (background "
         "exp(pol3)) to the expected spectrum of %.0f days:\n", measureDays);
  printf("   %-46s %8s %8s %10s %18s %14s %13s\n", "", "S/day", "B/day", "S/sqrt(B)", "L_alpha amount", "significance", "chi2/ndf");
  struct Set { const char* name; const Cuts* cuts; };
  std::vector<Set> sets;
  if (!gTOF) {   // the earlier sets cut on the middle-TOF bar and the time of flight, which the measured TOF replaces
    sets.push_back({"cuts of the 1 cm target", &cuts1cm});
    sets.push_back({"own cuts of muonic_xray_target_cuts.C", &cutsOwn});
  }
  sets.push_back({"the most efficient cuts (this optimisation)", &best});
  if (handChosen) sets.push_back({"the cuts of this analysis (chosen by hand)", &cuts});
  for (std::size_t k = 0; k < sets.size(); ++k) {
    double s, b;
    const double z = score(*sets[k].cuts, 2, &s, &b);
    TH1D* e = expected(spectrum(*sets[k].cuts, Form("set%zu", k)), true, Form("expectedSet%zu", k));
    const LineFit r = fitMeasurement(e, 1);
    const double c2 = chi2Of(e);
    printf("   %-46s %8.2f %8.1f %10.2f %10.1f +- %4.1f %8.1f sigma %7.1f/%d\n", sets[k].name, s, b, z, r.amount, r.error, r.z, c2,
           ndfFit);
    csv << "compare," << k << ",\"" << sets[k].name << "\"," << s << "," << b << "," << z << "," << r.z << "," << r.amount << ","
        << r.error << "\n";
    delete e;
  }
  // the null test, then the fits of this analysis (the known-energy one last: fKnown keeps it for the plots)
  LineFit fit[2], fitNull[2];
  for (int how : {0, 1}) fitNull[how] = fitMeasurement(expectNull, how);
  for (int how : {0, 1}) fit[how] = fitMeasurement(expect, how);
  gErrorIgnoreLevel = noteLevel;
  printf("\nthe fits to the expected spectrum of this analysis (%.0f-%.0f keV):\n   %-47s %-22s %24s %18s\n", fitLow, fitHigh, "",
         "", "L_alpha amount", "significance");
  for (int how : {0, 1})
    for (int n = 0; n < 2; ++n) {
      const LineFit& r = n ? fitNull[how] : fit[how];
      printf("   %-47s %-22s %14.1f +- %6.1f %12.1f sigma\n", modelNames[how], n ? "without muonic X-rays" : "expected", r.amount,
             r.error, r.z);
      csv << (n ? "null" : "expected") << "," << how << ",\"" << modelNames[how] << "\",,,," << r.z << "," << r.amount << ","
          << r.error << "\n";
    }
  const LineFit fitShown = fit[1];
  const double chi2Fit = chi2Of(expect);
  // what the simulation's own fluctuations add to it: their variance over that of measureDays days, bin by bin (a
  // simulated entry of the fast band counts 21 / 5.92 entries of a measurement, so its noise exceeds a measurement's)
  double chi2Mc = 0.;
  for (int b = expect->FindBin(fitLow + 1e-6); b <= expect->FindBin(fitHigh - 1e-6); ++b) {
    double var = 0.;
    for (int k = 0; k < 3; ++k) var += std::pow(measureDays * after[k]->GetBinError(b), 2);
    if (expect->GetBinContent(b) > 0.) chi2Mc += var / expect->GetBinContent(b);
  }
  printf("   the known-energy fit: L_beta %.1f +- %.1f; chi2/ndf %.1f/%d (Baker-Cousins, against %.0f days of statistics); the "
         "simulation's own fluctuations alone make about %.0f of it\n", fKnown->GetParameter(5), fKnown->GetParError(5), chi2Fit,
         ndfFit, measureDays, chi2Mc);

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
  auto heading = [](const char* text) {
    TLatex latex;
    latex.SetNDC(); latex.SetTextSize(0.045);
    latex.DrawLatex(0.12, 0.935, text);
  };
  auto frame = [](double x1, double y1, double x2, double y2, const char* titles) {
    auto f = gPad->DrawFrame(x1, y1, x2, y2, titles);
    f->GetXaxis()->SetTitleSize(0.045); f->GetYaxis()->SetTitleSize(0.045);
    f->GetXaxis()->SetTitleOffset(1.1); f->GetYaxis()->SetTitleOffset(1.25);
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
    std::vector<TH1D*> hS(nV), hB(nV);
    std::vector<std::vector<TH2D*>> hB2(nV, std::vector<TH2D*>(nV, nullptr));
    std::vector<std::vector<TGraph*>> gS(nV, std::vector<TGraph*>(nV, nullptr));
    for (int a = 0; a < nV; ++a) {
      const std::vector<double> ea = Edges(kVars[vars[a]]);
      hS[a] = new TH1D(Form("%s_s%d", name, a), "", ea.size() - 1, ea.data());
      hB[a] = new TH1D(Form("%s_b%d", name, a), "", ea.size() - 1, ea.data());
      for (int c = 0; c < a; ++c) {
        const std::vector<double> ec = Edges(kVars[vars[c]]);
        hB2[a][c] = new TH2D(Form("%s_b%d%d", name, a, c), "", ec.size() - 1, ec.data(), ea.size() - 1, ea.data());
        gS[a][c] = new TGraph();
      }
    }
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
        if (Passes(x, cuts, i)) (sig ? sKept : bKept)[a] += x.w;
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
    TCanvas c("corner", "corner plot", 1600, 1660);
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
    // a variable's split bin (e.g. beta >= 1): a grey line where it starts and its name in it, the variable along x or y
    auto splitMark = [](int v, bool alongX, double top = 0.97) {   // top: where a name along x hangs from
      const Variable& var = kVars[v];
      if (var.split <= 0.) return;
      gPad->Update();
      auto at = [](double f, bool x) {   // a fraction f of the frame along x or y, in the pad's coordinates
        const double u = x ? gPad->GetUxmin() + f * (gPad->GetUxmax() - gPad->GetUxmin())
                           : gPad->GetUymin() + f * (gPad->GetUymax() - gPad->GetUymin());
        return (x ? gPad->GetLogx() : gPad->GetLogy()) ? std::pow(10., u) : u;
      };
      const double centre = ShownAt(var.splitLow ? var.lo : var.hi, var);
      auto l = alongX ? new TLine(var.split, at(0., false), var.split, at(1., false)) : new TLine(at(0., true), var.split, at(1., true), var.split);
      l->SetLineColor(kGray + 2); l->SetLineStyle(2);
      l->Draw();
      TLatex text;
      text.SetTextFont(42); text.SetTextSize(0.06); text.SetTextColor(kGray + 2); text.SetTextAlign(32);
      text.SetTextAngle(alongX ? 90. : 0.);   // along x: hanging from the frame's top
      if (alongX) text.DrawLatex(centre, at(top, false), var.splitLabel);
      else text.DrawLatex(at(0.97, true), centre, var.splitLabel);
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
          TH1D *hb = hB[row], *hs = hS[row];
          for (auto h : {hb, hs}) if (h->Integral() > 0.) h->Scale(1. / h->Integral());
          const double top = std::max(hb->GetMaximum(), hs->GetMaximum());
          hb->SetMinimum(0.); hb->SetMaximum(1.45 * top);
          hb->SetFillColor(lightBlue); hb->SetLineColor(blue);
          hs->SetLineColor(red); hs->SetLineWidth(2);
          axisStyle(hb->GetXaxis(), bottomRow, kVars[i].name);
          axisStyle(hb->GetYaxis(), false, "");
          if (i == 0 && bottomRow) { hb->GetXaxis()->SetBinLabel(1, "none"); hb->GetXaxis()->SetBinLabel(2, "hit"); }
          hb->Draw("HIST");
          if (cuts[i].active) {
            for (auto [x1, x2] : {std::make_pair(kVars[i].lo, std::max(kVars[i].lo, cuts[i].lo)),
                                  std::make_pair(std::min(kVars[i].hi, cuts[i].hi), kVars[i].hi)}) {
              if (x2 <= x1) continue;
              auto box = new TBox(x1, 0., x2, top);
              box->SetFillStyle(3354); box->SetFillColor(kGray + 3); box->SetLineWidth(0);
              box->Draw();
            }
          }
          hs->Draw("HIST SAME");
          for (double x : bounds(i)) cutLine = line(x, 0., x, top);
          splitMark(i, true, 0.82);   // under the "cut keeps" line
          TLatex text;
          text.SetNDC(); text.SetTextSize(0.068);
          text.DrawLatex(0.23, 0.89, cuts[i].active ? Form("cut keeps S %.0f%%, B %.0f%%", 100. * sKept[row] / std::max(1e-12, sAll),
                                                           100. * bKept[row] / std::max(1e-12, bAll))
                                                    : "no cut");
          gPad->RedrawAxis();
        } else {            // background density and signal points: x the column's variable, y the row's
          TH2D* h = hB2[row][col];
          pad->SetLogz();
          axisStyle(h->GetXaxis(), bottomRow, kVars[j].name);
          axisStyle(h->GetYaxis(), col == 0, kVars[i].name);
          if (j == 0 && bottomRow) { h->GetXaxis()->SetBinLabel(1, "none"); h->GetXaxis()->SetBinLabel(2, "hit"); }
          h->SetMinimum(1e-3 * h->GetMaximum());   // three decades of background density
          h->Draw("COL");
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
    text.DrawLatex(0.02, 0.9, title);
    text.SetTextFont(42); text.SetTextSize(0.052);
    text.DrawLatex(0.02, 0.81, subtitle);
    auto marker = static_cast<TGraph*>(gS[1][0]->Clone("marker"));   // a legend marker big enough to see
    marker->SetMarkerSize(1.4);
    auto area = static_cast<TH1D*>(hB[0]->Clone("area"));
    auto legend = new TLegend(0.02, 0.46, 0.9, 0.77);
    legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(0.052);
    const char* unit = byEvent ? "events" : "entries";
    legend->AddEntry(marker, Form("%s with a muonic X-ray of the target (S)", unit), "p");
    legend->AddEntry(area, Form("other %s (B)", unit), "f");
    if (cutLine) legend->AddEntry(cutLine, cutsName, "l");
    legend->Draw();
    // as event_flow.C: expected in measureDays days, and the raw simulated numbers
    text.DrawLatex(0.02, 0.38, Form("in %.0f days: S %s, B %s", measureDays, thousands(measureDays * sAll).c_str(),
                                    thousands(measureDays * bAll).c_str()));
    text.SetTextSize(0.042);
    text.DrawLatex(0.02, 0.325, Form("simulated: S %s, B %s %s", thousands(nSig).c_str(), thousands(nBack).c_str(), unit));
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
    text.SetTextSize(0.045);
    for (std::size_t k = 0; k < drawn.size() && k < 4; ++k) text.DrawLatex(0.02, 0.27 - 0.055 * k, drawn[k].c_str());
    text.DrawLatex(0.02, 0.05, "diagonal: unit area, hatched = rejected; B log density below it");
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
  // the muon's cuts act on muonVars; its corner plots show muonShown, without the bottom-TOF hit and bar (0, 10)
  const std::vector<int> muonVars = usedOf({0, 10, 1, 11, 2, 3, 4});
  const std::vector<int> muonShown = usedOf({1, 11, 2, 3, 4}), muonHidden = usedOf({0});
  const std::vector<int> stackVars = usedOf(gCZT ? std::vector<int>{13, 5, 6, 9, 12, 7, 8} : std::vector<int>{13, 5, 6, 7, 8, 9});
  auto passesGroup = [&](const Entry& x, const std::vector<int>& group) {
    for (int i : group) if (!Passes(x, cuts, i)) return false;
    return true;
  };
  // the corner plots show the captured data, as event_flow.C counts it (asked for on 2026-11-02): the TOF-triggered
  // events with a CZT (stack) entry in time, at every energy (step 2a of event_flow.C), the cuts as chosen in the L_alpha
  // window drawn, each cut's share by itself (event_flow_cuts.png); one point per event, as event_flow.C counts them:
  // signal an event with a muonic X-ray of the target among its entries; the entry's variables those of its
  // representative entry, a muonic X-ray of the target first, then one that passes the entry's cuts (so that the event
  // passes them if any entry does), then one in the L_alpha window; else the first. Drawn once per group, under the names corner_muon_after and corner_stack_after
  // (since 2026-11-02: no plots after the other group's cuts, which no step of the flow shows)
  auto eventsOf = [&](std::function<bool(const Entry&)> keep) {
    std::vector<Entry> events;
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
      if (isNew) events.push_back(x);
      else if (better(x, events[it->second])) events[it->second] = x;
    }
    return events;
  };
  const std::vector<Entry> captured = eventsOf([](const Entry&) { return true; });
  const char* anEntry = gCZT ? "a CZT entry" : "a stack entry";
  const char* eventsName = gTOF ? "triggered events" : "events with a middle-TOF hit";
  const std::string capturedMuon = Form("%s with %s, all energies (step 2a)", eventsName, anEntry);
  const std::string capturedEntry = Form("%s with %s, all energies (step 2a; an entry each)", eventsName, anEntry);
  corner(captured, true, muonShown, [](const Entry&) { return true; }, "corner_muon_after", Form("The muon, %s", label),
         capturedMuon.c_str(), muonHidden);
  corner(captured, true, stackVars, [](const Entry&) { return true; }, "corner_stack_after", Form("The %s, %s", entryName, label),
         capturedEntry.c_str(), {});

  // 2. the cut flow
  {
    auto c = canvas("cutflow", true);
    c->SetBottomMargin(0.42);
    const int nSteps = flowS.size();
    auto axis = new TH1D("steps", ";;per day in the Al L-alpha window", nSteps, -0.5, nSteps - 0.5);
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
    note.DrawLatex(0.96, 0.8, Form("numbers: S/sqrt(B) in %.0f days after the step", measureDays));
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
  // the Al lines as simulated (the cascade model's energies; the measured Al K_alpha is 346.8 keV), named dx right of
  // them: L_alpha and L_beta, with kLines the K lines too (L_alpha and L_beta then named together)
  auto lineMarks = [&](double yLow, double yHigh, bool logy, bool kLines, double dx) {
    TLatex mark;
    mark.SetTextSize(0.03); mark.SetTextColor(kGray + 2); mark.SetTextAngle(90); mark.SetTextAlign(12);
    std::vector<std::pair<double, const char*>> marks = {{kAlpha, kLines ? "Al L_{#alpha}, L_{#beta}" : "Al L_{#alpha}"},
                                                         {kBeta, kLines ? "" : "Al L_{#beta}"}};
    if (kLines) { marks.push_back({354.9, "Al K_{#alpha}"}); marks.push_back({420.6, "Al K_{#beta}"}); }
    for (auto [eKeV, nameLine] : marks) {
      auto l = new TLine(eKeV, yLow, eKeV, yHigh);
      l->SetLineColor(kGray + 2); l->SetLineStyle(3);
      l->Draw();
      if (*nameLine) mark.DrawLatex(eKeV + dx, logy ? 1.3 * yHigh : yHigh, nameLine);
    }
  };
  // what a step holds, as event_flow.C counts it: each class's entries at eMin-eMax keV (the bins whose centre lies
  // there), in measureDays days and simulated, and the step's events at every energy (signal: a muonic X-ray of the
  // target among the entries that pass)
  struct StepCount { double total[3]; long raw[3]; double sDays, bDays; long sRaw, bRaw; };
  auto stepCount = [&](const std::array<TH1D*, 3>& h, const Cuts& cut, double eMin = 0., double eMax = 1e9) {
    StepCount n{};
    for (int k = 0; k < 3; ++k) {
      double sum = 0.;
      for (int b = 1; b <= h[k]->GetNbinsX(); ++b)
        if (h[k]->GetBinCenter(b) >= eMin && h[k]->GetBinCenter(b) < eMax) sum += h[k]->GetBinContent(b);
      n.total[k] = measureDays * sum;
    }
    std::map<long long, std::pair<double, bool>> events;   // per day, signal
    for (const auto& x : entries) {
      if (!PassesAll(x, cut)) continue;
      if (x.e >= eMin && x.e < eMax) ++n.raw[x.kind];
      auto& event = events[x.id];
      event.first = x.w;
      event.second = event.second || x.kind == 0;
    }
    for (const auto& [id, event] : events) {
      if (event.second) { n.sDays += measureDays * event.first; ++n.sRaw; }
      else { n.bDays += measureDays * event.first; ++n.bRaw; }
    }
    return n;
  };
  auto count = [&](double v) { return v < 100. ? std::string(Form("%.1f", v)) : thousands(v); };
  // the spectrum of the entries that pass the cut set c, per day and 2 keV at plotLow-plotHigh, as ROOT draws a histogram
  // (asked for on 2026-11-02): the simulated total as a line with its statistical error bars (the simulation's own), each
  // class as a thin line in its colour; with the step's events (at every energy) and each class's entries in the range
  // drawn, in measureDays days as event_flow.C counts them; mark: the L-line fit's range and its result
  auto spectrumPlot = [&](const std::array<TH1D*, 3>& h, const Cuts& cut, bool logy, bool mark, const char* name, const char* title) {
    auto c = canvas(name, logy);
    auto total = static_cast<TH1D*>(h[2]->Clone(Form("%s_total", name)));
    total->Add(h[1]); total->Add(h[0]);
    double top = 0., low = 1e30;
    for (int b = 1; b <= total->GetNbinsX(); ++b) {
      if (total->GetBinCenter(b) < plotLow || total->GetBinCenter(b) >= plotHigh) continue;
      top = std::max(top, total->GetBinContent(b) + total->GetBinError(b));
      for (int k = 0; k < 3; ++k) if (h[k]->GetBinContent(b) > 0.) low = std::min(low, h[k]->GetBinContent(b));
    }
    const StepCount n = stepCount(h, cut, plotLow, plotHigh);
    frame(plotLow, logy ? 0.5 * low : 0., plotHigh, logy ? 1000. * top : 1.6 * top,
          Form(";%s;entries per day per %.0f keV", energyTitle, total->GetBinWidth(1)));
    TBox* band = nullptr;
    if (mark) {   // the L-line fit's range, behind everything
      band = new TBox(fitLow, 0., fitHigh, 1.6 * top);
      band->SetFillColor(TColor::GetColor("#ececec")); band->SetLineColor(TColor::GetColor("#ececec"));
      band->Draw();
    }
    TH1D* line[3];
    for (int k = 2; k >= 0; --k) {
      line[k] = static_cast<TH1D*>(h[k]->Clone(Form("%s_line%d", name, k)));
      line[k]->SetLineColor(kindColour[k]); line[k]->SetLineWidth(1);
      line[k]->Draw("HIST SAME");
    }
    total->SetLineColor(kBlack); total->SetLineWidth(1); total->SetMarkerSize(0);   // thin: the classes show under it
    total->Draw("HIST SAME");
    total->Draw("E1 SAME");
    lineMarks(logy ? 0.5 * low : 0., logy ? 3. * top : 1.15 * top, logy, true, 0.01 * (plotHigh - plotLow));
    auto legend = legendAt(0.27, mark ? 0.5 : 0.6, 0.76, 0.89, 0.026);   // between the L lines and the K lines
    legend->SetHeader(Form("in %.0f days:", measureDays));
    legend->AddEntry((TObject*)nullptr, Form("events: S %s, B %s", count(n.sDays).c_str(), count(n.bDays).c_str()), "");
    legend->AddEntry((TObject*)nullptr, Form("(simulated: %s and %s)", thousands(n.sRaw).c_str(), thousands(n.bRaw).c_str()), "");
    legend->AddEntry((TObject*)nullptr, Form("entries at %.0f-%.0f keV (simulated):", plotLow, plotHigh), "");
    auto numbers = [&](double v, long r) { return count(v) + " (" + thousands(r) + ")"; };
    legend->AddEntry(total, Form("total: %s", numbers(n.total[0] + n.total[1] + n.total[2], n.raw[0] + n.raw[1] + n.raw[2]).c_str()), "le");
    for (int k = 0; k < 3; ++k) legend->AddEntry(line[k], Form("%s: %s", shortNames[k], numbers(n.total[k], n.raw[k]).c_str()), "l");
    if (mark) {
      legend->AddEntry(band, Form("fit of the L lines, %.0f-%.0f keV:", fitLow, fitHigh), "f");
      legend->AddEntry((TObject*)nullptr, Form("   L-alpha %.0f +/- %.0f entries, %.1f sigma", fitShown.amount, fitShown.error, fitShown.z), "");
    }
    legend->Draw();
    gPad->RedrawAxis();
    heading(title);
    save(c, name);
  };
  const auto beforeFull = spectrum(Cuts{}, "beforeFull", 500, 1000.), afterFull = spectrum(cuts, "afterFull", 500, 1000.);
  spectrumPlot(beforeFull, Cuts{}, true, false, "spectrum_before",
               Form("%s, before the cuts: %s", gTOF ? "Captured" : "Middle-TOF hit", label));
  spectrumPlot(afterFull, cuts, false, true, "spectrum_after", Form("After the corner-plot cuts: %s", label));
  {
    const StepCount all = stepCount(afterFull, cuts), drawn = stepCount(afterFull, cuts, plotLow, plotHigh);
    printf("\nafter the cuts, %.0f days: events S %.1f, B %.1f (simulated %ld, %ld)\n"
           "   entries at every energy: target muonic X-rays %.1f, other muonic X-rays %.1f, the rest %.1f (simulated %ld, %ld, %ld)\n"
           "   of these at %.0f-%.0f keV (drawn): %.1f, %.1f, %.1f (simulated %ld, %ld, %ld)\n", measureDays, all.sDays, all.bDays,
           all.sRaw, all.bRaw, all.total[0], all.total[1], all.total[2], all.raw[0], all.raw[1], all.raw[2], plotLow, plotHigh,
           drawn.total[0], drawn.total[1], drawn.total[2], drawn.raw[0], drawn.raw[1], drawn.raw[2]);
  }

  // the L-line fit of the expected spectrum after the cuts (step 3 of event_flow.C) over the fit's range, as ROOT draws a
  // fitted histogram (2026-11-02): the simulated total per day and 2 keV as a line with its statistical error bars (the
  // simulation's), the classes as thin lines, the fit as a red curve with its smooth background dashed, a box with the
  // fit's numbers (bigger than the other plots); in a panel below, the expectation minus the fitted smooth background,
  // with the fitted lines: the excess that the significance measures
  {
    const StepCount n = stepCount(after, cuts, fitLow, fitHigh);
    const int bLow = after[0]->FindBin(fitLow + 1e-6), bHigh = after[0]->FindBin(fitHigh - 1e-6);
    auto total = static_cast<TH1D*>(after[2]->Clone("fitTotal"));   // per day and 2 keV, with the simulation's errors
    total->Add(after[1]); total->Add(after[0]);
    auto c = new TCanvas("measurement_fit", "measurement_fit", 1600, 1600);
    auto upper = new TPad("upper", "", 0., 0.34, 1., 1.);
    upper->SetLeftMargin(0.12); upper->SetRightMargin(0.03); upper->SetTopMargin(0.12); upper->SetBottomMargin(0.015);
    auto lower = new TPad("lower", "", 0., 0., 1., 0.34);
    lower->SetLeftMargin(0.12); lower->SetRightMargin(0.03); lower->SetTopMargin(0.03); lower->SetBottomMargin(0.4);
    upper->Draw();
    lower->Draw();
    // the upper panel: the expectation, its classes, the fit
    upper->cd();
    double top = 0.;
    for (int b = bLow; b <= bHigh; ++b) top = std::max(top, total->GetBinContent(b) + total->GetBinError(b));
    auto background = new TF1("background", [=](double* x, double*) {
      const double u = (x[0] - 80.) / 50.;
      return std::exp(fKnown->GetParameter(0) + fKnown->GetParameter(1) * u + fKnown->GetParameter(2) * u * u +
                      fKnown->GetParameter(3) * u * u * u);
    }, fitLow, fitHigh, 0);
    auto perDay = [&](TF1* f, const char* name) {   // a curve of the fit (entries per 2 keV in measureDays days) per day
      return new TF1(name, [=](double* x, double*) { return f->Eval(x[0]) / measureDays; }, fitLow, fitHigh, 0);
    };
    TF1 *fitDay = perDay(fKnown, "fitPerDay"), *backgroundDay = perDay(background, "backgroundPerDay");
    const double curves = fitDay->GetMaximum(fitLow, fitHigh);
    auto upperFrame = gPad->DrawFrame(fitLow, 0., fitHigh, 2. * std::max(top, curves), ";;entries per day per 2 keV");
    upperFrame->GetXaxis()->SetLabelSize(0.);   // shared with the panel below
    upperFrame->GetYaxis()->SetTitleSize(0.052); upperFrame->GetYaxis()->SetLabelSize(0.045); upperFrame->GetYaxis()->SetTitleOffset(1.05);
    TH1D* line[3];
    for (int k = 2; k >= 0; --k) {
      line[k] = static_cast<TH1D*>(after[k]->Clone(Form("fitLine%d", k)));
      line[k]->SetLineColor(kindColour[k]); line[k]->SetLineWidth(2);
      line[k]->Draw("HIST SAME");
    }
    total->SetLineColor(kBlack); total->SetLineWidth(3); total->SetMarkerSize(0);
    total->Draw("HIST SAME");
    total->Draw("E1 SAME");
    drawDashed(backgroundDay, 14., 10., kRed, 4);
    fitDay->SetLineColor(kRed); fitDay->SetLineWidth(4); fitDay->SetNpx(2000);
    fitDay->Draw("SAME");
    lineMarks(0., 1.05 * std::max(top, curves), false, false, 2.);
    // the fit's numbers, in a box as ROOT shows a fit
    auto box = new TPaveText(0.62, 0.69, 0.96, 0.87, "NDC");
    box->SetBorderSize(1); box->SetFillColor(kWhite); box->SetTextAlign(12); box->SetTextFont(42); box->SetTextSize(0.03);
    box->AddText(Form("fit to the expectation, %.0f days", measureDays));
    box->AddText(Form("L_{#alpha}   %.1f #pm %.1f entries", fitShown.amount, fitShown.error));
    box->AddText(Form("L_{#beta}   %.1f #pm %.1f entries", fKnown->GetParameter(5), fKnown->GetParError(5)));
    box->AddText(Form("L_{#alpha} significance   %.1f #sigma", fitShown.z));
    box->AddText(Form("#chi^{2} / ndf   %.1f / %d", chi2Fit, ndfFit));
    box->Draw();
    auto dashed = [](int col, int width) { auto l = new TLine(); l->SetLineStyle(2); l->SetLineWidth(width); l->SetLineColor(col); return l; };
    auto legend = legendAt(0.57, 0.36, 0.96, 0.66, 0.03);
    legend->AddEntry(total, Form("simulated total: %s entries", count(n.total[0] + n.total[1] + n.total[2]).c_str()), "le");
    for (int k = 0; k < 3; ++k) legend->AddEntry(line[k], Form("%s: %s entries", shortNames[k], count(n.total[k]).c_str()), "l");
    legend->AddEntry(fitDay, "fit: L lines on a smooth background", "l");
    legend->AddEntry(dashed(kRed, 4), "its smooth background", "l");
    legend->Draw();
    TLatex text;
    text.SetNDC(); text.SetTextSize(0.032);
    text.DrawLatex(0.15, 0.84, Form("after the corner-plot cuts, %.0f days:", measureDays));
    text.DrawLatex(0.15, 0.795, Form("   %s target events, %s others", count(n.sDays).c_str(), count(n.bDays).c_str()));
    text.DrawLatex(0.15, 0.75, Form("   (their entries at %.0f-%.0f keV: the legend)", fitLow, fitHigh));
    gPad->RedrawAxis();
    // the lower panel: the expectation minus the fitted smooth background (the simulation's errors), and the fitted lines
    lower->cd();
    auto residual = new TH1D("fitResidual", "", bHigh - bLow + 1, after[0]->GetBinLowEdge(bLow), after[0]->GetBinLowEdge(bHigh + 1));
    residual->Sumw2();
    double rLow = 0., rHigh = 0.;
    for (int b = bLow; b <= bHigh; ++b) {
      const double r = total->GetBinContent(b) - backgroundDay->Eval(total->GetBinCenter(b)), error = total->GetBinError(b);
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
    auto lines = new TF1("lines", [=](double* x, double*) { return fitDay->Eval(x[0]) - backgroundDay->Eval(x[0]); }, fitLow, fitHigh, 0);
    lines->SetLineColor(kRed); lines->SetLineWidth(4); lines->SetNpx(2000);
    lines->Draw("SAME");
    TLatex panel;
    panel.SetNDC(); panel.SetTextSize(0.075);
    panel.DrawLatex(0.15, 0.85, "the expectation minus the fitted background; red: the fitted L lines");
    gPad->RedrawAxis();
    c->cd();
    TLatex title;
    title.SetNDC(); title.SetTextSize(0.036);
    title.DrawLatex(0.12, 0.945, Form("The Al L-line fit, %s", label));
    save(c, "measurement_fit", 0.0155);
  }
  printf("\nplots: %s/muonic_xray_select_{corner_muon_after,corner_stack_after,cutflow,spectrum_before,spectrum_after,"
         "measurement_fit}.png (and .pdf); table: %s/muonic_xray_select.csv\n", out.Data(), out.Data());
}
