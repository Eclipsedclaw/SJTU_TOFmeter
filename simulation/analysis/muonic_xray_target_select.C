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
// 1. Corner plots of ten candidate variables in the Al L_alpha window, five of the muon (bottom-TOF hit, middle-TOF
//    bar, middle-TOF energy, top-to-middle time, top-TOF energy) and five of the stack entry (highest and lowest layer
//    hit, hits in the entry, largest hit / entry, outermost hit), signal as red points, background as a blue density.
// 2. The most efficient cuts: greedy, one cut at a time the variable and threshold that raise S/sqrt(B) in the L_alpha
//    window the most (by at least 1%), then each threshold again with the others fixed; chosen on the even events and
//    scored on the odd ones and back, to see the optimisation's own gain. The cut flow, and the corner plots again,
//    each group after the other group's cuts.
// 3. The spectrum before and after the cuts by what made it, and simulated measurements of measureDays days drawn from
//    it (Poisson), fitted without truth: smooth background and Gaussian Al L_alpha and L_beta, either with free
//    position and width (background exp(pol2), as muonic_xray_target_truth.C fits the target's X-rays alone) or at the
//    known line energies with the stack's known width (background exp(pol3), only the amounts free), with the
//    likelihood-ratio significance of L_alpha; the median of `toys` measurements, one of them shown. A null test: the
//    same fits to measurements without any muonic X-ray, where an honest fit finds no L_alpha. The cut sets of the
//    1 cm target and the earlier own cuts of muonic_xray_target_cuts.C, compared with the known-energy fit.
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
// hit (10, 0 if none); with single pixels the pixel's number 1-16 (12). 10 and 12 are in the corner plots but not cut on:
// a set of bars or pixels is not a range, and their useful parts are variables 0 (no bottom hit) and 9 (central pixels).

#include <TParameter.h>

#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
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
// the pixel, single CZT pixels. 0-4, 10 and 11 are the muon's, the others the entry's
const int kNV = 13;
struct Variable { const char* name; int bins; double lo, hi; bool integer; };
bool Used(int i)
{
  if (i == 10) return gBottomBar;
  if (i == 11) return gTOF;
  if (i == 12) return gCZT && gPixel;
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
                             {"CZT pixel (1-16)", 16, 0.5, 16.5, true}};

struct Entry {
  double w;        // per day of its sample's exposure
  int kind, half;  // 0 muonic X-ray of the target (signal), 1 muonic X-ray of another volume, 2 the rest; event parity
  int sample;      // 0 the target's run (muons up to eSplit), 1 the 3-week run's faster muons
  bool trig;       // the TOF trigger fired (always true for skims without the measured TOF)
  bool timed;      // the top TOF timed the muon: variables 3 and 4 exist (and with the measured TOF the track, 1)
  double e;        // add-back energy with the response (keV)
  double v[kNV];
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
    default: break;   // 10 and 12, the bottom-TOF bar and the CZT pixel: shown, not cut (a set of them is no range)
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
    kVars[3] = {"TOF kinetic energy (MeV)", 40, 0., 200., false};
  }
  if (gCZT && gPixel) {   // a single pixel: the pixels hit in the event and their energy, where the pixel is
    energyTitle = "pixel energy with the response (keV)";
    kVars[5] = {"pixels hit in the event", 6, 0.5, 6.5, true};
    kVars[6] = {"other pixels' energy (MeV)", 40, 0., 10., false};
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
                            trig ? std::atan(std::hypot(slopeX, slopeY)) * 180. / M_PI : -1., pix + 1.}});
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
                          double(lowest), double(n), e > 0. ? emax / e : 0., rmax, 0., 0., 0.}});
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
                                                   "corner; shown, not cut)" : ""));
  }
  conditions.push_back(Form("  cuts: greedy on S/sqrt(B) in %g days in the L_alpha window, each step at least 1%% better and keeping "
                            "at least 30 signal entries (60 on all events), chosen on even and odd events alone to cross-validate",
                            measureDays));
  conditions.push_back(Form("  fits: %d simulated measurements of %g days (Poisson), %g-%g keV, seed %u", toys, measureDays, fitLow,
                            fitHigh, seed));
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
  std::vector<int> order0, order1, order;
  const Cuts cuts0 = optimise(0, order0), cuts1 = optimise(1, order1), cuts = optimise(2, order);
  const double zCv = 0.5 * (score(cuts0, 1) + score(cuts1, 0));
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
  list("all events:", cuts, order);
  list("even events:", cuts0, order0);
  list("odd events:", cuts1, order1);
  printf("   S/sqrt(B) on the events they were chosen on %.2f; cross-validated (even's cuts on the odd events and back) %.2f\n",
         score(cuts, 2), zCv);

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
  auto spectrum = [&](const Cuts& c, const char* tag) {
    std::array<TH1D*, 3> h;
    for (int k = 0; k < 3; ++k) {
      h[k] = new TH1D(Form("spectrum_%s_%d", tag, k), "", 100, 0., 200.);
      h[k]->Sumw2();
    }
    for (const auto& x : entries) if (PassesAll(x, c)) h[x.kind]->Fill(x.e, x.w);
    return h;
  };
  const auto before = spectrum(Cuts{}, "before"), after = spectrum(cuts, "after");

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

  // ---- simulated measurements: Poisson counts of measureDays days, fitted without truth, in two ways:
  //  free:  smooth background exp(pol2), Gaussian L_alpha and L_beta with free position and width (as
  //         muonic_xray_target_truth.C fits the target's X-rays alone);
  //  known: background exp(pol3), the Gaussian lines at their known energies with the stack's known width, only their
  //         amounts free (as a calibrated measurement fits them; an amount may come out negative);
  //  shape: as known, but the background's shape given (the simulated spectrum of the rest after the cuts), only its
  //         size free: the best case, of a background shape measured apart (e.g. a long run without the target);
  //  control: as shape, but the shape that of the same measurement's control region, the entries that pass every cut
  //         but the middle-TOF energy one and fail that one (nearly all fast muons): counted with the measurement,
  //         noise, signal leaking in and a shape a little off included.
  // The L_alpha significance: the likelihood ratio against the same fit without L_alpha, with the sign of its amount.
  // A null test: measurements without any muonic X-ray, where an honest fit finds no L_alpha. Fit range: fitLow-fitHigh.
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
  // the background shapes of the last two fits: the rest after the cut set's cuts per day, a measured control region
  Cuts controlCuts = cuts;
  controlCuts[2].active = true; controlCuts[2].lo = -1e9; controlCuts[2].hi = cuts[2].active ? cuts[2].lo : 1.75;
  const auto control = spectrum(controlCuts, "control");
  TH1D* shape = after[2];
  TH1D* controlToy = nullptr;
  TH1D** shapeNow = &shape;
  TH1D** controlNow = &controlToy;
  auto shapeModel = [&](const char* name, bool fromControl) {
    return new TF1(name, [=](double* x, double* p) {
      const double back = fromControl ? (*controlNow)->GetBinContent((*controlNow)->FindBin(x[0]))
                                      : measureDays * (*shapeNow)->GetBinContent((*shapeNow)->FindBin(x[0]));
      return p[0] * back + 2. * (p[1] * gaussian(x[0], kAlpha, SigmaKeV(kAlpha)) + p[2] * gaussian(x[0], kBeta, SigmaKeV(kBeta)));
    }, fitLow, fitHigh, 3);
  };
  TF1 *fFree = freeModel("free"), *fFreeNull = freeModel("freeNull"), *fKnown = knownModel("known"),
      *fKnownNull = knownModel("knownNull"), *fShape = shapeModel("shape", false), *fShapeNull = shapeModel("shapeNull", false),
      *fControl = shapeModel("control", true), *fControlNull = shapeModel("controlNull", true);
  fFreeNull->FixParameter(3, 0.); fFreeNull->FixParameter(4, kAlpha); fFreeNull->FixParameter(5, SigmaKeV(kAlpha));
  fKnownNull->FixParameter(4, 0.);
  fShapeNull->FixParameter(1, 0.);
  fControlNull->FixParameter(1, 0.);
  auto nll = [&](TH1D* h, TF1* f) {   // Poisson -log L of the fit range, at the bin centres as the fit
    double sum = 0.;
    for (int b = h->FindBin(fitLow + 1e-6); b <= h->FindBin(fitHigh - 1e-6); ++b) {
      const double mu = std::max(1e-12, f->Eval(h->GetBinCenter(b)));
      sum += mu - h->GetBinContent(b) * std::log(mu);
    }
    return sum;
  };
  struct LineFit { double z, amount, error, position, width; };
  auto fitMeasurement = [&](TH1D* h, int how) {   // 0 free, 1 known, 2 shape, 3 control
    double window = 0.;
    for (int b = h->FindBin(kAlpha - 2. * SigmaKeV(kAlpha)); b <= h->FindBin(kAlpha + 2. * SigmaKeV(kAlpha)); ++b) window += h->GetBinContent(b);
    const double back = std::log(std::max(1., h->GetBinContent(h->FindBin(80.))));
    if (how >= 2) {
      TF1 *f = how == 2 ? fShape : fControl, *f0 = how == 2 ? fShapeNull : fControlNull;
      const double ref = how == 2 ? measureDays * shape->GetBinContent(shape->FindBin(80.))
                                  : controlToy->GetBinContent(controlToy->FindBin(80.));
      f->SetParameters(h->GetBinContent(h->FindBin(80.)) / std::max(1e-9, ref), 0.1 * window, 0.02 * window);
      f->SetParLimits(1, -0.3 * window, 2. * window); f->SetParLimits(2, -0.3 * window, 2. * window);
      h->Fit(f, "LQNR");
      f0->SetParameter(0, f->GetParameter(0)); f0->SetParameter(2, f->GetParameter(2));
      f0->SetParLimits(2, -0.3 * window, 2. * window);
      h->Fit(f0, "LQNR");
      const double z = std::sqrt(std::max(0., 2. * (nll(h, f0) - nll(h, f))));
      return LineFit{f->GetParameter(1) < 0. ? -z : z, f->GetParameter(1), f->GetParError(1), kAlpha, SigmaKeV(kAlpha)};
    }
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
  TRandom3 toyRandom(seed);
  // a simulated measurement: Poisson counts of measureDays days from the spectrum h, with or without its muonic X-rays
  auto measurement = [&](const std::array<TH1D*, 3>& h, bool withXrays, const char* name) {
    auto m = new TH1D(name, "", 100, 0., 200.);
    for (int b = 1; b <= m->GetNbinsX(); ++b) {
      const double mu = measureDays * (h[2]->GetBinContent(b) + (withXrays ? h[0]->GetBinContent(b) + h[1]->GetBinContent(b) : 0.));
      m->SetBinContent(b, toyRandom.Poisson(mu));
      m->SetBinError(b, std::sqrt(m->GetBinContent(b)));
    }
    return m;
  };
  // over `toys` measurements: the median and 16-84% of the L_alpha significance, the fitted L_alpha amounts
  struct Toys { double z50, z16, z84, mean, rms; std::vector<double> amounts; };
  auto runToys = [&](const std::array<TH1D*, 3>& h, bool withXrays, int how, const char* tag) {
    Toys t;
    std::vector<double> z;
    for (int k = 0; k < toys; ++k) {
      TH1D* m = measurement(h, withXrays, Form("toy_%s_%d", tag, k));
      if (how == 3) controlToy = measurement(control, withXrays, Form("control_%s_%d", tag, k));
      const LineFit r = fitMeasurement(m, how);
      if (how == 3) { delete controlToy; controlToy = nullptr; }
      z.push_back(r.z);
      t.amounts.push_back(r.amount);
      delete m;
    }
    std::sort(z.begin(), z.end());
    t.z50 = z[toys / 2]; t.z16 = z[std::size_t(0.16 * toys)]; t.z84 = z[std::size_t(0.84 * toys)];
    double sum = 0., sum2 = 0.;
    for (double a : t.amounts) { sum += a; sum2 += a * a; }
    t.mean = sum / toys;
    t.rms = std::sqrt(std::max(0., sum2 / toys - t.mean * t.mean));
    return t;
  };
  double truthAlpha[3] = {}, truthFit[3] = {};   // muonic X-ray entries after the cuts in measureDays days: window, fit range
  for (const auto& x : entries) {
    if (!PassesAll(x, cuts)) continue;
    if (InAlpha(x.e)) truthAlpha[x.kind] += measureDays * x.w;
    if (x.e >= fitLow && x.e < fitHigh) truthFit[x.kind] += measureDays * x.w;
  }
  printf("\nfits without truth to %d simulated measurements of %.0f days each, after the most efficient cuts:\n", toys, measureDays);
  printf("   truth after the cuts: L_alpha window: target X-rays %.0f, other muonic X-rays %.0f, the rest %.0f; %.0f-%.0f keV: "
         "target X-rays %.0f, other muonic X-rays %.0f\n", truthAlpha[0], truthAlpha[1], truthAlpha[2], fitLow, fitHigh, truthFit[0],
         truthFit[1]);
  printf("   %-58s %22s %24s\n", "", "L_alpha significance", "fitted L_alpha amount");
  // the control region: what it holds, and how close its shape is to that of the background after the cuts
  {
    double inAlpha[3] = {};
    for (const auto& x : alpha) if (PassesAll(x, controlCuts)) inAlpha[x.kind] += x.w;
    double sa = 0., sb = 0.;
    for (int b = after[2]->FindBin(fitLow + 1e-6); b <= after[2]->FindBin(fitHigh - 1e-6); ++b) {
      sa += after[2]->GetBinContent(b);
      sb += control[0]->GetBinContent(b) + control[1]->GetBinContent(b) + control[2]->GetBinContent(b);
    }
    double chi2 = 0.;
    int bins = 0;
    for (int b = after[2]->FindBin(fitLow + 1e-6); b <= after[2]->FindBin(fitHigh - 1e-6); ++b, ++bins) {
      const double a = after[2]->GetBinContent(b) / sa, ea = after[2]->GetBinError(b) / sa;
      double c = 0., ec2 = 0.;
      for (int k = 0; k < 3; ++k) { c += control[k]->GetBinContent(b); ec2 += std::pow(control[k]->GetBinError(b), 2); }
      chi2 += std::pow(a - c / sb, 2) / (ea * ea + ec2 / (sb * sb));
    }
    printf("\nthe control region (%s): L_alpha window per day: target X-rays %.1f, other muonic X-rays %.1f, the rest %.1f;\n"
           "   its %.0f-%.0f keV shape against the background's after the cuts: chi2/ndf %.0f/%d\n",
           Form("middle-TOF energy < %.2f MeV, the other cuts as chosen", controlCuts[2].hi),
           inAlpha[0], inAlpha[1], inAlpha[2], fitLow, fitHigh, chi2, bins - 1);
  }
  const char* modelNames[4] = {"free position and width, background exp(pol2)", "known energies and width, background exp(pol3)",
                               "known energies and width, background shape known", "known energies and width, control-region shape"};
  Toys withLines[4], nullTest[4];
  for (int how = 0; how < 4; ++how) {
    withLines[how] = runToys(after, true, how, Form("lines%d", how));
    nullTest[how] = runToys(after, false, how, Form("null%d", how));
    for (int n = 0; n < 2; ++n) {
      const Toys& t = n ? nullTest[how] : withLines[how];
      printf("   %-47s %-10s %6.1f (%5.1f to %5.1f)   %7.0f +- %5.0f (spread)\n", modelNames[how], n ? "null test" : "all",
             t.z50, t.z16, t.z84, t.mean, t.rms);
      csv << (n ? "null" : "lines") << "," << how << ",\"" << modelNames[how] << "\",,,," << t.z50 << "," << t.z16 << ","
          << t.z84 << "\n";
    }
  }
  printf("\nthe cut sets compared: S/sqrt(B) in the L_alpha window, and the L_alpha significance of the fits with the known\n"
         "energies and width (median, 16-84%% of %d measurements of %.0f days), background exp(pol3) or of known shape:\n", toys, measureDays);
  printf("   %-46s %8s %8s %10s %22s %22s\n", "", "S/day", "B/day", "S/sqrt(B)", "exp(pol3)", "shape known");
  struct Set { const char* name; const Cuts* cuts; };
  const Set sets[3] = {{"cuts of the 1 cm target", &cuts1cm}, {"own cuts of muonic_xray_target_cuts.C", &cutsOwn},
                       {"the most efficient cuts (this optimisation)", &cuts}};
  // the earlier sets cut on the middle-TOF bar and the time of flight: with the measured TOF only this optimisation's
  for (int k = gTOF ? 2 : 0; k < 3; ++k) {
    double s, b;
    const double z = score(*sets[k].cuts, 2, &s, &b);
    const auto spec = spectrum(*sets[k].cuts, Form("set%d", k));
    const Toys t = k == 2 ? withLines[1] : runToys(spec, true, 1, Form("set%d", k));
    shape = spec[2];
    const Toys u = k == 2 ? withLines[2] : runToys(spec, true, 2, Form("shape%d", k));
    shape = after[2];
    printf("   %-46s %8.2f %8.1f %10.2f      %5.1f (%4.1f-%4.1f)      %5.1f (%4.1f-%4.1f)\n", sets[k].name, s, b, z, t.z50, t.z16,
           t.z84, u.z50, u.z16, u.z84);
    csv << "compare," << k << ",\"" << sets[k].name << "\"," << s << "," << b << "," << z << "," << t.z50 << "," << t.z16 << ","
        << t.z84 << "\n";
    csv << "compare_shape," << k << ",\"" << sets[k].name << "\"," << s << "," << b << "," << z << "," << u.z50 << "," << u.z16
        << "," << u.z84 << "\n";
  }
  // the measurement shown: the most efficient cuts, the fit with the known energies and width
  TH1D* shown = measurement(after, true, "measurement");
  TH1D* shownControl = measurement(control, true, "measurementControl");   // its control region, counted with it
  controlToy = shownControl;
  const LineFit fitShownControl = fitMeasurement(shown, 3);
  const LineFit fitShown = fitMeasurement(shown, 1);
  printf("\nthe measurement shown (%.0f days): L_alpha %.0f +- %.0f, %.1f sigma (background exp(pol3)); %.0f +- %.0f, %.1f sigma "
         "(control-region shape)\n", measureDays, fitShown.amount, fitShown.error, fitShown.z, fitShownControl.amount,
         fitShownControl.error, fitShownControl.z);

  // ---- plots, one per file (CLAUDE.md rule 8)
  gStyle->SetOptStat(0);
  const int red = kRed + 1, blue = TColor::GetColor("#2171b5"), lightBlue = TColor::GetColor("#9ecae1");
  const int kindColour[3] = {red, lightBlue, blue};
  const int kFitStyle = 11;
  gStyle->SetLineStyleString(kFitStyle, "40 20");   // the legend's sample of the dashed fit, as drawDashed draws it
  // the conditions in two lines at the foot of a plot (size: of the canvas height)
  auto footer = [&](TVirtualPad* c, double size, double gap) {
    c->cd();
    TLatex text;
    text.SetNDC(); text.SetTextFont(42); text.SetTextSize(size); text.SetTextColor(kGray + 2); text.SetTextAlign(11);
    text.DrawLatex(0.01, 0.006 + gap, stamp[0].c_str());
    text.DrawLatex(0.01, 0.006, stamp[1].c_str());
  };
  auto save = [&](TCanvas* c, const char* name) {
    footer(c, 0.02, 0.026);
    c->Print(out + Form("/muonic_xray_select_%s.pdf", name));
    c->SaveAs(out + Form("/muonic_xray_select_%s.png", name));
    delete c;
  };
  auto canvas = [](const char* name, bool logy) {
    auto c = new TCanvas(name, name, 900, 700);
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
    auto yPad = [](double y) { return gPad->GetLogy() ? std::log10(std::max(y, 1e-30)) : y; };
    TGraph* dash = nullptr;
    double length = 0., uLast = curve->GetXmin(), yLast = curve->Eval(uLast);
    for (double u = curve->GetXmin(); u <= curve->GetXmax() + 1e-9; u += 0.1) {
      const double y = curve->Eval(u);
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
  auto corner = [&](const std::vector<int>& vars, std::function<bool(const Entry&)> keep, const char* name, const char* title,
                    const char* subtitle) {
    const int nV = vars.size();
    std::vector<TH1D*> hS(nV), hB(nV);
    std::vector<std::vector<TH2D*>> hB2(nV, std::vector<TH2D*>(nV, nullptr));
    std::vector<std::vector<TGraph*>> gS(nV, std::vector<TGraph*>(nV, nullptr));
    for (int a = 0; a < nV; ++a) {
      const Variable& va = kVars[vars[a]];
      hS[a] = new TH1D(Form("%s_s%d", name, a), "", va.bins, va.lo, va.hi);
      hB[a] = new TH1D(Form("%s_b%d", name, a), "", va.bins, va.lo, va.hi);
      for (int c = 0; c < a; ++c) {
        const Variable& vc = kVars[vars[c]];
        hB2[a][c] = new TH2D(Form("%s_b%d%d", name, a, c), "", vc.bins, vc.lo, vc.hi, va.bins, va.lo, va.hi);
        gS[a][c] = new TGraph();
      }
    }
    double sAll = 0., bAll = 0.;
    std::vector<double> sKept(nV, 0.), bKept(nV, 0.);
    long nSig = 0, nBack = 0;
    TRandom3 spread(7);
    auto inRange = [](double v, const Variable& var) { return std::min(std::max(v, var.lo + 1e-6), var.hi - 1e-6); };
    for (const auto& x : alpha) {
      if (!keep(x)) continue;
      const bool sig = x.kind == 0;
      (sig ? sAll : bAll) += x.w;
      ++(sig ? nSig : nBack);
      for (int a = 0; a < nV; ++a) {
        const int i = vars[a];
        if (Passes(x, cuts, i)) (sig ? sKept : bKept)[a] += x.w;
        if (!Defined(x, i)) continue;
        const double va = inRange(x.v[i], kVars[i]);
        (sig ? hS : hB)[a]->Fill(va, x.w);
        for (int c = 0; c < a; ++c) {
          const int j = vars[c];
          if (!Defined(x, j)) continue;
          const double vc = inRange(x.v[j], kVars[j]);
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
    TLine* cutLine = nullptr;
    for (int row = 0; row < nV; ++row)
      for (int col = 0; col <= row; ++col) {
        c.cd();
        auto pad = new TPad(Form("pad%d%d", row, col), "", col * cell, 1. - (row + 1) * cellY, (col + 1) * cell, 1. - row * cellY);
        pad->SetLeftMargin(0.2); pad->SetBottomMargin(0.2); pad->SetRightMargin(0.03); pad->SetTopMargin(0.03);
        pad->Draw();
        pad->cd();
        const int i = vars[row], j = vars[col];
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
    legend->AddEntry(marker, Form("muonic X-rays of the target (S): %ld entries", nSig), "p");
    legend->AddEntry(area, Form("everything else (B): %ld entries", nBack), "f");
    if (cutLine) legend->AddEntry(cutLine, "the most efficient cuts", "l");
    legend->Draw();
    text.DrawLatex(0.02, 0.36, Form("per day: S %.1f, B %.0f", sAll, bAll));
    std::vector<std::string> drawn = {"cuts drawn:"};   // wrapped at the commas
    for (int i : vars) {
      if (!cuts[i].active) continue;
      const std::string next = Describe(cuts, i);
      if (drawn.back().size() + next.size() > 58) drawn.push_back(next);
      else drawn.back() += (drawn.back().back() == ':' ? " " : ", ") + next;
    }
    if (drawn.size() == 1 && drawn[0] == "cuts drawn:") drawn[0] += " none on these variables";
    text.SetTextSize(0.045);
    for (std::size_t k = 0; k < drawn.size() && k < 3; ++k) text.DrawLatex(0.02, 0.28 - 0.07 * k, drawn[k].c_str());
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
  const std::vector<int> muonVars = usedOf({0, 10, 1, 11, 2, 3, 4});
  const std::vector<int> stackVars = usedOf(gCZT ? std::vector<int>{5, 6, 9, 12, 7, 8} : std::vector<int>{5, 6, 7, 8, 9});
  const std::string beforeText = gTOF ? "Al L-alpha window, TOF trigger, before the cuts"
                                      : "Al L-alpha window, middle-TOF hit in coincidence, before the cuts";
  const char* beforeCuts = beforeText.c_str();
  auto passesGroup = [&](const Entry& x, const std::vector<int>& group) {
    for (int i : group) if (!Passes(x, cuts, i)) return false;
    return true;
  };
  corner(muonVars, [](const Entry&) { return true; }, "corner_muon",
         Form("The muon, %s", label), beforeCuts);
  corner(stackVars, [](const Entry&) { return true; }, "corner_stack",
         Form("The %s, %s", entryName, label), beforeCuts);
  corner(muonVars, [&](const Entry& x) { return passesGroup(x, stackVars); }, "corner_muon_after",
         Form("The muon after the %s cuts", gCZT ? "CZT" : "stack"), Form("Al L-alpha window, entries that pass the cuts on the %s", entryName));
  corner(stackVars, [&](const Entry& x) { return passesGroup(x, muonVars); }, "corner_stack_after",
         Form("The %s after the muon cuts", entryName), "Al L-alpha window, entries that pass the cuts on the muon");

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
    heading(Form("Cut flow of the most efficient cuts, %s", label));
    save(c, "cutflow");
  }

  // 3. the spectra by what made them: the rest (blue) at the bottom, other muonic X-rays (light blue), the target's red
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
      if (legend) legend->AddEntry(layer[k], Form("%s: %.0f entries", kindNames[k], h[k]->GetEntries()), "f");
    }
    return all;
  };
  auto lineMarks = [&](double yLow, double yHigh, bool logy) {
    TLatex mark;
    mark.SetTextSize(0.03); mark.SetTextColor(kGray + 2); mark.SetTextAngle(90); mark.SetTextAlign(12);
    for (auto [eKeV, nameLine] : {std::make_pair(kAlpha, "Al L_{#alpha}"), std::make_pair(kBeta, "Al L_{#beta}")}) {
      auto l = new TLine(eKeV, yLow, eKeV, yHigh);
      l->SetLineColor(kGray + 2); l->SetLineStyle(3);
      l->Draw();
      mark.DrawLatex(eKeV + 3., logy ? 1.3 * yHigh : yHigh, nameLine);
    }
  };
  auto spectrumPlot = [&](const std::array<TH1D*, 3>& h, bool logy, const char* name, const char* title) {
    auto c = canvas(name, logy);
    double top = 0., low = 1e30;
    for (int b = 1; b <= h[0]->GetNbinsX(); ++b) {
      const double sum = h[0]->GetBinContent(b) + h[1]->GetBinContent(b) + h[2]->GetBinContent(b);
      top = std::max(top, sum);
      if (h[0]->GetBinContent(b) > 0.) low = std::min(low, h[0]->GetBinContent(b));
    }
    frame(0., logy ? 0.5 * low : 0., 200., logy ? 60. * top : 1.45 * top, Form(";%s;entries per day per 2 keV", energyTitle));
    auto legend = legendAt(0.38, 0.72, 0.97, 0.89, 0.032);
    drawStack(h, 1., legend, name);
    lineMarks(logy ? 0.5 * low : 0., logy ? 3. * top : 1.15 * top, logy);
    legend->Draw();
    gPad->RedrawAxis();
    heading(title);
    save(c, name);
  };
  spectrumPlot(before, true, "spectrum_before", Form("%s, %s (stacked)", gTOF ? "TOF trigger, before the cuts" : "Before the cuts", label));
  spectrumPlot(after, false, "spectrum_after", Form("After the most efficient cuts, %s (stacked)", label));

  // a simulated measurement after the cuts, fitted without truth with the known energies and width
  {
    fitMeasurement(shown, 1);   // again, so that fKnown holds this measurement's fit
    auto c = canvas("measurement", false);
    double top = 0.;
    for (int b = 1; b <= shown->GetNbinsX(); ++b) top = std::max(top, shown->GetBinContent(b) + shown->GetBinError(b));
    frame(0., 0., 200., 1.85 * top, Form(";%s;entries per 2 keV in %.0f days", energyTitle, measureDays));
    drawStack(after, measureDays, nullptr, "expected");
    shown->SetMarkerStyle(20); shown->SetMarkerSize(0.8); shown->SetMarkerColor(kBlack); shown->SetLineColor(kBlack);
    shown->Draw("E1 SAME");
    auto background = new TF1("background", [=](double* x, double*) {
      const double u = (x[0] - 80.) / 50.;
      return std::exp(fKnown->GetParameter(0) + fKnown->GetParameter(1) * u + fKnown->GetParameter(2) * u * u +
                      fKnown->GetParameter(3) * u * u * u);
    }, fitLow, fitHigh, 0);
    drawDashed(fKnown, 12., 6., kBlack, 3);
    drawDashed(background, 3., 5., kGray + 2, 2);
    const LineFit fitShownShape = fitMeasurement(shown, 2);   // fShape: the fit with the background shape known
    drawDashed(fShape, 12., 6., kAzure + 2, 3);
    lineMarks(0., 1.05 * top, false);
    auto dashed = new TLine(); dashed->SetLineStyle(kFitStyle); dashed->SetLineWidth(3);
    auto dotted = new TLine(); dotted->SetLineStyle(3); dotted->SetLineWidth(2); dotted->SetLineColor(kGray + 2);
    auto fill = [](int col) { auto b = new TH1D(); b->SetFillColor(col); b->SetLineColor(col); return b; };
    auto legend = legendAt(0.57, 0.6, 0.97, 0.89, 0.028);
    legend->AddEntry(shown, Form("simulated measurement: %.0f entries", shown->Integral()), "pe");
    auto dashedShape = new TLine(); dashedShape->SetLineStyle(kFitStyle); dashedShape->SetLineWidth(3); dashedShape->SetLineColor(kAzure + 2);
    legend->AddEntry(dashed, "fit, smooth background", "l");
    legend->AddEntry(dotted, "its smooth background", "l");
    legend->AddEntry(dashedShape, "fit, background shape known", "l");
    legend->AddEntry(fill(red), "expected: muonic X-rays of the target", "f");
    legend->AddEntry(fill(lightBlue), "expected: other muonic X-rays", "f");
    legend->AddEntry(fill(blue), "expected: the rest", "f");
    legend->Draw();
    TLatex text;
    text.SetNDC(); text.SetTextSize(0.028);
    text.DrawLatex(0.14, 0.85, "Gaussian L lines of known energy and width");
    text.DrawLatex(0.14, 0.80, Form("smooth: L-alpha %.0f +/- %.0f, %.1f sigma", fitShown.amount, fitShown.error, fitShown.z));
    text.DrawLatex(0.14, 0.76, Form("   median of %d measurements %.1f sigma", toys, withLines[1].z50));
    text.SetTextColor(kAzure + 2);
    text.DrawLatex(0.14, 0.71, Form("shape known: L-alpha %.0f +/- %.0f, %.1f sigma", fitShownShape.amount, fitShownShape.error,
                                    fitShownShape.z));
    text.DrawLatex(0.14, 0.67, Form("   median of %d measurements %.1f sigma", toys, withLines[2].z50));
    gPad->RedrawAxis();
    heading(Form("A simulated %.0f-day measurement, %s", measureDays, label));
    save(c, "measurement");
  }

  // the null test: the fitted L_alpha amounts of measurements without any muonic X-ray, for both fits
  {
    double lo = 0., hi = 0.;
    for (int how = 0; how < 4; ++how) {
      for (double a : nullTest[how].amounts) { lo = std::min(lo, a); hi = std::max(hi, a); }
      hi = std::max(hi, withLines[how].mean);
    }
    const double margin = 0.1 * (hi - lo);
    lo -= margin; hi += margin;
    const int colours[4] = {kOrange + 7, kBlack, kAzure + 2, kGreen + 2};
    TH1D* h[4];
    double hMax = 0.;
    for (int how = 0; how < 4; ++how) {
      h[how] = new TH1D(Form("null%d", how), "", 50, lo, hi);
      for (double a : nullTest[how].amounts) h[how]->Fill(a);
      hMax = std::max(hMax, h[how]->GetMaximum());
    }
    auto c = canvas("null_test", false);
    frame(lo, 0., hi, 1.9 * hMax, Form(";fitted L-alpha amount (entries in %.0f days);simulated measurements", measureDays));
    auto legend = legendAt(0.14, 0.64, 0.97, 0.89, 0.027);
    for (int how = 0; how < 4; ++how) {
      h[how]->SetLineColor(colours[how]); h[how]->SetLineWidth(2);
      h[how]->Draw("HIST SAME");
      auto l = new TLine(withLines[how].mean, 0., withLines[how].mean, 1.1 * hMax);
      l->SetLineColor(colours[how]); l->SetLineStyle(2); l->SetLineWidth(2);
      l->Draw();
      legend->AddEntry(h[how], Form("%s: %.0f +/- %.0f", modelNames[how], nullTest[how].mean, nullTest[how].rms), "l");
    }
    auto sample = new TLine(); sample->SetLineStyle(2); sample->SetLineWidth(2);
    legend->AddEntry(sample, "dashed: the mean fitted amount with the muonic X-rays in", "l");
    legend->Draw();
    auto zero = new TLine(0., 0., 0., 1.1 * hMax);
    zero->SetLineColor(kGray + 1); zero->SetLineStyle(3);
    zero->Draw();
    heading(Form("Null test: %d measurements without muonic X-rays", toys));
    save(c, "null_test");
  }
  printf("\nplots: %s/muonic_xray_select_{corner_muon,corner_stack,corner_muon_after,corner_stack_after,cutflow,"
         "spectrum_before,spectrum_after,measurement,null_test}.png (and .pdf); table: %s/muonic_xray_select.csv\n", out.Data(), out.Data());
}
