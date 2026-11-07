// Cut optimisation for a cut profile (2026-11-05, asked for: the LYSO panel's own cuts in place of the CZT plane's): the
// counting significance of the target's L_alpha peak, Z_A = sqrt(2((S+B) ln(1+S/B) - S)), S the full-energy L_alpha
// entries (etrue within 0.3 keV of the line) in the L_alpha window (+-2 sigma of the skims' own resolution, their
// sigma1MeV), B every other entry there; counted in the slow band's days as muonic_xray_target_select.C and event_flow.C
// count them (the TOF trigger, a pixel entry within -5..+30 ns of the middle TOF, the slow band as simulated, the fast
// band scaled to its days and only for its background above 1 GeV).
// The lessons of the CZT optimisation (2026-11-02, the work folder's output/muonic_xray_czt/pixel_tof/cut_optimization.C):
// a fit's significance is no fair objective (its smooth background can pass under the peak), the counting one is; the
// simulation's own noise moves Z_A by about +-1 sigma, so a cut is only kept as far as both halves of the events agree.
// B's smooth part, the rest (entries without any muonic X-ray, from muons below and above 1 GeV), is not counted in the
// window but estimated from the whole fit range, 30-160 keV, times the window's share of the rest there before the cuts
// (the TOF trigger alone; smoothRest, the default since 2026-11-05): the fast band has few simulated entries, each of
// its days' weight, and counted in the window alone the search took cuts that let in more of them outside it (the
// first LYSO profile: the rest at 30-160 keV +70%, in the window +9%). The target X-rays' other deposits and the other
// volumes' muonic X-rays are peaked and counted in the window.
// The variables of a profile (analysis/CutProfile.h): track, zenith, topE, midE, tkin, otherE and npix on the grids
// below; entryE kept as in the start profile (a cut above the window does not change the peak's count).
// Search: coordinate descent (each variable in turn, its best window with the others fixed, until nothing changes) on
// every event from two starts, the TOF trigger alone and the start profile; then from the better result on the even
// events alone and on the odd ones alone (events by number), each scored on the other half. The profile written is the
// consensus of the two halves: for each variable the looser of their two choices (no cut where either has none).
// Each step keeps at least 30 simulated signal entries in the window.
// Writes in <out> (default <dir>): <name>.cuts, the new profile (copy it to analysis/cuts/ to use it by name), and
// <name>_optimisation.txt, the report (also printed): the starts, the optima, the cross-validation, the consensus and
// each variable's scan around it.
//   root -l -b -q 'analysis/cut_profile_optimize.C("output/lyso/analysis/si_prelim", "Si", "lyso", "czt")'
// smoothRest false: the rest counted in the window too (the CZT optimisation's way).

#include <TFile.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "CutProfile.h"

namespace {
const int kNV = 8;   // CutProfile::Variables(): track, zenith, topE, midE, tkin, entryE, otherE, npix
const double kOpen = 1e9;
struct Window { double lo = -kOpen, hi = kOpen; };
using Set = std::array<Window, kNV>;
// kind: 0 target X-ray, 1 other muonic X-ray, 2 the rest from the slow band, 3 from the fast band; window: in the L_alpha window
struct Entry { double w; int half, kind; bool full, window; double v[kNV]; };
struct Score { double zA = 0., s = 0., b = 0., bTarget = 0., bOther = 0., bSlow = 0., bFast = 0.; long sRaw = 0, bRaw = 0; };

double ZA(double s, double b) { return b > 0. ? std::sqrt(2. * ((s + b) * std::log(1. + s / b) - s)) : 0.; }
bool Passes(const Entry& x, const Set& c)
{
  for (int i = 0; i < kNV; ++i) if (x.v[i] < c[i].lo || x.v[i] > c[i].hi) return false;
  return true;
}
CutProfile ToProfile(const Set& c)
{
  CutProfile p;
  const auto& names = CutProfile::Variables();
  for (int i = 0; i < kNV; ++i)
    if (c[i].lo > -kOpen / 10. || c[i].hi < kOpen / 10.) p.range[names[i]] = {c[i].lo, c[i].hi};
  return p;
}
Set FromProfile(const CutProfile& p)
{
  Set c;
  const auto& names = CutProfile::Variables();
  for (int i = 0; i < kNV; ++i) if (p.Has(names[i])) c[i] = {p.Lo(names[i]), p.Hi(names[i])};
  return c;
}
std::string Describe(const Set& c) { const std::string d = ToProfile(c).Describe(); return d.empty() ? "no cuts" : d; }
}  // namespace

void cut_profile_optimize(const char* dir = "output/lyso/analysis/si_prelim", const char* target = "Si", const char* name = "lyso",
                          const char* start = "czt", const char* outDir = "",
                          const char* title = "the LYSO panel under the 4 cm target (12.69 x 12.69 cm: 20 x 20 pixels of 6.25 mm, each read out alone), with the measured TOF",
                          bool smoothRest = true)
{
  const double kFitLow = 30., kFitHigh = 160.;   // the fits' range (keV), as muonic_xray_target_select.C
  const std::string out = std::string(outDir).empty() ? dir : outDir;
  const double kAlpha = std::string(target) == "Si" ? 76.24 : 65.73;   // the target's L_alpha as simulated (keV)
  CutProfile startProfile;
  if (!startProfile.Load(start)) return;
  // ---- the window entries, from the skims
  double sigma1MeV = 0., days[2] = {};
  std::vector<Entry> entries;
  for (int band = 0; band < 2; ++band) {
    TFile file(Form("%s/xray_skim_%s.root", dir, band ? "fast" : "slow"));
    auto t = dynamic_cast<TTree*>(file.Get("cand"));
    auto info = dynamic_cast<TTree*>(file.Get("info"));
    if (!t || !info) { printf("no skim %s/xray_skim_%s.root\n", dir, band ? "fast" : "slow"); return; }
    if (band == 0) {
      auto p = dynamic_cast<TParameter<double>*>(file.Get("sigma1MeV"));
      if (!p) { printf("the skim %s/xray_skim_slow.root has no sigma1MeV\n", dir); return; }
      sigma1MeV = p->GetVal();
    }
    double lt = 0., live = 0.;   // summed: a merged skim has an info row per part
    info->SetBranchAddress("live_time", &lt);
    for (Long64_t i = 0; i < info->GetEntries(); ++i) { info->GetEntry(i); live += lt; }
    days[band] = live / 86400.;
    const double window = 2. * sigma1MeV * std::sqrt(kAlpha / 1000.);
    const double w = band ? days[0] / days[1] : 1.;   // entries of the slow band's days
    int event, src, trig, npix;
    float e, etrue, t0, tMidM, tkin, xTgt, yTgt, eOther, ekin, eMidM, eTopM, slopeX, slopeY;
    t->SetBranchAddress("event", &event); t->SetBranchAddress("src", &src); t->SetBranchAddress("trig", &trig);
    t->SetBranchAddress("npix", &npix); t->SetBranchAddress("e", &e); t->SetBranchAddress("etrue", &etrue);
    t->SetBranchAddress("t", &t0); t->SetBranchAddress("tMidM", &tMidM); t->SetBranchAddress("tkin", &tkin);
    t->SetBranchAddress("xTgt", &xTgt); t->SetBranchAddress("yTgt", &yTgt); t->SetBranchAddress("eOther", &eOther);
    t->SetBranchAddress("ekin", &ekin); t->SetBranchAddress("eMidM", &eMidM); t->SetBranchAddress("eTopM", &eTopM);
    t->SetBranchAddress("slopeX", &slopeX); t->SetBranchAddress("slopeY", &slopeY);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!trig || tMidM < 0. || t0 - tMidM < -5. || t0 - tMidM > 30.) continue;
      if (e < kFitLow || e >= kFitHigh) continue;
      const bool signal = src & 1;
      if (band == 1 && (signal || ekin <= 1000.)) continue;   // the fast band: its background above 1 GeV only
      const CutVars v = CutVars::FromSkim(xTgt, yTgt, slopeX, slopeY, eTopM, eMidM, tkin, e, eOther, npix);
      Entry x;
      x.w = w; x.half = event & 1; x.full = signal && std::abs(etrue - kAlpha) < 0.3; x.window = std::abs(e - kAlpha) <= window;
      x.kind = signal ? 0 : (src & 2) ? 1 : 2 + band;   // where B comes from, for the report
      const auto& names = CutProfile::Variables();
      for (int k = 0; k < kNV; ++k) x.v[k] = v.Get(names[k]);
      entries.push_back(x);
    }
  }
  const double window = 2. * sigma1MeV * std::sqrt(kAlpha / 1000.);
  // the window's share of the rest at 30-160 keV, after the trigger alone (and the start profile's entry-energy cut)
  double restAll = 0., restWindow = 0.;
  for (const auto& x : entries) {
    if (x.kind < 2) continue;
    restAll += x.w;
    if (x.window) restWindow += x.w;
  }
  const double share = restAll > 0. ? restWindow / restAll : 0.;
  auto score = [&](const Set& c, int half) {   // half -1: every event; 0, 1: the even or odd ones, weights doubled
    Score r;
    for (const auto& x : entries) {
      if (half >= 0 && x.half != half) continue;
      if (!Passes(x, c)) continue;
      const double w = (half < 0 ? 1. : 2.) * x.w;
      if (x.kind >= 2 && smoothRest) {   // the rest: its whole fit range, times the window's share
        (x.kind == 2 ? r.bSlow : r.bFast) += share * w;
        if (x.window) ++r.bRaw;
        continue;
      }
      if (!x.window) continue;
      if (x.full) { r.s += w; ++r.sRaw; continue; }
      ++r.bRaw;
      (x.kind == 0 ? r.bTarget : x.kind == 1 ? r.bOther : x.kind == 2 ? r.bSlow : r.bFast) += w;
    }
    r.b = r.bTarget + r.bOther + r.bSlow + r.bFast;
    r.zA = ZA(r.s, r.b);
    return r;
  };

  // ---- the grids (the CZT optimisation's, wider for the track and the other pixels)
  std::vector<std::vector<double>> loGrid(kNV, {-kOpen}), hiGrid(kNV, {kOpen});
  hiGrid[0] = {4., 5., 6., 7., 8., 9., 10., 11., 12., 14., 16., 20., kOpen};
  hiGrid[1] = {15., 20., 25., 30., 35., 40., 45., 50., 60., kOpen};
  loGrid[2] = {-kOpen, 1., 1.2, 1.4, 1.5, 1.6, 1.7, 1.8, 2., 2.2, 2.5}; hiGrid[2] = {3., 3.5, 4., 5., 6., 8., kOpen};
  loGrid[3] = {-kOpen, 1., 1.2, 1.4, 1.5, 1.6, 1.8, 2., 2.2, 2.5, 3.}; hiGrid[3] = {3., 4., 5., 6., 8., 10., kOpen};
  loGrid[4] = {-kOpen, 5., 10., 15., 20., 25., 30.}; hiGrid[4] = {40., 50., 60., 70., 80., 90., 100., 120., 150., 200., kOpen};
  hiGrid[6] = {0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1., 1.5, 2., 3., kOpen};
  hiGrid[7] = {1., 2., 3., 4., 5., 6., 8., 10., kOpen};
  const Set startSet = FromProfile(startProfile);
  loGrid[5] = {startSet[5].lo}; hiGrid[5] = {startSet[5].hi};   // the entry's energy: as the start profile
  const int kMinSignal = 30;   // simulated signal entries in the window, at least (on a half: 15)

  auto descend = [&](Set c, int half, int& evaluations) {
    double best = score(c, half).zA;
    for (int round = 0; round < 10; ++round) {
      bool changed = false;
      for (int i = 0; i < kNV; ++i) {
        Set trial = c, bestSet = c;
        for (double lo : loGrid[i])
          for (double hi : hiGrid[i]) {
            if (lo >= hi) continue;
            trial[i] = {lo, hi};
            ++evaluations;
            const Score r = score(trial, half);
            if (r.sRaw < (half < 0 ? kMinSignal : kMinSignal / 2)) continue;
            if (r.zA > best + 1e-3) { best = r.zA; bestSet = trial; changed = true; }
          }
        c = bestSet;
      }
      if (!changed) break;
    }
    return c;
  };

  std::ostringstream o;
  char when[32];
  const std::time_t now = std::time(nullptr);
  std::strftime(when, sizeof(when), "%Y-%m-%d %H:%M", std::localtime(&now));
  o << "Cut optimisation for the profile " << name << " (analysis/cut_profile_optimize.C, " << when << ")\n";
  o << Form("skims %s: slow band %.2f days, fast band %.2f days (scaled x %.3f, its background above 1 GeV only)\n", dir, days[0],
            days[1], days[0] / days[1]);
  o << Form("%s target, L_alpha %.2f keV, window +-%.2f keV (2 sigma of %.4g keV x sqrt(E / MeV), the skims' resolution); %zu "
            "entries at %.0f-%.0f keV after the trigger\n", target, kAlpha, window, sigma1MeV, entries.size(), kFitLow, kFitHigh);
  o << Form("objective: Z_A of S (full-energy L_alpha entries in the window) against B (every other entry there), in %.2f days\n",
            days[0]);
  if (smoothRest)
    o << Form("   B's rest (no muonic X-ray) from the whole %.0f-%.0f keV x %.4f, the window's share of it after the trigger alone\n\n",
              kFitLow, kFitHigh, share);
  else
    o << "   B counted in the window alone\n\n";
  auto row = [&](const char* label, const Set& c) {
    const Score a = score(c, -1), e = score(c, 0), d = score(c, 1);
    o << Form("   %-28s Z_A %5.2f (even %5.2f, odd %5.2f); S %6.1f (%ld), B %6.1f (%ld) = target %.1f | other muonic %.1f | "
              "slow %.1f | fast %.1f\n", label, a.zA, e.zA, d.zA, a.s, a.sRaw, a.b, a.bRaw, a.bTarget, a.bOther, a.bSlow, a.bFast);
    o << "      " << Describe(c) << "\n";
    return a;
  };
  o << "the starts (every event; in brackets the simulated entries; B by origin: the target's X-rays' other deposits, other\n"
       "volumes' muonic X-rays, the rest from muons below and above 1 GeV):\n";
  Set trigger;
  trigger[5] = startSet[5];
  row("TOF trigger only", trigger);
  row(Form("start profile %s", startProfile.name.c_str()), startSet);

  // ---- the search on every event, from both starts
  int evaluations = 0;
  const Set fromTrigger = descend(trigger, -1, evaluations), fromStart = descend(startSet, -1, evaluations);
  const Set best = score(fromTrigger, -1).zA > score(fromStart, -1).zA ? fromTrigger : fromStart;
  o << "\nthe optimum on every event (coordinate descent):\n";
  row("from the trigger only", fromTrigger);
  row(Form("from %s", startProfile.name.c_str()), fromStart);
  // ---- cross-validation: chosen on one half, scored on the other
  const Set onEven = descend(best, 0, evaluations), onOdd = descend(best, 1, evaluations);
  o << "\ncross-validation (from the better optimum, chosen on one half of the events alone):\n";
  const Score ee = score(onEven, 0), eo = score(onEven, 1), oo = score(onOdd, 1), oe = score(onOdd, 0);
  row("chosen on the even events", onEven);
  row("chosen on the odd events", onOdd);
  o << Form("   the even's choice: %.2f on the even events, %.2f on the odd; the odd's choice: %.2f on the odd, %.2f on the even\n",
            ee.zA, eo.zA, oo.zA, oe.zA);
  o << Form("   (the gain a half finds on itself but not on the other is the simulation's noise; %d evaluations)\n", evaluations);
  // ---- the consensus: for each variable the looser of the two halves' choices
  Set consensus;
  for (int i = 0; i < kNV; ++i) consensus[i] = {std::min(onEven[i].lo, onOdd[i].lo), std::max(onEven[i].hi, onOdd[i].hi)};
  o << "\nthe consensus (each variable the looser of the halves' choices), written as the profile:\n";
  const Score c = row(name, consensus);
  const Score s0 = score(startSet, -1);
  o << Form("   against the start profile %s: S %+.1f%%, B %+.1f%%, Z_A %+.2f\n", startProfile.name.c_str(),
            s0.s > 0. ? 100. * (c.s / s0.s - 1.) : 0., s0.b > 0. ? 100. * (c.b / s0.b - 1.) : 0., c.zA - s0.zA);
  // ---- each variable's scan around the consensus: its upper (or lower) side over its grid, the others kept
  o << "\neach variable's scan around the consensus (Z_A every event | even | odd; S; B):\n";
  const auto& names = CutProfile::Variables();
  for (int i = 0; i < kNV; ++i) {
    if (i == 5) continue;   // the entry's energy: not scanned
    for (int side = 0; side < 2; ++side) {
      const auto& grid = side ? hiGrid[i] : loGrid[i];
      if (grid.size() < 2) continue;
      o << Form("   %s, %s side:\n", names[i].c_str(), side ? "upper" : "lower");
      for (double v : grid) {
        Set trial = consensus;
        (side ? trial[i].hi : trial[i].lo) = v;
        if (trial[i].lo >= trial[i].hi) continue;
        const Score a = score(trial, -1), e = score(trial, 0), d = score(trial, 1);
        const bool chosen = v == (side ? consensus[i].hi : consensus[i].lo);
        o << Form("      %8s  %5.2f | %5.2f | %5.2f   S %6.1f   B %6.1f%s\n", std::abs(v) >= kOpen ? "open" : Form("%g", v), a.zA,
                  e.zA, d.zA, a.s, a.b, chosen ? "   <- the profile" : "");
      }
    }
  }
  printf("%s", o.str().c_str());
  gSystem->mkdir(out.c_str(), true);
  std::ofstream(Form("%s/%s_optimisation.txt", out.c_str(), name)) << o.str();

  // ---- the profile
  CutProfile profile = ToProfile(consensus);
  const std::string restNote = smoothRest ? Form(" (the rest from %.0f-%.0f keV x %.3f)", kFitLow, kFitHigh, share) : "";
  profile.notes = {
      Form("Cut profile %s: %s.", name, title),
      Form("Optimised on %s by analysis/cut_profile_optimize.C (report %s_optimisation.txt) on the skims of", when, name),
      Form("%s (%s target, L_alpha %.2f keV, window +-%.2f keV; slow band %.2f days, fast band %.2f days):", dir, target, kAlpha,
           window, days[0], days[1]),
      Form("the counting significance Z_A of the full-energy L_alpha entries in the window%s, %.2f in %.2f days (even events %.2f,",
           restNote.c_str(), c.zA, days[0], score(consensus, 0).zA),
      Form("odd %.2f; the start profile %s %.2f, the trigger alone %.2f); each variable the looser of the halves' choices.",
           score(consensus, 1).zA, startProfile.name.c_str(), s0.zA, score(trigger, -1).zA),
      "A cut keeps lo <= value <= hi; \"-\" is an open side (analysis/CutProfile.h).", "variable     lo       hi"};
  const std::string file = Form("%s/%s.cuts", out.c_str(), name);
  profile.Write(file);
  printf("\nwritten: %s and %s/%s_optimisation.txt\n", file.c_str(), out.c_str(), name);
}
