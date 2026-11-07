// The pixel panel's own radioactivity against the TOF trigger and the cuts (2026-11-03, for the Lu-176 of the LYSO panel,
// macro/muonic_xray_lyso.mac): how much of it each step of the CZT analysis lets through, next to the target's X-rays.
// Inputs:
//   <dir>/xray_skim_slow.root, xray_skim_fast.root: skims of a panel run with /tof/intrinsic/ (the decays within +-1 us of
//     every cosmic-ray event), analysis/muonic_xray_stack_skim.C with detector "czt" and the panel's sigma1MeV and pitch;
//   decayFiles: the decay-only runs (macro/lyso_decay_only.mac: one decay per event, live time = events / activity);
//   runDir: the panel run's folder, for the rates its job logs print (cosmic rays through the panel).
// The steps, for the panel's decays (truth: entries whose energy is all theirs, czt_src == 16) and for the target's
// X-rays (src bit 1):
//   1. the decays in the slow band's days (activity x time), 2. their single-pixel entries of 20-1000 keV with the
//   response (sigma = sigma1MeV x sqrt(E / MeV), 10 keV threshold; the decay-only runs), 3. those within -5..+30 ns of
//   the middle TOF of a triggered event (top and middle station within 20 ns, measured), 4. of these after the cuts of
//   the profile cutProfile (analysis/cuts/<name>.cuts, analysis/CutProfile.h; czt: the eight cuts of the CZT analysis,
//   event_flow.C), 5. of these in the target's L-alpha window (+-2 sigma of the panel).
// Steps 3-5 three ways: counted in the window; from the sidebands -990..-10 and +40..+990 ns (the same decays, 55 times
// as many, scaled by 35 / 1930; the one before the muon is what a measurement can use too, without truth); and, for
// step 3, expected as triggers x singles rate x 35 ns (the trigger rate of acceptance.mac, triggerPerMinute). Counts are those of the slow band's
// simulated days (the fast band scaled to them; as event_flow.C, the fast band only for its background above 1 GeV);
// the text also projects them to projectDays. Writes in <out> (default <dir>), one plot per file (PNG and PDF):
//   panel_intrinsic.txt              the numbers (also printed)
//   panel_intrinsic_singles          the decays' single-pixel spectrum, all pixels, per second and 2 keV (no trigger)
//   panel_intrinsic_time             the entries of triggered events against the middle-TOF time, by origin
//   panel_intrinsic_before           the spectrum within -5..+30 ns of triggered events, by origin, before the cuts
//   panel_intrinsic_after            the same after the cuts
//   panel_intrinsic_flow             (2026-11-05, to present the suppression) the steps as bars, from the decays to the
//                                    L-alpha window, each with the share it keeps, beside the target's X-rays
//   panel_intrinsic_suppression      (2026-11-05) the decays' spectrum with no trigger, in time with one and after the cuts
//                                    (the sidebands), beside the target's X-rays after the cuts; per 2 keV in 10 keV bins
//   root -l -b -q 'analysis/panel_intrinsic_check.C("output/lyso_concept", "output/lyso_concept/decay/decay_*.root", "output/lyso_concept/cosmic_al")'

#include <TArrow.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TBox.h>
#include <TParameter.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "CutProfile.h"

namespace {
// a count of entries: raw (simulated) and weighted to the slow band's days, with its error
struct Count {
  double raw = 0., w = 0., w2 = 0.;
  void Add(double weight) { raw += 1.; w += weight; w2 += weight * weight; }
  double Err() const { return std::sqrt(w2); }
};
std::string Show(const Count& c, double scale = 1.)
{
  return Form("%10.1f +- %6.1f (%6.0f simulated)", scale * c.w, scale * c.Err(), c.raw);
}
// the rates the panel run's job logs print: "  cosmic rays: 0.18 /cm2/min, through the panel's 400 cm2 1.2 /s: ..."
double PanelMuonRate(const char* runDir, const char* band)
{
  std::ifstream log(Form("%s/%s/log/job1.log", runDir, band));
  std::string line;
  while (std::getline(log, line)) {
    const auto at = line.find("cosmic rays: ");
    const auto through = line.find("through the panel's ");
    if (at == std::string::npos || through == std::string::npos) continue;
    std::istringstream words(line.substr(through + 20));
    double area = 0., rate = 0.;
    std::string unit;
    words >> area >> unit >> rate;
    return rate;
  }
  return -1.;
}
// the panel's outer width (cm) as the run's log prints it ("LYSO plane: ... frame (126.88 mm wide), ..."), 0 if not
double PanelWidth(const char* runDir)
{
  std::ifstream log(Form("%s/slow/log/job1.log", runDir));
  std::string line;
  for (int i = 0; i < 3000 && std::getline(log, line); ++i) {
    if (line.find(" plane: ") == std::string::npos) continue;
    const auto wide = line.find(" mm wide)"), open = line.rfind('(', line.find(" mm wide)"));
    if (wide != std::string::npos && open != std::string::npos) return std::stod(line.substr(open + 1, wide - open - 1)) / 10.;
    break;
  }
  return 0.;
}
void Save(TCanvas* c, const std::string& base)
{
  c->Print((base + ".png").c_str());
  c->Print((base + ".pdf").c_str());
}
}  // namespace

void panel_intrinsic_check(const char* dir = "output/lyso_concept", const char* decayFiles = "output/lyso_concept/decay/decay_*.root",
                           const char* runDir = "output/lyso_concept/cosmic_al", const char* target = "Al", double sigma1MeV = 26.85,
                           double triggerPerMinute = 493., double projectDays = 60., const char* outDir = "",
                           const char* cutProfile = "czt")
{
  CutProfile profile;
  if (!profile.Load(cutProfile)) return;
  const std::string out = std::string(outDir).empty() ? dir : outDir;
  const double line = std::string(target) == "Si" ? 76.24 : 65.73;   // the target's L-alpha as simulated (keV)
  const double halfWindow = 2. * sigma1MeV * std::sqrt(line / 1000.);   // +-2 sigma of the panel
  const double drawnHigh = 500.;   // both targets, as muonic_xray_target_select.C since 2026-11-05 (Al's ended at 450 before)
  const double tLo = -5., tHi = 30.;   // the window, ns from the middle TOF
  // the sidebands: before the muon, and after its prompt hits (there the truth tells the decays from its delayed ones)
  const double s1Lo = -990., s1Hi = -10., s2Lo = 40., s2Hi = 990.;
  const double sideScale = (tHi - tLo) / ((s1Hi - s1Lo) + (s2Hi - s2Lo));
  gStyle->SetOptStat(0);
  TRandom3 random(17601);

  // ---- 1-2. the decays alone: activity and single-pixel spectrum
  TChain decays("events"), decayRuns("run");
  decays.Add(decayFiles);
  decayRuns.Add(decayFiles);
  double decayLive = 0., lt = 0.;
  long long decayEvents = 0;
  int ne = 0;
  decayRuns.SetBranchAddress("live_time", &lt);
  decayRuns.SetBranchAddress("n_events", &ne);
  for (Long64_t i = 0; i < decayRuns.GetEntries(); ++i) { decayRuns.GetEntry(i); decayLive += lt; decayEvents += ne; }
  if (decayLive <= 0.) { printf("no decay-only runs %s\n", decayFiles); return; }
  const double activity = decayEvents / decayLive;
  auto singles = new TH1D("singles", "", 650, 0., 1300.);
  std::vector<double>* he = nullptr;
  decays.SetBranchAddress("czt_e", &he);
  double r10 = 0., rEntries = 0., rLine = 0.;   // entries per second: above 10 keV, 20-1000 keV, in the L-alpha window
  for (Long64_t i = 0; i < decays.GetEntries(); ++i) {
    decays.GetEntry(i);
    for (double e : *he) {
      const double keV = 1000. * (e + random.Gaus(0., sigma1MeV / 1000. * std::sqrt(e)));
      if (keV < 10.) continue;
      r10 += 1. / decayLive;
      singles->Fill(keV, 1. / decayLive);
      if (keV >= 20. && keV <= 1000.) rEntries += 1. / decayLive;
      if (std::abs(keV - line) <= halfWindow) rLine += 1. / decayLive;
    }
  }

  // ---- 3-5. the panel run: triggered events, by origin
  enum { kTarget, kIntrinsic, kPileUp, kOtherMu, kRest, kClasses };
  const char* className[kClasses] = {"muonic X-ray of the target", "panel's own decays (Lu-176) alone",
                                     "Lu-176 piled up on a muon's hit", "other muonic X-ray", "the rest (muon-made)"};
  const int colour[kClasses] = {TColor::GetColor("#d62728"), TColor::GetColor("#9467bd"), TColor::GetColor("#c5b0d5"),
                                TColor::GetColor("#ff7f0e"), TColor::GetColor("#1f77b4")};
  TH1D* before[kClasses];
  TH1D* after[kClasses];
  const int nBins = std::lround((drawnHigh - 20.) / 2.);
  for (int k = 0; k < kClasses; ++k) {
    before[k] = new TH1D(Form("before%d", k), "", nBins, 20., drawnHigh);
    after[k] = new TH1D(Form("after%d", k), "", nBins, 20., drawnHigh);
  }
  auto timeLu = new TH1D("timeLu", "", 200, -1000., 1000.);
  auto timeMu = new TH1D("timeMu", "", 200, -1000., 1000.);
  auto sideAfterSpectrum = new TH1D("sideAfterSpectrum", "", nBins, 20., drawnHigh);   // the sidebands' decays after the cuts
  Count inWindow[kClasses], afterCuts[kClasses], lineWindow[kClasses], sideBefore, sideAfter, sideLine, peak;
  double slowDays = 0., bandDays[2] = {0., 0.};
  long pixels = 0;   // the panel's pixels (the slow skim's pixelsPerSide squared)
  long long bandMuons[2] = {0, 0};
  const char* bands[2] = {"slow", "fast"};
  for (int b = 0; b < 2; ++b) {
    TFile file(Form("%s/xray_skim_%s.root", dir, bands[b]));
    auto t = dynamic_cast<TTree*>(file.Get("cand"));
    auto info = dynamic_cast<TTree*>(file.Get("info"));
    if (!t || !info) { printf("no skim %s/xray_skim_%s.root\n", dir, bands[b]); return; }
    double live = 0., l = 0.;
    long long mu = 0, m = 0;
    info->SetBranchAddress("live_time", &l);
    info->SetBranchAddress("muons", &m);
    for (Long64_t i = 0; i < info->GetEntries(); ++i) { info->GetEntry(i); live += l; mu += m; }
    bandDays[b] = live / 86400.;
    bandMuons[b] = mu;
    if (b == 0) slowDays = bandDays[0];
    if (auto n = dynamic_cast<TParameter<double>*>(file.Get("pixelsPerSide")); n && b == 0) pixels = std::lround(n->GetVal() * n->GetVal());
    const double w = slowDays / bandDays[b];
    int src, trig, npix;
    float e, etrue, t0, tMidM, tkin, xTgt, yTgt, eOther, ekin, eMidM, eTopM, slopeX, slopeY;
    t->SetBranchAddress("src", &src); t->SetBranchAddress("trig", &trig); t->SetBranchAddress("e", &e);
    t->SetBranchAddress("etrue", &etrue); t->SetBranchAddress("t", &t0); t->SetBranchAddress("tMidM", &tMidM);
    t->SetBranchAddress("tkin", &tkin); t->SetBranchAddress("xTgt", &xTgt); t->SetBranchAddress("yTgt", &yTgt);
    t->SetBranchAddress("eOther", &eOther); t->SetBranchAddress("ekin", &ekin); t->SetBranchAddress("eMidM", &eMidM);
    t->SetBranchAddress("eTopM", &eTopM); t->SetBranchAddress("slopeX", &slopeX); t->SetBranchAddress("slopeY", &slopeY);
    t->SetBranchAddress("npix", &npix);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!trig || tMidM < 0.) continue;
      if (b == 1 && ((src & 1) || ekin <= 1000.)) continue;   // the fast band: its background above 1 GeV only
      const int kind = (src & 1) ? kTarget : src == 16 ? kIntrinsic : (src & 16) ? kPileUp : (src & 2) ? kOtherMu : kRest;
      const double dt = t0 - tMidM;
      (kind == kIntrinsic ? timeLu : timeMu)->Fill(dt, w);
      // the cuts of the profile, each kept when lo <= v <= hi
      const bool cuts = profile.Pass(xTgt, yTgt, slopeX, slopeY, eTopM, eMidM, tkin, e, eOther, npix);
      const bool inLine = std::abs(e - line) <= halfWindow;
      if (kind == kIntrinsic && ((dt >= s1Lo && dt <= s1Hi) || (dt >= s2Lo && dt <= s2Hi))) {
        sideBefore.Add(w);
        if (cuts) { sideAfter.Add(w); sideAfterSpectrum->Fill(e, w); }
        if (cuts && inLine) sideLine.Add(w);
      }
      if (dt < tLo || dt > tHi) continue;
      inWindow[kind].Add(w);
      before[kind]->Fill(e, w);
      if (!cuts) continue;
      afterCuts[kind].Add(w);
      after[kind]->Fill(e, w);
      if (!inLine) continue;
      lineWindow[kind].Add(w);
      if (kind == kTarget && std::abs(etrue - line) < 0.3) peak.Add(w);
    }
  }

  // ---- the numbers
  const double T = slowDays * 86400.;   // s
  const double project = projectDays / slowDays;
  const double triggers = triggerPerMinute / 60. * T;
  const double expected = triggers * rEntries * (tHi - tLo) * 1e-9;
  const double muSlow = PanelMuonRate(runDir, "slow"), muFast = PanelMuonRate(runDir, "fast");
  std::ostringstream o;
  o << "The pixel panel's own radioactivity against the TOF trigger (" << target << " target, L-alpha " << line
    << " keV, window +-" << Form("%.2f", halfWindow) << " keV = 2 sigma of sigma1MeV " << sigma1MeV << " keV)\n";
  o << Form("panel run: slow band %.3f days (%lld muons), fast band %.3f days (%lld muons, scaled x %.3f to the slow band's days)\n",
            bandDays[0], bandMuons[0], bandDays[1], bandMuons[1], slowDays / bandDays[1]);
  o << Form("decay-only runs: %lld decays = %.3f s\n", decayEvents, decayLive);
  o << Form("cuts: profile %s (%s): %s\n\n", profile.name.c_str(), profile.path.c_str(), profile.Describe().c_str());
  o << "rates\n";
  if (pixels > 0)
    o << Form("   decays of the panel (activity)                      %10.0f /s   (%.1f /s per pixel of %ld)\n", activity,
              activity / pixels, pixels);
  else
    o << Form("   decays of the panel (activity)                      %10.0f /s\n", activity);
  o << Form("   single-pixel entries above 10 keV                   %10.0f /s\n", r10);
  o << Form("   single-pixel entries of 20-1000 keV                 %10.0f /s\n", rEntries);
  o << Form("   ... of these in the L-alpha window                  %10.1f /s\n", rLine);
  if (muSlow > 0. && muFast > 0.)
    o << Form("   cosmic muons through the panel (job logs)            %10.2f /s   (0-1 GeV %.2f, 1-1000 GeV %.2f): %.0f decays per muon\n",
              muSlow + muFast, muSlow, muFast, activity / (muSlow + muFast));
  o << Form("   TOF triggers (acceptance.mac)                        %10.2f /s   (%.0f decays per trigger)\n", triggerPerMinute / 60.,
            activity / (triggerPerMinute / 60.));
  o << Form("   chance of an entry of 20-1000 keV within the %.0f ns window of a trigger: %.2e\n\n", tHi - tLo,
            rEntries * (tHi - tLo) * 1e-9);
  o << Form("the steps, in the slow band's %.3f days [and projected to %.0f days]\n", slowDays, projectDays);
  o << Form("   1. decays of the panel                               %14.4g   [%.4g]\n", activity * T, activity * T * project);
  o << Form("   2. their entries of 20-1000 keV (singles)            %14.4g   [%.4g]\n", rEntries * T, rEntries * T * project);
  o << Form("      ... in the L-alpha window                         %14.4g   [%.4g]\n", rLine * T, rLine * T * project);
  o << "   3. in the window of a trigger (-5..+30 ns)\n";
  o << "      counted:   " << Show(inWindow[kIntrinsic]) << "   [" << Show(inWindow[kIntrinsic], project) << "]\n";
  o << "      sideband:  " << Show(sideBefore, sideScale) << "   [" << Show(sideBefore, sideScale * project) << "]\n";
  o << Form("      expected:  %10.1f (triggers %.4g x %.0f /s x %.0f ns)   [%.1f]\n", expected, triggers, rEntries, tHi - tLo,
            expected * project);
  o << "   4. after the cuts\n";
  o << "      counted:   " << Show(afterCuts[kIntrinsic]) << "   [" << Show(afterCuts[kIntrinsic], project) << "]\n";
  o << "      sideband:  " << Show(sideAfter, sideScale) << "   [" << Show(sideAfter, sideScale * project) << "]\n";
  o << "   5. of these in the L-alpha window\n";
  o << "      counted:   " << Show(lineWindow[kIntrinsic]) << "   [" << Show(lineWindow[kIntrinsic], project) << "]\n";
  o << "      sideband:  " << Show(sideLine, sideScale) << "   [" << Show(sideLine, sideScale * project) << "]\n\n";
  o << "every origin within -5..+30 ns of a trigger: before the cuts | after them | of these in the L-alpha window [projected]\n";
  for (int k = 0; k < kClasses; ++k)
    o << Form("   %-36s %10.1f | %8.1f | %7.1f   [%9.0f | %7.0f | %6.0f]\n", className[k], inWindow[k].w, afterCuts[k].w,
              lineWindow[k].w, inWindow[k].w * project, afterCuts[k].w * project, lineWindow[k].w * project);
  o << Form("   full-energy L-alpha photons of the target after the cuts: %.1f +- %.1f   [%.0f]\n", peak.w, peak.Err(),
            peak.w * project);
  o << Form("(sideband: the panel's decays alone at -990..-10 and +40..+990 ns of the middle TOF, times %.4f; its cuts on the\n"
            " other pixels see the hits around it, not the muon's prompt ones, so it may overestimate a little after the cuts)\n",
            sideScale);
  printf("%s", o.str().c_str());
  std::ofstream(out + "/panel_intrinsic.txt") << o.str();

  // ---- the plots, one per file
  const std::string days = Form("%.2f days of cosmic rays (the slow band's; the fast band scaled to them)", slowDays);
  // 1600 x 1100 since 2026-11-05 (asked for: the spectra got crowded), the texts as large as on the 1000 x 700 before
  const double ts = 0.64;
  auto bigCanvas = [](const char* name) {
    auto c = new TCanvas(name, "", 1600, 1100);
    c->SetLogy(); c->SetLeftMargin(0.08); c->SetRightMargin(0.03); c->SetTopMargin(0.05); c->SetBottomMargin(0.09);
    return c;
  };
  auto axes = [&](TAxis* x, TAxis* y) {
    for (TAxis* a : {x, y}) { a->SetTitleSize(0.04 * ts); a->SetLabelSize(0.035 * ts); }
    x->SetTitleOffset(1.0); y->SetTitleOffset(0.9);
  };
  TLatex note;
  note.SetNDC();
  note.SetTextFont(42);
  note.SetTextSize(0.03 * ts);
  {
    auto c = bigCanvas("c_singles");
    singles->SetTitle(pixels > 0 ? Form(";single-pixel energy with the response (keV);entries per second per 2 keV, all %ld pixels", pixels)
                                 : ";single-pixel energy with the response (keV);entries per second per 2 keV, all pixels");
    singles->SetLineColor(colour[kIntrinsic]); singles->SetLineWidth(2);
    singles->GetXaxis()->SetRangeUser(10., 1300.);
    axes(singles->GetXaxis(), singles->GetYaxis());
    singles->Draw("hist");
    c->Update();
    auto band = new TBox(line - halfWindow, std::pow(10., c->GetUymin()), line + halfWindow, std::pow(10., c->GetUymax()));
    band->SetFillColorAlpha(colour[kTarget], 0.25);
    band->Draw();
    singles->Draw("hist same");
    note.DrawLatex(0.52, 0.90, Form("the panel's Lu-176 alone, no trigger: %.0f decays/s", activity));
    note.DrawLatex(0.52, 0.868, Form("%.0f entries/s at 20-1000 keV, %.0f /s in the %s L_{#alpha} window", rEntries, rLine, target));
    note.DrawLatex(0.52, 0.836, Form("#sigma = %.2f keV #times #sqrt{E/MeV}, single pixels, 10 keV threshold", sigma1MeV));
    Save(c, out + "/panel_intrinsic_singles");
  }
  {
    auto c = bigCanvas("c_time");
    timeMu->SetTitle(";entry time - middle-TOF time (ns);entries per 10 ns");
    axes(timeMu->GetXaxis(), timeMu->GetYaxis());
    timeMu->SetLineColor(colour[kRest]); timeMu->SetLineWidth(2);
    timeLu->SetLineColor(colour[kIntrinsic]); timeLu->SetLineWidth(2);
    timeMu->SetMaximum(std::max(timeMu->GetMaximum(), timeLu->GetMaximum()) * 10.);
    timeMu->SetMinimum(0.5);
    timeMu->Draw("hist");
    timeLu->Draw("hist same");
    c->Update();
    for (double x : {tLo, tHi}) {
      auto l = new TLine(x, std::pow(10., c->GetUymin()), x, std::pow(10., c->GetUymax()));
      l->SetLineStyle(2); l->Draw();
    }
    auto legend = new TLegend(0.10, 0.86, 0.60, 0.94);
    legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(0.03 * ts);
    legend->AddEntry(timeLu, "the panel's Lu-176 alone", "l");
    legend->AddEntry(timeMu, "made by the muon (incl. its X-rays; Lu-176 piled up on them)", "l");
    legend->Draw();
    note.DrawLatex(0.105, 0.83, Form("triggered events, entries of 20-1000 keV, before the cuts; dashed: the %.0f..+%.0f ns window", tLo, tHi));
    note.DrawLatex(0.105, 0.798, days.c_str());
    Save(c, out + "/panel_intrinsic_time");
  }
  for (int s = 0; s < 2; ++s) {
    TH1D** h = s ? after : before;
    auto c = bigCanvas(Form("c_spec%d", s));
    auto stack = new THStack(Form("stack%d", s), Form(";single-pixel energy with the response (keV);entries per 2 keV"));
    auto legend = new TLegend(0.58, 0.76, 0.97, 0.94);
    legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(0.028 * ts);
    const Count* n = s ? afterCuts : inWindow;
    for (int k : {kRest, kOtherMu, kPileUp, kIntrinsic, kTarget}) {
      h[k]->SetFillColor(colour[k]); h[k]->SetLineColor(colour[k]);
      stack->Add(h[k]);
    }
    for (int k : {kTarget, kIntrinsic, kPileUp, kOtherMu, kRest})
      legend->AddEntry(h[k], Form("%s: %.1f", className[k], n[k].w), "f");
    double total = 0.;
    for (int k = 0; k < kClasses; ++k) total += h[k]->Integral();
    if (total <= 0.) c->SetLogy(0);   // nothing to draw on a log scale (a small run after the cuts)
    stack->SetMinimum(total > 0. ? 0.5 : 0.);
    stack->Draw("hist");
    stack->SetMaximum(total > 0. ? stack->GetMaximum() * 30. : 1.);
    stack->Draw("hist");
    axes(stack->GetXaxis(), stack->GetYaxis());
    legend->Draw();
    note.DrawLatex(0.11, 0.90, s ? "-5..+30 ns of a trigger, after the cuts" : "-5..+30 ns of a trigger, before the cuts");
    const double width = PanelWidth(runDir);
    note.DrawLatex(0.11, 0.868, width > 0. ? Form("4 cm %s target, LYSO panel %.4g x %.4g cm", target, width, width)
                                           : Form("4 cm %s target, LYSO panel", target));
    note.DrawLatex(0.11, 0.836, Form("%.2f days", slowDays));
    Save(c, out + (s ? "/panel_intrinsic_after" : "/panel_intrinsic_before"));
  }

  // ---- two plots to present the suppression (2026-11-05, the user asked for one): the steps as bars, and their spectra
  const int violet = colour[kIntrinsic], red = colour[kTarget], darkViolet = TColor::GetColor("#54278f");
  auto grouped = [](double v) {   // 27664.5 -> 27,665
    const std::string digits = Form("%.0f", std::round(v));
    std::string g;
    for (std::size_t i = 0; i < digits.size(); ++i) {
      if (i > 0 && (digits.size() - i) % 3 == 0) g += ',';
      g += digits[i];
    }
    return g;
  };
  auto rounded = [&](double v) {   // 93.67 billion, 27,665, 13.0
    if (v >= 1e9) return std::string(Form("%.2f billion", v / 1e9));
    if (v >= 1e7) return std::string(Form("%.1f million", v / 1e6));
    return v >= 100. ? grouped(v) : std::string(Form("%.1f", v));
  };
  auto kept = [&](double f) {   // the share a step keeps: 67%, 2.3%, 1 in 914, 1 in 3.4 million
    if (f >= 0.01) return std::string(Form(f >= 0.1 ? "%.0f%%" : "%.1f%%", 100. * f));
    return "1 in " + (1. / f >= 1e6 ? std::string(Form("%.1f million", 1e-6 / f)) : grouped(1. / f));
  };
  const double panelWidth = PanelWidth(runDir);
  const std::string panelName = panelWidth > 0. ? Form("LYSO panel %.4g x %.4g cm", panelWidth, panelWidth) : "LYSO panel";
  // the steps: no trigger (1-2a, the decay-only runs), then in time with a trigger, after the cuts, in the L-alpha window
  const double free3[3] = {activity * T, rEntries * T, rLine * T};
  const double lu3[3] = {inWindow[kIntrinsic].w, afterCuts[kIntrinsic].w, lineWindow[kIntrinsic].w};
  const double side3[3] = {sideScale * sideBefore.w, sideScale * sideAfter.w, sideScale * sideLine.w};
  const double x3[3] = {inWindow[kTarget].w, afterCuts[kTarget].w, lineWindow[kTarget].w};
  {
    // horizontal bars as event_flow.C draws its flow: rows 1-2a one bar, rows 3-5 three (counted, sidebands, X-rays)
    auto numberText = [&](double v, double before) {   // the number, and the share of the step before (none if 0)
      return (v > 0. ? rounded(v) : std::string("0")) + (v > 0. && before > 0. ? "  (" + kept(v / before) + ")" : std::string());
    };
    struct Series { std::vector<double> values; std::vector<std::string> text; int colour, fill; double offset, width; std::string name; };
    const std::vector<Series> all = {
      {{free3[0], free3[1], free3[2], 0., 0., 0.}, {numberText(free3[0], 0.), numberText(free3[1], 0.), numberText(free3[2], free3[1]), "", "", ""},
       violet, 1001, 0.31, 0.38, "the panel's Lu-176 alone"},
      {{0., 0., 0., lu3[0], lu3[1], lu3[2]}, {"", "", "", numberText(lu3[0], free3[1]), numberText(lu3[1], lu3[0]), numberText(lu3[2], lu3[1])},
       violet, 1001, 0.64, 0.26, ""},
      {{0., 0., 0., side3[0], side3[1], side3[2]}, {"", "", "", numberText(side3[0], 0.), numberText(side3[1], side3[0]), numberText(side3[2], side3[1])},
       violet, 3354, 0.37, 0.26, "Lu-176 from the sidebands"},
      {{0., 0., 0., x3[0], x3[1], x3[2]}, {"", "", "", numberText(x3[0], 0.), numberText(x3[1], x3[0]), numberText(x3[2], x3[1])},
       red, 1001, 0.10, 0.26, className[kTarget]}};
    const std::vector<std::string> labels = {Form("1. Lu-176 decays in the panel (%s /s)", grouped(activity).c_str()),
                                             "2. their single-pixel entries of 20-1000 keV",
                                             Form("2a. of these: in the %s L-alpha window", target),
                                             Form("3. in time with a TOF trigger (%.0f..+%.0f ns)", tLo, tHi),
                                             "4. after the corner-plot cuts",
                                             Form("5. of these: in the %s L-alpha window", target)};
    const std::string cuts = profile.Describe();
    std::size_t split = std::string::npos;   // the cuts in two lines, at the comma nearest the middle
    for (std::size_t at = cuts.find(", "); at != std::string::npos; at = cuts.find(", ", at + 1))
      if (split == std::string::npos || std::abs(at - 0.5 * cuts.size()) < std::abs(split - 0.5 * cuts.size())) split = at;
    const std::vector<std::string> footer = {
      Form("4 cm %s target, %s; 1-2a from %s simulated decays (%.2f s), scaled to %.2f days", target, panelName.c_str(),
           grouped(decayEvents).c_str(), decayLive, slowDays),
      Form("2: one pixel each, sigma = %.2f keV x sqrt(E/MeV), %.2f entries per decay; L-alpha window %.2f +- %.2f keV", sigma1MeV,
           rEntries / activity, line, halfWindow),
      Form("3: %.0f..+%.0f ns of the middle TOF of %.1f million triggers; expected: triggers x %s /s x %.0f ns = %s", tLo, tHi,
           triggers / 1e6, grouped(rEntries).c_str(), tHi - tLo, grouped(expected).c_str()),
      "4: " + (split == std::string::npos ? cuts : cuts.substr(0, split + 1)),
      "    " + (split == std::string::npos ? std::string() : cuts.substr(split + 2) + " ") + "(cut profile " + profile.name + ")",
      Form("hatched: Lu-176 at %.0f..%.0f and +%.0f..+%.0f ns, x %.4f (%.0fx the statistics); after the cuts an upper estimate", s1Lo,
           s1Hi, s2Lo, s2Hi, sideScale, 1. / sideScale),
      Form("not in the bars: Lu-176 piled up on a muon's hit, %.1f entries in time, %.1f after the cuts", inWindow[kPileUp].w,
           afterCuts[kPileUp].w),
      Form("simulated: muons of 0-1 GeV %.2f days, of 1-1000 GeV %.2f days (background only, scaled to %.2f days)", bandDays[0],
           bandDays[1], slowDays)};
    // sizes in pixels as on event_flow.png (its 1040-pixel canvas); rows of three bars taller
    const int n = labels.size();
    const std::vector<double> rowPx = {120., 120., 120., 185., 185., 185.};
    const double topPx = 177., bottomPx = 140. + 44. * footer.size();
    double H = topPx + bottomPx;
    for (double r : rowPx) H += r;
    const double scale = 1040. / H;
    std::vector<double> edges = {0.};   // the bins from the bottom row up: row r is bin n - r
    for (int r = n - 1; r >= 0; --r) edges.push_back(edges.back() + rowPx[r]);
    auto c = new TCanvas("c_flow", "", 1600, int(std::lround(H)));
    c->SetLeftMargin(0.37); c->SetRightMargin(0.04); c->SetTopMargin(topPx / H); c->SetBottomMargin(bottomPx / H);
    c->SetLogx();
    gStyle->SetHatchesLineWidth(2);
    auto legend = new TLegend(0.02, 1. - 161. / H, 0.98, 1. - 109. / H);
    legend->SetNColumns(3); legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextFont(42); legend->SetTextSize(0.032 * scale);
    TLatex number;
    number.SetTextFont(42); number.SetTextSize(0.036 * scale); number.SetTextAlign(12);
    const double xMin = 0.05, xMax = 2e14;   // room for the 0.7 of the sidebands and right of the billions for their number
    TH1D* frame = nullptr;
    for (std::size_t s = 0; s < all.size(); ++s) {
      auto h = new TH1D(Form("flow_%zu", s), "", n, edges.data());
      for (int r = 0; r < n; ++r) {
        h->SetBinContent(n - r, all[s].values[r]);
        if (s == 0) h->GetXaxis()->SetBinLabel(n - r, labels[r].c_str());
      }
      h->SetFillColor(all[s].colour); h->SetFillStyle(all[s].fill); h->SetLineColor(all[s].colour);
      h->SetBarWidth(all[s].width); h->SetBarOffset(all[s].offset);
      if (s == 0) {
        h->SetMinimum(xMin); h->SetMaximum(xMax);
        h->GetXaxis()->SetLabelSize(0.045 * scale);
        h->GetYaxis()->SetTitle(Form("entries in %.2f days (step 1: decays)", slowDays));
        h->GetYaxis()->SetTitleSize(0.045 * scale); h->GetYaxis()->SetLabelSize(0); h->GetYaxis()->CenterTitle();
        h->GetYaxis()->SetTitleOffset(1.1);
        h->Draw("HBAR");
        frame = h;
      } else {
        h->Draw("HBAR SAME");
      }
      if (!all[s].name.empty()) legend->AddEntry(h, all[s].name.c_str(), "f");
      number.SetTextColor(all[s].colour);
      for (int r = 0; r < n; ++r)
        if (!all[s].text[r].empty())
          number.DrawLatex(1.4 * std::max(all[s].values[r], xMin), edges[n - r - 1] + (all[s].offset + 0.5 * all[s].width) * rowPx[r],
                           all[s].text[r].c_str());
    }
    legend->Draw();
    TLatex tick;   // the value axis labelled every third decade: fifteen labels run into each other
    tick.SetTextFont(42); tick.SetTextSize(frame->GetYaxis()->GetTitleSize()); tick.SetTextAlign(23);
    for (int k = 0; k <= 12; k += 3) tick.DrawLatex(std::pow(10., k), -12., k == 0 ? "1" : Form("10^{%d}", k));
    TLatex words;
    words.SetNDC(); words.SetTextFont(62); words.SetTextSize(0.05 * scale);
    words.DrawLatex(0.02, 1. - 73. / H, Form("From Lu-176 decays to the %s L-alpha window in %.2f days", target, slowDays));
    words.SetTextFont(42); words.SetTextSize(0.032 * scale); words.SetTextColor(kGray + 2);
    for (std::size_t k = 0; k < footer.size(); ++k) words.DrawLatex(0.02, (21. + 44. * (footer.size() - 1 - k)) / H, footer[k].c_str());
    Save(c, out + "/panel_intrinsic_flow");
  }
  {
    // the spectra of the steps, 10 keV bins (no empty ones after the cuts) shown per 2 keV as the other spectra
    auto c = bigCanvas("c_suppression");
    auto noTrigger = new TH1D("noTrigger", "", nBins, 20., drawnHigh);
    for (int i = 1; i <= nBins; ++i) noTrigger->SetBinContent(i, T * singles->GetBinContent(singles->FindBin(noTrigger->GetBinCenter(i))));
    auto inTime = dynamic_cast<TH1D*>(before[kIntrinsic]->Clone("luInTime"));
    auto xrAfter = dynamic_cast<TH1D*>(after[kTarget]->Clone("xrAfter"));
    sideAfterSpectrum->Scale(sideScale);
    struct Curve { TH1D* h; int colour, style; std::string name; };
    const std::vector<Curve> curves = {
      {noTrigger, darkViolet, 1, "the panel's Lu-176 alone, no trigger: " + rounded(noTrigger->Integral())},
      {inTime, violet, 1, "... in time with a TOF trigger: " + rounded(inTime->Integral())},
      {sideAfterSpectrum, violet, 2, "... after the cuts, from the sidebands: " + rounded(sideAfterSpectrum->Integral()) + " (counted " +
                                         rounded(after[kIntrinsic]->Integral()) + ")"},
      {xrAfter, red, 1, std::string(className[kTarget]) + " after the cuts: " + rounded(xrAfter->Integral())}};
    const int merge = nBins % 5 == 0 ? 5 : 1;
    double yLow = 1e30, yHigh = 0.;
    for (const auto& k : curves) {
      k.h->Rebin(merge); k.h->Scale(1. / merge);
      k.h->SetFillStyle(0); k.h->SetLineColor(k.colour); k.h->SetLineStyle(k.style); k.h->SetLineWidth(3);
      for (int i = 1; i <= k.h->GetNbinsX(); ++i)
        if (k.h->GetBinContent(i) > 0.) { yLow = std::min(yLow, k.h->GetBinContent(i)); yHigh = std::max(yHigh, k.h->GetBinContent(i)); }
    }
    noTrigger->SetTitle(Form(";single-pixel energy with the response (keV);entries per 2 keV in %.2f days", slowDays));
    noTrigger->SetMinimum(0.3 * yLow); noTrigger->SetMaximum(3e3 * yHigh);
    axes(noTrigger->GetXaxis(), noTrigger->GetYaxis());
    c->SetLeftMargin(0.095); noTrigger->GetYaxis()->SetTitleOffset(1.15);   // clear of the labels up to 10^12
    noTrigger->Draw("hist");
    c->Update();
    auto band = new TBox(line - halfWindow, std::pow(10., c->GetUymin()), line + halfWindow, std::pow(10., c->GetUymax()));
    band->SetFillColorAlpha(red, 0.15);
    band->Draw();
    for (const auto& k : curves) k.h->Draw("hist same");
    auto legend = new TLegend(0.45, 0.44, 0.97, 0.66);
    legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(0.03 * ts);
    for (const auto& k : curves) legend->AddEntry(k.h, k.name.c_str(), "l");
    legend->Draw();
    note.DrawLatex(0.125, 0.90, Form("the panel's Lu-176 against the TOF trigger; entries of 20-%.0f keV", drawnHigh));
    note.DrawLatex(0.125, 0.868, Form("4 cm %s target, %s; %.2f days", target, panelName.c_str(), slowDays));
    note.DrawLatex(0.125, 0.836, Form("no trigger: %s simulated decays (%.2f s), scaled; %.0f keV bins; shaded: the %s L_{#alpha} window",
                                     grouped(decayEvents).c_str(), decayLive, 2. * merge, target));
    // arrows down from each step to the next, left of the L-alpha window, with the share each keeps
    auto mean = [](TH1D* h, double lo, double hi) {
      const int a = h->FindBin(lo + 1e-6), b = h->FindBin(hi - 1e-6);
      return h->Integral(a, b) / (b - a + 1);
    };
    const double xa = 30., yFree = mean(noTrigger, 20., 40.), yTime = mean(inTime, 20., 40.), yCut = mean(sideAfterSpectrum, 20., 40.);
    TLatex step;
    step.SetTextFont(42); step.SetTextSize(0.03 * ts); step.SetTextAlign(12); step.SetTextColor(kGray + 3);
    auto arrow = [&](double y1, double y2, const std::string& text) {
      auto a = new TArrow(xa, y1, xa, y2, 0.012, "|>");
      a->SetLineWidth(2); a->SetLineColor(kGray + 2); a->SetFillColor(kGray + 2);
      a->Draw();
      step.DrawLatex(xa + 6., std::sqrt(y1 * y2), text.c_str());
    };
    // the shares of the drawn curves (20 keV up to the plot's end, as in the legend), the bins merged into fewer, each per 2 keV
    const double drawn[3] = {noTrigger->Integral(), inTime->Integral(), sideAfterSpectrum->Integral()};
    arrow(0.5 * yFree, 2. * yTime, "kept " + kept(drawn[1] / drawn[0]) + Form(": in time with a TOF trigger (%.2f /s x %.0f ns)",
                                                                              triggerPerMinute / 60., tHi - tLo));
    arrow(0.5 * yTime, 2. * yCut, "kept " + kept(drawn[2] / drawn[1]) + ": the cuts (after them from the sidebands)");
    Save(c, out + "/panel_intrinsic_suppression");
  }
}
