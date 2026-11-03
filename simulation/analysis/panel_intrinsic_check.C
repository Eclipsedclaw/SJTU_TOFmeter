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
//   the middle TOF of a triggered event (top and middle station within 20 ns, measured), 4. of these after the eight
//   cuts of the CZT analysis (event_flow.C), 5. of these in the target's L-alpha window (+-2 sigma of the panel).
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
//   root -l -b -q 'analysis/panel_intrinsic_check.C("output/lyso_concept", "output/lyso_concept/decay/decay_*.root", "output/lyso_concept/cosmic_al")'

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
void Save(TCanvas* c, const std::string& base)
{
  c->Print((base + ".png").c_str());
  c->Print((base + ".pdf").c_str());
}
}  // namespace

void panel_intrinsic_check(const char* dir = "output/lyso_concept", const char* decayFiles = "output/lyso_concept/decay/decay_*.root",
                           const char* runDir = "output/lyso_concept/cosmic_al", const char* target = "Al", double sigma1MeV = 26.85,
                           double triggerPerMinute = 493., double projectDays = 60., const char* outDir = "")
{
  const std::string out = std::string(outDir).empty() ? dir : outDir;
  const double line = std::string(target) == "Si" ? 76.24 : 65.73;   // the target's L-alpha as simulated (keV)
  const double halfWindow = 2. * sigma1MeV * std::sqrt(line / 1000.);   // +-2 sigma of the panel
  const double drawnHigh = std::string(target) == "Si" ? 500. : 450.;
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
  Count inWindow[kClasses], afterCuts[kClasses], lineWindow[kClasses], sideBefore, sideAfter, sideLine, peak;
  double slowDays = 0., bandDays[2] = {0., 0.};
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
      // the eight cuts of the CZT analysis (event_flow.C), kept when lo <= v <= hi
      const bool cuts = std::max(std::abs(xTgt), std::abs(yTgt)) <= 8. &&
                        std::atan(std::hypot(slopeX, slopeY)) * 180. / M_PI <= 35. && eTopM >= 1.6 && eTopM <= 3.5 &&
                        eMidM >= 1.6 && tkin >= 10. && tkin <= 80. && e <= 600. && eOther / 1000. <= 0.6 && npix <= 4;
      const bool inLine = std::abs(e - line) <= halfWindow;
      if (kind == kIntrinsic && ((dt >= s1Lo && dt <= s1Hi) || (dt >= s2Lo && dt <= s2Hi))) {
        sideBefore.Add(w);
        if (cuts) sideAfter.Add(w);
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
  o << Form("decay-only runs: %lld decays = %.3f s\n\n", decayEvents, decayLive);
  o << "rates\n";
  o << Form("   decays of the panel (activity)                      %10.0f /s   (%.1f /s per pixel of 900)\n", activity, activity / 900.);
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
  o << "   4. after the eight cuts\n";
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
  TLatex note;
  note.SetNDC();
  note.SetTextFont(42);
  note.SetTextSize(0.03);
  {
    auto c = new TCanvas("c_singles", "", 1000, 700);
    c->SetLogy(); c->SetLeftMargin(0.12); c->SetRightMargin(0.04);
    singles->SetTitle(";single-pixel energy with the response (keV);entries per second per 2 keV, all 900 pixels");
    singles->SetLineColor(colour[kIntrinsic]); singles->SetLineWidth(2);
    singles->GetXaxis()->SetRangeUser(10., 1300.);
    singles->Draw("hist");
    c->Update();
    auto band = new TBox(line - halfWindow, std::pow(10., c->GetUymin()), line + halfWindow, std::pow(10., c->GetUymax()));
    band->SetFillColorAlpha(colour[kTarget], 0.25);
    band->Draw();
    singles->Draw("hist same");
    note.DrawLatex(0.42, 0.85, Form("the panel's Lu-176 alone, no trigger: %.0f decays/s", activity));
    note.DrawLatex(0.42, 0.80, Form("%.0f entries/s of 20-1000 keV, %.0f /s in the %s L_{#alpha} window (shaded)", rEntries, rLine, target));
    note.DrawLatex(0.42, 0.75, Form("#sigma = %.2f keV #times #sqrt{E/MeV}, single pixels, 10 keV threshold", sigma1MeV));
    Save(c, out + "/panel_intrinsic_singles");
  }
  {
    auto c = new TCanvas("c_time", "", 1000, 700);
    c->SetLogy(); c->SetLeftMargin(0.12); c->SetRightMargin(0.04);
    timeMu->SetTitle(";entry time - middle-TOF time (ns);entries per 10 ns");
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
    auto legend = new TLegend(0.14, 0.76, 0.70, 0.89);
    legend->SetBorderSize(0); legend->SetFillStyle(0);
    legend->AddEntry(timeLu, "the panel's Lu-176 alone", "l");
    legend->AddEntry(timeMu, "made by the muon (incl. its X-rays; Lu-176 piled up on them)", "l");
    legend->Draw();
    note.DrawLatex(0.14, 0.71, Form("triggered events, entries of 20-1000 keV, before the cuts; dashed: the %.0f..+%.0f ns window", tLo, tHi));
    note.DrawLatex(0.14, 0.66, days.c_str());
    Save(c, out + "/panel_intrinsic_time");
  }
  for (int s = 0; s < 2; ++s) {
    TH1D** h = s ? after : before;
    auto c = new TCanvas(Form("c_spec%d", s), "", 1000, 700);
    c->SetLogy(); c->SetLeftMargin(0.12); c->SetRightMargin(0.04);
    auto stack = new THStack(Form("stack%d", s), Form(";single-pixel energy with the response (keV);entries per 2 keV"));
    auto legend = new TLegend(0.50, 0.62, 0.95, 0.89);
    legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(0.028);
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
    legend->Draw();
    note.DrawLatex(0.14, 0.85, s ? "within -5..+30 ns of a trigger, after the eight cuts" : "within -5..+30 ns of a trigger, before the cuts");
    note.DrawLatex(0.14, 0.80, Form("4 cm %s target, LYSO panel", target));
    note.DrawLatex(0.14, 0.75, Form("%.2f days", slowDays));
    Save(c, out + (s ? "/panel_intrinsic_after" : "/panel_intrinsic_before"));
  }
}
