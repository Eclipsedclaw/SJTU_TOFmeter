// The TOF-trigger ON/OFF subtraction of the steel study (2026-11-20: U, Th and K at 50 ppm each in the steel; can
// "prompt minus random" show their muonic X-rays, with the room's background included?), per day on the CZT panel,
// from two files:
//   traceFile  beam_steel_trace_spectrum.root of analysis/beam_steel_trace_spectrum.C, run with the fractions wanted
//              (50 ppm each: output/czt/beam_steel_czt_trace/analysis_50ppm/): the beam day's spectra per element and
//              view, the plate's natural decays at those fractions, the emitted lines;
//   roomFile   room_background_spectrum.root of analysis/room_background_spectrum.C (the concrete floor, radon in the
//              castle's air, cosmic rays, the walls' estimate); "" leaves the room out.
// The two spectra of the subtraction, per day (s = rate x promptNs, the chance that a hit unrelated to a muon falls into
// that muon's window; 1.32e-3 at 1.32e5 mu-/s and 10 ns):
//   OFF = the free-running spectrum of the same beam day: every hit, the beam's (all muons, every delay), the plate's
//         natural decays and the room's;
//   ON  = the TOF-gated spectrum: the hits within promptNs of their own muon, plus s x OFF (the accidentals);
//   ON - s x OFF = the muon-correlated prompt spectrum alone: what each muon makes in its own window.
// For a real day both are counted (Poisson): Var(ON - s OFF) = ON + s^2 OFF.
// For each muonic X-ray line of U, Th and K (the emitted lines merged within one resolution, yield >= 1%, from 50 keV,
// the 8 strongest, as beam_steel_trace_spectrum.C) in +-1.5 sigma: S the element's own X-ray entries, B the rest, in the
// three spectra (OFF, ON, ON - s OFF); z = S / sqrt(Var) with B known exactly, S/B, and f3, the accuracy to which B must
// be known for 3 sigma (S = 3 sqrt(Var + (f3 B)^2); none where the counting alone keeps S below 3 sigma). Per element
// the windows' z added in quadrature (the strongest first, a window that overlaps one taken is left out).
// A simulated day: ON and OFF drawn (Poisson; Gaussian above a mean of 1000) from the expected spectra; ON - s OFF minus
// the expected spectrum without the U, Th, K X-rays (the steel's own spectrum known exactly) leaves their X-rays and a
// day's noise. The expected spectra carry the simulation's own noise (1e8 simulated muons for 1.14e10 a day: weight
// ~114 for the steel's elements, ~5.7 for U, Th, K); it is not added to a day's variance.
// fracU, fracTh, fracK only label the outputs: they must be the fractions the trace file was made with.
// Outputs in outDir (default: the trace file's folder): beam_steel_onoff.txt, .csv, .root (the spectra per day in
// 0.5 keV bins, the simulated day) and one plot per file:
//   onoff_off_spectrum    OFF by part: the beam, the plate's natural U, Th, K, the room's parts
//   onoff_on_spectrum     ON by part: the muons' own prompt hits (the plate's X-rays, the rest), the accidentals
//   onoff_subtracted      ON - s OFF with the X-rays of U, Th and K
//   onoff_residual_k      the simulated day minus the steel's spectrum at K L_alpha 141.8 and K K_alpha 712.7 keV
//                         (two pads of the same day)
//   onoff_residual_u_th   the same at U 295.4 / Th 282.6 keV and at the K_alpha lines of U and Th, 6.0-6.5 MeV
//   onoff_lines           every U, Th, K line: z per day and S/B in the three spectra (two pads)
// Only reads the two files; writes into outDir. Run in the work folder:
//   root -l -b -q 'analysis/beam_steel_onoff.C("output/czt/beam_steel_czt_trace/analysis_50ppm/beam_steel_trace_spectrum.root", "output/czt/beam_steel_czt_room/analysis/room_background_spectrum.root")'
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
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
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {
void Save(TCanvas* c, const std::string& base)
{
  c->Print((base + ".png").c_str());
  c->Print((base + ".pdf").c_str());
}

std::string Sci(double v, int digits = 2)
{
  if (v == 0.) return "0";
  const int e = int(std::floor(std::log10(std::fabs(v))));
  char s[64];
  snprintf(s, sizeof s, "%.*f#times10^{%d}", digits - 1, v / std::pow(10., e), e);
  return s;
}

double WindowSum(const TH1D* h, double lo, double hi)
{
  double sum = 0.;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    const double x = h->GetBinCenter(b);
    if (x >= lo && x <= hi) sum += h->GetBinContent(b);
  }
  return sum;
}

TH1D* Clone(const TH1D* h, const std::string& name, double scale = 1., int rebin = 1)
{
  auto c = dynamic_cast<TH1D*>(h->Clone(name.c_str()));
  c->SetDirectory(nullptr);
  if (scale != 1.) c->Scale(scale);
  if (rebin > 1) c->Rebin(rebin);
  return c;
}

struct XLine {   // a muonic X-ray line of U, Th or K and its numbers in the three spectra (OFF, ON, ON - s OFF)
  int z;
  std::string el;
  double e, y, lo, hi;
  double S[3], B[3], V[3], zs[3], f3[3];
};
}  // namespace

void beam_steel_onoff(const char* traceFile, const char* roomFile = "", const char* outDir = "", double fracU = 5e-5,
                      double fracTh = 5e-5, double fracK = 5e-5, int fakeSeed = 20261120)
{
  gStyle->SetOptStat(0);
  TFile tf(traceFile);
  if (tf.IsZombie()) { printf("cannot open %s\n", traceFile); return; }
  const std::string out = std::string(outDir).empty() ? std::string(gSystem->GetDirName(traceFile).Data()) : std::string(outDir);
  gSystem->mkdir(out.c_str(), true);
  double sigma1MeV = 6.36991, threshold = 20., promptNs = 10., rate = 1.32e5, seconds = 86400.;
  if (auto st = dynamic_cast<TNamed*>(tf.Get("settings")))
    sscanf(st->GetTitle(), "sigma1MeV %lf threshold %lf promptNs %lf rate %lf seconds %lf", &sigma1MeV, &threshold, &promptNs, &rate,
           &seconds);
  const double s = rate * promptNs * 1e-9;
  auto Sigma = [&](double kev) { return sigma1MeV * std::sqrt(kev / 1000.); };
  auto Get = [](TFile& f, const std::string& name) -> TH1D* {
    auto h = dynamic_cast<TH1D*>(f.Get(name.c_str()));
    if (!h) { printf("no histogram %s in %s\n", name.c_str(), f.GetName()); return nullptr; }
    auto c = dynamic_cast<TH1D*>(h->Clone((name + "_in").c_str()));
    c->SetDirectory(nullptr);
    return c;
  };

  // the beam day (trace file)
  const std::vector<std::pair<int, std::string>> elements = {{26, "Fe"}, {24, "Cr"}, {28, "Ni"}, {42, "Mo"}, {25, "Mn"},
                                                             {14, "Si"}, {7, "N"},   {6, "C"},   {15, "P"},  {16, "S"},
                                                             {19, "K"},  {90, "Th"}, {92, "U"}};
  std::map<int, TH1D*> prompt, promptX, onX;
  for (const auto& [z, el] : elements) {
    prompt[z] = Get(tf, "day_" + el + "_prompt");
    promptX[z] = Get(tf, "day_" + el + "_promptx");
    onX[z] = Get(tf, "day_" + el + "_onx");
    if (!prompt[z] || !promptX[z] || !onX[z]) return;
  }
  TH1D *otherPrompt = Get(tf, "day_other_prompt"), *viewOn = Get(tf, "view_on"), *viewPrompt = Get(tf, "view_prompt"),
       *natural = Get(tf, "natural");
  if (!otherPrompt || !viewOn || !viewPrompt || !natural) return;
  const int nBins = viewOn->GetNbinsX();

  // the room
  const bool withRoom = !std::string(roomFile).empty();
  TH1D *room = nullptr, *roomFloor = nullptr, *roomRadon = nullptr, *roomCosmic = nullptr, *roomWalls = nullptr;
  if (withRoom) {
    TFile rf(roomFile);
    if (rf.IsZombie()) { printf("cannot open %s\n", roomFile); return; }
    room = Get(rf, "room_total"); roomFloor = Get(rf, "room_floor"); roomRadon = Get(rf, "room_radon");
    roomCosmic = Get(rf, "room_cosmic"); roomWalls = Get(rf, "room_walls");
    if (!room || !roomFloor || !roomRadon || !roomCosmic || !roomWalls) return;
    if (room->GetNbinsX() != nBins || room->GetXaxis()->GetXmax() != viewOn->GetXaxis()->GetXmax()) {
      printf("the room's histograms have other bins than the trace file's: stopping\n");
      return;
    }
  } else {
    room = Clone(viewOn, "room_none");
    room->Reset();
  }

  // the spectra of the subtraction
  TH1D* own = Clone(otherPrompt, "own_prompt");   // the muons' own prompt hits
  for (auto& [z, h] : prompt) own->Add(h);
  TH1D* off = Clone(viewOn, "off");
  off->Add(room);
  TH1D* on = Clone(own, "on");
  on->Add(off, s);
  TH1D* sub = Clone(on, "on_minus_off");
  sub->Add(off, -s);
  TH1D* var = Clone(on, "var_on_minus_off");   // contents only: ON + s^2 OFF
  var->Add(off, s * s);
  double maxDevOn = 0., maxDevSub = 0.;   // ON against the trace's own gated view; ON - s OFF against the own prompt hits
  for (int b = 1; b <= nBins; ++b) {
    const double ref = viewPrompt->GetBinContent(b) + s * room->GetBinContent(b), o = own->GetBinContent(b);
    if (ref > 0.) maxDevOn = std::max(maxDevOn, std::fabs(on->GetBinContent(b) - ref) / ref);
    if (o > 0.) maxDevSub = std::max(maxDevSub, std::fabs(sub->GetBinContent(b) - o) / o);
  }
  TH1D *xU = promptX[92], *xTh = promptX[90], *xK = promptX[19];
  TH1D* xSteel = Clone(promptX[26], "x_steel");   // the steel's own elements' muonic X-rays (all but U, Th, K)
  for (const auto& [z, el] : elements)
    if (z != 26 && z != 19 && z != 90 && z != 92) xSteel->Add(promptX[z]);
  TH1D* restPrompt = Clone(own, "rest_prompt");   // prompt hits that are not muonic X-rays of the plate
  for (const auto& [z, el] : elements) restPrompt->Add(promptX[z], -1.);
  TH1D* model = Clone(sub, "steel_without_u_th_k");   // the expected ON - s OFF without the X-rays of U, Th and K
  for (TH1D* h : {xU, xTh, xK}) model->Add(h, -1.);
  TH1D* beamAll = Clone(viewOn, "beam_all");   // the beam's hits in OFF (every muon, every delay)
  beamAll->Add(natural, -1.);

  // the lines of U, Th and K
  std::vector<XLine> xlines;
  {
    auto lt = dynamic_cast<TTree*>(tf.Get("lines"));
    if (!lt) { printf("no lines tree in %s\n", traceFile); return; }
    int lz = 0;
    double le = 0., lpc = 0.;
    lt->SetBranchAddress("z", &lz);
    lt->SetBranchAddress("e_keV", &le);
    lt->SetBranchAddress("per_capture", &lpc);
    std::map<int, std::vector<std::pair<double, double>>> raw;
    for (Long64_t i = 0; i < lt->GetEntries(); ++i) {
      lt->GetEntry(i);
      if (lz == 19 || lz == 90 || lz == 92) raw[lz].push_back({le, lpc});
    }
    for (const auto& [z, el] : std::vector<std::pair<int, std::string>>{{19, "K"}, {90, "Th"}, {92, "U"}}) {
      auto& r = raw[z];
      std::sort(r.begin(), r.end());
      std::vector<std::pair<double, double>> merged;
      for (const auto& [e0, y] : r) {
        if (!merged.empty() && e0 - merged.back().first < Sigma(e0)) {
          auto& m = merged.back();
          m.first = (m.first * m.second + e0 * y) / (m.second + y);
          m.second += y;
        } else merged.push_back({e0, y});
      }
      std::vector<std::pair<double, double>> kept;
      for (const auto& m : merged) if (m.second >= 0.01 && m.first >= 50.) kept.push_back(m);
      std::sort(kept.begin(), kept.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
      if (kept.size() > 8) kept.resize(8);
      std::sort(kept.begin(), kept.end());
      for (const auto& [e0, y] : kept) {
        XLine x;
        x.z = z; x.el = el; x.e = e0; x.y = y; x.lo = e0 - 1.5 * Sigma(e0); x.hi = e0 + 1.5 * Sigma(e0);
        xlines.push_back(x);
      }
    }
  }
  // the three spectra: OFF (its own X-rays: the ungated view's), ON, ON - s OFF
  const char* viewNames[3] = {"OFF (free-running)", "ON (TOF-gated)", "ON - s OFF"};
  const char* viewTags[3] = {"off", "on", "on_minus_off"};
  TH1D* totals[3] = {off, on, sub};
  TH1D* vars[3] = {off, on, var};
  for (XLine& x : xlines)
    for (int v = 0; v < 3; ++v) {
      const double S = WindowSum(v == 0 ? onX[x.z] : promptX[x.z], x.lo, x.hi), T = WindowSum(totals[v], x.lo, x.hi),
                   V = WindowSum(vars[v], x.lo, x.hi);
      x.S[v] = S; x.B[v] = std::max(0., T - S); x.V[v] = V;
      x.zs[v] = V > 0. ? S / std::sqrt(V) : 0.;
      const double need = (S / 3.) * (S / 3.) - V;
      x.f3[v] = need > 0. && x.B[v] > 0. ? std::sqrt(need) / x.B[v] : 0.;
    }

  // a simulated day
  TRandom3 rnd(fakeSeed);
  auto Draw = [&](double mu) { return mu <= 0. ? 0. : mu < 1000. ? double(rnd.Poisson(mu)) : rnd.Gaus(mu, std::sqrt(mu)); };
  TH1D* onDay = Clone(on, "day_on");
  TH1D* offDay = Clone(off, "day_off");
  for (int b = 1; b <= nBins; ++b) {
    onDay->SetBinContent(b, Draw(on->GetBinContent(b)));
    offDay->SetBinContent(b, Draw(off->GetBinContent(b)));
  }
  TH1D* subDay = Clone(onDay, "day_on_minus_off");
  subDay->Add(offDay, -s);
  TH1D* resid = Clone(subDay, "day_residual");   // the day minus the steel's own spectrum
  resid->Add(model, -1.);
  for (int b = 1; b <= nBins; ++b) resid->SetBinError(b, std::sqrt(std::max(0., var->GetBinContent(b))));
  TH1D* expected = Clone(xU, "expected_u_th_k");
  expected->Add(xTh); expected->Add(xK);

  // the ROOT file
  {
    TFile of((out + "/beam_steel_onoff.root").c_str(), "RECREATE");
    for (TH1D* h : {off, on, sub, var, own, xSteel, restPrompt, model, beamAll, room, onDay, offDay, subDay, resid, expected}) h->Write();
    Clone(xU, "x_U")->Write(); Clone(xTh, "x_Th")->Write(); Clone(xK, "x_K")->Write();
    TNamed("settings", Form("s %g rate %g promptNs %g seconds %g sigma1MeV %g threshold %g fracU %g fracTh %g fracK %g trace %s room %s seed %d", s,
                            rate, promptNs, seconds, sigma1MeV, threshold, fracU, fracTh, fracK, traceFile, withRoom ? roomFile : "none",
                            fakeSeed)).Write();
  }

  // numbers, step by step
  std::ostringstream txt, csv;
  char line[800];
  const double lo0 = threshold, hi0 = 8000.;
  auto All = [&](const TH1D* h) { return WindowSum(h, lo0, hi0); };
  const double offT = All(off), onT = All(on), ownT = All(own);
  txt << "TOF-trigger ON/OFF subtraction (prompt minus random), one day of " << seconds << " s on the CZT 10 x 10 cm panel\n";
  txt << "  beam day: " << traceFile << "\n  room: " << (withRoom ? roomFile : "not included") << "\n";
  snprintf(line, sizeof line, "  steel with U %.3g, Th %.3g, K %.3g ppm (mass); %.3g mu-/s, prompt window %g ns: s = rate x window = %.4g, "
           "the chance that a hit unrelated to a muon falls into that muon's window\n", 1e6 * fracU, 1e6 * fracTh, 1e6 * fracK, rate, promptNs, s);
  txt << line;
  txt << "\nstep 1, the two spectra of the same beam day (entries per day, " << lo0 << "-" << hi0 << " keV):\n";
  snprintf(line, sizeof line, "  OFF, free-running: %.4g (%.4g /s)\n    the beam's hits (every muon, every delay): %.4g (%.3g%%)\n"
           "    the plate's natural U, Th, K: %.4g (%.3g%%)\n", offT, offT / seconds, All(beamAll), 100. * All(beamAll) / offT, All(natural),
           100. * All(natural) / offT);
  txt << line;
  if (withRoom) {
    snprintf(line, sizeof line, "    the room: %.4g (%.3g%%): floor %.4g, radon %.4g, cosmic rays %.4g, walls and ceiling %.4g\n", All(room),
             100. * All(room) / offT, All(roomFloor), All(roomRadon), All(roomCosmic), All(roomWalls));
    txt << line;
  }
  snprintf(line, sizeof line, "  ON, TOF-gated: %.4g\n    the muons' own prompt hits: %.4g (%.3g%%): the plate's muonic X-rays %.4g "
           "(U %.4g, Th %.4g, K %.4g), the rest %.4g\n    accidentals s x OFF: %.4g (%.3g%%): beam %.4g, natural %.4g, room %.4g\n",
           onT, ownT, 100. * ownT / onT, ownT - All(restPrompt), All(xU), All(xTh), All(xK), All(restPrompt), s * offT, 100. * s * offT / onT,
           s * All(beamAll), s * All(natural), s * All(room));
  txt << line;
  txt << "\nstep 2, the subtraction ON - s x OFF:\n";
  snprintf(line, sizeof line, "  ON agrees with the trace file's gated view (+ s x room) to %.2g per bin; ON - s OFF equals the muons' own "
           "prompt hits to %.2g per bin (rounding)\n", maxDevOn, maxDevSub);
  txt << line;
  snprintf(line, sizeof line, "  removed on average: %.4g entries (%.3g%% of ON), of them the room %.4g and the plate's natural decays %.4g\n",
           s * offT, 100. * s * offT / onT, s * All(room), s * All(natural));
  txt << line;
  snprintf(line, sizeof line, "  the price: the variance ON + s^2 OFF = %.6g x ON (the subtraction's own noise, s^2 OFF, is %.2g of it)\n",
           All(var) / onT, s * s * offT / All(var));
  txt << line;
  txt << "\nstep 3, the muonic X-ray lines of U, Th and K (+-1.5 sigma; S the element's own X-ray entries, B the rest; "
         "z = S / sqrt(Var) with B known exactly; f3: B known to this fraction for 3 sigma, '-' = not even with B exact):\n";
  txt << "  line (keV)  yield |  OFF: S, B, z, S/B         |  ON: S, B, z, S/B          |  ON - s OFF: S, B, z, S/B, f3\n";
  csv << "kind,element,line_keV,yield,view,window_lo,window_hi,S,B,var,z,S_over_B,f3\n";
  for (const XLine& x : xlines) {
    txt << Form("  %-2s %8.2f %6.3f", x.el.c_str(), x.e, x.y);
    for (int v = 0; v < 3; ++v) {
      txt << Form(" | %8.3g %8.3g %6.3f %8.2g", x.S[v], x.B[v], x.zs[v], x.B[v] > 0. ? x.S[v] / x.B[v] : 0.);
      csv << Form("line,%s,%.2f,%.4f,%s,%.2f,%.2f,%.6g,%.6g,%.6g,%.5g,%.5g,%.5g\n", x.el.c_str(), x.e, x.y, viewTags[v], x.lo, x.hi, x.S[v],
                  x.B[v], x.V[v], x.zs[v], x.B[v] > 0. ? x.S[v] / x.B[v] : 0., x.f3[v]);
    }
    txt << (x.f3[2] > 0. ? Form(" %8.2g", x.f3[2]) : "        -") << "\n";
  }
  txt << "\nstep 4, per element: the lines' z added in quadrature (non-overlapping windows), per day:\n";
  txt << "  element | OFF  | ON   | ON - s OFF | lines used (ON - s OFF)\n";
  for (const std::string el : {"K", "Th", "U"}) {
    double zc[3] = {0., 0., 0.};
    std::string used;
    for (int v = 0; v < 3; ++v) {
      std::vector<const XLine*> ls;
      for (const XLine& x : xlines) if (x.el == el) ls.push_back(&x);
      std::sort(ls.begin(), ls.end(), [v](const XLine* a, const XLine* b) { return a->zs[v] > b->zs[v]; });
      std::vector<const XLine*> taken;
      for (const XLine* x : ls) {
        bool overlap = false;
        for (const XLine* t : taken) if (x->lo < t->hi && x->hi > t->lo) overlap = true;
        if (overlap) continue;
        taken.push_back(x);
        zc[v] += x->zs[v] * x->zs[v];
        if (v == 2) used += Form("%s%.1f", used.empty() ? "" : " ", x->e);
      }
      zc[v] = std::sqrt(zc[v]);
    }
    txt << Form("  %-7s | %4.2f | %4.2f | %10.2f | %s\n", el.c_str(), zc[0], zc[1], zc[2], used.c_str());
    csv << Form("combined,%s,,,%s,,,,,,%.5g,,\ncombined,%s,,,%s,,,,,,%.5g,,\ncombined,%s,,,%s,,,,,,%.5g,,\n", el.c_str(), viewTags[0], zc[0],
                el.c_str(), viewTags[1], zc[1], el.c_str(), viewTags[2], zc[2]);
  }
  txt << "\nstep 5, one simulated day (seed " << fakeSeed << "): ON - s OFF minus the steel's own spectrum (known exactly), in each line's window:\n";
  txt << "  line (keV) | excess measured +- error | expected S (all three elements' X-rays) | z measured | a 0.1% error of the steel's spectrum there\n";
  for (const XLine& x : xlines) {
    const double ex = WindowSum(resid, x.lo, x.hi), sg = std::sqrt(WindowSum(var, x.lo, x.hi)), expS = WindowSum(expected, x.lo, x.hi),
                 m = 1e-3 * WindowSum(model, x.lo, x.hi);
    txt << Form("  %-2s %8.2f | %10.4g +- %9.3g | %10.4g | %6.2f | %9.3g (%.2g x S)\n", x.el.c_str(), x.e, ex, sg, expS, sg > 0. ? ex / sg : 0., m,
                expS > 0. ? m / expS : 0.);
    csv << Form("day,%s,%.2f,%.4f,on_minus_off,%.2f,%.2f,%.6g,%.6g,%.6g,%.5g,,\n", x.el.c_str(), x.e, x.y, x.lo, x.hi, ex, m, sg * sg,
                sg > 0. ? ex / sg : 0.);
  }
  if (withRoom) {
    txt << "\nstep 6, the room's part in each line's window (per day): in OFF | in ON (s x room) | in ON, share of B\n";
    for (const XLine& x : xlines) {
      const double r = WindowSum(room, x.lo, x.hi), n = WindowSum(natural, x.lo, x.hi);
      txt << Form("  %-2s %8.2f | room %9.3g = %.2g of OFF, natural %9.3g = %.2g of OFF | room %8.3g, natural %8.3g | %.2g, %.2g\n",
                  x.el.c_str(), x.e, r, r / WindowSum(off, x.lo, x.hi), n, n / WindowSum(off, x.lo, x.hi), s * r, s * n,
                  x.B[1] > 0. ? s * r / x.B[1] : 0., x.B[1] > 0. ? s * n / x.B[1] : 0.);
    }
  }
  // beam off, for comparison: the plate's own radioactivity at these fractions, without and with the room
  {
    struct NatLine { double e; const char* nuclide; const char* chain; };
    const std::vector<NatLine> nat = {{238.63, "Pb-212", "nat_90_232"}, {351.93, "Pb-214", "nat_92_238"}, {609.32, "Bi-214", "nat_92_238"},
                                      {911.20, "Ac-228", "nat_90_232"}, {1460.82, "K-40", "nat_19_40"},   {2614.51, "Tl-208", "nat_90_232"}};
    txt << "\nstep 7, beam off for comparison (one day): the plate's natural lines (+-1.5 sigma; S the line's chain, B the rest): "
           "without the room | with the room\n";
    for (const NatLine& nl : nat) {
      TH1D* chain = Get(tf, nl.chain);
      if (!chain) continue;
      const double sg = Sigma(nl.e), lo = nl.e - 1.5 * sg, hi = nl.e + 1.5 * sg;
      const double S = WindowSum(chain, lo, hi), B0 = std::max(0., WindowSum(natural, lo, hi) - S), R = WindowSum(room, lo, hi);
      txt << Form("  %8.2f %-7s | S %9.3g, B %9.3g, z %7.3g | room %9.3g: B %9.3g, z %7.3g, S/B %.2g\n", nl.e, nl.nuclide, S, B0,
                  S / std::sqrt(S + B0), R, B0 + R, S / std::sqrt(S + B0 + R), B0 + R > 0. ? S / (B0 + R) : 0.);
      csv << Form("beam_off,%s,%.2f,,with_room,%.2f,%.2f,%.6g,%.6g,%.6g,%.5g,%.5g,\n", nl.nuclide, nl.e, lo, hi, S, B0 + R, S + B0 + R,
                  S / std::sqrt(S + B0 + R), B0 + R > 0. ? S / (B0 + R) : 0.);
      delete chain;
    }
  }

  // plots
  const std::string beamText = Form("%s #mu^{-}/s", Sci(rate).c_str());
  const std::string ppmText = Form("U, Th, K %.3g, %.3g, %.3g ppm", 1e6 * fracU, 1e6 * fracTh, 1e6 * fracK);
  struct Part { TH1D* h; std::string label; int color, style; bool overTotal = false; };   // overTotal: drawn over the total, dashed
  auto SpectrumPlot = [&](const std::string& tag, const std::string& title, TH1D* total, const std::vector<Part>& parts) {
    const int rebin = 16;   // 8 keV
    auto c = new TCanvas(("c_" + tag).c_str(), "", 1100, 750);
    c->SetLogy(); c->SetLeftMargin(0.11); c->SetRightMargin(0.03); c->SetTopMargin(0.08);
    TH1D* t = Clone(total, tag + "_total", 1., rebin);
    t->SetTitle((title + ";energy in one pixel (keV);entries per 8 keV in one day").c_str());
    t->GetXaxis()->SetRangeUser(20., 7000.);
    t->SetLineColor(kBlack); t->SetLineWidth(2);
    const double top = t->GetMaximum();
    t->SetMaximum(top * 3000.); t->SetMinimum(top / 2e9);
    t->Draw("hist");
    auto leg = new TLegend(0.13, 0.72, 0.97, 0.91);
    leg->SetNColumns(2); leg->SetMargin(0.12); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.021);
    leg->AddEntry(t, ("all entries: " + Sci(WindowSum(total, 20., 7000.)) + " in 20-7000 keV").c_str(), "l");
    int i = 0;
    std::vector<TH1D*> over;
    for (const Part& p : parts) {
      TH1D* h = Clone(p.h, tag + "_p" + std::to_string(i++), 1., rebin);
      h->GetXaxis()->SetRangeUser(20., 7000.);
      h->SetLineColor(p.color); h->SetLineWidth(2); h->SetLineStyle(p.style);
      if (p.overTotal) over.push_back(h); else h->Draw("hist same");
      leg->AddEntry(h, (p.label + ": " + Sci(WindowSum(p.h, 20., 7000.))).c_str(), "l");
    }
    t->Draw("hist same");
    for (TH1D* h : over) h->Draw("hist same");
    leg->Draw();
    Save(c, out + "/" + tag);
  };
  {
    std::vector<Part> parts = {{beamAll, "the beam's hits (every muon, every delay; dashed over the total)", kOrange + 7, 2, true},
                               {natural, "the plate's natural U, Th, K decays", kYellow + 3, 1}};
    if (withRoom) {
      parts.push_back({roomFloor, "room: concrete floor (K, U, Th)", kCyan + 2, 1});
      parts.push_back({roomRadon, "room: radon in the castle's air", kGreen + 2, 1});
      parts.push_back({roomCosmic, "room: cosmic rays", kGray + 2, 1});
      parts.push_back({roomWalls, "room: walls and ceiling through the lead (estimate)", kPink + 1, 2});
    }
    parts.push_back({onX[92], "muonic X-rays of U", kViolet, 1});
    parts.push_back({onX[90], "muonic X-rays of Th", kAzure + 1, 1});
    parts.push_back({onX[19], "muonic X-rays of K", kOrange, 1});
    SpectrumPlot("onoff_off_spectrum", Form("OFF: the free-running spectrum of one day at %s, 1 cm steel with %s, CZT", beamText.c_str(), ppmText.c_str()),
                 off, parts);
  }
  {
    TH1D* accBeam = Clone(beamAll, "acc_beam", s);
    TH1D* accNat = Clone(natural, "acc_nat", s);
    std::vector<Part> parts = {{xSteel, "the muons' own prompt hits: muonic X-rays of the steel's elements", kRed, 1},
                               {restPrompt, "the muons' own prompt hits: the rest (capture #gamma, e^{-}, n)", kOrange + 7, 1},
                               {xU, "the muons' own prompt hits: muonic X-rays of U", kViolet, 1},
                               {xTh, "the muons' own prompt hits: muonic X-rays of Th", kAzure + 1, 1},
                               {xK, "the muons' own prompt hits: muonic X-rays of K", kOrange, 1},
                               {accBeam, Form("accidentals: other muons' hits (s = %.3g)", s), kGray + 1, 2},
                               {accNat, "accidentals: the plate's natural decays", kYellow + 3, 2}};
    if (withRoom) parts.push_back({Clone(room, "acc_room", s), "accidentals: the room", kCyan + 2, 2});
    SpectrumPlot("onoff_on_spectrum", Form("ON: the TOF-gated spectrum (%g ns) of the same day at %s, %s", promptNs, beamText.c_str(), ppmText.c_str()),
                 on, parts);
  }
  {
    std::vector<Part> parts = {{xSteel, "muonic X-rays of the steel's elements", kRed, 1},
                               {restPrompt, "the rest of the muons' prompt hits (capture #gamma, e^{-}, n)", kOrange + 7, 1},
                               {xU, "muonic X-rays of U", kViolet, 1},
                               {xTh, "muonic X-rays of Th", kAzure + 1, 1},
                               {xK, "muonic X-rays of K", kOrange, 1}};
    SpectrumPlot("onoff_subtracted", Form("ON - s OFF: what each muon makes within %g ns, one day at %s, %s", promptNs, beamText.c_str(),
                                          ppmText.c_str()),
                 sub, parts);
  }
  // the simulated day around the lines: two pads per plot
  auto ResidualPad = [&](const std::string& tag, double lo, double hi, int rebin, const std::vector<std::tuple<TH1D*, std::string, int>>& sig,
                         const std::vector<const XLine*>& marks) {
    TH1D* r = Clone(resid, tag + "_r", 1., rebin);
    TH1D* e = Clone(expected, tag + "_e", 1., rebin);
    TH1D* m = Clone(model, tag + "_m", 1e-3, rebin);   // 0.1% of the steel's spectrum
    r->GetXaxis()->SetRangeUser(lo, hi);
    double ymax = 0., ymin = 0.;   // from the day and the expected X-rays (the 0.1% curve may leave the frame)
    for (int b = r->FindBin(lo); b <= r->FindBin(hi); ++b) {
      ymax = std::max({ymax, r->GetBinContent(b) + r->GetBinError(b), e->GetBinContent(b)});
      ymin = std::min(ymin, r->GetBinContent(b) - r->GetBinError(b));
    }
    r->SetMaximum(ymin + (ymax - ymin) * 2.1); r->SetMinimum(ymin - 0.05 * (ymax - ymin));
    r->SetTitle(Form(";energy in one pixel (keV);entries per %g keV in one day", rebin * 0.5));
    r->SetMarkerStyle(20); r->SetMarkerSize(0.6); r->SetLineColor(kBlack); r->SetMarkerColor(kBlack);
    r->GetYaxis()->SetTitleOffset(1.25); r->GetYaxis()->SetMaxDigits(3);
    r->Draw("e1");
    m->SetLineColor(kGray + 1); m->SetLineStyle(2); m->SetLineWidth(2);
    m->GetXaxis()->SetRangeUser(lo, hi);
    m->Draw("hist same");
    auto leg = new TLegend(0.12, 0.70, 0.97, 0.90);
    leg->SetNColumns(2); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
    leg->AddEntry(r, "one day: (ON - s OFF) - the steel's own spectrum", "lep");
    leg->AddEntry(m, "0.1% of the steel's own spectrum", "l");
    for (const auto& [h, label, color] : sig) {
      TH1D* x = Clone(h, tag + "_" + label, 1., rebin);
      x->GetXaxis()->SetRangeUser(lo, hi);
      x->SetLineColor(color); x->SetLineWidth(2);
      x->Draw("hist same");
      leg->AddEntry(x, (label + " (expected)").c_str(), "l");
    }
    r->Draw("e1 same");
    TLine zero(lo, 0., hi, 0.);
    zero.SetLineColor(kGray + 2); zero.DrawClone();
    leg->Draw();
    TLatex tex;
    tex.SetNDC(); tex.SetTextSize(0.032); tex.SetTextAlign(13);
    double y = 0.66;
    for (const XLine* x : marks) {
      const double ex = WindowSum(resid, x->lo, x->hi), sg = std::sqrt(WindowSum(var, x->lo, x->hi));
      tex.DrawLatex(0.14, y, Form("%s %.1f keV: excess %s #pm %s (%.1f#sigma); expected %s (%s's own %s); B %s, S/B %.1e", x->el.c_str(),
                                  x->e, Sci(ex).c_str(), Sci(sg).c_str(), ex / sg, Sci(WindowSum(expected, x->lo, x->hi)).c_str(), x->el.c_str(),
                                  Sci(x->S[2]).c_str(), Sci(x->B[2]).c_str(), x->B[2] > 0. ? x->S[2] / x->B[2] : 0.));
      y -= 0.05;
    }
  };
  auto Find = [&](const std::string& el, double e) -> const XLine* {
    const XLine* best = nullptr;
    for (const XLine& x : xlines) if (x.el == el && (!best || std::fabs(x.e - e) < std::fabs(best->e - e))) best = &x;
    return best;
  };
  {
    auto c = new TCanvas("c_res_k", "", 1100, 1100);
    c->Divide(1, 2);
    c->cd(1); gPad->SetLeftMargin(0.11); gPad->SetRightMargin(0.03);
    ResidualPad("rk1", 120., 165., 2, {{xK, "K muonic X-rays", kOrange + 1}}, {Find("K", 141.8)});
    TLatex t; t.SetNDC(); t.SetTextSize(0.04); t.SetTextAlign(13);
    t.DrawLatex(0.11, 0.985, Form("A simulated day, ON - s OFF minus the steel's own spectrum (%s): K L_{#alpha} and K_{#alpha}", ppmText.c_str()));
    c->cd(2); gPad->SetLeftMargin(0.11); gPad->SetRightMargin(0.03);
    ResidualPad("rk2", 680., 745., 4, {{xK, "K muonic X-rays", kOrange + 1}}, {Find("K", 712.7)});
    Save(c, out + "/onoff_residual_k");
  }
  {
    auto c = new TCanvas("c_res_uth", "", 1100, 1100);
    c->Divide(1, 2);
    c->cd(1); gPad->SetLeftMargin(0.11); gPad->SetRightMargin(0.03);
    ResidualPad("ru1", 262., 316., 2, {{xTh, "Th muonic X-rays", kAzure + 1}, {xU, "U muonic X-rays", kViolet}}, {Find("Th", 282.6), Find("U", 295.4)});
    TLatex t; t.SetNDC(); t.SetTextSize(0.04); t.SetTextAlign(13);
    t.DrawLatex(0.11, 0.985, Form("A simulated day, ON - s OFF minus the steel's own spectrum (%s): U and Th", ppmText.c_str()));
    c->cd(2); gPad->SetLeftMargin(0.11); gPad->SetRightMargin(0.03);
    ResidualPad("ru2", 5980., 6560., 24, {{xTh, "Th muonic X-rays", kAzure + 1}, {xU, "U muonic X-rays", kViolet}}, {Find("Th", 6367.7), Find("U", 6420.)});
    Save(c, out + "/onoff_residual_u_th");
  }
  // every line: z and S/B in the three spectra
  {
    const int n = xlines.size();
    auto c = new TCanvas("c_lines", "", 1100, 900);
    c->Divide(1, 2);
    const int styles[3] = {24, 25, 20}, colors[3] = {kGray + 2, kBlue + 1, kBlack};
    for (int p = 0; p < 2; ++p) {
      c->cd(p + 1);
      gPad->SetLogy(); gPad->SetLeftMargin(0.09); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.19); gPad->SetTopMargin(p == 0 ? 0.11 : 0.05);
      gPad->SetGridy();
      auto frame = new TH1D(Form("frame%d", p), "", n, 0., n);
      frame->SetDirectory(nullptr);
      for (int i = 0; i < n; ++i) frame->GetXaxis()->SetBinLabel(i + 1, Form("%s %.1f", xlines[i].el.c_str(), xlines[i].e));
      frame->GetXaxis()->LabelsOption("v");
      frame->GetXaxis()->SetLabelSize(0.045);
      frame->SetMinimum(p == 0 ? 0.05 : 1e-6); frame->SetMaximum(p == 0 ? 30. : 0.1);
      frame->SetTitle(p == 0 ? ";;z = S / #sqrt{Var} (B known exactly)" : ";;S / B");   // the title over both pads
      frame->GetYaxis()->SetTitleSize(0.05); frame->GetYaxis()->SetTitleOffset(0.8);
      frame->Draw("axis");
      if (p == 0) {
        TLine three(0., 3., n, 3.);
        three.SetLineStyle(2); three.SetLineColor(kRed); three.DrawClone();
      }
      auto leg = new TLegend(0.10, p == 0 ? 0.78 : 0.80, 0.97, 0.90);
      leg->SetNColumns(3); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.045);
      for (int v = 0; v < 3; ++v) {
        auto g = new TGraph();
        for (int i = 0; i < n; ++i) {
          const double val = p == 0 ? xlines[i].zs[v] : (xlines[i].B[v] > 0. ? xlines[i].S[v] / xlines[i].B[v] : 0.);
          if (val > 0.) g->SetPoint(g->GetN(), i + 0.5 + 0.18 * (v - 1), val);
        }
        g->SetMarkerStyle(styles[v]); g->SetMarkerColor(colors[v]); g->SetMarkerSize(1.2);
        g->Draw("p same");
        leg->AddEntry(g, viewNames[v], "p");
      }
      leg->Draw();
    }
    c->cd();
    TLatex title;
    title.SetNDC(); title.SetTextSize(0.022); title.SetTextAlign(23);
    title.DrawLatex(0.5, 0.993, Form("The muonic X-ray lines of %s, one day at %s: z and S/B in OFF, ON and ON - s OFF", ppmText.c_str(),
                                     beamText.c_str()));
    Save(c, out + "/onoff_lines");
  }

  std::ofstream(out + "/beam_steel_onoff.txt") << txt.str();
  std::ofstream(out + "/beam_steel_onoff.csv") << csv.str();
  printf("%s", txt.str().c_str());
  printf("written into %s: beam_steel_onoff.txt, .csv, .root, onoff_off_spectrum, onoff_on_spectrum, onoff_subtracted, onoff_residual_k, "
         "onoff_residual_u_th, onoff_lines (.png, .pdf)\n", out.c_str());
}
