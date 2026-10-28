// Thicker Al targets, the selection taken apart. analysis/muonic_xray_target_compare.C puts every thickness through the
// cuts of the 1 cm target, and their efficiency falls with the depth of the stop in the target. From the same skims (the
// 1 cm target of the 3-week run, the runs of scripts/muonic_xray_target.sh):
//   root -l -b -q 'analysis/muonic_xray_target_cuts.C'
// 1. Which cut it is: the target's L_alpha entries of the thickest target by the depth of the stop, the share with a
//    middle-TOF hit in coincidence, and the share of those that pass each further cut on its own.
// 2. Why: the measured middle-TOF energy and top-to-middle time of the muon against that depth (a muon that stops deeper
//    was faster), beside those of the background of the fast muons.
// 3. Every thickness with its own slow-muon cuts: a scan of the cuts on the middle-TOF energy and on the top-to-middle
//    time, the other cuts as in muonic_xray_stack_cuts.C, for the best S/sqrt(B) in 21 days at Al L_alpha (with the
//    fast-muon background, primaries above eSplit, of the 1 cm run). Cross-validated: the best cuts of the even events
//    are scored on the odd ones and the other way round, so that a lucky fluctuation is not counted as a gain.
// 4. Whether the fast-muon background of the 1 cm run holds for a thick target: the background of the primaries of
//    0.5 GeV to eSplit, which cross every target as the faster ones do, per thickness.
// Writes <dir>/muonic_xray_target_cuts.png/.pdf/.csv.

#include <algorithm>
#include <fstream>
#include <sstream>

namespace {
// as analysis/muonic_xray_stack_cuts.C
double SigmaKeV(double eKeV) { return 40. * std::sqrt(eKeV / 1000.); }

struct Window { const char* name; double eKeV; };
const std::vector<Window> kWindows = {{"Al L_alpha", 65.73}, {"Al L_beta", 88.73}, {"Al K_alpha", 354.9}};

bool InWindow(double eKeV, const Window& w) { return std::abs(eKeV - w.eKeV) <= 2. * SigmaKeV(w.eKeV); }

// an entry in a window that passes the cuts not scanned: middle TOF in coincidence, no bottom TOF, middle bars 2-6,
// highest hit in the top two layers
struct Entry {
  int kind;           // 0 X-ray of the target; else 1 a primary up to eSplit, 2 a faster one
  int windows;        // bit w: in kWindows[w]
  int half;           // parity of the event, for the cross-validation
  bool top;           // the top TOF timed the muon
  float eMid, dt;     // measured middle-TOF energy (MeV), top-to-middle time (ns)
  float ekin, depth;  // truth: primary energy (MeV), depth of the stop below the top of the target (cm)
};

struct Sample {
  double cm = 1., top = 36.5, days = 0.;
  TString file;
  std::vector<Entry> entries;
  TH1D* depth[8] = {};   // the target's L_alpha entries by depth: all; coincidence; + each further cut alone; all cuts
};

double Quantile(std::vector<double> v, double q)
{
  if (v.empty()) return 0.;
  std::sort(v.begin(), v.end());
  return v[std::min<std::size_t>(v.size() - 1, std::size_t(q * v.size()))];
}
}  // namespace

void muonic_xray_target_cuts(const char* ref = "output/muonic_xray_stack/xray_skim_pbke.root",
                             const char* dir = "output/muonic_xray_target", const char* thicknesses = "2 3 4 5",
                             double eSplit = 1000., double sigmaT = 0.2, double sigmaE = 0.15, double eMidCut = 3.0,
                             double tofCut = 1.3)
{
  TH1::AddDirectory(false);
  std::vector<Sample> samples(1);
  samples[0].file = ref;
  std::istringstream list(thicknesses);
  for (double d; list >> d;) {
    Sample s;
    s.cm = d;
    s.top = std::max(36.5, 31.7 + d);   // as scripts/muonic_xray_target.sh
    s.file = Form("%s/t%gcm/xray_skim.root", dir, d);
    samples.push_back(s);
  }

  for (auto& s : samples) {
    TFile file(s.file);
    auto cand = dynamic_cast<TTree*>(file.Get("cand"));
    auto info = dynamic_cast<TTree*>(file.Get("info"));
    if (!cand || !info) {
      printf("no skim in %s (run analysis/muonic_xray_stack_skim.C)\n", s.file.Data());
      continue;
    }
    double lt = 0.;
    info->SetBranchAddress("live_time", &lt);
    for (Long64_t i = 0; i < info->GetEntries(); ++i) { info->GetEntry(i); s.days += lt / 86400.; }
    int proc, event, topLayer, src, tofTop, tofMid, barMid;
    float e, t, tTop, tMid, eMid, botFirst, ekin, endZ;
    cand->SetBranchStatus("*", 0);
    for (const char* b : {"proc", "event", "e", "topLayer", "t", "src", "ekin", "endZ", "tofTop", "tofMid", "tTop", "tMid",
                          "barMid", "eMid", "botFirst"})
      cand->SetBranchStatus(b, 1);
    cand->SetBranchAddress("proc", &proc); cand->SetBranchAddress("event", &event);
    cand->SetBranchAddress("e", &e); cand->SetBranchAddress("topLayer", &topLayer); cand->SetBranchAddress("t", &t);
    cand->SetBranchAddress("src", &src); cand->SetBranchAddress("ekin", &ekin); cand->SetBranchAddress("endZ", &endZ);
    cand->SetBranchAddress("tofTop", &tofTop); cand->SetBranchAddress("tofMid", &tofMid);
    cand->SetBranchAddress("tTop", &tTop); cand->SetBranchAddress("tMid", &tMid);
    cand->SetBranchAddress("barMid", &barMid); cand->SetBranchAddress("eMid", &eMid); cand->SetBranchAddress("botFirst", &botFirst);
    for (int k = 0; k < 8; ++k) s.depth[k] = new TH1D(Form("depth%g_%d", s.cm, k), "", 10, 0., 5.);
    TRandom3 random;
    for (Long64_t i = 0; i < cand->GetEntries(); ++i) {
      cand->GetEntry(i);
      const bool signal = src & 1, alpha = signal && InWindow(e, kWindows[0]);
      const double depth = s.top - endZ;
      if (alpha) s.depth[0]->Fill(depth);
      // the selection of muonic_xray_stack_cuts.C, with its TOF response and seeds
      if (!tofMid || t - tMid < -5. || t - tMid > 30.) continue;
      random.SetSeed(1000003ULL * proc + event + 1);
      const double eMidMeasured = eMid * (1. + sigmaE * random.Gaus());
      const double dtTof = tofTop ? (tMid + sigmaT * random.Gaus()) - (tTop + sigmaT * random.Gaus()) : -1.;
      const bool topBefore = tofTop && tTop < tMid + 5.;
      const bool bottom = botFirst >= 0. && botFirst - tMid > -5. && botFirst - tMid < 50.;
      const bool overTarget = barMid >= 2 && barMid <= 6, upper = topLayer >= 3;
      const bool slow = eMidMeasured > eMidCut, tofSlow = !topBefore || dtTof > tofCut;
      if (alpha) {
        const bool pass[5] = {!bottom, overTarget, slow, upper, tofSlow};
        s.depth[1]->Fill(depth);
        for (int c = 0; c < 5; ++c) if (pass[c]) s.depth[2 + c]->Fill(depth);
        if (pass[0] && pass[1] && pass[2] && pass[3] && pass[4]) s.depth[7]->Fill(depth);
      }
      if (bottom || !overTarget || !upper) continue;
      int windows = 0;
      for (int w = 0; w < 3; ++w) if (InWindow(e, kWindows[w])) windows |= 1 << w;
      if (windows)
        s.entries.push_back({signal ? 0 : ekin > eSplit ? 2 : 1, windows, event & 1, topBefore, float(eMidMeasured),
                             float(dtTof), ekin, float(depth)});
    }
  }
  samples.erase(std::remove_if(samples.begin(), samples.end(), [](const Sample& s) { return s.days <= 0.; }), samples.end());
  if (samples.empty() || samples[0].cm != 1.) {
    printf("the 1 cm skim %s is needed (its faster muons give the background of every thickness)\n", ref);
    return;
  }
  const int nS = samples.size(), kThick = nS - 1;
  const Sample& thick = samples[kThick];

  // ---- 3. the scan: the measured middle-TOF energy above a cut, and the top-to-middle time above a cut where the top
  // TOF timed the muon (the first time cut is none)
  std::vector<double> eCuts, tCuts = {-1e9};
  for (int i = 0; i <= 12; ++i) eCuts.push_back((4 + i) / 4.);   // 1-4 MeV
  for (int j = 8; j <= 18; ++j) tCuts.push_back(j / 10.);        // 0.8-1.8 ns
  const int nE = eCuts.size(), nT = tCuts.size(), nCells = 3 * nE * nT;
  auto cell = [&](int w, int i, int j) { return (w * nE + i) * nT + j; };
  int iRef = -1, jRef = -1;   // the cuts of the 1 cm target
  for (int i = 0; i < nE; ++i) if (std::abs(eCuts[i] - eMidCut) < 1e-6) iRef = i;
  for (int j = 0; j < nT; ++j) if (std::abs(tCuts[j] - tofCut) < 1e-6) jRef = j;
  if (iRef < 0 || jRef < 0) {
    printf("the cuts %.2f MeV, %.2f ns are not on the grid of the scan\n", eMidCut, tofCut);
    return;
  }
  // counts per sample and half (row 2k + half) and cell: the target's X-rays, the background of the primaries up to
  // eSplit, its part from 0.5 GeV; and the fast-muon background of the 1 cm run, per half
  std::vector<std::vector<double>> sig(2 * nS, std::vector<double>(nCells, 0.)), slowB = sig, midB = sig;
  std::vector<std::vector<double>> fastB(2, std::vector<double>(nCells, 0.));
  for (int k = 0; k < nS; ++k)
    for (const auto& x : samples[k].entries) {
      if (x.kind == 2 && k > 0) continue;   // the thicker runs have no primaries above eSplit
      for (int w = 0; w < 3; ++w) {
        if (!(x.windows >> w & 1)) continue;
        for (int i = 0; i < nE && x.eMid > eCuts[i]; ++i)   // the cuts rise: below one, below the rest
          for (int j = 0; j < nT; ++j) {
            if (x.top && x.dt <= tCuts[j]) continue;
            const int c = cell(w, i, j), r = 2 * k + x.half;
            if (x.kind == 0) sig[r][c] += 1.;
            else if (x.kind == 2) fastB[x.half][c] += 1.;
            else {
              slowB[r][c] += 1.;
              if (x.ekin > 500.) midB[r][c] += 1.;
            }
          }
      }
    }
  // counts of sample k at a cell, from half h (0, 1) or both (2)
  auto count = [&](const std::vector<std::vector<double>>& v, int k, int c, int h) {
    return h == 2 ? v[2 * k][c] + v[2 * k + 1][c] : v[2 * k + h][c];
  };
  // S/sqrt(B) in 21 days and its statistical error; the fast-muon background of the 1 cm run times fastScale
  auto score = [&](int k, int c, int h, double* err = nullptr, double fastScale = 1.) {
    const double share = h == 2 ? 1. : 0.5;   // of the exposure
    const double dK = share * samples[k].days, d1 = share * samples[0].days;
    const double nSig = count(sig, k, c, h), nSlow = count(slowB, k, c, h);
    const double nFast = fastScale * (h == 2 ? fastB[0][c] + fastB[1][c] : fastB[h][c]);
    const double sDay = nSig / dK, b = nSlow / dK + nFast / d1;
    const double z = b > 0. ? std::sqrt(21.) * sDay / std::sqrt(b) : 0.;
    if (err) *err = z * std::sqrt(1. / std::max(1., nSig) + 0.25 * (nSlow / (dK * dK) + fastScale * nFast / (d1 * d1)) / (b * b));
    return z;
  };
  // the cell of the best L_alpha S/sqrt(B) of sample k on half h (or both), among those with 10 target X-rays or more
  auto best = [&](int k, int h) {
    int bestCell = cell(0, iRef, jRef);
    double bestZ = -1.;
    for (int i = 0; i < nE; ++i)
      for (int j = 0; j < nT; ++j) {
        const int c = cell(0, i, j);
        if (count(sig, k, c, h) < 10.) continue;
        const double z = score(k, c, h);
        if (z > bestZ) { bestZ = z; bestCell = c; }
      }
    return bestCell;
  };
  struct Result { int i, j; double z, dz, zCv, zRef, dzRef; };
  std::vector<Result> results(nS);
  for (int k = 0; k < nS; ++k) {
    auto& r = results[k];
    const int c = best(k, 2), c0 = best(k, 0), c1 = best(k, 1);
    r.i = c / nT; r.j = c % nT;
    r.z = score(k, c, 2, &r.dz);
    r.zCv = 0.5 * (score(k, c0, 1) + score(k, c1, 0));
    r.zRef = score(k, cell(0, iRef, jRef), 2, &r.dzRef);
  }
  auto tName = [&](int j) { return j == 0 ? std::string("none") : std::string(Form("> %.1f ns", tCuts[j])); };
  auto perDay = [&](const std::vector<std::vector<double>>& v, int k, int c) { return count(v, k, c, 2) / samples[k].days; };
  auto fastDay = [&](int c) { return (fastB[0][c] + fastB[1][c]) / samples[0].days; };

  // ---- tables
  printf("\n=== Thicker Al targets, the cuts taken apart (TOF %.2f ns, %.0f%%; cuts of the 1 cm target > %.1f MeV, > %.1f ns) ===\n",
         sigmaT, 100. * sigmaE, eMidCut, tofCut);
  for (const auto& s : samples)
    printf("  %g cm: %s, %.2f days, %zu entries in a window after the cuts not scanned\n", s.cm, s.file.Data(), s.days,
           s.entries.size());
  printf("\n1. The target's L_alpha entries of the %g cm target by the depth of the stop: the share with a middle-TOF hit in\n"
         "   coincidence, of those the share that passes each further cut on its own, and the share after all cuts\n", thick.cm);
  printf("   depth (cm)  entries  coincid.  no bottom  bars 2-6  > %.1f MeV  top layers  > %.1f ns  all cuts\n", eMidCut, tofCut);
  for (int b = 1; b <= thick.depth[0]->GetNbinsX(); ++b) {
    const double n = thick.depth[0]->GetBinContent(b), n1 = thick.depth[1]->GetBinContent(b);
    if (n < 1.) continue;
    printf("   %3.1f-%3.1f    %6.0f    %5.3f", thick.depth[0]->GetBinLowEdge(b), thick.depth[0]->GetBinLowEdge(b + 1), n, n1 / n);
    for (int k = 2; k <= 6; ++k) printf("     %5.3f  ", n1 > 0. ? thick.depth[k]->GetBinContent(b) / n1 : 0.);
    printf("  %5.3f\n", thick.depth[7]->GetBinContent(b) / n);
  }

  // 2. the slow-muon variables by depth (thickest target, its L_alpha entries after the cuts not scanned) and of the
  // fast-muon background (1 cm run, L_alpha window, the same cuts)
  const int nDepth = 10;
  std::vector<std::vector<double>> eByDepth(nDepth), tByDepth(nDepth);
  for (const auto& x : thick.entries) {
    if (x.kind != 0 || !(x.windows & 1)) continue;
    const int b = std::min(nDepth - 1, std::max(0, int(x.depth / 0.5)));
    eByDepth[b].push_back(x.eMid);
    if (x.top) tByDepth[b].push_back(x.dt);
  }
  std::vector<double> eFast, tFast;
  for (const auto& x : samples[0].entries) {
    if (x.kind != 2 || !(x.windows & 1)) continue;
    eFast.push_back(x.eMid);
    if (x.top) tFast.push_back(x.dt);
  }
  printf("\n2. The slow-muon variables (median, 16%%-84%%), the %g cm target's L_alpha entries after the cuts not scanned:\n", thick.cm);
  printf("   depth (cm)  entries  middle-TOF energy (MeV)   timed  top-to-middle time (ns)\n");
  for (int b = 0; b < nDepth; ++b) {
    if (eByDepth[b].size() < 5) continue;
    printf("   %3.1f-%3.1f    %5zu     %4.2f  (%4.2f-%4.2f)      %5zu     %4.2f  (%4.2f-%4.2f)\n", 0.5 * b, 0.5 * (b + 1),
           eByDepth[b].size(), Quantile(eByDepth[b], 0.5), Quantile(eByDepth[b], 0.16), Quantile(eByDepth[b], 0.84),
           tByDepth[b].size(), Quantile(tByDepth[b], 0.5), Quantile(tByDepth[b], 0.16), Quantile(tByDepth[b], 0.84));
  }
  printf("   background of the primaries above %.0f MeV in the L_alpha window (1 cm run): %4.2f (%4.2f-%4.2f) MeV, %zu entries;\n"
         "   timed %zu: %4.2f (%4.2f-%4.2f) ns\n", eSplit, Quantile(eFast, 0.5), Quantile(eFast, 0.16), Quantile(eFast, 0.84),
         eFast.size(), tFast.size(), Quantile(tFast, 0.5), Quantile(tFast, 0.16), Quantile(tFast, 0.84));

  printf("\n3. Each thickness with its own slow-muon cuts (scan: middle TOF 1-4 MeV by 0.25, top-to-middle none or 0.8-1.8 ns\n"
         "   by 0.1; best Al L_alpha S/sqrt(B)); cross-validated: best of the even events scored on the odd ones and back\n");
  printf(" target  cuts: middle TOF  top-middle   S/day   B slow   B fast   S/sqrt(B) in 21 days       days to   S vs\n");
  printf("   cm                                                         same events   cross-val.   5 sigma   1 cm\n");
  std::ofstream csv(Form("%s/muonic_xray_target_cuts.csv", dir));
  csv << "target_cm,days,cuts,emid_cut_MeV,tof_cut_ns,window,S_day,B_slow_day,B_fast_day,S_sqrtB_21d,S_sqrtB_err,"
         "S_sqrtB_cv,days_to_5sigma\n";
  const double sRef1 = perDay(sig, 0, cell(0, iRef, jRef)), sOwn1 = perDay(sig, 0, cell(0, results[0].i, results[0].j));
  for (int k = 0; k < nS; ++k) {
    const auto& r = results[k];
    for (int own = 0; own < 2; ++own) {
      const int i = own ? r.i : iRef, j = own ? r.j : jRef, c = cell(0, i, j);
      const double z = own ? r.z : r.zRef, dz = own ? r.dz : r.dzRef, zFive = own ? r.zCv : r.zRef;
      printf("  %4.1f   %-5s > %4.2f MeV  %-9s  %6.1f   %6.1f   %6.1f   %5.2f +- %4.2f   %s    %6.1f    %4.2f\n", samples[k].cm,
             own ? "own" : "1 cm", eCuts[i], tName(j).c_str(), perDay(sig, k, c), perDay(slowB, k, c), fastDay(c), z, dz,
             own ? Form("%5.2f", r.zCv) : "  -  ", zFive > 0. ? 21. * std::pow(5. / zFive, 2) : 0.,
             perDay(sig, k, c) / (own ? sOwn1 : sRef1));
      for (int w = 0; w < 3; ++w) {
        const int cw = cell(w, i, j);
        double dzw = 0.;
        const double zw = score(k, cw, 2, &dzw);
        csv << samples[k].cm << "," << samples[k].days << "," << (own ? "own" : "1cm") << "," << eCuts[i] << ","
            << (j == 0 ? std::string("none") : std::string(Form("%.1f", tCuts[j]))) << "," << kWindows[w].name << "," << perDay(sig, k, cw) << "," << perDay(slowB, k, cw)
            << "," << fastDay(cw) << "," << zw << "," << dzw << "," << (own && w == 0 ? r.zCv : -1.) << ","
            << (zw > 0. ? 21. * std::pow(5. / zw, 2) : 0.) << "\n";
      }
    }
  }
  printf("   (days to 5 sigma: from the cross-validated S/sqrt(B) for the own cuts)\n");
  printf("   the own cuts of L_alpha at the other lines, S/sqrt(B) in 21 days (1 cm cuts in brackets):\n");
  for (int k = 0; k < nS; ++k) {
    printf("  %4.1f cm:", samples[k].cm);
    for (int w = 1; w < 3; ++w)
      printf("   %s %5.2f (%5.2f)", kWindows[w].name, score(k, cell(w, results[k].i, results[k].j), 2),
             score(k, cell(w, iRef, jRef), 2));
    printf("\n");
  }
  printf("   L_alpha S/sqrt(B) in 21 days against the middle-TOF cut (best time cut at each):\n   cut (MeV)");
  for (int i = 0; i < nE; ++i) printf(" %5.2f", eCuts[i]);
  printf("\n");
  std::vector<std::vector<double>> curve(nS, std::vector<double>(nE, 0.));
  for (int k = 0; k < nS; ++k) {
    printf("   %4.1f cm  ", samples[k].cm);
    for (int i = 0; i < nE; ++i) {
      for (int j = 0; j < nT; ++j)
        if (count(sig, k, cell(0, i, j), 2) >= 10.) curve[k][i] = std::max(curve[k][i], score(k, cell(0, i, j), 2));
      printf(" %5.2f", curve[k][i]);
    }
    printf("\n");
  }

  // 4. the background of the primaries of 0.5 GeV to eSplit (they cross every target) per day, against the 1 cm target
  printf("\n4. Background of the muons that cross the target, L_alpha window, per day (relative to 1 cm):\n");
  printf(" target   primaries of 500-%.0f MeV                            all primaries up to %.0f MeV   S/sqrt(B), own cuts,\n",
         eSplit, eSplit);
  printf("   cm     cuts of the 1 cm target      own cuts                cuts of the 1 cm target       fast bkg scaled so\n");
  std::vector<double> midRef(nS), midRefErr(nS), midOwn(nS), midOwnErr(nS), slowRef(nS), slowRefErr(nS), zScaled(nS);
  for (int k = 0; k < nS; ++k) {
    const int cR = cell(0, iRef, jRef), cO = cell(0, results[k].i, results[k].j);
    auto rel = [&](const std::vector<std::vector<double>>& v, int c, double& err) {
      const double n = count(v, k, c, 2), n1 = count(v, 0, c, 2);
      const double r = n1 > 0. ? (n / samples[k].days) / (n1 / samples[0].days) : 0.;
      err = r * std::sqrt(1. / std::max(1., n) + 1. / std::max(1., n1));
      return r;
    };
    midRef[k] = rel(midB, cR, midRefErr[k]);
    midOwn[k] = rel(midB, cO, midOwnErr[k]);
    slowRef[k] = rel(slowB, cR, slowRefErr[k]);
    zScaled[k] = score(k, cO, 2, nullptr, midOwn[k]);
    printf("  %4.1f    %5.1f (%4.2f +- %4.2f)      %5.1f (%4.2f +- %4.2f)    %5.1f (%4.2f +- %4.2f)          %5.2f\n",
           samples[k].cm, perDay(midB, k, cR), midRef[k], midRefErr[k], perDay(midB, k, cO), midOwn[k], midOwnErr[k],
           perDay(slowB, k, cR), slowRef[k], slowRefErr[k], zScaled[k]);
  }

  // ---- plots
  gStyle->SetOptStat(0);
  TCanvas c("c", "thicker targets: the cuts", 1800, 1000);
  c.Divide(3, 2);
  const int colour[7] = {kGray + 2, kAzure + 2, kGreen + 2, kOrange + 1, kRed + 1, kMagenta + 1, kCyan + 2};   // as the compare
  TLatex heading;
  heading.SetNDC(); heading.SetTextSize(0.05);
  auto pad = [&](int p) {
    c.cd(p);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.1); gPad->SetBottomMargin(0.12);
  };
  auto frame = [](double x1, double y1, double x2, double y2, const char* titles) {
    auto f = gPad->DrawFrame(x1, y1, x2, y2, titles);
    f->GetXaxis()->SetTitleSize(0.045); f->GetYaxis()->SetTitleSize(0.045);
    f->GetXaxis()->SetTitleOffset(1.15); f->GetYaxis()->SetTitleOffset(1.45);
    return f;
  };
  auto legendAt = [](double x1, double y1, double x2, double y2, double size) {
    auto l = new TLegend(x1, y1, x2, y2);
    l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(size);
    return l;
  };
  auto styled = [](TGraph* g, int col, int marker, int line) {
    g->SetMarkerStyle(marker); g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineStyle(line); g->SetLineWidth(2);
    g->SetMarkerSize(1.1);
    return g;
  };
  auto hLine = [](double x1, double x2, double y, int col, int style) {
    auto l = new TLine(x1, y, x2, y);
    l->SetLineColor(col); l->SetLineStyle(style); l->SetLineWidth(2);
    l->Draw();
    return l;
  };

  // 1. each cut by depth
  pad(1);
  frame(0., 0., thick.cm, 1.5, ";depth of the stop below the top of the target (cm);share of the target's L-alpha entries");
  const std::vector<std::string> cutNames = {"middle-TOF hit (of all)", "no bottom-TOF hit", "middle bars 2-6",
                                             Form("middle TOF > %.1f MeV", eMidCut), "top two stack layers",
                                             Form("top-to-middle > %.1f ns", tofCut), "all cuts (of all)"};
  const int cutColours[7] = {kBlack, kGray + 1, kCyan + 2, kRed + 1, kGreen + 2, kAzure + 2, kOrange + 1};
  const int cutMarkers[7] = {24, 25, 26, 20, 32, 21, 29};
  auto legend1 = legendAt(0.16, 0.68, 0.97, 0.89, 0.034);
  legend1->SetNColumns(2);
  for (int q = 0; q < 7; ++q) {
    TH1D* pass = thick.depth[q == 6 ? 7 : q + 1];
    TH1D* all = thick.depth[q == 0 || q == 6 ? 0 : 1];
    auto g = new TGraphErrors();
    for (int b = 1; b <= all->GetNbinsX(); ++b) {
      const double n = all->GetBinContent(b), f = n > 0. ? pass->GetBinContent(b) / n : 0.;
      if (n < 10.) continue;
      const int p = g->GetN();
      g->SetPoint(p, all->GetBinCenter(b) + 0.03 * (q - 3), f);
      g->SetPointError(p, 0., std::sqrt(f * (1. - f) / n));
    }
    styled(g, cutColours[q], cutMarkers[q], 1);
    if (q == 6) g->SetLineWidth(3);
    g->Draw("PL SAME");
    legend1->AddEntry(g, cutNames[q].c_str(), "lp");
  }
  legend1->Draw();
  heading.DrawLatex(0.14, 0.93, Form("Which cut loses the deep stops (%g cm target)", thick.cm));

  // 2. and 3. the slow-muon variables by depth, against the fast-muon background
  auto variable = [&](int p, const std::vector<std::vector<double>>& byDepth, const std::vector<double>& fast, double yMax,
                      double refCut, double ownCut, const char* titles, const char* title, const char* unit) {
    pad(p);
    frame(0., 0., thick.cm, yMax, titles);
    auto band = new TBox(0., Quantile(fast, 0.16), thick.cm, Quantile(fast, 0.84));
    band->SetFillColor(kGray); band->SetFillStyle(1001);
    band->Draw();
    hLine(0., thick.cm, Quantile(fast, 0.5), kGray + 2, 1);
    auto g = new TGraphAsymmErrors();
    for (int b = 0; b < nDepth; ++b) {
      if (byDepth[b].size() < 10 || 0.5 * b >= thick.cm) continue;
      const double m = Quantile(byDepth[b], 0.5);
      const int q = g->GetN();
      g->SetPoint(q, 0.5 * b + 0.25, m);
      g->SetPointError(q, 0., 0., m - Quantile(byDepth[b], 0.16), Quantile(byDepth[b], 0.84) - m);
    }
    styled(g, kRed + 1, 20, 1);
    g->Draw("P SAME");
    auto lRef = hLine(0., thick.cm, refCut, kBlack, 2);
    auto legend = legendAt(0.36, 0.66, 0.97, 0.89, 0.036);
    legend->AddEntry(g, "target X-rays (L-alpha): median, 16-84%", "lp");
    legend->AddEntry(band, Form("background, muons > %.0f GeV: 16-84%%", eSplit / 1000.), "f");
    legend->AddEntry(lRef, Form("cut of the 1 cm target: %.1f %s", refCut, unit), "l");
    if (ownCut > 0.) {
      auto lOwn = hLine(0., thick.cm, ownCut, kRed + 1, 3);
      legend->AddEntry(lOwn, Form("own cut of the %g cm target: %.2f %s", thick.cm, ownCut, unit), "l");
    }
    legend->Draw();
    gPad->RedrawAxis();
    heading.DrawLatex(0.14, 0.93, title);
  };
  variable(2, eByDepth, eFast, 9., eMidCut, eCuts[results[kThick].i],
           ";depth of the stop below the top of the target (cm);measured middle-TOF energy (MeV)",
           "Deeper stops come from faster muons", "MeV");
  variable(3, tByDepth, tFast, 3.6, tofCut, results[kThick].j > 0 ? tCuts[results[kThick].j] : -1.,
           ";depth of the stop below the top of the target (cm);measured top-to-middle time (ns)",
           "... and a shorter time of flight", "ns");

  // 4. S/sqrt(B) against the middle-TOF cut, the best time cut at each
  pad(4);
  double zMax = 0.;
  for (int k = 0; k < nS; ++k)
    zMax = std::max({zMax, results[k].z + results[k].dz, zScaled[k], *std::max_element(curve[k].begin(), curve[k].end())});
  frame(eCuts.front() - 0.1, 0., eCuts.back() + 0.1, 1.45 * zMax, ";cut: measured middle-TOF energy above (MeV);S/sqrt(B) at Al L-alpha in 21 days");
  auto legend4 = legendAt(0.16, 0.7, 0.97, 0.89, 0.036);
  legend4->SetNColumns(3);
  for (int k = 0; k < nS; ++k) {
    auto g = new TGraph();
    for (int i = 0; i < nE; ++i) if (curve[k][i] > 0.) g->SetPoint(g->GetN(), eCuts[i], curve[k][i]);
    styled(g, colour[k % 7], 20, 1);
    g->SetMarkerSize(0.8);
    g->Draw("PL SAME");
    auto chosen = new TGraph(1);
    chosen->SetPoint(0, eCuts[results[k].i], results[k].z);
    styled(chosen, colour[k % 7], 24, 1);
    chosen->SetMarkerSize(2.2);
    chosen->Draw("P SAME");
    legend4->AddEntry(g, Form("%g cm target", samples[k].cm), "lp");
  }
  auto vRef = new TLine(eMidCut, 0., eMidCut, 1.1 * zMax);
  vRef->SetLineStyle(2); vRef->SetLineWidth(2); vRef->Draw();
  legend4->AddEntry(vRef, "cut of the 1 cm target", "l");
  legend4->Draw();
  heading.DrawLatex(0.14, 0.93, "Each thickness has its own best cut");

  // 5. S/sqrt(B) against the thickness
  pad(5);
  const double xMax = samples.back().cm + 0.5;
  frame(0.5, 0., xMax, 1.45 * zMax, ";Al target thickness (cm);S/sqrt(B) at Al L-alpha in 21 days");
  auto gRef = new TGraphErrors(), gOwn = new TGraphErrors(), gCv = new TGraphErrors(), gScaled = new TGraphErrors();
  for (int k = 0; k < nS; ++k) {
    const auto& r = results[k];
    gRef->SetPoint(k, samples[k].cm - 0.06, r.zRef); gRef->SetPointError(k, 0., r.dzRef);
    gOwn->SetPoint(k, samples[k].cm, r.z); gOwn->SetPointError(k, 0., r.dz);
    gCv->SetPoint(k, samples[k].cm + 0.06, r.zCv); gCv->SetPointError(k, 0., r.dz);
    gScaled->SetPoint(k, samples[k].cm + 0.12, zScaled[k]); gScaled->SetPointError(k, 0., 0.);
  }
  styled(gRef, kGray + 2, 24, 2)->Draw("PL SAME");
  styled(gOwn, kRed + 1, 24, 3)->Draw("PL SAME");
  styled(gCv, kRed + 1, 20, 1)->Draw("PL SAME");
  styled(gScaled, kAzure + 2, 25, 1)->Draw("P SAME");
  auto legend5 = legendAt(0.16, 0.68, 0.97, 0.89, 0.036);
  legend5->AddEntry(gRef, "cuts of the 1 cm target", "lp");
  legend5->AddEntry(gOwn, "own cuts, on the events that chose them", "lp");
  legend5->AddEntry(gCv, "own cuts, cross-validated", "lp");
  legend5->AddEntry(gScaled, "own cuts, fast background as in panel 6", "p");
  legend5->Draw();
  heading.DrawLatex(0.14, 0.93, "Al L-alpha with each thickness's own cuts");

  // 6. the background of the muons that cross the target
  pad(6);
  double rMax = 0.;
  for (int k = 0; k < nS; ++k) rMax = std::max({rMax, midRef[k] + midRefErr[k], midOwn[k] + midOwnErr[k], slowRef[k] + slowRefErr[k]});
  frame(0.5, 0., xMax, 1.45 * rMax, ";Al target thickness (cm);background per day, relative to 1 cm");
  hLine(0.5, xMax, 1., kGray + 1, 2);
  auto gMidRef = new TGraphErrors(), gMidOwn = new TGraphErrors(), gSlow = new TGraphErrors();
  for (int k = 0; k < nS; ++k) {
    gMidRef->SetPoint(k, samples[k].cm - 0.06, midRef[k]); gMidRef->SetPointError(k, 0., midRefErr[k]);
    gMidOwn->SetPoint(k, samples[k].cm, midOwn[k]); gMidOwn->SetPointError(k, 0., midOwnErr[k]);
    gSlow->SetPoint(k, samples[k].cm + 0.06, slowRef[k]); gSlow->SetPointError(k, 0., slowRefErr[k]);
  }
  styled(gMidRef, kAzure + 2, 20, 1)->Draw("PL SAME");
  styled(gMidOwn, kAzure + 2, 24, 2)->Draw("PL SAME");
  styled(gSlow, kGray + 2, 25, 1)->Draw("PL SAME");
  auto legend6 = legendAt(0.16, 0.7, 0.97, 0.89, 0.034);
  legend6->AddEntry(gMidRef, Form("muons of 0.5-%.0f GeV, cuts of the 1 cm target", eSplit / 1000.), "lp");
  legend6->AddEntry(gMidOwn, Form("muons of 0.5-%.0f GeV, own cuts", eSplit / 1000.), "lp");
  legend6->AddEntry(gSlow, Form("all muons up to %.0f GeV, cuts of the 1 cm target", eSplit / 1000.), "lp");
  legend6->Draw();
  heading.DrawLatex(0.14, 0.93, "Background of the muons that cross the target");

  c.Print(Form("%s/muonic_xray_target_cuts.pdf", dir));
  c.SaveAs(Form("%s/muonic_xray_target_cuts.png", dir));
  printf("\nplots: %s/muonic_xray_target_cuts.png (and .pdf); table: %s/muonic_xray_target_cuts.csv\n", dir, dir);
}
