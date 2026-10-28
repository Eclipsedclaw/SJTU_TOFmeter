// Thicker Al targets: the X-ray selection of analysis/muonic_xray_stack_cuts.C on the 1 cm target of the 3-week run and
// on the runs of scripts/muonic_xray_target.sh (thicker targets, only the cosmic muons below 1 GeV), from their skims
// (analysis/muonic_xray_stack_skim.C, xray_skim.root in each folder):
//   root -l -b -q 'analysis/muonic_xray_target_compare.C'
// The target runs hold the signal and the backgrounds of the slow muons; the background of the faster muons (primaries
// above eSplit MeV) is that of the 3-week run, taken as the same for every thickness. Per thickness: mu- stops in the
// target with a prompt stack entry (truth) and the share whose muon came round the top lead layer (through a wall);
// entries with an X-ray of the target in +-2 sigma around Al L_alpha, L_beta and K_alpha before and after the cuts; the
// backgrounds after the cuts; S/sqrt(B) for 21 days and the days to 5 sigma; the efficiency of the cuts for the L_alpha
// entries versus the depth of the stop in the target. With the 1 cm run alone (no thicker skims yet) it gives the cut
// flow's last stage (S 412, B 3110 at L_alpha). Writes <dir>/muonic_xray_target_compare.png/.pdf/.csv.

#include <fstream>
#include <set>
#include <sstream>

namespace {
// as analysis/muonic_xray_stack_cuts.C
double SigmaKeV(double eKeV) { return 40. * std::sqrt(eKeV / 1000.); }

struct Window { const char* name; double eKeV; };
const std::vector<Window> kWindows = {{"Al L_alpha", 65.73}, {"Al L_beta", 88.73}, {"Al K_alpha", 354.9}};

bool InWindow(double eKeV, const Window& w) { return std::abs(eKeV - w.eKeV) <= 2. * SigmaKeV(w.eKeV); }

struct Sample {
  double cm = 1., top = 36.5, days = 0.;   // thickness, top of the target (cm), exposure
  TString file;
  bool hasPb = false;
  std::set<long long> stops, wallStops;   // mu- stopped in the target with a prompt entry; of those, round the top lead
  double sAll[3] = {}, sCut[3] = {}, bLow[3] = {}, bHigh[3] = {}, plates = 0., carbon = 0.;   // counts
  TH1D *depthAll = nullptr, *depthCut = nullptr, *signal = nullptr;
};
}  // namespace

void muonic_xray_target_compare(const char* ref = "output/muonic_xray_stack/xray_skim.root",
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
    int proc, event, topLayer, src, endVol, endProc, tofTop, tofMid, barMid;
    float e, t, tTop, tMid, eMid, botFirst, ekin, endT, endZ, pbKe = 0.f;
    cand->SetBranchAddress("proc", &proc); cand->SetBranchAddress("event", &event);
    cand->SetBranchAddress("e", &e); cand->SetBranchAddress("topLayer", &topLayer); cand->SetBranchAddress("t", &t);
    cand->SetBranchAddress("src", &src); cand->SetBranchAddress("ekin", &ekin);
    cand->SetBranchAddress("endVol", &endVol); cand->SetBranchAddress("endProc", &endProc);
    cand->SetBranchAddress("endT", &endT); cand->SetBranchAddress("endZ", &endZ);
    cand->SetBranchAddress("tofTop", &tofTop); cand->SetBranchAddress("tofMid", &tofMid);
    cand->SetBranchAddress("tTop", &tTop); cand->SetBranchAddress("tMid", &tMid);
    cand->SetBranchAddress("barMid", &barMid); cand->SetBranchAddress("eMid", &eMid); cand->SetBranchAddress("botFirst", &botFirst);
    s.hasPb = cand->GetBranch("pbKe") != nullptr;   // skims made before it was added have no pbKe
    if (s.hasPb) cand->SetBranchAddress("pbKe", &pbKe);
    s.depthAll = new TH1D(Form("depthAll%g", s.cm), "", 10, 0., 5.);
    s.depthCut = new TH1D(Form("depthCut%g", s.cm), "", 10, 0., 5.);
    s.signal = new TH1D(Form("signal%g", s.cm), "", 65, 30., 160.);
    TRandom3 random;
    for (Long64_t i = 0; i < cand->GetEntries(); ++i) {
      cand->GetEntry(i);
      const bool signal = src & 1, inAlpha = InWindow(e, kWindows[0]);
      if (endVol == 0 && endProc == 1 && t - endT > -5. && t - endT < 50.) {
        const long long id = 1000000000LL * proc + event;
        s.stops.insert(id);
        if (s.hasPb && pbKe < 0.) s.wallStops.insert(id);
      }
      if (signal) {
        for (int w = 0; w < 3; ++w)
          if (InWindow(e, kWindows[w])) s.sAll[w] += 1.;
        if (inAlpha) s.depthAll->Fill(s.top - endZ);
      }
      // the selection of muonic_xray_stack_cuts.C, with its TOF response and seeds
      if (!tofMid || t - tMid < -5. || t - tMid > 30.) continue;
      random.SetSeed(1000003ULL * proc + event + 1);
      const double eMidMeasured = eMid * (1. + sigmaE * random.Gaus());
      const double dtTof = tofTop ? (tMid + sigmaT * random.Gaus()) - (tTop + sigmaT * random.Gaus()) : -1.;
      const bool topBefore = tofTop && tTop < tMid + 5.;
      const bool bottom = botFirst >= 0. && botFirst - tMid > -5. && botFirst - tMid < 50.;
      if (bottom || barMid < 2 || barMid > 6 || eMidMeasured <= eMidCut || topLayer < 3 || (topBefore && dtTof <= tofCut))
        continue;
      if (signal) {
        if (inAlpha) s.depthCut->Fill(s.top - endZ);
        s.signal->Fill(e);
      }
      for (int w = 0; w < 3; ++w) {
        if (!InWindow(e, kWindows[w])) continue;
        if (signal) s.sCut[w] += 1.;
        else if (ekin > eSplit) s.bHigh[w] += 1.;
        else s.bLow[w] += 1.;
      }
      if (!signal && inAlpha && (src & 2) && ekin <= eSplit) {
        if (endVol == 1) s.plates += 1.;
        if (endVol == 4) s.carbon += 1.;
      }
    }
  }
  std::vector<Sample> good;
  for (auto& s : samples) if (s.days > 0.) good.push_back(s);
  if (good.empty() || good[0].cm != 1.) {
    printf("the 1 cm skim %s is needed (its faster muons give the background of every thickness)\n", ref);
    return;
  }
  const Sample& one = good[0];

  // ---- tables
  printf("\n=== Thicker Al targets: the selection of muonic_xray_stack_cuts.C (TOF %.2f ns, %.0f%%; > %.1f MeV; > %.1f ns) ===\n",
         sigmaT, 100. * sigmaE, eMidCut, tofCut);
  printf("1 cm: %s (all energies); thicker: %s/t<d>cm (muons below %.0f MeV); the background of the muons above is that of\n"
         "the 1 cm run for every thickness. Counts per day; S/sqrt(B) and days to 5 sigma for the +-2 sigma windows.\n\n",
         ref, dir, eSplit);
  printf(" target   days   mu- stops per day    through   muonic X-rays after the cuts\n");
  printf("   cm            (prompt entry)       a wall    of Al plates   of carbon (per day, L_alpha window)\n");
  for (const auto& s : good)
    printf("  %4.1f  %6.2f   %7.1f +- %4.1f     %s     %6.2f         %6.2f\n", s.cm, s.days, s.stops.size() / s.days,
           std::sqrt(double(s.stops.size())) / s.days,
           s.hasPb ? Form("%5.1f%%", 100. * s.wallStops.size() / std::max<double>(1., s.stops.size())) : "   - ",
           s.plates / s.days, s.carbon / s.days);
  gSystem->mkdir(dir, true);
  std::ofstream csv(Form("%s/muonic_xray_target_compare.csv", dir));
  csv << "target_cm,days,mu_stops_day,wall_share,window,S_before_day,S_day,S_err,efficiency,B_slow_day,B_fast_day,"
         "S_sqrtB_21d,S_sqrtB_err,days_to_5sigma\n";
  double zMax = 0.;
  TGraphErrors *gain[4], *signif[3];
  for (auto& g : gain) g = new TGraphErrors();
  for (auto& g : signif) g = new TGraphErrors();
  for (int w = 0; w < 3; ++w) {
    printf("\n%s (%.1f keV, window +-%.1f keV):\n", kWindows[w].name, kWindows[w].eKeV, 2. * SigmaKeV(kWindows[w].eKeV));
    printf(" target   S before   S after cuts    cut     B slow   B fast   S/sqrt(B)       days to   S vs 1 cm\n");
    printf("   cm     the cuts                   eff.   muons    muons    in 21 days      5 sigma\n");
    const double bFast = one.bHigh[w] / one.days, dbFast2 = one.bHigh[w] / (one.days * one.days);
    for (const auto& s : good) {
      const double sDay = s.sCut[w] / s.days, bSlow = s.bLow[w] / s.days, b = bSlow + bFast;
      const double z = b > 0. ? std::sqrt(21.) * sDay / std::sqrt(b) : 0.;
      const double rel = std::sqrt(1. / std::max(1., s.sCut[w]) + 0.25 * (s.bLow[w] / (s.days * s.days) + dbFast2) / (b * b));
      const double toFive = z > 0. ? 21. * std::pow(5. / z, 2) : 0.;
      const double ratio = sDay / (one.sCut[w] / one.days);
      const double dRatio = ratio * std::sqrt(1. / std::max(1., s.sCut[w]) + 1. / std::max(1., one.sCut[w]));
      printf("  %4.1f   %7.1f   %6.1f +- %4.1f   %5.3f   %6.1f   %6.1f   %5.2f +- %4.2f   %6.1f     %4.2f +- %4.2f\n", s.cm,
             s.sAll[w] / s.days, sDay, std::sqrt(s.sCut[w]) / s.days, s.sCut[w] / std::max(1., s.sAll[w]), bSlow, bFast, z,
             z * rel, toFive, ratio, dRatio);
      csv << s.cm << "," << s.days << "," << s.stops.size() / s.days << ","
          << (s.hasPb ? double(s.wallStops.size()) / std::max<double>(1., s.stops.size()) : -1.) << "," << kWindows[w].name << ","
          << s.sAll[w] / s.days << "," << sDay << "," << std::sqrt(s.sCut[w]) / s.days << ","
          << s.sCut[w] / std::max(1., s.sAll[w]) << "," << bSlow << "," << bFast << "," << z << "," << z * rel << "," << toFive << "\n";
      int p = gain[w]->GetN();
      gain[w]->SetPoint(p, s.cm, ratio);
      gain[w]->SetPointError(p, 0., dRatio);
      p = signif[w]->GetN();
      signif[w]->SetPoint(p, s.cm, z);
      signif[w]->SetPointError(p, 0., z * rel);
      zMax = std::max(zMax, z * (1. + rel));
    }
  }
  for (const auto& s : good) {
    const double ratio = (s.stops.size() / s.days) / (one.stops.size() / one.days);
    const int p = gain[3]->GetN();
    gain[3]->SetPoint(p, s.cm, ratio);
    gain[3]->SetPointError(p, 0., ratio * std::sqrt(1. / std::max<double>(1., s.stops.size()) + 1. / one.stops.size()));
  }

  // ---- plots
  gStyle->SetOptStat(0);
  TCanvas c("c", "thicker targets", 1400, 1000);
  c.Divide(2, 2);
  const int colours[3] = {kRed + 1, kOrange + 1, kAzure + 2};
  const char* lineNames[3] = {"Al L_{#alpha}", "Al L_{#beta}", "Al K_{#alpha}"};
  TLatex heading;
  heading.SetNDC();
  heading.SetTextSize(0.045);
  const double xMax = good.back().cm + 0.5;
  c.cd(1);
  gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.1);
  double gMax = 0.;
  for (auto g : gain) for (int p = 0; p < g->GetN(); ++p) gMax = std::max(gMax, g->GetY()[p] + g->GetEY()[p]);
  gPad->DrawFrame(0.5, 0., xMax, 1.2 * gMax, ";Al target thickness (cm);after the cuts, relative to 1 cm");
  auto legend = new TLegend(0.15, 0.68, 0.55, 0.88);
  legend->SetBorderSize(0); legend->SetFillStyle(0);
  gain[3]->SetMarkerStyle(24); gain[3]->SetMarkerColor(kGray + 2); gain[3]->SetLineColor(kGray + 2);
  gain[3]->Draw("PL SAME");
  legend->AddEntry(gain[3], "#mu^{-} stops with a prompt entry (no cuts)", "lp");
  for (int w = 0; w < 3; ++w) {
    gain[w]->SetMarkerStyle(20); gain[w]->SetMarkerColor(colours[w]); gain[w]->SetLineColor(colours[w]); gain[w]->SetLineWidth(2);
    gain[w]->Draw("PL SAME");
    legend->AddEntry(gain[w], Form("%s entries of the target", lineNames[w]), "lp");
  }
  legend->Draw();
  heading.DrawLatex(0.12, 0.93, "Signal after the cuts versus target thickness");
  c.cd(2);
  gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.1);
  gPad->DrawFrame(0.5, 0., xMax, 1.25 * zMax, ";Al target thickness (cm);S / #sqrt{B} in 21 days");
  auto legend2 = new TLegend(0.15, 0.72, 0.5, 0.88);
  legend2->SetBorderSize(0); legend2->SetFillStyle(0);
  for (int w = 0; w < 3; ++w) {
    signif[w]->SetMarkerStyle(20); signif[w]->SetMarkerColor(colours[w]); signif[w]->SetLineColor(colours[w]);
    signif[w]->SetLineWidth(2);
    signif[w]->Draw("PL SAME");
    legend2->AddEntry(signif[w], lineNames[w], "lp");
  }
  legend2->Draw();
  heading.DrawLatex(0.12, 0.93, "S/#sqrt{B} (fast-muon background of the 1 cm run)");
  const int sampleColours[] = {kGray + 2, kAzure + 2, kGreen + 2, kOrange + 1, kRed + 1, kMagenta + 1, kCyan + 2};
  c.cd(3);
  gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.1);
  gPad->DrawFrame(0., 0., 5., 1.05, ";depth of the stop below the top of the target (cm);efficiency of the cuts (L_{#alpha})");
  auto legend3 = new TLegend(0.65, 0.66, 0.95, 0.88);
  legend3->SetBorderSize(0); legend3->SetFillStyle(0);
  for (std::size_t k = 0; k < good.size(); ++k) {
    auto g = new TGraphErrors();
    for (int b = 1; b <= good[k].depthAll->GetNbinsX(); ++b) {
      const double n = good[k].depthAll->GetBinContent(b), m = good[k].depthCut->GetBinContent(b);
      if (n < 5.) continue;
      const double f = m / n;
      const int p = g->GetN();
      g->SetPoint(p, good[k].depthAll->GetBinCenter(b) + 0.04 * (k - 0.5 * (good.size() - 1.)), f);
      g->SetPointError(p, 0., std::sqrt(f * (1. - f) / n));
    }
    if (g->GetN() == 0) continue;
    g->SetMarkerStyle(20); g->SetMarkerColor(sampleColours[k % 7]); g->SetLineColor(sampleColours[k % 7]);
    g->Draw("P SAME");
    legend3->AddEntry(g, Form("%g cm target", good[k].cm), "lp");
  }
  legend3->Draw();
  heading.DrawLatex(0.12, 0.93, "Cut efficiency for the target's L_{#alpha} entries by stop depth");
  c.cd(4);
  gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.1);
  double sMax = 0.;
  for (auto& s : good) {
    s.signal->Scale(1. / s.days);
    sMax = std::max(sMax, s.signal->GetMaximum());
  }
  gPad->DrawFrame(30., 0., 160., 1.25 * sMax, ";add-back energy (keV);target X-ray entries after the cuts per day per 2 keV");
  auto legend4 = new TLegend(0.65, 0.66, 0.95, 0.88);
  legend4->SetBorderSize(0); legend4->SetFillStyle(0);
  for (std::size_t k = 0; k < good.size(); ++k) {
    good[k].signal->SetLineColor(sampleColours[k % 7]);
    good[k].signal->SetLineWidth(2);
    good[k].signal->Draw("HIST SAME");
    legend4->AddEntry(good[k].signal, Form("%g cm target", good[k].cm), "l");
  }
  legend4->Draw();
  heading.DrawLatex(0.12, 0.93, "The target's X-rays after the cuts");
  c.Print(Form("%s/muonic_xray_target_compare.pdf", dir));
  c.SaveAs(Form("%s/muonic_xray_target_compare.png", dir));
  printf("\nplots: %s/muonic_xray_target_compare.png (and .pdf); table: %s/muonic_xray_target_compare.csv\n", dir, dir);
}
