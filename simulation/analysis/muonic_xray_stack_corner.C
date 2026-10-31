// Corner plot of the cut variables of analysis/muonic_xray_stack_cuts.C: how the target's X-rays (S) and everything
// else (B) in the Al L_alpha window spread over the five variables the cuts use, each against each, for the entries
// with a middle-TOF hit in coincidence (the first stage of the cut flow), with the cuts drawn in. From a skim
// (analysis/muonic_xray_stack_skim.C):
//   root -l -b -q 'analysis/muonic_xray_stack_corner.C'                       the 1 cm target, with its cuts
//   root -l -b -q 'analysis/muonic_xray_stack_corner.C("output/muonic_xray_target/t4cm/xray_skim.root", 1.75, 1.1)'
// A thicker target's skim (scripts/muonic_xray_target.sh) holds only the muons up to eSplit: its background is
// completed with the faster muons of the 3-week skim ref, each part per day of its own exposure. The variables in the
// order of the cut flow, and what the cuts keep: a bottom-TOF hit within -5..+50 ns of the middle one (none), the
// middle-TOF bar (2-6, over the target), the measured middle-TOF energy (above eMidCut), the highest stack layer hit
// (3-4, the top YSO and Si) and the measured top-to-middle time (above tofCut where the top TOF timed the muon; the
// others pass). Diagonal: the shapes of S and B (unit area), what the cut rejects (hatched), and the shares of S and B
// it keeps on its own and after all the other cuts (N-1). Below it: B as a density (grey, per day, log scale) and S as
// points (red; spread a little on the integer variables). TOF response as the cuts macro (sigmaT per station, sigmaE,
// the same seeds). S counts every entry with a target X-ray (stk_src bit 1), as the cut studies do.
// Writes <skim dir>/muonic_xray_corner.png/.pdf.

namespace {
double SigmaKeV(double eKeV) { return 40. * std::sqrt(eKeV / 1000.); }
const double kAlphaKeV = 65.73;                // Al L_alpha, the window of muonic_xray_stack_cuts.C
const double kCutE1cm = 3.0, kCutT1cm = 1.3;   // the slow-muon cuts of the 1 cm target

// an entry of the L_alpha window with a middle-TOF hit in coincidence: weight (per day), target X-ray or not, whether
// the top TOF timed the muon, and the variables: bottom-TOF hit (0, 1), middle-TOF bar, measured middle-TOF energy
// (MeV), highest stack layer, measured top-to-middle time (ns)
struct Entry { double w; bool signal, timed; double v[5]; };
}  // namespace

void muonic_xray_stack_corner(const char* skim = "output/muonic_xray_stack/xray_skim_pbke.root", double eMidCut = 3.0,
                              double tofCut = 1.3, const char* ref = "output/muonic_xray_stack/xray_skim_pbke.root",
                              double eSplit = 1000., double sigmaT = 0.2, double sigmaE = 0.15)
{
  TH1::AddDirectory(false);
  std::vector<Entry> entries;
  // the entries of a skim; fastOnly: only the background of its primaries above eSplit
  auto load = [&](const char* path, bool fastOnly) {
    TFile file(path);
    auto cand = dynamic_cast<TTree*>(file.Get("cand"));
    auto info = dynamic_cast<TTree*>(file.Get("info"));
    if (!cand || !info) {
      printf("no skim in %s (run analysis/muonic_xray_stack_skim.C)\n", path);
      return 0.;
    }
    double lt = 0., days = 0.;
    info->SetBranchAddress("live_time", &lt);
    for (Long64_t i = 0; i < info->GetEntries(); ++i) { info->GetEntry(i); days += lt / 86400.; }
    int proc, event, topLayer, src, tofTop, tofMid, barMid;
    float e, t, tTop, tMid, eMid, botFirst, ekin;
    cand->SetBranchAddress("proc", &proc); cand->SetBranchAddress("event", &event);
    cand->SetBranchAddress("e", &e); cand->SetBranchAddress("topLayer", &topLayer); cand->SetBranchAddress("t", &t);
    cand->SetBranchAddress("src", &src); cand->SetBranchAddress("ekin", &ekin);
    cand->SetBranchAddress("tofTop", &tofTop); cand->SetBranchAddress("tofMid", &tofMid);
    cand->SetBranchAddress("tTop", &tTop); cand->SetBranchAddress("tMid", &tMid);
    cand->SetBranchAddress("barMid", &barMid); cand->SetBranchAddress("eMid", &eMid); cand->SetBranchAddress("botFirst", &botFirst);
    TRandom3 random;
    for (Long64_t i = 0; i < cand->GetEntries(); ++i) {
      cand->GetEntry(i);
      if (std::abs(e - kAlphaKeV) > 2. * SigmaKeV(kAlphaKeV)) continue;   // the L_alpha window
      if (!tofMid || t - tMid < -5. || t - tMid > 30.) continue;           // a middle-TOF hit in coincidence
      const bool signal = src & 1;
      if (fastOnly && (signal || ekin <= eSplit)) continue;
      // the TOF response of muonic_xray_stack_cuts.C, with its seeds
      random.SetSeed(1000003ULL * proc + event + 1);
      const double eMidMeasured = eMid * (1. + sigmaE * random.Gaus());
      const double dtTof = tofTop ? (tMid + sigmaT * random.Gaus()) - (tTop + sigmaT * random.Gaus()) : -1.;
      const bool timed = tofTop && tTop < tMid + 5.;
      const bool bottom = botFirst >= 0. && botFirst - tMid > -5. && botFirst - tMid < 50.;
      entries.push_back({1. / days, signal, timed, {bottom ? 1. : 0., double(barMid), eMidMeasured, double(topLayer), dtTof}});
    }
    return days;
  };
  const bool thick = std::string(skim) != ref;
  const double days = load(skim, false), refDays = thick ? load(ref, true) : days;
  if (days <= 0. || refDays <= 0.) return;
  TString label = "1 cm";
  TObjArray* match = TPRegexp("/t([0-9.]+)cm/").MatchS(skim);
  if (match->GetEntries() > 1) label = static_cast<TObjString*>(match->At(1))->GetString() + " cm";
  delete match;

  // ---- the cuts, alone and after all the others (N-1)
  const int nV = 5;
  auto passes = [&](const Entry& x, int i) {
    switch (i) {
      case 0: return x.v[0] < 0.5;                   // no bottom-TOF hit
      case 1: return x.v[1] >= 2. && x.v[1] <= 6.;   // a middle bar over the target
      case 2: return x.v[2] > eMidCut;               // a slow muon
      case 3: return x.v[3] >= 3.;                   // the highest hit in the top two layers
      default: return !x.timed || x.v[4] > tofCut;   // a slow muon, where the top TOF timed it
    }
  };
  double sAll = 0., bAll = 0., sPass = 0., bPass = 0.;
  double sAlone[nV] = {}, bAlone[nV] = {}, sOthers[nV] = {}, bOthers[nV] = {}, sNm1[nV] = {}, bNm1[nV] = {};
  for (const auto& x : entries) {
    bool pass[nV];
    int failed = 0, which = -1;
    for (int i = 0; i < nV; ++i) {
      pass[i] = passes(x, i);
      if (!pass[i]) { ++failed; which = i; }
    }
    (x.signal ? sAll : bAll) += x.w;
    if (failed == 0) (x.signal ? sPass : bPass) += x.w;
    for (int i = 0; i < nV; ++i) {
      if (pass[i]) (x.signal ? sAlone : bAlone)[i] += x.w;
      if (failed == 0 || (failed == 1 && which == i)) {   // passes every other cut
        (x.signal ? sOthers : bOthers)[i] += x.w;
        if (pass[i]) (x.signal ? sNm1 : bNm1)[i] += x.w;
      }
    }
  }
  const char* cutNames[nV] = {"no bottom-TOF hit within 50 ns", "middle-TOF bar 2-6", Form("middle-TOF energy > %.2f MeV", eMidCut),
                              "highest hit in the top two layers", Form("top-to-middle time > %.2f ns, if timed", tofCut)};
  printf("\n=== Cut variables in the Al L_alpha window, %s target: entries with a middle-TOF hit in coincidence ===\n", label.Data());
  printf("%s, %.2f days%s\n", skim, days, thick ? Form("; faster muons of %s, %.2f days", ref, refDays) : "");
  printf("per day before these cuts: S %.2f, B %.1f; after all of them: S %.2f, B %.1f, S/sqrt(B) in 21 days %.2f\n", sAll, bAll,
         sPass, bPass, bPass > 0. ? std::sqrt(21.) * sPass / std::sqrt(bPass) : 0.);
  printf("   cut                                         alone: S kept   B kept     after the others: S kept   B kept\n");
  for (int i = 0; i < nV; ++i)
    printf("   %-42s     %5.1f%%   %6.1f%%                      %5.1f%%   %6.1f%%\n", cutNames[i], 100. * sAlone[i] / sAll,
           100. * bAlone[i] / bAll, 100. * sNm1[i] / std::max(1e-12, sOthers[i]), 100. * bNm1[i] / std::max(1e-12, bOthers[i]));
  // the two slow-muon variables measure the same speed: their correlation where both are measured
  for (bool sig : {true, false}) {
    double n = 0., mx = 0., my = 0., xx = 0., yy = 0., xy = 0.;
    for (const auto& x : entries) {
      if (x.signal != sig || !x.timed) continue;
      n += x.w; mx += x.w * x.v[2]; my += x.w * x.v[4];
      xx += x.w * x.v[2] * x.v[2]; yy += x.w * x.v[4] * x.v[4]; xy += x.w * x.v[2] * x.v[4];
    }
    mx /= n; my /= n;
    const double r = (xy / n - mx * my) / std::sqrt((xx / n - mx * mx) * (yy / n - my * my));
    printf("   correlation of middle-TOF energy and top-to-middle time (timed entries), %s: %.2f\n", sig ? "S" : "B", r);
  }

  // ---- histograms: diagonal shapes, B densities and S points of each pair
  const char* names[nV] = {"bottom-TOF hit", "middle-TOF bar", "middle-TOF energy (MeV)", "highest stack layer",
                           "top-to-middle time (ns)"};
  const int nBins[nV] = {2, 9, 40, 5, 40};
  const double lo[nV] = {-0.5, -0.5, 0., -0.5, 0.}, hi[nV] = {1.5, 8.5, 10., 4.5, 4.};
  const bool integer[nV] = {true, true, false, true, false};
  auto defined = [](const Entry& x, int i) { return i != 4 || x.timed; };   // the time only where the top TOF timed the muon
  auto inRange = [&](double v, int i) { return std::min(std::max(v, lo[i] + 1e-6), hi[i] - 1e-6); };   // overflow: edge bins
  TH1D *hS[nV], *hB[nV];
  for (int i = 0; i < nV; ++i) {
    hS[i] = new TH1D(Form("s%d", i), "", nBins[i], lo[i], hi[i]);
    hB[i] = new TH1D(Form("b%d", i), "", nBins[i], lo[i], hi[i]);
  }
  TH2D* hB2[nV][nV] = {};
  TGraph* gS[nV][nV] = {};
  for (int i = 1; i < nV; ++i)
    for (int j = 0; j < i; ++j) {
      hB2[i][j] = new TH2D(Form("b%d%d", i, j), "", nBins[j], lo[j], hi[j], nBins[i], lo[i], hi[i]);
      gS[i][j] = new TGraph();
    }
  TRandom3 spread(7);
  for (const auto& x : entries) {
    double v[nV];
    for (int i = 0; i < nV; ++i) v[i] = inRange(x.v[i], i);
    for (int i = 0; i < nV; ++i) if (defined(x, i)) (x.signal ? hS : hB)[i]->Fill(v[i], x.w);
    for (int i = 1; i < nV; ++i)
      for (int j = 0; j < i; ++j) {
        if (!defined(x, i) || !defined(x, j)) continue;
        if (!x.signal) { hB2[i][j]->Fill(v[j], v[i], x.w); continue; }
        gS[i][j]->SetPoint(gS[i][j]->GetN(), v[j] + (integer[j] ? spread.Uniform(-0.3, 0.3) : 0.),
                           v[i] + (integer[i] ? spread.Uniform(-0.3, 0.3) : 0.));
      }
  }

  // ---- the corner plot
  gStyle->SetOptStat(0);
  {   // B density: nearly white to mid grey, so that the S points show everywhere
    double stops[2] = {0., 1.}, grey[2] = {0.97, 0.25};
    TColor::CreateGradientColorTable(2, stops, grey, grey, grey, 64);
  }
  TCanvas c("c", "corner plot of the cut variables", 1600, 1600);
  const double cell = 1. / nV;
  // where the cuts are, what they reject (diagonal), and the cuts of the 1 cm target where these differ
  const std::vector<double> bounds[nV] = {{0.5}, {1.5, 6.5}, {eMidCut}, {2.5}, {tofCut}};
  const std::vector<std::pair<double, double>> rejected[nV] = {{{0.5, 1.5}}, {{-0.5, 1.5}, {6.5, 8.5}}, {{0., eMidCut}},
                                                               {{-0.5, 2.5}}, {{0., tofCut}}};
  std::vector<double> refBounds[nV];
  if (std::abs(eMidCut - kCutE1cm) > 1e-6) refBounds[2] = {kCutE1cm};
  if (std::abs(tofCut - kCutT1cm) > 1e-6) refBounds[4] = {kCutT1cm};
  auto axisStyle = [](TAxis* a, bool shown, const char* title) {
    a->SetLabelSize(shown ? 0.075 : 0.);
    a->SetTitleSize(shown ? 0.08 : 0.);
    a->SetTitle(shown ? title : "");
    a->SetTitleOffset(1.15);
    a->SetNdivisions(505);
  };
  auto line = [](double x1, double y1, double x2, double y2, int col, int style) {
    auto l = new TLine(x1, y1, x2, y2);
    l->SetLineColor(col); l->SetLineStyle(style); l->SetLineWidth(2);
    l->Draw();
    return l;
  };
  TLine* cutLine = nullptr;
  TLine* refLine = nullptr;
  for (int r = 0; r < nV; ++r)
    for (int col = 0; col <= r; ++col) {
      c.cd();
      auto pad = new TPad(Form("pad%d%d", r, col), "", col * cell, 1. - (r + 1) * cell, (col + 1) * cell, 1. - r * cell);
      pad->SetLeftMargin(0.2); pad->SetBottomMargin(0.2); pad->SetRightMargin(0.03); pad->SetTopMargin(0.03);
      pad->Draw();
      pad->cd();
      const bool bottomRow = r == nV - 1;
      if (r == col) {   // the variable's shapes, unit area
        TH1D *b = hB[r], *s = hS[r];
        for (auto h : {b, s}) if (h->Integral() > 0.) h->Scale(1. / h->Integral());
        const double top = 1.45 * std::max(b->GetMaximum(), s->GetMaximum());
        b->SetMinimum(0.); b->SetMaximum(top);
        b->SetFillColor(kGray); b->SetLineColor(kGray + 2);
        s->SetLineColor(kRed + 1); s->SetLineWidth(2);
        axisStyle(b->GetXaxis(), bottomRow, names[r]);
        axisStyle(b->GetYaxis(), false, "");
        if (r == 0 && bottomRow) { b->GetXaxis()->SetBinLabel(1, "none"); b->GetXaxis()->SetBinLabel(2, "hit"); }
        b->Draw("HIST");
        for (const auto& [x1, x2] : rejected[r]) {
          auto box = new TBox(std::max(x1, lo[r]), 0., std::min(x2, hi[r]), top / 1.45);
          box->SetFillStyle(3354); box->SetFillColor(kGray + 3); box->SetLineWidth(0);
          box->Draw();
        }
        s->Draw("HIST SAME");
        for (double x : bounds[r]) cutLine = line(x, 0., x, top / 1.45, kBlack, 2);
        for (double x : refBounds[r]) refLine = line(x, 0., x, top / 1.45, kGray + 2, 3);
        TLatex text;
        text.SetNDC(); text.SetTextSize(0.068); text.SetTextColor(kBlack);
        text.DrawLatex(0.23, 0.89, Form("alone: S %.0f%%, B %.0f%%", 100. * sAlone[r] / sAll, 100. * bAlone[r] / bAll));
        text.DrawLatex(0.23, 0.81, Form("N-1: S %.0f%%, B %.0f%%", 100. * sNm1[r] / std::max(1e-12, sOthers[r]),
                                        100. * bNm1[r] / std::max(1e-12, bOthers[r])));
        gPad->RedrawAxis();
      } else {          // B density and S points of the pair: x the column's variable, y the row's
        TH2D* h = hB2[r][col];
        pad->SetLogz();
        axisStyle(h->GetXaxis(), bottomRow, names[col]);
        axisStyle(h->GetYaxis(), col == 0, names[r]);
        if (col == 0 && bottomRow) { h->GetXaxis()->SetBinLabel(1, "none"); h->GetXaxis()->SetBinLabel(2, "hit"); }
        h->SetMinimum(1e-3 * h->GetMaximum());   // three decades of B density
        h->Draw("COL");
        gS[r][col]->SetMarkerStyle(20); gS[r][col]->SetMarkerSize(0.45); gS[r][col]->SetMarkerColor(kRed + 1);
        if (gS[r][col]->GetN() > 0) gS[r][col]->Draw("P SAME");
        for (double x : bounds[col]) line(x, lo[r], x, hi[r], kBlack, 2);
        for (double y : bounds[r]) line(lo[col], y, hi[col], y, kBlack, 2);
        for (double x : refBounds[col]) line(x, lo[r], x, hi[r], kGray + 2, 3);
        for (double y : refBounds[r]) line(lo[col], y, hi[col], y, kGray + 2, 3);
        gPad->RedrawAxis();
      }
    }

  // the free upper right: what is shown (above), the numbers (below)
  c.cd();
  auto info = new TPad("info", "", 2. * cell + 0.03, 1. - 2. * cell + 0.03, 0.99, 0.99);
  info->Draw();
  info->cd();
  TLatex text;
  text.SetNDC(); text.SetTextFont(62); text.SetTextSize(0.075);
  text.DrawLatex(0.02, 0.9, Form("Cut variables, Al L-alpha window, %s target", label.Data()));
  text.SetTextFont(42); text.SetTextSize(0.052);
  text.DrawLatex(0.02, 0.81, "entries with a middle-TOF hit in coincidence (cut flow stage 1)");
  auto marker = static_cast<TGraph*>(gS[1][0]->Clone("marker"));   // a legend marker big enough to see
  marker->SetMarkerSize(1.4);
  auto legend = new TLegend(0.02, 0.44, 0.75, 0.77);
  legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(0.052);
  legend->AddEntry(marker, "target X-rays (S): points, line on the diagonal", "p");
  legend->AddEntry(hB[0], "everything else (B): log density, area on the diagonal", "f");
  legend->AddEntry(cutLine, "the cuts", "l");
  if (refLine) legend->AddEntry(refLine, Form("the cuts of the 1 cm target: %.1f MeV, %.1f ns", kCutE1cm, kCutT1cm), "l");
  legend->Draw();
  text.DrawLatex(0.02, 0.33, "kept: no bottom-TOF hit, middle bars 2-6 (over the target),");
  text.DrawLatex(0.02, 0.25, Form("middle-TOF energy > %.2f MeV, highest stack layer 3 or 4 (the top", eMidCut));
  text.DrawLatex(0.02, 0.17, Form("YSO and Si), top-to-middle time > %.2f ns where the top TOF timed it", tofCut));
  c.cd();
  auto numbers = new TPad("numbers", "", 3. * cell + 0.03, 1. - 3. * cell + 0.03, 0.99, 1. - 2. * cell + 0.01);
  numbers->Draw();
  numbers->cd();
  text.SetTextSize(0.1);
  text.DrawLatex(0.02, 0.85, Form("per day before these cuts: S %.1f, B %.0f", sAll, bAll));
  text.DrawLatex(0.02, 0.71, Form("after all of them: S %.1f, B %.0f", sPass, bPass));
  text.DrawLatex(0.02, 0.57, Form("S/sqrt(B) in 21 days: %.1f", bPass > 0. ? std::sqrt(21.) * sPass / std::sqrt(bPass) : 0.));
  text.DrawLatex(0.02, 0.4, "diagonal: unit area, hatched = rejected;");
  text.DrawLatex(0.02, 0.26, "kept by the cut alone, and N-1: after the others");
  text.DrawLatex(0.02, 0.1, Form("assumed TOF: %.1f ns per station, %.0f%% energy", sigmaT, 100. * sigmaE));
  const TString out = TString(gSystem->GetDirName(skim)) + "/muonic_xray_corner";
  c.Print(out + ".pdf");
  c.SaveAs(out + ".png");
  printf("plot: %s.png (and .pdf)\n", out.Data());
}
