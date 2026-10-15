// Muon spectrum before and after a lead slab, from scripts/lead_filter.sh (macro/lead_filter.mac):
// the EXPACS sea-level muons (from 9.25 MeV, zenith < 75 deg) on an 8 x 8 m lead slab in vacuum,
// files <dir>/lead_filter_t<cm>.root, or lead_filter_t<cm>_p<k>.root for a thickness run in parts
// (added up: the live times add).
//   root -l -b -q 'analysis/lead_filter_spectrum.C("output/lead_filter")'
//   root -l -b -q 'analysis/lead_filter_spectrum.C("output/lead_filter_below2GeV", 2)'
// before: kinetic energy when the muon reaches the slab (ekin; the source plane lies on it)
// after:  kinetic energy when it leaves the slab from below (end_ke: vacuum below, so it keeps it
//         to the end of the world); muons stopped in the lead have no "after"
// Both as the flux through a horizontal surface, per m2, s and GeV (each event weighs
// 1 / (live time x plane area)). "before" is the same for every slab, so it is averaged over them.
// below = 0: the whole spectrum was simulated (108 GeV); four panels: the spectra on log and
//   linear scales, the probability to get through versus the energy before, and the flux after /
//   before versus thickness, with a table of the stopped share, energy loss and slow muons.
// below = E (GeV): only the muons that can leave the lead below E were simulated
//   (scripts/lead_filter.sh --below E), so only the spectra below E are shown: linear in 20 MeV
//   bins, and log-log in 40 bins per decade from 1 MeV, with the flux of slow muons and its
//   statistical error.
// Writes <dir>/lead_filter_spectrum.pdf, .png and .csv (below E also _loglog.png and _loglog.csv).

namespace {
struct Slab {
  double t = 0., liveTime = 0., area = 0., eMax = 0.;   // cm, s, cm2, MeV
  std::vector<TString> files;
  TH1D *before = nullptr, *after = nullptr, *through = nullptr, *beforeLin = nullptr, *afterLin = nullptr;
  TH1D *beforeFine = nullptr, *afterFine = nullptr;   // 40 bins per decade from 1 MeV to the top of the linear range
  double nBefore = 0., nAfter = 0., nStopped = 0., nOther = 0.;   // per m2 per s
  double loss = 0.;   // sum of w x energy lost (MeV) by the muons that get through
};

TH1D* PerGeV(const TH1D* h, const char* name)
{
  auto d = static_cast<TH1D*>(h->Clone(name));
  for (int b = 1; b <= d->GetNbinsX(); ++b) {
    const double width = h->GetBinWidth(b) / 1000.;
    d->SetBinContent(b, h->GetBinContent(b) / width);
    d->SetBinError(b, h->GetBinError(b) / width);
  }
  return d;
}

// energy (MeV) below which a fraction q of the flux lies
double Quantile(const TH1D* h, double q)
{
  const double total = h->Integral();
  double sum = 0.;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    const double c = h->GetBinContent(b);
    if (c > 0. && sum + c >= q * total) {
      const double lo = h->GetBinLowEdge(b), hi = h->GetBinLowEdge(b + 1);
      return lo * std::pow(hi / lo, (q * total - sum) / c);
    }
    sum += c;
  }
  return 0.;
}
}  // namespace

void lead_filter_spectrum(const char* dir = "output/lead_filter", double below = 0.)
{
  TH1::AddDirectory(false);   // the histograms must outlive the input files

  // ---- the files, grouped by thickness: lead_filter_t<cm>.root or lead_filter_t<cm>_p<k>.root
  std::map<double, Slab> byThickness;
  if (void* d = gSystem->OpenDirectory(dir)) {
    while (const char* entry = gSystem->GetDirEntry(d)) {
      TString name = entry;
      if (!name.BeginsWith("lead_filter_t") || !name.EndsWith(".root")) continue;
      TString rest = name(13, name.Length() - 18);   // between "lead_filter_t" and ".root"
      if (rest.Index("_p") != kNPOS) rest = rest(0, rest.Index("_p"));
      Slab& s = byThickness[rest.Atof()];
      s.t = rest.Atof();
      s.files.push_back(Form("%s/%s", dir, entry));
    }
    gSystem->FreeDirectory(d);
  }
  if (byThickness.empty()) {
    printf("no lead_filter_t<cm>.root in %s (run sh scripts/lead_filter.sh)\n", dir);
    return;
  }

  std::vector<double> edges;   // 20 bins per decade, 5 MeV - 200 GeV
  for (double e = 5.; e < 2.0001e5; e *= std::pow(10., 0.05)) edges.push_back(e);
  const double linTop = below > 0. ? 1000. * below : 2000.;   // MeV
  const int linBins = int(std::lround(linTop / 20.));
  std::vector<double> fine;   // 1 MeV - linTop, 40 bins per decade, ending exactly on linTop
  const int nFine = int(std::ceil(40. * std::log10(linTop)));
  for (int k = 0; k <= nFine; ++k) fine.push_back(std::pow(linTop, double(k) / nFine));
  std::vector<Slab> slabs;
  long nEvents = 0;
  for (auto& [t, s] : byThickness) {
    // total live time of the parts first, so that every event weighs 1 / (sum of live times x area)
    std::vector<TString> good;
    for (const auto& path : s.files) {
      TFile file(path);
      auto run = dynamic_cast<TTree*>(file.Get("run"));
      if (!run || !file.Get("events") || run->GetEntries() != 1) {
        printf("skipping %s: incomplete\n", path.Data());
        continue;
      }
      double liveTime = 0., area = 0., eMax = 0.;
      run->SetBranchAddress("live_time", &liveTime);
      run->SetBranchAddress("area", &area);
      run->SetBranchAddress("e_max", &eMax);
      run->GetEntry(0);
      s.liveTime += liveTime;
      s.area = area;
      s.eMax = eMax;
      good.push_back(path);
    }
    if (good.empty() || s.liveTime <= 0.) continue;
    const double w = 1. / (s.liveTime * s.area * 1e-4);   // per m2 per s
    s.before = new TH1D(Form("before%g", t), "", int(edges.size()) - 1, edges.data());
    s.after = new TH1D(Form("after%g", t), "", int(edges.size()) - 1, edges.data());
    s.through = new TH1D(Form("through%g", t), "", int(edges.size()) - 1, edges.data());
    s.beforeLin = new TH1D(Form("beforeLin%g", t), "", linBins, 0., linTop);
    s.afterLin = new TH1D(Form("afterLin%g", t), "", linBins, 0., linTop);
    s.beforeFine = new TH1D(Form("beforeFine%g", t), "", nFine, fine.data());
    s.afterFine = new TH1D(Form("afterFine%g", t), "", nFine, fine.data());
    TH1D* hs[] = {s.before, s.after, s.through, s.beforeLin, s.afterLin, s.beforeFine, s.afterFine};
    for (TH1D* h : hs) h->Sumw2();
    for (const auto& path : good) {
      TFile file(path);
      auto events = dynamic_cast<TTree*>(file.Get("events"));
      double ekin = 0., endKe = 0., endZ = 0.;
      events->SetBranchStatus("*", false);
      const char* used[] = {"ekin", "end_ke", "end_z"};
      for (const char* b : used) events->SetBranchStatus(b, true);
      events->SetBranchAddress("ekin", &ekin);
      events->SetBranchAddress("end_ke", &endKe);
      events->SetBranchAddress("end_z", &endZ);
      const Long64_t n = events->GetEntries();
      for (Long64_t i = 0; i < n; ++i) {
        events->GetEntry(i);
        s.before->Fill(ekin, w);
        s.beforeLin->Fill(ekin, w);
        s.beforeFine->Fill(ekin, w);
        s.nBefore += w;
        if (endZ < -0.5 * t - 1e-3) {   // left the slab from below
          s.after->Fill(endKe, w);
          s.afterLin->Fill(endKe, w);
          s.afterFine->Fill(endKe, w);
          s.through->Fill(ekin, w);
          s.nAfter += w;
          s.loss += w * (ekin - endKe);
        } else if (endKe < 1e-3 && std::abs(endZ) <= 0.5 * t + 1e-3) {
          s.nStopped += w;
        } else {
          s.nOther += w;   // came back out of the top
        }
      }
      nEvents += n;
    }
    slabs.push_back(s);
  }
  if (slabs.empty()) return;

  // the spectrum before the lead, averaged over the slabs (the same muons arrive at every slab)
  auto before = static_cast<TH1D*>(slabs[0].before->Clone("beforeAll"));
  auto beforeLin = static_cast<TH1D*>(slabs[0].beforeLin->Clone("beforeLinAll"));
  auto beforeFine = static_cast<TH1D*>(slabs[0].beforeFine->Clone("beforeFineAll"));
  before->Reset();
  beforeLin->Reset();
  beforeFine->Reset();
  for (const auto& s : slabs) {
    before->Add(s.before, 1. / slabs.size());
    beforeLin->Add(s.beforeLin, 1. / slabs.size());
    beforeFine->Add(s.beforeFine, 1. / slabs.size());
  }
  double fluxBefore = 0.;
  for (const auto& s : slabs) fluxBefore += s.nBefore / slabs.size();
  auto upTo = [](const TH1D* h, double e, double& err) {
    return h->IntegralAndError(1, h->FindFixBin(e * 0.999999), err);
  };
  const int colours[] = {kBlue + 1, kAzure + 7, kGreen + 2, kOrange + 7, kRed + 1, kMagenta + 2, kGray + 2, kCyan + 3};
  gStyle->SetOptStat(0);
  const TString base = Form("%s/lead_filter_spectrum", dir);
  printf("\n=== Sea-level muons before and after a lead slab (%s): %ld muons, %zu thicknesses ===\n", dir, nEvents, slabs.size());

  if (below > 0.) {
    // ---- only the muons that can leave below E: spectra below E, linear
    printf("EXPACS mu+ and mu-, from 9.25 MeV, zenith < 75 deg; only the muons that can leave the lead below %g GeV\n", below);
    printf("Flux through a horizontal surface below a kinetic energy, per m2 per s (statistical error)\n");
    printf("%9s %11s %22s %22s %22s %22s\n", "lead (cm)", "up to (GeV)", "below 100 MeV", "below 300 MeV",
           "below 1 GeV", Form("below %g GeV", below));
    auto row = [&](const char* label, double eMax, const TH1D* h) {
      printf("%9s %11.3g", label, eMax / 1000.);
      const double limits[] = {100., 300., 1000., linTop};
      for (double l : limits) {
        double e = 0.;
        const double v = upTo(h, l, e);
        printf(" %13.4f +- %6.4f", v, e);
      }
      printf("\n");
    };
    row("before", 0., beforeLin);
    for (const auto& s : slabs) row(Form("%g", s.t), s.eMax, s.afterLin);
    const long perThickness = nEvents / long(slabs.size());
    printf("(about %ld simulated muons per thickness; \"up to\": top of the simulated spectrum)\n", perThickness);

    const TString csvName = base + ".csv";
    if (FILE* csv = fopen(csvName, "w")) {
      fprintf(csv, "E_low_MeV,E_high_MeV,before_m2sGeV,before_err");
      for (const auto& s : slabs) fprintf(csv, ",after_%gcm_m2sGeV,after_%gcm_err", s.t, s.t);
      fprintf(csv, "\n");
      for (int b = 1; b <= beforeLin->GetNbinsX(); ++b) {
        const double width = beforeLin->GetBinWidth(b) / 1000.;
        fprintf(csv, "%g,%g,%.6g,%.3g", beforeLin->GetBinLowEdge(b), beforeLin->GetBinLowEdge(b + 1),
                beforeLin->GetBinContent(b) / width, beforeLin->GetBinError(b) / width);
        for (const auto& s : slabs)
          fprintf(csv, ",%.6g,%.3g", s.afterLin->GetBinContent(b) / width, s.afterLin->GetBinError(b) / width);
        fprintf(csv, "\n");
      }
      fclose(csv);
    }

    TCanvas c("c", "lead filter, low energies", 1100, 800);
    gPad->SetGridy();
    gPad->SetBottomMargin(0.12);
    gPad->SetLeftMargin(0.11);
    TH1D* hb = PerGeV(beforeLin, "hbLin");
    hb->SetTitle(Form("muon spectrum before and after the lead, below %g GeV (%.2g muons per thickness);kinetic energy (MeV);dJ/dE (m^{-2} s^{-1} GeV^{-1})",
                      below, double(perThickness)));
    hb->SetLineColor(kBlack);
    hb->SetLineWidth(3);
    hb->SetMinimum(0.);
    double top = hb->GetMaximum();
    std::vector<TH1D*> curves;
    for (std::size_t k = 0; k < slabs.size(); ++k) {
      TH1D* h = PerGeV(slabs[k].afterLin, Form("haLin%zu", k));
      h->SetLineColor(colours[k % 8]);
      h->SetLineWidth(2);
      top = std::max(top, h->GetMaximum());
      curves.push_back(h);
    }
    hb->SetMaximum(1.18 * top);
    hb->GetXaxis()->SetTitleOffset(1.2);
    hb->Draw("HIST");
    for (TH1D* h : curves) h->Draw("HIST SAME");
    hb->Draw("HIST SAME");
    auto legend = new TLegend(0.62, 0.6, 0.89, 0.89);
    legend->SetTextSize(0.03);
    legend->AddEntry(hb, "before (sea level)", "l");
    for (std::size_t k = 0; k < slabs.size(); ++k) legend->AddEntry(curves[k], Form("after %g cm of lead", slabs[k].t), "l");
    legend->Draw();

    // the same on log-log axes, 40 bins per decade from 1 MeV
    TCanvas c2("c2", "lead filter, low energies, log-log", 1100, 800);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetBottomMargin(0.12);
    gPad->SetLeftMargin(0.11);
    TH1D* hbf = PerGeV(beforeFine, "hbFine");
    hbf->SetTitle(Form("muon spectrum before and after the lead, below %g GeV, log-log (%.2g muons per thickness);kinetic energy (MeV);dJ/dE (m^{-2} s^{-1} GeV^{-1})",
                       below, double(perThickness)));
    hbf->SetLineColor(kBlack);
    hbf->SetLineWidth(3);
    double fineMax = hbf->GetMaximum(), fineMin = 1e30;
    std::vector<TH1D*> fineCurves;
    for (std::size_t k = 0; k < slabs.size(); ++k) {
      TH1D* h = PerGeV(slabs[k].afterFine, Form("haFine%zu", k));
      h->SetLineColor(colours[k % 8]);
      h->SetLineWidth(2);
      fineMax = std::max(fineMax, h->GetMaximum());
      for (int b = 1; b <= h->GetNbinsX(); ++b)
        if (h->GetBinContent(b) > 0.) fineMin = std::min(fineMin, h->GetBinContent(b));
      fineCurves.push_back(h);
    }
    auto frame2 = gPad->DrawFrame(fine.front(), 0.5 * fineMin, fine.back(), 2.5 * fineMax);
    frame2->SetTitle(hbf->GetTitle());
    frame2->GetXaxis()->SetTitle("kinetic energy (MeV)");
    frame2->GetYaxis()->SetTitle("dJ/dE (m^{-2} s^{-1} GeV^{-1})");
    frame2->GetXaxis()->SetTitleOffset(1.2);
    for (TH1D* h : fineCurves) h->Draw("HIST SAME");
    // the simulated spectrum starts at 9.25 MeV: draw "before" from its first full bin, without the drop to 0
    int first = 1;
    while (first < hbf->GetNbinsX() && hbf->GetBinLowEdge(first) < 9.249) ++first;
    hbf->GetXaxis()->SetRange(first, hbf->GetNbinsX());
    hbf->Draw("HIST SAME ][");
    auto legend2 = new TLegend(0.52, 0.15, 0.89, 0.45);
    legend2->SetTextSize(0.03);
    legend2->AddEntry(hbf, "before (sea level, from 9.25 MeV)", "l");
    for (std::size_t k = 0; k < slabs.size(); ++k) legend2->AddEntry(fineCurves[k], Form("after %g cm of lead", slabs[k].t), "l");
    legend2->Draw();
    const TString csvLog = base + "_loglog.csv";
    if (FILE* csv = fopen(csvLog, "w")) {
      fprintf(csv, "E_low_MeV,E_high_MeV,before_m2sGeV,before_err");
      for (const auto& s : slabs) fprintf(csv, ",after_%gcm_m2sGeV,after_%gcm_err", s.t, s.t);
      fprintf(csv, "\n");
      for (int b = 1; b <= hbf->GetNbinsX(); ++b) {
        fprintf(csv, "%.6g,%.6g,%.6g,%.3g", hbf->GetBinLowEdge(b), hbf->GetBinLowEdge(b + 1), hbf->GetBinContent(b), hbf->GetBinError(b));
        for (TH1D* h : fineCurves) fprintf(csv, ",%.6g,%.3g", h->GetBinContent(b), h->GetBinError(b));
        fprintf(csv, "\n");
      }
      fclose(csv);
    }
    c.Print(base + ".pdf(");
    c2.Print(base + ".pdf)");
    c.SaveAs(base + ".png");
    c2.SaveAs(base + "_loglog.png");
    printf("plots: %s.pdf, %s.png, %s_loglog.png; tables: %s, %s\n", base.Data(), base.Data(), base.Data(),
           csvName.Data(), csvLog.Data());
    return;
  }

  // ---- the whole spectrum: table
  printf("EXPACS mu+ and mu-, 9.25 MeV - 108 GeV, zenith < 75 deg; flux through a horizontal surface, per m2 per s\n");
  printf("%9s %10s %9s %9s %12s %13s %13s %10s %10s %10s\n", "lead (cm)", "flux", "through", "stopped",
         "mean loss", "median E", "median E", "E<100 MeV", "E<300 MeV", "E<1 GeV");
  printf("%9s %10s %9s %9s %12s %13s %13s %10s %10s %10s\n", "", "after", "share", "share", "(MeV)", "before (GeV)",
         "after (GeV)", "after", "after", "after");
  double e = 0.;
  printf("%9s %10.2f %9s %9s %12s %13.3f %13s %10.2f %10.2f %10.2f\n", "before", fluxBefore, "", "", "",
         Quantile(before, 0.5) / 1000., "", upTo(beforeLin, 100., e), upTo(beforeLin, 300., e), upTo(beforeLin, 1000., e));
  for (const auto& s : slabs)
    printf("%9g %10.2f %9.3f %9.3f %12.0f %13.3f %13.3f %10.2f %10.2f %10.2f\n", s.t, s.nAfter, s.nAfter / s.nBefore,
           s.nStopped / s.nBefore, s.nAfter > 0. ? s.loss / s.nAfter : 0., Quantile(before, 0.5) / 1000.,
           Quantile(s.after, 0.5) / 1000., upTo(s.afterLin, 100., e), upTo(s.afterLin, 300., e), upTo(s.afterLin, 1000., e));
  printf("(\"E<x after\": flux of the muons that leave the lead with less than x, per m2 per s, from 20 MeV bins;\n"
         " the medians from 20 bins per decade. Muons coming back out of the top: below 0.1%% of the flux)\n");

  const TString csvName = base + ".csv";
  if (FILE* csv = fopen(csvName, "w")) {
    fprintf(csv, "E_low_MeV,E_high_MeV,before_m2s");
    for (const auto& s : slabs) fprintf(csv, ",after_%gcm_m2s,through_prob_%gcm", s.t, s.t);
    fprintf(csv, "\n");
    for (int b = 1; b <= before->GetNbinsX(); ++b) {
      fprintf(csv, "%.6g,%.6g,%.6g", before->GetBinLowEdge(b), before->GetBinLowEdge(b + 1), before->GetBinContent(b));
      for (const auto& s : slabs)
        fprintf(csv, ",%.6g,%.6g", s.after->GetBinContent(b),
                s.before->GetBinContent(b) > 0. ? s.through->GetBinContent(b) / s.before->GetBinContent(b) : 0.);
      fprintf(csv, "\n");
    }
    fclose(csv);
  }

  // ---- the whole spectrum: four panels
  TCanvas c("c", "lead filter", 1500, 1100);
  c.Divide(2, 2);
  std::vector<TH1D*> keep;

  // 1: log-log spectra
  c.cd(1);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetBottomMargin(0.12);
  TH1D* hb = PerGeV(before, "hbLog");
  hb->SetTitle("muon spectrum before and after the lead;kinetic energy (MeV);dJ/dE (m^{-2} s^{-1} GeV^{-1})");
  hb->SetLineColor(kBlack);
  hb->SetLineWidth(3);
  hb->GetXaxis()->SetRangeUser(5., 1e5);   // the simulated spectrum ends at 108 GeV
  hb->GetXaxis()->SetTitleOffset(1.2);
  hb->SetMinimum(1e-3);
  hb->SetMaximum(3. * hb->GetMaximum());
  hb->Draw("HIST");
  auto legend = new TLegend(0.14, 0.15, 0.48, 0.47);
  legend->SetTextSize(0.035);
  legend->AddEntry(hb, "before (sea level)", "l");
  for (std::size_t k = 0; k < slabs.size(); ++k) {
    TH1D* h = PerGeV(slabs[k].after, Form("haLog%zu", k));
    h->SetLineColor(colours[k % 8]);
    h->SetLineWidth(2);
    h->Draw("HIST SAME");
    legend->AddEntry(h, Form("after %g cm of lead", slabs[k].t), "l");
    keep.push_back(h);
  }
  legend->Draw();

  // 2: linear, below 2 GeV
  c.cd(2);
  gPad->SetGridy();
  gPad->SetBottomMargin(0.12);
  TH1D* hbl = PerGeV(beforeLin, "hbLin");
  hbl->SetTitle("below 2 GeV, linear scale;kinetic energy (MeV);dJ/dE (m^{-2} s^{-1} GeV^{-1})");
  hbl->SetLineColor(kBlack);
  hbl->SetLineWidth(3);
  hbl->SetMinimum(0.);
  hbl->SetMaximum(1.25 * hbl->GetMaximum());
  hbl->GetXaxis()->SetTitleOffset(1.2);
  hbl->Draw("HIST");
  for (std::size_t k = 0; k < slabs.size(); ++k) {
    TH1D* h = PerGeV(slabs[k].afterLin, Form("haLin%zu", k));
    h->SetLineColor(colours[k % 8]);
    h->SetLineWidth(2);
    h->Draw("HIST SAME");
    keep.push_back(h);
  }

  // 3: probability to get through versus the energy before
  c.cd(3);
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetBottomMargin(0.12);
  auto frame = gPad->DrawFrame(5., 0., 1e5, 1.08);
  frame->SetTitle("probability to get through the lead;kinetic energy before the lead (MeV);fraction leaving the slab below");
  frame->GetXaxis()->SetTitleOffset(1.2);
  for (std::size_t k = 0; k < slabs.size(); ++k) {
    auto p = static_cast<TH1D*>(slabs[k].through->Clone(Form("prob%zu", k)));
    p->Divide(p, slabs[k].before, 1., 1., "B");
    p->SetLineColor(colours[k % 8]);
    p->SetLineWidth(2);
    p->Draw("HIST SAME");
    keep.push_back(p);
  }
  auto half = new TLine(5., 0.5, 1e5, 0.5);
  half->SetLineStyle(2);
  half->Draw();

  // 4: flux after / before versus thickness, all and slow muons
  c.cd(4);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetBottomMargin(0.12);
  auto mg = new TMultiGraph("mg", "muon flux after the lead / before;lead thickness (cm);after / before");
  auto legend4 = new TLegend(0.14, 0.15, 0.62, 0.36);
  legend4->SetTextSize(0.035);
  const double limits[] = {1e30, 1000., 300., 100.};
  const char* names[] = {"all muons", "below 1 GeV", "below 300 MeV", "below 100 MeV"};
  const int lcol[] = {kBlack, kBlue + 1, kGreen + 2, kRed + 1};
  for (int l = 0; l < 4; ++l) {
    auto g = new TGraph();
    g->SetPoint(0, 0., 1.);
    for (const auto& s : slabs) {
      const double a = l == 0 ? s.nAfter : upTo(s.afterLin, limits[l], e);
      const double b = l == 0 ? s.nBefore : upTo(beforeLin, limits[l], e);
      if (b > 0.) g->SetPoint(g->GetN(), s.t, a / b);
    }
    g->SetMarkerStyle(20);
    g->SetMarkerColor(lcol[l]);
    g->SetLineColor(lcol[l]);
    g->SetLineWidth(2);
    mg->Add(g, "PL");
    legend4->AddEntry(g, names[l], "pl");
  }
  mg->Draw("A");
  mg->GetYaxis()->SetRangeUser(0., 1.3);
  mg->GetXaxis()->SetTitleOffset(1.2);
  legend4->Draw();

  c.Print(base + ".pdf");
  c.SaveAs(base + ".png");
  printf("plots: %s.pdf, %s.png; table: %s\n", base.Data(), base.Data(), csvName.Data());
}
