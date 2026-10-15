// Muon energy after the top lead: spectra of the top & middle & Compton coincidences (TMC, the
// selection of analysis/lead_energy_scan.C) versus lead thickness, from a lead x energy scan whose
// events have the columns pb_ke, top_ke and cmp_ke (macro/lead_thick_scan.mac):
//   root -l -b -q 'analysis/lead_after_spectrum.C("output/leadscan_thick")'
//   root -l -b -q 'analysis/lead_after_spectrum.C("output/leadscan_thick", 0.2, 0.1)'   (thresholds, MeV)
// Kinetic energies of the primary muon of each TMC event:
//   ekin    at the source plane (sea level)
//   pb_ke   where it leaves the lead layer (-1 if it went round the layer, or without lead; the
//           0 cm curve of that panel is the energy at sea level, where the lead would be)
//   top_ke  where it enters a bar of the top TOF (-1 if the top bar was fired by secondaries)
//   cmp_ke  where it enters the Compton detector (-1 if the coincidence came from secondaries)
// Mono-energetic points (leadscan_t<cm>_E<GeV>.root): each event weighs 3600 / live time of its
// run, so the sum over the points of a thickness gives rates per hour due to the sea-level muons
// of the whole band range. A point stands for its band E_low - E_high, so above 50 MeV each of its
// events is spread evenly over [x - (E - E_low), x + (E_high - E)] around its energy x, the band's
// muons arriving shifted by their offset from E; without this, thin lead gives a comb of the 60
// energies. The part of that interval below 50 MeV is folded back above 50 MeV, and events below
// 50 MeV stay at their own energies: there the muons stop or fail the selection depending on where
// they stop in the detector, which the selected events cannot extrapolate, and behind lead their
// own spread of angles and energy loss is already wide. So no rate moves across 50 MeV.
// The energy bins are the bands of the scan (the sea-level spectrum then
// has one point per bin), continued down to 1 MeV with the same width in log E. The spectrum
// check runs (leadspec_t<cm>.root, the continuous spectrum) give the same rates with fewer slow
// muons; the table and the open markers compare them.
// muonOnly = true: the TMC events must also have the muon itself (the primary) crossing the top and
// middle TOF (prim_st) and entering the Compton detector (cmp_ke >= 0); output names with "_muon".
// Writes <dir>/lead_after_spectrum.csv (rates per thickness and bin), <dir>/lead_after_spectrum.pdf
// and a .png per page (_spectra, _maps, _bythickness: the top-TOF spectrum of each thickness on its
// own, with its check run and the 0 cm spectrum).

namespace {
struct Spectra {
  TH1D *ini = nullptr, *pb = nullptr, *top = nullptr, *cmp = nullptr, *topLin = nullptr;
  double tmc = 0., tmc2 = 0.;   // rate per hour and its squared error
  double withPb = 0., withTop = 0., withCmp = 0.;
  int files = 0;
};

TH1D* MakeHist(const char* name, const std::vector<double>& edges)
{
  auto h = new TH1D(name, "", int(edges.size()) - 1, edges.data());
  h->Sumw2();
  return h;
}

void MakeSpectra(Spectra& s, const char* tag, const std::vector<double>& edges)
{
  s.ini = MakeHist(Form("ini_%s", tag), edges);
  s.pb = MakeHist(Form("pb_%s", tag), edges);
  s.top = MakeHist(Form("top_%s", tag), edges);
  s.cmp = MakeHist(Form("cmp_%s", tag), edges);
  s.topLin = new TH1D(Form("topLin_%s", tag), "", 50, 0., 1000.);
  s.topLin->Sumw2();
}

// w spread evenly over [a, b] (MeV, a kept above 0), in proportion to the overlap with each bin
void FillSpread(TH1D* h, double a, double b, double w)
{
  a = std::max(a, 1e-3);
  if (b <= a) {
    h->Fill(a, w);
    return;
  }
  const TAxis* axis = h->GetXaxis();
  const double xmin = axis->GetXmin(), xmax = axis->GetXmax();
  if (a < xmin) h->Fill(0.5 * xmin, w * (std::min(b, xmin) - a) / (b - a));       // underflow
  if (b > xmax) h->Fill(2. * xmax, w * (b - std::max(a, xmax)) / (b - a));        // overflow
  const int first = std::max(1, axis->FindFixBin(a)), last = std::min(axis->GetNbins(), axis->FindFixBin(b));
  for (int k = first; k <= last; ++k) {
    const double lo = std::max(a, axis->GetBinLowEdge(k)), hi = std::min(b, axis->GetBinUpEdge(k));
    if (hi > lo) h->Fill(axis->GetBinCenter(k), w * (hi - lo) / (b - a));
  }
}

// an event of a mono-energetic point at energy x (MeV) for its band, see above: spread evenly over
// [x - below, x + above], the part below kSpreadMin folded back above it; at x if x is below
const double kSpreadMin = 50.;   // MeV
void FillBandShift(TH1D* h, double x, double below, double above, double w)
{
  const double a = x - below, b = x + above;
  if (x < kSpreadMin || b <= a) {
    h->Fill(x, w);
    return;
  }
  if (a >= kSpreadMin) {
    FillSpread(h, a, b, w);
    return;
  }
  FillSpread(h, kSpreadMin, b, w * (b - kSpreadMin) / (b - a));
  FillSpread(h, kSpreadMin, 2. * kSpreadMin - a, w * (kSpreadMin - a) / (b - a));   // folded back
}

// energy (MeV) below which a fraction q of the rate lies, from the bin contents
double Quantile(const TH1D* h, double q)
{
  const double total = h->Integral();
  if (total <= 0.) return 0.;
  double sum = 0.;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    const double c = h->GetBinContent(b);
    if (sum + c >= q * total && c > 0.) {
      const double f = (q * total - sum) / c, lo = h->GetBinLowEdge(b), hi = h->GetBinLowEdge(b + 1);
      return lo * std::pow(hi / lo, f);
    }
    sum += c;
  }
  return h->GetBinLowEdge(h->GetNbinsX() + 1);
}

// rate per hour and GeV, for log-log plots
TH1D* PerGeV(const TH1D* h, const char* name)
{
  auto d = static_cast<TH1D*>(h->Clone(name));
  for (int b = 1; b <= d->GetNbinsX(); ++b) {
    const double w = d->GetBinWidth(b) / 1000.;
    d->SetBinContent(b, h->GetBinContent(b) / w);
    d->SetBinError(b, h->GetBinError(b) / w);
  }
  return d;
}
}  // namespace

void lead_after_spectrum(const char* dir = "output/leadscan_thick", double barThr = 0.2, double comptonThr = 0.1,
                         bool muonOnly = false)
{
  const TString suffix = muonOnly ? "_muon" : "";
  // ---- files and their run summaries
  struct FileInfo {
    TString path;
    bool check = false;
    double lead = 0., eMono = 0., eMin = 0., eMax = 0., liveTime = 0.;
  };
  std::vector<FileInfo> files;
  void* dirp = gSystem->OpenDirectory(dir);
  if (!dirp) {
    printf("cannot open %s\n", dir);
    return;
  }
  int noColumns = 0;
  while (const char* entry = gSystem->GetDirEntry(dirp)) {
    const TString name = entry;
    const bool isPoint = name.BeginsWith("leadscan_t"), isCheck = name.BeginsWith("leadspec_t");
    if (!(isPoint || isCheck) || !name.EndsWith(".root")) continue;
    FileInfo f;
    f.path = Form("%s/%s", dir, entry);
    f.check = isCheck;
    TFile file(f.path);
    auto run = dynamic_cast<TTree*>(file.Get("run"));
    auto events = dynamic_cast<TTree*>(file.Get("events"));
    if (!run || !events || run->GetEntries() != 1 || !run->GetBranch("lead_cm") || !run->GetBranch("e_mono")) continue;
    if (!events->GetBranch("pb_ke") || !events->GetBranch("top_ke") || (muonOnly && !events->GetBranch("prim_st"))) {
      ++noColumns;
      continue;
    }
    run->SetBranchAddress("lead_cm", &f.lead);
    run->SetBranchAddress("e_mono", &f.eMono);
    run->SetBranchAddress("e_min", &f.eMin);
    run->SetBranchAddress("e_max", &f.eMax);
    run->SetBranchAddress("live_time", &f.liveTime);
    run->GetEntry(0);
    if (f.liveTime <= 0. || isCheck != (f.eMono == 0.)) continue;
    files.push_back(f);
  }
  gSystem->FreeDirectory(dirp);
  if (noColumns > 0)
    printf("%d files in %s have no pb_ke / top_ke columns (made before they existed); skipped\n", noColumns, dir);
  std::vector<double> leads;
  std::vector<std::pair<double, double>> bands;   // (low, high) MeV of each energy point
  for (const auto& f : files) {
    if (std::find_if(leads.begin(), leads.end(), [&](double x) { return std::abs(x - f.lead) < 1e-6; }) == leads.end())
      leads.push_back(f.lead);
    if (!f.check && std::find_if(bands.begin(), bands.end(), [&](const std::pair<double, double>& b) {
                      return std::abs(b.first - f.eMin) < 1e-6 * f.eMin; }) == bands.end())
      bands.push_back({f.eMin, f.eMax});
  }
  if (bands.empty()) {
    printf("no mono-energetic points with pb_ke / top_ke in %s\n", dir);
    return;
  }
  std::sort(leads.begin(), leads.end());
  std::sort(bands.begin(), bands.end());
  const int nT = int(leads.size());

  // ---- energy bins: the bands if they form a log grid, continued down to 1 MeV; else 20 per decade
  std::vector<double> edges;
  const double ratio = bands[0].second / bands[0].first;
  bool logGrid = true;
  for (std::size_t i = 0; i < bands.size(); ++i) {
    if (std::abs(bands[i].second / bands[i].first - ratio) > 1e-3 * ratio) logGrid = false;
    if (i > 0 && std::abs(bands[i].first - bands[i - 1].second) > 1e-3 * bands[i].first) logGrid = false;
  }
  if (logGrid) {
    std::vector<double> below;
    for (double e = bands[0].first / ratio; e > 1.; e /= ratio) below.push_back(e);
    below.push_back(below.empty() ? bands[0].first / ratio : below.back() / ratio);
    edges.assign(below.rbegin(), below.rend());
    for (const auto& b : bands) edges.push_back(b.first);
    edges.push_back(bands.back().second);
  } else {
    const double top = bands.back().second;
    for (double e = 1.; e < top * 1.0001; e *= std::pow(10., 0.05)) edges.push_back(e);
    if (edges.back() < top) edges.push_back(top);
  }

  std::vector<Spectra> grid(nT), check(nT);
  for (int it = 0; it < nT; ++it) {
    MakeSpectra(grid[it], Form("g%d", it), edges);
    MakeSpectra(check[it], Form("c%d", it), edges);
  }

  // ---- TMC events of every file
  long nEvents = 0;
  for (const auto& f : files) {
    const int it = int(std::find_if(leads.begin(), leads.end(), [&](double x) { return std::abs(x - f.lead) < 1e-6; }) - leads.begin());
    Spectra& s = f.check ? check[it] : grid[it];
    ++s.files;
    const double w = 3600. / f.liveTime;
    TFile file(f.path);
    auto events = dynamic_cast<TTree*>(file.Get("events"));
    std::vector<int>* barSt = nullptr;
    std::vector<double>* barE = nullptr;
    double ekin = 0., pb = 0., top = 0., cmp = 0., stack0 = 0., stack1 = 0., ch0 = 0., ch1 = 0., ch2 = 0.;
    int primSt = 0;
    events->SetBranchStatus("*", false);
    const char* used[] = {"bar_st", "bar_e", "ekin", "pb_ke", "top_ke", "cmp_ke", "e_stack0", "e_stack1", "e_ch0", "e_ch1", "e_ch2",
                          "prim_st"};
    for (const char* b : used) events->SetBranchStatus(b, true);
    events->SetBranchAddress("bar_st", &barSt);
    events->SetBranchAddress("bar_e", &barE);
    events->SetBranchAddress("ekin", &ekin);
    events->SetBranchAddress("pb_ke", &pb);
    events->SetBranchAddress("top_ke", &top);
    events->SetBranchAddress("cmp_ke", &cmp);
    events->SetBranchAddress("e_stack0", &stack0);
    events->SetBranchAddress("e_stack1", &stack1);
    events->SetBranchAddress("e_ch0", &ch0);
    events->SetBranchAddress("e_ch1", &ch1);
    events->SetBranchAddress("e_ch2", &ch2);
    if (events->GetBranch("prim_st")) events->SetBranchAddress("prim_st", &primSt);
    const Long64_t n = events->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
      events->GetEntry(i);
      bool hitTop = false, hitMiddle = false;
      for (std::size_t k = 0; k < barSt->size(); ++k) {
        if ((*barE)[k] <= barThr) continue;
        hitTop = hitTop || (*barSt)[k] == 0;
        hitMiddle = hitMiddle || (*barSt)[k] == 1;
      }
      const bool compton = stack0 > comptonThr || stack1 > comptonThr || ch0 > comptonThr || ch1 > comptonThr || ch2 > comptonThr;
      if (!hitTop || !hitMiddle || !compton) continue;
      if (muonOnly && ((primSt & 3) != 3 || cmp < 0.)) continue;   // the muon itself, not only its secondaries
      ++nEvents;
      s.tmc += w;
      s.tmc2 += w * w;
      if (f.check) {
        s.ini->Fill(ekin, w);
        if (pb >= 0.) s.pb->Fill(pb, w);
        if (top >= 0.) { s.top->Fill(top, w); s.topLin->Fill(top, w); }
        if (cmp >= 0.) s.cmp->Fill(cmp, w);
      } else {   // the band of the point, around each energy
        const double below = f.eMono - f.eMin, above = f.eMax - f.eMono;
        FillSpread(s.ini, f.eMin, f.eMax, w);
        if (pb >= 0.) FillBandShift(s.pb, pb, below, above, w);
        if (top >= 0.) {
          FillBandShift(s.top, top, below, above, w);
          FillBandShift(s.topLin, top, below, above, w);
        }
        if (cmp >= 0.) FillBandShift(s.cmp, cmp, below, above, w);
      }
      if (pb >= 0.) s.withPb += w;
      if (top >= 0.) s.withTop += w;
      if (cmp >= 0.) s.withCmp += w;
    }
  }
  int nPoints = 0, nChecks = 0;
  for (int it = 0; it < nT; ++it) { nPoints += grid[it].files; nChecks += check[it].files; }
  printf("\n=== Muon energy after the lead, top & middle & Compton coincidences%s (%s) ===\n",
         muonOnly ? " made by the muon itself" : "", dir);
  printf("%d mono-energetic points and %d check runs, lead %g-%g cm, %d energy bands %.4g-%.4g MeV; "
         "%ld TMC events; TOF bar > %.2f MeV, Compton layer > %.2f MeV\n", nPoints, nChecks, leads.front(),
         leads.back(), int(bands.size()), bands.front().first, bands.back().second, nEvents, barThr, comptonThr);

  // ---- table per thickness
  printf("\n%5s %18s %12s %9s %9s %9s %11s %11s %9s %9s %9s\n", "lead", "TMC per hour", "check run", "pb_ke", "top_ke",
         "cmp_ke", "median top", "median cmp", "top<100", "top<300", "top<1000");
  printf("%5s %18s %12s %9s %9s %9s %11s %11s %9s %9s %9s\n", "cm", "(points)", "per hour", "share", "share", "share",
         "MeV", "MeV", "per h", "per h", "per h");
  for (int it = 0; it < nT; ++it) {
    const Spectra& g = grid[it];
    if (g.tmc <= 0.) continue;
    TString checkText = "-";
    if (check[it].tmc > 0.) checkText = Form("%.1f", check[it].tmc);
    // from the 20 MeV bins, whose edges fall on 100, 300 and 1000 MeV
    auto below = [&](double e) { return g.topLin->Integral(1, g.topLin->FindBin(e * 0.999999)); };
    printf("%5g %10.1f +- %5.1f %12s %9.3f %9.3f %9.3f %11.0f %11.0f %9.2f %9.2f %9.1f\n", leads[it], g.tmc,
           std::sqrt(g.tmc2), checkText.Data(), g.withPb / g.tmc, g.withTop / g.tmc, g.withCmp / g.tmc,
           Quantile(g.top, 0.5), Quantile(g.cmp, 0.5), below(100.), below(300.), below(1000.));
  }
  printf("(share: fraction of the TMC rate whose muon has that energy recorded; top<E: rate of TMC muons\n"
         " entering the top TOF below E MeV, from 20 MeV bins)\n");

  // ---- CSV: rates per hour per bin
  const TString csvName = Form("%s/lead_after_spectrum%s.csv", dir, suffix.Data());
  if (FILE* csv = fopen(csvName, "w")) {
    fprintf(csv, "lead_cm,E_low_MeV,E_high_MeV,rate_sea_level_h,rate_after_lead_h,rate_top_tof_h,rate_compton_h,"
                 "err_sea_level_h,err_after_lead_h,err_top_tof_h,err_compton_h,check_rate_top_tof_h\n");
    for (int it = 0; it < nT; ++it)
      for (int b = 1; b < int(edges.size()); ++b) {
        const Spectra& g = grid[it];
        fprintf(csv, "%g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", leads[it], edges[b - 1], edges[b],
                g.ini->GetBinContent(b), g.pb->GetBinContent(b), g.top->GetBinContent(b), g.cmp->GetBinContent(b),
                g.ini->GetBinError(b), g.pb->GetBinError(b), g.top->GetBinError(b), g.cmp->GetBinError(b),
                check[it].top->GetBinContent(b));
      }
    fclose(csv);
    printf("rates per thickness and bin: %s\n", csvName.Data());
  }

  // ---- plots
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  gStyle->SetNumberContours(64);
  std::vector<int> shown;
  const double wanted[] = {0., 5., 10., 20., 30., 40., 50.};
  for (double t : wanted)
    for (int it = 0; it < nT; ++it)
      if (std::abs(leads[it] - t) < 1e-6 && grid[it].tmc > 0.) shown.push_back(it);
  if (shown.size() < 3) {   // not those thicknesses: about seven evenly spaced
    shown.clear();
    for (int it = 0; it < nT; it += std::max(1, nT / 6)) shown.push_back(it);
  }
  const int colours[] = {kBlack, kBlue + 1, kAzure + 7, kGreen + 2, kOrange + 7, kRed + 1, kMagenta + 2, kGray + 2};

  // one legend for the four panels, in a strip on the left, so that it covers no curve
  TCanvas c1("c1", "after the lead", 1750, 1100);
  auto legendPad = new TPad("legendPad", "", 0., 0., 0.15, 1.);
  auto plotPad = new TPad("plotPad", "", 0.15, 0., 1., 1.);
  legendPad->Draw();
  plotPad->Draw();
  plotPad->Divide(2, 2);
  auto sharedLegend = new TLegend(0.04, 0.42, 0.98, 0.9);
  sharedLegend->SetTextFont(43);   // size in pixels: a narrow pad would scale relative sizes by its width
  sharedLegend->SetTextSize(19);
  sharedLegend->SetBorderSize(0);
  sharedLegend->SetHeader("lead on top", "L");
  TH1D* legendSea = nullptr;
  TH1D* legendCheck = nullptr;
  struct PadSpec { const char* title; int which; };
  const PadSpec pads[] = {
    {"TMC muons leaving the lead layer;kinetic energy after the lead (MeV);dR/dE (h^{-1} GeV^{-1})", 1},
    {"TMC muons entering the top TOF;kinetic energy at the top TOF (MeV);dR/dE (h^{-1} GeV^{-1})", 2},
    {"TMC muons entering the Compton detector;kinetic energy at the Compton detector (MeV);dR/dE (h^{-1} GeV^{-1})", 3}};
  for (int p = 0; p < 3; ++p) {
    plotPad->cd(p + 1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetBottomMargin(0.12);
    double yMax = 0., yMin = 1e30;
    std::vector<TH1D*> hs;
    TH1D* seaLevelCurve = nullptr;
    for (std::size_t k = 0; k < shown.size(); ++k) {
      const Spectra& g = grid[shown[k]];
      const TH1D* src = pads[p].which == 1 ? g.pb : (pads[p].which == 2 ? g.top : g.cmp);
      TString label = Form("%g cm lead", leads[shown[k]]);
      const bool atSeaLevel = pads[p].which == 1 && leads[shown[k]] == 0.;
      if (atSeaLevel) {   // no lead layer: its place is reached with the sea-level energy
        src = g.ini;
        label = "0 cm (no lead): at sea level";
      }
      if (src->Integral() <= 0.) continue;
      TH1D* h = PerGeV(src, Form("perGeV_%d_%zu", p, k));
      h->SetLineColor(colours[k % 8]);
      h->SetLineWidth(2);
      for (int b = 1; b <= h->GetNbinsX(); ++b)
        if (h->GetBinContent(b) > 0.) {
          yMax = std::max(yMax, h->GetBinContent(b));
          yMin = std::min(yMin, h->GetBinContent(b));
        }
      if (atSeaLevel) seaLevelCurve = h;
      hs.push_back(h);
      if (pads[p].which == 2) sharedLegend->AddEntry(h, leads[shown[k]] == 0. ? TString("0 cm (no lead)") : label, "l");
    }
    if (hs.empty()) continue;
    // the simulated sea-level spectrum starts at the lowest band: sea-level curves are drawn from
    // their first full bin, without a drop to zero below it
    auto fromFirstBand = [&](TH1D* h) {
      int first = 1;
      while (first < h->GetNbinsX() && h->GetBinLowEdge(first) < 0.999 * bands.front().first) ++first;
      h->GetXaxis()->SetRange(first, h->GetNbinsX());
    };
    auto frame = gPad->DrawFrame(1., std::max(0.3 * yMin, 1e-6 * yMax), edges.back(), 3. * yMax);
    frame->SetTitle(pads[p].title);
    frame->GetXaxis()->SetTitleOffset(1.25);
    for (TH1D* h : hs) {
      if (h == seaLevelCurve) {
        fromFirstBand(h);
        h->Draw("HIST SAME ][");
      } else {
        h->Draw("HIST SAME");
      }
    }
    if (pads[p].which == 2) {   // sea level, and the check runs as open markers
      TH1D* sea = PerGeV(grid[shown[0]].ini, "perGeV_sea");
      sea->SetLineColor(kGray + 1);
      sea->SetLineStyle(2);
      sea->SetLineWidth(2);
      fromFirstBand(sea);
      sea->Draw("HIST SAME ][");
      legendSea = sea;
      for (std::size_t k = 0; k < shown.size(); ++k) {
        const Spectra& c = check[shown[k]];
        if (c.top->Integral() <= 0.) continue;
        TH1D* h = PerGeV(c.top, Form("perGeVcheck_%zu", k));
        h->SetMarkerStyle(24);
        h->SetMarkerSize(0.6);
        h->SetMarkerColor(colours[k % 8]);
        h->SetLineColor(colours[k % 8]);
        h->Draw("E1 SAME");
        if (!legendCheck) legendCheck = h;
      }
    }
  }
  legendPad->cd();
  if (legendSea) sharedLegend->AddEntry(legendSea, "sea level, same events", "l");
  if (legendCheck) sharedLegend->AddEntry(legendCheck, "check runs", "pe");
  sharedLegend->Draw();
  auto legendNote = new TLatex();
  legendNote->SetNDC();
  legendNote->SetTextFont(43);
  legendNote->SetTextSize(16);
  const char* noteLines[] = {"Dashed grey and circles:", "top-TOF panel only.", "Top left: the 0 cm line is",
                             "the energy at sea level.", "Check runs: the continuous", "spectrum, 2x10^{6} muons."};
  for (int k = 0; k < 6; ++k) legendNote->DrawLatex(0.06, 0.38 - 0.028 * k, noteLines[k]);
  if (muonOnly) {
    legendNote->SetTextColor(kRed + 1);
    const char* muonLines[] = {"Muon only: the muon itself", "crossed the top and middle", "TOF and entered the", "Compton detector."};
    for (int k = 0; k < 4; ++k) legendNote->DrawLatex(0.06, 0.18 - 0.028 * k, muonLines[k]);
  }
  // slow muons at the top TOF, linear scale
  plotPad->cd(4);
  gPad->SetGridy();
  gPad->SetBottomMargin(0.12);
  double linMax = 0.;
  std::vector<TH1D*> lin;
  for (std::size_t k = 0; k < shown.size(); ++k) {
    auto h = static_cast<TH1D*>(grid[shown[k]].topLin->Clone(Form("lin_%zu", k)));
    h->Scale(1. / 20.);   // per MeV (20 MeV bins)
    h->SetLineColor(colours[k % 8]);
    h->SetLineWidth(2);
    linMax = std::max(linMax, h->GetMaximum());
    lin.push_back(h);
  }
  if (!lin.empty()) {
    lin[0]->SetTitle("slow TMC muons entering the top TOF;kinetic energy at the top TOF (MeV);dR/dE (h^{-1} MeV^{-1})");
    lin[0]->SetMaximum(1.15 * linMax);
    lin[0]->SetMinimum(0.);
    lin[0]->Draw("HIST");
    for (std::size_t k = 1; k < lin.size(); ++k) lin[k]->Draw("HIST SAME");
  }

  // maps: rate per hour per bin versus energy and thickness, and the slow-muon rates
  TCanvas c2("c2", "maps", 1500, 1100);
  c2.Divide(2, 2);
  std::vector<double> tEdges(nT + 1);
  for (int it = 0; it <= nT; ++it)
    tEdges[it] = it == 0 ? leads[0] - 0.5 * (nT > 1 ? leads[1] - leads[0] : 1.)
                         : (it == nT ? leads[nT - 1] + 0.5 * (nT > 1 ? leads[nT - 1] - leads[nT - 2] : 1.)
                                     : 0.5 * (leads[it - 1] + leads[it]));
  const char* mapTitles[] = {"TMC muons leaving the lead layer, per hour and bin", "TMC muons entering the top TOF, per hour and bin",
                             "TMC muons entering the Compton detector, per hour and bin"};
  for (int m = 0; m < 3; ++m) {
    auto h = new TH2D(Form("map%d", m), Form("%s;kinetic energy (MeV);lead thickness (cm)", mapTitles[m]),
                      int(edges.size()) - 1, edges.data(), nT, tEdges.data());
    for (int it = 0; it < nT; ++it) {
      const TH1D* src = m == 0 ? grid[it].pb : (m == 1 ? grid[it].top : grid[it].cmp);
      for (int b = 1; b <= src->GetNbinsX(); ++b) h->SetBinContent(b, it + 1, src->GetBinContent(b));
    }
    c2.cd(m + 1);
    gPad->SetLogx();
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.12);
    h->GetXaxis()->SetTitleOffset(1.25);
    h->GetXaxis()->SetRangeUser(1., edges.back());
    h->SetMinimum(std::max(1e-4 * h->GetMaximum(), 1e-3));
    h->Draw("COLZ");
  }
  c2.cd(4);
  gPad->SetLogy();
  gPad->SetGridy();
  auto mg = new TMultiGraph("mgSlow", "TMC muons per hour versus lead;lead thickness (cm);rate (h^{-1})");
  auto legend = new TLegend(0.5, 0.14, 0.89, 0.36);
  const double limits[] = {1e30, 1000., 300., 100.};
  const char* names[] = {"all TMC", "top TOF below 1 GeV", "top TOF below 300 MeV", "top TOF below 100 MeV"};
  const int limitColours[] = {kBlack, kBlue + 1, kGreen + 2, kRed + 1};
  for (int l = 0; l < 4; ++l) {
    auto g = new TGraphErrors();
    for (int it = 0; it < nT; ++it) {
      const Spectra& s = grid[it];
      if (s.tmc <= 0.) continue;
      double r = s.tmc, e = std::sqrt(s.tmc2);
      if (l > 0) {   // the 20 MeV bins, edges on the limits
        const int last = s.topLin->FindBin(limits[l] * 0.999999);
        r = s.topLin->IntegralAndError(1, last, e);
      }
      if (r <= 0.) continue;
      g->SetPoint(g->GetN(), leads[it], r);
      g->SetPointError(g->GetN() - 1, 0., e);
    }
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.7);
    g->SetMarkerColor(limitColours[l]);
    g->SetLineColor(limitColours[l]);
    if (g->GetN() > 0) {
      mg->Add(g, "PL");
      legend->AddEntry(g, names[l], "pl");
    }
  }
  mg->Draw("A");
  legend->Draw();

  // one panel per thickness: the top-TOF spectrum from the grid, its check run and the 0 cm spectrum
  std::vector<int> panels;
  for (int it : shown)
    if (leads[it] > 0.) panels.push_back(it);
  int it0 = -1;   // 0 cm, the reference
  for (int it = 0; it < nT; ++it)
    if (leads[it] == 0. && grid[it].tmc > 0.) it0 = it;
  TCanvas c3("c3", "top TOF spectrum per thickness", 1650, 1000);
  c3.Divide(3, (int(panels.size()) + 2) / 3);
  double pMax = 0., pMin = 1e30;   // the same axes for every panel
  for (int it : panels) {
    const TH1D* hs3[] = {grid[it].top, check[it].top};
    for (const TH1D* h : hs3)
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double v = h->GetBinContent(b) / (h->GetBinWidth(b) / 1000.);
        if (v > 0.) { pMax = std::max(pMax, v); pMin = std::min(pMin, v); }
      }
  }
  for (std::size_t k = 0; k < panels.size(); ++k) {
    const int it = panels[k];
    c3.cd(int(k) + 1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetBottomMargin(0.12);
    gPad->SetLeftMargin(0.13);
    auto frame = gPad->DrawFrame(1., std::max(0.3 * pMin, 1e-6 * pMax), edges.back(), 3. * pMax);
    frame->SetTitle(Form("TMC muons entering the top TOF, %g cm lead;kinetic energy at the top TOF (MeV);dR/dE (h^{-1} GeV^{-1})", leads[it]));
    frame->GetXaxis()->SetTitleOffset(1.25);
    const int colour = colours[(std::find(shown.begin(), shown.end(), it) - shown.begin()) % 8];
    auto legend3 = new TLegend(0.15, 0.15, 0.7, 0.36);
    legend3->SetTextSize(0.038);
    if (it0 >= 0) {
      TH1D* ref = PerGeV(grid[it0].top, Form("ref_%zu", k));
      ref->SetLineColor(kBlack);
      ref->SetLineWidth(2);
      ref->SetLineStyle(2);
      ref->Draw("HIST SAME ][");
      legend3->AddEntry(ref, "0 cm lead (grid)", "l");
    }
    TH1D* line = PerGeV(grid[it].top, Form("grid_%zu", k));
    line->SetLineColor(colour);
    line->SetLineWidth(3);
    line->Draw("HIST SAME ][");
    legend3->AddEntry(line, Form("%g cm lead (60-point grid)", leads[it]), "l");
    if (check[it].top->Integral() > 0.) {
      TH1D* pts = PerGeV(check[it].top, Form("check_%zu", k));
      for (int b = 1; b <= pts->GetNbinsX(); ++b)   // empty bins are not points
        if (check[it].top->GetBinContent(b) <= 0.) pts->SetBinError(b, 0.);
      pts->SetMarkerStyle(20);
      pts->SetMarkerSize(0.7);
      pts->SetMarkerColor(colour);
      pts->SetLineColor(colour);
      pts->Draw("E1 X0 SAME");
      const double nCheck = check[it].tmc2 > 0. ? check[it].tmc * check[it].tmc / check[it].tmc2 : 0.;   // equal weights
      legend3->AddEntry(pts, Form("%g cm lead, check run (%.0f TMC events)", leads[it], nCheck), "pe");
    }
    legend3->Draw();
  }

  const TString base = Form("%s/lead_after_spectrum%s", dir, suffix.Data());
  if (muonOnly) {   // mark every page of the muon-only selection
    TCanvas* pages[] = {&c2, &c3};
    for (TCanvas* page : pages) {
      page->cd(0);
      auto note = new TLatex(0.995, 0.997, "muon only: the muon itself crossed the top and middle TOF and entered the Compton detector");
      note->SetNDC();
      note->SetTextAlign(33);
      note->SetTextSize(0.016);
      note->SetTextColor(kRed + 1);
      note->Draw();
    }
  }
  c1.Print(base + ".pdf(");
  c2.Print(base + ".pdf");
  c3.Print(base + ".pdf)");
  c1.SaveAs(base + "_spectra.png");
  c2.SaveAs(base + "_maps.png");
  c3.SaveAs(base + "_bythickness.png");
  printf("plots: %s.pdf, %s_spectra.png, %s_maps.png, %s_bythickness.png\n", base.Data(), base.Data(), base.Data(),
         base.Data());
}
