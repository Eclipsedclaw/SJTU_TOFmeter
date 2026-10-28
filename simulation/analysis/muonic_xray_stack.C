// Muonic X-rays as the Si/YSO stack sees them, with a detector response: the add-back spectrum of the hits of
// scripts/muonic_xray_stack.sh (macro/muonic_xray_stack.mac, columns stk_*).
//   root -l -b -q 'analysis/muonic_xray_stack.C("output/muonic_xray_stack")'
// Every hit is smeared in energy with sigma/E = 4% at 1 MeV scaling as 1/sqrt(E), i.e. sigma = 0.04 MeV x
// sqrt(E / MeV) (the same in the Si and the YSO layers), and dropped below the threshold (10 keV, after the
// smearing). Add-back: the hits of an event within the prompt window (100 ns) of the first one, in any layer, are
// added up into one entry, so an X-ray that Compton-scatters in one layer and is absorbed in another counts with
// its full energy; capture gammas and decay electrons (~1 us later) make entries of their own. Positions are not
// used: summing every prompt hit of the stack does not depend on where they are (the 2 mm position resolution
// matters for a Compton reconstruction).
// Arguments: directory, resolution at 1 MeV (sigma/E), threshold (keV), prompt window (ns), K_alpha of the target
// (keV: 354.9 for the first-principles cascade, the default; 334.7 for Geant4's; measured 346.8), range and bin
// width (keV), seed of the smearing.
// The spectrum is per day of sea-level exposure (the live time of the finished processes), split by what made the
// hits of an entry (stk_src): X-rays of the Al target (alone or with other energy), muonic X-rays of other volumes
// (plates, stack), and everything else (the cosmic background, capture gammas). In a +-2 sigma window around the
// Al lines: the entries with a target X-ray (S), the others (B), and S / sqrt(B) for the exposure.
// Writes <dir>/muonic_xray_stack_addback.png/.pdf and _spectrum.csv.

namespace {
// what made the hits of an add-back entry: 0 X-rays of the target only, 1 target X-rays with other energy,
// 2 muonic X-rays of other volumes, 3 everything else
int Part(int src)
{
  if (src == 1) return 0;
  if (src & 1) return 1;
  if (src & 2) return 2;
  return 3;
}

// L lines (n -> 2) of the hydrogen-like muonic atom with the reduced mass (keV), as both cascade models emit them
double LLine(int n, int Z, double A)
{
  const double muMass = 105.6583755, eMass = 0.51099895;
  const double nucleus = A * 931.494 - Z * eMass;
  const double e1 = 13.6e-3 * Z * Z * (muMass * nucleus / (muMass + nucleus)) / eMass;
  return e1 * (0.25 - 1. / (n * n));
}
}  // namespace

void muonic_xray_stack(const char* dir = "output/muonic_xray_stack", double resolution = 0.04, double thresholdKeV = 10.,
                       double windowNs = 100., double kAlphaKeV = 354.9, double xMaxKeV = 500., double binKeV = 2.,
                       unsigned seed = 4357)
{
  std::vector<std::string> files;
  const bool withLog = !gSystem->AccessPathName(Form("%s/log", dir));   // AccessPathName: true if it does NOT exist
  if (void* dirp = gSystem->OpenDirectory(dir)) {
    while (const char* entry = gSystem->GetDirEntry(dirp)) {
      const std::string name = entry;
      if (name.rfind("stk_", 0) != 0 || name.size() <= 9 || name.substr(name.size() - 5) != ".root") continue;
      // with the log folder of scripts/muonic_xray_stack.sh: only the finished processes, never a file still written
      const std::string k = name.substr(4, name.size() - 9);
      if (withLog && gSystem->AccessPathName(Form("%s/log/job%s.done", dir, k.c_str()))) continue;
      files.push_back(std::string(dir) + "/" + name);
    }
    gSystem->FreeDirectory(dirp);
  }
  if (files.empty()) {
    printf("no stk_*.root in %s (run scripts/muonic_xray_stack.sh)\n", dir);
    return;
  }

  const int nX = int(std::lround(xMaxKeV / binKeV));
  const char* partNames[] = {"X-rays of the Al target", "target X-rays with other energy", "muonic X-rays of other volumes",
                             "everything else (cosmic background)"};
  TH1D* parts[4];
  for (int p = 0; p < 4; ++p) parts[p] = new TH1D(Form("part%d", p), "", nX, 0., xMaxKeV);   // keV, counts
  TRandom3 random(seed);
  double liveTime = 0.;
  long long muons = 0, entries = 0;
  int used = 0;
  for (const auto& path : files) {
    TFile file(path.c_str());
    auto run = dynamic_cast<TTree*>(file.Get("run"));
    auto events = dynamic_cast<TTree*>(file.Get("events"));
    // a process still running has no run tree yet
    if (!run || !events || run->GetEntries() != 1 || !events->GetBranch("stk_e")) continue;
    int n = 0;
    double lt = 0.;
    run->SetBranchAddress("n_events", &n);
    run->SetBranchAddress("live_time", &lt);
    run->GetEntry(0);
    muons += n;
    liveTime += lt;
    ++used;
    std::vector<int>* src = nullptr;
    std::vector<double>*e = nullptr, *t = nullptr;
    events->SetBranchStatus("*", 0);   // read only the columns used (much faster)
    for (const char* column : {"stk_src", "stk_e", "stk_t"}) events->SetBranchStatus(column, 1);
    events->SetBranchAddress("stk_src", &src);
    events->SetBranchAddress("stk_e", &e);
    events->SetBranchAddress("stk_t", &t);
    std::vector<std::tuple<double, double, int>> hits;   // time (ns), measured energy (MeV), source bits
    for (Long64_t i = 0; i < events->GetEntries(); ++i) {
      events->GetEntry(i);
      hits.clear();
      for (std::size_t k = 0; k < e->size(); ++k) {
        const double measured = (*e)[k] + random.Gaus(0., resolution * std::sqrt((*e)[k]));   // sigma in MeV
        if (measured * 1000. >= thresholdKeV) hits.emplace_back((*t)[k], measured, (*src)[k]);
      }
      std::sort(hits.begin(), hits.end());
      for (std::size_t k = 0; k < hits.size();) {   // add-back: the hits within the window of the first one
        const double start = std::get<0>(hits[k]);
        double sum = 0.;
        int bits = 0;
        for (; k < hits.size() && std::get<0>(hits[k]) < start + windowNs; ++k) {
          sum += std::get<1>(hits[k]);
          bits |= std::get<2>(hits[k]);
        }
        parts[Part(bits)]->Fill(sum * 1000.);
        ++entries;
      }
    }
  }
  if (used == 0) {
    printf("no finished process in %s yet\n", dir);
    return;
  }
  const bool cosmic = liveTime > 0.;
  const double perDay = cosmic ? 86400. / liveTime : 1.;   // without live time (a gun run): counts

  // the Al lines, +-2 sigma windows: entries with a target X-ray (S) and the others (B), counts of the exposure
  struct Window { const char* name; double eKeV; };
  const std::vector<Window> windows = {{"L_alpha", LLine(3, 13, 26.9815)}, {"L_beta", LLine(4, 13, 26.9815)},
                                       {"K_alpha", kAlphaKeV}};
  printf("\n=== Add-back spectrum of the Si/YSO stack (%s) ===\n", dir);
  if (cosmic)
    printf("%d finished processes, %lld muons, %.2f days of sea-level exposure; %lld add-back entries\n", used, muons,
           liveTime / 86400., entries);
  else
    printf("%d files, %lld primaries (no live time: counts); %lld add-back entries\n", used, muons, entries);
  printf("response: sigma/E = %.1f%% at 1 MeV x sqrt(1 MeV / E), threshold %g keV per hit, window %g ns\n",
         100. * resolution, thresholdKeV, windowNs);
  printf("%-8s %9s %9s %12s %12s %10s %14s\n", "line", "keV", "+-2 sigma", "S (target)", "B (others)", "S/sqrt(B)",
         cosmic ? "S, B per day" : "");
  std::vector<std::string> windowText;
  for (const auto& w : windows) {
    const double sigma = resolution * std::sqrt(w.eKeV / 1000.) * 1000.;   // keV
    const int b1 = parts[0]->FindBin(w.eKeV - 2. * sigma), b2 = parts[0]->FindBin(w.eKeV + 2. * sigma);
    const double s = parts[0]->Integral(b1, b2) + parts[1]->Integral(b1, b2);
    const double b = parts[2]->Integral(b1, b2) + parts[3]->Integral(b1, b2);
    const double significance = b > 0. ? s / std::sqrt(b) : 0.;
    printf("%-8s %9.1f %9.1f %12.0f %12.0f %10.2f", w.name, w.eKeV, 2. * sigma, s, b, significance);
    if (cosmic) printf("   %6.2f, %.1f", s * perDay, b * perDay);
    printf("\n");
    windowText.push_back(Form("%s %.0f keV: S = %.0f, B = %.0f, S/#sqrt{B} = %.1f", w.name, w.eKeV, s, b, significance));
  }

  TString base = Form("%s/muonic_xray_stack_addback", dir);
  std::ofstream csv((base + "_spectrum.csv").Data());
  csv << "e_low_keV,e_high_keV,target_xray,target_xray_with_other,other_xray,background,total"
      << (cosmic ? "  (per day of exposure)" : "  (counts)") << "\n";
  for (int b = 1; b <= nX; ++b) {
    double total = 0.;
    csv << parts[0]->GetBinLowEdge(b) << "," << parts[0]->GetBinLowEdge(b + 1);
    for (int p = 0; p < 4; ++p) {
      csv << "," << parts[p]->GetBinContent(b) * perDay;
      total += parts[p]->GetBinContent(b);
    }
    csv << "," << total * perDay << "\n";
  }

  // ---- plot: per day per bin, log scale; the total and its parts
  auto total = static_cast<TH1D*>(parts[0]->Clone("total"));
  for (int p = 1; p < 4; ++p) total->Add(parts[p]);
  auto targetAll = static_cast<TH1D*>(parts[0]->Clone("targetAll"));
  targetAll->Add(parts[1]);
  for (auto h : {total, targetAll, parts[2], parts[3]}) h->Scale(perDay);
  const int blue = TColor::GetColor("#2a78d6"), orange = TColor::GetColor("#eb6834"), aqua = TColor::GetColor("#1baf7a");
  const int grey = TColor::GetColor("#8c8a83"), inkSecondary = TColor::GetColor("#52514e");
  gStyle->SetOptStat(0);
  TCanvas c("c", "stack add-back", 1400, 850);
  auto pPlot = new TPad("pPlot", "", 0., 0.18, 1., 1.);
  auto pT = new TPad("pT", "", 0., 0., 1., 0.18);
  for (auto pad : {pPlot, pT}) pad->Draw();
  pPlot->cd();
  gPad->SetLogy();
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.03);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.11);
  double yMin = 1e30;
  for (auto h : {total, targetAll, parts[2]})
    for (int b = 1; b <= nX; ++b)
      if (h->GetBinContent(b) > 0.) yMin = std::min(yMin, h->GetBinContent(b));
  total->SetTitle(cosmic ? Form(";add-back energy (keV);entries per day per %g keV", binKeV)
                         : Form(";add-back energy (keV);entries per %g keV", binKeV));
  total->SetMinimum(yMin < 1e30 ? 0.5 * yMin : 1e-3);
  const double dataMax = total->GetMaximum();
  total->SetMaximum(dataMax * 300.);   // above the data: the line labels (~2 decades) and the legend
  for (auto axis : {total->GetXaxis(), total->GetYaxis()}) {
    axis->SetTitleSize(0.04);
    axis->SetLabelSize(0.036);
    axis->SetTitleColor(inkSecondary);
  }
  total->GetYaxis()->SetTitleOffset(0.9);
  total->SetLineColor(blue);
  total->SetLineWidth(3);
  parts[3]->SetLineColor(grey);
  parts[2]->SetLineColor(aqua);
  targetAll->SetLineColor(orange);
  for (auto h : {parts[3], parts[2], targetAll}) h->SetLineWidth(2);
  parts[3]->SetLineStyle(2);   // dashed: the background is nearly the whole total, so it lies on top of it
  total->Draw("HIST");
  parts[3]->Draw("HIST SAME");
  parts[2]->Draw("HIST SAME");
  targetAll->Draw("HIST SAME");
  // the Al lines: dashed markers, labels in the ink colour
  auto label = new TLatex();
  label->SetTextFont(42);
  label->SetTextSize(0.032);
  label->SetTextColor(inkSecondary);
  label->SetTextAlign(12);
  label->SetTextAngle(90);
  const char* latexNames[] = {"L_{#alpha}", "L_{#beta}", "K_{#alpha}"};
  for (std::size_t w = 0; w < windows.size(); ++w) {
    if (windows[w].eKeV > xMaxKeV) continue;
    // up to the space above the data (log scale: the frame goes to 30 x the highest bin)
    auto marker = new TLine(windows[w].eKeV, total->GetMinimum(), windows[w].eKeV, dataMax * 2.);
    marker->SetLineStyle(2);
    marker->SetLineColor(inkSecondary);
    marker->Draw();
    label->DrawLatex(windows[w].eKeV + 0.006 * xMaxKeV, dataMax * 2.5, Form("Al %s %.1f", latexNames[w], windows[w].eKeV));
  }
  auto title = new TLatex();
  title->SetNDC();
  title->SetTextFont(62);
  title->SetTextSize(0.04);
  title->DrawLatex(0.08, 0.935, "Add-back spectrum of the Si/YSO stack, sea-level muons through the setup with its Al target");
  auto legend = new TLegend(0.29, 0.7, 0.67, 0.88);   // between the L_beta and K_alpha markers
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.032);
  legend->AddEntry(total, "all entries (what the detector sees)", "l");
  legend->AddEntry(targetAll, "with an X-ray of the Al target", "l");
  legend->AddEntry(parts[2], partNames[2], "l");
  legend->AddEntry(parts[3], partNames[3], "l");
  legend->Draw();

  pT->cd();
  auto box = new TPaveText(0.01, 0.05, 0.99, 0.95, "NDC");
  box->SetFillColor(kWhite);
  box->SetBorderSize(1);
  box->SetTextAlign(12);
  box->SetTextFont(42);
  box->SetTextSize(0.13);
  box->SetTextColor(inkSecondary);
  box->AddText(cosmic ? Form("%.2f days of sea-level exposure (%lld muons, %d processes); %lld add-back entries.", liveTime / 86400.,
                             muons, used, entries)
                      : Form("%lld primaries (no live time: counts); %lld add-back entries.", muons, entries));
  box->AddText(Form("Response: sigma/E = %.0f%% at 1 MeV, scaling as 1/sqrt(E), in every layer; hits above %g keV; prompt window %g ns.",
                    100. * resolution, thresholdKeV, windowNs));
  box->AddText(("Counts of the exposure in +-2 sigma: " + windowText[0] + ";  " + windowText[2]).c_str());
  box->AddText(Form("Lines of the simulated cascade (K_{#alpha} %.1f keV; measured 346.8 keV). No TOF selection: every muon counts.",
                    kAlphaKeV));
  box->Draw();
  c.Print(base + ".pdf");
  c.SaveAs(base + ".png");
  printf("plots: %s.png, .pdf; spectrum: %s_spectrum.csv\n", base.Data(), base.Data());
}
