// Dummy muonic X-ray spectra of several materials (or cascade models) together: the emitted X-ray yield
// (muonic X-rays / stopped mu-) of mu- stopped in the target, from scripts/dummy_xray.sh
// (macro/dummy_xray.mac; one run per material and cascade model, in <dir>/<material>_<model>/).
//   root -l -b -q 'analysis/dummy_xray_compare.C("output/dummy_xray")'
// Arguments: directory, X-ray range and bin width (default 0-120 keV in 1 keV bins), y-axis top (1). Every folder
// <material>_<model> of the directory with mux*.root files is a run; for an element (G4_Al, G4_Si, ...)
// ROOT's element table gives Z and A, for the names of the lines. Selection: primary mu- that stopped
// and was captured at rest (end_proc muMinusCaptureAtRest) in the Target, and its muonic X-rays
// (mux_src 0, mux_parent 1) at their energy of creation: an ideal detector, no absorption, no
// resolution. Each photon is given to the line n -> n' of the hydrogen-like muonic atom (reduced mass)
// nearest to its energy, so the isotopes of an element (lines a few 0.01 keV apart) make one line. A
// cascade makes each transition, and reaches each level, at most once, so the uncertainty of a yield y
// of a line or a series from N stops is binomial, sqrt(y (1 - y) / N), and a yield is at most 1: a line or
// series above 1 is reported as an error. The y axis goes from 0 to 1, or to the last argument: a bin above
// it is drawn cut, with a triangle at the top, and its value is given under the plot.
// Writes <dir>/dummy_xray_compare.png/.pdf, _lines.csv (energy and yield of every line in every run)
// and _bins.csv (the histograms).

namespace {
using LineList = std::vector<std::tuple<double, int, int>>;   // energy (keV), n, n'

struct Run {
  std::string name, material, model, symbol, tag, label;   // e.g. G4_Si_firstPrinciples, G4_Si, firstPrinciples, Si
  int Z = 0;
  double A = 0.;
  LineList levels;
  std::vector<std::string> files;
  TH1D* hist = nullptr;
  int colour = 0, marker = 0;
  Long64_t generated = 0;
  long stops = 0, xrays = 0, xraysInRange = 0, unassigned = 0;
  std::map<std::pair<int, int>, std::pair<double, double>> lines;   // (n, n') -> photons, sum of their energies (keV)
  std::map<int, double> series;                                     // n' -> photons
};

double YieldError(double count, long stops)
{
  const double y = count / stops;
  return std::sqrt(std::max(0., y * (1. - y)) / stops);
}

// The lines n -> n' (n <= 14) of the hydrogen-like muonic atom with the reduced mass
LineList HydrogenLikeLines(int Z, double A)
{
  const double muMass = 105.6583755, eMass = 0.51099895;
  const double nucleus = A * 931.494 - Z * eMass;   // close enough for the reduced mass
  const double e1 = 13.6e-3 * Z * Z * (muMass * nucleus / (muMass + nucleus)) / eMass;   // keV
  LineList lines;
  for (int n = 2; n <= 14; ++n)
    for (int n2 = 1; n2 < n; ++n2) lines.emplace_back(e1 * (1. / (n2 * n2) - 1. / (n * n)), n, n2);
  return lines;
}

// The line nearest to eKeV, if within 0.06 keV, else {0, 0}
std::pair<int, int> Transition(double eKeV, const LineList& lines)
{
  std::pair<int, int> best = {0, 0};
  double bestDiff = 0.06;
  for (const auto& [e, n, n2] : lines)
    if (std::abs(e - eKeV) < bestDiff) {
      bestDiff = std::abs(e - eKeV);
      best = {n, n2};
    }
  return best;
}

// e.g. L_alpha for 3 -> 2 (ROOT TLatex with latex = true), n->n' beyond epsilon or Q
std::string LineName(int n, int n2, bool latex)
{
  if (n == 0) return "?";
  const char* series = "KLMNOPQ";
  const char* greek[] = {"alpha", "beta", "gamma", "delta", "epsilon"};
  if (n2 <= 7 && n - n2 <= 5) return Form(latex ? "%c_{#%s}" : "%c_%s", series[n2 - 1], greek[n - n2 - 1]);
  return Form(latex ? "%d#rightarrow%d" : "%d->%d", n, n2);
}

std::string ModelLabel(const std::string& model)
{
  if (model == "firstPrinciples") return "first-principles";
  if (model == "geant4") return "Geant4";
  return model;
}
}  // namespace

void dummy_xray_compare(const char* dir = "output/dummy_xray", double xMaxKeV = 120., double binKeV = 1.,
                        double yAxisMax = 1.)
{
  // ---- the runs: every folder <material>_<model> of dir with mux*.root files
  std::vector<Run> runs;
  TGeoElementTable elements(200);
  elements.BuildDefaultElements();
  if (void* dirp = gSystem->OpenDirectory(dir)) {
    while (const char* entry = gSystem->GetDirEntry(dirp)) {
      Run r;
      r.name = entry;
      if (r.name == "." || r.name == "..") continue;
      const std::string sub = std::string(dir) + "/" + r.name;
      if (void* subp = gSystem->OpenDirectory(sub.c_str())) {
        while (const char* file = gSystem->GetDirEntry(subp)) {
          const std::string f = file;
          if (f.rfind("mux", 0) == 0 && f.size() > 5 && f.substr(f.size() - 5) == ".root") r.files.push_back(sub + "/" + f);
        }
        gSystem->FreeDirectory(subp);
      }
      if (r.files.empty()) continue;
      const auto cut = r.name.rfind('_');
      r.material = cut == std::string::npos ? r.name : r.name.substr(0, cut);
      r.model = cut == std::string::npos ? "" : r.name.substr(cut + 1);
      r.symbol = r.material.rfind("G4_", 0) == 0 ? r.material.substr(3) : r.material;
      if (auto element = elements.FindElement(r.symbol.c_str())) {   // not for compounds: their lines stay unnamed
        r.Z = element->Z();
        r.A = element->A();
        r.levels = HydrogenLikeLines(r.Z, r.A);
      }
      runs.push_back(r);
    }
    gSystem->FreeDirectory(dirp);
  }
  if (runs.empty()) {
    printf("no <material>_<model>/mux*.root in %s (run scripts/dummy_xray.sh)\n", dir);
    return;
  }
  const std::vector<std::string> palette = {"#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"};
  const int markers[] = {20, 25, 22, 27, 21, 24, 23};
  if (runs.size() > palette.size()) {
    printf("%zu runs in %s: at most %zu fit in one plot\n", runs.size(), dir, palette.size());
    return;
  }
  std::sort(runs.begin(), runs.end(), [](const Run& a, const Run& b) { return a.Z != b.Z ? a.Z < b.Z : a.name < b.name; });
  std::vector<std::string> symbols, models;   // in the order of the runs
  for (const auto& r : runs) {
    if (std::find(symbols.begin(), symbols.end(), r.symbol) == symbols.end()) symbols.push_back(r.symbol);
    if (std::find(models.begin(), models.end(), r.model) == models.end()) models.push_back(r.model);
  }
  const int nX = int(std::lround(xMaxKeV / binKeV));
  for (std::size_t i = 0; i < runs.size(); ++i) {
    auto& r = runs[i];
    // what tells the runs apart: the element, the model, or both
    if (models.size() == 1) r.tag = r.symbol;
    else if (symbols.size() == 1) r.tag = ModelLabel(r.model);
    else r.tag = r.symbol + " " + ModelLabel(r.model);
    r.label = models.size() == 1 && r.Z > 0 ? r.tag + Form(" (Z = %d)", r.Z) : r.tag;
    r.colour = TColor::GetColor(palette[i].c_str());   // fixed order: the same run keeps its colour
    r.marker = markers[i];
    r.hist = new TH1D(("h_" + r.name).c_str(), "", nX, 0., xMaxKeV);   // keV
  }

  for (auto& r : runs) {
    for (const auto& path : r.files) {
      TFile file(path.c_str());
      auto run = dynamic_cast<TTree*>(file.Get("run"));
      auto events = dynamic_cast<TTree*>(file.Get("events"));
      if (!run || !events || run->GetEntries() != 1 || !events->GetBranch("mux_e")) {
        printf("skipping %s: incomplete or without the mux_* columns\n", path.c_str());
        continue;
      }
      int n = 0;
      run->SetBranchAddress("n_events", &n);
      run->GetEntry(0);
      r.generated += n;
      int pdg = 0;
      char endProc[128] = "", endVol[128] = "";
      std::vector<double>* e = nullptr;
      std::vector<int>*src = nullptr, *parent = nullptr;
      events->SetBranchAddress("pdg", &pdg);
      events->SetBranchAddress("end_proc", endProc);
      events->SetBranchAddress("end_vol", endVol);
      events->SetBranchAddress("mux_e", &e);
      events->SetBranchAddress("mux_src", &src);
      events->SetBranchAddress("mux_parent", &parent);
      for (Long64_t i = 0; i < events->GetEntries(); ++i) {
        events->GetEntry(i);
        if (pdg != 13 || std::string(endProc) != "muMinusCaptureAtRest" || std::string(endVol).rfind("Target", 0) != 0)
          continue;
        ++r.stops;
        for (std::size_t k = 0; k < e->size(); ++k) {
          if ((*parent)[k] != 1 || (*src)[k] != 0) continue;   // muonic X-rays of the primary only
          const double eKeV = (*e)[k] * 1000.;
          ++r.xrays;
          if (eKeV >= xMaxKeV) continue;
          ++r.xraysInRange;
          const auto t = Transition(eKeV, r.levels);
          if (t.first == 0) {
            ++r.unassigned;
            r.hist->Fill(eKeV);   // the others are filled below, at the energy of their line
            continue;
          }
          r.lines[t].first += 1.;
          r.lines[t].second += eKeV;
          r.series[t.second] += 1.;
        }
      }
    }
    if (r.stops == 0) {
      printf("no mu- stopped in the Target in %s/%s\n", dir, r.name.c_str());
      return;
    }
    // every line at the mean energy of its photons: the isotopes of an element (a few 0.01 keV apart)
    // never fall into two bins
    for (const auto& [t, sums] : r.lines) r.hist->Fill(sums.second / sums.first, sums.first);
    // per stopped mu-: the yield in each bin, with its binomial uncertainty
    for (int b = 1; b <= nX; ++b) {
      const double count = r.hist->GetBinContent(b);
      r.hist->SetBinContent(b, count / r.stops);
      r.hist->SetBinError(b, YieldError(count, r.stops));
    }
  }

  // ---- tables
  printf("\n=== Dummy muonic X-ray yields (%s), mu- stopped in the target, ideal detector ===\n", dir);
  printf("%-28s %4s %10s %10s %13s %15s\n", "run", "Z", "muons", "stopped", "X-rays/stop", Form("below %g keV", xMaxKeV));
  for (const auto& r : runs) {
    printf("%-28s %4d %10lld %10ld %13.3f %15.3f\n", r.name.c_str(), r.Z, r.generated, r.stops, double(r.xrays) / r.stops,
           double(r.xraysInRange) / r.stops);
    if (r.unassigned > 0)
      printf("   %.1f%% of its X-rays below %g keV match no hydrogen-like line: not in the line tables\n",
             100. * r.unassigned / r.xraysInRange, xMaxKeV);
  }
  std::string header;
  for (const auto& r : runs) header += Form(" %27s", r.tag.c_str());
  printf("\nseries below %g keV (lines n -> n' to one level n'; at most one per stop):\n%-10s%s\n", xMaxKeV, "series",
         header.c_str());
  std::map<int, bool> levels;
  for (const auto& r : runs)
    for (const auto& [n2, count] : r.series) levels[n2] = true;
  for (const auto& [n2, unused] : levels) {
    std::string row = n2 <= 7 ? Form("%c (n' = %d)", "KLMNOPQ"[n2 - 1], n2) : Form("n' = %d", n2);
    row.resize(10, ' ');
    for (const auto& r : runs) {
      const double c = r.series.count(n2) ? r.series.at(n2) : 0.;
      row += Form(" %18.4f +- %.4f", c / r.stops, YieldError(c, r.stops));
    }
    printf("%s\n", row.c_str());
  }

  std::map<std::pair<int, int>, bool> found;   // every line of any run
  for (const auto& r : runs)
    for (const auto& [t, sums] : r.lines) found[t] = true;
  std::vector<std::pair<int, int>> transitions;
  for (const auto& [t, unused] : found) transitions.push_back(t);
  std::sort(transitions.begin(), transitions.end(),   // by series (n'), then n
            [](const auto& a, const auto& b) { return a.second != b.second ? a.second < b.second : a.first < b.first; });
  TString base = Form("%s/dummy_xray_compare", dir);
  std::ofstream csv((base + "_lines.csv").Data());
  csv << "line,transition";
  for (const auto& r : runs) csv << ",energy_keV_" << r.name << ",yield_" << r.name << ",error_" << r.name;
  csv << "\n";
  header.clear();
  for (const auto& r : runs) header += Form(" %27s", (r.tag + ": keV, yield").c_str());
  printf("\nlines below %g keV with a yield of at least 1%% in one run (every line in the CSV):\n%-10s %7s%s\n", xMaxKeV,
         "name", "n->n'", header.c_str());
  struct Candidate { double eKeV, yield; std::string text; };
  std::vector<Candidate> candidates;   // for the labels of the plot
  for (const auto& t : transitions) {
    std::string row = Form("%-10s %7s", LineName(t.first, t.second, false).c_str(), Form("%d->%d", t.first, t.second));
    csv << LineName(t.first, t.second, false) << "," << t.first << "->" << t.second;
    double most = 0.;
    for (const auto& r : runs) {
      const auto it = r.lines.find(t);
      const double c = it == r.lines.end() ? 0. : it->second.first, y = c / r.stops;
      const double eKeV = c > 0. ? it->second.second / c : 0.;   // mean energy of its photons (isotopes)
      most = std::max(most, y);
      row += c > 0. ? Form(" %8.2f %8.4f +- %.4f", eKeV, y, YieldError(c, r.stops)) : Form(" %8s %8s    %6s", "-", "0", "");
      csv << "," << (c > 0. ? Form("%.2f", eKeV) : "") << "," << y << "," << YieldError(c, r.stops);
      // labels for the named lines (alpha to epsilon) only; below 5 keV they would be too crowded
      if (eKeV >= 5. && y >= 0.02 && t.first - t.second <= 5)
        candidates.push_back({eKeV, y, (symbols.size() > 1 ? r.symbol + " " : std::string()) + LineName(t.first, t.second, true) +
                                         Form(" %.1f", eKeV)});
    }
    csv << "\n";
    if (most >= 0.01) printf("%s\n", row.c_str());
  }
  std::ofstream bins((base + "_bins.csv").Data());
  bins << "e_low_keV,e_high_keV";
  for (const auto& r : runs) bins << ",yield_" << r.name << ",error_" << r.name;
  bins << "\n";
  for (int b = 1; b <= nX; ++b) {
    bins << runs[0].hist->GetBinLowEdge(b) << "," << runs[0].hist->GetBinLowEdge(b + 1);
    for (const auto& r : runs) bins << "," << r.hist->GetBinContent(b) << "," << r.hist->GetBinError(b);
    bins << "\n";
  }

  // a yield is at most 1: each transition happens, and each level is reached, at most once per stop
  for (const auto& r : runs) {
    for (const auto& [t, sums] : r.lines)
      if (sums.first > r.stops)
        printf("ERROR %s: line %d->%d has %g photons for %ld stops, a yield above 1, which is impossible: check the cascade\n",
               r.name.c_str(), t.first, t.second, sums.first, r.stops);
    for (const auto& [n2, count] : r.series)
      if (count > r.stops)
        printf("ERROR %s: the series to n' = %d has %g photons for %ld stops, a yield above 1, which is impossible\n",
               r.name.c_str(), n2, count, r.stops);
    for (int b = 1; b <= nX; ++b)
      if (r.hist->GetBinContent(b) > yAxisMax)
        printf("WARNING %s: %g-%g keV holds %.3f X-rays per stop%s, above the top of the y axis (%g): drawn cut there\n",
               r.name.c_str(), r.hist->GetBinLowEdge(b), r.hist->GetBinLowEdge(b + 1), r.hist->GetBinContent(b),
               r.hist->GetBinContent(b) > 1. ? " (several lines in one bin)" : "", yAxisMax);
  }

  // ---- plot: the yield per bin
  const int inkSecondary = TColor::GetColor("#52514e");
  gStyle->SetOptStat(0);
  TCanvas c("c", "dummy muonic X-rays", 1400, 850);
  auto pLin = new TPad("pLin", "", 0., 0.18, 1., 1.);
  auto pT = new TPad("pT", "", 0., 0., 1., 0.18);
  for (auto pad : {pLin, pT}) pad->Draw();
  double yMax = 0.;
  for (std::size_t i = 0; i < runs.size(); ++i) {
    auto& r = runs[i];
    yMax = std::max(yMax, r.hist->GetMaximum());
    r.hist->SetLineColor(r.colour);
    r.hist->SetLineWidth(i == 0 ? 2 : 3);   // the others are drawn first and wider: they show round the first run
    r.hist->SetMarkerColor(r.colour);       // for the legend; the markers are drawn by a TGraph
    r.hist->SetMarkerStyle(r.marker);
    r.hist->SetMarkerSize(1.1);
  }
  const double yTop = yAxisMax;   // 1 unless asked otherwise: a yield is at most 1 per stopped muon
  const double left = 0.08, right = 0.03, top = 0.1, bottom = 0.11;   // margins of the plot pad
  pLin->cd();
  gPad->SetLeftMargin(left);
  gPad->SetRightMargin(right);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  auto frame = new TH1D("frame", "", nX, 0., xMaxKeV);
  // "muon" in words: ROOT 6.32 leaves a gap before a Greek letter that follows other text
  const TString perBin = std::abs(binKeV - 1.) < 1e-9 ? TString("keV") : TString(Form("%g keV", binKeV));
  frame->SetTitle(Form(";photon energy at creation (keV);X-ray yield / stopped muon (per %s)", perBin.Data()));
  frame->SetMinimum(0.);
  frame->SetMaximum(yTop);
  for (auto axis : {frame->GetXaxis(), frame->GetYaxis()}) {
    axis->SetTitleSize(0.04);
    axis->SetLabelSize(0.036);
    axis->SetTitleColor(inkSecondary);
  }
  frame->GetYaxis()->SetTitleOffset(0.9);
  frame->Draw("AXIS");
  for (int i = int(runs.size()) - 1; i >= 0; --i) {
    runs[i].hist->Draw("HIST SAME");
    auto tops = new TGraph();   // a marker on top of each line, but not on the tiny ones
    auto cut = new TGraph();    // a triangle at the top of the axis for a bin above it
    for (int b = 1; b <= nX; ++b) {
      const double y = runs[i].hist->GetBinContent(b);
      if (y > yTop) cut->SetPoint(cut->GetN(), runs[i].hist->GetBinCenter(b), 0.97 * yTop);
      else if (y >= 0.01 * yMax) tops->SetPoint(tops->GetN(), runs[i].hist->GetBinCenter(b), y);
    }
    tops->SetMarkerStyle(runs[i].marker);
    tops->SetMarkerColor(runs[i].colour);
    tops->SetMarkerSize(1.1);
    if (tops->GetN() > 0) tops->Draw("P");
    cut->SetMarkerStyle(22);
    cut->SetMarkerColor(runs[i].colour);
    cut->SetMarkerSize(1.6);
    if (cut->GetN() > 0) cut->Draw("P");
  }
  // the legend where no line (or its label) reaches up to it: top left, top right, lower right, top middle
  auto keV = [&](double xNdc) { return (xNdc - left) / (1. - left - right) * xMaxKeV; };   // NDC of the pad -> axes
  auto yieldAt = [&](double yNdc) { return (yNdc - bottom) / (1. - top - bottom) * yTop; };
  const double legendWidth = 0.28, legendHeight = 0.047 * runs.size() + 0.024;
  const std::vector<std::pair<double, double>> corners = {{0.12, 0.86}, {0.67, 0.86}, {0.67, 0.62}, {0.4, 0.86}};
  double legendX = 0.67, legendTop = 0.86;
  for (const auto& corner : corners) {
    const int b1 = runs[0].hist->FindBin(std::max(0., keV(corner.first) - 1.));
    const int b2 = std::min(nX, runs[0].hist->FindBin(keV(corner.first + legendWidth) + 1.));
    bool free = true;
    for (const auto& r : runs)
      for (int b = b1; b <= b2; ++b)
        if (r.hist->GetBinContent(b) + 0.1 * yTop > yieldAt(corner.second - legendHeight)) free = false;
    if (free) {
      legendX = corner.first;
      legendTop = corner.second;
      break;
    }
  }
  auto legend = new TLegend(legendX, legendTop - legendHeight, legendX + legendWidth, legendTop);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.032);
  legend->SetMargin(0.15);
  for (const auto& r : runs) legend->AddEntry(r.hist, r.label.c_str(), "lp");
  legend->Draw();
  // labels of the strongest lines, in the ink colour (the material is in the text), each moved up until
  // it is clear of the other lines and labels
  pLin->Update();   // the axes of the pad, for the sizes of the labels
  struct Box { double x1, x2, y1, y2; };
  std::vector<Box> taken = {{keV(legendX), keV(legendX + legendWidth), yieldAt(legendTop - legendHeight), yieldAt(legendTop)}};
  std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) { return a.yield > b.yield; });
  auto label = new TLatex();
  label->SetTextFont(42);
  label->SetTextSize(0.03);
  label->SetTextColor(inkSecondary);
  std::vector<std::string> placed;   // a line of the same element in two models is labelled once
  for (const auto& cand : candidates) {
    if (placed.size() >= 10u) break;
    if (std::find(placed.begin(), placed.end(), cand.text) != placed.end()) continue;
    TLatex probe(0., 0., cand.text.c_str());
    probe.SetTextFont(42);
    probe.SetTextSize(0.03);
    const double w = probe.GetXsize(), h = probe.GetYsize();
    Box box = {cand.eKeV - w / 2., cand.eKeV + w / 2., cand.yield + 0.03 * yTop, 0.};
    int align = 21;   // centred on the line, but kept inside the frame near its ends
    if (box.x1 < 0.) { box.x1 = cand.eKeV; box.x2 = cand.eKeV + w; align = 11; }
    if (box.x2 > xMaxKeV) { box.x1 = cand.eKeV - w; box.x2 = cand.eKeV; align = 31; }
    const int b1 = std::max(1, runs[0].hist->FindBin(box.x1)), b2 = std::min(nX, runs[0].hist->FindBin(box.x2));
    for (int step = 0; step < 10; ++step) {   // up above every line under the label, then every label in the way
      double y1 = box.y1;
      for (const auto& r : runs)
        for (int b = b1; b <= b2; ++b) y1 = std::max(y1, r.hist->GetBinContent(b) + 0.03 * yTop);
      for (const auto& t : taken)
        if (box.x1 < t.x2 && t.x1 < box.x2 && y1 < t.y2 && t.y1 < y1 + h) y1 = t.y2 + 0.01 * yTop;
      if (y1 == box.y1) break;
      box.y1 = y1;
    }
    box.y2 = box.y1 + h;
    if (box.y2 > yTop) continue;   // no room left
    taken.push_back(box);
    label->SetTextAlign(align);
    label->DrawLatex(cand.eKeV, box.y1, cand.text.c_str());
    placed.push_back(cand.text);
  }
  std::string what;   // "Al and Si"
  for (std::size_t i = 0; i < symbols.size(); ++i) what += (i == 0 ? "" : i + 1 == symbols.size() ? " and " : ", ") + symbols[i];
  std::string modelNames;
  for (std::size_t i = 0; i < models.size(); ++i) modelNames += (i == 0 ? "" : " and ") + ModelLabel(models[i]);
  auto title = new TLatex();
  title->SetNDC();
  title->SetTextFont(62);
  title->SetTextSize(0.04);
  if (models.size() == 1) title->DrawLatex(0.08, 0.935, Form("Muonic X-ray yield of %s, 0-%g keV (%s cascade)", what.c_str(), xMaxKeV, modelNames.c_str()));
  else if (symbols.size() == 1) title->DrawLatex(0.08, 0.935, Form("Muonic X-ray yield of %s, 0-%g keV: %s cascades", what.c_str(), xMaxKeV, modelNames.c_str()));
  else title->DrawLatex(0.08, 0.935, Form("Muonic X-ray yield, 0-%g keV", xMaxKeV));

  pT->cd();
  auto box = new TPaveText(0.01, 0.05, 0.99, 0.95, "NDC");
  box->SetFillColor(kWhite);
  box->SetBorderSize(1);
  box->SetTextAlign(12);
  box->SetTextFont(42);
  box->SetTextColor(inkSecondary);
  std::string stopped, perStop, offScale;
  for (const auto& r : runs) {
    stopped += (stopped.empty() ? "" : ", ") + r.tag + Form(" %ld", r.stops);
    perStop += (perStop.empty() ? "" : ", ") + r.tag + Form(" %.3f", double(r.xraysInRange) / r.stops);
  }
  for (int b = 1; b <= nX; ++b) {   // the bins above the top of the axis, with their values
    std::string values;
    for (const auto& r : runs)
      if (r.hist->GetBinContent(b) > yTop) values += (values.empty() ? "" : ", ") + r.tag + Form(" %.3f", r.hist->GetBinContent(b));
    if (!values.empty())
      offScale += (offScale.empty() ? "" : "; ") + std::string(Form("%g-%g keV ", runs[0].hist->GetBinLowEdge(b),
                                                                       runs[0].hist->GetBinLowEdge(b + 1))) + values;
  }
  std::vector<std::string> notes;
  notes.push_back(Form("Muons stopped in the target: %s. X-rays per stop below %g keV: %s.", stopped.c_str(), xMaxKeV, perStop.c_str()));
  notes.push_back(Form("Ideal detector: energies at creation in %g keV bins, so a line's height is its yield per stopped muon"
                       " (lines closer than a bin add up).", binKeV));
  if (!offScale.empty())
    notes.push_back(Form("Above the top of the axis (drawn cut, triangle): %s.", offScale.c_str()));
  bool noK = true;   // the K lines (n -> 1) are where the two cascade models differ in energy
  for (const auto& r : runs)
    for (const auto& [eKeV, n, n2] : r.levels)
      if (n2 == 1 && eKeV < xMaxKeV) noK = false;
  // plain text: ROOT 6.32 spaces TLatex symbols and sub/superscripts after other text badly
  if (noK && symbols.size() == 1 && models.size() > 1)
    notes.push_back("The models put these lines at the same hydrogen-like energies (reduced mass); they differ only in the yields.");
  else if (noK && symbols.size() == 2 && runs.front().Z > 0 && runs.back().Z > 0)
    notes.push_back(Form("Line energies: hydrogen-like levels with the reduced mass, so they scale with Z squared (%s/%s: %d/%d = %.3f).",
                         symbols[1].c_str(), symbols[0].c_str(), runs.back().Z * runs.back().Z, runs.front().Z * runs.front().Z,
                         double(runs.back().Z * runs.back().Z) / (runs.front().Z * runs.front().Z)));
  else if (noK)
    notes.push_back("Line energies: hydrogen-like levels with the reduced mass, so they scale with Z squared.");
  notes.push_back("Lines below ~10 keV (unlabelled below 5 keV) would be absorbed in the target where they are made.");
  box->SetTextSize(notes.size() > 4 ? 0.11 : 0.13);   // five lines need smaller text to fit
  for (const auto& note : notes) box->AddText(note.c_str());
  box->Draw();
  c.Print(base + ".pdf");
  c.SaveAs(base + ".png");
  printf("\nplots: %s.png, .pdf; lines: %s_lines.csv; histograms: %s_bins.csv\n", base.Data(), base.Data(), base.Data());
}
