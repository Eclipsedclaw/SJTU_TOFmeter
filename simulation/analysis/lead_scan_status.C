// Chart of scripts/lead_scan_status.sh: which points of a lead x energy scan are finished, running
// or pending, from <dir>/scan_status.csv (lead_cm, E_GeV or "check", status 1 / 0.5 / 0).
//   sh scripts/lead_scan_status.sh [definition]      (runs this macro)
//   root -l -b -q 'analysis/lead_scan_status.C("output/leadscan_thick")'
// Writes <dir>/scan_status.png.

void lead_scan_status(const char* dir = "output/leadscan_thick")
{
  std::ifstream in(Form("%s/scan_status.csv", dir));
  if (!in) {
    printf("no %s/scan_status.csv: run sh scripts/lead_scan_status.sh first\n", dir);
    return;
  }
  struct Cell { double lead, energy, status; };
  std::vector<Cell> points;
  std::map<double, double> checks;   // lead -> status
  std::string line;
  std::getline(in, line);  // header
  while (std::getline(in, line)) {
    std::stringstream ss(line);
    std::string a, b, c;
    if (!std::getline(ss, a, ',') || !std::getline(ss, b, ',') || !std::getline(ss, c, ',')) continue;
    if (b == "check") checks[std::stod(a)] = std::stod(c);
    else points.push_back({std::stod(a), std::stod(b), std::stod(c)});
  }
  if (points.empty()) return;
  std::vector<double> leads, energies;
  for (const auto& p : points) {
    if (std::find(leads.begin(), leads.end(), p.lead) == leads.end()) leads.push_back(p.lead);
    if (std::find(energies.begin(), energies.end(), p.energy) == energies.end()) energies.push_back(p.energy);
  }
  std::sort(leads.begin(), leads.end());
  std::sort(energies.begin(), energies.end());
  const int nT = int(leads.size()), nE = int(energies.size());
  // bins: geometric midpoints between the energies, thickness +- half a step
  std::vector<double> eEdges(nE + 1), tEdges(nT + 1);
  for (int i = 1; i < nE; ++i) eEdges[i] = std::sqrt(energies[i - 1] * energies[i]);
  eEdges[0] = nE > 1 ? energies[0] * energies[0] / eEdges[1] : 0.5 * energies[0];
  eEdges[nE] = nE > 1 ? energies[nE - 1] * energies[nE - 1] / eEdges[nE - 1] : 2. * energies[0];
  for (int i = 1; i < nT; ++i) tEdges[i] = 0.5 * (leads[i - 1] + leads[i]);
  tEdges[0] = leads[0] - (nT > 1 ? 0.5 * (leads[1] - leads[0]) : 0.5);
  tEdges[nT] = leads[nT - 1] + (nT > 1 ? 0.5 * (leads[nT - 1] - leads[nT - 2]) : 0.5);

  // colour codes 1 pending, 2 running, 3 finished
  int nDone = 0, nRun = 0, cDone = 0;
  auto code = [](double s) { return s >= 0.75 ? 3. : (s >= 0.25 ? 2. : 1.); };
  auto map = new TH2D("status", "", nE, eEdges.data(), nT, tEdges.data());
  for (const auto& p : points) {
    map->Fill(p.energy, p.lead, code(p.status));
    if (p.status >= 0.75) ++nDone;
    else if (p.status >= 0.25) ++nRun;
  }
  auto check = new TH2D("check", "", 1, 0., 1., nT, tEdges.data());
  for (const auto& [lead, s] : checks) {
    check->Fill(0.5, lead, code(s));
    if (s >= 0.75) ++cDone;
  }

  int colours[3] = {TColor::GetColor("#d9d9d9"), TColor::GetColor("#f39c12"), TColor::GetColor("#27ae60")};
  gStyle->SetPalette(3, colours);
  gStyle->SetOptStat(0);
  for (TH2D* h : {map, check}) {
    h->SetMinimum(0.5);
    h->SetMaximum(3.5);
    h->SetContour(3);
  }
  TCanvas c("c", "scan status", 1500, 950);
  auto main = new TPad("main", "", 0., 0., 0.88, 1.);
  auto side = new TPad("side", "", 0.88, 0., 1., 1.);
  main->Draw();
  side->Draw();
  main->cd();
  gPad->SetLogx();
  gPad->SetLeftMargin(0.07);
  gPad->SetRightMargin(0.01);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.1);
  map->SetTitle(Form("%s: %d of %d points finished, %d running;muon kinetic energy (GeV);lead thickness (cm)",
                     dir, nDone, int(points.size()), nRun));
  map->GetXaxis()->SetMoreLogLabels();
  map->GetXaxis()->SetTitleOffset(1.2);
  map->Draw("COL");
  auto legend = new TLegend(0.07, 0.905, 0.5, 0.95);
  legend->SetNColumns(3);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  const char* names[3] = {"pending", "running", "finished"};
  for (int k = 0; k < 3; ++k) {
    auto box = new TH1D(Form("box%d", k), "", 1, 0., 1.);
    box->SetFillColor(colours[k]);
    box->SetLineColor(colours[k]);
    legend->AddEntry(box, names[k], "f");
  }
  legend->Draw();
  side->cd();
  gPad->SetLeftMargin(0.05);
  gPad->SetRightMargin(0.3);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.1);
  check->GetXaxis()->SetLabelSize(0);
  check->GetXaxis()->SetTickLength(0);
  check->GetYaxis()->SetLabelSize(0.12);
  check->Draw("COL Y+");
  auto label = new TLatex();
  label->SetNDC();
  label->SetTextSize(0.09);
  label->DrawLatex(0.05, 0.955, "check runs");
  label->DrawLatex(0.05, 0.915, Form("%d of %d", cDone, int(checks.size())));
  c.SaveAs(Form("%s/scan_status.png", dir));
}
