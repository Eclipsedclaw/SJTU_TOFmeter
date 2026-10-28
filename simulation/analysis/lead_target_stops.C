// Muon stops in the muonic X-ray target versus the thickness of the lead layer on top, from the lead-slab runs of
// scripts/lead_filter.sh (macro/lead_filter.mac: EXPACS sea-level muons through an 8 x 8 m lead slab in vacuum). Those
// runs sample the spectrum continuously, so they resolve the narrow window of energies that stop in a thin target (the
// mono-energetic grid of the lead x energy scans is too coarse for it). Every muon that leaves the slab from below enters
// the top face of the target at a random point, along the direction from its start to its end (so with its scattering
// in the lead), and is carried through what the default setup has between the lead and the target (steel plate, top
// TOF, Al1, middle TOF, Al2, air; each only where the track crosses it) with CSDA ranges, Bethe-Bloch with the density
// effect; it stops in the target if its range ends inside. Lead 0 cm: the muons of the thinnest slab before the slab.
// The slab stands for the top lead layer only: a track that meets z = 91 cm outside the inner faces of the walls
// (|x| or |y| > 45 cm) would come through a wall instead, so it is left out (its share is printed). Thicker targets keep
// the top at 36.5 cm and grow down towards the stack (at most 4.9 cm); for them also the share of the X-rays of a stop
// that reach the top face of the stack without interacting in the target (isotropic emission, NIST attenuation of Al).
// The ranges are checked against Geant4 on the slab itself (near-vertical muons through the lead).
//   root -l -b -q 'analysis/lead_target_stops.C("output/lead_filter_below2GeV")'
// Writes <out>/lead_target_stops.png/.pdf/.csv.

#include <fstream>
#include <sstream>

namespace {
constexpr double kMuon = 105.6584, kElectron = 0.51099895;   // MeV

double Interpolate(const std::vector<double>& x, const std::vector<double>& y, double v)
{
  std::size_t j = std::lower_bound(x.begin(), x.end(), v) - x.begin();
  j = std::clamp<std::size_t>(j, 1, x.size() - 1);
  return y[j - 1] + (y[j] - y[j - 1]) * (v - x[j - 1]) / (x[j] - x[j - 1]);
}

// a material for the muon's CSDA range: Bethe-Bloch with the density effect (Sternheimer parameters as in the PDG tables)
struct Medium {
  double za, iEv, rho, cBar, x0, x1, a, k, delta0;
  std::vector<double> logT, logR;
  double DEdx(double T) const   // MeV cm2/g
  {
    const double g = 1. + T / kMuon, bg2 = g * g - 1., b2 = bg2 / (g * g);
    const double wMax = 2. * kElectron * bg2 / (1. + 2. * g * kElectron / kMuon + std::pow(kElectron / kMuon, 2));
    const double x = 0.5 * std::log10(bg2), ln10 = std::log(10.);
    const double delta = x >= x1 ? 2. * ln10 * x - cBar
                         : x >= x0 ? 2. * ln10 * x - cBar + a * std::pow(x1 - x, k)
                                   : delta0 * std::pow(10., 2. * (x - x0));
    const double i = iEv * 1e-6;
    return 0.307075 * za / b2 * (0.5 * std::log(2. * kElectron * bg2 * wMax / (i * i)) - b2 - 0.5 * delta);
  }
  void Build()   // range table 1 MeV - 20 GeV; below 1 MeV (< 0.03 g/cm2) the range is taken as T / (2 dE/dx)
  {
    const int n = 4000;
    const double step = std::pow(2e4, 1. / n);
    double t = 1., r = 0.5 / DEdx(1.);
    for (int j = 0; j <= n; ++j, t *= step) {
      logT.push_back(std::log(t));
      logR.push_back(std::log(r));
      r += 0.5 * (1. / DEdx(t) + 1. / DEdx(t * step)) * t * (step - 1.);
    }
  }
  double Range(double T) const { return std::exp(Interpolate(logT, logR, std::log(T))); }   // g/cm2
  double After(double T, double x) const   // kinetic energy after x g/cm2, -1 if it stops in them
  {
    const double r = Range(T) - x;
    return r > 0. ? std::exp(Interpolate(logR, logT, std::log(r))) : -1.;
  }
};

struct Layer {
  const Medium* medium;
  double zTop, zBottom, halfX, halfY;   // cm
};
}  // namespace

void lead_target_stops(const char* dir = "output/lead_filter_below2GeV", const char* out = "output/lead_target_stops",
                       const char* thicknesses = "0 5 10 20 30 40 50", unsigned seed = 4357)
{
  TH1::AddDirectory(false);
  Medium steel{0.4658, 286.0, 8.00, 4.2911, -0.0012, 3.1531, 0.14680, 2.9632, 0.12};   // G4_STAINLESS-STEEL (Fe parameters)
  Medium plastic{0.54141, 64.7, 1.032, 3.1997, 0.1464, 2.4855, 0.16101, 3.2393, 0.};   // G4_PLASTIC_SC_VINYLTOLUENE
  Medium aluminium{0.48181, 166.0, 2.699, 4.2395, 0.1708, 3.0127, 0.08024, 3.6345, 0.12};
  Medium air{0.49919, 85.7, 1.205e-3, 10.5961, 1.7418, 4.2759, 0.10914, 3.3994, 0.};
  Medium lead{0.39575, 823.0, 11.35, 6.2018, 0.3776, 3.8073, 0.09359, 3.1608, 0.14};   // only for the check against Geant4
  for (Medium* m : {&steel, &plastic, &aluminium, &air, &lead}) m->Build();
  // between the lead and the target, from the top (README geometry table; the TOF stations are 9 bars of 40 x 4 cm)
  const std::vector<Layer> layers = {{&steel, 91., 90., 40., 40.},       {&plastic, 70.4, 69.8, 20., 18.},
                                     {&aluminium, 67.8, 67.5, 40., 31.}, {&plastic, 40.4, 39.8, 18., 20.},
                                     {&aluminium, 37.8, 37.5, 31., 40.}};
  const double airPath = 51.7;   // cm of air from the lead to the target top
  const double targetTop = 36.5, targetHalf = 10., stackTop = 31.6, stackHalf = 10., leadBottom = 91., wallInner = 45.;
  const std::vector<double> depths = {0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 4.9};   // target thicknesses (cm)
  const int nD = depths.size(), iOne = 1;
  // X-ray attenuation in Al (NIST XCOM, with coherent scattering, interpolated log-log): L_alpha 65.7, L_beta 88.7,
  // K_alpha 354.9 keV
  const double mu[3] = {0.2510 * 2.699, 0.1866 * 2.699, 0.09737 * 2.699};   // 1/cm
  const int nRays = 16;

  // the files of each thickness
  std::map<int, std::vector<TString>> files;
  if (void* d = gSystem->OpenDirectory(dir)) {
    while (const char* entry = gSystem->GetDirEntry(d)) {
      TString name = entry;
      if (!name.BeginsWith("lead_filter_t") || !name.EndsWith(".root")) continue;
      TString rest = name(13, name.Length() - 18);
      if (rest.Index("_p") != kNPOS) rest = rest(0, rest.Index("_p"));
      files[rest.Atoi()].push_back(Form("%s/%s", dir, entry));
    }
    gSystem->FreeDirectory(d);
  }
  if (files.empty()) { printf("no lead_filter_t<cm>.root in %s\n", dir); return; }

  struct Result {
    int t = 0;
    double days = 0.;                            // exposure of the files
    std::vector<double> stop[2], stop2[2];      // [mu-, mu+][target thickness]: rate per cm2 per s, its sum of w^2
    double wall[2] = {0., 0.}, reach[2] = {0., 0.};   // 1 cm stops through a wall; muons reaching the target face
    std::vector<double> xray[3];                 // mu- stops x share of the X-rays reaching the stack
    double rangeLoss = 0., rangePath = 0.;       // check against Geant4: lost CSDA range / straight path in the lead
  };
  std::vector<Result> results;
  TRandom3 random(seed);
  std::istringstream list(thicknesses);
  for (int t; list >> t;) {
    const bool before = t == 0;                  // no lead: the muons of the thinnest slab before it
    const int source = before ? files.begin()->first : t;
    if (!files.count(source)) { printf("no files for %d cm\n", t); continue; }
    Result res;
    res.t = t;
    for (int q = 0; q < 2; ++q) { res.stop[q].assign(nD, 0.); res.stop2[q].assign(nD, 0.); }
    for (auto& x : res.xray) x.assign(nD, 0.);
    double liveTime = 0., area = 0.;
    for (const auto& path : files[source]) {
      TFile file(path);
      auto run = dynamic_cast<TTree*>(file.Get("run"));
      if (!run || run->GetEntries() != 1) continue;
      double lt = 0., a = 0.;
      run->SetBranchAddress("live_time", &lt);
      run->SetBranchAddress("area", &a);
      run->GetEntry(0);
      liveTime += lt;
      area = a;
    }
    const double w = 1. / (liveTime * area);   // per cm2 per s
    res.days = liveTime / 86400.;
    for (const auto& path : files[source]) {
      TFile file(path);
      auto events = dynamic_cast<TTree*>(file.Get("events"));
      if (!events) continue;
      int pdg = 0;
      double ekin = 0., costh = 0., phi = 0., x0 = 0., y0 = 0., z0 = 0., endX = 0., endY = 0., endZ = 0., endKe = 0.;
      events->SetBranchStatus("*", false);
      for (const char* b : {"pdg", "ekin", "costh", "phi", "x0", "y0", "z0", "end_x", "end_y", "end_z", "end_ke"})
        events->SetBranchStatus(b, true);
      events->SetBranchAddress("pdg", &pdg); events->SetBranchAddress("ekin", &ekin);
      events->SetBranchAddress("costh", &costh); events->SetBranchAddress("phi", &phi);
      events->SetBranchAddress("x0", &x0); events->SetBranchAddress("y0", &y0); events->SetBranchAddress("z0", &z0);
      events->SetBranchAddress("end_x", &endX); events->SetBranchAddress("end_y", &endY);
      events->SetBranchAddress("end_z", &endZ); events->SetBranchAddress("end_ke", &endKe);
      const Long64_t n = events->GetEntries();
      for (Long64_t i = 0; i < n; ++i) {
        events->GetEntry(i);
        double T, ux, uy, uz;   // kinetic energy (MeV) and direction under the lead
        if (before) {
          const double s = std::sqrt(std::max(0., 1. - costh * costh));
          T = ekin; ux = s * std::cos(phi); uy = s * std::sin(phi); uz = -costh;
        } else {
          if (endZ >= -0.5 * t - 1e-3) continue;   // stopped in the slab
          const double dx = endX - x0, dy = endY - y0, dz = endZ - z0, norm = std::sqrt(dx * dx + dy * dy + dz * dz);
          T = endKe; ux = dx / norm; uy = dy / norm; uz = dz / norm;
          // check of the ranges: near-vertical muons well above the end of their range, in the lead
          const double path = t * lead.rho / costh;
          if (costh > 0.97 && ekin < 1000. && lead.After(ekin, path) > 100.) {
            res.rangeLoss += lead.Range(ekin) - lead.Range(endKe);
            res.rangePath += path;
          }
        }
        if (T > 400.) continue;   // too fast to stop in the target
        const double cz = -uz;
        const double xe = random.Uniform(-targetHalf, targetHalf), ye = random.Uniform(-targetHalf, targetHalf);
        auto at = [&](double z, double& x, double& y) { x = xe + ux / uz * (z - targetTop); y = ye + uy / uz * (z - targetTop); };
        double x, y;
        at(leadBottom, x, y);
        const bool fromTop = std::abs(x) < wallInner && std::abs(y) < wallInner;
        for (const auto& layer : layers) {
          at(0.5 * (layer.zTop + layer.zBottom), x, y);
          if (std::abs(x) < layer.halfX && std::abs(y) < layer.halfY)
            T = layer.medium->After(T, (layer.zTop - layer.zBottom) * layer.medium->rho / cz);
          if (T <= 0.) break;
        }
        if (T > 0.) T = air.After(T, airPath * air.rho / cz);
        if (T <= 0.) continue;   // stopped before the target
        const int q = pdg == 13 ? 0 : 1;
        if (fromTop) res.reach[q] += w;
        const double length = aluminium.Range(T) / aluminium.rho;   // cm, straight
        const double xs = xe + ux * length, ys = ye + uy * length, depth = cz * length;
        if (depth >= depths.back() || std::abs(xs) >= targetHalf || std::abs(ys) >= targetHalf) continue;   // not a stop
        if (!fromTop) {
          if (depth < depths[iOne]) res.wall[q] += w;
          continue;
        }
        for (int k = 0; k < nD; ++k)
          if (depth < depths[k]) { res.stop[q][k] += w; res.stop2[q][k] += w * w; }
        if (q != 0) continue;
        // its X-rays: isotropic, those towards the top face of the stack, attenuated on their way out of the target
        const double zs = targetTop - depth;
        for (int r = 0; r < nRays; ++r) {
          const double vz = random.Uniform(-1., 1.), a = random.Uniform(0., 2. * M_PI), s = std::sqrt(1. - vz * vz);
          const double vx = s * std::cos(a), vy = s * std::sin(a);
          if (vz >= 0.) continue;
          const double toStack = (stackTop - zs) / vz;
          if (std::abs(xs + vx * toStack) >= stackHalf || std::abs(ys + vy * toStack) >= stackHalf) continue;
          const double side = std::min(vx > 0. ? (targetHalf - xs) / vx : vx < 0. ? (-targetHalf - xs) / vx : 1e9,
                                       vy > 0. ? (targetHalf - ys) / vy : vy < 0. ? (-targetHalf - ys) / vy : 1e9);
          for (int k = 0; k < nD; ++k) {
            if (depth >= depths[k]) continue;
            const double inside = std::min(side, (targetTop - depths[k] - zs) / vz);
            for (int e = 0; e < 3; ++e) res.xray[e][k] += w * std::exp(-mu[e] * inside) / nRays;
          }
        }
      }
    }
    results.push_back(res);
    printf("lead %2d cm done (%s, %.2f days)\n", t, before ? "before the thinnest slab" : "after the slab", res.days);
  }
  if (results.empty()) return;

  // ---- tables
  const double perDay = 4. * targetHalf * targetHalf * 86400.;   // per cm2 per s -> per day on the target face
  const Result* ref = nullptr;
  for (const auto& r : results) if (r.t == 10) ref = &r;
  if (!ref) ref = &results.front();
  printf("\n=== Muons stopping in the 20 x 20 x 1 cm Al target versus the top lead layer (%s) ===\n", dir);
  printf("tracks through the top lead layer only; CSDA ranges through the default setup between the lead and the target\n\n");
  printf(" lead   mu- stops per day    vs %d cm   mu+ stops per day   through a wall    mu- reaching the     range check\n", ref->t);
  printf("   cm   (1 cm target)                                     (left out)        target face per day  G4 / CSDA\n");
  gSystem->mkdir(out, true);
  std::ofstream csv(Form("%s/lead_target_stops.csv", out));
  csv << "lead_cm,mu_minus_stops_day,err,relative,mu_plus_stops_day,err,wall_share,mu_minus_reaching_day,range_ratio\n";
  auto graph = new TGraphErrors();
  for (const auto& r : results) {
    const double s = r.stop[0][iOne] * perDay, ds = std::sqrt(r.stop2[0][iOne]) * perDay;
    const double sp = r.stop[1][iOne] * perDay, dsp = std::sqrt(r.stop2[1][iOne]) * perDay;
    const double wallShare = (r.wall[0] + r.wall[1]) / (r.wall[0] + r.wall[1] + r.stop[0][iOne] + r.stop[1][iOne]);
    const double rel = s / (ref->stop[0][iOne] * perDay);
    const double check = r.rangePath > 0. ? r.rangeLoss / r.rangePath : 0.;
    printf("   %2d   %7.1f +- %4.1f   %6.3f   %7.1f +- %4.1f     %5.1f%%        %9.0f          %s\n", r.t, s, ds, rel, sp, dsp,
           100. * wallShare, r.reach[0] * perDay, check > 0. ? Form("%.3f", check) : "-");
    csv << r.t << "," << s << "," << ds << "," << rel << "," << sp << "," << dsp << "," << wallShare << ","
        << r.reach[0] * perDay << "," << check << "\n";
    const int p = graph->GetN();
    graph->SetPoint(p, r.t, s);
    graph->SetPointError(p, 0., ds);
  }
  const auto best = std::max_element(results.begin(), results.end(),
                                     [](const Result& a, const Result& b) { return a.stop[0][1] < b.stop[0][1]; });
  printf("most mu- stops: %d cm of lead\n", best->t);

  printf("\nthicker targets with %d cm of lead (top at %.1f cm, growing down towards the stack):\n", ref->t, targetTop);
  printf("  target   mu- stops    vs 1 cm   share of the X-rays reaching the stack     X-rays at the stack vs 1 cm\n");
  printf("    cm      per day                %-12s %-12s %-12s      %-8s %-8s %-8s\n", "L_alpha", "L_beta", "K_alpha", "L_alpha", "L_beta", "K_alpha");
  csv << "\ntarget_cm,mu_minus_stops_day,relative,share_Lalpha,share_Lbeta,share_Kalpha,xray_rel_Lalpha,xray_rel_Lbeta,xray_rel_Kalpha\n";
  TGraph* gTarget[4];
  for (auto& g : gTarget) g = new TGraph();
  for (int k = 0; k < nD; ++k) {
    const double s = ref->stop[0][k] * perDay, rel = ref->stop[0][k] / ref->stop[0][iOne];
    double share[3], xrel[3];
    for (int e = 0; e < 3; ++e) {
      share[e] = ref->xray[e][k] / ref->stop[0][k];
      xrel[e] = ref->xray[e][k] / ref->xray[e][iOne];
    }
    printf("   %4.1f   %8.1f    %6.2f     %8.4f     %8.4f     %8.4f          %6.2f   %6.2f   %6.2f\n", depths[k], s, rel,
           share[0], share[1], share[2], xrel[0], xrel[1], xrel[2]);
    csv << depths[k] << "," << s << "," << rel << "," << share[0] << "," << share[1] << "," << share[2] << "," << xrel[0]
        << "," << xrel[1] << "," << xrel[2] << "\n";
    gTarget[0]->SetPoint(k, depths[k], rel);
    for (int e = 0; e < 3; ++e) gTarget[e + 1]->SetPoint(k, depths[k], xrel[e]);
  }
  printf("(share: of the photons of a stop, isotropic, the part that hits the top face of the stack without interacting\n"
         " in the target; Al attenuation %.3f, %.3f, %.3f /cm)\n", mu[0], mu[1], mu[2]);

  // ---- plots
  gStyle->SetOptStat(0);
  TCanvas c("c", "lead target stops", 1400, 600);
  c.Divide(2, 1);
  c.cd(1);
  gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.1);
  double yMax = 0.;
  for (int p = 0; p < graph->GetN(); ++p) yMax = std::max(yMax, graph->GetY()[p] + graph->GetEY()[p]);
  auto frame = gPad->DrawFrame(-2., 0., 52., 1.25 * yMax, ";top lead thickness (cm);#mu^{-} stops per day in the 1 cm target");
  frame->GetYaxis()->SetTitleOffset(1.3);
  graph->SetMarkerStyle(20); graph->SetMarkerColor(kAzure + 2); graph->SetLineColor(kAzure + 2); graph->SetLineWidth(2);
  graph->Draw("PL SAME");
  TLatex label;
  label.SetTextSize(0.035); label.SetTextAlign(21);
  for (int p = 0; p < graph->GetN(); ++p)
    label.DrawLatex(graph->GetX()[p], graph->GetY()[p] + graph->GetEY()[p] + 0.03 * yMax,
                    Form("%.2f", graph->GetY()[p] / (ref->stop[0][iOne] * perDay)));
  TLatex heading;
  heading.SetNDC(); heading.SetTextSize(0.04);
  heading.DrawLatex(0.13, 0.93, Form("Stops in the Al target, labels: relative to %d cm of lead", ref->t));
  c.cd(2);
  gPad->SetLeftMargin(0.11); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.1);
  double rMax = 0.;
  for (auto g : gTarget) for (int p = 0; p < g->GetN(); ++p) rMax = std::max(rMax, g->GetY()[p]);
  auto frame2 = gPad->DrawFrame(0., 0., 5., 1.15 * rMax, ";Al target thickness (cm);relative to the 1 cm target");
  frame2->GetYaxis()->SetTitleOffset(1.1);
  const int colours[4] = {kGray + 2, kRed + 1, kOrange + 1, kAzure + 2};
  const char* names[4] = {"#mu^{-} stops", "L_{#alpha} X-rays at the stack", "L_{#beta} X-rays at the stack",
                          "K_{#alpha} X-rays at the stack"};
  auto legend = new TLegend(0.14, 0.66, 0.6, 0.88);
  legend->SetBorderSize(0); legend->SetFillStyle(0);
  for (int g = 0; g < 4; ++g) {
    gTarget[g]->SetLineColor(colours[g]); gTarget[g]->SetMarkerColor(colours[g]);
    gTarget[g]->SetLineWidth(2); gTarget[g]->SetMarkerStyle(20);
    gTarget[g]->Draw("PL SAME");
    legend->AddEntry(gTarget[g], names[g], "lp");
  }
  legend->Draw();
  heading.DrawLatex(0.11, 0.93, Form("Thicker targets with %d cm of lead (top fixed, growing down)", ref->t));
  c.Print(Form("%s/lead_target_stops.pdf", out));
  c.SaveAs(Form("%s/lead_target_stops.png", out));
  printf("\nplots: %s/lead_target_stops.png (and .pdf); table: %s/lead_target_stops.csv\n", out, out);
}
