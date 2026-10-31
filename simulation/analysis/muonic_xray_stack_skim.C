// Skim of a run of scripts/muonic_xray_stack.sh for the X-ray selection studies (analysis/muonic_xray_stack_cuts.C):
// the detector response is applied once, the add-back entries are made, and only the candidates are kept, with what
// a selection can use. Processes first-last (finished ones, log/job<k>.done) into one file; several in parallel, then
// hadd:
//   root -l -b -q 'analysis/muonic_xray_stack_skim.C("output/muonic_xray_stack", 1, 39, "skim_1.root")'
// Response of the stack (as analysis/muonic_xray_stack.C): energy sigma = 0.04 MeV x sqrt(E / MeV) per hit, x and y
// sigma = 2 mm, hits below 10 keV dropped; add-back: the hits within 100 ns of the first. A TOF station is hit when a
// bar has more than 0.2 MeV; its time is the earliest such bar, its bar and energy those of the largest one (true
// values: a TOF response is applied in the cut study). Kept: entries of 20-1000 keV in events with a middle-TOF hit,
// or whose primary muon stopped (truth). Writes tree "cand" (one row per entry) and tree "info" (live time).
// detector "czt": the hits of the CZT plane (czt_*, scripts/muonic_xray_czt.sh) instead of the stack's, a pixel in place
// of a layer (layers: bit k a hit in pixel k; topLayer: the pixel number), e.g. for one band of that run, with a CZT's
// FWHM/E of 1.5% at 1 MeV (sigma 15 keV / 2.3548 = 6.37 keV at 1 MeV):
//   root -l -b -q 'analysis/muonic_xray_stack_skim.C("output/muonic_xray_czt/slow", 1, 8, "xray_skim.root", 4357, "czt", 6.37)'
// Each of the 16 pixels is read out alone: an entry per pixel hit (n = 1), no sum over pixels (a stopped muon's cascade
// sends several X-rays at once, e.g. Al L_alpha and K_alpha, which a sum would merge); pix is the pixel, npix the pixels
// with a hit within 100 ns of it (itself included) and eOther their energy without it. A CZT hit's x and y are its
// pixel's centre (2.5 cm pitch), not the smeared position of a stack hit.
// sigma1MeV: the energy sigma at 1 MeV (keV), 40 for the stack. etrue: the entry's energy before the response (keV).
// The TOF as measured, for both detectors, in branches of their own (the true values above stay for the cut studies):
// every bar's energy with tofFwhmE (FWHM / E, 10% by default), a bar hit above 0.2 MeV of it; a station's time that of
// its earliest hit bar with tofFwhmT (FWHM, 300 ps); both ends of every bar read out, so a station gives a point: across
// the bars the centre of its largest bar (4 cm pitch), along it the energy-weighted point with sigma vEff x sigma_t (the
// mean time of the two ends has sigma_t, half their difference times vEff the position; vEff 15 cm/ns by default).
// trig: the top and the middle station hit within 20 ns, the trigger. With it the track, the straight line through
// the two points, with its length between the planes, the speed beta = length / (c x time of flight) and the kinetic
// energy of a muon of that speed (tkin), and the point where the track crosses z = targetZ (the target's centre plane,
// 34.5 cm for the 4 cm target); the bottom station's largest bar and its point (barBotM, xBot, yBot). The TOF response
// has a random generator of its own, seeded per event, so the stack's and the CZT's hits are smeared as before. Every
// condition goes into the file: TNamed "conditions" (text, also printed), TNamed "detector" and TParameter<double>
// sigma1MeV, tofFwhmT, tofFwhmE, vEff and targetZ. Skims made in parts: merge them with analysis/skim_merge.C, which
// merges their conditions too (hadd keeps only the first part's text).

#include <TParameter.h>

#include <sstream>

namespace {
// where the primary ended: 0 target, 1 Al plates, 2 stack Si, 3 stack YSO, 4 TOF bars, 5 lead, 6 steel, 7 air or out,
// 8 CZT plane
int VolumeCode(const std::string& v)
{
  if (v.rfind("Target", 0) == 0) return 0;
  if (v.rfind("CZTPixel", 0) == 0) return 8;
  if (v.rfind("FrameworkAl", 0) == 0) return 1;
  if (v.rfind("SiliconLayer", 0) == 0) return 2;
  if (v.rfind("YSOLayer", 0) == 0) return 3;
  if (v.rfind("Scintillator", 0) == 0) return 4;
  if (v.rfind("Absorber", 0) == 0 || v.rfind("LeadWall", 0) == 0) return 5;
  if (v.rfind("FrameworkFe", 0) == 0 || v.rfind("SupportBeam", 0) == 0) return 6;
  return 7;
}
// TOF stations as DetectorConstruction.cc builds them: planes of the top, middle and bottom bars (cm), bars of 40 cm,
// 4 cm wide with their centres at -16 ... +16 cm across them, the middle ones along y and the others along x
const double kStationZ[3] = {70.1, 40.1, 2.6};
const double kBarHalfLength = 20.;
double BarCentre(int bar) { return -16. + 4. * bar; }
const double kMuonMass = 105.6584;   // MeV
const double kLightSpeed = 29.9792458;   // cm/ns
}  // namespace

void muonic_xray_stack_skim(const char* dir = "output/muonic_xray_stack", int first = 1, int last = 385,
                            const char* out = "skim.root", unsigned seed = 4357, const char* detector = "stack",
                            double sigma1MeV = 40., double tofFwhmT = 0.3, double tofFwhmE = 0.10, double vEff = 15.,
                            double targetZ = 34.5)
{
  const bool czt = std::string(detector) == "czt";
  const std::string hit = czt ? "czt_" : "stk_";   // the hit columns read
  const double fwhm = 2. * std::sqrt(2. * std::log(2.));   // 2.3548
  const double tofSigmaT = tofFwhmT / fwhm, tofSigmaE = tofFwhmE / fwhm, tofSigmaX = vEff * tofSigmaT;
  TFile output(Form("%s/%s", dir, out), "RECREATE");
  auto cand = new TTree("cand", "add-back X-ray candidates");
  auto info = new TTree("info", "exposure of the skimmed processes");
  int proc, event, n, layers, topLayer, src, pdg, endVol, endProc, tofTop, tofMid, barTop, barMid;
  float e, etrue, emax, rmax, t, tTop, tMid, eTop, eMid, botFirst, ekin, endT, endX, endY, endZ, pbKe;
  int pix, npix, trig, primSt, barBotM;
  float eOther, tTopM, tMidM, eTopM, eMidM, xTop, yTop, xMid, yMid, tof, pathL, beta, tkin, xTgt, yTgt, slopeX, slopeY,
      botM, xBot, yBot, keTop;
  cand->Branch("proc", &proc); cand->Branch("event", &event);
  cand->Branch("e", &e);               // add-back energy (keV, smeared); CZT: the pixel's energy
  cand->Branch("etrue", &etrue);       // the same hits before the response (keV)
  cand->Branch("n", &n);               // hits added up (CZT: 1)
  cand->Branch("layers", &layers);     // bit k: a hit in layer k (0, 2, 4 Si; 1, 3 YSO; from the bottom); CZT: pixel k
  cand->Branch("topLayer", &topLayer); // highest layer with a hit; CZT: the pixel
  cand->Branch("emax", &emax);         // largest hit (keV)
  cand->Branch("rmax", &rmax);         // largest |x| or |y| of its hits (cm, smeared): how close to the stack edge
  cand->Branch("t", &t);               // time of its first hit (ns)
  cand->Branch("src", &src);           // truth: stk_src bits of its hits (1 X-ray of the target)
  cand->Branch("pix", &pix);           // CZT: the pixel (copy number ix * 4 + iy); -1 for the stack
  cand->Branch("npix", &npix);         // CZT: pixels with a hit within 100 ns of this one, itself included; 0 stack
  cand->Branch("eOther", &eOther);     // CZT: their measured energy without this one (keV)
  cand->Branch("tofTop", &tofTop); cand->Branch("tTop", &tTop); cand->Branch("barTop", &barTop); cand->Branch("eTop", &eTop);
  cand->Branch("tofMid", &tofMid); cand->Branch("tMid", &tMid); cand->Branch("barMid", &barMid); cand->Branch("eMid", &eMid);
  cand->Branch("botFirst", &botFirst); // earliest bottom-TOF hit (ns), -1 if none
  // the TOF as measured, and what it reconstructs
  cand->Branch("trig", &trig);         // the trigger: top and middle stations hit within 20 ns
  cand->Branch("tTopM", &tTopM); cand->Branch("tMidM", &tMidM);   // station times (ns), -1 if not hit
  cand->Branch("eTopM", &eTopM); cand->Branch("eMidM", &eMidM);   // energy of the station's largest bar (MeV), 0 if not hit
  cand->Branch("xTop", &xTop); cand->Branch("yTop", &yTop);       // the stations' points (cm)
  cand->Branch("xMid", &xMid); cand->Branch("yMid", &yMid);
  cand->Branch("botM", &botM);         // time of the earliest bottom-TOF hit (ns), -1 if none
  cand->Branch("barBotM", &barBotM);   // the bottom station's largest bar (0-8), -1 if none
  cand->Branch("xBot", &xBot); cand->Branch("yBot", &yBot);       // and its point (cm)
  cand->Branch("tof", &tof);           // top-to-middle time of flight (ns), -1 without a trigger
  cand->Branch("pathL", &pathL);       // track length between the top and middle planes (cm)
  cand->Branch("beta", &beta);         // pathL / (c x tof)
  cand->Branch("tkin", &tkin);         // kinetic energy of a muon of that speed (MeV); 1e5 if beta >= 1 or no trigger
  cand->Branch("xTgt", &xTgt); cand->Branch("yTgt", &yTgt);       // the track at z = targetZ (cm); 1000 without trigger
  cand->Branch("slopeX", &slopeX); cand->Branch("slopeY", &slopeY); // its dx/dz and dy/dz
  cand->Branch("keTop", &keTop);       // truth: the primary's kinetic energy entering a top-TOF bar (MeV), -1 if never
  cand->Branch("primSt", &primSt);     // truth: the TOF stations the primary crossed (bit 0 top, 1 middle, 2 bottom)
  cand->Branch("pdg", &pdg); cand->Branch("ekin", &ekin);   // truth: the primary
  cand->Branch("endVol", &endVol);     // truth: where it ended (VolumeCode)
  cand->Branch("endProc", &endProc);   // truth: 1 mu- captured at rest, 2 decay at rest, 0 else
  cand->Branch("endT", &endT); cand->Branch("endX", &endX); cand->Branch("endY", &endY); cand->Branch("endZ", &endZ);
  cand->Branch("pbKe", &pbKe);         // truth: its kinetic energy leaving the top lead layer (MeV), -1 if it went round it
  double liveTime = 0.;
  long long muons = 0;
  int files = 0;
  info->Branch("live_time", &liveTime); info->Branch("muons", &muons); info->Branch("files", &files);

  TRandom3 random(seed + first);   // reproducible, different for every part
  TRandom3 tofRandom;              // the TOF response, seeded per event
  for (int k = first; k <= last; ++k) {
    if (gSystem->AccessPathName(Form("%s/log/job%d.done", dir, k))) continue;   // true: not finished
    TFile file(Form("%s/stk_%d.root", dir, k));
    auto run = dynamic_cast<TTree*>(file.Get("run"));
    auto events = dynamic_cast<TTree*>(file.Get("events"));
    if (!run || !events || run->GetEntries() != 1) continue;
    int nEvents = 0;
    double lt = 0.;
    run->SetBranchAddress("n_events", &nEvents);
    run->SetBranchAddress("live_time", &lt);
    run->GetEntry(0);
    liveTime += lt;
    muons += nEvents;
    ++files;
    int eventID = 0, primary = 0, primStations = 0;
    double kin = 0., tEnd = 0., xEnd = 0., yEnd = 0., zEnd = 0., keEnd = 0., kePb = -1., keTopBar = -1.;
    char vol[128] = "", process[128] = "";
    std::vector<int>*layer = nullptr, *hitSrc = nullptr, *barSt = nullptr, *barId = nullptr;
    std::vector<double>*he = nullptr, *ht = nullptr, *hx = nullptr, *hy = nullptr, *be = nullptr, *bt = nullptr,
    *bx = nullptr, *by = nullptr;
    const std::string cLayer = hit + (czt ? "pix" : "layer"), cSrc = hit + "src", cE = hit + "e", cT = hit + "t",
                      cX = hit + "x", cY = hit + "y";
    events->SetBranchStatus("*", 0);
    for (const char* c : {"eventID", "pdg", "ekin", "end_t", "end_x", "end_y", "end_z", "end_ke", "end_vol", "end_proc", "pb_ke", "top_ke",
                          "prim_st", cLayer.c_str(), cSrc.c_str(), cE.c_str(), cT.c_str(), cX.c_str(), cY.c_str(), "bar_st",
                          "bar_id", "bar_e", "bar_t", "bar_x", "bar_y"})
      events->SetBranchStatus(c, 1);
    events->SetBranchAddress("eventID", &eventID); events->SetBranchAddress("pdg", &primary);
    events->SetBranchAddress("ekin", &kin); events->SetBranchAddress("end_t", &tEnd); events->SetBranchAddress("end_z", &zEnd);
    events->SetBranchAddress("end_x", &xEnd); events->SetBranchAddress("end_y", &yEnd);
    events->SetBranchAddress("pb_ke", &kePb); events->SetBranchAddress("top_ke", &keTopBar);
    events->SetBranchAddress("prim_st", &primStations);
    events->SetBranchAddress("end_ke", &keEnd); events->SetBranchAddress("end_vol", vol); events->SetBranchAddress("end_proc", process);
    events->SetBranchAddress(cLayer.c_str(), &layer); events->SetBranchAddress(cSrc.c_str(), &hitSrc);
    events->SetBranchAddress(cE.c_str(), &he); events->SetBranchAddress(cT.c_str(), &ht);
    events->SetBranchAddress(cX.c_str(), &hx); events->SetBranchAddress(cY.c_str(), &hy);
    events->SetBranchAddress("bar_st", &barSt); events->SetBranchAddress("bar_id", &barId);
    events->SetBranchAddress("bar_e", &be); events->SetBranchAddress("bar_t", &bt);
    events->SetBranchAddress("bar_x", &bx); events->SetBranchAddress("bar_y", &by);
    struct Hit { double t, e, x, y; int layer, src; double etrue; };
    std::vector<Hit> hits;
    for (Long64_t i = 0; i < events->GetEntries(); ++i) {
      events->GetEntry(i);
      // TOF stations: hit above 0.2 MeV, earliest time, bar and energy of the largest bar
      tofTop = tofMid = 0; tTop = tMid = 1e9; barTop = barMid = -1; eTop = eMid = 0.; botFirst = -1.;
      for (std::size_t b = 0; b < be->size(); ++b) {
        if ((*be)[b] <= 0.2) continue;
        const int st = (*barSt)[b];
        const float tb = (*bt)[b], eb = (*be)[b];
        if (st == 0) { tofTop = 1; tTop = std::min(tTop, tb); if (eb > eTop) { eTop = eb; barTop = (*barId)[b]; } }
        if (st == 1) { tofMid = 1; tMid = std::min(tMid, tb); if (eb > eMid) { eMid = eb; barMid = (*barId)[b]; } }
        if (st == 2 && (botFirst < 0. || tb < botFirst)) botFirst = tb;
      }
      // the TOF as measured, seeded per event: the bar energies, then the stations' times and points
      tofRandom.SetSeed(1000003U * k + eventID + 1);
      bool sHit[3] = {false, false, false};
      double sT[3] = {1e9, 1e9, 1e9}, sE[3] = {0., 0., 0.}, sAlong[3] = {0., 0., 0.}, px[3] = {0., 0., 0.}, py[3] = {0., 0., 0.};
      int sBar[3] = {-1, -1, -1};
      for (std::size_t b = 0; b < be->size(); ++b) {
        const int st = (*barSt)[b];
        const double measured = (*be)[b] * (1. + tofSigmaE * tofRandom.Gaus());
        if (st < 0 || st > 2 || measured <= 0.2) continue;
        sHit[st] = true;
        sT[st] = std::min(sT[st], (*bt)[b]);
        if (measured > sE[st]) { sE[st] = measured; sBar[st] = (*barId)[b]; sAlong[st] = st == 1 ? (*by)[b] : (*bx)[b]; }
      }
      for (int st = 0; st < 3; ++st) {
        if (!sHit[st]) continue;
        sT[st] += tofSigmaT * tofRandom.Gaus();
        const double along = std::max(-kBarHalfLength, std::min(kBarHalfLength, sAlong[st] + tofSigmaX * tofRandom.Gaus()));
        px[st] = st == 1 ? BarCentre(sBar[st]) : along;
        py[st] = st == 1 ? along : BarCentre(sBar[st]);
      }
      trig = sHit[0] && sHit[1] && std::abs(sT[1] - sT[0]) < 20.;
      tTopM = sHit[0] ? sT[0] : -1.; tMidM = sHit[1] ? sT[1] : -1.; botM = sHit[2] ? sT[2] : -1.;
      eTopM = sE[0]; eMidM = sE[1];
      xTop = px[0]; yTop = py[0]; xMid = px[1]; yMid = py[1];
      barBotM = sHit[2] ? sBar[2] : -1; xBot = px[2]; yBot = py[2];
      tof = -1.; pathL = 0.; beta = 0.; tkin = 1e5; xTgt = yTgt = 1000.; slopeX = slopeY = 0.;
      if (trig) {
        const double dz = kStationZ[1] - kStationZ[0];
        slopeX = (px[1] - px[0]) / dz; slopeY = (py[1] - py[0]) / dz;
        pathL = std::sqrt(std::pow(px[1] - px[0], 2) + std::pow(py[1] - py[0], 2) + dz * dz);
        tof = sT[1] - sT[0];
        beta = tof > 0. ? pathL / (kLightSpeed * tof) : 1000.;
        if (beta < 1.) tkin = kMuonMass * (1. / std::sqrt(1. - beta * beta) - 1.);
        xTgt = px[1] + slopeX * (targetZ - kStationZ[1]);
        yTgt = py[1] + slopeY * (targetZ - kStationZ[1]);
      }
      const std::string procName = process;
      endProc = procName == "muMinusCaptureAtRest" ? 1 : (procName.rfind("Decay", 0) == 0 && keEnd <= 0.) ? 2 : 0;
      // a muon that stops in the target crosses the middle TOF (4 cm above it), often not the top one: keep every
      // event with a middle-TOF hit, and the stopped muons (truth)
      if (!tofMid && endProc == 0) continue;
      // the stack (CZT) hits with the response
      hits.clear();
      for (std::size_t h = 0; h < he->size(); ++h) {
        const double measured = (*he)[h] + random.Gaus(0., sigma1MeV / 1000. * std::sqrt((*he)[h]));   // MeV
        if (measured * 1000. < 10.) continue;
        // a CZT pixel tells only which pixel it is: its centre (pitch 2.5 cm, the plane centred on x = y = 0)
        const double x = czt ? (std::floor((*hx)[h] / 2.5) + 0.5) * 2.5 : (*hx)[h] + random.Gaus(0., 0.2);
        const double y = czt ? (std::floor((*hy)[h] / 2.5) + 0.5) * 2.5 : (*hy)[h] + random.Gaus(0., 0.2);
        hits.push_back({(*ht)[h], measured * 1000., x, y, (*layer)[h], (*hitSrc)[h], (*he)[h] * 1000.});
      }
      std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) { return a.t < b.t; });
      proc = k; event = eventID; pdg = primary; ekin = kin; endVol = VolumeCode(vol); endT = tEnd; endZ = zEnd; pbKe = kePb;
      endX = xEnd; endY = yEnd;
      keTop = keTopBar; primSt = primStations;
      if (czt) {   // each pixel read out alone: an entry per pixel hit
        for (const auto& h : hits) {
          t = h.t; e = h.e; etrue = h.etrue; n = 1; src = h.src;
          pix = h.layer; layers = 1 << pix; topLayer = pix; emax = e; rmax = std::max(std::abs(h.x), std::abs(h.y));
          npix = 0; eOther = -e;
          for (const auto& other : hits)
            if (std::abs(other.t - h.t) < 100.) { ++npix; eOther += other.e; }
          if (e >= 20. && e <= 1000.) cand->Fill();
        }
        continue;
      }
      pix = -1; npix = 0; eOther = 0.;
      for (std::size_t h = 0; h < hits.size();) {   // add-back entries
        t = hits[h].t;
        e = 0.; etrue = 0.; n = 0; layers = 0; topLayer = -1; emax = 0.; rmax = 0.; src = 0;
        for (; h < hits.size() && hits[h].t < t + 100.; ++h) {
          e += hits[h].e; etrue += hits[h].etrue; ++n;
          layers |= 1 << hits[h].layer;
          topLayer = std::max(topLayer, hits[h].layer);
          emax = std::max(emax, float(hits[h].e));
          rmax = std::max({rmax, float(std::abs(hits[h].x)), float(std::abs(hits[h].y))});
          src |= hits[h].src;
        }
        if (e >= 20. && e <= 1000.) cand->Fill();
      }
    }
  }
  info->Fill();

  // the conditions of this skim, in the file and printed
  std::ostringstream c;
  c << "skim: analysis/muonic_xray_stack_skim.C of " << dir << ", processes " << first << "-" << last << ": " << files
    << " finished, " << muons << " muons, live time " << liveTime / 86400. << " days; written " << TDatime().AsSQLString() << "\n";
  if (czt)
    c << "detector: the CZT plane (czt_* hits), each of its pixels read out alone: an entry per pixel hit, no sum over pixels; "
         "its x and y the pixel's centre (2.5 cm pitch); npix and eOther: the pixels with a hit within 100 ns\n";
  else
    c << "detector: the Si/YSO stack (stk_* hits): add-back of its hits within 100 ns of the first; x and y smeared by 2 mm\n";
  c << "energy response: sigma " << sigma1MeV << " keV x sqrt(E / MeV) per hit (FWHM/E " << fwhm * sigma1MeV / 10.
    << "% at 1 MeV), hits below 10 keV dropped; entries of 20-1000 keV\n";
  c << "TOF bar energy: sigma/E " << tofSigmaE << " (FWHM " << 100. * tofFwhmE << "%), a bar hit above 0.2 MeV of it\n";
  c << "TOF timing: a station's time that of its earliest hit bar, sigma " << tofSigmaT << " ns (FWHM " << 1000. * tofFwhmT
    << " ps); time of flight top to middle sigma " << std::sqrt(2.) * tofSigmaT << " ns\n";
  c << "TOF points: both bar ends read out; along the bar sigma vEff x sigma_t = " << tofSigmaX << " cm (vEff " << vEff
    << " cm/ns); across the bars the centre of the station's largest bar (4 cm pitch); planes at z = " << kStationZ[0] << ", "
    << kStationZ[1] << ", " << kStationZ[2] << " cm; the bottom station's bar and point in barBotM, xBot, yBot\n";
  c << "trigger (trig): the top and the middle station hit within 20 ns; track: the straight line through their points; "
       "beta = its length / (c x time of flight); tkin: kinetic energy of a muon (" << kMuonMass << " MeV) of that speed, "
       "1e5 if beta >= 1; xTgt, yTgt: the track at z = " << targetZ << " cm\n";
  c << "events kept: a middle-TOF hit (true energy above 0.2 MeV) or a primary that stopped (truth), as before the measured "
       "TOF; a trigger whose middle bar is above 0.2 MeV only as measured is missing (a bar at the threshold)\n";
  c << "seeds: hit response TRandom3(" << seed + first << "); TOF response TRandom3(1000003 x process + event + 1) per event\n";
  output.cd();   // the input files were the current directory
  TNamed("conditions", c.str().c_str()).Write();
  TNamed("detector", czt ? "czt" : "stack").Write();
  // merge mode 'f': merging skims made in parts (hadd, analysis/skim_merge.C) keeps the first value instead of adding
  // them up, and clears IsConst() if the parts differ
  TParameter<double>("sigma1MeV", sigma1MeV, 'f').Write();
  TParameter<double>("tofFwhmT", tofFwhmT, 'f').Write();
  TParameter<double>("tofFwhmE", tofFwhmE, 'f').Write();
  TParameter<double>("vEff", vEff, 'f').Write();
  TParameter<double>("targetZ", targetZ, 'f').Write();
  output.Write();
  printf("%s", c.str().c_str());
  printf("%s: %d processes, %.2f days, %lld candidates\n", out, files, liveTime / 86400., cand->GetEntries());
}
