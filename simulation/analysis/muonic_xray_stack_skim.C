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

namespace {
// where the primary ended: 0 target, 1 Al plates, 2 stack Si, 3 stack YSO, 4 TOF bars, 5 lead, 6 steel, 7 air or out
int VolumeCode(const std::string& v)
{
  if (v.rfind("Target", 0) == 0) return 0;
  if (v.rfind("FrameworkAl", 0) == 0) return 1;
  if (v.rfind("SiliconLayer", 0) == 0) return 2;
  if (v.rfind("YSOLayer", 0) == 0) return 3;
  if (v.rfind("Scintillator", 0) == 0) return 4;
  if (v.rfind("Absorber", 0) == 0 || v.rfind("LeadWall", 0) == 0) return 5;
  if (v.rfind("FrameworkFe", 0) == 0 || v.rfind("SupportBeam", 0) == 0) return 6;
  return 7;
}
}  // namespace

void muonic_xray_stack_skim(const char* dir = "output/muonic_xray_stack", int first = 1, int last = 385,
                            const char* out = "skim.root", unsigned seed = 4357)
{
  TFile output(Form("%s/%s", dir, out), "RECREATE");
  auto cand = new TTree("cand", "add-back X-ray candidates");
  auto info = new TTree("info", "exposure of the skimmed processes");
  int proc, event, n, layers, topLayer, src, pdg, endVol, endProc, tofTop, tofMid, barTop, barMid;
  float e, emax, rmax, t, tTop, tMid, eTop, eMid, botFirst, ekin, endT, endZ, pbKe;
  cand->Branch("proc", &proc); cand->Branch("event", &event);
  cand->Branch("e", &e);               // add-back energy (keV, smeared)
  cand->Branch("n", &n);               // hits added up
  cand->Branch("layers", &layers);     // bit k: a hit in layer k (0, 2, 4 Si; 1, 3 YSO; from the bottom)
  cand->Branch("topLayer", &topLayer); // highest layer with a hit
  cand->Branch("emax", &emax);         // largest hit (keV)
  cand->Branch("rmax", &rmax);         // largest |x| or |y| of its hits (cm, smeared): how close to the stack edge
  cand->Branch("t", &t);               // time of its first hit (ns)
  cand->Branch("src", &src);           // truth: stk_src bits of its hits (1 X-ray of the target)
  cand->Branch("tofTop", &tofTop); cand->Branch("tTop", &tTop); cand->Branch("barTop", &barTop); cand->Branch("eTop", &eTop);
  cand->Branch("tofMid", &tofMid); cand->Branch("tMid", &tMid); cand->Branch("barMid", &barMid); cand->Branch("eMid", &eMid);
  cand->Branch("botFirst", &botFirst); // earliest bottom-TOF hit (ns), -1 if none
  cand->Branch("pdg", &pdg); cand->Branch("ekin", &ekin);   // truth: the primary
  cand->Branch("endVol", &endVol);     // truth: where it ended (VolumeCode)
  cand->Branch("endProc", &endProc);   // truth: 1 mu- captured at rest, 2 decay at rest, 0 else
  cand->Branch("endT", &endT); cand->Branch("endZ", &endZ);
  cand->Branch("pbKe", &pbKe);         // truth: its kinetic energy leaving the top lead layer (MeV), -1 if it went round it
  double liveTime = 0.;
  long long muons = 0;
  int files = 0;
  info->Branch("live_time", &liveTime); info->Branch("muons", &muons); info->Branch("files", &files);

  TRandom3 random(seed + first);   // reproducible, different for every part
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
    int eventID = 0, primary = 0;
    double kin = 0., tEnd = 0., zEnd = 0., keEnd = 0., kePb = -1.;
    char vol[128] = "", process[128] = "";
    std::vector<int>*layer = nullptr, *hitSrc = nullptr, *barSt = nullptr, *barId = nullptr;
    std::vector<double>*he = nullptr, *ht = nullptr, *hx = nullptr, *hy = nullptr, *be = nullptr, *bt = nullptr;
    events->SetBranchStatus("*", 0);
    for (const char* c : {"eventID", "pdg", "ekin", "end_t", "end_z", "end_ke", "end_vol", "end_proc", "pb_ke", "stk_layer",
                          "stk_src", "stk_e", "stk_t", "stk_x", "stk_y", "bar_st", "bar_id", "bar_e", "bar_t"})
      events->SetBranchStatus(c, 1);
    events->SetBranchAddress("eventID", &eventID); events->SetBranchAddress("pdg", &primary);
    events->SetBranchAddress("ekin", &kin); events->SetBranchAddress("end_t", &tEnd); events->SetBranchAddress("end_z", &zEnd);
    events->SetBranchAddress("pb_ke", &kePb);
    events->SetBranchAddress("end_ke", &keEnd); events->SetBranchAddress("end_vol", vol); events->SetBranchAddress("end_proc", process);
    events->SetBranchAddress("stk_layer", &layer); events->SetBranchAddress("stk_src", &hitSrc);
    events->SetBranchAddress("stk_e", &he); events->SetBranchAddress("stk_t", &ht);
    events->SetBranchAddress("stk_x", &hx); events->SetBranchAddress("stk_y", &hy);
    events->SetBranchAddress("bar_st", &barSt); events->SetBranchAddress("bar_id", &barId);
    events->SetBranchAddress("bar_e", &be); events->SetBranchAddress("bar_t", &bt);
    struct Hit { double t, e, x, y; int layer, src; };
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
      const std::string procName = process;
      endProc = procName == "muMinusCaptureAtRest" ? 1 : (procName.rfind("Decay", 0) == 0 && keEnd <= 0.) ? 2 : 0;
      // a muon that stops in the target crosses the middle TOF (4 cm above it), often not the top one: keep every
      // event with a middle-TOF hit, and the stopped muons (truth)
      if (!tofMid && endProc == 0) continue;
      // the stack hits with the response
      hits.clear();
      for (std::size_t h = 0; h < he->size(); ++h) {
        const double measured = (*he)[h] + random.Gaus(0., 0.04 * std::sqrt((*he)[h]));   // MeV
        if (measured * 1000. < 10.) continue;
        hits.push_back({(*ht)[h], measured * 1000., (*hx)[h] + random.Gaus(0., 0.2), (*hy)[h] + random.Gaus(0., 0.2),
                        (*layer)[h], (*hitSrc)[h]});
      }
      std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) { return a.t < b.t; });
      proc = k; event = eventID; pdg = primary; ekin = kin; endVol = VolumeCode(vol); endT = tEnd; endZ = zEnd; pbKe = kePb;
      for (std::size_t h = 0; h < hits.size();) {   // add-back entries
        t = hits[h].t;
        e = 0.; n = 0; layers = 0; topLayer = -1; emax = 0.; rmax = 0.; src = 0;
        for (; h < hits.size() && hits[h].t < t + 100.; ++h) {
          e += hits[h].e; ++n;
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
  output.Write();
  printf("%s: %d processes, %.2f days, %lld candidates\n", out, files, liveTime / 86400., cand->GetEntries());
}
