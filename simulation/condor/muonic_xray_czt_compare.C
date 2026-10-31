// Compares two outputs of the same process of the CZT run (macro/muonic_xray_czt.mac with the same seed, run by
// scripts/muonic_xray_czt.sh or as an HTCondor job of condor/muonic_xray_czt_submit.sh), e.g. the cluster's test
// (2 x 10^5 muons) with the shared server's process 1 (5 x 10^6 muons, the same seed: its first 2 x 10^5 muons are the
// test's): the events of the muons both ran (eventID below the smaller n_events of the run trees), row by row: the
// muon's energy and end process and every CZT hit's pixel, energy and time. All identical: the same physics, bit for
// bit; another Geant4 version, compiler or CPU can change single events while the distributions still agree.
//   root -l -b -q 'condor/muonic_xray_czt_compare.C("output/muonic_xray_czt_condor_test/slow/stk_1.root", "server/slow/stk_1.root")'

#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>

#include <map>
#include <string>
#include <vector>

namespace {
struct Row {
  double ekin;
  std::string endProc;
  std::vector<int> pix;
  std::vector<double> e, t;
  bool operator==(const Row& o) const
  {
    return ekin == o.ekin && endProc == o.endProc && pix == o.pix && e == o.e && t == o.t;
  }
};

struct Run {
  long long muons = 0;
  long long seed = 0;
  std::map<int, Row> rows;  // by eventID
};

// the run tree's muons and seed and every event; false if the file lacks a tree or column, or an entry cannot be read
bool Read(const char* name, Run& r)
{
  TFile file(name);
  auto run = file.IsZombie() ? nullptr : dynamic_cast<TTree*>(file.Get("run"));
  auto events = run ? dynamic_cast<TTree*>(file.Get("events")) : nullptr;
  if (!events || run->GetEntries() < 1) {
    printf("%s: no run and events trees\n", name);
    return false;
  }
  run->GetEntry(0);
  r.muons = static_cast<long long>(run->GetLeaf("n_events")->GetValue());
  r.seed = static_cast<long long>(run->GetLeaf("seed")->GetValue());
  int id = 0;
  double ekin = 0.;
  char endProc[256] = "";  // a process name (TLeafC)
  std::vector<int>* pix = nullptr;
  std::vector<double>* e = nullptr;
  std::vector<double>* t = nullptr;
  auto proc = events->GetLeaf("end_proc");
  if (!proc || proc->GetMaximum() >= static_cast<int>(sizeof endProc) || events->SetBranchAddress("eventID", &id) < 0
      || events->SetBranchAddress("ekin", &ekin) < 0 || events->SetBranchAddress("end_proc", endProc) < 0
      || events->SetBranchAddress("czt_pix", &pix) < 0 || events->SetBranchAddress("czt_e", &e) < 0
      || events->SetBranchAddress("czt_t", &t) < 0) {
    printf("%s: the events tree lacks eventID, ekin, end_proc or the czt_* columns\n", name);
    return false;
  }
  for (Long64_t i = 0; i < events->GetEntries(); ++i) {
    if (events->GetEntry(i) <= 0) {
      printf("%s: cannot read entry %lld\n", name, i);
      return false;
    }
    r.rows[id] = Row{ekin, endProc, *pix, *e, *t};
  }
  return true;
}
}  // namespace

void muonic_xray_czt_compare(const char* a, const char* b)
{
  Run ra, rb;
  if (!Read(a, ra) || !Read(b, rb)) return;
  printf("a: %s, seed %lld, %lld muons, %zu events written\n", a, ra.seed, ra.muons, ra.rows.size());
  printf("b: %s, seed %lld, %lld muons, %zu events written\n", b, rb.seed, rb.muons, rb.rows.size());
  if (ra.seed != rb.seed) printf("the seeds differ: no event should be the same\n");
  const long long n = ra.muons < rb.muons ? ra.muons : rb.muons;
  long long same = 0, different = 0, onlyA = 0, onlyB = 0, inA = 0, inB = 0;
  int shown = 0;
  for (const auto& [id, row] : ra.rows) {
    if (id >= n) continue;
    ++inA;
    auto other = rb.rows.find(id);
    if (other == rb.rows.end()) {
      ++onlyA;
      continue;
    }
    if (row == other->second) {
      ++same;
      continue;
    }
    ++different;
    if (shown++ < 5)
      printf("  event %d differs: muon %.9g / %.9g MeV, end %s / %s, %zu / %zu CZT hits\n", id, row.ekin,
             other->second.ekin, row.endProc.c_str(), other->second.endProc.c_str(), row.e.size(), other->second.e.size());
  }
  for (const auto& [id, row] : rb.rows)
    if (id < n) {
      ++inB;
      if (!ra.rows.count(id)) ++onlyB;
    }
  printf("the %lld muons both ran: %lld and %lld events written; identical %lld, different %lld, only in a %lld, only in b %lld\n",
         n, inA, inB, same, different, onlyA, onlyB);
  if (inA + inB == 0) printf("no events to compare\n");
  else if (different == 0 && onlyA == 0 && onlyB == 0) printf("the same events, bit for bit\n");
  else printf("not the same events: compare their distributions instead\n");
}
