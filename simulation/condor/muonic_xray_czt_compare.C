// Compares two outputs of the same process of the CZT run (macro/muonic_xray_czt.mac with the same seed, run by
// scripts/muonic_xray_czt.sh or as an HTCondor job of condor/muonic_xray_czt_submit.sh), e.g. the cluster's test
// (2 x 10^5 muons) with the shared server's process 1 (5 x 10^6 muons, the same seed: its first 2 x 10^5 muons are the
// test's): the events of the muons both ran (eventID below the smaller n_events of the run trees), row by row: the
// muon's energy and end process and every CZT hit's pixel, energy and time. All identical: the same physics, bit for
// bit; another Geant4 version, compiler or CPU can change single events while the distributions still agree. Then the
// whole files' distributions are compared as two independent samples: the events written per muon (Poisson), the
// muon's end processes (chi2 of the table; the rare ones together), the muon's energy and the CZT energy of the event
// (Kolmogorov-Smirnov); consistent if all within 3 sigma or p > 0.01.
//   root -l -b -q 'condor/muonic_xray_czt_compare.C("output/muonic_xray_czt_condor_test/slow/stk_1.root", "server/slow/stk_1.root")'

#include <TFile.h>
#include <TLeaf.h>
#include <TMath.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
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

// Kolmogorov-Smirnov probability that the values of a quantity of the two runs' events come from one distribution
template <class F>
double KS(const Run& ra, const Run& rb, F value)
{
  std::vector<double> a, b;
  for (const auto& [id, row] : ra.rows) a.push_back(value(row));
  for (const auto& [id, row] : rb.rows) b.push_back(value(row));
  std::sort(a.begin(), a.end());
  std::sort(b.begin(), b.end());
  return TMath::KolmogorovTest(a.size(), a.data(), b.size(), b.data(), "");
}

// the two runs as independent samples; true if they agree within the statistics
bool CompareDistributions(const Run& ra, const Run& rb)
{
  const double na = ra.rows.size(), nb = rb.rows.size();
  if (na < 1 || nb < 1) {
    printf("no events in one of the files: no distributions to compare\n");
    return false;
  }
  const double fa = na / ra.muons, fb = nb / rb.muons;
  const double z = (fa - fb) / std::sqrt(fa / ra.muons + fb / rb.muons);
  printf("events written per 10^6 muons: a %.0f +- %.0f, b %.0f +- %.0f (%+.1f sigma)\n", 1e6 * fa,
         1e6 * std::sqrt(fa / ra.muons), 1e6 * fb, 1e6 * std::sqrt(fb / rb.muons), z);
  // end processes: the 2 x k table, a process expected fewer than 5 times in a run counted with the other rare ones
  std::map<std::string, std::pair<double, double>> n;
  for (const auto& [id, row] : ra.rows) n[row.endProc].first += 1.;
  for (const auto& [id, row] : rb.rows) n[row.endProc].second += 1.;
  std::vector<std::pair<std::string, std::pair<double, double>>> cells;
  std::pair<double, double> rare{0., 0.};
  for (const auto& [proc, c] : n) {
    const double share = (c.first + c.second) / (na + nb);
    if (share * std::min(na, nb) < 5.) {
      rare.first += c.first;
      rare.second += c.second;
    } else {
      cells.push_back({proc, c});
    }
  }
  if (rare.first + rare.second > 0.) cells.push_back({"(the rare ones)", rare});
  double chi2 = 0.;
  printf("the muon's end:            a        b\n");
  for (const auto& [proc, c] : cells) {
    const double share = (c.first + c.second) / (na + nb);
    chi2 += std::pow(c.first - share * na, 2) / (share * na) + std::pow(c.second - share * nb, 2) / (share * nb);
    printf("  %-22s %5.1f%%   %5.1f%%\n", proc.c_str(), 100. * c.first / na, 100. * c.second / nb);
  }
  const int ndf = static_cast<int>(cells.size()) - 1;
  const double pEnd = ndf > 0 ? TMath::Prob(chi2, ndf) : 1.;
  printf("  chi2 %.1f for %d degrees of freedom: p = %.3f\n", chi2, ndf, pEnd);
  const double pE = KS(ra, rb, [](const Row& r) { return r.ekin; });
  const double pCzt = KS(ra, rb, [](const Row& r) { return std::accumulate(r.e.begin(), r.e.end(), 0.); });
  printf("Kolmogorov-Smirnov: the muon's energy p = %.3f, the event's CZT energy p = %.3f\n", pE, pCzt);
  return std::fabs(z) < 3. && pEnd > 0.01 && pE > 0.01 && pCzt > 0.01;
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
  if (inA + inB > 0 && different == 0 && onlyA == 0 && onlyB == 0) {
    printf("the same events, bit for bit\n");
    return;
  }
  printf("%s: the distributions of the whole files as independent samples\n",
         inA + inB == 0 ? "no events to compare" : "not the same events");
  printf("%s\n", CompareDistributions(ra, rb) ? "consistent within the statistics" : "NOT consistent: look closer");
}
