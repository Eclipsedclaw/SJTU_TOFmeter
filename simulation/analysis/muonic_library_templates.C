// The muonic X-ray library's prompt spectra (2026-11-25, the user: unmix a compound's spectrum into the pure elements'
// spectra, as in fluorescence imaging; their choices: the prompt window, a phasor picture as a cross-check, the test
// samples Ni80Cr20, Fe50Ni50 and SS316). One pass over the jobs of one material, <band>/beam_<k>.root with
// log/job<k>.done (condor/muonic_library_submit.sh, or the steel trace run's beam band); analysis/muonic_unmixing.C fits.
//
// Entries: each pixel hit of the panel (czt_*) within promptNs (10 ns) of its muon (a tagged beam, one muon at a time:
// no accidentals), smeared by sigma1MeV x sqrt(E / MeV) (CZT 6.36991 keV: FWHM/E 1.5% at 1 MeV) and kept from
// threshold (20 keV) up, in 0.5 keV bins 0-7000 keV, as analysis/beam_steel_trace_spectrum.C. The muons an element in
// `drop` caught (e.g. "19,90,92": the steel run's K, Th, U) are left out, as if never simulated.
// Histograms (counts; the muons are in capture_z):
//   prompt                       every entry
//   prompt_x                     the entries the muonic X-rays made (czt_src bit 1)
//   prompt_half0, prompt_half1   the odd and the even jobs: two independent halves (the fit's self-test)
//   prompt_z<Z>, prompt_x_z<Z>   the entries of the muons element Z caught (cap_z; 0: none caught it)
//   capture_z, capture_z_half0, capture_z_half1   the muons per cap_z
// A mixture (a composition with more than one element) is also cut into measurements: its kept muons in job order,
// in slices of 1e5 (the first 1e7) and of 1e6 (up to 1e8):
//   slice5_<j>, slice6_<j>       their entries (complete slices only)
//   slice5_capz, slice6_capz     their muons (slice x cap_z)
// and a TNamed "settings" (the source, its composition, the muons, the smearing). Output: outFile, by default
// <band>/../analysis/templates/<band's name>.root. maxFiles > 0 reads only the first maxFiles jobs (a test).
// Run in the work folder, compiled once and then 10 at a time (rule 4), e.g. every library material and SS316 without
// its U, Th and K (each reads cap_z, czt_e, czt_t, czt_src only):
//   root -l -b -q -e '.L analysis/muonic_library_templates.C+'
//   ls -d output/spectrum_decode/muonic_library_czt/[0-9]*_* | xargs -P 10 -I{} root -l -b -q 'analysis/muonic_library_templates.C+("{}")'
//   root -l -b -q 'analysis/muonic_library_templates.C+("output/czt/beam_steel_czt_trace/beam", "output/spectrum_decode/muonic_library_czt/analysis/templates/SS316.root", "19,90,92")'
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace libtemplates {
const int kBins = 14000;   // 0.5 keV
const double kMaxKeV = 7000.;
const long kSlice5 = 100000, kSlice6 = 1000000;
const int kSlices5 = 100, kSlices6 = 100;

std::vector<int> FinishedJobs(const std::string& bandDir)   // as analysis/beam_steel_trace_spectrum.C
{
  std::vector<int> ks;
  if (void* d = gSystem->OpenDirectory(bandDir.c_str())) {
    while (const char* entry = gSystem->GetDirEntry(d)) {
      int k = 0;
      char tail[8] = "";
      if (sscanf(entry, "beam_%d.%5s", &k, tail) == 2 && std::string(tail) == "root" &&
          !gSystem->AccessPathName(Form("%s/log/job%d.done", bandDir.c_str(), k)))
        ks.push_back(k);
    }
    gSystem->FreeDirectory(d);
  }
  std::sort(ks.begin(), ks.end());
  return ks;
}

std::string MacroValue(const std::string& file, const std::string& key)   // a job macro's line "key value"
{
  std::ifstream in(file);
  std::string l;
  while (std::getline(in, l))
    if (l.rfind(key + " ", 0) == 0) return l.substr(key.size() + 1);
  return "";
}

TH1D* NewHist(const std::string& name, const std::string& what)
{
  auto h = new TH1D(name.c_str(), (what + ";energy (keV);entries / 0.5 keV").c_str(), kBins, 0., kMaxKeV);
  h->SetDirectory(nullptr);
  return h;
}
TH1D* NewCapture(const std::string& name)
{
  auto h = new TH1D(name.c_str(), "muons by the element that caught them;cap_z;muons", 100, -0.5, 99.5);
  h->SetDirectory(nullptr);
  return h;
}
}  // namespace libtemplates

void muonic_library_templates(const char* band, const char* outFile = "", const char* drop = "", int maxFiles = 0,
                              double sigma1MeV = 6.36991, double threshold = 20., double promptNs = 10.)
{
  using namespace libtemplates;
  std::string dir = band;
  while (dir.size() > 1 && dir.back() == '/') dir.pop_back();
  const size_t slash = dir.rfind('/');
  const std::string name = slash == std::string::npos ? dir : dir.substr(slash + 1);
  const std::string parent = slash == std::string::npos ? "." : dir.substr(0, slash);
  const std::string out = std::string(outFile).empty() ? parent + "/analysis/templates/" + name + ".root" : outFile;
  std::vector<int> jobs = FinishedJobs(dir);
  if (jobs.empty()) { printf("no finished jobs (beam_<k>.root with log/job<k>.done) in %s\n", dir.c_str()); return; }
  if (maxFiles > 0 && (int)jobs.size() > maxFiles) jobs.resize(maxFiles);
  const std::string composition = MacroValue(Form("%s/log/job%d.mac", dir.c_str(), jobs[0]), "/tof/target/composition");
  if (composition.empty()) { printf("no /tof/target/composition in %s/log/job%d.mac\n", dir.c_str(), jobs[0]); return; }
  const bool mixture = composition.find(',') != std::string::npos;
  std::set<int> dropZ;
  {
    std::stringstream s(drop);
    std::string z;
    while (std::getline(s, z, ','))
      if (!z.empty()) dropZ.insert(std::atoi(z.c_str()));
  }
  TRandom3 random(1 + std::hash<std::string>{}(name) % 900000000);   // a seed of its own per material

  TH1D* prompt = NewHist("prompt", name + ": prompt entries");
  TH1D* promptX = NewHist("prompt_x", name + ": prompt entries of the muonic X-rays");
  TH1D* half[2] = {NewHist("prompt_half0", name + ": prompt entries, odd jobs"),
                   NewHist("prompt_half1", name + ": prompt entries, even jobs")};
  TH1D* capture = NewCapture("capture_z");
  TH1D* captureHalf[2] = {NewCapture("capture_z_half0"), NewCapture("capture_z_half1")};
  std::map<int, std::pair<TH1D*, TH1D*>> byZ;   // cap_z -> prompt_z<Z>, prompt_x_z<Z>
  std::vector<TH1D*> s5, s6;
  TH2D *s5z = nullptr, *s6z = nullptr;
  if (mixture) {
    for (int j = 0; j < kSlices5; ++j) s5.push_back(NewHist(Form("slice5_%d", j), Form("%s: muons %ld-%ld", name.c_str(), j * kSlice5, (j + 1) * kSlice5)));
    for (int j = 0; j < kSlices6; ++j) s6.push_back(NewHist(Form("slice6_%d", j), Form("%s: muons %ld-%ld", name.c_str(), j * kSlice6, (j + 1) * kSlice6)));
    s5z = new TH2D("slice5_capz", "muons per slice of 1e5 by cap_z;slice;cap_z", kSlices5, -0.5, kSlices5 - 0.5, 100, -0.5, 99.5);
    s6z = new TH2D("slice6_capz", "muons per slice of 1e6 by cap_z;slice;cap_z", kSlices6, -0.5, kSlices6 - 0.5, 100, -0.5, 99.5);
    s5z->SetDirectory(nullptr);
    s6z->SetDirectory(nullptr);
  }
  long kept = 0, dropped = 0, entries = 0;
  TStopwatch watch;
  for (size_t f = 0; f < jobs.size(); ++f) {
    const int k = jobs[f];
    TFile file(Form("%s/beam_%d.root", dir.c_str(), k));
    auto run = dynamic_cast<TTree*>(file.Get("run"));
    auto ev = dynamic_cast<TTree*>(file.Get("events"));
    if (file.IsZombie() || !run || !ev) { printf("beam_%d.root: no run or events tree: stopping\n", k); return; }
    run->Draw("n_events:event_filter", "", "goff");
    double nRun = 0.;
    for (Long64_t i = 0; i < run->GetSelectedRows(); ++i) {
      nRun += run->GetV1()[i];
      if (run->GetV2()[i] != 0) { printf("beam_%d.root was run with /EventFilter other than 0: stopping\n", k); return; }
    }
    if (nRun != ev->GetEntries()) {   // every muon must be in the ntuple (the number of muons is the normalisation)
      printf("beam_%d.root: %.0f muons run, %lld in the events tree: stopping\n", k, nRun, ev->GetEntries());
      return;
    }
    int capZ = 0;
    std::vector<double>*e = nullptr, *t = nullptr;
    std::vector<int>* src = nullptr;
    ev->SetBranchStatus("*", 0);
    for (const char* b : {"cap_z", "czt_e", "czt_t", "czt_src"}) ev->SetBranchStatus(b, 1);
    ev->SetBranchAddress("cap_z", &capZ);
    ev->SetBranchAddress("czt_e", &e);
    ev->SetBranchAddress("czt_t", &t);
    ev->SetBranchAddress("czt_src", &src);
    const int h = (k - 1) % 2;
    for (Long64_t i = 0; i < ev->GetEntries(); ++i) {
      ev->GetEntry(i);
      if (dropZ.count(capZ)) { ++dropped; continue; }
      const long m = kept++;
      capture->Fill(capZ);
      captureHalf[h]->Fill(capZ);
      auto it = byZ.find(capZ);
      if (it == byZ.end())
        it = byZ.emplace(capZ, std::make_pair(NewHist(Form("prompt_z%d", capZ), Form("%s: prompt entries of the muons Z %d caught", name.c_str(), capZ)),
                                              NewHist(Form("prompt_x_z%d", capZ), Form("%s: their muonic X-rays' entries", name.c_str())))).first;
      TH1D *sl5 = nullptr, *sl6 = nullptr;
      if (mixture) {
        if (m < kSlice5 * kSlices5) { sl5 = s5[m / kSlice5]; s5z->Fill(m / kSlice5, capZ); }
        if (m < kSlice6 * kSlices6) { sl6 = s6[m / kSlice6]; s6z->Fill(m / kSlice6, capZ); }
      }
      for (size_t q = 0; q < e->size(); ++q) {
        const double eMeV = (*e)[q], tNs = std::max(0., (*t)[q]);
        if (eMeV <= 0. || tNs >= promptNs) continue;
        const double kev = 1000. * eMeV + random.Gaus(0., sigma1MeV * std::sqrt(eMeV));
        if (kev < threshold) continue;
        ++entries;
        prompt->Fill(kev);
        half[h]->Fill(kev);
        it->second.first->Fill(kev);
        if ((*src)[q] & 1) { promptX->Fill(kev); it->second.second->Fill(kev); }
        if (sl5) sl5->Fill(kev);
        if (sl6) sl6->Fill(kev);
      }
    }
    if ((f + 1) % 10 == 0 || f + 1 == jobs.size()) {
      printf("%s: %zu of %zu jobs, %ld muons kept, %ld dropped, %.3f prompt entries per muon, %.0f s\n", name.c_str(), f + 1,
             jobs.size(), kept, dropped, kept > 0 ? double(entries) / kept : 0., watch.RealTime());
      fflush(stdout);
      watch.Continue();
    }
  }
  if (kept == 0) { printf("%s: no muon kept\n", name.c_str()); return; }

  const size_t outSlash = out.rfind('/');
  if (outSlash != std::string::npos) gSystem->mkdir(out.substr(0, outSlash).c_str(), true);
  TFile fo(out.c_str(), "RECREATE");
  if (fo.IsZombie()) { printf("cannot write %s\n", out.c_str()); return; }
  for (TH1D* hh : {prompt, promptX, half[0], half[1], capture, captureHalf[0], captureHalf[1]}) hh->Write();
  for (auto& [z, p] : byZ) { p.first->Write(); p.second->Write(); }
  int n5 = 0, n6 = 0;
  if (mixture) {
    n5 = int(std::min<long>(kSlices5, kept / kSlice5));
    n6 = int(std::min<long>(kSlices6, kept / kSlice6));
    for (int j = 0; j < n5; ++j) s5[j]->Write();
    for (int j = 0; j < n6; ++j) s6[j]->Write();
    s5z->Write();
    s6z->Write();
  }
  TNamed("settings", Form("source %s composition %s drop %s jobs %zu muons %ld dropped %ld sigma1MeV %g threshold %g promptNs %g "
                          "slices5 %d slices6 %d", dir.c_str(), composition.c_str(), std::string(drop).empty() ? "-" : drop,
                          jobs.size(), kept, dropped, sigma1MeV, threshold, promptNs, n5, n6)).Write();
  fo.Close();
  printf("%s: %ld muons (%ld dropped), %ld prompt entries, %zu capturing elements%s -> %s\n", name.c_str(), kept, dropped,
         entries, byZ.size(), mixture ? Form(", %d + %d slices", n5, n6) : "", out.c_str());
}
