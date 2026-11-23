// Where the beam's mu- stop in the plate (2026-11-23, the muonic X-ray library: every material is tuned to the steel
// run's mean depth, 8.2424 mm). The depth is measured below the plate's top face, 33.5 cm (macro/beam_steel_lyso.mac:
// the 1 cm plate at 32.5-33.5 cm, the beam coming down from 33.55 cm). Only the primaries that came to rest in the
// plate count (end_ke 0, end_vol Target).
// Prints one line, which scripts/muonic_library_depth.sh reads:
//   depth <mean, mm> err <error of the mean> rms <mm> median <mm> stopped <muons at rest in the plate> events <all>
// files: one ROOT file or a pattern (TChain::Add), e.g. a library run's beam/beam_*.root. In the work folder:
//   root -l -b -q 'analysis/muonic_library_depth.C("output/muonic_library_czt/26_Fe/beam/beam_1.root")'
#include <TChain.h>
#include <TH1D.h>

#include <cmath>
#include <cstdio>

void muonic_library_depth(const char* files, double plateTopCm = 33.5)
{
  TChain chain("events");
  if (chain.Add(files) == 0) {
    printf("no file matches %s\n", files);
    return;
  }
  const Long64_t events = chain.GetEntries();
  TH1D h("depth", "", 10000, 0., 10.);   // mm, 1 um bins: the plate is 10 mm thick
  chain.Draw(Form("10. * (%.6f - end_z) >> depth", plateTopCm), "end_ke == 0 && end_vol == \"Target\"", "goff");
  const double n = h.GetEntries();
  if (n < 2) {
    printf("depth nan err nan rms nan median nan stopped %.0f events %lld\n", n, events);
    return;
  }
  double half = 0.5, median = 0.;
  h.GetQuantiles(1, &median, &half);
  printf("depth %.4f err %.4f rms %.4f median %.4f stopped %.0f events %lld\n", h.GetMean(), h.GetRMS() / std::sqrt(n),
         h.GetRMS(), median, n, events);
}
