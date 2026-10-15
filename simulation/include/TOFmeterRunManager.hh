//
// G4RunManager with a search mode for muon tridents: with /tof/tridents/search true, /run/beamOn N
// (also the Run button of the Qt toolbar, which sends /run/beamOn 1) does not stop after N muons
// but when N events with a trident have been kept for the viewer (/tof/tridents/keepForVis, see
// TOFmeterEventAction), or after /tof/tridents/searchLimit muons. /run/beamOn 0 and all runs with
// the search off are those of G4RunManager. Live time and output count the muons actually run.
//

#ifndef TOFmeterRunManager_h
#define TOFmeterRunManager_h 1

#include "G4RunManager.hh"

class TOFmeterRunManager : public G4RunManager
{
  public:
    void BeamOn(G4int nEvents, const char* macroFile = nullptr, G4int nSelect = -1) override;
};

#endif
