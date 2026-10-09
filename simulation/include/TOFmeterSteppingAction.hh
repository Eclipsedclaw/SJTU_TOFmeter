//
// Event-level scoring: energy deposits in the active volumes registered by
// Pb::DetectorConstruction (bars, camera layers, stack YSO, Pb absorber) are
// summed into TOFmeterEventAction. Independent of /TrackType and /TrackEdep,
// which only filter the step-level ASCII output of TOFmeterDetectorSD.
//

#ifndef TOFmeterSteppingAction_h
#define TOFmeterSteppingAction_h 1

#include "G4UserSteppingAction.hh"

namespace Pb { class DetectorConstruction; }
class TOFmeterEventAction;

class TOFmeterSteppingAction : public G4UserSteppingAction
{
  public:
    TOFmeterSteppingAction(const Pb::DetectorConstruction*, TOFmeterEventAction*);
    ~TOFmeterSteppingAction() override = default;

    void UserSteppingAction(const G4Step*) override;

  private:
    const Pb::DetectorConstruction* fDetector;
    TOFmeterEventAction* fEventAction;
};

#endif
