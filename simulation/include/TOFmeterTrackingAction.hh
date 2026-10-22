//
// Records where and how the primary particle ended (stop and capture/decay
// volume, or leaving the world) for the "events" ntuple, and the photons of
// mu- captures at rest (muonic X-rays and nuclear-capture gammas) when they start.
//

#ifndef TOFmeterTrackingAction_h
#define TOFmeterTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"

class TOFmeterEventAction;

class TOFmeterTrackingAction : public G4UserTrackingAction
{
  public:
    explicit TOFmeterTrackingAction(TOFmeterEventAction*);
    ~TOFmeterTrackingAction() override = default;

    void PreUserTrackingAction(const G4Track*) override;
    void PostUserTrackingAction(const G4Track*) override;

  private:
    TOFmeterEventAction* fEventAction;
    G4int fCascadeModelID;   // creator model of the muonic X-rays (atomic cascade)
    G4int fCaptureModelID;   // creator model of the nuclear-capture products
};

#endif
