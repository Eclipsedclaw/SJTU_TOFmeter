//
// Records where and how the primary particle ended (stop and capture/decay
// volume, or leaving the world) for the "events" ntuple.
//

#ifndef TOFmeterTrackingAction_h
#define TOFmeterTrackingAction_h 1

#include "G4UserTrackingAction.hh"

class TOFmeterEventAction;

class TOFmeterTrackingAction : public G4UserTrackingAction
{
  public:
    explicit TOFmeterTrackingAction(TOFmeterEventAction*);
    ~TOFmeterTrackingAction() override = default;

    void PostUserTrackingAction(const G4Track*) override;

  private:
    TOFmeterEventAction* fEventAction;
};

#endif
