#include "TOFmeterSteppingAction.hh"

#include "DetectorConstruction.hh"
#include "TOFmeterEventAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"

TOFmeterSteppingAction::TOFmeterSteppingAction(const Pb::DetectorConstruction* detector,
                                               TOFmeterEventAction* eventAction)
  : fDetector(detector), fEventAction(eventAction)
{}

void TOFmeterSteppingAction::UserSteppingAction(const G4Step* step)
{
  // Kinetic energy of the primary where it first enters a layer of the Compton detector
  const G4StepPoint* post = step->GetPostStepPoint();
  if (step->GetTrack()->GetTrackID() == 1 && post->GetStepStatus() == fGeomBoundary) {
    const G4VPhysicalVolume* next = post->GetPhysicalVolume();
    if (next && fDetector->IsComptonDetector(next->GetLogicalVolume()))
      fEventAction->SetComptonEntry(post->GetKineticEnergy());
  }

  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.) return;

  const G4StepPoint* pre = step->GetPreStepPoint();
  const G4VTouchable* touchable = pre->GetTouchable();
  const Pb::DetectorKind kind = fDetector->GetKind(touchable->GetVolume()->GetLogicalVolume());
  if (kind == Pb::DetectorKind::None) return;

  const G4ThreeVector pos = 0.5 * (pre->GetPosition() + step->GetPostStepPoint()->GetPosition());
  const G4double time = pre->GetGlobalTime();

  switch (kind) {
    case Pb::DetectorKind::Bar:  // bar copy number inside the station envelope
      fEventAction->AddBarHit(touchable->GetCopyNumber(1), touchable->GetCopyNumber(0), edep, time, pos);
      break;
    case Pb::DetectorKind::CameraLayer:  // copy number = CH index
      fEventAction->AddCameraHit(touchable->GetCopyNumber(0), edep, time, pos,
                                 step->GetTrack()->GetDefinition()->GetPDGEncoding(),
                                 step->GetTrack()->GetTrackID());
      break;
    case Pb::DetectorKind::StackYSO:
      fEventAction->AddStackHit(touchable->GetCopyNumber(0), edep);
      break;
    case Pb::DetectorKind::Absorber:
      fEventAction->AddAbsorberHit(edep);
      break;
    default:
      break;
  }
}
