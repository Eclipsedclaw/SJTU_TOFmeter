#include "TOFmeterSteppingAction.hh"

#include "DetectorConstruction.hh"
#include "TOFmeterEventAction.hh"

#include "G4EmProcessSubType.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"

TOFmeterSteppingAction::TOFmeterSteppingAction(const Pb::DetectorConstruction* detector,
                                               TOFmeterEventAction* eventAction)
  : fDetector(detector), fEventAction(eventAction)
{}

void TOFmeterSteppingAction::UserSteppingAction(const G4Step* step)
{
  // Kinetic energy of the primary where it first enters a layer of the Compton detector or a bar of
  // the top TOF (station copy number 0), and where it first leaves the lead layer
  const G4StepPoint* post = step->GetPostStepPoint();
  if (step->GetTrack()->GetTrackID() == 1 && post->GetStepStatus() == fGeomBoundary) {
    const G4VPhysicalVolume* next = post->GetPhysicalVolume();
    const G4LogicalVolume* nextLV = next ? next->GetLogicalVolume() : nullptr;
    const G4LogicalVolume* thisLV = step->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume();
    if (nextLV && fDetector->IsComptonDetector(nextLV)) fEventAction->SetComptonEntry(post->GetKineticEnergy());
    if (nextLV && fDetector->GetKind(nextLV) == Pb::DetectorKind::Bar && post->GetTouchable()->GetCopyNumber(1) == 0)
      fEventAction->SetTopTOFEntry(post->GetKineticEnergy());
    if (nextLV != thisLV && fDetector->GetKind(thisLV) == Pb::DetectorKind::Absorber)
      fEventAction->SetLeadExit(post->GetKineticEnergy());
  }

  // Photons of mu- captures: energy when they first leave the volume they were born in (the event
  // action keeps it only for the capture photons it knows)
  if (post->GetStepStatus() == fGeomBoundary && step->GetTrack()->GetDefinition() == G4Gamma::Gamma()
      && step->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()
           == step->GetTrack()->GetLogicalVolumeAtVertex())
    fEventAction->SetCapturePhotonExit(step->GetTrack()->GetTrackID(), post->GetKineticEnergy());

  // Muon tridents: mu+mu- pair production by the muon of this step (muToMuonPairProd)
  if (const G4VProcess* proc = post->GetProcessDefinedStep()) {
    if (proc->GetProcessType() == fElectromagnetic && proc->GetProcessSubType() == fMuonPairProdByCharged) {
      if (const auto* secondaries = step->GetSecondaryInCurrentStep()) fEventAction->AddTrident(*step, *secondaries);
    }
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
    case Pb::DetectorKind::Bar: {  // bar copy number inside the station envelope
      fEventAction->AddBarHit(touchable->GetCopyNumber(1), touchable->GetCopyNumber(0), edep, time, pos);
      const G4int pdg = step->GetTrack()->GetDefinition()->GetPDGEncoding();
      if (pdg == 13 || pdg == -13)  // which muon crossed which station (trident detection)
        fEventAction->MarkTOFTrack(step->GetTrack()->GetTrackID(), touchable->GetCopyNumber(1),
                                   touchable->GetCopyNumber(0));
      break;
    }
    case Pb::DetectorKind::CameraLayer:  // copy number = CH index
      fEventAction->AddCameraHit(touchable->GetCopyNumber(0), edep, time, pos,
                                 step->GetTrack()->GetDefinition()->GetPDGEncoding(),
                                 step->GetTrack()->GetTrackID());
      break;
    case Pb::DetectorKind::StackSi:   // copy number = layer of its kind, from the bottom
    case Pb::DetectorKind::StackYSO:
      fEventAction->AddStackHit(kind == Pb::DetectorKind::StackSi, touchable->GetCopyNumber(0), edep, time, pos,
                                step->GetTrack()->GetTrackID());
      break;
    case Pb::DetectorKind::Absorber:
      fEventAction->AddAbsorberHit(edep);
      break;
    default:
      break;
  }
}
