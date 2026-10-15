#include "TOFmeterTrackingAction.hh"

#include "TOFmeterEventAction.hh"

#include "G4EmProcessSubType.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

TOFmeterTrackingAction::TOFmeterTrackingAction(TOFmeterEventAction* eventAction)
  : fEventAction(eventAction)
{}

void TOFmeterTrackingAction::PostUserTrackingAction(const G4Track* track)
{
  // muons made by a trident or by a photon conversion to mu+mu-: energy at creation and at the end
  if (const G4VProcess* creator = track->GetCreatorProcess()) {
    const G4int sub = creator->GetProcessSubType();
    const G4int pdg = track->GetDefinition()->GetPDGEncoding();
    if (creator->GetProcessType() == fElectromagnetic && (pdg == 13 || pdg == -13)
        && (sub == fMuonPairProdByCharged || sub == fGammaConversionToMuMu)) {
      fEventAction->AddTridentMuon(pdg, sub, track->GetVertexKineticEnergy(), track->GetKineticEnergy(),
                                   track->GetTrackID(), track->GetParentID());
    }
  }

  if (track->GetTrackID() != 1) return;

  const G4StepPoint* last = track->GetStep()->GetPostStepPoint();
  const G4VProcess* process = last->GetProcessDefinedStep();
  // No volume after the last step means the primary left the world
  const G4VPhysicalVolume* volume = last->GetPhysicalVolume();
  fEventAction->SetPrimaryEnd(track->GetPosition(), track->GetGlobalTime(), track->GetKineticEnergy(),
                              process ? process->GetProcessName() : G4String("none"),
                              volume ? volume->GetName() : G4String("OutOfWorld"));
}
