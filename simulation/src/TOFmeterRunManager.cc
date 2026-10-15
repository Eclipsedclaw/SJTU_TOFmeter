#include "TOFmeterRunManager.hh"

#include "TOFmeterEventAction.hh"

void TOFmeterRunManager::BeamOn(G4int nEvents, const char* macroFile, G4int nSelect)
{
  auto eventAction = dynamic_cast<TOFmeterEventAction*>(userEventAction);
  if (nEvents <= 0 || !eventAction || !eventAction->TridentSearch()) {
    G4RunManager::BeamOn(nEvents, macroFile, nSelect);
    return;
  }
  // a search: the event action aborts the run when it has kept nEvents events
  const G4int muons = eventAction->BeginTridentSearch(nEvents);
  if (muons > 0) G4RunManager::BeamOn(muons, macroFile, nSelect);
  eventAction->EndTridentSearch(muons > 0 ? numberOfEventProcessed : 0);
}
