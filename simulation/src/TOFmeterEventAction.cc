//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: TOFmeterEventAction.cc,v 1.11 2006/06/29 17:48:05 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// @author Tsuguo Aramaki
// @date 2015 March 23
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "TOFmeterEventAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "global.h"
#include <unistd.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TOFmeterEventAction::TOFmeterEventAction()
{
  // "events" ntuple: one row per event (MeV, ns, cm). Scalar columns come
  // first and are filled in this order in EndOfEventAction.
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  fNtupleId = analysisManager->CreateNtuple("events", "TOF-meter event summary");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->CreateNtupleIColumn("pdg");        // primary
  analysisManager->CreateNtupleDColumn("ekin");       // primary kinetic energy
  analysisManager->CreateNtupleDColumn("costh");      // cos(zenith) of the primary, 1 = straight down
  analysisManager->CreateNtupleDColumn("phi");
  analysisManager->CreateNtupleDColumn("x0");
  analysisManager->CreateNtupleDColumn("y0");
  analysisManager->CreateNtupleDColumn("z0");
  analysisManager->CreateNtupleDColumn("end_x");      // where the primary stopped or left
  analysisManager->CreateNtupleDColumn("end_y");
  analysisManager->CreateNtupleDColumn("end_z");
  analysisManager->CreateNtupleDColumn("end_t");
  analysisManager->CreateNtupleDColumn("end_ke");
  analysisManager->CreateNtupleSColumn("end_proc");   // e.g. muMinusCaptureAtRest, Decay, Transportation
  analysisManager->CreateNtupleSColumn("end_vol");    // physical volume, OutOfWorld if it left
  analysisManager->CreateNtupleDColumn("e_ch0");      // Compton camera, summed per layer
  analysisManager->CreateNtupleDColumn("e_ch1");
  analysisManager->CreateNtupleDColumn("e_ch2");
  analysisManager->CreateNtupleDColumn("t_ch0");      // first deposit, -1 if none
  analysisManager->CreateNtupleDColumn("t_ch1");
  analysisManager->CreateNtupleDColumn("t_ch2");
  analysisManager->CreateNtupleIColumn("cam_prim");   // bit k: primary deposited in CHk
  analysisManager->CreateNtupleDColumn("e_stack0");   // YSO layers of the Si/YSO stack
  analysisManager->CreateNtupleDColumn("e_stack1");
  analysisManager->CreateNtupleDColumn("e_absorber");
  // bars with energy: station 0-2, bar 0-8, energy, first time, energy-weighted position
  analysisManager->CreateNtupleIColumn("bar_st", fBarStation);
  analysisManager->CreateNtupleIColumn("bar_id", fBarIndex);
  analysisManager->CreateNtupleDColumn("bar_e", fBarEdep);
  analysisManager->CreateNtupleDColumn("bar_t", fBarTime);
  analysisManager->CreateNtupleDColumn("bar_x", fBarX);
  analysisManager->CreateNtupleDColumn("bar_y", fBarY);
  analysisManager->CreateNtupleDColumn("bar_z", fBarZ);
  // every step with energy in the camera: layer (CH index), particle, track, energy, time, position
  analysisManager->CreateNtupleIColumn("cam_layer", fCamLayer);
  analysisManager->CreateNtupleIColumn("cam_pdg", fCamPdg);
  analysisManager->CreateNtupleIColumn("cam_trk", fCamTrack);
  analysisManager->CreateNtupleDColumn("cam_e", fCamE);
  analysisManager->CreateNtupleDColumn("cam_t", fCamT);
  analysisManager->CreateNtupleDColumn("cam_x", fCamX);
  analysisManager->CreateNtupleDColumn("cam_y", fCamY);
  analysisManager->CreateNtupleDColumn("cam_z", fCamZ);
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TOFmeterEventAction::~TOFmeterEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TOFmeterEventAction::BeginOfEventAction(const G4Event* evt)
{
  extern global_struct global;
  global.eventID = evt->GetEventID();

  fBars.clear();
  for (G4int i = 0; i < 3; ++i) { fCamEdep[i] = 0.; fCamTime[i] = -1.; }
  fCamPrimary = 0;
  fStackEdep[0] = fStackEdep[1] = 0.;
  fAbsorberEdep = 0.;
  fEndPos = G4ThreeVector();
  fEndTime = 0.;
  fEndKinE = 0.;
  fEndProcess = "none";
  fEndVolume = "none";
  fCamLayer.clear(); fCamPdg.clear(); fCamTrack.clear();
  fCamE.clear(); fCamT.clear(); fCamX.clear(); fCamY.clear(); fCamZ.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TOFmeterEventAction::AddBarHit(G4int station, G4int bar, G4double edep, G4double time,
                                    const G4ThreeVector& pos)
{
  BarSum& sum = fBars[100 * station + bar];
  sum.edep += edep;
  if (sum.time < 0. || time < sum.time) sum.time = time;
  sum.ex += edep * pos.x();
  sum.ey += edep * pos.y();
  sum.ez += edep * pos.z();
}

void TOFmeterEventAction::AddCameraHit(G4int layer, G4double edep, G4double time,
                                       const G4ThreeVector& pos, G4int pdg, G4int trackID)
{
  if (layer < 0 || layer > 2) return;
  fCamEdep[layer] += edep;
  if (fCamTime[layer] < 0. || time < fCamTime[layer]) fCamTime[layer] = time;
  if (trackID == 1) fCamPrimary |= (1 << layer);
  fCamLayer.push_back(layer);
  fCamPdg.push_back(pdg);
  fCamTrack.push_back(trackID);
  fCamE.push_back(edep / MeV);
  fCamT.push_back(time / ns);
  fCamX.push_back(pos.x() / cm);
  fCamY.push_back(pos.y() / cm);
  fCamZ.push_back(pos.z() / cm);
}

void TOFmeterEventAction::AddStackHit(G4int layer, G4double edep)
{
  if (layer >= 0 && layer < 2) fStackEdep[layer] += edep;
}

void TOFmeterEventAction::SetPrimaryEnd(const G4ThreeVector& pos, G4double time, G4double kinE,
                                        const G4String& process, const G4String& volume)
{
  fEndPos = pos;
  fEndTime = time;
  fEndKinE = kinE;
  fEndProcess = process;
  fEndVolume = volume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TOFmeterEventAction::EndOfEventAction(const G4Event* evt)
{
    G4int event_id = evt->GetEventID();
    if (event_id < 10 || event_id%1000 == 0) G4cout << ">>> Event " << evt->GetEventID() << G4endl;

    extern global_struct global;
    if (global.OutputFormat == 0) return; // ASCII output only

    const G4bool detectorEnergy = !fBars.empty() || fCamEdep[0] > 0. || fCamEdep[1] > 0. || fCamEdep[2] > 0.
                                  || fStackEdep[0] > 0. || fStackEdep[1] > 0.;
    if (global.EventFilter == 1) // skip muons that only crossed air
    {
      const G4bool stoppedInSetup = fEndVolume != "OutOfWorld" && fEndVolume != "WorldPhys";
      if (!detectorEnergy && fAbsorberEdep <= 0. && !stoppedInSetup) return;
    }
    if (global.EventFilter == 2 && !detectorEnergy) return; // keep only TOF or Compton-detector energy

    fBarStation.clear(); fBarIndex.clear(); fBarEdep.clear(); fBarTime.clear();
    fBarX.clear(); fBarY.clear(); fBarZ.clear();
    for (const auto& [key, sum] : fBars) {
      if (sum.edep <= 0.) continue;
      fBarStation.push_back(key / 100);
      fBarIndex.push_back(key % 100);
      fBarEdep.push_back(sum.edep / MeV);
      fBarTime.push_back(sum.time / ns);
      fBarX.push_back(sum.ex / sum.edep / cm);
      fBarY.push_back(sum.ey / sum.edep / cm);
      fBarZ.push_back(sum.ez / sum.edep / cm);
    }

    G4int pdg = 0;
    G4double ekin = 0., costh = 0., phi = 0.;
    G4ThreeVector vertex;
    if (const G4PrimaryVertex* pv = evt->GetPrimaryVertex(0)) {
      vertex = pv->GetPosition();
      if (const G4PrimaryParticle* primary = pv->GetPrimary(0)) {
        pdg = primary->GetPDGcode();
        ekin = primary->GetKineticEnergy();
        costh = -primary->GetMomentumDirection().z();
        phi = primary->GetMomentumDirection().phi();
      }
    }

    auto analysisManager = G4AnalysisManager::Instance();
    G4int col = 0;
    analysisManager->FillNtupleIColumn(fNtupleId, col++, event_id);
    analysisManager->FillNtupleIColumn(fNtupleId, col++, pdg);
    analysisManager->FillNtupleDColumn(fNtupleId, col++, ekin / MeV);
    analysisManager->FillNtupleDColumn(fNtupleId, col++, costh);
    analysisManager->FillNtupleDColumn(fNtupleId, col++, phi);
    analysisManager->FillNtupleDColumn(fNtupleId, col++, vertex.x() / cm);
    analysisManager->FillNtupleDColumn(fNtupleId, col++, vertex.y() / cm);
    analysisManager->FillNtupleDColumn(fNtupleId, col++, vertex.z() / cm);
    analysisManager->FillNtupleDColumn(fNtupleId, col++, fEndPos.x() / cm);
    analysisManager->FillNtupleDColumn(fNtupleId, col++, fEndPos.y() / cm);
    analysisManager->FillNtupleDColumn(fNtupleId, col++, fEndPos.z() / cm);
    analysisManager->FillNtupleDColumn(fNtupleId, col++, fEndTime / ns);
    analysisManager->FillNtupleDColumn(fNtupleId, col++, fEndKinE / MeV);
    analysisManager->FillNtupleSColumn(fNtupleId, col++, fEndProcess);
    analysisManager->FillNtupleSColumn(fNtupleId, col++, fEndVolume);
    for (G4int i = 0; i < 3; ++i) analysisManager->FillNtupleDColumn(fNtupleId, col++, fCamEdep[i] / MeV);
    for (G4int i = 0; i < 3; ++i)
      analysisManager->FillNtupleDColumn(fNtupleId, col++, fCamTime[i] < 0. ? -1. : fCamTime[i] / ns);
    analysisManager->FillNtupleIColumn(fNtupleId, col++, fCamPrimary);
    analysisManager->FillNtupleDColumn(fNtupleId, col++, fStackEdep[0] / MeV);
    analysisManager->FillNtupleDColumn(fNtupleId, col++, fStackEdep[1] / MeV);
    analysisManager->FillNtupleDColumn(fNtupleId, col++, fAbsorberEdep / MeV);
    analysisManager->AddNtupleRow(fNtupleId);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
