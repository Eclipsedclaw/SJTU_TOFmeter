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
// @author Jiancheng Zeng with Claude Opus 5.5
// @date October 2026
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "TOFmeterEventAction.hh"

#include "G4AnalysisManager.hh"
#include "G4EmProcessSubType.hh"
#include "G4GenericMessenger.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4ProcessTable.hh"
#include "G4LorentzVector.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "global.h"
#include <algorithm>
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
  analysisManager->CreateNtupleDColumn("cmp_ke");     // primary kinetic energy entering the Compton detector, -1 if never
  analysisManager->CreateNtupleDColumn("top_ke");     // ... entering a top-TOF bar, -1 if never
  analysisManager->CreateNtupleDColumn("pb_ke");      // ... leaving the lead layer, -1 if it stopped in it or missed it
  analysisManager->CreateNtupleIColumn("prim_st");    // TOF stations the primary crossed: bit 0 top, 1 middle, 2 bottom
  analysisManager->CreateNtupleIColumn("prim_top");   // first top-TOF bar of the primary, -1 if none
  analysisManager->CreateNtupleIColumn("prim_mid");   // first middle-TOF bar of the primary, -1 if none
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
  // muon tridents, one entry per interaction: vertex, kinetic energy of the muon before and after,
  // total energy of the produced pair, invariant mass of the two like-sign muons (the muon after
  // the interaction and the produced one of its charge), weight (1 / biasing factor), PDG code
  analysisManager->CreateNtupleDColumn("tri_x", fTriX);
  analysisManager->CreateNtupleDColumn("tri_y", fTriY);
  analysisManager->CreateNtupleDColumn("tri_z", fTriZ);
  analysisManager->CreateNtupleDColumn("tri_ein", fTriEin);
  analysisManager->CreateNtupleDColumn("tri_eout", fTriEout);
  analysisManager->CreateNtupleDColumn("tri_epair", fTriEpair);
  analysisManager->CreateNtupleDColumn("tri_mll", fTriMll);
  analysisManager->CreateNtupleDColumn("tri_w", fTriW);
  analysisManager->CreateNtupleIColumn("tri_pdg", fTriPdg);
  analysisManager->CreateNtupleIColumn("tri_parent", fTriParent);   // track ID of the muon, 1 = primary
  // muons made by a trident (tmu_proc 49) or by a photon conversion to mu+mu- (15, the second step
  // of a pseudotrident): PDG code, kinetic energy at creation and at the end of the track
  analysisManager->CreateNtupleIColumn("tmu_pdg", fTmuPdg);
  analysisManager->CreateNtupleIColumn("tmu_proc", fTmuProc);
  analysisManager->CreateNtupleDColumn("tmu_e0", fTmuE0);
  analysisManager->CreateNtupleDColumn("tmu_e1", fTmuE1);
  analysisManager->CreateNtupleIColumn("tmu_parent", fTmuParent);   // track ID of the muon that made it
  analysisManager->CreateNtupleIColumn("tmu_st", fTmuSt);           // TOF stations crossed, as prim_st
  analysisManager->CreateNtupleIColumn("tmu_top", fTmuTop);         // first top-TOF bar, -1 if none
  analysisManager->CreateNtupleIColumn("tmu_mid", fTmuMid);         // first middle-TOF bar, -1 if none
  // photons of mu- captures at rest, one entry each: energy at creation, energy when it first left the
  // volume it was born in (-1 if it never did, i.e. it was absorbed there), creation time and place,
  // source (0 muonic X-ray of the atomic cascade, 1 gamma of the nuclear capture), track ID of the muon
  analysisManager->CreateNtupleDColumn("mux_e", fMuxE);
  analysisManager->CreateNtupleDColumn("mux_exit", fMuxExit);
  analysisManager->CreateNtupleDColumn("mux_t", fMuxT);
  analysisManager->CreateNtupleDColumn("mux_x", fMuxX);
  analysisManager->CreateNtupleDColumn("mux_y", fMuxY);
  analysisManager->CreateNtupleDColumn("mux_z", fMuxZ);
  analysisManager->CreateNtupleIColumn("mux_src", fMuxSrc);
  analysisManager->CreateNtupleIColumn("mux_parent", fMuxParent);
  analysisManager->FinishNtuple();

  fMessenger = std::make_unique<G4GenericMessenger>(this, "/tof/tridents/", "Muon trident events");
  fMessenger->DeclareMethod("keepForVis", &TOFmeterEventAction::KeepTridents,
                            "Keep events for /vis/reviewKeptEvents: 0 none, 1 events with a trident, 2 events with "
                            "a detectable trident (its three muons all crossed the top and middle TOF); resets the count");
  fMessenger->DeclareProperty("maxKept", fMaxKept, "Stop the run when this many events have been kept (default 20)");
  fMessenger->DeclareProperty("search", fSearch,
                              "true: /run/beamOn N (also the Run button of the Qt toolbar) runs until N events are kept "
                              "as keepForVis says (1 if it is 0), at most searchLimit muons; false: normal runs");
  fMessenger->DeclareProperty("searchLimit", fSearchLimit, "Most muons one search may take (default 1000000)");
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
  fComptonEntryKE = -1.;
  fTopEntryKE = -1.;
  fLeadExitKE = -1.;
  fEndPos = G4ThreeVector();
  fEndTime = 0.;
  fEndKinE = 0.;
  fEndProcess = "none";
  fEndVolume = "none";
  fCamLayer.clear(); fCamPdg.clear(); fCamTrack.clear();
  fCamE.clear(); fCamT.clear(); fCamX.clear(); fCamY.clear(); fCamZ.clear();
  fTriX.clear(); fTriY.clear(); fTriZ.clear(); fTriEin.clear(); fTriEout.clear();
  fTriEpair.clear(); fTriMll.clear(); fTriW.clear(); fTriPdg.clear();
  fTmuPdg.clear(); fTmuProc.clear(); fTmuE0.clear(); fTmuE1.clear();
  fTriParent.clear(); fTmuParent.clear(); fTmuSt.clear(); fTmuTop.clear(); fTmuMid.clear();
  fTofTracks.clear();
  fCapturePhotons.clear();
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

void TOFmeterEventAction::AddTrident(const G4Step& step, const std::vector<const G4Track*>& secondaries)
{
  const G4StepPoint* post = step.GetPostStepPoint();
  const G4int pdg = step.GetTrack()->GetDefinition()->GetPDGEncoding();
  // the incoming muon after the interaction, plus the produced muon of the same charge
  G4LorentzVector likeSign(post->GetMomentum(), post->GetTotalEnergy());
  G4double pairEnergy = 0.;
  G4double weight = post->GetWeight();
  G4int nMuons = 0;
  for (const G4Track* t : secondaries) {
    const G4int q = t->GetDefinition()->GetPDGEncoding();
    if (q != 13 && q != -13) continue;
    ++nMuons;
    pairEnergy += t->GetTotalEnergy();
    weight = t->GetWeight();
    if (q == pdg) likeSign += G4LorentzVector(t->GetMomentum(), t->GetTotalEnergy());
  }
  if (nMuons == 0) return;  // interaction rejected by the integral method: nothing was produced

  fTriX.push_back(post->GetPosition().x() / cm);
  fTriY.push_back(post->GetPosition().y() / cm);
  fTriZ.push_back(post->GetPosition().z() / cm);
  fTriEin.push_back((post->GetKineticEnergy() + pairEnergy) / MeV);  // the nucleus takes no energy
  fTriEout.push_back(post->GetKineticEnergy() / MeV);
  fTriEpair.push_back(pairEnergy / MeV);
  fTriMll.push_back(likeSign.m() / MeV);
  fTriW.push_back(weight);
  fTriPdg.push_back(pdg);
  fTriParent.push_back(step.GetTrack()->GetTrackID());
}

void TOFmeterEventAction::AddTridentMuon(G4int pdg, G4int creatorSubType, G4double kinE0, G4double kinE1,
                                         G4int trackID, G4int parentID)
{
  const auto it = fTofTracks.find(trackID);
  const TOFTrack tof = (it == fTofTracks.end()) ? TOFTrack() : it->second;
  fTmuParent.push_back(parentID);
  fTmuSt.push_back(tof.mask);
  fTmuTop.push_back(tof.topBar);
  fTmuMid.push_back(tof.midBar);
  fTmuPdg.push_back(pdg);
  fTmuProc.push_back(creatorSubType);
  fTmuE0.push_back(kinE0 / MeV);
  fTmuE1.push_back(kinE1 / MeV);
}

void TOFmeterEventAction::AddCapturePhoton(const G4Track& track, G4int source)
{
  CapturePhoton& photon = fCapturePhotons[track.GetTrackID()];
  photon.energy = track.GetKineticEnergy();
  photon.time = track.GetGlobalTime();
  photon.position = track.GetPosition();
  photon.source = source;
  photon.parent = track.GetParentID();
}

void TOFmeterEventAction::SetCapturePhotonExit(G4int trackID, G4double kinE)
{
  const auto it = fCapturePhotons.find(trackID);
  if (it != fCapturePhotons.end() && it->second.exitEnergy < 0.) it->second.exitEnergy = kinE;
}

void TOFmeterEventAction::MarkTOFTrack(G4int trackID, G4int station, G4int bar)
{
  TOFTrack& t = fTofTracks[trackID];
  t.mask |= 1 << station;
  if (station == 0 && t.topBar < 0) t.topBar = bar;
  if (station == 1 && t.midBar < 0) t.midBar = bar;
}

void TOFmeterEventAction::KeepTridents(G4int mode)
{
  fKeepTridents = mode;
  fKept = 0;
}

// /run/beamOn n with /tof/tridents/search true (TOFmeterRunManager): keep events as keepForVis
// says (every trident if it is 0) and stop the run at the n-th kept one
G4int TOFmeterEventAction::BeginTridentSearch(G4int n)
{
  if (!G4ProcessTable::GetProcessTable()->FindProcess("muToMuonPairProd", "mu-")) {
    G4cout << "Trident search: no muon trident process (muToMuonPairProd); start TOFmeter with "
              "--preinit macro/physics_tridents.mac"
           << G4endl;
    return 0;
  }
  if (fKeepTridents <= 0) fKeepTridents = 1;
  fSearchTarget = n;
  fKept = 0;
  G4cout << "Trident search: looking for " << n << (n == 1 ? " event" : " events") << " with "
         << (fKeepTridents == 2 ? "a detectable trident" : "a trident") << " (at most " << fSearchLimit
         << " muons)" << G4endl;
  return fSearchLimit;
}

void TOFmeterEventAction::EndTridentSearch(G4int muons)
{
  if (fSearchTarget <= 0) return;
  G4cout << "Trident search: found " << fKept << " of " << fSearchTarget << " in " << muons << " muons";
  if (fKept < fSearchTarget) G4cout << " (/tof/tridents/searchLimit sets the most muons per search)";
  G4cout << G4endl;
  fSearchTarget = 0;
}

// A trident of the primary whose three muons (the primary after it and the produced pair) all
// crossed a bar of the top TOF and a bar of the middle TOF
G4bool TOFmeterEventAction::DetectableTrident() const
{
  const G4int topAndMiddle = 0x3;
  if ((fPrimSt & topAndMiddle) != topAndMiddle) return false;
  if (std::find(fTriParent.begin(), fTriParent.end(), 1) == fTriParent.end()) return false;
  G4int produced = 0;
  for (std::size_t i = 0; i < fTmuProc.size(); ++i) {
    if (fTmuProc[i] != fMuonPairProdByCharged || fTmuParent[i] != 1) continue;
    if ((fTmuSt[i] & topAndMiddle) != topAndMiddle) return false;
    ++produced;
  }
  return produced >= 2;
}

void TOFmeterEventAction::KeepForVisualisation(const G4Event* evt)
{
  const G4int maxKept = fSearchTarget > 0 ? fSearchTarget : fMaxKept;   // a search keeps its own number
  if (fKeepTridents <= 0 || fTriX.empty() || fKept >= maxKept) return;
  if (fKeepTridents == 2 && !DetectableTrident()) return;
  G4EventManager::GetEventManager()->KeepTheCurrentEvent();
  ++fKept;

  // one line per kept event, to go with its picture
  G4String top = std::to_string(fPrimTop), middle = std::to_string(fPrimMid);
  for (std::size_t i = 0; i < fTmuProc.size(); ++i) {
    if (fTmuProc[i] != fMuonPairProdByCharged) continue;
    top += " " + std::to_string(fTmuTop[i]);
    middle += " " + std::to_string(fTmuMid[i]);
  }
  G4cout << "Kept trident event " << fKept << "/" << maxKept << ": event " << evt->GetEventID() << ", muon "
         << fTriEin[0] / 1000. << " -> " << fTriEout[0] / 1000. << " GeV, pair " << fTriEpair[0] / 1000.
         << " GeV at (" << fTriX[0] << ", " << fTriY[0] << ", " << fTriZ[0] << ") cm, like-sign mass "
         << fTriMll[0] << " MeV; bars of the primary and the pair, top TOF: " << top << ", middle TOF: " << middle
         << G4endl;
  if (fKept >= maxKept) {
    // the run stops here; its label in the viewer (/vis/scene/add/eventID) then gives the muons used
    G4RunManager* runManager = G4RunManager::GetRunManager();
    if (G4Run* run = runManager->GetNonConstCurrentRun()) run->SetNumberOfEventToBeProcessed(evt->GetEventID() + 1);
    runManager->AbortRun(true);
  }
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

    // TOF stations and bars of the primary
    const auto primary = fTofTracks.find(1);
    fPrimSt = (primary == fTofTracks.end()) ? 0 : primary->second.mask;
    fPrimTop = (primary == fTofTracks.end()) ? -1 : primary->second.topBar;
    fPrimMid = (primary == fTofTracks.end()) ? -1 : primary->second.midBar;
    KeepForVisualisation(evt);

    extern global_struct global;
    if (global.OutputFormat == 0) return; // ASCII output only

    const G4bool detectorEnergy = !fBars.empty() || fCamEdep[0] > 0. || fCamEdep[1] > 0. || fCamEdep[2] > 0.
                                  || fStackEdep[0] > 0. || fStackEdep[1] > 0. || fComptonEntryKE >= 0.;
    if (global.EventFilter == 1) // skip muons that only crossed air
    {
      const G4bool stoppedInSetup = fEndVolume != "OutOfWorld" && fEndVolume != "WorldPhys";
      if (!detectorEnergy && fAbsorberEdep <= 0. && !stoppedInSetup) return;
    }
    if (global.EventFilter == 2 && !detectorEnergy) return; // keep only TOF or Compton-detector energy
    if (global.EventFilter == 3 && fTriX.empty() && fTmuPdg.empty()) return; // keep only muon-pair events
    if (global.EventFilter == 4) // keep only energy in both a top-TOF and a middle-TOF bar
    {
      G4bool top = false, middle = false;
      for (const auto& [key, sum] : fBars) {
        if (sum.edep <= 0.) continue;
        top = top || key / 100 == 0;
        middle = middle || key / 100 == 1;
      }
      if (!top || !middle) return;
    }
    if (global.EventFilter == 5 && fCapturePhotons.empty()) return; // keep only events with mu- capture photons

    fMuxE.clear(); fMuxExit.clear(); fMuxT.clear(); fMuxX.clear(); fMuxY.clear(); fMuxZ.clear();
    fMuxSrc.clear(); fMuxParent.clear();
    for (const auto& [trackID, photon] : fCapturePhotons) {
      fMuxE.push_back(photon.energy / MeV);
      fMuxExit.push_back(photon.exitEnergy < 0. ? -1. : photon.exitEnergy / MeV);
      fMuxT.push_back(photon.time / ns);
      fMuxX.push_back(photon.position.x() / cm);
      fMuxY.push_back(photon.position.y() / cm);
      fMuxZ.push_back(photon.position.z() / cm);
      fMuxSrc.push_back(photon.source);
      fMuxParent.push_back(photon.parent);
    }

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
    analysisManager->FillNtupleDColumn(fNtupleId, col++, fComptonEntryKE < 0. ? -1. : fComptonEntryKE / MeV);
    analysisManager->FillNtupleDColumn(fNtupleId, col++, fTopEntryKE < 0. ? -1. : fTopEntryKE / MeV);
    analysisManager->FillNtupleDColumn(fNtupleId, col++, fLeadExitKE < 0. ? -1. : fLeadExitKE / MeV);
    analysisManager->FillNtupleIColumn(fNtupleId, col++, fPrimSt);
    analysisManager->FillNtupleIColumn(fNtupleId, col++, fPrimTop);
    analysisManager->FillNtupleIColumn(fNtupleId, col++, fPrimMid);
    analysisManager->AddNtupleRow(fNtupleId);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
