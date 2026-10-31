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
// $Id: TOFmeterEventAction.hh,v 1.8 2006/06/29 17:47:35 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// @author Jiancheng Zeng with Claude Opus 5.5
// @date October 2026
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#ifndef TOFmeterEventAction_h
#define TOFmeterEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include <map>
#include <memory>
#include <unordered_map>
#include <vector>

class G4Event;
class G4GenericMessenger;
class G4Step;
class G4Track;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TOFmeterEventAction : public G4UserEventAction
{
  public:
    TOFmeterEventAction();
   ~TOFmeterEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

    // Event-level scoring, filled by TOFmeterSteppingAction / TOFmeterTrackingAction
    // and written to the "events" ntuple when /OutputFormat is 1 or 2
    void AddBarHit(G4int station, G4int bar, G4double edep, G4double time, const G4ThreeVector& pos);
    void AddCameraHit(G4int layer, G4double edep, G4double time, const G4ThreeVector& pos,
                      G4int pdg, G4int trackID);
    // energy in a layer of the Si/YSO stack: a Si layer (copy 0-2) or a YSO layer (copy 0-1), from the bottom
    void AddStackHit(G4bool silicon, G4int copy, G4double edep, G4double time, const G4ThreeVector& pos,
                     G4int trackID);
    // energy in a pixel of the CZT plane (/tof/czt/...), pixel = ix * pixels + iy
    void AddCZTHit(G4int pixel, G4double edep, G4double time, const G4ThreeVector& pos, G4int trackID);
    // where the energy of a track comes from, for the stack hits (stk_src): a muonic X-ray (source 0) or a
    // nuclear-capture gamma (source 1) of a mu- capture, else what its parent came from (tracking action)
    void SetTrackOrigin(const G4Track& track, G4int captureSource);
    static constexpr G4int kTargetXray = 1, kOtherXray = 2, kCaptureGamma = 4, kOtherOrigin = 8;
    void AddAbsorberHit(G4double edep) { fAbsorberEdep += edep; }
    // first entry of the primary into the Compton detector (later entries are ignored)
    void SetComptonEntry(G4double kinE) { if (fComptonEntryKE < 0.) fComptonEntryKE = kinE; }
    // first entry of the primary into a bar of the top TOF, and its first exit from the lead layer
    void SetTopTOFEntry(G4double kinE) { if (fTopEntryKE < 0.) fTopEntryKE = kinE; }
    void SetLeadExit(G4double kinE) { if (fLeadExitKE < 0.) fLeadExitKE = kinE; }
    // Muon tridents (/physics_lists/em/MuonToMuons): one entry per interaction (stepping action),
    // one per muon made by a trident or by a photon conversion to mu+mu- (tracking action)
    void AddTrident(const G4Step& step, const std::vector<const G4Track*>& secondaries);
    void AddTridentMuon(G4int pdg, G4int creatorSubType, G4double kinE0, G4double kinE1,
                        G4int trackID, G4int parentID);
    // A muon (track) deposited energy in a bar of a TOF station (0 top, 1 middle, 2 bottom)
    void MarkTOFTrack(G4int trackID, G4int station, G4int bar);
    // /tof/tridents/keepForVis: keep events for the visualisation (/vis/reviewKeptEvents),
    // 0 none, 1 with a trident, 2 with a detectable trident; resets the count of kept events
    void KeepTridents(G4int mode);
    // /tof/tridents/search: /run/beamOn N runs until N events are kept (TOFmeterRunManager)
    G4bool TridentSearch() const { return fSearch; }
    G4int BeginTridentSearch(G4int n);   // most muons to run, 0 without the trident process
    void EndTridentSearch(G4int muons);
    void SetPrimaryEnd(const G4ThreeVector& pos, G4double time, G4double kinE,
                       const G4String& process, const G4String& volume);
    // Photons of a mu- capture at rest (muMinusCaptureAtRest): muonic X-rays of the atomic cascade
    // (source 0) and gammas of the nuclear capture (source 1), recorded when they start (tracking
    // action), and the energy with which they first leave the volume they were born in (stepping action)
    void AddCapturePhoton(const G4Track& track, G4int source);
    void SetCapturePhotonExit(G4int trackID, G4double kinE);

  private:
    struct BarSum { G4double edep = 0., time = -1., ex = 0., ey = 0., ez = 0.; };
    struct TOFTrack { G4int mask = 0, topBar = -1, midBar = -1; };  // stations crossed, first bars
    struct CapturePhoton {
      G4double energy = 0., exitEnergy = -1., time = 0.;   // exitEnergy -1: never left its volume
      G4ThreeVector position;
      G4int source = 0, parent = 0;
    };
    struct StackDeposit {   // one step with energy in the stack (layer 0, 2, 4 Si and 1, 3 YSO, from the bottom) or in
      G4int layer = 0;      // the CZT plane (layer = pixel)
      G4double edep = 0., time = 0.;
      G4ThreeVector position;
      G4int origin = kOtherOrigin;
    };
    // hits of the stack or the CZT plane: per layer (pixel), the deposits within fStackHitWindow of the first one
    void MakeHits(std::vector<StackDeposit>& deposits, std::vector<G4int>& layer, std::vector<G4int>& src,
                  std::vector<G4double>& e, std::vector<G4double>& t, std::vector<G4double>& x,
                  std::vector<G4double>& y, std::vector<G4double>& z) const;

    G4bool DetectableTrident() const;
    void KeepForVisualisation(const G4Event*);

    G4int fNtupleId;
    std::map<G4int, BarSum> fBars;    // key 100 * station + bar
    G4double fCamEdep[3];
    G4double fCamTime[3];
    G4int fCamPrimary;                // bit k set if the primary deposited energy in CHk
    G4double fStackEdep[2];           // YSO layers of the stack
    G4double fSiEdep[3];              // Si layers of the stack
    std::vector<StackDeposit> fStackDeposits;
    std::vector<StackDeposit> fCZTDeposits;
    std::unordered_map<G4int, G4int> fOrigin;   // kTargetXray ... kOtherOrigin, by track ID
    G4double fStackHitWindow;         // /tof/stack/hitWindow
    G4double fAbsorberEdep;
    G4double fComptonEntryKE;         // -1 if the primary never entered the Compton detector
    G4double fTopEntryKE = -1.;       // -1 if it never entered a top-TOF bar
    G4double fLeadExitKE = -1.;       // -1 if it never left the lead layer (stopped in it, or missed it)
    G4ThreeVector fEndPos;
    G4double fEndTime;
    G4double fEndKinE;
    G4String fEndProcess;
    G4String fEndVolume;

    // vector columns of the "events" ntuple
    std::vector<G4int> fBarStation, fBarIndex;
    std::vector<G4double> fBarEdep, fBarTime, fBarX, fBarY, fBarZ;
    std::vector<G4int> fCamLayer, fCamPdg, fCamTrack;
    std::vector<G4double> fCamE, fCamT, fCamX, fCamY, fCamZ;
    std::vector<G4int> fStkLayer, fStkSrc;
    std::vector<G4double> fStkE, fStkT, fStkX, fStkY, fStkZ;
    std::vector<G4int> fCztPix, fCztSrc;
    std::vector<G4double> fCztE, fCztT, fCztX, fCztY, fCztZ;
    std::vector<G4double> fTriX, fTriY, fTriZ, fTriEin, fTriEout, fTriEpair, fTriMll, fTriW;
    std::vector<G4int> fTriPdg;
    std::vector<G4int> fTmuPdg, fTmuProc;
    std::vector<G4double> fTmuE0, fTmuE1;
    std::vector<G4int> fTriParent, fTmuParent, fTmuSt, fTmuTop, fTmuMid;
    std::map<G4int, TOFTrack> fTofTracks;   // muons that crossed a TOF bar, by track ID
    std::map<G4int, CapturePhoton> fCapturePhotons;   // photons of mu- captures, by track ID
    std::vector<G4double> fMuxE, fMuxExit, fMuxT, fMuxX, fMuxY, fMuxZ;
    std::vector<G4int> fMuxSrc, fMuxParent;
    G4int fPrimSt = 0, fPrimTop = -1, fPrimMid = -1;

    G4int fKeepTridents = 0;
    G4int fMaxKept = 20;
    G4int fKept = 0;
    G4bool fSearch = false;
    G4int fSearchLimit = 1000000;
    G4int fSearchTarget = 0;          // events to keep in the running search, 0 if none
    std::unique_ptr<G4GenericMessenger> fMessenger;
    std::unique_ptr<G4GenericMessenger> fStackMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
