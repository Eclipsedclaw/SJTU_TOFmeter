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
// @author Tsuguo Aramaki
// @date 2015 March 23
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#ifndef TOFmeterEventAction_h
#define TOFmeterEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include <map>
#include <vector>

class G4Event;

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
    void AddStackHit(G4int layer, G4double edep);
    void AddAbsorberHit(G4double edep) { fAbsorberEdep += edep; }
    // first entry of the primary into the Compton detector (later entries are ignored)
    void SetComptonEntry(G4double kinE) { if (fComptonEntryKE < 0.) fComptonEntryKE = kinE; }
    void SetPrimaryEnd(const G4ThreeVector& pos, G4double time, G4double kinE,
                       const G4String& process, const G4String& volume);

  private:
    struct BarSum { G4double edep = 0., time = -1., ex = 0., ey = 0., ez = 0.; };

    G4int fNtupleId;
    std::map<G4int, BarSum> fBars;    // key 100 * station + bar
    G4double fCamEdep[3];
    G4double fCamTime[3];
    G4int fCamPrimary;                // bit k set if the primary deposited energy in CHk
    G4double fStackEdep[2];
    G4double fAbsorberEdep;
    G4double fComptonEntryKE;         // -1 if the primary never entered the Compton detector
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
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
