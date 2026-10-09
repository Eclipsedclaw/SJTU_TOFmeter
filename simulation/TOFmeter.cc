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
// $Id: exampleN02.cc,v 1.12 2006/06/29 17:47:25 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// @author Tsuguo Aramaki
// @date 2015 March 23
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "TOFmeterDetectorMessenger.hh"
#include "TOFmeterPhysicsList.hh"
#include "TOFmeterPrimaryGeneratorAction.hh"
#include "TOFmeterRunAction.hh"
#include "TOFmeterEventAction.hh"
#include "TOFmeterSteppingAction.hh"
#include "TOFmeterTrackingAction.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "global.h"

#include <cstdlib>
#include <cstring>

#ifndef TOF_GEANT4_DATA_DIR
#define TOF_GEANT4_DATA_DIR ""
#endif

global_struct global;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Geant4 datasets: GEANT4_DATA_DIR from the environment wins; otherwise use the
  // directory given to CMake as TOF_GEANT4_DATA_DIR (if any)
  if (!std::getenv("GEANT4_DATA_DIR") && std::strlen(TOF_GEANT4_DATA_DIR) > 0)
    setenv("GEANT4_DATA_DIR", TOF_GEANT4_DATA_DIR, 0);

  // Run manager
  //
  G4RunManager * runManager = new G4RunManager;

  // User Initialization classes (mandatory)
  //
  auto detector = new Pb::DetectorConstruction;
  runManager->SetUserInitialization(detector);
  // /OutputDirectory, /OutputFile, /gun/seed, /update, ... (seed: /gun/seed or /random/setSeeds)
  auto detectorMessenger = new TOFmeterDetectorMessenger(detector);
  //

  // Reference physics list: QGSP_BERT unless the PHYSLIST environment variable
  // names another one, e.g. PHYSLIST=FTFP_BERT_EMZ (EM option 4).
  // Avoid *_HP lists and Shielding until the G4NDL dataset is installed.
  const char* listName = std::getenv("PHYSLIST");
  G4PhysListFactory factory;
  G4VModularPhysicsList* physics = factory.GetReferencePhysList(listName ? listName : "QGSP_BERT");
  if (!physics) return 1;

	physics->SetDefaultCutValue(1.0*mm) ;
	runManager->SetUserInitialization(physics);

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // User Action classes
  //
  auto gen_action = new TOFmeterPrimaryGeneratorAction();
  runManager->SetUserAction(gen_action);
  //
  auto event_action = new TOFmeterEventAction;
  runManager->SetUserAction(event_action);
  //
  G4UserRunAction* run_action = new TOFmeterRunAction(gen_action);
  runManager->SetUserAction(run_action);
  //
  runManager->SetUserAction(new TOFmeterSteppingAction(detector, event_action));
  runManager->SetUserAction(new TOFmeterTrackingAction(event_action));
  // Initialize G4 kernel
  //
  runManager->Initialize();
      
  // Get the pointer to the User Interface manager
  //
  G4UImanager * UI = G4UImanager::GetUIpointer();  
  G4UIExecutive* session = 0;
  if(argc==1)  // Define (G)UI terminal for interactive mode
  { 
    // G4UIterminal is a (dumb) terminal
    //
    session = new G4UIExecutive(argc, argv);
    UI->ApplyCommand("/control/execute");    
    session->SessionStart();
    delete session;
  }
  else   // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  delete visManager;

  delete detectorMessenger;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

