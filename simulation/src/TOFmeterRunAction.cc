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
// $Id: TOFmeterRunAction.cc,v 1.9 2006/06/29 17:48:16 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// @author Tsuguo Aramaki
// @date 2015 March 23 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TOFmeterRunAction.hh"
#include "TOFmeterPrimaryGeneratorAction.hh"
#include "CosmicRaySource.hh"
#include "DetectorConstruction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
/*
#include "G4ios.hh"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
*/
#include "global.h"
#include <cstdio>
#include <filesystem>
#include <fstream>

using namespace std;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TOFmeterRunAction::TOFmeterRunAction(TOFmeterPrimaryGeneratorAction* primaryGenerator)
: fPrimaryGenerator(primaryGenerator)
{
  // "run" ntuple: one row per run with the generator normalisation
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  fRunNtupleId = analysisManager->CreateNtuple("run", "Run summary and cosmic-muon normalisation");
  analysisManager->CreateNtupleIColumn("n_events");
  analysisManager->CreateNtupleIColumn("generator");   // /GeneralParticleSource: 0 gun, 1 GPS, 2 cosmic
  analysisManager->CreateNtupleIColumn("seed");
  analysisManager->CreateNtupleDColumn("flux");        // cosmic: muons /cm2/s through the plane
  analysisManager->CreateNtupleDColumn("area");        // cosmic: plane area, cm2
  analysisManager->CreateNtupleDColumn("live_time");   // cosmic: n_events / (flux area), s
  analysisManager->CreateNtupleDColumn("e_min");       // cosmic: kinetic energy range, MeV
  analysisManager->CreateNtupleDColumn("e_max");
  analysisManager->CreateNtupleDColumn("theta_max");   // cosmic: deg
  analysisManager->CreateNtupleDColumn("plane_z");     // cosmic: cm
  analysisManager->CreateNtupleSColumn("model");       // cosmic: expacs or guan
  analysisManager->CreateNtupleSColumn("particles");   // cosmic: species (expacs) or charge (guan)
  analysisManager->CreateNtupleDColumn("e_mono");      // cosmic: energy of every particle, MeV; 0 = spectrum
  analysisManager->CreateNtupleDColumn("lead_cm");     // lead layer on top, cm (0 = none)
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TOFmeterRunAction::~TOFmeterRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TOFmeterRunAction::BeginOfRunAction(const G4Run* aRun)
{
  extern global_struct global;
  char fname[1100];   // directory + file name + extension

  ((G4Run *)(aRun))->SetRunID(global.runnum);
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  if (global.outdir[0] == '\0') snprintf(global.outdir, sizeof(global.outdir), ".");
  if (global.outfile[0] == '\0') snprintf(global.outfile, sizeof(global.outfile), "test");
  std::error_code ec;
  std::filesystem::create_directories(global.outdir, ec);

  if (global.OutputFormat != 1) // ASCII step file
  {
    snprintf(fname, sizeof(fname), "%s/%s.dat", global.outdir, global.outfile );
    G4cout << "Output file: " << fname << G4endl;
    global.output.open (fname);
    if (!global.output.is_open())
    {
      G4ExceptionDescription ed;
      ed << "Cannot open " << fname << " for writing; no ASCII output this run";
      G4Exception("TOFmeterRunAction::BeginOfRunAction", "Run001", JustWarning, ed);
    }
  }
  if (global.OutputFormat != 0) // ROOT event ntuple
  {
    snprintf(fname, sizeof(fname), "%s/%s.root", global.outdir, global.outfile );
    G4cout << "Output file: " << fname << G4endl;
    G4AnalysisManager::Instance()->OpenFile(fname);
  }

  if (global.GPS == 2) fPrimaryGenerator->GetCosmicRaySource()->PrepareRun();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TOFmeterRunAction::EndOfRunAction(const G4Run* aRun)
{
	extern global_struct global;
	if (global.output.is_open()) global.output.close();

	const G4int nEvents = aRun->GetNumberOfEvent();
	const CosmicRaySource* cosmic = fPrimaryGenerator->GetCosmicRaySource();
	const G4bool isCosmic = (global.GPS == 2);
	const G4double liveTime = isCosmic ? cosmic->GetLiveTime(nEvents) : 0.;
	if (isCosmic)
	{
		G4cout << "Cosmic rays: " << nEvents << " events = " << liveTime/s << " s = "
		       << liveTime/hour << " h of equivalent live time" << G4endl;
	}

	if (global.OutputFormat != 0)
	{
		auto analysisManager = G4AnalysisManager::Instance();
		G4int col = 0;
		analysisManager->FillNtupleIColumn(fRunNtupleId, col++, nEvents);
		analysisManager->FillNtupleIColumn(fRunNtupleId, col++, global.GPS);
		analysisManager->FillNtupleIColumn(fRunNtupleId, col++, global.seed);
		analysisManager->FillNtupleDColumn(fRunNtupleId, col++, isCosmic ? cosmic->GetFlux()*cm2*s : 0.);
		analysisManager->FillNtupleDColumn(fRunNtupleId, col++, isCosmic ? cosmic->GetArea()/cm2 : 0.);
		analysisManager->FillNtupleDColumn(fRunNtupleId, col++, liveTime/s);
		analysisManager->FillNtupleDColumn(fRunNtupleId, col++, isCosmic ? cosmic->GetEnergyMin()/MeV : 0.);
		analysisManager->FillNtupleDColumn(fRunNtupleId, col++, isCosmic ? cosmic->GetEnergyMax()/MeV : 0.);
		analysisManager->FillNtupleDColumn(fRunNtupleId, col++, isCosmic ? cosmic->GetThetaMax()/deg : 0.);
		analysisManager->FillNtupleDColumn(fRunNtupleId, col++, isCosmic ? cosmic->GetPlaneCenter().z()/cm : 0.);
		analysisManager->FillNtupleSColumn(fRunNtupleId, col++, isCosmic ? cosmic->GetModel() : G4String("-"));
		analysisManager->FillNtupleSColumn(fRunNtupleId, col++, isCosmic ? cosmic->GetParticles() : G4String("-"));
		analysisManager->FillNtupleDColumn(fRunNtupleId, col++, isCosmic ? cosmic->GetMonoEnergy()/MeV : 0.);
		auto detector = dynamic_cast<const Pb::DetectorConstruction*>(
			G4RunManager::GetRunManager()->GetUserDetectorConstruction());
		analysisManager->FillNtupleDColumn(fRunNtupleId, col++, detector ? detector->GetLeadThickness()/cm : -1.);
		analysisManager->AddNtupleRow(fRunNtupleId);
		analysisManager->Write();
		analysisManager->CloseFile();
	}
	G4cout << "Run end  " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



