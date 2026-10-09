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
/// \file B2/B2a/include/DetectorConstruction.hh
/// \brief Definition of the B2a::DetectorConstruction class

#ifndef PbDetectorConstruction_h
#define PbDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "tls.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;

namespace Pb
{

class DetectorMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();  // 声明构造函数
    virtual ~DetectorConstruction();  // 声明析构函数

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // Set methods
    void SetMaxStep(G4double);
    void SetCheckOverlaps(G4bool);

private:
    void DefineMaterials();
    void SetVisualAttributes();  // 声明可视化属性设置方法

    // 逻辑体积指针
    G4LogicalVolume* logic_world;
    G4LogicalVolume* logic_Absorber;
    G4LogicalVolume* logic_framework_Fe;
    G4LogicalVolume* logic_framework_Al;
    G4LogicalVolume* logic_YSO;
    G4LogicalVolume* logic_Silicon;
    G4LogicalVolume* logic_CdZnTe;
    std::vector<G4LogicalVolume*> chamberComponents;

    // 材料指针
    G4Material* world_mat;
    G4Material* Absorber_mat;
    G4Material* framework_Fe_mat;
    G4Material* framework_Al_mat;
    G4Material* fChamberMaterial;
    G4Material* Silicon_mat;
    G4Material* fdetector_mat;
    G4Material* fdetector_mat1;
    G4int fNbOfChambers;
    G4Material* fTargetMaterial;

    G4UserLimits* fStepLimit;
    DetectorMessenger* fMessenger;
    G4bool fCheckOverlaps;
};

}
#endif