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
//
// TOF-meter geometry used by TOFmeter.cc.

#ifndef PbDetectorConstruction_h
#define PbDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include <memory>
#include <unordered_map>
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GenericMessenger;
class G4UImessenger;

namespace Pb
{

// Role of a volume, used for event-level scoring (TOFmeterSteppingAction); StackSi is passive
enum class DetectorKind { None, Bar, CameraLayer, StackYSO, StackSi, Absorber };

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();  // 声明构造函数
    virtual ~DetectorConstruction();  // 声明析构函数

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // Rebuild the geometry from the current parameters (/update)
    void UpdateGeometry();

    // Set methods
    void SetMaxStep(G4double);
    void SetCheckOverlaps(G4bool);

    DetectorKind GetKind(const G4LogicalVolume* lv) const;
    // Layers of the Compton detector: the Si/YSO stack and the camera
    G4bool IsComptonDetector(const G4LogicalVolume* lv) const;

    // Thickness of the lead layer on top, 0 when it is not built
    G4double GetLeadThickness() const { return fBuildAbsorber ? fLeadThickness : 0.; }

    // Passive boxes (/tof/blocks/add, /tof/blocks/clear); spec is
    // "name material x y z dx dy dz unit" with centre and full sizes
    void AddBlock(const G4String& spec);
    void ClearBlocks();

    // Every placed volume with its world-frame centre, axis-aligned size and
    // z range (/tof/det/print)
    void PrintGeometry();

private:
    void DefineMaterials();
    void DefineCommands();
    void SetVisualAttributes();  // 声明可视化属性设置方法
    void ConstructComptonCamera();
    void ConstructTargetAndDegrader();
    void ConstructBlocks();
    G4Material* FindMaterial(const G4String& name) const;

    G4VPhysicalVolume* fWorldPhys;

    // 逻辑体积指针
    G4LogicalVolume* logic_world;
    G4LogicalVolume* logic_Absorber;
    G4LogicalVolume* logic_framework_Fe;
    G4LogicalVolume* logic_framework_Al;
    G4LogicalVolume* logic_YSO;
    G4LogicalVolume* logic_Silicon;
    G4LogicalVolume* logic_CdZnTe;
    G4LogicalVolume* logic_CameraYSO;
    G4LogicalVolume* logic_CameraLYSO;
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
    G4Material* fLYSO_mat;
    G4int fNbOfChambers;
    G4Material* fTargetMaterial;

    G4UserLimits* fStepLimit;
    G4bool fCheckOverlaps;

    // Parts of the original setup that can be switched off (/tof/det/...)
    G4bool fBuildStations;
    G4bool fBuildFrame;
    G4bool fBuildStack;
    G4bool fBuildAbsorber;
    G4bool fBuildLeadWalls;
    G4double fLeadThickness;        // lead layer on the steel plate (/tof/det/leadThickness)
    G4String fWorldMaterialName;
    G4ThreeVector fWorldSize;

    // Compton camera (/tof/camera/...)
    G4bool fBuildCamera;
    G4ThreeVector fCameraPosition;  // centre of CH2, the entrance layer
    G4double fCameraRotationX;
    G4double fCameraSpacing21;      // CH2-CH1, centre to centre
    G4double fCameraSpacing10;      // CH1-CH0, centre to centre

    // Target and degrader (/tof/target/..., /tof/degrader/...)
    G4bool fBuildTarget;
    G4String fTargetMaterialName;
    G4ThreeVector fTargetSize;
    G4ThreeVector fTargetPosition;
    G4bool fBuildDegrader;
    G4String fDegraderMaterialName;
    G4ThreeVector fDegraderSize;
    G4ThreeVector fDegraderPosition;

    // Free-form boxes, e.g. lead shielding bricks (/tof/blocks/...)
    struct Block {
        G4String name;
        G4String material;
        G4ThreeVector position;
        G4ThreeVector size;
    };
    std::vector<Block> fBlocks;

    std::unordered_map<const G4LogicalVolume*, DetectorKind> fKinds;
    std::vector<std::unique_ptr<G4GenericMessenger>> fMessengers;
    std::unique_ptr<G4UImessenger> fBlockMessenger;
};

}
#endif
