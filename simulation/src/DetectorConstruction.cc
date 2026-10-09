#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "ChamberSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

#include "G4Element.hh"
#include "G4GenericMessenger.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4RunManager.hh"
#include "G4ThreeVector.hh"

namespace Pb
{

// 构造函数
DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction(),
      logic_world(nullptr),
      logic_Absorber(nullptr),
      logic_framework_Fe(nullptr),
      logic_framework_Al(nullptr),
            logic_YSO(nullptr),
      logic_Silicon(nullptr),
      logic_CdZnTe(nullptr),
      world_mat(nullptr),
      Absorber_mat(nullptr),
      framework_Fe_mat(nullptr),
      framework_Al_mat(nullptr),
      fChamberMaterial(nullptr),
      Silicon_mat(nullptr),
      fdetector_mat(nullptr),
            fdetector_mat1(nullptr),
      fNbOfChambers(0),
      fTargetMaterial(nullptr),
      fStepLimit(nullptr),
      fMessenger(nullptr),
      fCheckOverlaps(true)
{
    // 在构造函数中调用DefineMaterials确保材料被创建
    DefineMaterials();
    
    // 创建messenger
    fMessenger = new DetectorMessenger(this);
}

// 析构函数
DetectorConstruction::~DetectorConstruction()
{
    delete fMessenger;
}

void DetectorConstruction::DefineMaterials()
{
    G4NistManager* nist = G4NistManager::Instance();
    
    G4cout << "\n=== Defining Materials ===" << G4endl;
    
    // 使用FindOrBuildMaterial的返回值 - 真空材料
    world_mat = nist->FindOrBuildMaterial("G4_Galactic");
    if (!world_mat) {
        G4Exception("DetectorConstruction::DefineMaterials", 
                   "Material001", FatalException,
                   "Failed to create G4_Galactic material for world!");
    }
      
    fTargetMaterial = nist->FindOrBuildMaterial("G4_AIR");
    Absorber_mat = nist->FindOrBuildMaterial("G4_Pb");
    framework_Fe_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    framework_Al_mat = nist->FindOrBuildMaterial("G4_Al");
    fChamberMaterial = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    Silicon_mat = nist->FindOrBuildMaterial("G4_Si");
    
    // 验证标准材料是否创建成功
    G4cout << "Standard materials created:" << G4endl;
    G4cout << "  World material (AIR): " << (world_mat ? world_mat->GetName() : "NULL") << G4endl;
    G4cout << "  Absorber (Pb): " << (Absorber_mat ? Absorber_mat->GetName() : "NULL") << G4endl;
    G4cout << "  Silicon: " << (Silicon_mat ? Silicon_mat->GetName() : "NULL") << G4endl;

    // 修正：正确创建CdZnTe材料
    G4cout << "\nCreating CdZnTe material..." << G4endl;
    
    // 修正元素名称和比例
    G4Element* Zn = new G4Element("Zinc", "Zn", 30., 65.38 * g/mole);
    G4Element* Cd = new G4Element("Cadmium", "Cd", 48., 112.41 * g/mole);  // 修正名称
    G4Element* Te = new G4Element("Tellurium", "Te", 52., 127.60 * g/mole); // 修正名称
    
    G4cout << "  Elements created: Zn, Cd, Te" << G4endl;
    
    // CdZnTe通常的比例：Cd(0.9), Zn(0.1), Te(1.0)
    // 注意：元素比例总和应为1，这里是原子比例
    G4Material* CdZnTe = new G4Material("CdZnTe", 6.0 * g/cm3, 3);  // 实际密度约5.85 g/cm³
    
    // 使用原子数比例 (Cd₀.₉Zn₀.₁Te)
    CdZnTe->AddElement(Cd, 0.45);  // Cd的比例
    CdZnTe->AddElement(Zn, 0.05);  // Zn的比例
    CdZnTe->AddElement(Te, 0.50);  // Te的比例，总和为1.0
    
    fdetector_mat = CdZnTe;
    
    G4cout << "  CdZnTe material created: density = " 
           << fdetector_mat->GetDensity()/(g/cm3) << " g/cm3" << G4endl;


    // Create YSO (Y2SiO5) scintillator material
    G4cout << "\nCreating YSO material..." << G4endl;
    
    G4Element* Y = new G4Element("Yttrium", "Y", 39., 88.9059 * g/mole);
    G4Element* Si = new G4Element("Silicon_YSO", "Si", 14., 28.086 * g/mole);
    G4Element* O = new G4Element("Oxygen", "O", 8., 15.9994 * g/mole);
    
    G4cout << "  Elements created: Y, Si, O" << G4endl;
    
    // Y2SiO5: 2 Yttrium, 1 Silicon, 5 Oxygen atoms; density ~4.45 g/cm3
    G4Material* YSO = new G4Material("YSO", 4.45 * g/cm3, 3);
    
    // Use number-of-atoms interface for correct stoichiometry
    YSO->AddElement(Y, 2);   // 2 Y atoms
    YSO->AddElement(Si, 1);  // 1 Si atom
    YSO->AddElement(O, 5);   // 5 O atoms
    
    fdetector_mat1 = YSO;
    
    G4cout << "  YSO material created: density = " 
           << fdetector_mat1->GetDensity()/(g/cm3) << " g/cm3" << G4endl;
    // 输出所有材料信息用于调试
    G4cout << "\n=== All Materials ===" << G4endl;
    const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
    G4cout << "Total materials: " << materialTable->size() << G4endl;
    
    for (size_t i = 0; i < materialTable->size(); ++i) {
        G4Material* mat = (*materialTable)[i];
        G4cout << "  [" << i << "] " << mat->GetName() 
               << " (density: " << mat->GetDensity()/(g/cm3) << " g/cm3)";
        
        // 如果是自定义材料，显示组成
        if (mat == fdetector_mat) {
            G4cout << " [CUSTOM - CdZnTe]";
        }
        G4cout << G4endl;
    }
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4cout << "\n=== Constructing Geometry ===" << G4endl;
    chamberComponents.clear();
    
    // 确保材料已定义
    if (!world_mat) {
        G4Exception("DetectorConstruction::Construct",
                   "Geometry001", FatalException,
                   "World material is not defined!");
    }
    
    G4cout << "Using world material: " << world_mat->GetName() << G4endl;
    G4cout << "World material pointer valid: " << (world_mat != nullptr) << G4endl;
    
    // World
    G4double world_sizeX = 100 * cm, world_sizeY = 100 * cm, world_sizeZ = 160 * cm;
    auto solid_world = new G4Box("WorldSolid", 
                                 0.5 * world_sizeX, 
                                 0.5 * world_sizeY, 
                                 0.5 * world_sizeZ);
    
    G4bool checkOverlaps = true;
    
    // 创建世界逻辑体积
    logic_world = new G4LogicalVolume(solid_world, world_mat, "WorldLogic");
    
    // 验证逻辑体积的材料
    if (!logic_world->GetMaterial()) {
        G4Exception("DetectorConstruction::Construct",
                   "Geometry002", FatalException,
                   "World logical volume has no material!");
    }
    
    G4cout << "World logical volume created with material: " 
           << logic_world->GetMaterial()->GetName() << G4endl;
    
    auto physical_world = new G4PVPlacement(
        0,                                // 无旋转
        G4ThreeVector(0,0,0),             // 位置在原点
        logic_world,                      // 逻辑体积
        "WorldPhys",                      // 物理体积名称
        0,                                // 无母体积（世界是根体积）
        false,                            // 无布尔操作
        0,                                // 复制号
        checkOverlaps                     // 检查重叠
    );

    // ========== Chamber 组件开始 ==========
    
    // Absorber (铅)
    G4double Absorber_sizeX = 60 * cm, Absorber_sizeY = 60 * cm, Absorber_sizeZ = 10 * cm;
    auto solid_Absorber = new G4Box("AbsorberSolid", 
                                    0.5 * Absorber_sizeX, 
                                    0.5 * Absorber_sizeY, 
                                    0.5 * Absorber_sizeZ);
    
    logic_Absorber = new G4LogicalVolume(solid_Absorber, Absorber_mat, "AbsorberLogic");
    
    // 验证材料
    if (!logic_Absorber->GetMaterial()) {
        G4cout << "WARNING: Absorber has no material!" << G4endl;
    } else {
        G4cout << "Absorber material: " << logic_Absorber->GetMaterial()->GetName() << G4endl;
    }
    
    chamberComponents.push_back(logic_Absorber);

    new G4PVPlacement(
        0, 
        G4ThreeVector(0,0,-26*cm), 
        logic_Absorber, 
        "AbsorberPhys", 
        logic_world, 
        false, 
        0, 
        checkOverlaps
    );

    // Framework Fe
    G4double framework_Fe_sizeXY = 80 * cm, framework_Fe_sizeZ = 1 * cm;
    auto solid_framework_Fe = new G4Box("FrameworkFeSolid", 
                                        0.5 * framework_Fe_sizeXY, 
                                        0.5 * framework_Fe_sizeXY, 
                                        0.5 * framework_Fe_sizeZ);
    
    logic_framework_Fe = new G4LogicalVolume(solid_framework_Fe, framework_Fe_mat, "FrameworkFeLogic");
    chamberComponents.push_back(logic_framework_Fe);

    new G4PVPlacement(
        0, 
        G4ThreeVector(0,0,-20.5*cm), 
        logic_framework_Fe, 
        "FrameworkFePhys",
        logic_world, 
        false, 
        0, 
        checkOverlaps
    );

    // Framework Al
    G4double framework_Al_sizeLong  = 60.0 * cm;
    G4double framework_Al_sizeShort = 40.0 * cm;
    G4double framework_Al_sizeZ     = 3.0 * mm;

    auto solid_framework_Al = new G4Box(
    "FrameworkAlSolid",
    0.5 * framework_Al_sizeLong,   // local X = long side
    0.5 * framework_Al_sizeShort,  // local Y = short side
    0.5 * framework_Al_sizeZ
);

logic_framework_Al = new G4LogicalVolume(
    solid_framework_Al,
    framework_Al_mat,
    "FrameworkAlLogic"
);
chamberComponents.push_back(logic_framework_Al);

std::vector<G4double> alPositions = {
    2.35 * cm,   // Al1
    32.35 * cm,  // Al2
    52.35 * cm,  // Al3
    69.85 * cm   // Al4
};

// true:  long side along world X
// false: long side along world Y (rotate local long-X axis by 90 deg around Z)
std::vector<G4bool> alLongAlongX = {
    true,   // Al1: X
    false,  // Al2: Y
    true,   // Al3: X
    true    // Al4: X, perpendicular to Al2
};

for (size_t i = 0; i < alPositions.size(); ++i) {
    G4RotationMatrix* rot = nullptr;

    if (!alLongAlongX[i]) {
        rot = new G4RotationMatrix();
        rot->rotateZ(90.0 * deg);
    }

    new G4PVPlacement(
        rot,
        G4ThreeVector(0, 0, alPositions[i]),
        logic_framework_Al,
        "FrameworkAlPhys" + std::to_string(i + 1),
        logic_world,
        false,
        static_cast<G4int>(i),
        checkOverlaps
    );
}

 const G4double scintLong       = 40.0 * cm;
const G4double scintWidth      = 4.0  * cm;
const G4double scintThickness  = 0.6  * cm;
const G4double scintHalfZ      = 0.5 * scintThickness;
const G4double alHalfThick     = 0.5 * framework_Al_sizeZ;
const G4double scintAlGap      = 2.0 * cm;

// Pb is at negative Z, so the Pb-facing scintillator plane is below each Al plate.
auto scintillatorCenterZ = [&](G4double alCenterZ) -> G4double {
    return alCenterZ - alHalfThick - scintAlGap - scintHalfZ;
};

// Two bar geometries are used so that the long side directly follows world X/Y.
auto solid_ScintX = new G4Box(
    "ScintillatorSolidX",
    0.5 * scintLong,
    0.5 * scintWidth,
    scintHalfZ
);

auto solid_ScintY = new G4Box(
    "ScintillatorSolidY",
    0.5 * scintWidth,
    0.5 * scintLong,
    scintHalfZ
);

auto logic_ScintX = new G4LogicalVolume(
    solid_ScintX,
    fChamberMaterial,
    "ScintillatorLogicX"
);

auto logic_ScintY = new G4LogicalVolume(
    solid_ScintY,
    fChamberMaterial,
    "ScintillatorLogicY"
);

chamberComponents.push_back(logic_ScintX);
chamberComponents.push_back(logic_ScintY);

// Station-level copy numbers:
//   0 -> scintillators in front of Al1
//   1 -> scintillators in front of Al2
//   2 -> scintillators in front of Al4
//
// All 9 bars of one station deliberately share the same copy number so that
// the current ChamberSD/Analysis station semantics can remain layer-based.
// Individual bars remain identifiable by their physical-volume names.
struct ScintStation {
    G4int alIndex;        // 0-based Al index
    G4bool longAlongX;
    G4int chamberCopyNo;
    G4String tag;
};

std::vector<ScintStation> scintStations = {
    {0, true,  0, "Al1"},  // Al1: bars long along X
    {1, false, 1, "Al2"},  // Al2: bars long along Y
    {3, true,  2, "Al4"}   // Al4: bars long along X, perpendicular to Al2
};

for (const auto& station : scintStations) {
    const G4double zScint = scintillatorCenterZ(alPositions[station.alIndex]);

    for (G4int ibar = 0; ibar < 9; ++ibar) {
        // -16, -12, ..., +16 cm
        const G4double transverse =
            (-16.0 + 4.0 * static_cast<G4double>(ibar)) * cm;

        G4ThreeVector pos;
        G4LogicalVolume* scintLV = nullptr;

        if (station.longAlongX) {
            // long side along X; distribute bars along Y
            pos = G4ThreeVector(0.0, transverse, zScint);
            scintLV = logic_ScintX;
        } else {
            // long side along Y; distribute bars along X
            pos = G4ThreeVector(transverse, 0.0, zScint);
            scintLV = logic_ScintY;
        }

        new G4PVPlacement(
            nullptr,
            pos,
            scintLV,
            "ScintillatorPhys_" + station.tag + "_" + std::to_string(ibar + 1),
            logic_world,
            false,
            station.chamberCopyNo,
            checkOverlaps
        );
    }
}

    // Alternating Si / YSO detector stack between the Al frames at 32.35 cm and 52.35 cm
    G4double detectorLayerSizeXY = 20.0 * cm;
    G4double ysoLayerThickness = 3.0 * mm;       // YSO 保持 3mm
    G4double siliconLayerThickness = 6.0 * mm;   // Si 加倍到 6mm
    G4double detectorLayerGap = 3.0 * cm;
    G4double ysoLayerHalfZ = 0.5 * ysoLayerThickness;
    G4double siliconLayerHalfZ = 0.5 * siliconLayerThickness;
    G4double upperAlFrameZ = 52.35 * cm;
    G4double upperAlFrameHalfZ = 0.5 * framework_Al_sizeZ;
    G4double layerPitch = detectorLayerGap + ysoLayerThickness;  // 使用YSO厚度计算间距

    auto solid_Silicon = new G4Box("SiliconLayerSolid",
                                   0.5 * detectorLayerSizeXY,
                                   0.5 * detectorLayerSizeXY,
                                   siliconLayerHalfZ);
    auto solid_YSO = new G4Box("YSOLayerSolid",
                               0.5 * detectorLayerSizeXY,
                               0.5 * detectorLayerSizeXY,
                               ysoLayerHalfZ);

    logic_Silicon = new G4LogicalVolume(solid_Silicon, Silicon_mat, "SiliconLayerLogic");
    logic_YSO = new G4LogicalVolume(solid_YSO, fdetector_mat1, "YSOLayerLogic");
    logic_CdZnTe = nullptr;

    chamberComponents.push_back(logic_Silicon);
    chamberComponents.push_back(logic_YSO);

    // 为保证 5 层结构完整放置在两层 Al 框架之间，底部 Si 与 52.35 cm 铝框架朝向内部的一侧贴合
    G4double firstSiliconCenterZ = upperAlFrameZ - upperAlFrameHalfZ - siliconLayerHalfZ;
    std::vector<G4double> siliconLayerPositions = {
        firstSiliconCenterZ,
        firstSiliconCenterZ - 2.0 * layerPitch,
        firstSiliconCenterZ - 4.0 * layerPitch
    };
    std::vector<G4double> ysoLayerPositions = {
        firstSiliconCenterZ - layerPitch,
        firstSiliconCenterZ - 3.0 * layerPitch
    };

    for (size_t i = 0; i < siliconLayerPositions.size(); ++i) {
        new G4PVPlacement(
            0,
            G4ThreeVector(0, 0, siliconLayerPositions[i]),
            logic_Silicon,
            "SiliconLayerPhys" + std::to_string(i + 1),
            logic_world,
            false,
            i,
            checkOverlaps
        );
    }

    for (size_t i = 0; i < ysoLayerPositions.size(); ++i) {
        new G4PVPlacement(
            0,
            G4ThreeVector(0, 0, ysoLayerPositions[i]),
            logic_YSO,
            "YSOLayerPhys" + std::to_string(i + 1),
            logic_world,
            false,
            i,
            checkOverlaps
        );
       }

    G4cout << "\nAlternating Si / YSO detector stack:" << G4endl;
    G4cout << "  Si layer thickness: " << siliconLayerThickness/mm << " mm" << G4endl;
    G4cout << "  YSO layer thickness: " << ysoLayerThickness/mm << " mm" << G4endl;
    G4cout << "  Layer XY size: " << detectorLayerSizeXY/cm << " cm x "
           << detectorLayerSizeXY/cm << " cm" << G4endl;
    G4cout << "  Gap between adjacent layers: " << detectorLayerGap/cm << " cm" << G4endl;
    for (size_t i = 0; i < siliconLayerPositions.size(); ++i) {
        G4cout << "  Si layer " << (i + 1) << " center z = "
               << siliconLayerPositions[i]/cm << " cm" << G4endl;
    }
    for (size_t i = 0; i < ysoLayerPositions.size(); ++i) {
        G4cout << "  YSO layer " << (i + 1) << " center z = "
               << ysoLayerPositions[i]/cm << " cm" << G4endl;
    }

    // 设置可视化属性
    SetVisualAttributes();

    G4cout << "\n=== Geometry Construction Complete ===" << G4endl;
    G4cout << "Total chamber components: " << chamberComponents.size() << G4endl;
    G4cout << "World size: " << world_sizeX/cm << "x" << world_sizeY/cm 
           << "x" << world_sizeZ/cm << " cm" << G4endl;
    
    // 验证所有组件都有材料
    G4cout << "\n=== Material Validation ===" << G4endl;
    for (size_t i = 0; i < chamberComponents.size(); ++i) {
        G4LogicalVolume* lv = chamberComponents[i];
        if (lv && lv->GetMaterial()) {
            G4cout << "✓ " << lv->GetName() << ": " << lv->GetMaterial()->GetName() << G4endl;
        } else {
            G4cout << "✗ " << (lv ? lv->GetName() : "NULL") << ": NO MATERIAL!" << G4endl;
        }
    }
    
    return physical_world;
}

void DetectorConstruction::SetVisualAttributes()
{
    // 为不同组件设置不同的颜色
    
    // 世界 - 透明
    G4VisAttributes* worldVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.1)); // 半透明白色
    worldVisAtt->SetVisibility(true);
    worldVisAtt->SetForceWireframe(true);
    logic_world->SetVisAttributes(worldVisAtt);
    
    // 铅吸收体 - 灰色
    if (logic_Absorber) {
        G4VisAttributes* absorberVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)); // 灰色
        absorberVisAtt->SetVisibility(true);
        absorberVisAtt->SetForceSolid(true);
        logic_Absorber->SetVisAttributes(absorberVisAtt);
    }
    
    // 不锈钢框架 - 银色
    if (logic_framework_Fe) {
        G4VisAttributes* feVisAtt = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8)); // 银色
        feVisAtt->SetVisibility(true);
        feVisAtt->SetForceSolid(true);
        logic_framework_Fe->SetVisAttributes(feVisAtt);
    }
    
    // 铝框架 - 浅蓝色
    if (logic_framework_Al) {
        G4VisAttributes* alVisAtt = new G4VisAttributes(G4Colour(0.7, 0.7, 1.0)); // 浅蓝色
        alVisAtt->SetVisibility(true);
        alVisAtt->SetForceSolid(true);
        logic_framework_Al->SetVisAttributes(alVisAtt);
    }
    
    // 闪烁体 - 绿色
    for (size_t i = 0; i < chamberComponents.size(); ++i) {
        G4LogicalVolume* lv = chamberComponents[i];
        if (lv && lv->GetName().find("Scintillator") != std::string::npos) {
            G4VisAttributes* scintVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.3)); // 半透明绿色
            scintVisAtt->SetVisibility(true);
            scintVisAtt->SetForceSolid(true);
            lv->SetVisAttributes(scintVisAtt);
        }
    }
    
    // 硅探测器 - 红色
    if (logic_Silicon) {
        G4VisAttributes* siVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.5)); // 半透明红色
        siVisAtt->SetVisibility(true);
        siVisAtt->SetForceSolid(true);
        logic_Silicon->SetVisAttributes(siVisAtt);
    }

    // YSO 探测器 - 蓝绿色
    if (logic_YSO) {
        G4VisAttributes* ysoVisAtt = new G4VisAttributes(G4Colour(0.0, 0.8, 0.8, 0.5));
        ysoVisAtt->SetVisibility(true);
        ysoVisAtt->SetForceSolid(true);
        logic_YSO->SetVisAttributes(ysoVisAtt);
    }
    
    // CdZnTe探测器 - 黄色
    if (logic_CdZnTe) {
        G4VisAttributes* cdznteVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0)); // 黄色
        cdznteVisAtt->SetVisibility(true);
        cdznteVisAtt->SetForceSolid(true);
        logic_CdZnTe->SetVisAttributes(cdznteVisAtt);
    }
}
/*
void DetectorConstruction::ConstructSDandField()
{
    auto sdManager = G4SDManager::GetSDMpointer();
    
    // 创建统一的 Chamber 敏感探测器
    G4String chamberSDname = "/Chamber";
    auto chamberSD = new ChamberSD(chamberSDname, "ChamberHitsCollection");
    sdManager->AddNewDetector(chamberSD);

    G4cout << "\n=== Setting up Sensitive Detectors ===" << G4endl;
    G4cout << "Chamber SD name: " << chamberSDname << G4endl;
    
    // 为所有 chamber 组件设置相同的敏感探测器
    for (size_t i = 0; i < chamberComponents.size(); ++i) {
        if (chamberComponents[i]) {
            chamberComponents[i]->SetSensitiveDetector(chamberSD);
            G4cout << "  Set SD for: " << chamberComponents[i]->GetName() 
                   << " (material: " << chamberComponents[i]->GetMaterial()->GetName() << ")" << G4endl;
        }
    }
    
    G4cout << "Total components with sensitive detector: " << chamberComponents.size() << G4endl;
    G4cout << "=== SD Setup Complete ===" << G4endl;
}
*/


void DetectorConstruction::ConstructSDandField()
{
    auto sdManager = G4SDManager::GetSDMpointer();
    G4String chamberSDname = "/Chamber";
    auto chamberSD = new ChamberSD(chamberSDname, "ChamberHitsCollection");
    sdManager->AddNewDetector(chamberSD);
    for (auto lv : chamberComponents) {
        if (lv && (lv == logic_Absorber ||
                   lv->GetName().find("Scintillator") != std::string::npos ||
                   lv == logic_YSO)) {   // <--- 添加这一条件
            lv->SetSensitiveDetector(chamberSD);
        }
    }
    // 为 Absorber（铅砖）和所有 Scintillator 设置 SD
    //for (auto lv : chamberComponents) {
        //if (lv && (lv == logic_Absorber ||
                   //lv->GetName().find("Scintillator") != std::string::npos)) {
            //lv->SetSensitiveDetector(chamberSD);
        //}
    //}
}
void DetectorConstruction::SetMaxStep(G4double maxStep)
{
    if ((fStepLimit) && (maxStep > 0.)) {
        fStepLimit->SetMaxAllowedStep(maxStep);
    }
}

void DetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
    fCheckOverlaps = checkOverlaps;
}

}