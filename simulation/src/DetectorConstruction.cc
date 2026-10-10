#include "DetectorConstruction.hh"
#include "TOFmeterDetectorSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4ExtrudedSolid.hh"
#include "G4TwoVector.hh"
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
#include "G4Transform3D.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"
#include "G4UIparameter.hh"
#include "G4Point3D.hh"

#include <cfloat>
#include <functional>
#include <iomanip>
#include <sstream>

namespace
{
// /tof/blocks/add needs nine parameters; G4GenericMessenger would pass a method
// only the first word of the line, so it gets its own messenger
class BlockMessenger : public G4UImessenger
{
public:
    explicit BlockMessenger(Pb::DetectorConstruction* detector) : fDetector(detector)
    {
        fDirectory = new G4UIdirectory("/tof/blocks/");
        fDirectory->SetGuidance("Passive boxes such as shielding bricks; changes take effect with /update");

        fAddCmd = new G4UIcommand("/tof/blocks/add", this);
        fAddCmd->SetGuidance("Add a box: name material x y z dx dy dz [unit] (centre and full sizes)");
        for (const char* name : {"name", "material"}) {
            fAddCmd->SetParameter(new G4UIparameter(name, 's', false));
        }
        for (const char* name : {"x", "y", "z", "dx", "dy", "dz"}) {
            fAddCmd->SetParameter(new G4UIparameter(name, 'd', false));
        }
        auto unit = new G4UIparameter("unit", 's', true);
        unit->SetDefaultValue("cm");
        fAddCmd->SetParameter(unit);
        fAddCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

        fClearCmd = new G4UIcmdWithoutParameter("/tof/blocks/clear", this);
        fClearCmd->SetGuidance("Remove all boxes added with /tof/blocks/add");
        fClearCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    }

    ~BlockMessenger() override
    {
        delete fAddCmd;
        delete fClearCmd;
        delete fDirectory;
    }

    void SetNewValue(G4UIcommand* command, G4String value) override
    {
        if (command == fAddCmd) fDetector->AddBlock(value);
        else if (command == fClearCmd) fDetector->ClearBlocks();
    }

private:
    Pb::DetectorConstruction* fDetector;
    G4UIdirectory* fDirectory;
    G4UIcommand* fAddCmd;
    G4UIcmdWithoutParameter* fClearCmd;
};
}

namespace Pb
{

// 构造函数
DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction(),
      fWorldPhys(nullptr),
      logic_world(nullptr),
      logic_Absorber(nullptr),
      logic_framework_Fe(nullptr),
      logic_framework_Al(nullptr),
            logic_YSO(nullptr),
      logic_Silicon(nullptr),
      logic_CdZnTe(nullptr),
      logic_CameraYSO(nullptr),
      logic_CameraLYSO(nullptr),
      world_mat(nullptr),
      Absorber_mat(nullptr),
      framework_Fe_mat(nullptr),
      framework_Al_mat(nullptr),
      fChamberMaterial(nullptr),
      Silicon_mat(nullptr),
      fdetector_mat(nullptr),
            fdetector_mat1(nullptr),
      fLYSO_mat(nullptr),
      fNbOfChambers(0),
      fTargetMaterial(nullptr),
      fStepLimit(nullptr),
      fCheckOverlaps(true),
      fBuildStations(true),
      fBuildFrame(true),
      fBuildStack(true),
      fBuildAbsorber(true),
      fBuildLeadWalls(true),
      fLeadThickness(10.0 * cm),
      fWorldMaterialName("G4_AIR"),
      fWorldSize(4.0 * m, 4.0 * m, 3.0 * m),
      fBuildCamera(false),
      fCameraPosition(0., 0., 0.),
      fCameraRotationX(0.),
      fCameraSpacing21(25.0 * mm),
      fCameraSpacing10(35.0 * mm),
      fBuildTarget(false),
      fTargetMaterialName("G4_Fe"),
      fTargetSize(10.0 * cm, 10.0 * cm, 2.0 * cm),
      fTargetPosition(0., 0., 0.),
      fBuildDegrader(false),
      fDegraderMaterialName("G4_Pb"),
      fDegraderSize(40.0 * cm, 40.0 * cm, 5.0 * cm),
      fDegraderPosition(0., 0., 0.)
{
    // 在构造函数中调用DefineMaterials确保材料被创建
    DefineMaterials();

    // 几何参数的宏命令 (/tof/det, /tof/camera, /tof/target, /tof/degrader, /tof/blocks)
    DefineCommands();
}

// 析构函数
DetectorConstruction::~DetectorConstruction() = default;

void DetectorConstruction::DefineMaterials()
{
    G4NistManager* nist = G4NistManager::Instance();

    G4cout << "\n=== Defining Materials ===" << G4endl;

    // 世界材料: 默认空气 (/tof/det/worldMaterial 可改为 G4_Galactic)
    world_mat = nist->FindOrBuildMaterial("G4_AIR");
    nist->FindOrBuildMaterial("G4_Galactic");
    if (!world_mat) {
        G4Exception("DetectorConstruction::DefineMaterials",
                   "Material001", FatalException,
                   "Failed to create G4_AIR material for world!");
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
    // Note: AddElement(G4Element*, G4double) takes MASS fractions, so these
    // numbers are not the atom ratios of Cd0.9Zn0.1Te (material currently unused)
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

    // LYSO (Lu1.8Y0.2SiO5) for the CH0 absorber of the Compton camera;
    // density from Sun et al., arXiv:2608.21216 (~7.2 g/cm3)
    G4Element* Lu = nist->FindOrBuildElement("Lu");
    G4Material* LYSO = new G4Material("LYSO", 7.2 * g/cm3, 4);
    LYSO->AddElement(Lu, 18);
    LYSO->AddElement(Y, 2);
    LYSO->AddElement(Si, 10);
    LYSO->AddElement(O, 50);
    fLYSO_mat = LYSO;

    G4cout << "  LYSO material created: density = "
           << fLYSO_mat->GetDensity()/(g/cm3) << " g/cm3" << G4endl;
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

G4Material* DetectorConstruction::FindMaterial(const G4String& name) const
{
    // custom materials (YSO, LYSO, CdZnTe) first, then the NIST database
    G4Material* mat = G4Material::GetMaterial(name, false);
    if (!mat) mat = G4NistManager::Instance()->FindOrBuildMaterial(name);
    if (!mat) {
        G4ExceptionDescription ed;
        ed << "Unknown material \"" << name << "\" (use a NIST name such as G4_Fe, or YSO/LYSO/CdZnTe)";
        G4Exception("DetectorConstruction::FindMaterial", "Material002", FatalErrorInArgument, ed);
    }
    return mat;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4cout << "\n=== Constructing Geometry ===" << G4endl;
    chamberComponents.clear();
    fKinds.clear();

    // Volumes of a previous build were deleted by /update; parts that are
    // switched off must not leave dangling pointers behind
    logic_Absorber = nullptr;
    logic_framework_Fe = nullptr;
    logic_framework_Al = nullptr;
    logic_YSO = nullptr;
    logic_Silicon = nullptr;
    logic_CdZnTe = nullptr;
    logic_CameraYSO = nullptr;
    logic_CameraLYSO = nullptr;

    world_mat = FindMaterial(fWorldMaterialName);

    // 确保材料已定义
    if (!world_mat) {
        G4Exception("DetectorConstruction::Construct",
                   "Geometry001", FatalException,
                   "World material is not defined!");
    }

    G4cout << "Using world material: " << world_mat->GetName() << G4endl;
    G4cout << "World material pointer valid: " << (world_mat != nullptr) << G4endl;

    // World (large enough to hold the cosmic-muon source plane, /tof/det/worldSize)
    G4double world_sizeX = fWorldSize.x(), world_sizeY = fWorldSize.y(), world_sizeZ = fWorldSize.z();
    auto solid_world = new G4Box("WorldSolid",
                                 0.5 * world_sizeX,
                                 0.5 * world_sizeY,
                                 0.5 * world_sizeZ);

    G4bool checkOverlaps = fCheckOverlaps;

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

    // Copy number -1 (as in the original framework): in the ASCII output,
    // copyNb1 = -1 marks volumes placed directly in the world
    auto physical_world = new G4PVPlacement(
        0,                                // 无旋转
        G4ThreeVector(0,0,0),             // 位置在原点
        logic_world,                      // 逻辑体积
        "WorldPhys",                      // 物理体积名称
        0,                                // 无母体积（世界是根体积）
        false,                            // 无布尔操作
        -1,                               // 复制号
        checkOverlaps                     // 检查重叠
    );

    // ========== Chamber 组件开始 ==========

    // z points up: muons come from +z. From the top: lead layer, steel plate,
    // top TOF on Al1, middle TOF on Al2, Si/YSO stack (Compton detector) on Al3,
    // bottom TOF on Al4; four lead walls around everything.

    // Framework Fe: stainless-steel plate holding the lead layer
    G4double framework_Fe_sizeXY = 80 * cm, framework_Fe_sizeZ = 1 * cm;
    G4double framework_Fe_z = 20.5 * cm;
    const G4double frameTopZ = framework_Fe_z + 0.5 * framework_Fe_sizeZ;     // top of the steel plate

    // Stainless-steel L beams on the four corners of the steel plate: 10 cm legs, 1 cm thick
    const G4double beamLeg = 10.0 * cm, beamThick = 1.0 * cm;

    // Lead walls: 10 cm thick, 5 cm outside the edges of the steel plate
    const G4double leadWallThick = 10.0 * cm, leadWallGap = 5.0 * cm;
    const G4double leadWallInner = 0.5 * framework_Fe_sizeXY + leadWallGap;

    // Lead layer on top of the steel plate (/tof/det/leadThickness, 10 cm by default):
    // its bottom stays on the plate, so a thicker layer grows upwards
    if (fLeadThickness < 0.) {
        G4Exception("DetectorConstruction::Construct", "Geometry004", FatalErrorInArgument,
                    "/tof/det/leadThickness must not be negative");
    }
    const G4double leadLayerThick = fLeadThickness;
    const G4double leadLayerZ = frameTopZ + 0.5 * leadLayerThick;

    if (fBuildAbsorber && leadLayerThick > 0.) {
    // Absorber (铅), the lead layer on top. It reaches the inner faces of the lead
    // walls, whose tops are level with it, so the lead covers the whole top.
    G4double Absorber_sizeX = 2 * leadWallInner, Absorber_sizeY = 2 * leadWallInner, Absorber_sizeZ = leadLayerThick;
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
    fKinds[logic_Absorber] = DetectorKind::Absorber;

    new G4PVPlacement(
        0,
        G4ThreeVector(0,0,leadLayerZ),
        logic_Absorber,
        "AbsorberPhys",
        logic_world,
        false,
        0,
        checkOverlaps
    );
    }

    // Framework Al: the plates sit between the corner L beams. The long side spans the
    // distance between the inner faces of the beams (80 cm), the short side fits the
    // gap between the leg tips of two beams (80 + 2 x 1 - 2 x 10 = 62 cm).
    G4double framework_Al_sizeLong  = framework_Fe_sizeXY;
    G4double framework_Al_sizeShort = framework_Fe_sizeXY + 2 * beamThick - 2 * beamLeg;
    G4double framework_Al_sizeZ     = 3.0 * mm;

std::vector<G4double> alPositions = {
    -2.35 * cm,   // Al1: under the top TOF
    -32.35 * cm,  // Al2: under the middle TOF
    -52.35 * cm,  // Al3: under the Si/YSO stack
    -69.85 * cm   // Al4: under the bottom TOF
};

// true:  long side along world X
// false: long side along world Y (rotate local long-X axis by 90 deg around Z)
std::vector<G4bool> alLongAlongX = {
    true,   // Al1: X
    false,  // Al2: Y
    true,   // Al3: X
    true    // Al4: X, perpendicular to Al2
};

// Vertical extent of the frame and of the whole setup
const G4double setupBottomZ = alPositions[3] - 0.5 * framework_Al_sizeZ;   // bottom of Al4
const G4double leadTopZ = frameTopZ + leadLayerThick;                       // top of the lead layer

    if (fBuildFrame) {
    // Framework Fe, the stainless-steel plate under the lead layer
    auto solid_framework_Fe = new G4Box("FrameworkFeSolid",
                                        0.5 * framework_Fe_sizeXY,
                                        0.5 * framework_Fe_sizeXY,
                                        0.5 * framework_Fe_sizeZ);

    logic_framework_Fe = new G4LogicalVolume(solid_framework_Fe, framework_Fe_mat, "FrameworkFeLogic");
    chamberComponents.push_back(logic_framework_Fe);

    new G4PVPlacement(
        0,
        G4ThreeVector(0,0,framework_Fe_z),
        logic_framework_Fe,
        "FrameworkFePhys",
        logic_world,
        false,
        0,
        checkOverlaps
    );

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

    // Stainless-steel L beams (10 cm legs, 1 cm thick) on the four corners, from the
    // bottom of Al4 to the top of the steel plate. Each wraps a plate corner from
    // outside; the L is drawn for the (+x, +y) corner with its outer edge at the origin.
    const std::vector<G4TwoVector> beamProfile = {
        {0., 0.}, {0., -beamLeg}, {-beamThick, -beamLeg},
        {-beamThick, -beamThick}, {-beamLeg, -beamThick}, {-beamLeg, 0.}};
    auto solid_beam = new G4ExtrudedSolid("SupportBeamSolid", beamProfile,
                                          0.5 * (frameTopZ - setupBottomZ),
                                          G4TwoVector(), 1., G4TwoVector(), 1.);
    auto logic_beam = new G4LogicalVolume(solid_beam, framework_Fe_mat, "SupportBeamLogic");
    logic_beam->SetVisAttributes(G4VisAttributes(G4Colour(0.75, 0.75, 0.75)));
    chamberComponents.push_back(logic_beam);
    const G4double beamCorner = 0.5 * framework_Fe_sizeXY + beamThick;
    for (G4int k = 0; k < 4; ++k) {
        G4RotationMatrix rotation;
        rotation.rotateZ(90. * deg * k);
        const G4ThreeVector corner =
            rotation * G4ThreeVector(beamCorner, beamCorner, 0.5 * (frameTopZ + setupBottomZ));
        new G4PVPlacement(G4Transform3D(rotation, corner), logic_beam, "SupportBeamPhys",
                          logic_world, false, k, checkOverlaps);
    }
    }

    if (fBuildLeadWalls) {
    // Four 10 cm lead walls 5 cm outside the steel-plate edges, from the bottom of Al4
    // to the top of the lead layer (to the top of the steel plate when the layer is
    // 0 cm thick). The walls along x span the full outer width, the walls along y fit
    // between them.
    const G4double wallThick = leadWallThick;
    const G4double wallInner = leadWallInner;
    const G4double wallCentre = wallInner + 0.5 * wallThick;
    const G4double wallHalfZ = 0.5 * (leadTopZ - setupBottomZ);
    const G4double wallZ = 0.5 * (leadTopZ + setupBottomZ);
    auto solid_wallX = new G4Box("LeadWallXSolid", wallInner + wallThick, 0.5 * wallThick, wallHalfZ);
    auto solid_wallY = new G4Box("LeadWallYSolid", 0.5 * wallThick, wallInner, wallHalfZ);
    auto logic_wallX = new G4LogicalVolume(solid_wallX, Absorber_mat, "LeadWallXLogic");
    auto logic_wallY = new G4LogicalVolume(solid_wallY, Absorber_mat, "LeadWallYLogic");
    for (auto lv : {logic_wallX, logic_wallY}) {
        lv->SetVisAttributes(G4VisAttributes(G4Colour(0.4, 0.4, 0.45, 0.3)));
        chamberComponents.push_back(lv);
    }
    new G4PVPlacement(nullptr, G4ThreeVector(0., wallCentre, wallZ), logic_wallX, "LeadWallPosY",
                      logic_world, false, 0, checkOverlaps);
    new G4PVPlacement(nullptr, G4ThreeVector(0., -wallCentre, wallZ), logic_wallX, "LeadWallNegY",
                      logic_world, false, 1, checkOverlaps);
    new G4PVPlacement(nullptr, G4ThreeVector(wallCentre, 0., wallZ), logic_wallY, "LeadWallPosX",
                      logic_world, false, 2, checkOverlaps);
    new G4PVPlacement(nullptr, G4ThreeVector(-wallCentre, 0., wallZ), logic_wallY, "LeadWallNegX",
                      logic_world, false, 3, checkOverlaps);
    }

    if (fBuildStations) {
 const G4double scintLong       = 40.0 * cm;
const G4double scintWidth      = 4.0  * cm;
const G4double scintThickness  = 0.6  * cm;
const G4double scintHalfZ      = 0.5 * scintThickness;
const G4double alHalfThick     = 0.5 * framework_Al_sizeZ;
const G4double scintAlGap      = 2.0 * cm;
const G4int    nBarsPerStation = 9;

// Pb is on top (+z), so the Pb-facing scintillator plane is above each Al plate.
auto scintillatorCenterZ = [&](G4double alCenterZ) -> G4double {
    return alCenterZ + alHalfThick + scintAlGap + scintHalfZ;
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
fKinds[logic_ScintX] = DetectorKind::Bar;
fKinds[logic_ScintY] = DetectorKind::Bar;

// Station-level copy numbers:
//   0 -> top TOF, on Al1
//   1 -> middle TOF, on Al2
//   2 -> bottom TOF, on Al4
//
// The 9 bars of a station are placed in an air envelope that carries the
// station copy number, and each bar carries its own index 0-8. In the ASCII
// output this gives copyNb1 = station and copyNb = bar.
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

const G4double stationWidth = nBarsPerStation * scintWidth;

for (const auto& station : scintStations) {
    const G4double zScint = scintillatorCenterZ(alPositions[station.alIndex]);

    auto solid_Station = station.longAlongX
        ? new G4Box("StationSolid_" + station.tag, 0.5 * scintLong, 0.5 * stationWidth, scintHalfZ)
        : new G4Box("StationSolid_" + station.tag, 0.5 * stationWidth, 0.5 * scintLong, scintHalfZ);
    auto logic_Station = new G4LogicalVolume(solid_Station, world_mat, "StationLogic_" + station.tag);
    logic_Station->SetVisAttributes(G4VisAttributes::GetInvisible());

    for (G4int ibar = 0; ibar < nBarsPerStation; ++ibar) {
        // -16, -12, ..., +16 cm
        const G4double transverse =
            (-16.0 + 4.0 * static_cast<G4double>(ibar)) * cm;

        G4ThreeVector pos;
        G4LogicalVolume* scintLV = nullptr;

        if (station.longAlongX) {
            // long side along X; distribute bars along Y
            pos = G4ThreeVector(0.0, transverse, 0.0);
            scintLV = logic_ScintX;
        } else {
            // long side along Y; distribute bars along X
            pos = G4ThreeVector(transverse, 0.0, 0.0);
            scintLV = logic_ScintY;
        }

        new G4PVPlacement(
            nullptr,
            pos,
            scintLV,
            "ScintillatorPhys_" + station.tag + "_" + std::to_string(ibar + 1),
            logic_Station,
            false,
            ibar,
            checkOverlaps
        );
    }

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0.0, 0.0, zScint),
        logic_Station,
        "StationPhys_" + station.tag,
        logic_world,
        false,
        station.chamberCopyNo,
        checkOverlaps
    );
}
    }

    if (fBuildStack) {
    // Alternating Si / YSO detector stack (Compton detector) on the Al3 plate, below the middle TOF
    G4double detectorLayerSizeXY = 20.0 * cm;
    G4double ysoLayerThickness = 3.0 * mm;       // YSO 保持 3mm
    G4double siliconLayerThickness = 6.0 * mm;   // Si 加倍到 6mm
    G4double detectorLayerGap = 3.0 * cm;
    G4double ysoLayerHalfZ = 0.5 * ysoLayerThickness;
    G4double siliconLayerHalfZ = 0.5 * siliconLayerThickness;
    G4double stackAlFrameZ = alPositions[2];   // Al3
    G4double stackAlFrameHalfZ = 0.5 * framework_Al_sizeZ;
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
    fKinds[logic_YSO] = DetectorKind::StackYSO;
    fKinds[logic_Silicon] = DetectorKind::StackSi;

    // 底部 Si 贴合在 Al3 铝框架 (z = -52.35 cm) 的上表面, 其余各层向上排列
    G4double firstSiliconCenterZ = stackAlFrameZ + stackAlFrameHalfZ + siliconLayerHalfZ;
    std::vector<G4double> siliconLayerPositions = {
        firstSiliconCenterZ,
        firstSiliconCenterZ + 2.0 * layerPitch,
        firstSiliconCenterZ + 4.0 * layerPitch
    };
    std::vector<G4double> ysoLayerPositions = {
        firstSiliconCenterZ + layerPitch,
        firstSiliconCenterZ + 3.0 * layerPitch
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
    // layerPitch only adds the YSO thickness, so the face-to-face gap between a
    // 6 mm Si and a 3 mm YSO layer is smaller than detectorLayerGap
    G4cout << "  Gap between adjacent layers (face to face): "
           << (layerPitch - siliconLayerHalfZ - ysoLayerHalfZ)/cm << " cm" << G4endl;
    for (size_t i = 0; i < siliconLayerPositions.size(); ++i) {
        G4cout << "  Si layer " << (i + 1) << " center z = "
               << siliconLayerPositions[i]/cm << " cm" << G4endl;
    }
    for (size_t i = 0; i < ysoLayerPositions.size(); ++i) {
        G4cout << "  YSO layer " << (i + 1) << " center z = "
               << ysoLayerPositions[i]/cm << " cm" << G4endl;
    }
    }

    if (fBuildCamera) ConstructComptonCamera();
    ConstructTargetAndDegrader();
    ConstructBlocks();

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

    fWorldPhys = physical_world;
    return physical_world;
}

void DetectorConstruction::PrintGeometry()
{
    if (!fWorldPhys) {
        G4cout << "No geometry has been built yet" << G4endl;
        return;
    }
    // the table sets a fixed 2-digit format; restore the stream afterwards
    const std::ios_base::fmtflags oldFlags = G4cout.flags();
    const std::streamsize oldPrecision = G4cout.precision();
    G4cout << "\n=== Placed volumes (world frame, cm; size is the axis-aligned extent) ===\n"
           << "name:copy                       material                    centre x       y       z"
           << "     size x       y       z      z from      to" << G4endl;

    std::function<void(const G4VPhysicalVolume*, const G4Transform3D&, G4int)> print =
        [&](const G4VPhysicalVolume* pv, const G4Transform3D& toWorld, G4int depth) {
        G4ThreeVector lo, hi;
        pv->GetLogicalVolume()->GetSolid()->BoundingLimits(lo, hi);
        G4ThreeVector wMin(DBL_MAX, DBL_MAX, DBL_MAX), wMax(-DBL_MAX, -DBL_MAX, -DBL_MAX);
        for (G4int i = 0; i < 8; ++i) {
            const G4Point3D corner = toWorld * G4Point3D((i & 1) ? hi.x() : lo.x(),
                                                         (i & 2) ? hi.y() : lo.y(),
                                                         (i & 4) ? hi.z() : lo.z());
            wMin = G4ThreeVector(std::min(wMin.x(), corner.x()), std::min(wMin.y(), corner.y()),
                                 std::min(wMin.z(), corner.z()));
            wMax = G4ThreeVector(std::max(wMax.x(), corner.x()), std::max(wMax.y(), corner.y()),
                                 std::max(wMax.z(), corner.z()));
        }
        const G4ThreeVector centre = 0.5 * (wMin + wMax) / cm;
        const G4ThreeVector size = (wMax - wMin) / cm;
        const G4String name = G4String(2 * depth, ' ') + pv->GetName() + ":" + std::to_string(pv->GetCopyNo());
        G4cout << std::left << std::setw(32) << name << std::setw(26)
               << pv->GetLogicalVolume()->GetMaterial()->GetName() << std::right << std::fixed
               << std::setprecision(2) << std::setw(10) << centre.x() << std::setw(8) << centre.y()
               << std::setw(8) << centre.z() << std::setw(11) << size.x() << std::setw(8) << size.y()
               << std::setw(8) << size.z() << std::setw(12) << wMin.z() / cm << std::setw(8)
               << wMax.z() / cm << std::defaultfloat << G4endl;

        const G4LogicalVolume* lv = pv->GetLogicalVolume();
        for (size_t i = 0; i < lv->GetNoDaughters(); ++i) {
            const G4VPhysicalVolume* daughter = lv->GetDaughter(i);
            print(daughter,
                  toWorld * G4Transform3D(daughter->GetObjectRotationValue(), daughter->GetObjectTranslation()),
                  depth + 1);
        }
    };
    print(fWorldPhys, G4Transform3D(), 0);
    G4cout.flags(oldFlags);
    G4cout.precision(oldPrecision);
}

void DetectorConstruction::ConstructComptonCamera()
{
    // Three-layer Compton camera of Sun et al., arXiv:2608.21216, Sec. 2.1:
    // CH2 and CH1 are 3 mm YSO scatterers, CH0 is a 6 mm LYSO absorber, all
    // 47.04 x 47.04 mm. The paper gives the 25 mm (CH2-CH1) and 35 mm (CH1-CH0)
    // separations without saying whether they are centre to centre or face to
    // face; centre to centre is assumed (it matches the layer z in its Fig. 5a).
    // Pixels, BaSO4 wrapping, SiPMs, electronics and the Al cover are not modelled.
    const G4double side = 47.04 * mm;
    const G4double scatThickness = 3.0 * mm;
    const G4double absThickness = 6.0 * mm;
    const G4double margin = 1.0 * mm;

    // Layer centres in the module frame: CH2 at the origin, module axis along +z
    const G4double zCH2 = 0.;
    const G4double zCH1 = zCH2 - fCameraSpacing21;
    const G4double zCH0 = zCH1 - fCameraSpacing10;
    const G4double zTop = zCH2 + 0.5 * scatThickness + margin;
    const G4double zBottom = zCH0 - 0.5 * absThickness - margin;
    const G4double zEnvelope = 0.5 * (zTop + zBottom);

    auto solid_CameraEnvelope = new G4Box("CameraEnvelopeSolid",
                                          0.5 * side + margin,
                                          0.5 * side + margin,
                                          0.5 * (zTop - zBottom));
    auto logic_CameraEnvelope = new G4LogicalVolume(solid_CameraEnvelope, world_mat, "CameraEnvelopeLogic");
    logic_CameraEnvelope->SetVisAttributes(G4VisAttributes::GetInvisible());

    auto solid_CameraYSO = new G4Box("CameraYSOSolid", 0.5 * side, 0.5 * side, 0.5 * scatThickness);
    auto solid_CameraLYSO = new G4Box("CameraLYSOSolid", 0.5 * side, 0.5 * side, 0.5 * absThickness);
    logic_CameraYSO = new G4LogicalVolume(solid_CameraYSO, fdetector_mat1, "CameraYSOLogic");
    logic_CameraLYSO = new G4LogicalVolume(solid_CameraLYSO, fLYSO_mat, "CameraLYSOLogic");
    chamberComponents.push_back(logic_CameraYSO);
    chamberComponents.push_back(logic_CameraLYSO);
    fKinds[logic_CameraYSO] = DetectorKind::CameraLayer;
    fKinds[logic_CameraLYSO] = DetectorKind::CameraLayer;

    // Copy number = CH index of the paper
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zCH2 - zEnvelope), logic_CameraYSO,
                      "CameraCH2", logic_CameraEnvelope, false, 2, fCheckOverlaps);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zCH1 - zEnvelope), logic_CameraYSO,
                      "CameraCH1", logic_CameraEnvelope, false, 1, fCheckOverlaps);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zCH0 - zEnvelope), logic_CameraLYSO,
                      "CameraCH0", logic_CameraEnvelope, false, 0, fCheckOverlaps);

    // Rotate about x around the centre of CH2, then put CH2 at fCameraPosition.
    // Envelope copy number 10: camera hits have copyNb1 = 10 in the ASCII output.
    G4RotationMatrix rotation;
    rotation.rotateX(fCameraRotationX);
    const G4ThreeVector envelopeCentre = fCameraPosition + rotation * G4ThreeVector(0, 0, zEnvelope);
    new G4PVPlacement(G4Transform3D(rotation, envelopeCentre), logic_CameraEnvelope,
                      "CameraEnvelopePhys", logic_world, false, 10, fCheckOverlaps);

    G4cout << "\nCompton camera (CH2/CH1 YSO 3 mm, CH0 LYSO 6 mm, 47.04 mm square):" << G4endl;
    G4cout << "  CH2 centre: " << fCameraPosition/cm << " cm, rotation about x: "
           << fCameraRotationX/deg << " deg" << G4endl;
    G4cout << "  spacing CH2-CH1 / CH1-CH0 (centre to centre): " << fCameraSpacing21/mm
           << " / " << fCameraSpacing10/mm << " mm" << G4endl;
}

void DetectorConstruction::ConstructTargetAndDegrader()
{
    // 未知材料 (target) and 降速铅砖 (degrader); passive, off by default
    if (fBuildTarget) {
        auto solid_Target = new G4Box("TargetSolid", 0.5 * fTargetSize.x(), 0.5 * fTargetSize.y(), 0.5 * fTargetSize.z());
        auto logic_Target = new G4LogicalVolume(solid_Target, FindMaterial(fTargetMaterialName), "TargetLogic");
        logic_Target->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0.6, 0.0)));
        new G4PVPlacement(nullptr, fTargetPosition, logic_Target, "Target", logic_world, false, 0, fCheckOverlaps);
        G4cout << "\nTarget: " << fTargetMaterialName << ", size " << fTargetSize/cm
               << " cm at " << fTargetPosition/cm << " cm" << G4endl;
    }
    if (fBuildDegrader) {
        auto solid_Degrader = new G4Box("DegraderSolid", 0.5 * fDegraderSize.x(), 0.5 * fDegraderSize.y(), 0.5 * fDegraderSize.z());
        auto logic_Degrader = new G4LogicalVolume(solid_Degrader, FindMaterial(fDegraderMaterialName), "DegraderLogic");
        logic_Degrader->SetVisAttributes(G4VisAttributes(G4Colour(0.3, 0.3, 0.3)));
        new G4PVPlacement(nullptr, fDegraderPosition, logic_Degrader, "Degrader", logic_world, false, 0, fCheckOverlaps);
        G4cout << "\nDegrader: " << fDegraderMaterialName << ", size " << fDegraderSize/cm
               << " cm at " << fDegraderPosition/cm << " cm" << G4endl;
    }
}

void DetectorConstruction::ConstructBlocks()
{
    // 隔离铅砖 or any other passive box added with /tof/blocks/add
    for (size_t i = 0; i < fBlocks.size(); ++i) {
        const auto& block = fBlocks[i];
        auto solid = new G4Box(block.name + "Solid", 0.5 * block.size.x(), 0.5 * block.size.y(), 0.5 * block.size.z());
        auto logic = new G4LogicalVolume(solid, FindMaterial(block.material), block.name + "Logic");
        logic->SetVisAttributes(G4VisAttributes(G4Colour(0.45, 0.45, 0.55, 0.6)));
        new G4PVPlacement(nullptr, block.position, logic, block.name, logic_world, false,
                          static_cast<G4int>(i), fCheckOverlaps);
        G4cout << "Block " << block.name << ": " << block.material << ", size " << block.size/cm
               << " cm at " << block.position/cm << " cm" << G4endl;
    }
}

void DetectorConstruction::AddBlock(const G4String& spec)
{
    // "name material x y z dx dy dz unit": centre and full sizes
    std::istringstream is(spec);
    Block block;
    G4double x, y, z, dx, dy, dz;
    G4String unit;
    if (!(is >> block.name >> block.material >> x >> y >> z >> dx >> dy >> dz >> unit)) {
        G4Exception("DetectorConstruction::AddBlock", "Geometry003", JustWarning,
                    ("Expected \"name material x y z dx dy dz unit\", got \"" + spec + "\"; ignored").c_str());
        return;
    }
    const G4double u = G4UIcommand::ValueOf(unit);
    block.position = G4ThreeVector(x, y, z) * u;
    block.size = G4ThreeVector(dx, dy, dz) * u;
    fBlocks.push_back(block);
}

void DetectorConstruction::ClearBlocks()
{
    fBlocks.clear();
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

    // 康普顿相机 - YSO 散射层蓝色, LYSO 吸收层品红
    if (logic_CameraYSO) {
        G4VisAttributes* camYsoVisAtt = new G4VisAttributes(G4Colour(0.2, 0.4, 1.0));
        camYsoVisAtt->SetForceSolid(true);
        logic_CameraYSO->SetVisAttributes(camYsoVisAtt);
    }
    if (logic_CameraLYSO) {
        G4VisAttributes* camLysoVisAtt = new G4VisAttributes(G4Colour(0.9, 0.1, 0.9));
        camLysoVisAtt->SetForceSolid(true);
        logic_CameraLYSO->SetVisAttributes(camLysoVisAtt);
    }
}

void DetectorConstruction::ConstructSDandField()
{
    // One TOFmeterDetectorSD (step-level ASCII output, see README.txt) for all
    // active volumes. It is reused across /update so its hits collection stays valid.
    auto sdManager = G4SDManager::GetSDMpointer();
    G4VSensitiveDetector* chamberSD = sdManager->FindSensitiveDetector("DetectorSD", false);
    if (!chamberSD) {
        chamberSD = new TOFmeterDetectorSD("DetectorSD");
        sdManager->AddNewDetector(chamberSD);
    }
    for (auto lv : chamberComponents) {
        if (lv && (lv == logic_Absorber ||
                   lv->GetName().find("Scintillator") != std::string::npos ||
                   lv == logic_YSO ||   // <--- 添加这一条件
                   lv == logic_CameraYSO || lv == logic_CameraLYSO)) {
            lv->SetSensitiveDetector(chamberSD);
        }
    }
}

DetectorKind DetectorConstruction::GetKind(const G4LogicalVolume* lv) const
{
    auto it = fKinds.find(lv);
    return it == fKinds.end() ? DetectorKind::None : it->second;
}

G4bool DetectorConstruction::IsComptonDetector(const G4LogicalVolume* lv) const
{
    const DetectorKind kind = GetKind(lv);
    return kind == DetectorKind::StackYSO || kind == DetectorKind::StackSi || kind == DetectorKind::CameraLayer;
}

void DetectorConstruction::UpdateGeometry()
{
    // Destroy the old volumes and rebuild right away (Construct + ConstructSDandField),
    // so the new geometry can be drawn or checked before the next /run/beamOn
    auto runManager = G4RunManager::GetRunManager();
    runManager->ReinitializeGeometry(true, false);
    runManager->Initialize();
}

void DetectorConstruction::DefineCommands()
{
    auto det = std::make_unique<G4GenericMessenger>(this, "/tof/det/",
        "TOF-meter geometry; changes take effect with /update");
    det->DeclareProperty("worldMaterial", fWorldMaterialName, "World material (G4_AIR or G4_Galactic, ...)");
    det->DeclarePropertyWithUnit("worldSize", "cm", fWorldSize, "Full world size x y z");
    det->DeclareProperty("stations", fBuildStations, "Build the three 9-bar plastic stations");
    det->DeclareProperty("frame", fBuildFrame, "Build the Al plates, the steel plate and the corner beams");
    det->DeclareProperty("stack", fBuildStack, "Build the Si/YSO stack");
    det->DeclareProperty("absorber", fBuildAbsorber, "Build the lead layer on top (90x90 cm, out to the lead walls)");
    det->DeclarePropertyWithUnit("leadThickness", "cm", fLeadThickness,
        "Thickness of the lead layer on the steel plate (default 10 cm, 0 = no layer); its bottom "
        "stays on the plate and the lead walls end level with its top");
    det->DeclareProperty("leadWalls", fBuildLeadWalls, "Build the four 10 cm lead walls around the setup");
    det->DeclareProperty("camera", fBuildCamera, "Build the three-layer Compton camera");
    det->DeclareMethod("print", &DetectorConstruction::PrintGeometry,
        "Print every placed volume: material, world-frame centre, size and z range (cm)");
    fMessengers.push_back(std::move(det));

    auto camera = std::make_unique<G4GenericMessenger>(this, "/tof/camera/",
        "Compton camera placement; changes take effect with /update");
    camera->DeclarePropertyWithUnit("position", "cm", fCameraPosition, "Centre of CH2 (entrance layer)");
    camera->DeclarePropertyWithUnit("rotateX", "deg", fCameraRotationX,
        "Rotation about x around CH2; 0 = looking up (+z), 90 deg = looking towards -y");
    camera->DeclarePropertyWithUnit("spacing21", "mm", fCameraSpacing21, "CH2-CH1 spacing, centre to centre");
    camera->DeclarePropertyWithUnit("spacing10", "mm", fCameraSpacing10, "CH1-CH0 spacing, centre to centre");
    fMessengers.push_back(std::move(camera));

    auto target = std::make_unique<G4GenericMessenger>(this, "/tof/target/",
        "Target (unknown material); changes take effect with /update");
    target->DeclareProperty("enable", fBuildTarget, "Build the target box");
    target->DeclareProperty("material", fTargetMaterialName, "Target material (NIST or YSO/LYSO/CdZnTe)");
    target->DeclarePropertyWithUnit("size", "cm", fTargetSize, "Full target size x y z");
    target->DeclarePropertyWithUnit("position", "cm", fTargetPosition, "Target centre");
    fMessengers.push_back(std::move(target));

    auto degrader = std::make_unique<G4GenericMessenger>(this, "/tof/degrader/",
        "Degrader slab; changes take effect with /update");
    degrader->DeclareProperty("enable", fBuildDegrader, "Build the degrader box");
    degrader->DeclareProperty("material", fDegraderMaterialName, "Degrader material (NIST name)");
    degrader->DeclarePropertyWithUnit("size", "cm", fDegraderSize, "Full degrader size x y z");
    degrader->DeclarePropertyWithUnit("position", "cm", fDegraderPosition, "Degrader centre");
    fMessengers.push_back(std::move(degrader));

    fBlockMessenger = std::make_unique<BlockMessenger>(this);
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
