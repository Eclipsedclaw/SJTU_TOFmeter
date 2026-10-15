#include "TOFmeterPhysicsMessenger.hh"

#include "G4EmCalculator.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"

#include <cmath>
#include <cstdio>
#include <sstream>

TOFmeterPhysicsMessenger::TOFmeterPhysicsMessenger()
{
  fDirectory = new G4UIdirectory("/tof/physics/");
  fDirectory->SetGuidance("Physics information");

  fCrossSectionCmd = new G4UIcommand("/tof/physics/crossSection", this);
  fCrossSectionCmd->SetGuidance("Cross section of an EM process from its model: process particle material E1 E2 ...");
  fCrossSectionCmd->SetGuidance("Energies are kinetic, in GeV. Needs the physics tables (after /run/beamOn, 0 is enough).");
  for (const char* name : {"process", "particle", "material"}) {
    fCrossSectionCmd->SetParameter(new G4UIparameter(name, 's', false));
  }
  auto energies = new G4UIparameter("energies", 's', true);  // the last string takes the rest of the line
  energies->SetDefaultValue("1 2 5 10 20 50 100");
  fCrossSectionCmd->SetParameter(energies);
  fCrossSectionCmd->AvailableForStates(G4State_Idle);
}

TOFmeterPhysicsMessenger::~TOFmeterPhysicsMessenger()
{
  delete fCrossSectionCmd;
  delete fDirectory;
}

void TOFmeterPhysicsMessenger::SetNewValue(G4UIcommand* command, G4String value)
{
  if (command == fCrossSectionCmd) PrintCrossSection(value);
}

void TOFmeterPhysicsMessenger::PrintCrossSection(const G4String& arguments) const
{
  std::istringstream in(arguments);
  G4String processName, particleName, materialName;
  in >> processName >> particleName >> materialName;
  const G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
  const G4Material* material = G4Material::GetMaterial(materialName, false);
  if (!material) material = G4NistManager::Instance()->FindOrBuildMaterial(materialName);
  if (!particle || !material) {
    G4cout << "/tof/physics/crossSection: unknown particle " << particleName << " or material " << materialName
           << G4endl;
    return;
  }

  G4EmCalculator calculator;
  const G4double atoms = material->GetTotNbOfAtomsPerVolume();
  G4cout << "\nCross section of " << processName << " for " << particleName << " in " << materialName << " ("
         << material->GetDensity() / (g / cm3) << " g/cm3, " << material->GetNumberOfElements()
         << " element(s)), from the Geant4 model\n"
         << "  E_kin (GeV)   per atom (ub)   per volume (1/cm)   mean free path (m)   P(10 cm)" << G4endl;
  G4double energy = 0.;
  while (in >> energy) {
    const G4double perVolume =
      calculator.ComputeCrossSectionPerVolume(energy * GeV, particle, processName, material, 0.);
    const G4double perAtom = atoms > 0. ? perVolume / atoms : 0.;
    char line[160];
    std::snprintf(line, sizeof(line), "  %11.4g   %13.5g   %17.5g   %18.5g   %8.3g", energy, perAtom / microbarn,
                  perVolume * cm, perVolume > 0. ? 1. / perVolume / m : 0., 1. - std::exp(-perVolume * 10. * cm));
    G4cout << line << G4endl;
  }
}
