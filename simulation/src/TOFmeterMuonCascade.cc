#include "TOFmeterMuonCascade.hh"

#include "TOFmeterEmCaptureCascade.hh"

#include "G4EmCaptureCascade.hh"
#include "G4GenericMessenger.hh"
#include "G4HadronStoppingProcess.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4ProcessTable.hh"
#include "G4ios.hh"

TOFmeterMuonCascade::TOFmeterMuonCascade()
  : G4HadronicInteraction("TOFmeterMuonCascade"),
    fFirstPrinciples(new TOFmeterEmCaptureCascade()),
    fGeant4(new G4EmCaptureCascade())
{
  fMessenger = std::make_unique<G4GenericMessenger>(this, "/tof/physics/",
                                                    "Physics information and choices");
  fMessenger->DeclareMethod("muonCascade", &TOFmeterMuonCascade::SetModel,
                            "Atomic cascade of stopped mu- (muonic X-rays): firstPrinciples (TOFmeterEmCaptureCascade, "
                            "default) or geant4 (G4EmCaptureCascade)")
    .SetCandidates("firstPrinciples geant4");
  fMessenger->DeclareProperty("muonCascadeZmax", fZmax,
                              "Use the first-principles cascade only for elements up to this Z, the Geant4 one above "
                              "(default 92: every element)");
}

TOFmeterMuonCascade::~TOFmeterMuonCascade() = default;

TOFmeterMuonCascade* TOFmeterMuonCascade::Install()
{
  auto process = dynamic_cast<G4HadronStoppingProcess*>(
    G4ProcessTable::GetProcessTable()->FindProcess("muMinusCaptureAtRest", "mu-"));
  if (!process) {
    G4cout << "TOFmeterMuonCascade: the physics list has no muMinusCaptureAtRest; muonic cascade not changed"
           << G4endl;
    return nullptr;
  }
  auto cascade = new TOFmeterMuonCascade();   // owned by the hadronic interaction registry
  process->SetEmCascade(cascade);
  cascade->Print();
  return cascade;
}

G4HadFinalState* TOFmeterMuonCascade::ApplyYourself(const G4HadProjectile& projectile, G4Nucleus& nucleus)
{
  const G4bool firstPrinciples = fUseFirstPrinciples && nucleus.GetZ_asInt() <= fZmax;
  return (firstPrinciples ? fFirstPrinciples : fGeant4)->ApplyYourself(projectile, nucleus);
}

void TOFmeterMuonCascade::SetModel(G4String model)
{
  fUseFirstPrinciples = (model != "geant4");
  Print();
}

void TOFmeterMuonCascade::Print() const
{
  G4cout << "mu- atomic cascade (muonic X-rays): ";
  if (!fUseFirstPrinciples) G4cout << "Geant4 G4EmCaptureCascade";
  else if (fZmax >= 92) G4cout << "first-principles TOFmeterEmCaptureCascade";
  else G4cout << "first-principles TOFmeterEmCaptureCascade up to Z = " << fZmax << ", Geant4 above";
  G4cout << " (/tof/physics/muonCascade, muonCascadeZmax)" << G4endl;
}

void TOFmeterMuonCascade::ModelDescription(std::ostream& out) const
{
  out << "Selector of the atomic cascade of stopped mu-: the first-principles (n,l)-resolved cascade "
         "(TOFmeterEmCaptureCascade) or Geant4's G4EmCaptureCascade.\n";
}
