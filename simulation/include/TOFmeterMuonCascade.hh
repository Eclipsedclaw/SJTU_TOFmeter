//
// Atomic cascade of a mu- stopped and captured at rest (process muMinusCaptureAtRest): the model
// that emits the muonic X-rays and Auger electrons. Install() puts this selector into the mu-
// capture process of the physics list (G4HadronStoppingProcess::SetEmCascade), after the kernel
// is initialised. It passes each capture to one of two models:
//   firstPrinciples   TOFmeterEmCaptureCascade: (n,l)-resolved cascade from n = 14 with calculated
//                     radiative and K/L Auger rates (default)
//   geant4            G4EmCaptureCascade, the model of the Geant4 libraries
// /tof/physics/muonCascade firstPrinciples|geant4 chooses; /tof/physics/muonCascadeZmax Z uses the
// first-principles model only for elements up to Z and the Geant4 model above (default 92: all).
// Other negative particles stopped at rest keep their own processes and models.
//
// @author Jiancheng Zeng with Claude Opus 5.5
// @date October 2026
//

#ifndef TOFmeterMuonCascade_h
#define TOFmeterMuonCascade_h 1

#include "G4HadronicInteraction.hh"
#include "globals.hh"

#include <memory>

class G4GenericMessenger;

class TOFmeterMuonCascade : public G4HadronicInteraction
{
  public:
    TOFmeterMuonCascade();
    ~TOFmeterMuonCascade() override;

    // Puts a TOFmeterMuonCascade into the mu- capture process (after runManager->Initialize());
    // returns nullptr if the physics list has no such process
    static TOFmeterMuonCascade* Install();

    G4HadFinalState* ApplyYourself(const G4HadProjectile& projectile, G4Nucleus& nucleus) override;
    void ModelDescription(std::ostream& out) const override;

  private:
    void SetModel(G4String model);
    void Print() const;

    G4HadronicInteraction* fFirstPrinciples;   // both owned by the hadronic interaction registry
    G4HadronicInteraction* fGeant4;
    G4bool fUseFirstPrinciples = true;
    G4int fZmax = 92;
    std::unique_ptr<G4GenericMessenger> fMessenger;
};

#endif
