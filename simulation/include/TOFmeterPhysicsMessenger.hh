//
// Physics information (/tof/physics/...):
//   /tof/physics/crossSection <process> <particle> <material> <E1 E2 ... (kinetic, GeV)>
// prints the cross section of an electromagnetic process from its Geant4 model (G4EmCalculator):
// per volume, per atom (average for compounds), the mean free path and the probability over 10 cm.
// The process must be active and the physics tables built, i.e. after a /run/beamOn (0 events is
// enough).
//

#ifndef TOFmeterPhysicsMessenger_h
#define TOFmeterPhysicsMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4UIcommand;
class G4UIdirectory;

class TOFmeterPhysicsMessenger : public G4UImessenger
{
  public:
    TOFmeterPhysicsMessenger();
    ~TOFmeterPhysicsMessenger() override;

    void SetNewValue(G4UIcommand*, G4String) override;

  private:
    void PrintCrossSection(const G4String& arguments) const;

    G4UIdirectory* fDirectory;
    G4UIcommand* fCrossSectionCmd;
};

#endif
