//
// Sea-level cosmic rays for the TOF-meter simulation (/GeneralParticleSource 2,
// settings under /tof/cosmic/).
//
// Particles start on a horizontal rectangle above the setup and go downwards.
// Every event has weight 1, so a run of N particles corresponds to a live time
// N / (flux x area), with the flux through the plane integrated over the selected
// species, energy range and zenith range (GetLiveTime).
//
// Spectrum models (/tof/cosmic/model):
//  expacs  energy spectra of any species from an EXPACS table (PARMA model, T. Sato,
//          PLoS ONE 10(12): e0144679, 2015), converted with scripts/expacs_to_table.py.
//          EXPACS 3.x gives omnidirectional fluxes only, so the zenith dependence is a
//          model chosen per species (/tof/cosmic/angular): I ~ cos^n(theta), or for
//          muons the energy-dependent shape of Guan et al. The whole omnidirectional
//          flux is attributed to downward-going particles; nothing comes from below.
//  guan    muons only: Guan et al., arXiv:1509.06176 (Gaisser's formula extended to
//          low energies, curvature-corrected zenith angle of Chirkin, hep-ph/0407078),
//          mu+/mu- = /tof/cosmic/chargeRatio.
//

#ifndef CosmicRaySource_h
#define CosmicRaySource_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <vector>

class G4Event;
class G4GenericMessenger;
class G4ParticleDefinition;
class G4ParticleGun;

class CosmicRaySource
{
  public:
    CosmicRaySource();
    ~CosmicRaySource();

    void GeneratePrimaryVertex(G4Event*);

    // Check the settings and compute the flux normalisation (once per run)
    void PrepareRun();

    G4double GetFlux() const { return fFlux; }   // particles per unit area and time through the plane
    G4double GetArea() const { return 4. * fHalfX * fHalfY; }
    G4double GetLiveTime(G4double nEvents) const;

    const G4String& GetModel() const { return fModel; }
    G4String GetParticles() const;               // charge setting (guan) or species list (expacs)
    G4double GetEnergyMin() const { return fEnergyMin; }
    G4double GetEnergyMax() const { return fEnergyMax; }
    G4double GetThetaMax() const { return fThetaMax; }
    const G4ThreeVector& GetPlaneCenter() const { return fPlaneCenter; }

    // Guan et al. dI/dE in 1/(cm2 s sr GeV), total muon energy in GeV, true zenith cosine
    static G4double GuanIntensity(G4double totalEnergy, G4double cosTheta);

  private:
    // Guan model
    void PrepareGuan();
    void SampleGuan(G4ParticleDefinition*& particle, G4double& kineticEnergy, G4double& cosTheta) const;
    static G4double EnvelopeRatio(G4double totalEnergy, G4double cosTheta);
    G4double ChargeFraction() const;

    // EXPACS model
    void LoadExpacsTable();
    void PrepareExpacs();
    void SampleExpacs(G4ParticleDefinition*& particle, G4double& kineticEnergy, G4double& cosTheta) const;
    G4double SampleGuanCosTheta(G4double totalEnergyGeV) const;

    void CheckPlaneInsideWorld() const;

    G4ParticleGun* fGun;
    G4GenericMessenger* fMessenger;

    G4String fModel;          // expacs or guan
    G4String fExpacsFile;
    G4String fSpecies;        // expacs: "all" or a list of particle names
    G4String fAngular;        // expacs: pairs <species|all> <n|guan>
    G4ThreeVector fPlaneCenter;
    G4double fHalfX;
    G4double fHalfY;
    G4double fEnergyMin;      // kinetic (per nucleon for ions)
    G4double fEnergyMax;
    G4double fThetaMax;
    G4String fCharge;         // guan: both, mu-, mu+
    G4double fChargeRatio;    // guan: mu+/mu-

    G4double fFlux;           // for the settings used by PrepareRun
    G4double fCosMin;
    G4double fY1, fY2;        // guan: sampling variable (E + 3.64 GeV)^-1.7 at the energy limits
    G4bool fPrepared;

    // EXPACS table as read from fExpacsFile
    struct TableColumn {
      G4String name;
      std::vector<G4double> flux;   // omnidirectional, 1/(cm2 s MeV) per nucleon
    };
    G4String fLoadedFile;
    std::vector<G4String> fTableHeader;
    std::vector<G4double> fTableEnergy;   // kinetic energy per nucleon (MeV)
    std::vector<TableColumn> fTable;

    // Species used in this run: piecewise power-law plane flux J(E) = c E^g per segment
    struct Segment {
      G4double low, high, index;
    };
    struct ActiveSpecies {
      G4String name;
      G4ParticleDefinition* particle = nullptr;
      G4int nucleons = 1;
      G4bool guanAngular = false;
      G4double cosPower = 2.;
      std::vector<Segment> segments;
      std::vector<G4double> cumulative;     // over segments
      G4double planeFlux = 0.;              // 1/(cm2 s)
    };
    std::vector<ActiveSpecies> fActive;
    std::vector<G4double> fActiveCumulative;  // over species
};

#endif
