//
// Radioactivity of the pixel panel's own crystal (/tof/intrinsic/), e.g. the Lu-176 of LYSO (2026-11-03): decays of a
// nuclide at rest at random points of the panel's pixels (/tof/czt/...). With /GeneralParticleSource 2 every cosmic-ray
// event also gets the decays that fall into a time window around it (a Poisson number with mean activity x window, at
// times uniform in the window), so that random coincidences with the TOF trigger and pile-up come out at the true rate;
// with /GeneralParticleSource 3 every event is one decay at time 0 alone (the crystal's own spectrum; live time = events
// / activity). The activity follows from the pixels' mass, the natural isotopic abundance of the element (Geant4's NIST
// element) and the nuclide's half-life; /tof/intrinsic/activity overrides it. The nuclide decays by Geant4's radioactive
// decay (TOFmeter --preinit with /tof/physics/radioactiveDecay); its lifetime is set to 1 ps so that it decays where and
// when it is placed (Geant4 would sample its decay time from the 10^10-year lifetime, or take the nuclide as stable).
//

#ifndef IntrinsicSource_h
#define IntrinsicSource_h 1

#include "G4ThreeVector.hh"
#include "G4VUserPrimaryParticleInformation.hh"
#include "globals.hh"

#include <memory>

class G4Event;
class G4GenericMessenger;
class G4ParticleDefinition;

// marks the primaries of the intrinsic source (TOFmeterEventAction: source bit kIntrinsic of the hits)
class IntrinsicDecayTag : public G4VUserPrimaryParticleInformation
{
  public:
    void Print() const override {}
};

class IntrinsicSource
{
  public:
    IntrinsicSource();
    ~IntrinsicSource();

    G4bool IsEnabled() const { return fEnabled; }
    // the nuclide, the panel and the activity, printed beside the cosmic rays' flux (cosmicFlux: per area and time, 0
    // without them)
    void PrepareRun(G4double cosmicFlux);
    void AddDecays(G4Event* event);     // the decays in the window around a cosmic-ray event
    void AddOneDecay(G4Event* event);   // one decay at time 0 (/GeneralParticleSource 3)
    G4double GetActivity() const { return fActivity; }
    G4double GetLiveTime(G4double nEvents) const { return fActivity > 0. ? nEvents / fActivity : 0.; }

  private:
    void AddDecay(G4Event* event, G4double time);

    G4bool fEnabled;
    G4int fZ, fA;                 // the nuclide
    G4double fActivitySet;        // /tof/intrinsic/activity; 0: from the material
    G4double fWindowStart, fWindowEnd;
    G4double fActivity;           // decays per unit time in the pixels
    G4ParticleDefinition* fIon;
    // the panel's pixels
    G4int fPixels;
    G4double fPixelSize, fPitch, fThickness;
    G4ThreeVector fCentre;
    std::unique_ptr<G4GenericMessenger> fMessenger;
};

#endif
