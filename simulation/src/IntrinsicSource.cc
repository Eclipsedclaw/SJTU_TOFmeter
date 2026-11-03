//
// Radioactivity of the pixel panel's own crystal (/tof/intrinsic/): see IntrinsicSource.hh
//

#include "IntrinsicSource.hh"
#include "DetectorConstruction.hh"

#include "G4Element.hh"
#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4GenericMessenger.hh"
#include "G4IonTable.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4Poisson.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <map>

IntrinsicSource::IntrinsicSource()
  : fEnabled(false), fZ(71), fA(176), fActivitySet(0.), fWindowStart(-1000. * ns), fWindowEnd(1000. * ns), fActivity(0.),
    fIon(nullptr), fPixels(0), fPixelSize(0.), fPitch(0.), fThickness(0.)
{
  fMessenger = std::make_unique<G4GenericMessenger>(this, "/tof/intrinsic/",
    "Radioactivity of the pixel panel's crystal (/tof/czt/...), e.g. the Lu-176 of LYSO; needs Geant4's radioactive decay "
    "(TOFmeter --preinit with /tof/physics/radioactiveDecay)");
  fMessenger->DeclareProperty("enable", fEnabled,
    "Add the decays in the time window around every cosmic-ray event (/GeneralParticleSource 2); default false");
  fMessenger->DeclareProperty("Z", fZ, "Atomic number of the nuclide (default 71, Lu)");
  fMessenger->DeclareProperty("A", fA, "Mass number of the nuclide (default 176)");
  fMessenger->DeclarePropertyWithUnit("activity", "Bq", fActivitySet,
    "Decays per second in all pixels; 0 (default): from the pixels' mass, the natural abundance and the half-life");
  fMessenger->DeclarePropertyWithUnit("windowStart", "ns", fWindowStart,
    "Start of the window around a cosmic-ray event (its primary starts at time 0), default -1000 ns");
  fMessenger->DeclarePropertyWithUnit("windowEnd", "ns", fWindowEnd, "End of the window, default 1000 ns");
}

IntrinsicSource::~IntrinsicSource() = default;

void IntrinsicSource::PrepareRun(G4double cosmicFlux)
{
  const auto* detector = dynamic_cast<const Pb::DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  const auto panel = detector ? detector->GetPixelPanel() : Pb::DetectorConstruction::PixelPanel{};
  if (!detector || !panel.built || !panel.material) {
    G4Exception("IntrinsicSource::PrepareRun", "Intrinsic001", FatalErrorInArgument,
                "/tof/intrinsic/ needs the pixel panel (/tof/czt/enable true, then /update)");
    return;
  }
  fPixels = panel.pixels; fPixelSize = panel.pixelSize; fPitch = panel.pitch; fThickness = panel.thickness;
  fCentre = panel.centre;
  fIon = G4IonTable::GetIonTable()->GetIon(fZ, fA, 0.);
  if (!fIon) {
    G4ExceptionDescription ed;
    ed << "no nuclide Z = " << fZ << ", A = " << fA;
    G4Exception("IntrinsicSource::PrepareRun", "Intrinsic002", FatalErrorInArgument, ed);
    return;
  }
  // the nuclide's mean life as Geant4 knows it (G4NuclideTable), kept from before it is shortened below. Not
  // GetPDGStable(): G4IonTable flags every ion without a decay table as stable, yet keeps its lifetime, and
  // G4RadioactiveDecay looks only at the lifetime (it skips those above 1 year)
  static std::map<const G4ParticleDefinition*, G4double> meanLife;
  if (!meanLife.count(fIon)) meanLife[fIon] = fIon->GetPDGLifeTime() > 0. ? fIon->GetPDGLifeTime() : DBL_MAX;
  const G4double tau = meanLife[fIon];
  // atoms of the nuclide per volume of the pixels' material: each element of atomic number Z with its isotopes
  const G4Material* material = panel.material;
  G4double atoms = 0.;
  for (size_t i = 0; i < material->GetNumberOfElements(); ++i) {
    const G4Element* element = material->GetElement(i);
    if (element->GetZasInt() != fZ) continue;
    for (size_t k = 0; k < element->GetNumberOfIsotopes(); ++k)
      if (element->GetIsotope(k)->GetN() == fA)
        atoms += material->GetVecNbOfAtomsPerVolume()[i] * element->GetRelativeAbundanceVector()[k];
  }
  const G4double volume = fPixels * fPixels * fPixelSize * fPixelSize * fThickness;
  const G4double mass = volume * material->GetDensity();
  const G4double natural = (tau > 0. && tau < DBL_MAX) ? atoms * volume / tau : 0.;
  fActivity = fActivitySet > 0. ? fActivitySet : natural;
  // decay where and when placed: the lifetime shortened to 1 ps
  fIon->SetPDGStable(false);
  fIon->SetPDGLifeTime(1.e-3 * ns);

  const G4double window = fWindowEnd - fWindowStart;
  G4cout << "\nIntrinsic radioactivity of the " << panel.name << " panel (/tof/intrinsic/): " << fIon->GetParticleName()
         << ", half-life " << (tau < DBL_MAX ? tau * std::log(2.) / (365.25 * 86400. * s) : -1.) << " years; "
         << atoms * cm3 << " atoms/cm3 in " << material->GetName() << "; pixels " << fPixels << " x " << fPixels << " of "
         << fPixelSize / mm << " x " << fPixelSize / mm << " x " << fThickness / mm << " mm, " << volume / cm3 << " cm3, "
         << mass / g << " g\n  activity " << fActivity / becquerel << " Bq (" << fActivity / becquerel / (mass / g)
         << " Bq/g; from the material " << natural / becquerel << " Bq" << (fActivitySet > 0. ? ", set by hand" : "")
         << ")" << G4endl;
  if (cosmicFlux > 0.) {
    const G4double area = (fPixels * fPitch) * (fPixels * fPitch);
    const G4double muons = cosmicFlux * area;
    G4cout << "  cosmic rays: " << cosmicFlux * cm2 * minute << " /cm2/min, through the panel's " << area / cm2 << " cm2 "
           << muons * s << " /s: " << fActivity / muons << " decays per cosmic ray through the panel; per event "
           << fActivity * window << " decays on average in the window " << fWindowStart / ns << " to " << fWindowEnd / ns
           << " ns around it" << G4endl;
  }
}

void IntrinsicSource::AddDecays(G4Event* event)
{
  if (!fEnabled || fActivity <= 0. || !fIon) return;
  const G4double window = fWindowEnd - fWindowStart;
  const G4long n = G4Poisson(fActivity * window);
  for (G4long i = 0; i < n; ++i) AddDecay(event, fWindowStart + window * G4UniformRand());
}

void IntrinsicSource::AddOneDecay(G4Event* event)
{
  if (fIon) AddDecay(event, 0.);
}

void IntrinsicSource::AddDecay(G4Event* event, G4double time)
{
  // a pixel at random (all of the same volume), a point at random in it
  const G4int ix = std::min(fPixels - 1, static_cast<G4int>(fPixels * G4UniformRand()));
  const G4int iy = std::min(fPixels - 1, static_cast<G4int>(fPixels * G4UniformRand()));
  const G4ThreeVector position = fCentre + G4ThreeVector((ix - 0.5 * (fPixels - 1)) * fPitch + (G4UniformRand() - 0.5) * fPixelSize,
                                                         (iy - 0.5 * (fPixels - 1)) * fPitch + (G4UniformRand() - 0.5) * fPixelSize,
                                                         (G4UniformRand() - 0.5) * fThickness);
  auto vertex = new G4PrimaryVertex(position, time);
  auto nucleus = new G4PrimaryParticle(fIon, 0., 0., 0.);   // at rest
  nucleus->SetCharge(0.);
  nucleus->SetUserInformation(new IntrinsicDecayTag);
  vertex->SetPrimary(nucleus);
  event->AddPrimaryVertex(vertex);
}
