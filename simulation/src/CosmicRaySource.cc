#include "CosmicRaySource.hh"

#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "G4IonTable.hh"
#include "G4LogicalVolume.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4Navigator.hh"
#include "G4NistManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

namespace
{
// Guan et al., arXiv:1509.06176, Eq. (3): E in GeV, intensity in 1/(cm2 s sr GeV)
constexpr G4double kNorm = 0.14;
constexpr G4double kIndex = 2.7;
constexpr G4double kLowE = 3.64;
constexpr G4double kLowExp = 1.29;
constexpr G4double kSecondTermMax = 1.054;  // bound of the [1/(1+...) + 0.054/(1+...)] factor

// Curvature-corrected cos(theta), Guan et al. Eq. (2) (Chirkin's parametrisation)
constexpr G4double kP1 = 0.102573;
constexpr G4double kP2 = -0.068287;
constexpr G4double kP3 = 0.958633;
constexpr G4double kP4 = 0.0407253;
constexpr G4double kP5 = 0.817285;

G4double CosThetaStar(G4double c)
{
  const G4double num = c * c + kP1 * kP1 + kP2 * std::pow(c, kP3) + kP4 * std::pow(c, kP5);
  const G4double den = 1. + kP1 * kP1 + kP2 + kP4;
  return std::sqrt(num / den);
}

// Sampling variable y = (E + 3.64)^-1.7: the envelope 0.14 (E + 3.64)^-2.7 is flat in y
G4double ToY(G4double totalEnergy) { return std::pow(totalEnergy + kLowE, 1. - kIndex); }
G4double FromY(G4double y) { return std::pow(y, 1. / (1. - kIndex)) - kLowE; }

// Integral over [a, b] of j0 (E/e0)^g, written relative to e0 so that steep
// spectra cannot overflow
G4double SegmentIntegral(G4double j0, G4double e0, G4double g, G4double a, G4double b)
{
  if (std::abs(g + 1.) < 1e-9) return j0 * e0 * std::log(b / a);
  return j0 * e0 * (std::pow(b / e0, g + 1.) - std::pow(a / e0, g + 1.)) / (g + 1.);
}

// Sample of a power law E^g on [a, b]
G4double SamplePowerLaw(G4double g, G4double a, G4double b, G4double u)
{
  const G4double r = b / a;
  if (std::abs(g + 1.) < 1e-9) return a * std::pow(r, u);
  return a * std::pow(1. + u * (std::pow(r, g + 1.) - 1.), 1. / (g + 1.));
}

// Simpson's rule with n (even) intervals
template <typename F>
G4double Simpson(F f, G4double a, G4double b, G4int n)
{
  const G4double h = (b - a) / n;
  G4double sum = f(a) + f(b);
  for (G4int i = 1; i < n; ++i) sum += f(a + i * h) * ((i % 2) ? 4. : 2.);
  return sum * h / 3.;
}

// Particle of an EXPACS column: a Geant4 name, or "<symbol><A>" for an ion
G4ParticleDefinition* FindSpecies(const G4String& name, G4int& nucleons)
{
  if (auto particle = G4ParticleTable::GetParticleTable()->FindParticle(name)) {
    nucleons = std::max(1, particle->GetBaryonNumber());
    return particle;
  }
  std::size_t split = 0;
  while (split < name.size() && std::isalpha(static_cast<unsigned char>(name[split]))) ++split;
  if (split == 0 || split == name.size()) return nullptr;
  const G4int z = G4NistManager::Instance()->GetZ(name.substr(0, split));
  const G4int a = std::atoi(name.substr(split).c_str());
  if (z <= 0 || a < z) return nullptr;
  nucleons = a;
  return G4IonTable::GetIonTable()->GetIon(z, a, 0.);
}

G4bool IsMuon(const G4String& name) { return name == "mu+" || name == "mu-"; }
}  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CosmicRaySource::CosmicRaySource()
  : fGun(new G4ParticleGun(1)),
    fMessenger(nullptr),
    fModel("expacs"),
    fExpacsFile("data/expacs_spectrum.txt"),
    fSpecies("mu+ mu-"),
    fAngular("all 2 mu+ guan mu- guan"),
    fPlaneCenter(0., 0., 75. * cm),
    fHalfX(150. * cm),
    fHalfY(150. * cm),
    fEnergyMin(0.),
    fEnergyMax(1000. * GeV),
    fMonoEnergy(0.),
    fThetaMax(75. * deg),
    fCharge("both"),
    fChargeRatio(1.27),
    fFlux(0.),
    fCosMin(0.),
    fY1(0.),
    fY2(0.),
    fPrepared(false)
{
  fMessenger = new G4GenericMessenger(this, "/tof/cosmic/", "Cosmic rays, used with /GeneralParticleSource 2");
  fMessenger->DeclareProperty("model", fModel,
                              "Spectrum: expacs (any species, from an EXPACS table) or guan (muons, formula)")
    .SetCandidates("expacs guan");
  fMessenger->DeclareProperty("expacsFile", fExpacsFile,
                              "EXPACS table written by scripts/expacs_to_table.py");
  fMessenger->DeclareProperty("species", fSpecies,
                              "expacs: particles to generate, \"all\" or names such as mu+ mu- neutron proton "
                              "e- e+ gamma alpha");
  fMessenger->DeclareProperty("angular", fAngular,
                              "expacs: zenith laws as pairs <species|all> <n|guan>, I ~ cos^n(theta); "
                              "guan = energy-dependent muon shape of Guan et al.; all sets the default");
  fMessenger->DeclarePropertyWithUnit("planeCenter", "cm", fPlaneCenter,
                                      "Centre of the horizontal source plane");
  fMessenger->DeclarePropertyWithUnit("halfX", "cm", fHalfX, "Half length of the plane along x");
  fMessenger->DeclarePropertyWithUnit("halfY", "cm", fHalfY, "Half length of the plane along y");
  fMessenger->DeclarePropertyWithUnit("energyMin", "GeV", fEnergyMin,
                                      "Minimum kinetic energy (per nucleon for ions)");
  fMessenger->DeclarePropertyWithUnit("energyMax", "GeV", fEnergyMax,
                                      "Maximum kinetic energy (per nucleon for ions)");
  fMessenger->DeclarePropertyWithUnit("monoEnergy", "GeV", fMonoEnergy,
                                      "If > 0, kinetic energy (per nucleon for ions) of every particle, "
                                      "with the zenith distribution of that energy; the flux and live time "
                                      "still come from the spectrum between energyMin and energyMax, the "
                                      "band it stands for. 0 = sample the spectrum");
  fMessenger->DeclarePropertyWithUnit("thetaMax", "deg", fThetaMax, "Maximum zenith angle");
  fMessenger->DeclareProperty("charge", fCharge, "guan: muons to generate, both, mu- or mu+")
    .SetCandidates("both mu- mu+");
  fMessenger->DeclareProperty("chargeRatio", fChargeRatio, "guan: mu+/mu- ratio, used with charge = both");
}

CosmicRaySource::~CosmicRaySource()
{
  delete fMessenger;
  delete fGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double CosmicRaySource::GuanIntensity(G4double totalEnergy, G4double cosTheta)
{
  const G4double cs = CosThetaStar(cosTheta);
  const G4double e = totalEnergy;
  return kNorm * std::pow(e * (1. + kLowE / (e * std::pow(cs, kLowExp))), -kIndex)
         * (1. / (1. + 1.1 * e * cs / 115.) + 0.054 / (1. + 1.1 * e * cs / 850.));
}

G4double CosmicRaySource::EnvelopeRatio(G4double totalEnergy, G4double cosTheta)
{
  return GuanIntensity(totalEnergy, cosTheta) * cosTheta
         / (kNorm * std::pow(totalEnergy + kLowE, -kIndex));
}

G4double CosmicRaySource::ChargeFraction() const
{
  if (fCharge == "mu-") return 1. / (1. + fChargeRatio);
  if (fCharge == "mu+") return fChargeRatio / (1. + fChargeRatio);
  return 1.;
}

G4double CosmicRaySource::GetLiveTime(G4double nEvents) const
{
  return (fFlux > 0.) ? nEvents / (fFlux * GetArea()) : 0.;
}

G4String CosmicRaySource::GetParticles() const
{
  return fModel == "guan" ? fCharge : fSpecies;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CosmicRaySource::PrepareRun()
{
  if (fEnergyMin < 0. || fEnergyMax <= fEnergyMin || fThetaMax <= 0. || fThetaMax > 90. * deg
      || fHalfX <= 0. || fHalfY <= 0.)
  {
    G4Exception("CosmicRaySource::PrepareRun", "Cosmic001", FatalErrorInArgument,
                "Invalid /tof/cosmic settings (energy range, thetaMax or plane size)");
  }
  if (fMonoEnergy < 0. || (fMonoEnergy > 0. && (fMonoEnergy < fEnergyMin || fMonoEnergy > fEnergyMax))) {
    G4Exception("CosmicRaySource::PrepareRun", "Cosmic011", FatalErrorInArgument,
                "/tof/cosmic/monoEnergy must lie between energyMin and energyMax (0 = off)");
  }
  fCosMin = std::cos(fThetaMax);

  G4cout << "\n=== Cosmic rays (" << fModel << ") ===" << G4endl;
  if (fModel == "guan") PrepareGuan();
  else PrepareExpacs();
  CheckPlaneInsideWorld();

  const G4double fluxPerCm2s = fFlux * cm2 * s;
  G4cout << "  plane centre " << fPlaneCenter / cm << " cm, half sizes " << fHalfX / cm << " x "
         << fHalfY / cm << " cm, zenith < " << fThetaMax / deg << " deg" << G4endl
         << "  flux through the plane " << fluxPerCm2s << " /cm2/s = " << fluxPerCm2s * 60.
         << " /cm2/min; " << fluxPerCm2s * GetArea() / cm2 << " particles/s on the plane" << G4endl;
  if (fMonoEnergy > 0.) {
    G4cout << "  mono-energetic: every particle has " << fMonoEnergy / GeV
           << " GeV (per nucleon for ions), normalised to the flux above" << G4endl;
  }
  fPrepared = true;
}

void CosmicRaySource::PrepareGuan()
{
  if (fChargeRatio < 0.) {
    G4Exception("CosmicRaySource::PrepareGuan", "Cosmic003", FatalErrorInArgument,
                "/tof/cosmic/chargeRatio must not be negative");
  }
  const G4double mass = G4MuonMinus::Definition()->GetPDGMass() / GeV;
  fY1 = ToY(fEnergyMin / GeV + mass);
  fY2 = ToY(fEnergyMax / GeV + mass);

  // Flux through the plane = 2 pi Int dcos Int dE I(E, theta) cos(theta);
  // in y the integrand is EnvelopeRatio * 0.14/1.7, smooth enough for Simpson's rule
  const G4int nY = 400;
  const G4int nC = 200;
  G4double sum = 0.;
  for (G4int ic = 0; ic <= nC; ++ic) {
    const G4double c = fCosMin + (1. - fCosMin) * ic / nC;
    const G4double wc = (ic == 0 || ic == nC) ? 1. : ((ic % 2) ? 4. : 2.);
    for (G4int iy = 0; iy <= nY; ++iy) {
      const G4double y = fY2 + (fY1 - fY2) * iy / nY;
      const G4double wy = (iy == 0 || iy == nY) ? 1. : ((iy % 2) ? 4. : 2.);
      sum += wc * wy * EnvelopeRatio(FromY(y), c);
    }
  }
  const G4double integral = sum * ((1. - fCosMin) / nC / 3.) * ((fY1 - fY2) / nY / 3.);
  fFlux = twopi * kNorm / (kIndex - 1.) * integral * ChargeFraction() / (cm2 * s);

  G4cout << "  Guan et al. muons, kinetic energy " << fEnergyMin / GeV << " - " << fEnergyMax / GeV
         << " GeV, charge " << fCharge << " (mu+/mu- = " << fChargeRatio << ")" << G4endl;
}

void CosmicRaySource::LoadExpacsTable()
{
  if (fLoadedFile == fExpacsFile && !fTable.empty()) return;

  std::ifstream in(fExpacsFile);
  if (!in) {
    G4Exception("CosmicRaySource::LoadExpacsTable", "Cosmic004", FatalErrorInArgument,
                ("Cannot open the EXPACS table " + fExpacsFile
                 + " (run from the build directory or set /tof/cosmic/expacsFile)").c_str());
  }
  fTableHeader.clear();
  fTableEnergy.clear();
  fTable.clear();
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') {
      const std::size_t pos = line.find("columns:");
      if (pos == std::string::npos) {
        fTableHeader.push_back(line.substr(std::min<std::size_t>(2, line.size())));
        continue;
      }
      std::istringstream names(line.substr(pos + 8));
      std::string name;
      names >> name;  // energy column
      while (names >> name) fTable.push_back({name, {}});
      continue;
    }
    std::istringstream values(line);
    G4double energy = 0.;
    if (!(values >> energy)) continue;
    fTableEnergy.push_back(energy * MeV);
    for (auto& column : fTable) {
      G4double phi = 0.;
      values >> phi;
      column.flux.push_back(phi);
    }
  }
  if (fTable.empty() || fTableEnergy.size() < 2) {
    G4Exception("CosmicRaySource::LoadExpacsTable", "Cosmic005", FatalErrorInArgument,
                ("No '# columns:' line or no data in " + fExpacsFile).c_str());
  }
  fLoadedFile = fExpacsFile;
}

void CosmicRaySource::PrepareExpacs()
{
  LoadExpacsTable();

  // Species to generate
  std::set<G4String> wanted;
  {
    std::istringstream is(fSpecies);
    G4String name;
    while (is >> name) wanted.insert(name);
  }
  const G4bool all = wanted.count("all") > 0;
  for (const auto& name : wanted) {
    const G4bool known = std::any_of(fTable.begin(), fTable.end(),
                                     [&](const TableColumn& column) { return column.name == name; });
    if (name != "all" && !known) {
      G4Exception("CosmicRaySource::PrepareExpacs", "Cosmic006", FatalErrorInArgument,
                  ("Species " + name + " is not a column of " + fExpacsFile).c_str());
    }
  }

  // Zenith laws: "all n" sets the default, entries for single species override it
  G4double defaultPower = 2.;
  std::map<G4String, G4double> power;
  std::set<G4String> guanShape;
  {
    std::istringstream is(fAngular);
    G4String who, law;
    while (is >> who >> law) {
      if (law == "guan") {
        if (!IsMuon(who) && who != "all") {
          G4Exception("CosmicRaySource::PrepareExpacs", "Cosmic007", FatalErrorInArgument,
                      ("/tof/cosmic/angular: guan is for mu+ and mu- only, not " + who).c_str());
        }
        if (who == "all") guanShape.insert({"mu+", "mu-"});
        else guanShape.insert(who);
        continue;
      }
      char* end = nullptr;
      const G4double n = std::strtod(law.c_str(), &end);
      if (end == law.c_str() || *end != '\0' || n < 0.) {
        G4Exception("CosmicRaySource::PrepareExpacs", "Cosmic008", FatalErrorInArgument,
                    ("/tof/cosmic/angular: expected a power n >= 0 or guan, got " + law).c_str());
      }
      if (who == "all") defaultPower = n;
      else {
        power[who] = n;
        guanShape.erase(who);
      }
    }
  }

  fActive.clear();
  fActiveCumulative.clear();
  G4double total = 0.;
  for (const auto& column : fTable) {
    if (!all && !wanted.count(column.name)) continue;
    ActiveSpecies species;
    species.name = column.name;
    species.particle = FindSpecies(column.name, species.nucleons);
    if (!species.particle) {
      G4Exception("CosmicRaySource::PrepareExpacs", "Cosmic009", FatalErrorInArgument,
                  ("No Geant4 particle for the EXPACS column " + column.name).c_str());
    }
    species.guanAngular = IsMuon(species.name) && guanShape.count(species.name);
    species.cosPower = power.count(species.name) ? power[species.name] : defaultPower;

    // Fraction of the omnidirectional flux that crosses the plane downwards within thetaMax
    auto throughPlane = [&](G4double energyPerNucleon) -> G4double {
      if (species.guanAngular) {
        const G4double e = (energyPerNucleon + species.particle->GetPDGMass()) / GeV;
        const G4double down = Simpson([&](G4double c) { return GuanIntensity(e, c) * c; }, fCosMin, 1., 200);
        const G4double omni = Simpson([&](G4double c) { return GuanIntensity(e, c); }, 0., 1., 200);
        return omni > 0. ? down / omni : 0.;
      }
      const G4double n = species.cosPower;
      return (n + 1.) * (1. - std::pow(fCosMin, n + 2.)) / (n + 2.);
    };

    // Piecewise power law through the table points, cut to the energy range
    std::vector<G4double> plane(fTableEnergy.size(), 0.);
    for (std::size_t i = 0; i < fTableEnergy.size(); ++i) {
      if (column.flux[i] > 0.) plane[i] = column.flux[i] * throughPlane(fTableEnergy[i]);
    }
    G4double sum = 0.;
    for (std::size_t i = 0; i + 1 < fTableEnergy.size(); ++i) {
      if (plane[i] <= 0. || plane[i + 1] <= 0.) continue;
      const G4double a = std::max(fTableEnergy[i], fEnergyMin);
      const G4double b = std::min(fTableEnergy[i + 1], fEnergyMax);
      if (a >= b) continue;
      const G4double g = std::log(plane[i + 1] / plane[i]) / std::log(fTableEnergy[i + 1] / fTableEnergy[i]);
      sum += SegmentIntegral(plane[i], fTableEnergy[i], g, a, b);  // energies in MeV: 1/(cm2 s)
      species.segments.push_back({a, b, g});
      species.cumulative.push_back(sum);
    }
    species.planeFlux = sum;
    if (sum > 0.) {
      total += sum;
      fActive.push_back(species);
      fActiveCumulative.push_back(total);
    }
  }
  if (fActive.empty()) {
    G4Exception("CosmicRaySource::PrepareExpacs", "Cosmic010", FatalErrorInArgument,
                "No selected EXPACS species has flux in the energy range");
  }
  fFlux = total / (cm2 * s);

  G4cout << "  EXPACS table " << fExpacsFile << G4endl;
  for (const auto& line : fTableHeader) G4cout << "    " << line << G4endl;
  G4cout << "  kinetic energy " << fEnergyMin / MeV << " - " << fEnergyMax / MeV
         << " MeV (per nucleon for ions)" << G4endl
         << "  species     zenith law    flux through the plane (/cm2/s)   share" << G4endl;
  for (const auto& species : fActive) {
    char law[32];
    if (species.guanAngular) std::snprintf(law, sizeof(law), "guan");
    else std::snprintf(law, sizeof(law), "cos^%g", species.cosPower);
    char row[128];
    std::snprintf(row, sizeof(row), "  %-10s  %-12s  %13.5g                     %6.2f%%",
                  species.name.c_str(), law, species.planeFlux, 100. * species.planeFlux / total);
    G4cout << row << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CosmicRaySource::CheckPlaneInsideWorld() const
{
  auto navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  const G4VPhysicalVolume* world = navigator ? navigator->GetWorldVolume() : nullptr;
  if (!world) return;
  const G4VSolid* solid = world->GetLogicalVolume()->GetSolid();
  for (G4double sx : {-1., 1.}) {
    for (G4double sy : {-1., 1.}) {
      const G4ThreeVector corner = fPlaneCenter + G4ThreeVector(sx * fHalfX, sy * fHalfY, 0.);
      if (solid->Inside(corner) == kOutside) {
        G4Exception("CosmicRaySource::PrepareRun", "Cosmic002", FatalErrorInArgument,
                    "The cosmic source plane extends outside the world: enlarge /tof/det/worldSize "
                    "(then /update) or shrink /tof/cosmic/halfX, halfY");
      }
    }
  }
}

G4double CosmicRaySource::SampleGuanCosTheta(G4double totalEnergyGeV) const
{
  // Accept-reject in cos(theta) for the flux through the plane, I(E, theta) cos(theta),
  // under the bound 0.14 x 1.054 (E + 3.64)^-2.7 that holds for every angle
  const G4double bound = kNorm * kSecondTermMax * std::pow(totalEnergyGeV + kLowE, -kIndex);
  G4double cosTheta = 1.;
  do {
    cosTheta = fCosMin + (1. - fCosMin) * G4UniformRand();
  } while (G4UniformRand() * bound > GuanIntensity(totalEnergyGeV, cosTheta) * cosTheta);
  return cosTheta;
}

void CosmicRaySource::SampleGuan(G4ParticleDefinition*& particle, G4double& kineticEnergy,
                                 G4double& cosTheta) const
{
  const G4double mass = G4MuonMinus::Definition()->GetPDGMass();
  G4double energy = 0.;   // total, GeV
  if (fMonoEnergy > 0.) {
    energy = (fMonoEnergy + mass) / GeV;
    cosTheta = SampleGuanCosTheta(energy);
  }
  else {
    // Accept-reject against the envelope (flat in y and in cos(theta))
    do {
      energy = FromY(fY1 - G4UniformRand() * (fY1 - fY2));
      cosTheta = fCosMin + (1. - fCosMin) * G4UniformRand();
    } while (G4UniformRand() * kSecondTermMax > EnvelopeRatio(energy, cosTheta));
  }

  const G4bool positive =
    (fCharge == "mu+") || (fCharge == "both" && G4UniformRand() < fChargeRatio / (1. + fChargeRatio));
  particle = positive ? static_cast<G4ParticleDefinition*>(G4MuonPlus::Definition())
                      : G4MuonMinus::Definition();
  kineticEnergy = (fMonoEnergy > 0.) ? fMonoEnergy : std::max(0., energy * GeV - mass);
}

void CosmicRaySource::SampleExpacs(G4ParticleDefinition*& particle, G4double& kineticEnergy,
                                   G4double& cosTheta) const
{
  const G4double u = G4UniformRand() * fActiveCumulative.back();
  const std::size_t is = std::min<std::size_t>(
    fActive.size() - 1,
    std::upper_bound(fActiveCumulative.begin(), fActiveCumulative.end(), u) - fActiveCumulative.begin());
  const ActiveSpecies& species = fActive[is];

  G4double energyPerNucleon = fMonoEnergy;
  if (fMonoEnergy <= 0.) {
    const G4double v = G4UniformRand() * species.cumulative.back();
    const std::size_t k = std::min<std::size_t>(
      species.segments.size() - 1,
      std::upper_bound(species.cumulative.begin(), species.cumulative.end(), v) - species.cumulative.begin());
    const Segment& segment = species.segments[k];
    energyPerNucleon = SamplePowerLaw(segment.index, segment.low, segment.high, G4UniformRand());
  }

  if (species.guanAngular) {
    cosTheta = SampleGuanCosTheta((energyPerNucleon + species.particle->GetPDGMass()) / GeV);
  }
  else {
    // flux through the plane ~ cos^(n+1)(theta) dcos(theta)
    const G4double p = species.cosPower + 2.;
    const G4double lowest = std::pow(fCosMin, p);
    cosTheta = std::pow(lowest + G4UniformRand() * (1. - lowest), 1. / p);
  }
  particle = species.particle;
  kineticEnergy = energyPerNucleon * species.nucleons;
}

void CosmicRaySource::GeneratePrimaryVertex(G4Event* event)
{
  if (!fPrepared) PrepareRun();

  G4ParticleDefinition* particle = nullptr;
  G4double kineticEnergy = 0.;
  G4double cosTheta = 1.;
  if (fModel == "guan") SampleGuan(particle, kineticEnergy, cosTheta);
  else SampleExpacs(particle, kineticEnergy, cosTheta);

  const G4double sinTheta = std::sqrt(std::max(0., 1. - cosTheta * cosTheta));
  const G4double phi = twopi * G4UniformRand();
  fGun->SetParticleDefinition(particle);
  fGun->SetParticleEnergy(kineticEnergy);
  fGun->SetParticleMomentumDirection(
    G4ThreeVector(sinTheta * std::cos(phi), sinTheta * std::sin(phi), -cosTheta));
  fGun->SetParticlePosition(fPlaneCenter
                            + G4ThreeVector((2. * G4UniformRand() - 1.) * fHalfX,
                                            (2. * G4UniformRand() - 1.) * fHalfY, 0.));
  fGun->SetParticleTime(0.);
  fGun->GeneratePrimaryVertex(event);
}
