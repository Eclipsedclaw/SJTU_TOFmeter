//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
//
//-----------------------------------------------------------------------------
//
// GEANT4 Class file 
//
// File name:  TOFmeterEmCaptureCascade
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 22 April 2012 on base of G4MuMinusCaptureCascade
//
//
//-----------------------------------------------------------------------------
//
// Modifications: 
//  - October 2026: renamed from G4EmCaptureCascade for the TOFmeter application, because a
//    class of that name is part of the Geant4 libraries; installed in the mu- capture process
//    by TOFmeterMuonCascade (/tof/physics/muonCascade).
//
//-----------------------------------------------------------------------------

#include "TOFmeterEmCaptureCascade.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh" 
#include "G4MuonMinus.hh"
#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "G4AntiDeuteron.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4NucleiProperties.hh"

#include "G4PrimaryParticle.hh"
#include "G4Track.hh"

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <utility>
#include <vector>

std::ifstream pbarcas_si;
std::ifstream dbarcas_si;
std::ifstream pbarcas_al;
std::ifstream dbarcas_al;
std::ifstream pbarcas_c;
std::ifstream dbarcas_c;

G4double pbar_si[1000][4];
G4double dbar_si[2000][4];
G4double pbar_al[1000][4];
G4double dbar_al[2000][4];
G4double pbar_c[1000][4];
G4double dbar_c[2000][4];


namespace {

struct MuFPTransition {
  G4int n1 = 0;
  G4int l1 = 0;
  G4int n2 = 0;
  G4int l2 = 0;

  // Energy used in ALL THREE transition-rate calculations.
  // This comes from the K-shell table/interpolation, followed by E_n=E_1/n^2.
  G4double deltaERate = 0.0;

  // Energy assigned to the emitted secondary.
  // This comes from the pure hydrogenic muonic-atom levels with reduced mass.
  G4double deltaEOutput = 0.0;

  G4double gammaRad = 0.0;
  G4double gammaK = 0.0;
  G4double gammaL = 0.0;
};

enum class MuFPChannelType { kRadiative, kAugerK, kAugerL };

struct MuFPChannel {
  G4double rate = 0.0;
  G4int n2 = 0;
  G4int l2 = 0;
  G4double emittedEnergy = 0.0;
  MuFPChannelType type = MuFPChannelType::kRadiative;
};

class MuonFirstPrinciplesCascade {
public:
  MuonFirstPrinciplesCascade(G4int Z, G4int A, G4double muMass)
    : fZ(Z), fA(A), fMuonMass(muMass),
      fGammaFlu(1.0e15 * (1.0/second)), fAlpha(0.16), fInitialN(14)
  {
    InitializeKLevelEnergy();
    CalculateHydrogenicLevels();
    CalculateTransitionRates();
  }

  void Simulate(std::vector<G4DynamicParticle*>& secondaries,
                G4double& emittedEnergySum) const
  {
    secondaries.clear();
    emittedEnergySum = 0.0;

    // Modified statistical initial angular-momentum distribution:
    // P_l proportional to (2l+1) exp(alpha*l).
    std::vector<G4double> weights;
    weights.reserve(fInitialN);
    G4double norm = 0.0;
    for (G4int l = 0; l < fInitialN; ++l) {
      const G4double w = (2.0*l + 1.0) * std::exp(fAlpha*l);
      weights.push_back(w);
      norm += w;
    }

    G4double r = G4UniformRand() * norm;
    G4double accum = 0.0;
    G4int currentL = 0;
    for (G4int l = 0; l < fInitialN; ++l) {
      accum += weights[l];
      if (r <= accum) {
        currentL = l;
        break;
      }
    }

    G4int currentN = fInitialN;
    const G4int maxSteps = 100;
    G4int step = 0;

    // Electronic K/L ionisation energies used by the user's transition-rate
    // model.  They are needed here to give the emitted Auger electron its
    // physical kinetic energy T = DeltaE - E_ionisation.
    const G4double rydberg = 1.0967758e7 / m;
    const G4double EK = h_Planck * c_light * rydberg * fZ * fZ;
    const G4double EL = EK / 4.0;

    while (currentN > 1 && step < maxSteps) {
      ++step;
      std::vector<MuFPChannel> channels;

      for (const auto& t : fTransitions) {
        if (t.n1 != currentN || t.l1 != currentL) continue;

        // IMPORTANT:
        //   - the branching RATE was calculated with deltaERate (table/interpolated);
        //   - the displayed/generated secondary ENERGY uses deltaEOutput (hydrogenic).
        if (t.gammaRad > 0.0) {
          channels.push_back({t.gammaRad, t.n2, t.l2, t.deltaEOutput,
                              MuFPChannelType::kRadiative});
        }

        if (t.gammaK > 0.0) {
          const G4double kinetic = t.deltaEOutput - EK;
          if (kinetic > 0.0) {
            channels.push_back({t.gammaK, t.n2, t.l2, kinetic,
                                MuFPChannelType::kAugerK});
          }
        }

        // Gamma_L participates in the same Monte-Carlo competition.
        if (t.gammaL > 0.0) {
          const G4double kinetic = t.deltaEOutput - EL;
          if (kinetic > 0.0) {
            channels.push_back({t.gammaL, t.n2, t.l2, kinetic,
                                MuFPChannelType::kAugerL});
          }
        }
      }

      // With only Delta-l = +/-1 channels, a metastable 2s state can have no
      // outgoing transition.  This reproduces the current first-principles
      // model rather than silently inserting a standard-Geant4 transition.
      if (channels.empty()) break;

      G4double totalRate = 0.0;
      for (const auto& c : channels) totalRate += c.rate;
      if (totalRate <= 0.0) break;

      const G4double pick = G4UniformRand() * totalRate;
      G4double cumulative = 0.0;
      const MuFPChannel* selected = &channels.back();
      for (const auto& c : channels) {
        cumulative += c.rate;
        if (pick <= cumulative) {
          selected = &c;
          break;
        }
      }

      if (selected->emittedEnergy > 0.0) {
        G4ParticleDefinition* definition = nullptr;

        if (selected->type == MuFPChannelType::kRadiative) {
          definition = G4Gamma::Gamma();
        } else {
          definition = G4Electron::Electron();
        }

        secondaries.push_back(new G4DynamicParticle(
          definition, G4RandomDirection(), selected->emittedEnergy));
        emittedEnergySum += selected->emittedEnergy;
      }

      currentN = selected->n2;
      currentL = selected->l2;
    }
  }

private:
  static G4double Factorial(G4int n)
  {
    if (n < 0) return 0.0;
    G4double value = 1.0;
    for (G4int i = 2; i <= n; ++i) value *= i;
    return value;
  }

  static G4double EvalGenLaguerre(G4int n, G4int alpha, G4double x)
  {
    if (n < 0) return 0.0;
    if (n == 0) return 1.0;

    // No artificial x>30 cutoff.  Use long double internally to improve
    // stability of the three-term generalized-Laguerre recurrence.
    const long double xx = static_cast<long double>(x);
    const long double aa = static_cast<long double>(alpha);

    if (n == 1) {
      return static_cast<G4double>(1.0L + aa - xx);
    }

    long double previous = 1.0L;
    long double current = 1.0L + aa - xx;

    for (G4int k = 1; k < n; ++k) {
      const long double kk = static_cast<long double>(k);
      const long double next =
        ((2.0L*kk + 1.0L + aa - xx)*current
         - (kk + aa)*previous)
        / (kk + 1.0L);

      if (!std::isfinite(next)) return 0.0;

      previous = current;
      current = next;
    }

    return static_cast<G4double>(current);
  }

  G4double HydrogenicRadialWavefunction(G4int n, G4int l,
                                        G4int Z, G4double r,
                                        G4double massRatio) const
  {
    const G4double aEff = 1.0 / massRatio;
    const G4double term1 = (Z / (1.0*n*n)) / aEff;
    const G4double fact1 = Factorial(n - l - 1);
    const G4double fact2 = Factorial(n + l);
    if (fact2 <= 0.0) return 0.0;

    const G4double norm = std::sqrt(term1 * fact1 / fact2);
    const G4double rho = 2.0 * Z * r / (n * aEff);

    // No artificial rho>60 cutoff.  The physical exp(-rho/2) factor and
    // the finite integration range control the tail.
    const G4double laguerre = EvalGenLaguerre(n-l-1, 2*l+1, rho);
    G4double value = norm * std::exp(-rho/2.0)
                   * std::pow(rho, l+1) * laguerre;
    if (!std::isfinite(value)) return 0.0;

    const G4double threshold = 1.0e-80;
    if (std::abs(value) < threshold) {
      value = (value >= 0.0) ? threshold : -threshold;
    }
    return value;
  }

  G4double CalculateRadialMatrixElement(G4int n1, G4int l1,
                                        G4int n2, G4int l2)
  {
    if (std::abs(l1-l2) != 1) return 0.0;

    const std::uint64_t key =
      (static_cast<std::uint64_t>(n1) << 48) |
      (static_cast<std::uint64_t>(l1) << 40) |
      (static_cast<std::uint64_t>(n2) << 32) |
      (static_cast<std::uint64_t>(l2) << 24);

    auto found = fRadialCache.find(key);
    if (found != fRadialCache.end()) return found->second;

    const G4double massA = G4NucleiProperties::GetNuclearMass(fA, fZ);
    const G4double reducedMass = fMuonMass * massA / (fMuonMass + massA);
    const G4double massRatio = reducedMass / electron_mass_c2;
    const G4double aEff = 1.0 / massRatio;
    const G4int nMax = std::max(n1, n2);

    const G4double rMin = 1.0e-12;
    const G4double rMax = 5.0 * nMax*nMax * aEff / fZ;
    if (rMax <= rMin) return 0.0;

    // Keep the same log-grid trapezoidal integration strategy as the user's
    // CascadeCalculator implementation.
    const G4int nPoints = 20000;
    const G4double logMin = std::log(rMin);
    const G4double logMax = std::log(rMax);
    const G4double dlog = (logMax-logMin)/(nPoints-1);

    G4double integral = 0.0;
    G4double rPrev = std::exp(logMin);
    G4double fPrev =
      HydrogenicRadialWavefunction(n1, l1, fZ, rPrev, massRatio) *
      HydrogenicRadialWavefunction(n2, l2, fZ, rPrev, massRatio) * rPrev;

    for (G4int i = 1; i < nPoints; ++i) {
      const G4double rr = std::exp(logMin + i*dlog);
      const G4double ff =
        HydrogenicRadialWavefunction(n1, l1, fZ, rr, massRatio) *
        HydrogenicRadialWavefunction(n2, l2, fZ, rr, massRatio) * rr;
      integral += 0.5 * (fPrev + ff) * (rr-rPrev);
      rPrev = rr;
      fPrev = ff;
    }

    if (!std::isfinite(integral)) integral = 0.0;

    // IMPORTANT NORMALISATION:
    // The log-grid integration above is carried out with r measured in
    // Bohr-radius units and returns a quantity proportional to a_mu/a0,
    // where a_mu/a0 = 1/(lambda_red*Z).  The cascade-rate formulae of
    // Eisenberg/Kessler and Aramaki use the DIMENSIONLESS dipole integral
    //
    //     I = integral dx x^3 Rbar_f(x) Rbar_i(x),
    //     x = r / a_mu,  a_mu = a0/(lambda_red Z).
    //
    // Therefore convert the raw numerical integral to that dimensionless I:
    const G4double radialDimensionless = integral * massRatio * fZ;

    if (!std::isfinite(radialDimensionless)) {
      fRadialCache[key] = 0.0;
      return 0.0;
    }

    fRadialCache[key] = radialDimensionless;
    return radialDimensionless;
  }

  void InitializeKLevelEnergy()
  {
    static const G4int nlevels = 28;
    static const G4int listK[nlevels] = {
      1,2,4,6,8,11,14,17,18,21,24,26,29,32,
      38,40,41,44,49,53,55,60,65,70,75,81,85,92
    };
    static const G4double listKEnergy[nlevels] = {
      0.00275,0.011,0.043,0.098,0.173,0.326,0.524,
      0.765,0.853,1.146,1.472,1.708,2.081,2.475,
      3.323,3.627,3.779,4.237,5.016,5.647,5.966,
      6.793,7.602,8.421,9.249,10.222,10.923,11.984
    };

    fKLevelEnergy.fill(0.0);
    fKLevelEnergy[1] = listKEnergy[0] * MeV;
    G4int idx = 0;
    for (G4int i = 1; i < nlevels; ++i) {
      const G4int z1 = listK[idx];
      const G4int z2 = listK[i];
      if (z1+1 < z2) {
        const G4double dz = static_cast<G4double>(z2-z1);
        const G4double y1 = listKEnergy[idx]/static_cast<G4double>(z1*z1);
        const G4double y2 = listKEnergy[i]/static_cast<G4double>(z2*z2);
        for (G4int z = z1+1; z < z2; ++z) {
          fKLevelEnergy[z] =
            (y1 + (y2-y1)*(z-z1)/dz) * z*z * MeV;
        }
      }
      fKLevelEnergy[z2] = listKEnergy[i] * MeV;
      idx = i;
    }
  }

  void CalculateHydrogenicLevels()
  {
    // Two deliberately separate energy scales are kept.
    //
    // fLevelEnergy:
    //   Used ONLY for calculating Gamma_rad, Gamma_K and Gamma_L.
    //   E1 comes from the tabulated/interpolated muonic K-shell energy and
    //   higher levels follow E_n = E1/n^2.
    //
    // fHydrogenicLevelEnergy:
    //   Used ONLY for the energy assigned to generated gamma/Auger secondaries.
    //   It is the pure hydrogenic muonic-atom formula with nuclear reduced mass.
    fLevelEnergy.assign(fInitialN+2, 0.0);
    fHydrogenicLevelEnergy.assign(fInitialN+2, 0.0);

    const G4double massA = G4NucleiProperties::GetNuclearMass(fA, fZ);
    const G4double reducedMass =
      fMuonMass * massA / (fMuonMass + massA);

    // Rate-energy scale: K-shell table/interpolation.
    const G4double tableE1 = fKLevelEnergy[std::min(fZ,92)];

    // Output-energy scale: pure hydrogenic muonic atom.
    const G4double hydrogenicE1 =
      13.6 * eV * static_cast<G4double>(fZ*fZ)
      * reducedMass / electron_mass_c2;

    fLevelEnergy[0] = tableE1;
    fHydrogenicLevelEnergy[0] = hydrogenicE1;

    for (G4int i = 1; i <= fInitialN; ++i) {
      const G4int n = i+1;
      const G4double n2 = static_cast<G4double>(n*n);
      fLevelEnergy[i] = tableE1 / n2;
      fHydrogenicLevelEnergy[i] = hydrogenicE1 / n2;
    }
  }

  void CalculateTransitionRates()
  {
    fTransitions.clear();

    const G4double rydberg = 1.0967758e7 / m;

    // Use the muon-nucleus reduced mass consistently in the rate formulae.
    // This is the M* (or lambda=M*/m_e) appearing in the exotic-atom
    // radiative/Auger expressions and is the same mass scale used by the
    // hydrogenic radial wavefunctions above.
    const G4double massA = G4NucleiProperties::GetNuclearMass(fA, fZ);
    const G4double reducedMass = fMuonMass * massA / (fMuonMass + massA);
    const G4double massRatio = reducedMass / electron_mass_c2;
    const G4double a0OverM = Bohr_radius / massRatio;
    const G4double EK = h_Planck * c_light * rydberg * fZ*fZ;
    const G4double EL = EK/4.0;
    const G4double Zstar = fZ - 1.0;

    const G4double prefK =
      32.0*pi*fine_structure_const*c_light/Bohr_radius /
      (massRatio*massRatio) * (Zstar*Zstar)/(fZ*fZ);
    const G4double prefL =
      16.0*pi*fine_structure_const*c_light/Bohr_radius /
      (massRatio*massRatio) * (Zstar*Zstar)/(fZ*fZ);

    for (G4int n1 = 2; n1 <= fInitialN; ++n1) {
      for (G4int l1 = 0; l1 < n1; ++l1) {
        for (G4int n2 = 1; n2 < n1; ++n2) {
          const G4int candidates[2] = {l1-1, l1+1};
          for (G4int idx = 0; idx < 2; ++idx) {
            const G4int l2 = candidates[idx];
            if (l2 < 0 || l2 >= n2) continue;

            // deltaERate: table/interpolated energy difference used by ALL rates.
            const G4double deltaERate =
              fLevelEnergy[n2-1] - fLevelEnergy[n1-1];

            // deltaEOutput: pure hydrogenic energy difference used only for
            // the energy of the generated secondary particle.
            const G4double deltaEOutput =
              fHydrogenicLevelEnergy[n2-1]
              - fHydrogenicLevelEnergy[n1-1];

            if (deltaERate <= 0.0 || deltaEOutput <= 0.0) continue;

            // CalculateRadialMatrixElement now returns the dimensionless
            // literature integral I, not the unnormalised log-grid integral.
            const G4double radialI =
              CalculateRadialMatrixElement(n1,l1,n2,l2);
            if (radialI == 0.0) continue;

            const G4int lMax = std::max(l1,l2);
            const G4double radialSq = radialI*radialI;

            // Radiative transition rate.
            //
            // IMPORTANT: Gamma_rad now uses the SAME table/interpolated
            // deltaERate as Gamma_K and Gamma_L.  This is the explicit DeltaE^3
            // form that was already present as the alternative/commented formula
            // in the user's CascadeCalculator.cc:
            //
            // Gamma_rad = [4 e^2/(4 pi eps0)] * (a0/mu)^2 / Z^2
            //             * (DeltaE/hbar)^3
            //             * max(l1,l2)/(2l1+1)
            //             / [3 hbar c^3] * |I|^2 .
            const G4double radTerm1 =
              4.0 * eplus*eplus / (4.0*pi*epsilon0);
            const G4double radTerm2 = a0OverM*a0OverM;
            const G4double radTerm3 =
              std::pow(static_cast<G4double>(fZ),2);
            const G4double radTerm4 =
              std::pow(deltaERate/hbar_Planck,3);
            const G4double radTerm5 =
              static_cast<G4double>(lMax)/(2.0*l1+1.0);
            const G4double radTerm6 =
              3.0*hbar_Planck*std::pow(c_light,3);

            const G4double gammaRad =
              (radTerm1 * radTerm2 / radTerm3
               * radTerm4 * radTerm5 / radTerm6) * radialSq;

            G4double gammaK = 0.0;
            G4double gammaL = 0.0;

            const G4double TK = deltaERate-EK;
            if (TK > 0.0 && Zstar > 0.0) {
              const G4double y = Zstar*fine_structure_const /
                std::sqrt(std::pow(TK/electron_mass_c2,2)
                          + 2.0*TK/electron_mass_c2);
              const G4double a1 = static_cast<G4double>(lMax)/(3.0*(2.0*l1+1.0));
              const G4double a2 = y*y/(1.0+y*y);
              const G4double expFactor = std::exp(y*(4.0*std::atan(y)-pi));
              const G4double sinhFactor = std::sinh(pi*y);
              if (sinhFactor > 0.0 && std::isfinite(sinhFactor)) {
                const G4double raw = prefK*a1*a2*expFactor/sinhFactor*radialSq;
                if (raw > 0.0 && fGammaFlu > 0.0) {
                  gammaK = 1.0/(1.0/raw + 1.0/fGammaFlu);
                }
              }
            }

            const G4double TL = deltaERate-EL;
            if (TL > 0.0 && Zstar > 0.0) {
              const G4double y = Zstar*fine_structure_const /
                std::sqrt(std::pow(TL/electron_mass_c2,2)
                          + 2.0*TL/electron_mass_c2);
              const G4double a1 = static_cast<G4double>(lMax)/(3.0*(2.0*l1+1.0));
              const G4double a2 = y*y*(4.0+5.0*y*y)*(4.0+3.0*y*y)
                                / std::pow(4.0+y*y,3);
              // Eisenberg/Kessler L-shell form used in the user's code.
              const G4double expFactor =
                std::exp(y*(4.0*std::atan(y/2.0)-pi));
              const G4double sinhFactor = std::sinh(pi*y);
              if (sinhFactor > 0.0 && std::isfinite(sinhFactor)) {
                const G4double raw = prefL*a1*a2*expFactor/sinhFactor*radialSq;
                if (raw > 0.0) gammaL = raw;
              }
            }

            if (gammaRad > 0.0 || gammaK > 0.0 || gammaL > 0.0) {
              fTransitions.push_back({n1,l1,n2,l2,
                                      deltaERate,deltaEOutput,
                                      gammaRad,gammaK,gammaL});
            }
          }
        }
      }
    }
  }

  G4int fZ = 0;
  G4int fA = 0;
  G4double fMuonMass = 0.0;
  G4double fGammaFlu = 0.0;
  G4double fAlpha = 0.16;
  G4int fInitialN = 14;
  std::array<G4double,93> fKLevelEnergy{};
  std::vector<G4double> fLevelEnergy;
  std::vector<G4double> fHydrogenicLevelEnergy;
  std::vector<MuFPTransition> fTransitions;
  std::map<std::uint64_t,G4double> fRadialCache;
};

MuonFirstPrinciplesCascade& GetMuonFirstPrinciplesCascade(G4int Z, G4int A,
                                                           G4double muMass)
{
  // One cache per Geant4 worker thread.  Expensive radial integrals and rates
  // are calculated once for each target isotope, then reused event-by-event.
  using Key = std::pair<G4int,G4int>;
  thread_local std::map<Key,std::unique_ptr<MuonFirstPrinciplesCascade>> cache;

  const Key key(Z,A);
  auto it = cache.find(key);
  if (it == cache.end()) {
    std::unique_ptr<MuonFirstPrinciplesCascade> model(
      new MuonFirstPrinciplesCascade(Z,A,muMass));
    it = cache.emplace(key,std::move(model)).first;
  }
  return *(it->second);
}

} // namespace


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TOFmeterEmCaptureCascade::TOFmeterEmCaptureCascade()
  : G4HadronicInteraction("TOFmeterEmCaptureCascade")
{ 
  theElectron = G4Electron::Electron();
  theGamma = G4Gamma::Gamma();
  fMuMass = G4MuonMinus::MuonMinus()->GetPDGMass();
  fPbarMass = G4AntiProton::AntiProton()->GetPDGMass();
  fNbarMass = G4AntiNeutron::AntiNeutron()->GetPDGMass();
  fDbarMass = G4AntiDeuteron::AntiDeuteron()->GetPDGMass();

  fAntiParticleMass = 0.0;
  fAntiParticleType = 0;
  fTime = 0.0;

  // Calculate the Energy of K Mesoatom Level for this Element using
  // the Energy of Hydrogen Atom taken into account finite size of the
  // nucleus 
  static const G4int nlevels = 28;
  static const G4int listK[nlevels] = {
      1, 2,  4,  6,  8, 11, 14, 17, 18, 21, 24,
     26, 29, 32, 38, 40, 41, 44, 49, 53, 55,
     60, 65, 70, 75, 81, 85, 92};
  static const G4double listKEnergy[nlevels] = {
     0.00275, 0.011, 0.043, 0.098, 0.173, 0.326,
     0.524, 0.765, 0.853, 1.146, 1.472,
     1.708, 2.081, 2.475, 3.323, 3.627, 
     3.779, 4.237, 5.016, 5.647, 5.966,
     6.793, 7.602, 8.421, 9.249, 10.222,
    10.923,11.984};

  fKLevelEnergy[0] = 0.0;
  fKLevelEnergy[1] = listKEnergy[0];
  G4int idx = 1;
  for(G4int i=1; i<nlevels; ++i) {
    G4int z1 = listK[idx];
    G4int z2 = listK[i];
    if(z1+1 < z2) {
      G4double dz = G4double(z2 - z1);
      G4double y1 = listKEnergy[idx]/G4double(z1*z1);
      G4double y2 = listKEnergy[i]/G4double(z2*z2);
      for(G4int z=z1+1; z<z2; ++z) {
        fKLevelEnergy[z] = (y1 + (y2 - y1)*(z - z1)/dz)*z*z;
      }
    }
    fKLevelEnergy[z2] = listKEnergy[i];
    idx = i;  
  }
  //initial energy for all level, dbar=59, pbar=42, mu=14
  for(G4int i = 0; i<59; ++i) { fLevelEnergy[i] = 0.0; }

  // table for Silicon
  pbarcas_si.open("/home/mjxiao/G4data/cascade_pbarSi.data");
  for (G4int ii=0;ii<32;++ii)
  {
      for (G4int jj=0;jj<4;++jj)
      {
          pbarcas_si>>pbar_si[ii][jj];
      }
  }
  dbarcas_si.open("/home/mjxiao/G4data/cascade_dbarSi.data");
  for (G4int ii=0;ii<55;++ii)
  {
      for (G4int jj=0;jj<4;++jj)
      {
          dbarcas_si>>dbar_si[ii][jj];
      }
  }

  // tables for Al
  pbarcas_al.open("/home/mjxiao/G4data/cascade_pbarAl.data");
  for (G4int ii=0;ii<31;++ii)
  {
      for (G4int jj=0;jj<4;++jj)
      {
          pbarcas_al>>pbar_al[ii][jj];
      }
  }
  dbarcas_al.open("/home/mjxiao/G4data/cascade_dbarAl.data");
  for (G4int ii=0;ii<52;++ii)
  {
      for (G4int jj=0;jj<4;++jj)
      {
          dbarcas_al>>dbar_al[ii][jj];
      }
  }

  // tables for Carbon
  pbarcas_c.open("/home/mjxiao/G4data/cascade_pbarC.data");
  for (G4int ii=0;ii<27;++ii)
  {
      for (G4int jj=0;jj<4;++jj)
      {
          pbarcas_c>>pbar_c[ii][jj];
      }
  }
  dbarcas_c.open("/home/mjxiao/G4data/cascade_dbarC.data");
  for (G4int ii=0;ii<48;++ii)
  {
      for (G4int jj=0;jj<4;++jj)
      {
          dbarcas_c>>dbar_c[ii][jj];
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TOFmeterEmCaptureCascade::~TOFmeterEmCaptureCascade()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4HadFinalState* 
TOFmeterEmCaptureCascade::ApplyYourself(const G4HadProjectile& projectile, 
				  G4Nucleus& targetNucleus)
{
  result.Clear();
  result.SetStatusChange(isAlive);
  fTime = projectile.GetGlobalTime();
  fAntiParticleMass = projectile.GetDefinition()->GetPDGMass();
  fAntiParticleType = projectile.GetDefinition()->GetPDGEncoding();
  G4String fAntiParticleName = projectile.GetDefinition()->GetParticleName();
//  G4cout<<"Name: "<<fAntiParticleName<<" "<<"PDGCode: "<<fAntiParticleType<<" "<<"Mass: "<<fAntiParticleMass<<"MeV"<<G4endl;

  G4int Z = targetNucleus.GetZ_asInt(); 
  G4int A = targetNucleus.GetA_asInt(); 

  // ------------------------------------------------------------------
  // mu-: user's first-principles n=14 -> lower-n cascade
  // ------------------------------------------------------------------
  // The legacy/standard Geant4 muon probability block is NOT used.
  // For mu- we build the complete (n,l) transition map from n=14 using
  // calculated Gamma_rad, Gamma_K and Gamma_L, draw the initial l from
  // P(l) ~ (2l+1)exp(0.16*l), and select every step with
  // P_i = Gamma_i / sum_j Gamma_j.
  if (fAntiParticleType == 13) {
    MuonFirstPrinciplesCascade& muCascade =
      GetMuonFirstPrinciplesCascade(Z, A, fAntiParticleMass);

    std::vector<G4DynamicParticle*> muSecondaries;
    G4double muCascadeEnergy = 0.0;
    muCascade.Simulate(muSecondaries, muCascadeEnergy);

    for (auto* dp : muSecondaries) {
      if (!dp) continue;
      G4HadSecondary hs(dp);
      hs.SetTime(fTime);
      result.AddSecondary(hs);
    }

    // Keep the same energy-bookkeeping convention as the user's existing
    // CascadeCalculator / CaptureCascadeProcess implementation.
    result.SetLocalEnergyDeposit(muCascadeEnergy);
    return &result;
  }
  G4double massA = G4NucleiProperties::GetNuclearMass(A, Z);
  G4double mass = fAntiParticleMass * massA / (fAntiParticleMass + massA) ;
  G4double e = 13.6 * eV * (Z * Z) * mass/ electron_mass_c2;

  // precise corrections of energy only for K-shell
  fLevelEnergy[0] = fKLevelEnergy[std::min(Z, 92)];
  for(G4int i=1; i<60; ++i) {
    fLevelEnergy[i] = e/(G4double)((i+1)*(i+1));
  }
  fLevelEnergy[0] = fLevelEnergy[1];
  
  G4int nLevel = 0;
  if (fAntiParticleType==-2212) nLevel = 42;
  else if (fAntiParticleType==-1000010020) nLevel = 59;
  else nLevel = 14;

  // Advisor model: 42-th level for pbar, 60-th level for dbar.
  // The n=14 muon cascade has already been completed and returned above.
  G4double edep=0.;
  if (fAntiParticleType==-2212) edep = fLevelEnergy[41];
  else if (fAntiParticleType==-1000010020) edep = fLevelEnergy[58];
  else edep = fLevelEnergy[13];

  AddNewParticle(theElectron,edep);
  G4double deltaE = 0.;
  G4int iLevel = 0;
  
  if (fAntiParticleType==-2212)  // antiproton case, M.Xiao, May.31-2019
  {
    do {
         G4double var_pbar = 1.0*G4UniformRand();
	 if (Z==14) // use pbar_Si table
	 {
	   for (G4int k=1;k<nLevel;k++)
	   {
	     iLevel = k;
	     deltaE = fLevelEnergy[iLevel-1] -  fLevelEnergy[nLevel-1]; 
	     for (G4int kk=0;kk<32; kk++)
	     {
	       if (nLevel==pbar_si[kk][0]&&iLevel==pbar_si[kk][1]&&var_pbar<pbar_si[kk][2]) AddNewParticle(theGamma, deltaE);
	       else continue;
	    }
	  }
	}
	if (Z==13) //use pbar-Al table
	{
	  for (G4int k=1;k<nLevel;k++)
	  {
	    iLevel = k;
	    deltaE = fLevelEnergy[iLevel-1] -  fLevelEnergy[nLevel-1]; 
	    for (G4int kk=0;kk<31; kk++)
	    {
	      if (nLevel==pbar_al[kk][0]&&iLevel==pbar_al[kk][1]&&var_pbar<pbar_al[kk][2]) AddNewParticle(theGamma, deltaE);
	      else continue;
	    }
	  }
	}
	if (Z==6) //use pbar-Carbon table, added by M.Xiao, Jul.1-2019
	{
	  for (G4int k=1;k<nLevel;k++)
	  {
	    iLevel = k;
	    deltaE = fLevelEnergy[iLevel-1] -  fLevelEnergy[nLevel-1]; 
	    for (G4int kk=0;kk<27; kk++)
	    {
	      if (nLevel==pbar_c[kk][0]&&iLevel==pbar_c[kk][1]&&var_pbar<pbar_c[kk][2]) AddNewParticle(theGamma, deltaE);
	      else continue;
	    }
	  }
	}
		  
	--nLevel;
	edep += deltaE;
       } while( nLevel > 0 );
  }
  
  if (fAntiParticleType==-1000010020) // antideuteron case, M.Xiao, May.31-2019
  {
    do {
         G4double var_dbar = 1.0*G4UniformRand();
	 if (Z==14) // use dbar-Si table
	 {
	   for (G4int k=1;k<nLevel;k++)
	   {
	     iLevel = k;
	     deltaE = fLevelEnergy[iLevel-1] - fLevelEnergy[nLevel-1];
	     for (G4int kk=0; kk<55; kk++)
	     {
	       if (nLevel==dbar_si[kk][0]&&iLevel==dbar_si[kk][1]&&var_dbar<dbar_si[kk][2]) AddNewParticle(theGamma, deltaE);
	       else continue;
	     }	
	   }
	 }
	if (Z==13) // use dbr-Al table
	{
	  for (G4int k=1;k<nLevel;k++)
	  {
	    iLevel = k;
	    deltaE = fLevelEnergy[iLevel-1] - fLevelEnergy[nLevel-1];
	    for (G4int kk=0; kk<52; kk++)
	    {
	      if (nLevel==dbar_al[kk][0]&&iLevel==dbar_al[kk][1]&&var_dbar<dbar_al[kk][2]) AddNewParticle(theGamma, deltaE);
	      else continue;
	    }
	  }
	}
	 if (Z==6) // use dbar-Carbon table, added by M.Xiao, Jul.1-2019
	 {
	   for (G4int k=1;k<nLevel;k++)
	   {
	     iLevel = k;
	     deltaE = fLevelEnergy[iLevel-1] - fLevelEnergy[nLevel-1];
	     for (G4int kk=0; kk<48; kk++)
	     {
	       if (nLevel==dbar_c[kk][0]&&iLevel==dbar_c[kk][1]&&var_dbar<dbar_c[kk][2]) AddNewParticle(theGamma, deltaE);
	       else continue;
	     }	
	   }
	 }

	--nLevel;
	edep += deltaE; 
       } while(nLevel > 0);
  }

  // mu- does not reach this point: it is handled above by the
  // first-principles (n,l)-resolved cascade and returns immediately.
  // The old Geant4 simplified muon-cascade block has been removed.
  
  if (fAntiParticleType!=-2212&&fAntiParticleType!=-1000010020&&fAntiParticleType!=13) // for other antiparticles, no cascade x-rays, M.Xiao, May.31-2019
  {
    do {
         deltaE = 0.;
	 AddNewParticle(theGamma, deltaE);
		
	--nLevel;
      } while( nLevel > 0 );
  }

//  result.Clear(); // clear the above results
  result.SetLocalEnergyDeposit(edep);
  return &result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TOFmeterEmCaptureCascade::ModelDescription(std::ostream& outFile) const
{
  outFile << "Electromagnetic cascade model for stopped negative particles.\n"
          << "mu-: first-principles n=14 (n,l)-resolved cascade; all three "
          << "rates use table/interpolated DeltaE, while emitted energies use "
          << "pure hydrogenic levels.\n"
          << "anti-proton / anti-deuteron: advisor tabulated cascade model.\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....