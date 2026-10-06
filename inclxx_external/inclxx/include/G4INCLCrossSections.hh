#ifndef G4INCLCrossSections_hh
#define G4INCLCrossSections_hh 1

#include "G4INCLICrossSections.hh"
#include "G4INCLConfig.hh"

namespace G4INCL {
  namespace CrossSections {
      double elastic(Particle const * const p1, Particle const * const p2);
      double total(Particle const * const p1, Particle const * const p2);

      double NDeltaToNN(Particle const * const p1, Particle const * const p2);
      double NNToNDelta(Particle const * const p1, Particle const * const p2);
      double NNToxPiNN(const int xpi, Particle const * const p1, Particle const * const p2);
      double piNToDelta(Particle const * const p1, Particle const * const p2);
      double piNToxPiN(const int xpi, Particle const * const p1, Particle const * const p2);
      double piNToEtaN(Particle const * const p1, Particle const * const p2);
      double piNToOmegaN(Particle const * const p1, Particle const * const p2);
      double piNToEtaPrimeN(Particle const * const p1, Particle const * const p2);
	    double etaNToPiN(Particle const * const p1, Particle const * const p2);
	    double etaNToPiPiN(Particle const * const p1, Particle const * const p2);
      double etaNToLK(Particle const * const p1, Particle const * const p2);
      double etaNToSK(Particle const * const p1, Particle const * const p2);
      double omegaNToPiN(Particle const * const p1, Particle const * const p2);
      double omegaNToPiPiN(Particle const * const p1, Particle const * const p2);
      double omegaNToLK(Particle const * const p1, Particle const * const p2);
      double omegaNToSK(Particle const * const p1, Particle const * const p2);
      double etaPrimeNToPiN(Particle const * const p1, Particle const * const p2);

      double NNToNNEta(Particle const * const p1, Particle const * const p2);
      double NNToNNEtaExclu(Particle const * const p1, Particle const * const p2);
      double NNToNNEtaxPi(const int xpi, Particle const * const p1, Particle const * const p2);
      double NNToNDeltaEta(Particle const * const p1, Particle const * const p2);
      double NNToNNOmega(Particle const * const p1, Particle const * const p2);
      double NNToNNOmegaExclu(Particle const * const p1, Particle const * const p2);
      double NNToNNOmegaxPi(const int xpi, Particle const * const p1, Particle const * const p2);
      double NNToNDeltaOmega(Particle const * const p1, Particle const * const p2);
      
      /// \brief Strange cross sections
      double NNToNLK(Particle const * const p1, Particle const * const p2);
      double NNToNSK(Particle const * const p1, Particle const * const p2);
      double NNToNLKpi(Particle const * const p1, Particle const * const p2);
      double NNToNSKpi(Particle const * const p1, Particle const * const p2);
      double NNToNLK2pi(Particle const * const p1, Particle const * const p2);
      double NNToNSK2pi(Particle const * const p1, Particle const * const p2);
      double NNToNNKKb(Particle const * const p1, Particle const * const p2);
      double NNToMissingStrangeness(Particle const * const p1, Particle const * const p2);
	    double NDeltaToNLK(Particle const * const p1, Particle const * const p2);
	    double NDeltaToNSK(Particle const * const p1, Particle const * const p2);
	    double NDeltaToDeltaLK(Particle const * const p1, Particle const * const p2);
      double NDeltaToDeltaSK(Particle const * const p1, Particle const * const p2);
	    double NDeltaToNNKKb(Particle const * const p1, Particle const * const p2);
      double NpiToLK(Particle const * const p1, Particle const * const p2);
      double NpiToSK(Particle const * const p1, Particle const * const p2);
	    double p_pimToSzKz(Particle const * const p1, Particle const * const p2);
	    double p_pimToSmKp(Particle const * const p1, Particle const * const p2);
	    double p_pizToSzKp(Particle const * const p1, Particle const * const p2);
      double NpiToLKpi(Particle const * const p1, Particle const * const p2);
      double NpiToSKpi(Particle const * const p1, Particle const * const p2);
      double NpiToLK2pi(Particle const * const p1, Particle const * const p2);
      double NpiToSK2pi(Particle const * const p1, Particle const * const p2);
      double NpiToNKKb(Particle const * const p1, Particle const * const p2);
      double NpiToMissingStrangeness(Particle const * const p1, Particle const * const p2);
      double NLToNS(Particle const * const p1, Particle const * const p2);
      double NSToNL(Particle const * const p1, Particle const * const p2);
      double NSToNS(Particle const * const p1, Particle const * const p2);
      double NKToNK(Particle const * const p1, Particle const * const p2);
      double NKToNKpi(Particle const * const p1, Particle const * const p2);
      double NKToNK2pi(Particle const * const p1, Particle const * const p2);
      double NKbToNKb(Particle const * const p1, Particle const * const p2);
      double NKbToSpi(Particle const * const p1, Particle const * const p2);
      double NKbToLpi(Particle const * const p1, Particle const * const p2);
      double NKbToS2pi(Particle const * const p1, Particle const * const p2);
      double NKbToL2pi(Particle const * const p1, Particle const * const p2);
      double NKbToNKbpi(Particle const * const p1, Particle const * const p2);
      double NKbToNKb2pi(Particle const * const p1, Particle const * const p2);
      double NYelastic(Particle const * const p1, Particle const * const p2);
      double NKbelastic(Particle const * const p1, Particle const * const p2);
      double NKelastic(Particle const * const p1, Particle const * const p2);
      
      /// \brief antiparticle cross sections
      /// \brief Nucleon-AntiNucleon to Baryon-AntiBaryon cross sections
      double NNbarElastic(Particle const* const p1, Particle const* const p2);
      double NNbarCEX(Particle const* const p1, Particle const* const p2);

      double NNbarToLLbar(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-AntiNucleon to Nucleon-AntiNucleon + pions cross sections
      double NNbarToNNbarpi(Particle const* const p1, Particle const* const p2);
      double NNbarToNNbar2pi(Particle const* const p1, Particle const* const p2);
      double NNbarToNNbar3pi(Particle const* const p1, Particle const* const p2);
     
      /// \brief Nucleon-AntiNucleon total annihilation cross sections
      double NNbarToAnnihilation(Particle const* const p1, Particle const* const p2);     
      
      /** \brief Calculate the slope of the NN DDXS.
       *
       * \param energyCM energy in the CM frame, in MeV
       * \param iso total isospin of the system
       *
       * \return the slope of the angular distribution, in (GeV/c)^(-2)
       */
      double calculateNNAngularSlope(double energyCM, int iso);

      /** \brief Compute the "interaction distance".
       *
       * Defined on the basis of the average value of the N-N cross sections at
       * the given kinetic energy.
       *
       * \return the interaction distance
       */
      double interactionDistanceNN(const ParticleSpecies &aSpecies, const double kineticEnergy);
      double interactionDistanceNbarN(const ParticleSpecies &aSpecies, const double kineticEnergy);
      double interactionDistancenbarN(const ParticleSpecies &aSpecies, const double kineticEnergy);
      /** \brief Compute the "interaction distance".
       *
       * Defined on the basis of the average value of the pi-N cross sections at
       * the given kinetic energy.
       *
       * \return the interaction distance
       */
      double interactionDistancePiN(const double projectileKineticEnergy);

      /** \brief Compute the "interaction distance".
       *
       * Defined on the basis of the average value of the K-N cross sections at
       * the given kinetic energy.
       *
       * \return the interaction distance
       */
      double interactionDistanceKN(const double projectileKineticEnergy);

      /** \brief Compute the "interaction distance".
       *
       * Defined on the basis of the average value of the Kbar-N cross sections at
       * the given kinetic energy.
       *
       * \return the interaction distance
       */
      double interactionDistanceKbarN(const double projectileKineticEnergy);

      /** \brief Compute the "interaction distance".
       *
       * Defined on the basis of the average value of the Y-N cross sections at
       * the given kinetic energy.
       *
       * \return the interaction distance
       */
      double interactionDistanceYN(const double projectileKineticEnergy);

      void setCrossSections(ICrossSections *c);

      void deleteCrossSections();

      void initialize(Config const * const theConfig);

  }
}

#endif
