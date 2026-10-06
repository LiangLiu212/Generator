/** \file G4INCLParticleSampler.hh
 * \brief Class for sampling particles in a nucleus
 *
 * \date 18 July 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLPARTICLESAMPLER_HH_
#define G4INCLPARTICLESAMPLER_HH_

#include "G4INCLNuclearDensity.hh"
#include "G4INCLINuclearPotential.hh"
#include "G4INCLInterpolationTable.hh"

namespace G4INCL {

  class ParticleSampler {

    public:
      /** \brief Constructor.
       *
       * \param A the mass number
       * \param Z the charge number
       */
      ParticleSampler(const int A, const int Z, const int S);

      /// \brief Destructor
      ~ParticleSampler();

      /// \brief Getter for theDensity
      NuclearDensity const *getDensity() const { return theDensity; }

      /// \brief Getter for thePotential
      NuclearPotential::INuclearPotential const *getPotential() const { return thePotential; }

      /// \brief Getter for rpCorrelationCoefficient
      double getRPCorrelationCoefficient(const ParticleType t) const {
        assert(t==Proton || t==Neutron);
        return rpCorrelationCoefficient[t];
      }

      /// \brief Setter for theDensity
      void setDensity(NuclearDensity const * const d);

      /// \brief Setter for thePotential
      void setPotential(NuclearPotential::INuclearPotential const * const p);

      /// \brief Setter for rpCorrelationCoefficient
      void setRPCorrelationCoefficient(const ParticleType t, const double corrCoeff) {
        assert(t==Proton || t==Neutron);
        rpCorrelationCoefficient[t] = corrCoeff;
      }

      ParticleList sampleParticles(ThreeVector const &position);
      void sampleParticlesIntoList(ThreeVector const &position, ParticleList &theList);

    private:

      void updateSampleOneParticleMethods();

      typedef Particle *(ParticleSampler::*ParticleSamplerMethod)(const ParticleType t) const;

      /** \brief Sample a list of particles.
       *
       * This method is a pointer to the method that does the real work for protons.
       */
      ParticleSamplerMethod sampleOneProton;

      /** \brief Sample a list of particles.
       *
       * This method is a pointer to the method that does the real work for neutrons.
       */
      ParticleSamplerMethod sampleOneNeutron;

      /// \brief Sample one particle taking into account the rp-correlation
      Particle *sampleOneParticleWithRPCorrelation(const ParticleType t) const;

      /// \brief Sample one particle not taking into account the rp-correlation
      Particle *sampleOneParticleWithoutRPCorrelation(const ParticleType t) const;

      /// \brief Sample one particle with a fuzzy rp-correlation
      Particle *sampleOneParticleWithFuzzyRPCorrelation(const ParticleType t) const;

      /// \brief Mass number
      const int theA;

      /// \brief Charge number
      const int theZ;

      /// \brief Strangeness number
      const int theS;

      /// \brief Array of pointers to the r-space CDF table
      InterpolationTable const *theRCDFTable[UnknownParticle];

      /// \brief Array of pointers to the p-space CDF table
      InterpolationTable const *thePCDFTable[UnknownParticle];

      /// \brief Pointer to the Cluster's NuclearDensity
      NuclearDensity const *theDensity;

      /// \brief Pointer to the Cluster's NuclearPotential
      NuclearPotential::INuclearPotential const *thePotential;

      /// \brief Correlation coefficients for the r-p correlation
      double rpCorrelationCoefficient[UnknownParticle];
  };

}

#endif // G4INCLPARTICLESAMPLER_HH_
