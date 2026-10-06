#ifndef G4INCLPHASESPACERAUBOLDLYNCH_HH
#define G4INCLPHASESPACERAUBOLDLYNCH_HH

#include "G4INCLThreeVector.hh"
#include "G4INCLParticle.hh"
#include "G4INCLIPhaseSpaceGenerator.hh"
#include "G4INCLInterpolationTable.hh"
#include <vector>

namespace G4INCL {
  /// \brief Generate momenta using the RauboldLynch method
  class PhaseSpaceRauboldLynch : public IPhaseSpaceGenerator {

    public:
      PhaseSpaceRauboldLynch();
      virtual ~PhaseSpaceRauboldLynch();

      /// \brief Dummy copy constructor to silence Coverity warning
      PhaseSpaceRauboldLynch(PhaseSpaceRauboldLynch const &other);

      /// \brief Dummy assignment operator to silence Coverity warning
      PhaseSpaceRauboldLynch &operator=(PhaseSpaceRauboldLynch const &rhs);

      /** \brief Generate momenta according to a uniform, Lorentz-invariant phase-space model
       *
       * This function will assign momenta to the particles in the list that is
       * passed as an argument. The event is generated in the CM frame.
       *
       * \param sqrtS total centre-of-mass energy of the system
       * \param particles list of particles
       */
      void generate(const double sqrtS, ParticleList &particles);

      /// \brief Return the largest generated weight
      double getMaxGeneratedWeight() const;

    private:
      std::vector<double> masses;
      std::vector<double> sumMasses;
      std::vector<double> rnd;
      std::vector<double> invariantMasses;
      std::vector<double> momentaCM;
      size_t nParticles;
      double sqrtS;
      double availableEnergy;
      double maxGeneratedWeight;

      static const size_t nMasslessParticlesTable = 13;
      static const size_t wMaxNE = 30;
      static const double wMaxMasslessX[wMaxNE];
      static const double wMaxMasslessY[wMaxNE];
      static const double wMaxCorrectionX[wMaxNE];
      static const double wMaxCorrectionY[wMaxNE];
      static const double wMaxInterpolationMargin;

      InterpolationTable *wMaxMassless;
      InterpolationTable *wMaxCorrection;

      static const size_t wMaxNP = 20;

      /// \brief Precalculated coefficients: -ln(n)
      double prelog[wMaxNP];

      /// \brief Initialize internal structures (masses and sum of masses)
      void initialize(ParticleList &particles);

      /// \brief Compute the maximum possible weight using a naive algorithm
      double computeMaximumWeightNaive();

      /// \brief Compute the maximum possible weight using parametrizations
      double computeMaximumWeightParam();

      /// \brief Compute the maximum possible weight
      double computeWeight();

      /// \brief Generate an event
      void generateEvent(ParticleList &particles);
  };
}

#endif // G4INCLPHASESPACERAUBOLDLYNCH_HH
