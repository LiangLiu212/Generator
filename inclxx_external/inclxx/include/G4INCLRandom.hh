/*
 * G4INCLRandom.hh
 *
 *  \date 7 June 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLRANDOM_HH_
#define G4INCLRANDOM_HH_

#include <iostream>
#include <cmath>
#include <utility>
#include <limits>
#include "G4INCLIRandomGenerator.hh"
#include "G4INCLThreeVector.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLConfig.hh"
#ifdef INCLXX_IN_GEANT4_MODE
  #include "Randomize.hh"
#endif

namespace G4INCL {

  namespace Random {
    /**
     * Set the random number generator implementation to be used globally by INCL.
     *
     * @see G4INCL::IRandomGenerator
     */
    void setGenerator(G4INCL::IRandomGenerator *aGenerator);

    /**
     * Set the seeds of the current generator.
     *
     */
    void setSeeds(const SeedVector &sv);

    /**
     * Get the seeds of the current generator.
     *
     */
    SeedVector getSeeds();

    /**
     * Generate flat distribution of random numbers.
     */
    double shoot();

    /**
     * Return a random number in the ]0,1] interval
     */
    double shoot0();

    /**
     * Return a random number in the [0,1[ interval
     */
    double shoot1();

    /**
     * Return a random integer in the [0,n[ interval
     */
    template<typename T> T shootInteger(T n){
        return static_cast<T>(shoot1() * n);
    }

    /**
     * Generate random numbers using gaussian distribution.
     */
    double gauss(double sigma=1.);

#ifdef INCLXX_IN_GEANT4_MODE
    /**
     * Generate random numbers using gaussian distribution to be used only
     * for correlated pairs
     */
      G4double gaussWithMemory(G4double sigma=1.);
#endif

    /**
     * Generate isotropically-distributed ThreeVectors of given norm.
     */
    ThreeVector normVector(double norm=1.);

    /**
     * Generate ThreeVectors that are uniformly distributed in a sphere of
     * radius rmax.
     */
    ThreeVector sphereVector(double rmax=1.);

    /** \brief Generate Gaussianly-distributed ThreeVectors
     *
     * Generate ThreeVectors that are distributed as a three-dimensional
     * Gaussian of the given sigma.
     */
    ThreeVector gaussVector(double sigma=1.);

    /// \brief Generate pairs of correlated Gaussian random numbers
    std::pair<double,double> correlatedGaussian(const double corrCoeff, const double x0=0., const double sigma=1.);

    /// \brief Generate pairs of correlated uniform random numbers
    std::pair<double,double> correlatedUniform(const double corrCoeff);

    /**
     * Delete the generator
     */
    void deleteGenerator();

    /**
     * Check if the generator is initialized.
     */
    bool isInitialized();

#ifdef INCL_COUNT_RND_CALLS
    /// \brief Return the number of calls to the RNG
    unsigned long long getNumberOfCalls();
#endif

    /// \brief Save the status of the random-number generator
    void saveSeeds();

    /// \brief Get the saved status of the random-number generator
    SeedVector getSavedSeeds();

    /// \brief Initialize generator according to a Config object
#ifdef INCLXX_IN_GEANT4_MODE
    void initialize(Config const * const);
#else
    void initialize(Config const * const theConfig);
#endif

    class Adapter {
      public:
        using result_type = unsigned long long;

        static constexpr result_type min() {
          return std::numeric_limits<Adapter::result_type>::min();
        }

        static constexpr result_type max() {
          return std::numeric_limits<Adapter::result_type>::max();
        }

        result_type operator()() const {

          #ifdef INCLXX_IN_GEANT4_MODE
             return G4RandFlat::shootInt(INT_MAX);
          #else
             return shootInteger(max());
          #endif

        }

    };

    Adapter const &getAdapter();
  }

}

#endif /* G4INCLRANDOM_HH_ */
