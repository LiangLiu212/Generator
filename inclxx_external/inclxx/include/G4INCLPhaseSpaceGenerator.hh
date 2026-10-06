#ifndef G4INCLPHASESPACEGENERATOR_HH
#define G4INCLPHASESPACEGENERATOR_HH 1

#include "G4INCLIPhaseSpaceGenerator.hh"
#include "G4INCLConfig.hh"

namespace G4INCL {
  namespace PhaseSpaceGenerator {
    /// \brief Generate an event in the CM system
    void generate(const double sqrtS, ParticleList &particles);

    /** \brief Generate a biased event in the CM system
     *
     * This method first generates a "flat" event by calling generate(). The
     * particles are subsequently rotated in such a way that one of them
     * (identified by the parameter index) is biased towards the collisionAxis
     * with an exponential distribution of the form
     * \f[
     * \exp(B\cdot t)
     * \f]
     * where \f$t\f$ is the usual Mandelstam variable. The incoming momentum
     * is taken to be the momentum of particles[index] at the moment of the
     * call.
     *
     * \param sqrtS total energy in the centre of mass, in MeV
     * \param particles list of particles for which the event will be generated
     *        (modified on exit)
     * \param index index of the particle to be biased; all the other particles
     *        will follow
     * \param slope slope \f$B\f$ of the angular distribution: \f$\exp(Bt)\f$,
     *        in (GeV/c)^(-2)
     */
    void generateBiased(const double sqrtS, ParticleList &particles, const size_t index, const double slope);

    void setPhaseSpaceGenerator(IPhaseSpaceGenerator *g);

    IPhaseSpaceGenerator *getPhaseSpaceGenerator();

    void deletePhaseSpaceGenerator();

    void initialize(Config const * const theConfig);
  }
}

#endif
