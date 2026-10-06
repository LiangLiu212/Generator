#ifndef G4INCLPauliBlocking__hh
#define G4INCLPauliBlocking__hh 1

#include "G4INCLIPauli.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLConfig.hh"

namespace G4INCL {

  /// \brief Pauli blocking
  namespace Pauli {

    /** \brief Check Pauli blocking.
     *
     * Note: This is a "pure" function: it doesn't retain or modify
     * any state at all and thus only depends on its arguments.
     *
     * \param p list of modified and created particles
     * \param n the nucleus
     */
    bool isBlocked(ParticleList const &p, Nucleus const * const n);

    /** \brief Check CDPP blocking.
     *
     * Note: This is a "pure" function: it doesn't retain or modify
     * any state at all and thus only depends on its arguments.
     *
     * \param p list of created particles
     * \param n the nucleus
     */
    bool isCDPPBlocked(ParticleList const &p, Nucleus const * const n);

    /// \brief Get the Pauli blocker algorithm.
    IPauli * getBlocker();

    /// \brief Get the CDPP blocker algorithm.
    IPauli * getCDPP();

    /// \brief Set the Pauli blocker algorithm.
    void setBlocker(IPauli * const);

    /// \brief Set the CDPP blocker algorithm.
    void setCDPP(IPauli * const);

    /// \brief Delete blockers
    void deleteBlockers();

    /// \brief Initialise blockers according to a Config object
    void initialize(Config const * const aConfig);

  }

}

#endif
