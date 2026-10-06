/** \file G4INCLCrossSectionsTruncatedMultiPions.hh
 * \brief Truncated multipion cross sections
 *
 * \date 10th December 2014
 * \author Davide Mancusi
 */

#ifndef G4INCLCROSSSECTIONSTRUNCATEDMULTIPIONS_HH
#define G4INCLCROSSSECTIONSTRUNCATEDMULTIPIONS_HH

#include "G4INCLCrossSectionsMultiPions.hh"
#include <limits>

namespace G4INCL {
  /// \brief Truncated multipion cross sections

  class CrossSectionsTruncatedMultiPions : public CrossSectionsMultiPions {
    public:
      CrossSectionsTruncatedMultiPions(const int nPi=std::numeric_limits<int>::max());

      /** \brief Elastic particle-particle cross section
       *
       * The elastic pi-N cross section is calculated by difference:
       *   elastic = total - inelastic(pion production) - delta,
       * but we need to use the unmodified piNToDelta cross section, or else we
       * will violate the cross-section sum.
       *
       * This modification is only necessary for piN collisions.
       */
      virtual double elastic(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for Delta production - piN Channel
      virtual double piNToDelta(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for X pion production - piN Channel
      virtual double piNToxPiN(const int xpi, Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for X pion production - NN Channel
      virtual double NNToxPiNN(const int xpi, Particle const * const p1, Particle const * const p2);

    protected:

      /// \brief Maximum number of pions produced in TruncatedMultiPion collisions
      const int nMaxPi;

  };
}

#endif
