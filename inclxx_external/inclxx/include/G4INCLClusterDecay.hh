/** \file G4INCLClusterDecay.hh
 * \brief Static class for carrying out cluster decays
 *
 * \date 6th July 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLCLUSTERDECAY_HH
#define G4INCLCLUSTERDECAY_HH

#include "G4INCLCluster.hh"
#include "G4INCLParticleTable.hh"

namespace G4INCL {

  /// \brief Namespace for functions that handle decay of unstable clusters
  namespace ClusterDecay {

      /// \brief True if the cluster is stable.
      bool isStable(Cluster const * const c);

      /** \brief Carries out a cluster decay
       *
       * \param c cluster that should decay
       * \return list of decay products
       */
      ParticleList decay(Cluster * const c);

      enum ClusterDecayType {
        StableCluster,
        NeutronDecay,
        ProtonDecay,
        AlphaDecay,
        TwoProtonDecay,
        TwoNeutronDecay,
        ProtonUnbound,
        NeutronUnbound,
        LambdaUnbound,
        LambdaDecay
      };

      /** \brief Table for cluster decays
       *
       * Definition of "Stable": halflife > 1 ms
       *
       * These table includes decay data for clusters that INCL presently does
       * not produce. It can't hurt.
       *
       * Unphysical nuclides (A<Z) are marked as stable, but should never be
       * produced by INCL. If you find them in the output, something is fishy.
       */
      extern G4ThreadLocal ClusterDecayType clusterDecayMode[ParticleTable::clusterTableSSize][ParticleTable::clusterTableZSize][ParticleTable::clusterTableASize];

  }

}

#endif // G4INCLCLUSTERDECAY
