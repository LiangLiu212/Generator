/** \file G4INCLClustering.hh
 * \brief Static class for cluster formation
 *
 * \date 13th July 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLCLUSTERING_HH
#define G4INCLCLUSTERING_HH

#include "G4INCLIClusteringModel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"

namespace G4INCL {

  /// \brief Cluster formation
  namespace Clustering {
    /** \brief Call the clustering algorithm
     *
     * Choose a cluster candidate to be produced. At this point we
     * don't yet decide if it can pass through the Coulomb barrier or
     * not.
     */
    Cluster* getCluster(Nucleus *n, Particle *p);

    /// \brief Determine whether the cluster can escape or not
    bool clusterCanEscape(Nucleus const * const n, Cluster const * const c);

    /// \brief Get the clustering model
    IClusteringModel *getClusteringModel();

    /// \brief Set the clustering model
    void setClusteringModel(IClusteringModel * const model);

    /// \brief Delete the clustering model
    void deleteClusteringModel();

    /// \brief Initialize the clustering model based on the Config object
    void initialize(Config const * const theConfig);

  }
}

#endif
