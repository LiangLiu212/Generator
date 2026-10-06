#ifndef G4INCLIClusteringModel_hh
#define G4INCLIClusteringModel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLCluster.hh"
#include "G4INCLNucleus.hh"

namespace G4INCL {

  /**
   * An abstract interface to cluster formation model(s).
   */
  class IClusteringModel {
  public:
    IClusteringModel() {};
    virtual ~IClusteringModel() {};

    /**
     * Choose a cluster candidate to be produced. At this point we
     * don't yet decide if it can pass through the Coulomb barrier or
     * not.
     */
    virtual Cluster* getCluster(Nucleus*, Particle*) = 0;

    /**
     * Determine whether cluster can escape or not.
     */
    virtual bool clusterCanEscape(Nucleus const * const, Cluster const * const) = 0;
  };

}

#endif
