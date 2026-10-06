#ifndef G4INCLClusteringModelNone_hh
#define G4INCLClusteringModelNone_hh 1

#include "G4INCLIClusteringModel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLCluster.hh"
#include "G4INCLNucleus.hh"

namespace G4INCL {

  class ClusteringModelNone : public IClusteringModel {
  public:
    ClusteringModelNone() {}
    ~ClusteringModelNone() {}

    virtual Cluster* getCluster(Nucleus*, Particle*) { return 0; }
    virtual bool clusterCanEscape(Nucleus const * const, Cluster const * const) { return false; }

  };

}

#endif
