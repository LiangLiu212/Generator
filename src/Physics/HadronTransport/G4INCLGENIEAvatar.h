#ifndef G4INCLGENIEAVATAR_HH_
#define G4INCLGENIEAVATAR_HH_

#include "G4INCLInteractionAvatar.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLAllocationPool.hh"
#include "Physics/HadronTransport/G4INCLGENIEParticleRecord.h"

namespace G4INCL {

  /**
   * Decay avatar
   *
   * The reflection avatar is created when a particle reaches the boundary of the nucleus.
   * At this point it can either be reflected from the boundary or exit the nucleus.
   */
  class GENIEAvatar: public InteractionAvatar {
    public:
      GENIEAvatar(double time, Particle *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord);
      GENIEAvatar(double time, Cluster *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord);
      virtual ~GENIEAvatar();

      IChannel* getChannel();
      void fillFinalState(FinalState *fs);

      virtual void preInteraction();
      virtual void postInteraction(FinalState *fs);

      ParticleList getParticles() const {
	ParticleList theParticleList;
	theParticleList.push_back(particle1);
	return theParticleList;
      }

      std::string dump() const;
    private:
      Particle *particle1;
      Nucleus  *theNucleus;
      Cluster  *cluster;
      std::vector<GENIEParticleRecord> *genie_evtrec;

      INCL_DECLARE_ALLOCATION_POOL(GENIEAvatar)
  };

}

#endif /* G4INCLDECAYAVATAR_HH_ */
