#ifndef G4INCLDECAYAVATAR_HH_
#define G4INCLDECAYAVATAR_HH_

#include "G4INCLInteractionAvatar.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {

  /**
   * Decay avatar
   *
   * The reflection avatar is created when a particle reaches the boundary of the nucleus.
   * At this point it can either be reflected from the boundary or exit the nucleus.
   */
  class DecayAvatar: public InteractionAvatar {
  public:
    DecayAvatar(G4INCL::Particle *aParticle, double time, G4INCL::Nucleus *aNucleus, bool force=false);
    DecayAvatar(G4INCL::Particle *aParticle, G4INCL::Particle *bParticle, double time, G4INCL::Nucleus *aNucleus, bool force=false);
    virtual ~DecayAvatar();

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
    bool forced;
    ThreeVector const incidentDirection;

    INCL_DECLARE_ALLOCATION_POOL(DecayAvatar)
  };

}

#endif /* G4INCLDECAYAVATAR_HH_ */
