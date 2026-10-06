#ifndef G4INCLParticleEntryAvatar_hh
#define G4INCLParticleEntryAvatar_hh 1

#include "G4INCLIAvatar.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {

  enum EntryType{Default, APAR, ANAR, ADAR}; //D

  class ParticleEntryAvatar: public G4INCL::IAvatar {
  public:
    ParticleEntryAvatar(double, G4INCL::Nucleus*, G4INCL::Particle*, EntryType = Default); //D
    virtual ~ParticleEntryAvatar();
    virtual G4INCL::IChannel* getChannel();
    ParticleList getParticles() const {
      ParticleList theParticleList;
      theParticleList.push_back(theParticle);
      return theParticleList;
    };

    virtual void preInteraction() {};
    virtual void postInteraction(FinalState *);

    std::string dump() const;
   
  EntryType getEntryType() const { return theEType; }; //D
  
  private:
    Nucleus *theNucleus;
    Particle *theParticle;
    EntryType theEType; //D

    INCL_DECLARE_ALLOCATION_POOL(ParticleEntryAvatar)
  };
}

#endif
