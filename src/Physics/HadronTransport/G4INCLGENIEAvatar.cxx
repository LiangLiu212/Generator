#include "G4INCLGENIEAvatar.h"

#include "G4INCLGENIEQELChannel.h"
#include "G4INCLGENIEMECChannel.h"
#include "G4INCLGENIEDISChannel.h"
#include "G4INCLGENIERESChannel.h"
#include "G4INCLPauliBlocking.hh"
#include <sstream>
#include <string>
#include <cassert>

namespace G4INCL {

  GENIEAvatar::GENIEAvatar(double time, Particle *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord)
    : InteractionAvatar(time, n, p),
    particle1(p), theNucleus(n),
    genie_evtrec(eventRecord)
  {
    //    setType(GENIEAvatarType);
  }


  GENIEAvatar::GENIEAvatar(double time, Cluster *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord)
    : InteractionAvatar(time, n, nullptr),
    particle1(p), theNucleus(n), cluster(p),
    genie_evtrec(eventRecord)
  {
    //    setType(GENIEAvatarType);
  }



  GENIEAvatar::~GENIEAvatar() {

  }

  G4INCL::IChannel* GENIEAvatar::getChannel() {

    // TODO: using genie event record to fill final states of INCL
    // the initial states is hit nucleons, (for MEC channel, initial states
    // is NN cluster)
    // the final states should be hadron in nucleus
    switch((*genie_evtrec)[0].ScatteringType()){
      case 1: return new GENIEQELChannel(particle1, theNucleus, genie_evtrec);
      case 4: return new GENIERESChannel(particle1, theNucleus, genie_evtrec);
      case 3: return new GENIEDISChannel(particle1, theNucleus, genie_evtrec);
      case 10: return new GENIEMECChannel(cluster, theNucleus, genie_evtrec);
      default: return NULL;
    }

  }

  void GENIEAvatar::preInteraction() {
    // InteractionAvatar::preInteraction();
    return;
  }

  void GENIEAvatar::postInteraction(FinalState *fs) {
    // Make sure we have at least two particles in the final state
    // Removed because of neutral kaon decay

    // InteractionAvatar::postInteraction(fs);
    // TODO: modify the InteractionAvatar::postInteraction(FinalState *fs)
    // in here to handle the final state from primary neutrino interaction
    // QEL RES MEC DIS
    modified = fs->getModifiedParticles();
    created = fs->getCreatedParticles();
    Destroyed = fs->getDestroyedParticles();

    for(ParticleIter i=modified.begin(), e=modified.end(); i!=e; ++i ){
      (*i)->makeParticipant(); //FIXME
      theNucleus->getStore()->getBook().incrementCascading(); // FIXME
    }
    for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i ){
      (*i)->makeParticipant(); //FIXME
      //theNucleus->getStore()->getBook().incrementCascading(); // FIXME
    }
    for(ParticleIter i=Destroyed.begin(), e=Destroyed.end(); i!=e; ++i ){
      theNucleus->getStore()->getBook().incrementCascading(); // FIXME
    }
    return;
  }

  std::string GENIEAvatar::dump() const {
    std::stringstream ss;
    ss << "(avatar " << theTime << " 'decay" << '\n'
      << "(list " << '\n'
      << particle1->dump()
      << "))" << '\n';
    return ss.str();
  }
}
