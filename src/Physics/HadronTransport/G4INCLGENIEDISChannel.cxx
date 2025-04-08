#include "G4INCLGENIEDISChannel.h"
#include "G4INCLRandom.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLCrossSections.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  GENIEDISChannel::GENIEDISChannel(Particle *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord)
    :hitParticle(p), theNucleus(n),
    genie_evtrec(eventRecord)
  {
  }

  GENIEDISChannel::~GENIEDISChannel()
  {
  }

  void GENIEDISChannel::fillFinalState(FinalState *fs)
  {
    std::vector<GENIEParticleRecord>::iterator ip;

    for(ip = genie_evtrec->begin(); ip != genie_evtrec->end(); ip++){
      std::cout << "ip: " << ip->Pdg() << " " << ip->FirstMother() << "  " << ip->Status()  << "  " << hitParticle->getID() << std::endl;
      if(ip->Status() == 14){
	Particle *ihadron = new Particle(ip->Type(), ip->P3(), ip->X3());
	ip->setID(int(ihadron->getID()));
	ihadron->setMass(ip->Mass());
	ihadron->adjustEnergyFromMomentum();
	fs->addCreatedParticle(ihadron);
      }
    }
    fs->addDestroyedParticle(hitParticle);

    //exit(1);

  }

}
