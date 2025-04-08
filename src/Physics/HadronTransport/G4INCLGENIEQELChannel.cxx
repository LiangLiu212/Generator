#include "G4INCLGENIEQELChannel.h"
#include "G4INCLRandom.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLCrossSections.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  GENIEQELChannel::GENIEQELChannel(Particle *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord)
    :hitParticle(p), theNucleus(n),
    genie_evtrec(eventRecord)
  {
  }

  GENIEQELChannel::~GENIEQELChannel()
  {
  }

  void GENIEQELChannel::fillFinalState(FinalState *fs)
  {
    std::vector<GENIEParticleRecord>::iterator ip;

    for(ip = genie_evtrec->begin(); ip != genie_evtrec->end(); ip++){
      std::cout << "ip: " << ip->Pdg() << " " << ip->FirstMother() << "  " << ip->Status()  << "  " << hitParticle->getID() << std::endl;
      if(ip->Status() == 14){
	ip->setID(int(hitParticle->getID()));
	hitParticle->setType(ip->Type());
	hitParticle->setMomentum(ip->P3());
	hitParticle->setPosition(ip->X3());
	hitParticle->adjustEnergyFromMomentum();
        double energy = hitParticle->getEnergy() + hitParticle->getPotentialEnergy() - hitParticle->getEmissionQValueCorrection(theNucleus->getA(), theNucleus->getZ(), theNucleus->getS());  // FIXME
	hitParticle->setEnergy(energy); // FIXME
	hitParticle->adjustMomentumFromEnergy(); // FIXME
      }
    }
    fs->addModifiedParticle(hitParticle);
  }

}
