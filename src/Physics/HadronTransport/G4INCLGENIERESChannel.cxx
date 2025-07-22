#include "G4INCLGENIERESChannel.h"
#include "G4INCLRandom.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLCrossSections.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  GENIERESChannel::GENIERESChannel(Particle *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord)
    : hitParticle(p), theNucleus(n),
    genie_evtrec(eventRecord)
  {
  }

  GENIERESChannel::~GENIERESChannel()
  {
  }

  void GENIERESChannel::fillFinalState(FinalState *fs)
  {

    std::vector<GENIEParticleRecord>::iterator ip;
    for(ip = genie_evtrec->begin(); ip != genie_evtrec->end(); ip++){
      if(ip->Status() == 13){
        ip->setID(int(hitParticle->getID()));
        hitParticle->setType(ip->Type());
        hitParticle->setMass(ip->Mass());
        hitParticle->setMomentum(ip->P3());
        hitParticle->setPosition(ip->X3());
        hitParticle->adjustEnergyFromMomentum();
        fs->addModifiedParticle(hitParticle);
      }
      else if((ip->Status() == 14 && (ip->Type() == Proton || ip->Type() == Neutron))){
        ip->setID(int(hitParticle->getID()));
        hitParticle->setType(ip->Type());
        hitParticle->setMass(ip->Mass());
        hitParticle->setMomentum(ip->P3());
        hitParticle->setPosition(ip->X3());
        hitParticle->adjustEnergyFromMomentum();
        double energy = hitParticle->getEnergy() + hitParticle->getPotentialEnergy() - hitParticle->getEmissionQValueCorrection(theNucleus->getA(), theNucleus->getZ(), theNucleus->getS());  // FIXME
        hitParticle->setEnergy(energy); // FIXME
        hitParticle->adjustMomentumFromEnergy(); // FIXME
        fs->addModifiedParticle(hitParticle);
      }
      else if( (ip->Status() == 14) ){
        Particle *ihadron = new Particle(ip->Type(), ip->P3(), ip->X3());
        ip->setID(int(ihadron->getID()));
        ihadron->setMass(ip->Mass());
        ihadron->adjustEnergyFromMomentum();
        fs->addCreatedParticle(ihadron);
      }
    }
  }
}
