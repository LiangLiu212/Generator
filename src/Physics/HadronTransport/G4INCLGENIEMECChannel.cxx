#include "G4INCLGENIEMECChannel.h"
#include "G4INCLRandom.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLCrossSections.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  GENIEMECChannel::GENIEMECChannel(Cluster *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord): cluster(p), theNucleus(n),
    genie_evtrec(eventRecord)
  {
  }

  GENIEMECChannel::~GENIEMECChannel()
  {
  }

  void GENIEMECChannel::fillFinalState(FinalState *fs)
  {
    ParticleList const &particles = cluster->getParticles();
    int cluster_index = 0;
    std::vector<GENIEParticleRecord>::iterator ip;
    for(ip = genie_evtrec->begin(); ip != genie_evtrec->end(); ip++){
      std::cout << "ip: " << ip->Pdg() << " " << ip->FirstMother() << "  " << ip->Status()  << "  " <<  std::endl;
      if(ip->Status() == 14){
	ip->setID(int(particles.at(cluster_index)->getID()));
	particles.at(cluster_index)->setType(ip->Type());
	particles.at(cluster_index)->setMomentum(ip->P3());
	particles.at(cluster_index)->adjustEnergyFromMomentum();
//	double energy = particles.at(cluster_index)->getEnergy() +
//	  particles.at(cluster_index)->getPotentialEnergy() -
//	  particles.at(cluster_index)->getEmissionQValueCorrection(theNucleus->getA(), theNucleus->getZ(), theNucleus->getS()); // FIXME
//	particles.at(cluster_index)->setEnergy(energy); // FIXME
//	particles.at(cluster_index)->adjustMomentumFromEnergy(); // FIXME
	fs->addModifiedParticle(particles.at(cluster_index));
	cluster_index++;
      }
    }
    

    std::cout << "success!" << std::endl;
  }

}
