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
    delta_Z = 0;
    delta_Z += particle1->getZ();
  }


  GENIEAvatar::GENIEAvatar(double time, Cluster *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord)
    : InteractionAvatar(time, n, nullptr),
    theNucleus(n), cluster(p),
    genie_evtrec(eventRecord)
  {
    //    setType(GENIEAvatarType);
    delta_Z = 0;
    delta_Z += cluster->getZ();
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
    if((*genie_evtrec)[0].ScatteringType() != 10){
      int index = 0;
      double lepton_initial_energy = 0;
      ThreeVector leptonInitialMom;
      std::vector<GENIEParticleRecord>::iterator ip;
      for(ip = genie_evtrec->begin(); ip != genie_evtrec->end(); ip++){
        //std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << "  " << ip->ID() << " "
        //  << ip->Pdg() << " " << ip->Mass() << " " << std::sqrt(ip->P3().mag()*ip->P3().mag() + ip->Mass()*ip->Mass()) << std::endl;
        if(ip->RecordCode() == kProbe){
          lepton_initial_energy = std::sqrt(ip->P3().mag2() + ip->Mass()*ip->Mass());
          leptonInitialMom = ip->P3();
        }
        else if(ip->RecordCode() == kFinalStateLepton){
          leptonE = std::sqrt(ip->P3().mag2() + ip->Mass()*ip->Mass());
          leptonMom = ip->P3();
        }
        else if(ip->RecordCode() == kHitNucleon){
          if(ip->ScatteringType() != 10){
            ThreeVector n_mom = particle1->getMomentum();
            ip->setMomentum(n_mom);
            ip->setMass(particle1->getMass());
          }
        }
        index++;
      }
      oldTotalEnergy = lepton_initial_energy + particle1->getEnergy() - particle1->getPotentialEnergy();

      // transfrom the target nucleon to local energy frame
      KinematicsUtils::transformToLocalEnergyFrame(theNucleus, particle1);

      // make boost vector
      boostVector = (leptonInitialMom + particle1->getMomentum())/(lepton_initial_energy + particle1->getEnergy());
    }
    else{
      int index = 0;
      double lepton_initial_energy = 0;
      ThreeVector leptonInitialMom;
      std::vector<GENIEParticleRecord>::iterator ip;
      for(ip = genie_evtrec->begin(); ip != genie_evtrec->end(); ip++){
        if(ip->RecordCode() == kProbe){
          lepton_initial_energy = std::sqrt(ip->P3().mag2() + ip->Mass()*ip->Mass());
          leptonInitialMom = ip->P3();
        }
        else if(ip->RecordCode() == kFinalStateLepton){
          leptonE = std::sqrt(ip->P3().mag2() + ip->Mass()*ip->Mass());
          leptonMom = ip->P3();
        }
        else if(ip->RecordCode() == kHitNucleon){
        }
        index++;
      }

      G4INCL::ParticleList particles = cluster->getParticleList();
      oldTotalEnergy = lepton_initial_energy;

      ThreeVector local_mom = leptonInitialMom;
      double local_energy = lepton_initial_energy;
      for(G4INCL::ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
        oldTotalEnergy += (*i)->getEnergy() - (*i)->getPotentialEnergy();
        KinematicsUtils::transformToLocalEnergyFrame(theNucleus, (*i));
        local_mom += (*i)->getMomentum();
        local_energy += (*i)->getEnergy();
      }
      boostVector = local_mom / local_energy;
    }
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
    modifiedAndCreated = modified;
    modifiedAndCreated.insert(modifiedAndCreated.end(), created.begin(), created.end());


    bool success = enforceEnergyConservation(fs);
    if(!success){
      std::cout << "enforceEnergyConservation fs wrong!" << std::endl;
      fs->reset();
      fs->makeNoEnergyConservation();
      fs->setTotalEnergyBeforeInteraction(0.0);
      return; // Interaction is blocked. Return an empty final state.
    }

    for(ParticleIter i=modified.begin(), e=modified.end(); i!=e; ++i ){
      std::cout << "DEBUG " << (*i)->print() << std::endl;
    }

    /// update the event record after call enforceEnergyConservation;
    int index = 0;
    std::vector<GENIEParticleRecord>::iterator ip;
    ParticleIter imc = modifiedAndCreated.begin();
    for(ip = genie_evtrec->begin(); ip != genie_evtrec->end(); ip++){
      if(ip->RecordCode() == kFinalStateLepton){
        ip->setMomentum(leptonMom);
        ip->setMass(std::sqrt(std::max(leptonE*leptonE - leptonMom.mag2(), 0.)));
      }
      else if(ip->Status() == 14 || ip->Status() == 13){
        ThreeVector p_mom = (*imc)->getMomentum();
        ip->setMomentum(p_mom);
        ip->setMass((*imc)->getMass());
        imc++;
      }
      index++;
    }

    for(ParticleIter i=modified.begin(), e=modified.end(); i!=e; ++i ){
      (*i)->makeParticipant(); //FIXME
      delta_Z -= (*i)->getZ();
      theNucleus->getStore()->getBook().incrementCascading(); // FIXME
    }
    for(ParticleIter i=created.begin(), e=created.end(); i!=e; ++i ){
      (*i)->makeParticipant(); //FIXME
      delta_Z -= (*i)->getZ();
      //theNucleus->getStore()->getBook().incrementCascading(); // FIXME
    }
    for(ParticleIter i=Destroyed.begin(), e=Destroyed.end(); i!=e; ++i ){
      theNucleus->getStore()->getBook().incrementCascading(); // FIXME
    }
    theNucleus->setZ(theNucleus->getZ() - delta_Z);
    return;
  }

  bool GENIEAvatar::enforceEnergyConservation(FinalState * const fs){
    // Set up the violationE calculation
    violationEFunctor = new ViolationLeptonEMomentumFunctor(theNucleus, modifiedAndCreated, leptonE, leptonMom, boostVector, oldTotalEnergy, true);
    const RootFinder::Solution theSolution = RootFinder::solve(violationEFunctor, 1.0);
    if(theSolution.success) { // Apply the solution
      (*violationEFunctor)(theSolution.x);
    } else if(theNucleus){
      theNucleus->getStore()->getBook().incrementEnergyViolationInteraction();
    }
    delete violationEFunctor;
    violationEFunctor = NULL;
    return theSolution.success;

  }

  /* ***                                                      ***
   * *** InteractionAvatar::ViolationEMomentumFunctor methods ***
   * ***                                                      ***/

  GENIEAvatar::ViolationLeptonEMomentumFunctor::ViolationLeptonEMomentumFunctor(Nucleus * const nucleus, ParticleList const &modAndCre,
      double &leptonE, ThreeVector &leptonMom, ThreeVector const &boost,
      const double totalEnergyBeforeInteraction, const bool localE) :
    RootFunctor(0., 1E6),
    finalParticles(modAndCre),
    leptonEnergy(leptonE),
    leptonMomentum(leptonMom),
    initialEnergy(totalEnergyBeforeInteraction),
    theNucleus(nucleus),
    boostVector(boost),
    shouldUseLocalEnergy(localE)
  {
    // Store the particle momenta (necessary for the calls to
    // scaleParticleMomenta() to work)
    lepton = new Lepton(leptonMomentum, leptonEnergy);
    lepton->boost(boostVector);
    leptonMomentumCM = lepton->getMomentum();
    for(ParticleIter i=finalParticles.begin(), e=finalParticles.end(); i!=e; ++i) {
      (*i)->boost(boostVector);
      particleMomenta.push_back((*i)->getMomentum());
    }
  }

  GENIEAvatar::ViolationLeptonEMomentumFunctor::~ViolationLeptonEMomentumFunctor() {
    particleMomenta.clear();
    delete lepton;
  }


  double GENIEAvatar::ViolationLeptonEMomentumFunctor::operator()(const double alpha) const {
      //LOG("INCLCascadeIntranuke", pNOTICE) << "alpha : " << alpha;
    scaleParticleMomenta(alpha);

    double deltaE = 0.0;
    for(ParticleIter i=finalParticles.begin(), e=finalParticles.end(); i!=e; ++i){
      deltaE += (*i)->getEnergy() - (*i)->getPotentialEnergy();
    }
    deltaE += leptonEnergy;
    deltaE -= initialEnergy;
    return deltaE;
  }

  void GENIEAvatar::ViolationLeptonEMomentumFunctor::scaleParticleMomenta(const double alpha) const {

    lepton->setMomentum(leptonMomentumCM*alpha);
    lepton->adjustEnergyFromMomentum();
    lepton->boost(-boostVector);
    leptonEnergy = lepton->getEnergy();
    leptonMomentum = lepton->getMomentum();

    std::vector<ThreeVector>::const_iterator iP = particleMomenta.begin();
    for(ParticleIter i=finalParticles.begin(), e=finalParticles.end(); i!=e; ++i, ++iP) {
      (*i)->setMomentum((*iP)*alpha);
      (*i)->adjustEnergyFromMomentum();
      (*i)->rpCorrelate();
      (*i)->boost(-boostVector);
      if(theNucleus){
        theNucleus->updatePotentialEnergy(*i);
      } else {
        (*i)->setPotentialEnergy(0.);
      }

      //jcd      if(shouldUseLocalEnergy && !(*i)->isPion())  // This translates AECSVT's loops 1, 3 and 4
      if(shouldUseLocalEnergy && !(*i)->isPion() && !(*i)->isEta() && !(*i)->isOmega() &&
          !(*i)->isKaon() && !(*i)->isAntiKaon()  && !(*i)->isSigma() && !(*i)->isPhoton() && !(*i)->isLambda() && !(*i)->isAntiNucleon()) { // This translates AECSVT's loops 1, 3 and 4
        assert(theNucleus); // Local energy without a nucleus doesn't make sense
        const double energy = (*i)->getEnergy(); // Store the energy of the particle
        double locE = KinematicsUtils::getLocalEnergy(theNucleus, *i); // Initial value of local energy
        double locEOld;
        double deltaLocE = InteractionAvatar::locEAccuracy + 1E3;
        for(int iterLocE=0;
            deltaLocE>InteractionAvatar::locEAccuracy && iterLocE<InteractionAvatar::maxIterLocE;
            ++iterLocE) {
          locEOld = locE;
          (*i)->setEnergy(energy + locE); // Update the energy of the particle...
          (*i)->adjustMomentumFromEnergy();
          theNucleus->updatePotentialEnergy(*i); // ...update its potential energy...
          locE = KinematicsUtils::getLocalEnergy(theNucleus, *i); // ...and recompute locE.
          deltaLocE = std::abs(locE-locEOld);
        }
      }

      //jlrs  For lambdas and nuclei with masses higher than 19 also local energy
      if(shouldUseLocalEnergy && (*i)->isLambda() && theNucleus->getA()>19) {
        assert(theNucleus); // Local energy without a nucleus doesn't make sense
        const double energy = (*i)->getEnergy(); // Store the energy of the particle
        double locE = KinematicsUtils::getLocalEnergy(theNucleus, *i); // Initial value of local energy
        double locEOld;
        double deltaLocE = InteractionAvatar::locEAccuracy + 1E3;
        for(int iterLocE=0;
            deltaLocE>InteractionAvatar::locEAccuracy && iterLocE<InteractionAvatar::maxIterLocE;
            ++iterLocE) {
          locEOld = locE;
          (*i)->setEnergy(energy + locE); // Update the energy of the particle...
          (*i)->adjustMomentumFromEnergy();
          theNucleus->updatePotentialEnergy(*i); // ...update its potential energy...
          locE = KinematicsUtils::getLocalEnergy(theNucleus, *i); // ...and recompute locE.
          deltaLocE = std::abs(locE-locEOld);
        }
      }
    }
  }

  void GENIEAvatar::ViolationLeptonEMomentumFunctor::cleanUp(const bool success) const {
    if(!success)
      scaleParticleMomenta(1.);
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
