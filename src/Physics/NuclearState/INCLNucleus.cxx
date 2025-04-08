//____________________________________________________________________________
/*!

  \class    genie::INCLNucleus

  \brief    INCLXX nuclear model. Implements the NuclearModelI 
  interface.

  \ref      

  \author   Liang Liu, (liangliu@fnal.gov)

  \created  Oct. 2024

  \cpright  Copyright (c) 2003-2024, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#include <cassert>
#include <iostream>

#include <TSystem.h>
#include <TNtupleD.h>
#include <TTree.h>

#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/INCLNucleus.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"


#include "G4INCLCascade.hh"
#include "G4INCLRandom.hh"
#include "G4INCLStandardPropagationModel.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNuclearMassTable.hh"
#include "G4INCLGlobalInfo.hh"
#include "G4INCLNucleus.hh"

#include "G4INCLPauliBlocking.hh"

#include "G4INCLCrossSections.hh"

#include "G4INCLPhaseSpaceGenerator.hh"

#include "G4INCLLogger.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLNuclearDensityFactory.hh"

#include "G4INCLINuclearPotential.hh"

#include "G4INCLCoulombDistortion.hh"

#include "G4INCLClustering.hh"

#include "G4INCLIntersection.hh"

#include "G4INCLBinaryCollisionAvatar.hh"

#include "G4INCLCascadeAction.hh"
#include "G4INCLAvatarDumpAction.hh"

#include <cstring> 
#include <cstdlib>
#include <numeric>

#include "G4INCLPbarAtrestEntryChannel.hh"


#include "G4INCLGeant4Compat.hh"
#include "G4INCLCascade.hh"

#include "G4INCLClustering.hh"
#include "G4INCLParticle.hh"
#include "G4INCLIAvatar.hh"
#include "G4INCLIPropagationModel.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLStandardPropagationModel.hh"
#include "G4INCLRandom.hh"
#include "G4INCLRanecu.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLKinematicsUtils.hh"

// signal handler (for Linux and GCC)
#include "G4INCLSignalHandling.hh"

// For I/O
#include "IWriter.hh"
#include "ASCIIWriter.hh"
#include "ProtobufWriter.hh"
#include "INCLTree.hh"
#include "ROOTWriter.hh"
#include "HDF5Writer.hh"

// For configuration
#include "G4INCLConfig.hh"

// For logging
#include "G4INCLLogger.hh"

// Generic de-excitation interface
#include "G4INCLIDeExcitation.hh"

// ABLA v3p de-excitation
#ifdef INCL_DEEXCITATION_ABLAXX
#include "G4INCLAblaInterface.hh"
#endif

// ABLACXX de-excitation
#ifdef INCL_DEEXCITATION_ABLACXX
#include "G4INCLAblaxxInterface.hh"
#endif

// ABLA07 de-excitation
#ifdef INCL_DEEXCITATION_ABLA07
#include "G4INCLAbla07Interface.hh"
#endif

// SMM de-excitation
#ifdef INCL_DEEXCITATION_SMM
#include "G4INCLSMMInterface.hh"
#endif

// GEMINIXX de-excitation
#ifdef INCL_DEEXCITATION_GEMINIXX
#include "G4INCLGEMINIXXInterface.hh"
#endif



// INCL++
#include "G4INCLConfig.hh"
#include "G4INCLCascade.hh"
#include "G4INCLConfigEnums.hh"
#include "G4INCLParticle.hh"
// signal handler (for Linux and GCC)
#include "G4INCLSignalHandling.hh"

// Generic de-excitation interface
#include "G4INCLIDeExcitation.hh"

// ABLA v3p de-excitation
#ifdef INCL_DEEXCITATION_ABLAXX
#include "G4INCLAblaInterface.hh"
#endif

// ABLA07 de-excitation
#ifdef INCL_DEEXCITATION_ABLA07
#include "G4INCLAbla07Interface.hh"
#endif

// SMM de-excitation
#ifdef INCL_DEEXCITATION_SMM
#include "G4INCLSMMInterface.hh"
#endif

// GEMINIXX de-excitation
#ifdef INCL_DEEXCITATION_GEMINIXX
#include "G4INCLGEMINIXXInterface.hh"
#endif



using std::cout;
using std::endl;

using namespace genie;

//____________________________________________________________________________
INCLNucleus * INCLNucleus::fInstance = 0;

//____________________________________________________________________________
INCLNucleus::INCLNucleus():propagationModel_(0)
{
  //  this->Load();
  fInstance = 0;
  nucleon_index_ = -1;
  nucleus_ = nullptr;
  theConfig_ = nullptr;
  hitNucleon_ = nullptr;
}
//____________________________________________________________________________
INCLNucleus::~INCLNucleus()
{
  //  if(!gAbortingInErr) {
  //    cout << "INCLNucleus singleton dtor: Deleting inputs... " << endl;
  //  }
  //  delete fNuclSupprD2;
}
//____________________________________________________________________________
INCLNucleus * INCLNucleus::Instance()
{
  if(fInstance == 0) {
    LOG("INCLNucleus", pINFO) << "INCLNucleus late initialization";
    fInstance = new INCLNucleus;
    fInstance->theConfig_ = new G4INCL::Config();
    fInstance->nucleus_ = nullptr;
    fInstance->theDensityForLepton = nullptr;
  }
  return fInstance;
}

void INCLNucleus::configure(){

  theConfig_->init();
  theConfig_->setINCLXXDataFilePath(INCLXXDataFilePath_); // FIXME:: using config to set path
  theConfig_->setABLAXXDataFilePath(ablaxxDataFilePath_);
  theConfig_->setcABLA07DataFilePath(abla07DataFilePath_);
  theConfig_->setGEMINIXXDataFilePath(geminixxDataFilePath_);
  theConfig_->setDeExcitationType(deExcitationType_);

  //theConfig_->setRPCorrelationCoefficient(1.0); // Using r-p correlation without fuzzy

  // initialize INCL model
  G4INCL::Random::initialize(theConfig_);
  // Select the Pauli and CDPP blocking algorithms
  G4INCL::Pauli::initialize(theConfig_);
  // Set the cross-section set
  G4INCL::CrossSections::initialize(theConfig_);
  // Set the phase-space generator
  G4INCL::PhaseSpaceGenerator::initialize(theConfig_);
  // Initialize the INCL particle table:
  G4INCL::ParticleTable::initialize(theConfig_);
  // Select the Coulomb-distortion algorithm:
  G4INCL::CoulombDistortion::initialize(theConfig_);
  // Select the clustering algorithm:
  G4INCL::Clustering::initialize(theConfig_);
  // Initialize the value of cutNN in BinaryCollisionAvatar
  G4INCL::BinaryCollisionAvatar::setCutNN(theConfig_->getCutNN());
  // Initialize the value of strange cross section bias
  G4INCL::BinaryCollisionAvatar::setBias(theConfig_->getBias());



  //theConfig_->setLocalEnergyBBType(G4INCL::NeverLocalEnergy);
  //theConfig_->setLocalEnergyPiType(G4INCL::NeverLocalEnergy);


  // Propagation model is responsible for finding avatars and
  // transporting the particles. In principle this step is "hidden"
  // behind an abstract interface and the rest of the system does not
  // care how the transportation and avatar finding is done. This
  // should allow us to "easily" experiment with different avatar
  // finding schemes and even to support things like curved
  // trajectories in the future.
  propagationModel_ = new G4INCL::StandardPropagationModel(theConfig_->getLocalEnergyBBType(),theConfig_->getLocalEnergyPiType(),theConfig_->getHadronizationTime());
  if(theConfig_->getCascadeActionType() == G4INCL::AvatarDumpActionType)
    cascadeAction_ = new G4INCL::AvatarDumpAction();
  else
    cascadeAction_ = new G4INCL::CascadeAction();
//  cascadeAction_->beforeRunAction(theConfig_);
}

void INCLNucleus::initialize(const Target * tgt){
  // Skip to initialize a new nucleus if the INCL nucleus is not empty 
  // and the same species with GENIE target 
  

  // Initialize the INCL particle table:
  // FIXME: this is a temp set, 
  if(!theDensityForLepton){
  theConfig_->setRPCorrelationCoefficient(1.0);
  G4INCL::ParticleTable::initialize(theConfig_);
  // create the density that will use in neutrino interactiona
  theDensityForLepton = G4INCL::NuclearDensityFactoryHitByLepton::createDensity(tgt->A(), tgt->Z(), 0); // parameters are  A, Z, S = 0;
  theConfig_->setRPCorrelationCoefficient(G4INCL::Proton, 0.5); // FIXME, the default value for proton r-p correlation
  theConfig_->setRPCorrelationCoefficient(G4INCL::Neutron, 0.73); // FIXME, the default value for neutron r-p correlation
  G4INCL::ParticleTable::initialize(theConfig_);
  }


  if(nucleus_){
    if(!nucleus_->getStore()->getParticles().empty())
      if(nucleus_->getA() == tgt->A() && nucleus_->getZ() == tgt->Z())
	return ;
  }
  // initialize according process Event in INCL
  G4INCL::ParticleSpecies targetSpecies = G4INCL::ParticleSpecies(tgt->A(), tgt->Z());
  theConfig_->setTargetA(targetSpecies.theA);
  theConfig_->setTargetZ(targetSpecies.theZ);
  theConfig_->setTargetS(targetSpecies.theS);
  // define Nucleus and initialize it
  
  // ReInitialize the bias vector
  G4INCL::Particle::INCLBiasVector.clear();
  //Particle::INCLBiasVector.Clear();
  G4INCL::Particle::nextBiasedCollisionID = 0;

  // Set the target and the projectile 
  // implement prepare reaction
  
  // Reset the forced-transparent flag
  // forceTransparent = false; FIXME
  //
  // Initialise the maximum universe radius
  // INCL initialize universe radius according to particle species,
  // kenetic energy, and nucleus type. void INCL::initUniverseRadius(ParticleSpecies const &p, const double kineticEnergy, const int A, const int Z)
  // FIXME: I'm not sure if we need the interaction distance for neutrino scattering
  this->initUniverseRadius(targetSpecies.theA, targetSpecies.theZ);

  // FIXME the last two parameters need to be configed
  // theConfig_, G4INCL::ParticleTable::getMaximumNuclearRadius(G4INCL::Proton, targetSpecies.theA, targetSpecies.theZ)
  // G4INCL::NType
  if(nucleus_){
    delete nucleus_;
  }
  nucleus_ = new G4INCL::Nucleus(targetSpecies.theA, targetSpecies.theZ, targetSpecies.theS, 
  theConfig_, maxUniverseRadius_, G4INCL::Def);
  nucleus_->getStore()->getBook().reset();


  // sample the index of nucleon hitted by lepton
  int nucleon_pdg = tgt->HitNucPdg();
  LOG("INCLNucleus", pDEBUG) << "hit nucleon pdg : " << nucleon_pdg;

  if(nucleon_pdg == kPdgClusterNN){
    // randomly pick a neutron and then use the closest neutron to form the cluster NN
    clusterNN_ = getNNCluster(kPdgNeutron, kPdgNeutron);
  }
  else if(nucleon_pdg == kPdgClusterNP){
    // randomly pick a neutron and then use the closest proton to form the cluster NN
    clusterNN_ = getNNCluster(kPdgNeutron, kPdgProton);
  }
  else if(nucleon_pdg == kPdgClusterPP){
    // randomly pick a proton and then use the closest proton to form the cluster NN
    clusterNN_ = getNNCluster(kPdgProton, kPdgProton);
  }
  else if(pdg::IsProton(nucleon_pdg) || pdg::IsNeutron(nucleon_pdg)){
    hitNucleon_ = this->getNucleon(nucleon_pdg);
  }
  else{
    LOG("INCLNucleus", pFATAL) << "Can't get a valid nucleon! " << nucleon_pdg;
    exit(1);
  }

  propagationModel_->setNucleus(nucleus_);

  //  TODO :: delete the comment
  // set the index of nucleon hitted by lepton before initialize a nucleus
  // FIXME
//  nucleus_->setLeptonScatteringDensity(theDensityForLepton);
//  nucleus_->setLeptonHitNucleonIndex(nucleon_index_);
  //nucleus_->initializeParticles();



//  propagationModel_->setNucleus(nucleus_);
//  if(hitNucleon_) hitNucleon_ = nullptr;
//  hitNucleon_ = nucleus_->getStore()->getParticles().at(nucleon_index_);
//  LOG("INCLNucleus", pNOTICE) << hitNucleon_->print();

  // initialize max interaction distance
  // FIXME: in INCL, composite has non-zero max interaction distance.
  
  // maxInteractionDistance_ = 0;

  // set the min Remnant size to be 4
  // the min remnant is alpha particle
  // FIXME: it only works for nuclei with large A
  //
  minRemnantSize_ = 4;

  // cascade action is not related to simulation
  // it is just output the casade to file FIXME
  //cascadeAction_->beforeCascadeAction(propagationModel_);
  //
  // INCL need to decide whether the cascade can be ran or not
  // For genie, we need to run casecade for every events
  // const bool canRunCascade = preCascade(projectileSpecies, kineticEnergy);
  //
  LOG("INCLNucleus", pDEBUG) << nucleus_->getStore()->getParticles().at(2)->getPotentialEnergy() ;
  LOG("INCLNucleus", pDEBUG) << nucleus_->getStore()->getParticles().at(2)->getEnergy() -  nucleus_->getStore()->getParticles().at(2)->getPotentialEnergy()  ;
  LOG("INCLNucleus", pDEBUG) << nucleus_->getStore()->getParticles().at(2)->getMomentum().print() ;
  LOG("INCLNucleus", pDEBUG) << nucleus_->getStore()->getParticles().at(2)->getPosition().print() ;

}

void INCLNucleus::reset(const Target * tgt){
  // nucleus must exsit!
  if(!nucleus_)  LOG("INCLNucleus", pFATAL) << "nucleus doesn't exsit!";
  // can't reset a nucleus with different type
  if(!(nucleus_->getA() == tgt->A() && nucleus_->getZ() == tgt->Z())) 
    LOG("INCLNucleus", pFATAL) << "you are try to reset a nucleus with different type!";
  // reset the nucleus
  if(nucleus_){
    nucleus_->deleteParticles();
    nucleus_->getStore()->clear();
    nucleus_->getStore()->getBook().reset();

    // sample the index of nucleon hitted by lepton
    int nucleon_pdg = tgt->HitNucPdg();

    if(nucleon_pdg == kPdgClusterNN){
      // randomly pick a neutron and then use the closest neutron to form the cluster NN
      clusterNN_ = getNNCluster(kPdgNeutron, kPdgNeutron);
    }
    else if(nucleon_pdg == kPdgClusterNP){
      // randomly pick a neutron and then use the closest proton to form the cluster NN
      clusterNN_ = getNNCluster(kPdgNeutron, kPdgProton);
    }
    else if(nucleon_pdg == kPdgClusterPP){
      // randomly pick a proton and then use the closest proton to form the cluster NN
      clusterNN_ = getNNCluster(kPdgProton, kPdgProton);
    }
    else if(pdg::IsProton(nucleon_pdg) || pdg::IsNeutron(nucleon_pdg)){
      hitNucleon_ = this->getNucleon(nucleon_pdg);
    }
    else{
      LOG("INCLNucleus", pFATAL) << "Can't get a valid nucleon! " << nucleon_pdg;
      exit(1);
    }

    // set the index of nucleon hitted by lepton before initialize a nucleus
    LOG("INCLNucleus", pNOTICE) << nucleon_index_;
//    nucleus_->setLeptonScatteringDensity(theDensityForLepton);
//    nucleus_->setLeptonHitNucleonIndex(nucleon_index_);
//    nucleus_->initializeParticles();
    // reset the hit nucleon
//    if(hitNucleon_) hitNucleon_ = nullptr;
//    hitNucleon_ = nucleus_->getStore()->getParticles().at(nucleon_index_);
//    LOG("INCLNucleus", pNOTICE) << hitNucleon_->print();

    if(propagationModel_){
      delete propagationModel_;
    }
    propagationModel_ = new G4INCL::StandardPropagationModel(theConfig_->getLocalEnergyBBType(),theConfig_->getLocalEnergyPiType(),theConfig_->getHadronizationTime());
    propagationModel_->setNucleus(nucleus_);
  }
}

TVector3 INCLNucleus::getHitNucleonPosition(){
  if(!hitNucleon_){
	LOG("INCLNucleus", pFATAL) << "hit nucleon is not valid!";
	exit(1);
  }
  TVector3 v3(999999.,999999.,999999.);
  v3.SetXYZ(hitNucleon_->getPosition().getX(),
      hitNucleon_->getPosition().getY(),
      hitNucleon_->getPosition().getZ());
  return v3;
}
TVector3 INCLNucleus::getHitNucleonMomentum(){
  if(!hitNucleon_){
	LOG("INCLNucleus", pFATAL) << "hit nucleon is not valid!";
	exit(1);
  }
  TVector3 p3(999999.,999999.,999999.);
  p3.SetXYZ(hitNucleon_->getMomentum().getX(),
      hitNucleon_->getMomentum().getY(),
      hitNucleon_->getMomentum().getZ());
  return p3;
}
double INCLNucleus::getHitNucleonEnergy(){
  if(!hitNucleon_){
	LOG("INCLNucleus", pFATAL) << "hit nucleon is not valid!";
	exit(1);
  }
  return hitNucleon_->getEnergy();
}

double INCLNucleus::getHitNucleonMass(){
  if(!hitNucleon_){
	LOG("INCLNucleus", pFATAL) << "hit nucleon is not valid!";
	exit(1);
  }
  return hitNucleon_->getMass();
}

double INCLNucleus::getMass(){
  if(!nucleus_){
	LOG("INCLNucleus", pFATAL) << "nucleus is not valid!";
	exit(1);
  }
  return nucleus_->getMass();
}

G4INCL::Nucleus * INCLNucleus::getNuclues(){
  if(!nucleus_){
	LOG("INCLNucleus", pFATAL) << "nucleus is not valid!";
	exit(1);
  }
  return nucleus_;
}

G4INCL::Particle * INCLNucleus::getHitParticle(){
  if(!hitNucleon_){
    LOG("INCLNucleus", pFATAL) << "nucleus is not valid!";
    exit(1);
  }
  return hitNucleon_;
}
G4INCL::Cluster * INCLNucleus::getHitNNCluster(){
  if(!clusterNN_){
    LOG("INCLNucleus", pFATAL) << "cluster is not valid!";
    exit(1);
  }
  return clusterNN_;
}


G4INCL::StandardPropagationModel * INCLNucleus::getPropagationModel(){
  return propagationModel_;
}

double INCLNucleus::getRemovalEnergy(){
  if(!hitNucleon_){
	LOG("INCLNucleus", pFATAL) << "hit nucleon is not valid!";
	exit(1);
  }
  // FIXME: need to find the correct way for removal energy
 //   double removal_energy = 0;
 //   double nucleon_mass = hitNucleon_->getRealMass();
 //   double mag = hitNucleon_->getMomentum().mag();
 //   removal_energy = TMath::Sqrt(mag*mag + nucleon_mass*nucleon_mass) - hitNucleon_->getEnergy();
    // return removal_energy;
  //return hitNucleon_->getPotentialEnergy();
  //return nucleus_->getPotential()->getSeparationEnergy(hitNucleon_->getType());
  // FIXME: removal = potential - qvalue
  double removal = hitNucleon_->getPotentialEnergy() - hitNucleon_->getEmissionQValueCorrection(nucleus_->getA(),nucleus_->getZ(),nucleus_->getS());
  return removal;

}

void INCLNucleus::initUniverseRadius(const int A, const int Z){
  // This function is analogy to function in incl_physics/src/G4INCLCascade.cc
  // void INCL::initUniverseRadius(ParticleSpecies const &p, 
  //                const double kineticEnergy, const int A, 
  //                const int Z)
  double rMax = 0.0;
  // A should be large than 0
  // FIXME: 
  // 1. do we need to consider the isotopes?
  // 2. do we need to consider the extra-impact parameter for neutrino?
  // 	The xsec for neutrino-nucleus is ~10 fb
  // 	the xsec for hadron-nucleus is ~800 mb
  if(!(A > 0)) 
    throw std::runtime_error("Mass number A is not real!");
  const double pMaximumRadius = G4INCL::ParticleTable::getMaximumNuclearRadius(G4INCL::Proton,  A, Z);
  const double nMaximumRadius = G4INCL::ParticleTable::getMaximumNuclearRadius(G4INCL::Neutron, A, Z);
  const double maximumRadius = std::max(pMaximumRadius, nMaximumRadius);
  rMax = std::max(maximumRadius, rMax);
  maxUniverseRadius_ = rMax;
//  LOG("INCLNucleus", pINFO) << "max Universe Radius : " << maxUniverseRadius_; 
}

G4INCL::Cluster* INCLNucleus::getNNCluster(const int pdg1, const int pdg2){

  nucleus_->initializeParticles();
  LOG("INCLNucleus", pINFO) << "get cluster";
  cluster_index1_ = -1;
  cluster_index2_ = -1;
  RandomGen * rnd = RandomGen::Instance();
  G4INCL::Particle *cluster_N1;
  G4INCL::Particle *cluster_N2;
  if(pdg::IsProton(pdg1))
    cluster_index1_ = rnd->RndGen().Integer(nucleus_->getZ());
  else if(pdg::IsNeutron(pdg1))
    cluster_index1_ = rnd->RndGen().Integer(nucleus_->getA() - nucleus_->getZ()) + nucleus_->getZ();
  else {
    LOG("INCLNucleus", pFATAL) << "Can't get a valid nucleon! " << pdg1;
    exit(1);
  }

  cluster_N1 = nucleus_->getStore()->getParticles().at(cluster_index1_);
  G4INCL::ParticleList const &particles = nucleus_->getStore()->getParticles();
  if(pdg::IsProton(pdg2)){
    double size = 10e16;
    for(G4INCL::ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      if((*i)->getID() == cluster_N1->getID()) continue;
      if((*i)->getType() != G4INCL::Proton) continue;
      double space    = ((*i)->getPosition() - cluster_N1->getPosition()).mag2();
      double momentum = ((*i)->getMomentum() - cluster_N1->getMomentum()).mag2();
      double temp_size     = space;
      if(temp_size < size){
	size =  temp_size;
	cluster_N2 = (*i);
      }
    }
  }
  else if(pdg::IsNeutron(pdg2)){
    double size = 10e16;
    for(G4INCL::ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      if((*i)->getID() == cluster_N1->getID()) continue;
      if((*i)->getType() != G4INCL::Neutron) continue;
      double space    = ((*i)->getPosition() - cluster_N1->getPosition()).mag2();
      double momentum = ((*i)->getMomentum() - cluster_N1->getMomentum()).mag2();
      double temp_size     = space;
      if(temp_size < size){
	size =  temp_size;
	cluster_N2 = (*i);
      }
    }
  }
  else{
    LOG("INCLNucleus", pFATAL) << "Can't get a valid nucleon! " << pdg2;
    exit(1);
  }
  G4INCL::ParticleList selectedParticles;
  selectedParticles.push_back(cluster_N1);
  selectedParticles.push_back(cluster_N2);
  return (new G4INCL::Cluster(selectedParticles.begin(), selectedParticles.end()));

}

G4INCL::Particle * INCLNucleus::getNucleon(const int pdg){

  RandomGen * rnd = RandomGen::Instance();
  nucleon_index_ = -1;
  if(pdg::IsProton(pdg))
    nucleon_index_ = rnd->RndGen().Integer(nucleus_->getZ());
  else if(pdg::IsNeutron(pdg))
    nucleon_index_ = rnd->RndGen().Integer(nucleus_->getA() - nucleus_->getZ()) + nucleus_->getZ();
  else{
    LOG("INCLNucleus", pFATAL) << "Can't get a valid nucleon!";
    exit(1);
  }
    // set the index of nucleon hitted by lepton before initialize a nucleus
    LOG("INCLNucleus", pNOTICE) << nucleon_index_;
//    nucleus_->setLeptonScatteringDensity(theDensityForLepton);
//    nucleus_->setLeptonHitNucleonIndex(nucleon_index_);
    nucleus_->initializeParticles();
    // reset the hit nucleon
    return nucleus_->getStore()->getParticles().at(nucleon_index_);
}

#endif // __GENIE_INCL_ENABLED__

