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
INCLNucleus::INCLNucleus()
{
  //  this->Load();
  fInstance = 0;
  nucleon_index = -1;
  nucleus_ = nullptr;
  theConfig_ = nullptr;
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
    LOG("NuclData", pINFO) << "INCLNucleus late initialization";
    //    static INCLNucleus::Cleaner cleaner;
    //    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new INCLNucleus;
    fInstance->theConfig_ = new G4INCL::Config();
    fInstance->nucleus_ = nullptr;
  }
  return fInstance;
}


void INCLNucleus::init(){
//  G4INCL::Config *theConfig = new G4INCL::Config();
//  theConfig->init();
//  G4INCL::ParticleSpecies targetSpecies = G4INCL::ParticleSpecies("Ar40");
//
//  cout << targetSpecies.theA << endl;
//  cout << targetSpecies.theZ << endl;
//  //	cout << targetSpecies.theS << endl;
//  theConfig->setTargetA(targetSpecies.theA);
//  theConfig->setTargetZ(targetSpecies.theZ);
//  //	theConfig->setTargetS(targetSpecies.theS);
//  theConfig->setINCLXXDataFilePath("/root/inclxx/inclxx-v6.33.1-e5857a1/data");
//
//  //	cout << theConfig->summary() << endl;
//  //	cout << theConfig->getTargetA() << endl;
//  //	cout << theConfig->getProjectileType() << endl;
//  //	cout << theConfig->getINCLXXDataFilePath() << endl;
//
//
//  // initialize INCL model
//  G4INCL::Random::initialize(theConfig);
//  G4INCL::Clustering::initialize(theConfig);
//  G4INCL::ParticleTable::initialize(theConfig);
//  //	Cluster cc(18, 40, 0, true);
//  //	cc.initializeParticles();
//  //	cout << cc.print() << endl;
//
//  // define a new argon Nucleus and initialize it
//  G4INCL::Nucleus *p = new G4INCL::Nucleus(40, 18, 0, theConfig, G4INCL::ParticleTable::getMaximumNuclearRadius(G4INCL::Proton, 40, 18), G4INCL::NType);
//  p->initializeParticles();
//  cout << p->getStore()->getParticles().at(2)->getPotentialEnergy() << endl;
//  cout << p->getStore()->getParticles().at(2)->getEnergy() -  p->getStore()->getParticles().at(2)->getPotentialEnergy()  << endl;
//  cout << p->getStore()->getParticles().at(2)->getMomentum().print() << endl;
//  cout << p->getStore()->getParticles().at(2)->getPosition().print() << endl;
//
//  // reset the nucleus
//  p->deleteParticles();
//  p->getStore()->clear();
//  p->initializeParticles();
//
//  //	cout << p->print() << endl;
//  cout << p->getStore()->getParticles().at(2)->getPotentialEnergy() << endl;
//  cout << p->getStore()->getParticles().at(2)->getEnergy() -  p->getStore()->getParticles().at(2)->getPotentialEnergy()  << endl;
//  cout << p->getStore()->getParticles().at(2)->getMomentum().print() << endl;
//  cout << p->getStore()->getParticles().at(2)->getPosition().print() << endl;
//  cout << p->computeTotalEnergy() << endl;
//


}


void INCLNucleus::initialize(const GHepRecord * evrec){

  theConfig_->init();
  GHepParticle * nucleus = evrec->TargetNucleus();
  G4INCL::ParticleSpecies targetSpecies = G4INCL::ParticleSpecies(nucleus->Name());
  theConfig_->setTargetA(targetSpecies.theA);
  theConfig_->setTargetZ(targetSpecies.theZ);
  theConfig_->setTargetS(targetSpecies.theS);
  theConfig_->setINCLXXDataFilePath("/root/inclxx/inclxx-v6.33.1-e5857a1/data"); // FIXME:: using config to set path
 
  // initialize INCL model
  G4INCL::Random::initialize(theConfig_);
  G4INCL::Clustering::initialize(theConfig_);
  G4INCL::ParticleTable::initialize(theConfig_);
  //	Cluster cc(18, 40, 0, true);
  //	cc.initializeParticles();
  //	cout << cc.print() << endl;

  // define a new argon Nucleus and initialize it
 
  if(nucleus_){
    nucleus_->deleteParticles();
    nucleus_->getStore()->clear();
    delete nucleus_;
  }
  // FIXME the last two parameters need to be configed
  // theConfig_, G4INCL::ParticleTable::getMaximumNuclearRadius(G4INCL::Proton, targetSpecies.theA, targetSpecies.theZ)
  // G4INCL::NType
  this->nucleus_ = new G4INCL::Nucleus(targetSpecies.theA, 
      targetSpecies.theZ, 
      targetSpecies.theS, theConfig_, G4INCL::ParticleTable::getMaximumNuclearRadius(G4INCL::Proton, targetSpecies.theA, targetSpecies.theZ), G4INCL::NType);
  nucleus_->initializeParticles();

  LOG("INCLNucleus", pDEBUG) << nucleus_->getStore()->getParticles().at(2)->getPotentialEnergy() ;
  LOG("INCLNucleus", pDEBUG) << nucleus_->getStore()->getParticles().at(2)->getEnergy() -  nucleus_->getStore()->getParticles().at(2)->getPotentialEnergy()  ;
  LOG("INCLNucleus", pDEBUG) << nucleus_->getStore()->getParticles().at(2)->getMomentum().print() ;
  LOG("INCLNucleus", pDEBUG) << nucleus_->getStore()->getParticles().at(2)->getPosition().print() ;

  GHepParticle * nucleon = evrec->HitNucleon();
  LOG("INCLNucleus", pDEBUG) << "hit nucleon pdg : " << nucleon->Pdg();

  RandomGen * rnd = RandomGen::Instance();
  if(pdg::IsProton(nucleon->Pdg()))
    nucleon_index = rnd->RndGen().Integer(targetSpecies.theZ);
  else if(pdg::IsNeutron(nucleon->Pdg()))
    nucleon_index = rnd->RndGen().Integer(targetSpecies.theA - targetSpecies.theZ) + targetSpecies.theZ;
  else
    exit(1);
}

void INCLNucleus::reset(){
  // reset the nucleus
  if(nucleus_){
  	nucleus_->deleteParticles();
  	nucleus_->getStore()->clear();
  	nucleus_->initializeParticles();
  }
}
TVector3 INCLNucleus::getPosition(){
  if(nucleus_ && nucleon_index != -1){
    TVector3 v3_(999999.,999999.,999999.);
    v3_.SetXYZ(nucleus_->getStore()->getParticles().at(nucleon_index)->getPosition().getX(),
	nucleus_->getStore()->getParticles().at(nucleon_index)->getPosition().getY(),
	nucleus_->getStore()->getParticles().at(nucleon_index)->getPosition().getZ());
    return v3_;
  }
  else 
    exit(1);
}
TVector3 INCLNucleus::getMomentum(){
  if(nucleus_ && nucleon_index != -1){
    TVector3 v3_(999999.,999999.,999999.);
    v3_.SetXYZ(nucleus_->getStore()->getParticles().at(nucleon_index)->getMomentum().getX(),
	nucleus_->getStore()->getParticles().at(nucleon_index)->getMomentum().getY(),
	nucleus_->getStore()->getParticles().at(nucleon_index)->getMomentum().getZ());
    return v3_;
  }
  else 
    exit(1);

}
double INCLNucleus::getEnergy(){
  if(nucleus_ && nucleon_index != -1){
  return nucleus_->getStore()->getParticles().at(nucleon_index)->getEnergy() -  nucleus_->getStore()->getParticles().at(nucleon_index)->getPotentialEnergy();
  }
  else 
    exit(1);
}

#endif // __GENIE_INCL_ENABLED__

