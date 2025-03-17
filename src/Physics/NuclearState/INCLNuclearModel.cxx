 //____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

 @ Mar 18, 2016- Joe Johnston (SD)
   Update GenerateNucleon() and Prob() to accept a radius as the argument,
   and call the corresponding methods in the nuclear model with a radius.

*/
//____________________________________________________________________________

#include <sstream>

#include <TSystem.h>
#include <TNtupleD.h>
#include <TGraph2D.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Numerical/RandomGen.h"
#include "Physics/NuclearState/INCLNuclearModel.h"
#include "Physics/NuclearState/INCLNucleus.h"

using std::ostringstream;
using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________
INCLNuclearModel::INCLNuclearModel() :
NuclearModelI("genie::INCLNuclearModel")
{

}
//____________________________________________________________________________
INCLNuclearModel::INCLNuclearModel(string config) :
NuclearModelI("genie::INCLNuclearModel", config)
{

}
//____________________________________________________________________________
INCLNuclearModel::~INCLNuclearModel()
{

}
//____________________________________________________________________________
bool INCLNuclearModel::GenerateNucleon(const Target & target,
                                      double hitNucleonRadius) const
{

//  LOG("NucleusGenINCL", pINFO) << "Initialize a nucleus for a new event!";
  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  incl_nucleus->initialize(&target);
  incl_nucleus->reset(&target);
  incl_nucleus->initialize(&target);
  G4INCL::Particle *incl_nucleon = incl_nucleus->getHitParticle();

//  LOG("NucleusGenINCL", pINFO) << "Removal Energy " << incl_nucleus->getRemovalEnergy() / 1000;

  fCurrRemovalEnergy =  incl_nucleus->getRemovalEnergy() / 1000;
  TVector3 hit_mom = incl_nucleus->getHitNucleonMomentum();
  fCurrMomentum.SetXYZ(hit_mom.X()/1000, hit_mom.Y()/1000, hit_mom.Z()/1000);
  return true;
}
//____________________________________________________________________________
double INCLNuclearModel::Prob(double p, double w, const Target & target,
                             double hitNucRadius) const
{
  if(p < this->FermiMomentum(target, target.HitNucPdg()))
    return 1.0;
  else 
    return -1.0;
}
//____________________________________________________________________________
NuclearModel_t INCLNuclearModel::ModelType(const Target & target) const
{

}
//____________________________________________________________________________
double INCLNuclearModel::FermiMomentum( const Target &  t, int nucleon_pdg ) const {

  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  incl_nucleus->initialize(&t);
  G4INCL::Nucleus *incl_tgt = incl_nucleus->getNuclues();

  if(pdg::IsProton(nucleon_pdg))
    return incl_tgt->getPotential()->getFermiMomentum(G4INCL::Proton);
  else 
    return incl_tgt->getPotential()->getFermiMomentum(G4INCL::Neutron);
}
//____________________________________________________________________________
double INCLNuclearModel::LocalFermiMomentum( const Target & t, 
					    int nucleon_pdg, double /*radius*/ ) const {
  return this->FermiMomentum(t, nucleon_pdg);
}
//____________________________________________________________________________
void INCLNuclearModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void INCLNuclearModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void INCLNuclearModel::LoadConfig(void)
{
/*
  std::string inclxxpath;
  GetParamDef( "inclxx-data-dir", inclxxpath, std::string(""));
  LOG("NucleusGenINCL", pINFO) << this->expandEnvironmentPath(inclxxpath);

  std::string ablaxxpath;
  GetParamDef( "inclxx-ablaxx-dir", ablaxxpath, std::string(""));
  LOG("NucleusGenINCL", pINFO) << this->expandEnvironmentPath(ablaxxpath);

  std::string abla07path;
  GetParamDef( "inclxx-abla07-dir", abla07path, std::string(""));
  LOG("NucleusGenINCL", pINFO) << this->expandEnvironmentPath(abla07path);

  std::string geminixxpath;
  GetParamDef( "inclxx-geminixx-dir", geminixxpath, std::string(""));
  LOG("NucleusGenINCL", pINFO) << this->expandEnvironmentPath(geminixxpath);

  std::string deExType;
  G4INCL::DeExcitationType deExcitationType;
  GetParamDef( "inclxx-de-excitation", deExType, std::string(""));
  LOG("NucleusGenINCL", pINFO) << "inclxx-de-excitation : " << deExType;
  if(deExType.compare("ABLA07")){
    deExcitationType = G4INCL::DeExcitationABLA07;
  } else if(deExType.compare("ABLAXX")) {
    deExcitationType = G4INCL::DeExcitationABLAXX;
  } else if(deExType.compare("GEMINIXX")) {
    deExcitationType = G4INCL::DeExcitationGEMINIXX;
  } else {
    std::stringstream ss;
    ss<< "########################################################\n"
      << "###              WARNING WARNING WARNING             ###\n"
      << "###                                                  ###\n"
      << "### You are running the code without any coupling to ###\n"
      << "###              a de-excitation model!              ###\n"
      << "###    Results will be INCOMPLETE and UNPHYSICAL!    ###\n"
      << "###    Are you sure this is what you want to do?     ###\n"
      << "########################################################\n";
    LOG("NucleusGenINCL", pWARN) << '\n' << ss.str();
  }
  incl_nucleus->setINCLXXDataFilePath(this->expandEnvironmentPath(inclxxpath));
  incl_nucleus->setABLAXXDataFilePath(this->expandEnvironmentPath(ablaxxpath));
  incl_nucleus->setABLA07DataFilePath(this->expandEnvironmentPath(abla07path));
  incl_nucleus->setGEMINIXXDataFilePath(this->expandEnvironmentPath(geminixxpath));
  incl_nucleus->setDeExcitationType(deExcitationType);
  incl_nucleus->configure();
  */

 // INCLNucleus *incl_nucleus = INCLNucleus::Instance();

}

//____________________________________________________________________________
