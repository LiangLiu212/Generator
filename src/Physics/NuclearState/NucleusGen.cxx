//____________________________________________________________________________
/*!

\class    genie::NucleusGen

\brief    It visits the event record & computes a Fermi motion momentum for
          initial state nucleons bound in nuclei.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Liang Liu <liangliu \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  October 17, 2024

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__
#include <cstdlib>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TParticlePDG.h>
#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Physics/NuclearState/NucleusGen.h"

#include "Physics/NuclearState/NuclearModel.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
NucleusGen::NucleusGen() :
EventRecordVisitorI("genie::NucleusGen")
{

}
//___________________________________________________________________________
NucleusGen::NucleusGen(string config) :
EventRecordVisitorI("genie::NucleusGen", config)
{

}
//___________________________________________________________________________
NucleusGen::~NucleusGen()
{

}

//___________________________________________________________________________
void NucleusGen::ProcessEventRecord(GHepRecord * evrec) const
{
  // skip if not a nuclear target
  if(! evrec->Summary()->InitState().Tgt().IsNucleus()) return;

  // skip if no hit nucleon is set
  if(! evrec->HitNucleon()) return;


  // give hit nucleon a vertex
  this->setInitialStateVertex(evrec);
  // give hit nucleon a Fermi momentum
  this->setInitialStateMomentum(evrec);

  // handle the addition of the recoil nucleon
  // TODO:  INCL has it own SRC model
//  if ( fSecondEmitter ) fSecondEmitter -> ProcessEventRecord( evrec ) ;

  // add a recoiled nucleus remnant
  this->setTargetNucleusRemnant(evrec);
}

//___________________________________________________________________________
//  using INCL model to get the position and momentum of 
//  Hit  nucleon
void NucleusGen::setInitialStateVertex(GHepRecord * evrec) const{

}

void NucleusGen::setInitialStateMomentum(GHepRecord * evrec) const{

}

void NucleusGen::setTargetNucleusRemnant(GHepRecord * evrec)const{

}


//___________________________________________________________________________
void NucleusGen::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NucleusGen::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NucleusGen::LoadConfig(void)
{

//  RgKey nuclkey = "NuclearModel";
//  fNuclModel = 0;
//  fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
//  assert(fNuclModel);

//  this->GetParamDef("KeepHitNuclOnMassShell", fKeepNuclOnMassShell, false);
//
//  bool mom_dep_energy_removal_def = false;
//  this->GetParamDef("LFG-MomentumDependentErmv", mom_dep_energy_removal_def, false ) ;
//  // it defaults to whatever the nuclear model sets. Since only the LFG has this option
//  // this simple search is enough.
//
//  this->GetParamDef("MomentumDependentErmv", fMomDepErmv, mom_dep_energy_removal_def);
//
//  RgKey nuclearrecoilkey = "SecondNucleonEmitter" ;
//  fSecondEmitter = dynamic_cast<const SecondNucleonEmissionI *> (this->SubAlg(nuclearrecoilkey));

}
//____________________________________________________________________________


#endif // end  __GENIE_INCL_ENABLED__
