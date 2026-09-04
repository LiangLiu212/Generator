//____________________________________________________________________________
/*!

  \class    genie::NucleusGenerator

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

#include "Physics/NuclearState/NucleusGenI.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Interaction/Interaction.h"

namespace genie {

  NucleusGenI::NucleusGenI(string name ) :
    EventRecordVisitorI( name )
  {

  }
  //___________________________________________________________________________
  NucleusGenI::NucleusGenI(string name, string config) :
    EventRecordVisitorI( name, config)
  {

  }
  //___________________________________________________________________________
  NucleusGenI::~NucleusGenI()
  {

  }
  //____________________________________________________________________________
  void NucleusGenI::LoadConfig(void)
  {

  }
  //____________________________________________________________________________
  void NucleusGenI::SetRecordHitNucleon(GHepRecord * event_rec, const Interaction & interaction) const
  {
    GHepParticle * nucleon = event_rec->HitNucleon();
    if(!nucleon) return;
    const TLorentzVector p4 = interaction.InitState().Tgt().HitNucP4();
    nucleon->SetMomentum(p4);
    nucleon->SetRemovalEnergy(nucleon->Mass() - p4.E());
  }
  //____________________________________________________________________________
}

