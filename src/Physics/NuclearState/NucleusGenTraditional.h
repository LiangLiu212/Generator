//____________________________________________________________________________
/*!

  \class    genie::NucleusGenTraditional

  \brief    It visits the event record & combines the FermoMover and VertexGenerator
  computes a Fermi motion momentum and position for initial state nucleons 
  bound in nuclei. It is configured in NucleusGenerator. This new interface
  is compatible with previous implments.
  Is a concrete implementation of the EventRecordVisitorI interface.

  \author   Liang Liu <liangliu \at fnal.gov>
  Fermi National Accelerator Laboratory

  \created  Jan 17, 2025

  \cpright  Copyright (c) 2003-2024, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________


#ifndef _NUCLEUS_GEN_TRADITIONAL_H_
#define _NUCLEUS_GEN_TRADITIONAL_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Physics/NuclearState/NucleusGenI.h"
#include "Framework/GHEP/GHepParticle.h"

namespace genie {


  class NucleusGenTraditional : public NucleusGenI {

    public :
      NucleusGenTraditional();
      NucleusGenTraditional(string config);
      ~NucleusGenTraditional();

      //-- implement the EventRecordVisitorI interface
      void ProcessEventRecord(GHepRecord * event_rec) const;


      void GenerateVertex(GHepRecord *event_rec) const;
      void GenerateCluster(GHepRecord *event_rec) const;
      void setInitialStateVertex   (GHepRecord * evrec) const{
	// TODO: do noting for traditional nuclear model
      }
      void BindHitNucleon() const;
      void BindHitNucleon(Interaction& interaction, double& Eb, QELEvGen_BindingMode_t hitNucleonBindingMode) const;
      void GenerateNucleon(Interaction* interaction, ResamplingHitNucleon_t resampling_mode) const;
      bool isRPValid(double r, double p, const Target & tgt) const;
      void SetHitNucleonOnShellMom(TVector3 p3) const;

      //-- overload the Algorithm::Configure() methods to load private data
      //   members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);

    private:

      void LoadConfig (void);
      // 
      const EventRecordVisitorI *fFermiMover;
      const EventRecordVisitorI *fVertexGenerator;
  };

}      // genie namespace

#endif // _NUCLEUS_GEN_TRADITIONAL_H_
