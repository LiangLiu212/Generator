//____________________________________________________________________________
/*!

  \class    genie::NucleusGenI

  \brief    Interface of nucleus generator. It combines the 
            vertex generator and fermi mover for triditional 
	    nuclear model.
	    INCL nuclear model will provide a nucleon with both
	    position and momnetum

  \author   Liang Liu <liangliu \at fnal.gov>
  Fermi National Accelerator Laboratory

  \created  October 17, 2024

  \cpright  Copyright (c) 2003-2024, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#ifndef _NUCLEUS_GENERATOR_INTERFACE_H_
#define _NUCLEUS_GENERATOR_INTERFACE_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/Interaction/Target.h"

namespace genie {

  // A nucleus generator that will handle fermi motion and 
  // position simultaneously.
  //
  // 1. class NucleusGenTraditional combine the FermiMover 
  //    and VertexGenerator
  // 2. class NucleusGenINCL uses external INCLXX nuclear
  //    model.

  class NucleusGenI : public EventRecordVisitorI {

    public :
      virtual ~NucleusGenI();
      virtual const NuclearModelI* GetNuclearModel() const{
	return fNuclModel;
      }
      virtual void GenerateVertex(GHepRecord *event_rec) const  = 0;
      virtual void GenerateCluster(GHepRecord *event_rec) const  = 0;
      virtual void setInitialStateVertex   (GHepRecord * event_rec) const = 0; // ///< give hit nucleon a position

    protected:
      void LoadConfig (void);
      NucleusGenI(string name);
      NucleusGenI(string name, string config);
      // 
      const NuclearModelI *  fNuclModel;   ///< nuclear model
  };
}      // genie namespace
#endif // _NUCLEUS_GENERATOR_H_
