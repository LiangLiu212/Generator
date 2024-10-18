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

#ifndef _INCL_NUCLEUS_H_
#define _INCL_NUCLEUS_H_

// ROOT
#include "TVector3.h"


// INCL++
// For configuration
#include "G4INCLConfig.hh"
#include "G4INCLNucleus.hh"
// GENIE
#include "Framework/GHEP/GHepRecord.h"

namespace genie {

  class INCLNucleus {

    public: 
      static INCLNucleus * Instance (void);
      void init();
      void initialize(const GHepRecord * evrec);
      void reset();

      TVector3 getPosition();
      TVector3 getMomentum();
      double   getEnergy();

    private:
      INCLNucleus();
      ~INCLNucleus();

      static INCLNucleus *fInstance;

      TVector3 v3_; // position of initial nucleon 
      TVector3 p3_; // fermi momentum of initial nucleon
      double  energy_; // off-shell energy of initial nucleon
      G4INCL::Config *theConfig_;
      G4INCL::Nucleus *nucleus_;

      static int nucleon_index;


  };

}         // genie namespace
#endif    // _INCL_NUCLEUS_H_
#endif // __GENIE_INCL_ENABLED__
