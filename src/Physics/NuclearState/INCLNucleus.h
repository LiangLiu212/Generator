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
#include "Framework/Interaction/Target.h"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIPropagationModel.hh"
#include "G4INCLStandardPropagationModel.hh"
#include "G4INCLCascadeAction.hh"
#include "G4INCLEventInfo.hh"
#include "G4INCLGlobalInfo.hh"
#include "G4INCLLogger.hh"
#include "G4INCLConfig.hh"
#include "G4INCLRootFinder.hh"




namespace genie {

  class INCLNucleus {

    public: 
      static INCLNucleus * Instance (void);
      void initialize(const Target * tgt);
      void reset(const Target * tgt);
      void configure();

      TVector3 getHitNucleonPosition();
      TVector3 getHitNucleonMomentum();
      double   getHitNucleonEnergy();
      double   getHitNucleonMass();
      double   getMass();
      double   getRemovalEnergy();
      G4INCL::Nucleus * getNuclues();
      G4INCL::StandardPropagationModel * getPropagationModel();
      G4INCL::Particle *getHitParticle();
      G4INCL::Config *getConfig(){return theConfig_;}

      double getMaxUniverseRadius() {return maxUniverseRadius_;}

      void setINCLXXDataFilePath(std::string str){ INCLXXDataFilePath_ = str; }
      void setABLAXXDataFilePath(std::string str){ ablaxxDataFilePath_ = str; }
      void setABLA07DataFilePath(std::string str){ abla07DataFilePath_ = str; }
      void setGEMINIXXDataFilePath(std::string str){ geminixxDataFilePath_ = str; }
      void setDeExcitationType(G4INCL::DeExcitationType deExType){ deExcitationType_ = deExType; }

    private:
      INCLNucleus();
      ~INCLNucleus();
      void initUniverseRadius(const int A, const int Z);
      static INCLNucleus *fInstance;

      //TVector3 v3_; // position of initial nucleon 
      //TVector3 p3_; // fermi momentum of initial nucleon
      double  energy_; // off-shell energy of initial nucleon
      G4INCL::Config *theConfig_;
      G4INCL::Nucleus *nucleus_;
      G4INCL::StandardPropagationModel *propagationModel_;
      G4INCL::CascadeAction *cascadeAction_;
      G4INCL::Particle *hitNucleon_;

      int nucleon_index_;

      double maxUniverseRadius_;
      // double maxInteractionDistance_;
      double minRemnantSize_;

      std::string INCLXXDataFilePath_;
      std::string abla07DataFilePath_;
      std::string ablaxxDataFilePath_;
      std::string geminixxDataFilePath_;

      G4INCL::DeExcitationType deExcitationType_;

  };

}         // genie namespace
#endif    // _INCL_NUCLEUS_H_
#endif // __GENIE_INCL_ENABLED__
