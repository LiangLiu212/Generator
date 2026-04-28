/*
 * =====================================================================================
 *
 *       Filename:  G4INCLGENIECascadeAction.h
 *
 *    Description: : 
 *
 *        Version:  1.0
 *        Created:  04/22/2026 04:11:48 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Liang Liu (L. Liu), liangliu@fnal.gov
 *		    Fermi National Accelerator Laboratory
 *  Collaboration:  GENIE
 *
 * =====================================================================================
 */

#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#ifndef G4INCLGENIECASCADEACTION_HH
#define G4INCLGENIECASCADEACTION_HH 1

#include "G4INCLCascadeAction.hh"
#include <fstream>

#include "Framework/GHEP/GHepRecord.h"
#include "G4INCLGENIEParticleRecord.h"

namespace G4INCL {

  class GENIECascadeAction : public CascadeAction {

    // class CascadeAction make class INCL a friend because it needs to call private methods
    // but for GENIE interaction, we don't want to change the code inside G4INCL::INCL yet
    // Althrough GENIECascadeAction is inherited from CascadeAction, I only using the virtual 
    // functions.

    public:
      GENIECascadeAction();
      virtual ~GENIECascadeAction();

      virtual void beforeRunUserAction(Config const *);
      virtual void beforeCascadeUserAction(IPropagationModel *);
      virtual void beforePropagationUserAction(IPropagationModel *);
      virtual void beforeAvatarUserAction(IAvatar *, Nucleus *);
      virtual void beforeNPVAvatarUserAction();
      virtual void afterAvatarUserAction(IAvatar *, Nucleus *, FinalState *);
      virtual void afterNPVAvatarUserAction(IAvatar *, Nucleus *, FinalState *);
      virtual void afterPropagationUserAction(IPropagationModel *, IAvatar *);
      virtual void afterCascadeUserAction(Nucleus *);
      virtual void afterRunUserAction();

      void setGHepRecord(genie::GHepRecord *evr){
        evrec = evr;
      }

      // TODO: void file steps
      //
      // TODO: void fill final states

    private:
      //std::ofstream *oFile;
      long eventCounter;
      long stepCounter;

      genie::GHepRecord * evrec;

      ParticleList backup_mother;
      std::vector<G4INCL::GENIEParticleRecord> eventRecord;


      struct INCLRecord{
        int global_index;        // Each particles in INCLXX will have a unique ID, it is a global index for every simulation run.
        int pdgid;	       // PDG ID of particles in INCLXX
        int mother_index;        // mother index of particles in each event
        int local_index;         // local index of particles in each event
        TLorentzVector p4mom;
        TLorentzVector p4posi;
        G4INCL::ParticleType theType;   // INCL Particle type

        INCLRecord(int g_id, int p_id, int m_id, int l_id):
          global_index(g_id),
          pdgid(p_id),
          mother_index(m_id),
          local_index(l_id){}
        INCLRecord(int g_id, int p_id, int m_id, int l_id, TLorentzVector mom, TLorentzVector posi):
          global_index(g_id),
          pdgid(p_id),
          mother_index(m_id),
          local_index(l_id),
          p4mom(mom),
          p4posi(posi){}
        INCLRecord(int g_id, int p_id, int m_id, int l_id, TLorentzVector mom, TLorentzVector posi, G4INCL::ParticleType pType):
          global_index(g_id),
          pdgid(p_id),
          mother_index(m_id),
          local_index(l_id),
          p4mom(mom),
          p4posi(posi), 
          theType(pType){}
      };
      std::vector<INCLRecord> tempFinalState;


  };

}
#endif // G4INCLGENIECASCADEACTION_HH

#endif // __GENIE_INCL_ENABLED__
