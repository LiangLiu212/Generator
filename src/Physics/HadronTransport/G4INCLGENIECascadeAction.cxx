/*
 * =====================================================================================
 *
 *       Filename:  G4INCLGENIECascadeAction.cxx
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/22/2026 04:11:54 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Liang Liu (L. Liu), liangliu@fnal.gov
 *		    Fermi National Accelerator Laboratory
 *  Collaboration:  GENIE
 *
 * =====================================================================================
 */

#include "G4INCLGENIECascadeAction.h"
#include "G4INCLGENIEParticleRecord.h"
#include "Physics/NuclearState/INCLNucleus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/ProcessInfo.h"
#include <sstream>
#include <string>

namespace G4INCL {
  using namespace genie;

  GENIECascadeAction::GENIECascadeAction()
  {
  }

  GENIECascadeAction::~GENIECascadeAction() {}

  void GENIECascadeAction::beforeRunUserAction(Config const *){
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
  }

  void GENIECascadeAction::beforeCascadeUserAction(IPropagationModel * /*pm*/) {

  }

  void GENIECascadeAction::beforePropagationUserAction(IPropagationModel *){
  }
  void GENIECascadeAction::beforeAvatarUserAction(IAvatar *, Nucleus *) { }

  void GENIECascadeAction::beforeNPVAvatarUserAction() {
    const ProcessInfo & proc_info = evrec->Summary()->ProcInfo();
    INCLNucleus *incl_nucleus = INCLNucleus::Instance();

    // convert ghep event record to INCL Style.
    // G4INCL::GENIEParticleRecord is the bridge
    TObjArrayIter piter(evrec);
    GHepParticle * p = nullptr;
    eventRecord.clear();
    tempFinalState.clear();
    while ( (p = (GHepParticle *) piter.Next() ) ) {
      // the code of the particles in primary neutrino interaction
      G4INCL::GENIERecordCode recordCode;
      if(eventRecord.size() == evrec->ProbePosition())                         { recordCode = G4INCL::kProbe; }
      else if(eventRecord.size() == evrec->TargetNucleusPosition())            { recordCode = G4INCL::kTarget;}
      else if(eventRecord.size() == evrec->HitNucleonPosition())               { recordCode = G4INCL::kHitNucleon;}
      else if(eventRecord.size() == evrec->RemnantNucleusPosition())           { recordCode = G4INCL::kRemnant;}
      else if(eventRecord.size() == evrec->FinalStatePrimaryLeptonPosition())  { recordCode = G4INCL::kFinalStateLepton;}
      else { recordCode = G4INCL::kUnknown;}
      eventRecord.emplace_back(p, int(proc_info.ScatteringTypeId()), recordCode);
    }
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << "  " << eventRecord.size() << std::endl;

  }

  void GENIECascadeAction::afterAvatarUserAction(IAvatar *avatar, Nucleus *nucleus, FinalState *finalState) { }
  void GENIECascadeAction::afterNPVAvatarUserAction(IAvatar *avatar, Nucleus *nucleus, FinalState *finalState) { 
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;


  // update the event record after INCL postInteraction
  // INCL might rescale the four momentum of final states
  // particles, we update p4 in GHep event record

  TObjArrayIter piter(evrec);
  piter.Reset();   // rewind
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
  GHepParticle * p = nullptr;
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << "  " << eventRecord.size() << std::endl;
  auto er = eventRecord.begin();
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
  int idx =0;
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
  while ( (p = (GHepParticle *) piter.Next() ) ) {
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
    TLorentzVector *p4 = p->P4();
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << er->P3().print() << std::endl;
    p4->SetPx(er->P3().getX()/1000.);
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
    p4->SetPy(er->P3().getY()/1000.);
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
    p4->SetPz(er->P3().getZ()/1000.);
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
    p4->SetE(std::sqrt(er->P3().mag2() + er->Mass()*er->Mass())/1000.);
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
    tempFinalState.emplace_back(er->ID(), er->Pdg(), er->FirstMother(), idx++);
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
    er++;
  }
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;

  // std::vector<GENIEParticleRecord> *eventRecord = avatar->getEventRecord();
  ParticleList outgoing = finalState->getOutgoingParticles();
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
  for(ParticleIter iter=outgoing.begin(); iter!=outgoing.end(); ++iter){
    int outp_mother_idx = -1;
    int tmp_idx_ = -1;
    int pdg = 0;
    GHepParticle * p1 = nullptr;
    for(auto er = eventRecord.begin(); er != eventRecord.end(); er++){
      tmp_idx_++;
      if((*iter)->getID() == er->ID()){
        outp_mother_idx = tmp_idx_;
        pdg = er->Pdg();
      }
    }
    GHepParticle p(pdg, kIStStableFinalState, outp_mother_idx, -1, -1, -1, 
        TLorentzVector((*iter)->getMomentum().getX() / 1000,
          (*iter)->getMomentum().getY() / 1000,
          (*iter)->getMomentum().getZ() / 1000,
          (*iter)->getEnergy() / 1000),
        TLorentzVector((*iter)->getPosition().getX(),
          (*iter)->getPosition().getY(),
          (*iter)->getPosition().getZ(),
          0)
        );
    evrec->AddParticle(p);
    tempFinalState.emplace_back((*iter)->getID(), pdg, outp_mother_idx, idx++);
  }
  //evrec->Print(std::cout);
    
  }
  void GENIECascadeAction::afterPropagationUserAction(IPropagationModel *, IAvatar *) {
  }

  void GENIECascadeAction::afterCascadeUserAction(Nucleus * /*nucleus*/) {
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
    evrec->Print(std::cout);

  }
  void GENIECascadeAction::afterRunUserAction() {
    std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << std::endl;
    evrec->Print(std::cout);
  }

}
