#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

//---------------------
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

// GENIE
#include "INCLCascadeIntranuke.h"

// INCL++
#include "G4INCLConfig.hh"
#include "G4INCLCascade.hh"
#include "G4INCLConfigEnums.hh"
#include "G4INCLParticle.hh"
// signal handler (for Linux and GCC)
#include "G4INCLSignalHandling.hh"


#include "G4INCLPauliBlocking.hh"

#include "G4INCLCrossSections.hh"

#include "G4INCLPhaseSpaceGenerator.hh"

#include "G4INCLLogger.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLNuclearDensityFactory.hh"

#include "G4INCLINuclearPotential.hh"

#include "G4INCLCoulombDistortion.hh"

#include "G4INCLClustering.hh"

#include "G4INCLIntersection.hh"

#include "G4INCLBinaryCollisionAvatar.hh"

#include "G4INCLCascadeAction.hh"
#include "G4INCLAvatarDumpAction.hh"

#include "Framework/GHEP/GHepParticle.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"




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

//#ifdef HAS_BOOST_DATE_TIME
//#include <boost/date_time/posix_time/posix_time.hpp>
//namespace bpt = boost::posix_time;
//#endif

#ifdef HAS_BOOST_TIMER
#include <boost/timer/timer.hpp>
namespace bt = boost::timer;
#endif


// --------------------------------------Include for GENIE---------------------
// GENIE
#include "INCLConvertParticle.hh"
#include "INCLConfigParser.h"
// INCL++
//#include "ConfigParser.hh"

#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/SystemUtils.h"
#include "Physics/NuclearState/INCLNucleus.h"

// ROOT
#include "TSystem.h"

using namespace genie;
using namespace genie::utils;
using namespace G4INCL;
using std::ostringstream;
using namespace std;

INCLCascadeIntranuke::INCLCascadeIntranuke() :
EventRecordVisitorI("genie::INCLCascadeIntranuke"),
theINCLConfig(0), theINCLModel(0), theDeExcitation(0)
{
  LOG("INCLCascadeIntranuke", pDEBUG)
  << "default ctor";
}

//______________________________________________________________________________
INCLCascadeIntranuke::INCLCascadeIntranuke(string config) :
EventRecordVisitorI("genie::INCLCascadeIntranuke", config),
theINCLConfig(0), theINCLModel(0), theDeExcitation(0)
{
  LOG("INCLCascadeIntranuke", pDEBUG)
  << "ctor from config " << config;
}

//______________________________________________________________________________
INCLCascadeIntranuke::~INCLCascadeIntranuke()
{

  // Config is owned by model once handed over
  if ( theINCLConfig   ) { theINCLConfig=0;   }
  if ( theINCLModel    ) { delete theINCLModel;    theINCLModel=0;    }
  if ( theDeExcitation ) { delete theDeExcitation; theDeExcitation=0; }

}

//______________________________________________________________________________
void INCLCascadeIntranuke::LoadConfig(void)
{

}

//______________________________________________________________________________
bool INCLCascadeIntranuke::AddDataPathFlags(size_t& nflags, char** flags) {
  return false;

}

//______________________________________________________________________________

bool INCLCascadeIntranuke::LookForAndAddValidPath(std::vector<std::string>& datapaths,
    size_t defaultIndx,
 const char* optString,
 size_t& nflags, char** flags) {

  return false;
}

//______________________________________________________________________________
int INCLCascadeIntranuke::doCascade(GHepRecord * evrec) const {

  if ( ! theINCLConfig || ! theINCLModel ) return 0;
  return 0;
}

void INCLCascadeIntranuke::ProcessEventRecord(GHepRecord * evrec)  const {
    LOG("INCLCascadeIntranuke", pINFO) << "Start with this event";

    // FIXME: start with the NC process
    INCLNucleus *incl_nucleus = INCLNucleus::Instance();
    incl_target = incl_nucleus->getNuclues();
    hit_nucleon = incl_nucleus->getHitParticle();
    theConfig = incl_nucleus->getConfig();
    propagationModel =  incl_nucleus->getPropagationModel();
    incl_target->setParticleNucleusCollision();
    std::unique_ptr<G4INCL::FinalState> finalState(new FinalState);
 //   std::cout << hit_nucleon->print() << std::endl;
 //   std::cout << incl_target->print() << std::endl;
    prob = evrec->Probe();
    primarylepton = evrec->FinalStatePrimaryLepton();

    this->preCascade();

    double currentTime = 0.0;
    double temfin;
    double TLab;

    TObjArrayIter piter(evrec);
    GHepParticle * p = nullptr;

    tempFinalState.clear();

    stepFinalState.clear();
    istep = 0;

    int idx = 0;
    while ( (p = (GHepParticle *) piter.Next() ) ) {
      if(p->Status() == kIStHadronInTheNucleus){
	// FIXME: for NC QE, it is only neutron and proton
	temfin = 29.8 * std::pow(incl_target->getA(), 0.16);
	TLab = p->KinE() * 1000.;
        TLorentzVector *p4= p->P4();
        G4INCL::ThreeVector momentum (0.,0.,0.);
        momentum.setX(p4->Px()*1000);
        momentum.setY(p4->Py()*1000);
        momentum.setZ(p4->Pz()*1000);
        hit_nucleon->setMomentum(momentum);
        hit_nucleon->adjustEnergyFromMomentum();
        hit_nucleon->makeParticipant();
        finalState->addModifiedParticle(hit_nucleon);
	tempFinalState.emplace_back(hit_nucleon->getID(), p->Pdg(), p->FirstMother(), idx);
      }
      else{
	tempFinalState.emplace_back(-1, p->Pdg(), p->FirstMother(), idx);
      }
      idx++;
    }





//    const G4INCL::ThreeVector mom_lep(p4l->P4()->Px() * 1000.,
//	p4l->P4()->Py() * 1000.,
//	p4l->P4()->Pz() * 1000.);
//    const G4INCL::ThreeVector posi_lep(p4l->X4()->X(),
//	p4l->X4()->Y(),
//	p4l->X4()->Z());
//
//    G4INCL::ParticleType leptonType;
//    switch(p4l->Pdg()){
//      case  11: leptonType = G4INCL::Electron;  break;
//      case -11: leptonType = G4INCL::Positron;  break;
//      case  13: leptonType = G4INCL::Muon;      break;
//      case -13: leptonType = G4INCL::antiMuon;  break;
//      case  15: leptonType = G4INCL::Tau;       break;
//      case -15: leptonType = G4INCL::antiTau;   break;
//      case  12: leptonType = G4INCL::NuE;       break;
//      case -12: leptonType = G4INCL::antiNuE;   break;
//      case  14: leptonType = G4INCL::NuMu;      break;
//      case -14: leptonType = G4INCL::antiNuMu;  break;
//      case  16: leptonType = G4INCL::NuTau;     break;
//      case -16: leptonType = G4INCL::antiNuTau; break;
//      default:
//		exit(0);
//		break;
//    }
//
//    std::unique_ptr<G4INCL::Particle> lepton(new G4INCL::Particle(leptonType, mom_lep, posi_lep));
//    lepton->adjustEnergyFromMomentum();
//
//    LOG("INCLCascadeIntranuke", pINFO) << "lepton mass  : " << lepton->getMass();
//    LOG("INCLCascadeIntranuke", pINFO) << "lepton energy: " << lepton->getEnergy();
//
//    LOG("INCLCascadeIntranuke", pINFO) << "lepton " << lepton->print(); 

    //finalState->addOutgoingParticle(lepton.get());


    if(TLab > 2000.)
      temfin *= (5.8E4-TLab)/5.6E4;

    double maximumTime = temfin;

    propagationModel->setStoppingTime(maximumTime);
    propagationModel->setNucleus(incl_target);
    propagationModel->generateAllAvatars();

    // stopping time: 
    // INCL don't have the stopping time for neutrino.
    // we can calculate the longest stopping time for all the daughters from primary interaction.

    incl_target->applyFinalState(finalState.get());
    incl_target->getStore()->getBook().incrementCascading();
    incl_target->getStore()->getBook().incrementAcceptedCollisions();
    int step = 0;

 //   std::cout << finalState->print() << std::endl;



    while(step < 10000000 && continueCascade()){
      IAvatar *avatar = propagationModel->propagate(finalState.get());
      finalState->reset();
      if(avatar == 0) break; // No more avatars in the avatar list.
      LOG("INCLCascadeIntranuke", pINFO) << avatar->toString();
      avatar->fillFinalState(finalState.get());
 //     std::cout << finalState->print() << std::endl;
      incl_target->applyFinalState(finalState.get());
 //     std::cout << incl_target->print() << std::endl;
      LOG("INCLCascadeIntranuke", pINFO) << avatar->dump();
      G4INCL::ParticleList mother_list = avatar->getParticles();
      this->fillEventRecord(finalState.get(), mother_list, evrec, propagationModel->getCurrentTime());
      delete avatar;
      step++;
    }

    //finalState->reset();
    //finalState->addOutgoingParticle(lepton.get());
    //incl_target->applyFinalState(finalState.get());
    //finalState->reset();

 //   std::cout << incl_target->print() << std::endl;

    // put the nuclear remnant in the event record
 
    this->postCascade();


    LOG("INCLCascadeIntranuke", pINFO) << "nRemnants   : " << theEventInfo.nRemnants;
    LOG("INCLCascadeIntranuke", pINFO) << "nRemnants A : " << theEventInfo.ARem[0];
    LOG("INCLCascadeIntranuke", pINFO) << "nRemnants Z : " << theEventInfo.ZRem[0];
    LOG("INCLCascadeIntranuke", pINFO) << "nRemnants S : " << theEventInfo.SRem[0];

    LOG("INCLCascadeIntranuke", pINFO) << "nParticles   : " << theEventInfo.nParticles;

    const int n_outgoing = theEventInfo.nParticles;

    for(int i = 0; i < theEventInfo.nParticles; i++){
      LOG("INCLCascadeIntranuke", pINFO) << "Final State Particles PDG: " << theEventInfo.PDGCode[i];
      LOG("INCLCascadeIntranuke", pINFO) << "Final State Particles px : " << theEventInfo.px[i];
      LOG("INCLCascadeIntranuke", pINFO) << "Final State Particles py : " << theEventInfo.py[i];
      LOG("INCLCascadeIntranuke", pINFO) << "Final State Particles pz : " << theEventInfo.pz[i];
    }


    double Rem_p2 = theEventInfo.pxRem[0]*theEventInfo.pxRem[0]
      + theEventInfo.pyRem[0]*theEventInfo.pyRem[0]
      + theEventInfo.pzRem[0]*theEventInfo.pzRem[0];
    double Rem_Kin = theEventInfo.EKinRem[0];
    double Rem_M = (Rem_p2 - Rem_Kin*Rem_Kin) / ( 2.0 * Rem_Kin);
    double Rem_E =  std::sqrt(Rem_p2 + Rem_M*Rem_M);

    // GHepParticle *p4l = evrec->FinalStatePrimaryLepton();

    double Rem_px = theEventInfo.pxRem[0]/1000.;
    double Rem_py = theEventInfo.pyRem[0]/1000.;
    double Rem_pz = theEventInfo.pzRem[0]/1000.;

    TLorentzVector p4mom(Rem_px, Rem_py, Rem_pz, Rem_E);
    TLorentzVector p4posi(0,0,0,0);
    int A = incl_target->getA();
    int Z = incl_target->getZ();
    int pdg = genie::pdg::IonPdgCode( A , Z );

    evrec->AddParticle(pdg, kIStPreDeExNuclearRemnant, 3, -1, -1, -1, p4mom, p4posi);

    std::unique_ptr<G4INCL::IDeExcitation> theDeExcitation = std::make_unique<G4INCLAblaInterface>(theConfig);
    //theDeExcitation = new G4INCLAblaInterface(theConfig);
    theDeExcitation->deExcite(&theEventInfo);

 //   evrec->Print(std::cout);

    for(int i = 0; i < theEventInfo.nParticles; i++){
      LOG("INCLCascadeIntranuke", pINFO) << "Final State Particles PDG: " << theEventInfo.PDGCode[i];
      LOG("INCLCascadeIntranuke", pINFO) << "Final State Particles px : " << theEventInfo.px[i];
      LOG("INCLCascadeIntranuke", pINFO) << "Final State Particles py : " << theEventInfo.py[i];
      LOG("INCLCascadeIntranuke", pINFO) << "Final State Particles pz : " << theEventInfo.pz[i];
    }

    // particle from de-excitation
    const int remnant_id = evrec->GetEntries() - 1;
    for(int i = n_outgoing; i < theEventInfo.nParticles; i++){
      LOG("INCLCascadeIntranuke", pINFO) << "Final State Particles PDG: " << theEventInfo.PDGCode[i];
      LOG("INCLCascadeIntranuke", pINFO) << "Final State Particles PDG: " << this->INCLPDG_to_GHEPPDG(theEventInfo.PDGCode[i], theEventInfo.A[i], theEventInfo.Z[i]);
      int depdg = this->INCLPDG_to_GHEPPDG(theEventInfo.PDGCode[i], theEventInfo.A[i], theEventInfo.Z[i]);
      TParticlePDG * p = PDGLibrary::Instance()->Find(depdg);

      double p2 = theEventInfo.px[i]*theEventInfo.px[i]
        + theEventInfo.py[i]*theEventInfo.py[i]
        + theEventInfo.pz[i]*theEventInfo.pz[i];

      double M = p->Mass();

      double E = sqrt(p2/1000000. + M*M);


      TLorentzVector p4mom(theEventInfo.px[i] / 1000.,
	  theEventInfo.py[i] / 1000.,
	  theEventInfo.pz[i] / 1000., 
	  E);
      TLorentzVector p4posi(0,0,0,0);

      EGHepStatus ptype;
      if(theEventInfo.A[i] > 4){
	ptype = kIStFinalStateNuclearRemnant;
      } else {
	ptype = kIStDeExFinalState;
      }
      evrec->AddParticle(depdg, ptype, remnant_id, -1, -1, -1, p4mom, p4posi);
    }

 //   evrec->Print(std::cout);

    // print interaction in each step
    for(std::map<int, std::vector<INCLRecord>>::iterator it = stepFinalState.begin(); it != stepFinalState.end(); it++){
      LOG("INCLCascadeIntranuke", pINFO) << "step : " << it->first;
      for(auto ip = it->second.begin(); ip != it->second.end(); ip++){
	std::string fs_particle;
	switch (ip->mother_index) {
	  case -2: fs_particle = "Modified  particle "; break;
	  case -3: fs_particle = "Outgoing  particle "; break;
	  case -4: fs_particle = "Destoried particle "; break;
	  case -5: fs_particle = "Created   particle "; break;
	  case -6: fs_particle = "Mother    particle "; break;
	  default:;
	}
	LOG("INCLCascadeIntranuke", pINFO)  << "global id : " << ip->global_index << "  " 
	  << "pdg id : " << ip->pdgid << "  "
	  << "mother id : " << ip->mother_index << "  "
	  << "local id : " << ip->local_index;
      }
    }

    LOG("INCLCascadeIntranuke", pINFO) << "Done with this event";
}

bool INCLCascadeIntranuke::CanRescatter(const GHepParticle * p) const {

  return false;
}

void INCLCascadeIntranuke::TransportHadrons(GHepRecord * evrec) const {


}

//______________________________________________________________________________
int INCLCascadeIntranuke::pdgcpiontoA(int pdgc) const {

  if      ( pdgc == 2212 || pdgc == 2112 ) return 1;
  else if ( pdgc ==  211 || pdgc == -211 || pdgc == 111 ) return 0;
  return 0;  // return something

}

//______________________________________________________________________________
int INCLCascadeIntranuke::pdgcpiontoZ(int pdgc) const {

  if      ( pdgc == 2212 || pdgc == 211 ) return 1;
  else if ( pdgc == 2112 || pdgc == 111 ) return 0;
  else if ( pdgc == -211 ) return -1;
  return 0; // return something

}

//______________________________________________________________________________
bool INCLCascadeIntranuke::NeedsRescattering(const GHepParticle * p) const {

  // checks whether the particle should be rescattered
  assert(p);
  // attempt to rescatter anything marked as 'hadron in the nucleus'
  return ( p->Status() == kIStHadronInTheNucleus );

}

//______________________________________________________________________________
void INCLCascadeIntranuke::Configure(const Registry & config) {

  LOG("INCLCascadeIntranuke", pDEBUG)
    << "Configure from Registry: '" << config.Name() << "'\n"
    << config;

  Algorithm::Configure(config);
  this->LoadConfig();

}

//___________________________________________________________________________
void INCLCascadeIntranuke::Configure(string param_set) {

  LOG("INCLCascadeIntranuke", pDEBUG)
    << "Configure from param_set name: " << param_set;

  Algorithm::Configure(param_set);
  this->LoadConfig();

}

bool INCLCascadeIntranuke::continueCascade() const{

  bool continueCascade_ = true;

  LOG("INCLCascadeIntranuke", pINFO) << propagationModel->getCurrentTime() << "  "  <<  propagationModel->getStoppingTime() << "  " << continueCascade_;

  if(propagationModel->getCurrentTime() > propagationModel->getStoppingTime()){
    continueCascade_ = false;
  }
  LOG("INCLCascadeIntranuke", pINFO) << incl_target->getStore()->getBook().getCascading() << "  " << incl_target->getStore()->getIncomingParticles().empty() << "  " << continueCascade_;
  if(incl_target->getStore()->getBook().getCascading()==0 &&
      incl_target->getStore()->getIncomingParticles().empty()){
    continueCascade_ = false;
  }
  int minRemnantSize = 4;
  if(incl_target->getA() <= minRemnantSize) {
    continueCascade_ = false;
  }

  LOG("INCLCascadeIntranuke", pINFO) << incl_target->getTryCompoundNucleus()<<  "  " << continueCascade_;

  if(incl_target->getTryCompoundNucleus()) {
    continueCascade_ = false;
  }
  return continueCascade_;

}

void INCLCascadeIntranuke::fillEventRecord(G4INCL::FinalState *fs, G4INCL::ParticleList mother_list, GHepRecord * evrec, double time) const {
  std::vector<INCLRecord> stepParticleList;
  stepParticleList.clear();

  istep++;

  for(ParticleIter iter=mother_list.begin(); iter!=mother_list.end(); ++iter){
    this->fillStep(*iter, stepParticleList, -6, time);
    ParticleSpecies ptype = (*iter)->getSpecies();
    stepFinalState[istep].emplace_back((*iter)->getID(), ptype.getPDGCode(), -6, 0);
  }

  ParticleList modified = fs->getModifiedParticles();
  for(ParticleIter iter=modified.begin(); iter!=modified.end(); ++iter){
    this->fillStep(*iter, stepParticleList, -2, time);
    ParticleSpecies ptype = (*iter)->getSpecies();
    stepFinalState[istep].emplace_back((*iter)->getID(), ptype.getPDGCode(), -2, 0);
    //stepParticleList.emplace_back((*iter)->getID(), ptype.getPDGCode(), -2, 0, p4mom, p4posi);
  }
  ParticleList outgoing = fs->getOutgoingParticles();
  for(ParticleIter iter=outgoing.begin(); iter!=outgoing.end(); ++iter){
    this->fillStep(*iter, stepParticleList, -3, time);
    ParticleSpecies ptype = (*iter)->getSpecies();
    stepFinalState[istep].emplace_back((*iter)->getID(), ptype.getPDGCode(), -3, 0);
    //stepParticleList.emplace_back((*iter)->getID(), ptype.getPDGCode(), -3, 0);
  }
  ParticleList destroyed = fs->getDestroyedParticles();
  for(ParticleIter iter=destroyed.begin();  iter!=destroyed.end(); ++iter){
    this->fillStep(*iter, stepParticleList, -4, time);
    ParticleSpecies ptype = (*iter)->getSpecies();
    stepFinalState[istep].emplace_back((*iter)->getID(), ptype.getPDGCode(), -4, 0);
    //stepParticleList.emplace_back((*iter)->getID(), ptype.getPDGCode(), -4, 0);
  }
  ParticleList created = fs->getCreatedParticles();
  for(ParticleIter iter=created.begin(); iter!=created.end(); ++iter){
    this->fillStep(*iter, stepParticleList, -5, time);
    ParticleSpecies ptype = (*iter)->getSpecies();
    stepFinalState[istep].emplace_back((*iter)->getID(), ptype.getPDGCode(), -5, 0);
    //stepParticleList.emplace_back((*iter)->getID(), ptype.getPDGCode(), -5, 0);
  }

  LOG("INCLCascadeIntranuke", pINFO) << "1";

  int num_partiles = tempFinalState.size();
  int num_temp = stepParticleList.size();


  if(num_temp >= 1){
    std::map<int, int> mother_index;  // global id is the key; local id is the value
    mother_index.clear();
    int local_id = -1;
    for(auto it = stepParticleList.begin(); it != stepParticleList.end(); ++it){
      if(it->mother_index != -6) continue;
      for(auto ip = tempFinalState.begin(); ip != tempFinalState.end(); ++ip){
	LOG("INCLCascadeIntranuke", pINFO) << it->global_index << " " << ip->global_index;
	if(it->global_index == ip->global_index){
	  mother_index[ip->global_index] = ip->local_index;
	  local_id = ip->local_index;
	}
      }
    }

    int index_ = num_partiles;

    LOG("INCLCascadeIntranuke", pINFO) << "1";
    for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
    LOG("INCLCascadeIntranuke", pINFO) << "1";
      if(ip->mother_index == -4) continue;
    LOG("INCLCascadeIntranuke", pINFO) << "1";
      if(ip->mother_index == -6) continue;
    LOG("INCLCascadeIntranuke", pINFO) << " num_temp " << num_temp;
      if(num_temp == 2 && ip->mother_index == -3){
	tempFinalState.emplace_back(ip->global_index, ip->pdgid, local_id, index_, ip->p4mom, ip->p4posi);
	// pdg, status, first mother, second mother, first daughter, second daughter, mom, posi
	// FIXME: test code
	evrec->Particle(local_id)->SetRescatterCode(0);
	int pdg = ip->pdgid;

	LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << ip->pdgid;
	if(ip->theType == G4INCL::Composite){
	  if(pdg != 2212 && pdg != 2112 && pdg != 3122){
	    int A = pdg % 1000;
	    int Z = pdg / 1000;
	    pdg = genie::pdg::IonPdgCode( A , Z );
	  }
	}
	GHepParticle p(pdg, kIStStableFinalState, local_id, -1, -1, -1, ip->p4mom, ip->p4posi);
	evrec->AddParticle(p);
	index_++;
      }
      else if(num_temp > 2){
    LOG("INCLCascadeIntranuke", pINFO) << "1";
	tempFinalState.emplace_back(ip->global_index, ip->pdgid, local_id, index_, ip->p4mom, ip->p4posi);
    LOG("INCLCascadeIntranuke", pINFO) << "local_id " << local_id;
	// pdg, status, first mother, second mother, first daughter, second daughter, mom, posi
	// FIXME: test code
	evrec->Particle(local_id)->SetRescatterCode(1);
	int pdg = ip->pdgid;
	LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << ip->pdgid;
	if(ip->theType == G4INCL::Composite){
	  if(pdg != 2212 && pdg != 2112 && pdg != 3122){
	    int A = pdg % 1000;
	    int Z = pdg / 1000;
	    pdg = genie::pdg::IonPdgCode( A , Z );
	    LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << pdg;
	  }
	}
	GHepParticle p(pdg, kIStHadronInTheNucleus, local_id, -1, -1, -1, ip->p4mom, ip->p4posi);
	evrec->AddParticle(p);
	index_++;
      }
    }
  }

}

void INCLCascadeIntranuke::fillStep(G4INCL::Particle *par, std::vector<INCLRecord> &stepList, int type, double time) const {

  ParticleSpecies ptype = par->getSpecies();
  TLorentzVector p4mom;
  p4mom.SetPx(par->getMomentum().getX() / 1000.);
  p4mom.SetPy(par->getMomentum().getY() / 1000.);
  p4mom.SetPz(par->getMomentum().getZ() / 1000.);
  p4mom.SetE(par->getEnergy() / 1000.);
  TLorentzVector p4posi;
  p4posi.SetX(par->getPosition().getX());
  p4posi.SetY(par->getPosition().getY());
  p4posi.SetZ(par->getPosition().getZ());
  p4posi.SetT(time);
  stepList.emplace_back(par->getID(), ptype.getPDGCode(), type, 0, p4mom, p4posi, ptype.theType);
}

void INCLCascadeIntranuke::postCascade() const {
  theEventInfo.stoppingTime = propagationModel->getCurrentTime();
  // The event bias
  theEventInfo.eventBias = (Double_t) Particle::getTotalBias();

  // Check if the nucleus contains strange particles
  theEventInfo.sigmasInside = incl_target->containsSigma();
  theEventInfo.antikaonsInside = incl_target->containsAntiKaon();
  theEventInfo.lambdasInside = incl_target->containsLambda();
  theEventInfo.kaonsInside = incl_target->containsKaon();
  
  // Capture antiKaons and Sigmas and produce Lambda instead
  theEventInfo.absorbedStrangeParticle = incl_target->decayInsideStrangeParticles();
  
  // Emit strange particles still inside the nucleus
  incl_target->emitInsideStrangeParticles();
  theEventInfo.emitKaon = incl_target->emitInsideKaon();

  // Check if the nucleus contains deltas
  theEventInfo.deltasInside = incl_target->containsDeltas();

  // Take care of any remaining deltas
  theEventInfo.forcedDeltasOutside = incl_target->decayOutgoingDeltas();
  theEventInfo.forcedDeltasInside = incl_target->decayInsideDeltas();

  // Take care of any remaining etas, omegas, neutral Sigmas and/or neutral kaons
  double timeThreshold=theConfig->getDecayTimeThreshold();
  theEventInfo.forcedPionResonancesOutside = incl_target->decayOutgoingPionResonances(timeThreshold);
  incl_target->decayOutgoingSigmaZero(timeThreshold);
  incl_target->decayOutgoingNeutralKaon();

  // Apply Coulomb distortion, if appropriate
  // Note that this will apply Coulomb distortion also on pions emitted by
  // unphysical remnants (see decayInsideDeltas). This is at variance with
  // what INCL4.6 does, but these events are (should be!) so rare that
  // whatever we do doesn't (shouldn't!) make any noticeable difference.
  // CoulombDistortion::distortOut(incl_target->getStore()->getOutgoingParticles(), incl_target);

  // If the normal cascade predicted complete fusion, use the tabulated
  // masses to compute the excitation energy, the recoil, etc.
  if(incl_target->getStore()->getOutgoingParticles().size()==0
     && (!incl_target->getProjectileRemnant()
         || incl_target->getProjectileRemnant()->getParticles().size()==0)) {

    LOG("INCLCascadeIntranuke", pINFO) << "Cascade resulted in complete fusion, using realistic fusion kinematics";

    incl_target->useFusionKinematics();

    if(incl_target->getExcitationEnergy()<0.) {
      // Complete fusion is energetically impossible, return a transparent
      LOG("INCLCascadeIntranuke", pINFO) << "Complete-fusion kinematics yields negative excitation energy, returning a transparent!";
      theEventInfo.transparent = true;
      return;
    }

  }
  else {

    // Set the excitation energy
    incl_target->setExcitationEnergy(incl_target->computeExcitationEnergy());

    // Make a projectile pre-fragment out of the geometrical and dynamical
    // spectators
    //theEventInfo.nUnmergedSpectators = makeProjectileRemnant();

    // Compute recoil momentum, energy and spin of the nucleus
    int minRemnantSize = 4;
    if(incl_target->getA()==1 && minRemnantSize>1) {
      LOG("INCLCascadeIntranuke", pINFO) << "Computing one-nucleon recoil kinematics. We should never be here nowadays, cascade should stop earlier than this.";
    }
    incl_target->computeRecoilKinematics();

    // Make room for the remnant recoil by rescaling the energies of the
    // outgoing particles. FIXME
    if(incl_target->hasRemnant()) this->rescaleOutgoingForRecoil();
  }

  theEventInfo.clusterDecay = incl_target->decayOutgoingClusters() || incl_target->decayMe(); 
  incl_target->fillEventInfo(&theEventInfo);

}

bool INCLCascadeIntranuke::preCascade() const {
  // Reset theEventInfo
  theEventInfo.reset();
  EventInfo::eventNumber++;

  // Fill in the event information
  // neutrino projectile
  // theEventInfo.projectileType = projectileSpecies.theType;
  theEventInfo.Ap = 0;
  theEventInfo.Zp = 0;
  theEventInfo.Sp = 0;
  theEventInfo.Ep = prob->P4()->E() * 1000.;

  theEventInfo.At = incl_target->getA();
  theEventInfo.Zt = incl_target->getZ();
  theEventInfo.St = incl_target->getS();

  theEventInfo.transparent = false;
  theEventInfo.impactParameter = 0.;

  theEventInfo.effectiveImpactParameter = 0.;

  return true;

}

void INCLCascadeIntranuke::rescaleOutgoingForRecoil() const {
  TLorentzVector *p4prob = prob->P4();
  TLorentzVector *p4lep  = primarylepton->P4();
  G4INCL::ThreeVector transferQ((p4prob->Px() - p4lep->Px()) * 1000.,
      (p4prob->Py() - p4lep->Py()) * 1000.,
      (p4prob->Pz() - p4lep->Pz()) * 1000.);

  RecoilCMFunctor theRecoilFunctor(incl_target, theEventInfo, transferQ);

  // Apply the root-finding algorithm
  const G4INCL::RootFinder::Solution theSolution = G4INCL::RootFinder::solve(&theRecoilFunctor, 1.0);
  if(theSolution.success) {
    theRecoilFunctor(theSolution.x); // Apply the solution
  } else {
    LOG("INCLCascadeIntranuke", pINFO) << "Couldn't accommodate remnant recoil while satisfying energy conservation, root-finding algorithm failed.";
  }
}

int INCLCascadeIntranuke::INCLPDG_to_GHEPPDG(int pdg, int A, int Z) const{
  TParticlePDG * p = PDGLibrary::Instance()->Find(pdg);
  if(!p){
    if(A != 0){
    	return genie::pdg::IonPdgCode( A , Z );
    }
    else 
      return 0;
  }
  else
    return pdg;

}

#endif // __GENIE_INCL_ENABLED__
