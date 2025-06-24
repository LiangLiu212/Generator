#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

//---------------------
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

// GENIE
#include "INCLCascadeIntranuke.h"
#include "Framework/ParticleData/BaryonResUtils.h"

// INCL++
#include "G4INCLConfig.hh"
#include "G4INCLCascade.hh"
#include "G4INCLConfigEnums.hh"
#include "G4INCLParticle.hh"
// signal handler (for Linux and GCC)
#include "G4INCLSignalHandling.hh"


#include "G4INCLPauliBlocking.hh"

#include "G4INCLCrossSections.hh"
#include "G4INCLDecayAvatar.hh"

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

#include "G4INCLClusterDecay.hh"
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
#include "Physics/HadronTransport/G4INCLGENIEAvatar.h"
#include "Physics/HadronTransport/G4INCLGENIEParticleRecord.h"

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
  fResonanceDecayer = nullptr;
  fResonanceDecayer = dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("Decayer"));
  assert(fResonanceDecayer);
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

   evrec->Print(std::cout);
  // LOG("INCLCascadeIntranuke", pINFO) << evrec->Summary()->ProcInfo().ScatteringTypeAsString();
  // LOG("INCLCascadeIntranuke", pINFO) << evrec->Summary()->ProcInfo().ScatteringTypeId();
  // LOG("INCLCascadeIntranuke", pINFO) << evrec->Summary()->ProcInfo().InteractionTypeId();
  // LOG("INCLCascadeIntranuke", pINFO) << evrec->Summary()->ProcInfo().InteractionTypeAsString();
  // LOG("INCLCascadeIntranuke", pNOTICE) << "is resonacne : " << evrec->Summary()->ProcInfo().IsResonant();

  this->DecayResonance(evrec);

  // LOG("INCLCascadeIntranuke", pNOTICE) << "is resonacne : " << evrec->Summary()->ProcInfo().IsResonant();

  //  Get 'nuclear environment' at the beginning of hadron transport
  //  and keep track of the remnant nucleus A,Z

  //    GEvGenMode_t fGMode = evrec->EventGenerationMode();
  //    if(fGMode == kGMdLeptonNucleus ||
  //	  fGMode == kGMdDarkMatterNucleus ||
  //	  fGMode == kGMdNucleonDecay) {
  //      int inucl = evrec->RemnantNucleusPosition();
  //      evrec->Particle(inucl)->SetStatus(kIStIntermediateState);
  //    }

  // FIXME: start with the NC process
  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  incl_target = incl_nucleus->getNuclues();
  theConfig = incl_nucleus->getConfig();
  propagationModel =  incl_nucleus->getPropagationModel();
  incl_target->setParticleNucleusCollision();
  std::unique_ptr<G4INCL::FinalState> finalState(new FinalState);
  primarylepton = evrec->FinalStatePrimaryLepton();
  prob = evrec->Probe();
  target = evrec->TargetNucleus();

  this->preCascade();

  double currentTime = 0.0;
  //    double temfin;

  tempFinalState.clear();
  stepFinalState.clear();
  istep = 0;

  double maximumTime = 29.8 * std::pow(incl_target->getA(), 0.16);
  propagationModel->setStoppingTime(maximumTime);
  propagationModel->setNucleus(incl_target);
  propagationModel->generateAllAvatars();

  this->fillFinalState(evrec, finalState.get());

  /*
     int idx = 0;
     double TLab;
     temfin = 29.8 * std::pow(incl_target->getA(), 0.16);
     if(TLab > 2000.)
     temfin *= (5.8E4-TLab)/5.6E4;
     */

  //    double maximumTime = temfin;
  //    propagationModel->setStoppingTime(maximumTime);
  //    propagationModel->setNucleus(incl_target);
  //    propagationModel->generateAllAvatars();

  // stopping time: 
  // INCL don't have the stopping time for neutrino.
  // we can calculate the longest stopping time for all the daughters from primary interaction.

  incl_target->applyFinalState(finalState.get());
  //    incl_target->getStore()->getBook().incrementCascading();   // FIXME
  incl_target->getStore()->getBook().incrementAcceptedCollisions();
  int step = 0;

  while(step < 10000000 && continueCascade()){
    IAvatar *avatar = propagationModel->propagate(finalState.get());
    finalState->reset();
    if(avatar == 0) break; // No more avatars in the avatar list.
    G4INCL::ParticleList mother_list = avatar->getParticles();
    backup_mother.clear();
    for(G4INCL::ParticleIter imom =  mother_list.begin(); imom != mother_list.end(); imom++){
      this->fillStep(*imom, backup_mother, -1, -1);
    }
    avatar->fillFinalState(finalState.get());
    // Must fill event record before incl nucleus applyFinalState
    // applyFinalState will delete destroyed particles.
    this->fillEventRecord(finalState.get(), mother_list, evrec, propagationModel->getCurrentTime(), avatar->getType());
    incl_target->applyFinalState(finalState.get());
    delete avatar;
    step++;
  }

  // std::cout << incl_target->print() << std::endl;
  // put the nuclear remnant in the event record

  this->postCascade(evrec, finalState.get());

  TObjArrayIter piter(evrec);
  GHepParticle * p = nullptr;
  int stable_finalstate = 0;
  while ( (p = (GHepParticle *) piter.Next() ) ) {
    if(p->Status() == kIStStableFinalState){
      // FIXME: for NC QE, it is only neutron and proton
      stable_finalstate++;
    }
  }

  const int n_outgoing = theEventInfo.nParticles;

  // LOG("INCLCascadeIntranuke", pNOTICE) << "Final State Before De-Excitation : " << stable_finalstate - theEventInfo.nParticles;
  // LOG("INCLCascadeIntranuke", pNOTICE) << "nRemnants   : " << theEventInfo.nRemnants;
  // LOG("INCLCascadeIntranuke", pNOTICE) << "nRemnants A : " << theEventInfo.ARem[0];
  // LOG("INCLCascadeIntranuke", pNOTICE) << "nRemnants Z : " << theEventInfo.ZRem[0];
  // LOG("INCLCascadeIntranuke", pNOTICE) << "nRemnants S : " << theEventInfo.SRem[0];
  // LOG("INCLCascadeIntranuke", pNOTICE) << "nParticles   : " << theEventInfo.nParticles;
  // for(int i = 0; i < theEventInfo.nParticles; i++){
  //   LOG("INCLCascadeIntranuke", pNOTICE) << "Final State Particles PDG: " << theEventInfo.PDGCode[i];
  //   LOG("INCLCascadeIntranuke", pNOTICE) << "Final State Particles px : " << theEventInfo.px[i];
  //   LOG("INCLCascadeIntranuke", pNOTICE) << "Final State Particles py : " << theEventInfo.py[i];
  //   LOG("INCLCascadeIntranuke", pNOTICE) << "Final State Particles pz : " << theEventInfo.pz[i];
  // }

  // LOG("INCLCascadeIntranuke", pINFO) << incl_target->print();

  double Rem_p2 = theEventInfo.pxRem[0]*theEventInfo.pxRem[0]
    + theEventInfo.pyRem[0]*theEventInfo.pyRem[0]
    + theEventInfo.pzRem[0]*theEventInfo.pzRem[0];
  double Rem_Kin = theEventInfo.EKinRem[0];
  double Rem_M = (Rem_p2 - Rem_Kin*Rem_Kin) / ( 2.0 * Rem_Kin);

  // GHepParticle *p4l = evrec->FinalStatePrimaryLepton();

  double Rem_px = theEventInfo.pxRem[0]/1000.;
  double Rem_py = theEventInfo.pyRem[0]/1000.;
  double Rem_pz = theEventInfo.pzRem[0]/1000.;
  double Rem_E =  std::sqrt(Rem_p2 + Rem_M*Rem_M)/1000.;

  TLorentzVector p4mom(Rem_px, Rem_py, Rem_pz, Rem_E);

  TLorentzVector p4posi(0,0,0,0);
  int A = incl_target->getA();
  int Z = incl_target->getZ();
  int S = incl_target->getS(); // INCL and ABLA support hypernuclei
  int pdg = genie::pdg::IonPdgCode( A , Z, 0, 0 );

  // LOG("INCLCascadeIntranuke", pNOTICE) << pdg;
  TParticlePDG * prem = PDGLibrary::Instance()->Find(pdg);
  int PreDeExPDG = pdg;
  if(!prem) PreDeExPDG = kPdgHadronicBlob;
  evrec->AddParticle(PreDeExPDG, kIStPreDeExNuclearRemnant, 3, -1, -1, -1, p4mom, p4posi);

  switch(theConfig->getDeExcitationType()){
    case G4INCL::DeExcitationABLAXX:
      {
        std::unique_ptr<G4INCL::IDeExcitation> theDeExcitation = std::make_unique<G4INCLAblaInterface>(theConfig);
        theDeExcitation->deExcite(&theEventInfo);
        break;
      }
    case G4INCL::DeExcitationABLA07:
      {
        std::unique_ptr<G4INCL::IDeExcitation> theDeExcitation = std::make_unique<ABLA07CXX::Abla07Interface>(theConfig);
        theDeExcitation->deExcite(&theEventInfo);
        break;
      }
    case G4INCL::DeExcitationGEMINIXX:
      {
        std::unique_ptr<G4INCL::IDeExcitation> theDeExcitation = std::make_unique<G4INCLGEMINIXXInterface>(theConfig);
        theDeExcitation->deExcite(&theEventInfo);
        break;
      }
    default:
      {
        exit(1);
        break;
      }
  }

  // particle from de-excitation
  const int remnant_id = evrec->GetEntries() - 1;
  if(theEventInfo.nParticles > n_outgoing){
    for(int i = n_outgoing; i < theEventInfo.nParticles; i++){
      LOG("INCLCascadeIntranuke", pNOTICE) << "Final State Particles PDG: " << theEventInfo.PDGCode[i];
      LOG("INCLCascadeIntranuke", pNOTICE) << "Final State Particles PDG: " << this->INCLPDG_to_GHEPPDG(theEventInfo.PDGCode[i], theEventInfo.A[i], theEventInfo.Z[i], theEventInfo.S[i]);
      int depdg = this->INCLPDG_to_GHEPPDG(theEventInfo.PDGCode[i], theEventInfo.A[i], theEventInfo.Z[i], theEventInfo.S[i]);
      TParticlePDG * p = PDGLibrary::Instance()->Find(depdg);

      LOG("INCLCascadeIntranuke", pNOTICE) << depdg;

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
        ptype = kIStStableFinalState;
      }
      evrec->AddParticle(depdg, ptype, remnant_id, -1, -1, -1, p4mom, p4posi);
    }
  }
  else{
    TLorentzVector p4posi(0,0,0,0);
    int A = incl_target->getA();
    int Z = incl_target->getZ();
    int S = incl_target->getS();
    int pdg = genie::pdg::IonPdgCode( A , Z, std::abs(S), 0 );
    TParticlePDG * p = PDGLibrary::Instance()->Find(pdg);
    double M = p->Mass();
    double Rem_E = sqrt(Rem_p2/1000000. + M*M);
    TLorentzVector p4mom(Rem_px, Rem_py, Rem_pz, Rem_E);
    LOG("INCLCascadeIntranuke", pNOTICE) << pdg;
    evrec->AddParticle(pdg, kIStFinalStateNuclearRemnant, remnant_id, -1, -1, -1, p4mom, p4posi);
  }


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
      LOG("INCLCascadeIntranuke", pINFO)  << fs_particle  << "  global id : " << ip->global_index << "  " 
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

  if(propagationModel->getCurrentTime() > propagationModel->getStoppingTime()){
    continueCascade_ = false;
  }
  if(incl_target->getStore()->getBook().getCascading()==0 &&
      incl_target->getStore()->getIncomingParticles().empty()){
    continueCascade_ = false;
  }
  int minRemnantSize = 4;
  if(incl_target->getA() <= minRemnantSize) {
    continueCascade_ = false;
  }

  if(incl_target->getTryCompoundNucleus()) {
    continueCascade_ = false;
  }
  return continueCascade_;

}

void INCLCascadeIntranuke::fillEventRecord(G4INCL::FinalState *fs, G4INCL::ParticleList mother_list, GHepRecord * evrec, double time, G4INCL::AvatarType avaType) const {
  std::vector<INCLRecord> stepParticleList;
  stepParticleList.clear();
  istep++;
  LOG("INCLCascadeIntranuke", pDEBUG) << "istep : " << istep;
  for(ParticleIter iter=mother_list.begin(); iter!=mother_list.end(); ++iter){
    this->fillStep(*iter, stepParticleList, -6, time);
    ParticleSpecies ptype = (*iter)->getSpecies();
    stepFinalState[istep].emplace_back((*iter)->getID(), ptype.getPDGCode(), -6, 0);
    LOG("INCLCascadeIntranuke", pDEBUG) << "mother list ID : " << (*iter)->getID() << " pdg : " << ptype.getPDGCode(); 
  }
  ParticleList modified = fs->getModifiedParticles();
  for(ParticleIter iter=modified.begin(); iter!=modified.end(); ++iter){
    this->fillStep(*iter, stepParticleList, -2, time);
    ParticleSpecies ptype = (*iter)->getSpecies();
    stepFinalState[istep].emplace_back((*iter)->getID(), ptype.getPDGCode(), -2, 0);
    LOG("INCLCascadeIntranuke", pDEBUG) << "Modified ID : " << (*iter)->getID() << " pdg : " << ptype.getPDGCode(); 
  }
  ParticleList outgoing = fs->getOutgoingParticles();
  for(ParticleIter iter=outgoing.begin(); iter!=outgoing.end(); ++iter){
    this->fillStep(*iter, stepParticleList, -3, time);
    ParticleSpecies ptype = (*iter)->getSpecies();
    stepFinalState[istep].emplace_back((*iter)->getID(), ptype.getPDGCode(), -3, 0);
    LOG("INCLCascadeIntranuke", pDEBUG) << "Outgoing ID : " << (*iter)->getID() << " pdg : " << ptype.getPDGCode(); 
  }
  ParticleList destroyed = fs->getDestroyedParticles();
  for(ParticleIter iter=destroyed.begin();  iter!=destroyed.end(); ++iter){
    this->fillStep(*iter, stepParticleList, -4, time);
    ParticleSpecies ptype = (*iter)->getSpecies();
    stepFinalState[istep].emplace_back((*iter)->getID(), ptype.getPDGCode(), -4, 0);
    LOG("INCLCascadeIntranuke", pDEBUG) << "Destroyed ID : " << (*iter)->getID() << " pdg : " << ptype.getPDGCode(); 
  }
  ParticleList created = fs->getCreatedParticles();
  for(ParticleIter iter=created.begin(); iter!=created.end(); ++iter){
    this->fillStep(*iter, stepParticleList, -5, time);
    ParticleSpecies ptype = (*iter)->getSpecies();
    stepFinalState[istep].emplace_back((*iter)->getID(), ptype.getPDGCode(), -5, 0);
    LOG("INCLCascadeIntranuke", pDEBUG) << "Created ID : " << (*iter)->getID() << " pdg : " << ptype.getPDGCode(); 
  }

  int num_partiles = tempFinalState.size();
  int num_temp = stepParticleList.size();
  if(created.size() == 0 && outgoing.size() == 0 && modified.size() == 0) return;
  int index_ = num_partiles;
  if(mother_list.size() == 1){
    // FIXME: channel with mother size = 1 could be reflection, outgoing and decay
    // Try to match the mother in step to the event record history
    // if the mother of this step is not in the event record history and it is a reflection
    // avatar, we will not put this step into event record history
    int mother_position = -1;
    for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
      if(ip->mother_index != -6) continue;
      for(auto it = tempFinalState.begin(); it != tempFinalState.end(); ++it){
        if(it->global_index == ip->global_index){
          mother_position = it->local_index;
        }
      }
    }

    if(mother_position != -1){

      if(avaType == G4INCL::SurfaceAvatarType){
        //
        for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
          if(ip->mother_index == -6 || ip->mother_index == -4) continue;
          LOG("INCLCascadeIntranuke", pNOTICE) << "created: " << created.size();
          LOG("INCLCascadeIntranuke", pNOTICE) << "outgoing: " << outgoing.size();
          LOG("INCLCascadeIntranuke", pNOTICE) << "modified: " << modified.size();
          LOG("INCLCascadeIntranuke", pNOTICE) << "destroyed: " << destroyed.size();
          // put this particle in event record
          tempFinalState.emplace_back(ip->global_index, ip->pdgid, mother_position, index_++, ip->p4mom, ip->p4posi);
          evrec->Particle(mother_position)->SetRescatterCode(int(avaType));

          // get the pdg in genie style
          int pdg = ip->pdgid;
          LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << ip->pdgid;
          if(ip->theType == G4INCL::Composite){
            if(pdg != 2212 && pdg != 2112 && pdg != 3122){
              int S = pdg / 1000000;
              int pdg_no_s = pdg % 1000000;
              int A = pdg_no_s % 1000;
              int Z = pdg_no_s / 1000;
              pdg = genie::pdg::IonPdgCode( A , Z, S, 0);
            }
          }
          LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << pdg;
          // get the type of the particle
          EGHepStatus ptype;

          if(ip->mother_index == -3){
            //ptype = kIStPreDecayResonantState;
            ptype = kIStStableFinalState;
          }
          else{
            ptype = kIStHadronInTheNucleus;
          }

          GHepParticle p(pdg, ptype, mother_position, -1, -1, -1, ip->p4mom, ip->p4posi);
          evrec->AddParticle(p);
        }
      }
      else if(avaType == G4INCL::DecayAvatarType){

        /* If the decay was Pauli-blocked, make sure the propagation model
         * generates a new decay avatar on the next call to propagate().
         *
         * \bug{Note that we don't generate new decay avatars for deltas that
         * could not satisfy energy conservation. This is in keeping with
         * INCL4.6, but doesn't seem to make much sense to me (DM), as energy
         * conservation can be impossible to satisfy due to weird local-energy
         * conditions, for example, that evolve with time.}
         */

        /* decay channel without daughters, skip it
         * a decay must have more than 1 daughters
         */

        if((created.size() + outgoing.size() + modified.size()) != 1){
          evrec->Particle(mother_position)->SetRescatterCode(int(avaType));
          evrec->Particle(mother_position)->SetStatus(kIStDecayedState);
          for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
            if(ip->mother_index == -6 || ip->mother_index == -4) continue;
            LOG("INCLCascadeIntranuke", pNOTICE) << "created: " << created.size();
            LOG("INCLCascadeIntranuke", pNOTICE) << "outgoing: " << outgoing.size();
            LOG("INCLCascadeIntranuke", pNOTICE) << "modified: " << modified.size();
            LOG("INCLCascadeIntranuke", pNOTICE) << "destroyed: " << destroyed.size();
            // put this particle in event record
            tempFinalState.emplace_back(ip->global_index, ip->pdgid, mother_position, index_++, ip->p4mom, ip->p4posi);

            // get the pdg in genie style
            int pdg = ip->pdgid;
            LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << ip->pdgid;
            if(ip->theType == G4INCL::Composite){
              if(pdg != 2212 && pdg != 2112 && pdg != 3122){
                int S = pdg / 1000000;
                int pdg_no_s = pdg % 1000000;
                int A = pdg_no_s % 1000;
                int Z = pdg_no_s / 1000;
                pdg = genie::pdg::IonPdgCode( A , Z, S, 0);
              }
            }
            LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << pdg;
            // get the type of the particle
            EGHepStatus ptype;

            if(ip->mother_index == -3){
              //ptype = kIStPreDecayResonantState;
              ptype = kIStStableFinalState;
            }
            else{
              ptype = kIStHadronInTheNucleus;
            }

            GHepParticle p(pdg, ptype, mother_position, -1, -1, -1, ip->p4mom, ip->p4posi);
            evrec->AddParticle(p);
          }
        }
      }
      else if(avaType == G4INCL::UnknownParticle){

      }
    }
  }
  else if(mother_list.size() == 2){
    int first_mother_position = -1;
    int second_mother_position = -1;
    // find parents for nn-collision
    for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
      if(ip->mother_index != -6) continue;
      int mother_position = -1;
      bool is_spectator = true;
      for(auto it = tempFinalState.begin(); it != tempFinalState.end(); ++it){
        if(it->global_index == ip->global_index){
          mother_position = it->local_index;
          is_spectator = false;
        }
      }
      if(is_spectator){
        for(auto imom = backup_mother.begin(); imom != backup_mother.end(); ++imom){
          if(imom->global_index == ip->global_index){
            tempFinalState.emplace_back(ip->global_index, imom->pdgid, -1, index_++, imom->p4mom, imom->p4posi);
            GHepParticle p(imom->pdgid, kIStSpectator, -1, -1, -1, -1, imom->p4mom, imom->p4posi);
            evrec->AddParticle(p);
          }
        }
        mother_position = index_ - 1;
      }
      // assign the mother index
      if(first_mother_position == -1)
        first_mother_position = mother_position;
      else
        second_mother_position = mother_position;
    }

    if(first_mother_position > second_mother_position){
      int temp_position = first_mother_position;
      first_mother_position = second_mother_position;
      second_mother_position = temp_position;
    }
    LOG("INCLCascadeIntranuke", pNOTICE) << first_mother_position;
    LOG("INCLCascadeIntranuke", pNOTICE) << second_mother_position;
    evrec->Particle(first_mother_position)->SetRescatterCode(int(avaType));
    evrec->Particle(second_mother_position)->SetRescatterCode(int(avaType));

    for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
      if(ip->mother_index == -6 || ip->mother_index == -4) continue;
      LOG("INCLCascadeIntranuke", pNOTICE) << "created: " << created.size();
      LOG("INCLCascadeIntranuke", pNOTICE) << "outgoing: " << outgoing.size();
      LOG("INCLCascadeIntranuke", pNOTICE) << "modified: " << modified.size();
      LOG("INCLCascadeIntranuke", pNOTICE) << "destroyed: " << destroyed.size();
      // put this particle in event record
      tempFinalState.emplace_back(ip->global_index, ip->pdgid, first_mother_position, index_++, ip->p4mom, ip->p4posi);

      // get the pdg in genie style
      int pdg = ip->pdgid;
      LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << ip->pdgid;
      if(ip->theType == G4INCL::Composite){
        if(pdg != 2212 && pdg != 2112 && pdg != 3122){
          int S = pdg / 1000000;
          int pdg_no_s = pdg % 1000000;
          int A = pdg_no_s % 1000;
          int Z = pdg_no_s / 1000;
          pdg = genie::pdg::IonPdgCode( A , Z, S, 0);
        }
      }
      LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << pdg;
      // get the type of the particle
      EGHepStatus ptype;
      ptype = kIStHadronInTheNucleus;
      GHepParticle p(pdg, ptype, first_mother_position, second_mother_position, -1, -1, ip->p4mom, ip->p4posi);
      evrec->AddParticle(p);
    }
    // update the daughter's indices for spectator
    if(evrec->Particle(first_mother_position)->FirstDaughter() == -1){
      evrec->Particle(first_mother_position)->SetFirstDaughter(evrec->Particle(second_mother_position)->FirstDaughter());
      evrec->Particle(first_mother_position)->SetLastDaughter(evrec->Particle(second_mother_position)->LastDaughter());
    }
    else if(evrec->Particle(second_mother_position)->FirstDaughter() == -1){
      evrec->Particle(second_mother_position)->SetFirstDaughter(evrec->Particle(first_mother_position)->FirstDaughter());
      evrec->Particle(second_mother_position)->SetLastDaughter(evrec->Particle( first_mother_position)->LastDaughter());
    }

  }
  else{
    LOG("INCLCascadeIntranuke", pNOTICE) << "wrong mother size : " << mother_list.size();
    LOG("INCLCascadeIntranuke", pNOTICE) << "created: " << created.size();
    LOG("INCLCascadeIntranuke", pNOTICE) << "outgoing: " << outgoing.size();
    LOG("INCLCascadeIntranuke", pNOTICE) << "modified: " << modified.size();
    LOG("INCLCascadeIntranuke", pNOTICE) << "destroyed: " << destroyed.size();
    exit(1);
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

void INCLCascadeIntranuke::postCascade(GHepRecord * evrec, G4INCL::FinalState * finalState) const {
  // Fill in the event information
  theEventInfo.stoppingTime = propagationModel->getCurrentTime();

  // The event bias
  theEventInfo.eventBias = (Double_t) Particle::getTotalBias();

  // Check if the nucleus contains strange particles
  theEventInfo.sigmasInside = incl_target->containsSigma();
  theEventInfo.antikaonsInside = incl_target->containsAntiKaon();
  theEventInfo.lambdasInside = incl_target->containsLambda();
  theEventInfo.kaonsInside = incl_target->containsKaon();

  LOG("INCLCascadeIntranuke", pNOTICE) << "n Sigma and anti-Kaon " << incl_target->containsSigma() << "  " << incl_target->containsAntiKaon() << "  " << incl_target->containsLambda() << "  " << incl_target->containsKaon();

  // FIXME: It seems the INCL will transfer Sigma and anti-Kaon
  // into Lambda via interacting with the first nucleon in Store
  // list
  // - However, de-excitation
  // Capture antiKaons and Sigmas and produce Lambda instead
  theEventInfo.absorbedStrangeParticle = this->decayInsideStrangeParticles(evrec, finalState);

  // Emit strange particles still inside the nucleus
  this->emitInsideStrangeParticles(evrec, finalState);
  theEventInfo.emitKaon = this->emitInsideKaon(evrec, finalState);
  theEventInfo.emitLambda = this->emitInsideLambda(evrec, finalState);

  // Check if the nucleus contains deltas
  theEventInfo.deltasInside = incl_target->containsDeltas();

  // Take care of any remaining deltas

  if(theEventInfo.deltasInside){
    G4INCL::ParticleList const &inside = incl_target->getStore()->getOutgoingParticles();
  }
  //theEventInfo.forcedDeltasOutside = incl_target->decayOutgoingDeltas();   //FIXME: leave resonances to pythia
  theEventInfo.forcedDeltasInside = this->decayInsideDeltas(evrec, finalState);

  // Take care of any remaining etas, omegas, neutral Sigmas and/or neutral kaons
  double timeThreshold=theConfig->getDecayTimeThreshold();
  //theEventInfo.forcedPionResonancesOutside = incl_target->decayOutgoingPionResonances(timeThreshold);  //FIXME: leave resonances to pythia
  //incl_target->decayOutgoingSigmaZero(timeThreshold); //FIXME: leave resonances to pythia
  //incl_target->decayOutgoingNeutralKaon(); //FIXME: leave resonances to pythia

  //this->emitInsidePions(evrec, finalState); // FIXME: is it valid?

  // Apply Coulomb distortion, if appropriate
  // Note that this will apply Coulomb distortion also on pions emitted by
  // unphysical remnants (see decayInsideDeltas). This is at variance with
  // what INCL4.6 does, but these events are (should be!) so rare that
  // whatever we do doesn't (shouldn't!) make any noticeable difference.
  // CoulombDistortion::distortOut(incl_target->getStore()->getOutgoingParticles(), incl_target);

  // If the normal cascade predicted complete fusion, use the tabulated
  // masses to compute the excitation energy, the recoil, etc.
  //if(incl_target->getStore()->getOutgoingParticles().size()==0
  //    && (!incl_target->getProjectileRemnant()
  //	|| incl_target->getProjectileRemnant()->getParticles().size()==0)) {
  if( false && (!incl_target->getProjectileRemnant()
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

  theEventInfo.clusterDecay = this->decayOutgoingClusters(evrec, finalState) || this->decayMe(evrec, finalState); 
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

int INCLCascadeIntranuke::INCLPDG_to_GHEPPDG(int pdg, int A, int Z, int S) const{
  //  TParticlePDG * p = PDGLibrary::Instance()->Find(pdg);
  TDatabasePDG * fDatabasePDG = TDatabasePDG::Instance();
  TParticlePDG * p = fDatabasePDG->GetParticle(pdg);
  if(!p){
    if(A != 0){
      return genie::pdg::IonPdgCode( A , Z, std::abs(S), 0 );
    }
    else 
      return 0;
  }
  else
    return pdg;
}

void INCLCascadeIntranuke::fillFinalState(GHepRecord * evrec, G4INCL::FinalState * finalState) const{

  // neutrino interaction info
  const ProcessInfo & proc_info = evrec->Summary()->ProcInfo();

  // FIXME: start with the NC process
  INCLNucleus *incl_nucleus = INCLNucleus::Instance();

  double TLab;
  temfin = 29.8 * std::pow(incl_target->getA(), 0.16);

  TObjArrayIter piter(evrec);
  GHepParticle * p = nullptr;
  std::vector<G4INCL::GENIEParticleRecord> eventRecord;
  eventRecord.clear();
  while ( (p = (GHepParticle *) piter.Next() ) ) {
    // the code of the particles in primary neutrino interaction
    G4INCL::GENIERecordCode recordCode;
    if(eventRecord.size() == evrec->ProbePosition())
      recordCode = G4INCL::kProbe;
    else if(eventRecord.size() == evrec->HitNucleonPosition())
      recordCode = G4INCL::kHitNucleon;
    else if(eventRecord.size() == evrec->FinalStatePrimaryLeptonPosition())
      recordCode = G4INCL::kFinalStateLepton;
    else if(eventRecord.size() == evrec->TargetNucleusPosition())
      recordCode = G4INCL::kTarget;
    else if(eventRecord.size() == evrec->RemnantNucleusPosition())
      recordCode = G4INCL::kRemnant;
    else
      recordCode = G4INCL::kUnknown;

    eventRecord.emplace_back(p, int(proc_info.ScatteringTypeId()), recordCode);
  }


  G4INCL::IAvatar *avatar;
  if(proc_info.IsMEC()){
    avatar = new G4INCL::GENIEAvatar(0, incl_nucleus->getHitNNCluster(), incl_nucleus->getNuclues(), &eventRecord);
    avatar->fillFinalState(finalState);
  }
  else{
    avatar = new G4INCL::GENIEAvatar(0, incl_nucleus->getHitParticle(), incl_nucleus->getNuclues(), &eventRecord);
    avatar->fillFinalState(finalState);
  }

  // update the event record after INCL postInteraction
  //
  //

  TObjArrayIter piter1(evrec);
  GHepParticle * p1 = nullptr;
  auto er = eventRecord.begin();
  int idx =0;
  while ( (p1 = (GHepParticle *) piter1.Next() ) ) {
    TLorentzVector *p4 = p1->P4();
    p4->SetPx(er->P3().getX()/1000.);
    p4->SetPy(er->P3().getY()/1000.);
    p4->SetPz(er->P3().getZ()/1000.);
    p4->SetE(std::sqrt(er->P3().mag2() + er->Mass()*er->Mass())/1000.);
    tempFinalState.emplace_back(er->ID(), er->Pdg(), er->FirstMother(), idx++);
    er++;
  }
  //evrec->Print(std::cout);
  delete avatar;
  return;

}

G4INCL::ParticleType INCLCascadeIntranuke::PDG_to_INCLType(int pdg) const {
  switch(pdg){
    case 2212: return G4INCL::Proton;
    case 2112: return G4INCL::Neutron;
    case 211: return G4INCL::PiPlus;
    case -211: return G4INCL::PiMinus;
    case 111: return G4INCL::PiZero;
    case 2224: return G4INCL::DeltaPlusPlus;
    case 2214: return G4INCL::DeltaPlus;
    case 2114: return G4INCL::DeltaZero;
    case 1114: return G4INCL::DeltaMinus;
    case 221: return G4INCL::Eta;
    case 223: return G4INCL::Omega;
    case 331: return G4INCL::EtaPrime;
    case 22: return G4INCL::Photon;
    case 3122: return G4INCL::Lambda;
    case 3222: return G4INCL::SigmaPlus;
    case 3212: return G4INCL::SigmaZero;
    case 3112: return G4INCL::SigmaMinus;
    case -2212: return G4INCL::antiProton;
    case 3312: return G4INCL::XiMinus;
    case 3322: return G4INCL::XiZero;
    case -2112: return G4INCL::antiNeutron;
    case -3122: return G4INCL::antiLambda;
    case -3222: return G4INCL::antiSigmaPlus;
    case -3212: return G4INCL::antiSigmaZero;
    case -3112: return G4INCL::antiSigmaMinus;
    case -3312: return G4INCL::antiXiMinus;
    case -3322: return G4INCL::antiXiZero;
    case 321: return G4INCL::KPlus;
    case 311: return G4INCL::KZero;
    case -311: return G4INCL::KZeroBar;
    case -321: return G4INCL::KMinus;
    case 310: return G4INCL::KShort;
    case 130: return G4INCL::KLong;
    default:
              return G4INCL::UnknownParticle;
  }

}

void INCLCascadeIntranuke::DecayResonance(GHepRecord *evrec) const{
  if(evrec->Summary()->ProcInfo().IsResonant()){
    TObjArrayIter piter(evrec);
    GHepParticle * p = nullptr;
    bool decay_flag = false;
    while ( (p = (GHepParticle *) piter.Next() ) ) {
      if(p->Status() == kIStPreDecayResonantState && p->FirstDaughter() == -1 && PDG_to_INCLType(p->Pdg()) == G4INCL::UnknownParticle){
        decay_flag = true;
      }
    }
    if(decay_flag){
      fResonanceDecayer->ProcessEventRecord(evrec);
    }
  }
}

#include "INCLPostCascade.icc"

#endif // __GENIE_INCL_ENABLED__
