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

#include <cstdlib>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TParticlePDG.h>
#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Physics/NuclearState/NucleusGenTraditional.h"

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
#include "Physics/NuclearState/INCLNucleus.h"
#include "Physics/Common/VertexGenerator.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
NucleusGenTraditional::NucleusGenTraditional() :
NucleusGenI("genie::NucleusGenTraditional")
{

}
//___________________________________________________________________________
NucleusGenTraditional::NucleusGenTraditional(string config) :
NucleusGenI("genie::NucleusGenTraditional", config)
{

}
//___________________________________________________________________________
NucleusGenTraditional::~NucleusGenTraditional()
{

}

//___________________________________________________________________________
void NucleusGenTraditional::ProcessEventRecord(GHepRecord * evrec) const
{
  // skip if not a nuclear target
  if(! evrec->Summary()->InitState().Tgt().IsNucleus()) return;
  fVertexGenerator->ProcessEventRecord(evrec);
  fFermiMover->ProcessEventRecord(evrec);
}

void NucleusGenTraditional::GenerateVertex(GHepRecord * evrec) const{
  // This is function will be used in QEL-CC channel
  // skip if not a nuclear target
  LOG("NucleusGenTraditional", pNOTICE) << "new events";
  if(! evrec->Summary()->InitState().Tgt().IsNucleus()) return;
  fVertexGenerator->ProcessEventRecord(evrec);
}

void NucleusGenTraditional::GenerateCluster(GHepRecord * evrec) const{
  // This is function will be used in MEC channela
  // first, generate vertex for hit cluster
  if(! evrec->Summary()->InitState().Tgt().IsNucleus()) return;
  fVertexGenerator->ProcessEventRecord(evrec);

  GHepParticle * target_nucleus = evrec->TargetNucleus();
  assert(target_nucleus);
  GHepParticle * nucleon_cluster = evrec->HitNucleon();
  assert(nucleon_cluster);
  GHepParticle * remnant_nucleus = evrec->RemnantNucleus();
  assert(remnant_nucleus);

  // generate a Fermi momentum for each nucleon

  Target tgt(target_nucleus->Pdg());

  bool allowdup = true;
  PDGCodeList pdgv(allowdup);

  if(nucleon_cluster->Pdg() == kPdgClusterNN) {
     pdgv.push_back(kPdgNeutron);
     pdgv.push_back(kPdgNeutron);
  }
  else
  if(nucleon_cluster->Pdg() == kPdgClusterNP) {
     pdgv.push_back(kPdgNeutron);
     pdgv.push_back(kPdgProton);
  }
  else
  if(nucleon_cluster->Pdg() == kPdgClusterPP) {
     pdgv.push_back(kPdgProton);
     pdgv.push_back(kPdgProton);
  }
  else
  {
     LOG("NucleusGenTraditional", pERROR)
        << "Unknown di-nucleon cluster PDG code (" << nucleon_cluster->Pdg() << ")";
     exit(1);
  }

  assert(pdgv.size()==2);

  TVector3 p3a, p3b;
  double removalenergy1, removalenergy2;
  fNuclModel->GenerateCluster(tgt, pdgv, &p3a, &p3b, &removalenergy1, &removalenergy2);
  LOG("NucleusGenTraditional", pINFO)
    << "1st nucleon (code = " << pdgv[0] << ") generated momentum: ("
    << p3a.Px() << ", " << p3a.Py() << ", " << p3a.Pz() << "), "
    << "|p| = " << p3a.Mag();
  LOG("NucleusGenTraditional", pINFO)
    << "2nd nucleon (code = " << pdgv[1] << ") generated momentum: ("
    << p3b.Px() << ", " << p3b.Py() << ", " << p3b.Pz() << "), "
    << "|p| = " << p3b.Mag();

  TVector3 p3 = p3a + p3b;
  double M2n = PDGLibrary::Instance()->Find(nucleon_cluster->Pdg())-> Mass(); // nucleon cluster mass
  // nucleon cluster energy
  double EN = TMath::Sqrt(p3.Mag2() + M2n*M2n);

  TLorentzVector p4nclust   (   p3.Px(),    p3.Py(),    p3.Pz(),  EN   );
  nucleon_cluster->SetMomentum(p4nclust);

}
//___________________________________________________________________________

void NucleusGenTraditional::BindHitNucleon() const {
}

void NucleusGenTraditional::BindHitNucleon(Interaction& interaction, double& Eb, QELEvGen_BindingMode_t hitNucleonBindingMode) const {

  Target* tgt = interaction.InitState().TgtPtr();
  TLorentzVector* p4Ni = tgt->HitNucP4Ptr();

  // Initial nucleon 3-momentum (lab frame)
  TVector3 p3Ni = fNuclModel->Momentum3();

  // Look up the (on-shell) mass of the initial nucleon
  TDatabasePDG* tb = TDatabasePDG::Instance();
  double mNi = tb->GetParticle( tgt->HitNucPdg() )->Mass();

  // Set the (possibly off-shell) initial nucleon energy based on
  // the selected binding energy mode. Always put the initial nucleon
  // on shell if it is not part of a composite nucleus
  double ENi = 0.;
  if ( tgt->IsNucleus() && hitNucleonBindingMode != genie::kOnShell ) {

    // For a nuclear target with a bound initial struck nucleon, take binding
    // energy effects and Pauli blocking into account when computing QE
    // differential cross sections
    interaction.ResetBit( kIAssumeFreeNucleon );

    // Initial nucleus mass
    double Mi = tgt->Mass();

    // Final nucleus mass
    double Mf = 0.;

    // If we use the removal energy reported by the nuclear
    // model, then it implies a certain value for the final
    // nucleus mass
    if ( hitNucleonBindingMode == genie::kUseNuclearModel ) {
      Eb = fNuclModel->RemovalEnergy();
      // This equation is the definition that we assume
      // here for the "removal energy" (Eb) returned by the
      // nuclear model. It matches GENIE's convention for
      // the Bodek/Ritchie Fermi gas model.
      Mf = Mi + Eb - mNi;
    }
    // We can also assume that the final nucleus is in its
    // ground state. In this case, we can just look up its
    // mass from the standard table. This implies a particular
    // binding energy for the hit nucleon.
    else if ( hitNucleonBindingMode == genie::kUseGroundStateRemnant ) {
      // Determine the mass and proton numbers for the remnant nucleus
      int Af = tgt->A() - 1;
      int Zf = tgt->Z();
      if ( genie::pdg::IsProton( tgt->HitNucPdg()) ) --Zf;
      Mf = genie::PDGLibrary::Instance()->Find( genie::pdg::IonPdgCode(Af, Zf) )->Mass();

      // Deduce the binding energy from the final nucleus mass
      Eb = Mf - Mi + mNi;
    }
    // A third option is to assign the binding energy and final nuclear
    // mass in a way that is equivalent to the original Valencia model
    // treatment
    else if ( hitNucleonBindingMode == genie::kValenciaStyleQValue ) {
      // Compute the Q-value needed for the ground-state-to-ground-state
      // transition without nucleon removal (see Sec. III B of
      // https://arxiv.org/abs/nucl-th/0408005). Note that this will be
      // zero for NC and EM scattering since the proton and nucleon numbers
      // are unchanged. Also get Q_LFG, the difference in Fermi energies
      // between the final and initial struck nucleon species. This will
      // likewise be zero for NC and EM interactions.
      const genie::ProcessInfo& info = interaction.ProcInfo();
      double Qvalue = 0.;
      double Q_LFG = 0.;
      if ( info.IsWeakCC() ) {
        // Without nucleon removal, the final nucleon number remains the same
        int Af = tgt->A();
        // For CC interactions, the proton number will change by one. Choose
        // the right change by checking whether we're working with a neutrino
        // or an antineutrino.
        int Zf = tgt->Z();
        int probe_pdg = interaction.InitState().ProbePdg();
        if ( genie::pdg::IsAntiNeutrino(probe_pdg) ) {
          --Zf;
        }
        else if ( genie::pdg::IsNeutrino(probe_pdg) ) {
          ++Zf;
        }
        else {
          LOG( "QELEvent", pFATAL ) << "Unhandled probe PDG code " << probe_pdg
            << " encountered for a CC interaction in"
            << " genie::utils::BindHitNucleon()";
          gAbortingInErr = true;
          std::exit( 1 );
        }

        // We have what we need to get the Q-value. Get the final nuclear
        // mass (without nucleon removal)
        double mf_keep_nucleon = genie::PDGLibrary::Instance()
          ->Find( genie::pdg::IonPdgCode(Af, Zf) )->Mass();

        Qvalue = mf_keep_nucleon - Mi;

        // Get the Fermi energies for the initial and final nucleons. Include
        // the radial dependence if using the LFG.
        double hit_nucleon_radius = tgt->HitNucPosition();

        // Average of the proton and neutron masses. It may actually be better
        // to use the exact on-shell masses here. However, the original paper
        // uses the same value for protons and neutrons when computing the
        // difference in Fermi energies. For consistency with the Valencia
        // model publication, I'll do the same.
        const double mN = genie::constants::kNucleonMass;

        double kF_Ni = fNuclModel->LocalFermiMomentum( *tgt,
          tgt->HitNucPdg(), hit_nucleon_radius );
        double EFermi_Ni = std::sqrt( std::max(0., mN*mN + kF_Ni*kF_Ni) );

        double kF_Nf = fNuclModel->LocalFermiMomentum( *tgt,
          interaction.RecoilNucleonPdg(), hit_nucleon_radius );
        // (On-shell) final nucleon mass
        //double mNf = interaction.RecoilNucleon()->Mass();
        double EFermi_Nf = std::sqrt( std::max(0., mN*mN + kF_Nf*kF_Nf) );

        // The difference in Fermi energies is Q_LFG
        Q_LFG = EFermi_Nf - EFermi_Ni;
      }
      // Fail here for interactions that aren't one of EM/NC/CC. This is
      // intended to help avoid silently doing the wrong thing in the future.
      else if ( !info.IsEM() && !info.IsWeakNC() ) {
        LOG( "QELEvent", pFATAL ) << "Unhandled process type \"" << info
          << "\" encountered in genie::utils::BindHitNucleon()";
        gAbortingInErr = true;
        std::exit( 1 );
      }

      // On-shell total energy of the initial struck nucleon
      double ENi_OnShell = std::sqrt( std::max(0., mNi*mNi + p3Ni.Mag2()) );

      // Total energy of the remnant nucleus (with the hit nucleon removed)
      double Ef = Mi - ENi_OnShell + Qvalue - Q_LFG;

      // Mass and kinetic energy of the remnant nucleus (with the hit nucleon
      // removed)
      Mf = std::sqrt( std::max(0., Ef*Ef - p3Ni.Mag2()) );

      // Deduce the binding energy from the final nucleus mass
      Eb = Mf - Mi + mNi;

      LOG( "QELEvent", pDEBUG ) << "Qvalue = " << Qvalue
        << ", Q_LFG = " << Q_LFG << " at radius = " << tgt->HitNucPosition();
    }

    // The (lab-frame) off-shell initial nucleon energy is the difference
    // between the lab frame total energies of the initial and remnant nuclei
    ENi = Mi - std::sqrt( Mf*Mf + p3Ni.Mag2() );
  }
  else {
    // Keep the struck nucleon on shell either because
    // hitNucleonBindingMode == kOnShell or because
    // the target is a single nucleon
    ENi = std::sqrt( p3Ni.Mag2() + std::pow(mNi, 2) );
    Eb = 0.;

    // If we're dealing with a nuclear target but using the on-shell
    // binding mode, set the "assume free nucleon" flag. This turns off
    // Pauli blocking and the requirement that q0 > 0 in the QE cross section
    // models (an on-shell nucleon *is* a free nucleon)
    if ( tgt->IsNucleus() ) interaction.SetBit( kIAssumeFreeNucleon );
  }

  LOG( "QELEvent", pDEBUG ) << "Eb = " << Eb << ", pNi = " << p3Ni.Mag()
    << ", ENi = " << ENi;

  // Update the initial nucleon lab-frame 4-momentum in the interaction with
  // its current components
  p4Ni->SetVect( p3Ni );
  p4Ni->SetE( ENi );



} 

//___________________________________________________________________________

void NucleusGenTraditional::GenerateNucleon(Interaction* interaction, bool isRadius) const {
  // isRadius:
  // isRadius == 0 (false) put the nucleon in origin
  // isRadius == 1 (true) sampling the radius for nucleon
  Target* tgt = interaction->InitState().TgtPtr();
  if(isRadius){
    const VertexGenerator* vtx_gen = dynamic_cast<const VertexGenerator*>(fVertexGenerator);
    TVector3 vertex_pos = vtx_gen->GenerateVertex( interaction, tgt->A() );
    double radius = vertex_pos.Mag();
    tgt->SetHitNucPosition( radius );
    fNuclModel->GenerateNucleon(*tgt, radius);
  }
  else{
    tgt->SetHitNucPosition(0.);
    if ( tgt->IsNucleus() ){
      fNuclModel->GenerateNucleon(*tgt, 0.);
    }
    else{
      fNuclModel->SetRemovalEnergy(0.);
      interaction->SetBit( kIAssumeFreeNucleon );
    }
    fNuclModel->SetMomentum3( TVector3(0., 0., 0.) );
  }

}


//___________________________________________________________________________

// is function is a util function for QEL
bool NucleusGenTraditional::isRPValid(double r, double p, const Target & tgt) const {
  // TODO: document this, won't work for spectral functions
  double dummy_w = -1.;
  double prob = fNuclModel->Prob(p, dummy_w, tgt, r);
  return (prob > 0.);
}

//___________________________________________________________________________

void NucleusGenTraditional::SetHitNucleonOnShellMom(TVector3 p3) const {
  // Set the nucleon we're using to be upstream at max energy and unbound
  fNuclModel->SetMomentum3( p3 );
  fNuclModel->SetRemovalEnergy( 0. );
}

//___________________________________________________________________________
void NucleusGenTraditional::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NucleusGenTraditional::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NucleusGenTraditional::LoadConfig(void)
{
  bool is_configed = true;
  fFermiMover = nullptr;
  fVertexGenerator = nullptr;
  fFermiMover = dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("FermiMover"));
  if(!fFermiMover){
    is_configed = false;
    LOG("NucleusGenTraditional", pERROR) << "The SubAlg FermiMover is not configed";
  }
  fVertexGenerator = dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("VertexGenerator"));
  if(!fVertexGenerator){
    is_configed = false;
    LOG("NucleusGenTraditional", pERROR) << "The SubAlg VertexGenerator is not configed";
  }

  if(!is_configed){
    LOG("NucleusGenTraditional", pERROR) << "Configuration has failed.";
    exit(78);
  }

  RgKey nuclkey = "NuclearModel";
  fNuclModel = nullptr;
  fNuclModel = dynamic_cast<const NuclearModelI *>(fFermiMover->SubAlg(nuclkey));
  assert(fNuclModel);
}
//____________________________________________________________________________

