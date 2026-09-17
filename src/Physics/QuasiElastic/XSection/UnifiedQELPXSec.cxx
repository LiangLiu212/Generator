//____________________________________________________________________________
/*
 Copyright (c) 2003-2026, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Noah Steinberg <nsteinbe \at fnal.gov>
         Steven Gardiner <gardiner \at fnal.gov>
         Liang Liu <liangliu \at fnal.gov>
         Fermi National Acclerator Laboratory

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "TMath.h"
#include "TH2D.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

#include "UnifiedQELPXSec.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"
#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/NuclearState/SpectralFunc.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/Numerical/RandomGen.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//____________________________________________________________________________
UnifiedQELPXSec::UnifiedQELPXSec() :
XSecAlgorithmI("genie::UnifiedQELPXSec")
{

}
//____________________________________________________________________________
UnifiedQELPXSec::UnifiedQELPXSec(string config) :
XSecAlgorithmI("genie::UnifiedQELPXSec", config)
{

}
//____________________________________________________________________________
UnifiedQELPXSec::~UnifiedQELPXSec()
{

}
//____________________________________________________________________________
double UnifiedQELPXSec::XSec(const Interaction* interaction,
  KinePhaseSpace_t kps) const
{
  if ( !this->ValidProcess(interaction) ) return 0.;
  if ( !this->ValidKinematics(interaction) ) return 0.;

  // Get kinematics and init-state parameters
  const Kinematics&   kinematics = interaction->Kine();
  const InitialState& init_state = interaction->InitState();
  const Target& target = init_state.Tgt();

  TLorentzVector* temp_probeP4 = init_state.GetProbeP4( kRfLab );
  TLorentzVector probeP4 = *temp_probeP4;
  delete temp_probeP4;
  double E_probe = probeP4.E();

  TLorentzVector p4Ni = target.HitNucP4();
  double mNi = target.HitNucMass(); // on-shell initial hit nucleon mass
  double pNi = p4Ni.P();
  double E_NiOnShell = std::sqrt(pNi*pNi + mNi*mNi);
  TLorentzVector p4NiOnShell = TLorentzVector(p4Ni.Vect(), E_NiOnShell);
  double epsilon_B = E_NiOnShell - p4Ni.E();

  TLorentzVector lepP4 = kinematics.FSLeptonP4();
  double E_lep = lepP4.E();

  TLorentzVector p4Nf = kinematics.HadSystP4();
  double E_Nf = p4Nf.E();
  
  // Include phase space factors in cross section
  double xsec = fXSecScale / (E_lep * E_probe * E_NiOnShell * E_Nf) / 32. / kPi / kPi;

  // If we're dealing with a nuclear target, then apply Pauli blocking as
  // needed  
  if (fDoPauliBlocking && target.IsNucleus() && !interaction->TestBit(kIAssumeFreeNucleon) ) {
    double kF = fPauliBlocker->GetFermiMomentum(target, interaction->RecoilNucleonPdg(), target.HitNucPosition());
    if ( p4Nf.P() < kF ) {return 0.;}
  }

  // Scale cross section by the number of active nucleons
  int hit_nuc_pdg = target.HitNucPdg();
  // Number of active nucleons in the target
  int num_active = pdg::IsProton(hit_nuc_pdg) ? target.Z() : target.N();

  xsec *= num_active;

  // Spectator nucleons for the one-body / two-body current interference
  // (EM only for now)
  bool do_intf = fDoIntf && interaction->ProcInfo().IsEM();
  std::vector<IAOneTwoBodyInterferenceTensor::Spectator> spectators;
  if ( do_intf ) spectators = this->Spectators( target, p4Ni );
  if ( spectators.empty() ) do_intf = false;

  // Do we need to rotate so that \vec{q} || \vec{z}?
  if (fDoqAlongZ) {
    std::vector<TLorentzVector> otherp4 {p4Ni, p4NiOnShell, p4Nf};
    for ( const auto& sp : spectators ) otherp4.emplace_back( sp.p3, 0. );
    genie::utils::Rotate_qvec_alongZ(probeP4, lepP4, otherp4);
    p4Ni = otherp4[0];
    p4NiOnShell = otherp4[1];
    p4Nf = otherp4[2];
    for ( size_t s = 0; s < spectators.size(); ++s ) spectators[s].p3 = otherp4[3 + s].Vect();
  }

  // Compute form factors using Q2tilde (the effective Q2 value after
  // binding energy corrections)
  TLorentzVector qP4 = probeP4 - lepP4;
  TLorentzVector qTildeP4 = qP4;
  qTildeP4.SetE( qP4.E() - epsilon_B );

  double Q2 = -1. * qP4.M2();
  double Q2tilde = -1. * qTildeP4.M2();

  double Q2min = genie::controls::kMinQ2Limit; // CC/NC limit
  //std::cout << Q2min << "\n";
  if( interaction->ProcInfo().IsEM() ) Q2min = genie::utils::kinematics::electromagnetic::kMinQ2Limit; // EM limit
  // Make sure Q2 is physical
  if ( Q2 < Q2min ) {
    return 0.;
  }

  // Get the correct couplings and form factors model for the current
  // interaction
  double coupling_factor = 1.;
  const genie::ProcessInfo& proc_info = interaction->ProcInfo();
  if ( proc_info.IsWeakCC() ) {
    coupling_factor = kGF2 * fCos8c2;
    fFormFactors.SetModel( fCCFormFactorsModel );
  }
  else if ( proc_info.IsWeakNC() ) {
    coupling_factor = kGF2;
    fFormFactors.SetModel( fNCFormFactorsModel );
  }
  else if ( proc_info.IsEM() ) {
    coupling_factor = 32. * kPi * kPi * kAem2 / ( Q2*Q2 ) ;
    fFormFactors.SetModel( fEMFormFactorsModel );
  }
  else {
    LOG("UnifiedQE", pERROR) << "Unrecognized process type encountered"
      << " in genie::UnifiedQELPXSec::XSec()";
    return 0.;
  }
 
  // Apply the coupling factor to the differential cross section
  xsec *= coupling_factor;

  // For a bound hit nucleon, the corrected energy transfer qTilde0 needs to be
  // positive to enforce energy conservation
  if ( qTildeP4.E() <= 0. && interaction->InitState().Tgt().IsNucleus()
    && !interaction->TestBit(kIAssumeFreeNucleon) ) return 0.;

  // Set Q2 to Q2tilde (or keep the true Q2, see FormFactorsAtQ2Tilde) while
  // computing form factors
  double Q2ff = fFFAtQ2Tilde ? Q2tilde : Q2;
  interaction->KinePtr()->SetQ2( Q2ff );
  // Evaluate the form factors
  fFormFactors.Calculate( interaction );
  if ( do_intf ) fELFormFactors.Calculate( interaction );

  // Now that we've calculated them, store the true Q2 value
  interaction->KinePtr()->SetQ2( Q2 );

  // Compute the leptonic tensor
  // Note we have to pass a bool to LeptonTensor to use the SF conventions
  InteractionType_t type = interaction->ProcInfo().InteractionTypeId();
  LeptonTensor L_munu( probeP4, lepP4, init_state.ProbePdg(), type, true);

  // For CC use mean of proton/neutron mass
  // works for NC and EM as well 
  double xmn = ( mNi + interaction->RecoilNucleon()->Mass() ) / 2.;

  // Make a generic Rank2LorentzTensor 
  // object for the hadronic tensor
  std::shared_ptr<Rank2LorentzTensor> ATilde_munu;

  // call the hadron-tensor; now we have Noemi-hadron-tensor and Noemi-hadron-tensor-cc
  if(fTensorModel.find("Noemi-hadron-tensor") != std::string::npos) {
    ATilde_munu = std::make_shared<IASingleNucleonTensor>(qP4.E(), xmn, p4Ni, p4Nf, fFormFactors, fTensorModel);
  }
  else{
    LOG("UnifiedQE", pFATAL) << "Wrong setup of hadron tensor, the options should be Noemi-hadron-tensor or Noemi-hadron-tensor-cc";
    exit(1);
  }
      
  // Contract hadron and lepton tensors
  std::complex<double> contraction = L_munu * (*ATilde_munu);

  if ( std::abs(contraction.imag()) > kASmallNum ) {
    LOG("UnifiedQE", pWARN) << "Tensor contraction has nonvanishing imaginary part!";
  }
  
  double LA = contraction.real();

  // Add the interference between the one-body and the two-body currents. It
  // leads to the same single-nucleon knock-out final state, so it is simply
  // part of the quasielastic cross section.
  if ( do_intf ) {
    genie::twobody_currents_sf::ModelParams par = fIntfPar;
    par.xmn = xmn / genie::units::MeV;

    // N-Delta transition and pion form factors (couplings are in the
    // cross section prefactor, as for the one-body current)
    genie::twobody_currents_sf::FormFactors ff;
    ff.f1  = fFormFactors.F1V();
    ff.f2  = fFormFactors.xiF2V();
    ff.fa  = 0.;
    ff.fap = 0.;
    ff.fpiem = fELFormFactors.Gep() - fELFormFactors.Gen();
    double dipole2 = std::pow(1. + Q2ff / fMV2, 2);
    ff.cv3 = fCV3Norm / dipole2 / (1. + Q2ff / 4. / fMV2) * std::sqrt(1.5);
    ff.cv4 = fCV4Norm / dipole2 / (1. + Q2ff / 4. / fMV2) * std::sqrt(1.5);
    ff.cv5 = fCV5Norm / dipole2 / (1. + Q2ff / 0.776 / fMV2) * std::sqrt(1.5);
    ff.ca5 = 0.;

    IAOneTwoBodyInterferenceTensor A12_munu(par, ff, p4Ni, p4Nf, qP4,
      hit_nuc_pdg, interaction->RecoilNucleonPdg(), spectators, false, false);
    std::complex<double> contraction12 = L_munu * A12_munu;
    LA += contraction12.real();

    // The spectator sum is a Monte Carlo estimate, so an individual sample
    // can (rarely) drive the total below zero
    if ( LA < 0. ) {
      LOG("UnifiedQE", pINFO) << "One-body + interference contraction is"
        << " negative (" << LA << "), setting it to zero";
      LA = 0.;
    }
  }

  // Apply the tensor contraction to the cross section
  xsec *= LA;
  
  // Multiply by the analytic solution of the energy-conserving delta function
  // used by the kPSQELEvGen phase space. 
  xsec *= genie::utils::EnergyDeltaFunctionSolutionQEL( *interaction );

  // Check whether variable tranformation is needed
  if ( kps != kPSQELEvGen ) {
    // Compute the appropriate Jacobian for transformation to the requested
    // phase space
    double J = utils::kinematics::Jacobian(interaction, kPSQELEvGen, kps);
    xsec *= J;
  }

  return xsec;
}
//____________________________________________________________________________
const std::vector<IAOneTwoBodyInterferenceTensor::Spectator>&
UnifiedQELPXSec::Spectators(const Target& target, const TLorentzVector& p4Ni) const
{
  int tgt_pdg = target.Pdg();
  int hit_nuc_pdg = target.HitNucPdg();
  if ( tgt_pdg == fCachedTgtPdg && hit_nuc_pdg == fCachedHitNucPdg
    && p4Ni == fCachedP4Ni ) return fCachedSpectators;

  fCachedTgtPdg = tgt_pdg;
  fCachedHitNucPdg = hit_nuc_pdg;
  fCachedP4Ni = p4Ni;
  fCachedSpectators.clear();

  // The interference needs the mean-field part of the spectral function
  // of this nucleus. Without it there is no interference term.
  std::string tgt_str = std::to_string( tgt_pdg );
  const Registry& mf_config = fMFSpectralFunc->GetConfig();
  if ( !mf_config.Exists("SpectFuncTable@Pdg=" + tgt_str + "_" + std::to_string(kPdgProton))
    || !mf_config.Exists("SpectFuncTable@Pdg=" + tgt_str + "_" + std::to_string(kPdgNeutron)) )
  {
    return fCachedSpectators;
  }

  // The hit nucleon was sampled from the complete spectral function. Only
  // its mean-field part contributes: weight by S_MF / S_tot at the sampled
  // momentum and removal energy
  double pNi = p4Ni.P();
  double E_rmv = target.HitNucMass() - p4Ni.E();
  if ( E_rmv < 0. ) return fCachedSpectators;
  double prob_tot = fTotSpectralFunc->Prob( pNi, E_rmv, target );
  double prob_mf  = fMFSpectralFunc->Prob( pNi, E_rmv, target );
  if ( prob_tot <= 0. || prob_mf <= 0. ) return fCachedSpectators;
  double mf_fraction = prob_mf / prob_tot;

  // Inverse volume rho / A of nuclear matter at the Fermi momentum kF (MeV^3)
  double kF = 0.;
  RgKey kf_key = "TwoBodyIntf-FermiMomentum@Pdg=" + tgt_str;
  if ( this->GetConfig().Exists(kf_key) ) this->GetParam( kf_key, kF );
  else {
    std::string kf_table_name;
    this->GetParam( "FermiMomentumTable", kf_table_name );
    const FermiMomentumTable* kf_table
      = FermiMomentumTablePool::Instance()->GetTable( kf_table_name );
    kF = kf_table->FindClosestKF( tgt_pdg, kPdgProton );
  }
  kF /= genie::units::MeV;
  double rho_over_A = std::pow(kF, 3) / (1.5 * kPi2) / target.A();

  // Temporarily replace ROOT's gRandom RNG with GENIE's, as SpectralFunc does
  TRandom* old_gRandom = gRandom;
  RandomGen* rnd = RandomGen::Instance();
  gRandom = &rnd->RndGen();

  // If protons and neutrons share their mean-field table, one spectator
  // momentum serves both isospins (the two-body operators do not depend on
  // the spectator isospin, so this halves the cost). Otherwise each isospin
  // gets its own sample.
  std::string mf_table_p = mf_config.GetString(
    "SpectFuncTable@Pdg=" + tgt_str + "_" + std::to_string(kPdgProton) );
  std::string mf_table_n = mf_config.GetString(
    "SpectFuncTable@Pdg=" + tgt_str + "_" + std::to_string(kPdgNeutron) );
  bool shared_table = ( mf_table_p == mf_table_n );

  const int spect_pdgs[2] = { kPdgProton, kPdgNeutron };
  double weights[2];
  TH2D* sf_mf[2];
  for ( int is = 0; is < 2; ++is ) {
    Target spect_tgt( target );
    spect_tgt.SetHitNucPdg( spect_pdgs[is] );
    sf_mf[is] = fMFSpectralFunc->SelectSpectralFunction( spect_tgt );

    // Number of mean-field nucleons of this species (the histogram is
    // normalised to the mean-field fraction of one nucleon)
    int num_nuc = pdg::IsProton( spect_pdgs[is] ) ? target.Z() : target.N();
    double num_mf = num_nuc * sf_mf[is]->Integral();
    weights[is] = mf_fraction * rho_over_A * num_mf / fNumSpectators;
  }

  for ( int is = 0; is < (shared_table ? 1 : 2); ++is ) {
    for ( int s = 0; s < fNumSpectators; ++s ) {
      double p2 = 0., E2_rmv = 0.;
      sf_mf[is]->GetRandom2( p2, E2_rmv );
      double costheta = -1. + 2. * rnd->RndGen().Rndm();
      double sintheta = TMath::Sqrt(1. - costheta*costheta);
      double phi = 2. * kPi * rnd->RndGen().Rndm();

      IAOneTwoBodyInterferenceTensor::Spectator sp;
      sp.p3.SetXYZ( p2*sintheta*TMath::Cos(phi), p2*sintheta*TMath::Sin(phi), p2*costheta );
      sp.weight_p = ( shared_table || is == 0 ) ? weights[0] : 0.;
      sp.weight_n = ( shared_table || is == 1 ) ? weights[1] : 0.;
      fCachedSpectators.push_back( sp );
    }
  }

  gRandom = old_gRandom;

  return fCachedSpectators;
}
//____________________________________________________________________________
double UnifiedQELPXSec::Integral(const Interaction* in) const
{
  // Intended for use with genie::NewQELXSec, which is smart
  // enough to handle free nucleon vs. nuclear targets, different
  // nuclear models, etc.
  return fXSecIntegrator->Integrate(this, in);
}
//____________________________________________________________________________
bool UnifiedQELPXSec::ValidProcess(const Interaction* interaction) const
{
  if ( interaction->TestBit(kISkipProcessChk) ) return true;

  const InitialState& init_state = interaction->InitState();
  const ProcessInfo&  proc_info  = interaction->ProcInfo();

  // Calculation is only appropriate for complex nuclear targets,
  // not free nucleons.
  if ( !init_state.Tgt().IsNucleus() ) return false;

  if ( !proc_info.IsQuasiElastic() ) return false;

  if ( proc_info.IsEM() || proc_info.IsWeakNC() ) return true;

  // For weak CC interactions, check that the hit nucleon "works"
  else if ( proc_info.IsWeakCC() ) {
    int nucleon_pdg = init_state.Tgt().HitNucPdg();
    int probe_pdg = init_state.ProbePdg();

    if ( pdg::IsNeutron(nucleon_pdg) && pdg::IsNeutrino(probe_pdg) ) {
      return true;
    }
    else if ( pdg::IsProton(nucleon_pdg) && pdg::IsAntiNeutrino(probe_pdg) ) {
      return true;
    }
    else return false;
  }

  return false;
}
//____________________________________________________________________________
void UnifiedQELPXSec::Configure(const Registry& config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void UnifiedQELPXSec::Configure(std::string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void UnifiedQELPXSec::LoadConfig(void)
{
  double thc;
  GetParam( "CabibboAngle", thc ) ;
  fCos8c2 = TMath::Power(TMath::Cos(thc), 2);

  // cross section scaling factor
  GetParam( "QEL-CC-XSecScale", fXSecScale ) ;

  // load QEL form factors models
  fCCFormFactorsModel = dynamic_cast<const QELFormFactorsModelI*>(
    this->SubAlg("CCFormFactorsAlg") );
  assert( fCCFormFactorsModel );

  fNCFormFactorsModel = dynamic_cast<const QELFormFactorsModelI*>(
    this->SubAlg("NCFormFactorsAlg") );
  assert( fNCFormFactorsModel );

  fEMFormFactorsModel = dynamic_cast<const QELFormFactorsModelI*>(
    this->SubAlg("EMFormFactorsAlg") );
  assert( fEMFormFactorsModel );

  // Attach CC model for now. This will be updated later.
  fFormFactors.SetModel( fCCFormFactorsModel );

  // Pick a hadron tensor model
  GetParamDef("TensorModel", fTensorModel, std::string("Noemi-hadron-tensor"));

  // load xsec integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI*>(
    this->SubAlg("XSec-Integrator") );
  assert(fXSecIntegrator);

  // get nuclear model
  fNuclModel = dynamic_cast<const NuclearModelI*>(
    this->SubAlg("IntegralNuclearModel") );
  assert(fNuclModel);

  // get the algorithm ID for the PauliBlocker
  RgAlg pauliBlockID;
  GetParamDef( "PauliBlockerAlg", pauliBlockID,
    RgAlg("genie::PauliBlocker", "Default") );
  AlgId pbID = AlgId( pauliBlockID );

  AlgFactory* algf = AlgFactory::Instance();
  fPauliBlocker = dynamic_cast<const PauliBlocker*>(
    algf->GetAlgorithm(pbID) );
  assert( fPauliBlocker );

  // Decide whether or not it should be used in XSec()
  GetParamDef( "DoPauliBlocking", fDoPauliBlocking, true );
  
  // Decide whether or not we rotate so that q is along z
  GetParamDef( "DoRotate_qAlong_z", fDoqAlongZ, false );

  // Evaluate the form factors at Q2tilde (default) or at the true Q2
  GetParamDef( "FormFactorsAtQ2Tilde", fFFAtQ2Tilde, true );

  // Interference between the one-body and the two-body currents
  GetParamDef( "IncludeOneTwoBodyInterference", fDoIntf, false );
  fTotSpectralFunc = 0;
  fMFSpectralFunc = 0;
  fELFormFactorsModel = 0;
  fCachedTgtPdg = 0;
  fCachedHitNucPdg = 0;
  fCachedSpectators.clear();
  if ( fDoIntf ) {
    fTotSpectralFunc = dynamic_cast<const SpectralFunc*>( fNuclModel );
    fMFSpectralFunc = dynamic_cast<const SpectralFunc*>(
      this->SubAlg("TwoBodyIntf-MeanFieldSpectralFunc") );
    if ( !fTotSpectralFunc || !fMFSpectralFunc ) {
      LOG("UnifiedQE", pFATAL) << "The one-body / two-body current interference"
        << " needs genie::SpectralFunc algorithms for both IntegralNuclearModel"
        << " and TwoBodyIntf-MeanFieldSpectralFunc";
      exit(1);
    }

    fELFormFactorsModel = dynamic_cast<const ELFormFactorsModelI*>(
      this->SubAlg("ElasticFormFactorsModel") );
    assert( fELFormFactorsModel );
    fELFormFactors.SetModel( fELFormFactorsModel );

    // Couplings and cutoffs of the two-body currents. GENIE units (GeV) in
    // the configuration, MeV in twobody_currents_sf.
    double lpi, lpind;
    GetParamDef( "TwoBodyIntf-fPiNDelta", fIntfPar.fpind, 0.54 );
    GetParamDef( "TwoBodyIntf-fStar", fIntfPar.fstar, 2.13 );
    GetParamDef( "TwoBodyIntf-fPiNN2", fIntfPar.fpinn2, 0.081 * 4. * kPi );
    GetParamDef( "TwoBodyIntf-gA", fIntfPar.ga, 1.26 );
    GetParamDef( "TwoBodyIntf-LambdaPi", lpi, 1.300 );
    GetParamDef( "TwoBodyIntf-LambdaPiNDelta", lpind, 1.150 );
    fIntfPar.lpi = lpi / genie::units::MeV;
    fIntfPar.lpind = lpind / genie::units::MeV;

    PDGLibrary* pdglib = PDGLibrary::Instance();
    fIntfPar.xmd = pdglib->Find( kPdgP33m1232_DeltaP )->Mass() / genie::units::MeV;
    fIntfPar.xmpi = kPionMass / genie::units::MeV;
    fIntfPar.xmrho = pdglib->Find( kPdgRho0 )->Mass() / genie::units::MeV;
    fIntfPar.xmn = kNucleonMass / genie::units::MeV; // reset per event

    GetParamDef( "TwoBodyIntf-C3V", fCV3Norm, 2.13 );
    GetParamDef( "TwoBodyIntf-C4V", fCV4Norm, -1.15 );
    GetParamDef( "TwoBodyIntf-C5V", fCV5Norm, 0.48 );
    GetParamDef( "TwoBodyIntf-MV2", fMV2, 0.71 );

    GetParamDef( "TwoBodyIntf-NumSpectators", fNumSpectators, 1 );
    GetParamDef( "TwoBodyIntf-AchillesC4VC5V", fIntfPar.c4c5_broadcast, false );
  }

}

