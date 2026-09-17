//____________________________________________________________________________
/*!

\class    genie::UnifiedQELPXSec

\brief    Computes quasielastic neutrino-nucleus differential cross sections
          using a contraction of leptonic and hadronic tensors. Intended for
          use with a spectral function nuclear model.
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\author   Steven Gardiner <gardiner \at fnal.gov>
          Liang Liu <liangliu \at fnal.gov>
          Fermi National Acclerator Laboratory

\created  March 25, 2019

\cpright  Copyright (c) 2003-2026, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _CBF_SPECTRAL_FUNC_CROSS_SECTION_H_
#define _CBF_SPECTRAL_FUNC_CROSS_SECTION_H_

#include <string>
#include <complex>
#include <vector>

#include "Math/IFunction.h"

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/QuasiElastic/XSection/LeptonTensor.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/NuclearState/PauliBlocker.h"
#include "Physics/HadronTensors/Rank2LorentzTensor.h"
#include "Physics/HadronTensors/NucleonTensor.h"
#include "Physics/HadronTensors/IASingleNucleonTensor.h"
#include "Physics/HadronTensors/IAOneTwoBodyInterferenceTensor.h"
#include "Physics/QuasiElastic/XSection/ELFormFactors.h"
#include "TVector3.h"
#include "TLorentzVector.h"
namespace genie {

class QELFormFactorsModelI;
class QELFormFactors;
class ELFormFactorsModelI;
class SpectralFunc;
class XSecIntegratorI;

class UnifiedQELPXSec : public XSecAlgorithmI {

public:

  UnifiedQELPXSec();
  UnifiedQELPXSec(std::string config);
  virtual ~UnifiedQELPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction* i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction* i) const;
  bool ValidProcess      (const Interaction* i) const;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry& config);
  void Configure (std::string param_set);

private:
  void LoadConfig (void);

  /// Spectator nucleons (and their weights) used to evaluate the one-body /
  /// two-body current interference for the current hit nucleon. The sample is
  /// redrawn only when the hit nucleon changes, so that the cross section stays
  /// a smooth function of the lepton angles for a fixed hit nucleon (as the
  /// adaptive integration in genie::NewQELXSec requires).
  const std::vector<IAOneTwoBodyInterferenceTensor::Spectator>&
    Spectators(const Target& target, const TLorentzVector& p4Ni) const;

  const QELFormFactorsModelI* fCCFormFactorsModel;
  const QELFormFactorsModelI* fNCFormFactorsModel;
  const QELFormFactorsModelI* fEMFormFactorsModel;
  mutable QELFormFactors fFormFactors;

  const NuclearModelI* fNuclModel;
  const PauliBlocker* fPauliBlocker;
  const XSecIntegratorI* fXSecIntegrator;
  double fCos8c2; ///< cos^2(cabibbo angle)
  double fXSecScale; ///< external xsec scaling factor
  bool fDoPauliBlocking;
  bool fDoqAlongZ;
  std::string fTensorModel;

  /// Evaluate the form factors at Q2tilde (true) or at the true Q2 (false)
  bool fFFAtQ2Tilde;

  // One-body / two-body current interference leading to single-nucleon
  // knock-out (arXiv:2312.12545). EM only for now.
  bool fDoIntf;
  const SpectralFunc* fTotSpectralFunc; ///< complete spectral function
  const SpectralFunc* fMFSpectralFunc;  ///< its mean-field part
  const ELFormFactorsModelI* fELFormFactorsModel;
  mutable ELFormFactors fELFormFactors;
  genie::twobody_currents_sf::ModelParams fIntfPar; ///< MeV
  double fCV3Norm, fCV4Norm, fCV5Norm; ///< N-Delta vector form factors at Q2 = 0
  double fMV2;                         ///< their dipole mass squared (GeV^2)
  int fNumSpectators; ///< spectators sampled per hit nucleon and isospin

  mutable TLorentzVector fCachedP4Ni;
  mutable int fCachedTgtPdg;
  mutable int fCachedHitNucPdg;
  mutable std::vector<IAOneTwoBodyInterferenceTensor::Spectator> fCachedSpectators;
};

} // genie namespace

#endif
