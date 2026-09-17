//____________________________________________________________________________
/*

   \class    genie::IAOneTwoBodyInterferenceTensor

   \brief    Concrete implementation of HadronTensorI interface

   Impulse approximation response tensor for the interference between the
   one-body current and the two-body (meson-exchange + Delta) currents that
   leads to a single-nucleon knock-out final state, in the spectral-function
   formalism. The second nucleon is a spectator that stays in the Fermi sea;
   the caller supplies a Monte Carlo sample of spectators.

   The tensor uses the normalisation and index conventions of
   genie::IASingleNucleonTensor, so that the two can simply be added.

   Ref: A. Lovato, N. Rocco, N. Steinberg, arXiv:2312.12545

   \created  Sep 17, 2026
   \author   Liang Liu <liangliu \at fnal.gov>
   Fermi National Accelerator Laboratory

   \cpright  Copyright (c) 2003-2026, The GENIE Collaboration
   For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _IA_ONE_TWO_BODY_INTERFERENCE_TENSOR_H_
#define _IA_ONE_TWO_BODY_INTERFERENCE_TENSOR_H_

#include <vector>

#include <TLorentzVector.h>
#include <TVector3.h>

// GENIE includes
#include "Physics/HadronTensors/NucleonTensor.h"
#include "Physics/HadronTensors/twobody_currents_sf.h"

namespace genie {

  class IAOneTwoBodyInterferenceTensor : public NucleonTensor {

    public:

      // A spectator momentum and the (dimensionless) weights that the
      // contributions of a proton and of a neutron spectator of that momentum
      // to the tensor carry. A vanishing weight skips that isospin.
      struct Spectator {
        TVector3 p3;        // GeV
        double   weight_p;
        double   weight_n;
      };

      // p4Ni is the off-shell initial nucleon and q4 the true 4-momentum
      // transfer (GeV). par is in MeV, as twobody_currents_sf expects.
      IAOneTwoBodyInterferenceTensor(
          const genie::twobody_currents_sf::ModelParams& par,
          const genie::twobody_currents_sf::FormFactors& ff,
          const TLorentzVector& p4Ni, const TLorentzVector& p4Nf,
          const TLorentzVector& q4, int pdgNi, int pdgNf,
          const std::vector<Spectator>& spectators,
          bool has_axial, bool nc);

      // Overridden initialization function
      virtual void initialize_tensor(std::complex<double> (&hadron_tensor)[4][4]) const override;

      ~IAOneTwoBodyInterferenceTensor() {}

    protected:
      genie::twobody_currents_sf::ModelParams fPar;
      genie::twobody_currents_sf::FormFactors fFF;
      TLorentzVector fp4Ni;
      TLorentzVector fp4Nf;
      TLorentzVector fq4;
      int fPdgNi;
      int fPdgNf;
      std::vector<Spectator> fSpectators;
      bool fHasAxial;
      bool fNC;
  };

}
#endif
