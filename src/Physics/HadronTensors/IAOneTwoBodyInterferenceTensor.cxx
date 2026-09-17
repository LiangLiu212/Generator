//____________________________________________________________________________
/*

   Impulse approximation one-body / two-body current interference tensor

   Liang Liu <liangliu \at fnal.gov>
   Fermi National Accelerator Laboratory
   Copyright (c) 2003-2026, The GENIE Collaboration
   For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________
#include <cmath>

#include "Framework/Conventions/Units.h"
#include "Physics/HadronTensors/IAOneTwoBodyInterferenceTensor.h"

using namespace genie::twobody_currents_sf;

//____________________________________________________________________________
genie::IAOneTwoBodyInterferenceTensor::IAOneTwoBodyInterferenceTensor(
    const ModelParams& par, const FormFactors& ff,
    const TLorentzVector& p4Ni, const TLorentzVector& p4Nf,
    const TLorentzVector& q4, int pdgNi, int pdgNf,
    const std::vector<Spectator>& spectators, bool has_axial, bool nc) :
  fPar(par), fFF(ff), fp4Ni(p4Ni), fp4Nf(p4Nf), fq4(q4),
  fPdgNi(pdgNi), fPdgNf(pdgNf), fSpectators(spectators),
  fHasAxial(has_axial), fNC(nc)
{

}
//____________________________________________________________________________
void genie::IAOneTwoBodyInterferenceTensor::initialize_tensor(
    std::complex<double> (&hadron_tensor)[4][4]) const
{
  // twobody_currents_sf works in MeV with Lorentz index 0..3 = (t,x,y,z)
  const double to_MeV = 1. / genie::units::MeV;

  // Initial nucleon: off-shell energy (put on shell inside, which also
  // defines the shifted energy transfer). Final nucleon: on shell with the
  // same mass as used by the one-body tensor.
  const double p1[4] = { fp4Ni.E()*to_MeV, fp4Ni.X()*to_MeV, fp4Ni.Y()*to_MeV, fp4Ni.Z()*to_MeV };
  const double q[4]  = { fq4.E()*to_MeV,   fq4.X()*to_MeV,   fq4.Y()*to_MeV,   fq4.Z()*to_MeV };
  double pp1[4] = { 0., fp4Nf.X()*to_MeV, fp4Nf.Y()*to_MeV, fp4Nf.Z()*to_MeV };
  pp1[0] = std::sqrt( fPar.xmn*fPar.xmn + pp1[1]*pp1[1] + pp1[2]*pp1[2] + pp1[3]*pp1[3] );

  std::complex<double> sum[4][4];
  for (int mu = 0; mu < 4; ++mu)
    for (int nu = 0; nu < 4; ++nu) sum[mu][nu] = 0.;

  SpinCurrent J1b, J2b_p, J2b_n;
  std::complex<double> A[4][4];
  for ( const Spectator& sp : fSpectators ) {
    double p2[4] = { 0., sp.p3.X()*to_MeV, sp.p3.Y()*to_MeV, sp.p3.Z()*to_MeV };
    p2[0] = std::sqrt( fPar.xmn*fPar.xmn + p2[1]*p2[1] + p2[2]*p2[2] + p2[3]*p2[3] );

    if ( sp.weight_p != 0. && sp.weight_n != 0. ) {
      ComputeCurrents(fPar, fFF, p1, pp1, p2, q, fPdgNi, fPdgNf,
        fHasAxial, fNC, J1b, J2b_p, J2b_n);
    }
    else if ( sp.weight_p != 0. ) {
      ComputeCurrents(fPar, fFF, p1, pp1, p2, q, fPdgNi, fPdgNf, 2212,
        fHasAxial, fNC, J1b, J2b_p);
    }
    else if ( sp.weight_n != 0. ) {
      ComputeCurrents(fPar, fFF, p1, pp1, p2, q, fPdgNi, fPdgNf, 2112,
        fHasAxial, fNC, J1b, J2b_n);
    }

    if ( sp.weight_p != 0. ) {
      InterferenceHadronTensor(J1b, J2b_p, A);
      for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu) sum[mu][nu] += sp.weight_p * A[mu][nu];
    }
    if ( sp.weight_n != 0. ) {
      InterferenceHadronTensor(J1b, J2b_n, A);
      for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu) sum[mu][nu] += sp.weight_n * A[mu][nu];
    }
  }

  // Bring the tensor to the conventions of IASingleNucleonTensor:
  //  - average over the initial nucleon spin (1/2)
  //  - spinors normalised to ubar u = m instead of 2m (1/4), GeV^2 instead
  //    of MeV^2
  //  - Lorentz index order (x,y,z,t)
  const double norm = 0.5 * 0.25 * genie::units::MeV * genie::units::MeV;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      hadron_tensor[i][j] = norm * sum[(i + 1) % 4][(j + 1) % 4];
    }
  }
}
