/*
 * =====================================================================================
 *
 *       Filename:  twobody_currents_sf.h
 *
 *    Description:  One-body / two-body current interference leading to a
 *                  single-nucleon (1p1h) final state, spectral-function formalism.
 *
 *                  Native C++ port of N. Rocco's Fortran module dirac_matrices_intf
 *                  (ACHILLES, src/Achilles/fortran/currents_intf.f90 @ e02d266).
 *                  Reference: Lovato, Rocco, Steinberg, arXiv:2312.12545
 *
 *                  Everything in this module works in MeV, with Lorentz index
 *                  0..3 = (t, x, y, z), exactly as the Fortran does.
 *
 *         Author:  Liang Liu (L. Liu), liangliu@fnal.gov
 *		    Fermi National Accelerator Laboratory
 *  Collaboration:  GENIE
 *
 *  \cpright  Copyright (c) 2003-2026, The GENIE Collaboration
 *            For the full text of the license visit http://copyright.genie-mc.org
 *            or see $GENIE/LICENSE
 *
 * =====================================================================================
 */
#ifndef __TWOBODY_CURRENTS_SF_H__
#define __TWOBODY_CURRENTS_SF_H__

#include <complex>
#include <array>
#include <cmath>

namespace genie {
namespace twobody_currents_sf {

  // Masses and couplings of the meson-exchange / Delta currents (MeV)
  struct ModelParams {
    double xmd;      // Delta mass
    double xmn;      // nucleon mass
    double xmpi;     // pion mass
    double xmrho;    // rho mass
    double fpind;    // pi-N-Delta coupling entering the Delta width
    double fstar;    // pi-N-Delta coupling f*
    double fpinn2;   // f_piNN^2
    double ga;       // axial coupling (seagull / pion-pole axial pieces)
    double lpi;      // pi-N-N   cutoff
    double lpind;    // pi-N-Delta cutoff
    // ACHILLES @ e02d266 adds the C4V and C5V Delta-vertex terms to every
    // element of the Dirac matrix (Fortran scalar + array broadcast) instead
    // of multiplying the unit matrix. true reproduces that, for validation
    // against ACHILLES only.
    bool   c4c5_broadcast;
  };

  // Form factors, couplings stripped (GENIE keeps them in the xsec prefactor)
  struct FormFactors {
    double f1;     // struck-nucleon Dirac FF
    double f2;     // struck-nucleon Pauli FF (kappa included)
    double fa;     // axial
    double fap;    // pseudoscalar
    double fpiem;  // pion EM form factor (Gep - Gen)
    double cv3;    // N-Delta vector C3V
    double cv4;    // N-Delta vector C4V
    double cv5;    // N-Delta vector C5V
    double ca5;    // N-Delta axial  C5A
  };

  // J[f1][i1][mu]: final (f1) / initial (i1) struck-nucleon spin projections
  typedef std::array<std::array<std::array<std::complex<double>, 4>, 2>, 2> SpinCurrent;

  // One-body current J1b and two-body current J2b (spectator summed over spin,
  // direct minus exchange, divided by 2*E_spectator) for
  //   l + N(p1) -> l' + N(pp1), spectator nucleon N(p2) -> N(p2).
  // p1[0] is the OFF-shell struck-nucleon energy and q the true (omega, qvec);
  // the de Forest shift q0 -> q0 + p1[0] - E_onshell(p1) is applied inside.
  // Spinors are normalised to ubar u = 2 m.
  void ComputeCurrents(const ModelParams & par, const FormFactors & ff,
                       const double p1[4], const double pp1[4],
                       const double p2[4], const double q[4],
                       int pdg_in, int pdg_out, int pdg_spect,
                       bool has_axial, bool nc,
                       SpinCurrent & J1b, SpinCurrent & J2b);

  // Same, for a proton and a neutron spectator of the same momentum p2 (the
  // two-body operators do not depend on the spectator isospin, so this costs
  // little more than a single spectator)
  void ComputeCurrents(const ModelParams & par, const FormFactors & ff,
                       const double p1[4], const double pp1[4],
                       const double p2[4], const double q[4],
                       int pdg_in, int pdg_out, bool has_axial, bool nc,
                       SpinCurrent & J1b,
                       SpinCurrent & J2b_pspect, SpinCurrent & J2b_nspect);

  // R[mu][nu] = sum_spins ( J2b^mu conj(J1b^nu) + conj(J2b^mu) J1b^nu )
  // (N. Rocco's convention, used by the validation test)
  void InterferenceTensor(const SpinCurrent & J1b, const SpinCurrent & J2b,
                          std::complex<double> R[4][4]);

  // A[mu][nu] = sum_spins ( conj(J1b^mu) J2b^nu + conj(J2b^mu) J1b^nu )
  // Hermitian, same index/conjugation convention as the one-body tensor
  // sum_spins conj(J^mu) J^nu of onebody_currents_sf
  void InterferenceHadronTensor(const SpinCurrent & J1b, const SpinCurrent & J2b,
                                std::complex<double> A[4][4]);

  // sum_spins conj(J^mu) J^nu of a single current
  void SquaredHadronTensor(const SpinCurrent & J, std::complex<double> A[4][4]);

}  // namespace twobody_currents_sf
}  // namespace genie

#endif  // __TWOBODY_CURRENTS_SF_H__
