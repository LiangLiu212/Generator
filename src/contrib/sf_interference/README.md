# One-body / two-body current interference in `genie::UnifiedQELPXSec`

Interference between the one-body current and the two-body (pion-in-flight,
seagull, pion-pole, Delta) currents that leads to a **single-nucleon knock-out**
final state, in the spectral-function formalism of A. Lovato, N. Rocco and
N. Steinberg, arXiv:2312.12545. The second nucleon is a spectator that stays in
the Fermi sea. Since the final state is the quasielastic one, the term is added
to the hadron tensor of `UnifiedQELPXSec` (it is not a separate process, and no
negative-weight events are needed). **EM only for now.**

## Code

| file | role |
|------|------|
| `src/Physics/HadronTensors/twobody_currents_sf.{h,cxx}` | C++ port of N. Rocco's Fortran `dirac_matrices_intf` (ACHILLES `currents_intf.f90` @ e02d266): one-body current, two-body currents (direct - exchange, spectator spin summed) |
| `src/Physics/HadronTensors/IAOneTwoBodyInterferenceTensor.{h,cxx}` | the interference tensor in the conventions of `IASingleNucleonTensor` (so the two add) for a Monte Carlo sample of spectators |
| `src/Physics/QuasiElastic/XSection/UnifiedQELPXSec.{h,cxx}` | `IncludeOneTwoBodyInterference`, spectator sampling, `FormFactorsAtQ2Tilde` |
| `config/UnifiedQELPXSec.xml` | param set `Dipole-OneTwoBodyIntf`, documentation of the `TwoBodyIntf-*` parameters |
| `config/SpectralFunc.xml`, `data/evgen/nucl/spectral_functions/pke12_MF.data` | param set `MeanField`: mean-field part of the CBF 12C spectral function |

What is added to the one-body tensor, for a hit nucleon `(p1, E)` drawn from
the complete spectral function:

    A12 = S_MF(p1,E)/S_tot(p1,E) * rho/A * sum_tau2 N_MF(tau2) * < 1/2 sum_spins [ conj(J1b) J2b + conj(J2b) J1b ] >_p2

- only the mean-field part of the spectral function contributes (hole and
  spectator), hence the ratio `S_MF/S_tot` and `N_MF` = number of mean-field
  nucleons of isospin `tau2` (5.135 of 6 for 12C);
- `rho/A = kF^3/(1.5 pi^2)/A` is the inverse volume (`kF` = 225 MeV for 12C);
- `<>_p2` is the average over spectator momenta drawn from the mean-field
  momentum distribution, `J2b` includes `1/(2 E_2)`.

The spectators are drawn **once per hit nucleon** and cached: `genie::NewQELXSec`
integrates adaptively over the lepton angles for every nucleon throw, so the
cross section must be a smooth function of those angles for a fixed nucleon.
The Monte Carlo average over spectators is then part of the average over
nucleon throws (splines) or of the accept/reject loop (event generation).

If protons and neutrons share their mean-field table (12C), one spectator
momentum serves both spectator isospins: the two-body operators do not depend
on the spectator isospin. Cost: 32 us per hit nucleon for the two-body currents
(g++ -O2; one-body tensor 8 us), after the analytic contractions described in
the header of `twobody_currents_sf.cxx` (the literal port needed 236 us).

Only nuclei with a mean-field table in `SpectralFunc/MeanField` get the
interference; for all others the model is unchanged.

## Validation

1. `test_twobody_currents_sf.cxx` (standalone, only needs g++): the tensor
   reproduces N. Rocco's reference numbers of the ACHILLES unit test
   (`test/test_fortran_interference.f90`; EM, 12C, the four struck / spectator
   isospin combinations, 64 elements) to 1e-8 relative, and the one-body tensor
   of this module equals `onebody_currents_sf` to 1e-16 after the convention
   conversion applied by `IAOneTwoBodyInterferenceTensor`.
2. `gtest_qel_intf_window.cxx`: direct Monte Carlo integration of
   d2sigma/dOmega/domega in a scattering-angle window, several param sets
   evaluated on the same hit nucleons and lepton angles. With the
   validation-only sets of `config/` (put that directory first on `GXMLPATH`):

   e- 12C, 2.5 GeV, 14 < theta < 16 deg, 10^7 throws per hit nucleon, seed 20260917
   (GENIE_RC feature/qel-1b2b-interference, tune GEM26_22b_00_000, 2026-09-17):

   | param set | sigma(window) nb/atom | interference nb/atom |
   |-----------|----------------------|--------------|
   | ACHILLES v0.3.1 content (e02d266), QESpectral / Intf_Spectral_Func | 329.66 +/- 2.19 | 23.05 +/- 0.08 (+7.0%) |
   | `AchillesLike` (no Pauli blocking, FF at true Q2, Kelly) | 325.29 +/- 0.57 | |
   | `AchillesLike-OneTwoBodyIntf` (C4V/C5V as in ACHILLES) | 348.42 | 23.13 +/- 0.05 (+7.1%) |
   | `AchillesLike-OneTwoBodyIntf-C4C5unit` | 348.43 | 23.14 +/- 0.04 (+7.1%) |
   | `Dipole` (GENIE defaults) | 288.19 +/- 0.51 | |
   | `Dipole-OneTwoBodyIntf` | 309.37 | 21.17 +/- 0.04 (+7.3%) |

   (interference = paired difference on identical phase-space points, so its
   error is much smaller than that of either total.)
   Struck proton / neutron split of the interference: 12.48 / 10.65 nb
   (ACHILLES 12.39 / 10.66). The 15 MeV-bin omega spectra of `AchillesLike`
   and of its interference term agree with the ACHILLES event samples within
   the statistical errors of the latter over the whole quasielastic peak. With
   the GENIE defaults the interference is +11.7% of QE at omega = 0.10 GeV,
   +7.1% at the peak and +6.3% above it, and vanishes by omega ~ 0.5 GeV.
   No point in the window had a negative one-body + interference sum.

## Known differences with ACHILLES @ e02d266

- **C4V / C5V terms of the gamma-N-Delta vertex.** In `det_JaJb_JcJd` the
  Fortran adds them as *scalars* to a 4x4 array, which in Fortran adds the value
  to every element of the Dirac matrix; the axial `ca5` term next to them is
  multiplied by `id4`. The ACHILLES unit test has C4V = C5V = 0 and does not see
  this, the default `FormFactors.yml` has `cv4norm = -1.15`, `cv5norm = 0.48`.
  This port multiplies the unit matrix; `TwoBodyIntf-AchillesC4VC5V = true`
  reproduces ACHILLES (validation only). Numerically it hardly matters at the
  setting above: the paired difference is 0.015 +/- 0.032 nb.
- The pion-pole term (axial, unused for EM) uses the full `q^2` instead of
  `q0^2 - qz^2`, so it does not assume q along z.
- The in-medium Delta potential table (`rho_0p5.dat`) is not read: the Fortran
  looks the value up but never uses it.
- GENIE conventions that differ from ACHILLES `QESpectral` and are unrelated to
  the interference: Pauli blocking on, form factors evaluated at Q2tilde
  (`FormFactorsAtQ2Tilde`), BBA07 elastic form factors.

## To do

- Spline / event-generation closure test with `GEM26_22b_10_000` (gmkspl EMQE
  jobs launched 2026-09-17; GENIE_RC has no tunable EM Q2 cut, so gevgen around
  15 deg is inefficient there).

- CC and NC: the currents are ported (axial pieces, isospin raising structure)
  but not wired, not validated, and need `C5A`.
- Mean-field tables for other nuclei (ACHILLES ships 40Ar).
