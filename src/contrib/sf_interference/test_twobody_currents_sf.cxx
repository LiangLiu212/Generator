//____________________________________________________________________________
/*
  Standalone validation of genie::twobody_currents_sf against N. Rocco's
  reference one-body/two-body interference tensor. Kinematics, form factors and
  reference numbers are those of the ACHILLES unit test
  test/test_fortran_interference.f90 (@ e02d266): EM, 12C (kF = 225 MeV), all
  four struck-nucleon / spectator isospin combinations.

  No GENIE or ROOT libraries are needed:

    g++ -std=c++17 -O2 -I$GENIE/src \
        $GENIE/src/contrib/sf_interference/test_twobody_currents_sf.cxx \
        $GENIE/src/Physics/HadronTensors/twobody_currents_sf.cxx \
        $GENIE/src/Physics/HadronTensors/onebody_currents_sf.cxx \
        -o test_twobody_currents_sf && ./test_twobody_currents_sf

  Exit code 0 = all 64 tensor elements agree to the requested tolerance, and
  the one-body tensor built from this module's J1b reproduces
  genie::onebody_currents_sf once it is brought to that module's conventions
  (the same conversion IAOneTwoBodyInterferenceTensor applies).
*/
//____________________________________________________________________________
#include <cmath>
#include <complex>
#include <cstdio>

#include "Physics/HadronTensors/twobody_currents_sf.h"
#include "Physics/HadronTensors/onebody_currents_sf.h"

using namespace genie::twobody_currents_sf;

namespace {

  const double kRefNP[4][4] = {
    { 2.7740890726777619e-4, -3.7965422690241157e-4, -2.3231176568845108e-4,  4.1061096428384364e-4},
    {-4.3161000558052597e-4,  2.5764714497521469e-3,  2.0220656479648386e-4, -2.4244417210859243e-3},
    {-7.1436783649723686e-4,  1.6268240965361402e-4,  4.3096945527350689e-3,  1.0765531068087364e-4},
    { 4.0234418013681471e-4, -2.4126274543634810e-3,  2.2363340682766516e-4,  2.3185632439470393e-3}};

  const double kRefNN[4][4] = {
    {-3.3020181556906658e-6, -1.3001929821670444e-4, -5.9958819055489485e-4,  4.3404339239395113e-5},
    {-6.2171820196583785e-5,  4.0821888540744415e-4, -1.3414578647709242e-4, -4.0158278257070409e-4},
    { 1.4288406444307087e-5, -6.8658892892674422e-5, -2.7998584954695824e-4,  3.3064525168051652e-5},
    { 5.5817444445407680e-5, -4.1834036577303654e-4, -1.3053010694973955e-4,  3.7691312475286227e-4}};

  const double kRefPN[4][4] = {
    {-1.7885469087056062e-3, -7.8227095342742073e-4, -2.3389389889856866e-4,  2.5052082954706289e-5},
    {-2.0521391860787834e-4,  3.8534930856720922e-3,  2.7741445279762321e-4, -3.4703733478158249e-3},
    { 5.0733821810885205e-3,  8.5627806800367043e-4,  6.0672535084139809e-3,  1.7820579186363615e-3},
    { 1.2206241730577179e-4, -3.6152094763973706e-3,  3.5427731528894448e-4,  3.3032579363516269e-3}};

  const double kRefPP[4][4] = {
    { 3.5130429995087597e-4, -1.5649428658579787e-4, -9.0439040078329200e-4,  1.5847629614406479e-4},
    {-5.2010974451786391e-5,  6.0767008856932065e-4, -2.0035988588265645e-4, -5.8325143403928457e-4},
    { 3.8351027153701806e-5, -9.9801838469260902e-5, -4.1488281068150548e-4,  5.3449626473519358e-5},
    { 1.7778803543962618e-4, -6.0906130238444099e-4, -1.9776522199896026e-4,  5.8255693288547657e-4}};

  // returns the largest |computed - reference| / max|reference|
  double Check(const char * label, int pdg_in, int pdg_spect, const double ref[4][4])
  {
    const double kf = 225.0, A = 12.0, xmn = 938.91875;
    const double V = std::pow(kf, 3) / (1.5 * M_PI * M_PI) / A;

    ModelParams par;
    par.xmd   = 1232.25;
    par.xmn   = (938.27208816 + 939.56542054) / 2.0;   // ACHILLES Constant::mN
    par.xmpi  = 139.57018;
    par.xmrho = 775.8;
    par.fpind = 0.54;
    par.fstar = 2.13;
    par.fpinn2 = 0.081 * 4.0 * M_PI;
    par.ga    = 1.26;
    par.lpi   = 1300.0;
    par.lpind = 1150.0;
    par.c4c5_broadcast = false;

    FormFactors ff;
    if(pdg_in == 2112) { ff.f1 = -.017340697711910772;  ff.f2 = -0.99011610632251468; }
    else               { ff.f1 = 0.59137717615671226;   ff.f2 = 0.89838354961982203;  }
    ff.fa = 0.0;  ff.fap = 0.0;
    ff.fpiem = 0.47263563216390736;
    ff.cv3 = 1.2997163009564785;
    ff.cv4 = 0.0;  ff.cv5 = 0.0;  ff.ca5 = 0.0;

    const double p1[4]  = { 892.28663457549249, -77.036892302107034, -26.857649369116537, 67.942565322702237};
    const double p2[4]  = {1002.007566138, 96.6849111519814, -193.77595264909456, 274.8702302143906};
    const double pp1[4] = {1090.4949377893222, 281.3919769883637, -74.60793248850116, 470.7861894114025};
    const double q[4]   = { 198.20830321382988, 358.42886929047074, -47.750283119384619, 402.84362408870038};

    const double E_OS = std::sqrt(xmn*xmn + p1[1]*p1[1] + p1[2]*p1[2] + p1[3]*p1[3]);

    SpinCurrent J1b, J2b;
    ComputeCurrents(par, ff, p1, pp1, p2, q, pdg_in, pdg_in, pdg_spect, false, false, J1b, J2b);

    std::complex<double> R[4][4];
    InterferenceTensor(J1b, J2b, R);

    double ref_max = 0.0, dev_max = 0.0;
    for(int mu = 0; mu < 4; ++mu)
      for(int nu = 0; nu < 4; ++nu) ref_max = std::max(ref_max, std::abs(ref[mu][nu]));

    std::printf("\n%s  (struck %d, spectator %d)\n", label, pdg_in, pdg_spect);
    for(int mu = 0; mu < 4; ++mu)
    {
      for(int nu = 0; nu < 4; ++nu)
      {
        // flux, initial-spin average and final-nucleon phase space as in the Fortran test
        const double val = V * R[mu][nu].real() / 2.0 / (2.0 * E_OS) / (2.0 * pp1[0]);
        const double dev = std::abs(val - ref[mu][nu]);
        dev_max = std::max(dev_max, dev);
        std::printf("  R(%d,%d) = % .16e   ref % .16e   diff % .2e   imag % .1e\n",
            mu, nu, val, ref[mu][nu], val - ref[mu][nu], R[mu][nu].imag());
      }
    }
    std::printf("  max |diff| / max |ref| = %.3e\n", dev_max / ref_max);
    return dev_max / ref_max;
  }

  // One-body tensor: this module (MeV, ubar u = 2m, (t,x,y,z)) vs
  // onebody_currents_sf (GeV, ubar u = m, spin averaged, (x,y,z,t))
  double CheckOneBodyConventions()
  {
    ModelParams par;
    par.xmd = 1232.25;  par.xmn = 938.91875;  par.xmpi = 139.57018;  par.xmrho = 775.8;
    par.fpind = 0.54;  par.fstar = 2.13;  par.fpinn2 = 0.081 * 4.0 * M_PI;
    par.ga = 1.26;  par.lpi = 1300.0;  par.lpind = 1150.0;  par.c4c5_broadcast = false;

    FormFactors ff;
    ff.f1 = 0.59137717615671226;  ff.f2 = 0.89838354961982203;
    ff.fa = -0.7;  ff.fap = 0.3;   // arbitrary, to exercise the axial pieces too
    ff.fpiem = 0.47263563216390736;  ff.cv3 = 1.2997163009564785;
    ff.cv4 = 0.0;  ff.cv5 = 0.0;  ff.ca5 = 0.0;

    const double p1[4] = { 892.28663457549249, -77.036892302107034, -26.857649369116537, 67.942565322702237};
    const double p2[4] = {1002.007566138, 96.6849111519814, -193.77595264909456, 274.8702302143906};
    const double q[4]  = { 198.20830321382988, 358.42886929047074, -47.750283119384619, 402.84362408870038};
    double pp1[4] = {0., p1[1]+q[1], p1[2]+q[2], p1[3]+q[3]};
    pp1[0] = std::sqrt(par.xmn*par.xmn + pp1[1]*pp1[1] + pp1[2]*pp1[2] + pp1[3]*pp1[3]);

    SpinCurrent J1b, J2b;
    ComputeCurrents(par, ff, p1, pp1, p2, q, 2212, 2212, 2112, true, true, J1b, J2b);
    std::complex<double> A[4][4];
    SquaredHadronTensor(J1b, A);

    const double E_on = std::sqrt(par.xmn*par.xmn + p1[1]*p1[1] + p1[2]*p1[2] + p1[3]*p1[3]);
    const double wt   = q[0] + p1[0] - E_on;
    std::complex<double> H[4][4];
    genie::onebody_currents_sf::ComputeNucleonTensor(par.xmn*1e-3, q[0]*1e-3, wt*1e-3,
        p1[1]*1e-3, p1[2]*1e-3, p1[3]*1e-3, q[1]*1e-3, q[2]*1e-3, q[3]*1e-3,
        ff.f1, ff.f2, ff.fa, ff.fap, H);

    double ref_max = 0.0, dev_max = 0.0;
    for(int i = 0; i < 4; ++i)
    {
      for(int j = 0; j < 4; ++j)
      {
        const std::complex<double> mine = 0.5 * 0.25 * 1e-6 * A[(i+1)%4][(j+1)%4];
        ref_max = std::max(ref_max, std::abs(H[i][j]));
        dev_max = std::max(dev_max, std::abs(mine - H[i][j]));
      }
    }
    std::printf("\none-body tensor vs onebody_currents_sf: max |diff| / max |ref| = %.3e\n",
        dev_max / ref_max);
    return dev_max / ref_max;
  }

}  // anonymous namespace

int main()
{
  double worst = 0.0;
  worst = std::max(worst, Check("np", 2112, 2212, kRefNP));
  worst = std::max(worst, Check("nn", 2112, 2112, kRefNN));
  worst = std::max(worst, Check("pn", 2212, 2112, kRefPN));
  worst = std::max(worst, Check("pp", 2212, 2212, kRefPP));

  const double tol = 1e-6;
  std::printf("\nworst relative deviation = %.3e (tolerance %.0e): %s\n",
      worst, tol, worst < tol ? "PASS" : "FAIL");

  const double conv = CheckOneBodyConventions();
  std::printf("one-body convention check (tolerance 1e-10): %s\n", conv < 1e-10 ? "PASS" : "FAIL");

  return (worst < tol && conv < 1e-10) ? 0 : 1;
}
