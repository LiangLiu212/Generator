//____________________________________________________________________________
/*
  Inclusive (e,e') double-differential cross section of genie::UnifiedQELPXSec
  in a scattering-angle window, by direct Monte Carlo integration of the
  kPSQELEvGen differential cross section (hit nucleon from the nuclear model,
  lepton angles flat in the probe + nucleon COM frame). No splines, no event
  generation, no EM Q2 cut: meant for comparing param sets of the model with
  each other (e.g. with / without the one-body / two-body current interference)
  and with external calculations.

  All the listed UnifiedQELPXSec param sets are evaluated on the SAME hit
  nucleons and lepton angles, so their ratios have small statistical errors.
  The first one must be cheap (no interference): it alone is evaluated
  outside the angular window.

  Build (after sourcing the GENIE environment; the *_PKG_DIR variables are
  those of a spack-based setup, adapt the -L/-I paths to yours):
    g++ -std=c++17 -O2 gtest_qel_intf_window.cxx -o gtest_qel_intf_window \
        -I$GENIE/src $(root-config --cflags) $(genie-config --libs) \
        $(root-config --glibs) -lGeom -lMathMore \
        -I$LOG4CPP_PKG_DIR/include -L$LOG4CPP_PKG_DIR/lib -llog4cpp \
        -I$LIBXML2_PKG_DIR/include/libxml2 -L$LIBXML2_PKG_DIR/lib -lxml2 \
        -L$PYTHIA8_PKG_DIR/lib -lpythia8 -L$LHAPDF_PKG_DIR/lib -lLHAPDF \
        -L$GSL_PKG_DIR/lib -lgsl -lgslcblas

  Usage:
    gtest_qel_intf_window --tune GEM26_22b_00_000 -n 20000000 -e 2.5 \
        --theta-min 14 --theta-max 16 -t 1000060120 \
        -c Dipole,Dipole-OneTwoBodyIntf [-o out.txt] [--seed 1]
*/
//____________________________________________________________________________
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <TLorentzVector.h>
#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/StringUtils.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"

using namespace genie;

int main(int argc, char ** argv)
{
  CmdLnArgParser parser(argc, argv);
  long   n_throws  = parser.OptionExists('n') ? parser.ArgAsLong('n') : 1000000;
  double e_beam    = parser.OptionExists('e') ? parser.ArgAsDouble('e') : 2.5;
  int    tgt_pdg   = parser.OptionExists('t') ? parser.ArgAsInt('t') : 1000060120;
  double theta_min = parser.OptionExists("theta-min") ? parser.ArgAsDouble("theta-min") : 14.;
  double theta_max = parser.OptionExists("theta-max") ? parser.ArgAsDouble("theta-max") : 16.;
  long   seed      = parser.OptionExists("seed") ? parser.ArgAsLong("seed") : 20260917;
  std::string out  = parser.OptionExists('o') ? parser.ArgAsString('o') : "";
  std::string cfgs = parser.OptionExists('c') ? parser.ArgAsString('c') : "Dipole,Dipole-OneTwoBodyIntf";
  std::vector<std::string> configs = utils::str::Split(cfgs, ",");

  RunOpt::Instance()->ReadFromCommandLine(argc, argv);
  if ( !RunOpt::Instance()->Tune() ) {
    LOG("gtestQELIntf", pFATAL) << "No TuneId in RunOption";
    exit(-1);
  }
  RunOpt::Instance()->BuildTune();
  utils::app_init::MesgThresholds( RunOpt::Instance()->MesgThresholdFiles() );
  utils::app_init::RandGen(seed);

  AlgFactory * algf = AlgFactory::Instance();
  std::vector<const XSecAlgorithmI *> models;
  for ( const auto & c : configs ) {
    const XSecAlgorithmI * m = dynamic_cast<const XSecAlgorithmI *>(
      algf->GetAlgorithm("genie::UnifiedQELPXSec", c) );
    if ( !m ) { std::fprintf(stderr, "no UnifiedQELPXSec/%s\n", c.c_str()); return 1; }
    models.push_back(m);
  }
  const NuclearModelI * nucl_model = dynamic_cast<const NuclearModelI *>(
    algf->GetAlgorithm("genie::SpectralFunc", "Default") );

  const double omega_bin = 0.015;  // GeV
  const int    n_bins    = 60;
  const size_t n_models  = models.size();
  // [model][hit nucleon][bin]: sum of weights and of squared weights
  std::vector<std::vector<std::vector<double> > > sum(n_models,
    std::vector<std::vector<double> >(2, std::vector<double>(n_bins + 1, 0.)));
  std::vector<std::vector<std::vector<double> > > sum2 = sum;
  std::vector<long> n_negative(n_models, 0);
  // sum of w_m * w_k over all window points: errors of differences of totals
  std::vector<std::vector<double> > cov(n_models, std::vector<double>(n_models, 0.));
  std::vector<double> w_point(n_models, 0.);

  RandomGen * rnd = RandomGen::Instance();
  const int hit_pdgs[2] = { kPdgProton, kPdgNeutron };

  for ( int ih = 0; ih < 2; ++ih ) {
    Interaction * inter = Interaction::QELEM(tgt_pdg, hit_pdgs[ih], kPdgElectron, e_beam);
    Target * tgt = inter->InitStatePtr()->TgtPtr();

    for ( long i = 0; i < n_throws; ++i ) {
      nucl_model->GenerateNucleon(*tgt, 0.);
      tgt->SetHitNucPosition(0.);

      double cos_theta_0 = -1. + 2. * rnd->RndGen().Rndm();
      double phi_0 = 2. * constants::kPi * rnd->RndGen().Rndm();
      double Eb = 0.;

      // flat measure in (cos_theta_0, phi_0) -> weight 4 pi / N
      const double w_ps = 4. * constants::kPi / n_throws / units::nb;

      double xsec0 = utils::ComputeFullQELPXSec(inter, nucl_model, models[0],
        cos_theta_0, phi_0, Eb, kUseNuclearModel, 0., true);
      if ( xsec0 <= 0. ) continue;

      const TLorentzVector & lep = inter->Kine().FSLeptonP4();
      double theta = lep.Theta() * TMath::RadToDeg();
      if ( theta < theta_min || theta > theta_max ) continue;
      int bin = static_cast<int>( (e_beam - lep.E()) / omega_bin );
      if ( bin < 0 || bin > n_bins ) bin = n_bins;  // overflow

      for ( size_t m = 0; m < n_models; ++m ) {
        double xsec = (m == 0) ? xsec0 : utils::ComputeFullQELPXSec(inter, nucl_model,
          models[m], cos_theta_0, phi_0, Eb, kUseNuclearModel, 0., true);
        if ( xsec < xsec0 * 1e-12 && m > 0 ) ++n_negative[m];
        sum[m][ih][bin]  += xsec * w_ps;
        sum2[m][ih][bin] += xsec * w_ps * xsec * w_ps;
        w_point[m] = xsec * w_ps;
      }
      for ( size_t m = 0; m < n_models; ++m )
        for ( size_t k = 0; k < n_models; ++k ) cov[m][k] += w_point[m] * w_point[k];
    }
    delete inter;
  }

  FILE * f = out.empty() ? stdout : std::fopen(out.c_str(), "w");
  const double d_omega_sr = 2. * constants::kPi
    * ( std::cos(theta_min * TMath::DegToRad()) - std::cos(theta_max * TMath::DegToRad()) );

  std::fprintf(f, "# e- %d, E = %g GeV, theta in [%g, %g] deg (dOmega = %.6f sr), N = %ld throws per hit nucleon, seed %ld\n",
    tgt_pdg, e_beam, theta_min, theta_max, d_omega_sr, n_throws, seed);
  for ( size_t m = 0; m < n_models; ++m ) {
    double tot[2] = {0., 0.}, err2 = 0.;
    for ( int ih = 0; ih < 2; ++ih )
      for ( int b = 0; b <= n_bins; ++b ) { tot[ih] += sum[m][ih][b]; err2 += sum2[m][ih][b]; }
    std::fprintf(f, "# model %zu = UnifiedQELPXSec/%s : sigma(window) = %.4f +/- %.4f nb/atom (p %.4f + n %.4f), zero/negative in window: %ld\n",
      m, configs[m].c_str(), tot[0] + tot[1], std::sqrt(err2), tot[0], tot[1], n_negative[m]);
  }
  // differences of totals, with the error of the paired difference
  for ( size_t m = 1; m < n_models; ++m ) {
    for ( size_t k = 0; k < m; ++k ) {
      double tm = 0., tk = 0.;
      for ( int ih = 0; ih < 2; ++ih )
        for ( int b = 0; b <= n_bins; ++b ) { tm += sum[m][ih][b]; tk += sum[k][ih][b]; }
      double var = cov[m][m] + cov[k][k] - 2. * cov[m][k];
      std::fprintf(f, "# sigma[%zu] - sigma[%zu] = %.4f +/- %.4f nb/atom (%+.2f%%)\n",
        m, k, tm - tk, std::sqrt(std::max(0., var)), 100. * (tm - tk) / tk);
    }
  }
  std::fprintf(f, "# omega_lo omega_hi");
  for ( size_t m = 0; m < n_models; ++m ) std::fprintf(f, "  d2s[%zu] err[%zu]", m, m);
  std::fprintf(f, "   (nb/sr/GeV)\n");
  for ( int b = 0; b < n_bins; ++b ) {
    std::fprintf(f, "%.3f %.3f", b * omega_bin, (b + 1) * omega_bin);
    for ( size_t m = 0; m < n_models; ++m ) {
      double s = sum[m][0][b] + sum[m][1][b];
      double e = std::sqrt(sum2[m][0][b] + sum2[m][1][b]);
      std::fprintf(f, "  %.2f %.2f", s / d_omega_sr / omega_bin, e / d_omega_sr / omega_bin);
    }
    std::fprintf(f, "\n");
  }
  if ( !out.empty() ) std::fclose(f);
  return 0;
}
