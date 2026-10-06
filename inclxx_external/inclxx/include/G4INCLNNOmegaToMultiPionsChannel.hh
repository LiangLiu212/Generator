#ifndef G4INCLNNOmegaToMultiPionsChannel_hh
#define G4INCLNNOmegaToMultiPionsChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NNOmegaToMultiPionsChannel : public IChannel {
    public:
      NNOmegaToMultiPionsChannel(const int, Particle *, Particle *);
      virtual ~NNOmegaToMultiPionsChannel();

      void fillFinalState(FinalState *fs);

    private:
      int npion;
      int iso1; // like isosp, can be changed in isospinRepartition
      int iso2; // like isosp, can be changed in isospinRepartition
      int isosp[4];
      Particle *particle1, *particle2;

      static const double angularSlope;

      void inter2Part(const double p);
      void pn_ppPim();
      void pn_pnPi0();
      void pn_nnPip();
      void pp_nnPipPip();
      void nn_ppPimPim();
      void pn_pnPipPim();
      void pn_pnPi0Pi0();
      void pp_ppPipPim();
      void nn_nnPipPim();
      void pp_ppPi0Pi0();
      void nn_nnPi0Pi0();
      void pp_pnPipPi0();
      void pn_ppPimPi0();
      void pn_nnPipPi0();
      void nn_pnPimPi0();
      void pp_pnPipPi0Pi0();
      void nn_pnPimPi0Pi0();
      void pn_nnPipPi0Pi0();
      void pp_ppPipPimPi0();
      void nn_nnPipPimPi0();
      void pp_ppPi0Pi0Pi0();
      void nn_nnPi0Pi0Pi0();
      void pp_pnPipPipPim();
      void pp_nnPipPipPi0();
      void pn_ppPimPi0Pi0();
      void pn_ppPimPimPip();
      void pn_pnPi0PipPim();
      void pn_pnPi0Pi0Pi0();
      void pn_nnPipPipPim();
      void nn_pnPipPimPim();
      void nn_ppPimPimPi0();
      void pp_nnPipPipPi0Pi0();
      void pp_nnPipPipPipPim();
      void nn_ppPi0Pi0PimPim();
      void nn_ppPipPimPimPim();
      void pp_ppPi0Pi0Pi0Pi0();
      void nn_nnPi0Pi0Pi0Pi0();
      void pn_pnPi0Pi0Pi0Pi0();
      void pp_ppPipPi0Pi0Pim();
      void nn_nnPipPi0Pi0Pim();
      void pn_pnPipPi0Pi0Pim();
      void pp_ppPipPipPimPim();
      void nn_nnPipPipPimPim();
      void pn_pnPipPipPimPim();
      void pp_pnPipPi0Pi0Pi0();
      void pn_nnPipPi0Pi0Pi0();
      void pp_nnPipPi0Pi0Pi0();
      void pp_pnPipPipPi0Pim();
      void pn_nnPipPipPi0Pim();
      void pp_nnPipPipPi0Pim();
      void nn_pnPi0Pi0Pi0Pim();
      void pn_ppPi0Pi0Pi0Pim();
      void nn_pnPipPi0PimPim();
      void pn_ppPipPi0PimPim();
      void isospinRepartition();

      INCL_DECLARE_ALLOCATION_POOL(NNOmegaToMultiPionsChannel);
  };
}

#endif
