#ifndef G4INCLNNbarToAnnihilationChannel_hh
#define G4INCLNNbarToAnnihilationChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"
#include "G4INCLNucleus.hh"

namespace G4INCL {
  class NNbarToAnnihilationChannel : public IChannel {
    public:
      NNbarToAnnihilationChannel(Nucleus *, Particle *, Particle *);
      virtual ~NNbarToAnnihilationChannel();

      double read_file(std::string filename, std::vector<double>& probabilities, std::vector<std::vector<std::string>>& particle_types);
      int findStringNumber(double rdm, std::vector<double> yields);
      void fillFinalState(FinalState *fs);
      
      Nucleus *theNucleus;

    private:
      //Nucleus *theNucleus;
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NNbarToAnnihilationChannel);
  };
}

#endif
