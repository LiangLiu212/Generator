#ifndef G4INCLNNToNNOmegaChannel_hh
#define G4INCLNNToNNOmegaChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
	class NNToNNOmegaChannel : public IChannel {
	public:
		NNToNNOmegaChannel(Particle *, Particle *);
		virtual ~NNToNNOmegaChannel();
		
		void fillFinalState(FinalState *fs);
		
	private:
		int iso1; // like isosp, can be changed in isospinRepartition
		int iso2; // like isosp, can be changed in isospinRepartition
		Particle *particle1, *particle2;
		
		static const double angularSlope;
		
		
		INCL_DECLARE_ALLOCATION_POOL(NNToNNOmegaChannel);
	};
}

#endif
