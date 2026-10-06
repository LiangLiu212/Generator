#ifndef G4INCLNNToNNEtaChannel_hh
#define G4INCLNNToNNEtaChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
	class NNToNNEtaChannel : public IChannel {
	public:
		NNToNNEtaChannel(Particle *, Particle *);
		virtual ~NNToNNEtaChannel();
		
		void fillFinalState(FinalState *fs);
		
	private:
		int iso1; // like isosp, can be changed in isospinRepartition
		int iso2; // like isosp, can be changed in isospinRepartition
		Particle *particle1, *particle2;
		
		static const double angularSlope;
		
		
		INCL_DECLARE_ALLOCATION_POOL(NNToNNEtaChannel);
	};
}

#endif
