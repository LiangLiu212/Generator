
#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLAllocationPool.hh"
#include "G4INCLFinalState.hh"




#ifndef G4INCAntinucleiAtrestEntry_hh
#define G4INCLAntinucleiAtrestEntry_hh 1

namespace G4INCL{
	class FinalState;

	class AntinucleiAtrestEntryChannel: public IChannel{
	public :
	     AntinucleiAtrestEntryChannel(Nucleus *n, Cluster *ac, ThreeVector pos1, ThreeVector pos2);
	     AntinucleiAtrestEntryChannel(Nucleus *n, Particle *p);
	     virtual ~AntinucleiAtrestEntryChannel();
	     void fillFinalState(FinalState *fs);
         ThreeVector getAnnihilationPosition(ThreeVector nbarPos, ThreeVector pbarPos);
         ParticleList makeMesonStar();
         IAvatarList bringMesonStar(ParticleList const &pL, Nucleus * const n);

    private:
    	Nucleus *theNucleus;
    	Cluster *theantiComposite;
    	ThreeVector Posnbar; //Position of the annihilation from PbarAtrestEntryChannel
    	ThreeVector Pospbar; //Position of the annihilation from NbarAtrestEntryChannel
    	Particle *Meson; // For fillFinalState
    	int pbarListSize; //To know who is coming from pbar annihilation

    	INCL_DECLARE_ALLOCATION_POOL(AntinucleiAtrestEntryChannel)
	};
}

#endif