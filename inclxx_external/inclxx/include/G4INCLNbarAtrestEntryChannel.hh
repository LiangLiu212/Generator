/*
 * G4INCLNbarAtrestEntryChannel.cc
 *
 *  \date Aug 21, 2024
 * \author Olivier Lourgo
 */
#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLAllocationPool.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLICoulomb.hh"
#include <utility>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>


#ifndef G4INCLNbarAtrestEntry_hh
#define G4INCLNbarAtrestEntry_hh 1

namespace G4INCL{
	class FinalState;

	class NbarAtrestEntryChannel :public IChannel {
	public : 
	 NbarAtrestEntryChannel(Nucleus *n, Particle *p);
	 virtual ~NbarAtrestEntryChannel();

	 void fillFinalState(FinalState *fs);

	 ParticleList makeMesonStar();
	 IAvatarList bringMesonStar(ParticleList const &pL, Nucleus * const n);
	 bool ProtonIsTheVictim();
	 ThreeVector getAnnihilationPosition();

	 double Pabs(double x, double value);
	 double densityP();
	 double densityN();
	 double overlapP(double &x);
	 double overlapN(double &x);
	 double read_file(std::string filename, std::vector<double>& probabilities, std::vector<std::vector<std::string>>& particle_types);
     int findStringNumber(double rdm, std::vector<double> yields);

    private:
     Nucleus *theNucleus;
     Particle *theParticle;

     INCL_DECLARE_ALLOCATION_POOL(NbarAtrestEntryChannel)
	};
}
#endif