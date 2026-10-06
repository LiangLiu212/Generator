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


#ifndef G4INCLPbarAtrestEntry_hh
#define G4INCLPbarAtrestEntry_hh 1

namespace G4INCL {

  class FinalState;

  class PbarAtrestEntryChannel : public IChannel {
  public:
    PbarAtrestEntryChannel(Nucleus *n, Particle *p);
    virtual ~PbarAtrestEntryChannel();

    void fillFinalState(FinalState *fs);

    ParticleList makeMesonStar();
    double PbarCoulombicCascadeEnergy(int A, int Z);
    double n_annihilation(int A, int Z);
    IAvatarList bringMesonStar(ParticleList const &pL, Nucleus * const n);
    bool ProtonIsTheVictim();
    ThreeVector getAnnihilationPosition();
  
    double fctrl(const double arg1);
    double r1(const int n);
    double r2(const int n);
    double r3(double x, const int n);
    double r4(double x, const int n);
    double radial_wavefunction(double x, const int n);
    double densityP(double x);
    double densityN(double x);
    double overlapP(double &x, const int n);
    double overlapN(double &x, const int n);
    double read_file(std::string filename, std::vector<double>& probabilities, std::vector<std::vector<std::string>>& particle_types);
    int findStringNumber(double rdm, std::vector<double> yields);
    
  private:
    Nucleus *theNucleus;
    Particle *theParticle;


    INCL_DECLARE_ALLOCATION_POOL(PbarAtrestEntryChannel)
  };
}

#endif
