#ifndef G4INCLXXVInterfaceTally_hh
#define G4INCLXXVInterfaceTally_hh 1

#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4HadFinalState.hh"

class G4INCLXXVInterfaceTally {
  public:
    G4INCLXXVInterfaceTally() {}
    virtual ~G4INCLXXVInterfaceTally() {}

    virtual void Open() = 0;
    virtual void Close() = 0;
    virtual void Tally(G4HadProjectile const &aTrack, G4Nucleus const &theNucleus, G4HadFinalState const &result) = 0;
};

#endif

