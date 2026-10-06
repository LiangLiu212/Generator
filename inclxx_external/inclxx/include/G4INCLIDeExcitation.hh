#ifndef G4INCLIDeExcitation_hh
#define G4INCLIDeExcitation_hh 1

#include "G4INCLEventInfo.hh"
#include "G4INCLConfig.hh"

// FERMI_BREAKUP de-excitation
#ifdef INCL_DEEXCITATION_FERMI_BREAKUP
#include "G4INCLFermiBreakUpInterface.hh"
#endif

namespace G4INCL {
  class IDeExcitation {
  public:
    IDeExcitation(G4INCL::Config*);
    virtual ~IDeExcitation();

    virtual void deExcite(G4INCL::EventInfo*);
    virtual void deExciteRemnant(G4INCL::EventInfo*, const int) = 0;

  private:
    Config *theConfig;
#ifdef INCL_DEEXCITATION_FERMI_BREAKUP
    G4INCLFermiBreakUpInterface *theFermiBreakUp;
#endif
  };
}

#endif
