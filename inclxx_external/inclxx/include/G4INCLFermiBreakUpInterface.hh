#ifndef G4INCLFermiBreakUpInterface_hh
#define G4INCLFermiBreakUpInterface_hh 1

#include "G4INCLEventInfo.hh"

class G4INCLFermiBreakUp;

class G4INCLFermiBreakUpInterface {
public:
  G4INCLFermiBreakUpInterface();
  virtual ~G4INCLFermiBreakUpInterface();

  void deExciteRemnant(G4INCL::EventInfo *eventInfo, const int i);

private:
  G4INCLFermiBreakUp *theFermiBreakUpModel;
};

#endif
