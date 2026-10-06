
#pragma once

#include "G4Abla.hh"
#include "G4INCLConfig.hh"
#include "G4INCLEventInfo.hh"
#include "G4INCLIDeExcitation.hh"

class G4INCLAblaInterface : public G4INCL::IDeExcitation {
public:
  G4INCLAblaInterface(G4INCL::Config *);
  virtual ~G4INCLAblaInterface();

  virtual void deExciteRemnant(G4INCL::EventInfo *eventInfo, const int i);

private:
  G4INCL::Config *theConfig;
  G4AblaOutput *ablaResult;
  G4Abla *theABLAModel;
};
