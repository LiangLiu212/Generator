#ifndef G4INCLGEMINIXXINTERFACE_HH_
#define G4INCLGEMINIXXINTERFACE_HH_

#include "G4INCLIDeExcitation.hh"
#include "G4INCLEventInfo.hh"
#include "G4INCLConfig.hh"

double geminixxuserrng();

class G4INCLGEMINIXXInterface : public G4INCL::IDeExcitation {
  public:
    G4INCLGEMINIXXInterface(G4INCL::Config *theConfig);
    virtual ~G4INCLGEMINIXXInterface();

    virtual void deExciteRemnant(G4INCL::EventInfo *eventInfo, const int i);

  private:
    static const double theSpeedOfLight;
    char *theAssignment;
};

#endif // G4INCLGEMINIXXINTERFACE_HH_
