#ifndef G4INCLSMMINTERFACE_HH_
#define G4INCLSMMINTERFACE_HH_ 1

#include "G4INCLConfig.hh"
#include "G4INCLIDeExcitation.hh"
#include "G4INCLEventInfo.hh"

extern "C" {
  void smm_botvina_(int*, int*, float*, float*, float*, float*, int*, float*, int*);
  float smmuserrng_();
  void rndm1_(int*);

  extern struct
  {
    float acv[300];
    float zpcv[300];
    float pcv[300];
    float xcv[300];
    float ycv[300];
    float zcv[300];
    int iv;
  } volant_;

  extern struct
  {
    int nsmmfr;
  } smmfrg_;

}

namespace SMMCXX {

  class SMMInterface : public G4INCL::IDeExcitation {
  public:
    SMMInterface(G4INCL::Config*);
    virtual ~SMMInterface();

    virtual void deExciteRemnant(G4INCL::EventInfo *eventInfo, const int i);
  };
}

#endif // G4INCLSMMINTERFACE_HH_
