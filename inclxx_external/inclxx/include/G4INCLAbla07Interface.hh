#ifndef G4InclAbla07Interface_hh
#define G4InclAbla07Interface_hh 1

#include "G4INCLConfig.hh"
#include "G4INCLIDeExcitation.hh"
#include "G4INCLEventInfo.hh"
#include <string>

extern "C" {
  void ablantuple_(int*, int*);
  void rn07_(int*, int*, int*, double*, double*, double*, double*, double*);
  void abrainit_(char*, long int);
  double eflmac_(int*, int*, int*, int*);
  float abla07userrng_();

  extern struct
  {
    double ap;
    double zp;
    double at;
    double zt;
    double eap;
    double beta;
    double bmaxnuc;
    double crtot;
    double crnuc;
    double r_0;
    double r_p;
    double r_t;
    int imax;
    int irndm;
    double pi;
    double bfpro;
    double snpro;
    double shell;
  } abramain_;

  extern struct
  {
    double eap1;
    int ap1;
    int zp1;
  } proj_;

  const int nbiv = 300;
  extern struct
  {
    double ppx[nbiv];
    double ppy[nbiv];
    double ppz[nbiv];
    double eekin[nbiv];
    int ww[nbiv];
    int nna[nbiv];
    int nnz[nbiv];
    int nniv;
    int kprocout;
    int izf;
    int iaf;
    double esf;
  } zablaout_;

/*  extern struct
  {
    int isecdec;
    int imultbu;
    int kproc;
  } ablasmop_;*/

}

namespace ABLA07CXX {

  class Abla07Interface : public G4INCL::IDeExcitation {
  public:
    Abla07Interface(G4INCL::Config*);
    virtual ~Abla07Interface();

    virtual void deExciteRemnant(G4INCL::EventInfo *eventInfo, const int i);

  private:
    std::string historyIntToString(const int hist);
    class HistoryStepIntToChar {
      public:
        HistoryStepIntToChar();
        char operator()(const char hist);
        static const int historyDictionarySize = 10;
        static char historyDictionary[historyDictionarySize];
    };

    G4INCL::Config *theConfig;
    int verboseEvent;
    static const double fmp;
    static const double fmn;
  };
}

#endif
