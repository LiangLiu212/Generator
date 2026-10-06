/** \file G4INCLCrossSectionsStrangeness.hh
 * \brief Multipion, mesonic Resonances and strange cross sections
 *
 * \date 1st March 2016
 * \author Jason Hirtz
 */

#ifndef G4INCLCROSSSECTIONSSTRANGENESS_HH
#define G4INCLCROSSSECTIONSSTRANGENESS_HH

#include "G4INCLCrossSectionsMultiPionsAndResonances.hh"
#include "G4INCLConfig.hh"
//#include <limits>

namespace G4INCL {
  /// \brief Multipion, mesonic Resonances and strange cross sections

  class CrossSectionsStrangeness : public CrossSectionsMultiPionsAndResonances {
    public:
      CrossSectionsStrangeness();
	  
      /// \brief second new total particle-particle cross section
      virtual double total(Particle const * const p1, Particle const * const p2);
      
      /// \brief second new elastic particle-particle cross section
      virtual double elastic(Particle const * const p1, Particle const * const p2);
      
      /// \brief correction to old cross section
      virtual double piNToxPiN(const int xpi, Particle const * const p1, Particle const * const p2);
      virtual double NNToxPiNN(const int xpi, Particle const * const p1, Particle const * const p2);
	  
      /// \brief elastic scattering for Nucleon-Strange Particles cross sections
      virtual double NYelastic(Particle const * const p1, Particle const * const p2);
      virtual double NKbelastic(Particle const * const p1, Particle const * const p2);
      virtual double NKelastic(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-Nucleon to Stange particles cross sections
      virtual double NNToNLK(Particle const * const p1, Particle const * const p2);
      virtual double NNToNSK(Particle const * const p1, Particle const * const p2);
      virtual double NNToNLKpi(Particle const * const p1, Particle const * const p2);
      virtual double NNToNSKpi(Particle const * const p1, Particle const * const p2);
      virtual double NNToNLK2pi(Particle const * const p1, Particle const * const p2);
      virtual double NNToNSK2pi(Particle const * const p1, Particle const * const p2);
      virtual double NNToNNKKb(Particle const * const p1, Particle const * const p2);
      
      virtual double NNToMissingStrangeness(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-Delta to Stange particles cross sections
      virtual double NDeltaToNLK(Particle const * const p1, Particle const * const p2);
      virtual double NDeltaToNSK(Particle const * const p1, Particle const * const p2);
      virtual double NDeltaToDeltaLK(Particle const * const p1, Particle const * const p2);
      virtual double NDeltaToDeltaSK(Particle const * const p1, Particle const * const p2);
      
      virtual double NDeltaToNNKKb(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-Pion to Strange particles cross sections
      virtual double NpiToLK(Particle const * const p1, Particle const * const p2);
		  double p_pimToLK0(Particle const * const p1, Particle const * const p2);
      virtual double NpiToSK(Particle const * const p1, Particle const * const p2);
		  double p_pipToSpKp(Particle const * const p1, Particle const * const p2);
		  virtual double p_pimToSzKz(Particle const * const p1, Particle const * const p2);
		  virtual double p_pimToSmKp(Particle const * const p1, Particle const * const p2);
		  virtual double p_pizToSzKp(Particle const * const p1, Particle const * const p2);
      virtual double NpiToLKpi(Particle const * const p1, Particle const * const p2);
      virtual double NpiToSKpi(Particle const * const p1, Particle const * const p2);
      virtual double NpiToLK2pi(Particle const * const p1, Particle const * const p2);
      virtual double NpiToSK2pi(Particle const * const p1, Particle const * const p2);
      virtual double NpiToNKKb(Particle const * const p1, Particle const * const p2);
      
      virtual double NpiToMissingStrangeness(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-Hyperon quasi-elastic cross sections
      virtual double NLToNS(Particle const * const p1, Particle const * const p2);
      virtual double NSToNL(Particle const * const p1, Particle const * const p2);
      virtual double NSToNS(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-Kaon cross sections
      virtual double NKToNK(Particle const * const p1, Particle const * const p2);
      virtual double NKToNKpi(Particle const * const p1, Particle const * const p2);
      virtual double NKToNK2pi(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-antiKaon cross sections
      virtual double NKbToNKb(Particle const * const p1, Particle const * const p2);
      virtual double NKbToSpi(Particle const * const p1, Particle const * const p2);
      virtual double NKbToLpi(Particle const * const p1, Particle const * const p2);
		  virtual double p_kmToL_pz(Particle const * const p1, Particle const * const p2);
      virtual double NKbToS2pi(Particle const * const p1, Particle const * const p2);
      virtual double NKbToL2pi(Particle const * const p1, Particle const * const p2);
		  virtual double p_kmToL_pp_pm(Particle const * const p1, Particle const * const p2);
      virtual double NKbToNKbpi(Particle const * const p1, Particle const * const p2);
      virtual double NKbToNKb2pi(Particle const * const p1, Particle const * const p2);
      
      /// \brief eta-Nucleon cross sections
      virtual double etaNToLK(Particle const * const p1, Particle const * const p2);
      virtual double etaNToSK(Particle const * const p1, Particle const * const p2);
      
      /// \brief Omega-Nucleon cross sections
      virtual double omegaNToLK(Particle const * const p1, Particle const * const p2);
      virtual double omegaNToSK(Particle const * const p1, Particle const * const p2);
      virtual double omegaNToPiPiN(Particle const * const p1, Particle const * const p2);
   
  protected:
      /// \brief Maximum number of outgoing pions in NN collisions
      static const int nMaxPiNN;
	  
      /// \brief Maximum number of outgoing pions in piN collisions
      static const int nMaxPiPiN;
	  
      /// \brief Horner coefficients for s11pz
      const HornerC7 s11pzHC;
      /// \brief Horner coefficients for s01pp
      const HornerC8 s01ppHC;
      /// \brief Horner coefficients for s01pz
      const HornerC4 s01pzHC;
      /// \brief Horner coefficients for s11pm
      const HornerC4 s11pmHC;
      /// \brief Horner coefficients for s12pm
      const HornerC5 s12pmHC;
      /// \brief Horner coefficients for s12pp
      const HornerC3 s12ppHC;
      /// \brief Horner coefficients for s12zz
      const HornerC4 s12zzHC;
      /// \brief Horner coefficients for s02pz
      const HornerC4 s02pzHC;
      /// \brief Horner coefficients for s02pm
      const HornerC6 s02pmHC;
      /// \brief Horner coefficients for s12mz
      const HornerC4 s12mzHC;
      

  };
}

#endif
