/** \file G4INCLCrossSectionsAntiparticles.hh
 * \brief Multipion, mesonic Resonances, strange cross sections and antinucleon as projectile
 *
 * \date 31st March 2023
 * \author Demid Zharenov
 */

#ifndef G4INCLCROSSSECTIONSANTIPARTICLES_HH
#define G4INCLCROSSSECTIONSANTIPARTICLES_HH

#include "G4INCLCrossSectionsStrangeness.hh"
#include "G4INCLConfig.hh"
//#include <limits>

namespace G4INCL {
  /// \brief Multipion, mesonic Resonances and strange cross sections

//  class CrossSectionsAntiparticles : public CrossSectionsMultiPionsAndResonances {
  class CrossSectionsAntiparticles : public CrossSectionsStrangeness {
    public:
      CrossSectionsAntiparticles();
	  
      /// \brief second new total particle-particle cross section
      virtual double total(Particle const * const p1, Particle const * const p2);
     
      /// \brief old elastic particle-particle cross section
      virtual double elastic(Particle const * const p1, Particle const * const p2);
    
      /// \brief Nucleon-AntiNucleon to Nucleon-AntiNucleon cross sections
      virtual double NNbarElastic(Particle const* const p1, Particle const* const p2);
      virtual double NNbarCEX(Particle const* const p1, Particle const* const p2);

      virtual double NNbarToLLbar(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-AntiNucleon to Nucleon-AntiNucleon + pions cross sections
      virtual double NNbarToNNbarpi(Particle const* const p1, Particle const* const p2);
      virtual double NNbarToNNbar2pi(Particle const* const p1, Particle const* const p2);
      virtual double NNbarToNNbar3pi(Particle const* const p1, Particle const* const p2);
     
      /// \brief Nucleon-AntiNucleon total annihilation cross sections
      virtual double NNbarToAnnihilation(Particle const* const p1, Particle const* const p2);
         
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
