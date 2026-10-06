/** \file G4INCLCrossSectionsMultiPionsAndResonances.hh
 * \brief Multipion and mesonic Resonances cross sections
 *
 * \date 4th February 2014
 * \author Jean-Christophe David
 */

#ifndef G4INCLCROSSSECTIONSMULTIPIONSANDRESONANCES_HH
#define G4INCLCROSSSECTIONSMULTIPIONSANDRESONANCES_HH

#include "G4INCLCrossSectionsMultiPions.hh"
//#include <limits>

namespace G4INCL {
  /// \brief Multipion and mesonic Resonances cross sections

  class CrossSectionsMultiPionsAndResonances : public CrossSectionsMultiPions {
    public:
      CrossSectionsMultiPionsAndResonances();
	  
      /// \brief new elastic particle-particle cross section
      virtual double elastic(Particle const * const p1, Particle const * const p2);
	  
      /// \brief new total particle-particle cross section
      virtual double total(Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for X pion production - piN Channel (modified due to the mesonic resonances)
      virtual double piNToxPiN(const int xpi, Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross sections for mesonic resonance production - piN Channel
      virtual double piNToEtaN(Particle const * const p1, Particle const * const p2);
      virtual double piNToOmegaN(Particle const * const p1, Particle const * const p2);
      virtual double piNToEtaPrimeN(Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross sections for mesonic resonance absorption on nucleon - piN Channel
      virtual double etaNToPiN(Particle const * const p1, Particle const * const p2);
      virtual double omegaNToPiN(Particle const * const p1, Particle const * const p2);
      virtual double etaPrimeNToPiN(Particle const * const p1, Particle const * const p2);

			/// \brief Cross sections for mesonic resonance absorption on nucleon - pipiN Channel
			virtual double etaNToPiPiN(Particle const * const p1, Particle const * const p2);

      /// \brief Cross sections for eta-N --> LK/SK
      virtual double etaNToLK(Particle const * const p1, Particle const * const p2);
      virtual double etaNToSK(Particle const * const p1, Particle const * const p2);            

      /// \brief Cross section for Eta production (inclusive) - NN entrance channel
      virtual double NNToNNEta(Particle const * const particle1, Particle const * const particle2);
			
      /// \brief Cross section for Eta production  (exclusive) - NN entrance channel
      virtual double NNToNNEtaExclu(Particle const * const particle1, Particle const * const particle2);
			
      /// \brief Cross section for Omega production (inclusive) - NN entrance channel
      virtual double NNToNNOmega(Particle const * const particle1, Particle const * const particle2);
			
      /// \brief Cross section for Omega production  (exclusive) - NN entrance channel
      virtual double NNToNNOmegaExclu(Particle const * const particle1, Particle const * const particle2);
	  
      /// \brief Cross section for X pion production - NN Channel
      virtual double NNToxPiNN(const int xpi, Particle const * const p1, Particle const * const p2);
	  
		   	/// \brief Cross section for X pion production - NNEta Channel
		   	virtual double NNToNNEtaxPi(const int xpi, Particle const * const p1, Particle const * const p2);
	  
			   /// \brief Cross section for N-Delta-Eta production - NNEta Channel
			   virtual double NNToNDeltaEta(Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for X pion production - NNOmega Channel
      virtual double NNToNNOmegaxPi(const int xpi, Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for N-Delta-Eta production - NNOmega Channel
      virtual double NNToNDeltaOmega(Particle const * const p1, Particle const * const p2);
			
	  
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
	  
      /// \brief One over threshold for s11pz
      static const double s11pzOOT;
      /// \brief One over threshold for s01pp
      static const double s01ppOOT;
      /// \brief One over threshold for s01pz
      static const double s01pzOOT;
      /// \brief One over threshold for s11pm
      static const double s11pmOOT;
      /// \brief One over threshold for s12pm
      static const double s12pmOOT;
      /// \brief One over threshold for s12pp
      static const double s12ppOOT;
      /// \brief One over threshold for s12zz
      static const double s12zzOOT;
      /// \brief One over threshold for s02pz
      static const double s02pzOOT;
      /// \brief One over threshold for s02pm
      static const double s02pmOOT;
      /// \brief One over threshold for s12mz
      static const double s12mzOOT;
	  
	  
      /// \brief Internal function for pion cross sections
	  double piMinuspToEtaN(Particle const * const p1, Particle const * const p2);
	  double piMinuspToEtaN(const double ECM);
	  double piMinuspToOmegaN(Particle const * const p1, Particle const * const p2);
	  double piMinuspToOmegaN(const double ECM);
//      double piPluspOnePi(Particle const * const p1, Particle const * const p2);
//	  double piMinuspOnePi(Particle const * const p1, Particle const * const p2);
//      double piPluspTwoPi(Particle const * const p1, Particle const * const p2);
//	  double piMinuspTwoPi(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for One (more) pion production - piN entrance channel
//      virtual double piNOnePi(Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for Two (more) pion production - piN entrance channel
//      virtual double piNTwoPi(Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for Three (more) pion production - piN entrance channel
      ///virtual double piNThreePi(Particle const * const p1, Particle const * const p2);
	  
		   	/// \brief Isotopic Cross section for Eta production (inclusive) - NN entrance channel
		   	virtual double NNToNNEtaIso(const double ener, const int iso);
	  
		   /// \brief Isotopic Cross section for Eta production (exclusive) - NN entrance channel
		   	virtual double NNToNNEtaExcluIso(const double ener, const int iso);
			
		   	/// \brief Cross section for direct 1-pion production - NNEta  channel
			   virtual double NNToNNEtaOnePi(Particle const * const part1, Particle const * const part2);
						/// \brief Cross section for direct 1-pion production - NNEta  channel
		   	virtual double NNToNNEtaOnePiOrDelta(Particle const * const part1, Particle const * const part2);
			   /// \brief Cross section for direct 2-pion production - NNEta  channel
		   	virtual double NNToNNEtaTwoPi(Particle const * const part1, Particle const * const part2);
   			/// \brief Cross section for direct 3-pion production - NNEta  channel
	   		virtual double NNToNNEtaThreePi(Particle const * const part1, Particle const * const part2);
	   		/// \brief Cross section for direct 4-pion production - NNEta  channel
	   		virtual double NNToNNEtaFourPi(Particle const * const part1, Particle const * const part2);

	  
      /// \brief Isotopic Cross section for Omega production (inclusive) - NN entrance channel
      virtual double NNToNNOmegaIso(const double ener, const int iso);
	  
      /// \brief Isotopic Cross section for Omega production (exclusive) - NN entrance channel
      virtual double NNToNNOmegaExcluIso(const double ener, const int iso);
			
      /// \brief Cross section for direct 1-pion production - NNOmega  channel
      virtual double NNToNNOmegaOnePi(Particle const * const part1, Particle const * const part2);
      /// \brief Cross section for direct 1-pion production - NNOmega  channel
      virtual double NNToNNOmegaOnePiOrDelta(Particle const * const part1, Particle const * const part2);
      /// \brief Cross section for direct 2-pion production - NNOmega  channel
      virtual double NNToNNOmegaTwoPi(Particle const * const part1, Particle const * const part2);
      /// \brief Cross section for direct 3-pion production - NNOmega  channel
      virtual double NNToNNOmegaThreePi(Particle const * const part1, Particle const * const part2);
      /// \brief Cross section for direct 4-pion production - NNOmega  channel
      virtual double NNToNNOmegaFourPi(Particle const * const part1, Particle const * const part2);
   
			
		   	/// \brief Cross sections for mesonic resonance absorption on nucleon - elastic Channel
		   	virtual double etaNElastic(Particle const * const p1, Particle const * const p2);					
      virtual double omegaNElastic(Particle const * const p1, Particle const * const p2);					

			
      /// \brief Cross sections for mesonic resonance absorption on nucleon - inelastic Channel
      virtual double omegaNInelastic(Particle const * const p1, Particle const * const p2);					
			
      /// \brief Cross sections for omega-induced 2Pi emission on nucleon
      virtual double omegaNToPiPiN(Particle const * const p1, Particle const * const p2);     
      virtual double omegaNToLK(Particle const * const p1, Particle const * const p2);       
      virtual double omegaNToSK(Particle const * const p1, Particle const * const p2);  					
   
		};
}

#endif
