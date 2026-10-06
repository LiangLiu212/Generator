/** \file G4INCLCrossSectionsMultiPions.hh
 * \brief Cross sections used in INCL Multipions
 *
 * \date 26th November 2013
 * \author Jean-Christophe David
 */

#ifndef G4INCLCROSSSECTIONSMULTIPIONS_HH
#define G4INCLCROSSSECTIONSMULTIPIONS_HH

#include "G4INCLICrossSections.hh"
#include "G4INCLHornerFormEvaluator.hh"

namespace G4INCL {
  /// \brief Cross sections used in INCL Multipions

  class CrossSectionsMultiPions : public ICrossSections{
    public:
      CrossSectionsMultiPions();

      /// \brief Elastic particle-particle cross section
      virtual double elastic(Particle const * const p1, Particle const * const p2);

      /// \brief Total (elastic+inelastic) particle-particle cross section
      virtual double total(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for NDelta->NN
      virtual double NDeltaToNN(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for Delta production - NN Channel
      virtual double NNToNDelta(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for Delta production - piN Channel
      virtual double piNToDelta(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for X pion production - piN Channel
      virtual double piNToxPiN(const int xpi, Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for X pion production - NN Channel
      virtual double NNToxPiNN(const int xpi, Particle const * const p1, Particle const * const p2);

      /** \brief Calculate the slope of the NN DDXS.
       *
       * \param energyCM energy in the CM frame, in MeV
       * \param iso total isospin of the system
       *
       * \return the slope of the angular distribution
       */
      virtual double calculateNNAngularSlope(double energyCM, int iso);
	  
      /// \brief Cross sections for mesonic resonance production - piN Channel
      virtual double piNToEtaN(Particle const * const p1, Particle const * const p2);
      virtual double piNToOmegaN(Particle const * const p1, Particle const * const p2);
      virtual double piNToEtaPrimeN(Particle const * const p1, Particle const * const p2);
			
			/// \brief Cross sections for mesonic resonance absorption on nucleon - pipiN Channel
      virtual double etaNToPiPiN(Particle const * const p1, Particle const * const p2);			
      virtual double omegaNToPiPiN(Particle const * const p1, Particle const * const p2);	  
			
      /// \brief Cross sections for mesonic resonance absorption on nucleon - piN Channel
      virtual double etaNToPiN(Particle const * const p1, Particle const * const p2);
      virtual double omegaNToPiN(Particle const * const p1, Particle const * const p2);
      virtual double etaPrimeNToPiN(Particle const * const p1, Particle const * const p2); 
      
      /// \brief Cross sections for mesonic resonance absorption on nucleon - LK/SKChannel
      virtual double etaNToLK(Particle const * const p1, Particle const * const p2);     
      virtual double omegaNToLK(Particle const * const p1, Particle const * const p2);
      virtual double etaNToSK(Particle const * const p1, Particle const * const p2);     
      virtual double omegaNToSK(Particle const * const p1, Particle const * const p2);  
	  
	    /// \brief Cross section for Eta production - NN entrance channel
      virtual double NNToNNEta(Particle const * const particle1, Particle const * const particle2);
			
		  /// \brief Cross section for Eta production  (exclusive) - NN entrance channel
		  virtual double NNToNNEtaExclu(Particle const * const particle1, Particle const * const particle2);
	  
		  /// \brief Cross section for X pion production - NNEta Channel
			virtual double NNToNNEtaxPi(const int xpi, Particle const * const p1, Particle const * const p2);
	  
		  /// \brief Cross section for N-Delta-Eta production - NNEta Channel
		  virtual double NNToNDeltaEta(Particle const * const p1, Particle const * const p2);			
	  
      /// \brief Cross section for Eta production - NN entrance channel
      virtual double NNToNNOmega(Particle const * const particle1, Particle const * const particle2);
			
      /// \brief Cross section for Eta production  (exclusive) - NN entrance channel
      virtual double NNToNNOmegaExclu(Particle const * const particle1, Particle const * const particle2);
	  
      /// \brief Cross section for X pion production - NNEta Channel
      virtual double NNToNNOmegaxPi(const int xpi, Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for N-Delta-Eta production - NNEta Channel
      virtual double NNToNDeltaOmega(Particle const * const p1, Particle const * const p2);
      
      
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
      
      /// \brief Nucleon-Pion to Stange particles cross sections
      virtual double NpiToLK(Particle const * const p1, Particle const * const p2);
      virtual double NpiToSK(Particle const * const p1, Particle const * const p2);
		  virtual double p_pimToSzKz(Particle const * const p1, Particle const * const p2);
		  virtual double p_pimToSmKp(Particle const * const p1, Particle const * const p2);
		  virtual double p_pizToSzKp(Particle const * const p1, Particle const * const p2);
      virtual double NpiToLKpi(Particle const * const p1, Particle const * const p2);
      virtual double NpiToSKpi(Particle const * const p1, Particle const * const p2);
      virtual double NpiToLK2pi(Particle const * const p1, Particle const * const p2);
      virtual double NpiToSK2pi(Particle const * const p1, Particle const * const p2);
      virtual double NpiToNKKb(Particle const * const p1, Particle const * const p2);
      
      virtual double NpiToMissingStrangeness(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-Hyperon cross sections
      virtual double NLToNS(Particle const * const p1, Particle const * const p2);
      virtual double NSToNL(Particle const * const p1, Particle const * const p2);
      virtual double NSToNS(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-Kaon quasi-elastic and inelastic cross sections
      virtual double NKToNK(Particle const * const p1, Particle const * const p2);
      virtual double NKToNKpi(Particle const * const p1, Particle const * const p2);
      virtual double NKToNK2pi(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-antiKaon quasi-elastic and inelastic cross sections
      virtual double NKbToNKb(Particle const * const p1, Particle const * const p2);
      virtual double NKbToSpi(Particle const * const p1, Particle const * const p2);
      virtual double NKbToLpi(Particle const * const p1, Particle const * const p2);
      virtual double NKbToS2pi(Particle const * const p1, Particle const * const p2);
      virtual double NKbToL2pi(Particle const * const p1, Particle const * const p2);
      virtual double NKbToNKbpi(Particle const * const p1, Particle const * const p2);
      virtual double NKbToNKb2pi(Particle const * const p1, Particle const * const p2);
      
      /// \brief antiparticle cross sections
      /// \brief Nucleon-AntiNucleon to Baryon-AntiBaryon cross sections
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

      /// \brief Internal implementation of the NN elastic cross section
      double NNElastic(Particle const * const part1, Particle const * const part2);

      /// \brief Internal implementation of the NN elastic cross section with fixed isospin
      double NNElasticFixed(const double s, const int i);

      /// \brief Internal implementation of the NN total cross section
      double NNTot(Particle const * const part1, Particle const * const part2);

      /// \brief Internal implementation of the NN total cross section with fixed isospin
      double NNTotFixed(const double s, const int i);

      /// \brief Internal implementation of the isospin dependent NN reaction cross section
      double NNInelasticIso(const double ener, const int iso);

      /// \brief Cross section for direct 1-pion production + delta production - NN entrance channel
      virtual double NNOnePiOrDelta(const double ener, const int iso, const double xsiso);
      /// \brief Cross section for direct 2-pion production - NN entrance channel
      virtual double NNTwoPi(const double ener, const int iso, const double xsiso);
      /// \brief Cross section for direct 3-pion production - NN entrance channel
      virtual double NNThreePi(const double ener, const int iso, const double xsiso, const double xs1pi, const double xs2pi);

      /// \brief Cross section for direct 1-pion production - NN entrance channel
      virtual double NNOnePi(Particle const * const part1, Particle const * const part2);
      /// \brief Cross section for direct 1-pion production - NN entrance channel
      virtual double NNOnePiOrDelta(Particle const * const part1, Particle const * const part2);
      /// \brief Cross section for direct 2-pion production - NN entrance channel
      virtual double NNTwoPi(Particle const * const part1, Particle const * const part2);
      /// \brief Cross section for direct 3-pion production - NN entrance channel
      virtual double NNThreePi(Particle const * const part1, Particle const * const part2);
      /// \brief Cross section for direct 4-pion production - NN entrance channel
      virtual double NNFourPi(Particle const * const part1, Particle const * const part2);

      /// \brief Internal function for pion cross sections
      double spnPiPlusPHE(const double x);
      /// \brief Internal function for pion cross sections
      double spnPiMinusPHE(const double x);
      double piNIne(Particle const * const p1, Particle const * const p2);
      double piNTot(Particle const * const p1, Particle const * const p2);
      double piNTopiN(Particle const * const p1, Particle const * const p2);
      double piPluspIne(Particle const * const p1, Particle const * const p2);
	  double piMinuspIne(Particle const * const p1, Particle const * const p2);
      double piPluspOnePi(Particle const * const p1, Particle const * const p2);
	  double piMinuspOnePi(Particle const * const p1, Particle const * const p2);
      double piPluspTwoPi(Particle const * const p1, Particle const * const p2);
	  double piMinuspTwoPi(Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for One (more) pion production - piN entrance channel
      virtual double piNOnePi(Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for Two (more) pion production - piN entrance channel
      virtual double piNTwoPi(Particle const * const p1, Particle const * const p2);
	  
	  
  };
}

#endif
