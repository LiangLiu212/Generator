/** \file G4INCLCrossSectionsINCL46.hh
 * \brief Cross sections used in INCL4.6
 *
 * \date 25nd October 2013
 * \author Davide Mancusi
 */

#ifndef G4INCLCROSSSECTIONSINCL46_HH
#define G4INCLCROSSSECTIONSINCL46_HH

#include "G4INCLICrossSections.hh"

namespace G4INCL {
  /// \brief Cross sections used in INCL4.6

  class CrossSectionsINCL46 : public ICrossSections{
    public:

      /// \brief Elastic particle-particle cross section
      virtual double elastic(Particle const * const p1, Particle const * const p2);

      /// \brief Total (elastic+inelastic) particle-particle cross section
      virtual double total(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for NDelta->NN
      virtual double NDeltaToNN(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for NN->NDelta
      virtual double NNToNDelta(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for piN->NDelta
      virtual double piNToDelta(Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for X pion production - NN Channel
      virtual double NNToxPiNN(const int xpi, Particle const * const p1, Particle const * const p2);

      /// \brief Cross section for X pion production - piN Channel
      virtual double piNToxPiN(const int xpi, Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross sections for mesonic resonance production - piN Channel
      virtual double piNToEtaN(Particle const * const p1, Particle const * const p2);
      virtual double piNToOmegaN(Particle const * const p1, Particle const * const p2);
      virtual double piNToEtaPrimeN(Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross sections for mesonic resonance absorption on nucleons - piN Channel
      virtual double etaNToPiN(Particle const * const p1, Particle const * const p2);
      virtual double omegaNToPiN(Particle const * const p1, Particle const * const p2);
      virtual double etaPrimeNToPiN(Particle const * const p1, Particle const * const p2);
      
      /// \brief Cross sections for mesonic resonance absorption on nucleon - pipiN Channel
      virtual double etaNToPiPiN(Particle const * const p1, Particle const * const p2);     
      virtual double omegaNToPiPiN(Particle const * const p1, Particle const * const p2); 
      
      /// \brief Cross sections for mesonic resonance absorption on nucleon - LK/SK Channel
      virtual double etaNToLK(Particle const * const p1, Particle const * const p2);     
      virtual double omegaNToLK(Particle const * const p1, Particle const * const p2);
      virtual double etaNToSK(Particle const * const p1, Particle const * const p2);     
      virtual double omegaNToSK(Particle const * const p1, Particle const * const p2); 		
			
	    /// \brief Cross section for Eta production - NN entrance channel
      virtual double NNToNNEta(Particle const * const p1, Particle const * const p2);
			
			/// \brief Cross section for Eta production  (exclusive) - NN entrance channel
			virtual double NNToNNEtaExclu(Particle const * const p1, Particle const * const p2);
	  
			/// \brief Cross section for X pion production - NNEta Channel
		  virtual double NNToNNEtaxPi(const int xpi, Particle const * const p1, Particle const * const p2);

		  /// \brief Cross section for N-Delta-Eta production - NNEta Channel
			virtual double NNToNDeltaEta(Particle const * const p1, Particle const * const p2);
      
      
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
      
      /// \brief Nucleon-Kaon inelastic cross sections
      virtual double NKToNK(Particle const * const p1, Particle const * const p2);
      virtual double NKToNKpi(Particle const * const p1, Particle const * const p2);
      virtual double NKToNK2pi(Particle const * const p1, Particle const * const p2);
      
      /// \brief Nucleon-antiKaon inelastic cross sections
      virtual double NKbToNKb(Particle const * const p1, Particle const * const p2);
      virtual double NKbToSpi(Particle const * const p1, Particle const * const p2);
      virtual double NKbToLpi(Particle const * const p1, Particle const * const p2);
      virtual double NKbToS2pi(Particle const * const p1, Particle const * const p2);
      virtual double NKbToL2pi(Particle const * const p1, Particle const * const p2);
      virtual double NKbToNKbpi(Particle const * const p1, Particle const * const p2);
      virtual double NKbToNKb2pi(Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for Eta production - NN entrance channel
      virtual double NNToNNOmega(Particle const * const particle1, Particle const * const particle2);
			
      /// \brief Cross section for Eta production  (exclusive) - NN entrance channel
      virtual double NNToNNOmegaExclu(Particle const * const particle1, Particle const * const particle2);
	  
      /// \brief Cross section for X pion production - NNEta Channel
      virtual double NNToNNOmegaxPi(const int xpi, Particle const * const p1, Particle const * const p2);
	  
      /// \brief Cross section for N-Delta-Eta production - NNEta Channel
      virtual double NNToNDeltaOmega(Particle const * const p1, Particle const * const p2);
      
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
	  
      /** \brief Calculate the slope of the NN DDXS.
       *
       * \param energyCM energy in the CM frame, in MeV
       * \param iso total isospin of the system
       *
       * \return the slope of the angular distribution
       */
      virtual double calculateNNAngularSlope(double energyCM, int iso);

    protected:

      /// \brief Internal implementation of the elastic cross section
      double elasticNNLegacy(Particle const * const part1, Particle const * const part2);

      /// \brief Internal function for the delta-production cross section
      double deltaProduction(const int isospin, const double pLab);

      double spnPiPlusPHE(const double x);
      double spnPiMinusPHE(const double x);
	  
  };
}

#endif
