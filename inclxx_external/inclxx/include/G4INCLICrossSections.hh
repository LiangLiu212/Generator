/** \file G4INCLICrossSections.hh
 * \brief Abstract interface for the cross-section classes
 *
 * \date 25nd October 2013
 * \author Davide Mancusi
 */

#ifndef G4INCLICROSSSECTIONS_HH
#define G4INCLICROSSSECTIONS_HH

#include "G4INCLParticle.hh"

namespace G4INCL {
  /// \brief Abstract interface for the cross-section classes
  class ICrossSections {
    public:

      ICrossSections() {}
      virtual ~ICrossSections() {}

      /// \brief Elastic particle-particle cross section
      virtual double elastic(Particle const * const p1, Particle const * const p2) = 0;

      /// \brief Total (elastic+inelastic) particle-particle cross section
      virtual double total(Particle const * const p1, Particle const * const p2) = 0;

      /// \brief Cross section for NDelta->NN
      virtual double NDeltaToNN(Particle const * const p1, Particle const * const p2) = 0;

      /// \brief Cross section for NN->NDelta
      virtual double NNToNDelta(Particle const * const p1, Particle const * const p2) = 0;

      /// \brief Cross section for NN->xpiN
      virtual double NNToxPiNN(const int xpi, Particle const * const p1, Particle const * const p2) = 0;

      /// \brief Cross section for piN->NDelta
      virtual double piNToDelta(Particle const * const p1, Particle const * const p2) = 0;

      /// \brief Cross section for piN->piNpi
      virtual double piNToxPiN(const int xpi, Particle const * const p1, Particle const * const p2) = 0;
    
      /// \brief Cross section for PiN->EtaN
      virtual double piNToEtaN(Particle const * const p1, Particle const * const p2) = 0;
    
      /// \brief Cross section for PiN->OmegaN
      virtual double piNToOmegaN(Particle const * const p1, Particle const * const p2) = 0;
    
      /// \brief Cross section for PiN->EtaPrimeN
      virtual double piNToEtaPrimeN(Particle const * const p1, Particle const * const p2) = 0;
    
      /// \brief Cross section for EtaN->PiN
      virtual double etaNToPiN(Particle const * const p1, Particle const * const p2) = 0;
    
      /// \brief Cross section for EtaN->PiPiN
      virtual double etaNToPiPiN(Particle const * const p1, Particle const * const p2) = 0;
   
      /// \brief Cross section for EtaN->LK/SK
      virtual double etaNToLK(Particle const * const p1, Particle const * const p2) = 0;
      virtual double etaNToSK(Particle const * const p1, Particle const * const p2) = 0;
    
      /// \brief Cross section for OmegaN->PiN
      virtual double omegaNToPiN(Particle const * const p1, Particle const * const p2) = 0;
   
      /// \brief Cross section for OmegaN->PiPiN
      virtual double omegaNToPiPiN(Particle const * const p1, Particle const * const p2) = 0;
   
      /// \brief Cross section for OmegaN->LK/SK
      virtual double omegaNToLK(Particle const * const p1, Particle const * const p2) = 0;
      virtual double omegaNToSK(Particle const * const p1, Particle const * const p2) = 0;
   
      /// \brief Cross section for EtaPrimeN->PiN
      virtual double etaPrimeNToPiN(Particle const * const p1, Particle const * const p2) = 0;
    
      /// \brief Cross section for NN->NNEta (inclusive)
      virtual double NNToNNEta(Particle const * const p1, Particle const * const p2) = 0;
    
      /// \brief Cross section for NN->NNEta (exclusive)
      virtual double NNToNNEtaExclu(Particle const * const p1, Particle const * const p2) = 0;
      
      /// \brief Cross section for NN->NNEtaxPi
      virtual double NNToNNEtaxPi(const int xpi, Particle const * const p1, Particle const * const p2) = 0;
    
      /// \brief Cross section for N-Delta-Eta production - NNEta Channel
      virtual double NNToNDeltaEta(Particle const * const p1, Particle const * const p2) = 0;
         
      /// \brief Cross section for NN->NNEta (inclusive)
      virtual double NNToNNOmega(Particle const * const p1, Particle const * const p2) = 0;
    
      /// \brief Cross section for NN->NNEta (exclusive)
      virtual double NNToNNOmegaExclu(Particle const * const p1, Particle const * const p2) = 0;
      
      /// \brief Cross section for NN->NNEtaxPi
      virtual double NNToNNOmegaxPi(const int xpi, Particle const * const p1, Particle const * const p2) = 0;
    
      /// \brief Cross section for N-Delta-Eta production - NNEta Channel
      virtual double NNToNDeltaOmega(Particle const * const p1, Particle const * const p2) = 0;
      
      
      /// \brief elastic scattering for Nucleon-Strange Particles cross sections
      virtual double NYelastic(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NKbelastic(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NKelastic(Particle const * const p1, Particle const * const p2) = 0;
      
      /// \brief Nucleon-Nucleon to Stange particles cross sections
      virtual double NNToNLK(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NNToNSK(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NNToNLKpi(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NNToNSKpi(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NNToNLK2pi(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NNToNSK2pi(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NNToNNKKb(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NNToMissingStrangeness(Particle const * const p1, Particle const * const p2) = 0;
      
      /// \brief Nucleon-Delta to Stange particles cross sections
      virtual double NDeltaToNLK(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NDeltaToNSK(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NDeltaToDeltaLK(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NDeltaToDeltaSK(Particle const * const p1, Particle const * const p2) = 0;
      
      virtual double NDeltaToNNKKb(Particle const * const p1, Particle const * const p2) = 0;
         
      /// \brief Nucleon-Pion to Stange particles cross sections
      virtual double NpiToLK(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NpiToSK(Particle const * const p1, Particle const * const p2) = 0;
      virtual double p_pimToSzKz(Particle const * const p1, Particle const * const p2) = 0;
      virtual double p_pimToSmKp(Particle const * const p1, Particle const * const p2) = 0;
      virtual double p_pizToSzKp(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NpiToLKpi(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NpiToSKpi(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NpiToLK2pi(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NpiToSK2pi(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NpiToNKKb(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NpiToMissingStrangeness(Particle const * const p1, Particle const * const p2) = 0;
      
      /// \brief Nucleon-Hyperon cross sections
      virtual double NLToNS(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NSToNL(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NSToNS(Particle const * const p1, Particle const * const p2) = 0;
         
      /// \brief Nucleon-Kaon inelastic cross sections
      virtual double NKToNK(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NKToNKpi(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NKToNK2pi(Particle const * const p1, Particle const * const p2) = 0;
      
      /// \brief Nucleon-antiKaon inelastic cross sections
      virtual double NKbToNKb(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NKbToSpi(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NKbToLpi(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NKbToS2pi(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NKbToL2pi(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NKbToNKbpi(Particle const * const p1, Particle const * const p2) = 0;
      virtual double NKbToNKb2pi(Particle const * const p1, Particle const * const p2) = 0;
      
      /// \brief NNbar Particles cross sections
      /// \brief Nucleon-AntiNucleon to Baryon-AntiBaryon cross sections
      virtual double NNbarElastic(Particle const* const p1, Particle const* const p2) = 0;
      virtual double NNbarCEX(Particle const* const p1, Particle const* const p2) = 0;
      
      virtual double NNbarToLLbar(Particle const * const p1, Particle const * const p2) = 0;
      
      /// \brief Nucleon-AntiNucleon to Nucleon-AntiNucleon + pions cross sections
      virtual double NNbarToNNbarpi(Particle const* const p1, Particle const* const p2) = 0;
      virtual double NNbarToNNbar2pi(Particle const* const p1, Particle const* const p2) = 0;
      virtual double NNbarToNNbar3pi(Particle const* const p1, Particle const* const p2) = 0;
     
      /// \brief Nucleon-AntiNucleon total annihilation cross sections
      virtual double NNbarToAnnihilation(Particle const* const p1, Particle const* const p2) = 0;

      /** \brief Calculate the slope of the NN DDXS.
       *
       * \param energyCM energy in the CM frame, in MeV
       * \param iso total isospin of the system
       *
       * \return the slope of the angular distribution
       */
      virtual double calculateNNAngularSlope(double energyCM, int iso) = 0;

  };
}

#endif
