/** \file G4INCLGlobalInfo.hh
 * \brief Simple container for output of calculation-wide results.
 *
 * Contains the final results of an INCL calculation.
 *
 * \date 21 January 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLGLOBALINFO_HH
#define G4INCLGLOBALINFO_HH 1

#ifdef INCL_ROOT_USE
#include <Rtypes.h>
#endif

#include <string>
#include <vector>

namespace G4INCL {
#ifndef INCL_ROOT_USE
    typedef int Int_t;
    typedef short Short_t;
    typedef float Float_t;
#endif

    struct GlobalInfo {
      GlobalInfo() :
#ifdef INCL_ROOT_USE

#endif
        Ap(0),
        Zp(0),
        Sp(0),
        At(0),
        Zt(0),
        St(0),
        Ep((Float_t)0.0),
        nShots(0),
        geometricCrossSection((Float_t)0.0),
        biasFactor((Float_t)0.0),
        nTransparents(0),
        reactionCrossSection((Float_t)0.0),
        errorReactionCrossSection((Float_t)0.0),
        nNucleonAbsorptions(0),
        nucleonAbsorptionCrossSection((Float_t)0.0),
        nPionAbsorptions(0),
        pionAbsorptionCrossSection((Float_t)0.0),
        nForcedTransparents(0),
        nForcedCompoundNucleus(0),
        forcedCNCrossSection((Float_t)0.0),
        errorForcedCNCrossSection((Float_t)0.0),
        nCompleteFusion(0),
        completeFusionCrossSection((Float_t)0.0),
        errorCompleteFusionCrossSection((Float_t)0.0),
        nEnergyViolationInteraction(0),
        energyViolationInteractionCrossSection((Float_t)0.0)
      {
#ifdef INCL_ROOT_USE

#endif

      }
#ifdef INCL_ROOT_USE
      /** \brief Selection string for an abridged version of the ROOT tree */
      std::string rootSelection;
#endif
      /** \brief Name of the cascade model */
      std::string cascadeModel;
      /** \brief Name of the de-excitation model */
      std::string deexcitationModel;
      /** \brief Projectile mass number given as input */
      Short_t Ap;
      /** \brief Projectile charge number given as input */
      Short_t Zp;
      /** \brief Projectile strangeness number given as input */
      Short_t Sp;
      /** \brief Target mass number given as input */
      Short_t At;
      /** \brief Target charge number given as input */
      Short_t Zt;
      /** \brief Target strangeness number given as input */
      Short_t St;
      /** \brief Projectile kinetic energy given as input */
      Float_t Ep;
      /** \brief Number of shots */
      Int_t nShots;
      /** \brief Geometric cross section */
      Float_t geometricCrossSection;
      /** \brief Bias factor */
      Float_t biasFactor;
      /** \brief Number of transparent shots */
      Int_t nTransparents;
      /** \brief Calculated reaction cross section */
      Float_t reactionCrossSection;
      /** \brief Error on the calculated reaction cross section */
      Float_t errorReactionCrossSection;
      /** \brief Number of nucleon absorptions (no outcoming particles) */
      Int_t nNucleonAbsorptions;
      /** \brief Nucleon absorption cross section */
      Float_t nucleonAbsorptionCrossSection;
      /** \brief Number of nucleon absorptions (no outcoming pions) */
      Int_t nPionAbsorptions;
      /** \brief Pion absorption cross section */
      Float_t pionAbsorptionCrossSection;
      /** \brief Number of forced transparents */
      Int_t nForcedTransparents;
      /** \brief Number of forced compound-nucleus events */
      Int_t nForcedCompoundNucleus;
      /** \brief Calculated forced-compound-nucleus cross section */
      Float_t forcedCNCrossSection;
      /** \brief Error on the calculated forced-compound-nucleus cross section */
      Float_t errorForcedCNCrossSection;
      /** \brief Number of complete-fusion events (nParticles==0) */
      Int_t nCompleteFusion;
      /** \brief Calculated complete-fusion cross section (nParticles==0) */
      Float_t completeFusionCrossSection;
      /** \brief Error on the calculated complete-fusion cross section (nParticles==0) */
      Float_t errorCompleteFusionCrossSection;
      /** \brief Number of attempted collisions/decays for which the energy-conservation algorithm failed to find a solution. */
      Int_t nEnergyViolationInteraction;
      /** \brief Cross section for attempted collisions/decays for which the energy-conservation algorithm failed to find a solution. */
      Float_t energyViolationInteractionCrossSection;
      /** \brief Initial seeds for the pseudo-random-number generator */
      std::vector<Int_t> initialRandomSeeds;
      /** \brief Final seeds for the pseudo-random-number generator */
      std::vector<Int_t> finalRandomSeeds;
    };
}

#endif /* G4INCLGLOBALINFO_HH */
