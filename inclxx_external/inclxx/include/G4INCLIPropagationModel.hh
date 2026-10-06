/*
 * IPropagationModel.hh
 *
 *  \date 4 juin 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLIPropagationModel_hh
#define G4INCLIPropagationModel_hh

#include "G4INCLIAvatar.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLFinalState.hh"

namespace G4INCL {

    /**
     * Propagation model takes care of transporting the particles until something
     * interesting (i.e. an avatar) happens. This avatar is then returned back to the INCL
     * kernel for further processing.
     *
     * The propagation model idea abstracts the details of propagation. This allows us to
     * conveniently support multiple propagation models and to compare their results. Some
     * possible future propagation models are: straight line trajectories by using constant
     * time step and curved trajectories.
     */
    class IPropagationModel {
    public:
      IPropagationModel() {}
      virtual ~IPropagationModel() {}

      /**
       * Set the nucleus for the propagation model.
       *
       * @param nucleus Pointer to the nucleus
       */
      virtual void setNucleus(G4INCL::Nucleus *nucleus) = 0;

      /**
       * Get a pointer to the nucleus.
       *
       * @return G4INCL::Nuleus*
       */
      virtual G4INCL::Nucleus* getNucleus() = 0;

      virtual double shoot(ParticleSpecies const &projectileSpecies, const double kineticEnergy, const double impactParameter, const double phi) = 0;
    protected:
      virtual double shootAtrest(ParticleType const t, const double kineticEnergy) = 0; 
      virtual double shootParticle(ParticleType const t, const double kineticEnergy, const double impactParameter, const double phi) = 0;
      virtual double shootComposite(ParticleSpecies const &s, const double kineticEnergy, const double impactParameter, const double phi) = 0;
      virtual double shootCompositeAtrest(ParticleSpecies const &s, const double kineticEnergy) = 0;

    public:

      /**
       * Returns the current global time of the system.
       */
      virtual double getCurrentTime() = 0;

      /**
       * Set new stopping time to the propagation.
       */
      virtual void setStoppingTime(double) = 0;

      /**
       * Get the current stopping time.
       */
      virtual double getStoppingTime() = 0;

      /**
       * Propagate the particles and get the next avatar.
       *
       * @return G4INCL::IAvatar the next avatar
       */
      virtual G4INCL::IAvatar* propagate(FinalState const * const fs) = 0;
    };

}

#endif
