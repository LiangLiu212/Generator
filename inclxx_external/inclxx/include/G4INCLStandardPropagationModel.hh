/*
 * StandardPropagationModel.hh
 *
 *  \date 4 June 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLStandardPropagationModel_hh
#define G4INCLStandardPropagationModel_hh 1

#include "G4INCLNucleus.hh"
#include "G4INCLIPropagationModel.hh"
#include "G4INCLIAvatar.hh"
#include "G4INCLConfigEnums.hh"

#include <iterator>

namespace G4INCL {

    /**
     * Standard INCL4 particle propagation and avatar prediction
     *
     * This class implements the standard INCL4 avatar prediction and particle
     * propagation logic. The main idea is to predict all collisions between particles
     * and their reflections from the potential wall. After this we select the avatar
     * with the smallest time, propagate all particles to their positions at that time
     * and return the avatar to the INCL kernel @see G4INCL::Kernel.
     *
     * The particle trajectories in this propagation model are straight lines and all
     * particles are assumed to move with constant velocity.
     */
    class StandardPropagationModel: public G4INCL::IPropagationModel {
    public:
      StandardPropagationModel(LocalEnergyType localEnergyType, LocalEnergyType localEnergyDeltaType, const double hTime = 0.0);
      virtual ~StandardPropagationModel();

      double getCurrentTime();
      /**
       * Set the nucleus for this propagation model.
       */
      void setNucleus(G4INCL::Nucleus *nucleus);

      /**
       * Get the nucleus.
       */
      G4INCL::Nucleus* getNucleus();

      double shoot(ParticleSpecies const &projectileSpecies, const double kineticEnergy, const double impactParameter, const double phi);
      double shootParticle(ParticleType const t, const double kineticEnergy, const double impactParameter, const double phi);
      double shootComposite(ParticleSpecies const &s, const double kineticEnergy, const double impactParameter, const double phi);
      double shootAtrest(ParticleType const t, const double kineticEnergy);
      double shootCompositeAtrest(ParticleSpecies const &s, const double kineticEnergy); 
      
      /**
       * Set the stopping time of the simulation.
       */
      void setStoppingTime(double);

      /**
       * Get the current stopping time.
       */
      double getStoppingTime();

      /**
       * Add an avatar to the storage.
       */
      void registerAvatar(G4INCL::IAvatar *anAvatar);

      /** \brief Generate a two-particle avatar.
       *
       * Generate a two-particle avatar, if all the appropriate conditions are
       * met.
       */
      IAvatar *generateBinaryCollisionAvatar(Particle * const p1, Particle * const p2);

      /** \brief Get the reflection time.
       *
       * Returns the reflection time of a particle on the potential wall.
       *
       * \param aParticle pointer to the particle
       */
      double getReflectionTime(G4INCL::Particle const * const aParticle);

      /**
       * Get the predicted time of the collision between two particles.
       */
      double getTime(G4INCL::Particle const * const particleA,
		     G4INCL::Particle const * const particleB, double *minDistOfApproach) const;

      /** \brief Generate and register collisions between a list of updated particles and all the other particles.
       *
       * This method does not generate collisions among the particles in
       * updatedParticles; in other words, it generates a collision between one
       * of the updatedParticles and one of the particles ONLY IF the latter
       * does not belong to updatedParticles.
       *
       * If you intend to generate all possible collisions among particles in a
       * list, use generateCollisions().
       *
       * \param updatedParticles list of updated particles
       * \param particles list of particles
       */
      void generateUpdatedCollisions(const ParticleList &updatedParticles, const ParticleList &particles);

      /** \brief Generate and register collisions among particles in a list, except between those in another list.
       *
       * This method generates all possible collisions among the particles.
       * Each collision is generated only once.
       *
       * \param particles list of particles
       */
      void generateCollisions(const ParticleList &particles);

      /** \brief Generate and register collisions among particles in a list, except between those in another list.
       *
       * This method generates all possible collisions among the particles.
       * Each collision is generated only once. The collision is NOT generated
       * if BOTH collision partners belong to the except list.
       *
       * You should pass an empty list as the except parameter if you want to
       * generate all possible collisions among particles.
       *
       * \param particles list of particles
       * \param except list of excluded particles
       */
      void generateCollisions(const ParticleList &particles, const ParticleList &except);

      /** \brief Generate decays for particles that can decay.
       *
       * The list of particles given as an argument is allowed to contain also
       * stable particles.
       *
       * \param particles list of particles to (possibly) generate decays for
       */
      void generateDecays(const ParticleList &particles);

      /**
       * Update all avatars related to a particle.
       */
      void updateAvatars(const ParticleList &particles);

      /// \brief (Re)Generate all possible avatars.
      void generateAllAvatars();

#ifdef INCL_REGENERATE_AVATARS
      /** \brief (Re)Generate all possible avatars.
       *
       * This method excludes collision avatars between updated particles.
       */
      void generateAllAvatarsExceptUpdated(FinalState const * const fs);
#endif

      /**
       * Propagate all particles and return the first avatar.
       */
      G4INCL::IAvatar* propagate(FinalState const * const fs);

    private:
      G4INCL::Nucleus *theNucleus;
      double maximumTime;
      double currentTime;
      double hadronizationTime;
      bool firstAvatar;
      LocalEnergyType theLocalEnergyType, theLocalEnergyDeltaType;
      Particle backupParticle1, backupParticle2;
    };

}


#endif
