/* \file G4INCLInteractionAvatar.hh
 * \brief Virtual class for interaction avatars.
 *
 * This class is inherited by decay and collision avatars. The goal is to
 * provide a uniform treatment of common physics, such as Pauli blocking,
 * enforcement of energy conservation, etc.
 *
 *  \date Mar 1st, 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLINTERACTIONAVATAR_HH_
#define G4INCLINTERACTIONAVATAR_HH_

#include "G4INCLIAvatar.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLRootFinder.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {

  class InteractionAvatar : public G4INCL::IAvatar {
    public:
      InteractionAvatar(double, G4INCL::Nucleus*, G4INCL::Particle*);
      InteractionAvatar(double, G4INCL::Nucleus*, G4INCL::Particle*, G4INCL::Particle*);
      virtual ~InteractionAvatar();

      /// \brief Target accuracy in the determination of the local-energy Q-value
      static const double locEAccuracy;
      /// \brief Max number of iterations for the determination of the local-energy Q-value
      static const int maxIterLocE;

      /// \brief Release the memory allocated for the backup particles
      static void deleteBackupParticles();
      
    /**
     * static instance
     */
    static InteractionAvatar* Instance();
    
    void setSrcPartner(Particle *p /*, const ThreeVector m*/);
    
      /** \brief Apply local-energy transformation, if appropriate
       *
       * \param p particle to apply the transformation to
       */
      void preInteractionLocalEnergy(Particle * const p);
      
      ThreeVector getboostVector(){return boostVector;}
      
      void setboostVector(ThreeVector& v){boostVector = v;}

    protected:
      virtual G4INCL::IChannel* getChannel() = 0;

      bool bringParticleInside(Particle * const p);
      
      EventInfo theEventInfo;

      /** \brief Store the state of the particles before the interaction
       *
       * If the interaction cannot be realised for any reason, we will need to
       * restore the particle state as it was before. This is done by calling
       * the restoreParticles() method.
       */
      void preInteractionBlocking();

      void preInteraction();
      void postInteraction(FinalState *);

      /** \brief Restore the state of both particles.
       *
       * The state must first be stored by calling preInteractionBlocking().
       */
      void restoreParticles() const;
      
      void restoreSrcPartner(FinalState * fs);

      /// \brief true if the given avatar should use local energy
      bool shouldUseLocalEnergy() const;

      Nucleus *theNucleus;
      Particle *particle1, *particle2;
      static G4ThreadLocal Particle *backupParticle1, *backupParticle2;
      ThreeVector boostVector;
      double oldTotalEnergy, oldXSec;
      bool isPiN;
      double weight;

    private:
      static G4ThreadLocal InteractionAvatar* interactionAvatar;
      static G4ThreadLocal Particle *backupPartner;
      static ThreeVector mbackupPartner;
      
      /// \brief RootFunctor-derived object for enforcing energy conservation in N-N.
      class ViolationEMomentumFunctor : public RootFunctor {
        public:
          /** \brief Prepare for calling the () operator and scaleParticleMomenta
           *
           * The constructor sets the private class members.
           */
          ViolationEMomentumFunctor(Nucleus * const nucleus, ParticleList const &modAndCre, const double totalEnergyBeforeInteraction, ThreeVector const &boost, const bool localE);
          virtual ~ViolationEMomentumFunctor();

          /** \brief Compute the energy-conservation violation.
           *
           * \param x scale factor for the particle momenta
           * \return the energy-conservation violation
           */
          double operator()(const double x) const;

          /// \brief Clean up after root finding
          void cleanUp(const bool success) const;

        private:
          /// \brief List of final-state particles.
          ParticleList finalParticles;
          /// \brief CM particle momenta, as determined by the channel.
          std::vector<ThreeVector> particleMomenta;
          /// \brief Total energy before the interaction.
          double initialEnergy;
          /// \brief Pointer to the nucleus
          Nucleus *theNucleus;
          /// \brief Pointer to the boost vector
          ThreeVector const &boostVector;

          /// \brief True if we should use local energy
          const bool shouldUseLocalEnergy;

          /** \brief Scale the momenta of the modified and created particles.
           *
           * Set the momenta of the modified and created particles to alpha times
           * their original momenta (stored in particleMomenta). You must call
           * init() before using this method.
           *
           * \param alpha scale factor
           */
          void scaleParticleMomenta(const double alpha) const;

      };

      /// \brief RootFunctor-derived object for enforcing energy conservation in delta production
      class ViolationEEnergyFunctor : public RootFunctor {
        public:
          /** \brief Prepare for calling the () operator and setParticleEnergy
           *
           * The constructor sets the private class members.
           */
          ViolationEEnergyFunctor(Nucleus * const nucleus, Particle * const aParticle, const double totalEnergyBeforeInteraction, const bool localE);
          virtual ~ViolationEEnergyFunctor() {}

          /** \brief Compute the energy-conservation violation.
           *
           * \param x scale factor for the particle energy
           * \return the energy-conservation violation
           */
          double operator()(const double x) const;

          /// \brief Clean up after root finding
          void cleanUp(const bool success) const;

          /** \brief Set the energy of the particle.
           *
           * \param energy
           */
          void setParticleEnergy(const double energy) const;

        private:
          /// \brief Total energy before the interaction.
          double initialEnergy;
          /// \brief Pointer to the nucleus.
          Nucleus *theNucleus;
          /// \brief The final-state particle.
          Particle *theParticle;
          /// \brief The initial energy of the particle.
          double theEnergy;
          /// \brief The initial momentum of the particle.
          ThreeVector theMomentum;
          /** \brief Threshold for the energy of the particle
           *
           * The particle (a delta) cannot have less than this energy.
           */
          double energyThreshold;
          /// \brief Whether we should use local energy
          const bool shouldUseLocalEnergy;
      };

      RootFunctor *violationEFunctor;

    protected:
      /** \brief Enforce energy conservation.
       *
       * Final states generated by the channels might violate energy conservation
       * because of different reasons (energy-dependent potentials, local
       * energy...). This conservation law must therefore be enforced by hand. We
       * do so by rescaling the momenta of the final-state particles in the CM
       * frame. If this turns out to be impossible, this method returns false.
       *
       * \return true if the algorithm succeeded
       */
      bool enforceEnergyConservation(FinalState * const fs);

      ParticleList modified, created, modifiedAndCreated, Destroyed, ModifiedAndDestroyed;

      INCL_DECLARE_ALLOCATION_POOL(InteractionAvatar)
  };

}

#endif /* G4INCLINTERACTIONAVATAR_HH_ */
