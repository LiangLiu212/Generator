#ifndef G4INCLCascade_hh
#define G4INCLCascade_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIPropagationModel.hh"
#include "G4INCLCascadeAction.hh"
#include "G4INCLEventInfo.hh"
#include "G4INCLGlobalInfo.hh"
#include "G4INCLLogger.hh"
#include "G4INCLConfig.hh"
#include "G4INCLRootFinder.hh"
#ifndef INCLXX_IN_GEANT4_MODE
 #include "G4INCLGeant4Compat.hh"
#endif


namespace G4INCL {
  class INCL {
    public:
      INCL(Config const * const config);

      ~INCL();

      /// \brief Dummy copy constructor to silence Coverity warning
      INCL(const INCL &rhs);

      /// \brief Dummy assignment operator to silence Coverity warning
      INCL &operator=(const INCL &rhs);

      bool prepareReaction(const ParticleSpecies &projectileSpecies, const double kineticEnergy, const int A, const int Z, const int S);
      
      bool initializeTarget(const int A, const int Z, const int S, AnnihilationType theAType);
      double initUniverseRadiusForAntiprotonAtRest(const int A, const int Z);
      inline const EventInfo &processEvent() {
        return processEvent(
            theConfig->getProjectileSpecies(),
            theConfig->getProjectileKineticEnergy(),
            theConfig->getTargetA(),
            theConfig->getTargetZ(),
            theConfig->getTargetS()
            );
      }
      const EventInfo &processEvent(
          ParticleSpecies const &projectileSpecies,
          const double kineticEnergy,
          const int targetA,
          const int targetZ,
          const int targetS
          );

      void finalizeGlobalInfo(Random::SeedVector const &initialSeeds);
      const GlobalInfo &getGlobalInfo() const { return theGlobalInfo; }
      

    private:
      IPropagationModel *propagationModel;
      int theA, theZ, theS;
      bool targetInitSuccess;
      double maxImpactParameter;
      double maxUniverseRadius;
      double maxInteractionDistance;
      double fixedImpactParameter;
      CascadeAction *cascadeAction;
      Config const * const theConfig;
      Nucleus *nucleus;
      bool forceTransparent;
      
      EventInfo theEventInfo;
      GlobalInfo theGlobalInfo;

      /// \brief Remnant size below which cascade stops
      int minRemnantSize;

      /// \brief Class to adjust remnant recoil
      class RecoilFunctor : public RootFunctor {
        public:
          /** \brief Prepare for calling the () operator and scaleParticleEnergies
           *
           * The constructor sets the private class members.
           */
          RecoilFunctor(Nucleus * const n, const EventInfo &ei) :
            RootFunctor(0., 1E6),
            nucleus(n),
            outgoingParticles(n->getStore()->getOutgoingParticles()),
            theEventInfo(ei) {
              for(ParticleIter p=outgoingParticles.begin(), e=outgoingParticles.end(); p!=e; ++p) {
                particleMomenta.push_back((*p)->getMomentum());
                particleKineticEnergies.push_back((*p)->getKineticEnergy());
              }
              ProjectileRemnant * const aPR = n->getProjectileRemnant();
              if(aPR && (aPR->getA()>0 || aPR->getA()<0)) {
                particleMomenta.push_back(aPR->getMomentum());
                particleKineticEnergies.push_back(aPR->getKineticEnergy());
                outgoingParticles.push_back(aPR);
              }
            }
          virtual ~RecoilFunctor() {}

          /** \brief Compute the energy-conservation violation.
           *
           * \param x scale factor for the particle energies
           * \return the energy-conservation violation
           */
          double operator()(const double x) const {
            scaleParticleEnergies(x);
            return nucleus->getConservationBalance(theEventInfo,true).energy;
          }

          /// \brief Clean up after root finding
          void cleanUp(const bool success) const {
            if(!success)
              scaleParticleEnergies(1.);
          }

        private:
          /// \brief Pointer to the nucleus
          Nucleus *nucleus;
          /// \brief List of final-state particles.
          ParticleList outgoingParticles;
          // \brief Reference to the EventInfo object
          EventInfo const &theEventInfo;
          /// \brief Initial momenta of the outgoing particles
          std::list<ThreeVector> particleMomenta;
          /// \brief Initial kinetic energies of the outgoing particles
          std::list<double> particleKineticEnergies;

          /** \brief Scale the kinetic energies of the outgoing particles.
           *
           * \param rescale scale factor
           */
          void scaleParticleEnergies(const double rescale) const {
            // Rescale the energies (and the momenta) of the outgoing particles.
            ThreeVector pBalance = nucleus->getIncomingMomentum();
            std::list<ThreeVector>::const_iterator iP = particleMomenta.begin();
            std::list<double>::const_iterator iE = particleKineticEnergies.begin();
            for( ParticleIter i = outgoingParticles.begin(), e = outgoingParticles.end(); i!=e; ++i, ++iP, ++iE)
            {
              const double mass = (*i)->getMass();
              const double newKineticEnergy = (*iE) * rescale;

              (*i)->setMomentum(*iP);
              (*i)->setEnergy(mass + newKineticEnergy);
              (*i)->adjustMomentumFromEnergy();

              pBalance -= (*i)->getMomentum();
            }

            nucleus->setMomentum(pBalance);
            const double remnantMass = ParticleTable::getTableMass(nucleus->getA(),nucleus->getZ(),nucleus->getS()) + nucleus->getExcitationEnergy();
            const double pRem2 = pBalance.mag2();
            const double recoilEnergy = pRem2/
              (std::sqrt(pRem2+remnantMass*remnantMass) + remnantMass);
            nucleus->setEnergy(remnantMass + recoilEnergy);
          }
      };

      /// \brief Class to adjust remnant recoil in the reaction CM system
      class RecoilCMFunctor : public RootFunctor {
        public:
          /** \brief Prepare for calling the () operator and scaleParticleEnergies
           *
           * The constructor sets the private class members.
           */
          RecoilCMFunctor(Nucleus * const n, const EventInfo &ei) :
            RootFunctor(0., 1E6),
            nucleus(n),
            theIncomingMomentum(nucleus->getIncomingMomentum()),
            outgoingParticles(n->getStore()->getOutgoingParticles()),
            theEventInfo(ei) {
              if(theIncomingMomentum.mag() == 0.){ //change the condition
                thePTBoostVector = {0., 0., 0.};
                //INCL_WARN("PTBoostVector at rest is zero" << '\n');      
              }
              else{
                thePTBoostVector = nucleus->getIncomingMomentum()/(nucleus->getInitialEnergy()); //D
                //INCL_WARN("PTBoostVector" << '\n');
              }
              for(ParticleIter p=outgoingParticles.begin(), e=outgoingParticles.end(); p!=e; ++p) {
                (*p)->boost(thePTBoostVector);
                particleCMMomenta.push_back((*p)->getMomentum());
              }
              ProjectileRemnant * const aPR = n->getProjectileRemnant();
              if(aPR && aPR->getA()!=0) {
                aPR->boost(thePTBoostVector);
                particleCMMomenta.push_back(aPR->getMomentum());
                outgoingParticles.push_back(aPR);
              }
            }
          virtual ~RecoilCMFunctor() {}

          /** \brief Compute the energy-conservation violation.
           *
           * \param x scale factor for the particle energies
           * \return the energy-conservation violation
           */
          double operator()(const double x) const {
            scaleParticleCMMomenta(x);
            return nucleus->getConservationBalance(theEventInfo,true).energy;
          }

          /// \brief Clean up after root finding
          void cleanUp(const bool success) const {
            if(!success)
              scaleParticleCMMomenta(1.);
          }

        private:
          /// \brief Pointer to the nucleus
          Nucleus *nucleus;
          /// \brief Projectile-target CM boost vector
          ThreeVector thePTBoostVector;
          /// \brief Incoming momentum
          ThreeVector theIncomingMomentum;
          /// \brief List of final-state particles.
          ParticleList outgoingParticles;
          /// \brief Reference to the EventInfo object
          EventInfo const &theEventInfo;
          /// \brief Initial CM momenta of the outgoing particles
          std::list<ThreeVector> particleCMMomenta;

          /** \brief Scale the kinetic energies of the outgoing particles.
           *
           * \param rescale scale factor
           */
          void scaleParticleCMMomenta(const double rescale) const {
            // Rescale the CM momenta of the outgoing particles.
            ThreeVector remnantMomentum = theIncomingMomentum;
            std::list<ThreeVector>::const_iterator iP = particleCMMomenta.begin();
            for( ParticleIter i = outgoingParticles.begin(), e = outgoingParticles.end(); i!=e; ++i, ++iP)
            {
              (*i)->setMomentum(*iP * rescale);
              (*i)->adjustEnergyFromMomentum();
              (*i)->boost(-thePTBoostVector);

              remnantMomentum -= (*i)->getMomentum();
            }

            nucleus->setMomentum(remnantMomentum);
            const double remnantMass = ParticleTable::getTableMass(nucleus->getA(),nucleus->getZ(),nucleus->getS()) + nucleus->getExcitationEnergy();
            const double pRem2 = remnantMomentum.mag2();
            const double recoilEnergy = pRem2/
              (std::sqrt(pRem2+remnantMass*remnantMass) + remnantMass);
            nucleus->setEnergy(remnantMass + recoilEnergy);
          }
      };

      /** \brief Rescale the energies of the outgoing particles.
       *
       * Allow for the remnant recoil energy by rescaling the energy (and
       * momenta) of the outgoing particles.
       */
      void rescaleOutgoingForRecoil();

#ifndef INCLXX_IN_GEANT4_MODE
      /** \brief Run global conservation checks
       *
       * Check that energy and momentum are correctly conserved. If not, issue
       * a warning.
       *
       * Also feeds the balance variables in theEventInfo.
       *
       * \param afterRecoil whether to take into account nuclear recoil
       */
      void globalConservationChecks(bool afterRecoil);
#endif

      /** \brief Stopping criterion for the cascade
       *
       * Returns true if the cascade should continue, and false if any of the
       * stopping criteria is satisfied.
       */
      bool continueCascade();

      /** \brief Make a projectile pre-fragment out of geometrical spectators
       *
       * The projectile pre-fragment is assigned an excitation energy given
       * by \f$E_\mathrm{sp}-E_\mathrm{i,A}\f$, where \f$E_\mathrm{sp}\f$ is the
       * sum of the energies of the spectator particles, and \f$E_\mathrm{i,A}\f$
       * is the sum of the smallest \f$A\f$ particle energies initially present
       * in the projectile, \f$A\f$ being the mass of the projectile
       * pre-fragment. This is equivalent to assuming that the excitation
       * energy is given by the sum of the transitions of all excited
       * projectile components to the "holes" left by the participants.
       *
       * This method can modify the outgoing list and adds a projectile
       * pre-fragment.
       *
       * \return the number of dynamical spectators that were merged back in
       *         the projectile
       */
      int makeProjectileRemnant();

      /** \brief Make a compound nucleus
       *
       * Selects the projectile components that can actually enter their
       * potential and puts them into the target nucleus. If the CN excitation
       * energy turns out to be negative, the event is considered a
       * transparent. This method modifies theEventInfo and theGlobalInfo.
       */
      void makeCompoundNucleus();

      /// \brief Initialise the cascade
      bool preCascade(ParticleSpecies const &projectileSpecies, const double kineticEnergy);

      /// \brief The actual cascade loop
      void cascade();

      /// \brief Finalise the cascade and clean up
      void postCascade(const ParticleSpecies &projectileSpecies, const double kineticEnergy);

      /** \brief Initialise the maximum interaction distance.
       *
       * Used in forced CN events.
       */
      void initMaxInteractionDistance(ParticleSpecies const &p, const double kineticEnergy);

      /** \brief Initialize the universe radius.
       *
       * Used for determining the energy-dependent size of the volume particles
       * live in.
       */
      void initUniverseRadius(ParticleSpecies const &p, const double kineticEnergy, const int A, const int Z);

      /// \brief Update global counters and other members of theGlobalInfo object
      void updateGlobalInfo();

      /// \brief Fill probabilities and particle types from datafile and return probability sum for normalization
      double read_file(std::string filename, std::vector<double>& probabilities, 
                         std::vector<std::vector<G4String>>& particle_types);

      /// \brief This function will tell you the FS line number from the datafile based on your random value
      int findStringNumber(double rdm, std::vector<double> yields);

      /// \brief Initialise the "cascade" for pbar on H1
      void preCascade_pbarH1(ParticleSpecies const &projectileSpecies, const double kineticEnergy);

      /// \brief Initialise the "cascade" for pbar on H2
      void preCascade_pbarH2(ParticleSpecies const &projectileSpecies, const double kineticEnergy);

      /// \brief Initialise the "cascade" for nbar on H1
      void preCascade_nbarH1(ParticleSpecies const &projectileSpecies, const double kineticEnergy);

      /// \brief Initialise the "cascade" for nbar on H2
      void preCascade_nbarH2(ParticleSpecies const &projectileSpecies, const double kineticEnergy);

      /// \brief Finalise the "cascade" and clean up for pbar/nbar on H1
      void postCascade_pbarH1(ParticleList const &outgoingParticles);

      /// \brief Finalise the "cascade" and clean up for pbar/nbar on H2
      void postCascade_pbarH2(ParticleList const &outgoingParticles, ParticleList const &H2Particles);
  };
}

#endif
