/** \file G4INCLProjectileRemnant.hh
 * \brief Class for constructing a projectile-like remnant.
 *
 * \date 20 March 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLPROJECTILEREMNANT_HH_
#define G4INCLPROJECTILEREMNANT_HH_

#include "G4INCLCluster.hh"
#include "G4INCLRandom.hh"
#include "G4INCLAllocationPool.hh"
#include <vector>
#include <map>
#include <numeric>
#include <functional>

namespace G4INCL {

  class ProjectileRemnant : public Cluster {
    public:

    // typedefs for the calculation of the projectile excitation energy
    typedef std::vector<double> EnergyLevels;
    typedef std::map<long, double> EnergyLevelMap;

    ProjectileRemnant(ParticleSpecies const &species, const double kineticEnergy)
      : Cluster(species.theZ, species.theA, species.theS) {

      // Use the table mass
      setTableMass();

      // Set the kinematics
      const double projectileMass = getMass();
      const double energy = kineticEnergy + projectileMass;
      const double momentumZ = std::sqrt(energy*energy - projectileMass*projectileMass);

      // Initialise the particles
      initializeParticles();
      internalBoostToCM();
      putParticlesOffShell();

      // Store the energy levels of the ProjectileRemnant (used to compute its
      // excitation energy)
      storeEnergyLevels();

      // Boost the whole thing
      const ThreeVector aBoostVector = ThreeVector(0.0, 0.0, momentumZ / energy);
      boost(-aBoostVector);

      // Freeze the internal motion of the particles
      freezeInternalMotion();

      // Set as projectile spectator
      makeProjectileSpectator();
    }

    ~ProjectileRemnant() {
      deleteStoredComponents();
      // The ProjectileRemnant owns its particles
      deleteParticles();
      clearEnergyLevels();
    }

    /// \brief Reset the projectile remnant to the state at the beginning of the cascade
    void reset();

    /** \brief Remove a nucleon from the projectile remnant
     *
     * \param p particle to be removed
     * \param theProjectileCorrection correction to be given to the projectile total energy
     */
    void removeParticle(Particle * const p, const double theProjectileCorrection);

    /** \brief Add back dynamical spectators to the projectile remnant
     *
     * Try to add the dynamical spectators back to the projectile remnant.
     * Refuse to do so if this leads to a negative projectile excitation
     * energy.
     *
     * Return a list of rejected dynamical spectators.
     */
    ParticleList addDynamicalSpectators(ParticleList pL);

    /** \brief Add back dynamical spectators to the projectile remnant
     *
     * Try as hard as possible to add back all the dynamical spectators.  Don't
     * add spectators that lead to negative excitation energies. Start by
     * adding all of them, and repeatedly remove the most troublesome one until
     * the excitation energy becomes non-negative.
     *
     * Return a list of rejected dynamical spectators.
     */
    ParticleList addMostDynamicalSpectators(ParticleList pL);

    /** \brief Add back all dynamical spectators to the projectile remnant
     *
     * Return a list of rejected dynamical spectators.
     */
    ParticleList addAllDynamicalSpectators(ParticleList const &pL);

    /// \brief Clear the stored projectile components and delete the particles
    void deleteStoredComponents() {
      for(std::map<long,Particle*>::const_iterator p=storedComponents.begin(), e=storedComponents.end(); p!=e; ++p)
        delete p->second;
      clearStoredComponents();
    }

    /// \brief Clear the stored projectile components
    void clearStoredComponents() {
      storedComponents.clear();
    }

    /// \brief Clear the stored energy levels
    void clearEnergyLevels() {
      theInitialEnergyLevels.clear();
      theGroundStateEnergies.clear();
    }

    /** \brief Compute the excitation energy when a nucleon is removed
     *
     * Compute the excitation energy of the projectile-like remnant as the
     * difference between the initial and the present configuration. This
     * follows the algorithm proposed by A. Boudard in INCL4.2-HI, as
     * implemented in Geant4.
     *
     * \return the excitation energy
     */
    double computeExcitationEnergyExcept(const long exceptID) const;

    /** \brief Compute the excitation energy if some nucleons are put back
     *
     * \return the excitation energy
     */
    double computeExcitationEnergyWith(const ParticleList &pL) const;

    /// \brief Store the projectile components
    void storeComponents() {
      for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
        // Store the particles (needed for forced CN)
        storedComponents[(*p)->getID()]=new Particle(**p);
      }
    }

    /// \brief Get the number of the stored components
    int getNumberStoredComponents() const {
      return (int)storedComponents.size();
    }

    /// \brief Store the energy levels
    void storeEnergyLevels() {
      EnergyLevels energies;

      for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
        const double theCMEnergy = (*p)->getEnergy();
        // Store the CM energy in the EnergyLevels map
        theInitialEnergyLevels[(*p)->getID()] = theCMEnergy;
        energies.push_back(theCMEnergy);
      }

      std::sort(energies.begin(), energies.end());
      assert(energies.size()==(unsigned int)theA || energies.size()== -theA);
      theGroundStateEnergies.resize(energies.size());
      // Compute the partial sums of the CM energies -- they are our reference
      // ground-state energies for any number of nucleons
      std::partial_sum(energies.begin(), energies.end(), theGroundStateEnergies.begin());
    }

    EnergyLevels const &getGroundStateEnergies() const {
      return theGroundStateEnergies;
    }

    private:

    /** \brief Compute the excitation energy for a given configuration
     *
     * The function that does the real job of calculating the excitation energy
     * for a given configuration of energy levels.
     *
     * \param levels a configuration of energy levels
     * \return the excitation energy
     */
    double computeExcitationEnergy(const EnergyLevels &levels) const;

    EnergyLevels getPresentEnergyLevelsExcept(const long exceptID) const;

    EnergyLevels getPresentEnergyLevelsWith(const ParticleList &pL) const;

    /// \brief Shuffle the list of stored projectile components
    ParticleList shuffleStoredComponents() {
      ParticleList pL = getStoredComponents();
      std::shuffle(pL.begin(), pL.end(), Random::getAdapter());
      return pL;
    }

    ParticleList getStoredComponents() const {
      ParticleList pL;
      for(std::map<long,Particle*>::const_iterator p=storedComponents.begin(), e=storedComponents.end(); p!=e; ++p)
        pL.push_back(p->second);
      return pL;
    }

    /// \brief Return the stored momentum of a given projectile component
    ThreeVector const &getStoredMomentum(Particle const * const p) const {
      std::map<long,Particle*>::const_iterator i = storedComponents.find(p->getID());
      if(i==storedComponents.end()) {
        INCL_ERROR("Couldn't find particle " << p->getID() << " in the list of projectile components" << '\n');
        return p->getMomentum();
      } else {
        return i->second->getMomentum();
      }
    }

    /** \brief Add back a nucleon to the projectile remnant
     *
     * Try to add a dynamical spectator back to the projectile remnant. Refuse
     * to do so if this leads to a negative projectile excitation energy.
     * Return true on success, false on failure.
     */
    bool addDynamicalSpectator(Particle * const p);

    /// \brief Return the stored energy of a given projectile component
    /*    double getStoredEnergy(Particle const * const p) {
          std::map<long,Particle*>::const_iterator i = initialProjectileComponents.find(p->getID());
          if(i==initialProjectileComponents.end()) {
          INCL_ERROR("Couldn't find particle " << p->getID() << " in the list of projectile components" << '\n');
          return 0.;
          } else {
          return i->second->getEnergy();
          }
          }*/

    /** \brief Stored projectile components
     *
     * These particles are owned by the ProjectileRemnant.
     */
    std::map<long, Particle*> storedComponents;

    /// \brief Initial energy levels of the projectile
    EnergyLevelMap theInitialEnergyLevels;

    /// \brief Ground-state energies for any number of nucleons
    EnergyLevels theGroundStateEnergies;

    INCL_DECLARE_ALLOCATION_POOL(ProjectileRemnant)
  };
}

#endif // G4INCLPROJECTILEREMNANT_HH_

