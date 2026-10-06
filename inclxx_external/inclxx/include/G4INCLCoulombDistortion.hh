/** \file G4INCLCoulombDistortion.hh
 * \brief Static class for selecting Coulomb distortion.
 *
 * \date 14 February 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLCOULOMBDISTORTION_HH_
#define G4INCLCOULOMBDISTORTION_HH_

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLICoulomb.hh"
#include "G4INCLConfig.hh"

namespace G4INCL {

  /**
   * Coulomb distortion
   */
  namespace CoulombDistortion {
    /** \brief Modify the momentum of an incoming particle and position it on
     *         the surface of the nucleus.
     *
     * This method places Particle p on the surface of Nucleus n and modifies
     * the direction of its momentum to be tangent to the Coulomb trajectory in
     * that point.
     *
     * The input particle has to be prepared with its asymptotic momentum. Its
     * position is used only for the purpose of computing the asymptotic impact
     * parameter; in other words, this method only uses the components of the
     * position that are perpendicular to the momentum. The remaining component
     * is not used, and can be set to any value.
     *
     * This method returns a ParticleEntry avatar for the projectile.
     *
     * \param p incoming particle
     * \param n distorting nucleus
     * \return the ParticleEntryAvatar for the projectile particle
     **/
    ParticleEntryAvatar *bringToSurface(Particle *p, Nucleus * const n);

    ParticleEntryAvatar *bringToSurfaceAbar(Particle *p, Nucleus * const n);

    /** \brief Modify the momentum of an incoming cluster and position it on
     *         the surface of the target.
     *
     * Same as the Particle-based bringToSurface method, but for incoming heavy
     * ions.
     *
     * This method returns a list of ParticleEntry avatars for the participant
     * nucleons
     *
     * \param c incoming heavy ion
     * \param n distorting nucleus
     * \return a list of ParticleEntryAvatars
     **/
    IAvatarList bringToSurface(Cluster * const c, Nucleus * const n);

    /** \brief Modify the momentum of an outgoing particle. */
    void distortOut(ParticleList const &pL, Nucleus const * const n);

    /** \brief Return the maximum impact parameter for Coulomb-distorted
     *         trajectories. **/
    double maxImpactParameter(ParticleSpecies const &p, const double kinE, Nucleus const * const n);

    /** \brief Return the maximum impact parameter for Coulomb-distorted
     *         trajectories. **/
    double maxImpactParameter(Particle const * const p, Nucleus const * const n);

    /** \brief Set the Coulomb-distortion algorithm. */
    void setCoulomb(ICoulomb * const coulomb);

    /** \brief Delete the Coulomb-distortion object. */
    void deleteCoulomb();

    /// \brief Initialize the Coulomb-distortion algorithm
    void initialize(Config const * const theConfig);

  }
}

#endif /* G4INCLCOULOMBDISTORTION_HH_ */
