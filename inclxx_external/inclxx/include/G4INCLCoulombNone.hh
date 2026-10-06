/** \file G4INCLCoulombNone.hh
 * \brief Placeholder class for no Coulomb distortion.
 *
 * \date 14 February 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLCOULOMBNONE_HH_
#define G4INCLCOULOMBNONE_HH_

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLICoulomb.hh"
#include <utility>

namespace G4INCL {

  class CoulombNone : public ICoulomb {

    public:
    CoulombNone() {}
    virtual ~CoulombNone() {}

    /** \brief Position the particle on the surface of the nucleus.
     *
     * This method does not perform any distortion.
     *
     * \param p incoming particle
     * \param n distorting nucleus
     **/
    ParticleEntryAvatar *bringToSurface(Particle * const p, Nucleus * const n) const;

    /** \brief Position the particle on the surface of the nucleus.
     * ONLY FOR ANTIDEUTERON !!!
     * This method does not perform any distortion.
     *
     * \param p incoming particle
     * \param n distorting nucleus
     **/
    ParticleEntryAvatar *bringToSurfaceAbar(Particle * const p, Nucleus * const n) const;

    /** \brief Position the cluster on the surface of the nucleus.
     *
     * This method does not perform any distortion.
     *
     * \param c incoming cluster
     * \param n distorting nucleus
     **/
    IAvatarList bringToSurface(Cluster * const c, Nucleus * const n) const;

    /** \brief Modify the momenta of the outgoing particles.
     *
     * This method does not perform any distortion.
     */
    void distortOut(ParticleList const & /* pL */, Nucleus const * const /* n */) const {}

    /** \brief Return the maximum impact parameter for Coulomb-distorted
     *         trajectories. **/
    double maxImpactParameter(ParticleSpecies const &p, const double /*kinE*/, Nucleus const *
        const n) const {
      if(p.theType == Composite || p.theType == antiComposite)
        return 2.*ParticleTable::getLargestNuclearRadius(p.theA, p.theZ)
          + n->getUniverseRadius();
      else
        return n->getUniverseRadius();
    }

  };
}

#endif /* G4INCLCOULOMBNONE_HH_ */
