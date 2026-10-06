#ifndef KinematicsUtils_hh
#define KinematicsUtils_hh 1

#include "G4INCLThreeVector.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLParticleSpecies.hh"

namespace G4INCL {

  namespace KinematicsUtils {

    double fiveParFit(const double a, const double b, const double c, const double d, const double e, const double x);
    double compute_xs(const std::vector<double> coefficients, const double pLab);

    void transformToLocalEnergyFrame(Nucleus const * const n, Particle * const p);
    double getLocalEnergy(Nucleus const * const n, Particle * const p);

    ThreeVector makeBoostVector(Particle const * const p1, Particle const * const p2);
    double totalEnergyInCM(Particle const * const p1, Particle const * const p2);
    double squareTotalEnergyInCM(Particle const * const p1, Particle const * const p2);

    /** \brief gives the momentum in the CM frame of two particles.
     *
     * The formula is the following:
     * \f[ p_{CM}^2 = \frac{z^2 - m_1^2 m_2^2}{2 z + m_1^2 + m_2^2} \f]
     * where \f$z\f$ is the scalar product of the momentum four-vectors:
     * \f[ z = E_1 E_2 - \vec{p}_1\cdot\vec{p}_2 \f]
     *
     * \param p1 pointer to particle 1
     * \param p2 pointer to particle 2
     * \return the absolute value of the momentum of any of the two particles in
     * the CM frame, in MeV/c.
     */
    double momentumInCM(Particle const * const p1, Particle const * const p2);

    double momentumInCM(const double E, const double M1, const double M2);

    /** \brief gives the momentum in the lab frame of two particles.
     *
     * Assumes particle 1 carries all the momentum and particle 2 is at rest.
     *
     * The formula is the following:
     * \f[ p_{lab}^2 = \frac{s^2 - 2 s (m_1^2 + m_2^2) + {(m_1^2 - m_2^2)}^2}{4 m_2^2} \f]
     *
     * \param p1 pointer to particle 1
     * \param p2 pointer to particle 2
     * \return the absolute value of the momentum of particle 1 in the lab frame,
     * in MeV/c
     */
    double momentumInLab(Particle const * const p1, Particle const * const p2);
    double momentumInLab(const double s, const double m1, const double m2);
    double sumTotalEnergies(const ParticleList &);
    ThreeVector sumMomenta(const ParticleList &);
    double energy(const ThreeVector &p, const double m);
    double invariantMass(const double E, const ThreeVector & p);
    double squareInvariantMass(const double E, const ThreeVector & p);
    double gammaFromKineticEnergy(const ParticleSpecies &p, const double EKin);
  }
}

#endif
