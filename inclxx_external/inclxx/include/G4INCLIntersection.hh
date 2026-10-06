/** \file G4INCLIntersection.hh
 * \brief Simple class for computing intersections between a straight line and a sphere.
 *
 * \date 12 December 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLINTERSECTION_HH
#define G4INCLINTERSECTION_HH 1

#include "G4INCLThreeVector.hh"
#include <utility>

namespace G4INCL {

  /** \brief Intersection-point structure
   *
   * The structure contains the time and position of the intersection point
   * of a trajectory with a surface, if it exists.
   */
  struct Intersection {
    Intersection(const bool e, const double t, const ThreeVector &p) :
      exists(e), time(t), position(p) {}
    bool exists;
    double time;
    ThreeVector position;
  };

  namespace IntersectionFactory {

    /** \brief Compute the first intersection of a straight particle
     *         trajectory with a sphere.
     *
     * \param x0 the starting position of the trajectory
     * \param p the trajectory direction
     * \param r the radius of the sphere (centred in the origin)
     * \return an Intersection. The bool is true if an intersection exists,
     *         in which case its position is stored in the ThreeVector and
     *         its time in the double.
     */
    Intersection getEarlierTrajectoryIntersection(const ThreeVector &x0, const ThreeVector &p, const double r);
    /** \brief Compute the second intersection of a straight particle
     *         trajectory with a sphere.
     *
     * \param x0 the starting position of the trajectory
     * \param p the trajectory direction
     * \param r the radius of the sphere (centred in the origin)
     * \return an Intersection. The bool is true if an intersection exists,
     *         in which case its position is stored in the ThreeVector and
     *         its time in the double.
     */
    Intersection getLaterTrajectoryIntersection(const ThreeVector &x0, const ThreeVector &p, const double r);
    /** \brief Compute both intersections of a straight particle
     *         trajectory with a sphere.
     *
     * \param x0 the starting position of the trajectory
     * \param p the trajectory direction
     * \param r the radius of the sphere (centred in the origin)
     * \return an Intersection. The bool is true if an intersection exists,
     *         in which case its position is stored in the ThreeVector and
     *         its time in the double.
     */
    std::pair<Intersection,Intersection> getTrajectoryIntersections(const ThreeVector &x0, const ThreeVector &p, const double r);

    namespace {

      Intersection getTrajectoryIntersection(const ThreeVector &x0, const ThreeVector &v, const double r, const bool earliest) {
        const double scalarVelocity = v.mag();
        ThreeVector velocityUnitVector = v / scalarVelocity;

        ThreeVector positionTransverse = x0 - velocityUnitVector * x0.dot(velocityUnitVector);
        const double impactParameter = positionTransverse.mag();

        const double r2 = r*r;
        double distanceZ2 = r2 - impactParameter * impactParameter;
        if(distanceZ2 < 0.0)
          return Intersection(false, 0.0, ThreeVector());

        const double distanceZ = std::sqrt(distanceZ2);
        const ThreeVector position = positionTransverse + velocityUnitVector * (earliest ? -distanceZ : distanceZ);
        const double time = (position-x0).dot(velocityUnitVector)/scalarVelocity;
        return Intersection(true, time, position);
      }

    }

    inline Intersection getEarlierTrajectoryIntersection(const ThreeVector &x0, const ThreeVector &p, const double r) {
      return getTrajectoryIntersection(x0, p, r, true);
    }

    inline Intersection getLaterTrajectoryIntersection(const ThreeVector &x0, const ThreeVector &p, const double r) {
      return getTrajectoryIntersection(x0, p, r, false);
    }

    inline std::pair<Intersection,Intersection> getTrajectoryIntersections(const ThreeVector &x0, const ThreeVector &p, const double r) {
      return std::make_pair(
                            getTrajectoryIntersection(x0, p, r, true),
                            getTrajectoryIntersection(x0, p, r, false)
                           );
    }
  }

}

#endif /* G4INCLINTERSECTION_HH */
