/** \file G4INCLRootFinder.hh
 * \brief Static root-finder algorithm.
 *
 * Provides a (nearly) stateless root-finder algorithm.
 *
 * \date 2nd March 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLROOTFINDER_HH_
#define G4INCLROOTFINDER_HH_

#include <utility>
#include "G4INCLIFunction1D.hh"

namespace G4INCL {

  class RootFunctor : public IFunction1D {
  public:
    virtual void cleanUp(const bool success) const = 0;
    virtual ~RootFunctor() {};
  protected:
    RootFunctor(const double x0, const double x1) :
      IFunction1D(x0, x1)
    {};
  };

  namespace RootFinder {

    class Solution {
      public:
        Solution() :
          success(false),
          x(0.),
          y(0.)
      {}
        Solution( const double x0, const double y0) :
          success(true),
          x(x0),
          y(y0)
      {}
        ~Solution() {}

        bool success;
        double x;
        double y;
    };

    /** \brief Numerically solve a one-dimensional equation.
     *
     * Numerically solves the equation f(x)==0. This implementation uses the
     * false-position method.
     *
     * If a root is found, it can be retrieved using the getSolution() method,
     *
     * \param f pointer to a RootFunctor
     * \param x0 initial value of the function argument
     * \return a Solution object describing the root, if it was found
     */
    Solution solve(RootFunctor const * const f, const double x0);

  }
}
#endif /* G4INCLROOTFINDER_HH_ */
