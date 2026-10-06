/** \file G4INCLIFunction1D.hh
 * \brief Functor for 1-dimensional mathematical functions
 *
 * \date 16 July 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLIFUNCTION1D_HH_
#define G4INCLIFUNCTION1D_HH_ 1

#include <vector>

namespace G4INCL {

  // Forward declaration
  class InterpolationTable;

  /**
   * 1D function interface
   */
  class IFunction1D {
    public:
      IFunction1D() :
        xMin(0.),
        xMax(0.)
    {};
      IFunction1D(const double x0, const double x1) :
        xMin(x0),
        xMax(x1)
    {};

      virtual ~IFunction1D() {};

      /// \brief Return the minimum allowed value of the independent variable
      virtual inline double getXMinimum() const { return xMin; }

      /// \brief Return the maximum allowed value of the independent variable
      virtual inline double getXMaximum() const { return xMax; }

      /// \brief Compute the value of the function
      virtual double operator()(const double x) const = 0;

      /** \brief Integrate the function between two values
       *
       * \param x0 lower integration bound
       * \param x1 upper integration bound
       * \param step largest integration step size; if <0, 45 steps will be used
       * \return \f$\int_{x_0}^{x_1} f(x) dx\f$
       */
      virtual double integrate(const double x0, const double x1, const double step=-1.) const;

      /// \brief Return a pointer to the (numerical) primitive to this function
      IFunction1D *primitive() const;

      /// \brief Typedef to simplify the syntax of inverseCDFTable
      typedef double (* const ManipulatorFunc)(const double);

      /** \brief Return a pointer to the inverse of the CDF of this function
       *
       * The function parameter fWrap is wrapped around the return value of
       * operator(). If fWrap=NULL (default), fWrap=identity.
       */
      InterpolationTable *inverseCDFTable(ManipulatorFunc fWrap=0, const int nNodes=60) const;

    protected:
      /// \brief Minimum value of the independent variable
      double xMin;
      /// \brief Maximum value of the independent variable
      double xMax;

    private:
      /// \brief Coefficients for numerical integration
      static const double integrationCoefficients[];
  };

}

#endif // G4INCLIFUNCTION1D_HH_
