/** \file G4INCLNDFGaussian.hh
 * \brief Class for Gaussian density
 *
 * \date 16 July 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLNDFGAUSSIAN_HH_
#define G4INCLNDFGAUSSIAN_HH_

#include "G4INCLIFunction1D.hh"
#include <cmath>

namespace G4INCL {

  namespace NuclearDensityFunctions {

    class GaussianRP : public IFunction1D {
      public:
        GaussianRP(double maximumRadius, double standardDeviation) :
          IFunction1D(0., maximumRadius),
          theStandardDeviation(standardDeviation)
      {}

        inline double operator()(const double r) const {
          const double arg = std::pow((r/theStandardDeviation),2);
          return r*r*arg*std::exp(-arg/2.0);
        }

      protected:
        double theStandardDeviation;
    };

    class Gaussian : public IFunction1D {
      public:
        Gaussian(double maximumRadius, double standardDeviation) :
          IFunction1D(0., maximumRadius),
          theStandardDeviation(standardDeviation),
          normalisation(std::sqrt(2./Math::pi)/theStandardDeviation)
      {}

        inline double operator()(const double r) const {
          const double arg = std::pow((r/theStandardDeviation),2);
          return normalisation * arg * std::exp(-arg/2.0);
        }

      protected:
        double theStandardDeviation;
        const double normalisation;
    };

  }

}

#endif // G4INCLNDFGAUSSIAN_HH_

