/** \file G4INCLNDFWoodsSaxon.hh
 * \brief Class for Woods-Saxon density
 *
 * \date 16 July 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLNDFWOODSSAXON_HH_
#define G4INCLNDFWOODSSAXON_HH_

#include "G4INCLIFunction1D.hh"
#include <cmath>

namespace G4INCL {

  namespace NuclearDensityFunctions {

    class WoodsSaxonRP : public IFunction1D {
      public:
        WoodsSaxonRP(double radiusParameter, double maximumRadius, double diffusenessParameter) :
          IFunction1D(0., maximumRadius),
          theRadiusParameter(radiusParameter),
          theDiffusenessParameter(diffusenessParameter)
      {}

        inline double operator()(const double r) const {
          double wsax = std::pow(r,3)
            *std::exp((r - theRadiusParameter)/theDiffusenessParameter)
            /std::pow((1.0 + std::exp((r - theRadiusParameter)/theDiffusenessParameter)),2);
          return wsax/theDiffusenessParameter;
        }

        inline double getRadiusParameter() { return theRadiusParameter; };
        inline double getDiffusenessParameter() { return theDiffusenessParameter; };

        inline void setRadiusParameter(double r) { theRadiusParameter = r; };
        inline void setDiffusenessParameter(double a) { theDiffusenessParameter = a; };

      protected:
        double theRadiusParameter, theDiffusenessParameter;
    };

    class WoodsSaxon : public IFunction1D {
      public:
        WoodsSaxon(double radiusParameter, double maximumRadius, double diffusenessParameter) :
          IFunction1D(0., maximumRadius),
          theRadiusParameter(radiusParameter),
          theDiffusenessParameter(diffusenessParameter)
      {}

        inline double operator()(const double r) const {
          return r * r / (1.0 + std::exp((r - theRadiusParameter)/theDiffusenessParameter));
        }

        inline double getRadiusParameter() { return theRadiusParameter; };
        inline double getDiffusenessParameter() { return theDiffusenessParameter; };

        inline void setRadiusParameter(double r) { theRadiusParameter = r; };
        inline void setDiffusenessParameter(double a) { theDiffusenessParameter = a; };

      protected:
        double theRadiusParameter, theDiffusenessParameter;
    };

  }

}

#endif // G4INCLNDFWOODSSAXON_HH_
