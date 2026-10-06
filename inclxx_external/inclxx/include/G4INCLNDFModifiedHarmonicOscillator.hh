/** \file G4INCLNDFModifiedHarmonicOscillator.hh
 * \brief Class for modified harmonic oscillator density
 *
 * \date 16 July 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLNDFMODIFIEDHARMONICOSCILLATOR_HH_
#define G4INCLNDFMODIFIEDHARMONICOSCILLATOR_HH_

#include "G4INCLIFunction1D.hh"
#include <cmath>
#include <algorithm>

namespace G4INCL {

  namespace NuclearDensityFunctions {

    class ModifiedHarmonicOscillatorRP : public IFunction1D {
      public:
        ModifiedHarmonicOscillatorRP(double radiusParameter, double maximumRadius, double diffusenessParameter) :
          IFunction1D(0., maximumRadius),
          theRadiusParameter(radiusParameter),
          theDiffusenessParameter(diffusenessParameter)
      {}

        inline double operator()(const double r) const {
          const double arg = std::pow((r/theDiffusenessParameter),2);
          return std::max(0., -2.0* r*r *arg * (theRadiusParameter - 1.0 - theRadiusParameter*arg)*std::exp(-arg));
        }

        inline double getRadiusParameter() { return theRadiusParameter; };
        inline double getDiffusenessParameter() { return theDiffusenessParameter; };

        inline void setRadiusParameter(double r) { theRadiusParameter = r; };
        inline void setDiffusenessParameter(double a) { theDiffusenessParameter = a; };

      protected:
        double theRadiusParameter, theDiffusenessParameter;
    };

    class ModifiedHarmonicOscillator : public IFunction1D {
      public:
        ModifiedHarmonicOscillator(double radiusParameter, double maximumRadius, double diffusenessParameter) :
          IFunction1D(0., maximumRadius),
          theRadiusParameter(radiusParameter),
          theDiffusenessParameter(diffusenessParameter),
          normalisation(2./((theDiffusenessParameter+theRadiusParameter)*std::pow(theDiffusenessParameter,2.)))
      {}

        inline double operator()(const double r) const {
          const double arg = std::pow((r/theDiffusenessParameter),2);
          return r*r*( 1. + theRadiusParameter * arg) * std::exp(-arg);
        }

        inline double getRadiusParameter() { return theRadiusParameter; };
        inline double getDiffusenessParameter() { return theDiffusenessParameter; };

        inline void setRadiusParameter(double r) { theRadiusParameter = r; };
        inline void setDiffusenessParameter(double a) { theDiffusenessParameter = a; };

      protected:
        double theRadiusParameter, theDiffusenessParameter;
        const double normalisation;
    };

  }

}

#endif // G4INCLNDFMODIFIEDHARMONICOSCILLATOR_HH_

