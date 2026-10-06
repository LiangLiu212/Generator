/** \file G4INCLNDFParis.hh
 * \brief NDF* class for the deuteron density according to the Paris
 *        potential.
 *
 * \date 16 July 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLNDFPARIS_HH_
#define G4INCLNDFPARIS_HH_

#include "G4INCLIFunction1D.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLDeuteronDensity.hh"

namespace G4INCL {

  namespace NuclearDensityFunctions {

    class ParisR : public IFunction1D {
      public:
        ParisR() : // there are no free parameters in the Paris potential
          // We let the function go up to 17.42 fm
          IFunction1D(0., 17.42)
      {};

        inline double operator()(const double r) const {
          return DeuteronDensity::densityR(r);
        }
    };

    class ParisP : public IFunction1D {
      public:
        ParisP() : // there are no free parameters in the Paris potential
          // We let the function go up to 630 MeV/c
          IFunction1D(0., 630.)
      {};

        inline double operator()(const double p) const {
          return DeuteronDensity::densityP(p);
        }
    };

  }

}

#endif // G4INCLNDFPARIS_HH_

