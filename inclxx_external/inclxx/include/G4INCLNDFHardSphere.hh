/** \file G4INCLNDFHardSphere.hh
 * \brief NDF* class for the deuteron density according to the HardSphere
 *        potential.
 *
 * \date 16 July 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLNDFHARDSPHERE_HH_
#define G4INCLNDFHARDSPHERE_HH_

#include "G4INCLIFunction1D.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLDeuteronDensity.hh"

namespace G4INCL {

  namespace NuclearDensityFunctions {

    class HardSphere : public IFunction1D {
      public:
        HardSphere(const double rMax) :
          // We let the function go up to rMax. This can be the Fermi momentum
          // for a sphere in momentum space, or the surface radius for a
          // hard-sphere nucleus.
          IFunction1D(0., rMax),
          normalisation(3./std::pow(rMax,3.))
      {};

        inline double operator()(const double x) const {
          return x*x*normalisation;
        }
      protected:
        const double normalisation;
    };

  }

}

#endif // G4INCLNDFHARDSPHERE_HH_

