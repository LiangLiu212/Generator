/** \file G4INCLDeuteronDensity.hh
 * \brief Deuteron density in r and p according to the Paris potential.
 *
 * \date 6 March 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLDEUTERONDENSITY_HH_
#define G4INCLDEUTERONDENSITY_HH_ 1

#include "G4INCLIFunction1D.hh"

namespace G4INCL {

  /// \brief Static class for the deuteron density
  namespace DeuteronDensity {

    /** \brief PDF for a nucleon in r space
     *
     * The distribution is normalised to 1.
     *
     * \param r distance from the deuteron centre [fm]
     * \return 4 * pi * r^2 * |psi|^2
     */
    double densityR(const double r);

    /** \brief First derivative of the r-space density function
     *
     * \param r distance from the deuteron centre [fm]
     * \return d|psi|^2/dr
     */
    double derivDensityR(const double r);

    /** \brief PDF for a nucleon in p space
     *
     * The distribution is normalised to 1.
     *
     * \param p nucleon momentum [MeV/c]
     * \return 4 * pi * p^2 * d|psi|^2/dp
     */
    double densityP(const double p);

    double wavefunctionR(const int l, const double r);
    double wavefunctionP(const int l, const double p);
    double derivWavefunctionR(const int l, const double r);

  }

}

#endif // G4INCLDEUTERONDENSITY_HH_
