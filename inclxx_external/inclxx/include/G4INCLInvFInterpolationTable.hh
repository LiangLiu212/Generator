/** \file G4INCLInvFInterpolationTable.hh
 * \brief Simple interpolation table for the inverse of a IFunction1D functor
 *
 * \date 16 July 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLINVFINTERPOLATIONTABLE_HH_
#define G4INCLINVFINTERPOLATIONTABLE_HH_

#include "G4INCLInterpolationTable.hh"
#include <algorithm>
#include <functional>
#include <sstream>

namespace G4INCL {

  /// \brief Class for interpolating the inverse of a 1-dimensional function
  class InvFInterpolationTable : public InterpolationTable {
    public:
      InvFInterpolationTable(IFunction1D const &f, const unsigned int nNodes=30);
      virtual ~InvFInterpolationTable() {}
  };

}

#endif // G4INCLINVFINTERPOLATIONTABLE_HH_
