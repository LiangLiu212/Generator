/** \file G4INCLDeJongSpin.hh
 * \brief Simple class implementing De Jong's spin model for nucleus-nucleus
 *        collisions
 *
 * Reference: De Jong, Ignatyuk and Schmidt, Nucl. Phys. A613 (1997) 435-444.
 *
 * \date 2 April 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLDEJONGSPIN_HH_
#define G4INCLDEJONGSPIN_HH_

namespace G4INCL {
  namespace DeJongSpin {
    ThreeVector shoot(const int Ap, const int Af);
  }
}

#endif // G4INCLDEJONGSPIN_HH_
