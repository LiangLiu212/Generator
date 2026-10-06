/** \file G4INCLNuclearMassTable.hh
 * \brief Functions that encapsulate a mass table
 *
 * \date 22nd October 2013
 * \author Davide Mancusi
 */

#ifndef G4INCLNuclearMassTable_HH
#define G4INCLNuclearMassTable_HH

#ifndef INCLXX_IN_GEANT4_MODE

#include <map>
#include <string>

namespace G4INCL {

  namespace NuclearMassTable {
    void initialize(const std::string &path, const double pMass, const double nMass);
    double getMass(const int A, const int Z);
    double getMass(const int A, const int Z, const int S);
    void deleteTable();
  }

}

#endif // INCLXX_IN_GEANT4_MODE

#endif // G4INCLNuclearMassTable_HH
