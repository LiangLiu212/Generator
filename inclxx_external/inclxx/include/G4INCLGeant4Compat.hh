#ifndef G4INCLGeant4Compat_hh
#define G4INCLGeant4Compat_hh 1

// If we are not inside Geant4 we use our own built-in definitions of the G4 types
#ifndef INCLXX_IN_GEANT4_MODE
#include <iostream>
#include <string>
typedef int G4int;
typedef double G4double;
typedef bool G4bool;
typedef long G4long;
typedef std::string G4String;
#define G4cout std::cout;
#define G4endl std::endl;
#endif

#endif
