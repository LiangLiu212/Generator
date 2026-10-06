/*
 * \file G4INCLRandomSeedVector.hh
 *
 * \date 17 May 2013
 * \author Davide Mancusi
 */

#ifndef G4INCLRANDOMSEEDVECTOR_HH_
#define G4INCLRANDOMSEEDVECTOR_HH_

#include <vector>
#include <ostream>

namespace G4INCL {

  namespace Random {

    class SeedVector : private std::vector<int> {
      public:
        using std::vector<int>::at;
        using std::vector<int>::operator[];
        using std::vector<int>::size;
        using std::vector<int>::resize;
        using std::vector<int>::push_back;
        using std::vector<int>::clear;
        using std::vector<int>::begin;
        using std::vector<int>::end;

        friend std::ostream &operator<<(std::ostream &out, SeedVector const &sv);

    };

  }
}

#endif /* G4INCLIRANDOMSEEDVECTOR_HH_ */

