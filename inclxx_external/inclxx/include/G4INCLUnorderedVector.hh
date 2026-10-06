/*
 * \file G4INCLUnorderedVector.hh
 *
 * \date 2nd October 2014
 * \author Davide Mancusi
 */

#ifndef G4INCLUNORDEREDVECTOR_HH_
#define G4INCLUNORDEREDVECTOR_HH_

#include <vector>
#include <algorithm>

#if !defined(NDEBUG) && !defined(INCLXX_IN_GEANT4_MODE)
// Force instantiation of all the std::vector<Particle*> methods for debugging
// purposes
namespace G4INCL {
  class Particle;
}
template class std::vector<G4INCL::Particle*>;
#endif

namespace G4INCL {

  template<class T>
    class UnorderedVector : private std::vector<T> {
      public:
        UnorderedVector() {}
        using std::vector<T>::push_back;
        using std::vector<T>::pop_back;
        using std::vector<T>::size;
        using std::vector<T>::begin;
        using std::vector<T>::end;
        using std::vector<T>::rbegin;
        using std::vector<T>::rend;
        using std::vector<T>::front;
        using std::vector<T>::back;
        using std::vector<T>::clear;
        using std::vector<T>::empty;
        using std::vector<T>::insert;
        using std::vector<T>::erase;
        using std::vector<T>::operator[];
        using std::vector<T>::reserve;
        using std::vector<T>::resize;
        using std::vector<T>::at;
        using typename std::vector<T>::iterator;
        using typename std::vector<T>::reverse_iterator;
        using typename std::vector<T>::const_iterator;
        using typename std::vector<T>::const_reverse_iterator;

        void remove(const T &t) {
          const typename std::vector<T>::iterator removeMe = std::find(begin(), end(), t);
          assert(removeMe!=end());
          *removeMe = back();
          pop_back();
        }

        bool contains(const T &t) const {
          return (std::find(begin(), end(), t)!=end());
        }
    };

}

#endif // G4INCLUNORDEREDVECTOR_HH_
