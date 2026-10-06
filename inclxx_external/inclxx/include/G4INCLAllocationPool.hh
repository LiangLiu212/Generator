/** \file G4INCLAllocationPool.hh
 * \brief Singleton for recycling allocation of instances of a given class
 *
 * \date 2nd October 2014
 * \author Davide Mancusi
 */

#ifndef G4INCLALLOCATIONPOOL_HH
#define G4INCLALLOCATIONPOOL_HH

#if defined(INCL_USE_ALLOCATION_POOL) || defined(INCLXX_IN_GEANT4_MODE)

#include <stack>
#include <new>
#include <cstddef>

namespace G4INCL {

  template<typename T>
    class AllocationPool {
      public:
        static AllocationPool &getInstance() {
          if(!theInstance)
            theInstance = new AllocationPool<T>;
          return *theInstance;
        }

        T *getObject() {
          if(theStack.empty())
            return static_cast<T*>(::operator new(sizeof(T)));
          else {
            T *t = theStack.top();
            theStack.pop();
            return t;
          }
        }

        void recycleObject(T *t) {
          theStack.push(t);
        }

        void clear() {
          while(!theStack.empty()) { /* Loop checking, 10.07.2015, D.Mancusi */
            ::operator delete(theStack.top());
            theStack.pop();
          }
        }

      protected:
        AllocationPool() {}
        virtual ~AllocationPool() {
          clear();
        }

        static G4ThreadLocal AllocationPool *theInstance;

        std::stack<T*> theStack;

    };

  template<typename T>
    G4ThreadLocal AllocationPool<T> *AllocationPool<T>::theInstance = 0;

}

#define INCL_DECLARE_ALLOCATION_POOL(T) \
  public: \
    static void *operator new(size_t /* s */) { \
      ::G4INCL::AllocationPool<T> &allocator = ::G4INCL::AllocationPool<T>::getInstance(); \
      return allocator.getObject(); \
    } \
    static void operator delete(void *a, size_t /* s */) { \
      ::G4INCL::AllocationPool<T> &allocator = ::G4INCL::AllocationPool<T>::getInstance(); \
      allocator.recycleObject(static_cast<T *>(a)); \
    }

#else // defined(INCL_USE_ALLOCATION_POOL) || defined(INCLXX_IN_GEANT4_MODE)
#define INCL_DECLARE_ALLOCATION_POOL(T)
#endif

#endif // G4INCLALLOCATIONPOOL_HH
