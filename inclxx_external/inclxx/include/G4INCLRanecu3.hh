/*
 * \file G4INCLRanecu3.hh
 *
 *  \date 25 October 2013
 * \author Davide Mancusi
 */

#ifndef G4INCLRanecu3_HH_
#define G4INCLRanecu3_HH_

#include "G4INCLIRandomGenerator.hh"
#include <cassert>

namespace G4INCL {

  /** \brief Extended Ranecu-type RNG class
   *
   * This generator implements the C++ version of the generator suggested by
   * Badal and Sempau, Comp. Phys. Comm. 175 (2006) 440. It uses three 32-bit
   * seeds and has a period of ~5E27.
   */
  class Ranecu3 : public G4INCL::IRandomGenerator {
  public:
    Ranecu3();
    Ranecu3(const Random::SeedVector &sv);
    virtual ~Ranecu3();

    Random::SeedVector getSeeds() {
      Random::SeedVector sv;
      sv.push_back(iseed1);
      sv.push_back(iseed2);
      sv.push_back(iseed3);
      return sv;
    }

    void setSeeds(const Random::SeedVector &sv) {
      assert(sv.size()>=3);
      iseed1 = sv[0];
      iseed2 = sv[1];
      iseed3 = sv[2];
    }

    double flat();

  private:
    int iseed1, iseed2, iseed3;
    int i1, i2, i3, iz;
    const double uscale;
    const int m1, m2, m3;
    const int a1, a2, a3;
    const int q1, q2, q3;
    const int r1, r2, r3;
  };

}

#endif /* G4INCLRanecu3_HH_ */
