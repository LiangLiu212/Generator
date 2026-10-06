#ifndef G4INCLLorentzVector_hh
#define G4INCLLorentzVector_hh

#include "G4INCLThreeVector.hh"
#include "G4INCLKinematicsUtils.hh"

namespace G4INCL {

  class LorentzVector {
    public:
      LorentzVector() :
        ee(0.0)
    {}

      LorentzVector(const double aX, const double aY, const double aZ, const double aE) :
        ee(aE),
        p(aX,aY,aZ)
    {}

      LorentzVector(const ThreeVector &aP, const double aE) :
        ee(aE),
        p(aP)
    {}

      double e() const { return ee; }
      double t() const { return ee; }
      double x() const { return p.getX(); }
      double y() const { return p.getY(); }
      double z() const { return p.getZ(); }

      ThreeVector const &momentum() const { return p; }

      double m2() const { return KinematicsUtils::squareInvariantMass(ee,p); }
      double restMass2() const { return m2(); };
      double m() const {
        const double mass2=m2();
        if(mass2>0)
          return std::sqrt(mass2);
        else
          return -std::sqrt(-mass2);
      }
      double restMass() const { return m(); }

      void setE(const double aE) { ee=aE; }
      void setT(const double aE) { ee=aE; }
      void setX(const double aX) { p.setX(aX); }
      void setY(const double aY) { p.setY(aY); }
      void setZ(const double aZ) { p.setZ(aZ); }

      void setVect(const ThreeVector &aVect) { p=aVect; }

      LorentzVector &boost(double bx, double by, double bz){
        double b2 = bx*bx + by*by + bz*bz;
        double ggamma = 1.0 / std::sqrt(1.0 - b2);
        double bp = bx*x() + by*y() + bz*z();
        double gamma2 = b2 > 0 ? (ggamma - 1.0)/b2 : 0.0;

        setX(x() + gamma2*bp*bx + ggamma*bx*t());
        setY(y() + gamma2*bp*by + ggamma*by*t());
        setZ(z() + gamma2*bp*bz + ggamma*bz*t());
        setT(ggamma*(t() + bp));
        return *this;
      }

      LorentzVector &boost(const ThreeVector beta) {
        return boost(beta.getX(), beta.getY(), beta.getZ());
      }

      ThreeVector boostVector() const {
        if (ee == 0) {
          if (p.mag2() == 0) {
            return ThreeVector(0,0,0);
          } else {
            std::cerr << "LorentzVector::boostVector() - "
              << "boostVector computed for LorentzVector with t=0 -- infinite result"
              << std::endl;
            return p/ee;
          }
        }
        if (restMass2() <= 0) {
          std::cerr << "LorentzVector::boostVector() - "
            << "boostVector computed for a non-timelike LorentzVector " << std::endl;
          // result will make analytic sense but is physically meaningless
        }
        return p * (1./ee);
      }

    private:
      double ee;
      ThreeVector p;
  };

}

#endif // G4INCLLorentzVector_hh
