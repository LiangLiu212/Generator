/** \file G4INCLHornerFormEvaluator.hh
 * \brief Template-metaprogramming-based evaluator for polynomials in Horner form
 *
 * \date 3rd March 2014
 * \author Davide Mancusi
 */

#ifndef G4INCLHORNERFORMEVALUATOR_HH
#define G4INCLHORNERFORMEVALUATOR_HH

namespace G4INCL {

  template<int N>
    class HornerCoefficients {
      protected:
        double a[N];
      public:
        double &operator[](int i) { return a[i]; }
        const double &operator[](int i) const { return a[i]; }
    };

  struct HornerC1 : public HornerCoefficients<1> {
    HornerC1(
             const double a0
            ) {
      a[0] = a0;
    }
  };

  struct HornerC2 : public HornerCoefficients<2> {
    HornerC2(
             const double a0,
             const double a1
            ) {
      a[0] = a0;
      a[1] = a1;
    }
  };

  struct HornerC3 : public HornerCoefficients<3> {
    HornerC3(
             const double a0,
             const double a1,
             const double a2
            ) {
      a[0] = a0;
      a[1] = a1;
      a[2] = a2;
    }
  };

  struct HornerC4 : public HornerCoefficients<4> {
    HornerC4(
             const double a0,
             const double a1,
             const double a2,
             const double a3
            ) {
      a[0] = a0;
      a[1] = a1;
      a[2] = a2;
      a[3] = a3;
    }
  };

  struct HornerC5 : public HornerCoefficients<5> {
    HornerC5(
             const double a0,
             const double a1,
             const double a2,
             const double a3,
             const double a4
            ) {
      a[0] = a0;
      a[1] = a1;
      a[2] = a2;
      a[3] = a3;
      a[4] = a4;
    }
  };

  struct HornerC6 : public HornerCoefficients<6> {
    HornerC6(
             const double a0,
             const double a1,
             const double a2,
             const double a3,
             const double a4,
             const double a5
            ) {
      a[0] = a0;
      a[1] = a1;
      a[2] = a2;
      a[3] = a3;
      a[4] = a4;
      a[5] = a5;
    }
  };

  struct HornerC7 : public HornerCoefficients<7> {
    HornerC7(
             const double a0,
             const double a1,
             const double a2,
             const double a3,
             const double a4,
             const double a5,
             const double a6
            ) {
      a[0] = a0;
      a[1] = a1;
      a[2] = a2;
      a[3] = a3;
      a[4] = a4;
      a[5] = a5;
      a[6] = a6;
    }
  };

  struct HornerC8 : public HornerCoefficients<8> {
    HornerC8(
             const double a0,
             const double a1,
             const double a2,
             const double a3,
             const double a4,
             const double a5,
             const double a6,
             const double a7
            ) {
      a[0] = a0;
      a[1] = a1;
      a[2] = a2;
      a[3] = a3;
      a[4] = a4;
      a[5] = a5;
      a[6] = a6;
      a[7] = a7;
    }
  };

  template<int M>
    struct HornerEvaluator {
      template<int N>
        static double eval(const double x, HornerCoefficients<N> const &coeffs) {
          return coeffs[N-M] + x * HornerEvaluator<M-1>::eval(x, coeffs);
        }
    };

  template<>
    struct HornerEvaluator<1> {
      template<int N>
        static double eval(const double, HornerCoefficients<N> const &coeffs) {
          return coeffs[N-1];
        }
    };

}

#endif
