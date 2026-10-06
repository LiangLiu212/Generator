#ifndef G4INCLGlobals_hh
#define G4INCLGlobals_hh 1

#include <cmath>
#include <string>
#include <vector>
#include "G4INCLParticleType.hh"
#include <iostream>

namespace G4INCL {
  class Particle;

  namespace PhysicalConstants {
    /// \brief \f$\hbar c\f$ [MeV*fm]
    const double hc = 197.328;

    /// \brief \f$\hbar^2 c^2\f$ [MeV^2*fm^2]
    const double hcSquared = hc*hc;

    /// \brief Fermi momentum [MeV/c]
    const double Pf = 1.37*hc;
    //  const double Pf = 1.36828*hc;

    /** \brief Coulomb conversion factor [MeV*fm]
     *
     * \f[ e^2/(4 pi epsilon_0) \f]
     */
    const double eSquared = 1.439964;
  }

  namespace Math {
    const double pi = 3.14159265358979323846264338328;
    const double twoPi = 2.0 * pi;
    const double tenPi = 10.0 * pi;
    const double piOverTwo = 0.5 * pi;
    const double oneOverSqrtTwo = 1./std::sqrt((double)2.);
    const double oneOverSqrtThree = 1./std::sqrt((double)3.);
    const double oneThird = 1./3.;
    const double twoThirds = 2./3.;
    const double sqrtFiveThirds = std::sqrt(5./3.);
    const double sqrtThreeFifths = std::sqrt(3./5.);

    inline double toDegrees(double radians) {
      return radians * (180.0 / pi);
    }

    inline int heaviside(int n) {
      if(n < 0) return 0;
      else return 1;
    }

    inline double pow13(double x) {
      return std::pow(x, oneThird);
    }

    inline double powMinus13(double x) {
      return std::pow(x, -oneThird);
    }

    inline double pow23(double x) {
      return std::pow(x, twoThirds);
    }

    inline double aSinH(double x) {
      return std::log(x + std::sqrt(x*x+1.));
    }

    /**
     * A simple sign function that allows us to port fortran code to c++ more easily.
     */
    template <typename T> inline int sign(const T t) {
      return t > 0 ? 1: t < 0 ? -1 : 0;
    }

    /// brief Return the largest of the two arguments
    template <typename T> inline T max(const T t1, const T t2) {
      return t1 > t2 ? t1 : t2;
    }

    /// brief Return the smallest of the two arguments
    template <typename T> inline T min(const T t1, const T t2) {
      return t1 < t2 ? t1 : t2;
    }

    /** \brief Cumulative distribution function for Gaussian
     *
     * A public-domain approximation taken from Abramowitz and Stegun. Applies
     * to a Gaussian with mean=0 and sigma=1.
     *
     * \param x a Gaussian variable
     */
    double gaussianCDF(const double x);

    /** \brief Generic cumulative distribution function for Gaussian
     *
     * A public-domain approximation taken from Abramowitz and Stegun. Applies
     * to a generic Gaussian.
     *
     * \param x a Gaussian variable
     * \param x0 mean of the Gaussian
     * \param sigma standard deviation of the Gaussian
     */
    double gaussianCDF(const double x, const double x0, const double sigma);

    /** \brief Inverse cumulative distribution function for Gaussian
     *
     * A public-domain approximation taken from Abramowitz and Stegun. Applies
     * to a Gaussian with mean=0 and sigma=1.
     *
     * \param x a uniform variate
     * \return a Gaussian variate
     */
    double inverseGaussianCDF(const double x);

    /// \brief Calculates arcsin with some tolerance on illegal arguments
    double arcSin(const double x);

    /// \brief Calculates arccos with some tolerance on illegal arguments
    double arcCos(const double x);
  }

  namespace ParticleConfig {
    bool isPair(Particle const * const p1, Particle const * const p2, ParticleType t1, ParticleType t2);
  }

#ifndef INCLXX_IN_GEANT4_MODE
  namespace String {
    void wrap(std::string &str, const size_t lineLength=78, const std::string &separators=" \t");
    void replaceAll(std::string &str, const std::string &from, const std::string &to, const size_t maxPosition=std::string::npos);
    std::vector<std::string> tokenize(std::string const &str, const std::string &delimiters);
    bool isInteger(std::string const &str);
    std::string expandPath(std::string const &path);
  }
#endif
}
#endif
