/** \file G4INCLNaturalIsotopicDistributions.hh
 * \brief Classes that stores isotopic abundances
 *
 * \date 21st October 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLISOTOPICDISTRIBUTION_HH_
#define G4INCLISOTOPICDISTRIBUTION_HH_

#include <vector>
#include <map>

namespace G4INCL {

  /// \brief Holds an isotope and an abundance
  struct Isotope {
    Isotope(const int A, const double abundance);
    int theA;
    double theAbundance;
  };

  typedef std::vector<Isotope> IsotopeVector;
  typedef IsotopeVector::iterator IsotopeIter;

  /// \brief Class that stores isotopic abundances for a given element
  class IsotopicDistribution {
    public:

      /// \brief Constructor
      IsotopicDistribution(IsotopeVector const &aVector);

      /// \brief Draw a random isotope based on the abundance vector
      int drawRandomIsotope() const;

      /// \brief Get the isotope vector
      IsotopeVector const &getIsotopes() const;

    private:
      IsotopeVector theIsotopes;
  };

  /// \brief Class that stores isotopic abundances for a given element
  class NaturalIsotopicDistributions {
    public:

      /// \brief Constructor
      NaturalIsotopicDistributions();

      /** \brief Draw a random isotope
       *
       * \param Z the element number
       */
      int drawRandomIsotope(int const Z) const;

      /** \brief Get an isotopic distribution
       *
       * \param Z the element number
       */
      IsotopicDistribution const &getIsotopicDistribution(int const Z) const;

    private:
      std::map<int, IsotopicDistribution> theDistributions;
  };

}
#endif // G4INCLISOTOPICDISTRIBUTION_HH_
