/** \file G4INCLInterpolationTable.hh
 * \brief Simple interpolation table
 *
 * \date 30 January 2014
 * \author Davide Mancusi
 */

#ifndef G4INCLINTERPOLATIONTABLE_HH_
#define G4INCLINTERPOLATIONTABLE_HH_

#include "G4INCLIFunction1D.hh"
#include <algorithm>
#include <functional>
#include <sstream>

namespace G4INCL {

  /// \brief Interpolation node
  class InterpolationNode {
    public:
      InterpolationNode(const double x0, const double y0, const double yPrime0) :
        x(x0),
        y(y0),
        yPrime(yPrime0)
    {}

      virtual ~InterpolationNode() {}

      bool operator<(const InterpolationNode &rhs) const {
        return (x < rhs.x);
      }

      bool operator<=(const InterpolationNode &rhs) const {
        return (x <= rhs.x);
      }

      bool operator>(const InterpolationNode &rhs) const {
        return (x > rhs.x);
      }

      bool operator>=(const InterpolationNode &rhs) const {
        return (x >= rhs.x);
      }

      double getX() const { return x; }
      double getY() const { return y; }
      double getYPrime() const { return yPrime; }

      void setX(const double x0) { x=x0; }
      void setY(const double y0) { y=y0; }
      void setYPrime(const double yPrime0) { yPrime=yPrime0; }

      std::string print() const {
        std::stringstream message;
        message << "x, y, yPrime: " << x << '\t' << y << '\t' << yPrime << '\n';
        return message.str();
      }

    protected:
      /// \brief abscissa
      double x;
      /// \brief function value
      double y;
      /// \brief function derivative
      double yPrime;
  };

  /// \brief Class for interpolating the  of a 1-dimensional function
  class InterpolationTable : public IFunction1D {
    public:
      InterpolationTable(std::vector<double> const &x, std::vector<double> const &y);
      virtual ~InterpolationTable() {}

      std::size_t getNumberOfNodes() const { return nodes.size(); }

      std::vector<double> getNodeAbscissae() const;

      std::vector<double> getNodeValues() const;

      double operator()(const double x) const;

      std::string print() const;

    protected:
      InterpolationTable();

      /// \brief Initialise the values of the node derivatives
      void initDerivatives();

      /// \brief Interpolating nodes
      std::vector<InterpolationNode> nodes;

  };

}

#endif // G4INCLINTERPOLATIONTABLE_HH_
