/*
 * ThreeVector.hh
 *
 *  \date 4 June 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLThreeVector_hh
#define G4INCLThreeVector_hh 1

#include <string>
#include <sstream>
#include <cmath>

namespace G4INCL {

  class ThreeVector {
    public:
      ThreeVector()
        :x(0.0), y(0.0), z(0.0)
      {}

      ThreeVector(double ax, double ay, double az)
        :x(ax), y(ay), z(az)
      {}

      inline double getX() const { return x; }
      inline double getY() const { return y; }
      inline double getZ() const { return z; }

      inline double perp() const { return std::sqrt(x*x + y*y); }
      inline double perp2() const { return x*x + y*y; }
      /**
       * Get the length of the vector.
       */
      inline double mag() const { return std::sqrt(x*x + y*y + z*z); }

      /**
       * Get the square of the length.
       */
      inline double mag2() const { return (x*x + y*y + z*z); }

      /**
       * Theta angle
       */
      inline double theta() const {
        return x == 0.0 && y == 0.0 && z == 0.0 ? 0.0 : std::atan2(perp(),z);
      }

      /**
       * Phi angle
       */
      inline double phi() const {
        return x == 0.0 && y == 0.0 ? 0.0 : std::atan2(y,x);
      }

      /**
       * Dot product.
       */
      inline double dot(const ThreeVector &v) const {
        return (x*v.x + y*v.y + z*v.z);
      }

      /**
       * Vector product.
       */
      ThreeVector vector(const ThreeVector &v) const {
        return ThreeVector(
            y*v.z - z*v.y,
            z*v.x - x*v.z,
            x*v.y - y*v.x
            );
      }

      /// \brief Set the x coordinate
      inline void setX(double ax) { x =  ax; }

      /// \brief Set the y coordinate
      inline void setY(double ay) { y =  ay; }

      /// \brief Set the z coordinate
      inline void setZ(double az) { z =  az; }

      /// \brief Set all the coordinates
      inline void set(const double ax, const double ay, const double az) { x=ax; y=ay; z=az; }

      inline void operator+= (const ThreeVector &v) {
        x += v.x;
        y += v.y;
        z += v.z;
      }

      /// \brief Unary minus operator
      inline ThreeVector operator- () const {
        return ThreeVector(-x,-y,-z);
      }

      inline void operator-= (const ThreeVector &v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
      }

      template<typename T>
        inline void operator*= (const T &c) {
          x *= c;
          y *= c;
          z *= c;
        }

      template<typename T>
        inline void operator/= (const T &c) {
          const double oneOverC = 1./c;
          this->operator*=(oneOverC);
        }

      inline ThreeVector operator- (const ThreeVector &v) const {
        return ThreeVector(x-v.x, y-v.y, z-v.z);
      }

      inline ThreeVector operator+ (const ThreeVector &v) const {
        return ThreeVector(x+v.x, y+v.y, z+v.z);
      }

      /**
       * Divides all components of the vector with a constant number.
       */
      inline ThreeVector operator/ (const double C) const {
        const double oneOverC = 1./C;
        return ThreeVector(x*oneOverC, y*oneOverC, z*oneOverC);
      }

      inline ThreeVector operator* (const double C) const {
        return ThreeVector(x*C, y*C, z*C);
      }

      /** \brief Rotate the vector by a given angle around a given axis
       *
       * \param angle the rotation angle
       * \param axis the rotation axis, which must be a unit vector
       */
      inline void rotate(const double angle, const ThreeVector &axis) {
        // Use Rodrigues' formula
        const double cos = std::cos(angle);
        const double sin = std::sin(angle);
        (*this) = (*this) * cos + axis.vector(*this) * sin + axis * (axis.dot(*this)*(1.-cos));
      }

      /** \brief Return a vector orthogonal to this
       *
       * Simple algorithm from Hughes and Moeller, J. Graphics Tools 4 (1999)
       * 33.
       */
      ThreeVector anyOrthogonal() const {
        if(x<=y && x<=z)
          return ThreeVector(0., -z, y);
        else if(y<=x && y<=z)
          return ThreeVector(-z, 0., x);
        else
          return ThreeVector(-y, x, 0.);
      }

      std::string print() const {
        std::stringstream ss;
        ss <<"(x = " << x << "   y = " << y << "   z = " << z <<")";
        return ss.str();
      }

      std::string dump() const {
        std::stringstream ss;
        ss <<"(vector3 " << x << " " << y << " " << z << ")";
        return ss.str();
      }

    private:
      double x, y, z; //> Vector components
  };

}

#endif
