/** \file G4INCLINuclearPotential.hh
 * \brief Abstract interface to the nuclear potential.
 *
 * NuclearPotential-like classes should provide access to the value of the
 * potential of a particle in a particular context. For example, an instance of
 * a NuclearPotential class should be associated to every nucleus.
 *
 * \date 17 January 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLINUCLEARPOTENTIAL_HH
#define G4INCLINUCLEARPOTENTIAL_HH 1

#include "G4INCLParticle.hh"
#include "G4INCLRandom.hh"
#include "G4INCLDeuteronDensity.hh"
#include <map>
#include <cassert>

namespace G4INCL {

  namespace NuclearPotential {
    class INuclearPotential {
      public:
        INuclearPotential(const int A, const int Z, const bool pionPot) :
          theA(A),
          theZ(Z),
          pionPotential(pionPot)
        {
          if(pionPotential) {
            const double ZOverA = ((double) theZ) / ((double) theA);
            // As in INCL4.6, use the r0*A^(1/3) formula to estimate vc
            const double r = 1.12*Math::pow13((double)theA);

            const double xsi = 1. - 2.*ZOverA;
            const double vc = 1.25*PhysicalConstants::eSquared*theZ/r;
            vPiPlus = vPionDefault + 71.*xsi - vc;
            vPiZero = vPionDefault;
            vPiMinus = vPionDefault - 71.*xsi + vc;
            vKPlus = vKPlusDefault;
            vKZero = vKPlusDefault + 10.; // Hypothesis to be check
            vKMinus = vKMinusDefault;
            vKZeroBar = vKMinusDefault - 10.; // Hypothesis to be check
          } else {
            vPiPlus = 0.0;
            vPiZero = 0.0;
            vPiMinus = 0.0;
            vKPlus = 0.0;
            vKZero = 0.0;
            vKMinus = 0.0;
            vKZeroBar = 0.0;
          }
        }

        virtual ~INuclearPotential() {}

        /// \brief Do we have a pion potential?
        bool hasPionPotential() const { return pionPotential; }

        virtual double computePotentialEnergy(const Particle * const p) const = 0;

        /** \brief Return the Fermi energy for a particle.
         *
         * \param p pointer to a Particle
         * \return Fermi energy for that particle type
         **/
        inline double getFermiEnergy(const Particle * const p) const {
          std::map<ParticleType, double>::const_iterator i = fermiEnergy.find(p->getType());
          assert(i!=fermiEnergy.end());
          return i->second;
        }

        /** \brief Return the Fermi energy for a particle type.
         *
         * \param t particle type
         * \return Fermi energy for that particle type
         **/
        inline double getFermiEnergy(const ParticleType t) const {
          std::map<ParticleType, double>::const_iterator i = fermiEnergy.find(t);
          assert(i!=fermiEnergy.end());
          return i->second;
        }

        /** \brief Return the separation energy for a particle.
         *
         * \param p pointer to a Particle
         * \return separation energy for that particle type
         **/
        inline double getSeparationEnergy(const Particle * const p) const {
          std::map<ParticleType, double>::const_iterator i = separationEnergy.find(p->getType());
          assert(i!=separationEnergy.end());
          return i->second;
        }

        /** \brief Return the separation energy for a particle type.
         *
         * \param t particle type
         * \return separation energy for that particle type
         **/
        inline double getSeparationEnergy(const ParticleType t) const {
          std::map<ParticleType, double>::const_iterator i = separationEnergy.find(t);
          assert(i!=separationEnergy.end());
          return i->second;
        }

        /** \brief Return the Fermi momentum for a particle.
         *
         * \param p pointer to a Particle
         * \return Fermi momentum for that particle type
         **/
        inline double getFermiMomentum(const Particle * const p) const {
          if(p->isDelta()) {
            const double Tf = getFermiEnergy(p), mass = p->getMass();
            return std::sqrt(Tf*(Tf+2.*mass));
          } else {
            std::map<ParticleType, double>::const_iterator i = fermiMomentum.find(p->getType());
            assert(i!=fermiMomentum.end());
            return i->second;
          }
        }

        /** \brief Return the Fermi momentum for a particle type.
         *
         * \param t particle type
         * \return Fermi momentum for that particle type
         **/
        inline double getFermiMomentum(const ParticleType t) const {
          assert(t!=DeltaPlusPlus && t!=DeltaPlus && t!=DeltaZero && t!=DeltaMinus);
          std::map<ParticleType, double>::const_iterator i = fermiMomentum.find(t);
          return i->second;
        }

      protected:
        /// \brief Compute the potential energy for the given pion.
        double computePionPotentialEnergy(const Particle * const p) const {
          assert(p->getType()==PiPlus || p->getType()==PiZero || p->getType()==PiMinus);
          if(pionPotential && !p->isOutOfWell()) {
            switch( p->getType() ) {
              case PiPlus:
                return vPiPlus;
                break;
              case PiZero:
                return vPiZero;
                break;
              case PiMinus:
                return vPiMinus;
                break;
              default: // Pion potential is defined and non-zero only for pions
                return 0.0;
                break;
            }
          }
          else
            return 0.0;
        }

      protected:
        /// \brief Compute the potential energy for the given kaon.
        double computeKaonPotentialEnergy(const Particle * const p) const {
          assert(p->getType()==KPlus || p->getType()==KZero || p->getType()==KZeroBar || p->getType()==KMinus|| p->getType()==KShort|| p->getType()==KLong);
          if(pionPotential && !p->isOutOfWell()) { // if pionPotental false -> kaonPotential false
            switch( p->getType() ) {
              case KPlus:
                return vKPlus;
                break;
              case KZero:
                return vKZero;
                break;
              case KZeroBar:
                return vKZeroBar;
                break;
              case KShort:
              case KLong:
                 return 0.0; // Should never be in the nucleus
                 break;
               case KMinus:
                 return vKMinus;
                 break;
               default:
                 return 0.0;
                 break;
               }
            }
        else
          return 0.0;
        }

      protected:
        /// \brief Compute the potential energy for the given pion resonances (Eta, Omega and EtaPrime and Gamma also).
        double computePionResonancePotentialEnergy(const Particle * const p) const {
          assert(p->getType()==Eta || p->getType()==Omega || p->getType()==EtaPrime || p->getType()==Photon);
          if(pionPotential && !p->isOutOfWell()) {
            switch( p->getType() ) {
              case Eta:
//jcd             return vPiZero;
//jcd              return vPiZero*1.5;
                return 0.0; // (JCD: seems to give better results)
                break;
              case Omega:
                return 15.0; // S.Friedrich et al., Physics Letters B736(2014)26-32. (V. Metag in Hyperfine Interact (2015) 234:25-31 gives 29 MeV)
                break;
              case EtaPrime:
                return 37.0; // V. Metag in Hyperfine Interact (2015) 234:25-31
                break;
              case Photon:
                return 0.0;
                break;
              default: 
                return 0.0;
                break;
            }
          }
          else
            return 0.0;
        }

      protected:
        /// \brief The mass number of the nucleus
        const int theA;
        /// \brief The charge number of the nucleus
        const int theZ;
      private:
        const bool pionPotential;
        double vPiPlus, vPiZero, vPiMinus;
        static const double vPionDefault;
        double vKPlus, vKZero, vKZeroBar, vKMinus;
        static const double vKPlusDefault;
        static const double vKMinusDefault;
      protected:
        /* \brief map of Fermi energies per particle type */
        std::map<ParticleType,double> fermiEnergy;
        /* \brief map of Fermi momenta per particle type */
        std::map<ParticleType,double> fermiMomentum;
        /* \brief map of separation energies per particle type */
        std::map<ParticleType,double> separationEnergy;

    };



    /** \brief Create an INuclearPotential object
     *
     * This is the method that should be used to instantiate objects derived
     * from INuclearPotential. It uses a caching mechanism to minimise
     * thrashing and speed up the code.
     *
     * \param type the type of the potential to be created
     * \param theA mass number of the nucleus
     * \param theZ charge number of the nucleus
     * \param pionPotential whether pions should also feel the potential
     * \return a pointer to the nuclear potential
     */
    INuclearPotential const *createPotential(const PotentialType type, const int theA, const int theZ, const bool pionPotential);

    /// \brief Clear the INuclearPotential cache
    void clearCache();

  }

}

#endif /* G4INCLINUCLEARPOTENTIAL_HH_ */
