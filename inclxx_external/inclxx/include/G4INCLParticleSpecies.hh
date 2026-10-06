/*
 * G4INCLParticleSpecies.hh
 *
 *  \date Nov 25, 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLPARTICLESPECIES_HH_
#define G4INCLPARTICLESPECIES_HH_

#include "G4INCLParticleType.hh"
#include <string>

namespace G4INCL {

  class ParticleSpecies {
    public:
      /// \brief Convert a string to a particle species
      ParticleSpecies() :
        theType(UnknownParticle),
        theA(0),
        theZ(0),
        theS(0)
    {}
      ParticleSpecies(std::string const &pS);
      ParticleSpecies(ParticleType const t);
      ParticleSpecies(const int A, const int Z);
      ParticleSpecies(const int A, const int Z, const int S);
      // initialize particle with pdg code
      ParticleSpecies(int pdg);

      ParticleType theType;
      int theA, theZ, theS;
	  
      /** \brief Set a PDG Code (MONTE CARLO PARTICLE NUMBERING)
       *
       * \param 
       * \return integer (identifying number for each particle)
       */
      int getPDGCode() const;
	  
    private:
      /** \brief Parse a nuclide name
       *
       * Note: this function is UGLY. Look at it at your own peril.
       *
       * \param pS a normalised string (lowercase)
       */
      void parseNuclide(std::string const &pS);

      /** \brief Parse an element name
       *
       * Note: this function is UGLY. Look at it at your own peril.
       *
       * \param pS a normalised string (lowercase)
       * \return true if the parsing succeeded
       */
      bool parseElement(std::string const &pS);

      /** \brief Parse a IUPAC element name
       *
       * Note: this function is UGLY. Look at it at your own peril.
       *
       * \param s a normalised string (lowercase)
       * \return true if the parsing succeeded
       */
      bool parseIUPACElement(std::string const &s);
	  
  };

}

#endif /* G4INCLPARTICLESPECIES_HH_ */
