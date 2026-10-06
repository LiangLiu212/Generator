#ifndef G4INCLHFB_hh
#define G4INCLHFB_hh 1

#include <string>
#include <vector>
#include <cassert>

#include "G4INCLParticleType.hh"

namespace G4INCL {

    const int TableZSize = 128;
    const int TableASize = 300;

  namespace HFB {
#ifdef INCLXX_IN_GEANT4_MODE
    void initialize();
#else
    void initialize(const std::string &path);
#endif
    /// \brief Get the radius and diffuseness parameters from HFB calculations
    double getRadiusParameterHFB(const ParticleType t, const int A, const int Z);
    double getSurfaceDiffusenessHFB(const ParticleType t, const int A, const int Z);
  }
}

#endif
