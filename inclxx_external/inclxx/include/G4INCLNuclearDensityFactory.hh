#ifndef G4INCLNuclearDensityFactory_hh
#define G4INCLNuclearDensityFactory_hh 1

#include "G4INCLNuclearDensity.hh"
#include "G4INCLParticleSampler.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLInterpolationTable.hh"
#include <map>

namespace G4INCL {

  namespace NuclearDensityFactory {

    InterpolationTable *createRPCorrelationTable(const ParticleType t, const int A, const int Z);

    InterpolationTable *createRCDFTable(const ParticleType t, const int A, const int Z);

    InterpolationTable *createPCDFTable(const ParticleType t, const int A, const int Z);

    NuclearDensity const *createDensity(const int A, const int Z, const int S);

    void addRPCorrelationToCache(const int A, const int Z, const ParticleType t, InterpolationTable * const table);

    void addDensityToCache(const int A, const int Z, NuclearDensity * const density);

    void clearCache();

  }
}

#endif
