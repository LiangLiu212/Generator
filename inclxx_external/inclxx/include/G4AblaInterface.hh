#pragma once

#ifdef ABLAXX_IN_GEANT4_MODE

#include "G4Abla.hh"
#include "G4AblaFission.hh"
#include "G4Fragment.hh"
#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4ReactionProduct.hh"
#include "G4VPreCompoundModel.hh"

class G4AblaInterface : public G4VPreCompoundModel
{
  public:
    G4AblaInterface();
    virtual ~G4AblaInterface();

    virtual G4ReactionProductVector* DeExcite(G4Fragment& aFragment);

    virtual G4HadFinalState* ApplyYourself(G4HadProjectile const&, G4Nucleus&) { return NULL; }

    virtual void ModelDescription(std::ostream& outFile) const;
    virtual void DeExciteModelDescription(std::ostream& outFile) const;

  private:
    G4AblaOutput* ablaResult;
    G4Abla* theABLAModel;
    G4long eventNumber;

    /// \brief Convert an Abla particle to a G4ReactionProduct
    G4ReactionProduct* toG4Particle(G4int A, G4int Z, G4double kinE, G4double px, G4double py, G4double pz) const;

    /// \brief Convert A and Z to a G4ParticleDefinition
    G4ParticleDefinition* toG4ParticleDefinition(G4int A, G4int Z) const;
};

#endif // ABLAXX_IN_GEANT4_MODE
