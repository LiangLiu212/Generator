/*
 * G4INCLIChannel.h
 *
 *  \date Jun 5, 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLFINALSTATE_H_
#define G4INCLFINALSTATE_H_

#include "G4INCLParticle.hh"
#include <string>

namespace G4INCL {

  enum FinalStateValidity {
    ValidFS,
    PauliBlockedFS,
    NoEnergyConservationFS,
    ParticleBelowFermiFS,
    ParticleBelowZeroFS
  };

  /**
   * Final state of an interaction
   */
  class FinalState {
  public:
    FinalState();
    virtual ~FinalState();

    void reset();

    void setTotalEnergyBeforeInteraction(double E) { totalEnergyBeforeInteraction = E; };
    double getTotalEnergyBeforeInteraction() const { return totalEnergyBeforeInteraction; };

    void addModifiedParticle(Particle *p);
    void addOutgoingParticle(Particle *p);
    void addDestroyedParticle(Particle *p);
    void addCreatedParticle(Particle *p);
    void addEnteringParticle(Particle *p);

    ParticleList const &getModifiedParticles() const;
    ParticleList const &getOutgoingParticles() const;
    ParticleList const &getDestroyedParticles() const;
    ParticleList const &getCreatedParticles() const;
    ParticleList const &getEnteringParticles() const;
    
    ParticleList &getSrcModifiedParticles();

    FinalStateValidity getValidity() const { return validity; }
    void makeValid() { validity = ValidFS; }
    void makePauliBlocked() { validity = PauliBlockedFS; }
    void makeNoEnergyConservation() { validity = NoEnergyConservationFS; }
    void makeParticleBelowFermi() { validity = ParticleBelowFermiFS; }
    void makeParticleBelowZero() { validity = ParticleBelowZeroFS; }

    std::string print() const;

  private:
    ParticleList outgoing, created, destroyed, modified, entering;
    double totalEnergyBeforeInteraction;
    FinalStateValidity validity;
  };

}

#endif /* G4INCLFINALSTATE_H_ */
