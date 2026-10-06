/*
 * G4INCLSurfaceAvatar.hh
 *
 *  \date Jun 8, 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLSURFACEAVATAR_HH_
#define G4INCLSURFACEAVATAR_HH_

#include "G4INCLIAvatar.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {

  /**
   * Surface avatar
   *
   * The reflection avatar is created when a particle reaches the boundary of the nucleus.
   * At this point it can either be reflected from the boundary or exit the nucleus.
   */
  class SurfaceAvatar: public G4INCL::IAvatar {
  public:
    SurfaceAvatar(G4INCL::Particle *aParticle, double time, G4INCL::Nucleus *aNucleus);
    virtual ~SurfaceAvatar();

    IChannel* getChannel();
    void fillFinalState(FinalState *fs);

    virtual void preInteraction();
    virtual void postInteraction(FinalState *);

    ParticleList getParticles() const {
      ParticleList theParticleList;
      theParticleList.push_back(theParticle);
      return theParticleList;
    }

    std::string dump() const;

    /// \brief Calculate the transmission probability for the particle
    double getTransmissionProbability(Particle const * const particle);

    /// \brief Get the cosine of the refraction angle (precalculated by initializeRefractionVariables)
    double getCosRefractionAngle() const { return cosRefractionAngle; }

    /// \brief Get the outgoing momentum (precalculated by initializeRefractionVariables)
    double getOutgoingMomentum() const { return particlePOut; }

  private:
    void initializeRefractionVariables(Particle const * const particle);

  private:
    G4INCL::Particle *theParticle;
    G4INCL::Nucleus *theNucleus;

    double particlePIn;
    double particlePOut;
    double particleTOut;
    double TMinusV;
    double TMinusV2;
    double particleMass;
    double sinIncidentAngle;
    double cosIncidentAngle;
    double sinRefractionAngle;
    double cosRefractionAngle;
    double refractionIndexRatio;
    bool internalReflection;

    INCL_DECLARE_ALLOCATION_POOL(SurfaceAvatar)
  };

}

#endif /* G4INCLSURFACEAVATAR_HH_ */
