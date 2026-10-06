/*
 * G4INCLBinaryCollisionAvatar.hh
 *
 *  \date Jun 5, 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLBINARYCOLLISIONAVATAR_HH_
#define G4INCLBINARYCOLLISIONAVATAR_HH_

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLInteractionAvatar.hh"
#include "G4INCLAllocationPool.hh"
#include "G4INCLConfig.hh"

namespace G4INCL {

  class BinaryCollisionAvatar: public InteractionAvatar {
  public:
    BinaryCollisionAvatar(double, double, G4INCL::Nucleus*, G4INCL::Particle*, G4INCL::Particle*);
    virtual ~BinaryCollisionAvatar();
    G4INCL::IChannel* getChannel();
    ParticleList getParticles() const {
      ParticleList theParticleList;
      theParticleList.push_back(particle1);
      theParticleList.push_back(particle2);
      return theParticleList;
    };

    virtual void preInteraction();
    virtual void postInteraction(FinalState *);

    std::string dump() const;

    static void setCutNN(const double c) {
      cutNN = c;
      cutNNSquared = cutNN*cutNN;
    }

    static double getCutNN() { return cutNN; }

    static double getCutNNSquared() { return cutNNSquared; }
    
    /// \brief Get the global bias factor
    static double getBias() { return bias; }
    
    /// \brief Set the global bias factor
    static void setBias(const double b) { bias=b; }

  private:
    static G4ThreadLocal double cutNN;
    static G4ThreadLocal double cutNNSquared;
    static G4ThreadLocal double bias;

    double theCrossSection;
    bool isParticle1Spectator;
    bool isParticle2Spectator;
    bool isElastic;
    bool isStrangeProduction;
    
    void generateSrcPairsMethod(ParticleList &theList, const int then, const int thez, G4INCL::Particle*, G4INCL::Particle*);

    INCL_DECLARE_ALLOCATION_POOL(BinaryCollisionAvatar)
  };

}

#endif /* G4INCLBINARYCOLLISIONAVATAR_HH_ */
