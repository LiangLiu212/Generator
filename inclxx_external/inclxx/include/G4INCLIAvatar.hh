/*
 * IAvatar.hh
 *
 *  \date 4 June 2009
 * \author Pekka Kaitaniemi
 */

#ifndef IAVATAR_HH_
#define IAVATAR_HH_

#include "G4INCLIChannel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLUnorderedVector.hh"
#include "G4INCLAllocationPool.hh"
#include <string>

#if !defined(NDEBUG) && !defined(INCLXX_IN_GEANT4_MODE)
// Force instantiation of all the std::vector<IAvatar*> methods for debugging
// purposes
namespace G4INCL {
  class IAvatar;
}
template class std::vector<G4INCL::IAvatar*>;
#endif

namespace G4INCL {

  enum AvatarType {SurfaceAvatarType,
           CollisionAvatarType,
           DecayAvatarType,
           ParticleEntryAvatarType,
                   UnknownAvatarType};

  class IAvatar {
  public:
    IAvatar();
    IAvatar(double time);
    virtual ~IAvatar();

    virtual G4INCL::IChannel* getChannel() = 0;
    FinalState *getFinalState();
    void fillFinalState(FinalState *fs);
    virtual void preInteraction() = 0;
    virtual void postInteraction(FinalState *) = 0;

    double getTime() const { return theTime; };
    double setTime(double t) const { return t; };

    virtual ParticleList getParticles() const = 0;

    virtual std::string dump() const = 0;

    AvatarType getType() const { return type; };
    bool isACollision() const { return (type==CollisionAvatarType); };
    bool isADecay() const { return (type==DecayAvatarType); };
    void setType(AvatarType t) { type = t; };
    long getID() const { return ID; };

    std::string toString();
  private:
    long ID;
    AvatarType type;
    static G4ThreadLocal long nextID;
  protected:
    double theTime;

    INCL_DECLARE_ALLOCATION_POOL(IAvatar)
  };

  typedef UnorderedVector<IAvatar*> IAvatarList;
  typedef UnorderedVector<IAvatar*>::const_iterator IAvatarIter;
  typedef UnorderedVector<IAvatar*>::iterator IAvatarMutableIter;

}

#endif /* IAVATAR_HH_ */
