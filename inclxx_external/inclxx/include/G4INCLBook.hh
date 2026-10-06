#ifndef G4INCLBook_hh
#define G4INCLBook_hh 1

#include <map>
#include "G4INCLIAvatar.hh"

namespace G4INCL {
  class Book {
  public:
    Book() {
      reset();
    }
    ~Book() {};

    void reset() {
      nAcceptedCollisions = 0;
      nBlockedCollisions = 0;
      nAcceptedDecays = 0;
      nBlockedDecays = 0;
      nAcceptedSrc=0;
      nSrcPairs = 0;
      currentTime = 0.0;
      firstCollisionTime = 0.0;
      firstCollisionXSec = 0.0;
      firstCollisionSpectatorPosition = 0.0;
      firstCollisionSpectatorMomentum = 0.0;
      firstCollisionIsElastic = false;
      nAvatars[SurfaceAvatarType] = 0;
      nAvatars[CollisionAvatarType] = 0;
      nAvatars[DecayAvatarType] = 0;
      nAvatars[ParticleEntryAvatarType] = 0;
      nCascading = 0;
      nEmittedClusters = 0;
      nEnergyViolationInteraction = 0;
    };

    void incrementAcceptedSrcCollisions() { nAcceptedSrc++; };
    void incrementAcceptedCollisions() { nAcceptedCollisions++; };
    void incrementBlockedCollisions() { nBlockedCollisions++; };
    void incrementAcceptedDecays() { nAcceptedDecays++; };
    void incrementSrcPairs() { nSrcPairs++; };
    void incrementBlockedDecays() { nBlockedDecays++; };
    void incrementAvatars(AvatarType type) { nAvatars[type]++; };
    void incrementCascading() { nCascading++; }
    void decrementCascading() { nCascading--; }
    void incrementEmittedClusters() { nEmittedClusters++; }
    void incrementEnergyViolationInteraction() { nEnergyViolationInteraction++; }

    void setFirstCollisionTime(const double t) { firstCollisionTime = t; };
    double getFirstCollisionTime() const { return firstCollisionTime; };

    void setFirstCollisionXSec(const double x) { firstCollisionXSec = x; };
    double getFirstCollisionXSec() const { return firstCollisionXSec; };

    void setFirstCollisionSpectatorPosition(const double x) { firstCollisionSpectatorPosition = x; };
    double getFirstCollisionSpectatorPosition() const { return firstCollisionSpectatorPosition; };

    void setFirstCollisionSpectatorMomentum(const double x) { firstCollisionSpectatorMomentum = x; };
    double getFirstCollisionSpectatorMomentum() const { return firstCollisionSpectatorMomentum; };

    void setFirstCollisionIsElastic(const bool e) { firstCollisionIsElastic = e; };
    bool getFirstCollisionIsElastic() const { return firstCollisionIsElastic; };

    void setCurrentTime(double t) { currentTime = t; };
    double getCurrentTime() const { return currentTime; };
    
    void setSrcPairs(int n) { nSrcPairs=n; };
    int getSrcPairs() const { return nSrcPairs; };
    
    void setAcceptedSrcCollisions(int n) { nAcceptedSrc=n; };  
    int getAcceptedSrcCollisions() const { return nAcceptedSrc; };

    int getAcceptedCollisions() const { return nAcceptedCollisions; };
    int getBlockedCollisions() const {return nBlockedCollisions; };
    int getAcceptedDecays() const { return nAcceptedDecays; };
    int getBlockedDecays() const {return nBlockedDecays; };
    int getAvatars(AvatarType type) const { return nAvatars.find(type)->second; };
    int getCascading() const { return nCascading; };
    int getEmittedClusters() const { return nEmittedClusters; };
    int getEnergyViolationInteraction() const { return nEnergyViolationInteraction; };

  private:
    int nAcceptedSrc;
    int nAcceptedCollisions;
    int nBlockedCollisions;
    int nAcceptedDecays;
    int nBlockedDecays;
    int nSrcPairs;
    double currentTime;
    double firstCollisionTime;
    double firstCollisionXSec;
    double firstCollisionSpectatorPosition;
    double firstCollisionSpectatorMomentum;
    bool firstCollisionIsElastic;
    std::map<AvatarType,int> nAvatars;
    int nCascading;
    int nEmittedClusters;
    int nEnergyViolationInteraction;
  };
}

#endif
