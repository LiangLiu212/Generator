/*
 * G4INCLParticle.hh
 *
 *  \date Jun 5, 2009
 * \author Pekka Kaitaniemi
 */

#ifndef PARTICLE_HH_
#define PARTICLE_HH_

#include "G4INCLThreeVector.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLParticleType.hh"
#include "G4INCLParticleSpecies.hh"
#include "G4INCLLogger.hh"
#include "G4INCLUnorderedVector.hh"
#include "G4INCLAllocationPool.hh"
#include <sstream>
#include <string>

namespace G4INCL {

  class Particle;

  class ParticleList : public UnorderedVector<Particle*> {
    public:
      void rotatePositionAndMomentum(const double angle, const ThreeVector &axis) const;
      void rotatePosition(const double angle, const ThreeVector &axis) const;
      void rotateMomentum(const double angle, const ThreeVector &axis) const;
      void boost(const ThreeVector &b) const;
      double getParticleListBias() const;
      std::vector<int> getParticleListBiasVector() const;
  };

  typedef ParticleList::const_iterator ParticleIter;
  typedef ParticleList::iterator       ParticleMutableIter;

  class Particle {
  public:
    Particle();
    Particle(ParticleType t, double energy, ThreeVector const &momentum, ThreeVector const &position);
    Particle(ParticleType t, ThreeVector const &momentum, ThreeVector const &position);
    virtual ~Particle() {}

    /** \brief Copy constructor
     *
     * Does not copy the particle ID.
     */
    Particle(const Particle &rhs) :
      theZ(rhs.theZ),
      theA(rhs.theA),
      theS(rhs.theS),
      theParticipantType(rhs.theParticipantType),
      theType(rhs.theType),
      theEnergy(rhs.theEnergy),
      theFrozenEnergy(rhs.theFrozenEnergy),
      theMomentum(rhs.theMomentum),
      theFrozenMomentum(rhs.theFrozenMomentum),
      thePosition(rhs.thePosition),
      nCollisions(rhs.nCollisions),
      nDecays(rhs.nDecays),
      nSrcPair(rhs.nSrcPair),
      thePotentialEnergy(rhs.thePotentialEnergy),
      rpCorrelated(rhs.rpCorrelated),
      uncorrelatedMomentum(rhs.uncorrelatedMomentum),
      theParticleBias(rhs.theParticleBias),
      theNKaon(rhs.theNKaon),
#ifdef INCLXX_IN_GEANT4_MODE
      theParentResonancePDGCode(rhs.theParentResonancePDGCode),
      theParentResonanceID(rhs.theParentResonanceID),
#endif
      theHelicity(rhs.theHelicity),
      emissionTime(rhs.emissionTime),
      outOfWell(rhs.outOfWell),
      theSrcPartner(rhs.theSrcPartner),
      theMass(rhs.theMass)
      {
        if(rhs.thePropagationEnergy == &(rhs.theFrozenEnergy))
          thePropagationEnergy = &theFrozenEnergy;
        else
          thePropagationEnergy = &theEnergy;
        if(rhs.thePropagationMomentum == &(rhs.theFrozenMomentum))
          thePropagationMomentum = &theFrozenMomentum;
        else
          thePropagationMomentum = &theMomentum;
        // ID intentionally not copied
        ID = nextID++;
        
        theBiasCollisionVector = rhs.theBiasCollisionVector;
      }

  protected:
    /// \brief Helper method for the assignment operator
    void swap(Particle &rhs) {
      std::swap(theZ, rhs.theZ);
      std::swap(theA, rhs.theA);
      std::swap(theS, rhs.theS);
      std::swap(theParticipantType, rhs.theParticipantType);
      std::swap(theType, rhs.theType);
      if(rhs.thePropagationEnergy == &(rhs.theFrozenEnergy))
        thePropagationEnergy = &theFrozenEnergy;
      else
        thePropagationEnergy = &theEnergy;
      std::swap(theEnergy, rhs.theEnergy);
      std::swap(theFrozenEnergy, rhs.theFrozenEnergy);
      if(rhs.thePropagationMomentum == &(rhs.theFrozenMomentum))
        thePropagationMomentum = &theFrozenMomentum;
      else
        thePropagationMomentum = &theMomentum;
      std::swap(theMomentum, rhs.theMomentum);
      std::swap(theFrozenMomentum, rhs.theFrozenMomentum);
      std::swap(thePosition, rhs.thePosition);
      std::swap(nCollisions, rhs.nCollisions);
      std::swap(nDecays, rhs.nDecays);
      std::swap(nSrcPair, rhs.nSrcPair),
      std::swap(thePotentialEnergy, rhs.thePotentialEnergy);
      // ID intentionally not swapped

#ifdef INCLXX_IN_GEANT4_MODE
      std::swap(theParentResonancePDGCode, rhs.theParentResonancePDGCode);
      std::swap(theParentResonanceID, rhs.theParentResonanceID);
#endif

      std::swap(theHelicity, rhs.theHelicity);
      std::swap(emissionTime, rhs.emissionTime);
      std::swap(outOfWell, rhs.outOfWell);
      std::swap(theSrcPartner, rhs.theSrcPartner);

      std::swap(theMass, rhs.theMass);
      std::swap(rpCorrelated, rhs.rpCorrelated);
      std::swap(uncorrelatedMomentum, rhs.uncorrelatedMomentum);
      
      std::swap(theParticleBias, rhs.theParticleBias);
      std::swap(theBiasCollisionVector, rhs.theBiasCollisionVector);

    }

  public:

    /** \brief Assignment operator
     *
     * Does not copy the particle ID.
     */
    Particle &operator=(const Particle &rhs) {
      Particle temporaryParticle(rhs);
      swap(temporaryParticle);
      return *this;
    }

    /**
     * Get the particle type.
     * @see G4INCL::ParticleType
     */
    G4INCL::ParticleType getType() const {
      return theType;
    };

    /// \brief Get the particle species
    virtual G4INCL::ParticleSpecies getSpecies() const {
      return ParticleSpecies(theType);
    };

    void setType(ParticleType t) {
      theType = t;
      switch(theType)
      {
        case DeltaPlusPlus:
          theA = 1;
          theZ = 2;
          theS = 0;
          break;
        case Proton:
        case DeltaPlus:
          theA = 1;
          theZ = 1;
          theS = 0;
          break;
        case Neutron:
        case DeltaZero:
          theA = 1;
          theZ = 0;
          theS = 0;
          break;
        case DeltaMinus:
          theA = 1;
          theZ = -1;
          theS = 0;
          break;
        case PiPlus:
          theA = 0;
          theZ = 1;
          theS = 0;
          break;
        case PiZero:
        case Eta:
        case Omega:
        case EtaPrime:
        case Photon:
          theA = 0;
          theZ = 0;
          theS = 0;
          break;
        case PiMinus:
          theA = 0;
          theZ = -1;
          theS = 0;
          break;
        case Lambda:
          theA = 1;
          theZ = 0;
          theS = -1;
          break;
        case SigmaPlus:
          theA = 1;
          theZ = 1;
          theS = -1;
          break;
        case SigmaZero:
          theA = 1;
          theZ = 0;
          theS = -1;
          break;
        case SigmaMinus:
          theA = 1;
          theZ = -1;
          theS = -1;
          break;         
        case antiProton:
          theA = -1;
          theZ = -1;
          theS = 0;
          break;         
        case XiMinus:
          theA = 1;
          theZ = -1;
          theS = -2;
          break;
        case XiZero:
          theA = 1;
          theZ = 0;
          theS = -2;
          break;      
        case antiNeutron:
          theA = -1;
          theZ = 0;
          theS = 0;
          break;
        case antiLambda:
          theA = -1;
          theZ = 0;
          theS = 1;
          break;
        case antiSigmaMinus:
          theA = -1;
          theZ = 1;
          theS = 1;
          break;
        case antiSigmaPlus:
          theA = -1;
          theZ = -1;
          theS = 1;
          break;
        case antiSigmaZero:
          theA = -1;
          theZ = 0;
          theS = 1;
          break;
        case antiXiMinus:
          theA = -1;
          theZ = 1;
          theS = 2;
          break;
        case antiXiZero:
          theA = -1;
          theZ = 0;
          theS = 2;
          break;         
        case KPlus:
          theA = 0;
          theZ = 1;
          theS = 1;
          break;
        case KZero:
          theA = 0;
          theZ = 0;
          theS = 1;
          break;
        case KZeroBar:
          theA = 0;
          theZ = 0;
          theS = -1;
          break;
        case KShort:
          theA = 0;
          theZ = 0;
//        theS should not be defined
          break;
        case KLong:
          theA = 0;
          theZ = 0;
//        theS should not be defined
          break;
        case KMinus:
          theA = 0;
          theZ = -1;
          theS = -1;
          break;
        case Composite:
         // INCL_ERROR("Trying to set particle type to Composite! Construct a Cluster object instead" << '\n');
          theA = 0;
          theZ = 0;
          theS = 0;
          break;       
        case antiComposite:
          theA = 0;
          theZ = 0;
          theS = 0;
          break;
        case UnknownParticle:
          theA = 0;
          theZ = 0;
          theS = 0;
          INCL_ERROR("Trying to set particle type to Unknown!" << '\n');
          break;
      }

      if( !isResonance() && t!=Composite && t!=antiComposite )
        setINCLMass();
    }

    /**
     * Is this a nucleon?
     */
    bool isNucleon() const {
      if(theType == G4INCL::Proton || theType == G4INCL::Neutron)
    return true;
      else
    return false;
    };

    ParticipantType getParticipantType() const {
      return theParticipantType;
    }

    void setParticipantType(ParticipantType const p) {
      theParticipantType = p;
    }

    bool isParticipant() const {
      return (theParticipantType==Participant);
    }

    bool isTargetSpectator() const {
      return (theParticipantType==TargetSpectator);
    }

    bool isProjectileSpectator() const {
      return (theParticipantType==ProjectileSpectator);
    }

    virtual void makeParticipant() {
      theParticipantType = Participant;
    }

    virtual void makeTargetSpectator() {
      theParticipantType = TargetSpectator;
    }

    virtual void makeProjectileSpectator() {
      theParticipantType = ProjectileSpectator;
    }

    /** \brief Is this a pion? */
    bool isPion() const { return (theType == PiPlus || theType == PiZero || theType == PiMinus); }

    /** \brief Is this an eta? */
    bool isEta() const { return (theType == Eta); }

    /** \brief Is this an omega? */
    bool isOmega() const { return (theType == Omega); }

    /** \brief Is this an etaprime? */
    bool isEtaPrime() const { return (theType == EtaPrime); }

    /** \brief Is this a photon? */
    bool isPhoton() const { return (theType == Photon); }

    /** \brief Is it a resonance? */
    inline bool isResonance() const { return isDelta(); }

    /** \brief Is it a Delta? */
    inline bool isDelta() const {
      return (theType==DeltaPlusPlus || theType==DeltaPlus ||
          theType==DeltaZero || theType==DeltaMinus); }
    
    /** \brief Is this a Sigma? */
    bool isSigma() const { return (theType == SigmaPlus || theType == SigmaZero || theType == SigmaMinus); }     
    
    /** \brief Is this a Kaon? */
    bool isKaon() const { return (theType == KPlus || theType == KZero); } 
    
    /** \brief Is this an antiKaon? */
    bool isAntiKaon() const { return (theType == KZeroBar || theType == KMinus); }
    
    /** \brief Is this a Lambda? */
    bool isLambda() const { return (theType == Lambda); }

    /** \brief Is this a Nucleon or a Lambda? */
    bool isNucleonorLambda() const { return (isNucleon() || isLambda()); }
    
    /** \brief Is this an Hyperon? */
    bool isHyperon() const { return (isLambda() || isSigma() ); } //|| isXi()
    
    /** \brief Is this a Meson? */
    bool isMeson() const { return (isPion() || isKaon() || isAntiKaon() || isEta() || isEtaPrime() || isOmega()); }
    
    /** \brief Is this a Baryon? */
    bool isBaryon() const { return (isNucleon() || isResonance() || isHyperon()); }
    
    /** \brief Is this a Strange? */
    bool isStrange() const { return (isKaon() || isAntiKaon() || isHyperon()); }
    
    /** \brief Is this a Xi? */
    bool isXi() const { return (theType == XiZero || theType == XiMinus); } 
    
    /** \brief Is this an antinucleon? */
    bool isAntiNucleon() const { return (theType == antiProton || theType == antiNeutron); } 
     
    /** \brief Is this an antiSigma? */
    bool isAntiSigma() const { return (theType == antiSigmaPlus || theType == antiSigmaZero || theType == antiSigmaMinus); }     
    
    /** \brief Is this an antiXi? */
    bool isAntiXi() const { return (theType == antiXiZero || theType == antiXiMinus); } 
    
    /** \brief Is this an antiLambda? */
    bool isAntiLambda() const { return (theType == antiLambda); }
    
    /** \brief Is this an antiHyperon? */
    bool isAntiHyperon() const { return (isAntiLambda() || isAntiSigma() || isAntiXi()); }
    
    /** \brief Is this an antiBaryon? */
    bool isAntiBaryon() const { return (isAntiNucleon() || isAntiHyperon()); }
    
    /** \brief Is this an antiNucleon or an antiLambda? */
    bool isAntiNucleonorAntiLambda() const { return (isAntiNucleon() || isAntiLambda()); }

    /** \brief Returns the baryon number. */
    int getA() const { return theA; }

    /** \brief Returns the charge number. */
    int getZ() const { return theZ; }
    
    /** \brief Returns the strangeness number. */
    int getS() const { return theS; }
    
    /** \brief Returns the strangeness number. */
    int getSrcPair() const { return nSrcPair; }

    double getBeta() const {
      const double P = theMomentum.mag();
      return P/theEnergy;
    }

    /**
     * Returns a three vector we can give to the boost() -method.
     *
     * In order to go to the particle rest frame you need to multiply
     * the boost vector by -1.0.
     */
    ThreeVector boostVector() const {
      return theMomentum / theEnergy;
    }

    /**
     * Boost the particle using a boost vector.
     *
     * Example (go to the particle rest frame):
     * particle->boost(particle->boostVector());
     */
    void boost(const ThreeVector &aBoostVector) {
      const double beta2 = aBoostVector.mag2();
      const double gamma = 1.0 / std::sqrt(1.0 - beta2);
      const double bp = theMomentum.dot(aBoostVector);
      const double alpha = (gamma*gamma)/(1.0 + gamma);

      theMomentum = theMomentum + aBoostVector * (alpha * bp - gamma * theEnergy);
      theEnergy = gamma * (theEnergy - bp);
    }

    /** \brief Lorentz-contract the particle position around some center
     *
     * Apply Lorentz contraction to the position component along the
     * direction of the boost vector.
     *
     * \param aBoostVector the boost vector (velocity) [c]
     * \param refPos the reference position
     */
    void lorentzContract(const ThreeVector &aBoostVector, const ThreeVector &refPos) {
      const double beta2 = aBoostVector.mag2();
      const double gamma = 1.0 / std::sqrt(1.0 - beta2);
      const ThreeVector theRelativePosition = thePosition - refPos;
      const ThreeVector transversePosition = theRelativePosition - aBoostVector * (theRelativePosition.dot(aBoostVector) / aBoostVector.mag2());
      const ThreeVector longitudinalPosition = theRelativePosition - transversePosition;

      thePosition = refPos + transversePosition + longitudinalPosition / gamma;
    }

    /** \brief Get the cached particle mass. */
    inline double getMass() const { return theMass; }

    /** \brief Get the INCL particle mass. */
    inline double getINCLMass() const {
      switch(theType) {
        case Proton:
        case Neutron:
        case PiPlus:
        case PiMinus:
        case PiZero:
        case Lambda:
        case SigmaPlus:
        case SigmaZero:
        case SigmaMinus:       
        case antiProton: 
        case XiZero:
        case XiMinus:
        case antiNeutron:
        case antiLambda:
        case antiSigmaPlus:
        case antiSigmaZero:
        case antiSigmaMinus:
        case antiXiZero:
        case antiXiMinus:     
        case KPlus:
        case KZero:
        case KZeroBar:
        case KShort:
        case KLong:
        case KMinus:
        case Eta:
        case Omega:
        case EtaPrime:
        case Photon:                       
          return ParticleTable::getINCLMass(theType);
          break;

        case DeltaPlusPlus:
        case DeltaPlus:
        case DeltaZero:
        case DeltaMinus:
          return theMass;
          break;

        case Composite:
          return ParticleTable::getINCLMass(theA,theZ,theS);
          break;
        case antiComposite:
          return ParticleTable::getINCLMass(-theA,-theZ,theS);
          break;
        default:
          INCL_ERROR("Particle::getINCLMass: Unknown particle type." << '\n');
          return 0.0;
          break;
      }
    }

    /** \brief Get the tabulated particle mass. */
    inline virtual double getTableMass() const {
      switch(theType) {
        case Proton:
        case Neutron:
        case PiPlus:
        case PiMinus:
        case PiZero:
        case Lambda:
        case SigmaPlus:
        case SigmaZero:
        case SigmaMinus:       
        case antiProton:      
        case XiZero:
        case XiMinus:  
        case antiNeutron:
        case antiLambda:
        case antiSigmaPlus:
        case antiSigmaZero:
        case antiSigmaMinus:
        case antiXiZero:
        case antiXiMinus:  
        case KPlus:
        case KZero:
        case KZeroBar:
        case KShort:
        case KLong:
        case KMinus:
        case Eta:
        case Omega:
        case EtaPrime:
        case Photon:  
          return ParticleTable::getTableParticleMass(theType);
          break;

        case DeltaPlusPlus:
        case DeltaPlus:
        case DeltaZero:
        case DeltaMinus:
          return theMass;
          break;

        case Composite:
          return ParticleTable::getTableMass(theA,theZ,theS);
          break;
        case antiComposite:
          return ParticleTable::getTableMass(-theA,-theZ,theS);
          break;
        default:
          INCL_ERROR("Particle::getTableMass: Unknown particle type." << '\n');
          return 0.0;
          break;
      }
    }

    /** \brief Get the real particle mass. */
    inline double getRealMass() const {
      switch(theType) {
        case Proton:
        case Neutron:
        case PiPlus:
        case PiMinus:
        case PiZero:
        case Lambda:
        case SigmaPlus:
        case SigmaZero:
        case SigmaMinus:       
        case antiProton: 
        case XiZero:
        case XiMinus: 
        case antiNeutron:
        case antiLambda:
        case antiSigmaPlus:
        case antiSigmaZero:
        case antiSigmaMinus:
        case antiXiZero:
        case antiXiMinus:    
        case KPlus:
        case KZero:
        case KZeroBar:
        case KShort:
        case KLong:
        case KMinus:
        case Eta:
        case Omega:
        case EtaPrime:
        case Photon:    
          return ParticleTable::getRealMass(theType);
          break;

        case DeltaPlusPlus:
        case DeltaPlus:
        case DeltaZero:
        case DeltaMinus:
          return theMass;
          break;

        case Composite:
          return ParticleTable::getRealMass(theA,theZ,theS);
          break;
        case antiComposite:
          return ParticleTable::getRealMass(-theA,-theZ,theS);
          break;
        default:
          INCL_ERROR("Particle::getRealMass: Unknown particle type." << '\n');
          return 0.0;
          break;
      }
    }

    /// \brief Set the mass of the Particle to its real mass
    void setRealMass() { setMass(getRealMass()); }

    /// \brief Set the mass of the Particle to its table mass
    void setTableMass() { setMass(getTableMass()); }

    /// \brief Set the mass of the Particle to its table mass
    void setINCLMass() { setMass(getINCLMass()); }

    /**\brief Computes correction on the emission Q-value
     *
     * Computes the correction that must be applied to INCL particles in
     * order to obtain the correct Q-value for particle emission from a given
     * nucleus. For absorption, the correction is obviously equal to minus
     * the value returned by this function.
     *
     * \param AParent the mass number of the emitting nucleus
     * \param ZParent the charge number of the emitting nucleus
     * \return the correction
     */
    double getEmissionQValueCorrection(const int AParent, const int ZParent) const {
      const int SParent = 0;
      const int ADaughter = AParent - theA;
      const int ZDaughter = ZParent - theZ;
      const int SDaughter = 0;

      // Note the minus sign here
      double theQValue;
      if(isCluster())
        theQValue = -ParticleTable::getTableQValue(theA, theZ, theS, ADaughter, ZDaughter, SDaughter);
      else {
        const double massTableParent = ParticleTable::getTableMass(AParent,ZParent,SParent);
        const double massTableDaughter = ParticleTable::getTableMass(ADaughter,ZDaughter,SDaughter);
        const double massTableParticle = getTableMass();
        theQValue = massTableParent - massTableDaughter - massTableParticle;
      }

      const double massINCLParent = ParticleTable::getINCLMass(AParent,ZParent,SParent);
      const double massINCLDaughter = ParticleTable::getINCLMass(ADaughter,ZDaughter,SDaughter);
      const double massINCLParticle = getINCLMass();

      // The rhs corresponds to the INCL Q-value
      return theQValue - (massINCLParent-massINCLDaughter-massINCLParticle);
    }

    /**\brief Computes correction on the transfer Q-value
     *
     * Computes the correction that must be applied to INCL particles in
     * order to obtain the correct Q-value for particle transfer from a given
     * nucleus to another.
     *
     * Assumes that the receving nucleus is INCL's target nucleus, with the
     * INCL separation energy.
     *
     * \param AFrom the mass number of the donating nucleus
     * \param ZFrom the charge number of the donating nucleus
     * \param ATo the mass number of the receiving nucleus
     * \param ZTo the charge number of the receiving nucleus
     * \return the correction
     */
    double getTransferQValueCorrection(const int AFrom, const int ZFrom, const int ATo, const int ZTo) const {
      const int SFrom = 0;
      const int STo = 0;
      const int AFromDaughter = AFrom - theA;
      const int ZFromDaughter = ZFrom - theZ;
      const int SFromDaughter = 0;
      const int AToDaughter = ATo + theA;
      const int ZToDaughter = ZTo + theZ;
      const int SToDaughter = 0;
      const double theQValue = ParticleTable::getTableQValue(AToDaughter,ZToDaughter,SToDaughter,AFromDaughter,ZFromDaughter,SFromDaughter,AFrom,ZFrom,SFrom);

      const double massINCLTo = ParticleTable::getINCLMass(ATo,ZTo,STo);
      const double massINCLToDaughter = ParticleTable::getINCLMass(AToDaughter,ZToDaughter,SToDaughter);
      /* Note that here we have to use the table mass in the INCL Q-value. We
       * cannot use theMass, because at this stage the particle is probably
       * still off-shell; and we cannot use getINCLMass(), because it leads to
       * violations of global energy conservation.
       */
      const double massINCLParticle = getTableMass();

      // The rhs corresponds to the INCL Q-value for particle absorption
      return theQValue - (massINCLToDaughter-massINCLTo-massINCLParticle);
    }

    /**\brief Computes correction on the emission Q-value for hypernuclei
     *
     * Computes the correction that must be applied to INCL particles in
     * order to obtain the correct Q-value for particle emission from a given
     * nucleus. For absorption, the correction is obviously equal to minus
     * the value returned by this function.
     *
     * \param AParent the mass number of the emitting nucleus
     * \param ZParent the charge number of the emitting nucleus
     * \param SParent the strangess number of the emitting nucleus
     * \return the correction
     */
    double getEmissionQValueCorrection(const int AParent, const int ZParent, const int SParent) const {
      const int ADaughter = AParent - theA;
      const int ZDaughter = ZParent - theZ;
      const int SDaughter = SParent - theS;

      // Note the minus sign here
      double theQValue;
      if(isCluster())
        theQValue = -ParticleTable::getTableQValue(theA, theZ, theS, ADaughter, ZDaughter, SDaughter);
      else {
        const double massTableParent = ParticleTable::getTableMass(AParent,ZParent,SParent);
        const double massTableDaughter = ParticleTable::getTableMass(ADaughter,ZDaughter,SDaughter);
        const double massTableParticle = getTableMass();
        theQValue = massTableParent - massTableDaughter - massTableParticle;
      }

      const double massINCLParent = ParticleTable::getINCLMass(AParent,ZParent,SParent);
      const double massINCLDaughter = ParticleTable::getINCLMass(ADaughter,ZDaughter,SDaughter);
      const double massINCLParticle = getINCLMass();

      // The rhs corresponds to the INCL Q-value
      return theQValue - (massINCLParent-massINCLDaughter-massINCLParticle);
    }

    /**\brief Computes correction on the transfer Q-value for hypernuclei
     *
     * Computes the correction that must be applied to INCL particles in
     * order to obtain the correct Q-value for particle transfer from a given
     * nucleus to another.
     *
     * Assumes that the receving nucleus is INCL's target nucleus, with the
     * INCL separation energy.
     *
     * \param AFrom the mass number of the donating nucleus
     * \param ZFrom the charge number of the donating nucleus
     * \param SFrom the strangess number of the donating nucleus
     * \param ATo the mass number of the receiving nucleus
     * \param ZTo the charge number of the receiving nucleus
     * \param STo the strangess number of the receiving nucleus
     * \return the correction
     */
    double getTransferQValueCorrection(const int AFrom, const int ZFrom, const int SFrom, const int ATo, const int ZTo , const int STo) const {
      const int AFromDaughter = AFrom - theA;
      const int ZFromDaughter = ZFrom - theZ;
      const int SFromDaughter = SFrom - theS;
      const int AToDaughter = ATo + theA;
      const int ZToDaughter = ZTo + theZ;
      const int SToDaughter = STo + theS;
      const double theQValue = ParticleTable::getTableQValue(AToDaughter,ZToDaughter,SFromDaughter,AFromDaughter,ZFromDaughter,SToDaughter,AFrom,ZFrom,SFrom);

      const double massINCLTo = ParticleTable::getINCLMass(ATo,ZTo,STo);
      const double massINCLToDaughter = ParticleTable::getINCLMass(AToDaughter,ZToDaughter,SToDaughter);
      /* Note that here we have to use the table mass in the INCL Q-value. We
       * cannot use theMass, because at this stage the particle is probably
       * still off-shell; and we cannot use getINCLMass(), because it leads to
       * violations of global energy conservation.
       */
      const double massINCLParticle = getTableMass();

      // The rhs corresponds to the INCL Q-value for particle absorption
      return theQValue - (massINCLToDaughter-massINCLTo-massINCLParticle);
    }



    /** \brief Get the the particle invariant mass.
     *
     * Uses the relativistic invariant
     * \f[ m = \sqrt{E^2 - {\vec p}^2}\f]
     **/
    double getInvariantMass() const {
      const double mass = std::pow(theEnergy, 2) - theMomentum.dot(theMomentum);
      if(mass < 0.0) {
        INCL_ERROR("E*E - p*p is negative." << '\n');
        return 0.0;
      } else {
        return std::sqrt(mass);
      }
    };

    /// \brief Get the particle kinetic energy.
    inline double getKineticEnergy() const { return theEnergy - theMass; }

    /// \brief Get the particle potential energy.
    inline double getPotentialEnergy() const { return thePotentialEnergy; }

    /// \brief Set the particle potential energy.
    inline void setPotentialEnergy(double v) { thePotentialEnergy = v; }

    /**
     * Get the energy of the particle in MeV.
     */
    double getEnergy() const
    {
      return theEnergy;
    };

    /**
     * Set the mass of the particle in MeV/c^2.
     */
    void setMass(double mass)
    {
      this->theMass = mass;
    }

    /**
     * Set the energy of the particle in MeV.
     */
    void setEnergy(double energy)
    {
      this->theEnergy = energy;
    };

    /**
     * Get the momentum vector.
     */
    const G4INCL::ThreeVector &getMomentum() const
    {
      return theMomentum;
    };

    /** Get the angular momentum w.r.t. the origin */
    virtual G4INCL::ThreeVector getAngularMomentum() const
    {
      return thePosition.vector(theMomentum);
    };

    /**
     * Set the momentum vector.
     */
    virtual void setMomentum(const G4INCL::ThreeVector &momentum)
    {
      this->theMomentum = momentum;
    };

    /**
     * Set the position vector.
     */
    const G4INCL::ThreeVector &getPosition() const
    {
      return thePosition;
    };

    virtual void setPosition(const G4INCL::ThreeVector &position)
    {
      this->thePosition = position;
    };

    double getHelicity() { return theHelicity; };
    void setHelicity(double h) { theHelicity = h; };

    void propagate(double step) {
      thePosition += ((*thePropagationMomentum)*(step/(*thePropagationEnergy)));
    };

    /** \brief Return the number of collisions undergone by the particle. **/
    int getNumberOfCollisions() const { return nCollisions; }

    /** \brief Set the number of collisions undergone by the particle. **/
    void setNumberOfCollisions(int n) { nCollisions = n; }

    /** \brief Increment the number of collisions undergone by the particle. **/
    void incrementNumberOfCollisions() { nCollisions++; }

    /** \brief Return the number of decays undergone by the particle. **/
    int getNumberOfDecays() const { return nDecays; }

    /** \brief Set the number of decays undergone by the particle. **/
    void setNumberOfDecays(int n) { nDecays = n; }

    /** \brief Increment the number of decays undergone by the particle. **/
    void incrementNumberOfDecays() { nDecays++; }
    
    /** \brief Set the number of srcpairs. **/
    void setNumberOfSrcPair(int n) { nSrcPair = n; } 

    /** \brief Mark the particle as out of its potential well
     *
     * This flag is used to control pions created outside their potential well
     * in delta decay. The pion potential checks it and returns zero if it is
     * true (necessary in order to correctly enforce energy conservation). The
     * Nucleus::applyFinalState() method uses it to determine whether new
     * avatars should be generated for the particle.
     */
    void setOutOfWell() { outOfWell = true; }

    /// \brief Check if the particle is out of its potential well
    bool isOutOfWell() const { return outOfWell; }
    
    /// \brief Set and reset src partner 
    void setSrcPartner() { theSrcPartner = true; }  
    void resetSrcPartner() { theSrcPartner = false; nSrcPair=0; }

    /// \brief Check if the particle is a src partner    
    bool isSrcPartner() const { return theSrcPartner; }

    void setEmissionTime(double t) { emissionTime = t; }
    double getEmissionTime() { return emissionTime; };

    /** \brief Transverse component of the position w.r.t. the momentum. */
    ThreeVector getTransversePosition() const {
      return thePosition - getLongitudinalPosition();
    }

    /** \brief Longitudinal component of the position w.r.t. the momentum. */
    ThreeVector getLongitudinalPosition() const {
      return *thePropagationMomentum * (thePosition.dot(*thePropagationMomentum)/thePropagationMomentum->mag2());
    }

    /** \brief Rescale the momentum to match the total energy. */
    const ThreeVector &adjustMomentumFromEnergy();

    /** \brief Recompute the energy to match the momentum. */
    double adjustEnergyFromMomentum();

    bool isCluster() const {
      return ((theType == Composite || theType == antiComposite));
    }

    /// \brief Set the frozen particle momentum
    void setFrozenMomentum(const ThreeVector &momentum) { theFrozenMomentum = momentum; }

    /// \brief Set the frozen particle momentum
    void setFrozenEnergy(const double energy) { theFrozenEnergy = energy; }

    /// \brief Get the frozen particle momentum
    ThreeVector getFrozenMomentum() const { return theFrozenMomentum; }

    /// \brief Get the frozen particle momentum
    double getFrozenEnergy() const { return theFrozenEnergy; }

    /// \brief Get the propagation velocity of the particle
    ThreeVector getPropagationVelocity() const { return (*thePropagationMomentum)/(*thePropagationEnergy); }

    /** \brief Freeze particle propagation
     *
     * Make the particle use theFrozenMomentum and theFrozenEnergy for
     * propagation. The normal state can be restored by calling the
     * thawPropagation() method.
     */
    void freezePropagation() {
      thePropagationMomentum = &theFrozenMomentum;
      thePropagationEnergy = &theFrozenEnergy;
    }

    /** \brief Unfreeze particle propagation
     *
     * Make the particle use theMomentum and theEnergy for propagation. Call
     * this method to restore the normal propagation if the
     * freezePropagation() method has been called.
     */
    void thawPropagation() {
      thePropagationMomentum = &theMomentum;
      thePropagationEnergy = &theEnergy;
    }

    /** \brief Rotate the particle position and momentum
     *
     * \param angle the rotation angle
     * \param axis a unit vector representing the rotation axis
     */
    virtual void rotatePositionAndMomentum(const double angle, const ThreeVector &axis) {
      rotatePosition(angle, axis);
      rotateMomentum(angle, axis);
    }

    /** \brief Rotate the particle position
     *
     * \param angle the rotation angle
     * \param axis a unit vector representing the rotation axis
     */
    virtual void rotatePosition(const double angle, const ThreeVector &axis) {
      thePosition.rotate(angle, axis);
    }

    /** \brief Rotate the particle momentum
     *
     * \param angle the rotation angle
     * \param axis a unit vector representing the rotation axis
     */
    virtual void rotateMomentum(const double angle, const ThreeVector &axis) {
      theMomentum.rotate(angle, axis);
      theFrozenMomentum.rotate(angle, axis);
    }

    std::string print() const {
      std::stringstream ss;
      ss << "Particle (ID = " << ID << ") type = ";
      ss << ParticleTable::getName(theType);
      ss << ", SRC pair = " << nSrcPair;
      ss << ", Potential energy = " << thePotentialEnergy;
      ss << '\n'
        << "   energy = " << theEnergy << '\n'
        << "   momentum = "
        << theMomentum.print()
        << '\n'
        << "   position = "
        << thePosition.print()
        << '\n';
      return ss.str();
    };

    std::string dump() const {
      std::stringstream ss;
      ss << "(particle " << ID << " ";
      ss << ParticleTable::getName(theType);
      ss << nSrcPair << " ";
      ss << '\n'
        << thePosition.dump()
        << '\n'
        << theMomentum.dump()
        << '\n'
        << theEnergy << ")" << '\n';
      return ss.str();
    };

    long getID() const { return ID; };

    /**
     * Return a NULL pointer
     */
    ParticleList const *getParticles() const {
      INCL_WARN("Particle::getParticles() method was called on a Particle object" << '\n');
      return 0;
    }

    /** \brief Return the reflection momentum
     *
     * The reflection momentum is used by calls to getSurfaceRadius to compute
     * the radius of the sphere where the nucleon moves. It is necessary to
     * introduce fuzzy r-p correlations.
     */
    double getReflectionMomentum() const {
      if(rpCorrelated)
        return theMomentum.mag();
      else
        return uncorrelatedMomentum;
    }

    /// \brief Set the uncorrelated momentum
    void setUncorrelatedMomentum(const double p) { uncorrelatedMomentum = p; }

    /// \brief Make the particle follow a strict r-p correlation
    void rpCorrelate() { rpCorrelated = true; }

    /// \brief Make the particle not follow a strict r-p correlation
    void rpDecorrelate() { rpCorrelated = false; }

    /// \brief Get the cosine of the angle between position and momentum
    double getCosRPAngle() const {
      const double norm = thePosition.mag2()*thePropagationMomentum->mag2();
      if(norm>0.)
        return thePosition.dot(*thePropagationMomentum) / std::sqrt(norm);
      else
        return 1.;
    }

    /// \brief General bias vector function
    static double getTotalBias();
    static void setINCLBiasVector(std::vector<double> NewVector);
    static void FillINCLBiasVector(double newBias);
    static double getBiasFromVector(std::vector<int> VectorBias);

    static std::vector<int> MergeVectorBias(Particle const * const p1, Particle const * const p2);
    static std::vector<int> MergeVectorBias(std::vector<int> p1, Particle const * const p2);

    /// \brief Get the particle bias.
    double getParticleBias() const { return theParticleBias; };

    /// \brief Set the particle bias.
    void setParticleBias(double ParticleBias) { this->theParticleBias = ParticleBias; }

    /// \brief Get the vector list of biased vertices on the particle path.
    std::vector<int> getBiasCollisionVector() const { return theBiasCollisionVector; }

    /// \brief Set the vector list of biased vertices on the particle path.
    void setBiasCollisionVector(std::vector<int> BiasCollisionVector) {
	  this->theBiasCollisionVector = BiasCollisionVector;
	  this->setParticleBias(Particle::getBiasFromVector(BiasCollisionVector));
	  }
    
    /** \brief Number of Kaon inside de nucleus
     * 
     * Put in the Particle class in order to calculate the
     * "correct" mass of composit particle.
     * 
     */
     
    int getNumberOfKaon() const { return theNKaon; };
    void setNumberOfKaon(const int NK) { theNKaon = NK; }

#ifdef INCLXX_IN_GEANT4_MODE
    G4int getParentResonancePDGCode() const { return theParentResonancePDGCode; };
    void setParentResonancePDGCode(const G4int parentPDGCode) { theParentResonancePDGCode = parentPDGCode; };    
    G4int getParentResonanceID() const { return theParentResonanceID; };
    void setParentResonanceID(const G4int parentID) { theParentResonanceID = parentID; };
#endif
  
  public:
    /** \brief Time ordered vector of all bias applied
     * 
     * /!\ Caution /!\
     * methods Assotiated to G4VectorCache<T> are:
     * Push_back(…),
     * operator[],
     * Begin(),
     * End(),
     * Clear(),
     * Size() and 
     * Pop_back()
     * 
     */
#ifdef INCLXX_IN_GEANT4_MODE
      static std::vector<double> INCLBiasVector;
      //static G4VectorCache<double> INCLBiasVector;
#else
      static G4ThreadLocal std::vector<double> INCLBiasVector;
      //static G4VectorCache<double> INCLBiasVector;
#endif
    static G4ThreadLocal int nextBiasedCollisionID;
    
  protected:
    int theZ, theA, theS;
    ParticipantType theParticipantType;
    G4INCL::ParticleType theType;
    double theEnergy;
    double *thePropagationEnergy;
    double theFrozenEnergy;
    G4INCL::ThreeVector theMomentum;
    G4INCL::ThreeVector *thePropagationMomentum;
    G4INCL::ThreeVector theFrozenMomentum;
    G4INCL::ThreeVector thePosition;
    int nCollisions;
    int nDecays;
    int nSrcPair;
    double thePotentialEnergy;
    long ID;

    bool rpCorrelated;
    double uncorrelatedMomentum;
    
    double theParticleBias;
    /// \brief The number of Kaons inside the nucleus (update during the cascade)
    int theNKaon;

#ifdef INCLXX_IN_GEANT4_MODE
    G4int theParentResonancePDGCode;
    G4int theParentResonanceID;
#endif

  private:
    double theHelicity;
    double emissionTime;
    bool outOfWell;
    bool theSrcPartner;
    
    /// \brief Time ordered vector of all biased vertices on the particle path
    std::vector<int> theBiasCollisionVector;

    double theMass;
    static G4ThreadLocal long nextID;

    INCL_DECLARE_ALLOCATION_POOL(Particle)
  };
}

#endif /* PARTICLE_HH_ */
