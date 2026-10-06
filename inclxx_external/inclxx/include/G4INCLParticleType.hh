/*
 * G4INCLParticleType.hh
 *
 *  \date Feb 4, 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLPARTICLETYPE_HH_
#define G4INCLPARTICLETYPE_HH_

namespace G4INCL {

  enum ParticleType {
    Proton = 0,
    Neutron,
    PiPlus,
    PiMinus,
    PiZero,
    DeltaPlusPlus,
    DeltaPlus,
    DeltaZero,
    DeltaMinus,
    Composite,
    Eta,
    Omega,
    EtaPrime,
    Photon,
    Lambda,
    SigmaPlus,
    SigmaZero,
    SigmaMinus,  
    antiProton,
    XiMinus,
    XiZero,
    antiNeutron,
    antiLambda,
    antiSigmaPlus,
    antiSigmaZero,
    antiSigmaMinus,
    antiXiMinus,
    antiXiZero,
    KPlus,
    KZero,
    KZeroBar,
    KMinus,
    KShort,
    KLong, 
    antiComposite,
    // WARNING: if you add more particle types, you MUST add them BEFORE the
    // UnknownParticle type! This is because UnknownParticle is used as a
    // counter of the number of available particle types.
    UnknownParticle
  };

  enum ParticipantType {
    TargetSpectator,
    ProjectileSpectator,
    Participant
  };

}

#endif /* G4INCLPARTICLETYPE_HH_ */
