#ifndef G4INCLConfigEnums_hh
#define G4INCLConfigEnums_hh

namespace G4INCL {

  // Enumerator for Pauli-blocking algorithms
  enum PauliType {
    StatisticalPauli,
    StrictPauli,
    StrictStatisticalPauli,
    GlobalPauli,
    NoPauli
  };

  // Enumerator for Coulomb-distortion algorithms
  enum CoulombType {
    NonRelativisticCoulomb,
    NoCoulomb
  };

  // Enumerator for potential types
  enum PotentialType {
    IsospinEnergySmoothPotential,
    IsospinEnergyPotential,
    IsospinPotential,
    ConstantPotential
  };

  // Enumerator for local-energy types
  enum LocalEnergyType {
    AlwaysLocalEnergy,
    FirstCollisionLocalEnergy,
    NeverLocalEnergy
  };

  // Enumerator for de-excitation types
  enum DeExcitationType {
    DeExcitationNone
#ifdef INCL_DEEXCITATION_ABLAXX
    , DeExcitationABLAXX
#endif
#ifdef INCL_DEEXCITATION_ABLA07
    , DeExcitationABLA07
#endif
#ifdef INCL_DEEXCITATION_SMM
    , DeExcitationSMM
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
    , DeExcitationGEMINIXX
#endif
  };

  // Enumerator for cluster-algorithm types
  enum ClusterAlgorithmType {
    IntercomparisonClusterAlgorithm,
    NoClusterAlgorithm
  };

  // Enumerator for separation-energy types
  enum SeparationEnergyType {
    INCLSeparationEnergy,
    RealSeparationEnergy,
    RealForLightSeparationEnergy
  };

  // Enumerator for Fermi-momentum types
  enum FermiMomentumType {
    ConstantFermiMomentum,
    ConstantLightFermiMomentum,
    MassDependentFermiMomentum
  };

  // Enumerator for RNG
  enum RNGType {
    RanecuType,
    Ranecu3Type
  };

  // Enumerator for Cross-Section parametrizations
  enum CrossSectionsType {
    INCL46CrossSections,
    MultiPionsCrossSections,
	TruncatedMultiPionsCrossSections,
	MultiPionsAndResonancesCrossSections,
	StrangenessCrossSections,
  AntiparticlesCrossSections
  };

  // Enumerator for phase-space generator
  enum PhaseSpaceGeneratorType {
    KopylovType,
    RauboldLynchType
  };

  // Enumerator for cascade actions
  enum CascadeActionType {
    DefaultActionType,
    AvatarDumpActionType
  };

}

#endif
