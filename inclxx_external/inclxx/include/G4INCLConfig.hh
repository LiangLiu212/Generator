#ifndef G4INCLConfig_hh
#define G4INCLConfig_hh 1

#include "G4INCLParticleSpecies.hh"
#include "G4INCLConfigEnums.hh"
#include "G4INCLRandomSeedVector.hh"
#include <iostream>
#include <string>
#include <sstream>
#include <cassert>

class ConfigParser;

namespace G4INCL {

  /**
   * The INCL configuration object
   *
   * The Config object keeps track of various INCL physics options
   * (e.g. which Pauli blocking scheme to use, whether to use local
   * energy option or not, etc.
   */
  class Config {
  public:
    /// \brief Default constructor
    Config();

    /// \brief Default destructor
    ~Config();

    /// \brief Initialise the members
    void init();

    /// \brief Return a summary of the run configuration.
    std::string summary();

    /// \brief Get the verbosity.
    int getVerbosity() const { return verbosity; }

    /// \brief Get the run title.
    std::string const &getCalculationTitle() const { return title; }

    /// \brief Get the output file root.
    std::string const &getOutputFileRoot() const { return outputFileRoot; }

    /// \brief Get the number of shots.
    int getNumberOfShots() const { return nShots; }

    /// \brief Natural targets.
    bool isNaturalTarget() const { return naturalTarget; }

    /** \brief Get the target mass number.
     *
     * Note that A==0 means natural target. You should first check the
     * isNaturalTarget() method.
     */
    int getTargetA() const { return targetSpecies.theA; }

    /// \brief Get the target charge number.
    int getTargetZ() const { return targetSpecies.theZ; }

    /// \brief Get the target strangess number.
    int getTargetS() const { return targetSpecies.theS; }

    /// \brief Set target mass number
    void setTargetA(int A) { targetSpecies.theA = A; }

    /// \brief Set target charge number
    void setTargetZ(int Z) { targetSpecies.theZ = Z; }

    /// \brief Set target strangess number
    void setTargetS(int S) { targetSpecies.theS = S; }

    /// \brief Get the projectile type
    ParticleType getProjectileType() const { return projectileSpecies.theType; }

    /// \brief Get the projectile species
    ParticleSpecies getProjectileSpecies() const { return projectileSpecies; }

    /// \brief Set the projectile species
    void setProjectileSpecies(ParticleSpecies const &pars) { projectileSpecies=pars; }

    /// \brief Get the projectile kinetic energy.
    double getProjectileKineticEnergy() const { return projectileKineticEnergy; }

    /// \brief Set the projectile kinetic energy.
    void setProjectileKineticEnergy(double const kinE) { projectileKineticEnergy=kinE; }

    /// \brief Get the number of the verbose event.
    int getVerboseEvent() const { return verboseEvent; }

    /// \brief Get the INCL version ID.
    static std::string const getVersionID();

    /// \brief Get the INCL version hash.
    static std::string const getVersionHash();

    /// \brief Get the INCL version string.
    static std::string const getVersionString() {
      std::stringstream ss;
      ss << getVersionID() << "-" << getVersionHash();
      return ss.str();
    }

    /// \brief Get the seeds for the random-number generator.
    Random::SeedVector getRandomSeeds() const {
      return randomSeedVector;
    }

    /// \brief Get the Pauli-blocking algorithm.
    PauliType getPauliType() const { return pauliType; }
    void setPauliType(PauliType type) {pauliType = type;}
    void setPauliString(std::string str) {pauliString = str;}

    /// \brief Do we want CDPP?
    bool getCDPP() const { return CDPP; }

    /// \brief Get the Coulomb-distortion algorithm.
    CoulombType getCoulombType() const { return coulombType; }

    /// \brief Set the Coulomb-distortion algorithm.
    void setCoulombType(CoulombType const c) { coulombType = c; }

    /// \brief Get the type of the potential for nucleons.
    PotentialType getPotentialType() const { return potentialType; }

    /// \brief Set the type of the potential for nucleons.
    void setPotentialType(PotentialType type) { potentialType = type; }

    /// \brief Do we want the pion potential?
    bool getPionPotential() const { return pionPotential; }

    /// \brief Set the type of the potential for nucleons.
    void setPionPotential(const bool pionPot) { pionPotential = pionPot; }

    /// \brief Get the type of local energy for N-N avatars.
    LocalEnergyType getLocalEnergyBBType() const { return localEnergyBBType; }

    /// \brief Set the type of local energy for N-N avatars.
    void setLocalEnergyBBType(const LocalEnergyType t) { localEnergyBBType=t; }

    /// \brief Get the type of local energy for pi-N and decay avatars.
    LocalEnergyType getLocalEnergyPiType() const { return localEnergyPiType; }

    /// \brief Set the type of local energy for N-N avatars.
    void setLocalEnergyPiType(const LocalEnergyType t) { localEnergyPiType=t; }

    /// \brief Get the log file name.
    std::string const &getLogFileName() const { return logFileName; }

    /// \brief Get the de-excitation model.
    void setDeExcitationType(DeExcitationType t) { deExcitationType = t; }
    DeExcitationType getDeExcitationType() const { return deExcitationType; }

    /// \brief Get the de-excitation string.
    std::string getDeExcitationString() const { return deExcitationString; }

    /// \brief Get the clustering algorithm.
    ClusterAlgorithmType getClusterAlgorithm() const { return clusterAlgorithmType; }

    /// \brief Set the clustering algorithm.
    void setClusterAlgorithm(ClusterAlgorithmType const c) { clusterAlgorithmType = c; }
    void setClusterAlgorithmString(std::string const str) { clusterAlgorithmString = str; }

    /// \brief Get the maximum mass for production of clusters.
    int getClusterMaxMass() const { return clusterMaxMass; }

    /// \brief Set the maximum mass for production of clusters.
    void setClusterMaxMass(const int clm){ clusterMaxMass=clm; }

    /// \brief Get back-to-spectator
    bool getBackToSpectator() const { return backToSpectator; }

    /// \brief Set back-to-spectator
    void setBackToSpectator(const bool b) { backToSpectator = b; }

    /// \brief Whether to use real masses
    bool getUseRealMasses() const { return useRealMasses; }

    /// \brief Set whether to use real masses
    void setUseRealMasses(bool use) { useRealMasses = use; }

    /// \brief Set the INCLXX datafile path
    void setINCLXXDataFilePath(std::string const &path) { INCLXXDataFilePath=path; }
    
    /// \brief Set the ABLAXX datafile path
#ifdef INCL_DEEXCITATION_ABLAXX
    void setABLAXXDataFilePath(std::string const &path) { ablaxxDataFilePath=path; }
    double getFissionDispCoeff() const { return fissdisscoeff; }
    double getLevDensAv() const { return levdensav; }
    double getLevDensAs() const { return levdensas; }
    double getLevDensAk() const { return levdensak; }
    double getTempFreezeOut() const { return tfreezeout; }
    float getEmissionBarrierFactorH2() const { return barrierfactorh2; }
    float getEmissionBarrierFactorH3() const { return barrierfactorh3; }
    float getEmissionBarrierFactorHe3() const { return barrierfactorhe3; }
    float getEmissionBarrierFactorHe4() const { return barrierfactorhe4; }
    float getEmissionBarrierFactorHe6() const { return barrierfactorhe6; }
#endif    

    std::string const &getINCLXXDataFilePath() const {
      return INCLXXDataFilePath;
    }

#ifdef INCL_DEEXCITATION_ABLAXX
    std::string const &getABLAXXDataFilePath() const {
      return ablaxxDataFilePath;
    }
#endif

#ifdef INCL_DEEXCITATION_ABLA07
    void setcABLA07DataFilePath(std::string const &path) { abla07DataFilePath=path; }
    std::string const &getABLA07DataFilePath() const {
      return abla07DataFilePath;
    }
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
    void setGEMINIXXDataFilePath(std::string const &path) { geminixxDataFilePath=path; }
    std::string const &getGEMINIXXDataFilePath() const {
      return geminixxDataFilePath;
    }
#endif

    double getImpactParameter() const { return impactParameter; }

    /// \brief Get the separation-energy type
    SeparationEnergyType getSeparationEnergyType() const { return separationEnergyType; }

    /// \brief Get the Fermi-momentum type
    FermiMomentumType getFermiMomentumType() const { return fermiMomentumType; }

    /// \brief Set the Fermi-momentum type
    void setFermiMomentumType(FermiMomentumType const f) { fermiMomentumType=f; }

    /// \brief Get the Fermi momentum
    double getFermiMomentum() const { return fermiMomentum; }

    /// \brief Set the Fermi momentum
    void setFermiMomentum(const double p) { fermiMomentum = p; }

    double getCutNN() const { return cutNN; }

#ifdef INCL_ROOT_USE
    std::string const &getROOTSelectionString() const {
      return rootSelectionString;
    }
#endif

#ifdef INCL_DEEXCITATION_FERMI_BREAKUP
    int getMaxMassFermiBreakUp() const {
      return maxMassFermiBreakUp;
    }

    int getMaxChargeFermiBreakUp() const {
      return maxChargeFermiBreakUp;
    }
#endif

    /// \brief Get the r-p correlation coefficient
    double getRPCorrelationCoefficient(const ParticleType t) const {
      assert(t==Proton || t==Neutron);
      return ((t==Proton) ? rpCorrelationCoefficientProton : rpCorrelationCoefficientNeutron);
    }

    /// \brief Set the r-p correlation coefficient
    void setRPCorrelationCoefficient(const ParticleType t, const double corrCoeff) {
      assert(t==Proton || t==Neutron);
      if(t==Proton)
        rpCorrelationCoefficientProton=corrCoeff;
      else
        rpCorrelationCoefficientNeutron=corrCoeff;
    }

    /// \brief Set the r-p correlation coefficient
    void setRPCorrelationCoefficient(const double corrCoeff) {
      setRPCorrelationCoefficient(Proton,corrCoeff);
      setRPCorrelationCoefficient(Neutron,corrCoeff);
    }

    /// \brief Get the neutron-skin thickness
    double getNeutronSkin() const { return neutronSkin; }

    /// \brief Set the neutron-skin thickness
    void setNeutronSkin(const double d) { neutronSkin=d; }

    /// \brief Get the neutron-halo size
    double getNeutronHalo() const { return neutronHalo; }

    /// \brief Set the neutron-halo size
    void setNeutronHalo(const double d) { neutronHalo=d; }

    /// \brief True if we should use refraction
    bool getRefraction() const { return refraction; }

    /// \brief Set the refraction variable
    void setRefraction(const bool r) { refraction = r; }

    /// \brief Get the RNG type
    RNGType getRNGType() const { return rngType; }

    /// \brief Set the RNG type
    void setRNGType(RNGType const r) { rngType=r; }

    /// \brief Get the phase-space-generator type
    PhaseSpaceGeneratorType getPhaseSpaceGeneratorType() const { return phaseSpaceGeneratorType; }

    /// \brief Set the phase-space-generator type
    void setPhaseSpaceGeneratorType(PhaseSpaceGeneratorType const p) { phaseSpaceGeneratorType=p; }

    /// \brief Get the cascade-action type
    CascadeActionType getCascadeActionType() const { return cascadeActionType; }

    /// \brief Set the cascade-action type
    void setCascadeActionType(CascadeActionType const c) { cascadeActionType=c; }

    /// \brief Get the autosave frequency
    unsigned int getAutosaveFrequency() const { return autosaveFrequency; }

    /// \brief Set the autosave frequency
    void setAutosaveFrequency(const unsigned int f) { autosaveFrequency=f; }

    /// \brief Get the Cross Section type
    CrossSectionsType getCrossSectionsType() const { return crossSectionsType; }

    /// \brief Get the maximum number of pions for multipion collisions
    int getMaxNumberMultipions() const { return maxNumberMultipions; }

    /// \brief Set the maximum number of pions for multipion collisions
    void setMaxNumberMultipions(const int n) { maxNumberMultipions=n; }

    /// \brief Set the Cross Section type
    void setCrossSectionsType(CrossSectionsType const c) { crossSectionsType=c; }

    /// \brief Get the hadronization time
    double getHadronizationTime() const { return hadronizationTime; }

    /// \brief Set the hadronization time
    void setHadronizationTime(const double t) { hadronizationTime=t; }

#ifdef INCL_ROOT_USE
    bool getConciseROOTTree() const { return conciseROOTTree; }
#endif

    bool getInverseKinematics() const { return inverseKinematics; }
    
    bool getsrcPairConfig() const { return srcPairCorrelations; }
    
    float getsrcPairDist() const { return srcPairDistance; }

    /// \brief Get the decay time threshold time
    double getDecayTimeThreshold() const { return decayTimeThreshold; }

    /// \brief Set decay time threshold time
    void setDecayTimeThreshold(const double t) { decayTimeThreshold=t; }

    /// \brief Get the bias
    double getBias() const { return bias; }

    /// \brief Get the pbar at rest annihilation threshold
    double getAtrestThreshold() const { return atrestThreshold; }

    /// \brief Set the pbar at rest annihilation threshold
    void setAtrestThreshold(const double t) { atrestThreshold=t;}

    /// \brief Get the nbar at rest annihilation threshold
    double getnbAtrestThreshold() const {return nbatrestThreshold;} 
    
    /// \brief Set the nbar at rest annihilation threshold
    void setnbAtrestThreshold(const double t){ nbatrestThreshold=t;}

    /// \brief Get the dbar at rest annihilation threshold
    double getdbAtrestThreshold() const {return dbatrestThreshold;}

    /// \brief Set the dbar at rest annihilation threshold
    void setdbAtrestThreshold(const double t){ dbatrestThreshold=t;} 


  private:

    int verbosity;
    std::string inputFileName;
    std::string title;
    std::string outputFileRoot;
    std::string fileSuffix;
    std::string logFileName;

    int nShots;

    std::string targetString;
    ParticleSpecies targetSpecies;
    bool naturalTarget;

    std::string projectileString;
    ParticleSpecies projectileSpecies;
    double projectileKineticEnergy;

    int verboseEvent;

    std::string randomSeeds;
    Random::SeedVector randomSeedVector;

    std::string pauliString;
    PauliType pauliType;
    bool CDPP;

    std::string coulombString;
    CoulombType coulombType;

    std::string potentialString;
    PotentialType potentialType;
    bool pionPotential;

    std::string localEnergyBBString;
    LocalEnergyType localEnergyBBType;

    std::string localEnergyPiString;
    LocalEnergyType localEnergyPiType;

    std::string deExcitationModelList;
    std::string deExcitationOptionDescription;
    std::string deExcitationString;
    DeExcitationType deExcitationType;
#ifdef INCL_DEEXCITATION_ABLAXX
    std::string ablaxxDataFilePath;
    double fissdisscoeff;
    double levdensav;
    double levdensas;
    double levdensak;
    double tfreezeout;
    float barrierfactorh2;
    float barrierfactorh3;
    float barrierfactorhe3;
    float barrierfactorhe4;
    float barrierfactorhe6;
#endif
#ifdef INCL_DEEXCITATION_ABLA07
    std::string abla07DataFilePath;
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
    std::string geminixxDataFilePath;
#endif
    std::string INCLXXDataFilePath;

    std::string clusterAlgorithmString;
    ClusterAlgorithmType clusterAlgorithmType;

    int clusterMaxMass;

    bool backToSpectator;

    bool useRealMasses;

    double impactParameter;

    std::string separationEnergyString;
    SeparationEnergyType separationEnergyType;

    std::string fermiMomentumString;
    FermiMomentumType fermiMomentumType;

    double fermiMomentum;

    double cutNN;

    //bool ann;
    
    double bias;

    double atrestThreshold;

    double nbatrestThreshold;

    double dbatrestThreshold;

#ifdef INCL_ROOT_USE
    std::string rootSelectionString;
#endif

#ifdef INCL_DEEXCITATION_FERMI_BREAKUP
    int maxMassFermiBreakUp;
    int maxChargeFermiBreakUp;
#endif

    double rpCorrelationCoefficient;
    double rpCorrelationCoefficientProton;
    double rpCorrelationCoefficientNeutron;

    double neutronSkin;
    double neutronHalo;

    bool refraction;

    std::string randomNumberGenerator;
    RNGType rngType;

    std::string phaseSpaceGenerator;
    PhaseSpaceGeneratorType phaseSpaceGeneratorType;

    unsigned int autosaveFrequency;

    std::string crossSectionsString;
    CrossSectionsType crossSectionsType;
    int maxNumberMultipions;

    std::string cascadeAction;
    CascadeActionType cascadeActionType;

    double hadronizationTime;

#ifdef INCL_ROOT_USE
    bool conciseROOTTree;
#endif

    bool inverseKinematics;
    
    bool srcPairCorrelations;
    
    float srcPairDistance;

    double decayTimeThreshold;

    friend class ::ConfigParser;
  };

}

#endif
