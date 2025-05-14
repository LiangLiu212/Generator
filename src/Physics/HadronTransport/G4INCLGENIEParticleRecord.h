#ifndef G4INCLGENIEPARTICLERECORD_H_
#define G4INCLGENIEPARTICLERECORD_H_

#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "G4INCLThreeVector.hh"
#include "G4INCLParticleType.hh"

namespace G4INCL {
  /*
   * Convert GENIE Event Record to INCL Style
   *
   */

  enum GENIERecordCode {
    kUnknown = -1,
    kProbe,
    kTarget,
    kHitNucleon,
    kRemnant,
    kFinalStateLepton,
    kHadron
  };

  class GENIEParticleRecord {
    public:
      GENIEParticleRecord(genie::GHepParticle* p, int scType, GENIERecordCode recordCode);
      ~GENIEParticleRecord();
      // Basic properties
      int           Pdg            (void) const { return  fPdgCode;            }
      int           ID             (void) const { return  fID;                 }
      int           Status         (void) const { return  fStatus;             }
      int           FirstMother    (void) const { return  fFirstMother;        }
      int           LastMother     (void) const { return  fLastMother;         }
      int           FirstDaughter  (void) const { return  fFirstDaughter;      }
      int           LastDaughter   (void) const { return  fLastDaughter;       }
      int           ScatteringType (void) const { return  fScatteringType;     }
      double        Mass           (void) const { return  fMass;               }
      ThreeVector   P3             (void) const { return  fP3;                 }
      ThreeVector   X3             (void) const { return  fX3;                 }
      ParticleType  Type           (void) const { return  fPType;              } 
      void          setID          (int id);

      void          setMomentum    (ThreeVector &mom) {  fP3 = mom;            }
      void          setMass        (double       m)   {  fMass = m;            }


      GENIERecordCode RecordCode   (void) const { return  fGenieRecordCode;    }

    private:
      int           fPdgCode;        ///< particle PDG code
      int           fID;             ///< particle INCL ID
      int           fStatus;         ///< particle status (Convert GENIE status code to INCL FSI code 0: no rescattering, 1: do rescattering)
      int           fFirstMother;    ///< first mother idx
      int           fLastMother;     ///< last mother idx
      int           fFirstDaughter;  ///< first daughter idx
      int           fLastDaughter;   ///< last daughter idx
      int           fScatteringType; ///< scattering type  (QEL, RES, DIS, ...)
      double        fMass;           ///< invariant mass (MeV)
      ThreeVector   fP3;             ///< momentum 3-vector (MeV)
      ThreeVector   fX3;             ///< position 3-vector (in the target nucleus coordinate system / x,y,z in fm / t from the moment of the primary interaction in ys(yocto second = 10^-24 s)
      ParticleType fPType;           ///< INCL particle type
      GENIERecordCode fGenieRecordCode;



      ParticleType PDG_to_INCLType(int pdg) const {
	switch(pdg){
	  case 2212: return G4INCL::Proton;
	  case 2112: return G4INCL::Neutron;
	  case 211: return G4INCL::PiPlus;
	  case -211: return G4INCL::PiMinus;
	  case 111: return G4INCL::PiZero;
	  case 2224: return G4INCL::DeltaPlusPlus;
	  case 2214: return G4INCL::DeltaPlus;
	  case 2114: return G4INCL::DeltaZero;
	  case 1114: return G4INCL::DeltaMinus;
	  case 221: return G4INCL::Eta;
	  case 223: return G4INCL::Omega;
	  case 331: return G4INCL::EtaPrime;
	  case 22: return G4INCL::Photon;
	  case 3122: return G4INCL::Lambda;
	  case 3222: return G4INCL::SigmaPlus;
	  case 3212: return G4INCL::SigmaZero;
	  case 3112: return G4INCL::SigmaMinus;
	  case -2212: return G4INCL::antiProton;
	  case 3312: return G4INCL::XiMinus;
	  case 3322: return G4INCL::XiZero;
	  case -2112: return G4INCL::antiNeutron;
	  case -3122: return G4INCL::antiLambda;
	  case -3222: return G4INCL::antiSigmaPlus;
	  case -3212: return G4INCL::antiSigmaZero;
	  case -3112: return G4INCL::antiSigmaMinus;
	  case -3312: return G4INCL::antiXiMinus;
	  case -3322: return G4INCL::antiXiZero;
	  case 321: return G4INCL::KPlus;
	  case 311: return G4INCL::KZero;
	  case -311: return G4INCL::KZeroBar;
	  case -321: return G4INCL::KMinus;
	  case 310: return G4INCL::KShort;
	  case 130: return G4INCL::KLong;
	  default:
		    return G4INCL::UnknownParticle;
	}
      }
  };
}

#endif /*  G4INCLPARTICLERECORD_H_  */
