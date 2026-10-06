//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4INCLVFermiFragment.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)
//
// Modifications:
// 01.04.2011 General cleanup by V.Ivanchenko

#ifndef G4INCLVFermiFragment_h
#define G4INCLVFermiFragment_h 1

#include "G4INCLCluster.hh"
#include "G4INCLParticleVector.hh"
#include "G4INCLLorentzVector.hh"

class G4INCLVFermiFragment 
{
public:

  G4INCLVFermiFragment(int anA, int aZ, int Pol, double ExE);

  virtual ~G4INCLVFermiFragment();
  
private:

  G4INCLVFermiFragment(const G4INCLVFermiFragment &right);  
  const G4INCLVFermiFragment & operator=(const G4INCLVFermiFragment &right);
  bool operator==(const G4INCLVFermiFragment &right) const;
  bool operator!=(const G4INCLVFermiFragment &right) const;
  
public:

  virtual G4INCL::ParticleVector * GetFragment(const G4INCL::LorentzVector & aMomentum) const = 0;

  inline int GetA(void) const 
  {
    return A;
  }
  
  inline int GetZ(void) const 
  {
    return Z;
  }
  
  inline int GetPolarization(void) const 
  {
    return Polarization;
  }

  inline double GetExcitationEnergy(void) const 
  {
    return ExcitEnergy;
  }

  inline double GetFragmentMass(void) const
  {
    return fragmentMass;
  }

  inline double GetTotalEnergy(void) const
  {
    return (GetFragmentMass() + GetExcitationEnergy());
  }

  inline bool IsStable() const
  {
    return isStable;
  }

protected:

  bool isStable;

  int A;

  int Z;

  int Polarization;

  double ExcitEnergy;

  double fragmentMass;

};


#endif


