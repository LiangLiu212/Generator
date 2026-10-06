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
// $Id: G4INCLStableFermiFragment.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)
//
// Modifications:
// 01.04.2011 General cleanup by V.Ivanchenko

#ifndef G4INCLStableFermiFragment_h
#define G4INCLStableFermiFragment_h 1

#include "G4INCLVFermiFragment.hh"
#include "G4INCLParticleVector.hh"
#include "G4INCLFermiPhaseSpaceDecay.hh"
#include <vector>

class G4INCLStableFermiFragment : public G4INCLVFermiFragment 
{
public:

  G4INCLStableFermiFragment(int anA, int aZ, int Pol, double ExE);

  virtual ~G4INCLStableFermiFragment();

  virtual G4INCL::ParticleVector * GetFragment(const G4INCL::LorentzVector & aMomentum) const;
  
private:

  G4INCLStableFermiFragment(const G4INCLStableFermiFragment &right);
  
  const G4INCLStableFermiFragment & operator=(const G4INCLStableFermiFragment &right);
  bool operator==(const G4INCLStableFermiFragment &right) const;
  bool operator!=(const G4INCLStableFermiFragment &right) const;

  G4INCLFermiPhaseSpaceDecay thePhaseSpace;
  std::vector<double> decayMasses;
};


#endif


