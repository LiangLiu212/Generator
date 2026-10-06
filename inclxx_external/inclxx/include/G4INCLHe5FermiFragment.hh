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
// $Id: G4INCLHe5FermiFragment.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)
//
// Modifications:
// 01.04.2011 General cleanup by V.Ivanchenko: integer Z and A, constructor

#ifndef G4INCLHe5FermiFragment_h
#define G4INCLHe5FermiFragment_h 1

#include "G4INCLUnstableFermiFragment.hh"

class G4INCLHe5FermiFragment : public G4INCLUnstableFermiFragment
{
public:

  G4INCLHe5FermiFragment(int anA, int aZ, int Pol, double ExE); 

  ~G4INCLHe5FermiFragment();
  
private:

  G4INCLHe5FermiFragment(const G4INCLHe5FermiFragment &right);
  const G4INCLHe5FermiFragment & operator=(const G4INCLHe5FermiFragment &right);
  bool operator==(const G4INCLHe5FermiFragment &right) const;
  bool operator!=(const G4INCLHe5FermiFragment &right) const;
  
};


#endif


