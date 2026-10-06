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
// $Id: G4INCLVFermiBreakUp.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)
//
// Modifications:
// 01.04.2011 General cleanup by V.Ivanchenko

#ifndef G4INCLVFermiBreakUp_h
#define G4INCLVFermiBreakUp_h 1

#include "G4INCLCluster.hh"
#include "G4INCLParticleVector.hh"

class G4INCLVFermiBreakUp 
{
public:

  G4INCLVFermiBreakUp();
  virtual ~G4INCLVFermiBreakUp();

  virtual G4INCL::ParticleVector * BreakItUp(const G4INCL::Cluster &theNucleus) = 0;
  
private:

  G4INCLVFermiBreakUp(const G4INCLVFermiBreakUp &right);  
  const G4INCLVFermiBreakUp & operator=(const G4INCLVFermiBreakUp &right);
  bool operator==(const G4INCLVFermiBreakUp &right) const;
  bool operator!=(const G4INCLVFermiBreakUp &right) const;
};

#endif


