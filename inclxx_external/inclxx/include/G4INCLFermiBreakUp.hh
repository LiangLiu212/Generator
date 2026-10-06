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
// $Id: G4INCLFermiBreakUp.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4INCLFermiBreakUp_h
#define G4INCLFermiBreakUp_h 1

#include "G4INCLVFermiBreakUp.hh"
#include "G4INCLFermiConfigurationList.hh"

class G4INCLFermiBreakUp : public G4INCLVFermiBreakUp 
{
public:

  G4INCLFermiBreakUp();
  virtual ~G4INCLFermiBreakUp();

  G4INCL::ParticleVector * BreakItUp(const G4INCL::Cluster &theNucleus);
  
private:

  G4INCLFermiBreakUp(const G4INCLFermiBreakUp &right);  
  const G4INCLFermiBreakUp & operator=(const G4INCLFermiBreakUp &right);
  bool operator==(const G4INCLFermiBreakUp &right) const;
  bool operator!=(const G4INCLFermiBreakUp &right) const;

  G4INCLFermiConfigurationList theConfigurationList;
};


#endif


