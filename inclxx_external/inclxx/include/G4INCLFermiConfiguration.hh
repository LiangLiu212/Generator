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
// $Id: G4INCLFermiConfiguration.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)
// 23.04.2011 V.Ivanchenko: make this class to be a simple container and no physics

#ifndef G4INCLFermiConfiguration_h
#define G4INCLFermiConfiguration_h 1

#include "G4INCLVFermiFragment.hh"
#include "G4INCLParticleVector.hh"
#include "G4INCLCluster.hh"

class G4INCLFermiConfiguration 
{
public:
  // Constructors
  G4INCLFermiConfiguration(const std::vector<const G4INCLVFermiFragment*>&);

  ~G4INCLFermiConfiguration();

  G4INCL::ParticleVector* GetFragments(const G4INCL::Cluster & theNucleus);

  inline int GetA() const;
  inline int GetZ() const;
  inline double GetMass() const;
  
  inline const std::vector<const G4INCLVFermiFragment*>& GetFragmentList();

private:

  inline G4INCLFermiConfiguration(const G4INCLFermiConfiguration &);
  inline const G4INCLFermiConfiguration & operator=(const G4INCLFermiConfiguration &);
  inline bool operator==(const G4INCLFermiConfiguration &) const;
  inline bool operator!=(const G4INCLFermiConfiguration &) const;
  
  int    totalZ;
  int    totalA;

  double totalMass;

  std::vector<const G4INCLVFermiFragment*> Configuration;

};

inline int G4INCLFermiConfiguration::GetA() const
{
  return totalA;
}

inline int G4INCLFermiConfiguration::GetZ() const
{
  return totalZ;
}

inline double G4INCLFermiConfiguration::GetMass() const
{
  return totalMass;
}

inline const std::vector<const G4INCLVFermiFragment*>& 
G4INCLFermiConfiguration::GetFragmentList()
{
  return Configuration;
}

#endif


