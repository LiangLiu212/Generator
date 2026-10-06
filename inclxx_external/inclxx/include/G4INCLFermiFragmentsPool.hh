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
// $Id: G4INCLVFermiBreakUp.cc,v 1.5 2006-06-29 20:13:13 gunter Exp $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modifications:
// 01.04.2011 General cleanup by V.Ivanchenko - more clean usage of static


#ifndef G4INCLFermiFragmentsPool_hh 
#define G4INCLFermiFragmentsPool_hh 1

#include "G4INCLVFermiFragment.hh"
#include "G4INCLFermiConfiguration.hh"
#include <vector>

class G4INCLFermiFragmentsPool
{
public:

  static G4INCLFermiFragmentsPool* Instance();

  ~G4INCLFermiFragmentsPool();

  const std::vector<G4INCLFermiConfiguration*>* 
  GetConfigurationList(int Z, int A, double mass);

  const G4INCLVFermiFragment* GetFragment(int Z, int A);

  inline int GetMaxZ() const;

  inline int GetMaxA() const;
  
private:

  G4INCLFermiFragmentsPool();

  void Initialise();

  bool IsExist(int Z, int A, std::vector<const G4INCLVFermiFragment*>&);

  inline bool IsAvailable(int Z, int A);

  static G4ThreadLocal G4INCLFermiFragmentsPool* theInstance;

  std::vector<const G4INCLVFermiFragment*> fragment_pool;

  int maxZ;
  int maxA;
  int verbose;
 
  // list of configuration sorted by A for 1, 2, 3, 4 final fragments
  std::vector<G4INCLFermiConfiguration*> list1[17]; 
  std::vector<G4INCLFermiConfiguration*> list2[17]; 
  std::vector<G4INCLFermiConfiguration*> list3[17];
  std::vector<G4INCLFermiConfiguration*> list4[17];
  // list of exotic configurations
  std::vector<G4INCLFermiConfiguration*> listextra;
};

inline bool G4INCLFermiFragmentsPool::IsAvailable(int Z, int A)
{
  bool res = true;
  if     (2 == Z && 5 == A) { res = false; }
  else if(3 == Z && 5 == A) { res = false; }
  else if(4 == Z && 8 == A) { res = false; }
  else if(5 == Z && 9 == A) { res = false; }
  return res;
}

inline int G4INCLFermiFragmentsPool::GetMaxZ() const
{
  return maxZ;
}

inline int G4INCLFermiFragmentsPool::GetMaxA() const
{
  return maxA;
}

#endif

