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
// $Id: G4INCLFermiConfigurationList.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)
//
// 23.04.2011 V.Ivanchenko: make this class to be responsible for
//            selection of decay channel and decay
//

#ifndef G4INCLFermiConfigurationList_h
#define G4INCLFermiConfigurationList_h 1

#include "G4INCLFermiConfiguration.hh"
#include "G4INCLVFermiFragment.hh"
#include "G4INCLFermiPhaseSpaceDecay.hh"
#include "G4INCLParticleVector.hh"
#include <vector>

class G4INCLFermiFragmentsPool;
class G4Pow;

class G4INCLFermiConfigurationList 
{
public:

  G4INCLFermiConfigurationList();

  ~G4INCLFermiConfigurationList();

  G4INCL::ParticleVector* GetFragments(const G4INCL::Cluster& theNucleus);
  
private:

  double CoulombBarrier(const std::vector<const G4INCLVFermiFragment*>& v);

  double DecayProbability(int A, double TotalE, G4INCLFermiConfiguration*);

  const std::vector<const G4INCLVFermiFragment*>*
  SelectConfiguration(int Z, int A, double mass);

  G4INCLFermiConfigurationList(const G4INCLFermiConfigurationList &right);  
  const G4INCLFermiConfigurationList & operator=(const G4INCLFermiConfigurationList &right);
  bool operator==(const G4INCLFermiConfigurationList &right) const;
  bool operator!=(const G4INCLFermiConfigurationList &right) const;
  
  G4INCLFermiFragmentsPool* thePool;

  std::vector<double> NormalizedWeights;
  
  // Kappa = V/V_0 it is used in calculation of Coulomb energy
  static const double Kappa;
  
  // Nuclear radius r0 (is a model parameter)
  static const double r0;

  double Coef;
  double ConstCoeff;
  size_t   nmax;
  
  G4Pow* g4pow;

  G4INCLFermiPhaseSpaceDecay thePhaseSpace;

};


#endif


