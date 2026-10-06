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
// $Id: G4ExcitationHandler.hh,v 1.13 2010-11-17 16:20:31 vnivanch Exp $
//
// Hadronic Process: Phase space decay for the Fermi BreakUp model
// by V. Lara
//
// Modifications:
// 01.04.2011 General cleanup by V.Ivanchenko: 
//          - IsotropicVector is inlined
//          - Momentum computation return zero or positive value
//          - DumpProblem method is added providing more information
//          - Reduced usage of exotic std functions  
//

#ifndef G4INCLFermiPhaseSpaceDecay_hh
#define G4INCLFermiPhaseSpaceDecay_hh 1

#include <vector>
#include <cmath>

#include "G4INCLLorentzVector.hh"
#include "G4INCLThreeVector.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"

class G4INCLFermiPhaseSpaceDecay
{
public:

  G4INCLFermiPhaseSpaceDecay();
  ~G4INCLFermiPhaseSpaceDecay();
  
  inline std::vector<G4INCL::LorentzVector*> * 
  Decay(const double,  const std::vector<double>&) const;

private:

  inline double PtwoBody(double E, double P1, double P2) const;
  
  inline G4INCL::ThreeVector IsotropicVector(const double Magnitude = 1.0) const;

  inline double BetaKopylov(int) const; 

  std::vector<G4INCL::LorentzVector*> * 
  TwoBodyDecay(double, const std::vector<double>&) const;

  std::vector<G4INCL::LorentzVector*> * 
  NBodyDecay(double, const std::vector<double>&) const;

  std::vector<G4INCL::LorentzVector*> * 
  KopylovNBodyDecay(double, const std::vector<double>&) const;

  void DumpProblem(double E, double P1, double P2, double P) const;

  G4INCLFermiPhaseSpaceDecay(const G4INCLFermiPhaseSpaceDecay&);
  const G4INCLFermiPhaseSpaceDecay & operator=(const G4INCLFermiPhaseSpaceDecay &); 
  bool operator==(const G4INCLFermiPhaseSpaceDecay&);
  bool operator!=(const G4INCLFermiPhaseSpaceDecay&);

};

inline double 
G4INCLFermiPhaseSpaceDecay::PtwoBody(double E, double P1, double P2) const
{
  double res = 0.0;
  double P = (E+P1+P2)*(E+P1-P2)*(E-P1+P2)*(E-P1-P2)/(4.0*E*E);
  if (P>0.0) { res = std::sqrt(P); }
  else { DumpProblem(E,P1,P2,P); }
  return res;
}

inline std::vector<G4INCL::LorentzVector*> * G4INCLFermiPhaseSpaceDecay::
Decay(double parent_mass_parameter, const std::vector<double>& fragment_masses) const
{
  return KopylovNBodyDecay(parent_mass_parameter,fragment_masses);
}

inline double G4INCLFermiPhaseSpaceDecay::BetaKopylov(int K) const
{
  int N = 3*K - 5;
  double xN = double(N);
  double F;
  //double Fmax=std::pow((3.*K-5.)/(3.*K-4.),(3.*K-5.)/2.)*std::sqrt(1-((3.*K-5.)/(3.*K-4.))); 
  // VI variant
  double Fmax = std::sqrt(std::pow(xN/(xN + 1),N)/(xN + 1)); 
  double chi;
  do {
    chi = G4INCL::Random::shoot();
    F = std::sqrt(std::pow(chi,N)*(1-chi));      
  } while ( Fmax*G4INCL::Random::shoot() > F);  
  return chi;
}

inline G4INCL::ThreeVector 
G4INCLFermiPhaseSpaceDecay::IsotropicVector(double Magnitude) const
  // Samples a isotropic random vectorwith a magnitud given by Magnitude.
  // By default Magnitude = 1.0
{
  double CosTheta = 2.0*G4INCL::Random::shoot() - 1.0;
  double SinTheta = std::sqrt((1. - CosTheta)*(1. + CosTheta));
  double Phi = G4INCL::Math::twoPi*G4INCL::Random::shoot();
  G4INCL::ThreeVector Vector(Magnitude*std::cos(Phi)*SinTheta,
		       Magnitude*std::sin(Phi)*SinTheta,
		       Magnitude*CosTheta);
  return Vector;
}

#endif
