//____________________________________________________________________________
/*!

\class    genie::INCLNuclearModel

\brief    This class is a hook for  nuclear models and allows associating each
          one of them with specific nuclei.
          Is a concrete implementation of the NuclearModelI interface.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  May 07, 2004

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#ifndef _INCL_NUCLEAR_MODEL_H_
#define _INCL_NUCLEAR_MODEL_H_

#include "Physics/NuclearState/NuclearModelI.h"

namespace genie {

class INCLNuclearModel : public NuclearModelI {

public:
  INCLNuclearModel();
  INCLNuclearModel(string config);
  virtual ~INCLNuclearModel();

  using NuclearModelI::GenerateNucleon;  // inherit versions not overridden here
  using NuclearModelI::Prob;

  //-- Allow GenerateNucleon to be called with a radius
  virtual bool   GenerateNucleon (const Target & t,
                                  double hitNucleonRadius) const;
  virtual double  Prob           (double p, double w, const Target & t,
                                  double hitNucleonRadius) const;

  //-- implement the NuclearModelI interface
  bool GenerateNucleon (const Target & t) const {
    return GenerateNucleon(t,0.0);
  }
  double Prob (double p, double w, const Target & t) const {
    return Prob(p,w,t,0.0);
  }
  NuclearModel_t ModelType       (const Target & t) const;

  virtual double FermiMomentum( const Target & t, int nucleon_pdg ) const ;
  virtual double LocalFermiMomentum( const Target & t, int nucleon_pdg, double radius ) const ;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string config);

private:
  void LoadConfig(void);

};

}      // genie namespace
#endif // _INCL_NUCLEAR_MODEL_H_
