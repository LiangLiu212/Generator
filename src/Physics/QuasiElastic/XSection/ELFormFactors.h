//____________________________________________________________________________
/*!

\class    genie::ELFormFactors

\brief    A class holding the Elastic Form Factors Ge,Gm

          This class is using the \b Strategy Pattern. \n

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  Apr 20, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _EL_FORM_FACTORS_H_
#define _EL_FORM_FACTORS_H_

#include <iostream>

#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"

using std::ostream;

namespace genie {

class ELFormFactors;
class Interaction;

ostream & operator << (ostream & stream, const ELFormFactors & ff);

class ELFormFactors {

public:
  ELFormFactors();
  ELFormFactors(const ELFormFactors & form_factors);
  virtual ~ELFormFactors() { }

  //! Attach an algorithm
  void   SetModel  (const ELFormFactorsModelI * model);

  //! Calculate the form factors for the input interaction using the attached algorithm
  void   Calculate (const Interaction * interaction);

  //! Get the computed form factor Gep
  double Gep (void) const { return fGep; }

  //! Get the computed form factor Gmp
  double Gmp (void) const { return fGmp; }

  //! Get the computed form factor Gen
  double Gen (void) const { return fGen; }

  //! Get the computed form factor Gmn
  double Gmn (void) const { return fGmn; }

  //! Get the attached model
  const ELFormFactorsModelI * Model (void) const {return fModel;}

  void   Reset    (Option_t * opt="");
  void   Copy     (const ELFormFactors & ff);
  bool   Compare  (const ELFormFactors & ff) const;
  void   Print    (ostream & stream) const;

  bool             operator == (const ELFormFactors & ff) const;
  ELFormFactors &  operator =  (const ELFormFactors & ff);
  friend ostream & operator << (ostream & stream, const ELFormFactors & ff);

private:

  double fGep;
  double fGmp;
  double fGen;
  double fGmn;

  const ELFormFactorsModelI * fModel;
};

}        // genie namespace

#endif   // _QEL_FORM_FACTORS_H_
