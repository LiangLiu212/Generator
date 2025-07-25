//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include "Framework/Registry/RegistryItemTypeDef.h"

using std::endl;

//____________________________________________________________________________
RgAlg::RgAlg()
{

}
//____________________________________________________________________________
RgAlg::RgAlg(string n, string c) :
name(n),
config(c)
{

}
//____________________________________________________________________________
RgAlg::~RgAlg()
{

}
//____________________________________________________________________________
ostream & operator << (ostream & stream, const RgAlg & alg)
{
  stream << alg.name << "/" << alg.config;
  return stream;
}
//____________________________________________________________________________
RgAlg & RgAlg::operator = (const RgAlg & alg)
{
  this->name   = alg.name;
  this->config = alg.config;
  return (*this);
}
//____________________________________________________________________________
