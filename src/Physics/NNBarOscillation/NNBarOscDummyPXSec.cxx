//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Jeremy Hewes, Georgia Karagiorgi
 University of Manchester
*/
//____________________________________________________________________________

#include "Physics/NNBarOscillation/NNBarOscDummyPXSec.h"

using namespace genie;

//____________________________________________________________________________
NNBarOscDummyPXSec::NNBarOscDummyPXSec() :
XSecAlgorithmI("genie::NNBarOscDummyPXSec")
{

}
//____________________________________________________________________________
NNBarOscDummyPXSec::NNBarOscDummyPXSec(string config) :
XSecAlgorithmI("genie::NNBarOscDummyPXSec", config)
{

}
//____________________________________________________________________________
NNBarOscDummyPXSec::~NNBarOscDummyPXSec()
{

}
//____________________________________________________________________________
double NNBarOscDummyPXSec::XSec(const Interaction * , KinePhaseSpace_t ) const
{
  return 0;
}
//____________________________________________________________________________
double NNBarOscDummyPXSec::Integral(const Interaction * ) const
{
  return 0;
}
//____________________________________________________________________________
bool NNBarOscDummyPXSec::ValidProcess(const Interaction * ) const
{
  return true;
}
//____________________________________________________________________________
