//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steve Dytman <dytman \at pitt.edu>
 Univ. of Pittsburgh

 Joe Johnston <jpj13 \at pitt.edu>
 Univ. of Pittsburgh
*/
//____________________________________________________________________________

#include "Physics/QuasiElastic/XSection/NievesQELException.h"
#include "Framework/Messenger/Messenger.h"

using std::endl;
using namespace genie::exceptions;

//___________________________________________________________________________
namespace genie {
 namespace exceptions {
  ostream & operator<< (ostream& stream, const NievesQELException & exc)
  {
   exc.Print(stream);
   return stream;
  }
 }
}
//___________________________________________________________________________
NievesQELException::NievesQELException()
{
  this->Init();
}
//___________________________________________________________________________
NievesQELException::NievesQELException(const NievesQELException & exc)
{
  this->Copy(exc);
}
//___________________________________________________________________________
NievesQELException::~NievesQELException()
{

}
//___________________________________________________________________________
void NievesQELException::Init(void)
{
  fReason = "";
}
//___________________________________________________________________________
void NievesQELException::Copy(const NievesQELException & exc)
{
  fReason = exc.fReason;
}
//___________________________________________________________________________
void NievesQELException::Print(ostream & stream) const
{
  stream << "**EXCEPTION Reason: " << this->ShowReason() << endl;
}
//___________________________________________________________________________
