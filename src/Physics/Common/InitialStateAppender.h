//____________________________________________________________________________
/*!

\class    genie::InitialStateAppender

\brief    Appends the initial state information to the event record.
          Is a concerete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  October 04, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _INITIAL_STATE_APPENDER_H_
#define _INITIAL_STATE_APPENDER_H_

#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {

class InitialStateAppender : public EventRecordVisitorI {

public :

  InitialStateAppender();
  InitialStateAppender(string config);
  ~InitialStateAppender();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;

private :

  void AddNeutrino       (GHepRecord * event_rec) const;
  void AddNucleus        (GHepRecord * event_rec) const;
  void AddStruckParticle (GHepRecord * event_rec) const;

};

}      // genie namespace

#endif // _INITIAL_STATE_APPENDER_H_
