//____________________________________________________________________________
/*!

\class   genie::EventRecordVisitorI

\brief   Defines the EventRecordVisitorI interface.
         Concrete implementations of this interface use the 'Visitor' Design
         Pattern to perform an operation on an EventRecord.

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

\created October 04, 2004

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org       
*/
//____________________________________________________________________________

#ifndef _EVENT_RECORD_VISITOR_I_H_
#define _EVENT_RECORD_VISITOR_I_H_

#include "Framework/Algorithm/Algorithm.h"

namespace genie {

class GHepRecord;

class EventRecordVisitorI : public Algorithm {

public :

  virtual ~EventRecordVisitorI();

  //-- define the EventRecordVisitorI interface

  virtual void ProcessEventRecord(GHepRecord * event_rec) const = 0;

protected :

  EventRecordVisitorI();
  EventRecordVisitorI(string name);
  EventRecordVisitorI(string name, string config);
};

}      // genie namespace

#endif // _EVENT_RECORD_VISITOR_I_H_
