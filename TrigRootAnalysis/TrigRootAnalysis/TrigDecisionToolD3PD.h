// Dear emacs, this is -*- c++ -*-
// $Id: TrigDecisionToolD3PD.h 465546 2011-10-31 14:25:20Z krasznaa $
#ifndef TRIGROOTANALYSIS_TRIGDECISIONTOOLD3PD_H
#define TRIGROOTANALYSIS_TRIGDECISIONTOOLD3PD_H

// ROOT include(s):
#include <TNamed.h>
#include <TString.h>

// Local include(s):
#include "TrigConfigSvcD3PD.h"
#include "DataAccess.h"
#include "ConfigAccess.h"
#include "ChainGroupHandling.h"
#include "Conditions.h"

// Forward declaration(s):
class TTree;

namespace D3PD {

   /**
    *  @short Tool for interpreting trigger information in the D3PDs
    *
    *         This tool can be used to interpret the trigger information from the
    *         D3PDs, which is stored in an encoded way. The way in which the data
    *         is stored is explained in
    *
    *         https://twiki.cern.ch/twiki/bin/view/Atlas/TriggerD3PDMaker
    *
    *         Currently the only implemented feature is to ask the tool if a trigger
    *         item passed various conditions.
    *
    * @author Attila Krasznahorkay <Attila.Krasznahorkay@cern.ch>
    *
    * $Revision: 465546 $
    * $Date: 2011-10-31 14:25:20 +0000 (Mon, 31 Oct 2011) $
    */
   class TrigDecisionToolD3PD : public virtual TNamed,
                                public D3PD::Trig::DataAccess,
                                public D3PD::Trig::ConfigAccess,
                                public D3PD::Trig::ChainGroupHandling {

   public:
      /// Constructor specifying the D3PD trees that the tool should act on
      TrigDecisionToolD3PD( ::TTree* eventTree = 0, ::TTree* confTree = 0,
                            const ::TString& prefix = "trig_" );

      /// Find out whether a given LVL1 item / HLT chain was passed in the event
      ::Bool_t IsPassed( const std::string& name,
                         D3PD::TrigDefs::DecisionTypes type = D3PD::TrigDefs::Physics );

      /// Get the list of triggers from this chain group that passed the current event
      std::vector< std::string >
      GetPassedTriggers( D3PD::TrigDefs::DecisionTypes type = D3PD::TrigDefs::Physics );
      /// Get the list of LVL1 triggers from the chain group that passed the current event
      std::vector< std::string >
      GetPassedL1Triggers( D3PD::TrigDefs::DecisionTypes type = D3PD::TrigDefs::Physics );
      /// Get the list of LVL2 triggers from the chain group that passed the current event
      std::vector< std::string >
      GetPassedL2Triggers( D3PD::TrigDefs::DecisionTypes type = D3PD::TrigDefs::Physics );
      /// Get the list of EF triggers from the chain group that passed the current event
      std::vector< std::string >
      GetPassedEFTriggers( D3PD::TrigDefs::DecisionTypes type = D3PD::TrigDefs::Physics );

      ClassDef( D3PD::TrigDecisionToolD3PD, 0 )

   }; // class TrigDecisionToolD3PD

} // namespace D3PD

#endif // TRIGROOTANALYSIS_TRIGDECISIONTOOLD3PD_H
