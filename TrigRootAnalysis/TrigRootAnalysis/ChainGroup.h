// Dear emacs, this is -*- c++ -*-
// $Id: ChainGroup.h 465546 2011-10-31 14:25:20Z krasznaa $
#ifndef TRIGROOTANALYSIS_CHAINGROUP_H
#define TRIGROOTANALYSIS_CHAINGROUP_H

// STL include(s):
#include <vector>
#include <string>

// ROOT include(s):
#include <TNamed.h>

// Local include(s):
#include "Conditions.h"

namespace D3PD {

   // Forward declaration(s):
   class TrigConfigSvcD3PD;
   namespace Trig {
      class IDataAccess;
      class ChainGroupHandling;
   }

   /**
    *  @short Class implementing the chain group functionality
    *
    *         This class acts in pretty much the same way as the ChainGroup class in
    *         the TrigDecisionTool package. Such objects can only be produced by
    *         TrigDecisionToolD3PD, it has no public constructor.
    *
    * @author Attila Krasznahorkay <Attila.Krasznahorkay@cern.ch>
    *
    * $Revision: 465546 $
    * $Date: 2011-10-31 14:25:20 +0000 (Mon, 31 Oct 2011) $
    */
   class ChainGroup : public ::TNamed {

      /// Allow the TDT to create these objects
      friend class D3PD::Trig::ChainGroupHandling;

   protected:
      /// Constructor receiving all the needed information
      ChainGroup( const std::vector< std::string >& triggerNames,
                  const D3PD::Trig::IDataAccess& parent,
                  D3PD::TrigConfigSvcD3PD& svc );

   public:
      /// Find out if the chain group passed the selection in the event
      ::Bool_t IsPassed( D3PD::TrigDefs::DecisionTypes type = D3PD::TrigDefs::Physics );
      /// Get the overall prescale of the chain group
      ::Float_t GetPrescale();
      /// Get the list of triggers matching the selection
      const std::vector< std::string >& GetListOfTriggers();

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
      

   private:
      /// Update the object using the trigger configuration
      ::Bool_t Update();
      /// Function splitting a comma separated list into a vector
      static std::vector< std::string > ToVector( const std::string& names );

      const D3PD::Trig::IDataAccess& m_parent; ///< Interface for accessing the trigger data
      D3PD::TrigConfigSvcD3PD& m_configSvc; ///< Reference to the configuration service

      const std::vector< std::string > m_triggerNames; ///< The names given by the user

      ::Int_t m_smk; ///< The last SMK that was used to update the object
      std::vector< std::string > m_existingTriggers; ///< Names of the triggers
      std::vector< Int_t >       m_existingIDs;      ///< IDs of the triggers

      ClassDef( D3PD::ChainGroup, 0 )

   }; // class ChainGroup

} // namespace D3PD

#endif // TRIGROOTANALYSIS_CHAINGROUP_H
