// Dear emacs, this is -*- c++ -*-
// $Id: PyTrigDecisionToolD3PD.h 465546 2011-10-31 14:25:20Z krasznaa $
#ifndef TRIGROOTANALYSIS_PYTRIGDECISIONTOOLD3PD_H
#define TRIGROOTANALYSIS_PYTRIGDECISIONTOOLD3PD_H

// ROOT include(s):
#include <TNamed.h>
#include <TString.h>

// Local include(s):
#include "Conditions.h"
#include "TrigDecisionToolD3PD.h"

// Forward declaration(s):
class TTree;

namespace D3PD {

   // Forward declaration(s):
   class TrigConfigSvcD3PD;

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
    *         This is a simplified version of the tool for PyROOT, as Python
    *         has trouble understanding the inheritance structure of
    *         <code>D3PD::TrigDecisionToolD3PD</code>.
    *
    * @author Doug Schouten <doug.schouten@triumf.ca>
    * @author Attila Krasznahorkay <Attila.Krasznahorkay@cern.ch>
    *
    * $Revision: 465546 $
    * $Date: 2011-10-31 14:25:20 +0000 (Mon, 31 Oct 2011) $
    */
   class PyTrigDecisionToolD3PD : public ::TNamed {

   public:
      /// Constructor specifying the D3PD trees that the tool should act on
      PyTrigDecisionToolD3PD( ::TTree* eventTree = 0, ::TTree* confTree = 0,
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

      /// Get the pointer to the event-wise tree currently being used
      ::TTree* GetEventTree() const;
      /// Set the pointer to the event-wise tree
      ::Bool_t SetEventTree( ::TTree* eventTree );
      /// Write the branches to an output tree
      ::Bool_t WriteToTree( ::TTree* outputTree, ::Int_t detailLevel = 10 );

      /// Load a new entry (new event)
      ::Int_t GetEntry( ::Long64_t entry, ::Int_t getall = 0 );

      /// Get the detail level that the D3PD was produced with
      ::Int_t GetDetailLevel() const;

      /// Create a chain group from a single pattern
      D3PD::ChainGroup GetChainGroup( const std::string& pattern );

      /// Return a reference to the trigger configuration service
      D3PD::TrigConfigSvcD3PD& GetConfigSvc();

   private:
      /// The "c++ version" of the tool
      TrigDecisionToolD3PD m_tdt;

      ClassDef( D3PD::PyTrigDecisionToolD3PD, 0 )

   }; // class PyTrigDecisionToolD3PD

} // namespace D3PD

#endif // TRIGROOTANALYSIS_PYTRIGDECISIONTOOLD3PD_H
