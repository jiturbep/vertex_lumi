// $Id: TrigDecisionToolD3PD.cxx 469469 2011-11-21 09:37:13Z krasznaa $

// STL include(s):
#include <algorithm>

// Local include(s):
#include "TrigRootAnalysis/TrigDecisionToolD3PD.h"

ClassImp( D3PD::TrigDecisionToolD3PD )

namespace D3PD {

   /**
    * To read any encoded trigger information from a D3PD, the file has to contain
    * a metadata tree called "TrigConfTree". This is located in the file in a
    * subdirectory called "<d3pd tree name>Meta". For instance in the egammaD3PDs
    * the tree is under "egammaMeta/TrigConfTree". If you don't have such a tree
    * in your file, the tool will not be usable.
    *
    * The tool of course also needs to access the event-level TTree holding the
    * event-wise trigger information.
    *
    * @param eventTree Pointer to the event-wise D3PD TTree
    * @param confTree  Pointer to "TrigConfTree"
    * @param prefix    Prefix used by the trigger tools. "trig_" by default.
    */
   TrigDecisionToolD3PD::TrigDecisionToolD3PD( TTree* eventTree, TTree* confTree,
                                               const TString& prefix )
      : TNamed( "TrigDecisionToolD3PD", "Trigger decision tool" ),
        Trig::DataAccess( eventTree, prefix ),
        Trig::ConfigAccess( confTree ),
        Trig::ChainGroupHandling() {

   }

   /**
    * Perhaps the most important function that people need. It can tell the user
    * whether a given LVL1 item or HLT chain passed in the current event.
    *
    * By default the function returns the "physics decision", but it can be asked
    * to give more detailed results as well. The function tries to be smart, so if
    * the user asks for information that is not available in the given D3PD, it
    * should just print a warning. (And return kFALSE.)
    *
    * For a description of the types of trigger results, have a look at the
    * documentation of the TrigDefs::DecisionTypes enumeration.
    *
    * @param name The name of the trigger we're interested in
    * @param type The type of the trigger decision to look up
    * @returns <code>kTRUE</code> if the trigger was passed in the current event,
    *          <code>kFALSE</code> otherwise
    */
   Bool_t TrigDecisionToolD3PD::IsPassed( const std::string& name,
                                          TrigDefs::DecisionTypes type ) {

      // Check if this trigger is configured for this event:
      if( ! GetConfigSvc().IsConfigured( name ) ) {
         Warning( "IsPassed", "Trigger not configured: %s",
                  name.c_str() );
         return kFALSE;
      }

      Int_t id = 0; // Helper variable

      //
      // If we're asking about a LVL1 item:
      //
      if( name.find( "L1_" ) == 0 ) {

         // Find out if the LVL1 item exists in the current configuration:
         if( ( id = GetConfigSvc().GetCTPId( name ) ) < 0 ) {
            Warning( "IsPassed", "LVL1 item not known: %s",
                     name.c_str() );
            return kFALSE;
         }

         // Variables helping in looking up the correct bit in the result:
         const Int_t word = id / 32;
         const Int_t bit  = id % 32;

         //
         // Give the result based on a few conditions:
         //
         switch( type ) {

         case TrigDefs::Physics:
            if( GetDetailLevel() >= 0 ) {
               return ( ( ( * GetL1Result( Trig::DataAccess::TAV ) )[ word ] >> bit ) & 0x1 );
            } else {
               Warning( "IsPassed",
                        "The requested information is not available in the input" );
               return kFALSE;
            }
            break;
         case TrigDefs::L1_TBP:
            if( GetDetailLevel() >= 1 ) {
               return ( ( ( * GetL1Result( Trig::DataAccess::TBP ) )[ word ] >> bit ) & 0x1 );
            } else {
               Warning( "IsPassed",
                        "The requested information is not available in the input" );
               return kFALSE;
            }
            break;
         case TrigDefs::L1_TAP:
            if( GetDetailLevel() >= 1 ) {
               return ( ( ( * GetL1Result( Trig::DataAccess::TAP ) )[ word ] >> bit ) & 0x1 );
            } else {
               Warning( "IsPassed",
                        "The requested information is not available in the input" );
               return kFALSE;
            }
            break;
         default:
            Warning( "IsPassed", "Decision type not recognised. Returning physics decision." );
            if( GetDetailLevel() >= 0 ) {
               return ( ( ( * GetL1Result( Trig::DataAccess::TAV ) )[ word ] >> bit ) & 0x1 );
            } else {
               Warning( "IsPassed",
                        "The requested information is not available in the input" );
               return kFALSE;
            }
            break;
         }
      }

      // Else we have to be talking about an HLT chain:
      if( ( id = GetConfigSvc().GetChainId( name ) ) < 0 ) {
         Warning( "IsPassed", "HLT chain not known: %s",
                  name.c_str() );
         return kFALSE;
      }

      //
      // If this is a LVL2 chain:
      //
      if( name.find( "L2_" ) == 0 ) {

         //
         // Decide on which kind of result the user wants:
         //
         std::vector< short >::const_iterator begin, end;
         switch( type ) {

         case TrigDefs::Physics:
            if( GetDetailLevel() >= 0 ) {
               begin = GetL2Result( Trig::DataAccess::Physics )->begin();
               end   = GetL2Result( Trig::DataAccess::Physics )->end();
            } else {
               Warning( "IsPassed",
                        "The requested information is not available in the input" );
               return kFALSE;
            }
            break;
         case TrigDefs::requireDecision:
            if( GetDetailLevel() >= 1 ) {
               begin = GetL2Result( Trig::DataAccess::Raw )->begin();
               end   = GetL2Result( Trig::DataAccess::Raw )->end();
            } else {
               Warning( "IsPassed",
                        "The requested information is not available in the input" );
               return kFALSE;
            }
            break;
         case TrigDefs::allowResurrected:
            if( GetDetailLevel() >= 2 ) {
               begin = GetL2Result( Trig::DataAccess::Resurrected )->begin();
               end   = GetL2Result( Trig::DataAccess::Resurrected )->end();
            } else {
               Warning( "IsPassed",
                        "The requested information is not available in the input" );
               return kFALSE;
            }
            break;
         case TrigDefs::passedThrough:
            if( GetDetailLevel() >= 2 ) {
               begin = GetL2Result( Trig::DataAccess::PassedThrough )->begin();
               end   = GetL2Result( Trig::DataAccess::PassedThrough )->end();
            } else {
               Warning( "IsPassed",
                        "The requested information is not available in the input" );
               return kFALSE;
            }
            break;
         default:
            Warning( "IsPassed",
                     "Decision type not recognised. Returning physics decision." );
            if( GetDetailLevel() >= 0 ) {
               begin = GetL2Result( Trig::DataAccess::Physics )->begin();
               end   = GetL2Result( Trig::DataAccess::Physics )->end();
            } else {
               Warning( "IsPassed",
                        "The requested information is not available in the input" );
               return kFALSE;
            }
            break;
         }

         //
         // See if the trigger chain met this condition: 
         //
         if( std::find( begin, end, id ) != end ) {
            return kTRUE;
         } else {
            return kFALSE;
         }
      }

      //
      // If this is an EF chain:
      //
      if( name.find( "EF_" ) == 0 ) {

         //
         // Decide on which kind of result the user wants:
         //
         std::vector< short >::const_iterator begin, end;
         switch( type ) {

         case TrigDefs::Physics:
            if( GetDetailLevel() >= 0 ) {
               begin = GetEFResult( Trig::DataAccess::Physics )->begin();
               end   = GetEFResult( Trig::DataAccess::Physics )->end();
            } else {
               Warning( "IsPassed",
                        "The requested information is not available in the input" );
               return kFALSE;
            }
            break;
         case TrigDefs::requireDecision:
            if( GetDetailLevel() >= 1 ) {
               begin = GetEFResult( Trig::DataAccess::Raw )->begin();
               end   = GetEFResult( Trig::DataAccess::Raw )->end();
            } else {
               Warning( "IsPassed",
                        "The requested information is not available in the input" );
               return kFALSE;
            }
            break;
         case TrigDefs::allowResurrected:
            if( GetDetailLevel() >= 2 ) {
               begin = GetEFResult( Trig::DataAccess::Resurrected )->begin();
               end   = GetEFResult( Trig::DataAccess::Resurrected )->end();
            } else {
               Warning( "IsPassed",
                        "The requested information is not available in the input" );
               return kFALSE;
            }
            break;
         case TrigDefs::passedThrough:
            if( GetDetailLevel() >= 2 ) {
               begin = GetEFResult( Trig::DataAccess::PassedThrough )->begin();
               end   = GetEFResult( Trig::DataAccess::PassedThrough )->end();
            } else {
               Warning( "IsPassed",
                        "The requested information is not available in the input" );
               return kFALSE;
            }
            break;
         default:
            Warning( "IsPassed",
                     "Decision type not recognised. Returning physics decision." );
            if( GetDetailLevel() >= 0 ) {
               begin = GetEFResult( Trig::DataAccess::Physics )->begin();
               end   = GetEFResult( Trig::DataAccess::Physics )->end();
            } else {
               Warning( "IsPassed",
                        "The requested information is not available in the input" );
               return kFALSE;
            }
            break;
         }

         //
         // See if the trigger chain met this condition: 
         //
         if( std::find( begin, end, id ) != end ) {
            return kTRUE;
         } else {
            return kFALSE;
         }
      }

      //
      // If everything went well, we should never reach this point. If we do, then
      // some of the assumptions about the trigger configuration have to be wrong.
      //
      Error( "IsPassed", "We've arrived at a serious programming bug!" );
      return kFALSE;
   }

   /**
    * This is a very slow function. It uses the <code>GetPassedL1Triggers(...)</code>,
    * <code>GetPassedL2Triggers(...)</code> and <code>GetPassedEFTriggers(...)</code>
    * functions, and merges their results.
    *
    * So if you know what kind of triggers you're looking for, you're better off
    * using one of the specific functions.
    *
    * @param type The trigger decision type you're interested in
    * @returns The names of all the triggers that fired in the current event
    */
   std::vector< std::string >
   TrigDecisionToolD3PD::GetPassedTriggers( TrigDefs::DecisionTypes type ) {

      // Get all the fired trigger from the 3 trigger levels:
      std::vector< std::string > result   = GetPassedL1Triggers( type );
      const std::vector< std::string > l2 = GetPassedL2Triggers( type );
      const std::vector< std::string > ef = GetPassedEFTriggers( type );

      // Now merge them:
      result.insert( result.end(), l2.begin(), l2.end() );
      result.insert( result.end(), ef.begin(), ef.end() );

      return result;
   }

   /**
    * Get a list of all the LVL1 triggers that passed in the current event.
    * Keep in mind that the function is quite slow. It has to do multiple lookups
    * in an std::map, that I don't see avoiding. It also returns the name list
    * by value, which can take a lot of CPU cycles. (Creating the object on the
    * heap and returning a pointer to it would be terribly prone to memory leaks.)
    *
    * Use this function only in interactive analysis, or just in the final steps
    * of an analysis. (On a "small" number of events.)
    *
    * @param type The trigger decision type you're interested in
    * @returns The names of all the LVL1 items that fired in the current event
    */
   std::vector< std::string >
   TrigDecisionToolD3PD::GetPassedL1Triggers( TrigDefs::DecisionTypes type ) {

      // The result object:
      std::vector< std::string > result;

      // Select the LVL1 result to investigate:
      const std::vector< unsigned int >* l1result = 0;
      switch( type ) {

      case TrigDefs::Physics:
         if( GetDetailLevel() >= 0 ) {
            l1result = GetL1Result( Trig::DataAccess::TAV );
         } else {
            Warning( "GetPassedL1Triggers",
                     "The requested information is not available in the input" );
            return result;
         }
         break;
      case TrigDefs::L1_TBP:
         if( GetDetailLevel() >= 1 ) {
            l1result = GetL1Result( Trig::DataAccess::TBP );
         } else {
            Warning( "GetPassedL1Triggers",
                     "The requested information is not available in the input" );
            return result;
         }
         break;
      case TrigDefs::L1_TAP:
         if( GetDetailLevel() >= 1 ) {
            l1result = GetL1Result( Trig::DataAccess::TAP );
         } else {
            Warning( "GetPassedL1Triggers",
                     "The requested information is not available in the input" );
            return result;
         }
         break;
      default:
         Warning( "GetPassedL1Triggers",
                  "Decision type not recognised. Returning physics decision." );
         if( GetDetailLevel() >= 0 ) {
            l1result = GetL1Result( Trig::DataAccess::TAV );
         } else {
            Warning( "GetPassedL1Triggers",
                     "The requested information is not available in the input" );
            return result;
         }
         break;
      }

      // Loop over all bits of the LVL1 result, and decode which LVL1 items have passed:
      for( Int_t id = 0; id < 256; ++id ) {

         // Variables helping in looking up the correct bit in the result:
         const Int_t word = id / 32;
         const Int_t bit  = id % 32;

         // If this bit is on, add the corresponding trigger to the result list:
         if( ( ( *l1result )[ word ] >> bit ) & 0x1 ) {
            result.push_back( GetConfigSvc().GetNameFromCTPId( id ) );
         }
      }

      return result;
   }

   /**
    * Get a list of all the LVL2 triggers that passed in the current event.
    * Keep in mind that the function is quite slow. It has to do multiple lookups
    * in an std::map, that I don't see avoiding. It also returns the name list
    * by value, which can take a lot of CPU cycles. (Creating the object on the
    * heap and returning a pointer to it would be terribly prone to memory leaks.)
    *
    * Use this function only in interactive analysis, or just in the final steps
    * of an analysis. (On a "small" number of events.)
    *
    * @param type The trigger decision type you're interested in
    * @returns The names of all the LVL2 items that fired in the current event
    */
   std::vector< std::string >
   TrigDecisionToolD3PD::GetPassedL2Triggers( D3PD::TrigDefs::DecisionTypes type ) {

      // The result object:
      std::vector< std::string > result;

      // Select the LVL1 result to investigate:
      const std::vector< short >* l2result = 0;
      switch( type ) {

      case TrigDefs::Physics:
         if( GetDetailLevel() >= 0 ) {
            l2result = GetL2Result( Trig::DataAccess::Physics );
         } else {
            Warning( "GetPassedL2Triggers",
                     "The requested information is not available in the input" );
            return result;
         }
         break;
      case TrigDefs::requireDecision:
         if( GetDetailLevel() >= 1 ) {
            l2result = GetL2Result( Trig::DataAccess::Raw );
         } else {
            Warning( "GetPassedL2Triggers",
                     "The requested information is not available in the input" );
            return result;
         }
         break;
      case TrigDefs::allowResurrected:
         if( GetDetailLevel() >= 2 ) {
            l2result = GetL2Result( Trig::DataAccess::Resurrected );
         } else {
            Warning( "GetPassedL2Triggers",
                     "The requested information is not available in the input" );
            return result;
         }
         break;
      case TrigDefs::passedThrough:
         if( GetDetailLevel() >= 2 ) {
            l2result = GetL2Result( Trig::DataAccess::PassedThrough );
         } else {
            Warning( "GetPassedL2Triggers",
                     "The requested information is not available in the input" );
            return result;
         }
         break;
      default:
         Warning( "GetPassedL2Triggers",
                  "Decision type not recognised. Returning physics decision." );
         if( GetDetailLevel() >= 0 ) {
            l2result = GetL2Result( Trig::DataAccess::Physics );
         } else {
            Warning( "GetPassedL2Triggers",
                     "The requested information is not available in the input" );
            return result;
         }
         break;
      }

      // Loop over all chain IDs and add all names to the result list:
      std::vector< short >::const_iterator itr = l2result->begin();
      std::vector< short >::const_iterator end = l2result->end();
      for( ; itr != end; ++itr ) {
         result.push_back( GetConfigSvc().GetL2NameFromChainId( *itr ) );
      }

      return result;
   }

   /**
    * Get a list of all the EF triggers that passed in the current event.
    * Keep in mind that the function is quite slow. It has to do multiple lookups
    * in an std::map, that I don't see avoiding. It also returns the name list
    * by value, which can take a lot of CPU cycles. (Creating the object on the
    * heap and returning a pointer to it would be terribly prone to memory leaks.)
    *
    * Use this function only in interactive analysis, or just in the final steps
    * of an analysis. (On a "small" number of events.)
    *
    * @param type The trigger decision type you're interested in
    * @returns The names of all the EF items that fired in the current event
    */
   std::vector< std::string >
   TrigDecisionToolD3PD::GetPassedEFTriggers( D3PD::TrigDefs::DecisionTypes type ) {

      // The result object:
      std::vector< std::string > result;

      // Select the LVL1 result to investigate:
      const std::vector< short >* efresult = 0;
      switch( type ) {

      case TrigDefs::Physics:
         if( GetDetailLevel() >= 0 ) {
            efresult = GetEFResult( Trig::DataAccess::Physics );
         } else {
            Warning( "GetPassedEFTriggers",
                     "The requested information is not available in the input" );
            return result;
         }
         break;
      case TrigDefs::requireDecision:
         if( GetDetailLevel() >= 1 ) {
            efresult = GetEFResult( Trig::DataAccess::Raw );
         } else {
            Warning( "GetPassedEFTriggers",
                     "The requested information is not available in the input" );
            return result;
         }
         break;
      case TrigDefs::allowResurrected:
         if( GetDetailLevel() >= 2 ) {
            efresult = GetEFResult( Trig::DataAccess::Resurrected );
         } else {
            Warning( "GetPassedEFTriggers",
                     "The requested information is not available in the input" );
            return result;
         }
         break;
      case TrigDefs::passedThrough:
         if( GetDetailLevel() >= 2 ) {
            efresult = GetEFResult( Trig::DataAccess::PassedThrough );
         } else {
            Warning( "GetPassedEFTriggers",
                     "The requested information is not available in the input" );
            return result;
         }
         break;
      default:
         Warning( "GetPassedEFTriggers",
                  "Decision type not recognised. Returning physics decision." );
         if( GetDetailLevel() >= 0 ) {
            efresult = GetEFResult( Trig::DataAccess::Physics );
         } else {
            Warning( "GetPassedEFTriggers",
                     "The requested information is not available in the input" );
            return result;
         }
         break;
      }

      // Loop over all chain IDs and add all names to the result list:
      std::vector< short >::const_iterator itr = efresult->begin();
      std::vector< short >::const_iterator end = efresult->end();
      for( ; itr != end; ++itr ) {
         result.push_back( GetConfigSvc().GetEFNameFromChainId( *itr ) );
      }

      return result;
   }

} // namespace D3PD
