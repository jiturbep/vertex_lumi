// Dear emacs, this is -*- c++ -*-
// $Id: TrigConfigSvcD3PD.h 465546 2011-10-31 14:25:20Z krasznaa $
#ifndef TRIGROOTANALYSIS_TRIGCONFIGSVCD3PD_H
#define TRIGROOTANALYSIS_TRIGCONFIGSVCD3PD_H

// STL include(s):
#include <map>
#include <string>
#include <vector>

// ROOT include(s):
#include <TNamed.h>

// Local include(s):
#include "VarHandle.h"

// Forward declaration(s):
class TTree;
class TBranch;

namespace D3PD {

   /**
    *  @short Class acting as a trigger configuration service on D3PDs
    *
    *         This class can be used to ease reading the information out of the
    *         trigger configuration metadata trees of D3PDs.
    *
    *         The user is required to get hold of the trigger configuration
    *         tree himself/herself, but after that this class takes care of loading
    *         the correct entries from the TTree.
    *
    * @author Attila Krasznahrokay <Attila.Krasznahorkay@cern.ch>
    *
    * $Revision: 465546 $
    * $Date: 2011-10-31 14:25:20 +0000 (Mon, 31 Oct 2011) $
    */
   class TrigConfigSvcD3PD : public ::TNamed {

      /// A convenience definition for an internal type
      typedef std::pair< ::Int_t, std::pair< ::Int_t, ::Int_t > > DBKeys_t;

   public:
      /// Constructor with a pointer to the D3PD trigger configuration tree
      TrigConfigSvcD3PD( ::TTree* confTree = 0 );
      /// Destructor
      ~TrigConfigSvcD3PD();

      /// Set the object to use a new configuration tree
      ::Bool_t SetConfigTree( ::TTree* confTree );

      /// Load the configuration belonging to a DB entry
      ::Bool_t Load( ::Int_t SMK, ::Int_t L1PSK, ::Int_t HLTPSK,
                     ::Bool_t silent = kFALSE );
      /// Check if a given LVL1 item / HLT chain is configured
      ::Bool_t IsConfigured( const std::string& name ) const;

      /// Get the CTP ID of a LVL1 item
      ::Int_t GetCTPId( const std::string& name ) const;
      /// Get the chain ID of an HLT chain
      ::Int_t GetChainId( const std::string& name ) const;
      /// Get the prescale of a trigger (LVL1 or HLT)
      ::Float_t GetPrescale( const std::string& name ) const;
      /// Get the rerun prescale value of an HLT chain
      ::Float_t GetRerunPrescale( const std::string& name ) const;
      /// Get the passthrough value of an HLT chain
      ::Float_t GetPassthrough( const std::string& name ) const;
      /// Get the lower chain name of an HLT chain
      const std::string& GetLowerChainName( const std::string& name ) const;

      /// Get the item name corrsponding to a given CTP ID
      const std::string& GetNameFromCTPId( ::Int_t ctpid ) const;
      /// Get the LVL2 chain name corrsponding to a given chain ID
      const std::string& GetL2NameFromChainId( ::Int_t chainid ) const;
      /// Get the EF chain name corrsponding to a given chain ID
      const std::string& GetEFNameFromChainId( ::Int_t chainid ) const;

      /// Get all the configured LVL1 item names
      const std::vector< std::string >& GetLVL1Items() const;
      /// Get all the configured HLT chain names
      const std::vector< std::string >& GetHLTChains() const;

      /// Print a summary of all the available configurations
      void PrintSummary();
      /// Print details about the currently loaded configuration
      void PrintConfig();

   private:
      /// Initialize the object from the specified TTree
      ::Bool_t Initialize();
      /// Reset the object
      void Reset();
      /// Extract the configured LVL1 item and HLT chain names
      ::Bool_t UpdateTriggerNames();
      /// Create the key object describing a configuration
      static DBKeys_t MakeDBKeys( ::Int_t SMK, ::Int_t L1PSK, ::Int_t HLTPSK );

      ::TTree* m_tree; ///< Pointer to the configuration tree

      ::Bool_t m_initialized; ///< Internal status flag

      //
      // Variables for reading the information from the configuration
      // TTree:
      //
      ::Long64_t m_entry;
#ifndef __MAKECINT__
      D3PDReaderPriv::VarHandle< ::Int_t > m_smk;
      D3PDReaderPriv::VarHandle< ::Int_t > m_l1psk;
      D3PDReaderPriv::VarHandle< ::Int_t > m_hltpsk;

      D3PDReaderPriv::VarHandle< std::map< std::string, int >* >   m_lvl1NameMap;
      D3PDReaderPriv::VarHandle< std::map< std::string, float >* > m_lvl1PrescaleMap;

      D3PDReaderPriv::VarHandle< std::map< std::string, int >* >   m_hltNameMap;
      D3PDReaderPriv::VarHandle< std::map< std::string, float >* > m_hltPrescaleMap;
      D3PDReaderPriv::VarHandle< std::map< std::string, float >* > m_hltRerunPrescaleMap;
      D3PDReaderPriv::VarHandle< std::map< std::string, float >* > m_hltPassthroughMap;

      D3PDReaderPriv::VarHandle< std::map< std::string, std::string >* > m_lowerChainNameMap;
#endif // __MAKECINT__

      std::map< int, std::string > m_lvl1IDMap; ///< Reverse of m_lvl1NameMap
      std::map< int, std::string > m_lvl2IDMap; ///< Reverse of m_hltNameMap (for LVL2 chains)
      std::map< int, std::string > m_efIDMap;   ///< Reverse of m_hltNameMap (for EF chains)

      std::map< DBKeys_t, ::Int_t > m_configKeyMap; ///< Internal cache of the configuration keys

      std::vector< std::string > m_confLVL1Items; ///< Internal cache of the LVL1 item names
      std::vector< std::string > m_confHLTChains; ///< Internal cache of the HLT chain names

      std::string m_dummy; ///< Dummy string needed for technical reasons

      ClassDef( D3PD::TrigConfigSvcD3PD, 0 )

   }; // class TrigConfigSvcD3PD

} // namespace D3PD

#endif // TRIGROOTANALYSIS_TRIGCONFIGSVCD3PD_H
