// $Id: TrigConfigSvcD3PD.cxx 465546 2011-10-31 14:25:20Z krasznaa $

// ROOT include(s):
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>

// Local include(s):
#include "TrigRootAnalysis/TrigConfigSvcD3PD.h"

ClassImp( D3PD::TrigConfigSvcD3PD )

namespace D3PD {

   /**
    * The class allows itself to be instantiated without a proper configuration tree.
    * This is supposed to allow one to create the configuration object only once, and
    * then use it to read the configuration from all the files which are opened in
    * the analysis job.
    *
    * @param confTree Pointer to a trigger configuration metadata TTree
    */
   TrigConfigSvcD3PD::TrigConfigSvcD3PD( ::TTree* confTree )
      : TNamed( "TrigConfigSvcD3PD", "Trigger configuration service" ),
        m_tree( confTree ), m_initialized( kFALSE ), m_entry( 0 ),
        m_smk( this, "SMK", &m_entry ),
        m_l1psk( this, "L1PSK", &m_entry ),
        m_hltpsk( this, "HLTPSK", &m_entry ),
        m_lvl1NameMap( this, "LVL1NameMap", &m_entry ),
        m_lvl1PrescaleMap( this, "LVL1PrescaleMap", &m_entry ),
        m_hltNameMap( this, "HLTNameMap", &m_entry ),
        m_hltPrescaleMap( this, "HLTPrescaleMap", &m_entry ),
        m_hltRerunPrescaleMap( this, "HLTRerunPrescaleMap", &m_entry ),
        m_hltPassthroughMap( this, "HLTPassthroughMap", &m_entry ),
        m_lowerChainNameMap( this, "HLTLowerChainNameMap", &m_entry ),
        m_dummy( "" ) {

      if( m_tree ) Initialize();
   }

   /**
    * The destructor tried to clean up after the object as good as it can.
    */
   TrigConfigSvcD3PD::~TrigConfigSvcD3PD() {

      Reset();
   }

   /**
    * This function can be used to set the object to use a new trigger configuration
    * metadata TTree. If you give it a null-pointer, the object resets itself, and lets
    * go of the previous TTree without attaching itself to a new TTree yet.
    *
    * Note that unfortunately we can't just check if the passed TTree pointer
    * is the same as the cached one. ROOT has the habit of putting the same sort
    * of TTree into the same memory space when the second file is opened for instance.
    * I've seen some weird problems because of this.
    *
    * @param confTree Pointer to a trigger configuration metadata TTree
    * @returns <code>kTRUE</code> if the operation was successful,
    *          <code>kFALSE</code> otherwise
    */
   ::Bool_t TrigConfigSvcD3PD::SetConfigTree( ::TTree* confTree ) {

      //
      // If the user gave a null-pointer, reset the object and return gracefully:
      //
      if( ! confTree ) {
         Reset();
         return kTRUE;
      }

      //
      // Reset the object and start using a new TTree:
      //
      Reset();
      m_tree = confTree;
      return Initialize();
   }

   /**
    * This function can be used to load a new particular configuration in memory.
    * The object can of course only load configurations which exist in the input
    * file, and complains if the user wants to load something non-available.
    *
    * The function is fairly quick, so it can in principle be used in every event
    * in an event loop.
    *
    * @param SMK Super Master Key
    * @param L1PSK LVL1 Prescale Key
    * @param HLTPSK HLT Prescale Key
    * @param silent Set to <code>kTRUE</code> if you don't want info messages
    *               from the function
    * @returns <code>kTRUE</code> if the operation was successful,
    *          <code>kFALSE</code> otherwise
    */
   ::Bool_t TrigConfigSvcD3PD::Load( ::Int_t SMK, ::Int_t L1PSK, ::Int_t HLTPSK,
                                     ::Bool_t silent ) {

      //
      // Make sure that the object is initialized:
      //
      if( ! Initialize() ) {
         Error( "Load", "Couldn't load configuration" );
         return kFALSE;
      }

      //
      // Make quick work if the correct entry is already loaded:
      //
      if( ( SMK == m_smk() ) && ( L1PSK == m_l1psk() ) && ( HLTPSK == m_hltpsk() ) ) {
         return kTRUE;
      }

      //
      // Find which entry corresponds to this configuration:
      //
      DBKeys_t keys = MakeDBKeys( SMK, L1PSK, HLTPSK );
      std::map< DBKeys_t, Int_t >::const_iterator itr;
      if( ( itr = m_configKeyMap.find( keys ) ) == m_configKeyMap.end() ) {
         Error( "Load", "Configuration with the following keys not available:" );
         Error( "Load", "  SMK: %i, L1PSK: %i, HLTPSK: %i",
                SMK, L1PSK, HLTPSK );
         Error( "Load", "Available configurations are:" );
         std::map< DBKeys_t, Int_t >::const_iterator db_itr = m_configKeyMap.begin();
         std::map< DBKeys_t, Int_t >::const_iterator db_end = m_configKeyMap.end();
         for( ; db_itr != db_end; ++db_itr ) {
            Error( "Load", "  SMK: %i, L1PSK: %i, HLTPSK: %i",
                   db_itr->first.first, db_itr->first.second.first,
                   db_itr->first.second.second );
         }
         return kFALSE;
      }

      //
      // Load the requested entry:
      //
      m_entry = itr->second;
      TChain* chain = dynamic_cast< TChain* >( m_tree );
      if( chain ) {
         m_entry = chain->LoadTree( itr->second );
      }

      //
      // Update the cached trigger names:
      //
      if( ! UpdateTriggerNames() ) {
         Error( "Load", "Failed to cache the configured trigger names" );
         return kFALSE;
      }

      if( ! silent ) {
         Info( "Load", "Loaded configuration: SMK: %i, L1PSK: %i, HLTPSK: %i",
               SMK, L1PSK, HLTPSK );
      }
      return kTRUE;
   }

   /**
    * This function checks whether a particular LVL1 item or HLT chain is available
    * in the currently loaded configuration.
    *
    * @param name Name of the LVL1 item or HLT chain
    * @returns <code>kTRUE</code> if the trigger is available in the configuration,
    *          <code>kFALSE</code> otherwise
    */
   ::Bool_t TrigConfigSvcD3PD::IsConfigured( const std::string& name ) const {

      // Do a small sanity check:
      if( ( ! m_hltNameMap.IsAvailable() ) || ( ! m_lvl1NameMap.IsAvailable() ) ) {
         Error( "IsConfigured", "Configuration meta-tree not initialized yet!" );
         return kFALSE;
      }

      //
      // Try to find it as an HLT chain (I guess this is a more common query):
      //
      if( m_hltNameMap()->find( name ) != m_hltNameMap()->end() ) {
         return kTRUE;
      }

      //
      // Now try to find it as a LVL1 item:
      //
      if( m_lvl1NameMap()->find( name ) != m_lvl1NameMap()->end() ) {
         return kTRUE;
      }

      // If it hasn't been found, then it wasn't configured:
      return kFALSE;
   }

   /**
    * This function can return the CTP ID of a LVL1 item. The CTP ID is an ID from 0
    * to 255 that specifies which bit in the CTP result shows the result of the
    * particular LVL1 item.
    *
    * @param name Name of the LVL1 item
    * @returns The CTP ID belonging to the specified LVL1 item
    */
   ::Int_t TrigConfigSvcD3PD::GetCTPId( const std::string& name ) const {

      // Do a small sanity check:
      if( ! m_lvl1NameMap.IsAvailable() ) {
         Error( "GetCTPId", "Configuration meta-tree not initialized yet!" );
         return kFALSE;
      }

      //
      // Find the ID of the requested item:
      //
      std::map< std::string, int >::const_iterator itr;
      if( ( itr = m_lvl1NameMap()->find( name ) ) == m_lvl1NameMap()->end() ) {
         Warning( "GetCTPId",  "LVL1 item not known: %s",
                  name.c_str() );
         return -1;
      }

      return itr->second;
   }

   /**
    * This function can return the chain ID of an HLT chain. The chain ID is an ID from 0
    * to 8192, specifying which bit in the HLT result shows the result of the particular
    * HLT chain.
    *
    * @param name Name of the HLT chain
    * @returns The chain ID belonging to the specified HLT chain
    */
   ::Int_t TrigConfigSvcD3PD::GetChainId( const std::string& name ) const {

      // Do a small sanity check:
      if( ! m_hltNameMap.IsAvailable() ) {
         Error( "GetChainId", "Configuration meta-tree not initialized yet!" );
         return kFALSE;
      }

      //
      // Find the ID of the requested chain:
      //
      std::map< std::string, int >::const_iterator itr;
      if( ( itr = m_hltNameMap()->find( name ) ) == m_hltNameMap()->end() ) {
         Warning( "GetChainId", "HLT chain not known: %s",
                  name.c_str() );
         return -1;
      }

      return itr->second;
   }

   /**
    * This function can return the prescale value of a LVL1 item or HLT chain.
    * Note that the function doesn't return the cumulative prescale! So for an EF trigger
    * chain it only returns the prescale value that was used at the EF level, not
    * taking the LVL1 and LVL2 prescales into account.
    *
    * @param name Name of the LVL1 item or HLT chain name
    * @returns The prescale value of the specified trigger
    */
   ::Float_t TrigConfigSvcD3PD::GetPrescale( const std::string& name ) const {

      // Do a small sanity check:
      if( ( ! m_hltPrescaleMap.IsAvailable() ) || ( ! m_lvl1PrescaleMap.IsAvailable() ) ) {
         Error( "GetPrescale", "Configuration meta-tree not initialized yet!" );
         return kFALSE;
      }

      //
      // Check if this is an HLT chain:
      //
      std::map< std::string, float >::const_iterator itr;
      if( ( itr = m_hltPrescaleMap()->find( name ) ) != m_hltPrescaleMap()->end() ) {
         return itr->second;
      }

      //
      // If not, then it should be a LVL1 item:
      //
      if( ( itr = m_lvl1PrescaleMap()->find( name ) ) != m_lvl1PrescaleMap()->end() ) {
         return itr->second;
      }

      Warning( "GetPrescale", "Trigger not known: %s",
               name.c_str() );
      return -1.0;
   }

   /**
    * This function can return the rerun prescale of a given HLT chain. This is
    * usually only needed for in-depth trigger analyses.
    *
    * @param name The name of the HLT chain
    * @returns The rerun prescale value belonging to the specified HLT chain
    */
   ::Float_t TrigConfigSvcD3PD::GetRerunPrescale( const std::string& name ) const {

      //
      // Check if this information is available:
      //
      if( ! m_hltRerunPrescaleMap.IsAvailable() ) {
         Warning( "GetRerunPrescale", "The requested information is not available" );
         return -1.0;
      }

      //
      // Find the rerun prescale of the requested chain:
      //
      std::map< std::string, float >::const_iterator itr;
      if( ( itr = m_hltRerunPrescaleMap()->find( name ) ) != m_hltRerunPrescaleMap()->end() ) {
         return itr->second;
      }

      Warning( "GetRerunPrescale", "HLT chain now known: %s",
               name.c_str() );
      return -1.0;
   }

   /**
    * This function can return the passthrough value set for a given HLT chain.
    * This is usually only needed for in-depth trigger analyses.
    *
    * @param name The name of the HLT chain
    * @returns The passthrough value belonging to the specified HLT chain
    */
   ::Float_t TrigConfigSvcD3PD::GetPassthrough( const std::string& name ) const {

      //
      // Check if this information is available:
      //
      if( ! m_hltPassthroughMap.IsAvailable() ) {
         Warning( "GetPassthrough", "The requested information is not available" );
         return -1.0;
      }

      //
      // Find the passthrough value of the requested chain:
      //
      std::map< std::string, float >::const_iterator itr;
      if( ( itr = m_hltPassthroughMap()->find( name ) ) != m_hltPassthroughMap()->end() ) {
         return itr->second;
      }

      Warning( "GetPassthrough", "HLT chain now known: %s",
               name.c_str() );
      return -1.0;
   }

   /**
    * This function can be used to find the name of the HLT chain / LVL1 item that seeds
    * another HLT chain. For instance "EF_mu20" is usually seeded by "L2_mu20".
    *
    * Note: This information may not be available in some older D3PDs.
    *
    * @param name The name of the HLT chain
    * @returns The name of the item/chain seeding the specified HLT chain
    */
   const std::string& TrigConfigSvcD3PD::GetLowerChainName( const std::string& name ) const {

      //
      // Check if this information is available:
      //
      if( ! m_lowerChainNameMap.IsAvailable() ) {
         Warning( "GetLowerChainName", "The requested information is not available" );
         return m_dummy;
      }

      //
      // Find the lower chain name of the requested chain:
      //
      std::map< std::string, std::string >::const_iterator itr;
      if( ( itr = m_lowerChainNameMap()->find( name ) ) != m_lowerChainNameMap()->end() ) {
         return itr->second;
      }

      Warning( "GetLowerChainName", "HLT chain now known: %s",
               name.c_str() );
      return m_dummy;
   }

   /**
    * This function can be used to do a "reverse lookup" in the LVL1 configuration.
    * Meaning, that if you know that a certain CTP Id fired in the current event,
    * you can use this function to check which LVL1 item that corresponds to.
    *
    * @param ctpid The CTP Id of the item in question
    * @returns The name of the LVL1 item, or an empty string if the lookup is
    *          not successful
    */
   const std::string& TrigConfigSvcD3PD::GetNameFromCTPId( ::Int_t ctpid ) const {

      //
      // Find this CTP Id in the internal reverse list:
      //
      std::map< int, std::string >::const_iterator itr = m_lvl1IDMap.find( ctpid );
      if( itr != m_lvl1IDMap.end() ) {
         return itr->second;
      }

      Warning( "GetNameFromCTPId", "Couldn't find LVL1 item for CTP Id: %i",
               ctpid );
      return m_dummy;
   }

   /**
    * This function can be used to do a "reverse lookup" in the LVL2 configuration.
    * Meaning, that if you know that a certain Chain Id fired in the current event,
    * you can use this function to check which LVL2 chain that corresponds to.
    *
    * @param chainid The Chain Id of the chain in question
    * @returns The name of the LVL2 chain, or an empty string if the lookup is
    *          not successful
    */
   const std::string& TrigConfigSvcD3PD::GetL2NameFromChainId( ::Int_t chainid ) const {

      //
      // Find this Chain Id in the internal reverse list:
      //
      std::map< int, std::string >::const_iterator itr = m_lvl2IDMap.find( chainid );
      if( itr != m_lvl2IDMap.end() ) {
         return itr->second;
      }

      Warning( "GetL2NameFromChainId", "Couldn't find LVL2 chain for Chain Id: %i",
               chainid );
      return m_dummy;
   }

   /**
    * This function can be used to do a "reverse lookup" in the EF configuration.
    * Meaning, that if you know that a certain Chain Id fired in the current event,
    * you can use this function to check which EF chain that corresponds to.
    *
    * @param chainid The Chain Id of the chain in question
    * @returns The name of the EF chain, or an empty string if the lookup is
    *          not successful
    */
   const std::string& TrigConfigSvcD3PD::GetEFNameFromChainId( ::Int_t chainid ) const {

      //
      // Find this Chain Id in the internal reverse list:
      //
      std::map< int, std::string >::const_iterator itr = m_efIDMap.find( chainid );
      if( itr != m_efIDMap.end() ) {
         return itr->second;
      }

      Warning( "GetEFNameFromChainId", "Couldn't find EF chain for Chain Id: %i",
               chainid );
      return m_dummy;
   }

   /**
    * This function can be used to get a list of all the configured LVL1 trigger items
    * for the current configuration.
    *
    * @returns A list of all the configured LVL1 trigger item names
    */
   const std::vector< std::string >& TrigConfigSvcD3PD::GetLVL1Items() const {

      return m_confLVL1Items;
   }

   /**
    * This function can be used to get a list of all the configured HLT trigger chains
    * for the current configuration.
    *
    * @returns A list of all the configured HLT trigger chain names
    */
   const std::vector< std::string >& TrigConfigSvcD3PD::GetHLTChains() const {

      return m_confHLTChains;
   }

   /**
    * This function can be used to print some short summary about the trigger
    * configurations stored in a given D3PD file. Mostly just for debugging.
    */
   void TrigConfigSvcD3PD::PrintSummary() {

      //
      // Make sure that the object is initialized:
      //
      if( ! Initialize() ) {
         Error( "PrintSummary", "Couldn't load configuration" );
         return;
      }

      //
      // Print some minimal information about the available configurations:
      //
      Info( "PrintSummary", "The following configurations are available:" );
      std::map< DBKeys_t, Int_t >::const_iterator db_itr = m_configKeyMap.begin();
      std::map< DBKeys_t, Int_t >::const_iterator db_end = m_configKeyMap.end();
      for( ; db_itr != db_end; ++db_itr ) {
         if( ! Load( db_itr->first.first, db_itr->first.second.first,
                     db_itr->first.second.second, kTRUE ) ) {
            Fatal( "PrintSummary", "Couldn't load an apparently existing configuration!" );
            Fatal( "PrintSummary", "There is some serious problem with the code!" );
            return;
         }

         Info( "PrintSummary", "  SMK: %i, L1PSK: %i, HLTPSK: %i, "
               "#%lu LVL1 items, #%lu HLT chains",
               db_itr->first.first,
               db_itr->first.second.first,
               db_itr->first.second.second,
               static_cast< long unsigned int >( m_lvl1NameMap()->size() ),
               static_cast< long unsigned int >( m_hltNameMap()->size() ) );
      }

      return;
   }

   /**
    * This function can be used to print quite detailed information about the currently
    * loaded trigger configuration. Mostly used for debugging.
    */
   void TrigConfigSvcD3PD::PrintConfig() {

      //
      // Make sure that the object is initialized:
      //
      if( ! Initialize() ) {
         Error( "PrintSummary", "Couldn't load configuration" );
         return;
      }

      //
      // Print the information about the loaded LVL1 items:
      //
      Info( "PrintConfig", "*** LVL1 configuration ***" );
      std::map< std::string, int >::const_iterator name_itr = m_lvl1NameMap()->begin();
      std::map< std::string, int >::const_iterator name_end = m_lvl1NameMap()->end();
      std::map< std::string, float >::const_iterator prescale_itr;
      for( ; name_itr != name_end; ++name_itr ) {
         if( ( prescale_itr = m_lvl1PrescaleMap()->find( name_itr->first ) ) ==
             m_lvl1PrescaleMap()->end() ) {
            Fatal( "PrintConfig",
                   "There is an internal inconsistency in the LVL1 configuration!" );
            return;
         }
         Info( "PrintConfig", "Name: %s, \tCTP ID: %i, \tPrescale: %f",
               name_itr->first.c_str(), name_itr->second, prescale_itr->second );
      }

      //
      // Print the information about the loaded HLT chains:
      //
      Info( "PrintConfig", "*** HLT configuration ***" );
      name_itr = m_hltNameMap()->begin();
      name_end = m_hltNameMap()->end();
      std::map< std::string, float >::const_iterator rprescale_itr;
      std::map< std::string, float >::const_iterator passthrough_itr;
      std::map< std::string, std::string >::const_iterator lowerChain_itr;
      for( ; name_itr != name_end; ++name_itr ) {
         if( ( ( prescale_itr = m_hltPrescaleMap()->find( name_itr->first ) ) ==
               m_hltPrescaleMap()->end() ) ||
             ( ( rprescale_itr = m_hltRerunPrescaleMap()->find( name_itr->first ) ) ==
               m_hltRerunPrescaleMap()->end() ) ||
             ( ( passthrough_itr = m_hltPassthroughMap()->find( name_itr->first ) ) ==
               m_hltPassthroughMap()->end() ) ||
             ( ( lowerChain_itr = m_lowerChainNameMap()->find( name_itr->first ) ) ==
               m_lowerChainNameMap()->end() ) ) {
            Fatal( "PrintConfig",
                   "There is an internal inconsistency in the HLT configuration!" );
            return;
         }
         Info( "PrintConfig", "Name: %s, \tChain ID: %i, \tPrescale: %f, "
               "\tRerun-prescale: %f, \tPassthrough: %f, "
               "\tLower chain name: %s",
               name_itr->first.c_str(), name_itr->second,
               prescale_itr->second, rprescale_itr->second,
               passthrough_itr->second,
               lowerChain_itr->second.c_str() );
      }

      return;
   }

   /**
    * This function is used internally to attach the different branches of the trigger
    * configuration metadata TTree to the internal pointers of the object. It also
    * maps the available configurations for faster access. In the end it sets the object
    * to describe the first configuration available in the metadata TTree.
    *
    * @returns <code>kTRUE</code> if successful,
    *          <code>kFALSE</code> otherwise
    */
   Bool_t TrigConfigSvcD3PD::Initialize() {

      if( m_initialized ) return kTRUE;

      m_configKeyMap.clear();

      // Check if a metadata tree has already been defined:
      if( ! m_tree ) {
         Error( "Initialize", "No trigger configuration metadata TTree specified!" );
         return kFALSE;
      }

      // Check if we actually received a TChain, not a TTree:
      TChain* chain = dynamic_cast< TChain* >( m_tree );
      if( chain ) {
         // Just to make sure that the first TTree is loaded...
         if( ! chain->GetTree() ) {
            chain->LoadTree( 0 );
         }
      }

      //
      // Connect the most basic variables:
      //
      m_smk.ReadFrom( m_tree );
      m_l1psk.ReadFrom( m_tree );
      m_hltpsk.ReadFrom( m_tree );

      //
      // Check whether these exist in the TTree:
      //
      if( ( ! m_smk.IsAvailable() ) ||
          ( ! m_l1psk.IsAvailable() ) ||
          ( ! m_hltpsk.IsAvailable() ) ) {
         Error( "Initialize", "The supplied TTree doesn't seem to be a trigger "
                "configuration metadata tree!" );
         return kFALSE;
      }

      //
      // Index the configuration entries:
      //
      for( Long64_t i = 0; i < m_tree->GetEntries(); ++i ) {

         // If the user specified a TChain as the metadata tree, then we have to be
         // a bit smarter:
         if( chain ) {
            m_entry = chain->LoadTree( i );
            if( m_entry < 0 ) break;
         }
         // In case it's a simple TTree, then we don't have to be all that careful:
         else {
            m_entry = i;
         }

         DBKeys_t keys = MakeDBKeys( m_smk(), m_l1psk(), m_hltpsk() );
         m_configKeyMap[ keys ] = i;
      }

      //
      // Reset the branches to the first entry:
      //
      if( chain ) {
         chain->LoadTree( 0 );
      }
      m_entry = 0;

      //
      // Connect the LVL1 branches:
      //
      m_lvl1NameMap.ReadFrom( m_tree );
      m_lvl1PrescaleMap.ReadFrom( m_tree );

      //
      // Check whether the LVL1 branches exist:
      //
      if( ( ! m_lvl1NameMap.IsAvailable() ) ||
          ( ! m_lvl1PrescaleMap.IsAvailable() ) ) {
         Error( "Initialize", "The supplied TTree doesn't hold the LVL1 branches" );
         Error( "Initialize", "This is VERY strange" );
         return kFALSE;
      }

      //
      // Connect the basic HLT branches:
      //
      m_hltNameMap.ReadFrom( m_tree );
      m_hltPrescaleMap.ReadFrom( m_tree );

      //
      // Check whether the basic HLT branches exist:
      //
      if( ( ! m_hltNameMap.IsAvailable() ) ||
          ( ! m_hltPrescaleMap.IsAvailable() ) ) {
         Error( "Initialize", "The supplied TTree doesn't hold the HLT branches" );
         Error( "Initialize", "This is VERY strange" );
         return kFALSE;
      }

      //
      // Connect the extra HLT branches:
      //
      m_hltRerunPrescaleMap.ReadFrom( m_tree );
      m_hltPassthroughMap.ReadFrom( m_tree );
      m_lowerChainNameMap.ReadFrom( m_tree );

      //
      // Check whether these extra branches exist:
      //
      if( ( ! m_hltRerunPrescaleMap.IsAvailable() ) ||
          ( ! m_hltPassthroughMap.IsAvailable() ) ||
          ( ! m_lowerChainNameMap.IsAvailable() ) ) {
         Warning( "Initialize", "The supplied TTree doesn't have some HLT branches" );
         Warning( "Initialize", "This is expected for slightly older D3PD files" );
      }

      //
      // Update the cached trigger names:
      //
      if( ! UpdateTriggerNames() ) {
         Error( "Initialize", "Failed to cache the configured trigger names" );
         return kFALSE;
      }

      Info( "Initialize",
            "Loaded configuration: SMK: %i, L1PSK: %i, HLTPSK: %i",
            m_smk(), m_l1psk(), m_hltpsk() );

      m_initialized = kTRUE;
      return kTRUE;
   }

   /**
    * This function is used internally to clear the memory allocated by the object.
    */
   void TrigConfigSvcD3PD::Reset() {

      m_initialized = kFALSE;

      m_smk.ReadFrom( 0 );
      m_l1psk.ReadFrom( 0 );
      m_hltpsk.ReadFrom( 0 );

      m_lvl1NameMap.ReadFrom( 0 );
      m_lvl1PrescaleMap.ReadFrom( 0 );

      m_hltNameMap.ReadFrom( 0 );
      m_hltPrescaleMap.ReadFrom( 0 );
      m_hltRerunPrescaleMap.ReadFrom( 0 );
      m_hltPassthroughMap.ReadFrom( 0 );

      m_lowerChainNameMap.ReadFrom( 0 );

      return;
   }

   /**
    * This function is used internally to update the variables caching the names of
    * all the LVL1 trigger items and HLT trigger chains.
    *
    * @returns <code>kTRUE</code> if the operation was successful, or
    *          <code>kFALSE</code> otherwise
    */
   Bool_t TrigConfigSvcD3PD::UpdateTriggerNames() {

      // Do a small sanity check:
      if( ( ! m_lvl1NameMap.IsAvailable() ) || ( ! m_hltNameMap.IsAvailable() ) ) {
         Error( "UpdateTriggerNames", "Basic trigger config branches not configured yet!" );
         return kFALSE;
      }

      // Clear the currently cached information:
      m_confLVL1Items.clear();
      m_confHLTChains.clear();
      m_lvl1IDMap.clear();
      m_lvl2IDMap.clear();
      m_efIDMap.clear();

      // Cache the names of the LVL1 items:
      std::map< std::string, int >::const_iterator itr = m_lvl1NameMap()->begin();
      std::map< std::string, int >::const_iterator end = m_lvl1NameMap()->end();
      for( ; itr != end; ++itr ) {
         m_confLVL1Items.push_back( itr->first );
         m_lvl1IDMap[ itr->second ] = itr->first;
      }

      // Cache the names of the HLT chains:
      itr = m_hltNameMap()->begin();
      end = m_hltNameMap()->end();
      for( ; itr != end; ++itr ) {
         m_confHLTChains.push_back( itr->first );
         if( itr->first.find( "L2_" ) == 0 ) {
            m_lvl2IDMap[ itr->second ] = itr->first;
         } else if( itr->first.find( "EF_" ) == 0 ) {
            m_efIDMap[ itr->second ] = itr->first;
         } else {
            Error( "UpdateTriggerNames",
                   "Can't decide if chain \"%s\" is a LVL2 or EF chain",
                   itr->first.c_str() );
            return kFALSE;
         }
      }

      return kTRUE;
   }

   /**
    * This is just a convenience function for generating DBKeys_t objects quickly.
    *
    * @param SMK Super Master Key
    * @param L1PSK LVL1 Prescale Key
    * @param HLTPSK HLT Prescale Key
    * @returns A DBKeys_t object describing the specified configuration
    */
   TrigConfigSvcD3PD::DBKeys_t
   TrigConfigSvcD3PD::MakeDBKeys( Int_t SMK, Int_t L1PSK, Int_t HLTPSK ) {

      return std::make_pair( SMK, std::make_pair( L1PSK, HLTPSK ) );
   }

} // namespace D3PD
