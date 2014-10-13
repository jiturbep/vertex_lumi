// Dear emacs, this is -*- c++ -*-
// $Id: ChainGroupHandling.h 457181 2011-09-05 13:16:56Z krasznaa $
#ifndef TRIGROOTANALYSIS_CHAINGROUPHANDLING_H
#define TRIGROOTANALYSIS_CHAINGROUPHANDLING_H

// STL include(s):
#include <vector>
#include <string>

// ROOT include(s):
#include <TNamed.h>

// Local include(s):
#include "IDataAccess.h"
#include "IConfigAccess.h"
#include "ChainGroup.h"

namespace D3PD {

   namespace Trig {

      /**
       *  @short Class implementing the functions that handle ChainGroup-s
       *
       *         The functionality of the TDT is split up in many pieces. This
       *         piece takes care of creating ChainGroup objects in the correct
       *         way.
       *
       * @author Attila Krasznahorkay <Attila.Krasznahorkay@cern.ch>
       *
       * $Revision: 457181 $
       * $Date: 2011-09-05 14:16:56 +0100 (Mon, 05 Sep 2011) $
       */
      class ChainGroupHandling : public virtual ::TNamed,
                                 public virtual D3PD::Trig::IDataAccess,
                                 public virtual D3PD::Trig::IConfigAccess {

      public:
         /// Default constructor
         ChainGroupHandling();

         /// Create a chain group from a list of patterns
         D3PD::ChainGroup GetChainGroup( const std::vector< std::string >& patterns );
         /// Create a chain group from a single pattern
         D3PD::ChainGroup GetChainGroup( const std::string& pattern );

         ClassDef( D3PD::Trig::ChainGroupHandling, 0 )

      }; // class ChainGroupHandling

   } // namespace Trig

} // namespace D3PD

#endif // TRIGROOTANALYSIS_CHAINGROUPHANDLING_H
