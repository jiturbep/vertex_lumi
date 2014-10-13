// Dear emacs, this is -*- c++ -*-
// $Id$
#ifndef TRIGROOTANALYSIS_CONFIGACCESS_H
#define TRIGROOTANALYSIS_CONFIGACCESS_H

// ROOT include(s):
#include <TNamed.h>

// Local include(s):
#include "IConfigAccess.h"
#include "IDataAccess.h"
#include "TrigConfigSvcD3PD.h"

// Forward declaration(s):
class TTree;

namespace D3PD {

   namespace Trig {

      /**
       *  @short Implementation of the IConfigAccess interface
       *
       *         This class implements a functionality of the TDT that's used by
       *         multiple other components of the class.
       *
       *         Note that the implementation needs to access the
       *         IDataAccess interface, which makes the class structure a
       *         bit complicated...
       *
       * @author Attila Krasznahorkay <Attila.Krasznahorkay@cern.ch>
       *
       * $Revision$
       * $Date$
       */
      class ConfigAccess : public virtual ::TNamed,
                           public virtual D3PD::Trig::IConfigAccess,
                           public virtual D3PD::Trig::IDataAccess {

      public:
         /// Constructor with the configuration tree
         ConfigAccess( ::TTree* confTree );

         /// Return a reference to the trigger configuration service
         virtual D3PD::TrigConfigSvcD3PD& GetConfigSvc( ::Bool_t update = kTRUE );

      private:
         D3PD::TrigConfigSvcD3PD m_configSvc; ///< The trigger configuration service

         ClassDef( D3PD::Trig::ConfigAccess, 0 )

      }; // class ConfigAccess

   } // namespace Trig

} // namespace D3PD

#endif // TRIGROOTANALYSIS_CONFIGACCESS_H
