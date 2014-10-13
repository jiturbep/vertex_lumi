// Dear emacs, this is -*- c++ -*-
// $Id: IConfigAccess.h 469517 2011-11-21 13:06:16Z krasznaa $
#ifndef TRIGROOTANALYSIS_ICONFIGACCESS_H
#define TRIGROOTANALYSIS_ICONFIGACCESS_H

namespace D3PD {

   // Forward declaration(s):
   class TrigConfigSvcD3PD;

   namespace Trig {

      /**
       *  @short Interface providing access to the configuration service
       *
       *         In order to be able to split up the functionality of the TDT into
       *         multiple parts, the separate parts communicate with each other
       *         through such interfaces. This interface gives access to all the other
       *         components to the trigger configuration service object.
       *
       * @author Attila Krasznahorkay <Attila.Krasznahorkay@cern.ch>
       *
       * $Revision: 469517 $
       * $Date: 2011-11-21 13:06:16 +0000 (Mon, 21 Nov 2011) $
       */
      class IConfigAccess {

      public:
         /// Virtual destructor to make vtable happy
         virtual ~IConfigAccess() {}

         /// Function giving access to the trigger configuration service object
         virtual D3PD::TrigConfigSvcD3PD& GetConfigSvc( ::Bool_t update = kTRUE ) = 0;

      }; // class IConfigAccess

   } // namespace Trig

} // namespace D3PD

#endif // TRIGROOTANALYSIS_ICONFIGACCESS_H
