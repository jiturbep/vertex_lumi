#include "GoodRunsLists/TGoodRunsList.h"
#include "TString.h"

namespace DQ
{
  void SetXMLFile(const TString& filename);
  Root::TGoodRunsList& GetGRL();
  bool PassRunLB(int runnr, int lbnr);
}

