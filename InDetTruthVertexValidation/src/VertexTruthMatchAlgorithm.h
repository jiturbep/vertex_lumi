#ifndef VTXTMEXAMPLE
#define VTXTMEXAMPLE

// header for AthAlgorithm
#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "InDetTruthVertexValidation/IVertexTruthMatchAthAlgTool.h"

class ITHistSvc;

class VertexTruthMatchAlgorithm : public AthAlgorithm {

public:
  // Default constructor format for AthAlgorithm
  VertexTruthMatchAlgorithm(const std::string& name, ISvcLocator* pSvcLocator);
  virtual ~VertexTruthMatchAlgorithm();

  // Required methods for all Athena algs
  virtual StatusCode initialize();
  virtual StatusCode execute();
  virtual StatusCode finalize();

private:
  ITHistSvc *m_tHistSvc;

  ToolHandle< IVertexTruthMatchAthAlgTool > m_VtxTmHandle;

  //debug histos
  TH1F * h_TM_goodUnmatched_pt;
  TH1F * h_TM_goodUnmatched_eta;  
  TH2F * h_TM_goodUnmatched_pteta;
  TH1F * h_TM_highestW;
  TH1F * h_TM_2ndHighestW;
  TH1F * h_TM_OtherHighestW;
  TH1F * h_TM_numGenMatched;
  TH1F * h_TM_numGenMatchedAll;
  TH1F * h_TM_class;
  TH1F * h_TM_fakeRelWeight;  
  TH1F * h_TM_splitDeltaZ; ///< Delta-Z between main and split vertices
  TH1F * h_TM_splitDeltaZSig; ///< Delta-Z/sigma(Delta-Z) between main and split vertices


};


#endif
