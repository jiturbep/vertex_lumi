/** @file Vertex truth matching algorithm */

#ifndef InDetVertexTruthMatch_h
#define InDetVertexTruthMatch_h

#include "VertexMatchInfo.h"

#include <vector>
#include <string>

#include "TDatabasePDG.h"
#include "TH1F.h"
#include "TH2F.h"

/** Class for vertex truth matching.
 * Match reconstructed vertices to truth level interactions vertices
 * through the chain: track -> particle -> genEvent -> genVertex
 * Categorize reconstructed vertices depending on their composition.
 */
class InDetVertexTruthMatch {

 public:

  InDetVertexTruthMatch();
  ~InDetVertexTruthMatch();

  //---------------------
  // --- Settings
  //---------------------

  /// Set matching hit-based track matching probability to be used
  void SetTrackMatchProbability(float trkMatchProb);

  /// Set weight threshold for single vertex-interaction matching
  void SetVtxMatchWeight(float vtxWTh);

  /// Set weight threshold for storing interaction informations
  void SetVtxStoreWeight(float vtxWTh);

  /** Set requirements on generated vertices to be considered.
   * Depending on version the defined requirements are:
   * 1: InTimePileup (default)
   * 2: InTimePileup, 2 charged particles with |eta| < 2.5, pT > 400 MeV 
   * 3: InTimePileup, 5 charged particles with |eta| < 2.5, pT > 400 MeV 
   */
  void SetGenVertexRequirementVersion(Int_t version);
  //or set individual cuts
  void SetGenVertexReqNPart( int npart ) { req_ngen = npart; }
  void SetGenVertexReqPt( double pt ) { req_genpt = pt; }
  void SetGenVertexReqEta( double eta ) { req_geneta = eta; }
  void SetGenVertexReqNTracks( int ntrks ) { req_ntrks = ntrks; }

  /** Set debugging printout.
   * >=1 enable debug histograms
   * >=2 print results
   * >=3 print more info..
   * >=4,5 print even more info!
   */
  void SetDebug(int dbg);

  //--------------------------
  // --- Set input information
  //--------------------------
  ///Reconstructed vertices inputs (D3PD)
  void SetRecoVtxInfo(int p_vx_n, 
		      std::vector<int> *p_vx_trk_n, std::vector<std::vector<float> > *p_vx_trk_weight, std::vector<std::vector<int> > *p_vx_trk_index);//,
// 		      std::vector<float> *p_vx_x, std::vector<float> *p_vx_y, std::vector<float> *p_vx_z,
// 		      std::vector<float> *p_vx_cov_x, std::vector<float> *p_vx_cov_y, std::vector<float> *p_vx_cov_z);
  
  ///Reconstructed traks inputs (D3PD)
  void SetRecoTrkInfo(int p_trk_n,
		      std::vector<float> *p_trk_mc_probability, std::vector<int> *p_trk_mc_index,
		      std::vector<float> *p_trk_pt); //, std::vector<float> *p_trk_eta, std::vector<float> *p_trk_d0_wrtPV, std::vector<float> *p_trk_z0_wrtPV);
  
  ///Generated particles inputs (D3PD)
  void SetGenPartInfo(int p_mcpart_n,
		      std::vector<int> *p_mcpart_mcevt_index,// std::vector<int> *p_mcpart_mcprodvtx_index,
		      std::vector<int> *p_mcpart_type,
		      std::vector<float> *p_mcpart_pt, std::vector<float> *p_mcpart_eta, std::vector<int> *p_mcpart_barcode, std::vector<int> *p_mcpart_status);
  
  ///Generated vertices inputs (D3PD)
  void SetGenVtxInfo(int p_mcvtx_n,
		     std::vector<int> *p_mcvtx_mcevt_index);//,
  //		     std::vector<float> *p_mcvtx_x, std::vector<float> *p_mcvtx_y, std::vector<float> *p_mcvtx_z);
  
  ///Generated events inputs (D3PD)
  void SetGenEventsInfo(int p_mcevt_n, std::vector<short> *p_mcevt_pileUpType, std::vector<int> *p_mcevt_nparticle);
  
  //---------------------
  // --- Run algorithm
  //---------------------
  //providing run/evt number will cache good gen vertex info
  void MatchVertices(int runNumber=-1, long evtNumber=-1);

  //---------------------
  // --- Retrieve results
  //---------------------

  /// Return the mother vertex of a split
  int getMotherOfSplit(int idxRecoVtx);
  // --- Store results: matched vertices (one day will be protected)
  std::vector<VertexMatchInfo> matchedVtx;

  VertexMatchInfo GetRecoVtxInfo(int recoIdx);

  //get truth match info for reco vertex this is matched or merged hard scatter vertex
  VertexMatchInfo GetHardScatterBestMatch();
  std::vector<VertexMatchInfo> GetHardScatterInfo();




  // --- Debug plots

//   TH1F * h_TM_goodUnmatched_pt;
//   TH1F * h_TM_goodUnmatched_eta;  
//   TH2F * h_TM_goodUnmatched_pteta;
//   TH1F * h_TM_highestW;
//   TH1F * h_TM_2ndHighestW;
//   TH1F * h_TM_OtherHighestW;
//   TH1F * h_TM_numGenMatched;
//   TH1F * h_TM_numGenMatchedAll;
//   TH1F * h_TM_class;
//   TH1F * h_TM_fakeRelWeight;  
//   TH1F * h_TM_splitDeltaZ; ///< Delta-Z between main and split vertices
//   TH1F * h_TM_splitDeltaZSig; ///< Delta-Z/sigma(Delta-Z) between main and split vertices

  /// Define debug histograms, if a valid pointer is not already set
//   void InitDebugHistograms(std::string prefixName="");

//   /// Write into current directory debug histograms
//   void WriteDebugHisto();

  // --- Utility functions
  /// Apply "quality" requirement to generated vertices to consider
  //  if given run and event number will cache for multiple calls
  int isGoodGenVertex(int vertexToCheck, int runNumber =-1, long evtNumber=-1); //return != 0 if gen vertex is good

 protected:
  // --- Settings
  float vtxTMWThreshold; //to define 
  float vtxTMWThresholdStore;
  float vtxTMTrackMatchProb;

  /// Set requirements on "good" generated vertices
  Int_t genVertexRequirementVersion;
  int req_ngen; //good generated parts
  float req_genpt; //pt requirement on gen parts
  float req_geneta; //eta requirement on gen parts
  int req_ntrks; //require that gen vertex also matches reco tracks


  int debugTM;

  // --- Helpers
  /// PDG Table
  TDatabasePDG *m_PDG;

  // --- Caching system for IsGoodGenVertex
  // Store info on all the vertices for a given run/evt number
  bool cacheValid; //override flag for use in match vertices when no run/evt number given
  int m_runNumber;
  int m_evtNumber;
  std::vector<bool> m_isGenVtxGood;


  // --- Store all needed information
  // reconstructed vertices
  int vx_n;
//   std::vector<float> *vx_x;
//   std::vector<float> *vx_y;
//   std::vector<float> *vx_z;
//   std::vector<float> *vx_cov_x;
//   std::vector<float> *vx_cov_y;
//   std::vector<float> *vx_cov_z;
  std::vector<int> *vx_trk_n;
  std::vector<std::vector<float> > *vx_trk_weight;
  std::vector<std::vector<int> > *vx_trk_index;
  
  // reconstructed tracks
  int             trk_n;
  std::vector<float>   *trk_pt;
//   std::vector<float>   *trk_eta;
//   std::vector<float>   *trk_d0_wrtPV;
//   std::vector<float>   *trk_z0_wrtPV;
  std::vector<float>   *trk_mc_probability;
  std::vector<int>     *trk_mc_index;

  // particles
  int             mcpart_n;
  std::vector<float>   *mcpart_pt;
  std::vector<float>   *mcpart_eta;
  std::vector<int>     *mcpart_mcevt_index;
  //  std::vector<int>     *mcpart_mcprodvtx_index;
  std::vector<int>     *mcpart_type;  
  std::vector<int>     *mcpart_barcode;
  std::vector<int>     *mcpart_status;

  // genVertices
  int             mcvtx_n;
//   std::vector<float>   *mcvtx_x;
//   std::vector<float>   *mcvtx_y;
//   std::vector<float>   *mcvtx_z;
  std::vector<int>     *mcvtx_mcevt_index;

  // genEvents
  int              mcevt_n;
  std::vector<short>   *mcevt_pileUpType;
  std::vector<int>     *mcevt_nparticle;



};


#endif
