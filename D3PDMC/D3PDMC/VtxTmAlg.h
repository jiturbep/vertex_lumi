/** @file Vertex truth matching algorithm */

#ifndef VtxTmAlg_h
#define VtxTmAlg_h

#include "VertexTruthMatch.h"

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
class VtxTmAlg {
  protected:
    // --- Settings
    float vtxTMWThreshold; //to define
    float vtxTMWThresholdStore;
    float vtxTMTrackMatchProb;

    /// Set requirements on "good" generated vertices
    Int_t genVertexRequirementVersion;

    int debugTM;

    // --- Helpers
    /// PDG Table
    TDatabasePDG *m_PDG;

    // --- Caching system for IsGoodGenVertex
    // Store info on all the vertices for a given run/evt number
    int m_runNumber;
    int m_runNumber_tight;
    int m_evtNumber;
    int m_evtNumber_tight;
    std::vector<bool> m_isGenVtxGood;
    std::vector<bool> m_isGenInteraction;
    std::vector<int> m_NTracksGenVtx;


    // --- Store all needed information
    // reconstructed vertices
    int vx_n;
    std::vector<float> *vx_x;
    std::vector<float> *vx_y;
    std::vector<float> *vx_z;
    std::vector<float> *vx_cov_x;
    std::vector<float> *vx_cov_y;
    std::vector<float> *vx_cov_z;
    std::vector<int> *vx_trk_n;
    std::vector<std::vector<float> > *vx_trk_weight;
    std::vector<std::vector<int> > *vx_trk_index;

    // reconstructed tracks
    int             trk_n;
    std::vector<float>   *trk_pt;
    std::vector<float>   *trk_eta;
    std::vector<float>   *trk_d0_wrtPV;
    std::vector<float>   *trk_z0_wrtPV;
    std::vector<float>   *trk_mc_probability;
    std::vector<int>     *trk_mc_index;

    // particles
    int             mcpart_n;
    std::vector<float>   *mcpart_pt;
    std::vector<float>   *mcpart_eta;
    std::vector<int>     *mcpart_mcevt_index;
    std::vector<int>     *mcpart_mcprodvtx_index;
    std::vector<int>     *mcpart_type;
    std::vector<int>     *mcpart_barcode;
    std::vector<int>     *mcpart_status;

    // genVertices
    int             mcvtx_n;
    std::vector<float>   *mcvtx_x;
    std::vector<float>   *mcvtx_y;
    std::vector<float>   *mcvtx_z;
    std::vector<int>     *mcvtx_mcevt_index;

    // genEvents
    int              mcevt_n;
    std::vector<short>   *mcevt_pileUpType;
    std::vector<int>     *mcevt_nparticle;

    // --- Debug plots
  public:
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

    /// Define debug histograms, if a valid pointer is not already set
    void InitDebugHistograms(std::string prefixName="");

  public:
    // --- Store results: matched vertices (one day will be protected)
    std::vector<VertexTruthMatch> matchedVtx;

  public:
    VtxTmAlg();
    ~VtxTmAlg();

    // --- Settings

    /// Set matching hit-based track matching probability to be used
    void SetTrackMatchProbability(float trkMatchProb);

    /// Set weight threshold for single vertex-interaction matching
    void SetVtxMatchWeight(float vtxWTh);

    /// Set weight threshold for storing interaction informations
    void SetVtxStoreWeight(float vtxWTh);

    /** Set requirements on generated vertices to be considered.
     * Depending on version the defined requirements are:
     * 1: InTimePileup
     * 2: InTimePileup, 2 charged particles with |eta| < 2.4, pT > 400 MeV (default)
     */
    void SetGenVertexRequirementVersion(Int_t version);

    /** Set debugging printout.
     * >=1 enable debug histograms
     * >=2 print results
     * >=3 print more info..
     * >=4,5 print even more info!
     */
    void SetDebug(int dbg);

    // --- Set input information
    ///Reconstructed vertices inputs (D3PD)
    void SetRecoVtxInfo(int p_vx_n,
                        std::vector<int> *p_vx_trk_n, std::vector<std::vector<float> > *p_vx_trk_weight, std::vector<std::vector<int> > *p_vx_trk_index,
                        std::vector<float> *p_vx_x, std::vector<float> *p_vx_y, std::vector<float> *p_vx_z,
                        std::vector<float> *p_vx_cov_x, std::vector<float> *p_vx_cov_y, std::vector<float> *p_vx_cov_z);

    ///Reconstructed traks inputs (D3PD)
    void SetRecoTrkInfo(int p_trk_n,
                        std::vector<float> *p_trk_mc_probability, std::vector<int> *p_trk_mc_index,
                        std::vector<float> *p_trk_pt, std::vector<float> *p_trk_eta, std::vector<float> *p_trk_d0_wrtPV, std::vector<float> *p_trk_z0_wrtPV);

    ///Generated particles inputs (D3PD)
    void SetGenPartInfo(int p_mcpart_n,
                        std::vector<int> *p_mcpart_mcevt_index, std::vector<int> *p_mcpart_mcprodvtx_index, std::vector<int> *p_mcpart_type,
                        std::vector<float> *p_mcpart_pt, std::vector<float> *p_mcpart_eta, std::vector<int> *p_mcpart_barcode, std::vector<int> *p_mcpart_status);

    ///Generated vertices inputs (D3PD)
    void SetGenVtxInfo(int p_mcvtx_n,
                       std::vector<int> *p_mcvtx_mcevt_index,
                       std::vector<float> *p_mcvtx_x, std::vector<float> *p_mcvtx_y, std::vector<float> *p_mcvtx_z);

    ///Generated events inputs (D3PD)
    void SetGenEventsInfo(int p_mcevt_n, std::vector<short> *p_mcevt_pileUpType, std::vector<int> *p_mcevt_nparticle);

    // --- Run algorithm
    void MatchVertices(int runNumber, long evtNumber);

    // --- Retrieve results
    // up to now direct access to result vector matchedVtx

    /// Return the mother vertex of a split
    int getMotherOfSplit(int idxRecoVtx);

    /// Write into current directory debug histograms
    void WriteDebugHisto();

    // --- Utility functions
    /// Apply "quality" requirement to generated vertices to consider
    int isGoodGenVertex(int vertexToCheck, int runNumber, long evtNumber); //return != 0 if gen vertex is good
    int GetNTracksGenVtx(int vertexToCheck, int runNumber, long evtNumber);
    int isGenInteraction(int vertexToCheck, int runNumber, long evtNumber);
};

#endif
