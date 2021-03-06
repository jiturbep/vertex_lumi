#ifndef AnaVtxTree_h
#define AnaVtxTree_h

// Remark A. Wildauer: when recreating the base class with root MakeClass all you need to do is to add
// using namespace std;
// to the header.

//Special flag to include MC truth
#include "VtxTree.h"

//Include useful definition to be used also in other macros (ex. drawing)
#include "CommonTypes.h"

//Include truth-matching algorithm
#include "VtxTmAlg.h"

//C++ includes
#include <fstream>
#include <map>
#include <limits>
#include <string>
#include <sstream>

//ROOT includes
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"
#include "TChain.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStopwatch.h"

//Trigger Decision Tool
#include "TrigRootAnalysis/TrigDecisionToolD3PD.h"

//My includes
#include "ETypes.h"

class AnaVtxTree : public VtxTree {
  public:
    // -- Const, enum
    // - Categorization based on distance reco-generated
    enum TagVtxEfficiencyType {
      TagVtxCorrect=0, // only one vertex within threshold, correctly tagged
      TagVtxSplitCorrect=1, // >1 vertices within threshold, right one tagged
      TagVtxSplitWrong=2, // >1 vertices within threshold, wrong one tagged
      TagVtxWrong=3, // >=1 vertex within threshold, but not the tagged one
      TagVtxWrongFake=4, // no vertexes within threshold
      TagVtxNumTypes=5
    };

    enum VtxZMatch {
      VtxZM_Match, //1-1 association
      VtxZM_Merge, //2-1
      VtxZM_Split, //1-2
      VtxZM_Fake,  //0-1
      VtxZM_Ineff, //1-0
      VtxZM_Match_2, //2-2 association
      VtxZM_Others, //other combinations
      VtxZM_NMatch //Number of possibilities
    };

    //Vertex distance metric
    enum VtxDistanceMetric {
      VtxDistM_deltaZ, // Delta Z
      VtxDistM_deltaZsig, // Delta Z / sigma(Delta Z)
      VtxDistM_3Dsig // Quadratic sum of Delta x_i/Sigma(Delta x_i) (i=x,y,z)
    };

    enum VtxType {
      VtxSimulated, // MC vertex
      VtxReconstructedNoBC, // non-BC reco vertex
      VtxReconstructedBC, // BC reco vertex
      VtxReconstructedSplit // Split vertex
    };
    float cfgMaxDistance;
    VtxDistanceMetric cfgMetric;

    enum VtxTmResult {
      TmReal,
      TmFake, // Vertex doesn't have at least (NTrkCut) tracks from a single interaction
      TmSplit, // Vertex shares dominant generated interaction with a higher SumPt2 vertex
      TmReject, // Vertex doesn't pass selection
      TmUnknown
    };

  protected:

    // -- Histograms
    static const Int_t NBinsZAxis = 600;
    static const Float_t LowZAxis = -300.;
    static const Float_t HighZAxis = 300.;
    static const Int_t NBinsXAxis = 400;
    static const Float_t LowXAxis = -2.;
    static const Float_t HighXAxis = 2.;
    static const Int_t NBinsYAxis = 400;
    static const Float_t LowYAxis = -1.;
    static const Float_t HighYAxis = 3.;
    static const Int_t MaxNGenInt = 120; //ngenint goes to up to 114, actualIntPerXing goes up to 71
    static const Int_t MaxNVtxRecon = 100;
    vector<Int_t> nTrkCuts;

    TH1F* h_NTrig_NGenInt;
    std::map<Int_t, TH1F*> h_NTrig_NVtxRecon;
    TH2F* h_actualInt_NGenInt;
    TH2F* h_actualInt_mcvtxn;

    // Fake correction
    std::map<Int_t, TH1F*> h_fakes_NGenInt;
    std::map<Int_t, TH1F*> h_truefakes_NGenInt;
    std::map<Int_t, TH1F*> h_splits_NGenInt;
    std::map<Int_t, TH1F*> h_real_NGenInt;
    std::map<Int_t, TH1F*> h_all_NGenInt;
    std::map<Int_t, TH1F*> h_all_actualInt;

    std::map<Int_t, TH1F*> h_fakes_NVtxRecon;
    std::map<Int_t, TH1F*> h_truefakes_NVtxRecon;
    std::map<Int_t, TH1F*> h_splits_NVtxRecon;
    std::map<Int_t, TH1F*> h_real_NVtxRecon;
    std::map<Int_t, TH1F*> h_all_NVtxRecon;

    std::map<Int_t, TH1F*> h_fake_events_NGenInt;
    std::map<Int_t, TH1F*> h_all_events_NGenInt;
    std::map<Int_t, TH1F*> h_fake_events_NVtxRecon;
    std::map<Int_t, TH1F*> h_all_events_NVtxRecon;

    std::map<Int_t, TH1F*> h_event_contains_fake_NGenInt;
    std::map<Int_t, TH1F*> h_event_contains_any_NGenInt;
    std::map<Int_t, TH1F*> h_event_contains_fake_NVtxRecon;
    std::map<Int_t, TH1F*> h_event_contains_any_NVtxRecon;

    std::map<Int_t, TProfile*> h_GoodVertices_actualInt;
    std::map<Int_t, TH2F*> h_GoodVertices_actualInt_2D;
    std::map<Int_t, TProfile*> h_GoodVertices_NGenInt;
    std::map<Int_t, TH2F*> h_GoodVertices_NGenInt_2D;
    std::map<Int_t, TProfile*> h_GoodVertices_mcvtxn;
    std::map<Int_t, TH2F*> h_GoodVertices_mcvtxn_2D;
    std::map<Int_t, TProfile*> h_NVtxRecon_actualInt;
    std::map<Int_t, TProfile*> h_NVtxRecon_NGenInt;
    std::map<Int_t, TProfile*> h_NReal_actualInt;
    std::map<Int_t, TProfile*> h_NReal_NGenInt;
    std::map<Int_t, TProfile*> h_NReal_NVtxRecon;


    // Masking correction
    std::map<Int_t, TH2F*> h_all_dz_NGenInt;
    std::map<Int_t, TH2F*> h_all_dz_NVtxRecon;
    TH2F *h_privtx_z_mu; // Primary vertex only: x axis = mu, y axis = PosZ
    TH2F *h_privtx_z_NGenInt;
    TH2F *h_privtx_z_NVtxRecon;

    // Split studies
    std::map<Int_t, TH2F*> h_split_dz_NGenInt;
    std::map<Int_t, TH2F*> h_split_dzsig_NGenInt;
    std::map<Int_t, TH2F*> h_split_dz_dzsig;
    std::map<Int_t, TH2F*> h_real_dz_NGenInt;
    std::map<Int_t, TH2F*> h_real_dzsig_NGenInt;
    std::map<Int_t, TH2F*> h_real_dz_dzsig;

    std::map<Int_t, TH2F*> h_split_dz_NVtxRecon;
    std::map<Int_t, TH2F*> h_split_dzsig_NVtxRecon;
    std::map<Int_t, TH2F*> h_real_dz_NVtxRecon;
    std::map<Int_t, TH2F*> h_real_dzsig_NVtxRecon;

    // Debug plots
    std::map<Int_t, TH1F*> h_fakes_x;
    std::map<Int_t, TH1F*> h_fakes_y;
    std::map<Int_t, TH1F*> h_fakes_z;
    std::map<Int_t, TH1F*> h_fakes_cov_xx;
    std::map<Int_t, TH1F*> h_fakes_cov_yy;
    std::map<Int_t, TH1F*> h_fakes_cov_zz;
    std::map<Int_t, TH1F*> h_fakes_cov_xy;
    std::map<Int_t, TH1F*> h_fakes_cov_xz;
    std::map<Int_t, TH1F*> h_fakes_cov_yz;
    std::map<Int_t, TH1F*> h_fakes_chi2ndf;
    std::map<Int_t, TH1F*> h_fakes_chi2;
    std::map<Int_t, TH1F*> h_fakes_ndf;
    std::map<Int_t, TH1F*> h_fakes_ntrk;
    std::map<Int_t, TH1F*> h_fakes_sumPt;

    std::map<Int_t, TH1F*> h_splits_x;
    std::map<Int_t, TH1F*> h_splits_y;
    std::map<Int_t, TH1F*> h_splits_z;
    std::map<Int_t, TH1F*> h_splits_cov_xx;
    std::map<Int_t, TH1F*> h_splits_cov_yy;
    std::map<Int_t, TH1F*> h_splits_cov_zz;
    std::map<Int_t, TH1F*> h_splits_cov_xy;
    std::map<Int_t, TH1F*> h_splits_cov_xz;
    std::map<Int_t, TH1F*> h_splits_cov_yz;
    std::map<Int_t, TH1F*> h_splits_chi2ndf;
    std::map<Int_t, TH1F*> h_splits_chi2;
    std::map<Int_t, TH1F*> h_splits_ndf;
    std::map<Int_t, TH1F*> h_splits_ntrk;
    std::map<Int_t, TH1F*> h_splits_sumPt;

    std::map<Int_t, TH1F*> h_real_x;
    std::map<Int_t, TH1F*> h_real_y;
    std::map<Int_t, TH1F*> h_real_z;
    std::map<Int_t, TH1F*> h_real_cov_xx;
    std::map<Int_t, TH1F*> h_real_cov_yy;
    std::map<Int_t, TH1F*> h_real_cov_zz;
    std::map<Int_t, TH1F*> h_real_cov_xy;
    std::map<Int_t, TH1F*> h_real_cov_xz;
    std::map<Int_t, TH1F*> h_real_cov_yz;
    std::map<Int_t, TH1F*> h_real_chi2ndf;
    std::map<Int_t, TH1F*> h_real_chi2;
    std::map<Int_t, TH1F*> h_real_ndf;
    std::map<Int_t, TH1F*> h_real_ntrk;
    std::map<Int_t, TH1F*> h_real_sumPt;

    // Truth info
    TH1F* h_truth_x;
    TH1F* h_truth_y;
    TH1F* h_truth_z;

    //Truth-Matching
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

    TH1F *h_MCVrt_RecoSimu_VtxM[VtxZM_NMatch];
    TH1F *h_MCVrt_Num_SimuVtxGood; //events as function of simulated good vertices
    TH1F *h_MCVrt_Num_RecoVtxGood; // events as function of reconstructed good vertices
    TH1F *h_MCVrt_RecoSimu_VtxZM_Dist_Match;
    TH1F *h_MCVrt_RecoSimu_VtxZM_Dist_Split_gen_1;
    TH1F *h_MCVrt_RecoSimu_VtxZM_Dist_Split_gen_2;
    TH1F *h_MCVrt_RecoSimu_VtxZM_Dist_Split_rec;
    TH1F *h_MCVrt_RecoSimu_VtxZM_Dist_Merge_firstGen;
    TH1F *h_MCVrt_RecoSimu_VtxZM_Dist_Merge_secondGen;
    TH1F *h_MCVrt_RecoSimu_VtxZM_Dist_Fake;
    TH1F *h_MCVrt_RecoSimu_VtxZM_Dist_Ineff;

    //Truth-matching studies
    TH1F *h_TM_splitDeltaZ; ///< Delta-Z between main and split vertices
    TH1F *h_TM_splitDeltaZSig; ///< Delta-Z/sigma(Delta-Z) between main and split vertices
    /** Study hard-scattering interaction (hs).
     * Bin 0: total number of events
     * Bin 1: events where a vertex is matched to hs
     * Bin 2: events where the vertex matched to hs is merged
     * Bin 3: Number of split vertices which originate (main contribution) from hs
     * Bin 4: events where the matched vertex to hs (bin 1) is tagged as the PV (first in the collection, higher sumpt2)
     */
    TH1F *h_TM_class_hs;
    /** Study pile-up interactions (pu).
     * Bin 0: total number of pile-up interactions
     * Bin 1: pu vertices with a vertex matched
     * Bin 2: pu vertices that are a merge of pu interactions
     * Bin 3: pu vertices that are a split from a pu interaction
     */
    TH1F *h_TM_class_pu;


    // -- counters
    Int_t m_TotEvents; // Total events
    Int_t m_TotRawEvents; //before GRL requirement
    Int_t m_VtxRec; //Number of Vx reconstructed
    Int_t m_VtxRecPri; //Number of Primary Vx reconstructed
    Int_t m_VtxRecSec; //Number of Primary Vx reconstructed
    Int_t m_VtxRecPU; //Number of Primary Vx reconstructed
    Int_t m_VtxSim; //Number of Vx simulated
    Int_t m_VtxGoodSim; //Number of good Vx simulated
    Int_t m_NoVtxButTracks; //Events with no Tagged Vertex but with selected tracks
    Int_t m_VtxTagRec; //Number of tagged Vx reconstructed
    Int_t m_VtxTagZSel; //Number of tagged Vx Z-selected
    Int_t m_VtxTagPurPos; //Number of tagged Vx with positive purity
    Int_t m_VtxAtLease2Trk; // Number of events with at least 2 tracks good for Vx finding
    Int_t m_VtxOnlyTagRecSim1;
    Int_t m_VtxTagRecSim1;
    Int_t m_VtxTagZSel_wrtAnyGenVtx;
    Int_t m_VtxTagThreshold[TagVtxNumTypes];
    Int_t m_VtxIDEfficiency; ///< Identification of PV efficiency
    // Monitoring / text output
    std::map<Int_t, Float_t> m_vtx_total;
    std::map<Int_t, Float_t> m_vtx_fake, m_vtx_real, m_vtx_split;
    Float_t m_ngenint;
    std::map<Int_t, Float_t> m_efficiency; // m_vtx_total[*nTrkCut] / m_ngenint


    // -- Compute numeric results
    void computeStats();
    void computeMCOnlyStats();
    void computeMeanRms(TH1F *h_mean, TH1F *h_rms, TH1F *h_n);
    void computeMeanRms(TH1F *h_mean, TH1F *h_rms, TH2F *h_2dplot);
    void computeMeanRmsReverseAxis(TH1F *h_mean, TH1F *h_rms, TH2F *h_2dplot);
    void computeMeanRmsFit(TH1F *h_mean, TH1F *h_rms, TH2F *h_2dplot);
    Double_t m_AvgVtxRec;
    Double_t m_AvgVtxSim;
    Double_t m_VtxTagRecoEff;
    Double_t m_VtxTagRecoEffErr;
    Double_t m_VtxRecoEffPU;
    Double_t m_VtxRecoEffPUErr;
    Double_t m_VtxTagZSelEff;
    Double_t m_VtxTagZSelEffErr;
    Double_t m_VtxTagZEff;
    Double_t m_VtxTagZEffErr;
    Double_t m_VtxTagPurSelEff;
    Double_t m_VtxTagPurSelEffErr;
    Double_t m_VtxFakeReco;
    Double_t m_VtxFakeRecoErr;
    Double_t m_VtxTagZSelGenAnyEff;
    Double_t m_VtxTagZSelGenAnyEffErr;
    ETypes::EDouble m_VtxTagTSelEff[TagVtxNumTypes];
    //Categorization based on vertex->tracks->truth-particle->gen-vertex matching
    VtxTmAlg TmVtxNBc; ///< Truth-matching of non beam-constrained vertices

    //other results
    Double_t m_VtxResX; // X resolution
    Double_t m_VtxResY; // X resolution
    Double_t m_VtxResZ; // X resolution
    Double_t m_VtxResXErr; // X resolution
    Double_t m_VtxResYErr; // X resolution
    Double_t m_VtxResZErr; // X resolution

    // -- temporary quick results Canvases
    TCanvas * c_resCanvas;
    TCanvas * c_resCanvas2;
    TLegend * l_rC3a;
    TLegend * l_rC3b;

    /// -- Fit functions
    TF1 *vtxFitd0phiFcn;

    // -- Output
    TFile *outputNtp;
    TString outputFileName;
    ofstream outputTxt;

    // -- Class Parameters
    Float_t VtxEffCriteriaZ;
    Float_t VtxEffCriteriaT;
    bool performD0PhiFit;
    /** Debug Level.
     * Set higher positive numbers for increasing debug.
     * Set specific negative number for enabling a particular sub-system debug:
     * -12,-13: z-distance truth vertex matching
     * -11: comupte RMS/mean
     * -14: z-distance class matching
     * -15: new truth-matching (based on track matching)
     * -17: write out vertex info to TTree if dz < 5mm.
     */
    Int_t debugLevel;
    Int_t qualityVertexVersion;
    Int_t qualityVertexModifier;
    Int_t qualityVertexParameter;
    //vertex truth-matching parameters
    float vtxTMWThreshold;
    float vtxTMWThresholdStore;
    float vtxTMTrackMatchProb;

    // -- Helper functions
    /*TH1F* defineHistogram(TString name,
                          unsigned int bins, Double_t lbin, Double_t ubin,
                          TString xaxistitle, TString yaxistitle,
                          bool save=true, TString save_name="");
    TH1I* defineIHistogram(TString name,
                           unsigned int bins, Int_t lbin, Int_t ubin,
                           TString xaxistitle, TString yaxistitle,
                           bool save=true, TString save_name="");
    TH2F* define2DHistogram(TString name, unsigned
                            int xbins, Double_t lxbin, Double_t uxbin,
                            unsigned int ybins, Double_t lybin, Double_t uybin,
                            TString xaxistitle, TString yaxistitle,
                            bool save=true, TString save_name="");
    TH2F* define2DHistogram(TString name, unsigned
                            int xbins, Double_t lxbin, Double_t uxbin,
                            unsigned int ybins, Double_t lybin, Double_t uybin,
                            TString xaxistitle, TString yaxistitle,
                            TH1F * &h_mean, TH1F * &h_rms, //created to store mean and RMS of the 2D histogram
                            bool save=true, TString save_name="");
    TH3F* define3DHistogram(TString name,
                            unsigned int xbins, Double_t lxbin, Double_t uxbin,
                            unsigned int ybins, Double_t lybin, Double_t uybin,
                            unsigned int zbins, Double_t lzbin, Double_t uzbin,
                            TString xaxistitle, TString yaxistitle, TString zaxistitle,
                            bool save, TString save_name);*/

    bool trigger(TString triggerName);
    int isGoodVertex(int vertexToCheck=-1); //return index of tagged vertex
    int calculateEfficiency(Double_t &eff, Double_t &effErr, const Int_t num, const Int_t den);
    int calculateEfficiencyAsymErrors(Double_t &eff, Double_t &effErr_p, Double_t &effErr_m, const Int_t num, const Int_t den);
    TagVtxEfficiencyType getTagVertexClass(Int_t TagVtx, Int_t simVtx, bool forceDebug=false);
    Double_t getTDistance(Int_t VtxIdx, Int_t simVtxIdx);
    Int_t countNGenTracks(Float_t pEtaMax = 2.5, Float_t pPtMin = 150.);
    Float_t GetVtxDistance(Int_t VxIndex1, VtxType VxType1, Int_t VxIndex2, VtxType VxType2, VtxDistanceMetric cfgMetric);
    TString MakeHistoNamesCOORD(TString t, int ia); ///< Replace work 'COORD' with 'X','Y','Z' depending on ia = 0,1,2

    // -- Map to store histograms to be saved: [name, [type, TH1-object]]
    //map<std::string, pair<string, TH1* > > storageHistQueue;
    void saveHistograms();


  public:
    AnaVtxTree();
    ~AnaVtxTree() {};

    virtual Bool_t Process(Long64_t entry);
    virtual void SlaveBegin(TTree *tree);
    virtual void SlaveTerminate();

    Float_t SetEffCriteriaZ(Float_t criteria) {
      return (VtxEffCriteriaZ = criteria);
    }
    Float_t SetEffCriteriaT(Float_t criteria) {
      return (VtxEffCriteriaT = criteria);
    }
    TString SetOutputFileName(TString f) {
      return (outputFileName = f);
    };
    //Set debugLevels:
    // 1..50: increasing verbosity
    // -11: debug ComputeMeanRms function
    // -12: debug getTagVertexClass function
    // -13: debug getTagVertexClass function - Wrong vertices
    Int_t   SetDebugLevel(Int_t debug);
    //Set/Get qualityVertexVersion: x.y.z
    // rev = 0: major version (x)
    // rev = 1: modifier (y)
    // rev = 2: parameter (z)
    TString SetQualityVertexVersion(TString version);
    Int_t   GetQualityVertexVersion(Int_t rev=0);

    //Trigger filtering, need to provide also trigger meta-data tree
    void SetTriggerFilter(TString p_triggerName, TChain *p_trigMetaDataTree=0);

    //Optional filter on actual mu values
    void SetMuFilter(float p_muMin, float p_muMax);

    //Optional filter on chi2/ndf
    void SetMaxChi2Ndf(Float_t p_max_chi2ndf);

    //Set hard-scatter skip
    void SetHS(bool p_skip_hs);
    bool skip_hs;

    //Stopwatch for performance study
    void DoTiming();
    bool do_timing;

  protected:
    //Trigger decision tool
    D3PD::TrigDecisionToolD3PD* triggerTool;
    TChain *trigMetaDataTree;
    TString triggerName; ///< Trigger filtering, if not empty

    // Mu filtering
    float muFilter_min;
    float muFilter_max;

    // Cut on dz significance
    float dzsig_max;
    Float_t max_chi2ndf;
    
  private:
    
    // Helper functions to concatenate strings and doubles
    std::string str( const std::string &inStr, const double &dbl ) {
      std::stringstream combine;
      combine << inStr << dbl;
      return combine.str();
    }

    TString cstr( const std::string &inStr, const double &dbl ) {
      return str(inStr,dbl).c_str();
    }

    std::string str( const std::string &str1, const std::string &str2 ) {
      std::stringstream combine;
      combine << str1 << str2;
      return combine.str();
    }

    TString cstr( const std::string &str1, const std::string &str2 ) {
      return str(str1,str2).c_str();
    }
    
    
};

#endif
