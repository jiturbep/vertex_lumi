/**
  *  Class LumiVtx
  *  Author: David R. Yu (dryu@berkeley.edu)
  * Input: path to D3PD-based ntuple, path to timestamps
  *  Output: luminosity
*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <stdlib.h>

#include "TROOT.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TLeafI.h"
#include "TLeafD.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"

#include "GlobalSettings/GlobalSettings.h"

#include "PileupCorrections/PileupMaskingCorrection.h"
#include "PileupCorrections/FakeCorrection.h"

#ifdef USE_UNFOLD
#include "VtxLumiUnfold/VertexLumiUnfold.h"
#endif

using namespace std;

class LumiVtx {
  public:
    //Basic methods
    LumiVtx(); ///< Constructor.
    ~LumiVtx(); ///< Destructor

    void LoadTTree(TString path, TString branch_nvtx, TString branch_nevt, TString branch_ntrig, std::vector<std::pair<TString, Int_t> > *selection, bool do_physics_run = false ); /// Load LB:NVtx and LB:NEvt maps from TTrees, produced earlier from Simone's vertexing D3PD
    void LoadLivefraction(TString path);
    void LoadLivefraction();
    void LoadPrescale(TString path);
    void LoadPrescale();
    void SetSigmaVisNVtx(Float_t value);
    void SetSigmaVisNEvt(Float_t value);
    void SetSigmaVis(GlobalSettings::AnalysisMethods method, Float_t value, Float_t error=-1);
    void LoadTimestamps(TString path);
    void LoadDurations(TString path);
    void SetAllDurations(const unsigned int &LBlow, const unsigned int &LBhigh);
    void LoadDurationsFromLumiNtuple(TString path);
    void CorrectDataRate(TString triggerType);

    // - Corrections
    void InitializeMaskingCorrection(TString p_energy, TString p_settings, Int_t p_ntrkcut, TH1D *h_z, TString save_name);
    void InitializeMaskingCorrection(TString p_energy, TString p_settings, Int_t p_ntrkcut, TString save_name);
    void InitializeFakeCorrection(TString p_energy, TString p_settings,Int_t p_nTrkCut);
    void SetMuScale(Float_t p_mu_scale);
    bool correct_masking, correct_fakes;
    void CorrectPileupEffects(TH2D *h_z_plb, TString p_energy, TString p_settings, Int_t p_ntrkcut, TString p_run);
    void CorrectMasking();
    void ApplyMcCorrection(TString path, TString objectname);
    PileupMaskingCorrection *pmc;
    FakeCorrection *fc;

    // - unfolding-related
    /** Set response matrix file name and histogram base name
     * (histogram name = histoBaseName + "_NTrk" + ntrkcut
     */
    void SetUnfResponseMatrix(TString fileName, TString histoBaseName);

    /** Load vertex histograms and unfold them.
     * Store results in the appropriate members (n_unf*).
     */
#ifdef USE_UNFOLD
    void LoadVertexCountsUnfold(TString p_path, Int_t p_bcid, Int_t nTrkCut, TString triggerType);
#endif

  protected:
    // - unfolding-related

    TString responseMatrixFileName;
    TString responseMatrixHName;

#ifdef USE_UNFOLD
    /// Unfold input histogram and fills with unfolded results
    void UnfoldRawVertexCounts(TH1D *h_nVtx, Int_t p_bcid, Int_t nTrkCut, Int_t pLB, TString triggerType);

    /// Unfolding class interface
    VertexLumiUnfold *m_unfold;
#endif

  public:

    void CalculateRates(); /// Make LB:rVtx map, using LB durations.
    void CalculateAllLbLum(); /// Calculates the average-inst and total luminosity for each lumiblock.
    void CalculateAllRunLum(); /// Calculates the average-inst and total luminosity for the whole run.

    //Methods to set members manually
    void SetNBunches(Int_t p_n_bcids);

    //Methods to get members
    ////Vertex counting luminosity
    Double_t GetLbInstLum(Int_t n, GlobalSettings::AnalysisMethods method=GlobalSettings::kVtxC);
    Double_t GetLbInstLumError(Int_t n, GlobalSettings::AnalysisMethods method=GlobalSettings::kVtxC);
    Double_t GetLbTotalLum(Int_t n, GlobalSettings::AnalysisMethods method=GlobalSettings::kVtxC);
    Double_t GetLbTotalLumError(Int_t n, GlobalSettings::AnalysisMethods method=GlobalSettings::kVtxC);
    Double_t GetRunTotalLum(GlobalSettings::AnalysisMethods method=GlobalSettings::kVtxC);
    Double_t GetRunTotalLumError(GlobalSettings::AnalysisMethods method=GlobalSettings::kVtxC);

    ////Event counting luminosity
    Double_t GetEvtLbInstLum(Int_t n);
    std::pair<Double_t, Double_t> GetEvtLbInstLumError(Int_t n);
    Double_t GetEvtLbTotalLum(Int_t n);
    std::pair<Double_t, Double_t> GetEvtLbTotalLumError(Int_t n);
    Double_t GetEvtRunTotalLum();
    std::pair<Double_t, Double_t> GetEvtRunTotalLumError();

    ////Pileup masking
    Double_t GetNVtxMasked(Int_t n);
    Double_t GetNVtx(Int_t n);
    Double_t GetNEvt(Int_t n);
    Double_t GetNVtxFromNEvt(Int_t n);

    Int_t GetFirstLb();
    Int_t GetLastLb();
    Int_t GetNLb();
    Int_t GetLbFromDist(Int_t n);
    Double_t GetLbDuration(Int_t n);
    Double_t GetTotalDuration();
    Double_t GetTotalNVtxRaw();
    Double_t GetTotalNVtx();
    Double_t GetLbNVtxRaw(Int_t n);
    Double_t GetLbNVtx(Int_t n);
    Double_t GetTotalNEvtRaw(); /// Returns the raw number of NVtxTight events observed
    Double_t GetTotalNEvt(); /// Returns the scaled-up number of NVtxTight expected
    Double_t GetLbNEvtRaw(Int_t n);
    Double_t GetLbNEvt(Int_t n);
    Double_t GetTotalNTotalRaw();
    Double_t GetTotalNTotal();
    Double_t GetLbNTotalRaw(Int_t n);
    Double_t GetLbNTotal(Int_t n);
    Double_t GetTotalNVtxEvt(); /// Returns the total number of vertices as determined from analytically correcting event counts.
    Double_t GetLbNVtxEvt(Int_t n); /// Returns the total number of vertices for lumiblock LB[n], as determined from analytically correcting event counts.

    Double_t GetOflLum();
    Double_t GetOflLum(Int_t n);

    Double_t GetMaskingCorrection(Int_t n);
    Double_t GetFakeMu(Int_t, TString p_method);


    std::vector<Int_t> lumiblocks;
    std::vector<Int_t> lumiblocks_tag;
    std::vector<Int_t> lumiblocks_metadata;

    //Raw counts
    // -- n_vtx = number of tight vertices
    // -- n_evt = number of event with >=1 tight vertices
    // -- n_trig = number of events recorded
    // -- n_unf = average number of vertices before trigger/unfolding corrections
    //    0: very raw
    //    1: after trigger correction
    std::map<Int_t, Double_t> n_vtx_raw;
    std::map<Int_t, Double_t> n_vtx_raw_err;
    Double_t total_n_vtx_raw;
    Double_t total_n_vtx_raw_err;
    std::map<Int_t, Double_t> n_evt_raw;
    std::map<Int_t, Double_t> n_evt_raw_err;
    Double_t total_n_evt_raw;
    Double_t total_n_evt_raw_err;
    std::map<Int_t, Double_t> n_trig_raw;
    Double_t total_n_trig_raw;
    std::map<Int_t, Double_t> n_unf_raw[2], n_unf_raw_err[2];


    //Counts scaled up to the full values (un-prescaled, un-random-triggered, etc.).
    std::map<Int_t, Double_t> n_vtx;
    std::map<Int_t, Double_t> n_vtx_err;
    Double_t total_n_vtx, total_n_vtx_err;
    std::map<Int_t, Double_t> n_evt; // Not sure whether this will be used, the "corrected event count" doesn't really mean anything. I think.
    std::map<Int_t, Double_t> n_evt_err;
    Double_t total_n_evt;
    std::map<Int_t, Double_t> n_trig; // Number of bunch crossing per lumiblock.
    Double_t total_n_trig; // Silly name, but this is the total number of bunch crossings in the run.
    std::map<Int_t, Double_t> n_unf, n_unf_err;

    //Luminosity from vertex counting
    std::map<Int_t, Double_t> r_vtx;
    std::map<Int_t, Double_t> r_vtx_err;
    Double_t total_r_vtx;
    Double_t total_r_vtx_err;
    std::map<Int_t, Double_t> lb_inst_lum; // LB average inst lum
    std::map<Int_t, Double_t> lb_inst_lum_err;
    std::map<Int_t, Double_t> lb_total_lum;
    std::map<Int_t, Double_t> lb_total_lum_err;
    Double_t run_trig_lum; // Run total luminosity
    Double_t run_trig_lum_err;
    Double_t total_ofl_lum;

    //Luminosity from event counting
    std::map<Int_t, Double_t> r_vtx_evt;
    std::map<Int_t, Double_t> r_vtx_evt_err_up;
    std::map<Int_t, Double_t> r_vtx_evt_err_down;
    Double_t total_n_vtx_evt;
    Double_t total_n_vtx_evt_err_up;
    Double_t total_n_vtx_evt_err_down;
    //    Double_t total_r_vtx_evt;
    //      Double_t total_r_vtx_evt_err_up;
    //      Double_t total_r_vtx_evt_err_down;
    std::map<Int_t, Double_t> lb_inst_lum_evt; // LB average inst lum
    std::map<Int_t, Double_t> lb_inst_lum_evt_err_up;
    std::map<Int_t, Double_t> lb_inst_lum_evt_err_down;
    std::map<Int_t, Double_t> lb_total_lum_evt;
    std::map<Int_t, Double_t> lb_total_lum_evt_err_up;
    std::map<Int_t, Double_t> lb_total_lum_evt_err_down;
    Double_t run_trig_lum_evt; // Run total luminosity
    Double_t run_trig_lum_evt_err_up;
    Double_t run_trig_lum_evt_err_down;

    //Luminosity from unfolding
    std::map<Int_t, Double_t> lb_inst_lum_unf, lb_inst_lum_unf_err;
    std::map<Int_t, Double_t> lb_total_lum_unf, lb_total_lum_unf_err;
    Double_t run_total_lum_unf, run_total_lum_unf_err;
    Double_t run_total_ofl_lum_unf, run_total_ofl_lum_unf_err;

    //General
    std::map<Int_t, Double_t> ofl_lum; /// Offline preferred luminosity
    std::map<Int_t, Double_t> n_evt_rec, n_evt_trig; /// Number of events recorded/triggered. Not used, for now; the deadtime derives from these numbers.
    std::map<Int_t, std::pair<Double_t, Double_t> > pLB_timestamps;
    std::map<Int_t, Double_t> lb_duration; /// Lumiblock duration vs. LB
    Double_t total_duration;
    std::map<Int_t, Double_t> live_fraction; /// Live fractions vs. LB
    std::map<Int_t, Float_t> prescale; /// Prescale vs. LB

    // Corrections
    std::map<Int_t, Float_t> masking_corrections;
    std::map<Int_t, Float_t> mu_fake_list_nvtx;
    std::map<Int_t, Float_t> mu_fake_list_nevt;
    /*
    std::map<Int_t, Float_t> fake_corrections_nvtx;
    std::map<Int_t, Float_t> fake_corrections_nevt;
    */
    //Settings and stuff
    std::map<GlobalSettings::AnalysisMethods, Double_t> sigmavis, sigmavis_err;

    //Stuff for counting bunches. Note: it doesn't actually work, since it counts BCIDs recorded in the TAG, and hence might miss some!
    std::vector<Int_t> bcid_list;

    //Pileup correction
    std::map<Int_t, Double_t> n_vtx_masked;
    /* OLD STUFF
    std::map<Int_t, TGraph*> pu_corr_map;
    Float_t PileupCorrect(Float_t n_vtx_obs);
    Double_t sz; // sigma_z of the vertex distribution
    Double_t p_mask;
    */

    //Trigger mask, if a specific trigger is neeeded from the TAG.
    //trig_mask = L1, L2, or EF, trig_item will be taken %32.
    //void AddTrigMask(TString trig_type, Int_t trig_item);
    //std::vector< std::pair<TString, Int_t> > trig_mask_list;
    //bool apply_trig_mask;

    //Verbosity
    void SetVerbose(Int_t level);

    // Systematic flags
    void SetSystematicUncertaintyFlag(TString p_systematic_name);

  private:
    TFile *f_tag;
    TString ttree_path, rootfile_path, path_pu_corr;
    Int_t verbose;
    Int_t pLB_min, pLB_max, n_pLBs;
    Int_t n_bcids;
    std::map<TString, bool> systematic_uncertainty_list;
};


