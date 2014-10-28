#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <stdlib.h>


//Compile as standalone executable
//adds specific includes here
#include "TROOT.h"
#include "TF1.h"
#include "TString.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TLeafI.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TRandom.h"
#include "TFitResult.h"
#include "TMatrixT.h"
#include "TMatrixD.h"
#include "TMatrixTSym.h"
#include "TVectorD.h"
#include "TVector.h"
#include "TSpline.h"

#include "VdM/vdMScanData.h"
#include "PileupCorrections/PileupMaskingCorrection.h"
#include "PileupCorrections/FakeCorrection.h"
#ifdef USE_UNFOLD
#include "VtxLumiUnfold/VertexLumiUnfold.h"
#endif

class VanDerMeerAnalysis {
  public:

    void ConvertToRandom(Int_t bunches, Float_t rate);

    VanDerMeerAnalysis(TString p_vtx_method, Int_t p_ntrkcut, TString p_save_tag);
    ~VanDerMeerAnalysis();
    TString vtx_method;
    Int_t ntrkcut;
    //Int_t prescale;
    TString trigger_type;
    TString save_tag;
	
    void SetScanEndpoints(Float_t p_low, Float_t p_high);
    //void SetPrescale(Int_t p_ps);
    void LoadVertexCounts(TString p_path, Int_t p_bcid, TString p_trigger_type);
    void LoadPlbTimestamps(TString p_path, TString axis);
    void LoadPlbPrescales(TString p_path);
    void LoadDeadtime(TString p_path, TString p_path_tmp = "");
    void LoadBunchIntensities(TString p_path, Int_t bcid);
    void LoadNominalSeparations(TString p_path);

    //unfolding-related
    /** Set response matrix file name and histogram base name
     * (histogram name = histoBaseName + "_NTrk" + ntrkcut
     */
    void SetUnfResponseMatrix(TString fileName, TString histoBaseName);

  protected:
    TString responseMatrixFileName;
    TString responseMatrixHName;

#ifdef USE_UNFOLD
    /// Unfold input histograms and fills nvtx_pLB with unfolded results
    void LoadVertexCountsUnfold(TString p_path, Int_t p_bcid);
    void UnfoldRawVertexCounts(TH1D *h_nVtx, Int_t p_bcid, Int_t nTrkCut, Int_t pLB);

    VertexLumiUnfold *m_unfold;
#endif

  public:

    std::map<Int_t, Double_t> nvtx_pLB;
    std::map<Int_t, Double_t> nvtx_err_pLB;
    std::map<Int_t, std::pair<Double_t, Double_t> > pLB_timestamps_x, pLB_timestamps_y, pLB_timestamps;
    std::vector<Int_t> plb_list_x, plb_list_y, plb_list;
    std::map<Int_t, Double_t> live_fractions, prescale;
    std::map<Int_t, Double_t> bunch_intensities_1, bunch_intensities_2;
    std::map<Int_t, Double_t> nominal_separation;

    void CalculateMuPlb(); // Convert vertex counts to average mu per pLB, using pLB timestamps and deadtimes.
    std::map<Int_t, Double_t> mu_pLB;
    std::map<Int_t, Double_t> mu_err_pLB;
    std::map<Int_t, Double_t> mu_raw_pLB;
    std::map<Int_t, Double_t> mu_raw_err_pLB;
    std::map<Int_t, Double_t> mu_real_pLB;
    std::map<Int_t, Double_t> mu_real_err_pLB;

    void InitializeFakeCorrection(TString p_energy, TString p_settings, Int_t p_ntrkcut);
    void CorrectPileupEffects(TH2D *h_z_plb, TString p_energy, TString p_settings, Int_t p_ntrkcut, TString p_run, Int_t p_bcid);
#ifdef REMOVED_051612
    void InitializeMaskingCorrection(TString p_energy, TString p_settings, Int_t p_ntrkcut, TH1D *h_z, TString save_tag);
    void InitializeMaskingCorrection(TString p_energy, TString p_settings, Int_t p_ntrkcut, TString save_tag);
    void CorrectPileupEffects();
    PileupMaskingCorrection *pmc;
#endif
    FakeCorrection *fc;

    void CalculateMuSpecPlb(); // Divide by bunch intensities
    std::map<Int_t, Double_t> musp_pLB;
    std::map<Int_t, Double_t> musp_err_pLB;

    // Variety of fit functions
    std::vector<TString> fit_functions;
    TGraphErrors* tg_musp_x;
    TGraphErrors* tg_musp_x_ini;
    TGraphErrors* tg_musp_y;
    TGraphErrors* tg_musp_y_ini;
    std::map<TString, TF1*> f_x, f_y;
    TSpline *fit_spline_x;
    TSpline *fit_spline_y;
    std::map<TString, TFitResult> fr_x_actual, fr_y_actual;
    std::map<TString, TFitResultPtr> fr_x, fr_y;
    std::map<TString, TFitResultPtr> fr_x_ini, fr_y_ini;
    std::map<TString, Int_t> fit_status;

    // Parameters needed for luminosity calculation
    std::map<TString, Double_t> mu_max_x, mu_max_y;
    std::map<TString, Double_t> Sigma_x, Sigma_y;
    std::map<TString, Double_t> c_x, c_y;
    std::map<TString, Double_t> r_x, r_y;

    std::map<TString, Double_t> mu_max_x_err, mu_max_y_err;
    std::map<TString, Double_t> Sigma_x_err, Sigma_y_err;
    std::map<TString, Double_t> c_x_err, c_y_err;


    // Stuff needed for error estimation
    std::map<TString, TMatrixDSym*> fit_cov_x, fit_cov_y;
    TMatrixDSym *fit_cov_x_dg, *fit_cov_y_dg;
    TMatrixDSym *fit_cov_x_dgc, *fit_cov_y_dgc;
    TMatrixDSym *fit_cov_x_dgcl, *fit_cov_y_dgcl;
    std::map<TString, TMatrixD*> v_dsigma_dpx, v_dsigma_dpy;
    std::map<TString, TMatrixD*> v_dsigma_dpx_T, v_dsigma_dpy_T;

    std::map<TString, Double_t> lumi_sp;
    std::map<TString, Double_t> lumi_sp_err;
    std::map<TString, Double_t> lumi;
    std::map<TString, Double_t> lumi_err;
    std::map<TString, Double_t> sigma_vis;
    std::map<TString, Double_t> sigma_vis_err;

    // Pileup corrections
    std::map<Int_t, Double_t> mu_fake_list, mu_fake_uncertainty_list, masking_correction_factors;

    void FitVdmCurves();

    void Finalize(); // Calculate sigma_vis and related.

    // -- Get methods
    Double_t GetChi2Ndf(TString p_fit_name, TString axis);
    Double_t GetFitParameter(TString p_fit_name, TString axis, TString p_fit_parameter);
    Double_t GetFitParameterError(TString p_fit_name, TString axis, TString p_fit_parameter);
    Double_t GetVdmParameter(TString p_fit_name, TString p_variable_name);
    Double_t GetVdmParameterError(TString p_fit_name, TString p_variable_name);
    Double_t GetLumiSp(TString p_fit_name);
    Double_t GetLumiSpErr(TString p_fit_name);
    Double_t GetLumi(TString p_fit_name);
    Double_t GetLumiErr(TString p_fit_name);
    Double_t GetSigmaVis(TString p_fit_name);
    Double_t GetSigmaVisErr(TString p_fit_name);
    TGraphErrors *GetTGraphX();
    TGraphErrors *GetTGraphY();

    void SetSystematicUncertaintyFlag(TString p_systematic_name);
    std::map<TString, bool> systematic_uncertainty_list;

    void DebugPlots(TString p_path, TString p_tag);


  private:
    bool init_vertex_counts, init_deadtime, init_bunch_intensities;
    std::map<TString, bool> init_plb_timestamps;
    bool init_fc, init_pmc;

    Float_t low_displacement, high_displacement;
	Double_t mean_x, mean_y, max_x, max_y, myRMS_x, myRMS_y, sumw_x, sumw_y, sumdev2_x, sumdev2_y;

    Double_t CalculateTSplineIntegral(TSpline *ts, Double_t xmin, Double_t xmax, Int_t npx);
    Double_t CalculateTSplineMax(TSpline *ts, Double_t xmin, Double_t xmax);


};

