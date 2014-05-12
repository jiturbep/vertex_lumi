/**
  *  Author: David R. Yu (dryu@berkeley.edu)
  *  Class to calculate pileup correction lookup TGraphs.
  *  Input: TH1D of pairwise dz distribution
  *  Output: multiple things related to the pileup correction.
  *
  *  How to run:
  *    1. Initialize PileupMaskingCorrection(TH1D* h1, TString string1) with h1 = dz distribution, string1 = prefix (for every member's name). \
  *    2. Alter settings if needed (see the constructor for default settings).
  *    3. Do void PileupMaskingCorrection::RunAllCalculations() to run the whole thing.
  *    4. The relevant members can be accessed/modified directly, for plotting and further analysis, or obtained through Get... functions.
  *
*/

#ifndef PileupMaskingCorrection_h
#define PileupMaskingCorrection_h

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
#include "TH3D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TLeafI.h"
#include "TLine.h"
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
#include "atlasstyle/AtlasStyle.h"
#include "atlasstyle/AtlasUtils.h"
#include "atlasstyle/AtlasLabels.h"
#include "GlobalSettings/GlobalSettings.h"

class PileupMaskingCorrection {
  public:
    PileupMaskingCorrection(TH1D *h1, TString string1="");
#ifdef OLD
    PileupMaskingCorrection(TString f_cache, TString pmask_name, TString string1);
#endif
#ifdef REMOVED_061412
    PileupMaskingCorrection(TString p_source, TString p_energy, TString p_settings, TString p_ntrkcut);
#endif
    PileupMaskingCorrection(TString p_tag, Int_t p_ntrkcut);
    ~PileupMaskingCorrection(); // destructor
    
    //Objects
    TH1D *h_dz;
    TH1D *h_dz_rebinned;
    TH1D *h_dz_random;
    TH1D *h_dz_expected;
    TF1 *f_dz_excluded;
    TF1 *f_dz_full;
    TF1 *f_dz_full_normalized;
    TH1D *h_dz_rebinned_residuals;
    TH1D *h_pmask_dz;
    Double_t pmask;
    TGraph *tg_pmask_cumulative;
    Double_t pmask_err;
    TGraphErrors *tg_mu_obs_vs_mu_actual;
    TGraphErrors *tg_pileup_correction;

    //Settings
    TString tag;
    Int_t rebin_factor;
    Float_t exclude_dz;
    bool low_stats;
    Float_t max_dz;
    bool do_ll;
    bool redo_cache;
    Double_t MaxNGenInt;

    //Methods
    bool pmask_available;
    void GenerateDzDistribution(TH1D *h_z, TString tag = "");
    void LoadDzDistribution(TH1D *h_dz_in);
    void GenerateNewPmask();
    void FitExcluded();
    void FitExcludedResiduals();
    void MakeFullGaussian();
    void MakeDifferentialPmask();

#ifdef REMOVED_061412
    void GenerateCorrection(Float_t sigma_z);
    void CalculateTotalPmask(Float_t sigma_z);
#endif
    void LoadCorrection(TGraphErrors *tg1, TGraphErrors *tg2);
    void GenerateCorrection(TH1D *h_dz_input);
    void CalculateTotalPmask(TH1D *h_dz_input);
    void MakePuCorrTGraphs();

    Double_t GetCorrectionFactor(Double_t mu_in); // Main function to get the correction.
    Double_t GetCorrectionFactorHigh(Double_t mu_in, Double_t mu_in_err);
    Double_t GetCorrectionFactorLow(Double_t mu_in, Double_t mu_in_err);
    Double_t GetGaussianSigma();

    TH1D *GetExpectedDzDistribution();
    TF1 *GetExcludedDzFitFunc();
    TF1 *GetFullDzFitFunc();
    TH1D *GetDifferentialPmask();
    Double_t GetTotalPmask();
    Double_t GetTotalPmaskError();
    TGraphErrors *GetMuMap();
    TGraphErrors *GetMuCorrection();

    void Save(TString path_rootfile = "include/pmask_cache.root", TString suffix = "", bool new_file = true);
    void SaveDebugHistograms(TString path_rootfile, TString p_suffix, bool new_file);
    void SetPmaskScale(float p_pmask_scale);

    /* Fit with gaussian */
    Double_t fitFunc_gaussian(Double_t *x, Double_t *par);
    Double_t fitFunc_gaussian_exclude(Double_t *x, Double_t *par);//

    /* Fit with a dz template */
    Double_t fitFunc_generic(Double_t *x, Double_t *par);
    Double_t fitFunc_generic_exclude(Double_t *x, Double_t *par);
    TH1D *dz_template;//

    Double_t meanNObs(Double_t p_mask, Double_t n_gen);
    Double_t meanMuObs(Double_t mu_actual, Double_t p_mask);
    bool reject;
    bool is_MC;

  private:
    bool finished;
    float pmask_scale;
};

#endif
