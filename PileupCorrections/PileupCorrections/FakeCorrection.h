/**
  *  Author: David R. Yu (dryu@berkeley.edu)
  *  Class to correct for vertex fakes.
  *
  * Should know about a histogram with (fake fraction) vs. (mu).
  *
*/

#ifndef FakeCorrection_h
#define FakeCorrection_h

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

class FakeCorrection {
  public:
    //Settings
    TString tag; // Tag attached to the name of objects in the cache root file

    //Objects
    //TGraphErrors *tg_fakefraction_NVtxRecon;
    //TGraphErrors *tg_fakefraction_NGenInt;
    //TGraphErrors *tg_fakefraction_mu;
    TGraphErrors *tg_fakefraction_MuReconMC;
    TGraphErrors *tg_fakemu_MuReconMC;
    TGraphErrors *tg_evt_mufake_murecon;
    TGraphErrors *tg_nevt_mu_real_vs_mu_recon;

    //TGraphErrors *tg_fakemu_NGenInt;
    //TGraphErrors *tg_fakemu_mu;
    //TGraphErrors *tg_fakemu_NVtxTightAMC;

    //Methods: three possible inputs,
    // - Mu, from BCM (easiest)
    // - Mu_vis, from vertexing (need to make a map of sorts...)
    // - Mu_vis after masking correction.
    FakeCorrection(TString p_tag, Int_t p_ntrkcut);
    Double_t GetFakeFractionFromMuReconMC(Double_t mu_in);
    Double_t GetFakeMuFromMuReconMC(Double_t mu_in);
    Double_t GetFakeMuUncertaintyFromMuReconMC(Double_t mu_in); // Return statistical uncertainty
    Double_t GetFakeMuHighFromMuReconMC(Double_t mu_in);
    Double_t GetFakeMuLowFromMuReconMC(Double_t mu_in);
    //Double_t GetFakeEventProbabilityFromMuReconMC(Double_t mu_in);
    //Double_t GetMuFakeFromMuReconEvt(Double_t mu_in);
    //Double_t GetMuFakeUncertaintyFromMuReconEvt(Double_t mu_in);

    // -- Event counting
    Double_t GetMuFakeFromMuReconNEvt(Double_t mu_in);
    Double_t GetMuFakeUncertaintyFromMuReconNEvt(Double_t mu_in);


    TGraphErrors* GetFakeCorrectionTGraph();
    void SetMuScale(float p_mu_scale);
    void SetMuScaleUncertainty(float p_mu_scale_uncertainty);

    GlobalSettings gs;
  private:
    TGraph *ConvertTH1ToTGraph(TH1D *h_in);
    float mu_scale;
    float mu_scale_uncertainty;

};


#endif
