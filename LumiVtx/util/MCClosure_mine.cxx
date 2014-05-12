/**
  *  Monte Carlo closure test.
  *  Author: Julia Iturbe 
  *
  ** Input: -i (input MC sample) -t (corrections MC) -o (output location)
  *
  **/

//#define DEBUG
#define NOSPLIT
#define INCLUDE_HACKS

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <stdlib.h>
#include <string>

#include "TROOT.h"
#include "TString.h"
#include "TH1D.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TLeafI.h"
#include "TLeafD.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TPaveText.h"

#include "LumiVtx/LumiVtx.h"
#include "VdM/VtxCalibration.h"
#include "atlasstyle/AtlasUtils.h"
#include "atlasstyle/AtlasLabels.h"

using namespace std;

TGraphErrors* GetPercentDifference(TGraphErrors *tg_in, Double_t value);
void ApplyXshift(TGraphErrors *tg_in, Float_t shift);

// --- Usage
void usage(const char *pname);

int main(int argc, char **argv) {

  SetAtlasStyle();
//  gStyle->SetOptFit(111);
  gStyle->SetOptFit(000);
  gStyle->SetStatY(1);
	gStyle->SetStatX(1);
	gStyle->SetStatW(0.1);
	gStyle->SetStatH(0.1); 
  // --- Scan for command line parameters

  //Input/output detection
  TString input_tag = "";
  TString corrections_tag = "";
  TString output_prefix = "";
  Float_t mu_scale = -1.;

  int c;
  extern int optind;
  extern char* optarg;
  while (1) {
    int option_index = 0;
    static struct option long_options[] = {
      {"input",  required_argument,  0,  'i'},
      {"tag", required_argument, 0, 't'},
      {"outPrefix",  required_argument,  0,  'o'},
      {"help",  no_argument,    0,  'h'},
      {0, 0, 0, 0}
    };
    c=getopt_long(argc, argv, "a:i:t:o:",long_options,&option_index);
    if (c == -1) {
      break;
    }
    switch(c) {
      case 'a':
        input_tag = optarg;
        corrections_tag = optarg;
        output_prefix = optarg;
        break;
      case 'i': {
        input_tag = optarg;
        break;
      }
      case 't': {
        corrections_tag = optarg;
        break;
      }
      case 'o': {
        output_prefix = optarg;
        break;
      }
      case 'h': {
        return 0;
        break;
      }
    }
  } // End input parameter scanning.

  if (input_tag == "" || corrections_tag == "" || output_prefix == "") {
    cout << "Please specify input, corrections, and output stuff" << endl;
    return 0;
  }

  cout << "Loading TTree" << endl;
  
  TFile *f_histograms = new TFile(TString("/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/") + input_tag + TString("/17.2-VtxLumi/InDetTrackD3PD_results_actualXing_both.root"), "READ");
  //TFile *f_histograms = new TFile(TString("/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/") + input_tag + TString("/17.2-VtxLumi/InDetTrackD3PD_results_NGenInt_both.root"), "READ");

  TFile *f_counts = new TFile(TString("/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/") + input_tag + TString("/17.2-VtxLumi/InDetTrackD3PD_results_actualXing_both_tree.root"), "READ");
  //TFile *f_counts = new TFile(TString("/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/") + input_tag + TString("/17.2-VtxLumi/InDetTrackD3PD_results_NGenInt_both_tree.root"), "READ");

  TFile *f_raw = new TFile(TString("/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/") + input_tag + TString("/17.2-VtxLumi/MCClosureTest_mumax20sample_NGenIntVsNVertices.root"), "READ");
  TFile *f_raw_1 = new TFile(TString("/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/") + input_tag + TString("/17.2-VtxLumi/MCClosureTest_mumax75sample_NGenIntVsNVertices.root"), "READ");
  //TFile *f_raw = new TFile(TString("/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/") + input_tag + TString("/17.2-VtxLumi/MCClosureTest_both.root"), "READ");

  TTree *t = (TTree*)f_counts->Get("t_vdm");
  Long64_t n_entries = 100;

  std::vector<Int_t> nTrkCuts;
  //nTrkCuts.push_back(3);
  //nTrkCuts.push_back(4);
  nTrkCuts.push_back(5);

  std::map<Int_t, Int_t> marker_styles;
  marker_styles[3] = 20;
  marker_styles[4] = 21;
  marker_styles[5] = 22;

  std::map<Int_t, Int_t> marker_colors;
  marker_colors[3] = kRed+1;
  marker_colors[4] = kBlue+1;
  marker_colors[5] = kBlack;

  std::map<Int_t, Float_t> x_shifts;
  //x_shifts[5] = -0.05;

  std::vector<TGraphErrors*> graphs_muvis_ngenint, graphs_murec_ngenint;

  std::map<Int_t, TGraphErrors*> tg_nvtxrecon_ngenint; // Total number of vertices reconstructed vs. number of generated interactions
  std::map<Int_t, TGraphErrors*> tg_murec_ngenint, tg_mureal_ngenint, tg_muvis_ngenint, tg_ratio_murec_muvis;
  std::map<Int_t, TGraph*> tg_comparison, tg_murec_ngenint_res, tg_muvis_ngenint_res, tg_ratio_data_fit;

  std::map<Int_t, std::map<Int_t, Float_t> > m_nvtxrecon_ngenint, m_murec_ngenint, m_nvtxrecon_err_ngenint, m_murec_err_ngenint;

  TLegend *l_nolines = new TLegend(0.85, 0.75, 0.95, 0.85);
  l_nolines->SetFillColor(0);
  l_nolines->SetBorderSize(1);

  TCanvas *c_nvtxrecon_ngenint = new TCanvas("c_nvtxrecon_ngenint", "c_nvtxrecon_ngenint", 1200, 800);
  TCanvas *c_murec_ngenint = new TCanvas("c_murec_ngenint", "c_murec_ngenint", 1200, 800);
  TCanvas *c_mureal_ngenint = new TCanvas("c_mureal_ngenint", "c_mureal_ngenint", 1200, 800);
  TCanvas *c_muvis_ngenint = new TCanvas("c_muvis_ngenint", "c_muvis_ngenint", 1200, 800);
  TCanvas *c_ratio_murec_muvis = new TCanvas("c_ratio_murec_muvis", "c_ratio_murec_muvis", 1200, 800);
  TCanvas *c_comparison = new TCanvas("c_comparison", "c_comparison", 1200, 800);
  TCanvas *c_ratio_data_fit = new TCanvas("c_ratio_data_fit", "c_ratio_data_fit", 1200, 800);

  bool draw_first = true;
  bool TwoFits = false;

  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {

    cout << "On NTrk " << *nTrkCut << endl;

    TString hname = "h_NGenInt_NVtx";
    //TString hname = "h_mu_NVtx";
    hname += *nTrkCut;
    TH1D *hist20 = (TH1D*)f_raw->Get(hname);
    TH1D *hist75 = (TH1D*)f_raw_1->Get(hname);
    TH1D *hist = (TH1D*)hist20->Clone();
    hist->Add(hist75);
    n_entries = hist->GetNbinsX();
    n_entries = 55;
    cout << "Number of entries = " << n_entries << endl;
    

    TString tgname = "tg_nvtxrecon_ngenint_NTrk";
    tgname += *nTrkCut;
    tg_nvtxrecon_ngenint[*nTrkCut] = new TGraphErrors(n_entries);
    tg_nvtxrecon_ngenint[*nTrkCut]->SetName(tgname);

    tgname = "tg_murec_ngenint_NTrk";
    tgname += *nTrkCut;
    tg_murec_ngenint[*nTrkCut] = new TGraphErrors(n_entries);
    tg_murec_ngenint[*nTrkCut]->SetName(tgname);

    tgname = "tg_murec_ngenint_res_NTrk";
    tgname += *nTrkCut;
    tg_murec_ngenint_res[*nTrkCut] = new TGraphErrors(n_entries);
    tg_murec_ngenint_res[*nTrkCut]->SetName(tgname);

    tgname = "tg_mureal_ngenint_NTrk";
    tgname += *nTrkCut;
    tg_mureal_ngenint[*nTrkCut] = new TGraphErrors(n_entries);
    tg_mureal_ngenint[*nTrkCut]->SetName(tgname);

    tgname = "tg_muvis_ngenint_NTrk";
    tgname += *nTrkCut;
    tg_muvis_ngenint[*nTrkCut] = new TGraphErrors(n_entries);
    tg_muvis_ngenint[*nTrkCut]->SetName(tgname);

    tgname = "tg_muvis_ngenint_res_NTrk";
    tgname += *nTrkCut;
    tg_muvis_ngenint_res[*nTrkCut] = new TGraphErrors(n_entries);
    tg_muvis_ngenint_res[*nTrkCut]->SetName(tgname);

    tgname = "tg_ratio_murec_muvis_NTrk";
    tgname += *nTrkCut;
    tg_ratio_murec_muvis[*nTrkCut] = new TGraphErrors(n_entries);
    tg_ratio_murec_muvis[*nTrkCut]->SetName(tgname);

    tgname = "tg_comparison_NTrk";
    tgname += *nTrkCut;
    tg_comparison[*nTrkCut] = new TGraph(n_entries);
    tg_comparison[*nTrkCut]->SetName(tgname);

    tgname = "tg_ratio_data_fit_NTrk";
    tgname += *nTrkCut;
    tg_ratio_data_fit[*nTrkCut] = new TGraph(n_entries);
    tg_ratio_data_fit[*nTrkCut]->SetName(tgname);

    TF1 *line0 = new TF1("line0","[0]*x", 0, 22);
    TF1 *line0high = new TF1("line0high","[0]*x", 38, 71);

    Double_t chi2, ndf, slope, slope_err, data, data_err, fiteval, res, b, b_err, ratio;

    TString lineforbox;

    Double_t FIGURE2_RATIO = 0.35;
    Double_t SUBFIGURE_MARGIN = 0.05;

    TPaveText *statsbox = new TPaveText(0.17,0.75,0.4,0.9, "NDC");
    statsbox->SetTextSize(0.03);
    statsbox->SetFillColor(0);
    statsbox->SetTextAlign(12);

    TPaveText *statsbox1 = new TPaveText(0.17,0.75,0.4,0.9, "NDC");
    statsbox1->SetTextSize(0.03);
    statsbox1->SetFillColor(0);
    statsbox1->SetTextAlign(12);

    //Bin 0 is underflow, bin 1 corresponds to poisson_mean = 0 and there is nothing there, bin 2 corresponds to poisson_mean = 1 

    for (int i = 0; i <= 20; i++){
      //Starting at bin 3 because value of poisson mean = 1 is not good for the fits and doesn't provide useful information
      cout << i << ", " << hist->GetBinCenter(i+3) << ", " << hist->GetBinContent(i+3) << endl;
      tg_murec_ngenint[*nTrkCut]->SetPoint(i, hist->GetBinCenter(i+3), hist->GetBinContent(i+3));
      tg_murec_ngenint[*nTrkCut]->SetPointError(i, 0, hist->GetBinError(i+3));
    }

    for (int i = 0; i <= 33; i++){
      //Starting at bin 3 because value of poisson mean = 1 is not good for the fits and doesn't provide useful information
      cout << i << ", " << hist->GetBinCenter(i+39) << ", " << hist->GetBinContent(i+39) << endl;
      tg_murec_ngenint[*nTrkCut]->SetPoint(i+21, hist->GetBinCenter(i+39), hist->GetBinContent(i+39));
      tg_murec_ngenint[*nTrkCut]->SetPointError(i+21, 0, hist->GetBinError(i+39));
    }

    cout << "Number of points in tg_murec_ngenint: " << tg_murec_ngenint[*nTrkCut]->GetN() << endl;

    if (TwoFits == 0){
      cout << "Fitting One Line" << endl;
      tg_murec_ngenint[*nTrkCut]->Fit("line0","SEMR");
      TF1 *fit = tg_murec_ngenint[*nTrkCut]->GetFunction("line0");
      fit->SetLineColor(kBlue);
      fit->SetLineWidth(1);
      Double_t chi2 = fit->GetChisquare();
      Double_t ndf = fit->GetNDF();
      Double_t slope = fit->GetParameter(0);
      Double_t slope_err = fit->GetParError(0);
      Double_t data, data_err, fiteval, res;
      for (int i=0; i < tg_murec_ngenint[*nTrkCut]->GetN(); i++){
        data = tg_murec_ngenint[*nTrkCut]->GetY()[i];
        data_err = tg_murec_ngenint[*nTrkCut]->GetEY()[i];
        fiteval = fit->Eval(tg_murec_ngenint[*nTrkCut]->GetX()[i]);
        if (fit->Eval(0)!=0){cout << "Fit is not passing through zero!!"<<endl;}
        res = (data-fiteval)/data_err;
        if(data_err>0){
          tg_murec_ngenint_res[*nTrkCut]->SetPoint(i, tg_murec_ngenint[*nTrkCut]->GetX()[i], res);        
        }
      }
      std::ostringstream s1, s2;
      s1 << "Fit #chi^{2}/ndf = " << chi2 << "/" << ndf << " = " << chi2/ndf;
      lineforbox = s1.str(); 
      statsbox->AddText(lineforbox);
      s2 << "p0 = " << slope << ", err = " << slope_err;
      lineforbox = s2.str(); 
      statsbox->AddText(lineforbox);
      
      cout << "Number of points in tg_murec_ngenint_res: " << tg_murec_ngenint_res[*nTrkCut]->GetN() << endl;
      cout << endl;
      cout << "///////////// Fit straight line to mu_raw tgraph, chi2/ndf = " << chi2 << "/" << ndf << "=" << chi2/ndf << endl;
      cout << endl;
    } else{
      cout << "Fitting two lines" << endl;
      
      tg_murec_ngenint[*nTrkCut]->Fit("line0","SEMR+");
      TF1 *fit = tg_murec_ngenint[*nTrkCut]->GetFunction("line0");
      fit->SetLineColor(kBlue);
      fit->SetLineWidth(1);
      chi2 = fit->GetChisquare();
      ndf = fit->GetNDF();
      slope = fit->GetParameter(0);
      slope_err = fit->GetParError(0);
      for (Int_t i=0; i < 21; i++){
        data = tg_murec_ngenint[*nTrkCut]->GetY()[i];
        data_err = tg_murec_ngenint[*nTrkCut]->GetEY()[i];
        fiteval = fit->Eval(tg_murec_ngenint[*nTrkCut]->GetX()[i]);
        if (fit->Eval(0)!=0){cout << "Fit is not passing through zero!!"<<endl;}
        res = (data-fiteval)/data_err;
        if (data_err>0){
          tg_murec_ngenint_res[*nTrkCut]->SetPoint(i, tg_murec_ngenint[*nTrkCut]->GetX()[i], res);
        }
      }
      std::ostringstream s1, s2, s3, s4;
      s1 << "Fit 0-22 #chi^{2}/ndf = " << chi2 << "/" << ndf << " = " << chi2/ndf;
      lineforbox = s1.str(); 
      statsbox->AddText(lineforbox);
      s2 << "p0 = " << slope << ", err = " << slope_err;
      lineforbox = s2.str(); 
      statsbox->AddText(lineforbox);
      cout << "///////////// Fit straight line to mu_rec tgraph from 0 to 22, chi2/ndf = " << chi2 << "/" << ndf << "=" << chi2/ndf << endl;

      tg_murec_ngenint[*nTrkCut]->Fit("line0high","SEMR+");
      TF1 *fithigh = tg_murec_ngenint[*nTrkCut]->GetFunction("line0high");
      fithigh->SetLineColor(kGreen);
      fithigh->SetLineWidth(1);
      chi2 = fithigh->GetChisquare();
      ndf = fithigh->GetNDF();
      slope = fithigh->GetParameter(0);
      slope_err = fithigh->GetParError(0);
      b = fithigh->GetParameter(1);
      b_err = fithigh->GetParError(1);
      for (Int_t i=0; i < 34; i++){
        data = tg_murec_ngenint[*nTrkCut]->GetY()[i+21];
        data_err = tg_murec_ngenint[*nTrkCut]->GetEY()[i+21];
        fiteval = fithigh->Eval(tg_murec_ngenint[*nTrkCut]->GetX()[i+21]);
        if (fithigh->Eval(0)!=0){cout << "Fit is not passing through zero!!"<<endl;}
        res = (data-fiteval)/data_err;
        if (data_err>0){
          tg_murec_ngenint_res[*nTrkCut]->SetPoint(i+21, tg_murec_ngenint[*nTrkCut]->GetX()[i+21], res);
        }
      }
      s3 << "Fit 38-71 #chi^{2}/ndf = " << chi2 << "/" << ndf << " = " << chi2/ndf;
      lineforbox = s3.str(); 
      statsbox->AddText(lineforbox);
      s4 << "p0 = " << slope << ", err = " << slope_err;
      lineforbox = s4.str(); 
      statsbox->AddText(lineforbox);

      cout << "///////////// Fit straight line to mu_rec tgraph from 38 to 71, chi2/ndf = " << chi2 << "/" << ndf << "=" << chi2/ndf<< endl;
      cout << "Number of points in tg_murec_ngenint_res: " << tg_murec_ngenint_res[*nTrkCut]->GetN() << endl;
      cout << endl;
    }

    TString legend_entry = "NTrk";
    legend_entry += *nTrkCut;
    l_nolines->AddEntry(tg_murec_ngenint[*nTrkCut], legend_entry, "p");

    TString bname;
    bname = "NVtx";
    bname += *nTrkCut;
    bname += "_BSC";
    for (int entry = 1; entry <= 100; entry++) {
      
      t->GetEntry(entry);

      Int_t c_ngenint = t->GetLeaf("pLB")->GetValue(0);
      Float_t c_nvtxrecon = t->GetLeaf(bname)->GetValue(0);
      Float_t c_ntrig = (Float_t)t->GetLeaf("NTrig")->GetValue(0);

      tg_nvtxrecon_ngenint[*nTrkCut]->SetPoint(entry-1, c_ngenint, c_nvtxrecon);
      m_nvtxrecon_ngenint[*nTrkCut][c_ngenint] = c_nvtxrecon;
      m_nvtxrecon_err_ngenint[*nTrkCut][c_ngenint] = TMath::Sqrt(c_nvtxrecon);

      if (t->GetLeaf("NTrig")->GetValue(0) > 0 && c_nvtxrecon > 0) {
        m_murec_ngenint[*nTrkCut][c_ngenint] = 1. * c_nvtxrecon / c_ntrig;
        m_murec_err_ngenint[*nTrkCut][c_ngenint] = 1. * TMath::Sqrt(c_nvtxrecon) / c_ntrig;
      } 
    }
    cout << "////////////// Comparison plot " << endl;
    for (int k = 2; k <= 22; k++){
      cout << k-2 << ", " << hist->GetBinCenter(k+1) << ", " << m_murec_ngenint[*nTrkCut][k] << ", " << hist->GetBinContent(k+1) << ", " << m_murec_ngenint[*nTrkCut][k]/hist->GetBinContent(k+1) << endl;
      //if (m_murec_ngenint[*nTrkCut][k] > 0){
      tg_comparison[*nTrkCut]->SetPoint(k-2,hist->GetBinCenter(k+1), m_murec_ngenint[*nTrkCut][k]/hist->GetBinContent(k+1));
      //}
    }
    for (int k = 0; k <= 33; k++){
      cout << k << ", " << hist->GetBinCenter(k+39) << ", " << m_murec_ngenint[*nTrkCut][k+38] << ", " << hist->GetBinContent(k+39) << ", " << m_murec_ngenint[*nTrkCut][k+38]/hist->GetBinContent(k+39) << endl;
      //if (m_murec_ngenint[*nTrkCut][k] > 0){
      tg_comparison[*nTrkCut]->SetPoint(k+21,hist->GetBinCenter(k+39), m_murec_ngenint[*nTrkCut][k+38]/hist->GetBinContent(k+39));
      //}
    }
    cout << "Number of points in tg_comparison: " << tg_comparison[*nTrkCut]->GetN() << endl;

    // Apply pileup corrections
    
    cout << "[MCClosure] Creating new FakeCorrection(" << corrections_tag << "," << *nTrkCut << ")" << endl;
    FakeCorrection *fc = new FakeCorrection(corrections_tag, *nTrkCut);

    cout << "[MCClosure] Creating new PileupMaskingCorrection(" << corrections_tag << "," << *nTrkCut << ")" << endl;
		PileupMaskingCorrection *mc = new PileupMaskingCorrection(corrections_tag, *nTrkCut);
		mc->is_MC = kTRUE;
    hname = "hist/PriVtxZpLB_BCID0";
    TH2D *h_z_plb = (TH2D*)f_histograms->Get(hname);
    TH1D *h_z = (TH1D*)h_z_plb->ProjectionX();
    hname = "h_z_NTrk";
    hname += *nTrkCut;
    h_z->SetName(hname);
    TCanvas *canvas = new TCanvas("h_z","h_z",800,600);
    //h_z->GetXaxis()->SetRangeUser(-200,200);
    h_z->Draw();
    canvas->Print("h_z_from_MC_Closure_Test.pdf");

    TString pmc_tag = corrections_tag;
    pmc_tag += *nTrkCut;
    cout << "[MCClosure] Next step: GenerateDzDistribution" << endl;
    mc->GenerateDzDistribution(h_z, pmc_tag);
    cout << "[MCClosure] Next step: GenerateCorrection(mc->GetExpectedDzDistribution())" << endl;
    mc->GenerateCorrection(mc->GetExpectedDzDistribution());

    Int_t count = 0;
    Int_t count_r = 0;
    for (int i=0; i < tg_murec_ngenint[*nTrkCut]->GetN(); i++) {
      Double_t mu_raw = tg_murec_ngenint[*nTrkCut]->GetY()[i];
      Double_t mu_raw_err = tg_murec_ngenint[*nTrkCut]->GetEY()[i];
      Double_t initial_masking_correction_factor = mc->GetCorrectionFactor(mu_raw);
      if (initial_masking_correction_factor < 1.) {
        initial_masking_correction_factor = 1.;
      }
      //cout << initial_masking_correction_factor << endl;
      /* Subtraction method*/
      Double_t mu_fake = fc->GetFakeMuFromMuReconMC(mu_raw * initial_masking_correction_factor);
      Double_t mu_real = mu_raw - mu_fake;
      Double_t final_masking_correction_factor = mc->GetCorrectionFactor(mu_real);
      //cout << final_masking_correction_factor << endl;
      Double_t mu_vis = mu_real * final_masking_correction_factor;
      Double_t mu_vis_err = tg_murec_ngenint[*nTrkCut]->GetEY()[i] * final_masking_correction_factor;

      cout << i << ", mu_raw=" << mu_raw << ", imcf=" << initial_masking_correction_factor << ", mu_fake=" << mu_fake << ", mu_real=" << mu_real << " fmcf=" << final_masking_correction_factor << ", mu_vis=" << mu_vis << endl;

      /* Multiplication method */
      /*
      Double_t fake_fraction = fc->GetFakeFractionFromMuReconMC(mu_raw * initial_masking_correction_factor);
      Double_t mu_vis = mu_raw * (1. - fake_fraction) * initial_masking_correction_factor;
      */

      if (mu_vis<=0){ cout << "For " << i << ", muvis is smaller than 0." << endl;}
      
      if (mu_vis>0){  
        tg_muvis_ngenint[*nTrkCut]->SetPoint(i, tg_murec_ngenint[*nTrkCut]->GetX()[i], mu_vis);
        tg_muvis_ngenint[*nTrkCut]->SetPointError(i, 0., mu_vis_err);
        
        tg_ratio_murec_muvis[*nTrkCut]->SetPoint(i, tg_murec_ngenint[*nTrkCut]->GetX()[i], mu_raw/mu_vis);
        tg_ratio_murec_muvis[*nTrkCut]->SetPointError(i, 0, (mu_raw/mu_vis)*TMath::Sqrt( (mu_raw_err/mu_raw)*(mu_raw_err/mu_raw) + (mu_vis_err/mu_vis)*(mu_vis_err/mu_vis)) );

        //cout << count << " ,Ratio muraw/muvis (%): " << ((mu_raw/mu_vis)-1)*100 << endl;

        count += 1;
      }
      if (mu_real>0){
        tg_mureal_ngenint[*nTrkCut]->SetPoint(count_r, tg_murec_ngenint[*nTrkCut]->GetX()[i], mu_real);
        tg_mureal_ngenint[*nTrkCut]->SetPointError(count_r, 0., 0.);
        count_r+=1;
      }
    }

    cout << "Number of points in tg_muvis_ngenint: " << tg_murec_ngenint[*nTrkCut]->GetN() << endl;
    cout << "Number of points in tg_ratio_murec_muvis: " << tg_ratio_murec_muvis[*nTrkCut]->GetN() << endl;

    if (TwoFits == 0){
      cout << "Fitting One Line" << endl;
      tg_muvis_ngenint[*nTrkCut]->Fit("line0","SEMR");
      TF1 *fit1 = tg_muvis_ngenint[*nTrkCut]->GetFunction("line0");
      fit1->SetLineColor(kBlue);
      fit1->SetLineWidth(1);
      chi2 = fit1->GetChisquare();
      ndf = fit1->GetNDF();
      slope = fit1->GetParameter(0);
      slope_err = fit1->GetParError(0);
      for (Int_t i=0; i < tg_muvis_ngenint[*nTrkCut]->GetN(); i++){
        data = tg_muvis_ngenint[*nTrkCut]->GetY()[i];
        data_err = tg_muvis_ngenint[*nTrkCut]->GetEY()[i];
        fiteval = fit1->Eval(tg_muvis_ngenint[*nTrkCut]->GetX()[i]);
        if (fit1->Eval(0)!=0){cout << "Fit is not passing through zero!!"<<endl;}
        res = (data-fiteval)/data_err;
        ratio = (1 - (data/fiteval))*100;
        if(data_err>0){
          tg_muvis_ngenint_res[*nTrkCut]->SetPoint(i, tg_muvis_ngenint[*nTrkCut]->GetX()[i], res);
          tg_ratio_data_fit[*nTrkCut]->SetPoint(i, tg_muvis_ngenint[*nTrkCut]->GetX()[i], ratio);        
        }
      }
      std::ostringstream s5, s6;
      s5 << "Fit #chi^{2}/ndf = " << chi2 << "/" << ndf << " = " << chi2/ndf;
      lineforbox = s5.str(); 
      statsbox1->AddText(lineforbox);
      s6 << "p0 = " << slope << ", err = " << slope_err;
      lineforbox = s6.str(); 
      statsbox1->AddText(lineforbox);
      cout << endl;
      cout << "///////////// Fit straight line to mu_vis tgraph from 0 to 22, chi2/ndf = " << chi2 << "/" << ndf << "=" << chi2/ndf << endl;
      cout << endl;
    } else{
      cout << "Fitting with two lines" << endl;
      tg_muvis_ngenint[*nTrkCut]->Fit("line0","SEMR+");
      TF1 *fit1 = tg_muvis_ngenint[*nTrkCut]->GetFunction("line0");
      fit1->SetLineColor(kBlue);
      fit1->SetLineWidth(1);
      chi2 = fit1->GetChisquare();
      ndf = fit1->GetNDF();
      slope = fit1->GetParameter(0);
      slope_err = fit1->GetParError(0);
      for (Int_t i=0; i < 21; i++){
        data = tg_muvis_ngenint[*nTrkCut]->GetY()[i];
        data_err = tg_muvis_ngenint[*nTrkCut]->GetEY()[i];
        fiteval = fit1->Eval(tg_muvis_ngenint[*nTrkCut]->GetX()[i]);
        if (fit1->Eval(0)!=0){cout << "Fit is not passing through zero!!"<<endl;}
        res = (data-fiteval)/data_err;
        ratio = (1 - (data/fiteval))*100;
        cout << "For: " << tg_muvis_ngenint[*nTrkCut]->GetX()[i] << ", Data: " << data << ", Fit: " << fiteval << "DataErr: " << data_err << ", Res: " << res << endl;
        if (data_err>0){
          tg_muvis_ngenint_res[*nTrkCut]->SetPoint(i, tg_muvis_ngenint[*nTrkCut]->GetX()[i], res);
          tg_ratio_data_fit[*nTrkCut]->SetPoint(i, tg_muvis_ngenint[*nTrkCut]->GetX()[i], ratio);
        }
      }
      std::ostringstream s5, s6, s7, s8;
      s5 << "Fit 0-22 #chi^{2}/ndf = " << chi2 << "/" << ndf << " = " << chi2/ndf;
      lineforbox = s5.str(); 
      statsbox1->AddText(lineforbox);
      s6 << "p0 = " << slope << ", err = " << slope_err;
      lineforbox = s6.str(); 
      statsbox1->AddText(lineforbox);
      cout << endl;
      cout << "///////////// Fit straight line to mu_vis tgraph from 0 to 22, chi2/ndf = " << chi2 << "/" << ndf << "=" << chi2/ndf << endl;
      cout << endl;

      tg_muvis_ngenint[*nTrkCut]->Fit("line0high","SEMR+");
      TF1 *fithigh1 = tg_muvis_ngenint[*nTrkCut]->GetFunction("line0high");
      fithigh1->SetLineColor(kGreen);
      fithigh1->SetLineWidth(1);
      chi2 = fithigh1->GetChisquare();
      ndf = fithigh1->GetNDF();
      slope = fithigh1->GetParameter(0);
      slope_err = fithigh1->GetParError(0);
      for (Int_t i=0; i < 34; i++){
        data = tg_muvis_ngenint[*nTrkCut]->GetY()[i+21];
        data_err = tg_muvis_ngenint[*nTrkCut]->GetEY()[i+21];
        fiteval = fithigh1->Eval(tg_muvis_ngenint[*nTrkCut]->GetX()[i+21]);
        if (fithigh1->Eval(0)!=0){cout << "Fit is not passing through zero!!"<<endl;}
        res = (data-fiteval)/data_err;
        ratio = (1 - (data/fiteval))*100;
        cout << "For: " << tg_muvis_ngenint[*nTrkCut]->GetX()[i] << ", Data: " << data << ", Fit: " << fiteval << "DataErr: " << data_err << ", Res: " << res << endl;
        //cout << i << " PM: " << tg_muvis_ngenint[*nTrkCut]->GetX()[i+21] << ", Data: " << data << ", Fit: " << fiteval << ", Res: " << res << ", Ratio(%): " << ratio << endl;
        if (data_err>0){
          tg_muvis_ngenint_res[*nTrkCut]->SetPoint(i+21, tg_muvis_ngenint[*nTrkCut]->GetX()[i+21], res);
          tg_ratio_data_fit[*nTrkCut]->SetPoint(i+21, tg_muvis_ngenint[*nTrkCut]->GetX()[i+21], ratio);
        }
      }
      s7 << "Fit 38-71 #chi^{2}/ndf = " << chi2 << "/" << ndf << " = " << chi2/ndf;
      lineforbox = s7.str(); 
      statsbox1->AddText(lineforbox);
      s8 << "p0 = " << slope << ", err = " << slope_err;
      lineforbox = s8.str(); 
      statsbox1->AddText(lineforbox);
      cout << "Number of points in tg_muvis_ngenint_res: " << tg_muvis_ngenint_res[*nTrkCut]->GetN() << endl;
      cout << endl;
      cout << "///////////// Fit straight line to mu_vis tgraph from 38 to 71, chi2/ndf = " << chi2 << "/" << ndf << "=" << chi2/ndf << endl;
      cout << endl;
    }

    TString draw_options;
    if (draw_first) {
      draw_first = false;
      draw_options = "ap";
    } else {
      draw_options = "p";
    }

    // Optionally, draw a label on the plot specifying which sample is used.
    if (draw_options == "ap") {
      myText(0.01, 0.01, kBlack, input_tag, 0.5);
    }

    c_nvtxrecon_ngenint->cd();
    //tg_nvtxrecon_ngenint[*nTrkCut]->SetMinimum(0.2);
    //tg_nvtxrecon_ngenint[*nTrkCut]->SetMaximum(0.9);
    tg_nvtxrecon_ngenint[*nTrkCut]->SetMarkerStyle(marker_styles[*nTrkCut]);
    tg_nvtxrecon_ngenint[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_nvtxrecon_ngenint[*nTrkCut]->GetXaxis()->SetTitle("#mu");
    tg_nvtxrecon_ngenint[*nTrkCut]->GetYaxis()->SetTitle("Total number of vertices reconstructed");
    tg_nvtxrecon_ngenint[*nTrkCut]->Draw(draw_options);

    c_comparison->cd();
    tg_comparison[*nTrkCut]->SetMarkerStyle(marker_styles[*nTrkCut]);
    tg_comparison[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_comparison[*nTrkCut]->GetXaxis()->SetTitle("#mu");
    tg_comparison[*nTrkCut]->GetYaxis()->SetTitle("Ratio of the two definitions");
    tg_comparison[*nTrkCut]->Draw(draw_options);

    c_ratio_murec_muvis->cd();
    tg_ratio_murec_muvis[*nTrkCut]->SetMarkerStyle(marker_styles[*nTrkCut]);
    tg_ratio_murec_muvis[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_ratio_murec_muvis[*nTrkCut]->GetXaxis()->SetTitle("#mu");
    tg_ratio_murec_muvis[*nTrkCut]->GetYaxis()->SetTitle("Ratio of murec/muvis");
    tg_ratio_murec_muvis[*nTrkCut]->Draw(draw_options);

    c_ratio_data_fit->cd();
    tg_ratio_data_fit[*nTrkCut]->SetMarkerStyle(marker_styles[*nTrkCut]);
    tg_ratio_data_fit[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_ratio_data_fit[*nTrkCut]->GetXaxis()->SetTitle("Poisson mean");
    tg_ratio_data_fit[*nTrkCut]->GetYaxis()->SetTitle("1 - #frac{<Nvertices reconstructed>_{corrected}}{Fit} (%)");
    //tg_ratio_data_fit[*nTrkCut]->Print();
    tg_ratio_data_fit[*nTrkCut]->Draw(draw_options);
    TLine *lone = new TLine(tg_ratio_data_fit[*nTrkCut]->GetXaxis()->GetXmin(), 1, tg_ratio_data_fit[*nTrkCut]->GetXaxis()->GetXmax(), 1);
    lone->SetLineColor(kBlue);
    lone->Draw("same");
    TLine *lmone = new TLine(tg_ratio_data_fit[*nTrkCut]->GetXaxis()->GetXmin(), -1, tg_ratio_data_fit[*nTrkCut]->GetXaxis()->GetXmax(), -1);
    lmone->SetLineColor(kBlue);
    lmone->Draw("same");
    TLine *lzero = new TLine(tg_ratio_data_fit[*nTrkCut]->GetXaxis()->GetXmin(), 0, tg_ratio_data_fit[*nTrkCut]->GetXaxis()->GetXmax(), 0);
    lzero->SetLineColor(kRed);
    lzero->Draw("same");

    TLegend *l_muinsteps = new TLegend(0.85, 0.20, 0.95, 0.30);
    l_muinsteps->SetFillColor(0);
    l_muinsteps->SetBorderSize(1);

    c_murec_ngenint->cd();    
    //tg_murec_ngenint[*nTrkCut]->SetMinimum(0.2);
    //tg_murec_ngenint[*nTrkCut]->SetMaximum(0.9);
    tg_murec_ngenint[*nTrkCut]->SetMarkerStyle(24);
    tg_murec_ngenint[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_murec_ngenint[*nTrkCut]->GetXaxis()->SetTitleOffset(5);
    tg_murec_ngenint[*nTrkCut]->GetXaxis()->SetLabelSize(0.0);
    tg_murec_ngenint[*nTrkCut]->GetXaxis()->SetTitleSize(0.0);
    tg_murec_ngenint[*nTrkCut]->GetYaxis()->SetTitle("<Nvertices reconstructed>");
    tg_murec_ngenint[*nTrkCut]->GetYaxis()->SetLabelSize(0.03);
    tg_murec_ngenint[*nTrkCut]->GetYaxis()->SetTitleSize(0.03);
    Double_t xmin = tg_murec_ngenint[*nTrkCut]->GetXaxis()->GetXmin();
    Double_t xmax = tg_murec_ngenint[*nTrkCut]->GetXaxis()->GetXmax();
    tg_murec_ngenint_res[*nTrkCut]->SetMarkerStyle(24);
    tg_murec_ngenint_res[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_murec_ngenint_res[*nTrkCut]->GetXaxis()->SetTitle("Poisson mean");
    tg_murec_ngenint_res[*nTrkCut]->GetXaxis()->SetLabelSize(0.03);
    tg_murec_ngenint_res[*nTrkCut]->GetXaxis()->SetTitleSize(0.03);
    tg_murec_ngenint_res[*nTrkCut]->GetYaxis()->SetTitle("data-fit/#sigma_{data}");
    tg_murec_ngenint_res[*nTrkCut]->GetYaxis()->SetLabelSize(0.02);
    tg_murec_ngenint_res[*nTrkCut]->GetYaxis()->SetTitleSize(0.03);
    tg_murec_ngenint_res[*nTrkCut]->SetMaximum(5);
    tg_murec_ngenint_res[*nTrkCut]->SetMinimum(-5);
    tg_murec_ngenint_res[*nTrkCut]->GetXaxis()->SetRangeUser(xmin,xmax);
    tg_murec_ngenint[*nTrkCut]->Draw(draw_options);
    l_muinsteps->AddEntry(tg_murec_ngenint[*nTrkCut],"#mu_{raw}","P");
    statsbox->Draw("same");
    TLine *line1 = new TLine(xmin,1,xmax,1);
    line1->SetLineColor(kRed);
    c_murec_ngenint->SetBottomMargin(FIGURE2_RATIO);
    TPad *Pad = new TPad( "p_test", "p_test", 0, 0, 1, 1.0-SUBFIGURE_MARGIN, 0, 0, 0); // create new pad, fullsize to have equal font-sizes in both plots                                   
    Pad->SetTopMargin(1-FIGURE2_RATIO);  // top-boundary (should be 1 - thePad->GetBottomMargin() )                                 
    Pad->SetFillStyle(0);  // needs to be transparent                                                                             
    Pad->Draw();
    Pad->cd();
    tg_murec_ngenint_res[*nTrkCut]->Draw(draw_options);
    line1->Draw("same");
    graphs_murec_ngenint.push_back(tg_murec_ngenint[*nTrkCut]);

    c_mureal_ngenint->cd();
    //tg_muvis_ngenint[*nTrkCut]->SetMinimum(0.2);
    //tg_muvis_ngenint[*nTrkCut]->SetMaximum(0.9);
    tg_mureal_ngenint[*nTrkCut]->SetMarkerStyle(25);
    tg_mureal_ngenint[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_mureal_ngenint[*nTrkCut]->GetXaxis()->SetTitle("#mu reconstructed");
    tg_mureal_ngenint[*nTrkCut]->GetYaxis()->SetTitle("<Nvertices reconstructed> after fake correction");
    tg_mureal_ngenint[*nTrkCut]->Draw(draw_options);
    l_muinsteps->AddEntry(tg_mureal_ngenint[*nTrkCut],"#mu_{real}","P");

    c_muvis_ngenint->cd();
    //tg_muvis_ngenint[*nTrkCut]->SetMinimum(0.2);
    //tg_muvis_ngenint[*nTrkCut]->SetMaximum(0.9);
    tg_muvis_ngenint[*nTrkCut]->SetMarkerStyle(26);
    tg_muvis_ngenint[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    TString title = "#mu in steps, nTrkCut";
    title+=*nTrkCut;
    tg_muvis_ngenint[*nTrkCut]->SetTitle(title);
    tg_muvis_ngenint[*nTrkCut]->GetXaxis()->SetTitleOffset(5);
    tg_muvis_ngenint[*nTrkCut]->GetXaxis()->SetLabelSize(0.0);
    tg_muvis_ngenint[*nTrkCut]->GetXaxis()->SetTitleSize(0.0);
    tg_muvis_ngenint[*nTrkCut]->GetYaxis()->SetTitle("<Nvertices reconstructed>  after corrections");
    tg_muvis_ngenint[*nTrkCut]->GetYaxis()->SetLabelSize(0.03);
    tg_muvis_ngenint[*nTrkCut]->GetYaxis()->SetTitleSize(0.03);
    xmin = tg_muvis_ngenint[*nTrkCut]->GetXaxis()->GetXmin();
    xmax = tg_muvis_ngenint[*nTrkCut]->GetXaxis()->GetXmax();
    tg_muvis_ngenint_res[*nTrkCut]->SetMarkerStyle(26);
    tg_muvis_ngenint_res[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_muvis_ngenint_res[*nTrkCut]->GetXaxis()->SetTitle("Poisson mean");
    tg_muvis_ngenint_res[*nTrkCut]->GetXaxis()->SetLabelSize(0.03);
    tg_muvis_ngenint_res[*nTrkCut]->GetXaxis()->SetTitleSize(0.03);
    tg_muvis_ngenint_res[*nTrkCut]->GetYaxis()->SetTitle("data-fit/#sigma_{data}");
    tg_muvis_ngenint_res[*nTrkCut]->GetYaxis()->SetLabelSize(0.02);
    tg_muvis_ngenint_res[*nTrkCut]->GetYaxis()->SetTitleSize(0.03);
    tg_muvis_ngenint_res[*nTrkCut]->SetMaximum(5);
    tg_muvis_ngenint_res[*nTrkCut]->SetMinimum(-5);
    tg_muvis_ngenint_res[*nTrkCut]->GetXaxis()->SetRangeUser(xmin, xmax);
    tg_muvis_ngenint[*nTrkCut]->Draw(draw_options);
    statsbox1->Draw("same");
    TLine *line11 = new TLine(xmin,1,xmax,1);
    line11->SetLineColor(kRed);
    c_muvis_ngenint->SetBottomMargin(FIGURE2_RATIO);
    TPad *Pad1 = new TPad( "p_test1", "p_test1", 0, 0, 1, 1.0-SUBFIGURE_MARGIN, 0, 0, 0); // create new pad, fullsize to have equal font-sizes in both plots                                   
    Pad1->SetTopMargin(1-FIGURE2_RATIO);  // top-boundary (should be 1 - thePad->GetBottomMargin() )                                 
    Pad1->SetFillStyle(0);  // needs to be transparent                                                                             
    Pad1->Draw();
    Pad1->cd();
    tg_muvis_ngenint_res[*nTrkCut]->Draw(draw_options);
    line11->Draw("same");
    graphs_muvis_ngenint.push_back(tg_muvis_ngenint[*nTrkCut]);
    l_muinsteps->AddEntry(tg_muvis_ngenint[*nTrkCut],"#mu_{vis}","P");

    TCanvas *c_muinsteps_ngenint = new TCanvas("c_muinsteps_ngenint", "c_muinsteps_ngenint", 1200, 800);
    TLegend *l_trkcut = new TLegend(0.85, 0.30, 0.95, 0.40);
    l_trkcut->SetFillColor(0);
    l_trkcut->SetBorderSize(1);
    c_muinsteps_ngenint->cd();
    //c_muinsteps_ngenint->SetLogy();
    tg_muvis_ngenint[*nTrkCut]->Draw("ap");
    tg_murec_ngenint[*nTrkCut]->Draw("p");
    tg_mureal_ngenint[*nTrkCut]->Draw("p");
    //TPaveStats *st = (TPaveStats*)tg_muvis_ngenint[*nTrkCut]->GetListOfFunctions()->FindObject("stats");
   	//st->SetX1NDC(0.9); //new x start position
   	//st->SetX2NDC(0.9); //new x end position
    TString legend_justtrakcut = "nTrkCut"; legend_justtrakcut+=*nTrkCut;
    l_trkcut->AddEntry(tg_muvis_ngenint[*nTrkCut],legend_justtrakcut,"P");
    l_muinsteps->Draw("same");
    l_trkcut->Draw("same");
    TString name = TString("") + output_prefix + TString("/Mine/c_muinsteps_ngenint_nTrkCut") ;
    name += *nTrkCut;
    c_muinsteps_ngenint->SaveAs(name+TString(".pdf"));

  }
  cout << "All done." << endl;
  cout << "output_prefix " << output_prefix << endl;

  c_nvtxrecon_ngenint->cd();
  l_nolines->Draw();
  c_nvtxrecon_ngenint->SaveAs(TString("") + output_prefix + TString("/Mine/c_nvtxrecon_poissonmean.pdf"));

  c_comparison->cd();
  TPad *c_Pad = new TPad("grid","",0,0,1,1);
  c_Pad->SetTicks(1,1);
  c_Pad->SetGrid(1,1);
  //l_nolines->Draw();
  c_comparison->SaveAs(TString("") + output_prefix + TString("/Mine/c_comparison.pdf"));

  c_murec_ngenint->cd();
  l_nolines->Draw();
  c_murec_ngenint->SaveAs(TString("") + output_prefix + TString("/Mine/c_muraw_poissonmean.pdf"));

  c_muvis_ngenint->cd();
  l_nolines->Draw();
  c_muvis_ngenint->SaveAs(TString("") + output_prefix + TString("/Mine/c_muvis_poissonmean.pdf"));

  c_mureal_ngenint->cd();
  l_nolines->Draw();
  c_mureal_ngenint->SaveAs(TString("") + output_prefix + TString("/Mine/c_mureal_poissonmean.pdf"));

  c_ratio_murec_muvis->cd();
  l_nolines->Draw();
  c_ratio_murec_muvis->SaveAs(TString("") + output_prefix + TString("/Mine/c_ratio_murec_muvis.pdf"));

  c_ratio_data_fit->cd();
  //l_nolines->Draw();
  c_ratio_data_fit->SaveAs(TString("") + output_prefix + TString("/Mine/c_ratio_data_fit.pdf"));


  TFile *outputFile_muvis = new TFile("muvis_ngenint.root","RECREATE");
  for (vector<TGraphErrors*>::iterator tgraph = graphs_muvis_ngenint.begin(); tgraph != graphs_muvis_ngenint.end(); ++tgraph) {
    (*tgraph)->Write();
  }
  outputFile_muvis->Close();

  TFile *outputFile_murec = new TFile("murec_ngenint.root","RECREATE");
  for (vector<TGraphErrors*>::iterator tgraph = graphs_murec_ngenint.begin(); tgraph != graphs_murec_ngenint.end(); ++tgraph) {
    (*tgraph)->Write();
  }
  outputFile_murec->Close();

}

TGraphErrors* GetPercentDifference(TGraphErrors *tg_in, Double_t value) {

  cout << "[MCClosure] DEBUG : GetPercentDifference(" << tg_in->GetName() << ", " << value << ")" << endl;

  TGraphErrors *tg_out = new TGraphErrors(tg_in->GetN());
  tg_out->SetName(tg_in->GetName() + TString("_difference"));

  for (int i = 0; i < tg_in->GetN(); i++) {

    tg_out->SetPoint(i, tg_in->GetX()[i], (tg_in->GetY()[i] - value) / value * 100.);
    tg_out->SetPointError(i, 0., tg_in->GetEY()[i] / value * 100.);

  }

  return tg_out;

}

void ApplyXshift(TGraphErrors *tg_in, Float_t shift) {

  cout << "[McClosure] INFO : ApplyXshift()..." << std::flush;
  for (int i = 0; i < tg_in->GetN(); i++) {
    tg_in->SetPoint(i, tg_in->GetX()[i] + shift, tg_in->GetY()[i]);
    tg_in->SetPointError(i, tg_in->GetEX()[i], tg_in->GetEY()[i]);
  }
  cout << "done" << endl;

}
