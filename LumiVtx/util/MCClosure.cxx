/**
  *  Monte Carlo closure test.
  *  Author: David R. Yu
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
  gStyle->SetOptFit(111);
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
  //TFile *f_histograms = new TFile(TString("/eliza18/atlas/dryu/Luminosity/VertexCounts/MC/") + input_tag + TString("/v9/InDetTrackD3PD_results.root"), "READ");
  //TFile *f_histograms = new TFile(TString("/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/") + input_tag + TString("/17.2-VtxLumi/InDetTrackD3PD_results_actualXing_both.root"), "READ");
  TFile *f_histograms = new TFile(TString("/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/") + input_tag + TString("/17.2-VtxLumi/InDetTrackD3PD_results_NGenInt_both.root"), "READ");
  //TFile *f_counts = new TFile(TString("/eliza18/atlas/dryu/Luminosity/VertexCounts/MC/") + input_tag + TString("/v9/InDetTrackD3PD_results_tree.root"), "READ");
  //TFile *f_counts = new TFile(TString("/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/") + input_tag + TString("/17.2-VtxLumi/InDetTrackD3PD_results_actualXing_both_tree.root"), "READ");
  TFile *f_counts = new TFile(TString("/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/") + input_tag + TString("/17.2-VtxLumi/InDetTrackD3PD_results_NGenInt_both_tree.root"), "READ");
  TTree *t = (TTree*)f_counts->Get("t_vdm");
  Long64_t n_entries = 110;

  std::vector<Int_t> nTrkCuts;
  nTrkCuts.push_back(3);
  nTrkCuts.push_back(4);
  nTrkCuts.push_back(5);
  //nTrkCuts.push_back(6);
  //nTrkCuts.push_back(7);
  //nTrkCuts.push_back(8);
  //nTrkCuts.push_back(10);

  //nTrkCuts.push_back(5);

  std::map<Int_t, Int_t> marker_styles;
  marker_styles[3] = 20;
  marker_styles[4] = 21;
  marker_styles[5] = 22;
  marker_styles[6] = 20;
  marker_styles[7] = 21;
  marker_styles[8] = 22;
  marker_styles[10] = 20;

  std::map<Int_t, Int_t> marker_colors;
  marker_colors[3] = kRed+1;
  marker_colors[4] = kMagenta+1;
  marker_colors[5] = kBlack;
  marker_colors[6] = kBlue+1;
  marker_colors[7] = kCyan+2;
  marker_colors[8] = kOrange+7;
 	marker_colors[10] = kGreen+2;

  std::map<Int_t, Float_t> x_shifts;
  //x_shifts[5] = -0.05;

  std::vector<TGraphErrors*> graphs_muvis_ngenint, graphs_murec_ngenint;

  std::map<Int_t, TGraphErrors*> tg_nvtxrecon_ngenint; // Total number of vertices reconstructed vs. number of generated interactions
  std::map<Int_t, TGraphErrors*> tg_murec_ngenint, tg_mureal_ngenint, tg_muvis_ngenint;
  std::map<Int_t, TGraphErrors*> tg_vtx_efficiency_ngenint;
  std::map<Int_t, TF1*> f_vtx_efficiency_linear;
  std::map<Int_t, TGraphErrors*> tg_diff_vtx_efficiency_ngenint;
  std::map<Int_t, TLine*> tl_average_vtx_efficiency_ngenint;

  std::map<Int_t, std::map<Int_t, Float_t> > m_nvtxrecon_ngenint, m_murec_ngenint, m_nvtxrecon_err_ngenint, m_murec_err_ngenint;
  std::map<Int_t, TGraphErrors*> tg_murec_muinel, tg_mureal_muinel, tg_muvis_muinel;
  std::map<Int_t, TGraphErrors*> tg_vtx_efficiency_muinel;
  std::map<Int_t, TGraphErrors*> tg_diff_vtx_efficiency_muinel;
  std::map<Int_t, TLine*> tl_average_vtx_efficiency_muinel;

  std::map<Int_t, Float_t> vtx_efficiency_ngenint1;
  std::map<Int_t, TLine*> tl_vtx_efficiency_ngenint1_x55, tl_vtx_efficiency_ngenint1_x22;

  TLegend *l = new TLegend(0.55, 0.72, 0.95, 0.95);
  l->SetFillColor(0);
  l->SetBorderSize(1);
  l->SetNColumns(2);

  TCanvas *c_vtx_efficiency_ngenint = new TCanvas("c_vtx_efficiency_ngenint", "c_vtx_efficiency_ngenint", 1000, 1200);
  TPad *p_vtx_efficiency_ngenint_top = new TPad("p_vtx_efficiency_ngenint_top", "p_vtx_efficiency_ngenint_top", 0., 0.333, 1., 1.);
  p_vtx_efficiency_ngenint_top->SetBottomMargin(0.02);
  p_vtx_efficiency_ngenint_top->Draw();
  c_vtx_efficiency_ngenint->cd();
  TPad *p_vtx_efficiency_ngenint_bottom = new TPad("p_vtx_efficiency_ngenint_bottom", "p_vtx_efficiency_ngenint_bottom", 0., 0., 1., 0.333);
  p_vtx_efficiency_ngenint_bottom->SetTopMargin(0.02);
  p_vtx_efficiency_ngenint_bottom->SetBottomMargin(0.2);
  p_vtx_efficiency_ngenint_bottom->Draw();

  TCanvas *c_vtx_efficiency_muinel = new TCanvas("c_vtx_efficiency_muinel", "c_vtx_efficiency_muinel", 1000, 1200);
  TPad *p_vtx_efficiency_muinel_top = new TPad("p_vtx_efficiency_muinel_top", "p_vtx_efficiency_muinel_top", 0., 0.333, 1., 1.);
  p_vtx_efficiency_muinel_top->SetBottomMargin(0.02);
  p_vtx_efficiency_muinel_top->Draw();
  c_vtx_efficiency_muinel->cd();
  TPad *p_vtx_efficiency_muinel_bottom = new TPad("p_vtx_efficiency_muinel_bottom", "p_vtx_efficiency_muinel_bottom", 0., 0., 1., 0.333);
  p_vtx_efficiency_muinel_bottom->SetTopMargin(0.01);
  p_vtx_efficiency_muinel_bottom->SetBottomMargin(0.15);
  p_vtx_efficiency_muinel_bottom->Draw();

  TLegend *l_nolines = new TLegend(0.85, 0.75, 0.95, 0.85);
  l_nolines->SetFillColor(0);
  l_nolines->SetBorderSize(1);

  std::map<Int_t, TGraphErrors*> tg_raw_vtx_efficiency_ngenint, tg_raw_vtx_efficiency_muinel;
  TCanvas *c_raw_vtx_efficiency_ngenint = new TCanvas("c_raw_vtx_efficiency_ngenint", "c_raw_vtx_efficiency_ngenint", 1200, 800);
  TCanvas *c_raw_vtx_efficiency_muinel = new TCanvas("c_raw_vtx_efficiency_muinel", "c_raw_vtx_efficiency_muinel", 1200, 800);
  TCanvas *c_nvtxrecon_ngenint = new TCanvas("c_nvtxrecon_ngenint", "c_nvtxrecon_ngenint", 1200, 800);
  TCanvas *c_murec_ngenint = new TCanvas("c_murec_ngenint", "c_murec_ngenint", 1200, 800);
  TCanvas *c_muvis_ngenint = new TCanvas("c_muvis_ngenint", "c_muvis_ngenint", 1200, 800);
  TCanvas *c_mureal_ngenint = new TCanvas("c_mureal_ngenint", "c_mureal_ngenint", 1200, 800);
  TCanvas *c_murec_muinel = new TCanvas("c_murec_muinel", "c_murec_muinel", 1200, 800);
  TCanvas *c_muvis_muinel = new TCanvas("c_muvis_muinel", "c_muvis_muinel", 1200, 800);
  TCanvas *c_mureal_muinel = new TCanvas("c_mureal_muinel", "c_mureal_muinel", 1200, 800);

  bool draw_first = true;

  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {

    cout << "On NTrk " << *nTrkCut << endl;

    TString tgname = "tg_nvtxrecon_ngenint_NTrk";
    tgname += *nTrkCut;
    tg_nvtxrecon_ngenint[*nTrkCut] = new TGraphErrors(n_entries);
    tg_nvtxrecon_ngenint[*nTrkCut]->SetName(tgname);

    tgname = "tg_murec_ngenint_NTrk";
    tgname += *nTrkCut;
    tg_murec_ngenint[*nTrkCut] = new TGraphErrors(n_entries);
    tg_murec_ngenint[*nTrkCut]->SetName(tgname);

    TString bname;
    bname = "NVtx";
    bname += *nTrkCut;
    bname += "_BSC";
    for (int entry = 0; entry < n_entries; entry++) {
      t->GetEntry(entry);

      Int_t c_ngenint = t->GetLeaf("pLB")->GetValue(0);
      Float_t c_nvtxrecon = t->GetLeaf(bname)->GetValue(0);
      Float_t c_ntrig = (Float_t)t->GetLeaf("NTrig")->GetValue(0);

      tg_nvtxrecon_ngenint[*nTrkCut]->SetPoint(entry, c_ngenint, c_nvtxrecon);
      m_nvtxrecon_ngenint[*nTrkCut][c_ngenint] = c_nvtxrecon;
      m_nvtxrecon_err_ngenint[*nTrkCut][c_ngenint] = TMath::Sqrt(c_nvtxrecon);

      cout << "c_ngenint = " << c_ngenint << ", c_nvtxrecon = " << c_nvtxrecon << ", c_ntrig = " << c_ntrig << endl;

      if (t->GetLeaf("NTrig")->GetValue(0) > 0 && c_nvtxrecon > 0) {
        tg_murec_ngenint[*nTrkCut]->SetPoint(entry, c_ngenint, c_nvtxrecon / c_ntrig);
        tg_murec_ngenint[*nTrkCut]->SetPointError(entry, 0., TMath::Sqrt(c_nvtxrecon) / c_ntrig);

        m_murec_ngenint[*nTrkCut][c_ngenint] = 1. * c_nvtxrecon / c_ntrig;
        m_murec_err_ngenint[*nTrkCut][c_ngenint] = 1. * TMath::Sqrt(c_nvtxrecon) / c_ntrig;

      } //else {
        //tg_murec_ngenint[*nTrkCut]->SetPoint(entry, c_ngenint, 0.);
        //m_murec_ngenint[*nTrkCut][c_ngenint] = 0.;
      //}
      //cout << "pLB " << c_ngenint << " : NVtx = " << c_nvtxrecon << endl;
      //cout << "\tmu_rec = " << m_murec_ngenint[*nTrkCut][c_ngenint] << endl;

      // Get reconstruction efficiency without pileup
      if (c_ngenint == 1) {
        vtx_efficiency_ngenint1[*nTrkCut] = c_nvtxrecon / c_ntrig;
      }

    }

    tgname= "tg_murec_muinel_NTrk";
    tgname += *nTrkCut;
    tg_murec_muinel[*nTrkCut] = new TGraphErrors(n_entries);
    tg_murec_muinel[*nTrkCut]->SetName(tgname);

    tgname= "tg_mureal_muinel_NTrk";
    tgname += *nTrkCut;
    tg_mureal_muinel[*nTrkCut] = new TGraphErrors(n_entries);
    tg_mureal_muinel[*nTrkCut]->SetName(tgname);

    tgname= "tg_muvis_muinel_NTrk";
    tgname += *nTrkCut;
    tg_muvis_muinel[*nTrkCut] = new TGraphErrors(n_entries);
    tg_muvis_muinel[*nTrkCut]->SetName(tgname);

    Int_t point = 0;
    for (Float_t c_muinel = 1; c_muinel <= n_entries; c_muinel++) {

      Float_t average_murec = 0.;
      Float_t average_murec_err = 0.;
      Float_t weight = TMath::Exp(-1. * c_muinel); // Include the zero bin in the weight.

      for (Int_t c_ngenint = 1; c_ngenint < n_entries; c_ngenint++) {
        if (m_murec_ngenint[*nTrkCut][c_ngenint] == 0) {
          continue;
        }

        Float_t poisson_factor = TMath::Exp(-1. * c_muinel) * TMath::Power(c_muinel, c_ngenint) / TMath::Factorial(c_ngenint);
        average_murec += poisson_factor * m_murec_ngenint[*nTrkCut][c_ngenint];
        average_murec_err += TMath::Power(poisson_factor * m_murec_err_ngenint[*nTrkCut][c_ngenint], 2);

        weight += poisson_factor;
      }
      average_murec = average_murec / weight;
      average_murec_err = TMath::Sqrt(average_murec_err) / weight;
      tg_murec_muinel[*nTrkCut]->SetPoint(point, c_muinel, average_murec);
      tg_murec_muinel[*nTrkCut]->SetPointError(point, 0., average_murec_err);

      //if(m_murec_ngenint[*nTrkCut][c_muinel]>0){
      //tg_murec_muinel[*nTrkCut]->SetPoint(point, c_muinel, m_murec_ngenint[*nTrkCut][c_muinel]);
      //tg_murec_muinel[*nTrkCut]->SetPointError(point, 0., m_murec_err_ngenint[*nTrkCut][c_muinel]);}
      //cout << "m_murec_ngenint[*nTrkCut]["<<c_muinel<<"] = " << m_murec_ngenint[*nTrkCut][c_muinel] << endl;

      point++;
    }

    // Apply pileup corrections
    tgname = "tg_mureal_ngenint_NTrk";
    tgname += *nTrkCut;
    tg_mureal_ngenint[*nTrkCut] = new TGraphErrors(n_entries);
    tg_mureal_ngenint[*nTrkCut]->SetName(tgname);

    tgname = "tg_muvis_ngenint_NTrk";
    tgname += *nTrkCut;
    tg_muvis_ngenint[*nTrkCut] = new TGraphErrors(n_entries);
    tg_muvis_ngenint[*nTrkCut]->SetName(tgname);

    //FakeCorrection *fc = new FakeCorrection("mc_8TeV_17.2_default", *nTrkCut);
    cout << "[MCClosure] Creating new FakeCorrection(" << corrections_tag << "," << *nTrkCut << ")" << endl;
    FakeCorrection *fc = new FakeCorrection(corrections_tag, *nTrkCut);
		//TString masking_correction_tag = "data_8TeV_17.2-VtxLumi_207216";
		TString masking_correction_tag = "mc_8TeV_17.2_VtxLumi_2newsets";
    cout << "[MCClosure] Creating new PileupMaskingCorrection(" << masking_correction_tag << "," << *nTrkCut << ")" << endl;
		PileupMaskingCorrection *mc = new PileupMaskingCorrection(masking_correction_tag, *nTrkCut);
    TString hname = "hist/PriVtxZpLB_BCID0";
    TH2D *h_z_plb = (TH2D*)f_histograms->Get(hname);
    TH1D *h_z = (TH1D*)h_z_plb->ProjectionX();
    hname = "h_z_NTrk";
    hname += *nTrkCut;
    h_z->SetName(hname);

    TString pmc_tag = corrections_tag;
    pmc_tag += *nTrkCut;
    cout << "[MCClosure] Next step: GenerateDzDistribution" << endl;
    mc->GenerateDzDistribution(h_z, pmc_tag);
    cout << "[MCClosure] Next step: GenerateCorrection(mc->GetExpectedDzDistribution())" << endl;
    mc->GenerateCorrection(mc->GetExpectedDzDistribution());

    for (int i=0; i < tg_murec_ngenint[*nTrkCut]->GetN(); i++) {
      Double_t mu_raw = tg_murec_ngenint[*nTrkCut]->GetY()[i];
      Double_t initial_masking_correction_factor = mc->GetCorrectionFactor(mu_raw);
      if (initial_masking_correction_factor < 1.) {
        initial_masking_correction_factor = 1.;
      }
      cout << i << endl;
      cout << "For mu_raw = " << mu_raw << ", initial_masking_correction_factor = " << initial_masking_correction_factor << endl;
      /* Subtraction method*/

      Double_t mu_fake = fc->GetFakeMuFromMuReconMC(mu_raw * initial_masking_correction_factor);
      //cout << "For mu_raw*initial_masking_correction_factor = " << mu_raw*initial_masking_correction_factor << ", mu_fake = " << mu_fake << endl;

      Double_t mu_real = mu_raw - mu_fake;
      Double_t final_masking_correction_factor = mc->GetCorrectionFactor(mu_real);
      //cout << "For mu_real = mu_raw - mu_fake = " << mu_real << ", final_masking_correction_factor = " << final_masking_correction_factor << endl;
      Double_t mu_vis = mu_real * final_masking_correction_factor;
      cout << "Finally, mu_vis = mu_real*final_masking_correction_factor = " << mu_real << endl;

      Double_t mu_vis_err = tg_murec_ngenint[*nTrkCut]->GetEY()[i] * final_masking_correction_factor;

      /* Multiplication method */
      /*
      Double_t fake_fraction = fc->GetFakeFractionFromMuReconMC(mu_raw * initial_masking_correction_factor);
      Double_t mu_vis = mu_raw * (1. - fake_fraction) * initial_masking_correction_factor;
      */
      if (mu_vis>0){
      tg_muvis_ngenint[*nTrkCut]->SetPoint(i, tg_murec_ngenint[*nTrkCut]->GetX()[i], mu_vis);
      tg_muvis_ngenint[*nTrkCut]->SetPointError(i, 0., mu_vis_err);}

      if (mu_real>0){
      tg_mureal_ngenint[*nTrkCut]->SetPoint(i, tg_murec_ngenint[*nTrkCut]->GetX()[i], mu_real);
      tg_mureal_ngenint[*nTrkCut]->SetPointError(i, 0., 0.);}

    }

    tgname = "tg_mureal_muinel_NTrk";
    tgname += *nTrkCut;
    tg_mureal_muinel[*nTrkCut] = new TGraphErrors(80);
    tg_mureal_muinel[*nTrkCut]->SetName(tgname);

    tgname = "tg_muvis_muinel_NTrk";
    tgname += *nTrkCut;
    tg_muvis_muinel[*nTrkCut] = new TGraphErrors(80);
    tg_muvis_muinel[*nTrkCut]->SetName(tgname);

    for (int i = 0; i < n_entries; i++) {
      
      Float_t mu_raw = tg_murec_muinel[*nTrkCut]->GetY()[i];
      Double_t initial_masking_correction_factor = mc->GetCorrectionFactor(mu_raw);
      if (initial_masking_correction_factor < 1.) {
        initial_masking_correction_factor = 1.;
      }
      /* Subtraction method*/

      Double_t mu_fake = fc->GetFakeMuFromMuReconMC(mu_raw * initial_masking_correction_factor);

      Double_t mu_real = mu_raw - mu_fake;
      Double_t final_masking_correction_factor = mc->GetCorrectionFactor(mu_real);
      Double_t mu_vis = mu_real * final_masking_correction_factor;

      Double_t mu_vis_err = tg_murec_ngenint[*nTrkCut]->GetEY()[i] * final_masking_correction_factor;

      tg_mureal_muinel[*nTrkCut]->SetPoint(i, tg_murec_muinel[*nTrkCut]->GetX()[i], mu_real);
      tg_mureal_muinel[*nTrkCut]->SetPointError(i, 0., 0);

      tg_muvis_muinel[*nTrkCut]->SetPoint(i, tg_murec_muinel[*nTrkCut]->GetX()[i], mu_vis);
      tg_muvis_muinel[*nTrkCut]->SetPointError(i, 0., mu_vis_err);
    }

    tgname = "tg_vtx_efficiency_ngenint_NTrk";
    tgname += *nTrkCut;
    tg_vtx_efficiency_ngenint[*nTrkCut] = new TGraphErrors(n_entries);
    tg_vtx_efficiency_ngenint[*nTrkCut]->SetName(tgname);

    for (int i=0; i < tg_muvis_ngenint[*nTrkCut]->GetN(); i++) {
      if (tg_muvis_ngenint[*nTrkCut]->GetY()[i] == 0) {
        continue;
      }
      tg_vtx_efficiency_ngenint[*nTrkCut]->SetPoint(i, tg_muvis_ngenint[*nTrkCut]->GetX()[i], tg_muvis_ngenint[*nTrkCut]->GetY()[i] / tg_muvis_ngenint[*nTrkCut]->GetX()[i] * 100.);
      tg_vtx_efficiency_ngenint[*nTrkCut]->SetPointError(i, 0., tg_muvis_ngenint[*nTrkCut]->GetEY()[i] / tg_muvis_ngenint[*nTrkCut]->GetX()[i] * 100.);
    }

    tgname = "tg_vtx_efficiency_muinel_NTrk";
    tgname += *nTrkCut;
    tg_vtx_efficiency_muinel[*nTrkCut] = new TGraphErrors(n_entries);
    tg_vtx_efficiency_muinel[*nTrkCut]->SetName(tgname);

    for (int i=0; i < tg_muvis_muinel[*nTrkCut]->GetN(); i++) {
      tg_vtx_efficiency_muinel[*nTrkCut]->SetPoint(i, tg_muvis_muinel[*nTrkCut]->GetX()[i], tg_muvis_muinel[*nTrkCut]->GetY()[i] / tg_muvis_muinel[*nTrkCut]->GetX()[i]);
      tg_vtx_efficiency_muinel[*nTrkCut]->SetPointError(i, 0., tg_muvis_muinel[*nTrkCut]->GetEY()[i] / tg_muvis_muinel[*nTrkCut]->GetX()[i]);
    }

    // -- Calculate average
    Double_t average = 0.;
    Int_t points = 0;
    for (int i = 0; i < tg_vtx_efficiency_ngenint[*nTrkCut]->GetN(); i++) {
      if (tg_vtx_efficiency_ngenint[*nTrkCut]->GetY()[i] > 0) {
        if (tg_vtx_efficiency_ngenint[*nTrkCut]->GetY()[i] <= 0) {
          continue;
        }
        average += tg_vtx_efficiency_ngenint[*nTrkCut]->GetY()[i];
        points++;
      }
    }
    average = average / points;
    //cout << "Average = " << average << endl;
    tl_average_vtx_efficiency_ngenint[*nTrkCut] = new TLine(0, average, n_entries, average);
    tl_average_vtx_efficiency_ngenint[*nTrkCut]->SetLineColor(marker_colors[*nTrkCut]);
    tl_average_vtx_efficiency_ngenint[*nTrkCut]->SetLineStyle(2);

    TString fname = "f_vtx_efficiency_linear_NTrk";
    fname += *nTrkCut;
    f_vtx_efficiency_linear[*nTrkCut] = new TF1(fname, "[0] + [1] * x", 0., n_entries);
    f_vtx_efficiency_linear[*nTrkCut]->SetParameter(0, average);
    f_vtx_efficiency_linear[*nTrkCut]->SetParameter(1, 0.);
    tg_vtx_efficiency_ngenint[*nTrkCut]->Fit(f_vtx_efficiency_linear[*nTrkCut], "QR0");
    cout << "Fit results for NTrk" << *nTrkCut << ":" << endl;
    cout << "\teff(mu) = [0] + [1] * x" << endl;
    cout << "\t [0] = " << f_vtx_efficiency_linear[*nTrkCut]->GetParameter(0) << " +/- " << f_vtx_efficiency_linear[*nTrkCut]->GetParError(0) << endl;
    cout << "\t [1] = " << f_vtx_efficiency_linear[*nTrkCut]->GetParameter(1) << " +/- " << f_vtx_efficiency_linear[*nTrkCut]->GetParError(1) << endl;
    cout << "\t chi2/ndf = " << f_vtx_efficiency_linear[*nTrkCut]->GetChisquare() << " / " << f_vtx_efficiency_linear[*nTrkCut]->GetNDF() << " = " << f_vtx_efficiency_linear[*nTrkCut]->GetChisquare() / f_vtx_efficiency_linear[*nTrkCut]->GetNDF() << endl;


    tg_diff_vtx_efficiency_ngenint[*nTrkCut] = GetPercentDifference(tg_vtx_efficiency_ngenint[*nTrkCut], average);
    //ApplyXshift(tg_diff_vtx_efficiency_ngenint[*nTrkCut], x_shifts[*nTrkCut]);

    average = 0.;
    points = 0;
    for (int i = 0; i < tg_vtx_efficiency_muinel[*nTrkCut]->GetN(); i++) {
      if (tg_vtx_efficiency_muinel[*nTrkCut]->GetY()[i] > 0) {
        average += tg_vtx_efficiency_muinel[*nTrkCut]->GetY()[i];
        points++;
      }
    }
    average = average / points;
    cout << "Average = " << average << endl;
    tl_average_vtx_efficiency_muinel[*nTrkCut] = new TLine(0, average, 22, average);
    tl_average_vtx_efficiency_muinel[*nTrkCut]->SetLineColor(marker_colors[*nTrkCut]);
    tl_average_vtx_efficiency_muinel[*nTrkCut]->SetLineStyle(2);

    tg_diff_vtx_efficiency_muinel[*nTrkCut] = GetPercentDifference(tg_vtx_efficiency_muinel[*nTrkCut], average);

    TString draw_options;
    if (draw_first) {
      draw_first = false;
      draw_options = "ap";
    } else {
      draw_options = "p";
    }

    // -- Also, make a TLine showing the reconstruction efficiency at NGenInt=1
    tl_vtx_efficiency_ngenint1_x55[*nTrkCut] = new TLine(0, vtx_efficiency_ngenint1[*nTrkCut], 55, vtx_efficiency_ngenint1[*nTrkCut]);
    tl_vtx_efficiency_ngenint1_x55[*nTrkCut]->SetLineColor(marker_colors[*nTrkCut]);
    tl_vtx_efficiency_ngenint1_x55[*nTrkCut]->SetLineStyle(1);
    tl_vtx_efficiency_ngenint1_x22[*nTrkCut] = new TLine(0, vtx_efficiency_ngenint1[*nTrkCut], 22, vtx_efficiency_ngenint1[*nTrkCut]);
    tl_vtx_efficiency_ngenint1_x22[*nTrkCut]->SetLineColor(marker_colors[*nTrkCut]);
    tl_vtx_efficiency_ngenint1_x22[*nTrkCut]->SetLineStyle(1);

    c_vtx_efficiency_ngenint->cd();
    p_vtx_efficiency_ngenint_top->cd();
    tg_vtx_efficiency_ngenint[*nTrkCut]->SetMinimum(0);
    tg_vtx_efficiency_ngenint[*nTrkCut]->SetMaximum(100);
    tg_vtx_efficiency_ngenint[*nTrkCut]->SetMarkerStyle(marker_styles[*nTrkCut]);
    tg_vtx_efficiency_ngenint[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_vtx_efficiency_ngenint[*nTrkCut]->GetXaxis()->SetTitle("NGenInt");
    tg_vtx_efficiency_ngenint[*nTrkCut]->GetXaxis()->SetLabelSize(0.);
    tg_vtx_efficiency_ngenint[*nTrkCut]->GetXaxis()->SetTitleSize(0.);
    tg_vtx_efficiency_ngenint[*nTrkCut]->GetYaxis()->SetTitle("Vertex Efficiency (%)");
    tg_vtx_efficiency_ngenint[*nTrkCut]->GetXaxis()->Set(100, 0., n_entries);
    tg_vtx_efficiency_ngenint[*nTrkCut]->Draw(draw_options);

    f_vtx_efficiency_linear[*nTrkCut]->SetLineColor(marker_colors[*nTrkCut]+1);
    f_vtx_efficiency_linear[*nTrkCut]->SetLineStyle(1);
    f_vtx_efficiency_linear[*nTrkCut]->SetLineWidth(1);
    f_vtx_efficiency_linear[*nTrkCut]->Draw("same");

    //tl_average_vtx_efficiency_ngenint[*nTrkCut]->Draw("same");
    //tl_vtx_efficiency_ngenint1_x55[*nTrkCut]->Draw("same");

    c_vtx_efficiency_ngenint->cd();
    p_vtx_efficiency_ngenint_bottom->cd();
    tg_diff_vtx_efficiency_ngenint[*nTrkCut]->SetMinimum(-6.);
    tg_diff_vtx_efficiency_ngenint[*nTrkCut]->SetMaximum(6.);
    tg_diff_vtx_efficiency_ngenint[*nTrkCut]->SetMarkerStyle(marker_styles[*nTrkCut]);
    tg_diff_vtx_efficiency_ngenint[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_diff_vtx_efficiency_ngenint[*nTrkCut]->SetLineColor(marker_colors[*nTrkCut]);
    tg_diff_vtx_efficiency_ngenint[*nTrkCut]->GetXaxis()->SetTitle("#mu_{inel}");
    tg_diff_vtx_efficiency_ngenint[*nTrkCut]->GetYaxis()->SetTitle("% from average");
    tg_diff_vtx_efficiency_ngenint[*nTrkCut]->GetYaxis()->SetTitleSize(0.08);
    tg_diff_vtx_efficiency_ngenint[*nTrkCut]->GetXaxis()->SetTitleSize(0.08);
    tg_diff_vtx_efficiency_ngenint[*nTrkCut]->GetYaxis()->SetTitleOffset(0.7);
    tg_diff_vtx_efficiency_ngenint[*nTrkCut]->GetXaxis()->SetTitleOffset(0.7);
    tg_diff_vtx_efficiency_ngenint[*nTrkCut]->GetXaxis()->Set(100, 0., n_entries);
    tg_diff_vtx_efficiency_ngenint[*nTrkCut]->Draw(draw_options);

    // Draw -1, 0, 1% lines
    if (draw_options == "ap") {
      TLine *tl_zero = new TLine(0., 0., tg_diff_vtx_efficiency_ngenint[*nTrkCut]->GetXaxis()->GetXmax(), 0.);
      tl_zero->SetLineStyle(1);
      tl_zero->SetLineColor(kGray);
      tl_zero->Draw("same");

      TLine *tl_plus = new TLine(0., 1., tg_diff_vtx_efficiency_ngenint[*nTrkCut]->GetXaxis()->GetXmax(), 1.);
      tl_plus->SetLineStyle(1);
      tl_plus->SetLineColor(kGray+3);
      tl_plus->Draw("same");

      TLine *tl_minus = new TLine(0., -1., tg_diff_vtx_efficiency_ngenint[*nTrkCut]->GetXaxis()->GetXmax(), -1.);
      tl_minus->SetLineStyle(1);
      tl_minus->SetLineColor(kGray+3);
      tl_minus->Draw("same");
    }

    // Optionally, draw a label on the plot specifying which sample is used.
    if (draw_options == "ap") {
      myText(0.01, 0.01, kBlack, input_tag, 0.5);
    }

    // Now vs. artificially generated mu_inel

    c_vtx_efficiency_muinel->cd();
    p_vtx_efficiency_muinel_top->cd();
    tg_vtx_efficiency_muinel[*nTrkCut]->SetMinimum(0.0);
    tg_vtx_efficiency_muinel[*nTrkCut]->SetMaximum(0.5);
    tg_vtx_efficiency_muinel[*nTrkCut]->SetMarkerStyle(marker_styles[*nTrkCut]);
    tg_vtx_efficiency_muinel[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_vtx_efficiency_muinel[*nTrkCut]->GetXaxis()->SetLabelSize(0.);
    tg_vtx_efficiency_muinel[*nTrkCut]->GetXaxis()->SetTitleSize(0.);
    tg_vtx_efficiency_muinel[*nTrkCut]->GetYaxis()->SetTitle("Vertex Efficiency");
    tg_vtx_efficiency_muinel[*nTrkCut]->Draw(draw_options);
    tl_average_vtx_efficiency_muinel[*nTrkCut]->Draw("same");
    tl_vtx_efficiency_ngenint1_x22[*nTrkCut]->Draw("same");

    p_vtx_efficiency_muinel_bottom->cd();
    tg_diff_vtx_efficiency_muinel[*nTrkCut]->SetMinimum(-5.);
    tg_diff_vtx_efficiency_muinel[*nTrkCut]->SetMaximum(5.);
    tg_diff_vtx_efficiency_muinel[*nTrkCut]->SetMarkerStyle(marker_styles[*nTrkCut]);
    tg_diff_vtx_efficiency_muinel[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_diff_vtx_efficiency_muinel[*nTrkCut]->GetXaxis()->SetTitle("#mu_{inel}");
    tg_diff_vtx_efficiency_muinel[*nTrkCut]->GetYaxis()->SetTitle("% from average");
    tg_diff_vtx_efficiency_muinel[*nTrkCut]->Draw(draw_options);

    // Draw -1, 0, 1% lines
    if (draw_options == "ap") {
      TLine *tl_zero = new TLine(0., 0., tg_diff_vtx_efficiency_muinel[*nTrkCut]->GetXaxis()->GetXmax(), 0.);
      tl_zero->SetLineStyle(1);
      tl_zero->SetLineColor(kGray);
      tl_zero->Draw("same");

      TLine *tl_plus = new TLine(0., 1., tg_diff_vtx_efficiency_muinel[*nTrkCut]->GetXaxis()->GetXmax(), 1.);
      tl_plus->SetLineStyle(1);
      tl_plus->SetLineColor(kGray+3);
      tl_plus->Draw("same");

      TLine *tl_minus = new TLine(0., -1., tg_diff_vtx_efficiency_muinel[*nTrkCut]->GetXaxis()->GetXmax(), -1.);
      tl_minus->SetLineStyle(1);
      tl_minus->SetLineColor(kGray+3);
      tl_minus->Draw("same");
    }


    TString legend_entry = "NTrk";
    legend_entry += *nTrkCut;
    l->AddEntry(tg_vtx_efficiency_muinel[*nTrkCut], legend_entry, "p");
    l_nolines->AddEntry(tg_vtx_efficiency_muinel[*nTrkCut], legend_entry, "p");
    legend_entry = "NTrk";
    legend_entry += *nTrkCut;
    legend_entry += " fit";
    l->AddEntry(f_vtx_efficiency_linear[*nTrkCut], legend_entry, "l");
    //legend_entry = "NTrk"; legend_entry += *nTrkCut; legend_entry += " mean";
    //l->AddEntry(tl_average_vtx_efficiency_muinel[*nTrkCut], legend_entry, "l");
    //legend_entry = "NTrk"; legend_entry += *nTrkCut; legend_entry += "w/o pileup";
    //l->AddEntry(tl_vtx_efficiency_ngenint1_x22[*nTrkCut], legend_entry, "l");

    // -- Raw vertex efficiencies, for comparison
    cout << "[McClosure] INFO : Raw vertex efficiencies" << endl;
    c_raw_vtx_efficiency_ngenint->cd();
    tg_raw_vtx_efficiency_ngenint[*nTrkCut] = new TGraphErrors(tg_murec_ngenint[*nTrkCut]->GetN());
    tgname = "tg_raw_vtx_efficiency_ngenint_NTrk";
    tgname += *nTrkCut;
    tg_raw_vtx_efficiency_ngenint[*nTrkCut]->SetName(tgname);

    for (int i = 0; i < tg_murec_ngenint[*nTrkCut]->GetN(); i++) {
      tg_raw_vtx_efficiency_ngenint[*nTrkCut]->SetPoint(i, tg_murec_ngenint[*nTrkCut]->GetX()[i], tg_murec_ngenint[*nTrkCut]->GetY()[i] / tg_murec_ngenint[*nTrkCut]->GetX()[i]);
      tg_raw_vtx_efficiency_ngenint[*nTrkCut]->SetPointError(i, 0., tg_murec_ngenint[*nTrkCut]->GetEY()[i] / tg_murec_ngenint[*nTrkCut]->GetX()[i]);
    }

    tg_raw_vtx_efficiency_ngenint[*nTrkCut]->SetMinimum(0.0);
    tg_raw_vtx_efficiency_ngenint[*nTrkCut]->SetMaximum(0.5);
    tg_raw_vtx_efficiency_ngenint[*nTrkCut]->SetMarkerStyle(marker_styles[*nTrkCut]);
    tg_raw_vtx_efficiency_ngenint[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_raw_vtx_efficiency_ngenint[*nTrkCut]->GetXaxis()->SetTitle("NGenInt");
    tg_raw_vtx_efficiency_ngenint[*nTrkCut]->GetYaxis()->SetTitle("Raw Vertex Efficiency");
    tg_raw_vtx_efficiency_ngenint[*nTrkCut]->Draw(draw_options);

    c_raw_vtx_efficiency_muinel->cd();
    tg_raw_vtx_efficiency_muinel[*nTrkCut] = new TGraphErrors(tg_murec_muinel[*nTrkCut]->GetN());
    tgname = "tg_raw_vtx_efficiency_muinel_NTrk";
    tgname += *nTrkCut;
    tg_raw_vtx_efficiency_muinel[*nTrkCut]->SetName(tgname);

    for (int i = 0; i < tg_murec_muinel[*nTrkCut]->GetN(); i++) {
      tg_raw_vtx_efficiency_muinel[*nTrkCut]->SetPoint(i, tg_murec_muinel[*nTrkCut]->GetX()[i], tg_murec_muinel[*nTrkCut]->GetY()[i] / tg_murec_muinel[*nTrkCut]->GetX()[i]);
      tg_raw_vtx_efficiency_muinel[*nTrkCut]->SetPointError(i, 0., tg_murec_muinel[*nTrkCut]->GetEY()[i] / tg_murec_muinel[*nTrkCut]->GetX()[i]);
    }

    tg_raw_vtx_efficiency_muinel[*nTrkCut]->SetMinimum(0.0);
    tg_raw_vtx_efficiency_muinel[*nTrkCut]->SetMaximum(0.5);
    tg_raw_vtx_efficiency_muinel[*nTrkCut]->SetMarkerStyle(marker_styles[*nTrkCut]);
    tg_raw_vtx_efficiency_muinel[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_raw_vtx_efficiency_muinel[*nTrkCut]->GetXaxis()->SetTitle("#mu_{inel}");
    tg_raw_vtx_efficiency_muinel[*nTrkCut]->GetYaxis()->SetTitle("Raw Vertex Efficiency");
    tg_raw_vtx_efficiency_muinel[*nTrkCut]->Draw(draw_options);

    c_nvtxrecon_ngenint->cd();
    //tg_nvtxrecon_ngenint[*nTrkCut]->SetMinimum(0.2);
    //tg_nvtxrecon_ngenint[*nTrkCut]->SetMaximum(0.9);
    tg_nvtxrecon_ngenint[*nTrkCut]->SetMarkerStyle(marker_styles[*nTrkCut]);
    tg_nvtxrecon_ngenint[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_nvtxrecon_ngenint[*nTrkCut]->GetXaxis()->SetTitle("#mu");
    tg_nvtxrecon_ngenint[*nTrkCut]->GetYaxis()->SetTitle("Total number of vertices reconstructed");
    tg_nvtxrecon_ngenint[*nTrkCut]->Draw(draw_options);

    TLegend *l_muinsteps = new TLegend(0.85, 0.20, 0.95, 0.30);
    l_muinsteps->SetFillColor(0);
    l_muinsteps->SetBorderSize(1);

    TF1 *line0 = new TF1("line0","[0]*x", 0, n_entries);

    c_murec_ngenint->cd();    
    //tg_murec_ngenint[*nTrkCut]->SetMinimum(0.2);
    //tg_murec_ngenint[*nTrkCut]->SetMaximum(0.9);
    tg_murec_ngenint[*nTrkCut]->SetMarkerStyle(24);
    tg_murec_ngenint[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_murec_ngenint[*nTrkCut]->GetXaxis()->SetTitle("NGenInt");
    tg_murec_ngenint[*nTrkCut]->GetYaxis()->SetTitle("Number of vertices reconstructed");
    tg_murec_ngenint[*nTrkCut]->Draw(draw_options);
    l_muinsteps->AddEntry(tg_murec_ngenint[*nTrkCut],"#mu_{raw}","P");
    graphs_murec_ngenint.push_back(tg_murec_ngenint[*nTrkCut]);

    c_mureal_ngenint->cd();
    //tg_muvis_ngenint[*nTrkCut]->SetMinimum(0.2);
    //tg_muvis_ngenint[*nTrkCut]->SetMaximum(0.9);
    tg_mureal_ngenint[*nTrkCut]->SetMarkerStyle(25);
    tg_mureal_ngenint[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_mureal_ngenint[*nTrkCut]->GetXaxis()->SetTitle("NGenInt");
    tg_mureal_ngenint[*nTrkCut]->GetYaxis()->SetTitle("Number of vertices after fake correction");
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
    tg_muvis_ngenint[*nTrkCut]->GetXaxis()->SetTitle("NGenInt");
    tg_muvis_ngenint[*nTrkCut]->GetYaxis()->SetTitle("Number of vertices after corrections");
    tg_muvis_ngenint[*nTrkCut]->Draw(draw_options);
    graphs_muvis_ngenint.push_back(tg_muvis_ngenint[*nTrkCut]);
    cout << "///////////// Performing straight line fit ///////////////" << endl;
    tg_muvis_ngenint[*nTrkCut]->Fit("line0","SEMR");
    TF1 *fit = tg_muvis_ngenint[*nTrkCut]->GetFunction("line0");
    fit->SetLineColor(kBlue);
    fit->SetLineWidth(1);
    Double_t chi2 = fit->GetChisquare();
    Double_t ndf = fit->GetNDF();
    Double_t slope = fit->GetParameter(0);
    Double_t slope_err = fit->GetParError(0);
    cout << "///////////// Fit straight line chi2/ndf = " << chi2 << "/" << ndf << "=" << chi2/ndf << endl;
    l_muinsteps->AddEntry(tg_muvis_ngenint[*nTrkCut],"#mu_{vis}","P");

    TCanvas *c_muinsteps_ngenint = new TCanvas("c_muinsteps_ngenint", "c_muinsteps_ngenint", 1200, 800);
    TLegend *l_trkcut = new TLegend(0.85, 0.30, 0.95, 0.40);
    l_trkcut->SetFillColor(0);
    l_trkcut->SetBorderSize(1);
    c_muinsteps_ngenint->cd();
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
    TString name = TString("") + output_prefix + TString("/c_muinsteps_ngenint_nTrkCut") ;
    name += *nTrkCut;
    c_muinsteps_ngenint->SaveAs(name+TString(".pdf"));

    c_murec_muinel->cd();    
    tg_murec_muinel[*nTrkCut]->SetMarkerStyle(24);
    tg_murec_muinel[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_murec_muinel[*nTrkCut]->GetXaxis()->SetTitle("#mu_{inel}");
    tg_murec_muinel[*nTrkCut]->GetYaxis()->SetTitle("Number of vertices reconstructed");
    tg_murec_muinel[*nTrkCut]->Draw(draw_options);

    c_mureal_muinel->cd();
    tg_mureal_muinel[*nTrkCut]->SetMarkerStyle(25);
    tg_mureal_muinel[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_mureal_muinel[*nTrkCut]->GetXaxis()->SetTitle("#mu_{inel}");
    tg_mureal_muinel[*nTrkCut]->GetYaxis()->SetTitle("Number of vertices after fake correction");
    tg_mureal_muinel[*nTrkCut]->Draw(draw_options);

    c_muvis_muinel->cd();
    tg_muvis_muinel[*nTrkCut]->SetMarkerStyle(26);
    tg_muvis_muinel[*nTrkCut]->SetMarkerColor(marker_colors[*nTrkCut]);
    tg_muvis_muinel[*nTrkCut]->SetTitle(title);
    tg_muvis_muinel[*nTrkCut]->GetXaxis()->SetTitle("#mu_{inel}");
    tg_muvis_muinel[*nTrkCut]->GetYaxis()->SetTitle("Number of vertices after corrections");
    tg_muvis_muinel[*nTrkCut]->Draw(draw_options);
    cout << "///////////// Performing straight line fit ///////////////" << endl;
    tg_muvis_muinel[*nTrkCut]->Fit("line0","SEMR");
    TF1 *fit2 = tg_muvis_muinel[*nTrkCut]->GetFunction("line0");
    fit2->SetLineColor(kBlue);
    fit2->SetLineWidth(1);
    Double_t chi22 = fit2->GetChisquare();
    Double_t ndf2 = fit2->GetNDF();
    Double_t slope2= fit2->GetParameter(0);
    Double_t slope_err2 = fit2->GetParError(0);
    cout << "///////////// Fit straight line chi2/ndf = " << chi22<< "/" << ndf2 << "=" << chi22/ndf2 << endl;

    TCanvas *c_muinsteps_muinel = new TCanvas("c_muinsteps_muinel", "c_muinsteps_muinel", 1200, 800);
    c_muinsteps_muinel->cd();
    tg_muvis_muinel[*nTrkCut]->Draw("ap");
    tg_murec_muinel[*nTrkCut]->Draw("p");
    tg_mureal_muinel[*nTrkCut]->Draw("p");
    l_muinsteps->Draw("same");
    l_trkcut->Draw("same");
    TString name_muinel = TString("") + output_prefix + TString("/c_muinsteps_muinel_nTrkCut") ;
    name_muinel += *nTrkCut;
    c_muinsteps_muinel->SaveAs(name_muinel+TString(".pdf"));

  }
  cout << "All done." << endl;
  cout << "output_prefix " << output_prefix << endl;
  c_vtx_efficiency_ngenint->cd();
  l->Draw();
  c_vtx_efficiency_ngenint->SaveAs(TString("") + output_prefix + TString("/c_vtx_efficiency_ngenint.pdf"));

  c_vtx_efficiency_muinel->cd();
  l->Draw();
  c_vtx_efficiency_muinel->SaveAs(TString("") + output_prefix + TString("/c_vtx_efficiency_muinel.pdf"));

  c_raw_vtx_efficiency_ngenint->cd();
  l_nolines->Draw();
  c_raw_vtx_efficiency_ngenint->SaveAs(TString("") + output_prefix + TString("/c_raw_vtx_efficiency_ngenint.pdf"));

  c_raw_vtx_efficiency_muinel->cd();
  l_nolines->Draw();
  c_raw_vtx_efficiency_muinel->SaveAs(TString("") + output_prefix + TString("/c_raw_vtx_efficiency_muinel.pdf"));

  c_nvtxrecon_ngenint->cd();
  l_nolines->Draw();
  c_nvtxrecon_ngenint->SaveAs(TString("") + output_prefix + TString("/c_nvtxrecon_ngenint.pdf"));

  c_murec_ngenint->cd();
  l_nolines->Draw();
  c_murec_ngenint->SaveAs(TString("") + output_prefix + TString("/c_muraw_ngenint.pdf"));

  c_muvis_ngenint->cd();
  l_nolines->Draw();
  c_muvis_ngenint->SaveAs(TString("") + output_prefix + TString("/c_muvis_ngenint.pdf"));

  c_mureal_ngenint->cd();
  l_nolines->Draw();
  c_mureal_ngenint->SaveAs(TString("") + output_prefix + TString("/c_mureal_ngenint.pdf"));

  c_murec_muinel->cd();
  l_nolines->Draw();
  c_murec_muinel->SaveAs(TString("") + output_prefix + TString("/c_muraw_muinel.pdf"));

  c_muvis_muinel->cd();
  l_nolines->Draw();
  c_muvis_muinel->SaveAs(TString("") + output_prefix + TString("/c_muvis_muinel.pdf"));

  c_mureal_muinel->cd();
  l_nolines->Draw();
  c_mureal_muinel->SaveAs(TString("") + output_prefix + TString("/c_mureal_muinel.pdf"));

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
