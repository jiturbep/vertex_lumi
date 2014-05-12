/**
  *  Author: David R. Yu (dryu@berkeley.edu)
  *  Make fake fraction TGraphs and save to include/fake_cache.root
  *
*/

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

#include "GlobalSettings/GlobalSettings.h"

#include "atlasstyle/AtlasStyle.h"
#include "atlasstyle/AtlasUtils.h"
#include "atlasstyle/AtlasLabels.h"

#include "PileupCorrections/PileupMaskingCorrection.h"

using namespace std;

TCanvas* PlotNTrkGraphs(std::map<Int_t, TGraphErrors*> *graphs, TString legend_prefix, bool save_pdf, TString save_path, TString x_title, TString y_title);
TCanvas* PlotNTrkGraphs(std::map<Int_t, TH1D*> *graphs, TString legend_prefix, bool save_pdf, TString save_path);

int main(int argc, char **argv) {

  SetAtlasStyle();

  TString tag, output_prefix;
  bool deleteOldHistograms = false;
  bool doEventBased = false;
  bool redo_pileup_corrections = false;
  // --- Scan for command line parameters
  int c;
  extern int optind;
  extern char* optarg;
  while (1) {
    int option_index = 0;
    static struct option long_options[] = { {0,0,0,0} };
    c=getopt_long(argc, argv, "t:o:",long_options,&option_index);
    if (c == -1) break;
    switch(c) {
      case 't':
        {
          tag = optarg;
          break;
        }
      case 'o': 
        {
          output_prefix = optarg;
          break;
        }
    }
  }

  TString path_vertex_histograms;

  std::vector<TString> mc_7TeV_samples;
  mc_7TeV_samples.push_back("mc_7TeV_16.X_normal_pythia6_pu");
  mc_7TeV_samples.push_back("mc_7TeV_16.X_hybrid_pythia6_pu");
  mc_7TeV_samples.push_back("mc_7TeV_17.2_normal_pythia6_pu");
  mc_7TeV_samples.push_back("mc_7TeV_17.2_VtxLumi_pythia6_pu");
  mc_7TeV_samples.push_back("mc_7TeV_17.2_normal_pythia8_pu");
  mc_7TeV_samples.push_back("mc_7TeV_17.2_VtxLumi_pythia8_pu");
  mc_7TeV_samples.push_back("mc_7TeV_17.2_normal_pythia6_pu_noslim");
  mc_7TeV_samples.push_back("mc_7TeV_17.2_normal_pythia8DL_pu");
  mc_7TeV_samples.push_back("mc_7TeV_17.2_normal_pythia8DL_pu_noslim");
  mc_7TeV_samples.push_back("mc_7TeV_17.2_normal_pythia8_pu_bs45");
  mc_7TeV_samples.push_back("mc_7TeV_17.2_VtxLumi_pythia8_pu_bs45");
  mc_7TeV_samples.push_back("mc_7TeV_17.2_normal_pythia8_pu_bs55");

  std::vector<TString> mc_8TeV_samples;
  mc_8TeV_samples.push_back("mc_8TeV_17.2_normal_pythia8_pu");
  mc_8TeV_samples.push_back("mc_8TeV_17.2_VtxLumi_pythia8_pu");
  mc_8TeV_samples.push_back("mc_8TeV_17.2_VtxLumi_2newsets");
  mc_8TeV_samples.push_back("mc_8TeV_17.2_VtxLumi_mumax20");
  mc_8TeV_samples.push_back("mc_8TeV_17.2_VtxLumi_mumax75");

  if ((find(mc_7TeV_samples.begin(), mc_7TeV_samples.end(), tag) == mc_7TeV_samples.end()) && (find(mc_8TeV_samples.begin(), mc_8TeV_samples.end(), tag) == mc_8TeV_samples.end())) {
    cerr << "[initializeFakeCorrection] ERROR : Unknown sample specified: " << tag << ". Exiting..." << endl;
    exit(1);
  } else {
    path_vertex_histograms = GlobalSettings::path_D3PDMCResults;
    path_vertex_histograms += "/";
    path_vertex_histograms += tag;
    path_vertex_histograms += "/";
    path_vertex_histograms += GlobalSettings::path_D3PDMCResults_v;
    path_vertex_histograms += "InDetTrackD3PD_results_bothsamples.root";
  }

  cout << "[initializeFakeCorr] path_vertex_histograms = " << path_vertex_histograms << endl;

  if (output_prefix == "") {
    output_prefix = tag;
  }

  TFile *f_in = new TFile(path_vertex_histograms, "READ");
  
  std::vector<Int_t> nTrkCuts;
  //nTrkCuts.push_back(2);
  //nTrkCuts.push_back(3);
  //nTrkCuts.push_back(4);
  nTrkCuts.push_back(5);
  //nTrkCuts.push_back(6);
  //nTrkCuts.push_back(7);
  //nTrkCuts.push_back(8);
  //nTrkCuts.push_back(10);
  
  std::map<Int_t, Color_t> nTrkColors;
  //nTrkColors[2] = kYellow;
  nTrkColors[3] = kRed;
  nTrkColors[4] = kBlue;
  nTrkColors[5] = kBlack;
  nTrkColors[6] = kGreen;
  nTrkColors[7] = kGreen+3;
  nTrkColors[8] = kCyan+2;
  nTrkColors[10] = kBlue+3;
  

  Int_t mu_points = 2000;
  Float_t mu_max = 71.; //maximun value of ei_actualIntPerXing variable in the high mu sample

  std::map<Int_t, TH1D*> h_NTrig_NGenInt, h_z;
  std::map<Int_t, TH1D*> h_fakes_NGenInt, h_recon_NGenInt;
  std::map<Int_t, TH1D*> h_mu_fake_NGenInt, h_mu_recon_NGenInt;
  std::map<Int_t, TGraphErrors*> tg_mu_fake_mu, tg_mu_recon_mu;
  std::map<Int_t, TGraphErrors*> tg_mu_fake_MuReconMC, tg_fake_fraction_MuReconMC;

  std::map<Int_t, PileupMaskingCorrection*> pmc;

  bool draw_first = true;

  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {

    cout << "On NTrk " << *nTrkCut << endl;

    TString hname, tgname;

    cout << "Loading histograms from D3PD_results" << endl;
    hname = "hist/h_fakes_NGenInt_NTrk"; hname += *nTrkCut;
    h_fakes_NGenInt[*nTrkCut] = (TH1D*)f_in->Get(hname);

    hname = "hist/h_all_NGenInt_NTrk"; hname += *nTrkCut;
    h_recon_NGenInt[*nTrkCut] = (TH1D*)f_in->Get(hname);

    hname = "hist/h_NTrig_NGenInt";
    h_NTrig_NGenInt[*nTrkCut] = (TH1D*)f_in->Get(hname);
    h_NTrig_NGenInt[*nTrkCut]->GetXaxis()->SetTitle("NGenInt");
    h_NTrig_NGenInt[*nTrkCut]->GetYaxis()->SetTitle("Simulated Events");

    hname = "h_mu_fake_NGenInt_NTrk"; hname += *nTrkCut;
    h_mu_fake_NGenInt[*nTrkCut] = (TH1D*)h_fakes_NGenInt[*nTrkCut]->Clone();
    h_mu_fake_NGenInt[*nTrkCut]->SetName(hname);
    h_mu_fake_NGenInt[*nTrkCut]->GetXaxis()->SetTitle("NGenInt");
    h_mu_fake_NGenInt[*nTrkCut]->GetYaxis()->SetTitle("#mu_{fake}");
    h_mu_fake_NGenInt[*nTrkCut]->Divide(h_NTrig_NGenInt[*nTrkCut]);
    for (int bin = 1; bin <= h_mu_fake_NGenInt[*nTrkCut]->GetNbinsX(); bin++) {
      if (h_NTrig_NGenInt[*nTrkCut]->GetBinContent(bin) > 0) {
        h_mu_fake_NGenInt[*nTrkCut]->SetBinError(bin, TMath::Sqrt(h_fakes_NGenInt[*nTrkCut]->GetBinContent(bin)) / h_NTrig_NGenInt[*nTrkCut]->GetBinContent(bin));
      } else {
        h_mu_fake_NGenInt[*nTrkCut]->SetBinError(bin, 0.);
      }
    }

    // Set maximum mu value for x-axes
    Float_t mu_limit = 0.;
    for (int bin = 1; bin <= h_NTrig_NGenInt[*nTrkCut]->GetNbinsX(); bin++) {
      if (h_NTrig_NGenInt[*nTrkCut]->GetBinContent(bin) < 50) {
        mu_limit = h_NTrig_NGenInt[*nTrkCut]->GetBinCenter(bin);
        break;
      } else {
        if (tag == "mc_8TeV_17.2_VtxLumi_2newsets") {mu_limit = 80;}
        else{ mu_limit = 22; }
      }
    }
    //mu_limit = 80;
    cout << "[initializeFakeCorr] DEBUG: Restricting mu" << endl;
    cout << "mu_limit = " << mu_limit << ", mu_max = " << mu_max << endl;
    cout << "h_mu_fake_NGenInt[*nTrkCut]->GetNbinsX() = " << h_mu_fake_NGenInt[*nTrkCut]->GetNbinsX() << endl;

    hname = "h_mu_recon_NGenInt_NTrk"; hname += *nTrkCut;
    h_mu_recon_NGenInt[*nTrkCut] = (TH1D*)h_recon_NGenInt[*nTrkCut]->Clone();
    h_mu_recon_NGenInt[*nTrkCut]->SetName(hname);
    h_mu_recon_NGenInt[*nTrkCut]->GetXaxis()->SetTitle("NGenInt");
    h_mu_recon_NGenInt[*nTrkCut]->GetYaxis()->SetTitle("#mu_{rec}");
    h_mu_recon_NGenInt[*nTrkCut]->Divide(h_NTrig_NGenInt[*nTrkCut]);
    for (int bin = 1; bin <= h_mu_recon_NGenInt[*nTrkCut]->GetNbinsX(); bin++) {
      if (h_NTrig_NGenInt[*nTrkCut]->GetBinContent(bin) > 0.) {
        h_mu_recon_NGenInt[*nTrkCut]->SetBinError(bin, TMath::Sqrt(h_recon_NGenInt[*nTrkCut]->GetBinContent(bin)) / h_NTrig_NGenInt[*nTrkCut]->GetBinContent(bin));
      } else {
        h_mu_recon_NGenInt[*nTrkCut]->SetBinError(bin, 0.);  
      }
    }

    cout << "Smearing with poisson distributions" << endl;

    tgname = "tg_mu_fake_mu_NTrk"; tgname += *nTrkCut;
    tg_mu_fake_mu[*nTrkCut] = new TGraphErrors(mu_points);
    tg_mu_fake_mu[*nTrkCut]->SetName(tgname);
    tg_mu_fake_mu[*nTrkCut]->GetHistogram()->GetXaxis()->SetTitle("#mu_{inel}");
    tg_mu_fake_mu[*nTrkCut]->GetHistogram()->GetYaxis()->SetTitle("#mu_{fake}");

    for (int i = 1; i <= mu_points; i++) {
      Float_t current_mu = mu_max * i / mu_points;
      Float_t average_mu_fake = 0;
      Float_t average_mu_fake_err = 0;

      if (current_mu > mu_limit) break;

      for (int j = 1; j <= h_mu_fake_NGenInt[*nTrkCut]->GetNbinsX(); j++) {
      //for (int j = 1; j <= 10; j++) {
        Int_t current_NGenInt = TMath::Nint(h_mu_fake_NGenInt[*nTrkCut]->GetBinCenter(j));
        Float_t poisson_factor = TMath::Exp(-1. * current_mu) * TMath::Power(current_mu, current_NGenInt) / TMath::Factorial(current_NGenInt);

        Float_t current_mu_fake = h_mu_fake_NGenInt[*nTrkCut]->GetBinContent(j);
        Float_t current_mu_fake_err = h_mu_fake_NGenInt[*nTrkCut]->GetBinError(j);

        if (current_mu_fake >= 0) {
          average_mu_fake += poisson_factor * current_mu_fake;
          //cout << "Line 239, average_mu_fake = "<< average_mu_fake << " = " << poisson_factor<<" * "<<current_mu_fake<< endl;
          average_mu_fake_err += TMath::Power(poisson_factor * current_mu_fake_err, 2);
        }
      }
      average_mu_fake_err = TMath::Sqrt(average_mu_fake_err);

      tg_mu_fake_mu[*nTrkCut]->SetPoint(i-1, current_mu, average_mu_fake);
      tg_mu_fake_mu[*nTrkCut]->SetPointError(i-1, 0., average_mu_fake_err);
    }

    tgname = "tg_mu_recon_mu_NTrk"; tgname += *nTrkCut;
    tg_mu_recon_mu[*nTrkCut] = new TGraphErrors(mu_points);
    tg_mu_recon_mu[*nTrkCut]->SetName(tgname);
    tg_mu_recon_mu[*nTrkCut]->GetHistogram()->GetXaxis()->SetTitle("#mu_{inel}");
    tg_mu_recon_mu[*nTrkCut]->GetHistogram()->GetYaxis()->SetTitle("#mu_{rec}");

    for (int i = 1; i <= mu_points; i++) {
      Float_t current_mu = mu_max * i / mu_points;
      Float_t average_mu_recon = 0;
      Float_t average_mu_recon_err = 0;

      if (current_mu > mu_limit) break;
 
      for (int j = 1; j <= h_mu_recon_NGenInt[*nTrkCut]->GetNbinsX(); j++) {
        Int_t current_NGenInt = TMath::Nint(h_mu_recon_NGenInt[*nTrkCut]->GetBinCenter(j));
        Float_t poisson_factor = TMath::Exp(-1. * current_mu) * TMath::Power(current_mu, current_NGenInt) / TMath::Factorial(current_NGenInt);

        Float_t current_mu_recon = h_mu_recon_NGenInt[*nTrkCut]->GetBinContent(j);
        Float_t current_mu_recon_err = h_mu_recon_NGenInt[*nTrkCut]->GetBinError(j);

        if (current_mu_recon >= 0) {
          average_mu_recon += poisson_factor * current_mu_recon;
          average_mu_recon_err += TMath::Power(poisson_factor * current_mu_recon_err, 2);
        }
      }
      average_mu_recon_err = TMath::Sqrt(average_mu_recon_err);

      tg_mu_recon_mu[*nTrkCut]->SetPoint(i-1, current_mu, average_mu_recon);
      tg_mu_recon_mu[*nTrkCut]->SetPointError(i-1, 0., average_mu_recon_err);
    }
    
    // -- Final TGraphs: x-axis is mu_recon after masking correction applied.
    tgname = "tg_mu_fake_MuReconMC_NTrk"; tgname += *nTrkCut;
    tg_mu_fake_MuReconMC[*nTrkCut] = new TGraphErrors(mu_points);
    tg_mu_fake_MuReconMC[*nTrkCut]->SetName(tgname);
    tg_mu_fake_MuReconMC[*nTrkCut]->GetHistogram()->GetXaxis()->SetTitle("#mu_{rec} * f_{mask}(#mu_{rec})");
    tg_mu_fake_MuReconMC[*nTrkCut]->GetHistogram()->GetYaxis()->SetTitle("#mu_{fake}");

    if (*nTrkCut != 2) {
      pmc[*nTrkCut] = new PileupMaskingCorrection(tag, *nTrkCut);
      pmc[*nTrkCut]->is_MC = kTRUE;
    }
    TH2D *h_z_mu = (TH2D*)f_in->Get("hist/h_privtx_z_mu")->Clone();
    TString name = "h_z_NTrk"; name += *nTrkCut;
    h_z[*nTrkCut] = (TH1D*)h_z_mu->ProjectionX(name);

    if (*nTrkCut != 2) {
      pmc[*nTrkCut]->GenerateDzDistribution(h_z[*nTrkCut]);
      pmc[*nTrkCut]->GenerateCorrection(pmc[*nTrkCut]->GetExpectedDzDistribution());
    }

    cout << "[initializeFakeCorrection] INFO: tg_mu_fake_mu[*nTrkCut]->GetXaxis()->GetXmax() " << tg_mu_fake_mu[*nTrkCut]->GetXaxis()->GetXmax() << endl;

    for (int i = 1; i <= mu_points; i++) {

      Float_t current_mu_inel = tg_mu_fake_mu[*nTrkCut]->GetX()[i-1];
      if (current_mu_inel > mu_limit) break;

      Float_t current_mu_fake = tg_mu_fake_mu[*nTrkCut]->GetY()[i-1];
      Float_t current_mu_fake_err = tg_mu_fake_mu[*nTrkCut]->GetEY()[i-1];

      Float_t current_mu_recon = tg_mu_recon_mu[*nTrkCut]->GetY()[i-1];
      Float_t current_mu_recon_err = tg_mu_recon_mu[*nTrkCut]->GetEY()[i-1];

      Float_t masking_correction_factor = 1;
      if (*nTrkCut != 2) {
        masking_correction_factor = pmc[*nTrkCut]->GetCorrectionFactor(current_mu_recon);
      }
      if (masking_correction_factor < 1.) masking_correction_factor = 1.;
      Float_t current_MuReconMC = current_mu_recon * masking_correction_factor;
      Float_t current_MuReconMC_err = current_mu_recon_err * masking_correction_factor;

      cout << i << ", mu_inel=" << current_mu_inel << ", mu_fake=" << current_mu_fake << ", mu_recon=" << current_mu_recon << ", mcf=" << masking_correction_factor << ", MuReconMC=" << current_MuReconMC << endl;

      tg_mu_fake_MuReconMC[*nTrkCut]->SetPoint(i-1, current_MuReconMC, current_mu_fake);
      tg_mu_fake_MuReconMC[*nTrkCut]->SetPointError(i-1, current_MuReconMC_err, current_mu_fake_err);
    }

    tgname = "tg_fake_fraction_MuReconMC_NTrk"; tgname += *nTrkCut;
    tg_fake_fraction_MuReconMC[*nTrkCut] = new TGraphErrors(mu_points);
    tg_fake_fraction_MuReconMC[*nTrkCut]->SetName(tgname);
    tg_fake_fraction_MuReconMC[*nTrkCut]->GetXaxis()->SetTitle("#mu_{rec} * f_{mask}(#mu_{rec})");
    tg_fake_fraction_MuReconMC[*nTrkCut]->GetYaxis()->SetTitle("Fake fraction");

    for (int i = 1; i <= mu_points; i++) {

      Float_t current_mu_inel = tg_mu_fake_mu[*nTrkCut]->GetX()[i-1];
      if (current_mu_inel > mu_limit) break;

      Float_t current_mu_fake = tg_mu_fake_mu[*nTrkCut]->GetY()[i-1];
      Float_t current_mu_fake_err = tg_mu_fake_mu[*nTrkCut]->GetEY()[i-1];

      Float_t current_mu_recon = tg_mu_recon_mu[*nTrkCut]->GetY()[i-1];
      Float_t current_mu_recon_err = tg_mu_recon_mu[*nTrkCut]->GetEY()[i-1];

      Float_t masking_correction_factor = 1;
      if (*nTrkCut != 2) {
        masking_correction_factor = pmc[*nTrkCut]->GetCorrectionFactor(current_mu_recon);
      }
      if (masking_correction_factor < 1.) masking_correction_factor = 1.;
      Float_t current_MuReconMC = current_mu_recon * masking_correction_factor;
      //cout << "Line 351 current_MuReconMC = " << current_MuReconMC << " = " << current_mu_recon<< " * " << masking_correction_factor << endl;
      Float_t current_MuReconMC_err = current_mu_recon_err * masking_correction_factor;

      if (current_mu_recon != 0) {
        tg_fake_fraction_MuReconMC[*nTrkCut]->SetPoint(i-1, current_MuReconMC, current_mu_fake / current_mu_recon);
        tg_fake_fraction_MuReconMC[*nTrkCut]->SetPointError(i-1, current_MuReconMC_err, current_mu_fake / current_mu_recon * TMath::Sqrt(TMath::Power(current_mu_fake_err / current_mu_fake, 2) + TMath::Power(current_mu_recon_err / current_mu_recon, 2)));
      } else {
        tg_fake_fraction_MuReconMC[*nTrkCut]->SetPoint(i-1, current_MuReconMC, 0.);
      }
    }
  }

  /**
    *  Event counting method
    *  - Method:   For each bin of NGenInt, count <average # real vertices> and <average # fake vertices>.
    *         The probability of a fake event is P(#real = 0) x P(#fake > 0)
    */

  cout << "Starting event counting corrections." << endl;

  // -- 1. Fake event probability vs. NGenInt. 
  std::map<Int_t, TH1D*> h_fake_event_probability_NGenInt, h_recon_event_probability_NGenInt;
  std::map<Int_t, TGraphErrors*> tg_fake_event_probability_mu, tg_recon_event_probability_mu, tg_pfake_vs_precon;
  std::map<Int_t, TGraphErrors*> tg_nevt_mu_real_vs_mu_recon, tg_nevt_realfraction_vs_mu_recon;

  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    
    //cout << "p_all_fake(NGenInt) for NTrk" << *nTrkCut << ":" << endl;
    
    TString hname = "h_fake_event_probability_NGenInt_NTrk"; hname += *nTrkCut;
    h_fake_event_probability_NGenInt[*nTrkCut] = new TH1D(hname, hname, h_fakes_NGenInt[*nTrkCut]->GetNbinsX(), h_fakes_NGenInt[*nTrkCut]->GetXaxis()->GetXmin(), h_fakes_NGenInt[*nTrkCut]->GetXaxis()->GetXmax());
    h_fake_event_probability_NGenInt[*nTrkCut]->GetXaxis()->SetTitle("NGenInt");
    h_fake_event_probability_NGenInt[*nTrkCut]->GetYaxis()->SetTitle("Fake evt probability");

    hname = "h_recon_event_probability_NGenInt_NTrk"; hname += *nTrkCut;
    h_recon_event_probability_NGenInt[*nTrkCut] = new TH1D(hname, hname, h_fakes_NGenInt[*nTrkCut]->GetNbinsX(), h_fakes_NGenInt[*nTrkCut]->GetXaxis()->GetXmin(), h_fakes_NGenInt[*nTrkCut]->GetXaxis()->GetXmax());

    for (int i = 1; i <= h_fakes_NGenInt[*nTrkCut]->GetNbinsX(); i++) {

      Int_t c_ngenint = h_fakes_NGenInt[*nTrkCut]->GetBinCenter(i);
      Float_t p_fake;
      Float_t p_real;
      Float_t p_recon;
      Float_t p_fake_err;
      Float_t p_real_err;
      Float_t p_recon_err;
      if (h_NTrig_NGenInt[*nTrkCut]->GetBinContent(i) * c_ngenint > 0) {
        p_fake = h_fakes_NGenInt[*nTrkCut]->GetBinContent(i) / (h_NTrig_NGenInt[*nTrkCut]->GetBinContent(i) * c_ngenint);
        p_real = (h_recon_NGenInt[*nTrkCut]->GetBinContent(i) - h_fakes_NGenInt[*nTrkCut]->GetBinContent(i)) / (h_NTrig_NGenInt[*nTrkCut]->GetBinContent(i) * c_ngenint);
        p_recon = h_recon_NGenInt[*nTrkCut]->GetBinContent(i) / (h_NTrig_NGenInt[*nTrkCut]->GetBinContent(i) * c_ngenint);

        p_fake_err = TMath::Sqrt(p_fake * (1. - p_fake) / (h_NTrig_NGenInt[*nTrkCut]->GetBinContent(i) * c_ngenint));
        p_real_err = TMath::Sqrt(p_real * (1. - p_real) / (h_NTrig_NGenInt[*nTrkCut]->GetBinContent(i) * c_ngenint));
        p_recon_err = TMath::Sqrt(p_recon * (1. - p_recon) / (h_NTrig_NGenInt[*nTrkCut]->GetBinContent(i) * c_ngenint));

      } else {
        p_fake = 0.;
        p_real = 0.;
        p_recon = 0.;

        p_fake_err = 0.;
        p_real_err = 0.;
        p_recon_err = 0.;
      }

      Float_t p_no_real = TMath::Power(1. - p_real, c_ngenint);
      Float_t p_no_fake = TMath::Power(1. - p_fake, c_ngenint);
      Float_t p_no_recon = TMath::Power(1. - p_recon, c_ngenint);

      Float_t p_no_real_err = p_real_err * c_ngenint * TMath::Power(1. - p_real, c_ngenint - 1);
      Float_t p_no_fake_err = p_fake_err * c_ngenint * TMath::Power(1. - p_fake, c_ngenint - 1);
      Float_t p_no_recon_err = p_recon_err * c_ngenint * TMath::Power(1. - p_recon, c_ngenint - 1);

      Float_t p_all_fake = p_no_real * (1. - p_no_fake);
      Float_t p_all_fake_err = TMath::Sqrt(TMath::Power(p_no_real_err * (1. - p_no_fake), 2) + TMath::Power(p_no_real * p_no_fake_err, 2));

      h_fake_event_probability_NGenInt[*nTrkCut]->SetBinContent(i, p_all_fake);
      h_fake_event_probability_NGenInt[*nTrkCut]->SetBinError(i, p_all_fake_err);

      h_recon_event_probability_NGenInt[*nTrkCut]->SetBinContent(i, 1. - p_no_recon);
      h_recon_event_probability_NGenInt[*nTrkCut]->SetBinError(i, p_no_recon_err);

      //cout << "p_all_fake(" << c_ngenint << ") = " << p_all_fake << endl;

    }

    // -- 2. Poisson distributions

    TString tgname = "tg_fake_event_probability_mu_NTrk"; tgname += *nTrkCut;
    tg_fake_event_probability_mu[*nTrkCut] = new TGraphErrors(mu_points);
    tg_fake_event_probability_mu[*nTrkCut]->SetName(tgname);

    tgname = "tg_recon_event_probability_mu_NTrk"; tgname += *nTrkCut;
    tg_recon_event_probability_mu[*nTrkCut] = new TGraphErrors(mu_points);
    tg_recon_event_probability_mu[*nTrkCut]->SetName(tgname);

    tgname = "tg_pfake_vs_precon_NTrk"; tgname += *nTrkCut;
    tg_pfake_vs_precon[*nTrkCut] = new TGraphErrors(mu_points);
    tg_pfake_vs_precon[*nTrkCut]->SetName(tgname);

    //cout << "p_fake(p_recon) for NTrk" << *nTrkCut << ":" << endl;

    for (int i = 0; i < mu_points; i++) {

      Float_t c_mu = mu_max * i / mu_points;

      Float_t weight = 0.;
      Float_t p_fake_event = 0.;
      Float_t p_recon_event = 0.;
      Float_t p_fake_event_err = 0.;
      Float_t p_recon_event_err = 0.;

      for (int j = 1; j <= h_fake_event_probability_NGenInt[*nTrkCut]->GetNbinsX(); j++) {

        Int_t c_ngenint = TMath::Nint(h_fake_event_probability_NGenInt[*nTrkCut]->GetBinCenter(j));
        Float_t poisson_factor = TMath::Exp(-1. * c_mu) * TMath::Power(c_mu, c_ngenint) / TMath::Factorial(c_ngenint);
        weight += poisson_factor;

        p_fake_event += poisson_factor * h_fake_event_probability_NGenInt[*nTrkCut]->GetBinContent(j);
        p_recon_event += poisson_factor * h_recon_event_probability_NGenInt[*nTrkCut]->GetBinContent(j);

        p_fake_event_err += TMath::Power(poisson_factor * h_fake_event_probability_NGenInt[*nTrkCut]->GetBinError(j), 2);
        p_recon_event_err += TMath::Power(poisson_factor * h_recon_event_probability_NGenInt[*nTrkCut]->GetBinError(j), 2);

      }
      if (i % 300 == 0) {
        //cout << "[debug] weight = " << weight << endl;
        //cout << "[debug] p_fake_event = " << p_fake_event << endl;
      }
      p_fake_event = p_fake_event / weight;
      p_recon_event = p_recon_event / weight;
      p_fake_event_err = TMath::Sqrt(p_fake_event_err) / weight;
      p_recon_event_err = TMath::Sqrt(p_recon_event_err) / weight;

      tg_fake_event_probability_mu[*nTrkCut]->SetPoint(i, c_mu, p_fake_event);
      tg_fake_event_probability_mu[*nTrkCut]->SetPointError(i, 0., p_fake_event_err);

      tg_recon_event_probability_mu[*nTrkCut]->SetPoint(i, c_mu, p_recon_event);
      tg_recon_event_probability_mu[*nTrkCut]->SetPointError(i, 0., p_recon_event_err);

      tg_pfake_vs_precon[*nTrkCut]->SetPoint(i, p_recon_event, p_fake_event);
      tg_pfake_vs_precon[*nTrkCut]->SetPointError(i, p_recon_event_err, p_fake_event_err);

      //if (i % 300 == 0) cout << "\tp_fake(" << p_recon_event << ") = " << p_fake_event << endl;
    }

    // -- Convert probabilities to mus
    tgname = "tg_nevt_mu_real_vs_mu_recon_NTrk"; tgname += *nTrkCut;
    tg_nevt_mu_real_vs_mu_recon[*nTrkCut] = new TGraphErrors(mu_points);
    tg_nevt_mu_real_vs_mu_recon[*nTrkCut]->SetName(tgname);

    tgname = "tg_nevt_realfraction_vs_mu_recon_NTrk"; tgname += *nTrkCut;
    tg_nevt_realfraction_vs_mu_recon[*nTrkCut] = new TGraphErrors(mu_points);
    tg_nevt_realfraction_vs_mu_recon[*nTrkCut]->SetName(tgname);

    for (int i = 0; i < mu_points; i++) {
      Float_t c_precon = tg_pfake_vs_precon[*nTrkCut]->GetX()[i];
      Float_t c_pfake = tg_pfake_vs_precon[*nTrkCut]->GetY()[i];
      Float_t c_precon_err = tg_pfake_vs_precon[*nTrkCut]->GetEX()[i];
      Float_t c_pfake_err = tg_pfake_vs_precon[*nTrkCut]->GetEY()[i];

      Float_t c_mu_recon = -1. * TMath::Log(1. - c_precon);
      Float_t c_mu_real = -1. * TMath::Log(1. - (c_precon - c_pfake));
      Float_t c_mu_recon_err = 1. / (1. - c_precon) * c_precon_err;
      Float_t c_mu_real_err = TMath::Sqrt(TMath::Power(1. / (1. - (c_precon - c_pfake)) * c_precon_err, 2) + TMath::Power(1. / (1. - (c_precon - c_pfake)) * c_pfake_err, 2));
    
      tg_nevt_mu_real_vs_mu_recon[*nTrkCut]->SetPoint(i, c_mu_recon, c_mu_real);
      tg_nevt_mu_real_vs_mu_recon[*nTrkCut]->SetPointError(i, c_mu_recon_err, c_mu_real_err);

      tg_nevt_realfraction_vs_mu_recon[*nTrkCut]->SetPoint(i, c_mu_recon, c_mu_real / c_mu_recon);
      tg_nevt_realfraction_vs_mu_recon[*nTrkCut]->SetPointError(i, c_mu_recon_err, c_mu_real / c_mu_recon * TMath::Sqrt(TMath::Power(c_mu_real_err / c_mu_real, 2) + TMath::Power(c_mu_recon_err / c_mu_recon, 2)));

    }
  }

#ifdef REMOVED_061812
  /* Legacy stuff for Simone */
  std::map<Int_t, TH1D*> h_fake_events_NGenInt;
  std::map<Int_t, TH1D*> h_fake_events_NVtxRecon;
  std::map<Int_t, TH1D*> h_fake_event_probability_NGenInt;
  std::map<Int_t, TH1D*> h_fake_event_probability_NVtxRecon;
  std::map<Int_t, TH1D*> h_NTrig_NVtxRecon;

  std::map<Int_t, TH1D*> h_event_contains_fake_NGenInt, h_event_contains_fake_NVtxRecon;
  std::map<Int_t, TH1D*> h_event_contains_any_NGenInt, h_event_contains_any_NVtxRecon;
  std::map<Int_t, TH1D*> h_p_evt_recon_NGenInt, h_p_evt_fake_NGenInt;
  std::map<Int_t, TGraphErrors*> tg_evt_mufake_murecon, tg_evt_mufake_mu, tg_evt_murecon_mu;

  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    cout << "On NTrk " << *nTrkCut << endl;

    TString hname, tgname;

    cout << "Loading histograms from D3PD_results" << endl;
    hname = "hist/h_NTrig_NVtxRecon_NTrk"; hname += *nTrkCut;
    h_NTrig_NVtxRecon[*nTrkCut] = (TH1D*)f_in->Get(hname);

    hname = "hist/h_fake_events_NGenInt_NTrk"; hname += *nTrkCut;
    h_fake_events_NGenInt[*nTrkCut] = (TH1D*)f_in->Get(hname);

    hname = "hist/h_fake_events_NVtxRecon_NTrk"; hname += *nTrkCut;
    h_fake_events_NVtxRecon[*nTrkCut] = (TH1D*)f_in->Get(hname);

    hname = "h_fake_event_probability_NGenInt_NTrk"; hname += *nTrkCut;
    h_fake_event_probability_NGenInt[*nTrkCut] = (TH1D*)h_fake_events_NGenInt[*nTrkCut]->Clone();
    h_fake_event_probability_NGenInt[*nTrkCut]->SetName(hname);
    h_fake_event_probability_NGenInt[*nTrkCut]->Divide(h_NTrig_NGenInt[*nTrkCut]);

    hname = "h_fake_event_probability_NVtxRecon_NTrk"; hname += *nTrkCut;
    h_fake_event_probability_NVtxRecon[*nTrkCut] = (TH1D*)h_fake_events_NVtxRecon[*nTrkCut]->Clone();
    h_fake_event_probability_NVtxRecon[*nTrkCut]->SetName(hname);
    h_fake_event_probability_NVtxRecon[*nTrkCut]->Divide(h_NTrig_NVtxRecon[*nTrkCut]);

    hname = "hist/h_event_contains_any_NGenInt_NTrk"; hname += *nTrkCut;
    h_event_contains_any_NGenInt[*nTrkCut] = (TH1D*)f_in->Get(hname);

    hname = "hist/h_event_contains_fake_NGenInt_NTrk"; hname += *nTrkCut;
    h_event_contains_fake_NGenInt[*nTrkCut] = (TH1D*)f_in->Get(hname);

    cout << "Calculating event probabilities" << endl;

    hname = "h_p_evt_recon_NGenInt_NTrk"; hname += *nTrkCut;
    h_p_evt_recon_NGenInt[*nTrkCut] = (TH1D*)h_event_contains_any_NGenInt[*nTrkCut]->Clone();
    h_p_evt_recon_NGenInt[*nTrkCut]->SetName(hname);
    h_p_evt_recon_NGenInt[*nTrkCut]->Divide(h_NTrig_NGenInt[*nTrkCut]);
    for (int bin = 1; bin <= h_p_evt_recon_NGenInt[*nTrkCut]->GetNbinsX(); bin++) {
      if (h_NTrig_NGenInt[*nTrkCut]->GetBinContent(bin) > 0) {
        Float_t p = h_p_evt_recon_NGenInt[*nTrkCut]->GetBinContent(bin) / h_NTrig_NGenInt[*nTrkCut]->GetBinContent(bin);
        h_p_evt_recon_NGenInt[*nTrkCut]->SetBinError(bin, TMath::Sqrt(p * (1. - p) / h_NTrig_NGenInt[*nTrkCut]->GetBinContent(bin)));
      } else {
        h_p_evt_recon_NGenInt[*nTrkCut]->SetBinError(bin, 0.);
      }
    }

    hname = "h_p_evt_fake_NGenInt_NTrk"; hname += *nTrkCut;
    h_p_evt_fake_NGenInt[*nTrkCut] = (TH1D*)h_event_contains_fake_NGenInt[*nTrkCut]->Clone();
    h_p_evt_fake_NGenInt[*nTrkCut]->SetName(hname);
    h_p_evt_fake_NGenInt[*nTrkCut]->Divide(h_NTrig_NGenInt[*nTrkCut]);
    for (int bin = 1; bin <= h_p_evt_fake_NGenInt[*nTrkCut]->GetNbinsX(); bin++) {
      if (h_NTrig_NGenInt[*nTrkCut]->GetBinContent(bin) > 0.) {
        Float_t p = h_p_evt_fake_NGenInt[*nTrkCut]->GetBinContent(bin) / h_NTrig_NGenInt[*nTrkCut]->GetBinContent(bin);
        h_p_evt_fake_NGenInt[*nTrkCut]->SetBinError(bin, TMath::Sqrt(p * (1. - p) / h_NTrig_NGenInt[*nTrkCut]->GetBinContent(bin)));
      } else {
        h_p_evt_fake_NGenInt[*nTrkCut]->SetBinError(bin, 0.);
      }
    }

    cout << "Smearing with a Poisson distribution" << endl;

    tgname = "tg_evt_mufake_murecon_NTrk"; tgname += *nTrkCut;
    tg_evt_mufake_murecon[*nTrkCut] = new TGraphErrors(mu_points);
    tg_evt_mufake_murecon[*nTrkCut]->SetName(tgname);

    tgname = "tg_evt_mufake_mu_NTrk"; tgname += *nTrkCut;
    tg_evt_mufake_mu[*nTrkCut] = new TGraphErrors(mu_points);
    tg_evt_mufake_mu[*nTrkCut]->SetName(tgname);

    tgname = "tg_evt_murecon_mu_NTrk"; tgname += *nTrkCut;
    tg_evt_murecon_mu[*nTrkCut] = new TGraphErrors(mu_points);
    tg_evt_murecon_mu[*nTrkCut]->SetName(tgname);


    for (int i = 1; i <= mu_points; i++) {
      Float_t current_mu = mu_max * i / mu_points;

      Float_t average_p_fake = 0.;
      Float_t average_p_recon = 0.;
      Float_t average_p_fake_err = 0.;
      Float_t average_p_recon_err = 0.;

      for (int bin = 1; bin <= h_p_evt_fake_NGenInt[*nTrkCut]->GetNbinsX(); bin++) {
        Int_t current_NGenInt = TMath::Nint(h_p_evt_fake_NGenInt[*nTrkCut]->GetBinCenter(bin));
        Float_t poisson_factor = TMath::Exp(-1. * current_mu) * TMath::Power(current_mu, current_NGenInt) / TMath::Factorial(current_NGenInt);

        Float_t current_p_fake = h_p_evt_fake_NGenInt[*nTrkCut]->GetBinContent(bin);
        Float_t current_p_fake_err = h_p_evt_fake_NGenInt[*nTrkCut]->GetBinError(bin);
        if (current_p_fake >= 0) {
          average_p_fake += current_p_fake * poisson_factor;
          average_p_fake_err += TMath::Power(poisson_factor * current_p_fake_err, 2);
        }

        Float_t current_p_recon = h_p_evt_recon_NGenInt[*nTrkCut]->GetBinContent(bin);
        Float_t current_p_recon_err = h_p_evt_recon_NGenInt[*nTrkCut]->GetBinError(bin);
        if (current_p_recon >= 0) {
          average_p_recon += current_p_recon * poisson_factor;
          average_p_recon_err += TMath::Power(poisson_factor * current_p_recon_err, 2);
        }
      }

      average_p_fake_err = TMath::Sqrt(average_p_fake_err);
      average_p_recon_err = TMath::Sqrt(average_p_recon_err);

      Float_t average_mu_fake = -1. * TMath::Log(1. - average_p_fake);
      Float_t average_mu_fake_err = -1. * TMath::Log(1. - (average_p_fake + average_p_fake_err)) - average_mu_fake;

      Float_t average_mu_recon = -1. * TMath::Log(1. - average_p_recon);
      Float_t average_mu_recon_err = -1. * TMath::Log(1. - (average_p_recon + average_p_recon_err)) - average_mu_recon;

      tg_evt_mufake_mu[*nTrkCut]->SetPoint(i-1, current_mu, average_mu_fake);
      tg_evt_mufake_mu[*nTrkCut]->SetPointError(i-1, 0., average_mu_fake_err);

      tg_evt_murecon_mu[*nTrkCut]->SetPoint(i-1, current_mu, average_mu_recon);
      tg_evt_murecon_mu[*nTrkCut]->SetPointError(i-1, 0., average_mu_recon_err);

      tg_evt_mufake_murecon[*nTrkCut]->SetPoint(i-1, average_mu_recon, average_mu_fake);
      tg_evt_mufake_murecon[*nTrkCut]->SetPointError(i-1, average_mu_recon_err, average_mu_fake_err);

      if (i % 100 == 0) {
        cout << "mu = " << current_mu << endl;
        cout << "average_p_recon = " << average_p_recon << " +/- " << average_p_recon_err << endl;
        cout << "average_mu_recon = " << average_mu_recon << " + " << average_mu_recon_err << " - " << average_mu_recon - -1. * TMath::Log(1. - (average_p_recon - average_p_recon_err)) << endl;
      }
    }
  }
#endif

  /* True fakes */
  
  //Added the VtxLumi ones
  std::vector<TString> realfake_samples;
  realfake_samples.push_back("mc_7TeV_17.2_normal_pythia8_pu_bs45");
  realfake_samples.push_back("mc_7TeV_17.2_normal_pythia8_pu_bs55");
  realfake_samples.push_back("mc_7TeV_17.2_VtxLumi_pythia8_pu_bs45");
  realfake_samples.push_back("mc_7TeV_17.2_VtxLumi_pythia8_pu");
  realfake_samples.push_back("mc_8TeV_17.2_normal_pythia8_pu");
  realfake_samples.push_back("mc_8TeV_17.2_VtxLumi_2newsets");

  if (find(realfake_samples.begin(), realfake_samples.end(), tag) != realfake_samples.end()) {

    std::map<Int_t, TH1D*> h_truefakes_NGenInt;
    std::map<Int_t, TGraphErrors*> tg_truefakemu_NGenInt;

    for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {

      TString name;

      name = "hist/h_truefakes_NGenInt_NTrk"; name += *nTrkCut;
      h_truefakes_NGenInt[*nTrkCut] = (TH1D*)f_in->Get(name);

      tg_truefakemu_NGenInt[*nTrkCut] = new TGraphErrors(h_truefakes_NGenInt[*nTrkCut]->GetNbinsX());
      name = "tg_truefakemu_NGenInt_NTrk"; name += *nTrkCut;
      tg_truefakemu_NGenInt[*nTrkCut]->SetName(name);

      for (int bin = 1; bin <= h_truefakes_NGenInt[*nTrkCut]->GetNbinsX(); bin++) {
        Float_t c_ngenint = h_truefakes_NGenInt[*nTrkCut]->GetBinCenter(bin);
        Int_t c_truefakes = h_truefakes_NGenInt[*nTrkCut]->GetBinContent(bin);
        Int_t c_ntrig = h_NTrig_NGenInt[*nTrkCut]->GetBinContent(bin);
        Float_t c_truefakemu = 1. * c_truefakes / c_ntrig;

        Float_t c_truefakemu_err = TMath::Sqrt(c_truefakes) / c_ntrig;

        tg_truefakemu_NGenInt[*nTrkCut]->SetPoint(bin-1, c_ngenint, c_truefakemu);
        tg_truefakemu_NGenInt[*nTrkCut]->SetPointError(bin-1, 0., c_truefakemu_err);

      }
    }
    PlotNTrkGraphs(&tg_truefakemu_NGenInt, "NTrk", true, GlobalSettings::path_fakeCorrection + TString("/") + tag, "Number of generated interactions", "#mu_{true fake}");
  }

  // Drawing
  cout << "Drawing..." << endl;

  TCanvas *c_fakefraction_MuReconMC = new TCanvas("c_fakefraction_MuReconMC", "c_fakefraction_MuReconMC", 800, 800);
  c_fakefraction_MuReconMC->SetRightMargin(0.2);
  TLegend *l_fakefraction = new TLegend(0.81, 0.3, 0.99, 0.7);
  l_fakefraction->SetFillColor(0);
  l_fakefraction->SetBorderSize(1);

  draw_first = true;

  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    TString draw_options;
    if (draw_first) {
      draw_first = false;
      draw_options = "ap";
    } else {
      draw_options = "p";
    }

    c_fakefraction_MuReconMC->cd();
    tg_fake_fraction_MuReconMC[*nTrkCut]->SetMarkerColor(nTrkColors[*nTrkCut]);
    tg_fake_fraction_MuReconMC[*nTrkCut]->SetMarkerStyle(20);
    tg_fake_fraction_MuReconMC[*nTrkCut]->GetXaxis()->SetTitle("#mu_{vis} (after masking correction)");
    tg_fake_fraction_MuReconMC[*nTrkCut]->GetYaxis()->SetTitle("Fake fraction");
    tg_fake_fraction_MuReconMC[*nTrkCut]->Draw(draw_options);

    TString legend_entry = "NTrk"; legend_entry += *nTrkCut;
    l_fakefraction->AddEntry(tg_fake_fraction_MuReconMC[*nTrkCut], legend_entry, "p");
  }

  c_fakefraction_MuReconMC->cd();
  l_fakefraction->Draw();
  c_fakefraction_MuReconMC->SaveAs(GlobalSettings::path_fakeCorrection + TString("/") + tag + "/" + c_fakefraction_MuReconMC->GetName() + TString(".eps"));
  c_fakefraction_MuReconMC->SaveAs(GlobalSettings::path_fakeCorrection + TString("/") + tag + "/" + c_fakefraction_MuReconMC->GetName() + TString(".pdf"));


  // Make canvases
  PlotNTrkGraphs(&h_mu_fake_NGenInt, "NTrk", true, GlobalSettings::path_fakeCorrection + TString("/") + tag);
  PlotNTrkGraphs(&tg_mu_fake_mu, "NTrk", true, GlobalSettings::path_fakeCorrection + TString("/") + tag, "#mu_{inel}", "#mu_{fake}");
  PlotNTrkGraphs(&tg_fake_fraction_MuReconMC, "NTrk", true, GlobalSettings::path_fakeCorrection + TString("/") + tag, "#mu_{rec} * f_{mask}(#mu_{rec})", "Fake fraction");
  PlotNTrkGraphs(&tg_mu_fake_MuReconMC, "NTrk", true, GlobalSettings::path_fakeCorrection + TString("/") + tag, "#mu_{rec} * f_{mask}(#mu_{rec})", "#mu_{fake}");
  // Write TGraphs to cache

  TString path_out = GlobalSettings::path_fakeCorrection; path_out += "/"; path_out += tag; path_out += "/fakerates.root";
  TFile *f_out = new TFile(path_out, "RECREATE");
  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    h_z[*nTrkCut]->Write();
    h_NTrig_NGenInt[*nTrkCut]->Write();
    h_fakes_NGenInt[*nTrkCut]->Write();
    h_recon_NGenInt[*nTrkCut]->Write();
    h_mu_fake_NGenInt[*nTrkCut]->Write();
    h_mu_recon_NGenInt[*nTrkCut]->Write();
    tg_mu_fake_mu[*nTrkCut]->Write();
    tg_mu_recon_mu[*nTrkCut]->Write();
    tg_mu_fake_MuReconMC[*nTrkCut]->Write();
    tg_fake_fraction_MuReconMC[*nTrkCut]->Write();
    h_fake_event_probability_NGenInt[*nTrkCut]->Write();
    h_recon_event_probability_NGenInt[*nTrkCut]->Write();
    tg_fake_event_probability_mu[*nTrkCut]->Write();
    tg_recon_event_probability_mu[*nTrkCut]->Write();
    tg_pfake_vs_precon[*nTrkCut]->Write();
    tg_nevt_mu_real_vs_mu_recon[*nTrkCut]->Write();
    tg_nevt_realfraction_vs_mu_recon[*nTrkCut]->Write();
    /*
    h_fakes_NGenInt[*nTrkCut]->Write();
    h_all_NGenInt[*nTrkCut]->Write();
    h_NTrig_NGenInt[*nTrkCut]->Write();
    tg_fakefraction_NGenInt[*nTrkCut]->Write();
    tg_fakemu_NGenInt[*nTrkCut]->Write();
    tg_fakefraction_mu[*nTrkCut]->Write();
    tg_fakemu_mu[*nTrkCut]->Write();
    h_fakes_NVtxRecon[*nTrkCut]->Write();
    h_all_NVtxRecon[*nTrkCut]->Write();
    h_NTrig_NVtxRecon[*nTrkCut]->Write();
    tg_fakefraction_NVtxRecon[*nTrkCut]->Write();
    tg_fakemu_NVtxRecon[*nTrkCut]->Write();
    tg_fakefraction_MuRecon[*nTrkCut]->Write();
    tg_fakemu_MuRecon[*nTrkCut]->Write();
    tg_fakefraction_MuReconMC[*nTrkCut]->Write();
    tg_fakemu_MuReconMC[*nTrkCut]->Write();
    h_fake_events_NGenInt[*nTrkCut]->Write();
     h_fake_events_NVtxRecon[*nTrkCut]->Write();
    h_fake_event_probability_NGenInt[*nTrkCut]->Write();
    h_fake_event_probability_NVtxRecon[*nTrkCut]->Write();
    h_p_evt_fake_NGenInt[*nTrkCut]->Write();
    h_p_evt_recon_NGenInt[*nTrkCut]->Write();
    tg_evt_mufake_murecon[*nTrkCut]->Write();
    tg_evt_mufake_mu[*nTrkCut]->Write();
    tg_evt_murecon_mu[*nTrkCut]->Write();
    */
    //tg_fake_event_probability_mu[*nTrkCut]->Write();
    //tg_fake_event_probability_MuRecon[*nTrkCut]->Write();
    //tg_fake_event_probability_MuReconMC[*nTrkCut]->Write();
    //tg_fake_vs_total_event_probability[*nTrkCut]->Write();

  }
  f_out->Close();

}

TCanvas* PlotNTrkGraphs(std::map<Int_t, TGraphErrors*> *graphs, TString legend_prefix, bool save_pdf, TString save_path, TString x_title, TString y_title) {

  cout << "[debug] PlotNTrkGraphs() (TGraphErrors)" << endl;

  // Name canvas based on the first TGraph name
  TCanvas *c;
  bool first = true;
  //TString x_title;
  //TString y_title;

  Double_t xmin = 10000000;
  Double_t xmax = -10000000;
  Double_t ymin = 10000000;
  Double_t ymax = -10000000;

  std::map<Int_t, Int_t> colors;
  colors[0] = kBlue+3;
	colors[1] = kRed;
  colors[2] = kBlue;
  colors[3] = kBlack;
  colors[4] = kGreen;
  colors[5] = kGreen+3;
  colors[6] = kCyan+2;
  colors[7] = kYellow;


  std::map<Int_t, Int_t> markers;
  markers[0] = 20;
  markers[1] = 21;
  markers[2] = 22;
  markers[3] = 23;
  markers[4] = 24;
  markers[5] = 25;
  markers[6] = 26;
  markers[7] = 27;
  /*markers[0] = 7;
  markers[1] = 7;
  markers[2] = 7;*/

  for (std::map<Int_t, TGraphErrors*>::iterator graph = (*graphs).begin(); graph != (*graphs).end(); ++graph) {

    if (first) {
      first = false;
      TString cname = (*graph).second->GetName();
      if (cname.BeginsWith("tg_")) {
        cname.Remove(0, 3);
        //cname.ReplaceAll("_NTrk2", "");
				cname.ReplaceAll("_NTrk3", "");
				cname.ReplaceAll("_NTrk4", "");
				cname.ReplaceAll("_NTrk5", "");
        cname.ReplaceAll("_NTrk6", "");
        cname.ReplaceAll("_NTrk7", "");
        cname.ReplaceAll("_NTrk8", "");
        cname.ReplaceAll("_NTrk10", "");
      }
      cname.Prepend("c_");
      c = new TCanvas(cname, cname, 800, 800);

      //x_title = (*graph).second->GetXaxis()->GetTitle();
      //y_title = (*graph).second->GetYaxis()->GetTitle();
    }

    for (int i = 0; i < (*graph).second->GetN(); i++) {
      Double_t x = (*graph).second->GetX()[i];
      Double_t y = (*graph).second->GetY()[i];
      Double_t ex = (*graph).second->GetEX()[i];
      Double_t ey = (*graph).second->GetEY()[i];

      if (x + ex > xmax) xmax = x + ex;
      if (x - ex < xmin) xmin = x - ex;
      if (y + ey > ymax) ymax = y + ey;
      if (y - ey < ymin) ymin = y - ey;
    }
  }

  TH1F *frame = new TH1F("frame", "frame", 100, xmin - 0.05 * (xmax - xmin), xmax + 0.4 * (xmax - xmin));
  frame->SetMinimum(ymin - 0.05 * (ymax - ymin));
  cout << "Setting minimum to " << ymin - 0.05 * (ymax - ymin) << endl;
  frame->SetMaximum(ymax + 0.05 * (ymax - ymin));
  cout << "Setting maximum to " << ymax + 0.05 * (ymax - ymin) << endl;
  frame->GetXaxis()->SetTitle(x_title);
  frame->GetXaxis()->SetLabelSize(0.035);
  frame->GetXaxis()->SetTitleSize(0.035);
  frame->GetYaxis()->SetTitle(y_title);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetTitleSize(0.035);
  frame->GetYaxis()->SetTitleOffset(1.6);
  frame->Draw("axis");

  TLegend *l = new TLegend(0.75, 0.2, 0.9, 0.35);
  l->SetFillColor(0);
  l->SetBorderSize(1);

  Int_t counter = -1;
  for (std::map<Int_t, TGraphErrors*>::iterator graph = (*graphs).begin(); graph != (*graphs).end(); ++graph) {
    counter = counter + 1;
    (*graph).second->SetMarkerColor(colors[counter]);
    (*graph).second->SetMarkerStyle(markers[counter]);
    (*graph).second->Draw("p");
    TString legend_entry = legend_prefix; legend_entry += " "; legend_entry += (*graph).first;
    l->AddEntry((*graph).second, legend_entry, "p");
  }
  l->Draw();
  if (save_pdf) {
    c->SaveAs(save_path + "/" + c->GetName() + TString(".pdf"));
  }
  return c;

}

TCanvas* PlotNTrkGraphs(std::map<Int_t, TH1D*> *graphs, TString legend_prefix, bool save_pdf, TString save_path) {

  cout << "[debug] PlotNTrkGraphs called (TH1D)" << endl;

  // Name canvas based on the first TGraph name
  TCanvas *c;
  bool first = true;
  TString x_title;
  TString y_title;

  Double_t xmin = 10000000;
  Double_t xmax = -10000000;
  Double_t ymin = 10000000;
  Double_t ymax = -10000000;

  std::map<Int_t, Int_t> colors;
  colors[0] = kBlue+3;
  colors[1] = kRed;
  colors[2] = kBlue;
  colors[3] = kBlack;
  colors[4] = kGreen;
  colors[5] = kGreen+3;
  colors[6] = kCyan+2;
  colors[7] = kYellow;
  /*colors[0] = kRed;
  colors[1] = kBlue;
  colors[2] = kBlack;*/


  std::map<Int_t, Int_t> markers;
  markers[0] = 20;
  markers[1] = 21;
  markers[2] = 22;
  markers[3] = 23;
  markers[4] = 24;
  markers[5] = 25;
  markers[6] = 26;
  markers[7] = 27;
  /*markers[0] = 8;
  markers[1] = 8;
  markers[2] = 8;*/

  for (std::map<Int_t, TH1D*>::iterator graph = (*graphs).begin(); graph != (*graphs).end(); ++graph) {

    if (first) {
      first = false;
      TString cname = (*graph).second->GetName();
      if (cname.BeginsWith("h_")) {
        cname.Remove(0, 2);
        //cname.ReplaceAll("_NTrk2", "");
				cname.ReplaceAll("_NTrk3", "");
				cname.ReplaceAll("_NTrk4", "");
        cname.ReplaceAll("_NTrk5", "");
        cname.ReplaceAll("_NTrk6", "");
        cname.ReplaceAll("_NTrk7", "");
        cname.ReplaceAll("_NTrk8", "");
        cname.ReplaceAll("_NTrk10", "");
      }
      cname.Prepend("c_");
      c = new TCanvas(cname, cname, 800, 800);

      x_title = (*graph).second->GetXaxis()->GetTitle();
      y_title = (*graph).second->GetYaxis()->GetTitle();
    }

    for (int i = 1; i <= (*graph).second->GetNbinsX(); i++) {
      Double_t c_xmin = (*graph).second->GetXaxis()->GetXmin();
      Double_t c_xmax = (*graph).second->GetXaxis()->GetXmax();
      Double_t y = (*graph).second->GetBinContent(i);
      Double_t ey = (*graph).second->GetBinError(i);

      if (c_xmax > xmax) xmax = c_xmax;
      if (c_xmin < xmin) xmin = c_xmin;
      if (y + ey > ymax) ymax = y + ey;
      if (y - ey < ymin) ymin = y - ey;
    }
  }

  TH1F *frame = new TH1F("frame", "frame", 100, xmin - 0.05 * (xmax - xmin), xmax + 0.4 * (xmax - xmin));
  frame->SetMinimum(ymin - 0.05 * (ymax - ymin));
  cout << "Setting minimum to " << ymin - 0.05 * (ymax - ymin) << endl;
  frame->SetMaximum(ymax + 0.05 * (ymax - ymin));
  cout << "Setting maximum to " << ymax + 0.05 * (ymax - ymin) << endl;
  frame->GetXaxis()->SetTitle(x_title);
  frame->GetXaxis()->SetLabelSize(0.035);
  frame->GetXaxis()->SetTitleSize(0.035);
  frame->GetYaxis()->SetTitle(y_title);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetTitleSize(0.035);
  frame->GetYaxis()->SetTitleOffset(1.6);
  frame->Draw("axis");

  TLegend *l = new TLegend(0.75, 0.2, 0.9, 0.35);
  l->SetFillColor(0);
  l->SetBorderSize(1);

  Int_t counter = -1;
  for (std::map<Int_t, TH1D*>::iterator graph = (*graphs).begin(); graph != (*graphs).end(); ++graph) {
    counter = counter + 1;
    (*graph).second->SetMarkerColor(colors[counter]);
    (*graph).second->SetMarkerStyle(markers[counter]);
    (*graph).second->Draw("p same");
    TString legend_entry = legend_prefix; legend_entry += " "; legend_entry += (*graph).first;
    l->AddEntry((*graph).second, legend_entry, "p");
  }
  l->Draw();
  if (save_pdf) {
    c->SaveAs(save_path + "/" + c->GetName() + TString(".pdf"));
  }
  return c;

}

