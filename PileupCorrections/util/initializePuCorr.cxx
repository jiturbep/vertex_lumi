#define DEBUG

/**
  *  Pileup masking standalone executable: generate a new pmask_dz histogram using either the May VdM scan, or low-pileup MC.
  *  Author: David R. Yu (dryu@lbl.gov)
  *  January 30, 2012
  *  
  **/

#define REDO_DZ
  
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
#include "TList.h"
#include "TObject.h"
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

#include "PileupCorrections/PileupMaskingCorrection.h"
#include "GlobalSettings/GlobalSettings.h"

using namespace std;
  
int main(int argc, char **argv) {

  SetAtlasStyle();

  cout << "[initializePuCorr] INFO : Start masking correction initialization." << endl;
  #ifdef DEBUG
  cout << "[initializePuCorr] DEBUG: Enabling debug." << endl;
  #endif
  // --- Scan for command line parameters
    
  //Input/output detection
  TString tag;
  TString output_prefix;
  bool set_verbose = false;
  
  int c;
  extern char* optarg;
  while (1) {
    int option_index = 0;
    static struct option long_options[] = { 
      {"tag", required_argument, 0, 't'},
      {"outputPrefix",  required_argument,  0,  'o'},
      {"verbose",  no_argument,    0,  'v'},
      {0, 0, 0, 0}
    };
    c = getopt_long(argc, argv, "t:o:v",long_options,&option_index);
    if (c == -1) break;
    switch(c) {
      case 't':
        {
          tag = optarg;
          cout << "[initializePuCorr] INFO : Read in tag " << tag << endl;
          break;
        }
      case 'o':
        {
          output_prefix = optarg;
          break;
        }
      case 'v':
        {
          set_verbose = true;
          break;
        }
    }
  } // End input parameter scanning.

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
  mc_8TeV_samples.push_back("mc_8TeV_17.2_VtxLumi_BothSamples"); //Added 03aug2014
  mc_8TeV_samples.push_back("mc_8TeV_17.2_VtxLumi_LowMuSample"); //''
  mc_8TeV_samples.push_back("mc_8TeV_17.2_VtxLumi_HighMuSample"); //''

  std::vector<TString> data_samples;
  data_samples.push_back("data_7TeV_17.2_normal");
  data_samples.push_back("data_8TeV_17.2_normal");
  data_samples.push_back("data_7TeV_17.2_VtxLumi");
  data_samples.push_back("data_8TeV_17.2_VtxLumi");
  data_samples.push_back("data_8TeV_17.2_VtxLumi_201351");
  data_samples.push_back("data_8TeV_17.2_VtxLumi_207216");
  data_samples.push_back("data_8TeV_17.2_VtxLumi_207219");
  data_samples.push_back("data_8TeV_17.2_VtxLumi_214984");
  data_samples.push_back("data_8TeV_17.2_VtxLumi_215021");
  

  TString input_type;
  TString input_path;
  TString output_path;
  TString data_run;
  TString mc_energy;

  if (find(mc_7TeV_samples.begin(), mc_7TeV_samples.end(), tag) != mc_7TeV_samples.end()) {
    input_type = "mc";
    input_path = GlobalSettings::path_D3PDMCResults; input_path += tag; input_path += "/"; input_path += GlobalSettings::path_D3PDMCResults_v; input_path += "/InDetTrackD3PD_results.root";
    mc_energy = "7";
  } else if (find(mc_8TeV_samples.begin(), mc_8TeV_samples.end(), tag) != mc_8TeV_samples.end()) {
    input_type = "mc";
    input_path = GlobalSettings::path_D3PDMCResults; input_path += tag; input_path += "/"; input_path += GlobalSettings::path_D3PDMCResults_v; input_path += "/InDetTrackD3PD_results.root";
    mc_energy = "8";
  } else if (tag == "data_7TeV_17.2_normal") {
    input_type = "data";
    data_run = "182013";
    input_path = GlobalSettings::path_inputRawCount; input_path += GlobalSettings::path_VdM_prefix; input_path += "182013/17.2-normal/"; input_path += GlobalSettings::path_inputRawCount_v; input_path += "InDetTrackD3PD_results.root";
  } else if (tag == "data_7TeV_17.2_VtxLumi") {
    input_type = "data";
    data_run = "182013";
    input_path = GlobalSettings::path_inputRawCount; input_path += GlobalSettings::path_VdM_prefix; input_path += "182013/17.2-VtxLumi/"; input_path += GlobalSettings::path_inputRawCount_v; input_path += "InDetTrackD3PD_results.root";
  } else if (tag == "data_8TeV_17.2_normal") {
    input_type = "data";
    data_run = "201351";
    input_path = GlobalSettings::path_inputRawCount; input_path += "/VdMScan-201351/17.2-normal/"; input_path += GlobalSettings::path_inputRawCount_v; input_path += "/InDetTrackD3PD_results.root";
  } else if (tag == "data_8TeV_17.2_VtxLumi_201351") {
    input_type = "data";
    data_run = "201351";
    //input_path = GlobalSettings::path_inputRawCount; input_path += "/VdMScan-201351/17.2-VtxLumi/"; input_path += GlobalSettings::path_inputRawCount_v; input_path += "/InDetTrackD3PD_results_21042014_NormalD3PD.root";
    input_path = GlobalSettings::path_inputRawCount; input_path += "/VdMScan-201351/17.2-VtxLumi/"; input_path += GlobalSettings::path_inputRawCount_v; input_path += "/InDetTrackD3PD_results_NormalD3PD_04Aug2014.root";
  } else if (tag == "data_8TeV_17.2_VtxLumi_207216") {
    input_type = "data";
    data_run = "207216";
    //input_path = GlobalSettings::path_inputRawCount; input_path += "/VdMScan-207216/17.2-VtxLumi/"; input_path += GlobalSettings::path_inputRawCount_v; input_path += "/InDetTrackD3PD_results_21042014_NormalD3PD.root";
    input_path = GlobalSettings::path_inputRawCount; input_path += "/VdMScan-207216/17.2-VtxLumi/"; input_path += GlobalSettings::path_inputRawCount_v; input_path += "/InDetTrackD3PD_results_NormalD3PD_05Aug2014.root";
  } else if (tag == "data_8TeV_17.2_VtxLumi_207219") {
    input_type = "data";
    data_run = "207219";
    //input_path = GlobalSettings::path_inputRawCount; input_path += "/VdMScan-207219/17.2-VtxLumi/"; input_path += GlobalSettings::path_inputRawCount_v; input_path += "/InDetTrackD3PD_results_21042014_NormalD3PD.root";
    input_path = GlobalSettings::path_inputRawCount; input_path += "/VdMScan-207219/17.2-VtxLumi/"; input_path += GlobalSettings::path_inputRawCount_v; input_path += "/InDetTrackD3PD_results_NormalD3PD_05Aug2014.root";
  } else if (tag == "data_8TeV_17.2_VtxLumi_214984") {
    input_type = "data";
    data_run = "214984";
    //input_path = GlobalSettings::path_inputRawCount; input_path += "/VdMScan-214984/17.2-VtxLumi/"; input_path += GlobalSettings::path_inputRawCount_v; input_path += "/InDetTrackD3PD_results_21042014_NormalD3PD.root";
    input_path = GlobalSettings::path_inputRawCount; input_path += "/VdMScan-214984/17.2-VtxLumi/"; input_path += GlobalSettings::path_inputRawCount_v; input_path += "/InDetTrackD3PD_results_NormalD3PD_05Aug2014.root";
  } else if (tag == "data_8TeV_17.2_VtxLumi_215021") {
    input_type = "data";
    data_run = "215021";
    //input_path = GlobalSettings::path_inputRawCount; input_path += "/VdMScan-215021/17.2-VtxLumi/"; input_path += GlobalSettings::path_inputRawCount_v; input_path += "/InDetTrackD3PD_results_21042014_NormalD3PD.root";
    input_path = GlobalSettings::path_inputRawCount; input_path += "/VdMScan-215021/17.2-VtxLumi/"; input_path += GlobalSettings::path_inputRawCount_v; input_path += "/InDetTrackD3PD_results.root";
  } else {
    cerr << "Input not recognized: " << tag << ". Exiting..." << endl;
    exit(1);
  }

  cout << "[initializePuCorr] Input path = "  << input_path << endl;

  if (output_prefix == "") {
    output_prefix = tag;
  }

  cout << "[initializePuCorr] INFO : Saving results with prefix " << output_prefix << endl;

  TFile *f_in = new TFile(input_path, "READ");

  #ifdef DEBUG
  TString debug_path = GlobalSettings::path_maskingCorrection; debug_path += "/"; debug_path += output_prefix; debug_path += "/debug.root";
  TFile *f_debug = new TFile(debug_path, "RECREATE");
  cout << "[initializePuCorr] DEBUG : Writing debug histograms to " << debug_path << endl;
  #endif

  std::vector<Int_t> nTrkCuts;
  //nTrkCuts.push_back(2); 
  nTrkCuts.push_back(3);
  nTrkCuts.push_back(4);
  nTrkCuts.push_back(5);
  //nTrkCuts.push_back(6);
  //nTrkCuts.push_back(7);
  //nTrkCuts.push_back(8);
  //nTrkCuts.push_back(10);

  // Fork code into data and MC here
  if (input_type == "data") {
    //Some run-specific stuff
    std::map<Int_t, Color_t> marker_colors;
    std::vector<Int_t> bcids;

    if (data_run == "182013") {
      bcids.push_back(81);
      bcids.push_back(867);
      bcids.push_back(2752);

      marker_colors[81] = kRed - 3;
      marker_colors[867] = kGreen - 3;
      marker_colors[2752] = kBlue + 3;
    } else if (data_run == "200805") {
      bcids.push_back(1);
      bcids.push_back(149);
      bcids.push_back(334);

      marker_colors[1] = kRed - 3;
      marker_colors[149] = kGreen - 3;
      marker_colors[334] = kBlue + 3;
    } //Added VdMScan 201351
      else if (data_run == "201351") {
      bcids.push_back(1);
      bcids.push_back(241);
      bcids.push_back(2881);
      bcids.push_back(3121);

      marker_colors[1] = kRed - 3;
      marker_colors[241] = kGreen - 3;
      marker_colors[2881] = kBlue + 3;
      marker_colors[3121] = kBlack;
    } //Added July 2012 VdMScan 207216
      else if (data_run == "207216") {
      bcids.push_back(1);
      bcids.push_back(721);
      bcids.push_back(1821);

      marker_colors[1] = kRed - 3;
      marker_colors[721] = kGreen - 3;
      marker_colors[1821] = kBlue + 3;
    } //Added July 2012 VdMScan 207219
      else if (data_run == "207219") {
      bcids.push_back(1);
      bcids.push_back(721);
      bcids.push_back(1821);

      marker_colors[1] = kRed - 3;
      marker_colors[721] = kGreen - 3;
      marker_colors[1821] = kBlue + 3;
    } //Added November 2012 VdMScan 214984
      else if (data_run == "214984") {
      bcids.push_back(1);
      bcids.push_back(2361);
      bcids.push_back(2881);

      marker_colors[1] = kRed - 3;
      marker_colors[2361] = kGreen - 3;
      marker_colors[2881] = kBlue + 3;
    } //Added November 2012 VdMScan 215021
      else if (data_run == "215021") {
      bcids.push_back(1);
      bcids.push_back(2361);
      bcids.push_back(2881);

      marker_colors[1] = kRed - 3;
      marker_colors[2361] = kGreen - 3;
      marker_colors[2881] = kBlue + 3;
    }

    std::vector<std::pair<Int_t, Int_t> > good_pLBs;

    if (data_run == "182013") {
      good_pLBs.push_back(std::make_pair<Int_t, Int_t>(38, 66));
      good_pLBs.push_back(std::make_pair<Int_t, Int_t>(116, 152));
      good_pLBs.push_back(std::make_pair<Int_t, Int_t>(202, 223));
      good_pLBs.push_back(std::make_pair<Int_t, Int_t>(273, 293));
      good_pLBs.push_back(std::make_pair<Int_t, Int_t>(343, 357));
      good_pLBs.push_back(std::make_pair<Int_t, Int_t>(469, 481));
      good_pLBs.push_back(std::make_pair<Int_t, Int_t>(531, 700));
    } else if (data_run == "200805") {
      good_pLBs.push_back(std::make_pair<Int_t,Int_t>(179, 211));
    } else if (data_run == "201351") {
      good_pLBs.push_back(std::make_pair<Int_t,Int_t>(18, 133)); //Scan1
      good_pLBs.push_back(std::make_pair<Int_t,Int_t>(140, 251)); //Scan2
      good_pLBs.push_back(std::make_pair<Int_t,Int_t>(354, 464)); //Scan3
      //Low <mu> (from run query)
      //good_pLBs.push_back(std::make_pair<Int_t,Int_t>(18, 64)); 
      //good_pLBs.push_back(std::make_pair<Int_t,Int_t>(78, 133));
      //good_pLBs.push_back(std::make_pair<Int_t,Int_t>(140, 164));
      //good_pLBs.push_back(std::make_pair<Int_t,Int_t>(173, 251));
      //good_pLBs.push_back(std::make_pair<Int_t,Int_t>(354, 442));
      //good_pLBs.push_back(std::make_pair<Int_t,Int_t>(454, 457));
      //good_pLBs.push_back(std::make_pair<Int_t,Int_t>(114, 131)); //Test
    } else if (data_run == "207216") {
      good_pLBs.push_back(std::make_pair<Int_t,Int_t>(55, 172)); //Scan4  ?? It was wrong, should be fine now
      good_pLBs.push_back(std::make_pair<Int_t,Int_t>(290, 401)); //Scan5 
      good_pLBs.push_back(std::make_pair<Int_t,Int_t>(409, 516)); //Scan6
      //good_pLBs.push_back(std::make_pair<Int_t,Int_t>(74, 91)); //Test
    } else if (data_run == "207219") {
      good_pLBs.push_back(std::make_pair<Int_t,Int_t>(86, 118)); //Scan8
    } else if (data_run == "214984") {
      good_pLBs.push_back(std::make_pair<Int_t,Int_t>(10, 68)); //Scan10
      good_pLBs.push_back(std::make_pair<Int_t,Int_t>(75, 197)); //Scan11
      good_pLBs.push_back(std::make_pair<Int_t,Int_t>(451, 561)); //Scan14
    } else if (data_run == "215021") {
      good_pLBs.push_back(std::make_pair<Int_t,Int_t>(4, 114)); //Scan15
    }

    std::vector<Int_t> pLB_list;
    
    for (vector<std::pair<Int_t, Int_t> >::iterator pLB_interval = good_pLBs.begin(); pLB_interval != good_pLBs.end(); ++pLB_interval) {
      for (Int_t current_pLB = (*pLB_interval).first; current_pLB <= (*pLB_interval).second; current_pLB++) {
        pLB_list.push_back(current_pLB);
      }
    }

    Int_t n_pLBs = pLB_list.size();

    // - Make maps of the beamspot sigma_z values. Maps are BCID : pLB : sigma_z.
    std::map<Int_t, std::map<Int_t, Double_t> > sigma_z;
    std::map<Int_t, std::map<Int_t, bool> > low_stats;
    std::map<Int_t, TGraphErrors*> tg_sigma_z;
    std::map<Int_t, TGraphErrors*> tg_z;
    std::map<Int_t, TGraph*> tg_chi2ndf;
    
    for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
      TString name;
      
      name = "tg_sigma_z_BCID"; name += *bcid;
      tg_sigma_z[*bcid] = new TGraphErrors(n_pLBs);
      tg_sigma_z[*bcid]->SetName(name);
      
      name = "tg_z_BCID"; name += *bcid;
      tg_z[*bcid] = new TGraphErrors(n_pLBs);
      tg_z[*bcid]->SetName(name);
      
      name = "tg_chi2ndf_BCID"; name += *bcid;
      tg_chi2ndf[*bcid] = new TGraphErrors(n_pLBs);
      tg_chi2ndf[*bcid]->SetName(name);
      
      name = "hist/PriVtxZpLB_BCID"; name += *bcid;
      TH2D *h_z_plb_tmp = (TH2D*)f_in->Get(name);

      for (vector<Int_t>::iterator current_pLB = pLB_list.begin(); current_pLB != pLB_list.end(); ++current_pLB) {     
        Int_t current_point = distance(pLB_list.begin(), current_pLB);
        //if ((current_point % 20) == 0) {
          cout << "On point " << current_point << " / " << n_pLBs << endl;
        //}
      
        TString tag = "BCID"; tag += *bcid; tag += "_pLB"; tag += *current_pLB;
        
        Int_t ybin_z = h_z_plb_tmp->GetYaxis()->FindBin(*current_pLB);
        TH1D *h_z_tmp = (TH1D*)h_z_plb_tmp->ProjectionX("_x", ybin_z, ybin_z);
        h_z_tmp->SetName(TString("h_z_") + tag);
        h_z_tmp->Sumw2();
        
        if (h_z_tmp->Integral() < 500) {
          low_stats[*bcid][*current_pLB] = true;
          continue;
        }
        
        low_stats[*bcid][*current_pLB] = false;
        
        TF1 *f_gaus = new TF1("f_gaus", "gaus(0)", -200., 200.);
        f_gaus->SetParameter(0, h_z_tmp->Integral());
        f_gaus->SetParameter(1, 0.);
        f_gaus->SetParameter(2, 50.);
        h_z_tmp->Fit(f_gaus, "QRE+");
        
        sigma_z[*bcid][*current_pLB] = f_gaus->GetParameter(2);
        cout << "sigma_z["<<*bcid<<"]["<<*current_pLB<<"] = " << sigma_z[*bcid][*current_pLB] << endl;
        
        tg_sigma_z[*bcid]->SetPoint(current_point, *current_pLB, f_gaus->GetParameter(2));
        tg_sigma_z[*bcid]->SetPointError(current_point, 0, f_gaus->GetParError(2));
        tg_z[*bcid]->SetPoint(current_point, *current_pLB, f_gaus->GetParameter(1));
        tg_z[*bcid]->SetPointError(current_point, 0, f_gaus->GetParError(1));
        tg_chi2ndf[*bcid]->SetPoint(current_point, *current_pLB, f_gaus->GetChisquare()/f_gaus->GetNDF());
      
        #ifdef DEBUG
        f_debug->cd();
        h_z_tmp->Write();
        #endif
      }
      #ifdef DEBUG
      f_debug->cd();
      tg_sigma_z[*bcid]->Write();
      tg_z[*bcid]->Write();
      tg_chi2ndf[*bcid]->Write();
      #endif
    }
    #ifdef DEBUG
    f_debug->Close();
    #endif
    
    // Make whole-BCID z-distributions (outside of scans)
    std::map<Int_t, TH1D*> h_z;
    for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
      TString name;
      name = "hist/PriVtxZpLB_BCID"; name += *bcid;
      TH2D *h_z_plb_tmp = (TH2D*)f_in->Get(name);

      name = "h_z_BCID"; name += *bcid;
      h_z[*bcid] = new TH1D(name, name, h_z_plb_tmp->GetXaxis()->GetNbins(), h_z_plb_tmp->GetXaxis()->GetXmin(),h_z_plb_tmp->GetXaxis()->GetXmax());
    
      //Add up slices of the TH2D
      cout << "[initializePuCorr] Add up slices of the TH2D to make z distribution histogram" << endl;
      for (vector<std::pair<Int_t, Int_t> >::iterator pLB_interval = good_pLBs.begin(); pLB_interval != good_pLBs.end(); ++pLB_interval) {
        Int_t bin1 = h_z_plb_tmp->GetYaxis()->FindBin((*pLB_interval).first);
        Int_t bin2 = h_z_plb_tmp->GetYaxis()->FindBin((*pLB_interval).second);
        h_z[*bcid]->Add(h_z_plb_tmp->ProjectionX("tmp", bin1+10, bin2-10)); // Hack for April
        cout << "bin1 = " << bin1 << ", bin2 = " << bin2 << "h_zRMS = " << h_z[*bcid]->GetRMS() << endl;
      }
      h_z[*bcid]->Sumw2();

      TString path_z = "data_z_distribution_BCID"; path_z+=*bcid; path_z+=".root";
      TFile *f_z;
      f_z = new TFile(path_z,"RECREATE");
      h_z[*bcid]->Write();

      /*TFile *f_z = new TFile(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/z_distribution.root"), "UPDATE");
      h_z_copy[*bcid] = new TH1D(name, name, h_z_plb_tmp->GetXaxis()->GetNbins(), h_z_plb_tmp->GetXaxis()->GetXmin(),h_z_plb_tmp->GetXaxis()->GetXmax());
      h_z_copy[*bcid] = h_z[*bcid]->Clone();
      h_z_copy[*bcid]->Scale(1/h_z_copy[*bcid]); 
      TString nameh = "Z distribution, BCID"; nameh+=*bcid;
      h_z_copy[*bcid]->SetTitle(nameh);
      h_z_copy[*bcid]->GetXaxis()->SetTitle("Z [mm]");
      h_z_copy[*bcid]->Write();
      f_z->Close();*/
    }
      
    //Proceed to PMCs
    
    std::map<Int_t, std::map<Int_t, PileupMaskingCorrection*> > pmc;
    std::map<Int_t, std::map<Int_t, TH1D*> > h_pmask_dz, h_dz_rebinned;
    std::map<Int_t, TH1D*> h_pmask_dz_avg;
    std::map<Int_t, TCanvas*> c_pmask_dz;
    std::map<Int_t, TLegend*> l_pmask_dz;
    std::map<Int_t, std::map<Int_t, TCanvas*> > c_dz;
    std::map<Int_t, std::map<Int_t, TLegend*> > l_dz;
    std::map<Int_t, std::map<Int_t, TGraphErrors*> > tg_masking_correction;

    bool first = true;
    
    for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
      //if (*nTrkCut != 5) continue;
      cout << "[initializePuCorr] ************************************************************************" << endl;
      cout << "[initializePuCorr] *** NTrk = " << *nTrkCut << " ***" << endl;
      cout << "[initializePuCorr] ************************************************************************" << endl;
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        //if (*bcid == 867 || *bcid == 2752) continue;
        cout << "[initializePuCorr] ************************************************************************" << endl;
        cout << "[initializePuCorr] *** BCID = " << *bcid << "***" << endl;
        cout << "[initializePuCorr] ************************************************************************" << endl;

        TString tag1 = "BCID"; tag1 += *bcid; tag1 += "_NTrkCut"; tag1 += *nTrkCut;
        
        // - Make DZ histogram
        
        TString name = "hist/VtxDzTightTight_pLB_"; name += tag1;
        TH2D *h2 = (TH2D*)f_in->Get(name);
        cout << "[initializePuCorr] VtxDzTightTight has " << h2->GetEntries() << "entries." << endl;
                    
        TString tag = "BCID"; tag += *bcid; tag += "_NTrkCut"; tag += *nTrkCut;
      
        TH1D *h_dz = new TH1D(TString("h_dz_") + tag, TString("h_dz_") + tag, h2->GetXaxis()->GetNbins(), h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax());      
        TH1D *h_dz_copy;    
        
        for (vector<std::pair<Int_t, Int_t> >::iterator pLB_interval = good_pLBs.begin(); pLB_interval != good_pLBs.end(); ++pLB_interval) {
          Int_t bin1 = h2->GetYaxis()->FindBin((*pLB_interval).first);
          Int_t bin2 = h2->GetYaxis()->FindBin((*pLB_interval).second);
          h_dz->Add(h2->ProjectionX("tmp2", bin1, bin2));
        }
        h_dz->Sumw2();
        //Creating debug file
        TString path_deltaz = "data_deltaz_NTrkCut"; path_deltaz+=*nTrkCut; path_deltaz+=".root";
        TFile *f_deltaz;
        f_deltaz = new TFile(path_deltaz,"RECREATE");
        h_dz->Write();

        /*TFile *f_dz = new TFile(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/dz_distribution.root"), "UPDATE");
        h_dz_copy = (TH1D*)h_dz->Clone();
        //h_dz_copy->Scale(1/h_dz_copy->Integral());
        h_dz_copy->Rebin(20);
      	TString nameh = "deltaZ distribution"; nameh+=tag;
      	h_dz_copy->SetTitle(nameh);
      	h_dz_copy->GetXaxis()->SetTitle("deltaZ [mm]");
      	h_dz_copy->Write();
      	f_dz->Close();*/
        
        //Initialize PileupMaskingCorrection
        cout << "[initializePuCorr] INFO: Initializing PileupCorrection " << endl;
        pmc[*nTrkCut][*bcid] = new PileupMaskingCorrection(h_dz, tag);
        
        // - Expected delta Z distribution
        #ifdef REDO_DZ
        cout << "[initializePuCorr] INFO: REDO_DZ " << endl;
        pmc[*nTrkCut][*bcid]->GenerateDzDistribution(h_z[*bcid]);
        TFile *f_dz_expected = new TFile(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/dz_expected.root"), "UPDATE");
        TH1D *h_tmp = (TH1D*)pmc[*nTrkCut][*bcid]->GetExpectedDzDistribution();
        name = "h_dz_expected_BCID"; name += *bcid; name += "_NTrk"; name += *nTrkCut;
        h_tmp->SetName(name);
        h_tmp->Write();
        f_dz_expected->Close();
        #else
        cout << "[initializePuCorr] INFO: Loading dz distribution " << endl;
        name = "h_dz_expected_BCID"; name += *bcid; name += "_NTrk"; name += *nTrkCut;
        TH1D *h_tmp = (TH1D*)f_dz_expected->Get(name);
        pmc[*nTrkCut][*bcid]->LoadDzDistribution(h_tmp);
        #endif
        
        // - Make p_mask vs. delta Z histogram
        pmc[*nTrkCut][*bcid]->GenerateNewPmask();
        
        // - Calculate total p_mask
        #ifdef DEBUG
        cout << "[initializePuCorr] DEBUG : Expected Delta Z distribution:" << endl;
        cout << "[initializePuCorr] DEBUG : \tBins: " << pmc[*nTrkCut][*bcid]->GetExpectedDzDistribution()->GetNbinsX() << endl;
        cout << "[initializePuCorr] DEBUG : \tRange: " << pmc[*nTrkCut][*bcid]->GetExpectedDzDistribution()->GetXaxis()->GetXmin() << ", " << pmc[*nTrkCut][*bcid]->GetExpectedDzDistribution()->GetXaxis()->GetXmax() << endl;
        cout << "[initializePuCorr] DEBUG : p_mask dz distribution:" << endl;
        cout << "[initializePuCorr] DEBUG : \tBins: " << pmc[*nTrkCut][*bcid]->GetDifferentialPmask()->GetNbinsX() << endl;
        cout << "[initializePuCorr] DEBUG : \tRange : " << pmc[*nTrkCut][*bcid]->GetDifferentialPmask()->GetXaxis()->GetXmin() << ", " << pmc[*nTrkCut][*bcid]->GetDifferentialPmask()->GetXaxis()->GetXmax() << endl;
        #endif

        pmc[*nTrkCut][*bcid]->GenerateCorrection(pmc[*nTrkCut][*bcid]->GetExpectedDzDistribution());
        cout << "[initializePuCorr] INFO : Total p_mask = " << pmc[*nTrkCut][*bcid]->GetTotalPmask() << " +/- " << pmc[*nTrkCut][*bcid]->GetTotalPmaskError() << endl;
        h_pmask_dz[*nTrkCut][*bcid] = (TH1D*)pmc[*nTrkCut][*bcid]->GetDifferentialPmask()->Clone();

        tg_masking_correction[*nTrkCut][*bcid] = (TGraphErrors*)pmc[*nTrkCut][*bcid]->GetMuCorrection()->Clone();
        TString tgname = "tg_masking_correction_NTrk"; tgname += *nTrkCut; tgname += "_BCID"; tgname += *bcid;
        tg_masking_correction[*nTrkCut][*bcid]->SetName(tgname);
        
        if (first) {
          first = false;
          pmc[*nTrkCut][*bcid]->Save(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/pmask_cache.root"), tag, true);
        } else {
          pmc[*nTrkCut][*bcid]->Save(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/pmask_cache.root"), tag, false);
        }
      }

      // Make an average p_mask over BCIDs
      TString name = "h_pmask_dz_NTrk"; name += *nTrkCut;
      bool first_histogram = true;
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        if (first_histogram) {
          h_pmask_dz_avg[*nTrkCut] = (TH1D*)h_pmask_dz[*nTrkCut][*bcid]->Clone();
          h_pmask_dz_avg[*nTrkCut]->SetName(name);
          h_pmask_dz_avg[*nTrkCut]->SetTitle(name);
          first_histogram = false;
        } else {
          h_pmask_dz_avg[*nTrkCut]->Add(h_pmask_dz[*nTrkCut][*bcid]);
        }
      }
      h_pmask_dz_avg[*nTrkCut]->Scale(1./bcids.size());
      TFile *f_out = new TFile(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/pmask_cache.root"), "UPDATE");
      h_pmask_dz_avg[*nTrkCut]->Write();

      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        tg_masking_correction[*nTrkCut][*bcid]->Write();
      }
      //Draw all three p_mask_dz functions + average
      name = "c_pmask_dz_NTrk"; name += *nTrkCut;
      
      c_pmask_dz[*nTrkCut] = new TCanvas(name, name, 1200, 800);
      l_pmask_dz[*nTrkCut] = new TLegend(0.75,0.75,0.94,0.94);
      l_pmask_dz[*nTrkCut]->SetFillColor(0);
      l_pmask_dz[*nTrkCut]->SetBorderSize(1);
      
      h_pmask_dz_avg[*nTrkCut]->SetLineColor(kBlack);
      h_pmask_dz_avg[*nTrkCut]->SetMinimum(-0.1);
      h_pmask_dz_avg[*nTrkCut]->SetMaximum(1.1);
      h_pmask_dz_avg[*nTrkCut]->Rebin(2);
      h_pmask_dz_avg[*nTrkCut]->Scale(0.5);
      h_pmask_dz_avg[*nTrkCut]->GetXaxis()->SetRange(h_pmask_dz_avg[*nTrkCut]->GetXaxis()->FindBin(-50), h_pmask_dz_avg[*nTrkCut]->GetXaxis()->FindBin(50));
      h_pmask_dz_avg[*nTrkCut]->GetXaxis()->SetTitle("#Deltaz (mm)");
      h_pmask_dz_avg[*nTrkCut]->GetYaxis()->SetTitle("p_{mask}(#Deltaz)");
      h_pmask_dz_avg[*nTrkCut]->Draw("hist");
      l_pmask_dz[*nTrkCut]->AddEntry(h_pmask_dz_avg[*nTrkCut], "Avg.", "l");
      
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        h_pmask_dz[*nTrkCut][*bcid]->SetMarkerColor(marker_colors[*bcid]);
        h_pmask_dz[*nTrkCut][*bcid]->SetMarkerStyle(20);
        h_pmask_dz[*nTrkCut][*bcid]->Rebin(2);
        h_pmask_dz[*nTrkCut][*bcid]->Scale(0.5);
        h_pmask_dz[*nTrkCut][*bcid]->Draw("same");
        TString legend_entry = "BCID "; legend_entry += *bcid;
        l_pmask_dz[*nTrkCut]->AddEntry(h_pmask_dz[*nTrkCut][*bcid], legend_entry, "p");
      }
      
      l_pmask_dz[*nTrkCut]->Draw();
      
      c_pmask_dz[*nTrkCut]->Write();
      c_pmask_dz[*nTrkCut]->SaveAs(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/") + c_pmask_dz[*nTrkCut]->GetName() + TString(".png"));
      c_pmask_dz[*nTrkCut]->SaveAs(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/") + c_pmask_dz[*nTrkCut]->GetName() + TString(".eps"));
      
      // - Draw expected dz fit
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {

        name = "c_dz_NTrk"; name += *nTrkCut; name += "_BCID"; name += *bcid;
        c_dz[*nTrkCut][*bcid] = new TCanvas(name, name, 1200, 800);
        l_dz[*nTrkCut][*bcid] = new TLegend(0.7, 0.7, 0.95, 0.95);
        l_dz[*nTrkCut][*bcid]->SetFillColor(0);
        l_dz[*nTrkCut][*bcid]->SetBorderSize(1);
        h_dz_rebinned[*nTrkCut][*bcid] = (TH1D*)pmc[*nTrkCut][*bcid]->h_dz_rebinned;
        h_dz_rebinned[*nTrkCut][*bcid]->GetXaxis()->SetTitle("#Deltaz (mm)");
        h_dz_rebinned[*nTrkCut][*bcid]->SetMaximum(h_dz_rebinned[*nTrkCut][*bcid]->GetMaximum() * 1.15);
        h_dz_rebinned[*nTrkCut][*bcid]->GetXaxis()->SetRangeUser(-300,300);
        h_dz_rebinned[*nTrkCut][*bcid]->Draw("hist func");
        name = "fit_gaussian_excluded_BCID"; name += *bcid; name += "_NTrkCut"; name += *nTrkCut;
        TF1 *fit = (TF1*)h_dz_rebinned[*nTrkCut][*bcid]->GetFunction(name);
        l_dz[*nTrkCut][*bcid]->AddEntry(h_dz_rebinned[*nTrkCut][*bcid], "Observed #Deltaz", "l");
        if( fit != 0 ) {
          fit->SetLineColor(2);
          fit->SetLineWidth(1);
          fit->Draw("same");
          l_dz[*nTrkCut][*bcid]->AddEntry(fit, "Fit to sides", "l");
        }
        l_dz[*nTrkCut][*bcid]->Draw();
        c_dz[*nTrkCut][*bcid]->SaveAs(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/") + c_dz[*nTrkCut][*bcid]->GetName() + TString(".png"));
        c_dz[*nTrkCut][*bcid]->SaveAs(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/") + c_dz[*nTrkCut][*bcid]->GetName() + TString(".eps"));
      }

      f_out->Close();
      
    }
    
    f_in->Close();  
    
  } else if (input_type == "mc") {
    
    cout << "Initializing pileup masking correction for MC." << endl;
    
    std::map<Int_t, PileupMaskingCorrection*> pmc;
    std::map<Int_t, TH1D*> h_pmask_dz, h_dz_rebinned;
    std::map<Int_t, TCanvas*> c_pmask_dz, c_dz;
    std::map<Int_t, TLegend*> l_pmask_dz, l_dz;

    bool first = true;
    
    for (std::vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    
      cout << endl << "////////////////// [initializePuCorr] INFO : Computing masking correction for NTrk " << *nTrkCut << endl;
      cout << endl;
      cout << endl;
      TString name;
      
      // Z and DZ histograms
      //name = "hist/h_dz_NGenInt_NTrk"; name += *nTrkCut;
      name = "hist/h_all_dz_NGenInt_NTrk"; name += *nTrkCut;
      f_in->cd();
      TH2D *h_dz_NGenInt = (TH2D*)f_in->Get(name);
      cout << "[initializePuCorr] INFO: h_dz_NGenInt->GetNbinsX() " << h_dz_NGenInt->GetNbinsX() << endl;
      cout << "[initializePuCorr] INFO: h_dz_NGenInt->GetNbinsY() " << h_dz_NGenInt->GetNbinsY() << endl;
      TH1D *h_dz;
      TH1D *h_tmp_scaled;
      TH1D *h_tmp_copy;
      Double_t weight = 0.;
      //Double_t MaxNGenInt = h_dz_NGenInt->GetNbinsY();
      //Double_t MaxNGenInt = 40; // MC Low Mu sample
      Double_t MaxNGenInt = 100; // MC High Mu sample

      
      Float_t mu;
      if (mc_energy == "7") {
        mu = 1.97;
      } else if (mc_energy == "8") {
        mu = 5; // HACK: change this to 0.5 once you have enough statistics at low mu to deal with it!
      }
      TString path_deltaz = "MC_deltaz_NTrkCut"; path_deltaz+=*nTrkCut; path_deltaz+=".root";
      TFile *f_deltaz;
      f_deltaz = new TFile(path_deltaz,"RECREATE");
      cout << "[initializePuCorr] DEBUG: MaxNGenInt = " << MaxNGenInt << endl;
      for (Int_t current_ngenint = 2; current_ngenint <11; current_ngenint++) {
        Int_t bin = h_dz_NGenInt->GetYaxis()->FindBin(current_ngenint);
        TString tmpname = "h_dz_NGenInt"; tmpname += current_ngenint; tmpname += "_NTrk"; tmpname += *nTrkCut;
        TH1D *h_tmp = (TH1D*)h_dz_NGenInt->ProjectionX(tmpname, bin, bin);
        h_tmp->Sumw2();
        if (current_ngenint < MaxNGenInt){
          h_tmp_copy = (TH1D*)h_tmp->Clone();
          h_tmp_copy->Rebin(10);
          //if (h_tmp_copy->Integral() != 0) h_tmp_copy->Scale(poisson_factor / h_tmp_copy->Integral());
          h_tmp_copy->GetXaxis()->SetLimits(-200.0,200.0);
          h_tmp_copy->Write();
        }

        if (current_ngenint == 2) {
          h_dz = (TH1D*)h_tmp->Clone();
          TString name3 = "h_dz_NTrk"; name3 += *nTrkCut;
          h_dz->SetName(name3);
          //h_dz->Scale(h_tmp->Integral());
          weight += h_tmp->Integral();
        } else {
          h_dz->Add(h_tmp,1);
          //h_dz->Add(h_tmp, h_tmp->Integral());
          weight+=h_tmp->Integral();
        }
      }

      //h_dz->Scale(1.0/weight);
      h_dz->Scale(1./h_dz->Integral());
      h_dz->Write();
      TH1D *h_dz_copy = (TH1D*)h_dz->Clone();
      h_dz_copy->SetName((h_dz->GetName())+TString("_copy"));
      h_dz_copy->Rebin(10);
      h_dz_copy->Write();
      cout << "weight = " << weight << endl;      

      cout << "[initializePuCorr] DEBUG : h_dz has " << h_dz->Integral() << " entries" << endl;
      cout << "[initializePuCorr] DEBUG : h_dz has " << h_dz->GetNbinsX() << " X bins" << endl;

      /*//The old version of the code, however, created the delta z distribution reweighting with a poisson distribution
      for (Int_t current_ngenint = 2; current_ngenint < MaxNGenInt; current_ngenint++) {
        Float_t poisson_factor = TMath::Exp(-1. * mu) * TMath::Power(mu, current_ngenint) / TMath::Factorial(current_ngenint);
        Int_t bin = h_dz_NGenInt->GetYaxis()->FindBin(current_ngenint);
        TString tmpname = "h_dz_NGenInt"; tmpname += current_ngenint; tmpname += "_NTrk"; tmpname += *nTrkCut;
        TH1D *h_tmp = (TH1D*)h_dz_NGenInt->ProjectionX(tmpname, bin, bin);
        h_tmp->Sumw2();
        if (current_ngenint < MaxNGenInt){
          h_tmp_copy = (TH1D*)h_tmp->Clone();
          h_tmp_copy->Rebin(10);
          if (h_tmp_copy->Integral() != 0) h_tmp_copy->Scale(poisson_factor / h_tmp_copy->Integral());
          h_tmp_copy->GetXaxis()->SetLimits(-200.0,200.0);
          h_tmp_copy->Write();
        }

        if (h_tmp->Integral() != 0) h_tmp->Scale(poisson_factor / h_tmp->Integral());

        if (current_ngenint == 2) {
          h_dz = (TH1D*)h_tmp->Clone();
          TString name3 = "h_dz_NTrk"; name3 += *nTrkCut;
          h_dz->SetName(name3);
        } else {
          h_dz->Add(h_tmp, 1);
        }
      }

      h_dz->Scale(1./h_dz->Integral());
      h_dz->Write();
      TH1D *h_dz_copy = (TH1D*)h_dz->Clone();
      h_dz_copy->SetName((h_dz->GetName())+TString("_copy"));
      h_dz_copy->Rebin(10);
      h_dz_copy->Write();

      cout << "[initializePuCorr] DEBUG : h_dz has " << h_dz->Integral() << " entries" << endl;
      cout << "[initializePuCorr] DEBUG : h_dz has " << h_dz->GetNbinsX() << " X bins" << endl;*/
    
      ////////////// Z distribution  
      
      TH2D *h_z_NGenInt = (TH2D*)f_in->Get("hist/h_privtx_z_NGenInt");
      TH1D *h_z;
      TH1D *h_tmp_copy2;
      Double_t weight1 = 0.;

      TString path_z = "MC_z_distribution_NTrkCut"; path_z+=*nTrkCut; path_z+=".root";
      TFile *f_z;
      f_z = new TFile(path_z,"RECREATE");
      
      for (Int_t current_ngenint = 1; current_ngenint < 11; current_ngenint++) {
        Float_t poisson_factor = TMath::Exp(-1. * mu) * TMath::Power(mu, current_ngenint) / TMath::Factorial(current_ngenint);
        //cout << "[INFO] Line 724, Originally, poisson factor is "<< poisson_factor << ", but I am setting it to 1" << endl;
        //poisson_factor=1;

        Int_t bin = h_z_NGenInt->GetYaxis()->FindBin(current_ngenint);
        TString tmpname = "h_z_NGenInt"; tmpname += current_ngenint; tmpname += "_NTrk"; tmpname += *nTrkCut;
        TH1D *h_tmp = (TH1D*)h_z_NGenInt->ProjectionX(tmpname, bin, bin);
        //h_tmp->Print("all");
        h_tmp->Sumw2();
        if (current_ngenint < 11){
          h_tmp_copy2 = (TH1D*)h_tmp->Clone();
          h_tmp_copy2->Rebin(10);
          //if (h_tmp_copy2->Integral() != 0) h_tmp_copy2->Scale(poisson_factor / h_tmp_copy2->Integral());
          h_tmp_copy2->GetXaxis()->SetLimits(-200.0,200.0);
          h_tmp_copy2->Write();
        }

        //if (h_tmp->Integral() != 0) h_tmp->Scale(poisson_factor / h_tmp->Integral());

        if (current_ngenint == 1) {
          TString name3 = "h_z_NTrk"; name3 += *nTrkCut;
          h_z = (TH1D*)h_tmp->Clone();
          h_z->SetName(name3);
          //h_z->Scale(h_tmp->Integral());
          //weight1+=h_tmp->Integral();
        } else {
          h_z->Add(h_tmp,1);
          //h_z->Add(h_tmp,h_tmp->Integral());
          //weight1+=h_tmp->Integral();
        }
      }
      //h_z->Scale(1./weight1);
      //cout << "[initializePuCorr] DEBUG : h_z has " << h_z->Integral() << " entries" << endl;
      h_z->Scale(1./h_z->Integral());
      //h_z->Rebin(10);
      h_z->Write();
      TH1D *h_z_copy = (TH1D*)h_z->Clone();
      h_z_copy->SetName((h_z->GetName())+TString("_copy"));
      //h_z_copy->Scale(1./h_z_copy->Integral());
      h_z_copy->Rebin(20);
      h_z_copy->Write();
      //cout << "weight1 = " << weight1 << endl;
      //f_z->Close();

      cout << "[initializePuCorr] DEBUG : h_z has " << h_z->Integral() << " entries" << endl;

      //Pileup masking correction
      
      TString tag = "NTrk"; tag += *nTrkCut;
      pmc[*nTrkCut] = new PileupMaskingCorrection(h_dz, tag);
      pmc[*nTrkCut]->low_stats = false;
      pmc[*nTrkCut]->rebin_factor = 10;
      pmc[*nTrkCut]->is_MC = kTRUE;
      pmc[*nTrkCut]->MaxNGenInt = MaxNGenInt;
      
      // - Generate expected delta Z distribution
      pmc[*nTrkCut]->GenerateDzDistribution(h_z,tag);
      // - Make p_mask vs. delta Z histogram
      pmc[*nTrkCut]->GenerateNewPmask();
      
      // - Calculate total p_mask
      pmc[*nTrkCut]->GenerateCorrection(pmc[*nTrkCut]->GetExpectedDzDistribution());
      cout << "\t\t Total p_mask = " << pmc[*nTrkCut]->GetTotalPmask() << " +/- " << pmc[*nTrkCut]->GetTotalPmaskError() << endl;
      h_pmask_dz[*nTrkCut] = (TH1D*)pmc[*nTrkCut]->GetDifferentialPmask()->Clone();
      
      #ifdef DEBUG
      cout << "[initializePuCorr] DEBUG : Expected Delta Z distribution:" << endl;
      cout << "[initializePuCorr] DEBUG : \tBins: " << pmc[*nTrkCut]->GetExpectedDzDistribution()->GetNbinsX() << endl;
      cout << "[initializePuCorr] DEBUG : \tRange: " << pmc[*nTrkCut]->GetExpectedDzDistribution()->GetXaxis()->GetXmin() << ", " << pmc[*nTrkCut]->GetExpectedDzDistribution()->GetXaxis()->GetXmax() << endl;
      cout << "[initializePuCorr] DEBUG : p_mask dz distribution:" << endl;
      cout << "[initializePuCorr] DEBUG : \tBins: " << pmc[*nTrkCut]->GetDifferentialPmask()->GetNbinsX() << endl;
      cout << "[initializePuCorr] DEBUG : \tRange : " << pmc[*nTrkCut]->GetDifferentialPmask()->GetXaxis()->GetXmin() << ", " << pmc[*nTrkCut]->GetDifferentialPmask()->GetXaxis()->GetXmax() << endl;
      #endif
      
      // - Start saving
      
      if (first) {
        first = false;
        pmc[*nTrkCut]->Save(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/pmask_cache.root"), tag, true);
        if (set_verbose) {
          pmc[*nTrkCut]->SaveDebugHistograms(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/debug.root"), "", true);
        }
      } else {
        pmc[*nTrkCut]->Save(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/pmask_cache.root"), tag, false);
        if (set_verbose) {
          pmc[*nTrkCut]->SaveDebugHistograms(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/debug.root"), "", false);
        }
      }


      TFile *f_out = new TFile(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/pmask_cache.root"), "UPDATE");
      h_pmask_dz[*nTrkCut]->Write();
      
      // - Draw pmask_dz
      
      name = "c_pmask_dz_NTrk"; name += *nTrkCut;
      c_pmask_dz[*nTrkCut] = new TCanvas(name, name, 1200, 800);
      name = "p_{mask} vs. #Deltaz, NTrk"; name += *nTrkCut;
      h_pmask_dz[*nTrkCut]->SetTitle(name);
      h_pmask_dz[*nTrkCut]->Rebin(8);
      h_pmask_dz[*nTrkCut]->Scale(1./8.);
      h_pmask_dz[*nTrkCut]->GetXaxis()->SetTitle("#Deltaz (mm)");
      h_pmask_dz[*nTrkCut]->GetYaxis()->SetTitle("p_{mask}");
      h_pmask_dz[*nTrkCut]->GetXaxis()->SetRange(h_pmask_dz[*nTrkCut]->GetXaxis()->FindBin(-50), h_pmask_dz[*nTrkCut]->GetXaxis()->FindBin(50));
      h_pmask_dz[*nTrkCut]->SetMinimum(-0.1);
      h_pmask_dz[*nTrkCut]->SetMaximum(1.1);

      h_pmask_dz[*nTrkCut]->Draw("hist");
      
      c_pmask_dz[*nTrkCut]->SaveAs(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/") + c_pmask_dz[*nTrkCut]->GetName() + TString(".png"));
      c_pmask_dz[*nTrkCut]->Write();

      // - Draw expected dz fit
      name = "c_dz_NTrk"; name += *nTrkCut;
      c_dz[*nTrkCut] = new TCanvas(name, name, 1200, 800);
      l_dz[*nTrkCut] = new TLegend(0.7, 0.7, 0.95, 0.95);
      l_dz[*nTrkCut]->SetFillColor(0);
      l_dz[*nTrkCut]->SetBorderSize(1);
      h_dz_rebinned[*nTrkCut] = (TH1D*)pmc[*nTrkCut]->h_dz_rebinned;
      h_dz_rebinned[*nTrkCut]->GetXaxis()->SetTitle("#Deltaz (mm)");
      h_dz_rebinned[*nTrkCut]->GetXaxis()->SetRangeUser(-300,300);
      h_dz_rebinned[*nTrkCut]->Draw("hist func");
      name = "fit_gaussian_excluded_NTrk"; name += *nTrkCut;
      TF1 *fit = (TF1*)h_dz_rebinned[*nTrkCut]->GetFunction(name);
      l_dz[*nTrkCut]->AddEntry(h_dz_rebinned[*nTrkCut], "Observed #Deltaz", "l");
      if( fit != 0 ) {
        fit->SetLineColor(2);
        fit->SetLineWidth(1);
        fit->Draw("same");
        l_dz[*nTrkCut]->AddEntry(fit, "Fit to sides", "l");
      }
      l_dz[*nTrkCut]->Draw();
      c_dz[*nTrkCut]->SaveAs(GlobalSettings::path_maskingCorrection + TString("/") + output_prefix + TString("/") + c_dz[*nTrkCut]->GetName() + TString(".png"));
      
      f_out->Close();
    }
  }
  cout << "Done!" << endl;
}
