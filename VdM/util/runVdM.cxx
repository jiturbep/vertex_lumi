// Executable for VdM calibration of vertex-based luminosity measurement.
//TESTMESSAGE
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <stdlib.h>

#include "TROOT.h"
#include "TF1.h"
#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TLeafI.h"
#include "TLeafF.h"
#include "TStopwatch.h"
#include "TMath.h"

#include "GlobalSettings/GlobalSettings.h"
#include "VdM/VanDerMeerAnalysis.h"

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<double>+;
#else
template class vector<double>;
#endif

//---helpers
void usage();
//--- Main function
int main(int argc, char **argv) {
  TString run, settings;
  bool deleteOldHistograms = false;
  int doAnalysisMethod = GlobalSettings::kVtxC; //default: vertex counting
  //bool redo_pileup_corrections = false;
  bool verbose = false;

  /**
    *  Systematic uncertainty evaluation: each of these options will apply some change to the method, e.g. vary the mu scaling parameter, and write out the resulting visible cross section to <settings>-<uncertainty name>
    */
  std::map<TString, bool> systematic_uncertainty_list;
  systematic_uncertainty_list["fake_high"] = false; // Highest fake mu from mu scaling
  systematic_uncertainty_list["fake_low"] = false; // Lowest fake mu from mu scaling
  systematic_uncertainty_list["masking_toy_scaling"] = false; // Scale p_mask up by 2%, as suggested by the simulation.
  systematic_uncertainty_list["bs-45mm"] = false;
  systematic_uncertainty_list["bs-55mm"] = false;

  // --- Scan for command line parameters
  int c;
  //extern int optind;
  extern char* optarg;
  while (1) {
    int option_index = 0;
    static struct option long_options[] = { {0,0,0,0} };
    c=getopt_long(argc, argv, "veu:nr:s:fh",long_options,&option_index);
    if (c == -1) {
      break;
    }
    switch(c) {
      case 'r': {
        run = optarg;
        break;
      }
      case 's': {
        settings = optarg;
        break;
      }
      case 'n': {
        cout << "Creating new histogram file" << endl;
        deleteOldHistograms = true;
        break;
      }
      case 'e': {
        cout << "Doing event-based analysis" << endl;
        doAnalysisMethod = GlobalSettings::kEvtC;
        break;
      }
      case 'u': {
        if (systematic_uncertainty_list.find(optarg) == systematic_uncertainty_list.end()) {
          cerr << "Unknown systematic uncertainty specified: " << optarg << ". Exiting..." << endl;
          exit(1);
        } else {
          cout << "Systematic uncertainty mode: " << optarg << endl;
          systematic_uncertainty_list[optarg] = true;
        }
        break;
      }
      case 'v': {
        cout << "Writing out debugging histograms." << endl;
        verbose = true;
        break;
      }
      case 'f': {
        cout << "Doing Unfolding-based method."<< endl;
        doAnalysisMethod = GlobalSettings::kUnfC;
        break;
      }
      case 'h':
      default:
        usage();
        return 0;
    }
  }

  if (run == "") {
    cerr << "Error: must specify run. Exiting..." << endl;
    usage();
    exit(1);
  } else if ((run != "178013") && (run != "178064") && (run != "182013") && (run != "201351") && (run != "207216") && (run != "207219") && (run != "214984") && (run != "215021")) {
    cerr << "Error: run not recognized: " << run << ". Exiting..." << endl;
    usage();
    exit(1);
  }

  if (settings == "") {
    cerr << "Error: must specify settings. Exiting..." << endl;
    usage();
    exit(1);
  } else if ((settings != "16.X-normal") && (settings != "17.2-normal") && (settings != "17.2-VtxLumi")) {
    cerr << "Error: settings not recognized: " << settings << ". Exiting..." << endl;
    usage();
    exit(1);
  }

  //Load global settings
  GlobalSettings gs; //load defaults
  //std::string prefix( run == "182013" ? GlobalSettings::path_VdM_prefix : GlobalSettings::path_mu_prefix );

  vector<Int_t> nTrkCuts;
  //if (doAnalysisMethod == GlobalSettings::kEvtC) {
  //nTrkCuts.push_back(2);
  //}
  //nTrkCuts.push_back(2);
  nTrkCuts.push_back(3);
  nTrkCuts.push_back(4);
  nTrkCuts.push_back(5);
  //nTrkCuts.push_back(6);
  //nTrkCuts.push_back(7);
  //nTrkCuts.push_back(8);
  //nTrkCuts.push_back(10);

  vector<TString> fit_functions;
  //fit_functions.push_back("sgcl");
  fit_functions.push_back("sgc");
  //fit_functions.push_back("sg");
  //fit_functions.push_back("dg");
  fit_functions.push_back("dgc");
  //fit_functions.push_back("dgcl");
  //fit_functions.push_back("spline");

  // Input and output
  stringstream s_path_vertex_counts, s_path_vertex_histograms;
  s_path_vertex_counts << gs.path_inputRawCount << "/" << "VdMScan-" << run << "/" << settings << "/" << gs.path_inputRawCount_v << gs.fname_inputRawCount_Tree;
  s_path_vertex_histograms << gs.path_inputRawCount << "/" << "VdMScan-" << run << "/" << settings << "/" << gs.path_inputRawCount_v << gs.fname_inputRawCount_Histo;

  stringstream s_path_results;
  stringstream s_path_debug;
  s_path_results << gs.path_outputVdM << "/" << "VdMScan-" << run << "/" << settings;
  for (std::map<TString, bool>::iterator it = systematic_uncertainty_list.begin(); it != systematic_uncertainty_list.end(); ++it) {
    if ((*it).second) {
      s_path_results << "_" << (*it).first;
    }
  }
  s_path_results << "/";
  if (doAnalysisMethod == GlobalSettings::kVtxC) {
    s_path_results << "NVtx/";
  } else if (doAnalysisMethod == GlobalSettings::kEvtC) {
    s_path_results << "NEvt/";
  } else {
    s_path_results << "NUnf/";
  }
  s_path_debug << s_path_results.rdbuf(); //the same up to the very last
  s_path_results << gs.fname_outputVdM;
  s_path_debug << gs.fname_outputDebugVdM;

  // -- Clear old files
  TString ts_fileToClear;
  TFile *f_fileToClear;
  if (deleteOldHistograms) {
    f_fileToClear = new TFile(TString(s_path_results.str()),"RECREATE");
    f_fileToClear->Close();
    f_fileToClear = new TFile(TString(s_path_debug.str()), "RECREATE");
    f_fileToClear->Close();
  }


  // Run-specific things
  vector<Int_t> bcidList;
  vector<Int_t> scans;
  std::map<Int_t, TString> x_timestamps, y_timestamps;
  TString trigger_type;
  vector<TString> path_scan_ntuple;
  TString file_deadtime;
  TString file_deadtime_2;
  TString energy;

  if (run == "182013") {
    cout << "[RunVdM] INFO : Starting VdM analysis of run 182013" << endl;

    bcidList.push_back(81);
    bcidList.push_back(867);
    bcidList.push_back(2752);

    scans.push_back(1);
    scans.push_back(2);

    x_timestamps[1] = gs.path_timestamps + "/182013/scan/x1_timestamps.dat";
    x_timestamps[2] = gs.path_timestamps + "/182013/scan/x2_timestamps.dat";
    y_timestamps[1] = gs.path_timestamps + "/182013/scan/y1_timestamps.dat";
    y_timestamps[2] = gs.path_timestamps + "/182013/scan/y2_timestamps.dat";
    
    file_deadtime = gs.path_deadtime + "/dt_182013.root";
    path_scan_ntuple.push_back(gs.path_lumiNtuples + "/2011MayScan1Raw-v8.root");

    trigger_type = "MBTS";

    energy = "7";

  } else if (run == "201351") {
    cout << "[RunVdM] INFO : Starting VdM analysis of run 201351" << endl;

    bcidList.push_back(1);
    bcidList.push_back(241);
    bcidList.push_back(2881);
    bcidList.push_back(3121);

    scans.push_back(1);
    scans.push_back(2);
    scans.push_back(3);

    x_timestamps[1] = gs.path_timestamps + "/201351/scan/x1_timestamps.dat";
    x_timestamps[2] = gs.path_timestamps + "/201351/scan/x2_timestamps.dat";
    x_timestamps[3] = gs.path_timestamps + "/201351/scan/x3_timestamps.dat";
    y_timestamps[1] = gs.path_timestamps + "/201351/scan/y1_timestamps.dat";
    y_timestamps[2] = gs.path_timestamps + "/201351/scan/y2_timestamps.dat";
    y_timestamps[3] = gs.path_timestamps + "/201351/scan/y3_timestamps.dat";
    
    //file_deadtime = gs.path_deadtime + "/Julia/201351_1_trig_L1_MBTS_2_BGRP7.root";
    //file_deadtime_2 = gs.path_deadtime + "/Julia/201351_2_trig_L1_MBTS_2_BGRP7.root";
    file_deadtime = gs.path_deadtime + "/James/201351_trig_L1_MBTS_2_BGRP7.root";
    path_scan_ntuple.push_back(gs.path_lumiNtuples + "/2012AprScan1Raw-v13.root");
    path_scan_ntuple.push_back(gs.path_lumiNtuples + "/2012AprScan2Raw-v13.root");

    //TODO: Deadtime corrections (SP)

    //trigger_type = "BGRP7";
	  trigger_type = "MBTS";

    energy = "8";
  }
  
  else if (run == "207216") {
    cout << "[RunVdM] INFO : Starting VdM analysis of run 207216" << endl;

    bcidList.push_back(1);
    bcidList.push_back(721);
    bcidList.push_back(1821);

    scans.push_back(4);
    scans.push_back(5);
    scans.push_back(6);

    x_timestamps[4] = gs.path_timestamps + "/207216/scan/x4_timestamps.dat";
    x_timestamps[5] = gs.path_timestamps + "/207216/scan/x5_timestamps.dat";
    x_timestamps[6] = gs.path_timestamps + "/207216/scan/x6_timestamps.dat";
    y_timestamps[4] = gs.path_timestamps + "/207216/scan/y4_timestamps.dat";
    y_timestamps[5] = gs.path_timestamps + "/207216/scan/y5_timestamps.dat";
    y_timestamps[6] = gs.path_timestamps + "/207216/scan/y6_timestamps.dat";
    
    file_deadtime = gs.path_deadtime + "/Julia/207216_trig_L1_MBTS_2_BGRP7.root";
    path_scan_ntuple.push_back(gs.path_lumiNtuples + "/2012JulScan1Raw-v13.root");

    //TODO: Deadtime corrections (SP)

    trigger_type = "MBTS";

    energy = "8";
  }
  
    else if (run == "207219") {
    cout << "[RunVdM] INFO : Starting VdM analysis of run 207219" << endl;

    bcidList.push_back(1);
    bcidList.push_back(721);
    bcidList.push_back(1821);

    scans.push_back(8);

    x_timestamps[8] = gs.path_timestamps + "/207219/scan/x8_timestamps.dat";
    y_timestamps[8] = gs.path_timestamps + "/207219/scan/y8_timestamps.dat";
    
    file_deadtime = gs.path_deadtime + "/Julia/207219_trig_L1_MBTS_2_BGRP7.root";
    path_scan_ntuple.push_back(gs.path_lumiNtuples + "/2012JulScan2Raw-v13.root");

    trigger_type = "MBTS";

    energy = "8";
  }
  
  else if (run == "214984") {
    cout << "[RunVdM] INFO : Starting VdM analysis of run 214984" << endl;

    bcidList.push_back(1);
    bcidList.push_back(2361);
    bcidList.push_back(2881);

    scans.push_back(10);
    scans.push_back(11);
    scans.push_back(14);
    x_timestamps[10] = gs.path_timestamps + "/214984/scan/x10_timestamps.dat";
    y_timestamps[10] = gs.path_timestamps + "/214984/scan/y10_timestamps.dat";
    x_timestamps[11] = gs.path_timestamps + "/214984/scan/x11_timestamps.dat";
    y_timestamps[11] = gs.path_timestamps + "/214984/scan/y11_timestamps.dat";
    x_timestamps[14] = gs.path_timestamps + "/214984/scan/x14_timestamps.dat";
    y_timestamps[14] = gs.path_timestamps + "/214984/scan/y14_timestamps.dat";

    file_deadtime = gs.path_deadtime + "/214984_trig_L1_MBTS_2_BGRP7.root";
    path_scan_ntuple.push_back(gs.path_lumiNtuples + "/2012NovScan1Raw-v13.root");
    path_scan_ntuple.push_back(gs.path_lumiNtuples + "/2012NovScan2Raw-v13.root");

    trigger_type = "MBTS";

    energy = "8";
  }
  
  else if (run == "215021") {
    cout << "[RunVdM] INFO : Starting VdM analysis of run 215021" << endl;

    bcidList.push_back(1);
    bcidList.push_back(2361);
    bcidList.push_back(2881);

    scans.push_back(15);

    x_timestamps[15] = gs.path_timestamps + "/215021/scan/x15_timestamps.dat";
    y_timestamps[15] = gs.path_timestamps + "/215021/scan/y15_timestamps.dat";

    file_deadtime = gs.path_deadtime + "/Julia/215021_trig_L1_MBTS_2_BGRP7.root";
    path_scan_ntuple.push_back(gs.path_lumiNtuples + "/2012NovScan3Raw-v13.root");

    trigger_type = "MBTS";

    energy = "8";
  }

  // -- Output TTree
  TFile *f_out = new TFile(TString(s_path_results.str()), "RECREATE");

  Int_t current_scan, current_bcid;
  Int_t current_ntrkcut;

  std::vector<TString> fit_methods;
  std::vector<double> current_sigma_vis;
  std::vector<double> current_sigma_vis_err;
  std::vector<double> current_lumi_sp;
  std::vector<double> current_lumi_sp_err;
  std::vector<double> current_lumi;
  std::vector<double> current_lumi_err;

  std::vector<TString>*p_fit_methods = &fit_methods;
  std::vector<double> *p_current_sigma_vis = &current_sigma_vis;
  std::vector<double> *p_current_sigma_vis_err = &current_sigma_vis_err;
  std::vector<double> *p_current_lumi_sp = &current_lumi_sp;
  std::vector<double> *p_current_lumi_sp_err = &current_lumi_sp_err;
  std::vector<double> *p_current_lumi = &current_lumi;

  std::vector<double> current_sigma_x;
  std::vector<double> current_sigma_y;
  std::vector<double> current_mu_max_x;
  std::vector<double> current_mu_max_y;
  std::vector<double> current_c_x;
  std::vector<double> current_c_y;
  std::vector<double> current_r_x;
  std::vector<double> current_r_y;
  std::vector<double> current_chi2ndf_x;
  std::vector<double> current_chi2ndf_y;

  std::vector<double> *p_current_sigma_x = &current_sigma_x;
  std::vector<double> *p_current_sigma_y = &current_sigma_y;
  std::vector<double> *p_current_mu_max_x = &current_mu_max_x;
  std::vector<double> *p_current_mu_max_y = &current_mu_max_y;
  std::vector<double> *p_current_c_x = &current_c_x;
  std::vector<double> *p_current_c_y = &current_c_y;
  std::vector<double> *p_current_r_x = &current_r_x;
  std::vector<double> *p_current_r_y = &current_r_y;
  std::vector<double> *p_current_chi2ndf_x = &current_chi2ndf_x;
  std::vector<double> *p_current_chi2ndf_y = &current_chi2ndf_y;

  std::vector<double> current_sigma_x_err;
  std::vector<double> current_sigma_y_err;
  std::vector<double> current_mu_max_x_err;
  std::vector<double> current_mu_max_y_err;
  std::vector<double> current_c_x_err;
  std::vector<double> current_c_y_err;

  std::vector<double> *p_current_sigma_x_err = &current_sigma_x_err;
  std::vector<double> *p_current_sigma_y_err = &current_sigma_y_err;
  std::vector<double> *p_current_mu_max_x_err = &current_mu_max_x_err;
  std::vector<double> *p_current_mu_max_y_err = &current_mu_max_y_err;
  std::vector<double> *p_current_c_x_err = &current_c_x_err;
  std::vector<double> *p_current_c_y_err = &current_c_y_err;


  TString title = "VdM results for run ";
  title += run;
  TTree *t_vdm = new TTree("VdmResults", title);
  t_vdm->Branch("Scan", &current_scan, "Scan/I");
  t_vdm->Branch("BCID", &current_bcid, "BCID/I");
  t_vdm->Branch("NTrkCut", &current_ntrkcut, "NTrkCut/I");
  /*
    TBronch *b_sigma_vis = t_vdm->Bronch("SigmaVis", "vector<float>", &p_current_sigma_vis);
    TBronch *b_sigma_vis_err = t_vdm->Bronch("SigmaVisErr", "vector<float>", &p_current_sigma_vis_err);
    TBronch *b_lumi_sp = t_vdm->Bronch("LumiSp", "vector<float>", &p_current_lumi_sp);
    TBronch *b_lumi_sp_err = t_vdm->Bronch("LumiSpErr", "vector<float>", &p_current_lumi_sp_err);
  */

  t_vdm->Branch("FitMethods", &p_fit_methods);
  t_vdm->Branch("SigmaVis", &p_current_sigma_vis);
  t_vdm->Branch("SigmaVisErr", &p_current_sigma_vis_err);
  t_vdm->Branch("LumiSp", &p_current_lumi_sp);
  t_vdm->Branch("LumiSpErr", &p_current_lumi_sp_err);
  t_vdm->Branch("Lumi", &p_current_lumi_sp);
  t_vdm->Branch("LumiErr", &p_current_lumi_sp_err);

  t_vdm->Branch("SigmaX", &p_current_sigma_x);
  t_vdm->Branch("SigmaY", &p_current_sigma_y);
  t_vdm->Branch("MuMaxX", &p_current_mu_max_x);
  t_vdm->Branch("MuMaxY", &p_current_mu_max_y);
  t_vdm->Branch("CX", &p_current_c_x);
  t_vdm->Branch("CY", &p_current_c_y);
  t_vdm->Branch("RX", &p_current_r_x);
  t_vdm->Branch("RY", &p_current_r_y);
  t_vdm->Branch("Chi2NdfX", &p_current_chi2ndf_x);
  t_vdm->Branch("Chi2NdfY", &p_current_chi2ndf_y);

  t_vdm->Branch("SigmaXErr", &p_current_sigma_x_err);
  t_vdm->Branch("SigmaYErr", &p_current_sigma_y_err);
  t_vdm->Branch("MuMaxXErr", &p_current_mu_max_x_err);
  t_vdm->Branch("MuMaxYErr", &p_current_mu_max_y_err);
  t_vdm->Branch("CXErr", &p_current_c_x_err);
  t_vdm->Branch("CYErr", &p_current_c_y_err);



  for (vector<Int_t>::iterator bcid = bcidList.begin(); bcid != bcidList.end(); ++bcid) {
    current_bcid = *bcid;
    for (vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {
      current_scan = *scan;
      for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
        current_ntrkcut = *nTrkCut;

        cout << endl;
        cout << "[RunVdM] INFO : On BCID " << *bcid << " / Scan " << *scan << " / NTrk " << *nTrkCut << endl;
        cout << endl;

        TString method;
        if (doAnalysisMethod == GlobalSettings::kVtxC) {
          method = "Vtx";
        } else if (doAnalysisMethod == GlobalSettings::kEvtC) {
          method = "Evt";
        } else if (doAnalysisMethod == GlobalSettings::kUnfC) {
          method = "Unf";
        }

        TString save_tag = "BCID";
        save_tag += *bcid;
        save_tag += "_scan";
        save_tag += *scan;
        save_tag += "_NTrk";
        save_tag += *nTrkCut;

        VanDerMeerAnalysis *vdma = new VanDerMeerAnalysis(method, *nTrkCut, save_tag);
        for (map<TString, bool>::iterator it = systematic_uncertainty_list.begin(); it != systematic_uncertainty_list.end(); ++it) {
          if ((*it).second) {
            vdma->SetSystematicUncertaintyFlag((*it).first);
          }
        }

        vdma->LoadPlbTimestamps(x_timestamps[*scan], "x");
        vdma->LoadPlbTimestamps(y_timestamps[*scan], "y");
        vdma->LoadPlbPrescales(file_deadtime);
        //if (run == "182013") {
        //  vdma->prescale = 2; prescales[pLB ]
        //}
        //if (run == "201351") vdma->LoadDeadtime(file_deadtime,file_deadtime_2);
        vdma->LoadDeadtime(file_deadtime);
        for (vector<TString>::iterator it = path_scan_ntuple.begin(); it != path_scan_ntuple.end(); ++it) {
          vdma->LoadNominalSeparations(*it);
          vdma->LoadBunchIntensities(*it, *bcid);
        }

        if (doAnalysisMethod == GlobalSettings::kUnfC) {
          vdma->SetUnfResponseMatrix((gs.path_unf_ResponseMatrix+"/"+gs.fname_unf_ResponseMatrix).c_str(),
          gs.other_unf_ResponseMatrixHNameBase.c_str());
          vdma->LoadVertexCounts(TString(s_path_vertex_histograms.str()), *bcid, trigger_type);
        } else {
          vdma->LoadVertexCounts(TString(s_path_vertex_counts.str()), *bcid, trigger_type);
        }

        vdma->CalculateMuPlb();

        #ifdef RANDOM_TRIGGER
        vdma->ConvertToRandom(2, 9000.);
        #endif

        // -- Pileup corrections
        if (doAnalysisMethod != GlobalSettings::kUnfC) {
          vdma->InitializeFakeCorrection(energy, settings, *nTrkCut);
          TFile *f_histograms = new TFile(TString(s_path_vertex_histograms.str()), "READ");
          TString hname = "hist/PriVtxZpLB_BCID";
          hname += *bcid;
          cout << "[RunVdM] INFO : h_z_plb for masking correction is coming from "<< s_path_vertex_histograms.str() << hname << endl;
          TH2D *h_z_plb = (TH2D*)f_histograms->Get(hname);

          vdma->CorrectPileupEffects(h_z_plb, energy, settings, *nTrkCut, run, *bcid);
        }


        #ifdef REMOVED_051612
        /*
          You need to initialize the masking correction for each pLB individually! If you project over the whole run, you get a wider z distribution!
        */
        if ((*nTrkCut != 2) && (doAnalysisMethod == GlobalSettings::kVtxC)) {
          TFile *f_cache = new TFile("include/cache.root", "READ");
          TString save_tag = "BCID";
          save_tag += *bcid;
          save_tag += "_scan";
          save_tag += *scan;
          save_tag += "_NTrk";
          save_tag += *nTrkCut;
          save_tag += "_";
          save_tag += energy;
          save_tag += "_";
          save_tag += settings;
          TString h_dz_name = "h_dz_expected_";
          h_dz_name += save_tag;

          if (f_cache->Get(h_dz_name) && (!redo_pileup_corrections)) {
            vdma->InitializeMaskingCorrection(energy, settings, *nTrkCut, save_tag);
          } else {
            TFile *f_histograms = new TFile(TString(s_path_vertex_histograms.str()), "READ");
            TString hname = "hist/PriVtxZpLB_BCID";
            hname += *bcid;
            TH2D *h_z_pLB = (TH2D*)f_histograms->Get(hname);
            TH1D *h_z = (TH1D*)h_z_pLB->ProjectionX();
            vdma->InitializeMaskingCorrection(energy, settings, *nTrkCut, h_z, save_tag);
          }
        }
        vdma->CorrectPileupEffects();

        #endif

        vdma->CalculateMuSpecPlb();
        vdma->FitVdmCurves();
	      for (std::vector<TString>::iterator fit_function = fit_functions.begin(); fit_function != fit_functions.end(); ++fit_function){
          if (*scan == 11 || *scan == 5 || *scan == 10){
            vdma->r_x[*fit_function]=1;
            cout << "[runVdM] Scan = " << *scan << ", therefore changed r_x to 1" << endl;
	        } 
	      }
        vdma->Finalize();

        if (verbose) {
          TString save_tag = "BCID";
          save_tag += *bcid;
          save_tag += "_scan";
          save_tag += *scan;
          save_tag += "_NTrk";
          save_tag += *nTrkCut;
          vdma->DebugPlots(TString(s_path_debug.str()), save_tag);
        }

        fit_methods.clear();
        current_sigma_vis.clear();
        current_sigma_vis_err.clear();
        current_lumi_sp.clear();
        current_lumi_sp_err.clear();

        current_sigma_x.clear();
        current_sigma_y.clear();
        current_mu_max_x.clear();
        current_mu_max_y.clear();
        current_c_x.clear();
        current_c_y.clear();
        current_r_x.clear();
        current_r_y.clear();
        current_chi2ndf_x.clear();
        current_chi2ndf_y.clear();

        current_sigma_x_err.clear();
        current_sigma_y_err.clear();
        current_mu_max_x_err.clear();
        current_mu_max_y_err.clear();
        current_c_x_err.clear();
        current_c_y_err.clear();

        for (vector<TString>::iterator it = fit_functions.begin(); it != fit_functions.end(); ++it) {

          fit_methods.push_back(*it);

          current_sigma_vis.push_back(vdma->GetSigmaVis(*it));
          current_sigma_vis_err.push_back(vdma->GetSigmaVisErr(*it));
          current_lumi_sp.push_back(vdma->GetLumiSp(*it));
          current_lumi_sp_err.push_back(vdma->GetLumiSpErr(*it));
          current_lumi.push_back(vdma->GetLumi(*it));
          current_lumi_err.push_back(vdma->GetLumiErr(*it));

          current_sigma_x.push_back(vdma->GetVdmParameter(*it, "Sigma_x"));
          current_sigma_y.push_back(vdma->GetVdmParameter(*it, "Sigma_y"));
          current_mu_max_x.push_back(vdma->GetVdmParameter(*it, "mu_max_x"));
          current_mu_max_y.push_back(vdma->GetVdmParameter(*it, "mu_max_y"));
          current_c_x.push_back(vdma->GetVdmParameter(*it, "c_x"));
          current_c_y.push_back(vdma->GetVdmParameter(*it, "c_y"));
          current_r_x.push_back(vdma->r_x[*it]);
          current_r_y.push_back(vdma->r_y[*it]);
          current_chi2ndf_x.push_back(vdma->GetChi2Ndf(*it, "x"));
          current_chi2ndf_y.push_back(vdma->GetChi2Ndf(*it, "y"));

          current_sigma_x_err.push_back(vdma->GetVdmParameterError(*it, "Sigma_x"));
          current_sigma_y_err.push_back(vdma->GetVdmParameterError(*it, "Sigma_y"));
          current_mu_max_x_err.push_back(vdma->GetVdmParameterError(*it, "mu_max_x"));
          current_mu_max_y_err.push_back(vdma->GetVdmParameterError(*it, "mu_max_y"));
          current_c_x_err.push_back(vdma->GetVdmParameterError(*it, "c_x"));
          current_c_y_err.push_back(vdma->GetVdmParameterError(*it, "c_y"));
        }
        t_vdm->Fill();

        f_out->cd();
        TGraphErrors *tg_x = (TGraphErrors*)vdma->GetTGraphX()->Clone();
        TString name = "tg_x_BCID";
        name += *bcid;
        name += "_scan";
        name += *scan;
        name += "_NTrkCut";
        name += *nTrkCut;
        tg_x->SetName(name);
        tg_x->Write();
        TGraphErrors *tg_y = (TGraphErrors*)vdma->GetTGraphY()->Clone();
        name = "tg_y_BCID";
        name += *bcid;
        name += "_scan";
        name += *scan;
        name += "_NTrkCut";
        name += *nTrkCut;
        tg_y->SetName(name);
        tg_y->Write();

        delete vdma;

      }
    }
  }
  f_out->cd();
  t_vdm->Write();
  f_out->Close();


}


void usage() {
  cerr << "runVdM -r run -s settings [-e | -u] [-n] [-u systName] [-v]" << endl;
  cerr << "r: Specify run number (178013, 178064, 183013, 201351)" << endl;
  cerr << "s: Specify tracking settings (16.X-normal, 17.2-normal, 17.2-VtxLumi)" << endl;
  cerr << "n: Delete previous histograms" << endl;
  cerr << "u: Enable given systematic uncertainty" << endl;
  cerr << "v: Write out debug histograms" << endl;
  cerr << endl << "By default VtxC analysis is run, otherwise:" << endl;
  cerr << "e: Doing event-based analysis EvC (default: VtxC)" << endl;
  cerr << "f: Doing unfolding-based analysis UnfC (default: VtxC)" << endl;
}
