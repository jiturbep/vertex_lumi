/**TEST
  *  Interface for class LumiVtx to calculate vertex-based luminosities.
  *  Author: David R. Yu
  *
  ** Input: -r <run number> -s <"run" or "scan"> -o <output prefix> [-v]
  *
  ** Output: root file one TTrees, with pLB : <lumi from all algs>.
  *
  *  Example: ./GetLuminosity.exe -r 182013 -s scan
*/

//#define DEBUG
#define NOSPLIT

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

#include "GlobalSettings/GlobalSettings.h"
#include "LumiVtx/LumiVtx.h"
#include "VdM/VtxCalibration.h"

using namespace std;



// --- Usage
void usage(const char *pname);

// --- Helpers
bool fileExists(TString path);

int main(int argc, char **argv) {

#ifdef INCLUDE_HACKS
  cout << "[GetLuminosity] WARNING : Flag INCLUDE_HACKS defined. Do not use final values for anything serious!!" << endl;
#endif

  GlobalSettings gs;

  std::vector<Int_t> known_runs;
  known_runs.push_back(182013);
  known_runs.push_back(188949);
  known_runs.push_back(188951);
  known_runs.push_back(191373);
  known_runs.push_back(200805);
  known_runs.push_back(216399);
  known_runs.push_back(216416);
  known_runs.push_back(216432);

  //Input/output detection
  TString run;
  TString settings;
  Int_t single_NTrkCut = -1;
  Int_t single_BCID = -1;
  TString output_prefix = "";
  TString output_path = gs.path_lumiVtx;//in alternative to the full prefix
  bool set_verbose = false;
  bool redo_masking_correction = false;
  Float_t mu_scale = -1.;

  std::map<TString, bool> systematic_uncertainty_list;
  systematic_uncertainty_list["fake_low"] = false;
  systematic_uncertainty_list["fake_high"] = false;
  systematic_uncertainty_list["masking_toy_scaling"] = false;
  systematic_uncertainty_list["bs-45mm"] = false;
  systematic_uncertainty_list["bs-55mm"] = false;

  std::vector<TString> vtx_settings;
  vtx_settings.push_back("NVtx");
//  vtx_settings.push_back("NEvt");
//  vtx_settings.push_back("NUnf");
  
  std::vector<Int_t> nTrkCuts;
 	//nTrkCuts.push_back(2);
  //nTrkCuts.push_back(3);
  //nTrkCuts.push_back(4);
  nTrkCuts.push_back(5);
  //nTrkCuts.push_back(6);
  //nTrkCuts.push_back(7);
  //nTrkCuts.push_back(8);
  //nTrkCuts.push_back(10);

  // --- Scan for command line parameters
  int c;
  extern int optind;
  extern char* optarg;
  while (1) {
    int option_index = 0;
    static struct option long_options[] = {
      {"run",  required_argument,  0,  'r'},
      {"settings", required_argument, 0, 's'},
      {"outPrefix",  required_argument,  0,  'o'},
      {"help",  no_argument,    0,  'h'},
      {"verbose",  no_argument,    0,  'v'},
      {0, 0, 0, 0}
    };
    c=getopt_long(argc, argv, "mr:s:o:vhf:u:n:i:p:",long_options,&option_index);
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
      case 'o': {
        output_prefix = optarg;
        break;
      }
      case 'p': {
        output_path = optarg;
        break;
      }
      case 'v': {
        set_verbose = true;
        break;
      }
      case 'h': {
        usage((const char*)argv[0]);
        return 0;
        break;
      }
      case 'm': {
        cout << "[GetLuminosity] INFO : Forcing remake of pileup correction." << endl;
        redo_masking_correction = true;
        break;
      }
      case 'f': {
        cout << "[GetLuminosity] INFO : Setting fake correction mu scale to " << optarg << endl;
        stringstream ss;
        ss.str(optarg);
        try {
          ss >> mu_scale;
        } catch (const char *) {
          cerr << "[GetLuminosity] ERROR : Invalid mu scale specified : " << optarg << endl;
          return(1);
          break;
        }
        break;
      }
      case 'u': {
        if (systematic_uncertainty_list.find(optarg) == systematic_uncertainty_list.end()) {
          cerr << "[GetLuminosity] ERROR : Unrecognized systematic uncertainty test specified: " << optarg << ". Exiting..." << endl;
          exit(1);
        }
        systematic_uncertainty_list[optarg] = true;
        break;
      }
      case 'n': {
        cout << "[GetLuminosity] INFO : Running on a single NTrkCut: " << optarg << endl;
        stringstream ss;
        ss.str(optarg);
        try {
          ss >> single_NTrkCut;
        } catch (const char *) {
          cerr << "[GetLuminosity] ERROR : Invalid single NTrkCut specified: " << optarg << endl;
          return(1);
          break;
        }
        break;
      }
      case 'i': {
        cout << "[GetLuminosity] INFO : Running on a single BCID: " << optarg << endl;
        stringstream ss;
        ss.str(optarg);
        try {
          ss >> single_BCID;;
        } catch (const char *) {
          cerr << "[GetLuminosity] ERROR : Invalid single BCID specified: " << optarg << endl;
          exit(1);
          break;
        }
        break;
      }
    }
  } // End input parameter scanning.

  //Check parameters
  if (run == "") {
    cerr << "Input error: must specify a run number (or all)." << endl;
    usage((const char*)argv[0]);
    exit(1);
  }

  if (settings == "") {
    cerr << "Input error: must specify a setting" << endl;
    usage((const char*)argv[0]);
    exit(1);
  }

  if (output_prefix == "") {
    stringstream ss_folder;
    ss_folder << output_path << gs.getPrefix( run ) << run << "/" << settings;
    //ss_folder << output_path << gs.getPrefix( run ) << run ;
    for (map<TString, bool>::iterator it = systematic_uncertainty_list.begin(); it != systematic_uncertainty_list.end(); ++it) {
      if ((*it).second) {
        ss_folder << "_" << (*it).first;
      }
    }
    //boost::filesystem::create_directories(ss_folder.str());
    output_prefix = ss_folder.str() ; output_prefix += "/lumi_r" ; output_prefix += run;
    if (single_NTrkCut > 0) {
      output_prefix += "_NTrk" ; output_prefix += single_NTrkCut ;
    }
    if (single_BCID > 0) {
      output_prefix += "_BCID" ; output_prefix += single_BCID;
    }
  }
  output_prefix += "WitCal_2newMCs";
  //output_prefix += "_JulyCal_2newMCs";
  //output_prefix += "_NovCal_2newMCs";

  if (!settings) {
    cerr << "Input error: must specify -s settings." << endl;
    usage((const char*)argv[0]);
    exit(1);
  }

  //Done checking input

  /**********************************************************************/
  // Complete list of settings and input files
  /**********************************************************************/


  if (single_NTrkCut > 0) {
    if (find(nTrkCuts.begin(), nTrkCuts.end(), single_NTrkCut) == nTrkCuts.end()) {
      cerr << "[GetLuminosity] ERROR : If specifying a single NTrkCut, it must be 5, 7 or 10. Exiting..." << endl;
      return(1);
    } else {
      nTrkCuts.clear();
      nTrkCuts.push_back(single_NTrkCut);
    }
  }

  //Inputs from D3PD processing
  stringstream ss_path_vertex_counts;
  ss_path_vertex_counts << gs.path_inputRawCount << "/"
                        << gs.getPrefix( run ) << run << "/" << settings << "/"
                        << gs.path_inputRawCount_v << gs.fname_inputRawCount_Tree;
                        cout<<"[GetLuminosity] Line 257, ss_path_vertex_counts = "<< ss_path_vertex_counts.str() <<endl;
  TString path_vertex_counts = ss_path_vertex_counts.str();
  stringstream ss_path_vertex_histograms;
  ss_path_vertex_histograms << gs.path_inputRawCount << "/"
                            << gs.getPrefix( run ) << run << "/" << settings << "/"
                            << gs.path_inputRawCount_v << gs.fname_inputRawCount_Histo;
                            cout<<"[GetLuminosity] Line 263, ss_path_vertex_histograms = "<< ss_path_vertex_histograms.str() <<endl;
  TString path_vertex_histograms = ss_path_vertex_histograms.str();
  /*
    if (run == "188951") {
      //folder = "~/Luminosity/LumiVtx/D3PD/results/188951/v7/";
      //folder = "~/Luminosity/LumiVtx/D3PD_muscan/results/188951/blzphi/";
      folder = "/eliza18/atlas/dryu/Luminosity/VertexCounts/188951/17.2-normal/";
    } else if (run == "182013") {
      //folder = "~/Luminosity/vdm/D3PD/results/182013/"; folder += subrun; folder += "/v6/";
      folder = "~/Luminosity/LumiVtx/D3PD/results/182013/v7/";
    } else if (run == "188949") {
      folder = "~/Luminosity/vdm/D3PD/results/188949/scan/v6/";
    } else if (run == "191373") {
      folder = "~/Luminosity/LumiVtx/D3PD/results/191373/hacktrigger/";
    } else if (run == "200805") {
      folder = "/eliza18/atlas/dryu/Luminosity/VertexCounts/200805/17.2-normal/v9/";
    }
  */

  // Other stuff comes from different place :( So, we have to do some special detection.
  vector<TString> timestamp_paths, duration_paths, lumi_ntuple_paths;
  vector<Int_t> bcidList;
  TString livefraction_path(""), plb_skiplist_path("");
  vector<Int_t> plb_skiplist;

  bool do_physics_run = false;

  if (run == "182013") {
    timestamp_paths.push_back(gs.path_timestamps + "/182013/run/all_timestamps.dat");
    bcidList.push_back(81);
    bcidList.push_back(867);
    bcidList.push_back(2752);

    livefraction_path = gs.path_deadtime + "/dt_182013.root";
    plb_skiplist_path = gs.path_plb_skiplist + "/plb_skiplist_182013.txt";
  } else if (run == "191373") {
    timestamp_paths.push_back(gs.path_timestamps + "/191373/run/191373_timestamps.dat");
    bcidList.push_back(1);
    livefraction_path = gs.path_deadtime + "/dt_191373.root";
  } else if (run == "188949") {
    timestamp_paths.push_back(gs.path_timestamps + "/188949/scan/timestamps1.dat");
    timestamp_paths.push_back(gs.path_timestamps + "/188949/scan/timestamps2.dat");

    bcidList.push_back(200);
    bcidList.push_back(999);

    plb_skiplist_path = gs.path_plb_skiplist + "/plb_skiplist_188949.txt";
  } else if (run == "188951") {
    timestamp_paths.push_back(gs.path_timestamps + "/188951/scan/timestamps3.dat");

    bcidList.push_back(200);
    bcidList.push_back(999);

    plb_skiplist_path = gs.path_plb_skiplist + "/plb_skiplist_188951.txt";

  } else if (run == "200805") {
    duration_paths.push_back(gs.path_timestamps + "/200805/200805_durations.dat");

    bcidList.push_back(1);
    bcidList.push_back(149);
    bcidList.push_back(334);
  } else if (run == "201351") {
    timestamp_paths.push_back(gs.path_timestamps + "/201351/run/all_timestamps.dat");

    bcidList.push_back(1);
    bcidList.push_back(241);
    bcidList.push_back(2881);
    bcidList.push_back(3121);
  /*} else if (run == "216399") {
    stringstream ss;
    ss << gs.path_PhysRun_lumiNtuples + "/r" << run << ".root";
    lumi_ntuple_paths.push_back(TString(ss.str()));
    bcidList.push_back(1);
    for (int i = 301; i<=315 ; i++){
      bcidList.push_back(i);
    }
    for (int j = 1195; j<=1209 ; j++){
      bcidList.push_back(j);
    } 
  } else if (run == "216416") {
    stringstream ss;
    ss << gs.path_PhysRun_lumiNtuples + "/r" << run << ".root";
    lumi_ntuple_paths.push_back(TString(ss.str()));
    for (int i = 301; i<=315 ; i++){
      bcidList.push_back(i);
    }
    for (int j = 1195; j<=1209 ; j++){
      bcidList.push_back(j);
    }
    for (int k = 2089; k<=2103 ; k++){
      bcidList.push_back(k);
    }
    for (int l = 2977; l<=2991 ; l++){
      bcidList.push_back(l);
    }
  } else if (run == "216432") {
    stringstream ss;
    ss << gs.path_PhysRun_lumiNtuples + "/r" << run << ".root";
    lumi_ntuple_paths.push_back(TString(ss.str()));
    for (int i = 301; i<=315 ; i++){
      bcidList.push_back(i);
    }
    for (int j = 1195; j<=1209 ; j++){
      bcidList.push_back(j);
    }
    for (int k = 2089; k<=2103 ; k++){
      bcidList.push_back(k);
    }
    for (int l = 2977; l<=2991 ; l++){
      bcidList.push_back(l);
    }*/
  } else {
    // Physics run
    do_physics_run = true;
    stringstream ss;
    ss << gs.path_PhysRun_lumiNtuples + "/r" << run << ".root";
    lumi_ntuple_paths.push_back(TString(ss.str()));
    bcidList.push_back(0);
  }

  // Check single BCID specification
  if (single_BCID > 0) {
    if (find(bcidList.begin(), bcidList.end(), single_BCID) == bcidList.end()) {
      cerr << "[GetLuminosity] ERROR : Specified single BCID " << single_BCID
           << " does not belong to run " << run << ". Exiting..." << endl;
      return(1);
    } else {
      bcidList.clear();
      bcidList.push_back(single_BCID);
    }
  }

  //pLBs to skip because of moving beam
  if (plb_skiplist_path != "") {
    ifstream f_plb_skiplist(plb_skiplist_path.Data());
    if (not f_plb_skiplist.is_open()) {
      cerr << "Unable to open pLB skip list: " << plb_skiplist_path << endl;
      exit(1);
    }
    Int_t plb_to_skip;
    while (f_plb_skiplist.good()) {
      string line;
      getline(f_plb_skiplist, line);
      if (line.empty()) {
        continue;
      }
      stringstream istr(line);
      istr >> plb_to_skip;
      plb_skiplist.push_back(plb_to_skip);
    }
    sort(plb_skiplist.begin(), plb_skiplist.end());
  }


  //Energies
  TString energy;
  if (atoi(run) < 200000) {
    energy = "7";
  } else {
    energy = "8";
//    energy = "7";
  }

  /**********************************************************************/
  // Output TTree
  /**********************************************************************/

  //Declare TTrees here.
  // - Per-LB values
  // -- Containers for filling

  Int_t current_run;
  Int_t current_bcid;
  Int_t current_lb;
  Float_t current_lb_duration;
////////////////////////////////////////////
/*CHANGED THE SIZES OF VECTORS FROM 4 TO 6*/
  Float_t current_lb_inst_luminosity_nvtx[6];
  Float_t current_lb_inst_luminosity_nvtx_err[6];

  Float_t current_lb_inst_luminosity_nevt[6];
  Float_t current_lb_inst_luminosity_nevt_errup[6];
  Float_t current_lb_inst_luminosity_nevt_errdown[6];

  Float_t current_lb_inst_luminosity_nunf[6];
  Float_t current_lb_inst_luminosity_nunf_err[6];
  Float_t current_lb_total_luminosity_nunf[6];
  Float_t current_lb_total_luminosity_nunf_err[6];

  Float_t current_lb_nvtx[6];
  Float_t current_lb_nvtx_masked[6];
  Float_t current_lb_nevt[6];
  Float_t current_lb_nvtx_from_nevt[6];

  Float_t current_lb_masking_correction[6];
  Float_t current_lb_fake_mu_nvtx[6];
  Float_t current_lb_fake_mu_nevt[6];
///////////////////////////////////////////////

  Float_t current_lb_inst_luminosity_bcmhor, current_lb_inst_luminosity_bcmhand, current_lb_inst_luminosity_bcmvor, current_lb_inst_luminosity_bcmvand, current_lb_inst_luminosity_lucidor, current_lb_inst_luminosity_lucidand, current_lb_inst_luminosity_lucidhitor, current_lb_inst_luminosity_lucidhitand;
  Float_t current_lb_inst_luminosity_bcmhor_err, current_lb_inst_luminosity_bcmhand_err, current_lb_inst_luminosity_bcmvor_err, current_lb_inst_luminosity_bcmvand_err, current_lb_inst_luminosity_lucidor_err, current_lb_inst_luminosity_lucidand_err, current_lb_inst_luminosity_lucidhitor_err, current_lb_inst_luminosity_lucidhitand_err;

  std::map<TString, Float_t> sigma_vis;
  std::map<TString, Float_t> sigma_vis_err;

  TTree *t_lb = new TTree("LumiVtx", "Vertex-based luminosity vs. pLB");
  TBranch *b_run2 = t_lb->Branch("Run", &current_run, "Run/I");
  TBranch *b_bcid = t_lb->Branch("BCID", &current_bcid, "BCID/I");
  TBranch *b_lb = t_lb->Branch("LB", &current_lb, "LB/I");
  TBranch *b_lb_duration = t_lb->Branch("Duration", &current_lb_duration, "Duration/F");

////////////////////////////////////////////
/*CHANGED THE SIZES OF VECTORS FROM 4 TO 6*/
  TBranch *b_inst_luminosity_nvtx[6];
  TBranch *b_inst_luminosity_nvtx_err[6];
  TBranch *b_masking_correction[6];
  TBranch *b_fake_mu_nvtx[6];
  TBranch *b_fake_mu_nevt[6];
  TBranch *b_inst_luminosity_nunf[6];
  TBranch *b_inst_luminosity_nunf_err[6];
  TBranch *b_total_luminosity_nunf[6];
  TBranch *b_total_luminosity_nunf_err[6];
////////////////////////////////////////////

  for (int i=0; i<nTrkCuts.size(); i++) {
    TString name = "Linst_NVtx"; name += nTrkCuts.at(i);
    b_inst_luminosity_nvtx[i] = t_lb->Branch( name, &(current_lb_inst_luminosity_nvtx[i]), (name+"/F") );
    b_inst_luminosity_nvtx_err[i] = t_lb->Branch( (name+"_err"), &(current_lb_inst_luminosity_nvtx_err[i]), (name+"_err/F") );

    name = "MaskingCorrection_NVtx"; name += nTrkCuts.at(i);
    b_masking_correction[i] = t_lb->Branch( name, &(current_lb_masking_correction[i]), (name+"/F") );

    name = "FakeMu_NVtx"; name += nTrkCuts.at(i);
    b_fake_mu_nvtx[i] = t_lb->Branch( name, &(current_lb_fake_mu_nvtx[i]), (name+"/F") );

    name = "FakeMu_NEvt"; name += nTrkCuts.at(i);
    b_fake_mu_nevt[i] = t_lb->Branch( name, &(current_lb_fake_mu_nevt[i]), (name+"/F") );

    name = "Linst_NUnf"; name += nTrkCuts.at(i);
    b_inst_luminosity_nunf[i] = t_lb->Branch( name, &(current_lb_inst_luminosity_nunf[i]), (name+"/F") );
    b_inst_luminosity_nunf[i] = t_lb->Branch( (name+"_err"), &(current_lb_inst_luminosity_nunf_err[i]), (name+"_err/F") );

    name = "Ltotal_NUnf"; name += nTrkCuts.at(i);
    b_total_luminosity_nunf[i] = t_lb->Branch( name, &(current_lb_total_luminosity_nunf[i]), (name+"/F") );
    b_total_luminosity_nunf[i] = t_lb->Branch( (name+"_err"), &(current_lb_total_luminosity_nunf_err[i]), (name+"_err/F") );

  }

  TBranch *b_inst_luminosity_nevt[4];
  TBranch *b_inst_luminosity_nevt_errup[4];
  TBranch *b_inst_luminosity_nevt_errdown[4];

  for (int i=0; i<nTrkCuts.size(); i++) {
    TString name = "Linst_NEvt"; name += nTrkCuts.at(i);
    TString leafname = name;
    leafname += "/F";
    b_inst_luminosity_nevt[i] = t_lb->Branch(name, &(current_lb_inst_luminosity_nevt[i]), leafname);

    name = "Linst_NEvt"; name += nTrkCuts.at(i);
    name += "_errup";
    leafname = name;
    leafname += "/F";
    b_inst_luminosity_nevt_errup[i] = t_lb->Branch(name, &(current_lb_inst_luminosity_nevt_errup[i]), leafname);


    name = "Linst_NEvt"; name += nTrkCuts.at(i);
    name += "_errdown";
    leafname = name;
    leafname += "/F";
    b_inst_luminosity_nevt_errdown[i] = t_lb->Branch(name, &(current_lb_inst_luminosity_nevt_errdown[i]), leafname);
  }

  TBranch *b_inst_luminosity_bcmhor = t_lb->Branch("Linst_bcmhor", &current_lb_inst_luminosity_bcmhor, "Linst_bcmhor/F");
  TBranch *b_inst_luminosity_bcmhand = t_lb->Branch("Linst_bcmhand", &current_lb_inst_luminosity_bcmhand, "Linst_bcmhand/F");
  TBranch *b_inst_luminosity_bcmvor = t_lb->Branch("Linst_bcmvor", &current_lb_inst_luminosity_bcmvor, "Linst_bcmvor/F");
  TBranch *b_inst_luminosity_bcmvand = t_lb->Branch("Linst_bcmvand", &current_lb_inst_luminosity_bcmvand, "Linst_bcmvand/F");
  TBranch *b_inst_luminosity_lucidor = t_lb->Branch("Linst_lucidor", &current_lb_inst_luminosity_lucidor, "Linst_lucidor/F");
  TBranch *b_inst_luminosity_lucidand = t_lb->Branch("Linst_lucidand", &current_lb_inst_luminosity_lucidand, "Linst_lucidand/F");
  TBranch *b_inst_luminosity_lucidhitor = t_lb->Branch("Linst_lucidhitor", &current_lb_inst_luminosity_lucidhitor, "Linst_lucidhitor/F");
  TBranch *b_inst_luminosity_lucidhitand = t_lb->Branch("Linst_lucidhitand", &current_lb_inst_luminosity_lucidhitand, "Linst_lucidhitand/F");

  TBranch *b_inst_luminosity_bcmhor_err = t_lb->Branch("Linst_bcmhor_err", &current_lb_inst_luminosity_bcmhor_err, "Linst_bcmhor_err/F");
  TBranch *b_inst_luminosity_bcmhand_err = t_lb->Branch("Linst_bcmhand_err", &current_lb_inst_luminosity_bcmhand_err, "Linst_bcmhand_err/F");
  TBranch *b_inst_luminosity_bcmvor_err = t_lb->Branch("Linst_bcmvor_err", &current_lb_inst_luminosity_bcmvor_err, "Linst_bcmvor_err/F");
  TBranch *b_inst_luminosity_bcmvand_err = t_lb->Branch("Linst_bcmvand_err", &current_lb_inst_luminosity_bcmvand_err, "Linst_bcmvand_err/F");
  TBranch *b_inst_luminosity_lucidor_err = t_lb->Branch("Linst_lucidor_err", &current_lb_inst_luminosity_lucidor_err, "Linst_lucidor_err/F");
  TBranch *b_inst_luminosity_lucidand_err = t_lb->Branch("Linst_lucidand_err", &current_lb_inst_luminosity_lucidand_err, "Linst_lucidand_err/F");
  TBranch *b_inst_luminosity_lucidhitor_err = t_lb->Branch("Linst_lucidhitor_err", &current_lb_inst_luminosity_lucidhitor_err, "Linst_lucidhitor_err/F");
  TBranch *b_inst_luminosity_lucidhitand_err = t_lb->Branch("Linst_lucidhitand_err", &current_lb_inst_luminosity_lucidhitand_err, "Linst_lucidhitand_err/F");


////////////////////////////////////////////
/*CHANGED THE SIZES OF VECTORS FROM 4 TO 6*/
  TBranch *b_nvtx[6];
  TBranch *b_nvtx_masked[6];
  TBranch *b_nevt[6];
  TBranch *b_nvtx_from_nevt[6];
  for (int i=0; i<nTrkCuts.size(); i++) {
    TString name = "NVtx_NVtx";
    name += (Int_t)(*(nTrkCuts.begin() + i));
    TString leafname = name;
    leafname += "/F";
    b_nvtx[i] = t_lb->Branch(name, &(current_lb_nvtx[i]), leafname);

    name = "NVtxMasked_NVtx";
    name += (Int_t)(*(nTrkCuts.begin() + i));
    leafname = name;
    leafname += "/F";
    b_nvtx_masked[i] = t_lb->Branch(name, &(current_lb_nvtx_masked[i]), leafname);

    name = "NEvt_NEvt";
    name += (Int_t)(*(nTrkCuts.begin() + i));
    leafname = name;
    leafname += "/F";
    b_nevt[i] = t_lb->Branch(name, &(current_lb_nevt[i]), leafname);

    name = "NVtx_NEvt";
    name += (Int_t)(*(nTrkCuts.begin() + i));
    leafname = name;
    leafname += "/F";
    b_nvtx_from_nevt[i] = t_lb->Branch(name, &(current_lb_nvtx_from_nevt[i]), leafname);
  }

  TBranch *b_sigma_vis_nvtx[6];
  TBranch *b_sigma_vis_nevt[6];
  TBranch *b_sigma_vis_nunf[6];
  for (int i=0; i<nTrkCuts.size(); i++) {
    TString name, leafname, algname;
    name = "SigmaVis_NVtx";
    name += (Int_t)(nTrkCuts[i]);
    leafname = name;
    leafname += "/F";
    algname = "NVtx";
    algname += (Int_t)(nTrkCuts[i]);
    sigma_vis[algname] = 0.0;
    sigma_vis_err[algname] = 0.0;
    b_sigma_vis_nvtx[i] = t_lb->Branch(name, &(sigma_vis[algname]), leafname);

    name = "SigmaVis_NEvt";
    name += (Int_t)(nTrkCuts[i]);
    leafname = name;
    leafname += "/F";
    algname = "NEvt";
    algname += (Int_t)(nTrkCuts[i]);
    sigma_vis[algname] = 0.0;
    sigma_vis_err[algname] = 0.0;
    b_sigma_vis_nevt[i] = t_lb->Branch(name, &(sigma_vis[algname]), leafname);

    name = "SigmaVis_NUnf";
    name += (Int_t)(nTrkCuts[i]);
    leafname = name;
    leafname += "/F";
    algname = "NUnf";
    algname += (Int_t)(nTrkCuts[i]);
    sigma_vis[algname] = 0.0;
    sigma_vis_err[algname] = 0.0;
    b_sigma_vis_nunf[i] = t_lb->Branch(name, &(sigma_vis[algname]), leafname);
  }

  TBranch *b_sigma_vis_bcmhor = t_lb->Branch("SigmaVis_bcmhor",   &(sigma_vis["bcmhor"]), "SigmaVis_bcmhor/F");
  TBranch *b_sigma_vis_bcmhand = t_lb->Branch("SigmaVis_bcmhand",  &(sigma_vis["bcmhand"]), "SigmaVis_bcmhand/F");
  TBranch *b_sigma_vis_bcmvor = t_lb->Branch("SigmaVis_bcmvor",   &(sigma_vis["bcmvor"]), "SigmaVis_bcmvor/F");
  TBranch *b_sigma_vis_bcmvand = t_lb->Branch("SigmaVis_bcmvand",  &(sigma_vis["bcmvand"]), "SigmaVis_bcmvand/F");
  TBranch *b_sigma_vis_lucidor = t_lb->Branch("SigmaVis_lucidor",  &(sigma_vis["lucidor"]), "SigmaVis_lucidor/F");
  TBranch *b_sigma_vis_lucidand = t_lb->Branch("SigmaVis_lucidand",  &(sigma_vis["lucidand"]), "SigmaVis_lucidand/F");
  TBranch *b_sigma_vis_lucidhitor = t_lb->Branch("SigmaVis_lucidhitor",  &(sigma_vis["lucidhitor"]), "SigmaVis_lucidhitor/F");
  TBranch *b_sigma_vis_lucidhitand = t_lb->Branch("SigmaVis_lucidhitand", &(sigma_vis["lucidhitand"]), "SigmaVis_lucidhitand/F");

  /**********************************************************************/
  // Luminosity from other algorithms
  /**********************************************************************/
  std::map<Int_t, std::map<Int_t, Float_t> > lumi_bcmhor, lumi_bcmhand, lumi_bcmvor, lumi_bcmvand, lumi_lucidor, lumi_lucidand, lumi_lucidhitor, lumi_lucidhitand; // lumi_bcm is from Stefan, all others from Benedetto. Hopefully lumi_bcm == lumi_bcmhor!
  std::map<Int_t, std::map<Int_t, Float_t> > lumi_bcmhor_err, lumi_bcmhand_err, lumi_bcmvor_err, lumi_bcmvand_err, lumi_lucidor_err, lumi_lucidand_err, lumi_lucidhitor_err, lumi_lucidhitand_err; // lumi_bcm is from Stefan_err, all others from Benedetto. Hopefully lumi_bcm == lumi_bcmhor!

  if (atoi(run) < 200000) {
    sigma_vis["bcmhor"] = 4689;
    sigma_vis["bcmvor"] = 4736;
    sigma_vis["lucidor"] = 42510;
    sigma_vis["lucidand"] = 13440;
    sigma_vis["bcmhand"] = 135.9;
    sigma_vis["bcmvand"] = 138.7;
    sigma_vis["lucidhitor"] = 42510 / 26;
    sigma_vis["lucidhitand"] = 13440 / 26;
  } else {
  //From OLCOnlineConfiguration 6 June 2012
    /*
  //Values I calculated 01/August/2013
    sigma_vis["bcmhor"] = 4974.5; //4860.5OLCO
    sigma_vis["bcmvor"] = 4934;
    sigma_vis["lucidor"] = 34908.6; //36575;OLCO
    sigma_vis["lucidand"] = 9550; //OLCO
    sigma_vis["bcmhand"] = 149.95; //OLCO
    sigma_vis["bcmvand"] = 148.7; //OLCO
    sigma_vis["lucidhitor"] = 2628.6; //3286;OLCO  
    sigma_vis["lucidhitand"] = 13440 / 26; //2011
    */
  //Values I calculated and averaged over July scans
    /*sigma_vis["bcmhor"] = 4868.821097; 
    sigma_vis["bcmvor"] = 4888.050745;
    sigma_vis["lucidor"] = 34570.42712; 
    sigma_vis["lucidand"] = 9550; //OLCO
    sigma_vis["bcmhand"] = 149.95; //OLCO
    sigma_vis["bcmvand"] = 148.7; //OLCO
    sigma_vis["lucidhitor"] = 2646.376093; 
    sigma_vis["lucidhitand"] = 13440 / 26; *///2011
    
 //Values from Witold's presentation of the 7 of November 2013
    sigma_vis["bcmhor"] = 4976.5;
    sigma_vis["bcmvor"] = 4937.3;
    sigma_vis["lucidor"] = 34709.0;
    sigma_vis["lucidhitor"] = 2609.5;

 
//    exit(1);
  }

  //Load Stefan's hard-to-process text file
  if (run == "182013") {
    TString bcm_filename = gs.path_lumiTxt + "/" + gs.fname_lumiBcmH;
    ifstream mu_bcmh_file(bcm_filename.Data());
    if (not mu_bcmh_file.is_open()) {
      cerr << "Unable to open Stefan's mu file: " << bcm_filename.Data() << endl;
      exit(1);
    }

    TString next_line = "separator";
    Int_t c_plb, c_bcid;
    Float_t c_mu_gp0b, c_mu_gp0L, c_mu_spline;
    TString blank;

    while (mu_bcmh_file.good()) {
      string line;
      getline(mu_bcmh_file, line);
      TString ts_line = line;
      if (ts_line.Contains("********") || next_line == "separator") {
        next_line = "pLB_info";
        continue;
      } else if (next_line == "pLB_info") {
        stringstream istr(line);
        istr >> c_plb >> blank;
        next_line = "skip1";
        continue;
      } else if (next_line == "skip1") {
        next_line = "skip2";
        continue;
      } else if (next_line == "skip2") {
        next_line = "skip3";
        continue;
      } else if (next_line == "skip3") {
        next_line = "skip4";
      } else if (next_line == "skip4") {
        next_line = "mu_line";
      } else if (next_line == "mu_line") {
        stringstream istr(line);
        istr >> c_bcid >> c_mu_gp0b >> c_mu_gp0L >> c_mu_spline;
        lumi_bcmhor[c_bcid][c_plb] = c_mu_gp0L * 11245.5 / 71500. * 4.6915 / 4.711;
        if (c_bcid == 2752) {
          next_line = "separator";
        } else {
          next_line = "mu_line";
        }
        continue;
      }
    }

    bcm_filename = gs.path_lumiTxt + "/" + gs.fname_lumiBcmV;
    ifstream mu_bcmv_file(bcm_filename.Data());
    if (not mu_bcmv_file.is_open()) {
      cerr << "Unable to open Stefan's mu file: " << bcm_filename.Data() << endl;
      exit(1);
    }

    next_line = "separator";

    while (mu_bcmv_file.good()) {
      string line;
      getline(mu_bcmv_file, line);
      TString ts_line = line;
      if (ts_line.Contains("********") || next_line == "separator") {
        next_line = "pLB_info";
        continue;
      } else if (next_line == "pLB_info") {
        stringstream istr(line);
        istr >> c_plb >> blank;
        next_line = "skip1";
        continue;
      } else if (next_line == "skip1") {
        next_line = "skip2";
        continue;
      } else if (next_line == "skip2") {
        next_line = "skip3";
        continue;
      } else if (next_line == "skip3") {
        next_line = "skip4";
      } else if (next_line == "skip4") {
        next_line = "mu_line";
      } else if (next_line == "mu_line") {
        stringstream istr(line);
        istr >> c_bcid >> c_mu_gp0b >> c_mu_gp0L >> c_mu_spline;
        lumi_bcmvor[c_bcid][c_plb] = c_mu_gp0L * 11245.5 / 71500. * 4.7365 / 4.757;
        if (c_bcid == 2752) {
          next_line = "separator";
        } else {
          next_line = "mu_line";
        }
        continue;
      }
    }

  } //end of if run = 182013

  //Load Benedetto's text file
  if (run == "188951") {
    TString allmu_filename = gs.path_lumiTxt + "/" + gs.fname_muScanAll;
    ifstream mu_file(allmu_filename.Data());
    if (not mu_file.is_open()) {
      cerr << "Unable to open Benedetto's mu file: " << allmu_filename.Data() << endl;
      exit(1);
    }

    Int_t c_plb, c_bcid;
    Float_t c_bcmhor, c_bcmhand, c_bcmvor, c_bcmvand, c_lucidor, c_lucidand;

    while (mu_file.good()) {
      string line;
      getline(mu_file, line);
      if (line.empty()) {
        continue;
      }
      stringstream istr(line);
      istr >> c_plb >> c_bcid >> c_bcmhor >> c_bcmhand >> c_bcmvor >> c_bcmvand >> c_lucidor >> c_lucidand;
      //cout << "values from Benedetto's file:  " << c_bcmhor << " " << c_bcmhand << endl;
	    //cout << "plb:  " << c_plb << ", bcid " << c_bcid << endl;
	    if ((c_bcid != 200) && (c_bcid != 999)) {
        continue;
      }
      lumi_bcmhor[c_bcid][c_plb] = c_bcmhor * 11245.5 / 71500. * 4.6915 / 4.711;
      lumi_bcmhand[c_bcid][c_plb] = c_bcmhand * 11245.5 / 71500.;
      lumi_bcmvor[c_bcid][c_plb] = c_bcmvor * 11245.5 / 71500. * 4.7365 / 4.757;
      lumi_bcmvand[c_bcid][c_plb] = c_bcmvand * 11245.5 / 71500.;
      lumi_lucidor[c_bcid][c_plb] = c_lucidor * 11245.5 / 71500.;
      lumi_lucidand[c_bcid][c_plb] = c_lucidand * 11245.5 / 71500.;
    }
  } else if (run == "191373") {
    TString allmu_filename = gs.path_lumiTxt + "/" + gs.fname_r191373;
    ifstream mu_file(allmu_filename.Data());
    if (not mu_file.is_open()) {
      cerr << "Unable to open Benedetto's mu file: " << allmu_filename.Data() << endl;
      exit(1);
    }

    Int_t c_run, c_plb, c_bcid;
    Float_t c_bcmhor, c_bcmhand, c_bcmvor, c_bcmvand, c_lucidor, c_lucidand, c_lucidhitor, c_lucidhitand, c_lbdur, c_nbcids, c_xsec;
    Float_t c_bcmhor_err, c_bcmvor_err, c_lucidor_err, c_lucidhitor_err;
    while (mu_file.good()) {
      c_bcid = 1;
      string line;
      getline(mu_file, line);
      if (line.empty()) {
        continue;
      }
      stringstream istr(line);
      istr >> c_run >> c_plb >> c_bcmhor >> c_bcmhor_err >> c_bcmhand >> c_bcmvor >> c_bcmvor_err >> c_bcmvand >> c_lucidor >> c_lucidor_err >> c_lucidand >> c_lucidhitor >> c_lucidhitor_err >> c_lucidhitand >> c_lbdur >> c_nbcids >> c_xsec;
      if (c_bcid != 1) {
        continue;
      }
      lumi_bcmhor[c_bcid][c_plb] = c_bcmhor * 11245.5 / 71500.;
      lumi_bcmhor_err[c_bcid][c_plb] = c_bcmhor_err * 11245.5 / 71500.;
      lumi_bcmhand[c_bcid][c_plb] = c_bcmhand * 11245.5 / 71500.;
      lumi_bcmvor[c_bcid][c_plb] = c_bcmvor * 11245.5 / 71500.;
      lumi_bcmvor_err[c_bcid][c_plb] = c_bcmvor_err * 11245.5 / 71500.;
      lumi_bcmvand[c_bcid][c_plb] = c_bcmvand * 11245.5 / 71500.;
      lumi_lucidor[c_bcid][c_plb] = c_lucidor * 11245.5 / 71500.;
      lumi_lucidor_err[c_bcid][c_plb] = c_lucidor_err * 11245.5 / 71500.;
      lumi_lucidand[c_bcid][c_plb] = c_lucidand * 11245.5 / 71500.;
      lumi_lucidhitor[c_bcid][c_plb] = c_lucidhitor * 11245.5 / 71500.;
      lumi_lucidhitor_err[c_bcid][c_plb] = c_lucidhitor_err * 11245.5 / 71500.;
      lumi_lucidhitand[c_bcid][c_plb] = c_lucidhitand * 11245.5 / 71500.;

    }
  }

  // Load Eric's lumi ntuple
  std::map<Int_t, Int_t> n_bcids;
  Int_t n_bcids_run = 0;
  Int_t n_files = lumi_ntuple_paths.size();
  Int_t counter_files = 0;
  if (lumi_ntuple_paths.size() > 0) {
    for (vector<TString>::iterator itf = lumi_ntuple_paths.begin(); itf != lumi_ntuple_paths.end(); ++itf) {
      counter_files++;
      if (counter_files > n_files) {
        break;
      }

      TFile *f_other = new TFile(*itf, "READ");
      TTree *t_other = (TTree*)f_other->Get("lumiData");

      ULong64_t       LBDATA_StartTime;
      ULong64_t       LBDATA_EndTime;
      UInt_t          LBDATA_Run;
      UInt_t          LBDATA_LB;
      UInt_t          LBDATA_stable;
      Float_t         AvgBeam1;
      Float_t         AvgBeam2;
      Float_t         pmtA;
      Float_t         pmtC;
      Float_t         muToLumi;
      UInt_t          Status[3564];
      Float_t         RawLucOR[3564];
      Float_t         RawLucA[3564];
      Float_t         RawLucC[3564];
      Float_t         RawLucAND[3564];
      Float_t         RawLucHitOR[3564];
      Float_t         RawLucHitAND[3564];
      Float_t         RawBcmHOR[3564];
      Float_t         RawBcmHAND[3564];
      Float_t         RawBcmHXORA[3564];
      Float_t         RawBcmHXORC[3564];
      Float_t         RawBcmVOR[3564];
      Float_t         RawBcmVAND[3564];
      Float_t         RawBcmVXORA[3564];
      Float_t         RawBcmVXORC[3564];
      Float_t         RawZdcOR[3564];
      Float_t         RawZdcORA[3564];
      Float_t         RawZdcORC[3564];
      Float_t         RawZdcAND[3564];
      Float_t         MuLucOR[3564];
      Float_t         MuLucA[3564];
      Float_t         MuLucC[3564];
      Float_t         MuLucAND[3564];
      Float_t         MuLucHitOR[3564];
      Float_t         MuLucHitAND[3564];
      Float_t         MuBcmHOR[3564];
      Float_t         MuBcmHAND[3564];
      Float_t         MuBcmHXORA[3564];
      Float_t         MuBcmHXORC[3564];
      Float_t         MuBcmVOR[3564];
      Float_t         MuBcmVAND[3564];
      Float_t         MuBcmVXORA[3564];
      Float_t         MuBcmVXORC[3564];
      Float_t         MuZdcOR[3564];
      Float_t         MuZdcORA[3564];
      Float_t         MuZdcORC[3564];
      Float_t         MuZdcAND[3564];
      Float_t         Beam1[3564];
      Float_t         Beam2[3564];

      // List of branches
      TBranch        *b_LBDATA;   //!
      TBranch        *b_AvgBeam1;   //!
      TBranch        *b_AvgBeam2;   //!
      TBranch        *b_pmtA;   //!
      TBranch        *b_pmtC;   //!
      TBranch        *b_muToLumi;   //!
      TBranch        *b_Status;   //!
      TBranch        *b_RawLucOR;   //!
      TBranch        *b_RawLucA;   //!
      TBranch        *b_RawLucC;   //!
      TBranch        *b_RawLucAND;   //!
      TBranch        *b_RawLucHitOR;   //!
      TBranch        *b_RawLucHitAND;   //!
      TBranch        *b_RawBcmHOR;   //!
      TBranch        *b_RawBcmHAND;   //!
      TBranch        *b_RawBcmHXORA;   //!
      TBranch        *b_RawBcmHXORC;   //!
      TBranch        *b_RawBcmVOR;   //!
      TBranch        *b_RawBcmVAND;   //!
      TBranch        *b_RawBcmVXORA;   //!
      TBranch        *b_RawBcmVXORC;   //!
      TBranch        *b_RawZdcOR;   //!
      TBranch        *b_RawZdcORA;   //!
      TBranch        *b_RawZdcORC;   //!
      TBranch        *b_RawZdcAND;   //!
      TBranch        *b_MuLucOR;   //!
      TBranch        *b_MuLucA;   //!
      TBranch        *b_MuLucC;   //!
      TBranch        *b_MuLucAND;   //!
      TBranch        *b_MuLucHitOR;   //!
      TBranch        *b_MuLucHitAND;   //!
      TBranch        *b_MuBcmHOR;   //!
      TBranch        *b_MuBcmHAND;   //!
      TBranch        *b_MuBcmHXORA;   //!
      TBranch        *b_MuBcmHXORC;   //!
      TBranch        *b_MuBcmVOR;   //!
      TBranch        *b_MuBcmVAND;   //!
      TBranch        *b_MuBcmVXORA;   //!
      TBranch        *b_MuBcmVXORC;   //!
      TBranch        *b_MuZdcOR;   //!
      TBranch        *b_MuZdcORA;   //!
      TBranch        *b_MuZdcORC;   //!
      TBranch        *b_MuZdcAND;   //!
      TBranch        *b_Beam1;   //!
      TBranch        *b_Beam2;   //!

      t_other->SetMakeClass(1);

      t_other->SetBranchAddress("LBDATA", &LBDATA_StartTime, &b_LBDATA);
      t_other->SetBranchAddress("Status", Status, &b_Status);
      t_other->SetBranchAddress("RawLucOR", RawLucOR, &b_RawLucOR);
      t_other->SetBranchAddress("RawLucA", RawLucA, &b_RawLucA);
      t_other->SetBranchAddress("RawLucC", RawLucC, &b_RawLucC);
      t_other->SetBranchAddress("RawLucAND", RawLucAND, &b_RawLucAND);
      t_other->SetBranchAddress("RawLucHitOR", RawLucHitOR, &b_RawLucHitOR);
      t_other->SetBranchAddress("RawLucHitAND", RawLucHitAND, &b_RawLucHitAND);
      t_other->SetBranchAddress("RawBcmHOR", RawBcmHOR, &b_RawBcmHOR);
      t_other->SetBranchAddress("RawBcmHAND", RawBcmHAND, &b_RawBcmHAND);
      t_other->SetBranchAddress("RawBcmHXORA", RawBcmHXORA, &b_RawBcmHXORA);
      t_other->SetBranchAddress("RawBcmHXORC", RawBcmHXORC, &b_RawBcmHXORC);
      t_other->SetBranchAddress("RawBcmVOR", RawBcmVOR, &b_RawBcmVOR);
      t_other->SetBranchAddress("RawBcmVAND", RawBcmVAND, &b_RawBcmVAND);
      t_other->SetBranchAddress("RawBcmVXORA", RawBcmVXORA, &b_RawBcmVXORA);
      t_other->SetBranchAddress("RawBcmVXORC", RawBcmVXORC, &b_RawBcmVXORC);


      /*
            UInt_t       Status[3564];
            Float_t         MuLucOR[3564];
            Float_t         MuLucAND[3564];
            Float_t         MuLucHitOR[3564];
            Float_t         MuLucHitAND[3564];
            Float_t         MuBcmHOR[3564];
            Float_t         MuBcmHAND[3564];
            Float_t         MuBcmHORC[3564];
            Float_t         MuBcmHAND25[3564];
            Float_t         MuBcmVOR[3564];
            Float_t         MuBcmVAND[3564];
            Float_t         MuBcmVORC[3564];
            Float_t         MuBcmVAND25[3564];
            TBranch      *b_Status;
            TBranch        *b_MuLucOR;   //!
            TBranch        *b_MuLucAND;   //!
            TBranch        *b_MuLucHitOR;   //!
            TBranch        *b_MuLucHitAND;   //!
            TBranch        *b_MuBcmHOR;   //!
            TBranch        *b_MuBcmHAND;   //!
            TBranch        *b_MuBcmVOR;   //!
            TBranch        *b_MuBcmVAND;   //!
            t_other->SetBranchAddress("Status", Status, &b_Status);
            t_other->SetBranchAddress("MuLucOR", MuLucOR, &b_MuLucOR);
            t_other->SetBranchAddress("MuLucAND", MuLucAND, &b_MuLucAND);
            t_other->SetBranchAddress("MuLucHitOR", MuLucHitOR, &b_MuLucHitOR);
            t_other->SetBranchAddress("MuLucHitAND", MuLucHitAND, &b_MuLucHitAND);
            t_other->SetBranchAddress("MuBcmHOR", MuBcmHOR, &b_MuBcmHOR);
            t_other->SetBranchAddress("MuBcmHAND", MuBcmHAND, &b_MuBcmHAND);
            t_other->SetBranchAddress("MuBcmVOR", MuBcmVOR, &b_MuBcmVOR);
            t_other->SetBranchAddress("MuBcmVAND", MuBcmVAND, &b_MuBcmVAND);
      */

      Long64_t entries = t_other->GetEntriesFast();
      for (int i=0; i<entries; i++) {
        t_other->GetEntry(i);

        Int_t current_LB = t_other->GetLeaf("LB")->GetValue(0);
        UInt_t current_status = t_other->GetLeaf("Status")->GetValue(0);

        // This part is different depending on whether we want a per-BCID measurement. If pLBs are available, we can do BCID-by-BCID:
        if (!do_physics_run) {
          n_bcids[current_LB] = 1;
          for (vector<Int_t>::iterator bcid = bcidList.begin(); bcid != bcidList.end(); bcid++) {
            lumi_bcmhor[*bcid][current_LB] = -1. * TMath::Log(1. - RawBcmHOR[*bcid]) * 11245.5 / sigma_vis["bcmhor"];
            lumi_bcmhand[*bcid][current_LB] = -1. * TMath::Log(1. - RawBcmHAND[*bcid]) * 11245.5 / sigma_vis["bcmhand"];
            lumi_bcmvor[*bcid][current_LB] = -1. * TMath::Log(1. - RawBcmVOR[*bcid]) * 11245.5 / sigma_vis["bcmvor"];
            lumi_bcmvand[*bcid][current_LB] = -1. * TMath::Log(1. - RawBcmVAND[*bcid]) * 11245.5 / sigma_vis["bcmvand"];
            lumi_lucidor[*bcid][current_LB] = -1. * TMath::Log(1. - RawLucOR[*bcid]) * 11245.5 / sigma_vis["lucidor"];
            lumi_lucidand[*bcid][current_LB] = -1. * TMath::Log(1. - RawLucAND[*bcid]) * 11245.5 / sigma_vis["lucidand"];
            lumi_lucidhitor[*bcid][current_LB] = -1. * TMath::Log(1. - RawLucHitOR[*bcid]) * 11245.5 / sigma_vis["lucidhitor"];
            lumi_lucidhitand[*bcid][current_LB] = -1. * TMath::Log(1. - RawLucHitAND[*bcid]) * 11245.5 / sigma_vis["lucidhitand"];

            /*
            lumi_bcmhor[*bcid][current_LB] = MuBcmHOR[*bcid] * 11245.5 / 71500.;
            lumi_bcmhand[*bcid][current_LB] = MuBcmHAND[*bcid] * 11245.5 / 71500.;
            lumi_bcmvor[*bcid][current_LB] = MuBcmVOR[*bcid] * 11245.5 / 71500.;
            lumi_bcmvand[*bcid][current_LB] = MuBcmVAND[*bcid] * 11245.5 / 71500.;
            lumi_lucidor[*bcid][current_LB] = MuLucOR[*bcid] * 11245.5 / 71500.;
            lumi_lucidand[*bcid][current_LB] = MuLucAND[*bcid] * 11245.5 / 71500.;
            lumi_lucidhitor[*bcid][current_LB] = MuLucHitOR[*bcid] * 11245.5 / 71500.;
            lumi_lucidhitand[*bcid][current_LB] = MuLucHitAND[*bcid] * 11245.5 / 71500.;
            */
          }
        } else {
          // Otherwise, we need to take an average. This will maybe take a while, but nothing to do about it.
          n_bcids[current_LB] = 0;
          for (int j = 0; j < 3563; j++) { 
            if (Status[j] == 0) {
              continue;
            }

            n_bcids[current_LB]++;
            /*
              Recompute the luminosity using the raw counting rates.
              -Log(1-raw_rate) = mu_vis
              Use most recent cross sections.
            */


            if (!lumi_bcmhor[0][current_LB]) {
              lumi_bcmhor[0][current_LB] = -1. * TMath::Log(1. - RawBcmHOR[j]) * 11245.5 / sigma_vis["bcmhor"];
              lumi_bcmhand[0][current_LB] = -1. * TMath::Log(1. - RawBcmHAND[j]) * 11245.5 / sigma_vis["bcmhand"];
              lumi_bcmvor[0][current_LB] = -1. * TMath::Log(1. - RawBcmVOR[j]) * 11245.5 / sigma_vis["bcmvor"];
              lumi_bcmvand[0][current_LB] = -1. * TMath::Log(1. - RawBcmVAND[j]) * 11245.5 / sigma_vis["bcmvand"];
              lumi_lucidor[0][current_LB] = -1. * TMath::Log(1. - RawLucOR[j]) * 11245.5 / sigma_vis["lucidor"];
              lumi_lucidand[0][current_LB] = -1. * TMath::Log(1. - RawLucAND[j]) * 11245.5 / sigma_vis["lucidand"];
              lumi_lucidhitor[0][current_LB] = -1. * TMath::Log(1. - RawLucHitOR[j]) * 11245.5 / sigma_vis["lucidhitor"];
              lumi_lucidhitand[0][current_LB] = -1. * TMath::Log(1. - RawLucHitAND[j]) * 11245.5 / sigma_vis["lucidhitand"];

              /*
              lumi_bcmhor[0][current_LB] = MuBcmHOR[j] * 11245.5 / 71500.;
              lumi_bcmhand[0][current_LB] = MuBcmHAND[j] * 11245.5 / 71500.;
              lumi_bcmvor[0][current_LB] = MuBcmVOR[j] * 11245.5 / 71500.;
              lumi_bcmvand[0][current_LB] = MuBcmVAND[j] * 11245.5 / 71500.;
              lumi_lucidor[0][current_LB] = MuLucOR[j] * 11245.5 / 71500.;
              lumi_lucidand[0][current_LB] = MuLucAND[j] * 11245.5 / 71500.;
              lumi_lucidhitor[0][current_LB] = MuLucHitOR[j] * 11245.5 / 71500.;
              lumi_lucidhitand[0][current_LB] = MuLucHitAND[j] * 11245.5 / 71500.;
              */
            } else {
              lumi_bcmhor[0][current_LB] += -1. * TMath::Log(1. - RawBcmHOR[j]) * 11245.5 / sigma_vis["bcmhor"];
              lumi_bcmhand[0][current_LB] += -1. * TMath::Log(1. - RawBcmHAND[j]) * 11245.5 / sigma_vis["bcmhand"];
              lumi_bcmvor[0][current_LB] += -1. * TMath::Log(1. - RawBcmVOR[j]) * 11245.5 / sigma_vis["bcmvor"];
	       		  lumi_bcmvand[0][current_LB] += -1. * TMath::Log(1. - RawBcmVAND[j]) * 11245.5 / sigma_vis["bcmvand"];
              lumi_lucidor[0][current_LB] += -1. * TMath::Log(1. - RawLucOR[j]) * 11245.5 / sigma_vis["lucidor"];
              lumi_lucidand[0][current_LB] += -1. * TMath::Log(1. - RawLucAND[j]) * 11245.5 / sigma_vis["lucidand"];
              lumi_lucidhitor[0][current_LB] += -1. * TMath::Log(1. - RawLucHitOR[j]) * 11245.5 / sigma_vis["lucidhitor"];
              lumi_lucidhitand[0][current_LB] += -1. * TMath::Log(1. - RawLucHitAND[j]) * 11245.5 / sigma_vis["lucidhitand"];

              /*
              lumi_bcmhor[0][current_LB] += MuBcmHOR[j] * 11245.5 / 71500.;
              lumi_bcmhand[0][current_LB] += MuBcmHAND[j] * 11245.5 / 71500.;
              lumi_bcmvor[0][current_LB] += MuBcmVOR[j] * 11245.5 / 71500.;
              lumi_bcmvand[0][current_LB] += MuBcmVAND[j] * 11245.5 / 71500.;
              lumi_lucidor[0][current_LB] += MuLucOR[j] * 11245.5 / 71500.;
              lumi_lucidand[0][current_LB] += MuLucAND[j] * 11245.5 / 71500.;
              lumi_lucidhitor[0][current_LB] += MuLucHitOR[j] * 11245.5 / 71500.;
              lumi_lucidhitand[0][current_LB] += MuLucHitAND[j] * 11245.5 / 71500.;
              */
            }
          } // end BCID loop
          
          // Divide by nBCIDS - this gives an average lumi over all BCIDs
          //lumi_bcmhor[0][current_LB] /= n_bcids[current_LB];
          //lumi_bcmhand[0][current_LB] /= n_bcids[current_LB];
          //lumi_bcmvor[0][current_LB] /= n_bcids[current_LB];
          //cout << "[GetLuminosity] Line 1047 n_bcids[current_LB] = " << n_bcids[current_LB] << endl;
          //lumi_bcmvand[0][current_LB] /= n_bcids[current_LB];
          //lumi_lucidor[0][current_LB] /= n_bcids[current_LB];
          //lumi_lucidand[0][current_LB] /= n_bcids[current_LB];
          //lumi_lucidhitor[0][current_LB] /= n_bcids[current_LB];
          //lumi_lucidhitand[0][current_LB] /= n_bcids[current_LB];

          // Keep a persistent copy of the number of bunches.
          if (n_bcids[current_LB] > 0) {
            if (n_bcids_run == 0) {
              n_bcids_run = n_bcids[current_LB];
            } else {
              if (n_bcids[current_LB] != n_bcids_run) {
                cout << "[GetLuminosity] WARNING : LB " << current_LB << " has " << n_bcids[current_LB] << " bunches, while the number saved is " << n_bcids_run << endl;
              }
            }
          }

        } // end if/else for physics/special run
      }
    }
  }

  /**********************************************************************/
  // Calculations
  /**********************************************************************/
  //Visible cross sections
  std::map<Int_t, Float_t> sigmavis_nvtx;
  std::map<Int_t, Float_t> sigmavis_nevt;

  /* R16, 7 TeV
    sigmavis_nvtx[2] = 51620;
    sigmavis_nvtx[5] = 38570;
    sigmavis_nvtx[7] = 32190;
    sigmavis_nvtx[10] = 25280;

    sigmavis_nevt[2] = 51160;
    sigmavis_nevt[5] = 38590;
    sigmavis_nevt[7] = 32210;
    sigmavis_nevt[10] = 25310;

  */
  /* R17, 7 TeV*//*
  sigmavis_nvtx[2] = 1000000;
  sigmavis_nvtx[3] = 1000000;
  sigmavis_nvtx[4] = 1000000;
  sigmavis_nvtx[5] = 38016.6;
  sigmavis_nvtx[7] = 31602.2;
  sigmavis_nvtx[10] = 24672.8;

  sigmavis_nevt[2] = 50995;
  sigmavis_nevt[3] = 50995;
  sigmavis_nevt[4] = 50995;
  sigmavis_nevt[5] = 38015.3;
  sigmavis_nevt[7] = 31614.7;
  sigmavis_nevt[10] = 24693.2;
*/

  /* R17, 8 TeV from MC */
  /*
    sigmavis_nvtx[2] = 100000000;
    sigmavis_nvtx[5] = 39240.9;
    sigmavis_nvtx[7] = 32800.3;
    sigmavis_nvtx[10] = 25815.22;

    sigmavis_nevt[2] = 10000000.7;
    sigmavis_nevt[5] = 39239.8;
    sigmavis_nevt[7] = 32813.1;
    sigmavis_nevt[10] = 25836.6;
  */
  /* R17.2, default settings, 8 TeV from VdM scan */
  /*
    sigmavis_nvtx[2] = 100000000;
    sigmavis_nvtx[5] = 38471.1;
    sigmavis_nvtx[7] = 32244.7;
    sigmavis_nvtx[10] = 25427.6;

    sigmavis_nevt[2] = 49960.0;
    sigmavis_nevt[5] = 38463.9;
    sigmavis_nevt[7] = 32244.7;
    sigmavis_nevt[10] = 25449.5;
  */

  //July Scans
  /*sigmavis_nvtx[3]= 24561.8344976;
	sigmavis_nvtx[4]= 20030.5455829;
	sigmavis_nvtx[5]= 16488.4822889;
	sigmavis_nvtx[6]= 13711.0841629;
	sigmavis_nvtx[7]= 11513.065241;
	sigmavis_nvtx[8]= 9631.63929099;
	sigmavis_nvtx[10]= 6677.63028981;*/
	//November Scans
	/*sigmavis_nvtx[3]= 25253.5936583;
	sigmavis_nvtx[4]= 20568.9567171;
	sigmavis_nvtx[5]= 16938.1358851;
	sigmavis_nvtx[6]= 14114.5976408;
	sigmavis_nvtx[7]= 11824.8763754;
	sigmavis_nvtx[8]= 9888.40189616;
	sigmavis_nvtx[10]= 6852.46564036;*/
	//Witold Presentation: various corrections applied
	sigmavis_nvtx[5] = 16879.0;
      				

  VtxCalibration *vc = new VtxCalibration(gs.path_outputVdM);
  vc->Initialize();
  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    for( vector<TString>::const_iterator itVtxSetting = vtx_settings.begin(); itVtxSetting != vtx_settings.end(); ++itVtxSetting ) {
    
      TString algname( *itVtxSetting );
      TString detailed_settings = settings;
      /*
       *
       * NVtx counting
       *
       */    
      if( algname == "NVtx" ) {
        algname += *nTrkCut;

        for (map<TString, bool>::iterator it = systematic_uncertainty_list.begin(); it != systematic_uncertainty_list.end(); ++it) {
          if ((*it).second) {
            detailed_settings += "_";
            detailed_settings += (*it).first;
          }
        }
        sigma_vis[algname] = sigmavis_nvtx[*nTrkCut];  
        //sigma_vis[algname] = 1000. * vc->GetCrossSection(atof(energy), *nTrkCut, detailed_settings, "NVtx");
        cout << "[GetLuminosity] INFO : sigma_vis (NVtx, track cut = "<< *nTrkCut <<") = " << sigma_vis[algname] << endl;
      }

      /*
       *
       * NEvt counting
       *
       */    
      else if( algname == "NEvt" ) {
        algname += *nTrkCut;
        sigma_vis[algname] = 1000. * vc->GetCrossSection(atof(energy), *nTrkCut, detailed_settings, "NEvt");
        cout << "[GetLuminosity] INFO : sigma_vis (NEvt) = " << sigma_vis[algname] << endl;
      }
      
      /*
       *
       * NUnf counting
       *
       */    
      else if( algname == "NUnf" ) {
        algname += *nTrkCut;
        #ifdef USE_UNFOLD
        //if not USE_UNFOLD it's set already to 0.0
        sigma_vis[algname] = 1000. * vc->GetCrossSection(atof(energy), *nTrkCut, detailed_settings, "NUnf");
        sigma_vis_err[algname] = 1000. * vc->GetCrossSectionError(atof(energy), *nTrkCut, detailed_settings, "NUnf");
        #endif
        cout << "[GetLuminosity] INFO : sigma_vis (" << algname << ") = "
             << sigma_vis[algname] << " +/- " << sigma_vis_err[algname]<< endl;
      }
    }
  }

  std::map<Int_t, std::map<Int_t, LumiVtx*> > luminosityContainer;

  for (vector<Int_t>::iterator bcid = bcidList.begin(); bcid != bcidList.end(); ++bcid) {
    for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {

      TString detailed_settings = settings;
      for (map<TString, bool>::iterator it = systematic_uncertainty_list.begin(); it != systematic_uncertainty_list.end(); ++it) {
        if ((*it).second) {
          detailed_settings += "_";
          detailed_settings += (*it).first;
        }
      }

      luminosityContainer[*nTrkCut][*bcid] = new LumiVtx();

      for( vector<TString>::const_iterator itVtxSetting = vtx_settings.begin(); itVtxSetting != vtx_settings.end(); ++itVtxSetting ) {    
        TString algname( *itVtxSetting );

        if( algname = "NVtx" ) {
          algname += *nTrkCut;
          cout << "[GetLuminosity] Line 1193 sigma_vis["<<algname<<"] = " << sigma_vis[algname] << endl;
          luminosityContainer[*nTrkCut][*bcid]->SetSigmaVisNVtx(sigma_vis[algname]);
        } else if ( algname == "NEvt" ) {
          algname += *nTrkCut;
          luminosityContainer[*nTrkCut][*bcid]->SetSigmaVisNEvt(sigma_vis[algname]);
        } else if ( algname == "NUnf" ) {
          algname += *nTrkCut;
          luminosityContainer[*nTrkCut][*bcid]->SetSigmaVis(GlobalSettings::kUnfC, sigma_vis[algname], sigma_vis_err[algname]);
        }
      }

      for (std::map<TString, bool>::iterator it = systematic_uncertainty_list.begin(); it != systematic_uncertainty_list.end(); ++it) {
       	if ((*it).second) {
          luminosityContainer[*nTrkCut][*bcid]->SetSystematicUncertaintyFlag((*it).first);
        }
      }

      if (do_physics_run) { 
				n_bcids_run = 1; //changed this from 1 to 1368 (colliding bunches from ATLAS Run Query) back to 1
        luminosityContainer[*nTrkCut][*bcid]->SetNBunches(n_bcids_run);
      }

      std::vector<std::pair<TString, Int_t> > selection;
      if (!do_physics_run) {
        selection.push_back(std::make_pair<TString, Int_t>("BunchId", *bcid));
      }

      TString b_nvtx_name = "NVtx";
      b_nvtx_name += *nTrkCut;
      TString b_nevt_name = "NEvt";
      b_nevt_name += *nTrkCut;
      #ifdef NOSPLIT
      b_nvtx_name += "_BSC";
      b_nevt_name += "_BSC";
      #else
      b_nvtx_name += "_ASC";
      b_nevt_name += "_ASC";
      #endif
            
      std::cout << "[GetLuminosity] INFO : Loading luminosity container for nTrkCut="<<*nTrkCut<<", BCID="<<*bcid << std::endl;      
      luminosityContainer[*nTrkCut][*bcid]->LoadTTree(path_vertex_counts, b_nvtx_name, b_nevt_name, "NTrig", &selection, do_physics_run);
      std::cout << "[GetLuminosity] INFO : Finished loading luminosity container for nTrkCut="<<*nTrkCut<<", BCID="<<*bcid << std::endl; 
     
      if (timestamp_paths.size() > 0) {
        for (vector<TString>::iterator itf = timestamp_paths.begin(); itf != timestamp_paths.end(); ++itf) {
          luminosityContainer[*nTrkCut][*bcid]->LoadTimestamps(*itf);
        }
      } else if (duration_paths.size() > 0) {
        for (vector<TString>::iterator itf = duration_paths.begin(); itf != duration_paths.end(); ++itf) {
          luminosityContainer[*nTrkCut][*bcid]->LoadDurations(*itf);
        }
      //} else if ( do_physics_run ) {
      //  // Don't care about instantaneous luminosities - set all lumiblock durations to 1
      //  //luminosityContainer[*nTrkCut][*bcid]->SetAllDurations(0, 10000);
      } else {
        // Assume that the lumi ntuple is available
        for (vector<TString>::iterator itf = lumi_ntuple_paths.begin(); itf != lumi_ntuple_paths.end(); ++itf) {
          luminosityContainer[*nTrkCut][*bcid]->LoadDurationsFromLumiNtuple(*itf);
        }
      }

      if ((run == "182013") || (run == "191373")) {
        luminosityContainer[*nTrkCut][*bcid]->LoadLivefraction(livefraction_path);
        luminosityContainer[*nTrkCut][*bcid]->LoadPrescale(livefraction_path);
      } else {
        luminosityContainer[*nTrkCut][*bcid]->LoadLivefraction();
        luminosityContainer[*nTrkCut][*bcid]->LoadPrescale();
      }

      TString triggerType;
      if ((run == "182013") || (run == "191373")) {
        luminosityContainer[*nTrkCut][*bcid]->CorrectDataRate("full");
        triggerType = "MBTS_2";
      } else if (do_physics_run){
        triggerType = "data";
				luminosityContainer[*nTrkCut][*bcid]->CorrectDataRate("data");
      } else{
        luminosityContainer[*nTrkCut][*bcid]->CorrectDataRate("fixedrate");
        triggerType = "BGRP7"; //we're assume random trigger basically (e.g. EF_rd0_Filled_NoAlg)
      }

      #ifdef USE_UNFOLD
      //we load and compute unfolding results here, in one step, since we need all corrections already loaded
      luminosityContainer[*nTrkCut][*bcid]->SetUnfResponseMatrix((gs.path_unf_ResponseMatrix+"/"+gs.fname_unf_ResponseMatrix).c_str(), gs.other_unf_ResponseMatrixHNameBase.c_str());
      luminosityContainer[*nTrkCut][*bcid]->LoadVertexCountsUnfold(path_vertex_histograms, *bcid, *nTrkCut, triggerType);
      #endif

      luminosityContainer[*nTrkCut][*bcid]->CalculateRates();

      luminosityContainer[*nTrkCut][*bcid]->InitializeFakeCorrection(energy, settings, *nTrkCut);

      if (mu_scale > 0) {
        luminosityContainer[*nTrkCut][*bcid]->SetMuScale(mu_scale);
      }

      TFile *f_histograms = new TFile(path_vertex_histograms, "READ");
      cout << "Line 1288 path_vertex_histograms = " << path_vertex_histograms << endl;
      TString hname = "hist/PriVtxZpLB_BCID";
      hname += *bcid;
      TH2D *h_z_plb = (TH2D*)f_histograms->Get(hname);

      std::cout << "[GetLuminosity] INFO : Correcting pileup effects for nTrkCut="<<*nTrkCut<<", BCID="<<*bcid << std::endl;      
      luminosityContainer[*nTrkCut][*bcid]->CorrectPileupEffects(h_z_plb, energy, settings, *nTrkCut, "207216");
      std::cout << "[GetLuminosity] INFO : Finished correcting pileup effects for nTrkCut="<<*nTrkCut<<", BCID="<<*bcid << std::endl;      


      #ifdef REMOVED_052612
      // Removed: I think it's better to calculate masking corrections individually in as small increments as possible, e.g. single pLBs, if possible.
      if (*nTrkCut != 2) {
        //Check for cached masking correction
        TString save_name( run ); save_name += "_BCID"; save_name += *bcid ; save_name += "_NTrk" ; save_name += *nTrkCut ; save_name += "_" ; save_name += energy ; save_name += "_" ; save_name += settings ;
        TString cache_path = "cache/";
        cache_path += save_name;
        cache_path += ".root";

        if (fileExists(cache_path) && !redo_masking_correction) {
          // Cache exists
          cout << "Loading cached masking correction" << endl;
          TFile *f_cache = new TFile(cache_path, "READ");
          luminosityContainer[*nTrkCut][*bcid]->InitializeMaskingCorrection(energy, settings, *nTrkCut, save_name);
        } else {
          // Make masking correction from scratch
          cout << "Generating new masking correction" << endl;

          // - Find z-distribution
          TString name_z = "hist/PriVtxZpLB_BCID";
          name_z += *bcid;
          TFile *f_z = new TFile(path_vertex_histograms, "READ");
          TH2D *h_z_pLB = (TH2D*)f_z->Get(name_z);
          h_z_pLB->SetDirectory(0);
          TH1D *h_z = (TH1D*)h_z_pLB->ProjectionX();

          // - Initialize masking correction
          luminosityContainer[*nTrkCut][*bcid]->InitializeMaskingCorrection(energy, settings, *nTrkCut, h_z, save_name);
          f_z->Close();
        }
      }
      #endif

//#ifdef REMOVED_053012
//      TString objectname = "tg_MC_truthmap_mu_NTrk";
//      objectname += *nTrkCut;
//      luminosityContainer[*nTrkCut][*bcid]->ApplyMcCorrection("/u/dryu/Luminosity/muScan/VertexStudies/mc_correction/D3PD/results///MC_truthmap.root", objectname);
//#endif

      luminosityContainer[*nTrkCut][*bcid]->CalculateAllLbLum();
      luminosityContainer[*nTrkCut][*bcid]->CalculateAllRunLum();
    } // End nTrkCut loop

    //Load Stefan's TGraph
    /*
        TGraph *tg_bcm;
        if (run == "182013") {
          TString bcm_filename = "input/lumi_bcm_"; bcm_filename += run; bcm_filename += ".root";
          TFile *f_bcm = new TFile(bcm_filename, "READ");
          TString bcm_name = "BCMH_OR_LumiCorr_"; bcm_name += *bcid;
          tg_bcm = (TGraph*)f_bcm->Get(bcm_name);
          if (tg_bcm) {
            for (int i=0; i<tg_bcm->GetN(); i++) {
              lumi_bcm[*bcid][tg_bcm->GetX()[i]] = tg_bcm->GetY()[i];
            }
          }
        }
    */


    /**********************************************************************/
    //  Write to output tree
    /**********************************************************************/

    cout << "[GetLuminosity] Calculation done." << endl;
    cout << "Total luminosity summary, BCID = " << *bcid << endl;
    for( vector<TString>::const_iterator itVtxSetting = vtx_settings.begin(); itVtxSetting != vtx_settings.end(); ++itVtxSetting ) {
      if( *itVtxSetting == "NVtx" ) {
        cout << "- Vertex counting:"  << endl;
        for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
          cout << " NTrkCut: " << *nTrkCut << "; L = "
               << luminosityContainer[*nTrkCut][*bcid]->GetRunTotalLum(GlobalSettings::kVtxC) << endl;
        }
      } else if ( *itVtxSetting == "NEvt" ) {
        cout << "- Event counting:"  << endl;
        for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
          cout << " NTrkCut: " << *nTrkCut << "; L = "
             << luminosityContainer[*nTrkCut][*bcid]->GetRunTotalLum(GlobalSettings::kEvtC) << endl;
        }
      } else if ( *itVtxSetting == "NUnf" ) {
        cout << "- Unfolding counting:"  << endl;
        for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
          cout << " NTrkCut: " << *nTrkCut << "; L = "
             << luminosityContainer[*nTrkCut][*bcid]->GetRunTotalLum(GlobalSettings::kUnfC) << endl;
        }
      }
    }

    cout << "\n" << "Creating output Tree for saving results." << endl;
    Int_t n_pLBs = luminosityContainer[nTrkCuts[0]][*bcid]->GetNLb();
    current_run = atoi(run);
    
    for (Int_t i = 0; i < n_pLBs; i++) {
      current_bcid = *bcid;
/*Test*/
      current_lb = luminosityContainer[nTrkCuts[0]][*bcid]->GetLbFromDist(i);
	  const int current_lb2 = luminosityContainer[nTrkCuts[0]][*bcid]->GetLbFromDist(i);
      if (find(plb_skiplist.begin(), plb_skiplist.end(), current_lb2) != plb_skiplist.end()) {
        continue;
      }
      current_lb_duration = luminosityContainer[nTrkCuts[0]][*bcid]->GetLbDuration(i);

      for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
        int j = distance(nTrkCuts.begin(), nTrkCut);

        current_lb_inst_luminosity_nvtx[j] = luminosityContainer[*nTrkCut][*bcid]->GetLbInstLum(i);
        cout << "current_lb_inst_luminosity_nvtx[" << j << "] = " << current_lb_inst_luminosity_nvtx[j] << endl;
        current_lb_inst_luminosity_nvtx_err[j] = luminosityContainer[*nTrkCut][*bcid]->GetLbInstLumError(i);

        current_lb_inst_luminosity_nevt[j] = luminosityContainer[*nTrkCut][*bcid]->GetEvtLbInstLum(i);
        current_lb_inst_luminosity_nevt_errdown[j] = (luminosityContainer[*nTrkCut][*bcid]->GetEvtLbInstLumError(i)).first;
        current_lb_inst_luminosity_nevt_errup[j] = (luminosityContainer[*nTrkCut][*bcid]->GetEvtLbInstLumError(i)).second;

        current_lb_nvtx[j] = luminosityContainer[*nTrkCut][*bcid]->GetNVtx(i);

        if (luminosityContainer[*nTrkCut][*bcid]->GetNVtxMasked(i)) {
          current_lb_nvtx_masked[j] = luminosityContainer[*nTrkCut][*bcid]->GetNVtxMasked(i);
        }

        current_lb_nevt[j] = luminosityContainer[*nTrkCut][*bcid]->GetNEvt(i);

        current_lb_nvtx_from_nevt[j] = luminosityContainer[*nTrkCut][*bcid]->GetNVtxFromNEvt(i);

        current_lb_masking_correction[j] = luminosityContainer[*nTrkCut][*bcid]->GetMaskingCorrection(i);
        current_lb_fake_mu_nvtx[j] = luminosityContainer[*nTrkCut][*bcid]->GetFakeMu(i, "NVtx");
        current_lb_fake_mu_nevt[j] = luminosityContainer[*nTrkCut][*bcid]->GetFakeMu(i, "NEvt");

        current_lb_inst_luminosity_nunf[j] = luminosityContainer[*nTrkCut][*bcid]->GetLbInstLum(i, GlobalSettings::kUnfC);
        current_lb_inst_luminosity_nunf_err[j] = luminosityContainer[*nTrkCut][*bcid]->GetLbInstLumError(i, GlobalSettings::kUnfC);
        current_lb_total_luminosity_nunf[j] = luminosityContainer[*nTrkCut][*bcid]->GetLbTotalLum(i, GlobalSettings::kUnfC);
        current_lb_total_luminosity_nunf_err[j] = luminosityContainer[*nTrkCut][*bcid]->GetLbTotalLumError(i, GlobalSettings::kUnfC);
      }
      if (lumi_bcmhor[*bcid].count(current_lb2) > 0) {
        current_lb_inst_luminosity_bcmhor = lumi_bcmhor[*bcid][current_lb2];
      } else {
        current_lb_inst_luminosity_bcmhor = 0.;
      }
      if (lumi_bcmhand[*bcid].count(current_lb2) > 0) {
        current_lb_inst_luminosity_bcmhand = lumi_bcmhand[*bcid][current_lb2];
      } else {
        current_lb_inst_luminosity_bcmhand = 0.;
      }
      if (lumi_bcmvor[*bcid].count(current_lb2) > 0) {
        current_lb_inst_luminosity_bcmvor = lumi_bcmvor[*bcid][current_lb2];
      } else {
        current_lb_inst_luminosity_bcmvor = 0.;
      }
      if (lumi_bcmvand[*bcid].count(current_lb2) > 0) {
        current_lb_inst_luminosity_bcmvand = lumi_bcmvand[*bcid][current_lb2];
      } else {
        current_lb_inst_luminosity_bcmvand = 0.;
      }

      if (lumi_lucidor[*bcid].count(current_lb2) > 0) {
        current_lb_inst_luminosity_lucidor = lumi_lucidor[*bcid][current_lb2];
      } else {
        current_lb_inst_luminosity_lucidor = 0.;
      }

      if (lumi_lucidand[*bcid].count(current_lb2) > 0) {
        current_lb_inst_luminosity_lucidand = lumi_lucidand[*bcid][current_lb2];
      } else {
        current_lb_inst_luminosity_lucidand = 0.;
      }
  
      if (lumi_lucidhitor[*bcid].count(current_lb2) > 0) {
        current_lb_inst_luminosity_lucidhitor = lumi_lucidhitor[*bcid][current_lb2];
      } else {
        current_lb_inst_luminosity_lucidhitor = 0.;
      }

      if (lumi_lucidhitand[*bcid].count(current_lb2) > 0) {
        current_lb_inst_luminosity_lucidhitand = lumi_lucidhitand[*bcid][current_lb2];
      } else {
        current_lb_inst_luminosity_lucidhitand = 0.;
      }

      // -- Errors
	 
      if (lumi_bcmhor_err[*bcid].count(current_lb2) > 0) {
        current_lb_inst_luminosity_bcmhor_err = lumi_bcmhor_err[*bcid][current_lb2];
      } else {
        current_lb_inst_luminosity_bcmhor_err = 0.;
      }

      if (lumi_bcmhand_err[*bcid].count(current_lb2) > 0) {
        current_lb_inst_luminosity_bcmhand_err = lumi_bcmhand_err[*bcid][current_lb2];
      } else {
        current_lb_inst_luminosity_bcmhand_err = 0.;
      }

      if (lumi_bcmvor_err[*bcid].count(current_lb2) > 0) {
        current_lb_inst_luminosity_bcmvor_err = lumi_bcmvor_err[*bcid][current_lb2];
      } else {
        current_lb_inst_luminosity_bcmvor_err = 0.;
      }

      if (lumi_bcmvand_err[*bcid].count(current_lb2) > 0) {
        current_lb_inst_luminosity_bcmvand_err = lumi_bcmvand_err[*bcid][current_lb2];
      } else {
        current_lb_inst_luminosity_bcmvand_err = 0.;
      }

      if (lumi_lucidor_err[*bcid].count(current_lb2) > 0) {
        current_lb_inst_luminosity_lucidor_err = lumi_lucidor_err[*bcid][current_lb2];
      } else {
        current_lb_inst_luminosity_lucidor_err = 0.;
      }

      if (lumi_lucidand_err[*bcid].count(current_lb2) > 0) {
        current_lb_inst_luminosity_lucidand_err = lumi_lucidand_err[*bcid][current_lb2];
      } else {
        current_lb_inst_luminosity_lucidand_err = 0.;
      }

      if (lumi_lucidhitor_err[*bcid].count(current_lb2) > 0) {
        current_lb_inst_luminosity_lucidhitor_err = lumi_lucidhitor_err[*bcid][current_lb2];
      } else {
        current_lb_inst_luminosity_lucidhitor_err = 0.;
      }

      if (lumi_lucidhitand_err[*bcid].count(current_lb2) > 0) {
        current_lb_inst_luminosity_lucidhitand_err = lumi_lucidhitand_err[*bcid][current_lb2];
      } else {
        current_lb_inst_luminosity_lucidhitand_err = 0.;
      }

      t_lb->Fill();
    }
  }


  //Save TTrees
  TString path_output = output_prefix;
  path_output += ".root";
  cout << "[GetLuminosity] INFO : Saving output TTree to: " << path_output << endl;
  TFile *f_out = new TFile(path_output,"RECREATE");
  t_lb->Write();
  f_out->Close();

  cout << "All done." << endl;

}

void usage(const char *pname) {
  std::cout << pname << " -r run -s settings [options]" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "\t -r, --run Run: Set run number to Run" << endl;
  std::cout << "\t -s, --settings Settings: Set reconstruction setup Settings" << endl;
  //std::cout << "\t -f, --inTXT InputTextFile: Path to tab-separated text file containing TAG/metadata paths." << std::endl;
  std::cout << "\t -p path: Specify output path, outPrefix will be biult as path/run/settings[/syst]/lumi_r" << std::endl;
  std::cout << "\t -o, --outPrefix OutPrefix: Specify full output prefix for output files. Directories must already exist!" << std::endl;
  std::cout << "\t -v, --verbose: Enable verbose printout." << std::endl;
  std::cout << "\t -m: Force re-making pileup correction." << std::endl;
  std::cout << "\t -f Scale: Enable mu-scaling in fake correction of Scale" << std::endl;
  std::cout << "\t -u Syst: Enable systematics Syst" << std::endl;
  std::cout << "\t -n nTrkCut: Run only with nTrkCut track cut" << std::endl;
  std::cout << "\t -i bcid: Run only on the given bcid" << std::endl;
  std::cout << "\t -h, --help: Print this screen." << std::endl;
}


bool fileExists(TString path) {

  ifstream ifile(path);
  if (ifile.is_open()) {
    ifile.close();
    return true;
  } else {
    return false;
  }
}
