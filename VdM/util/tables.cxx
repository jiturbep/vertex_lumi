#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <stdlib.h>
#include <string>

#include "GlobalSettings/GlobalSettings.h"

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
#include "TLegend.h"

#include "atlasstyle/AtlasStyle.h"
#include "atlasstyle/AtlasUtils.h"
#include "atlasstyle/AtlasLabels.h"

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<double>+;
#else
template class vector<double>;
#endif

int main(int argc, char **argv) {
  SetAtlasStyle();

  TString vtx_method;
  TString run;
  TString settings;
  // --- Scan for command line parameters
  int c;
  extern int optind;
  extern char* optarg;
  while (1) {
    int option_index = 0;
    static struct option long_options[] = { {0,0,0,0} };
    c=getopt_long(argc, argv, "r:s:m:",long_options,&option_index);
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
      case 'm': {
        vtx_method = optarg;
        break;
      }
    }
  }

  stringstream ss_folder;
  ss_folder << GlobalSettings::path_outputVdM << GlobalSettings::path_VdM_prefix << run << "/" << settings << "/" << vtx_method << "/";
  TString filename(ss_folder.str());
  filename += "vdm_results.root";

  cout << "Opening file " << filename << endl;
  TFile *f_in = new TFile(filename, "READ");
  if (!f_in->IsOpen()) {
    cerr << "Unable to open input file " << filename << ". Exiting..." << endl;
    exit(1);
  }

  TString default_fit;
  if (run == "182013") {
    default_fit = "sgcl";
  } else if (run == "201351") {
    default_fit = "dgc";
  } else if (run == "207216") {
    default_fit = "dgc";
  } else if (run == "207219") {
    default_fit = "dgc";
  } else if (run == "214984") {
    default_fit = "sgc";
  } else if (run == "215021") {
    default_fit = "sgc";
  }

  vector<TString> fit_functions;
  //fit_functions.push_back("sgcl");
  fit_functions.push_back("sgc");
  //fit_functions.push_back("sg");
  //fit_functions.push_back("dg");
  fit_functions.push_back("dgc");
  //fit_functions.push_back("dgcl");
  //fit_functions.push_back("spline");

  vector<Int_t> bcids;
  if (run == "182013") {
    bcids.push_back(81);
    bcids.push_back(867);
    bcids.push_back(2752);
  } else if (run == "201351") {
    bcids.push_back(1);
    bcids.push_back(241);
    bcids.push_back(2881);
    bcids.push_back(3121);
  }
  else if (run == "207216") {
    bcids.push_back(1);
    bcids.push_back(721);
    bcids.push_back(1821);
  }
  else if (run == "207219") {
    bcids.push_back(1);
    bcids.push_back(721);
    bcids.push_back(1821);
  }
  else if (run == "214984") {
    bcids.push_back(1);
    bcids.push_back(2361);
    bcids.push_back(2881);
  }
  else if (run == "215021") {
    bcids.push_back(1);
    bcids.push_back(2361);
    bcids.push_back(2881);
  }

  vector<Int_t> scans;
  if (run == "182013") {
    scans.push_back(1);
    scans.push_back(2);
  } else if (run == "201351") {
    scans.push_back(1);
    scans.push_back(2);
    scans.push_back(3);
  }
  else if (run == "207216") {
    scans.push_back(4);
    scans.push_back(5);
    scans.push_back(6);
  }
  else if (run == "207219") {
    scans.push_back(8);
  }
  else if (run == "214984") {
    scans.push_back(10);
    scans.push_back(11);
    scans.push_back(14);
  }
  else if (run == "215021") {
    scans.push_back(15);
  }

  vector<Int_t> nTrkCuts;
  //if (vtx_method == "NEvt") nTrkCuts.push_back(2);
  //nTrkCuts.push_back(2);
  //nTrkCuts.push_back(3);
  //nTrkCuts.push_back(4);
  nTrkCuts.push_back(5);
  //nTrkCuts.push_back(6);
  //nTrkCuts.push_back(7);
  //nTrkCuts.push_back(8);
  //nTrkCuts.push_back(10);

  // -- Read TTree into maps
  std::map<TString, Float_t> m_sigma_vis, m_sigma_vis_err, m_lumi_sp, m_lumi_sp_err, m_lumi, m_lumi_err;
  std::map<TString, Float_t> m_Sigma_x, m_Sigma_x_err, m_Sigma_y, m_Sigma_y_err, m_mu_max_x, m_mu_max_x_err, m_mu_max_y, m_mu_max_y_err, m_c_x, m_c_x_err, m_c_y, m_c_y_err, m_r_x, m_r_y;
  std::map<TString, Float_t> m_chi2ndf_x, m_chi2ndf_y;

  TTree *t_vdm = (TTree*)f_in->Get("VdmResults");
  if (!t_vdm) {
    cerr << "TTree VdMResults not found in input file! Exiting..." << endl;
    exit(1);
  }

  // Declaration of leaf types
  Int_t           Scan;
  Int_t           BCID;
  Int_t           NTrkCut;
  vector<double>  *SigmaVis;
  vector<double>  *SigmaVisErr;
  vector<double>  *LumiSp;
  vector<double>  *LumiSpErr;
  vector<double>  *Lumi;
  vector<double>  *LumiErr;
  vector<double>  *SigmaX;
  vector<double>  *SigmaY;
  vector<double>  *MuMaxX;
  vector<double>  *MuMaxY;
  vector<double>  *CX;
  vector<double>  *CY;
  vector<double>  *RX;
  vector<double>  *RY;
  vector<double>  *Chi2NdfX;
  vector<double>  *Chi2NdfY;
  vector<double>  *SigmaXErr;
  vector<double>  *SigmaYErr;
  vector<double>  *MuMaxXErr;
  vector<double>  *MuMaxYErr;
  vector<double>  *CXErr;
  vector<double>  *CYErr;

  // List of branches
  TBranch        *b_Scan;   //!
  TBranch        *b_BCID;   //!
  TBranch        *b_NTrkCut;   //!
  TBranch        *b_SigmaVis;   //!
  TBranch        *b_SigmaVisErr;   //!
  TBranch        *b_LumiSp;   //!
  TBranch        *b_LumiSpErr;   //!
  TBranch        *b_Lumi;   //!
  TBranch        *b_LumiErr;   //!
  TBranch        *b_SigmaX;   //!
  TBranch        *b_SigmaY;   //!
  TBranch        *b_MuMaxX;   //!
  TBranch        *b_MuMaxY;   //!
  TBranch        *b_CX;   //!
  TBranch        *b_CY;   //!
  TBranch        *b_RX;   //!
  TBranch        *b_RY;   //!
  TBranch        *b_Chi2NdfX;   //!
  TBranch        *b_Chi2NdfY;   //!
  TBranch        *b_SigmaXErr;   //!
  TBranch        *b_SigmaYErr;   //!
  TBranch        *b_MuMaxXErr;   //!
  TBranch        *b_MuMaxYErr;   //!
  TBranch        *b_CXErr;   //!
  TBranch        *b_CYErr;   //!

  // Set object pointer
  SigmaVis = 0;
  SigmaVisErr = 0;
  LumiSp = 0;
  LumiSpErr = 0;
  Lumi = 0;
  LumiErr = 0;
  SigmaX = 0;
  SigmaY = 0;
  MuMaxX = 0;
  MuMaxY = 0;
  CX = 0;
  CY = 0;
  RX = 0;
  RY = 0;
  Chi2NdfX = 0;
  Chi2NdfY = 0;
  SigmaXErr = 0;
  SigmaYErr = 0;
  MuMaxXErr = 0;
  MuMaxYErr = 0;
  CXErr = 0;
  CYErr = 0;

  t_vdm->SetMakeClass(1);

  t_vdm->SetBranchAddress("Scan", &Scan, &b_Scan);
  t_vdm->SetBranchAddress("BCID", &BCID, &b_BCID);
  t_vdm->SetBranchAddress("NTrkCut", &NTrkCut, &b_NTrkCut);
  t_vdm->SetBranchAddress("SigmaVis", &SigmaVis, &b_SigmaVis);
  t_vdm->SetBranchAddress("SigmaVisErr", &SigmaVisErr, &b_SigmaVisErr);
  t_vdm->SetBranchAddress("LumiSp", &LumiSp, &b_LumiSp);
  t_vdm->SetBranchAddress("LumiSpErr", &LumiSpErr, &b_LumiSpErr);
  t_vdm->SetBranchAddress("Lumi", &Lumi, &b_Lumi);
  t_vdm->SetBranchAddress("LumiErr", &LumiErr, &b_LumiErr);
  t_vdm->SetBranchAddress("SigmaX", &SigmaX, &b_SigmaX);
  t_vdm->SetBranchAddress("SigmaY", &SigmaY, &b_SigmaY);
  t_vdm->SetBranchAddress("MuMaxX", &MuMaxX, &b_MuMaxX);
  t_vdm->SetBranchAddress("MuMaxY", &MuMaxY, &b_MuMaxY);
  t_vdm->SetBranchAddress("CX", &CX, &b_CX);
  t_vdm->SetBranchAddress("CY", &CY, &b_CY);
  t_vdm->SetBranchAddress("RX", &RX, &b_RX);
  t_vdm->SetBranchAddress("RY", &RY, &b_RY);
  t_vdm->SetBranchAddress("Chi2NdfX", &Chi2NdfX, &b_Chi2NdfX);
  t_vdm->SetBranchAddress("Chi2NdfY", &Chi2NdfY, &b_Chi2NdfY);
  t_vdm->SetBranchAddress("SigmaXErr", &SigmaXErr, &b_SigmaXErr);
  t_vdm->SetBranchAddress("SigmaYErr", &SigmaYErr, &b_SigmaYErr);
  t_vdm->SetBranchAddress("MuMaxXErr", &MuMaxXErr, &b_MuMaxXErr);
  t_vdm->SetBranchAddress("MuMaxYErr", &MuMaxYErr, &b_MuMaxYErr);
  t_vdm->SetBranchAddress("CXErr", &CXErr, &b_CXErr);
  t_vdm->SetBranchAddress("CYErr", &CYErr, &b_CYErr);

  Long64_t entries = t_vdm->GetEntriesFast();

  for (int i=0; i<entries; i++) {
    t_vdm->GetEntry(i);

    TString tag_base = vtx_method;
    tag_base += NTrkCut;
    tag_base += "_BCID";
    tag_base += BCID;
    tag_base += "_scan";
    tag_base += Scan;

    for (vector<TString>::iterator fit_function = fit_functions.begin(); fit_function != fit_functions.end(); ++fit_function) {

      TString tag = tag_base;
      tag += "_";
      tag += *fit_function;
      Int_t vector_idx = distance(fit_functions.begin(), fit_function);

      m_sigma_vis[tag] = SigmaVis->at(vector_idx); /// (1. - 0.0072);
      m_sigma_vis_err[tag] = SigmaVisErr->at(vector_idx); /// (1. - 0.0072);
      m_lumi_sp[tag] = LumiSp->at(vector_idx);
      m_lumi_err[tag] = LumiErr->at(vector_idx);
      m_lumi[tag] = Lumi->at(vector_idx); //* (1. - 0.0072);
      m_lumi_sp_err[tag] = LumiSpErr->at(vector_idx); //* (1. - 0.0072);
      m_Sigma_x[tag] = SigmaX->at(vector_idx);
      m_Sigma_x_err[tag] = SigmaXErr->at(vector_idx);
      m_Sigma_y[tag] = SigmaY->at(vector_idx);
      m_Sigma_y_err[tag] = SigmaYErr->at(vector_idx);
      m_mu_max_x[tag] = MuMaxX->at(vector_idx);
      m_mu_max_x_err[tag] = MuMaxXErr->at(vector_idx);
      m_mu_max_y[tag] = MuMaxY->at(vector_idx);
      m_mu_max_y_err[tag] = MuMaxYErr->at(vector_idx);
      m_c_x[tag] = CX->at(vector_idx);
      m_c_x_err[tag] = CXErr->at(vector_idx);
      m_c_y[tag] = CY->at(vector_idx);
      m_c_y_err[tag] = CYErr->at(vector_idx);
      m_r_x[tag] = RX->at(vector_idx);
      m_r_y[tag] = RY->at(vector_idx);
      m_chi2ndf_x[tag] = Chi2NdfX->at(vector_idx);
      m_chi2ndf_y[tag] = Chi2NdfY->at(vector_idx);

    }
  }

  // -- Lumi-weighted averages
  // -- Get sigma_vis-weighted average of sigma_vis for default fit
  // -- Get lumi-weighted average of lumi_sp for default fit
  std::map<Int_t, Float_t> sigma_vis_avg, sigma_vis_err, sigma_vis_avg2, sigma_vis_rms;
  std::map<Int_t, std::map<Int_t, Float_t> > sigma_vis_scan_avg, sigma_vis_scan_err, sigma_vis_scan_avg2, sigma_vis_scan_rms;
  std::map<Int_t, Float_t> lumi_sp_avg, lumi_sp_err, lumi_sp_avg2, lumi_sp_rms;
  std::map<Int_t, std::map<Int_t, Float_t> > lumi_sp_scan_avg, lumi_sp_scan_err, lumi_sp_scan_avg2, lumi_sp_scan_rms;
  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    Float_t sv_weight = 0.;
    Float_t sv_avg_sum = 0.;
    Float_t sv_avg_sum2 = 0.;
    Float_t sv_err_sum2 = 0.;
    Float_t ls_weight = 0.;
    Float_t ls_avg_sum = 0.;
    Float_t ls_avg_sum2 = 0.;
    Float_t ls_err_sum2 = 0.;
    for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {

      Float_t sv_scan_weight = 0.;
      Float_t sv_scan_avg_sum =0.;
      Float_t sv_scan_avg_sum2 = 0.;
      Float_t sv_scan_err_sum2 = 0.;
      Float_t ls_scan_weight = 0.;
      Float_t ls_scan_avg_sum =0.;
      Float_t ls_scan_avg_sum2 = 0.;
      Float_t ls_scan_err_sum2 = 0.;

      for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {

        TString tag = vtx_method;
        tag += *nTrkCut;
        tag += "_BCID";
        tag += *bcid;
        tag += "_scan";
        tag += *scan;
        tag += "_";
        tag += default_fit;
        //Float_t current_weight = m_lumi[tag];
        Float_t sv_current_weight = m_sigma_vis_err[tag]; //Now weighting by sigma_vis error instead of luminosity
        Float_t ls_current_weight = m_lumi_sp_err[tag];
        sv_weight += sv_current_weight;
        ls_weight += ls_current_weight;
        sv_avg_sum += m_sigma_vis[tag] * sv_current_weight;
        ls_avg_sum += m_lumi_sp[tag] * ls_current_weight;
        sv_avg_sum2 += sv_current_weight * TMath::Power(m_sigma_vis[tag], 2);
        ls_avg_sum2 += ls_current_weight * TMath::Power(m_lumi_sp[tag], 2);
        sv_err_sum2 += TMath::Power(sv_current_weight * m_sigma_vis_err[tag], 2);
        ls_err_sum2 += TMath::Power(ls_current_weight * m_lumi_sp_err[tag], 2);

        sv_scan_weight += sv_current_weight;
        sv_scan_avg_sum += m_sigma_vis[tag] * sv_current_weight;
        sv_scan_avg_sum2 += sv_current_weight * TMath::Power(m_sigma_vis[tag], 2);
        sv_scan_err_sum2 += TMath::Power(sv_current_weight * m_sigma_vis_err[tag], 2);
        ls_scan_weight += ls_current_weight;
        ls_scan_avg_sum += m_lumi_sp[tag] * ls_current_weight;
        ls_scan_avg_sum2 += ls_current_weight * TMath::Power(m_lumi_sp[tag], 2);
        ls_scan_err_sum2 += TMath::Power(ls_current_weight * m_lumi_sp_err[tag], 2);

      }
      sigma_vis_scan_avg[*nTrkCut][*scan] = sv_scan_avg_sum / sv_scan_weight;
      sigma_vis_scan_err[*nTrkCut][*scan] = TMath::Sqrt(sv_scan_err_sum2) / sv_scan_weight;
      sigma_vis_scan_avg2[*nTrkCut][*scan] = sv_scan_avg_sum2 / sv_scan_weight;
      sigma_vis_scan_rms[*nTrkCut][*scan] = TMath::Sqrt(sigma_vis_scan_avg2[*nTrkCut][*scan] - TMath::Power(sigma_vis_scan_avg[*nTrkCut][*scan], 2));
      lumi_sp_scan_avg[*nTrkCut][*scan] = ls_scan_avg_sum / ls_scan_weight;
      lumi_sp_scan_err[*nTrkCut][*scan] = TMath::Sqrt(ls_scan_err_sum2) / ls_scan_weight;
      lumi_sp_scan_avg2[*nTrkCut][*scan] = ls_scan_avg_sum2 / ls_scan_weight;
      lumi_sp_scan_rms[*nTrkCut][*scan] = TMath::Sqrt(lumi_sp_scan_avg2[*nTrkCut][*scan] - TMath::Power(lumi_sp_scan_avg[*nTrkCut][*scan], 2));


    }
    sigma_vis_avg[*nTrkCut] = sv_avg_sum / sv_weight;
    sigma_vis_err[*nTrkCut] = TMath::Sqrt(sv_err_sum2) / sv_weight;
    sigma_vis_avg2[*nTrkCut] = sv_avg_sum2 / sv_weight;
    sigma_vis_rms[*nTrkCut] = TMath::Sqrt(sigma_vis_avg2[*nTrkCut] - (sigma_vis_avg[*nTrkCut] * sigma_vis_avg[*nTrkCut]));
    lumi_sp_avg[*nTrkCut] = ls_avg_sum / ls_weight;
    lumi_sp_err[*nTrkCut] = TMath::Sqrt(ls_err_sum2) / ls_weight;
    lumi_sp_avg2[*nTrkCut] = ls_avg_sum2 / ls_weight;
    lumi_sp_rms[*nTrkCut] = TMath::Sqrt(lumi_sp_avg2[*nTrkCut] - (lumi_sp_avg[*nTrkCut] * lumi_sp_avg[*nTrkCut]));
  }
	
	//Text file with sigma_vis value
	ofstream myfile;
	for (vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {
    for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    	stringstream ss;
			ss << "/afs/cern.ch/user/j/jiturbep/InfoFromVdMScans/vtx_counting/results_sigma_vis_value_scan" << *scan << "_NTrkCut" << *nTrkCut <<".txt";
			string filename = ss.str();
			myfile.open(filename.c_str());
  		myfile.precision(11);
    	myfile << sigma_vis_scan_avg[*nTrkCut][*scan] << " " << sigma_vis_scan_err[*nTrkCut][*scan];
    	myfile.close();
    	cout << "Created text file: "<< filename << endl;
  	}
	}
	
	//Text file with lumi_sp value
	ofstream myfile2;
	for (vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {
    for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    	stringstream ss;
			ss << "/afs/cern.ch/user/j/jiturbep/InfoFromVdMScans/vtx_counting/results_lumi_sp_value_scan" << *scan << "_NTrkCut" << *nTrkCut <<".txt";
			string filename = ss.str();
			myfile2.open(filename.c_str());
  		myfile2.precision(11);
    	myfile2 << lumi_sp_scan_avg[*nTrkCut][*scan] << " " << lumi_sp_scan_err[*nTrkCut][*scan];
    	myfile2.close();
    	cout << "Created text file: "<< filename << endl;
  	}
	}
	
	//Text file with fit parameters
	ofstream myfile3;
	for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {
		for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
		  stringstream ss;
			ss << "/afs/cern.ch/user/j/jiturbep/InfoFromVdMScans/vtx_counting/results_parameters_scan" << *scan << "_NTrkCut" << *nTrkCut <<".txt";
			string filename = ss.str();
			myfile3.open(filename.c_str());
  		myfile3.precision(11);
  	  for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString tag = vtx_method;
        tag += *nTrkCut;
        tag += "_BCID";
        tag += *bcid;
        tag += "_scan";
        tag += *scan;
        tag += "_";
        tag += default_fit;
        myfile3 << *bcid << " " << m_sigma_vis[tag] << " " <<  m_sigma_vis_err[tag] << " " << m_Sigma_x[tag] << " " << m_Sigma_x_err[tag] << " " << m_Sigma_y[tag] << " " << m_Sigma_y_err[tag] << " " << m_mu_max_x[tag] << " " << m_mu_max_x_err[tag] << " " << m_mu_max_y[tag] << " " << m_mu_max_y_err[tag] << " " << m_chi2ndf_x[tag] << " " << m_chi2ndf_y[tag] << " " << m_lumi_sp[tag] << " " << m_lumi_sp_err[tag] << " " << m_r_x[tag] << " " << m_r_y[tag] << " " << m_c_x[tag] << " " << m_c_x_err[tag] << " " << m_c_y[tag] << " " << m_c_y_err[tag] << endl; 
      }
     myfile3.close();
    }
  }
	
	for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {
    for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
      cout << "Scan" << *scan << " / BCID" << *bcid << endl;
      for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
        TString tag = vtx_method;
        tag += *nTrkCut;
        tag += "_BCID";
        tag += *bcid;
        tag += "_scan";
        tag += *scan;
        tag += "_";
        tag += default_fit;
        cout << " nTrkCut " << *nTrkCut << endl;
        cout.precision(11);
        cout << " sigma_vis = " << m_sigma_vis[tag] << endl;
        cout << " error = " << m_sigma_vis_err[tag];
      }
    }
  }
	
  cout.precision();
  cout << "\\begin{table}" << endl;
  cout << "\t\\scriptsize" << endl;
  cout << "\t\\begin{tabular}{|c|c|c|c|";
  //if (vtx_method == "NEvt") cout << "c|";
  cout << "}" << endl;
  cout << "\t\t\\hline" << endl;
  cout << "\t\tScan/BCID\\ \\ NTrk ";
  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    cout << "& " << *nTrkCut;
  }
  cout << " \\\\" << endl;
  cout << "\t\t\\hline" << endl;

  for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {
    for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
      cout << "\t\t" << *scan << " / " << *bcid;
      for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
        TString tag = vtx_method;
        tag += *nTrkCut;
        tag += "_BCID";
        tag += *bcid;
        tag += "_scan";
        tag += *scan;
        tag += "_";
        tag += default_fit;
        cout.precision(5);
        cout << " & " << m_sigma_vis[tag] << " $\\pm$ ";
        cout.precision(2);
        cout << m_sigma_vis_err[tag];
      }
      cout << "\\\\" << endl;
      cout << "\t\t\\hline" << endl;
    }
  }

  for (vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {

    cout << "\t\t\\textcolor{complement}{Scan " << *scan << " Avg}";
    for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
      cout.precision(5);
      cout << " & \\textcolor{complement}{" << sigma_vis_scan_avg[*nTrkCut][*scan];
      cout.precision(2);
      cout << " $\\pm$ " << sigma_vis_scan_err[*nTrkCut][*scan] << "}";
    }
    cout << "\\\\" << endl;
    cout << "\t\t\\hline" << endl;
    cout << "\t\tScan " << *scan << " RMS";
    for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
      cout.precision(2);
      cout << " & " << sigma_vis_scan_rms[*nTrkCut][*scan];
    }
    cout << "\\\\" << endl;
    cout << "\t\t\\hline" << endl;
  }

  cout << "\t\t\\textcolor{complement}{Total Avg}";
  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    cout.precision(5);
    cout << " & \\textcolor{complement}{" << sigma_vis_avg[*nTrkCut];
    cout.precision(2);
    cout << " $\\pm$ " << sigma_vis_err[*nTrkCut] << "}";
  }
  cout << "\\\\" << endl;
  cout << "\t\t\\hline" << endl;
  cout << "\t\tTotal RMS";
  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    cout.precision(3);
    cout << " & " << sigma_vis_rms[*nTrkCut];
  }
  cout << "\\\\" << endl;
  cout << "\t\t\\hline" << endl;
  cout << "\t\\end{tabular}" << endl;
  cout << "\\end{table}" << endl;

}
