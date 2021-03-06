#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <stdlib.h>

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
  fit_functions.push_back("sgcl");
  fit_functions.push_back("sgc");
  fit_functions.push_back("sg");
  fit_functions.push_back("dg");
  fit_functions.push_back("dgc");
  fit_functions.push_back("dgcl");
  fit_functions.push_back("spline");

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
  } else if (run == "207216") {
    bcids.push_back(1);
    bcids.push_back(721);
    bcids.push_back(1821);
  } else if (run == "207219") {
    bcids.push_back(1);
    bcids.push_back(721);
    bcids.push_back(1821);
  } else if (run == "214984") {
    bcids.push_back(1);
    bcids.push_back(2361);
    bcids.push_back(2881);
  } else if (run == "215021") {
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
  } else if (run == "207216") {
    scans.push_back(4);
    scans.push_back(5);
    scans.push_back(6);
  } else if (run == "207219") {
    scans.push_back(1);
  } else if (run == "214984") {
    scans.push_back(11);
    scans.push_back(14);
  } else if (run == "215021") {
    scans.push_back(1);
  }

  vector<Int_t> nTrkCuts;
  //if (vtx_method == "NEvt") nTrkCuts.push_back(2);
  nTrkCuts.push_back(2);
  nTrkCuts.push_back(3);
  nTrkCuts.push_back(4);
  nTrkCuts.push_back(5);
  nTrkCuts.push_back(7);
  nTrkCuts.push_back(10);

  // -- Read TTree into maps
  std::map<TString, Float_t> m_Sigma_x, m_Sigma_x_err, m_Sigma_y, m_Sigma_y_err, m_MuMax_x, m_MuMax_x_err, m_MuMax_y, m_MuMax_y_err, m_lumi_sp, m_lumi_sp_err;

  TTree *t_vdm = (TTree*)f_in->Get("VdmResults");
  if (!t_vdm) {
    cerr << "TTree VdMResults not found in input file! Exiting..." << endl;
    exit(1);
  }

  // Declaration of leaf types
  Int_t           Scan;
  Int_t           BCID;
  Int_t           NTrkCut;
	vector<double>  *LumiSp;
  vector<double>  *LumiSpErr;
  vector<double>  *SigmaX;
  vector<double>  *SigmaY;
  vector<double>  *SigmaXErr;
  vector<double>  *SigmaYErr;
  vector<double>  *MuMaxX;
  vector<double>  *MuMaxY;
  vector<double>  *MuMaxXErr;
  vector<double>  *MuMaxYErr;

  // List of branches
  TBranch        *b_Scan;   //!
  TBranch        *b_BCID;   //!
  TBranch        *b_NTrkCut;   //!
  TBranch        *b_LumiSp;   //!
  TBranch        *b_LumiSpErr;   //!
  TBranch        *b_SigmaX;   //!
  TBranch        *b_SigmaY;   //!
  TBranch        *b_SigmaXErr;   //!
  TBranch        *b_SigmaYErr;   //!
  TBranch        *b_MuMaxX;   //!
  TBranch        *b_MuMaxY;   //!
  TBranch        *b_MuMaxXErr;   //!
  TBranch        *b_MuMaxYErr;   //!

  // Set object pointer
  SigmaX = 0;
  SigmaY = 0;
  MuMaxX = 0;
  MuMaxY = 0;
  LumiSp = 0;
  LumiSpErr = 0;
  SigmaXErr = 0;
  SigmaYErr = 0;
  MuMaxXErr = 0;
  MuMaxYErr = 0;

  t_vdm->SetMakeClass(1);

  t_vdm->SetBranchAddress("Scan", &Scan, &b_Scan);
  t_vdm->SetBranchAddress("BCID", &BCID, &b_BCID);
  t_vdm->SetBranchAddress("NTrkCut", &NTrkCut, &b_NTrkCut);
  t_vdm->SetBranchAddress("LumiSp", &LumiSp, &b_LumiSp);
  t_vdm->SetBranchAddress("LumiSpErr", &LumiSpErr, &b_LumiSpErr);
  t_vdm->SetBranchAddress("SigmaX", &SigmaX, &b_SigmaX);
  t_vdm->SetBranchAddress("SigmaY", &SigmaY, &b_SigmaY);
  t_vdm->SetBranchAddress("MuMaxX", &MuMaxX, &b_MuMaxX);
  t_vdm->SetBranchAddress("MuMaxY", &MuMaxY, &b_MuMaxY);
  t_vdm->SetBranchAddress("SigmaXErr", &SigmaXErr, &b_SigmaXErr);
  t_vdm->SetBranchAddress("SigmaYErr", &SigmaYErr, &b_SigmaYErr);
  t_vdm->SetBranchAddress("MuMaxXErr", &MuMaxXErr, &b_MuMaxXErr);
  t_vdm->SetBranchAddress("MuMaxYErr", &MuMaxYErr, &b_MuMaxYErr);

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

      
      m_lumi_sp[tag] = LumiSp->at(vector_idx);
      m_lumi_sp_err[tag] = LumiSpErr->at(vector_idx); //* (1. - 0.0072);
	  m_Sigma_x[tag] = SigmaX->at(vector_idx);
      m_Sigma_x_err[tag] = SigmaXErr->at(vector_idx);
      m_Sigma_y[tag] = SigmaY->at(vector_idx);
      m_Sigma_y_err[tag] = SigmaYErr->at(vector_idx);
      m_MuMax_x[tag] = MuMaxX->at(vector_idx);
      m_MuMax_x_err[tag] = MuMaxXErr->at(vector_idx);
      m_MuMax_y[tag] = MuMaxY->at(vector_idx);
      m_MuMax_y_err[tag] = MuMaxYErr->at(vector_idx);

    }
  }

  cout.precision(6);
  cout << "\\begin{table}" << endl;
  cout << "\t\\scriptsize" << endl;
  cout << "\t\\begin{tabular}{|c|c|c|c|c|c|";
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
      cout << "\t\t /" << *scan << " / " << *bcid;
      for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
        TString tag = vtx_method;
        tag += *nTrkCut;
        tag += "_BCID";
        tag += *bcid;
        tag += "_scan";
        tag += *scan;
        tag += "_";
        tag += default_fit;
        cout.precision(6);
        cout << " & " << m_Sigma_x[tag] << " $\\pm$ ";
        cout.precision(6);
        cout << m_Sigma_x_err[tag];
      }
      cout << "\\\\" << endl;
      cout << "\t\t\\hline" << endl;
    }
  }

for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {
    for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
      cout << "\t\t " << *scan << " / " << *bcid;
      for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
        TString tag = vtx_method;
        tag += *nTrkCut;
        tag += "_BCID";
        tag += *bcid;
        tag += "_scan";
        tag += *scan;
        tag += "_";
        tag += default_fit;
        cout.precision(6);
        cout << " & " << m_Sigma_y[tag] << " $\\pm$ ";
        cout.precision(6);
        cout << m_Sigma_y_err[tag];
      }
      cout << "\\\\" << endl;
      cout << "\t\t\\hline" << endl;
    }
  }
  
    for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {
    for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
      cout << "\t\t /" << *scan << " / " << *bcid;
      for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
        TString tag = vtx_method;
        tag += *nTrkCut;
        tag += "_BCID";
        tag += *bcid;
        tag += "_scan";
        tag += *scan;
        tag += "_";
        tag += default_fit;
        cout.precision(6);
        cout << " & " << m_MuMax_x[tag] << " $\\pm$ ";
        cout.precision(4);
        cout << m_MuMax_x_err[tag];
      }
      cout << "\\\\" << endl;
      cout << "\t\t\\hline" << endl;
    }
  }

for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {
    for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
      cout << "\t\t " << *scan << " / " << *bcid;
      for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
        TString tag = vtx_method;
        tag += *nTrkCut;
        tag += "_BCID";
        tag += *bcid;
        tag += "_scan";
        tag += *scan;
        tag += "_";
        tag += default_fit;
        cout.precision(6);
        cout << " & " << m_MuMax_y[tag] << " $\\pm$ ";
        cout.precision(4);
        cout << m_MuMax_y_err[tag];
      }
      cout << "\\\\" << endl;
      cout << "\t\t\\hline" << endl;
    }
  }

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
        cout.precision(6);
        cout << " & " << m_lumi_sp[tag] << " $\\pm$ ";
        cout.precision(4);
        cout << m_lumi_sp_err[tag];
      }
      cout << "\\\\" << endl;
      cout << "\t\t\\hline" << endl;
    }
  }

  cout << "\\\\" << endl;
  cout << "\t\t\\hline" << endl;
  cout << "\t\\end{tabular}" << endl;
  cout << "\\end{table}" << endl;

cout << "Sigma_x = " << endl;
for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {
  cout << "Scan " << *scan << endl;
    for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString tag = vtx_method;
        tag += 5;
        tag += "_BCID";
        tag += *bcid;
        tag += "_scan";
        tag += *scan;
        tag += "_";
        tag += default_fit;
        cout.precision(9); 
        cout << m_Sigma_x[tag] << endl;       
    }
  }
  cout << "Sigma_y = " << endl;
  for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {
    cout << "Scan " << *scan << endl;
    for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString tag = vtx_method;
        tag += 5;
        tag += "_BCID";
        tag += *bcid;
        tag += "_scan";
        tag += *scan;
        tag += "_";
        tag += default_fit;
        cout.precision(9); 
        cout << m_Sigma_y[tag] << endl;       
    }
  }
  cout << "Sigma_x_err = " << endl;
  for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {
    cout << "Scan " << *scan << endl;
    for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString tag = vtx_method;
        tag += 5;
        tag += "_BCID";
        tag += *bcid;
        tag += "_scan";
        tag += *scan;
        tag += "_";
        tag += default_fit;
        cout.precision(6); 
        cout << m_Sigma_x_err[tag] << endl;       
    }
  }
  cout << "Sigma_y_err = " << endl; 
  for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {
    cout << "Scan " << *scan << endl;
    for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString tag = vtx_method;
        tag += 5;
        tag += "_BCID";
        tag += *bcid;
        tag += "_scan";
        tag += *scan;
        tag += "_";
        tag += default_fit;
        cout.precision(6);
        cout << m_Sigma_y_err[tag] << endl;       
    }
  }



}
