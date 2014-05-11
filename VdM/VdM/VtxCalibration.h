/**
  *  Class VtxCalibration
  *  Container for vertexing visible cross sections
  *  Author: David Yu (dryu@lbl.gov)
  */

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

using namespace std;

class VtxCalibration {

  public:
    VtxCalibration(TString p_path);

    ~VtxCalibration();

    void Initialize();
    void PrintLatexTable();

    double GetCrossSection(Float_t p_energy, Int_t p_ntrkcut, TString p_settings, TString p_vtx_method);
    double GetCrossSectionError(Float_t p_energy, Int_t p_ntrkcut, TString p_settings, TString p_vtx_method);
  private:
    // Map is energy : ntrkcut : settings : vertex method
    std::map<Float_t, std::map<Int_t, std::map<TString, std::map<TString, double> > > > xsec, xsec_err;
    std::map<TString, Float_t> run_to_energy;
    std::map<TString, Int_t> run_to_default_fit;
    std::map<TString, vector<Int_t> > run_to_bcids;

  public:
    TString path_base;
    //Make public to allow user-settings
    std::vector<TString> runs, reco_settings, vtx_settings;
    vector<Int_t> nTrkCuts;


};
