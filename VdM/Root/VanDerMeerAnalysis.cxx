#define RANDOM_TRIGGER
#include "VdM/VanDerMeerAnalysis.h"

#include "TKey.h"
#include "TCollection.h"
#include "TDirectory.h"

using namespace std;

Double_t fit_func_sg(Double_t *x, Double_t *par);
Double_t fit_func_sgc(Double_t *x, Double_t *par);
Double_t fit_func_sgcl(Double_t *x, Double_t *par);
Double_t fit_func_dg(Double_t *x, Double_t *par);
Double_t fit_func_dgc(Double_t *x, Double_t *par);
Double_t fit_func_dgcl(Double_t *x, Double_t *par);
Double_t fit_func_spline(Double_t *x, Double_t *par);

VanDerMeerAnalysis::VanDerMeerAnalysis(TString p_vtx_method, Int_t p_ntrkcut, TString p_save_tag) {

  vtx_method = p_vtx_method;
  ntrkcut = p_ntrkcut;
  //prescale = 1;
  save_tag = p_save_tag;

  if (p_vtx_method != "Evt" && p_vtx_method != "Vtx" && p_vtx_method != "Unf") {
    cerr << "[VanDerMeerAnalysis] ERROR : Must specify argt 1 = Vtx or Evt or Unf. Exiting..." << endl;
    exit(1);
  }

  if (p_ntrkcut != 2 && p_ntrkcut != 5 && p_ntrkcut != 7 && p_ntrkcut != 10 && p_ntrkcut != 3 && p_ntrkcut != 4 && p_ntrkcut != 6 && p_ntrkcut != 8 ) {
    cerr << "[VanDerMeerAnalysis] ERROR : Invalid NTrkCut specified in inialization: " << p_ntrkcut << ". Exiting..." << endl;
    exit(1);
  }

  /*if (p_ntrkcut == 2 && (p_vtx_method == "Vtx" || p_vtx_method == "Unf")) {
    cerr << "[VanDerMeerAnalysis] ERROR: NTrk=2 not supported in Vtx and Unf methods." << endl;
    exit(1);
  }*/

  init_vertex_counts = false;
  init_plb_timestamps["x"] = false;
  init_plb_timestamps["y"] = false;
  init_deadtime = false;
  init_bunch_intensities = false;
  init_fc = false;
  init_pmc = false;

  low_displacement = 10000000.;
  high_displacement = -10000000.;

  myRMS_x = 0;
  myRMS_y = 0;
  mean_x = 0;
  mean_y = 0;
  max_x = 0;
  max_y = 0;
  sumw_x = 0;
  sumw_y = 0;
  sumdev2_x = 0;
  sumdev2_y = 0;

  systematic_uncertainty_list["fake_high"] = false; // Highest fake mu from mu scaling
  systematic_uncertainty_list["fake_low"] = false; // Lowest fake mu from mu scaling
  systematic_uncertainty_list["masking_toy_scaling"] = false; // Scale p_mask up by 2%, as suggested by the simulation.
  systematic_uncertainty_list["bs-45mm"] = false;
  systematic_uncertainty_list["bs-55mm"] = false;

  TH1::SetDefaultSumw2(); //in case it's not already set

	#ifdef USE_UNFOLD
  m_unfold = 0;
	#endif
  fc = 0;
  tg_musp_x = 0;
  tg_musp_y = 0;
  tg_musp_x_ini = 0;
  tg_musp_y_ini = 0;
  fit_spline_x = 0;
  fit_spline_y = 0;
  fit_cov_x_dg = 0;
  fit_cov_y_dg = 0;

}

VanDerMeerAnalysis::~VanDerMeerAnalysis() {
#ifdef USE_UNFOLD
  if (m_unfold) {
    delete m_unfold;
  }
#endif
  if (fc) {
    delete fc;
  }
  if (tg_musp_x) {
    delete tg_musp_x;
  }
  if (tg_musp_x_ini) {
    delete tg_musp_x_ini;
  }
  if (tg_musp_y) {
    delete tg_musp_y;
  }
  if (tg_musp_y_ini) {
    delete tg_musp_y_ini;
  }
  if (fit_spline_x) {
    delete fit_spline_x;
  }
  if (fit_spline_y) {
    delete fit_spline_y;
  }
  if (fit_cov_x_dg) {
    delete fit_cov_x_dg;
  }
  if (fit_cov_y_dg) {
    delete fit_cov_y_dg;
  }


  //now clear maps
  for (std::map<TString, TF1*>::iterator itf = f_x.begin();
       itf != f_x.end(); ++ itf) {
    if (itf->second) {
      delete itf->second;
    }
  }
  f_x.clear();

  for (std::map<TString, TF1*>::iterator itf = f_y.begin();
       itf != f_y.end(); ++ itf) {
    if (itf->second) {
      delete itf->second;
    }
  }
  f_y.clear();

  for (std::map<TString, TMatrixDSym*>::iterator itM = fit_cov_x.begin();
       itM != fit_cov_x.end(); ++itM) {
    if (itM->second) {
      delete itM->second;
    }
  }
  fit_cov_x.clear();

  for (std::map<TString, TMatrixDSym*>::iterator itM = fit_cov_y.begin();
       itM != fit_cov_y.end(); ++itM) {
    if (itM->second) {
      delete itM->second;
    }
  }
  fit_cov_y.clear();

  for (std::map<TString, TMatrixD*>::iterator itM = v_dsigma_dpx.begin();
       itM != v_dsigma_dpx.end(); ++itM) {
    if (itM->second) {
      delete itM->second;
    }
  }
  v_dsigma_dpx.clear();

  for (std::map<TString, TMatrixD*>::iterator itM = v_dsigma_dpy.begin();
       itM != v_dsigma_dpy.end(); ++itM) {
    if (itM->second) {
      delete itM->second;
    }
  }
  v_dsigma_dpy.clear();

  for (std::map<TString, TMatrixD*>::iterator itM = v_dsigma_dpx_T.begin();
       itM != v_dsigma_dpx_T.end(); ++itM) {
    if (itM->second) {
      delete itM->second;
    }
  }
  v_dsigma_dpx_T.clear();

  for (std::map<TString, TMatrixD*>::iterator itM = v_dsigma_dpy_T.begin();
       itM != v_dsigma_dpy_T.end(); ++itM) {
    if (itM->second) {
      delete itM->second;
    }
  }
  v_dsigma_dpy_T.clear();

}

void VanDerMeerAnalysis::SetScanEndpoints(Float_t p_low, Float_t p_high) {

  low_displacement = p_low;
  high_displacement = p_high;
}

void VanDerMeerAnalysis::LoadVertexCounts(TString p_path, Int_t p_bcid, TString p_trigger_type) {

  cout << "[VanDerMeerAnalysis] INFO : Loading vertex counts from " << p_path << endl;

  trigger_type = p_trigger_type;
  cout << "[VDMA] Line 189 trigger_type = " << trigger_type << endl;

  if (vtx_method == "Unf") {
    #ifdef USE_UNFOLD
    LoadVertexCountsUnfold(p_path, p_bcid);
    #else
    cerr << "Unfolding not supported. Please compile with USE_UNFOLD in cmt/Makefile.RootCore" << endl;
    exit(1);
    #endif
  } else {
    //load from tree file
    TFile *f_in = new TFile(p_path, "READ");
    TTree *t_vdm = (TTree*)f_in->Get("t_vdm");

    Long64_t n_entries = t_vdm->GetEntriesFast();

    stringstream ss_vtx_leaf;
    ss_vtx_leaf << "N" << vtx_method << ntrkcut << "_BSC";
    TString vtx_leaf(ss_vtx_leaf.str());

    for (Int_t entry = 0; entry < n_entries; entry++) {

      t_vdm->GetEntry(entry);
      Int_t current_bcid = int(t_vdm->GetLeaf("BunchId")->GetValue(0));
      if (current_bcid != p_bcid) {
        continue;
      }

      Int_t current_pLB = int(t_vdm->GetLeaf("pLB")->GetValue(0));
      //t_vdm->Print();
      //std::cout << "vtx_leaf: " << vtx_leaf << std::endl;
      Float_t current_nvtx = t_vdm->GetLeaf(vtx_leaf)->GetValue(0);
      nvtx_pLB[current_pLB] = current_nvtx * prescale[current_pLB];
      //cout << "Prescale[" << current_pLB << "] = " << prescale[current_pLB] << endl;

      if (vtx_method == "Vtx") {
        nvtx_err_pLB[current_pLB] = (0.5 + TMath::Sqrt(current_nvtx + 0.25)) * prescale[current_pLB];
      } else if (vtx_method == "Unf") {
        //load error from tree -- not used anymore
        TLeaf *l_unf_err = t_vdm->GetLeaf(vtx_leaf+"_Err");
        if (!l_unf_err) {
          cerr << "WARNING: Cannot load errors from t_vdm tree for Unf method." << endl;
        } else {
          //This is not a per-BCID prescale
          nvtx_err_pLB[current_pLB] = l_unf_err->GetValue(0) * prescale[current_pLB];
        }
      }

      if (p_trigger_type == "BGRP7") {

        Float_t current_ntrig = t_vdm->GetLeaf("NTrig")->GetValue(0);
        Float_t current_duration = pLB_timestamps[current_pLB].second - pLB_timestamps[current_pLB].first;

        nvtx_pLB[current_pLB] *= 11245.5 * current_duration / current_ntrig;
        if (vtx_method == "Vtx") {
          nvtx_err_pLB[current_pLB] *= 11245.5 * current_duration / current_ntrig;
        }
      }

    }
    f_in->Close();
  }

  #ifdef DEBUG_VDM
  cout << "[VanDerMeerAnalysis] DEBUG : Vertex counts" << endl;
  for (vector<Int_t>::iterator it = plb_list_x.begin(); it != plb_list_x.end(); it++) {
    cout << "[VanDerMeerAnalysis] DEBUG : NVtx[" << *it << "] = " << nvtx_pLB[*it] << endl;
  }
  for (vector<Int_t>::iterator it = plb_list_y.begin(); it != plb_list_y.end(); it++) {
    cout << "[VanDerMeerAnalysis] DEBUG : NVtx[" << *it << "] = " << nvtx_pLB[*it] << endl;
  }
  #endif

  init_vertex_counts = true;

}


void VanDerMeerAnalysis::LoadPlbTimestamps(TString p_path, TString axis) {

  cout << "[VanDerMeerAnalysis] INFO : Loading pLB timestamps for " << axis << " scan, from " << p_path << endl;

  ifstream tsfile(p_path.Data());
  if (!tsfile.is_open()) {
    cerr << "[VanDerMeerAnalysis] ERROR : Unable to open requested timestamp file: " << p_path << ". Exiting..." << endl;
    exit(1);
  }

  Int_t current_pLB;
  Double_t ts_start, ts_end;

  bool first = true;
  Int_t first_pLB;

  while (tsfile.good()) {
    string line;
    getline(tsfile, line);
    if (line.empty()) {
      continue;
    }
    stringstream istr(line);
    if (line.find_first_of("PseudoLB") != string::npos) {
      continue;  // header
    }

    istr >> current_pLB >> ts_start >> ts_end;
    std::pair<Double_t, Double_t> ts = make_pair<Double_t, Double_t>(ts_start, ts_end);

    if (first) {
      first = false;
      first_pLB = current_pLB;
    } else {
      if ((current_pLB - first_pLB) % 2 != 0) {
        continue;
      }
    }
    plb_list.push_back(current_pLB);
    pLB_timestamps[current_pLB] = ts;
    if (axis == "x") {
      pLB_timestamps_x[current_pLB] = ts;
      plb_list_x.push_back(current_pLB);
    } else {
      pLB_timestamps_y[current_pLB] = ts;
      plb_list_y.push_back(current_pLB);
    }
  }

  tsfile.close();

  #ifdef DEBUG_VDM
  cout << "pLB timestamps:" << endl;
  cout.precision(10);
  for (vector<Int_t>::iterator it = plb_list_x.begin(); it != plb_list_x.end(); it++) {
    cout << *it << " : " << pLB_timestamps_x[*it].first << " - " << pLB_timestamps_x[*it].second << endl;
  }
  for (vector<Int_t>::iterator it = plb_list_y.begin(); it != plb_list_y.end(); it++) {
    cout << *it << " : " << pLB_timestamps_y[*it].first << " - " << pLB_timestamps_y[*it].second << endl;
  }
  #endif


  init_plb_timestamps[axis] = true;
}

void VanDerMeerAnalysis::LoadPlbPrescales(TString p_path) {

  cout << "[VanDerMeerAnalysis] INFO : Loading prescales from " << p_path << endl;

  Int_t current_pLB;
  Double_t current_prescale;

  if (p_path != "") {
    TFile *f_dt = new TFile(p_path, "READ");
    TH1F *h_dt = (TH1F*)f_dt->Get("PlbPrescale");
    for (int i = 1; i <= h_dt->GetNbinsX(); i++) {
      current_pLB = TMath::Nint(h_dt->GetBinCenter(i));
      current_prescale = h_dt->GetBinContent(i);
      prescale[current_pLB] = current_prescale;
    }

    f_dt->Close();
  } else {
    for (vector<Int_t>::iterator it = plb_list.begin(); it != plb_list.end(); ++it) {
      prescale[*it] = 1.;
    }
  }

#ifdef DEBUG_VDM
  cout << "[VanDerMeerAnalysis] DEBUG : Prescales" << endl;
  for (vector<Int_t>::iterator it = plb_list_x.begin(); it != plb_list_x.end(); it++) {
    cout << "[VanDerMeerAnalysis] DEBUG : f[" << *it << "] = " << prescale[*it] << endl;
  }
  for (vector<Int_t>::iterator it = plb_list_y.begin(); it != plb_list_y.end(); it++) {
    cout << "[VanDerMeerAnalysis] DEBUG : f[" << *it << "] = " << prescale[*it] << endl;
  }

#endif

  //init_deadtime = true;
}


void VanDerMeerAnalysis::LoadDeadtime(TString p_path, TString p_path_tmp) {

  cout << "[VanDerMeerAnalysis] INFO : Loading deadtime from " << p_path << endl;

  Int_t current_pLB;
  Double_t current_livefraction;

  if (p_path != "") {
    TFile *f_dt = new TFile(p_path, "READ");
    TH1F *h_dt = (TH1F*)f_dt->Get("PlbLiveFractions");
    for (int i = 1; i <= h_dt->GetNbinsX(); i++) {
      current_pLB = TMath::Nint(h_dt->GetBinCenter(i));
      current_livefraction = h_dt->GetBinContent(i);
      if (current_livefraction <= 1.) {
        live_fractions[current_pLB] = current_livefraction;
      } else {
        live_fractions[current_pLB] = 1.;
      }

    }
    
    if (p_path_tmp != "") {
    TFile *f_dt_2 = new TFile(p_path_tmp, "READ");
    TH1F *h_dt_2 = (TH1F*)f_dt_2->Get("PlbLiveFractions");
    for (int i = 1; i <= h_dt_2->GetNbinsX(); i++) {
      current_pLB = TMath::Nint(h_dt_2->GetBinCenter(i));
      current_livefraction = h_dt_2->GetBinContent(i);
      if (current_livefraction <= 1.) {
        live_fractions[current_pLB] = current_livefraction;
      } else {
        live_fractions[current_pLB] = 1.;
      }

    	}
    }

    f_dt->Close();
  } else {
    for (vector<Int_t>::iterator it = plb_list.begin(); it != plb_list.end(); ++it) {
      live_fractions[*it] = 1.;
    }
  }

#ifdef DEBUG_VDM
  cout << "[VanDerMeerAnalysis] DEBUG : Live Fractions" << endl;
  for (vector<Int_t>::iterator it = plb_list_x.begin(); it != plb_list_x.end(); it++) {
    cout << "[VanDerMeerAnalysis] DEBUG : f[" << *it << "] = " << live_fractions[*it] << endl;
  }
  for (vector<Int_t>::iterator it = plb_list_y.begin(); it != plb_list_y.end(); it++) {
    cout << "[VanDerMeerAnalysis] DEBUG : f[" << *it << "] = " << live_fractions[*it] << endl;
  }

#endif

  init_deadtime = true;
}

void VanDerMeerAnalysis::LoadBunchIntensities(TString p_path, Int_t bcid) {
  cout << "[VanDerMeerAnalysis] INFO : Loading bunch intensities from " << p_path << endl;

  TFile *f_scan = new TFile(p_path, "READ");
  TTree *t_scan = (TTree*)f_scan->Get("vdMScanData");

  t_scan->SetBranchStatus("*", 0);
  t_scan->SetBranchStatus("SCANDATA", 1);
  t_scan->SetBranchStatus("BCT_B1BunchIntensity", 1);
  t_scan->SetBranchStatus("BCT_B2BunchIntensity", 1);
  t_scan->SetBranchStatus("BCT_B1BCID", 1);
  t_scan->SetBranchStatus("BCT_B2BCID", 1);

  Int_t bcid1_tmp[2400];
  Int_t bcid2_tmp[2400];
  Float_t n1_tmp[2400];
  Float_t n2_tmp[2400];
  UInt_t current_pLB;

  t_scan->SetBranchAddress("BCT_B1BunchIntensity", &n1_tmp);
  t_scan->SetBranchAddress("BCT_B2BunchIntensity", &n2_tmp);
  t_scan->SetBranchAddress("BCT_B1BCID", &bcid1_tmp);
  t_scan->SetBranchAddress("BCT_B2BCID", &bcid2_tmp);

  Long64_t n_entries = t_scan->GetEntriesFast();

  for (int i = 0; i < n_entries; i++) {
    t_scan->GetEntry(i);

    // pLB
    if (p_path == "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/LumiVtxExtras/CoolScanNtuple//2012NovScan1Raw-v13.root") {
        current_pLB = int( t_scan->GetLeaf("ScanLB")->GetValue(0) ) + 1;
    }
    else{current_pLB = int( t_scan->GetLeaf("ScanLB")->GetValue(0) );}

    // Find index of BCID
    int beam1_index = -1;
    int beam2_index = -1;
    for (int j = 0; j < 2400; j++) {
      if (bcid1_tmp[j] == bcid) {
        beam1_index = j;
        break;
      }
    }
    for (int j = 0; j < 2400; j++) {
      if (bcid2_tmp[j] == bcid) {
        beam2_index = j;
        break;
      }
    }
    if ((beam1_index >= 0) && (beam2_index >= 0)) {
 	  //Double_t ghost_charge_scaling = (1 - 0.0072);
 	  Double_t ghost_charge_scaling = 1;
 	 // cout << "[VanDerMeerAnalysis] Applying ghost charge scaling, changed to be 1" << endl;
      bunch_intensities_1[current_pLB] = n1_tmp[beam1_index] * TMath::Power(10., 11) * TMath::Sqrt(ghost_charge_scaling);
      bunch_intensities_2[current_pLB] = n2_tmp[beam2_index] * TMath::Power(10., 11) * TMath::Sqrt(ghost_charge_scaling);			 
    } else {
      cerr << "[VanDerMeerAnalysis] ERROR : Did not find BCID " << bcid << " in Eric's scan ntuple. Exiting..." << endl;
      exit(1);
    }
  }

#ifdef DEBUG_VDM
  cout << "[VanDerMeerAnalysis] DEBUG : Bunch Intensities" << endl;
  for (std::vector<Int_t>::iterator current_pLB = plb_list_x.begin(); current_pLB != plb_list_x.end(); current_pLB++) {
    cout << "[VanDerMeerAnalysis] DEBUG : N1[" << *current_pLB << "] = " << bunch_intensities_1[*current_pLB] << endl;
    cout << "[VanDerMeerAnalysis] DEBUG : N2[" << *current_pLB << "] = " << bunch_intensities_2[*current_pLB] << endl;
  }

#endif
  init_bunch_intensities = true;

  f_scan->Close();
}

void VanDerMeerAnalysis::LoadNominalSeparations(TString p_path) {

  cout << "[VanDerMeerAnalysis] INFO : Loading nominal separations from " << p_path << endl;

  TFile *f_scan = new TFile(p_path, "READ");
  TTree *t_scan = (TTree*)f_scan->Get("vdMScanData");

  t_scan->SetBranchStatus("*", 0);
  t_scan->SetBranchStatus("SCANDATA", 1);

  Int_t current_pLB;
  //Float_t current_NominalSeparation;

  Long64_t n_entries = t_scan->GetEntriesFast();

  for (int i = 0; i < n_entries; i++) {
    t_scan->GetEntry(i);
    if (t_scan->GetLeaf("AcquisitionFlag")->GetValue(0) != 1) {
      continue;
    }
    if (p_path == "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/LumiVtxExtras/CoolScanNtuple//2012NovScan1Raw-v13.root") {
        current_pLB = int( t_scan->GetLeaf("ScanLB")->GetValue(0) ) + 1;
    }
    else {current_pLB = int( t_scan->GetLeaf("ScanLB")->GetValue(0) );}

    nominal_separation[current_pLB] = t_scan->GetLeaf("NominalSeparation")->GetValue(0) * 1000.;
    cout << "[VDMA] pLB " << current_pLB << " has separation " << nominal_separation[current_pLB] << endl;

    if (nominal_separation[current_pLB] > high_displacement) {
      high_displacement = nominal_separation[current_pLB];
      //cout << "[VDMA] Line 473 high_displacement = " << high_displacement << endl;
    }
    if (nominal_separation[current_pLB] < low_displacement) {
      low_displacement = nominal_separation[current_pLB];
      //cout << "[VDMA] Line 477 low_displacement = " << low_displacement << endl;
    }
  }

  f_scan->Close();

}

void VanDerMeerAnalysis::InitializeFakeCorrection(TString p_energy, TString p_settings, Int_t p_ntrkcut) {

  cout << "[VanDerMeerAnalysis] INFO : Initializing fake correction" << endl;

  // Old fake correction : fc = new FakeCorrection(p_energy, p_settings, p_ntrkcut);
  if (systematic_uncertainty_list["bs-45mm"]) {
    fc = new FakeCorrection("mc_7TeV_17.2_normal_pythia8_pu_bs45", p_ntrkcut);
  } else if (systematic_uncertainty_list["bs-55mm"]) {
    fc = new FakeCorrection("mc_7TeV_17.2_normal_pythia8_pu_bs55", p_ntrkcut); 
    //fc = new FakeCorrection("mc_7TeV_17.2_VtxLumi_pythia8_pu", p_ntrkcut); 
  } else if (p_settings == "17.2-normal") {
    fc = new FakeCorrection("mc_8TeV_17.2_normal_pythia8_pu", p_ntrkcut); //CHANGED this to 7, it was 8 //changed back to 8
  } else if (p_settings = "17.2-VtxLumi") {
    //fc = new FakeCorrection("mc_8TeV_17.2_VtxLumi_2newsets", p_ntrkcut); //CHANGED this to 7, it was 8 // changed back to 8
    fc = new FakeCorrection("mc_8TeV_17.2_VtxLumi_BothSamples", p_ntrkcut);
  } else {
    cerr << "Fake correction only available for r17.2 (you specified " << p_settings << "). Exiting..." << endl;
    exit(1);
  }
  init_fc = true;
}

void VanDerMeerAnalysis::CalculateMuPlb() {

  // Vertex counting method (both reco or unfolded)
  if (vtx_method == "Vtx") {

    for (vector<Int_t>::iterator it = plb_list_x.begin(); it != plb_list_x.end(); ++it) {

      Double_t current_duration = pLB_timestamps_x[*it].second - pLB_timestamps_x[*it].first;
      mu_pLB[*it] = nvtx_pLB[*it] / (current_duration * 11245.5);
      mu_err_pLB[*it] = nvtx_err_pLB[*it] / (current_duration * 11245.5);
      //cout << "[VDMA] Line 516, for x mu_pLB["<<*it<<"] = " << mu_pLB[*it] << endl;

      if (init_deadtime) {
        mu_pLB[*it] = mu_pLB[*it] / live_fractions[*it];
        mu_err_pLB[*it] = mu_err_pLB[*it] / live_fractions[*it];
        //cout << "[VDMA] Line 529, for x, deadtime, mu_pLB["<<*it<<"] = " << mu_pLB[*it] << endl;
        //cout << "[VDMA] Line 530, for x, deadtime, live_fractions["<<*it<<"] = " << live_fractions[*it] << endl;
      }

    }

    for (vector<Int_t>::iterator it = plb_list_y.begin(); it != plb_list_y.end(); ++it) {

      Double_t current_duration = pLB_timestamps_y[*it].second - pLB_timestamps_y[*it].first;
      mu_pLB[*it] = nvtx_pLB[*it] / (current_duration * 11245.5);
      mu_err_pLB[*it] = nvtx_err_pLB[*it] / (current_duration * 11245.5);
      //cout << "[VDMA] Line 531, for y mu_pLB["<<*it<<"] = " << mu_pLB[*it] << endl;

      if (init_deadtime) {
        mu_pLB[*it] = mu_pLB[*it] / live_fractions[*it];
        mu_err_pLB[*it] = mu_err_pLB[*it] / live_fractions[*it];
        //cout << "[VDMA] Line 545, for y, deadtime, mu_pLB["<<*it<<"] = " << mu_pLB[*it] << endl;
        //cout << "[VDMA] Line 546, for y, deadtime, live_fractions["<<*it<<"] = " << live_fractions[*it] << endl;
      }

    }


  } else if (vtx_method == "Evt") {
    for (vector<Int_t>::iterator it = plb_list_x.begin(); it != plb_list_x.end(); ++it) {

      Double_t current_duration = pLB_timestamps_x[*it].second - pLB_timestamps_x[*it].first;
      Double_t current_nevt = nvtx_pLB[*it];
      if (init_deadtime) {
        current_nevt = current_nevt / live_fractions[*it];
      }
      Double_t p = current_nevt / (current_duration* 11245.5);

      //Double_t p_fake = fc->GetFakeEventProbabilityFromMuReconMC(p);

      mu_pLB[*it] = -1. * TMath::Log(1. - (current_nevt / current_duration) / 11245.5);

      // -- Statistical uncertainty
      if (p > 1.) {
        cout << "[VanDerMeerAnalysis] WARNING : p = " << p << " exceeds 1!" << endl;
        mu_err_pLB[*it] = 100.;
      } else {
        Double_t dp = TMath::Sqrt(p * (1. - p) / (current_duration * 11245.5 * live_fractions[*it]));
        mu_err_pLB[*it] = dp / (1. - p);
      }
    }

    for (vector<Int_t>::iterator it = plb_list_y.begin(); it != plb_list_y.end(); ++it) {

      Double_t current_duration = pLB_timestamps_y[*it].second - pLB_timestamps_y[*it].first;
      Double_t current_nevt = nvtx_pLB[*it];
      if (init_deadtime) {
        current_nevt = current_nevt / live_fractions[*it];
      }

      mu_pLB[*it] = -1. * TMath::Log(1. - (current_nevt / current_duration) / 11245.5);

      // -- Statistical uncertainty
      Double_t p = nvtx_pLB[*it] / (current_duration * 11245.5 * live_fractions[*it]);
      if (p > 1.) {
        cout << "[VanDerMeerAnalysis] WARNING : p = " << p << " exceeds 1!" << endl;
        mu_err_pLB[*it] = 100.;
      } else {
        Double_t dp = TMath::Sqrt(p * (1. - p) / (current_duration * 11245.5 * live_fractions[*it]));
        mu_err_pLB[*it] = dp / (1. - p);
      }
    }

  } else if (vtx_method == "Unf") {
    //numbers are already mu values
    vector<Int_t> allpLB;
    allpLB.insert(allpLB.end(), plb_list_x.begin(), plb_list_x.end());
    allpLB.insert(allpLB.end(), plb_list_y.begin(), plb_list_y.end());
    for (vector<Int_t>::iterator it = allpLB.begin(); it != allpLB.end(); ++it) {
      mu_pLB[*it] = nvtx_pLB[*it];
      mu_err_pLB[*it] = nvtx_err_pLB[*it];
    }
  }

}

void VanDerMeerAnalysis::SetUnfResponseMatrix(TString fileName, TString histoBaseName) {
  responseMatrixFileName = fileName;
  responseMatrixHName = histoBaseName;
  cout << "[VanDerMeerAnalysis] Selected Response Matrix " << responseMatrixHName
       << " from file " << responseMatrixFileName << endl;
}

#ifdef USE_UNFOLD
void VanDerMeerAnalysis::LoadVertexCountsUnfold(TString p_path, Int_t p_bcid) {
  //First setup unfolding, if not done yet
  //init unfolding machinery
  if (!m_unfold) {
    cout << "[VanDerMeerAnalysis] Setting up unfolding " << endl;
    m_unfold = new VertexLumiUnfold();
    TString tf_name;
    tf_name = responseMatrixHName;
    tf_name += "_NTrk";
    tf_name += ntrkcut;
    if (m_unfold->SetResponseMatrix(responseMatrixFileName, tf_name) != 0) {
      cerr << "Error setupping unfolding with Matrix TF " << tf_name
           << " from file " << responseMatrixHName << endl;
      exit(2);
    }
    m_unfold->SetUnfoldMethod(RooUnfold::kBayes, 4);
    //m_unfold->SetUnfoldMethod(RooUnfold::kSVD, 4);
  }


  //Then load histograms
  TFile *f_in = new TFile(p_path, "READ");
  if (!f_in) {
    cerr << "Cannot open " << p_path << endl;
    exit(1);
  }

  TDirectory *f_dir = f_in;
  if (f_in->FindKey("hist") != 0) {
    //Assume it's a directory (should be), so cd in
    f_in->cd("hist");
    f_dir = gDirectory;
  }
  TIter nextH(f_dir->GetListOfKeys());
  TKey *keyH;
  while ((keyH = (TKey*)nextH())) {
    TClass *cl = gROOT->GetClass(keyH->GetClassName());
    //cout << "Pre-proc " << cl->GetName() << " , " << keyTest->GetName() << endl;
    if (!cl->InheritsFrom("TH2F")) {
      continue;  //expect our interesting histograms to be TH2F
    }

    //try to load reco observable
    TH2F *h_nVtxpLB = (TH2F*)keyH->ReadObj();
    if (!h_nVtxpLB) {
      cerr << "WARNING: Cannot load histogram from file " << p_path
           << ".Very dangerous, do not ignore this error if you don't know why." << endl;
      h_nVtxpLB = 0;
      //will skip this one.
    } else {
      h_nVtxpLB = (TH2F*) h_nVtxpLB->Clone();
      h_nVtxpLB->SetDirectory(0);
    }
    TString h_nVtxpLB_name = h_nVtxpLB->GetName();
    if (!h_nVtxpLB_name.BeginsWith("NVtxpLB_")) {
      delete h_nVtxpLB;
      continue; //not interesting for us
    }

    //now we have the histogram, get its properties
    Int_t current_bcid(-1);
    Int_t nTrk(-1);
    TObjArray *h_name_tokens = h_nVtxpLB_name.Tokenize("_");
    if (h_name_tokens->GetEntries() != 3) {
      cerr << "WARNING: Unrecognized histogram name: " << h_nVtxpLB_name << endl;
      delete h_nVtxpLB;
      continue; //skip it
    }
    TObjString *tok_bcid = (TObjString*)h_name_tokens->At(1);
    TObjString *tok_nTrk = (TObjString*)h_name_tokens->At(2);
    TString str_bcid;
    TString str_nTrk;
    str_bcid = tok_bcid->String().Remove(0, TString("BCID").Length());
    str_nTrk = tok_nTrk->String().Remove(0, TString("NTrkCut").Length());
    if (!str_bcid.IsDigit()) {
      cerr << "Malformed histogram name (BCID): " << h_nVtxpLB_name << endl;
      delete h_nVtxpLB;
      continue;
    }
    if (!str_nTrk.IsDigit()) {
      cerr << "Malformed histogram name (NTrkCut): " << h_nVtxpLB_name << endl;
      delete h_nVtxpLB;
      continue;
    }
    current_bcid = str_bcid.Atoi();
    nTrk = str_nTrk.Atoi();

    //select bcid and nTrk cuts
    if (current_bcid != p_bcid) {
      delete h_nVtxpLB;
      continue;
    }
    if (nTrk != ntrkcut) {
      delete h_nVtxpLB;
      continue;
    }

    //now loop over defined pLB, first x-scans, then y-scans
    vector<Int_t> allpLB;
    allpLB.insert(allpLB.end(), plb_list_x.begin(), plb_list_x.end());
    allpLB.insert(allpLB.end(), plb_list_y.begin(), plb_list_y.end());
    for (vector<Int_t>::iterator it = allpLB.begin(); it != allpLB.end(); ++it) {

      TString h_nVtxName;
      h_nVtxName = "h_nVtx_BCID";
      h_nVtxName += current_bcid;
      h_nVtxName += "_NTrkCut";
      h_nVtxName += ntrkcut;
      h_nVtxName += "_pLB";
      h_nVtxName += *it;

      Int_t h_bin_pLB = h_nVtxpLB->GetYaxis()->FindBin(*it);
      TH1D *h_nVtx = h_nVtxpLB->ProjectionX(h_nVtxName, h_bin_pLB, h_bin_pLB);
      if (!h_nVtx) {
        cerr << "WARNING: pLB " << *it << " not defined in input histogram: " << h_nVtxpLB->GetName() << endl;
        continue;
      }

      //apply corrections, unfold and store results
      UnfoldRawVertexCounts(h_nVtx, current_bcid, ntrkcut, *it);
      delete h_nVtx;
    }
    delete h_nVtxpLB;

  } //end loop over keys of the input file

  f_in->Close();
}

#ifdef DEBUG_VDM
// Quyick porting for debugging
Double_t debugMeanError( TH1D *m_h_trueVtx, Double_t *error) {
  Double_t x, w, err;

  Double_t stats[4] = {0.0, 0.0, 0.0, 0.0}; //sum of: weights, weights^2, x*weight, x*x*weight
  for (Int_t binx = 1; binx <= m_h_trueVtx->GetNbinsX(); ++binx) {
    x   = m_h_trueVtx->GetXaxis()->GetBinLowEdge(binx); //assume bins a left-edge centered
    w   = TMath::Abs(m_h_trueVtx->GetBinContent(binx));
    err = TMath::Abs(m_h_trueVtx->GetBinError(binx));
    stats[0] += w;
    stats[1] += err*err;
    stats[2] += w*x;
    stats[3] += w*x*x;
  }

  /*
  if (error) {
    Double_t effEntries = stats[1] ? stats[0]*stats[0]/stats[1] : TMath::Abs(stats[0])
    *error = effEntries > 0 ? TMath::Sqrt((stats[3]-stats[2]*stats[2]) / effEntries) : 0;
  }
  */
  if (error) {
    *error = m_h_trueVtx->GetMeanError();
  }

  if (stats[0] > 0) {
    return stats[2] / stats[0];
  }
  return -1; //no data?
}
#endif

void VanDerMeerAnalysis::UnfoldRawVertexCounts(TH1D *h_nVtx, Int_t p_bcid, Int_t nTrkCut, Int_t pLB) {
  //Duration of the current pLB
  Float_t current_duration = pLB_timestamps[pLB].second - pLB_timestamps[pLB].first;

  //Total number of bunch-crossings in the pLB = Delta_T * f_r
  Long64_t nBC = current_duration * 11245.5;

  //Total number of triggerable events after prescale (i.e. 100% trigger efficiency after prescale)
  Long64_t nBC_prescaled = nBC / prescale;

  //Number of triggered events that end in the plot
  // note: if does not enter the plot, we should consider it a 0-vertex anyway (just a check now)
  Long64_t nTrig = h_nVtx->GetEntries();

  //This is the number of non-triggered events, excluding prescale effects
  Long64_t nMissedEvents = nBC_prescaled - nTrig;

  //Live-fraction accounting for dead-time
  //(note that here we're only interested in the livefraction of this pLB,
  //  it's irrelevant what other pLB do, since we get mu_vis right away)
  Double_t liveFraction(1.0);
  if (init_deadtime) {
    liveFraction = live_fractions[pLB];
  }

  cout << "Unfolding " << h_nVtx->GetName() << endl;
  cout << " Trigger correction (0-vertex bin): " << nMissedEvents << " * " << liveFraction << "(lf)"
       << " = " << nMissedEvents*liveFraction << "(" << nMissedEvents*liveFraction / nTrig << ")" << endl;
#ifdef DEBUG_VDM
  cout << " Number of BC (prescaled BC)= " << nBC << "(" << nBC_prescaled << ")"
       << ", livetime: " << liveFraction << ", duration: " << current_duration << endl;
  cout << " Triggered events: " << nTrig << endl;
  cout << " (Old) Raw mean: " << h_nVtx->GetMean() << " +/- " << h_nVtx->GetMeanError() << endl;
  Double_t debug_mean, debug_error;
  debug_mean = debugMeanError(h_nVtx, &debug_error);
  cout << " New mean: " << debug_mean << " +/- " << debug_error << endl;
#endif

  //heavy-debug //////////////////
  /*
  cout << "=========================" << endl;
  cout << "HEAVY-DEBUG: Original h_nVtx:" << endl;
  h_nVtx->Print("all");
  cout << "=========================" << endl;
  */
  /////////////////////////////////

  //apply corrections
  if (trigger_type != "BGRP7") { //this is a random trigger, typically EF_rd0_Filled_NoAlg
    //If not random trigger, we need to correct for this
    //Assume 100% trigger efficiency. Inject 0-vertex events as needed: nMissedEvents*liveFraction
    // Also assumes livefraction constant in pLB
    // Note: This is also the place where a trigger efficiency can be applied
    Int_t zeroBin = h_nVtx->FindBin(0.0);
    double val, err;
    val = h_nVtx->GetBinContent(zeroBin)+nMissedEvents*liveFraction;
    err = TMath::Power(h_nVtx->GetBinError(zeroBin),2) + nMissedEvents*liveFraction;
    err = TMath::Sqrt(err);
    h_nVtx->SetBinContent(zeroBin, val);
    h_nVtx->SetBinError(zeroBin, err);
  } else {
    //this is supposed to be random trigger... check this at 0.1% level
    if (nMissedEvents*liveFraction / nTrig > 0.001)
      cerr << "WARNING: Missing events in a random-trigger:" << nMissedEvents
           << "( / " << nTrig << " ~(liveFraction) " << nMissedEvents*liveFraction / nTrig << ")"
           << " for historgam " << h_nVtx->GetName() << endl;
  }

  //heavy-debug ///////////////////
  /*
  cout << "=========================" << endl;
  cout << "HEAVY-DEBUG: Corrected h_nVtx:" << endl;
  cout << " Raw mean: " << h_nVtx->GetMean() << " +/- " << h_nVtx->GetMeanError() << endl;
  debug_mean = debugMeanError(h_nVtx, &debug_error);
  cout << " Debug mean: " << debug_mean << "+ +/- " << debug_error << endl;
  h_nVtx->Print("all");
  cout << "=========================" << endl;
  */
  /////////////////////////////////

  //and unfold it!
  m_unfold->VertexUnfold(h_nVtx);
  TH1D* h_nvtx_unfold = (TH1D*)m_unfold->GetInteractionsSpectrum();
  TString recoUnfoldSpectrumName = "unfolded_";
  recoUnfoldSpectrumName += h_nVtx->GetName();
  h_nvtx_unfold->SetName(recoUnfoldSpectrumName.Data());
  recoUnfoldSpectrumName = "Unfolded spectrum ";
  recoUnfoldSpectrumName += h_nVtx->GetName();
  h_nvtx_unfold->SetTitle(recoUnfoldSpectrumName.Data());

  //finally store results (this is actually already mu_vis)
  Double_t errMu;
  nvtx_pLB[pLB] = m_unfold->GetMeanInteractions(&errMu);//same as h_nvtx_unfold->GetMean()
  nvtx_err_pLB[pLB] = errMu;//h_nvtx_unfold->GetMeanError();

  //heavy-debug ///////////////////
  /*
  cout << "=========================" << endl;
  cout << "HEAVY-DEBUG: Unfolded h_nVtx:" << endl;
  cout << " Raw mean: " << h_nvtx_unfold->GetMean() << " +/- " << h_nvtx_unfold->GetMeanError() << endl;
  debug_mean = debugMeanError(h_nvtx_unfold, &debug_error);
  cout << " Debug mean: " << debug_mean << "+ +/- " << debug_error << endl;
  h_nvtx_unfold->Print("all");
  cout << "=========================" << endl;
  */
  /////////////////////////////////


  //Debug
  cout << "Ufolded " << h_nVtx->GetName() << "; "
       << "Unfolded mu: " << nvtx_pLB[pLB] << " +/- " << nvtx_err_pLB[pLB] << endl;
#ifdef DEBUG_VDM
  cout << " Old Mean mu_rec: " << h_nVtx->GetMean() << " +/- " << h_nVtx->GetMeanError() << endl;
#endif

}
#endif

void VanDerMeerAnalysis::ConvertToRandom(Int_t bunches, Float_t rate) {

  cout << "[VanDerMeerAnalysis] WARNING : Converting uncertainties to those from a random trigger at " << rate << " Hz, " << bunches << " bunches." << endl;

  Float_t recorded_fraction = rate / bunches / 11245.5;

  for (vector<Int_t>::iterator it = plb_list_x.begin(); it != plb_list_x.end(); ++it) {

    Double_t current_duration = pLB_timestamps_x[*it].second - pLB_timestamps_x[*it].first;
    Float_t counts = mu_pLB[*it] * 11245.5 * current_duration;
    Float_t counts_recorded = counts * recorded_fraction;
    Float_t counts_recorded_err = (0.5 + TMath::Sqrt(0.25 + counts_recorded)) / recorded_fraction;

    mu_err_pLB[*it] = counts_recorded_err / (11245.5 * current_duration);
  }

  for (vector<Int_t>::iterator it = plb_list_y.begin(); it != plb_list_y.end(); ++it) {

    Double_t current_duration = pLB_timestamps_y[*it].second - pLB_timestamps_y[*it].first;
    Float_t counts = mu_pLB[*it] * 11245.5 * current_duration;
    Float_t counts_recorded = counts * recorded_fraction;
    Float_t counts_recorded_err = (0.5 + TMath::Sqrt(0.25 + counts_recorded)) / recorded_fraction;

    mu_err_pLB[*it] = counts_recorded_err / (11245.5 * current_duration);
  }


}

#ifdef REMOVED_051612
void VanDerMeerAnalysis::InitializeMaskingCorrection(TString p_energy, TString p_settings, Int_t p_ntrkcut, TH1D *h_z, TString save_tag) {
  /**
    *  h_z is used to generate an expected dz distribution.
    *  Expected dz histogram is saved to include/cache.root : h_dz_expected_<save_tag>.
    */

  cout << "[VanDerMeerAnalysis] INFO : Initializing masking correction." << endl;

  TString s_tag = "data_";
  s_tag += p_energy;
  s_tag += "TeV_";
  s_tag += p_settings;
  pmc = new PileupMaskingCorrection(s_tag, p_ntrkcut);

  pmc->GenerateDzDistribution(h_z);
  pmc->GenerateCorrection(pmc->GetExpectedDzDistribution());

  TFile *f_cache = new TFile("include/cache.root", "UPDATE");

  pmc->GetExpectedDzDistribution()->Write(TString("h_dz_expected_") + save_tag, TObject::kOverwrite);
  pmc->GetMuMap()->Write(TString("tg_mu_obs_vs_mu_actual_") + save_tag, TObject::kOverwrite);
  pmc->GetMuCorrection()->Write(TString("tg_pileup_correction_") + save_tag, TObject::kOverwrite);

  f_cache->Close();

  init_pmc = true;
}

void VanDerMeerAnalysis::InitializeMaskingCorrection(TString p_energy, TString p_settings, Int_t p_ntrkcut, TString save_tag) {
  /**
    *  Load some time-intensive histograms from include/cache.root.
    *  dz histograms saved/loaded as h_dz_expected_<save_tag>
    */

  cout << "[VanDerMeerAnalysis] INFO : Initializing masking and fake correction, using cached masking correction." << endl;

  TString s_ntrkcut = "";
  s_ntrkcut += p_ntrkcut;
  TString s_tag = "data_";
  s_tag += p_energy;
  s_tag += "TeV_";
  s_tag += p_settings;
  pmc = new PileupMaskingCorrection(s_tag, p_ntrkcut);

  TFile *f_cache = new TFile("include/cache.root", "READ");
  TH1F *h_dz_expected = (TH1F*)f_cache->Get(TString("h_dz_expected_") + save_tag);
  TGraphErrors *tg_mu_obs_vs_mu_actual = (TGraphErrors*)f_cache->Get(TString("tg_mu_obs_vs_mu_actual_") + save_tag);
  TGraphErrors *tg_pileup_correction = (TGraphErrors*)f_cache->Get(TString("tg_pileup_correction_") + save_tag);
  pmc->LoadCorrection(tg_mu_obs_vs_mu_actual, tg_pileup_correction);
  f_cache->Close();

  init_pmc = true;
}
#endif


void VanDerMeerAnalysis::CorrectPileupEffects(TH2D *h_z_plb, TString p_energy, TString p_settings, Int_t p_ntrkcut, TString p_run, Int_t p_bcid) {

  /**
    *  New method as of 5-16-12: need to calculate a unique masking correction for each plb.
    *  So, input is the z distribution for the plb, plus a tag and ntrkcut to choose the pmask_dz distribution from cache.
    *  Order of things:
    *    1. Do an initial pileup correction.
    *    2. Look up mu_fake with the pileup-corrected mu_raw
    *    3. Subtract mu_fake to get mu_real = mu_raw - mu_fake
    *    4. Look up and apply another pileup correction, this time with only mu_real.
    */

  cout << "[VanDerMeerAnalysis] INFO : Correcting pileup effects" << endl;

  for (vector<Int_t>::iterator it = plb_list.begin(); it != plb_list.end(); ++it) {

    cout << endl;
    cout << endl;
    cout << "///////////////////////////////////////////////////////////" << endl;
    cout << "[VanDerMeerAnalysis] INFO : On " << p_energy << " TeV, " << p_settings << ", NTrk" << p_ntrkcut << ", pLB " << *it << endl;
    cout << "///////////////////////////////////////////////////////////" << endl;
    cout << endl;
    cout << endl;

    Double_t mu_raw = mu_pLB[*it];
    Double_t initial_masking_correction_factor = 1.;
    PileupMaskingCorrection *mc;
    mu_raw_pLB[*it] = mu_pLB[*it];

    if (vtx_method == "Vtx") {
      TString p_tag = "data_";
      p_tag += p_energy;
      p_tag += "TeV_";
      p_tag += p_settings;
      p_tag += "_";
      p_tag += p_run;

      Int_t bin = h_z_plb->GetYaxis()->FindBin(*it);
      TString hname = "h_z_plb";
      hname += *it;
      TH1D *h_z = (TH1D*)h_z_plb->ProjectionX(hname, bin, bin);
      cout << "[VanDerMeerAnalysis] INFO : Writting h_z_plb histogram into root file" << endl;
      TFile *f_z_plb = new TFile(hname+".root", "new");
      h_z->Write();
      f_z_plb->Close();

      cout << "[VanDerMeerAnalysis] INFO : h_z RMS = " << h_z->GetRMS() << endl;

      //mc = new PileupMaskingCorrection(p_tag, p_ntrkcut);
      mc = new PileupMaskingCorrection(p_tag, p_ntrkcut, p_bcid, true);
      if (systematic_uncertainty_list["masking_toy_scaling"]) {
        mc->SetPmaskScale(0.98);
      }
      TString pmc_tag = "";
      pmc_tag += p_energy;
      pmc_tag += "_";
      pmc_tag += p_settings;
      pmc_tag += "_";
      pmc_tag += p_ntrkcut;
      pmc_tag += "_";
      pmc_tag += *it;
      //mc->GenerateDzDistribution(h_z, pmc_tag);
      //mc->GenerateCorrection(mc->GetExpectedDzDistribution());
      initial_masking_correction_factor = mc->GetCorrectionFactor(mu_raw);
      cout << "[VanDerMeerAnalysis] INFO : mu_raw = " << mu_raw << endl;
      cout << "[VanDerMeerAnalysis] INFO : imcf = " << initial_masking_correction_factor << endl;

      // // Load PileupMaskingCorrection factors from cached file
      // PileupMaskingCorrection mc(p_tag, p_ntrkcut);
      // initial_masking_correction_factor = mc.GetCorrectionFactor(mu_raw);
      
      if (initial_masking_correction_factor < 1.) {
        cout << "[VanDerMeerAnalysis] WARNING : pLB " << *it << " has initial masking correction factor " << initial_masking_correction_factor << " < 1. Setting to 1." << endl;
        initial_masking_correction_factor = 1.;
      }

      Double_t mu_fake;
      if (systematic_uncertainty_list["fake_low"]) {
        mu_fake = fc->GetFakeMuLowFromMuReconMC(mu_raw * initial_masking_correction_factor);
      } else if (systematic_uncertainty_list["fake_high"]) {
        mu_fake = fc->GetFakeMuHighFromMuReconMC(mu_raw * initial_masking_correction_factor);
      } else {
        mu_fake = fc->GetFakeMuFromMuReconMC(mu_raw * initial_masking_correction_factor);
      }

      Double_t mu_fake_uncertainty = fc->GetFakeMuUncertaintyFromMuReconMC(mu_raw * initial_masking_correction_factor);
      if (mu_fake < 0. || mu_fake > mu_raw) {
        cout << "[VanDerMeerAnalysis] WARNING : pLB " << *it << " has invalid mu_fake " << mu_fake << ".  Setting to 0, but needs debugging!!" << endl;
        mu_fake = 0.;
      }

      Double_t mu_real = mu_raw - mu_fake;
      mu_real_pLB[*it] = mu_real;

      cout << "[VanDerMeerAnalysis] INFO : mu_fake = " << mu_fake << endl;
      cout << "[VanDerMeerAnalysis] INFO : mu_real = " << mu_real << endl;

      Double_t final_masking_correction_factor = 1.;
      if (vtx_method == "Vtx") {
        final_masking_correction_factor = mc->GetCorrectionFactor(mu_real);
        //final_masking_correction_factor = mc.GetCorrectionFactor(mu_real);
        if (final_masking_correction_factor < 1.) {
          cout << "[VanDerMeerAnalysis] WARNING : pLB " << *it << " has final masking correction factor " << final_masking_correction_factor << " < 1. Setting to 1." << endl;
          final_masking_correction_factor = 1.;
        }
      }

      mu_pLB[*it] = mu_real * final_masking_correction_factor;
      mu_err_pLB[*it] =  TMath::Sqrt(TMath::Power(mu_err_pLB[*it], 2) + TMath::Power(mu_fake_uncertainty, 2)) * final_masking_correction_factor;

      mu_fake_list[*it] = mu_fake;
      masking_correction_factors[*it] = final_masking_correction_factor;

      cout << "[VanDerMeerAnalysis] INFO : fmcf = " << final_masking_correction_factor << endl;
      cout << "[VanDerMeerAnalysis] INFO : mu_vis = " << mu_pLB[*it] << endl;

    } else if (vtx_method == "NEvt") {

      Float_t mu_fake = fc->GetMuFakeFromMuReconNEvt(mu_pLB[*it]);
      Float_t mu_fake_err = fc->GetMuFakeUncertaintyFromMuReconNEvt(mu_pLB[*it]);
      if (mu_fake < 0.) {
        cout << "[VanDerMeerAnalysis] WARNING : pLB " << *it << " has mu_fake = " << mu_fake << ". Setting to 0." << endl;
        mu_fake = 0.;
        mu_fake_err = 0.;
      }
      mu_pLB[*it] = mu_pLB[*it] - mu_fake;
      mu_err_pLB[*it] = TMath::Sqrt(TMath::Power(mu_err_pLB[*it], 2) + TMath::Power(mu_fake_err, 2));
      mu_fake_list[*it] = mu_fake;

    }
  }

  #ifdef REMOVED_051612
  Double_t old_mu = mu_pLB[*it];
  //if (old_mu < 0.05) continue; // I don't think our corrections are well-understood in this region.

  Double_t masking_correction_factor = (init_pmc ? pmc->GetCorrectionFactor(old_mu) : 1.);

  if (masking_correction_factor < 1.) {
    masking_correction_factor = 1.;
  }

  // Apply corrections
  Double_t old_mu_AMC = (init_pmc ? old_mu * pmc->GetCorrectionFactor(old_mu) : old_mu);
  Double_t fake_fraction = (init_fc ? fc->GetFakeFractionFromMuVisAMC(old_mu_AMC) : 0.);
  mu_pLB[*it] = mu_pLB[*it] * (1. - fake_fraction) * masking_correction_factor;
  mu_err_pLB[*it] = mu_err_pLB[*it] * (1. - fake_fraction) * masking_correction_factor;

  // Save correction factors
  fake_fractions[*it] = fake_fraction;
  masking_correction_factors[*it] = masking_correction_factor;

  }
  #endif

}

void VanDerMeerAnalysis::CalculateMuSpecPlb() {

  cout << "[VanDerMeerAnalysis] INFO : Calculating mu_sp" << endl;

  for (vector<Int_t>::iterator it = plb_list_x.begin(); it != plb_list_x.end(); ++it) {
    musp_pLB[*it] = mu_pLB[*it] / (bunch_intensities_1[*it] * bunch_intensities_2[*it]);
    musp_err_pLB[*it] = mu_err_pLB[*it] / (bunch_intensities_1[*it] * bunch_intensities_2[*it]);
  }

  for (vector<Int_t>::iterator it = plb_list_y.begin(); it != plb_list_y.end(); ++it) {
    musp_pLB[*it] = mu_pLB[*it] / (bunch_intensities_1[*it] * bunch_intensities_2[*it]);
    musp_err_pLB[*it] = mu_err_pLB[*it] / (bunch_intensities_1[*it] * bunch_intensities_2[*it]);
  }
}

/*void VanDerMeerAnalysis::FitVdmCurves() {
  cout << "[VanDerMeerAnalysis] INFO : Fitting VdM curves" << endl;

  // -- Make TGraphs
  // First, check that there are an equal number of x and y points
  if (plb_list_x.size() != plb_list_y.size()) {
    cerr << "[VanDerMeerAnalysis] ERROR : Different number of scan points in X and Y scans! Exiting..." << endl;
    exit(1);
  }

  tg_musp_x = new TGraphErrors(plb_list_x.size());
  tg_musp_x->SetName(TString("tg_musp_x_") + save_tag);
  tg_musp_y = new TGraphErrors(plb_list_y.size());
  tg_musp_y->SetName(TString("tg_musp_y_") + save_tag);

  Double_t mu_max_guess_x, mu_max_guess_y;

  for (unsigned int i=0; i<plb_list_x.size(); i++) {

    tg_musp_x->SetPoint(i, nominal_separation[plb_list_x[i]], musp_pLB[plb_list_x[i]]);
    tg_musp_x->SetPointError(i, 0., musp_err_pLB[plb_list_x[i]]);
    tg_musp_y->SetPoint(i, nominal_separation[plb_list_x[i]], musp_pLB[plb_list_y[i]]);
    tg_musp_y->SetPointError(i, 0., musp_err_pLB[plb_list_y[i]]);

    if (i == 12) {
      mu_max_guess_x = musp_pLB[plb_list_x[i]];
      mu_max_guess_y = musp_pLB[plb_list_y[i]];
    }
  }

  fit_functions.push_back("sgcl");
  cout << "Line 1152, low_displacement = " << low_displacement << ", high_displacement = " << high_displacement << endl;
  f_x["sgcl"] = new TF1("f_x_sgcl", fit_func_sgc, low_displacement, high_displacement, 4);
  f_x["sgcl"]->SetNpx(10000);

  f_x["sgcl"]->SetParName(0, "mu_max");
  f_x["sgcl"]->SetParName(1, "sigma_a");
  f_x["sgcl"]->SetParName(2, "h0");
  f_x["sgcl"]->SetParName(3, "c");

  f_x["sgcl"]->SetParameter(0, mu_max_guess_x);
  f_x["sgcl"]->SetParameter(1, 55.);
  f_x["sgcl"]->SetParameter(2, 0.);
  f_x["sgcl"]->SetParameter(3, mu_max_guess_x / 1000.);

  f_x["sgcl"]->SetParLimits(0, mu_max_guess_x / 100., mu_max_guess_x * 100.);
  f_x["sgcl"]->SetParLimits(1, 10., 500.);
  f_x["sgcl"]->SetParLimits(2, -100., 100.);
  f_x["sgcl"]->SetParLimits(3, 0., mu_max_guess_x / 10.);

  f_y["sgcl"] = new TF1("f_y_sgcl", fit_func_sgc, low_displacement, high_displacement, 4);
  f_y["sgcl"]->SetNpx(10000);

  f_y["sgcl"]->SetParName(0, "mu_max");
  f_y["sgcl"]->SetParName(1, "sigma_a");
  f_y["sgcl"]->SetParName(2, "h0");
  f_y["sgcl"]->SetParName(3, "c");

  f_y["sgcl"]->SetParameter(0, mu_max_guess_y);
  f_y["sgcl"]->SetParameter(1, 55.);
  f_y["sgcl"]->SetParameter(2, 0.);
  f_y["sgcl"]->SetParameter(3, mu_max_guess_y / 1000.);

  f_y["sgcl"]->SetParLimits(0, mu_max_guess_y / 100., mu_max_guess_y * 100.);
  f_y["sgcl"]->SetParLimits(1, 10., 500.);
  f_y["sgcl"]->SetParLimits(2, -100., 100.);
  f_y["sgcl"]->SetParLimits(3, 0., mu_max_guess_y / 10.);


  fit_functions.push_back("sgc");
  f_x["sgc"] = new TF1("f_x_sgc", fit_func_sgc, low_displacement, high_displacement, 4);
  f_x["sgc"]->SetNpx(10000);

  f_x["sgc"]->SetParName(0, "mu_max");
  f_x["sgc"]->SetParName(1, "sigma_a");
  f_x["sgc"]->SetParName(2, "h0");
  f_x["sgc"]->SetParName(3, "c");

  f_x["sgc"]->SetParameter(0, mu_max_guess_x);
  f_x["sgc"]->SetParameter(1, 55.);
  f_x["sgc"]->SetParameter(2, 0.);
  f_x["sgc"]->SetParameter(3, mu_max_guess_x / 1000.);

  f_x["sgc"]->SetParLimits(0, mu_max_guess_x / 100., mu_max_guess_x * 100.);
  f_x["sgc"]->SetParLimits(1, 10., 500.);
  f_x["sgc"]->SetParLimits(2, -100., 100.);
  f_x["sgc"]->SetParLimits(3, 0., mu_max_guess_x / 10.);

  f_y["sgc"] = new TF1("f_y_sgc", fit_func_sgc, low_displacement, high_displacement, 4);
  f_y["sgc"]->SetNpx(10000);

  f_y["sgc"]->SetParName(0, "mu_max");
  f_y["sgc"]->SetParName(1, "sigma_a");
  f_y["sgc"]->SetParName(2, "h0");
  f_y["sgc"]->SetParName(3, "c");

  f_y["sgc"]->SetParameter(0, mu_max_guess_y);
  f_y["sgc"]->SetParameter(1, 55.);
  f_y["sgc"]->SetParameter(2, 0.);
  f_y["sgc"]->SetParameter(3, mu_max_guess_y / 1000.);

  f_y["sgc"]->SetParLimits(0, mu_max_guess_y / 100., mu_max_guess_y * 100.);
  f_y["sgc"]->SetParLimits(1, 10., 500.);
  f_y["sgc"]->SetParLimits(2, -100., 100.);
  f_y["sgc"]->SetParLimits(3, 0., mu_max_guess_y / 10.);


  fit_functions.push_back("sg");
  f_x["sg"] = new TF1("f_x_sg", fit_func_sgc, low_displacement, high_displacement, 4);
  f_x["sg"]->SetNpx(10000);

  f_x["sg"]->SetParName(0, "mu_max");
  f_x["sg"]->SetParName(1, "sigma_a");
  f_x["sg"]->SetParName(2, "h0");
  f_x["sg"]->SetParName(3, "c");

  f_x["sg"]->SetParameter(0, mu_max_guess_x);
  f_x["sg"]->SetParameter(1, 55.);
  f_x["sg"]->SetParameter(2, 0.);
  f_x["sg"]->FixParameter(3, 0.);

  f_x["sg"]->SetParLimits(0, mu_max_guess_x / 100., mu_max_guess_x * 100.);
  f_x["sg"]->SetParLimits(1, 10., 500.);
  f_x["sg"]->SetParLimits(2, -100., 100.);

  f_y["sg"] = new TF1("f_y_sg", fit_func_sgc, low_displacement, high_displacement, 4);
  f_y["sg"]->SetNpx(10000);

  f_y["sg"]->SetParName(0, "mu_max");
  f_y["sg"]->SetParName(1, "sigma_a");
  f_y["sg"]->SetParName(2, "h0");
  f_y["sg"]->SetParName(3, "c");

  f_y["sg"]->SetParameter(0, mu_max_guess_y);
  f_y["sg"]->SetParameter(1, 55.);
  f_y["sg"]->SetParameter(2, 0.);
  f_y["sg"]->FixParameter(3, 0.);

  f_y["sg"]->SetParLimits(0, mu_max_guess_y / 100., mu_max_guess_y * 100.);
  f_y["sg"]->SetParLimits(1, 10., 500.);
  f_y["sg"]->SetParLimits(2, -100., 100.);

  ///// Double Gaussian /////
  fit_functions.push_back("dg");
  //f_x["dg"] = new TF1("f_x_dg", fit_func_dg, low_displacement, high_displacement, 6);
  f_x["dg"] = new TF1("f_x_dg", fit_func_dg, low_displacement, high_displacement, 6);
  f_x["dg"]->SetNpx(10000);

  f_x["dg"]->SetParName(0, "mu_max");
  f_x["dg"]->SetParName(1, "sigma_a");
  f_x["dg"]->SetParName(2, "h0");
  f_x["dg"]->SetParName(3, "c");
  f_x["dg"]->SetParName(4, "sigma_b");
  f_x["dg"]->SetParName(5, "f_a");

  f_x["dg"]->SetParameter(0, mu_max_guess_x);
  f_x["dg"]->SetParameter(1, 20.);
  f_x["dg"]->SetParameter(2, 0.);
  f_x["dg"]->FixParameter(3, 0.);
  f_x["dg"]->SetParameter(4, 40.);
  f_x["dg"]->SetParameter(5, 0.75);

  f_x["dg"]->SetParLimits(0, mu_max_guess_x / 100., mu_max_guess_x * 100.);
  f_x["dg"]->SetParLimits(1, 5., 1000.);
  f_x["dg"]->SetParLimits(2, -100., 100.);
  f_x["dg"]->SetParLimits(4, 5., 1000.);
  f_x["dg"]->SetParLimits(5, 0., 1.);

  //f_y["dg"] = new TF1("f_y_dg", fit_func_dg, low_displacement, high_displacement, 6);
  f_y["dg"] = new TF1("f_y_dg", fit_func_dg, low_displacement, high_displacement, 6);
  f_y["dg"]->SetNpx(10000);

  f_y["dg"]->SetParName(0, "mu_max");
  f_y["dg"]->SetParName(1, "sigma_a");
  f_y["dg"]->SetParName(2, "h0");
  f_y["dg"]->SetParName(3, "c");
  f_y["dg"]->SetParName(4, "sigma_b");
  f_y["dg"]->SetParName(5, "f_a");

  f_y["dg"]->SetParameter(0, mu_max_guess_y);
  f_y["dg"]->SetParameter(1, 20.);
  f_y["dg"]->SetParameter(2, 0.);
  f_y["dg"]->FixParameter(3, 0.);
  f_y["dg"]->SetParameter(4, 40.);
  f_y["dg"]->SetParameter(5, 0.75);

  f_y["dg"]->SetParLimits(0, mu_max_guess_y / 100., mu_max_guess_y * 100.);
  f_y["dg"]->SetParLimits(1, 5.,1000.);
  f_y["dg"]->SetParLimits(2, -100., 100.);
  f_y["dg"]->SetParLimits(4, 5., 1000.);
  f_y["dg"]->SetParLimits(5, 0., 1.);

  ///// Double Gaussian plus constant (not luminosity) /////
  fit_functions.push_back("dgc");
  f_x["dgc"] = new TF1("f_x_dgc", fit_func_dgc, low_displacement, high_displacement, 6);
  f_x["dgc"]->SetNpx(10000);

  f_x["dgc"]->SetParName(0, "mu_max");
  f_x["dgc"]->SetParName(1, "sigma_a");
  f_x["dgc"]->SetParName(2, "h0");
  f_x["dgc"]->SetParName(3, "c");
  f_x["dgc"]->SetParName(4, "Sigma_x");
  f_x["dgc"]->SetParName(5, "f_a");

  f_x["dgc"]->SetParameter(0, 3e-22);
  f_x["dgc"]->SetParameter(1, 20.);
  f_x["dgc"]->SetParameter(2, 0.);
  f_x["dgc"]->SetParameter(3, 3.5e-36);
  f_x["dgc"]->SetParameter(4, 26.);
  f_x["dgc"]->SetParameter(5, 0.25);

  f_x["dgc"]->SetParLimits(0, 3.0e-24, 3.0e-20);
  f_x["dgc"]->SetParLimits(1, 5., 1000.);
  f_x["dgc"]->SetParLimits(2, -100., 100.);
  f_x["dgc"]->SetParLimits(3, 0., 3.e-34);
  f_x["dgc"]->SetParLimits(4, 5., 1000.);
  f_x["dgc"]->SetParLimits(5, 0., 1.);
  
  cout << "mu_max_guess_x = " << mu_max_guess_x << endl;

  f_y["dgc"] = new TF1("f_y_dgc", fit_func_dgc, low_displacement, high_displacement, 6);
  f_y["dgc"]->SetNpx(10000);

  f_y["dgc"]->SetParName(0, "mu_max");
  f_y["dgc"]->SetParName(1, "sigma_a");
  f_y["dgc"]->SetParName(2, "h0");
  f_y["dgc"]->SetParName(3, "c");
  f_y["dgc"]->SetParName(4, "Sigma_y");
  f_y["dgc"]->SetParName(5, "f_a");

  f_y["dgc"]->SetParameter(0, 3.2e-22);
  f_y["dgc"]->SetParameter(1, 30.);
  f_y["dgc"]->SetParameter(2, 0.);
  f_y["dgc"]->SetParameter(3, 1.e-26);
  f_y["dgc"]->SetParameter(4, 33.);
  f_y["dgc"]->SetParameter(5, 0.32);

  f_y["dgc"]->SetParLimits(0, 3.2e-24, 3.2e-20);
  f_y["dgc"]->SetParLimits(1, 5.,1000.);
  f_y["dgc"]->SetParLimits(2, -100., 100.);
  f_y["dgc"]->SetParLimits(3, 0., 2.e-24);
  f_y["dgc"]->SetParLimits(4, 5., 1000.);
  f_y["dgc"]->SetParLimits(5, 0., 1.);
  
  cout << "mu_max_guess_y = " << mu_max_guess_y << endl;

  ///// Double Gaussian plus constant (luminosity)/////
  fit_functions.push_back("dgcl");
  f_x["dgcl"] = new TF1("f_x_dgcl", fit_func_dgcl, low_displacement, high_displacement, 6);
  f_x["dgcl"]->SetNpx(10000);

  f_x["dgcl"]->SetParName(0, "mu_max");
  f_x["dgcl"]->SetParName(1, "sigma_a");
  f_x["dgcl"]->SetParName(2, "h0");
  f_x["dgcl"]->SetParName(3, "c");
  f_x["dgcl"]->SetParName(4, "sigma_b");
  f_x["dgcl"]->SetParName(5, "f_a");

  f_x["dgcl"]->SetParameter(0, mu_max_guess_x);
  f_x["dgcl"]->SetParameter(1, 20.);
  f_x["dgcl"]->SetParameter(2, 0.);
  f_x["dgcl"]->SetParameter(3, mu_max_guess_x / 1000.);
  f_x["dgcl"]->SetParameter(4, 40.);
  f_x["dgcl"]->SetParameter(5, 0.75);

  f_x["dgcl"]->SetParLimits(0, mu_max_guess_x / 100., mu_max_guess_x * 100.);
  f_x["dgcl"]->SetParLimits(1, 5., 1000.);
  f_x["dgcl"]->SetParLimits(2, -100., 100.);
  f_x["dgcl"]->SetParLimits(3, 0., mu_max_guess_x * 100.);
  f_x["dgcl"]->SetParLimits(4, 5., 1000.);
  f_x["dgcl"]->SetParLimits(5, 0., 1.);

  f_y["dgcl"] = new TF1("f_y_dgcl", fit_func_dgcl, low_displacement, high_displacement, 6);
  f_y["dgcl"]->SetNpx(10000);

  f_y["dgcl"]->SetParName(0, "mu_max");
  f_y["dgcl"]->SetParName(1, "sigma_a");
  f_y["dgcl"]->SetParName(2, "h0");
  f_y["dgcl"]->SetParName(3, "c");
  f_y["dgcl"]->SetParName(4, "sigma_b");
  f_y["dgcl"]->SetParName(5, "f_a");

  f_y["dgcl"]->SetParameter(0, mu_max_guess_y);
  f_y["dgcl"]->SetParameter(1, 20.);
  f_y["dgcl"]->SetParameter(2, 0.);
  f_y["dgcl"]->SetParameter(3, mu_max_guess_y / 1000.);
  f_y["dgcl"]->SetParameter(4, 40.);
  f_y["dgcl"]->SetParameter(5, 0.75);

  f_y["dgcl"]->SetParLimits(0, mu_max_guess_y / 100., mu_max_guess_y * 100.);
  f_y["dgcl"]->SetParLimits(1, 5.,1000.);
  f_y["dgcl"]->SetParLimits(2, -100., 100.);
  f_y["dgcl"]->SetParLimits(3, 0., mu_max_guess_y * 100.);
  f_y["dgcl"]->SetParLimits(4, 5., 1000.);
  f_y["dgcl"]->SetParLimits(5, 0., 1.);

  fit_functions.push_back("spline");
  fit_status["spline"] = true;
  fit_spline_x = new TSpline3("fit_spline_x", tg_musp_x);
  fit_spline_y = new TSpline3("fit_spline_y", tg_musp_y);

  // -- Perform fits
  for (std::vector<TString>::iterator fit_function = fit_functions.begin(); fit_function != fit_functions.end(); ++fit_function) {
    if (*fit_function != "spline") {
      //fr_x_actual[*fit_function] = *(tg_musp_x->Fit(f_x[*fit_function], "SQER0+"));
      //fr_y_actual[*fit_function] = *(tg_musp_y->Fit(f_y[*fit_function], "SQER0+"));
      //fr_x[*fit_function] = &fr_x_actual[*fit_function];
      //fr_y[*fit_function] = &fr_y_actual[*fit_function];
      TString fit_options = "SERM+";
      fr_x[*fit_function] = tg_musp_x->Fit(f_x[*fit_function], fit_options);
      fr_y[*fit_function] = tg_musp_y->Fit(f_y[*fit_function], fit_options);
      //if (!(int)fr_x[*fit_function] && !(int)fr_y[*fit_function]) {
        fit_status[*fit_function] = true;
      //} else {
      //  fit_status[*fit_function] = false;
      //}
    }
  }
}*/

///Attempt to integrate the centring correction

void VanDerMeerAnalysis::FitVdmCurves() {
  cout << "[VanDerMeerAnalysis] INFO : Fitting VdM curves" << endl;

  // -- Make TGraphs
  // First, check that there are an equal number of x and y points
  if (plb_list_x.size() != plb_list_y.size()) {
    cerr << "[VanDerMeerAnalysis] ERROR : Different number of scan points in X and Y scans! Exiting..." << endl;
    exit(1);
  }

  tg_musp_x = new TGraphErrors(plb_list_x.size());
  tg_musp_x->SetName(TString("tg_musp_x_") + save_tag);
  tg_musp_y = new TGraphErrors(plb_list_y.size());
  tg_musp_y->SetName(TString("tg_musp_y_") + save_tag);

  for (unsigned int i=0; i<plb_list_x.size(); i++) {

    tg_musp_x->SetPoint(i, nominal_separation[plb_list_x[i]], musp_pLB[plb_list_x[i]]);
    tg_musp_x->SetPointError(i, 0., musp_err_pLB[plb_list_x[i]]);
    tg_musp_y->SetPoint(i, nominal_separation[plb_list_y[i]], musp_pLB[plb_list_y[i]]);
    tg_musp_y->SetPointError(i, 0., musp_err_pLB[plb_list_y[i]]);

  }

  max_x = 0;
  for( int iPoint = 0; iPoint < tg_musp_x->GetN(); ++iPoint ) {
    double x(0), y(0);
    tg_musp_x->GetPoint( iPoint, x, y );
    if( y > max_x ) { max_x = y; }
  }

  max_y = 0;
  for( int iPoint = 0; iPoint < tg_musp_y->GetN(); ++iPoint ) {
    double x(0), y(0);
    tg_musp_y->GetPoint( iPoint, x, y );
    if( y > max_y ) { max_y = y; }
  }

  double sumwx(0), sumpointsx(0);
  for( int iPoint = 0; iPoint < tg_musp_x->GetN(); ++iPoint ) {
    double x(0), w(0);
    tg_musp_x->GetPoint( iPoint, x, w );
    sumwx += w;
    sumpointsx += x*w;
  }
  mean_x =  (sumpointsx / sumwx);

  double sumwy(0), sumpointsy(0);
  for( int iPoint = 0; iPoint < tg_musp_y->GetN(); ++iPoint ) {
    double x(0), w(0);
    tg_musp_y->GetPoint( iPoint, x, w );
    sumwy += w;
    sumpointsy += x*w;
  }
  mean_y =  (sumpointsy / sumwy);

  //mean_x = tg_musp_x->GetMean();
  //max_x = TMath::MaxElement(tg_musp_x->GetN(),tg_musp_x->GetY());
  for( int iPoint = 0; iPoint < tg_musp_x->GetN(); ++iPoint ) {
    double x(0), w(0);
    tg_musp_x->GetPoint( iPoint, x, w );
    sumw_x += w;
    sumdev2_x += (x-mean_x)*(x-mean_x)*w;
  }
  myRMS_x = ( sqrt(sumdev2_x / sumw_x) );

  //mean_y = tg_musp_y->GetMean();
  //max_y = TMath::MaxElement(tg_musp_y->GetN(),tg_musp_y->GetY());
  for( int iPoint = 0; iPoint < tg_musp_y->GetN(); ++iPoint ) {
    double x(0), w(0);
    tg_musp_y->GetPoint( iPoint, x, w );
    sumw_y += w;
    sumdev2_y += (x-mean_y)*(x-mean_y)*w;
  }
  myRMS_y = ( sqrt(sumdev2_y / sumw_y) );

  cout << "tg_musp_x, myRMS = " << myRMS_x << ", GetRMS() = " << tg_musp_x->GetRMS() << endl;
  cout << "mean_x = " << mean_x << ", GetMean() = " << tg_musp_x->GetMean() << endl;
  cout << "max_x = " << max_x<< ", MaxElement = " << TMath::MaxElement(tg_musp_x->GetN(),tg_musp_x->GetY()) << endl;
  cout << "tg_musp_y, myRMS = " << myRMS_y << ", GetRMS() = " << tg_musp_y->GetRMS() << endl;
  cout << "mean_y = " << mean_y << ", GetMean() = " << tg_musp_y->GetMean() << endl;
  cout << "max_y = " << max_y << ", MaxElement = " << TMath::MaxElement(tg_musp_y->GetN(),tg_musp_y->GetY()) << endl;
 
  ////////////////Hack for April vdM scans, excluding outermost points from the fit
  low_displacement = -127.0;
  high_displacement = 125.0;

  cout << "low_displacement = " << low_displacement << ", high_displacement = " << high_displacement << endl;

  //fit_functions.push_back("sgcl");
  f_x["sgcl"] = new TF1("f_x_sgcl", fit_func_sgc, low_displacement, high_displacement, 4);
  f_x["sgcl"]->SetNpx(10000);

  f_x["sgcl"]->SetParName(0, "mu_max");
  f_x["sgcl"]->SetParName(1, "sigma_a");
  f_x["sgcl"]->SetParName(2, "h0");
  f_x["sgcl"]->SetParName(3, "c");

  f_x["sgcl"]->SetParameter(0, TMath::MaxElement(tg_musp_x->GetN(),tg_musp_x->GetY()));
  f_x["sgcl"]->SetParameter(1, 0.37*tg_musp_x->GetRMS());
  f_x["sgcl"]->SetParameter(2, tg_musp_x->GetMean());
  f_x["sgcl"]->SetParameter(3, 0.);

  //f_x["sgcl"]->SetParLimits(0, mu_max_guess_x / 100., mu_max_guess_x * 100.);
  //f_x["sgcl"]->SetParLimits(1, 10., 500.);
  //f_x["sgcl"]->SetParLimits(2, -100., 100.);
  //f_x["sgcl"]->SetParLimits(3, 0., mu_max_guess_x / 10.);

  f_y["sgcl"] = new TF1("f_y_sgcl", fit_func_sgc, low_displacement, high_displacement, 4);
  f_y["sgcl"]->SetNpx(10000);

  f_y["sgcl"]->SetParName(0, "mu_max");
  f_y["sgcl"]->SetParName(1, "sigma_a");
  f_y["sgcl"]->SetParName(2, "h0");
  f_y["sgcl"]->SetParName(3, "c");

  f_y["sgcl"]->SetParameter(0, TMath::MaxElement(tg_musp_y->GetN(),tg_musp_y->GetY()));
  f_y["sgcl"]->SetParameter(1, 0.37*tg_musp_y->GetRMS());
  f_y["sgcl"]->SetParameter(2, tg_musp_y->GetMean());
  f_y["sgcl"]->SetParameter(3, 0.);
	
  //f_y["sgcl"]->SetParLimits(0, mu_max_guess_y / 100., mu_max_guess_y * 100.);
  //f_y["sgcl"]->SetParLimits(1, 10., 500.);
  //f_y["sgcl"]->SetParLimits(2, -100., 100.);
  //f_y["sgcl"]->SetParLimits(3, 0., mu_max_guess_y / 10.);


  fit_functions.push_back("sgc");
  f_x["sgc"] = new TF1("f_x_sgc", fit_func_sgc, low_displacement, high_displacement, 4);
  f_x["sgc"]->SetNpx(10000);

  f_x["sgc"]->SetParName(0, "mu_max");
  f_x["sgc"]->SetParName(1, "sigma_a");
  f_x["sgc"]->SetParName(2, "h0");
  f_x["sgc"]->SetParName(3, "c");

  f_x["sgc"]->SetParameter(0, max_x);
  f_x["sgc"]->SetParameter(1, myRMS_x);
  f_x["sgc"]->SetParameter(2, mean_x);
  f_x["sgc"]->SetParameter(3, max_x/20.0);

  f_x["sgc"]->SetParLimits(0, 0.7*max_x, 1.3*max_x);
  f_x["sgc"]->SetParLimits(1, 0.7*myRMS_x, 1.3*myRMS_x);
  f_x["sgc"]->SetParLimits(2, mean_x-0.5*myRMS_x, mean_x+0.5*myRMS_x);
  f_x["sgc"]->SetParLimits(3, 0., max_x / 10.);

  f_y["sgc"] = new TF1("f_y_sgc", fit_func_sgc, low_displacement, high_displacement, 4);
  f_y["sgc"]->SetNpx(10000);

  f_y["sgc"]->SetParName(0, "mu_max");
  f_y["sgc"]->SetParName(1, "sigma_a");
  f_y["sgc"]->SetParName(2, "h0");
  f_y["sgc"]->SetParName(3, "c");

  f_y["sgc"]->SetParameter(0, max_y);
  f_y["sgc"]->SetParameter(1, myRMS_y);
  f_y["sgc"]->SetParameter(2, mean_y);
  f_y["sgc"]->SetParameter(3, max_y/20.0);

  f_y["sgc"]->SetParLimits(0, 0.7*max_y, 1.3*max_y);
  f_y["sgc"]->SetParLimits(1, 0.7*myRMS_y, 1.3*myRMS_y);
  f_y["sgc"]->SetParLimits(2, mean_y-0.5*myRMS_y, mean_y+0.5*myRMS_y);
  f_y["sgc"]->SetParLimits(3, 0., max_y / 10.);


  //fit_functions.push_back("sg");
  f_x["sg"] = new TF1("f_x_sg", fit_func_sgc, low_displacement, high_displacement, 4);
  f_x["sg"]->SetNpx(10000);

  f_x["sg"]->SetParName(0, "mu_max");
  f_x["sg"]->SetParName(1, "sigma_a");
  f_x["sg"]->SetParName(2, "h0");
  f_x["sg"]->SetParName(3, "c");

  f_x["sg"]->SetParameter(0, TMath::MaxElement(tg_musp_x->GetN(),tg_musp_x->GetY()));
  f_x["sg"]->SetParameter(1, 0.37*tg_musp_x->GetRMS());
  f_x["sg"]->SetParameter(2, tg_musp_x->GetMean());
  f_x["sg"]->FixParameter(3, 0.);

  //f_x["sg"]->SetParLimits(0, mu_max_guess_x / 100., mu_max_guess_x * 100.);
  //f_x["sg"]->SetParLimits(1, 10., 500.);
  //f_x["sg"]->SetParLimits(2, -100., 100.);

  f_y["sg"] = new TF1("f_y_sg", fit_func_sgc, low_displacement, high_displacement, 4);
  f_y["sg"]->SetNpx(10000);

  f_y["sg"]->SetParName(0, "mu_max");
  f_y["sg"]->SetParName(1, "sigma_a");
  f_y["sg"]->SetParName(2, "h0");
  f_y["sg"]->SetParName(3, "c");

  f_y["sg"]->SetParameter(0, TMath::MaxElement(tg_musp_y->GetN(),tg_musp_y->GetY()));
  f_y["sg"]->SetParameter(1, 0.37*tg_musp_y->GetRMS());
  f_y["sg"]->SetParameter(2, tg_musp_y->GetMean());
  f_y["sg"]->FixParameter(3, 0.);
  
    //f_y["sg"]->SetParLimits(0, mu_max_guess_y / 100., mu_max_guess_y * 100.);
  //f_y["sg"]->SetParLimits(1, 10., 500.);
  //f_y["sg"]->SetParLimits(2, -100., 100.);

  ///// Double Gaussian /////
  //fit_functions.push_back("dg");
  //f_x["dg"] = new TF1("f_x_dg", fit_func_dg, low_displacement, high_displacement, 6);
  f_x["dg"] = new TF1("f_x_dg", fit_func_dg, low_displacement, high_displacement, 6);
  f_x["dg"]->SetNpx(10000);

  f_x["dg"]->SetParName(0, "mu_max");
  f_x["dg"]->SetParName(1, "sigma_a");
  f_x["dg"]->SetParName(2, "h0");
  f_x["dg"]->SetParName(3, "c");
  f_x["dg"]->SetParName(4, "sigma_b");
  f_x["dg"]->SetParName(5, "f_a");

  f_x["dg"]->SetParameter(0, TMath::MaxElement(tg_musp_x->GetN(),tg_musp_x->GetY()));
  f_x["dg"]->SetParameter(1, 0.37*tg_musp_x->GetRMS());
  f_x["dg"]->SetParameter(2, tg_musp_x->GetMean());
  f_x["dg"]->FixParameter(3, 0.);
  f_x["dg"]->SetParameter(4, 0.37*tg_musp_x->GetRMS());
  f_x["dg"]->SetParameter(5, 0.75);

  //f_x["dg"]->SetParLimits(0, mu_max_guess_x / 100., mu_max_guess_x * 100.);
  //f_x["dg"]->SetParLimits(1, 5., 1000.);
  //f_x["dg"]->SetParLimits(2, -100., 100.);
  //f_x["dg"]->SetParLimits(4, 5., 1000.);
  //f_x["dg"]->SetParLimits(5, 0., 1.);

  //f_y["dg"] = new TF1("f_y_dg", fit_func_dg, low_displacement, high_displacement, 6);
  f_y["dg"] = new TF1("f_y_dg", fit_func_dg, low_displacement, high_displacement, 6);
  f_y["dg"]->SetNpx(10000);

  f_y["dg"]->SetParName(0, "mu_max");
  f_y["dg"]->SetParName(1, "sigma_a");
  f_y["dg"]->SetParName(2, "h0");
  f_y["dg"]->SetParName(3, "c");
  f_y["dg"]->SetParName(4, "sigma_b");
  f_y["dg"]->SetParName(5, "f_a");

  f_y["dg"]->SetParameter(0, TMath::MaxElement(tg_musp_y->GetN(),tg_musp_y->GetY()));
  f_y["dg"]->SetParameter(1, 0.37*tg_musp_y->GetRMS()); 
  f_y["dg"]->SetParameter(2, tg_musp_y->GetMean());
  f_y["dg"]->FixParameter(3, 0.);
  f_y["dg"]->SetParameter(4, 0.37*tg_musp_y->GetRMS()); 
  f_y["dg"]->SetParameter(5, 0.75);

  //f_y["dg"]->SetParLimits(0, mu_max_guess_y / 100., mu_max_guess_y * 100.);
  //f_y["dg"]->SetParLimits(1, 5.,1000.);
  //f_y["dg"]->SetParLimits(2, -100., 100.);
  //f_y["dg"]->SetParLimits(4, 5., 1000.);
  //f_y["dg"]->SetParLimits(5, 0., 1.);

  ///// Double Gaussian plus constant (not luminosity) /////
  fit_functions.push_back("dgc");
  f_x["dgc"] = new TF1("f_x_dgc", fit_func_dgc, low_displacement, high_displacement, 6);
  f_x["dgc"]->SetNpx(10000);

  f_x["dgc"]->SetParName(0, "mu_max");
  //f_x["dgc"]->SetParName(1, "alpha");
  f_x["dgc"]->SetParName(1, "sigma_a");
  f_x["dgc"]->SetParName(2, "h0");
  f_x["dgc"]->SetParName(3, "c");
  f_x["dgc"]->SetParName(4, "Sigma_x");
  f_x["dgc"]->SetParName(5, "f_a");
  
  //James
  f_x["dgc"]->SetParameter(0, max_x);
  //f_x["dgc"]->SetParameter(1, 0.5);
  //f_x["dgc"]->SetParameter(1, myRMS_x*2.0);
  f_x["dgc"]->SetParameter(1, myRMS_x);
  f_x["dgc"]->SetParameter(2, mean_x);
  //f_x["dgc"]->SetParameter(3, max_x/20.0);
  f_x["dgc"]->SetParameter(3, 0.0);
  f_x["dgc"]->SetParameter(4, myRMS_x);
  //f_x["dgc"]->SetParameter(5, 0.2);
  f_x["dgc"]->SetParameter(5, 0.5);

  //OLD
  /*f_x["dgc"]->SetParameter(0, TMath::MaxElement(tg_musp_x->GetN(),tg_musp_x->GetY()));
  f_x["dgc"]->SetParameter(1, 0.37*tg_musp_x->GetRMS());
  f_x["dgc"]->SetParameter(2, tg_musp_x->GetMean());
  f_x["dgc"]->SetParameter(3, 0.);
  f_x["dgc"]->SetParameter(4, 0.37*tg_musp_x->GetRMS());
  f_x["dgc"]->SetParameter(5, 0.75);

  //James
  /*f_x["dgc"]->SetParLimits(0, 0.7*max_x, 1.3*max_x);
  //f_x["dgc"]->SetParLimits(1, 0.0, 1.0);
  //f_x["dgc"]->SetParLimits(1, 0.7*myRMS_x, 14.0*myRMS_x);
  f_x["dgc"]->SetParLimits(1, 0.7*myRMS_x, 1.3*myRMS_x);
  f_x["dgc"]->SetParLimits(2, mean_x-0.5*myRMS_x, mean_x+0.5*myRMS_x);
  f_x["dgc"]->SetParLimits(3, -1.*max_x/10.0, max_x/10.0);
  f_x["dgc"]->SetParLimits(4, 0.7*myRMS_x, 1.3*myRMS_x);
  f_x["dgc"]->SetParLimits(5, 0.0, 1.0);*/

  //OLD
  /*f_x["dgc"]->SetParLimits(1, 0., 10*tg_musp_x->GetRMS());
  //f_x["dgc"]->SetParLimits(1, 25/100, 25*100);
  //f_x["dgc"]->SetParLimits(2, -100., 100.);
  //f_x["dgc"]->SetParLimits(3, 0., 3.e-34);
  f_x["dgc"]->SetParLimits(4, 0., 10*tg_musp_x->GetRMS());
  //f_x["dgc"]->SetParLimits(4, 25/100, 25*100);
  f_x["dgc"]->SetParLimits(5, 0.001, 0.999);*/

  f_x["dgc"]->SetParLimits(0, 0.2*max_x, 5.0*max_x);
  //f_x["dgc"]->SetParLimits(1, 0.0, 1.0);
  //f_x["dgc"]->SetParLimits(1, 0.7*myRMS_x, 14.0*myRMS_x);
  f_x["dgc"]->SetParLimits(1, 0.1*myRMS_x, 10.0*myRMS_x);
  f_x["dgc"]->SetParLimits(2, mean_x-0.5*myRMS_x, mean_x+0.5*myRMS_x);
  f_x["dgc"]->SetParLimits(3, -1.*max_x/10.0, max_x/10.0);
  f_x["dgc"]->SetParLimits(4, 0.1*myRMS_x, 10.0*myRMS_x);
  f_x["dgc"]->SetParLimits(5, 0.0, 1.0);

  f_y["dgc"] = new TF1("f_y_dgc", fit_func_dgc, low_displacement, high_displacement, 6);
  f_y["dgc"]->SetNpx(10000);

  f_y["dgc"]->SetParName(0, "mu_max");
  //f_y["dgc"]->SetParName(1, "alpha");
  f_y["dgc"]->SetParName(1, "sigma_a");
  f_y["dgc"]->SetParName(2, "h0");
  f_y["dgc"]->SetParName(3, "c");
  f_y["dgc"]->SetParName(4, "Sigma_y");
  f_y["dgc"]->SetParName(5, "f_a");

  //James
  f_y["dgc"]->SetParameter(0, max_y);
  //f_y["dgc"]->SetParameter(1, 0.5);
  //f_y["dgc"]->SetParameter(1, myRMS_y*2.0);
  f_y["dgc"]->SetParameter(1, myRMS_y);
  f_y["dgc"]->SetParameter(2, mean_y);
  //f_y["dgc"]->SetParameter(3, max_y/20.0);
  f_y["dgc"]->SetParameter(3, 0.0);
  f_y["dgc"]->SetParameter(4, myRMS_y);
  //f_y["dgc"]->SetParameter(5, 0.2);
  f_y["dgc"]->SetParameter(5, 0.5);

  //OLD
  /*f_y["dgc"]->SetParameter(0, TMath::MaxElement(tg_musp_y->GetN(),tg_musp_y->GetY()));
  f_y["dgc"]->SetParameter(1, 0.37*tg_musp_y->GetRMS());
  f_y["dgc"]->SetParameter(2, tg_musp_y->GetMean());
  f_y["dgc"]->SetParameter(3, 0.);
  f_y["dgc"]->SetParameter(4, 0.37*tg_musp_y->GetRMS());
  f_y["dgc"]->SetParameter(5, 0.75);*/

  //James
  /*f_y["dgc"]->SetParLimits(0, 0.7*max_y, 1.3*max_y);
  //f_y["dgc"]->SetParLimits(1, 0.0, 1.0);
  f_y["dgc"]->SetParLimits(1, 0.7*myRMS_y, 14.0*myRMS_y);
  f_y["dgc"]->SetParLimits(2, mean_y-0.5*myRMS_y, mean_y+0.5*myRMS_y);
  f_y["dgc"]->SetParLimits(3, -1.*max_y/10.0, max_y/10.0);
  f_y["dgc"]->SetParLimits(4, 0.7*myRMS_y, 1.3*myRMS_y);
  f_y["dgc"]->SetParLimits(5, 0.0, 1.0);*/

  //OLD
  /*f_y["dgc"]->SetParLimits(1, 0., 10*tg_musp_y->GetRMS());
  //f_y["dgc"]->SetParLimits(1, 25/100, 25*100);
  //f_y["dgc"]->SetParLimits(2, -100., 100.);
  //f_y["dgc"]->SetParLimits(3, 0., 3.e-34);
  f_y["dgc"]->SetParLimits(4, 0., 10*tg_musp_y->GetRMS());
  //f_y["dgc"]->SetParLimits(4, 25/100, 25*100);
  f_y["dgc"]->SetParLimits(5, 0.001, 0.999);*/

  f_y["dgc"]->SetParLimits(0, 0.2*max_y, 5.0*max_y);
  //f_y["dgc"]->SetParLimits(1, 0.0, 1.0);
  //f_y["dgc"]->SetParLimits(1, 0.7*myRMS_y, 14.0*myRMS_y);
  f_y["dgc"]->SetParLimits(1, 0.1*myRMS_y, 10.0*myRMS_y);
  f_y["dgc"]->SetParLimits(2, mean_y-0.5*myRMS_y, mean_y+0.5*myRMS_y);
  f_y["dgc"]->SetParLimits(3, -1.*max_y/10.0, max_y/10.0);
  f_y["dgc"]->SetParLimits(4, 0.1*myRMS_y, 10.0*myRMS_y);
  f_y["dgc"]->SetParLimits(5, 0.0, 1.0);

  ///// Double Gaussian plus constant (luminosity)/////
  //fit_functions.push_back("dgcl");
  f_x["dgcl"] = new TF1("f_x_dgcl", fit_func_dgcl, low_displacement, high_displacement, 6);
  f_x["dgcl"]->SetNpx(10000);

  f_x["dgcl"]->SetParName(0, "mu_max");
  f_x["dgcl"]->SetParName(1, "sigma_a");
  f_x["dgcl"]->SetParName(2, "h0");
  f_x["dgcl"]->SetParName(3, "c");
  f_x["dgcl"]->SetParName(4, "sigma_b");
  f_x["dgcl"]->SetParName(5, "f_a");

  f_x["dgcl"]->SetParameter(0, TMath::MaxElement(tg_musp_x->GetN(),tg_musp_x->GetY()));
  f_x["dgcl"]->SetParameter(1, 0.37*tg_musp_x->GetRMS());
  f_x["dgcl"]->SetParameter(2, tg_musp_x->GetMean());
  f_x["dgcl"]->SetParameter(3, 0.);
  f_x["dgcl"]->SetParameter(4, 0.37*tg_musp_x->GetRMS());
  f_x["dgcl"]->SetParameter(5, 0.75);

  //f_x["dgcl"]->SetParLimits(0, mu_max_guess_x / 100., mu_max_guess_x * 100.);
  //f_x["dgcl"]->SetParLimits(1, 5., 1000.);
  //f_x["dgcl"]->SetParLimits(2, -100., 100.);
  //f_x["dgcl"]->SetParLimits(3, 0., mu_max_guess_x * 100.);
  //f_x["dgcl"]->SetParLimits(4, 5., 1000.);
  //f_x["dgcl"]->SetParLimits(5, 0., 1.);

  f_y["dgcl"] = new TF1("f_y_dgcl", fit_func_dgcl, low_displacement, high_displacement, 6);
  f_y["dgcl"]->SetNpx(10000);

  f_y["dgcl"]->SetParName(0, "mu_max");
  f_y["dgcl"]->SetParName(1, "sigma_a");
  f_y["dgcl"]->SetParName(2, "h0");
  f_y["dgcl"]->SetParName(3, "c");
  f_y["dgcl"]->SetParName(4, "sigma_b");
  f_y["dgcl"]->SetParName(5, "f_a");

  f_y["dgcl"]->SetParameter(0, TMath::MaxElement(tg_musp_y->GetN(),tg_musp_y->GetY()));
  f_y["dgcl"]->SetParameter(1, 0.37*tg_musp_y->GetRMS());
  f_y["dgcl"]->SetParameter(2, tg_musp_y->GetMean());
  f_y["dgcl"]->SetParameter(3, 0.);
  f_y["dgcl"]->SetParameter(4, 0.37*tg_musp_y->GetRMS());
  f_y["dgcl"]->SetParameter(5, 0.75);

  //f_y["dgcl"]->SetParLimits(0, mu_max_guess_y / 100., mu_max_guess_y * 100.);
  //f_y["dgcl"]->SetParLimits(1, 5.,1000.);
  //f_y["dgcl"]->SetParLimits(2, -100., 100.);
  //f_y["dgcl"]->SetParLimits(3, 0., mu_max_guess_y * 100.);
  //f_y["dgcl"]->SetParLimits(4, 5., 1000.);
  //f_y["dgcl"]->SetParLimits(5, 0., 1.);

  //fit_functions.push_back("spline");
  fit_status["spline"] = true;
  fit_spline_x = new TSpline3("fit_spline_x", tg_musp_x);
  fit_spline_y = new TSpline3("fit_spline_y", tg_musp_y);

  // -- Perform fits
  for (std::vector<TString>::iterator fit_function = fit_functions.begin(); fit_function != fit_functions.end(); ++fit_function) {
    if (*fit_function != "spline") {
      //fr_x_actual[*fit_function] = *(tg_musp_x->Fit(f_x[*fit_function], "SQER0+"));
      //fr_y_actual[*fit_function] = *(tg_musp_y->Fit(f_y[*fit_function], "SQER0+"));
      //fr_x[*fit_function] = &fr_x_actual[*fit_function];
      //fr_y[*fit_function] = &fr_y_actual[*fit_function];
      TString fit_options = "SEMR+";
      cout << "[VDMA] Performing fit on x with fit function " << *fit_function << endl;
      fr_x[*fit_function] = tg_musp_x->Fit(f_x[*fit_function], fit_options);
      cout << "[VDMA] Performing fit on y with fit function " << *fit_function << endl;
      fr_y[*fit_function] = tg_musp_y->Fit(f_y[*fit_function], fit_options);
      //if (!(int)fr_x[*fit_function] && !(int)fr_y[*fit_function]) {
      fit_status[*fit_function] = true;
      //} else {
      //  fit_status[*fit_function] = false;
      //}
      TString name_x = f_x[*fit_function]->GetName();
      TString name_y = f_y[*fit_function]->GetName();
      r_x[*fit_function] = tg_musp_x->GetFunction(name_x)->Eval(fr_x[*fit_function]->Parameter(2)) / tg_musp_x->GetFunction(name_x)->Eval(0.);
      cout << "[VDMA] r_x[" << *fit_function<< "] = " << r_x[*fit_function] << endl;
      r_y[*fit_function] = tg_musp_y->GetFunction(name_y)->Eval(fr_y[*fit_function]->Parameter(2)) / tg_musp_y->GetFunction(name_y)->Eval(0.);
      cout << "[VDMA] r_y[" << *fit_function<< "] = " << r_y[*fit_function] << endl;
    }
  }
}


void VanDerMeerAnalysis::Finalize() {

  cout << "[VanDerMeerAnalysis] INFO : Calculating sigma_vis values" << endl;

  // -- Single gaussian + constant(L)
//  mu_max_x["sgcl"] = fr_x["sgcl"]->Parameter(0) + fr_x["sgcl"]->Parameter(3);
//  Sigma_x["sgcl"] = fr_x["sgcl"]->Parameter(1);
//  c_x["sgcl"] = fr_x["sgcl"]->Parameter(3);
//  mu_max_y["sgcl"] = fr_y["sgcl"]->Parameter(0) + fr_y["sgcl"]->Parameter(3);
//  Sigma_y["sgcl"] = fr_y["sgcl"]->Parameter(1);
//  c_y["sgcl"] = fr_y["sgcl"]->Parameter(3);//

//  mu_max_x_err["sgcl"] = TMath::Sqrt(TMath::Power(fr_x["sgcl"]->ParError(0), 2) + TMath::Power(fr_x["sgcl"]->ParError(3), 2));
//  Sigma_x_err["sgcl"] = fr_x["sgcl"]->ParError(1);
//  c_x_err["sgcl"] = fr_x["sgcl"]->ParError(3);
//  mu_max_y_err["sgcl"] = TMath::Sqrt(TMath::Power(fr_y["sgcl"]->ParError(0), 2) + TMath::Power(fr_y["sgcl"]->ParError(3), 2));
//  Sigma_y_err["sgcl"] = fr_y["sgcl"]->ParError(1);
//  c_y_err["sgcl"] = fr_y["sgcl"]->ParError(3);

  // -- Single gaussian + constant(B)
  mu_max_x["sgc"] = fr_x["sgc"]->Parameter(0);
  Sigma_x["sgc"] = fr_x["sgc"]->Parameter(1);
  c_x["sgc"] = fr_x["sgc"]->Parameter(3);
  mu_max_y["sgc"] = fr_y["sgc"]->Parameter(0);
  Sigma_y["sgc"] = fr_y["sgc"]->Parameter(1);
  c_y["sgc"] = fr_y["sgc"]->Parameter(3);

  mu_max_x_err["sgc"] = fr_x["sgc"]->ParError(0);
  Sigma_x_err["sgc"] = fr_x["sgc"]->ParError(1);
  c_x_err["sgc"] = fr_x["sgc"]->ParError(3);
  mu_max_y_err["sgc"] = fr_y["sgc"]->ParError(0);
  Sigma_y_err["sgc"] = fr_y["sgc"]->ParError(1);
  c_y_err["sgc"] = fr_y["sgc"]->ParError(3);

  // -- Single gaussian
//  mu_max_x["sg"] = fr_x["sg"]->Parameter(0);
//  Sigma_x["sg"] = fr_x["sg"]->Parameter(1);
//  c_x["sg"] = 0.;
//  mu_max_y["sg"] = fr_y["sg"]->Parameter(0);
//  Sigma_y["sg"] = fr_y["sg"]->Parameter(1);
//  c_y["sg"] = 0.;//

//  mu_max_x_err["sg"] = fr_x["sg"]->ParError(0);
//  Sigma_x_err["sg"] = fr_x["sg"]->ParError(1);
//  c_x_err["sg"] = 0.;
//  mu_max_y_err["sg"] = fr_y["sg"]->ParError(0);
//  Sigma_y_err["sg"] = fr_y["sg"]->ParError(1);
//  c_y_err["sg"] = 0.;

  // -- Double gaussian
//  mu_max_x["dg"] = fr_x["dg"]->Parameter(0);
//  Double_t f_a_x = fr_x["dg"]->Parameter(5);
//  Sigma_x["dg"] = 1. / (f_a_x / fr_x["dg"]->Parameter(1) + (1. - f_a_x) / fr_x["dg"]->Parameter(4));
//  c_x["dg"] = 0.;
//  mu_max_y["dg"] = fr_y["dg"]->Parameter(0);
//  Double_t f_a_y = fr_y["dg"]->Parameter(5);
//  Sigma_y["dg"] = 1. / (f_a_y / fr_y["dg"]->Parameter(1) + (1. - f_a_y) / fr_y["dg"]->Parameter(4));
//  c_y["dg"] = 0.;//

//  mu_max_x_err["dg"] = fr_x["dg"]->ParError(0);
//  Sigma_x_err["dg"] = Sigma_x["dg"] * TMath::Sqrt(
//                        TMath::Power(f_a_x / TMath::Power(fr_x["dg"]->Parameter(1), 2) * fr_x["dg"]->ParError(1), 2) +
//                        TMath::Power((1. - f_a_x) / TMath::Power(fr_x["dg"]->Parameter(4), 2) * fr_x["dg"]->ParError(4), 2) +
//                        TMath::Power((1./fr_x["dg"]->Parameter(1) - 1./fr_x["dg"]->Parameter(4)) * fr_x["dg"]->ParError(5), 2));
//  c_x_err["dg"] = 0.;
//  mu_max_y_err["dg"] = fr_y["dg"]->ParError(0);
//  Sigma_y_err["dg"] = Sigma_y["dg"] * TMath::Sqrt(
//                        TMath::Power(f_a_y / TMath::Power(fr_y["dg"]->Parameter(1), 2) * fr_y["dg"]->ParError(1), 2) +
//                        TMath::Power((1. - f_a_y) / TMath::Power(fr_y["dg"]->Parameter(4), 2) * fr_y["dg"]->ParError(4), 2) +
//                        TMath::Power((1./fr_y["dg"]->Parameter(1) - 1./fr_y["dg"]->Parameter(4)) * fr_y["dg"]->ParError(5), 2));
//  c_y_err["dg"] = 0.;

  // -- Double gaussian plus constant (not luminosity)
  mu_max_x["dgc"] = fr_x["dgc"]->Parameter(0);
  Double_t f_a_x = fr_x["dgc"]->Parameter(5);
  Sigma_x["dgc"] = fr_x["dgc"]->Parameter(4);
  c_x["dgc"] = fr_x["dgc"]->Parameter(3);
  mu_max_y["dgc"] = fr_y["dgc"]->Parameter(0);
  Double_t f_a_y = fr_y["dgc"]->Parameter(5);
  Sigma_y["dgc"] = fr_y["dgc"]->Parameter(4);
  c_y["dgc"] = fr_y["dgc"]->Parameter(3);

  mu_max_x_err["dgc"] = fr_x["dgc"]->ParError(0);
  Sigma_x_err["dgc"] = fr_x["dgc"]->ParError(4);
  c_x_err["dgc"] = fr_x["dgc"]->ParError(3);
  mu_max_y_err["dgc"] = fr_y["dgc"]->ParError(0);
  Sigma_y_err["dgc"] = fr_y["dgc"]->ParError(4);
  c_y_err["dgc"] = fr_y["dgc"]->ParError(3);

//  // -- Double gaussian plus constant (luminosity)
//  mu_max_x["dgcl"] = fr_x["dgcl"]->Parameter(0) + fr_x["dgcl"]->Parameter(3);
//  f_a_x = fr_x["dgcl"]->Parameter(5);
//  Sigma_x["dgcl"] = 1. / (f_a_x / fr_x["dgcl"]->Parameter(1) + (1. - f_a_x) / fr_x["dgcl"]->Parameter(4));
//  c_x["dgcl"] = fr_x["dgcl"]->Parameter(3);
//  mu_max_y["dgcl"] = fr_y["dgcl"]->Parameter(0) + fr_y["dgcl"]->Parameter(3);
//  f_a_y = fr_y["dgcl"]->Parameter(5);
//  Sigma_y["dgcl"] = 1. / (f_a_y / fr_y["dgcl"]->Parameter(1) + (1. - f_a_y) / fr_y["dgcl"]->Parameter(4));
//  c_y["dgcl"] = fr_y["dgcl"]->Parameter(3);//

//  mu_max_x_err["dgcl"] = TMath::Sqrt(TMath::Power(fr_x["dgcl"]->ParError(0), 2) + TMath::Power(fr_x["dgcl"]->ParError(3), 2));
//  Sigma_x_err["dgcl"] = Sigma_x["dgcl"] * TMath::Sqrt(
//                          TMath::Power(f_a_x / TMath::Power(fr_x["dgcl"]->Parameter(1), 2) * fr_x["dgcl"]->ParError(1), 2) +
//                          TMath::Power((1. - f_a_x) / TMath::Power(fr_x["dgcl"]->Parameter(4), 2) * fr_x["dgcl"]->ParError(4), 2) +
//                          TMath::Power((1./fr_x["dgcl"]->Parameter(1) - 1./fr_x["dgcl"]->Parameter(4)) * fr_x["dgcl"]->ParError(5), 2));
//  c_x_err["dgcl"] = fr_x["dgcl"]->ParError(3);
//  mu_max_y_err["dgcl"] = TMath::Sqrt(TMath::Power(fr_y["dgcl"]->ParError(0), 2) + TMath::Power(fr_y["dgcl"]->ParError(3), 2));
//  Sigma_y_err["dgcl"] = Sigma_y["dgcl"] * TMath::Sqrt(
//                          TMath::Power(f_a_y / TMath::Power(fr_y["dgcl"]->Parameter(1), 2) * fr_y["dgcl"]->ParError(1), 2) +
//                          TMath::Power((1. - f_a_y) / TMath::Power(fr_y["dgcl"]->Parameter(4), 2) * fr_y["dgcl"]->ParError(4), 2) +
//                          TMath::Power((1./fr_y["dgcl"]->Parameter(1) - 1./fr_y["dgcl"]->Parameter(4)) * fr_y["dgcl"]->ParError(5), 2));
//  c_y_err["dgcl"] = fr_y["dgcl"]->ParError(3);//
//

//  // -- Spline
//  Double_t spline_integral_x = CalculateTSplineIntegral(fit_spline_x, low_displacement, high_displacement, 10000);
//  mu_max_x["spline"] = CalculateTSplineMax(fit_spline_x, low_displacement, high_displacement);
//  Sigma_x["spline"] = spline_integral_x / TMath::Sqrt(2 * TMath::Pi()) / mu_max_x["spline"];
//  c_x["spline"] = 0.;//

//  Double_t spline_integral_y = CalculateTSplineIntegral(fit_spline_y, low_displacement, high_displacement, 10000);
//  mu_max_y["spline"] = CalculateTSplineMax(fit_spline_y, low_displacement, high_displacement);
//  Sigma_y["spline"] = spline_integral_y / TMath::Sqrt(2 * TMath::Pi()) / mu_max_y["spline"];
//  c_y["spline"] = 0.;//

//  mu_max_x_err["spline"] = 0.;
//  Sigma_x_err["spline"] = 0.;
//  c_x_err["spline"] = 0.;
//  mu_max_y_err["spline"] = 0.;
//  Sigma_y_err["spline"] = 0.;
//  c_y_err["spline"] = 0.;

  for (vector<TString>::iterator fit_function = fit_functions.begin(); fit_function != fit_functions.end(); ++fit_function) {
    lumi_sp[*fit_function] = 11245.5 / (2 * TMath::Pi() * Sigma_x[*fit_function] * Sigma_y[*fit_function]) * TMath::Power(10., 22) * TMath::Power(10, -19);
    sigma_vis[*fit_function] = TMath::Pi() * (mu_max_x[*fit_function]*r_y[*fit_function] + mu_max_y[*fit_function]*r_x[*fit_function]) * Sigma_x[*fit_function] * Sigma_y[*fit_function] * TMath::Power(10, 19);
    lumi[*fit_function] = lumi_sp[*fit_function] * (bunch_intensities_1[plb_list_x[12]] + bunch_intensities_1[plb_list_y[12]]) * (bunch_intensities_2[plb_list_x[12]] + bunch_intensities_2[plb_list_y[12]]) / 4.;

    cout << endl;
    cout << "[VanDerMeerAnalysis] INFO : Results for fit function " << *fit_function << ", N" << vtx_method << ntrkcut << ":" << endl;
    if (fit_status[*fit_function]) {
      cout << "[VanDerMeerAnalysis] WARNING : Fit failed!!!" << endl;
    }
    if (*fit_function != "spline") {
      cout << "[VanDerMeerAnalysis] INFO : X chi2/ndf = " << fr_x[*fit_function]->Chi2() / fr_x[*fit_function]->Ndf() << endl;
    }
    cout << "[VanDerMeerAnalysis] INFO : mu_max_x = " << mu_max_x[*fit_function] << endl;
    cout << "[VanDerMeerAnalysis] INFO : Sigma_x = " << Sigma_x[*fit_function] << endl;
    cout << "[VanDerMeerAnalysis] INFO : c_x = " << c_x[*fit_function] << endl;
    if (*fit_function == "dg" ) {
      cout << "[VanDerMeerAnalysis] INFO : f_a_x = " << f_a_x << endl;
      cout << "[VanDerMeerAnalysis] INFO : sigma_a = " << fr_x["dg"]->Parameter(1) << endl;
      cout << "[VanDerMeerAnalysis] INFO : sigma_b = " << fr_x["dg"]->Parameter(4) << endl;
    }
    if (*fit_function == "dgc" ) {
      cout << "[VanDerMeerAnalysis] INFO : mu_max_x_err = " << mu_max_x_err["dgc"] << endl;
      cout << "[VanDerMeerAnalysis] INFO : Sigma_x_err = " << Sigma_x_err["dgc"] << endl;
    }
    if (*fit_function != "spline") {
      cout << "[VanDerMeerAnalysis] INFO : Y chi2/ndf = " << fr_y[*fit_function]->Chi2() / fr_y[*fit_function]->Ndf() << endl;
    }
    cout << "[VanDerMeerAnalysis] INFO : mu_max_y = " << mu_max_y[*fit_function] << endl;
    cout << "[VanDerMeerAnalysis] INFO : Sigma_y = " << Sigma_y[*fit_function] << endl;
    cout << "[VanDerMeerAnalysis] INFO : c_y = " << c_y[*fit_function] << endl;
    if (*fit_function == "dg") {
      cout << "[VanDerMeerAnalysis] INFO : f_a_y = " << f_a_y << endl;
      cout << "[VanDerMeerAnalysis] INFO : sigma_a = " << fr_y["dg"]->Parameter(1) << endl;
      cout << "[VanDerMeerAnalysis] INFO : sigma_b = " << fr_y["dg"]->Parameter(4) << endl;
    }
    if (*fit_function == "dgc" ) {
      cout << "[VanDerMeerAnalysis] INFO : mu_max_y_err = " << mu_max_y_err["dgc"] << endl;
      cout << "[VanDerMeerAnalysis] INFO : Sigma_y_err = " << Sigma_y_err["dgc"] << endl;
    }

    cout << "[VanDerMeerAnalysis] INFO : ----------------------------------------" << endl;
    cout << "[VanDerMeerAnalysis] INFO : sigma_vis = " << sigma_vis[*fit_function] << endl;
    cout << "[VanDerMeerAnalysis] INFO : lumi_sp = " << lumi_sp[*fit_function] << endl;
    cout << endl;
  }

  // -- ERRORS UGH
  
  //NOTE: After applying the centring correction, I only applied this to the two functions I'm currently interested in dgc and sgc
  for (vector<TString>::iterator fit_function = fit_functions.begin(); fit_function != fit_functions.end(); ++fit_function) {
    cout << "Calculating errors for fit function: " << *fit_function << endl;
    if (*fit_function == "spline") { 
      continue;
    }
    if (!fit_status[*fit_function]) {
      continue;
    }

    Int_t npar;
    if (*fit_function == "sgcl" || *fit_function == "sgc" || *fit_function == "sg") {
      npar = 4;
    } else if (*fit_function == "dgcl" || *fit_function == "dgc" || *fit_function == "dg") {
      npar = 6;
    }
    if (*fit_function != "dg") {
      cout << "GetCovarianceMatrix for fit function: " << *fit_function << ", parameters: " << npar << endl;
      fit_cov_x[*fit_function] = new TMatrixDSym(npar);
      *fit_cov_x[*fit_function] = fr_x[*fit_function]->GetCovarianceMatrix();

      fit_cov_y[*fit_function] = new TMatrixDSym(npar);
      *fit_cov_y[*fit_function] = fr_y[*fit_function]->GetCovarianceMatrix();
    }
  }
//  // -- SGCL
//  Double_t dsigma_dpx_4[] = {
//    TMath::Pi() * Sigma_x["sgcl"] * Sigma_y["sgcl"], //mu
//    TMath::Pi() * (mu_max_x["sgcl"] + mu_max_y["sgcl"]) * Sigma_y["sgcl"], //sigma
//    0., //center
//    TMath::Pi() * Sigma_x["sgcl"] * Sigma_y["sgcl"]//constant
//  };
//  Double_t dsigma_dpy_4[] = {
//    TMath::Pi() * Sigma_x["sgcl"] * Sigma_y["sgcl"], //mu
//    TMath::Pi() * (mu_max_x["sgcl"] + mu_max_y["sgcl"]) * Sigma_x["sgcl"], //sigma
//    0., //center
//    TMath::Pi() * Sigma_x["sgcl"] * Sigma_y["sgcl"]//constant
//  };
//  Double_t dsigma_dpx_T_4[] = {
//    TMath::Pi() * Sigma_x["sgcl"] * Sigma_y["sgcl"], //mu
//    TMath::Pi() * (mu_max_x["sgcl"] + mu_max_y["sgcl"]) * Sigma_y["sgcl"], //sigma
//    0., //center
//    TMath::Pi() * Sigma_x["sgcl"] * Sigma_y["sgcl"]//constant
//  };
//  Double_t dsigma_dpy_T_4[] = {
//    TMath::Pi() * Sigma_x["sgcl"] * Sigma_y["sgcl"], //mu
//    TMath::Pi() * (mu_max_x["sgcl"] + mu_max_y["sgcl"]) * Sigma_x["sgcl"], //sigma
//    0., //center
//    TMath::Pi() * Sigma_x["sgcl"] * Sigma_y["sgcl"]//constant
//  };//

//  if (fit_status["sgcl"]) {
//    v_dsigma_dpx["sgcl"] = new TMatrixD(4, 1);
//    v_dsigma_dpx_T["sgcl"] = new TMatrixD(1, 4);
//    v_dsigma_dpx["sgcl"]->Use(4, 1, dsigma_dpx_4);
//    v_dsigma_dpx_T["sgcl"]->Use(1, 4, dsigma_dpx_T_4);
//    v_dsigma_dpy["sgcl"] = new TMatrixD(4, 1);
//    v_dsigma_dpy_T["sgcl"] = new TMatrixD(1, 4);
//    v_dsigma_dpy["sgcl"]->Use(4, 1, dsigma_dpy_4);
//    v_dsigma_dpy_T["sgcl"]->Use(1, 4, dsigma_dpy_T_4);//

//    *v_dsigma_dpx["sgcl"] = *fit_cov_x["sgcl"] * *v_dsigma_dpx["sgcl"];
//    *v_dsigma_dpy["sgcl"] = *fit_cov_y["sgcl"] * *v_dsigma_dpy["sgcl"];//

//    sigma_vis_err["sgcl"] = TMath::Power(10, 19) * TMath::Sqrt(
//                              ((*v_dsigma_dpx_T["sgcl"]) * (*v_dsigma_dpx["sgcl"]))(0, 0) +
//                              ((*v_dsigma_dpy_T["sgcl"]) * (*v_dsigma_dpy["sgcl"]))(0, 0)
//                            );//

//  } else {
//    sigma_vis_err["sgcl"] = 0.;
//  }

  // -- SGC
  if (fit_status["sgc"]) {
    Double_t dsigma_dpx_4[] = {}; 
    Double_t dsigma_dpy_4[] = {}; 
    dsigma_dpx_4[0] = TMath::Pi() * Sigma_x["sgc"] * Sigma_y["sgc"] * r_y["sgc"];
    dsigma_dpx_4[1] = TMath::Pi() * (mu_max_x["sgc"]*r_y["sgc"] + mu_max_y["sgc"]*r_x["sgc"]) * Sigma_y["sgc"];
    dsigma_dpx_4[2] = 0.;
    dsigma_dpx_4[3] = 0.;
    dsigma_dpy_4[0] = TMath::Pi() * Sigma_x["sgc"] * Sigma_y["sgc"] * r_x["sgc"];
    dsigma_dpy_4[1] = TMath::Pi() * (mu_max_x["sgc"]*r_y["sgc"] + mu_max_y["sgc"]*r_x["sgc"]) * Sigma_x["sgc"];
    dsigma_dpy_4[2] = 0.;
    dsigma_dpy_4[3] = 0.;

    Double_t dsigma_dpx_T_4[] = {};
    Double_t dsigma_dpy_T_4[] = {};
    dsigma_dpx_T_4[0] = TMath::Pi() * Sigma_x["sgc"] * Sigma_y["sgc"] * r_y["sgc"];
    dsigma_dpx_T_4[1] = TMath::Pi() * (mu_max_x["sgc"]*r_y["sgc"] + mu_max_y["sgc"]*r_x["sgc"]) * Sigma_y["sgc"];
    dsigma_dpx_T_4[2] = 0.;
    dsigma_dpx_T_4[3] = 0.;
    dsigma_dpy_T_4[0] = TMath::Pi() * Sigma_x["sgc"] * Sigma_y["sgc"] * r_x["sgc"];
    dsigma_dpy_T_4[1] = TMath::Pi() * (mu_max_x["sgc"]*r_y["sgc"] + mu_max_y["sgc"]*r_x["sgc"]) * Sigma_x["sgc"];
    dsigma_dpy_T_4[2] = 0.;
    dsigma_dpy_T_4[3] = 0.;

    v_dsigma_dpx["sgc"] = new TMatrixD(4, 1);
    v_dsigma_dpx_T["sgc"] = new TMatrixD(1, 4);
    v_dsigma_dpx["sgc"]->Use(4, 1, dsigma_dpx_4);
    v_dsigma_dpx_T["sgc"]->Use(1, 4, dsigma_dpx_T_4);
    v_dsigma_dpy["sgc"] = new TMatrixD(4, 1);
    v_dsigma_dpy_T["sgc"] = new TMatrixD(1, 4);
    v_dsigma_dpy["sgc"]->Use(4, 1, dsigma_dpy_4);
    v_dsigma_dpy_T["sgc"]->Use(1, 4, dsigma_dpy_T_4);

    *v_dsigma_dpx["sgc"] = *fit_cov_x["sgc"] * *(v_dsigma_dpx["sgc"]);
    *v_dsigma_dpy["sgc"] = *fit_cov_y["sgc"] * *(v_dsigma_dpy["sgc"]);

    sigma_vis_err["sgc"] = TMath::Power(10, 19) * TMath::Sqrt(
                             ((*v_dsigma_dpx_T["sgc"]) * (*v_dsigma_dpx["sgc"]))(0, 0) +
                             ((*v_dsigma_dpy_T["sgc"]) * (*v_dsigma_dpy["sgc"]))(0, 0)
                           );
  } else {
    sigma_vis_err["sgc"] = 0.;
  }

//  // -- SG
//  if (fit_status["sg"]) {
//    dsigma_dpx_4[0] = TMath::Pi() * Sigma_x["sg"] * Sigma_y["sg"];
//    dsigma_dpx_4[1] = TMath::Pi() * (mu_max_x["sg"] + mu_max_y["sg"]) * Sigma_y["sg"];
//    dsigma_dpx_4[2] = 0.;
//    dsigma_dpx_4[3] = -1. * TMath::Pi() * Sigma_x["sg"] * Sigma_y["sg"];
//    dsigma_dpy_4[0] = TMath::Pi() * Sigma_x["sg"] * Sigma_y["sg"];
//    dsigma_dpy_4[1] = TMath::Pi() * (mu_max_x["sg"] + mu_max_y["sg"]) * Sigma_x["sg"];
//    dsigma_dpy_4[2] = 0.;
//    dsigma_dpy_4[3] = -1. * TMath::Pi() * Sigma_x["sg"] * Sigma_y["sg"];//

//    dsigma_dpx_T_4[0] = TMath::Pi() * Sigma_x["sg"] * Sigma_y["sg"];
//    dsigma_dpx_T_4[1] = TMath::Pi() * (mu_max_x["sg"] + mu_max_y["sg"]) * Sigma_y["sg"];
//    dsigma_dpx_T_4[2] = 0.;
//    dsigma_dpx_T_4[3] = -1. * TMath::Pi() * Sigma_x["sg"] * Sigma_y["sg"];
//    dsigma_dpy_T_4[0] = TMath::Pi() * Sigma_x["sg"] * Sigma_y["sg"];
//    dsigma_dpy_T_4[1] = TMath::Pi() * (mu_max_x["sg"] + mu_max_y["sg"]) * Sigma_x["sg"];
//    dsigma_dpy_T_4[2] = 0.;
//    dsigma_dpy_T_4[3] = -1. * TMath::Pi() * Sigma_x["sg"] * Sigma_y["sg"];//

//    v_dsigma_dpx["sg"] = new TMatrixD(4, 1);
//    v_dsigma_dpx_T["sg"] = new TMatrixD(1, 4);
//    v_dsigma_dpx["sg"]->Use(4, 1, dsigma_dpx_4);
//    v_dsigma_dpx_T["sg"]->Use(1, 4, dsigma_dpx_T_4);
//    v_dsigma_dpy["sg"] = new TMatrixD(4, 1);
//    v_dsigma_dpy_T["sg"] = new TMatrixD(1, 4);
//    v_dsigma_dpy["sg"]->Use(4, 1, dsigma_dpy_4);
//    v_dsigma_dpy_T["sg"]->Use(1, 4, dsigma_dpy_T_4);//

//    *v_dsigma_dpx["sg"] = *fit_cov_x["sg"] * *v_dsigma_dpx["sg"];
//    *v_dsigma_dpy["sg"] = *fit_cov_y["sg"] * *v_dsigma_dpy["sg"];//

//    sigma_vis_err["sg"] = TMath::Power(10, 19) * TMath::Sqrt(
//                            ((*v_dsigma_dpx_T["sg"]) * (*v_dsigma_dpx["sg"]))(0, 0) +
//                            ((*v_dsigma_dpy_T["sg"]) * (*v_dsigma_dpy["sg"]))(0, 0)
//                          );
//  } else {
//    sigma_vis_err["sg"] = 0.;
//  }//

//  // -- DG
//  Double_t dsigma_dpx_6[] = {
//    TMath::Pi() * Sigma_x["dg"] * Sigma_y["dg"],
//    TMath::Pi() * (mu_max_x["dg"] + mu_max_y["dg"]) * Sigma_y["dg"] * TMath::Power(Sigma_x["dg"], 2) * fr_x["dg"]->Parameter(5) / TMath::Power(fr_x["dg"]->Parameter(1), 2),
//    0.,
//    TMath::Pi() * Sigma_x["dg"] * Sigma_y["dg"],
//    TMath::Pi() * (mu_max_x["dg"] + mu_max_y["dg"]) * Sigma_y["dg"] * TMath::Power(Sigma_x["dg"], 2) * (1. - fr_x["dg"]->Parameter(5)) / TMath::Power(fr_x["dg"]->Parameter(1), 4),
//    TMath::Pi() * (mu_max_x["dg"] + mu_max_y["dg"]) * Sigma_y["dg"] * TMath::Power(Sigma_x["dg"], 2) * (1./fr_x["dg"]->Parameter(1) - 1./fr_x["dg"]->Parameter(4)),
//  };
//  Double_t dsigma_dpy_6[] = {
//    TMath::Pi() * Sigma_x["dg"] * Sigma_y["dg"],
//    TMath::Pi() * (mu_max_x["dg"] + mu_max_y["dg"]) * Sigma_x["dg"] * TMath::Power(Sigma_y["dg"], 2) * fr_y["dg"]->Parameter(5) / TMath::Power(fr_y["dg"]->Parameter(1), 2),
//    0.,
//    TMath::Pi() * Sigma_x["dg"] * Sigma_y["dg"],
//    TMath::Pi() * (mu_max_x["dg"] + mu_max_y["dg"]) * Sigma_x["dg"] * TMath::Power(Sigma_y["dg"], 2) * (1. - fr_y["dg"]->Parameter(5)) / TMath::Power(fr_y["dg"]->Parameter(1), 4),
//    TMath::Pi() * (mu_max_x["dg"] + mu_max_y["dg"]) * Sigma_x["dg"] * TMath::Power(Sigma_y["dg"], 2) * (1./fr_y["dg"]->Parameter(1) - 1./fr_y["dg"]->Parameter(4)),
//  };//

//  Double_t dsigma_dpx_T_6[] = {
//    TMath::Pi() * Sigma_x["dg"] * Sigma_y["dg"],
//    TMath::Pi() * (mu_max_x["dg"] + mu_max_y["dg"]) * Sigma_y["dg"] * TMath::Power(Sigma_x["dg"], 2) * fr_x["dg"]->Parameter(5) / TMath::Power(fr_x["dg"]->Parameter(1), 2),
//    0.,
//    TMath::Pi() * Sigma_x["dg"] * Sigma_y["dg"],
//    TMath::Pi() * (mu_max_x["dg"] + mu_max_y["dg"]) * Sigma_y["dg"] * TMath::Power(Sigma_x["dg"], 2) * (1. - fr_x["dg"]->Parameter(5)) / TMath::Power(fr_x["dg"]->Parameter(1), 4),
//    TMath::Pi() * (mu_max_x["dg"] + mu_max_y["dg"]) * Sigma_y["dg"] * TMath::Power(Sigma_x["dg"], 2) * (1./fr_x["dg"]->Parameter(1) - 1./fr_x["dg"]->Parameter(4)),
//  };
//  Double_t dsigma_dpy_T_6[] = {
//    TMath::Pi() * Sigma_x["dg"] * Sigma_y["dg"],
//    TMath::Pi() * (mu_max_x["dg"] + mu_max_y["dg"]) * Sigma_x["dg"] * TMath::Power(Sigma_y["dg"], 2) * fr_y["dg"]->Parameter(5) / TMath::Power(fr_y["dg"]->Parameter(1), 2),
//    0.,
//    TMath::Pi() * Sigma_x["dg"] * Sigma_y["dg"],
//    TMath::Pi() * (mu_max_x["dg"] + mu_max_y["dg"]) * Sigma_x["dg"] * TMath::Power(Sigma_y["dg"], 2) * (1. - fr_y["dg"]->Parameter(5)) / TMath::Power(fr_y["dg"]->Parameter(1), 4),
//    TMath::Pi() * (mu_max_x["dg"] + mu_max_y["dg"]) * Sigma_x["dg"] * TMath::Power(Sigma_y["dg"], 2) * (1./fr_y["dg"]->Parameter(1) - 1./fr_y["dg"]->Parameter(4)),
//  };//
//

//  if (fit_status["dg"]) {
//    fit_cov_x_dg = new TMatrixDSym(6);
//    fit_cov_y_dg = new TMatrixDSym(6);//

//    *fit_cov_x_dg = fr_x["dg"]->GetCovarianceMatrix();
//    *fit_cov_y_dg = fr_y["dg"]->GetCovarianceMatrix();//

//    TMatrixD *v_dsigma_dpx_dg = new TMatrixD(6, 1);
//    TMatrixD *v_dsigma_dpx_T_dg = new TMatrixD(1, 6);
//    v_dsigma_dpx_dg->Use(6, 1, dsigma_dpx_6);
//    v_dsigma_dpx_T_dg->Use(1, 6, dsigma_dpx_T_6);
//    TMatrixD *v_dsigma_dpy_dg = new TMatrixD(6, 1);
//    TMatrixD *v_dsigma_dpy_T_dg = new TMatrixD(1, 6);
//    v_dsigma_dpy_dg->Use(6, 1, dsigma_dpy_6);
//    v_dsigma_dpy_T_dg->Use(1, 6, dsigma_dpy_T_6);//

//    *v_dsigma_dpx_dg = (*fit_cov_x_dg) * (*v_dsigma_dpx_dg);
//    *v_dsigma_dpy_dg = (*fit_cov_y_dg) * (*v_dsigma_dpy_dg);//

//    sigma_vis_err["dg"] = TMath::Power(10, 19) * TMath::Sqrt(
//                            ((*v_dsigma_dpx_T_dg) * (*v_dsigma_dpx_dg))(0, 0) +
//                            ((*v_dsigma_dpy_T_dg) * (*v_dsigma_dpy_dg))(0, 0)
//                          );
//    cout << "sigma_vis_err[dg] = " << sigma_vis_err["dg"] << endl;
//    
//    delete v_dsigma_dpx_dg;
//    delete v_dsigma_dpx_T_dg;
//    delete v_dsigma_dpy_dg;
//    delete v_dsigma_dpy_T_dg;
//  } else {
//    sigma_vis_err["dg"] = 0.;
//  }

  // -- DGC
  /*Double_t dsigma_dpx_6_dgc[] = {
    TMath::Pi() * Sigma_x["dgc"] * Sigma_y["dgc"],
    TMath::Pi() * (mu_max_x["dgc"] + mu_max_y["dgc"]) * Sigma_y["dgc"] * TMath::Power(Sigma_x["dgc"], 2) * fr_x["dgc"]->Parameter(5) / TMath::Power(fr_x["dgc"]->Parameter(1), 2),
    0.,
    TMath::Pi() * Sigma_x["dgc"] * Sigma_y["dgc"],
    TMath::Pi() * (mu_max_x["dgc"] + mu_max_y["dgc"]) * Sigma_y["dgc"] * TMath::Power(Sigma_x["dgc"], 2) * (1. - fr_x["dgc"]->Parameter(5)) / TMath::Power(fr_x["dgc"]->Parameter(1), 4),
    TMath::Pi() * (mu_max_x["dgc"] + mu_max_y["dgc"]) * Sigma_y["dgc"] * TMath::Power(Sigma_x["dgc"], 2) * (1./fr_x["dgc"]->Parameter(1) - 1./fr_x["dgc"]->Parameter(4)),
  };
  Double_t dsigma_dpy_6_dgc[] = {
    TMath::Pi() * Sigma_x["dgc"] * Sigma_y["dgc"],
    TMath::Pi() * (mu_max_x["dgc"] + mu_max_y["dgc"]) * Sigma_x["dgc"] * TMath::Power(Sigma_y["dgc"], 2) * fr_y["dgc"]->Parameter(5) / TMath::Power(fr_y["dgc"]->Parameter(1), 2),
    0.,
    TMath::Pi() * Sigma_x["dgc"] * Sigma_y["dgc"],
    TMath::Pi() * (mu_max_x["dgc"] + mu_max_y["dgc"]) * Sigma_x["dgc"] * TMath::Power(Sigma_y["dgc"], 2) * (1. - fr_y["dgc"]->Parameter(5)) / TMath::Power(fr_y["dgc"]->Parameter(1), 4),
    TMath::Pi() * (mu_max_x["dgc"] + mu_max_y["dgc"]) * Sigma_x["dgc"] * TMath::Power(Sigma_y["dgc"], 2) * (1./fr_y["dgc"]->Parameter(1) - 1./fr_y["dgc"]->Parameter(4)),
  };

  Double_t dsigma_dpx_T_6_dgc[] = {
    TMath::Pi() * Sigma_x["dgc"] * Sigma_y["dgc"],
    TMath::Pi() * (mu_max_x["dgc"] + mu_max_y["dgc"]) * Sigma_y["dgc"] * TMath::Power(Sigma_x["dgc"], 2) * fr_x["dgc"]->Parameter(5) / TMath::Power(fr_x["dgc"]->Parameter(1), 2),
    0.,
    TMath::Pi() * Sigma_x["dgc"] * Sigma_y["dgc"],
    TMath::Pi() * (mu_max_x["dgc"] + mu_max_y["dgc"]) * Sigma_y["dgc"] * TMath::Power(Sigma_x["dgc"], 2) * (1. - fr_x["dgc"]->Parameter(5)) / TMath::Power(fr_x["dgc"]->Parameter(1), 4),
    TMath::Pi() * (mu_max_x["dgc"] + mu_max_y["dgc"]) * Sigma_y["dgc"] * TMath::Power(Sigma_x["dgc"], 2) * (1./fr_x["dgc"]->Parameter(1) - 1./fr_x["dgc"]->Parameter(4)),
  };
  Double_t dsigma_dpy_T_6_dgc[] = {
    TMath::Pi() * Sigma_x["dgc"] * Sigma_y["dgc"],
    TMath::Pi() * (mu_max_x["dgc"] + mu_max_y["dgc"]) * Sigma_x["dgc"] * TMath::Power(Sigma_y["dgc"], 2) * fr_y["dgc"]->Parameter(5) / TMath::Power(fr_y["dgc"]->Parameter(1), 2),
    0.,
    TMath::Pi() * Sigma_x["dgc"] * Sigma_y["dgc"],
    TMath::Pi() * (mu_max_x["dgc"] + mu_max_y["dgc"]) * Sigma_x["dgc"] * TMath::Power(Sigma_y["dgc"], 2) * (1. - fr_y["dgc"]->Parameter(5)) / TMath::Power(fr_y["dgc"]->Parameter(1), 4),
    TMath::Pi() * (mu_max_x["dgc"] + mu_max_y["dgc"]) * Sigma_x["dgc"] * TMath::Power(Sigma_y["dgc"], 2) * (1./fr_y["dgc"]->Parameter(1) - 1./fr_y["dgc"]->Parameter(4)),
  };


  if (fit_status["dgc"]) {
    fit_cov_x_dgc = new TMatrixDSym(6);
    fit_cov_y_dgc = new TMatrixDSym(6);

    *fit_cov_x_dgc = fr_x["dgc"]->GetCovarianceMatrix();
    *fit_cov_y_dgc = fr_y["dgc"]->GetCovarianceMatrix();

    TMatrixD *v_dsigma_dpx_dgc = new TMatrixD(6, 1);
    TMatrixD *v_dsigma_dpx_T_dgc = new TMatrixD(1, 6);
    v_dsigma_dpx_dgc->Use(6, 1, dsigma_dpx_6_dgc);
    v_dsigma_dpx_T_dgc->Use(1, 6, dsigma_dpx_T_6_dgc);
    TMatrixD *v_dsigma_dpy_dgc = new TMatrixD(6, 1);
    TMatrixD *v_dsigma_dpy_T_dgc = new TMatrixD(1, 6);
    v_dsigma_dpy_dgc->Use(6, 1, dsigma_dpy_6_dgc);
    v_dsigma_dpy_T_dgc->Use(1, 6, dsigma_dpy_T_6_dgc);

    *v_dsigma_dpx_dgc = (*fit_cov_x_dgc) * (*v_dsigma_dpx_dgc);
    *v_dsigma_dpy_dgc = (*fit_cov_y_dgc) * (*v_dsigma_dpy_dgc);

    sigma_vis_err["dgc"] = TMath::Power(10, 19) * TMath::Sqrt(
                             ((*v_dsigma_dpx_T_dgc) * (*v_dsigma_dpx_dgc))(0, 0) +
                             ((*v_dsigma_dpy_T_dgc) * (*v_dsigma_dpy_dgc))(0, 0)
                           );
    delete v_dsigma_dpx_dgc;
    delete v_dsigma_dpx_T_dgc;
    delete v_dsigma_dpy_dgc;
    delete v_dsigma_dpy_T_dgc;
  } else {
    sigma_vis_err["dgc"] = 0.;
  }*/

// -- DGC
  Double_t dsigma_dpx_6_dgc[] = {
    TMath::Pi() * Sigma_x["dgc"] * Sigma_y["dgc"] * TMath::Power(10, 19) * r_y["dgc"],
    0.0,
    0.0,
    0.0,
    TMath::Pi() * (mu_max_x["dgc"]*r_y["dgc"] + mu_max_y["dgc"]*r_x["dgc"]) * Sigma_y["dgc"] * TMath::Power(10, 19),
    0.0,
  };
  Double_t dsigma_dpy_6_dgc[] = {
    TMath::Pi() * Sigma_x["dgc"] * Sigma_y["dgc"] * TMath::Power(10, 19) * r_x["dgc"],
    0.0,
    0.0,
    0.0,
    TMath::Pi() * (mu_max_x["dgc"]*r_y["dgc"] + mu_max_y["dgc"]*r_x["dgc"]) * Sigma_x["dgc"] * TMath::Power(10, 19),
    0.0,
  };

  Double_t dsigma_dpx_T_6_dgc[] = {
    TMath::Pi() * Sigma_x["dgc"] * Sigma_y["dgc"] * TMath::Power(10, 19) * r_y["dgc"],
    0.0,
    0.0,
    0.0,
    TMath::Pi() * (mu_max_x["dgc"]*r_y["dgc"] + mu_max_y["dgc"]*r_x["dgc"]) * Sigma_y["dgc"] * TMath::Power(10, 19),
    0.0,
  };
  Double_t dsigma_dpy_T_6_dgc[] = {
    TMath::Pi() * Sigma_x["dgc"] * Sigma_y["dgc"] * TMath::Power(10, 19) * r_x["dgc"],
    0.0,
    0.0,
    0.0,
    TMath::Pi() * (mu_max_x["dgc"]*r_y["dgc"] + mu_max_y["dgc"]*r_x["dgc"]) * Sigma_x["dgc"] * TMath::Power(10, 19),
    0.0,
  };


  if (fit_status["dgc"]) {
    fit_cov_x_dgc = new TMatrixDSym(6);
    fit_cov_y_dgc = new TMatrixDSym(6);

    *fit_cov_x_dgc = fr_x["dgc"]->GetCovarianceMatrix();
    *fit_cov_y_dgc = fr_y["dgc"]->GetCovarianceMatrix();

    TMatrixD *v_dsigma_dpx_dgc = new TMatrixD(6, 1);
    TMatrixD *v_dsigma_dpx_T_dgc = new TMatrixD(1, 6);
    v_dsigma_dpx_dgc->Use(6, 1, dsigma_dpx_6_dgc);
    v_dsigma_dpx_T_dgc->Use(1, 6, dsigma_dpx_T_6_dgc);
    TMatrixD *v_dsigma_dpy_dgc = new TMatrixD(6, 1);
    TMatrixD *v_dsigma_dpy_T_dgc = new TMatrixD(1, 6);
    v_dsigma_dpy_dgc->Use(6, 1, dsigma_dpy_6_dgc);
    v_dsigma_dpy_T_dgc->Use(1, 6, dsigma_dpy_T_6_dgc);

    *v_dsigma_dpx_dgc = (*fit_cov_x_dgc) * (*v_dsigma_dpx_dgc);
    *v_dsigma_dpy_dgc = (*fit_cov_y_dgc) * (*v_dsigma_dpy_dgc);

    sigma_vis_err["dgc"] = TMath::Sqrt(
                             ((*v_dsigma_dpx_T_dgc) * (*v_dsigma_dpx_dgc))(0, 0) +
                             ((*v_dsigma_dpy_T_dgc) * (*v_dsigma_dpy_dgc))(0, 0)
                           );
    
    delete v_dsigma_dpx_dgc;
    delete v_dsigma_dpx_T_dgc;
    delete v_dsigma_dpy_dgc;
    delete v_dsigma_dpy_T_dgc;
  } else {
    sigma_vis_err["dgc"] = 0.;
  }
	cout << "double gaussian repar, sigma_vis_err = " << sigma_vis_err["dgc"] << endl;
  //// -- DGCL
//  Double_t dsigma_dpx_6_dgcl[] = {
//    TMath::Pi() * Sigma_x["dgcl"] * Sigma_y["dgcl"],
//    TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_y["dgcl"] * TMath::Power(Sigma_x["dgcl"], 2) * fr_x["dgcl"]->Parameter(5) / TMath::Power(fr_x["dgcl"]->Parameter(1), 2),
//    0.,
//    TMath::Pi() * Sigma_x["dgcl"] * Sigma_y["dgcl"],
//    //TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_y["dgcl"] * TMath::Power(Sigma_x["dgcl"], 2) * (1. - fr_x["dgcl"]->Parameter(5)) / TMath::Power(fr_x["dgcl"]->Parameter(1), 4),
//    TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_y["dgcl"] * TMath::Power(Sigma_x["dgcl"], 2) * (1. - fr_x["dgcl"]->Parameter(5)) / TMath::Power(fr_x["dgcl"]->Parameter(4), 2),
//    //TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_y["dgcl"] * TMath::Power(Sigma_x["dgcl"], 2) * (1./fr_x["dgcl"]->Parameter(1) - 1./fr_x["dgcl"]->Parameter(4)),
//    TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_y["dgcl"] * TMath::Power(Sigma_x["dgcl"], 2) * (1./fr_x["dgcl"]->Parameter(4) - 1./fr_x["dgcl"]->Parameter(1)),
//  };
//  Double_t dsigma_dpy_6_dgcl[] = {
//    TMath::Pi() * Sigma_x["dgcl"] * Sigma_y["dgcl"],
//    TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_x["dgcl"] * TMath::Power(Sigma_y["dgcl"], 2) * fr_y["dgcl"]->Parameter(5) / TMath::Power(fr_y["dgcl"]->Parameter(1), 2),
//    0.,
//    TMath::Pi() * Sigma_x["dgcl"] * Sigma_y["dgcl"],
//    //TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_x["dgcl"] * TMath::Power(Sigma_y["dgcl"], 2) * (1. - fr_y["dgcl"]->Parameter(5)) / TMath::Power(fr_y["dgcl"]->Parameter(1), 4),
//    TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_x["dgcl"] * TMath::Power(Sigma_y["dgcl"], 2) * (1. - fr_y["dgcl"]->Parameter(5)) / TMath::Power(fr_y["dgcl"]->Parameter(4), 2),
//    //TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_x["dgcl"] * TMath::Power(Sigma_y["dgcl"], 2) * (1./fr_y["dgcl"]->Parameter(1) - 1./fr_y["dgcl"]->Parameter(4)),
//    TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_x["dgcl"] * TMath::Power(Sigma_y["dgcl"], 2) * (1./fr_y["dgcl"]->Parameter(4) - 1./fr_y["dgcl"]->Parameter(1)),
//  };//

//  Double_t dsigma_dpx_T_6_dgcl[] = {
//    TMath::Pi() * Sigma_x["dgcl"] * Sigma_y["dgcl"],
//    TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_y["dgcl"] * TMath::Power(Sigma_x["dgcl"], 2) * fr_x["dgcl"]->Parameter(5) / TMath::Power(fr_x["dgcl"]->Parameter(1), 2),
//    0.,
//    TMath::Pi() * Sigma_x["dgcl"] * Sigma_y["dgcl"],
//    //TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_y["dgcl"] * TMath::Power(Sigma_x["dgcl"], 2) * (1. - fr_x["dgcl"]->Parameter(5)) / TMath::Power(fr_x["dgcl"]->Parameter(1), 4),
//    TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_y["dgcl"] * TMath::Power(Sigma_x["dgcl"], 2) * (1. - fr_x["dgcl"]->Parameter(5)) / TMath::Power(fr_x["dgcl"]->Parameter(4), 2),
//    //TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_y["dgcl"] * TMath::Power(Sigma_x["dgcl"], 2) * (1./fr_x["dgcl"]->Parameter(1) - 1./fr_x["dgcl"]->Parameter(4)),
//    TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_y["dgcl"] * TMath::Power(Sigma_x["dgcl"], 2) * (1./fr_x["dgcl"]->Parameter(4) - 1./fr_x["dgcl"]->Parameter(1)),
//  };
//  Double_t dsigma_dpy_T_6_dgcl[] = {
//    TMath::Pi() * Sigma_x["dgcl"] * Sigma_y["dgcl"],
//    TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_x["dgcl"] * TMath::Power(Sigma_y["dgcl"], 2) * fr_y["dgcl"]->Parameter(5) / TMath::Power(fr_y["dgcl"]->Parameter(1), 2),
//    0.,
//    TMath::Pi() * Sigma_x["dgcl"] * Sigma_y["dgcl"],
//    //TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_x["dgcl"] * TMath::Power(Sigma_y["dgcl"], 2) * (1. - fr_y["dgcl"]->Parameter(5)) / TMath::Power(fr_y["dgcl"]->Parameter(1), 4),
//    TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_x["dgcl"] * TMath::Power(Sigma_y["dgcl"], 2) * (1. - fr_y["dgcl"]->Parameter(5)) / TMath::Power(fr_y["dgcl"]->Parameter(4), 2),
//    //TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_x["dgcl"] * TMath::Power(Sigma_y["dgcl"], 2) * (1./fr_y["dgcl"]->Parameter(1) - 1./fr_y["dgcl"]->Parameter(4)),
//    TMath::Pi() * (mu_max_x["dgcl"] + mu_max_y["dgcl"]) * Sigma_x["dgcl"] * TMath::Power(Sigma_y["dgcl"], 2) * (1./fr_y["dgcl"]->Parameter(4) - 1./fr_y["dgcl"]->Parameter(1)),
//  };//
//

//  if (fit_status["dgcl"]) {
//    fit_cov_x_dgcl = new TMatrixDSym(6);
//    fit_cov_y_dgcl = new TMatrixDSym(6);//

//    *fit_cov_x_dgcl = fr_x["dgcl"]->GetCovarianceMatrix();
//    *fit_cov_y_dgcl = fr_y["dgcl"]->GetCovarianceMatrix();//

//    TMatrixD *v_dsigma_dpx_dgcl = new TMatrixD(6, 1);
//    TMatrixD *v_dsigma_dpx_T_dgcl = new TMatrixD(1, 6);
//    v_dsigma_dpx_dgcl->Use(6, 1, dsigma_dpx_6_dgcl);
//    v_dsigma_dpx_T_dgcl->Use(1, 6, dsigma_dpx_T_6_dgcl);
//    TMatrixD *v_dsigma_dpy_dgcl = new TMatrixD(6, 1);
//    TMatrixD *v_dsigma_dpy_T_dgcl = new TMatrixD(1, 6);
//    v_dsigma_dpy_dgcl->Use(6, 1, dsigma_dpy_6_dgcl);
//    v_dsigma_dpy_T_dgcl->Use(1, 6, dsigma_dpy_T_6_dgcl);//

//    *v_dsigma_dpx_dgcl = (*fit_cov_x_dgcl) * (*v_dsigma_dpx_dgcl);
//    *v_dsigma_dpy_dgcl = (*fit_cov_y_dgcl) * (*v_dsigma_dpy_dgcl);//

//    sigma_vis_err["dgcl"] = TMath::Power(10, 19) * TMath::Sqrt(
//                              ((*v_dsigma_dpx_T_dgcl) * (*v_dsigma_dpx_dgcl))(0, 0) +
//                              ((*v_dsigma_dpy_T_dgcl) * (*v_dsigma_dpy_dgcl))(0, 0)
//                            );
//                            
//    cout << "sigma_vis_err[dgcl] = " << sigma_vis_err["dgcl"] << endl;
//    delete v_dsigma_dpx_dgcl;
//    delete v_dsigma_dpx_T_dgcl;
//    delete v_dsigma_dpy_dgcl;
//    delete v_dsigma_dpy_T_dgcl;
//  } else {
//    sigma_vis_err["dgcl"] = 0.;
//  }//
//

//  sigma_vis_err["spline"] = 0.;

  // -- L_sp errors
  //newerlumi_sp_err["sgcl"] = lumi_sp["sgcl"] * TMath::Sqrt(
                         // TMath::Power(fr_x["sgcl"]->ParError(1) / fr_x["sgcl"]->Parameter(1), 2) +
                          //TMath::Power(fr_y["sgcl"]->ParError(1) / fr_x["sgcl"]->Parameter(1), 2));

  //lumi_sp_err["sgc"] = lumi_sp["sgc"] * TMath::Sqrt(
                         //TMath::Power(fr_x["sgc"]->ParError(1) / fr_x["sgc"]->Parameter(1), 2) +
                         //TMath::Power(fr_y["sgc"]->ParError(1) / fr_x["sgc"]->Parameter(1), 2));
                         
  lumi_sp_err["sgc"] = TMath::Sqrt(TMath::Power( (lumi_sp["sgc"]*fr_x["sgc"]->ParError(1)) / fr_x["sgc"]->Parameter(1), 2) + 
                                   TMath::Power( (lumi_sp["sgc"]*fr_y["sgc"]->ParError(1)) / fr_y["sgc"]->Parameter(1), 2) );

  //lumi_sp_err["sg"] = lumi_sp["sg"] * TMath::Sqrt(
                        //TMath::Power(fr_x["sg"]->ParError(1) / fr_x["sg"]->Parameter(1), 2) +
                        //TMath::Power(fr_y["sg"]->ParError(1) / fr_x["sg"]->Parameter(1), 2));

  //lumi_sp_err["dg"] = lumi_sp["dg"] * TMath::Sqrt(
                        //TMath::Power(fr_x["dg"]->ParError(1) * fr_x["dg"]->Parameter(5) / TMath::Power(fr_x["dg"]->Parameter(1), 2), 2) +
                        //TMath::Power(fr_x["dg"]->ParError(4) * (1. - fr_x["dg"]->Parameter(5)) / TMath::Power(fr_x["dg"]->Parameter(4), 2), 2) +
                        //TMath::Power(fr_x["dg"]->ParError(5) * (1./fr_x["dg"]->Parameter(1) - 1./fr_x["dg"]->Parameter(4)), 2) +
                        //TMath::Power(fr_y["dg"]->ParError(1) * fr_y["dg"]->Parameter(5) / TMath::Power(fr_y["dg"]->Parameter(1), 2), 2) +
                        //TMath::Power(fr_y["dg"]->ParError(4) * (1. - fr_y["dg"]->Parameter(5)) / TMath::Power(fr_y["dg"]->Parameter(4), 2), 2) +
                        //TMath::Power(fr_y["dg"]->ParError(5) * (1./fr_y["dg"]->Parameter(1) - 1./fr_y["dg"]->Parameter(4)), 2));
                        
  lumi_sp_err["dgc"] = TMath::Sqrt(TMath::Power( (lumi_sp["dgc"]*fr_x["dgc"]->ParError(4)) / fr_x["dgc"]->Parameter(4), 2) + 
                                   TMath::Power( (lumi_sp["dgc"]*fr_y["dgc"]->ParError(4)) / fr_y["dgc"]->Parameter(4), 2) );

  //lumi_sp_err["dgc"] = lumi_sp["dgc"] * TMath::Sqrt(
                         //TMath::Power(fr_x["dgc"]->ParError(1) * fr_x["dgc"]->Parameter(5) / TMath::Power(fr_x["dgc"]->Parameter(1), 2), 2) +
                         //TMath::Power(fr_x["dgc"]->ParError(4) * (1. - fr_x["dgc"]->Parameter(5)) / TMath::Power(fr_x["dgc"]->Parameter(4), 2), 2) +
                         //TMath::Power(fr_x["dgc"]->ParError(5) * (1./fr_x["dgc"]->Parameter(1) - 1./fr_x["dgc"]->Parameter(4)), 2) +
                         //TMath::Power(fr_y["dgc"]->ParError(1) * fr_y["dgc"]->Parameter(5) / TMath::Power(fr_y["dgc"]->Parameter(1), 2), 2) +
                         //TMath::Power(fr_y["dgc"]->ParError(4) * (1. - fr_y["dgc"]->Parameter(5)) / TMath::Power(fr_y["dgc"]->Parameter(4), 2), 2) +
                         //TMath::Power(fr_y["dgc"]->ParError(5) * (1./fr_y["dgc"]->Parameter(1) - 1./fr_y["dgc"]->Parameter(4)), 2));

  //lumi_sp_err["dgcl"] = lumi_sp["dgcl"] * TMath::Sqrt(
                          //TMath::Power(fr_x["dgcl"]->ParError(1) * fr_x["dgcl"]->Parameter(5) / TMath::Power(fr_x["dgcl"]->Parameter(1), 2), 2) +
                          //TMath::Power(fr_x["dgcl"]->ParError(4) * (1. - fr_x["dgcl"]->Parameter(5)) / TMath::Power(fr_x["dgcl"]->Parameter(4), 2), 2) +
                          //TMath::Power(fr_x["dgcl"]->ParError(5) * (1./fr_x["dgcl"]->Parameter(1) - 1./fr_x["dgcl"]->Parameter(4)), 2) +
                          //TMath::Power(fr_y["dgcl"]->ParError(1) * fr_y["dgcl"]->Parameter(5) / TMath::Power(fr_y["dgcl"]->Parameter(1), 2), 2) +
                          //TMath::Power(fr_y["dgcl"]->ParError(4) * (1. - fr_y["dgcl"]->Parameter(5)) / TMath::Power(fr_y["dgcl"]->Parameter(4), 2), 2) +
                          //TMath::Power(fr_y["dgcl"]->ParError(5) * (1./fr_y["dgcl"]->Parameter(1) - 1./fr_y["dgcl"]->Parameter(4)), 2));

  //lumi_sp_err["spline"] = 0.;

  for (vector<TString>::iterator fit_function = fit_functions.begin(); fit_function != fit_functions.end(); ++fit_function) {
    lumi_err[*fit_function] = lumi_sp_err[*fit_function] * (bunch_intensities_1[plb_list_x[12]] + bunch_intensities_1[plb_list_y[12]]) * (bunch_intensities_2[plb_list_x[12]] + bunch_intensities_2[plb_list_y[12]]) / 4.;
  }
}

Double_t VanDerMeerAnalysis::GetFitParameter(TString p_fit_name, TString axis, TString p_fit_parameter) {

  Double_t val;

  if (axis == "x") {
    Int_t idx = fr_x[p_fit_name]->Index(p_fit_parameter.Data());
    if (idx != -1) {
      val = fr_x[p_fit_name]->Parameter(idx);
    } else {
      cout << "[VanDerMeerAnalysis] WARNING : requested fit " << p_fit_name << " and fit parameter " << p_fit_parameter << " could not be found. Returning -1000." << endl;
      val = -1000;
    }
  } else if (axis == "y") {
    Int_t idx = fr_y[p_fit_name]->Index(p_fit_parameter.Data());
    if (idx != -1) {
      val = fr_y[p_fit_name]->Parameter(idx);
    } else {
      cout << "[VanDerMeerAnalysis] WARNING : requested fit " << p_fit_name << " and fit parameter " << p_fit_parameter << " could not be found. Returning -1000." << endl;
      val = -1000;
    }
  } else {
    cout << "[VanDerMeerAnalysis] WARNING : must specify axis when requesting fit parameters" << endl;
  }

  return val;

}

Double_t VanDerMeerAnalysis::GetFitParameterError(TString p_fit_name, TString axis, TString p_fit_parameter) {

  Double_t val;

  if (axis == "x") {
    Int_t idx = fr_x[p_fit_name]->Index(p_fit_parameter.Data());
    if (idx != -1) {
      val = fr_x[p_fit_name]->ParError(idx);
    } else {
      cout << "[VanDerMeerAnalysis] WARNING : requested fit " << p_fit_name << " and fit parameter " << p_fit_parameter << " could not be found. Returning -1000." << endl;
      val = -1000;
    }
  } else if (axis == "y") {
    Int_t idx = fr_y[p_fit_name]->Index(p_fit_parameter.Data());
    if (idx != -1) {
      val = fr_y[p_fit_name]->ParError(idx);
    } else {
      cout << "[VanDerMeerAnalysis] WARNING : requested fit " << p_fit_name << " and fit parameter " << p_fit_parameter << " could not be found. Returning -1000." << endl;
      val = -1000;
    }
  } else {
    cout << "[VanDerMeerAnalysis] WARNING : must specify axis when requesting fit parameters" << endl;
  }

  return val;

}

Double_t VanDerMeerAnalysis::GetVdmParameter(TString p_fit_name, TString p_variable_name) {

  Double_t val;

  std::map<TString, int> pars;
  pars["Sigma_x"] = 1;
  pars["Sigma_y"] = 2;
  pars["mu_max_x"] = 3;
  pars["mu_max_y"] = 4;
  pars["c_x"] = 5;
  pars["c_y"] = 6;

  switch(pars[p_variable_name]) {

    case 1:
      val = Sigma_x[p_fit_name];
      break;
    case 2:
      val = Sigma_y[p_fit_name];
      break;
    case 3:
      val = mu_max_x[p_fit_name];
      break;
    case 4:
      val = mu_max_y[p_fit_name];
      break;
    case 5:
      val = c_x[p_fit_name];
      break;
    case 6:
      val = c_y[p_fit_name];
      break;
  }

  return val;

}

Double_t VanDerMeerAnalysis::GetVdmParameterError(TString p_fit_name, TString p_variable_name) {

  Double_t val;

  if (p_fit_name == "spline") {
    return 0.;
  }

  std::map<TString, int> pars;
  pars["Sigma_x"] = 1;
  pars["Sigma_y"] = 2;
  pars["mu_max_x"] = 3;
  pars["mu_max_y"] = 4;
  pars["c_x"] = 5;
  pars["c_y"] = 6;

  switch(pars[p_variable_name]) {

    case 1:
      val = Sigma_x_err[p_fit_name];
      break;
    case 2:
      val = Sigma_y_err[p_fit_name];
      break;
    case 3:
      val = mu_max_x_err[p_fit_name];
      break;
    case 4:
      val = mu_max_y_err[p_fit_name];
      break;
    case 5:
      val = c_x_err[p_fit_name];
      break;
    case 6:
      val = c_y_err[p_fit_name];
      break;
  }

  return val;

}


Double_t VanDerMeerAnalysis::GetChi2Ndf(TString p_fit_name, TString axis) {

  Double_t val;

  if (p_fit_name == "spline") {
    val = 0.;
  } else {
    if (axis == "x") {
      val = fr_x[p_fit_name]->Chi2() / fr_x[p_fit_name]->Ndf();
    } else if (axis == "y") {
      val = fr_y[p_fit_name]->Chi2() / fr_y[p_fit_name]->Ndf();
    }
  }

  return val;

}

Double_t VanDerMeerAnalysis::GetLumiSp(TString p_fit_name) {

  if (lumi_sp.find(p_fit_name) == lumi_sp.end()) {
    cout << "[VanDerMeerAnalysis] WARNING : In GetLumiSp, requested fit method " << p_fit_name << " was not found. Returning -1." << endl;
    return -1;
  }

  return lumi_sp[p_fit_name];

}

Double_t VanDerMeerAnalysis::GetLumiSpErr(TString p_fit_name) {

  if (lumi_sp_err.find(p_fit_name) == lumi_sp_err.end()) {
    cout << "[VanDerMeerAnalysis] WARNING : In GetLumiSpErr, requested fit method " << p_fit_name << " was not found. Returning -1." << endl;
    return -1;
  }

  return lumi_sp_err[p_fit_name];

}

Double_t VanDerMeerAnalysis::GetLumi(TString p_fit_name) {
  if (lumi.find(p_fit_name) == lumi.end()) {
    cout << "[VanDerMeerAnalysis] WARNING : In GetLumi, requested fit method " << p_fit_name << " was not found. Returning -1." << endl;
    return -1;
  }

  return lumi[p_fit_name];

}

Double_t VanDerMeerAnalysis::GetLumiErr(TString p_fit_name) {

  if (lumi_err.find(p_fit_name) == lumi_err.end()) {
    cout << "[VanDerMeerAnalysis] WARNING : In GetLumiErr, requested fit method " << p_fit_name << " was not found. Returning -1." << endl;
    return -1;
  }

  return lumi_err[p_fit_name];

}


Double_t VanDerMeerAnalysis::GetSigmaVis(TString p_fit_name) {

  if (sigma_vis.find(p_fit_name) == sigma_vis.end()) {
    cout << "[VanDerMeerAnalysis] WARNING : In GetSigmaVis, requested fit method " << p_fit_name << " was not found. Returning -1." << endl;
    return -1;
  }

  return sigma_vis[p_fit_name];

}

Double_t VanDerMeerAnalysis::GetSigmaVisErr(TString p_fit_name) {

  if (sigma_vis_err.find(p_fit_name) == sigma_vis_err.end()) {
    cout << "[VanDerMeerAnalysis] WARNING : In GetSigmaVisErr, requested fit method " << p_fit_name << " was not found. Returning -1." << endl;
    return -1;
  }

  cout << "[VanDerMeerAnalysis] INFO : In GetSigmaVisErr, requested fit method " << p_fit_name << " has error " << sigma_vis_err[p_fit_name] << endl;

  return sigma_vis_err[p_fit_name];

}

TGraphErrors* VanDerMeerAnalysis::GetTGraphX() {

  return tg_musp_x;

}

TGraphErrors* VanDerMeerAnalysis::GetTGraphY() {

  return tg_musp_y;

}

Double_t VanDerMeerAnalysis::CalculateTSplineIntegral(TSpline *ts, Double_t xmin, Double_t xmax, Int_t npx) {
  Double_t current_sum=0.;
  for (int i=0; i<=npx; i++) {
    Double_t x = xmin + (xmax-xmin)*i/npx;
    current_sum += ts->Eval(x) * (xmax-xmin)/npx;
  }
  return current_sum;
}

Double_t VanDerMeerAnalysis::CalculateTSplineMax(TSpline *ts, Double_t xmin, Double_t xmax) {
  Double_t current_ymax = 0.;//Value of maximum
  Double_t current_xmax = 0.;//x-coordinate of maximum
  Int_t npx = 10000;
  for (int i=0; i<=npx; i++) {
    Double_t x = xmin + i*(xmax-xmin)/npx;
    Double_t current_value = ts->Eval(x);
    if (current_value > current_ymax) {
      current_ymax = current_value;
      current_xmax = x;
    }
  }
  Double_t xmin2 = current_xmax - (xmax-xmin)/npx;
  Double_t xmax2 = current_xmax + (xmax-xmin)/npx;
  for (int j=0; j<=npx; j++) {
    Double_t x = xmin2 + j*(xmax2-xmin2)/npx;
    Double_t current_value = ts->Eval(x);
    if (current_value > current_ymax) {
      current_ymax = current_value;
      current_xmax = x;
    }
  }
  //Two iterations enough?
  return current_ymax;
}

void VanDerMeerAnalysis::SetSystematicUncertaintyFlag(TString p_systematic_name) {

  systematic_uncertainty_list[p_systematic_name] = true;

}


void VanDerMeerAnalysis::DebugPlots(TString p_path, TString p_tag) {


  /**
    *  Write out some debug info:
    *  - TGraphs of {n_vtx, mu, mu_sp, fake_fraction, pileup_correction_factor} vs. separation with SGCL
    */

  /*
    tg_musp_x = new TGraphErrors(plb_list_x.size());
    tg_musp_x->SetName("tg_musp_x");
    tg_musp_y = new TGraphErrors(plb_list_y.size());
    tg_musp_y->SetName("tg_musp_y");

    Double_t mu_max_guess_x, mu_max_guess_y;

    for (int i=0; i<plb_list_x.size(); i++) {

      tg_musp_x->SetPoint(i, plb_list_x[i], musp_pLB[plb_list_x[i]]);
      tg_musp_x->SetPoint(i, plb_list_x[i], musp_err_pLB[plb_list_x[i]]);
      tg_musp_y->SetPoint(i, plb_list_y[i], musp_pLB[plb_list_y[i]]);
      tg_musp_y->SetPoint(i, plb_list_y[i], musp_err_pLB[plb_list_y[i]]);

      if (i == 12) {
        mu_max_guess_x = musp_pLB[plb_list_x[i]];
        mu_max_guess_y = musp_pLB[plb_list_y[i]];
      }
    }

  */

  cout << "[VanDerMeerAnalysis] INFO : Saving debug histograms to " << p_path << " with save tag " << p_tag << endl;

  TFile *f_debug = new TFile(p_path, "UPDATE");

  std::map<TString, TGraphErrors*> tg_n_vtx, tg_mu, tg_mu_fake, tg_masking_correction_factors, tg_mu_corr_fake, tg_mu_corr_mask, tg_mu_corr_total, tg_mask_corr;

  TString name;

  // x
  name = "tg_n_vtx_x";
  name += p_tag;
  tg_n_vtx["x"] = new TGraphErrors(plb_list_x.size());
  tg_n_vtx["x"]->SetName(name);

  name = "tg_mu_x";
  name += p_tag;
  tg_mu["x"] = new TGraphErrors(plb_list_x.size());
  tg_mu["x"]->SetName(name);

  name = "tg_mu_fake_x";
  name += p_tag;
  tg_mu_fake["x"] = new TGraphErrors(plb_list_x.size());
  tg_mu_fake["x"]->SetName(name);

  name = "tg_masking_correction_factors_x";
  name += p_tag;
  tg_masking_correction_factors["x"] = new TGraphErrors(plb_list_x.size());
  tg_masking_correction_factors["x"]->SetName(name);

  name = "tg_mu_corr_total_x";
  name += p_tag;
  tg_mu_corr_total["x"] = new TGraphErrors(plb_list_x.size());
  tg_mu_corr_total["x"]->SetName(name);

  name = "tg_mu_corr_fake_x";
  name += p_tag;
  tg_mu_corr_fake["x"] = new TGraphErrors(plb_list_x.size());
  tg_mu_corr_fake["x"]->SetName(name);

  name = "tg_mu_corr_mask_x";
  name += p_tag;
  tg_mu_corr_mask["x"] = new TGraphErrors(plb_list_x.size());
  tg_mu_corr_mask["x"]->SetName(name);

  name = "tg_mask_corr_x";
  name += p_tag;
  tg_mask_corr["x"] = new TGraphErrors(plb_list_x.size());
  tg_mask_corr["x"]->SetName(name);

  for (unsigned int i=0; i<plb_list_x.size(); i++) {
    Int_t current_pLB = plb_list_x[i];

    tg_n_vtx["x"]->SetPoint(i, nominal_separation[current_pLB], nvtx_pLB[current_pLB]);
    tg_n_vtx["x"]->SetPointError(i, 0., nvtx_err_pLB[current_pLB]);

    tg_mu["x"]->SetPoint(i, nominal_separation[current_pLB], mu_pLB[current_pLB]);
    tg_mu["x"]->SetPointError(i, 0., mu_err_pLB[current_pLB]);

    tg_mu_fake["x"]->SetPoint(i, nominal_separation[current_pLB], mu_fake_list[current_pLB]);
    tg_mu_fake["x"]->SetPointError(i, 0., mu_fake_uncertainty_list[current_pLB]);

    tg_mask_corr["x"]->SetPoint(i, mu_raw_pLB[current_pLB], masking_correction_factors[current_pLB]);

    tg_masking_correction_factors["x"]->SetPoint(i, nominal_separation[current_pLB], masking_correction_factors[current_pLB]);
    cout << "mu_raw = " << mu_raw_pLB[current_pLB] << ", fake correction (%) = " << (mu_fake_list[current_pLB]/mu_raw_pLB[current_pLB])*100 << endl;
    tg_mu_corr_fake["x"]->SetPoint(i,mu_raw_pLB[current_pLB],(mu_fake_list[current_pLB]/mu_raw_pLB[current_pLB])*100);
    cout << "mu_real = " << mu_real_pLB[current_pLB] << ", masking correction (%) = " << (1-masking_correction_factors[current_pLB])*100 << endl;
    tg_mu_corr_mask["x"]->SetPoint(i,mu_real_pLB[current_pLB],(1-masking_correction_factors[current_pLB])*100);
    cout << "Total correction for this mu_raw is (%) = " << ((mu_raw_pLB[current_pLB]-mu_pLB[current_pLB])/mu_raw_pLB[current_pLB])*100 << endl;
    tg_mu_corr_total["x"]->SetPoint(i, mu_raw_pLB[current_pLB], ((mu_raw_pLB[current_pLB]-mu_pLB[current_pLB])/mu_raw_pLB[current_pLB])*100);
  }

  tg_n_vtx["x"]->Write();
  tg_mu["x"]->Write();
  tg_musp_x->Write();
  tg_mu_fake["x"]->Write();
  tg_masking_correction_factors["x"]->Write();
  tg_mu_corr_fake["x"]->Write();
  tg_mu_corr_mask["x"]->Write();
  tg_mu_corr_total["x"]->Write();
  tg_mask_corr["x"]->Write();

  // y
  name = "tg_n_vtx_y";
  name += p_tag;
  tg_n_vtx["y"] = new TGraphErrors(plb_list_y.size());
  tg_n_vtx["y"]->SetName(name);

  name = "tg_mu_y";
  name += p_tag;
  tg_mu["y"] = new TGraphErrors(plb_list_y.size());
  tg_mu["y"]->SetName(name);

  name = "tg_mu_fake_y";
  name += p_tag;
  tg_mu_fake["y"] = new TGraphErrors(plb_list_y.size());
  tg_mu_fake["y"]->SetName(name);

  name = "tg_masking_correction_factors_y";
  name += p_tag;
  tg_masking_correction_factors["y"] = new TGraphErrors(plb_list_y.size());
  tg_masking_correction_factors["y"]->SetName(name);

  for (unsigned int i=0; i<plb_list_y.size(); i++) {
    Int_t current_pLB = plb_list_y[i];

    tg_n_vtx["y"]->SetPoint(i, nominal_separation[current_pLB], nvtx_pLB[current_pLB]);
    tg_n_vtx["y"]->SetPointError(i, 0., nvtx_err_pLB[current_pLB]);

    tg_mu["y"]->SetPoint(i, nominal_separation[current_pLB], mu_pLB[current_pLB]);
    tg_mu["y"]->SetPointError(i, 0., mu_err_pLB[current_pLB]);

    tg_mu_fake["y"]->SetPoint(i, nominal_separation[current_pLB], mu_fake_list[current_pLB]);
    tg_mu_fake["y"]->SetPointError(i, nominal_separation[current_pLB], mu_fake_uncertainty_list[current_pLB]);

    tg_masking_correction_factors["y"]->SetPoint(i, nominal_separation[current_pLB], masking_correction_factors[current_pLB]);

  }

  tg_n_vtx["y"]->Write();
  tg_mu["y"]->Write();
  tg_musp_y->Write();
  tg_mu_fake["y"]->Write();
  tg_masking_correction_factors["y"]->Write();

  //fc->GetFakeCorrectionTGraph()->Write(TString("tg_fake_correction_") + p_tag);
  //pmc->GetDifferentialPmask()->Write(TString("h_pmask_dz_") + p_tag);
  f_debug->Close();

  for (map<TString, TGraphErrors*>::iterator itG = tg_n_vtx.begin();
       itG != tg_n_vtx.end(); ++itG) {
    if (itG->second) {
      delete itG->second;
    }
  }
  tg_n_vtx.clear();

  for (map<TString, TGraphErrors*>::iterator itG = tg_mu.begin();
       itG != tg_mu.end(); ++itG) {
    if (itG->second) {
      delete itG->second;
    }
  }
  tg_mu.clear();

  for (map<TString, TGraphErrors*>::iterator itG = tg_mu_fake.begin();
       itG != tg_mu_fake.end(); ++itG) {
    if (itG->second) {
      delete itG->second;
    }
  }
  tg_mu_fake.clear();

  for (map<TString, TGraphErrors*>::iterator itG = tg_masking_correction_factors.begin();
       itG != tg_masking_correction_factors.end(); ++itG) {
    if (itG->second) {
      delete itG->second;
    }
  }
  tg_masking_correction_factors.clear();

}

Double_t fit_func_sg(Double_t *x, Double_t *par) {

  Double_t mu_max=par[0];
  Double_t h0 = par[2];
  return mu_max * TMath::Exp(-TMath::Power((x[0]-h0),2)/(2*TMath::Power(par[1],2)));

}
Double_t fit_func_dg(Double_t *x, Double_t *par) {

  /*
  par[0]=mu_max (fixed)
  par[1]=sigma_a
  par[2]=h0 (fixed)
  par[3]=c (fixed)
  par[4]=sigma_b
  par[5]=f_a
  */
  //Fixed parameters
  //Double_t mu_max = par[0];
  //Double_t h0 = par[2];
  //Double_t sigma_b = par[4] * par[1] * (1-par[5]) / ( par[1] - par[5]*par[4] );
  Double_t Sigma_h = 1./(par[5]/par[1] + (1.-par[5])/par[4]);

  return par[3] + par[0] * Sigma_h * (
           par[5]/par[1] * TMath::Exp(- TMath::Power(x[0]-par[2],2)/(2*TMath::Power(par[1],2))) +
           (1.-par[5])/par[4] * TMath::Exp(-TMath::Power((x[0]-par[2]),2)/(2*TMath::Power(par[4],2))) );


}

//Double_t fit_func_dgc(Double_t *x, Double_t *par) {

  /*
  par[0]=mu_max (fixed)
  par[1]=sigma_a
  par[2]=h0 (fixed)
  par[3]=c
  par[4]=sigma_b
  par[5]=f_a
  */
  //Fixed parameters
  //Double_t mu_max = par[0];
  //Double_t h0 = par[2];
  //Double_t sigma_b = par[4] * par[1] * (1-par[5]) / ( par[1] - par[5]*par[4] );
  //Double_t Sigma_h = 1./(par[5]/par[1] + (1.-par[5])/par[4]);

 // return par[3] + par[0] * Sigma_h * (
           //par[5]/par[1] * TMath::Exp(- TMath::Power(x[0]-par[2],2)/(2*TMath::Power(par[1],2))) +
          //(1.-par[5])/par[4] * TMath::Exp(-TMath::Power((x[0]-par[2]),2)/(2*TMath::Power(par[4],2))) );


//}

/*Double_t fit_func_dgc(Double_t *x, Double_t *par) {

  //par[0]=mu_max 
	//par[1]=alpha
	//par[2]=h0 
	//par[3]=c
	//par[4]=Sigma_x
	//par[5]=f_a
							
  //Double_t sigma_a = ((1.-par[5])*par[4])/(par[1]) + (par[5]*par[4]);
  //Double_t sigma_b = par[1]*sigma_a;
  Double_t sigma_a = par[1];
  Double_t sigma_b = ((1-par[5])*par[4]*sigma_a)/(sigma_a-(par[5]*par[4]));
	return par[3] + par[0] * par[4] * (
	       par[5]/sigma_a * TMath::Exp(- TMath::Power(x[0]-par[2],2)/(2*TMath::Power(sigma_a,2))) +
	       (1.-par[5])/sigma_b * TMath::Exp(-TMath::Power((x[0]-par[2]),2)/(2*TMath::Power(sigma_b,2))) );
}*/

Double_t fit_func_dgc(Double_t *x, Double_t *par) {

  // [0] = mu_max
  // [1] = sigma_a
  // [2] = h0
  // [3] = c
  // [4] = Sigma_x
  // [5] = fa

  return par[0] * TMath::Abs(par[4]) * ( par[5]/TMath::Abs(par[1]) * TMath::Exp( -0.5*TMath::Power( ( (x[0]-par[2])/TMath::Abs(par[1]) ), 2 ) ) + (1.0-par[5]) / ((1.0-par[5])*TMath::Abs(par[1])*TMath::Abs(par[4])/(TMath::Abs(par[1])-par[5]*TMath::Abs(par[4]))) * TMath::Exp( -0.5*( TMath::Power( ( (x[0]-par[2])/((1.0-par[5])*TMath::Abs(par[1])*TMath::Abs(par[4])/(TMath::Abs(par[1])-par[5]*TMath::Abs(par[4]))) ), 2) ) ) ) + par[3];
}

Double_t fit_func_dgcl(Double_t *x, Double_t *par) {

  /*
  par[0]=mu_max (fixed)
  par[1]=sigma_a
  par[2]=h0 (fixed)
  par[3]=c
  par[4]=sigma_b
  par[5]=f_a
  */
  //Fixed parameters
  //Double_t mu_max = par[0];
  //Double_t h0 = par[2];
  //Double_t sigma_b = par[4] * par[1] * (1-par[5]) / ( par[1] - par[5]*par[4] );
  Double_t Sigma_h = 1./( (par[5]/par[1]) + ( (1.- par[5])/par[4] ) );

  return par[3] + par[0] * Sigma_h * ( par[5]/par[1] * TMath::Exp(-1.*TMath::Power(x[0]-par[2],2)/(2*TMath::Power(par[1],2))) +
           (1.-par[5])/par[4] * TMath::Exp(-1.*TMath::Power((x[0]-par[2]),2)/(2*TMath::Power(par[4],2))) );


}
Double_t fit_func_sgc(Double_t *x, Double_t *par) {

  /*
  par[0]=mu_max (fixed)
  par[1]=sigma_a
  par[2]=h0 (fixed)
  par[3]=c
  */
  Double_t mu_max=par[0];
  Double_t h0 = par[2];
  return par[3] + mu_max * TMath::Exp(-TMath::Power((x[0]-h0),2)/(2*TMath::Power(par[1],2)));


}


