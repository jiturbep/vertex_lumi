#define DEBUG_LUMIVTX

#include "LumiVtx/LumiVtx.h"

#include "TKey.h"
#include "TCollection.h"
#include "TDirectory.h"

Double_t fitFuncSimpleGaussian(Double_t *x, Double_t *par);

using namespace std;

LumiVtx::LumiVtx() : fc(0) {
  verbose = 0;
  pLB_min = 10000;
  pLB_max = 0;

  n_bcids = 1;

  correct_masking = false;
  correct_fakes = false;

  systematic_uncertainty_list["fake_low"] = false;
  systematic_uncertainty_list["fake_high"] = false;
  systematic_uncertainty_list["masking_toy_scaling"] = false;
  //Changed to true
  systematic_uncertainty_list["bs-45mm"] = false;
  systematic_uncertainty_list["bs-55mm"] = false;

#ifdef USE_UNFOLD
  m_unfold = 0;
#endif


}

LumiVtx::~LumiVtx() {
#ifdef USE_UNFOLD
  if (m_unfold) {
    delete m_unfold;
  }
#endif

}

void LumiVtx::LoadTTree(TString path, TString branch_nvtx, TString branch_nevt, TString branch_ntrig, std::vector<std::pair<TString, Int_t> > *selection, bool do_physics_run) {

  // SP, Note: We'll not attempt to load unfolding results here.
  //  Corrections are not loaded yet, since they need pLB which are loaded here (??).
  // Call specific method for loading unfolding results LoadVertexCountsUnfold(...)

  cout << "[LumiVtx] INFO : Loading TTree... at " << path << std::endl;
  //Note: and also determine the pLB range! So, this should be called as early as possible, because later things will use this info!


  TFile *f_in = new TFile(path, "READ");
  TTree *t_vdm = (TTree*)f_in->Get("t_vdm");

  t_vdm->SetBranchStatus("*", 0);
  t_vdm->SetBranchStatus(branch_nvtx, 1);
  t_vdm->SetBranchStatus(branch_nevt, 1);
  t_vdm->SetBranchStatus(branch_ntrig, 1);
  t_vdm->SetBranchStatus("pLB", 1);
  for (vector<std::pair<TString, Int_t> >::iterator it = selection->begin(); it != selection->end(); ++it) {
    t_vdm->SetBranchStatus((*it).first, 1);
  }
  Float_t current_nvtx;
  Float_t current_nevt;
  Float_t current_ntrig;
  Int_t current_pLB;
  t_vdm->SetBranchAddress(branch_nvtx, &current_nvtx);
  t_vdm->SetBranchAddress(branch_nevt, &current_nevt);
  t_vdm->SetBranchAddress("pLB", &current_pLB);
  t_vdm->SetBranchAddress(branch_ntrig, &current_ntrig);


  Long64_t n_entries = t_vdm->GetEntriesFast();
  for (int i=0; i < n_entries; i++) {
    t_vdm->GetEntry(i);

    //Impose selection
    bool good = true;
    for (vector<std::pair<TString, Int_t> >::iterator it = selection->begin(); it != selection->end(); ++it) {
      Int_t value = t_vdm->GetLeaf((*it).first)->GetValue(0);
      if (value != (*it).second) {
        good = false;
        break;
      }
    }
    if (!good) {
      continue;
    }

    //Fill maps
    n_vtx_raw[current_pLB] = current_nvtx;
    n_evt_raw[current_pLB] = current_nevt;
    n_trig_raw[current_pLB] = current_ntrig;
    n_vtx_raw_err[current_pLB] = TMath::Sqrt(current_nvtx);
    n_evt_raw_err[current_pLB] = 0.; //Can't calculate this quite yet, need deadtime.

#ifdef DEBUG_LUMIVTX
    cout << "[LumiVtx] DEBUG : pLB " << current_pLB << " : NVtx = " << current_nvtx << " / NEvt = " << current_nevt << " / NTrig =" << current_ntrig << endl;
    cout << "[LumiVtx] DEBUG : n_vtx_raw_err[" << current_pLB << "] = " << n_vtx_raw_err[current_pLB] << endl;
#endif

		/*if ( do_physics_run ) {
			if ( current_ntrig != 0 ) {
				if (current_pLB < pLB_min) pLB_min = current_pLB;
				if (current_pLB > pLB_max) pLB_max = current_pLB;
				lumiblocks.push_back(current_pLB);
				lb_duration[current_pLB] = 1.0;
			} 
		}*/  

    //Do this from the timestamp file -- a little more controllable.
    //Max/min pLB
    //if (current_pLB < pLB_min) pLB_min = current_pLB;
    //if (current_pLB > pLB_max) pLB_max = current_pLB;
    //lumiblocks.push_back(current_pLB);

    

  }
  cout << "done" << endl;
}

void LumiVtx::LoadLivefraction(TString path) {
  cout << "Importing deadtime..." << endl;

  #ifdef DEBUG_LUMIVTX
  cout << "Path: " << path << endl;
  #endif

  TFile *f_in = new TFile(path, "READ");
  TH1F *h_livefraction = (TH1F*)f_in->Get("PlbLiveFractions");
  for (int i=1; i<h_livefraction->GetNbinsX(); i++) {
    live_fraction[TMath::Nint(h_livefraction->GetBinCenter(i))] = h_livefraction->GetBinContent(i);
    #ifdef DEBUG_LUMIVTX
    if (i % 50 == 0) {
      cout << "pLB " << h_livefraction->GetBinCenter(i) << " : live fraction " << h_livefraction->GetBinContent(i) << endl;
    }
    #endif
  }
  /*
  for (Int_t current_pLB = pLB_min; current_pLB < pLB_max; current_pLB++) {
    live_fraction[current_pLB] = h_livefraction->GetBinContent(h_livefraction->FindBin(current_pLB));
    #ifdef DEBUG_LUMIVTX
    cout << "pLB " << current_pLB << ", live fraction = " << live_fraction[current_pLB] << endl;
    #endif
  }
  */
}

void LumiVtx::LoadLivefraction() {
  //If no argument is given, set all live fractions to 1.
  cout << "[LumiVtx] Setting all live fraction corrections to 1." << endl;
  for (Int_t current_pLB = pLB_min; current_pLB < pLB_max; current_pLB++) {
    live_fraction[current_pLB] = 1.0;
  }
}

void LumiVtx::LoadPrescale(TString path) {
  cout << "Importing prescales..." << endl;

  #ifdef DEBUG_LUMIVTX
  cout << "Path: " << path << endl;
  #endif

  TFile *f_in = new TFile(path, "READ");
  TH1F *h_prescale = (TH1F*)f_in->Get("PlbPrescale");
  for (int i=1; i<h_prescale->GetNbinsX(); i++) {

    prescale[TMath::Nint(h_prescale->GetBinCenter(i))] = (Float_t)h_prescale->GetBinContent(i);

#ifdef DEBUG_LUMIVTX
    if (i % 50 == 0) {
      cout << "pLB " << h_prescale->GetBinCenter(i) << " : prescale " << h_prescale->GetBinContent(i) << endl;
    }
#endif
  }
}

void LumiVtx::LoadPrescale() {
  //If no argument is given, set all prescales to 1.
  cout << "[LumiVtx] Setting all prescales to 1." << endl;
  for (Int_t current_pLB = pLB_min; current_pLB < pLB_max; current_pLB++) {
    prescale[current_pLB] = 1.;
  }
}

void LumiVtx::SetSigmaVisNVtx(Float_t value) {
  sigmavis[GlobalSettings::kVtxC] = value;
  cout << "[LumiVtx] INFO : Using sigmavis_nvtx = " << sigmavis[GlobalSettings::kVtxC] << endl;
}

void LumiVtx::SetSigmaVisNEvt(Float_t value) {
  sigmavis[GlobalSettings::kEvtC] = value;
  cout << "[LumiVtx] INFO : Using sigmavis_nevt = " << sigmavis[GlobalSettings::kEvtC] << endl;
}

void LumiVtx::SetSigmaVis(GlobalSettings::AnalysisMethods method, Float_t value, Float_t error) {
  switch (method) {
    case GlobalSettings::kVtxC:
      SetSigmaVisNVtx(value);
      cout << "[LumiVtx] Line 205 value = " << value << endl; 
      break;
    case GlobalSettings::kEvtC:
      SetSigmaVisNEvt(value);
      break;
    case GlobalSettings::kUnfC:
      sigmavis[GlobalSettings::kUnfC] = value;
      sigmavis_err[GlobalSettings::kUnfC] = error;
      cout << "[LumiVtx] INFO: Using sigmavis for UnfC = " << sigmavis[GlobalSettings::kUnfC]
           << " +/- " << sigmavis_err[GlobalSettings::kUnfC] << endl;
      break;
    default:
      cerr << "[LumiVtx::SetSigmaVis] Requested invalid method: " << method << endl;
  }
}

void LumiVtx::LoadTimestamps(TString path) {

  cout << "[LumiVtx] INFO : Loading timestamps from " << path << endl;

  ifstream f_in(path.Data());
  if (not f_in.is_open()) {
    cerr << "Unable to open timestamps text file: " << path << endl;
    exit(1);
  }

  Int_t current_pLB;
  Double_t ts_start, ts_end;
  TString everything_else;

  while (f_in.good()) {
    string line;
    getline(f_in, line);
    if (line.empty()) {
      continue;
    }
    stringstream istr(line);
    istr >> current_pLB >> ts_start >> ts_end >> everything_else;

    pLB_timestamps[current_pLB] = std::make_pair<Double_t, Double_t>(ts_start, ts_end);
    lb_duration[current_pLB] = ts_end - ts_start;

    if (current_pLB < pLB_min) {
      pLB_min = current_pLB;
    }
    if (current_pLB > pLB_max) {
      pLB_max = current_pLB;
    }
    lumiblocks.push_back(current_pLB);


#ifdef DEBUG_LUMIVTX
    cout << "Set timestamp " << current_pLB << " : ("<< ts_start << ", " << ts_end << "), lb_duration[" << current_pLB << "] = " << lb_duration[current_pLB] << endl;
#endif
  }
  sort(lumiblocks.begin(), lumiblocks.end());
  n_pLBs = pLB_max - pLB_min + 1;
#ifdef DEBUG_LUMIVTX
  cout << "pLB range : " << pLB_min << " - " << pLB_max << endl;;
#endif

}

void LumiVtx::LoadDurations(TString path) {

  cout << "[LumiVtx] INFO : Loading lumiblock durations from " << path << endl;

  ifstream f_in(path.Data());
  if (not f_in.is_open()) {
    cerr << "Unable to open timestamps text file: " << path << endl;
    exit(1);
  }

  Int_t current_pLB;
  Double_t current_duration;
  TString everything_else;

  while (f_in.good()) {
    string line;
    getline(f_in, line);
    if (line.empty()) {
      continue;
    }
    stringstream istr(line);
    istr >> current_pLB >> current_duration >> everything_else;

    lb_duration[current_pLB] = current_duration;

    if (current_pLB < pLB_min) {
      pLB_min = current_pLB;
    }
    if (current_pLB > pLB_max) {
      pLB_max = current_pLB;
    }
    lumiblocks.push_back(current_pLB);


#ifdef DEBUG_LUMIVTX
    cout << "Set lb_duration[" << current_pLB << "] = " << lb_duration[current_pLB] << endl;
#endif
  }
  sort(lumiblocks.begin(), lumiblocks.end());
  n_pLBs = pLB_max - pLB_min + 1;
#ifdef DEBUG_LUMIVTX
  cout << "pLB range : " << pLB_min << " - " << pLB_max << endl;;
#endif

}

void LumiVtx::SetAllDurations(const unsigned int &LBlow, const unsigned int &LBhigh) {
  pLB_min = LBlow;
  pLB_max = LBhigh;
  for (unsigned int iLB = pLB_min; iLB <= pLB_max; ++iLB) {
    lb_duration[iLB] = 1.0;
    lumiblocks.push_back(iLB);
  }
}

void LumiVtx::LoadDurationsFromLumiNtuple(TString path) {

  cout << "[LumiVtx] INFO : Loading duration from Eric's lumi ntuple at " << path << endl;

  TFile *f_in = TFile::Open(path, "READ");
  TTree *t_lumi = (TTree*)f_in->Get("lumiData");

  t_lumi->SetBranchStatus("*", 0);
  t_lumi->SetBranchStatus("LBDATA", 1);

  Long64_t entries = t_lumi->GetEntriesFast();
  for (int i = 0; i < entries; i++) {
    t_lumi->GetEntry(i);
    Int_t current_lb = t_lumi->GetLeaf("LB")->GetValue(0);
    ULong64_t current_ts_start = t_lumi->GetLeaf("StartTime")->GetValue(0);
    ULong64_t current_ts_end = t_lumi->GetLeaf("EndTime")->GetValue(0);
    lb_duration[current_lb] = (current_ts_end - current_ts_start) / TMath::Power(10, 9);

    if (current_lb < pLB_min) {
      pLB_min = current_lb;
    }
    if (current_lb > pLB_max) {
      pLB_max = current_lb;
    }
    lumiblocks.push_back(current_lb);

  }

  sort(lumiblocks.begin(), lumiblocks.end());
  n_pLBs = pLB_max - pLB_min + 1; 
#ifdef DEBUG_LUMIVTX
  cout << "[LumiVtx] DEBUG : Found " << n_pLBs << " lumiblocks." << endl;
  for (std::map<Int_t, Double_t>::iterator it = lb_duration.begin(); it != lb_duration.end(); ++it) {

    cout << "[LumiVtx] DEBUG : LB " << (*it).first << " : " << (*it).second << endl;
  }
#endif
}

#ifdef USE_UNFOLD
void LumiVtx::LoadVertexCountsUnfold(TString p_path, Int_t p_bcid, Int_t nTrkCut, TString triggerType) {
  //First setup unfolding, if not done yet
  //init unfolding machinery
  if (!m_unfold) {
    cout << "[VanDerMeerAnalysis] Setting up unfolding " << endl;
    m_unfold = new VertexLumiUnfold();
    TString tf_name;
    tf_name = responseMatrixHName;
    tf_name += "_NTrk";
    tf_name += nTrkCut;
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
    if (p_bcid > 0 && current_bcid != p_bcid) {
      delete h_nVtxpLB;
      continue;
    }
    if (nTrk != nTrkCut) {
      delete h_nVtxpLB;
      continue;
    }

    //now loop over defined pLB
    for (Int_t h_bin_pLB = 1; h_bin_pLB < h_nVtxpLB->GetYaxis()->GetNbins(); h_bin_pLB++) {
      Int_t pLB = h_nVtxpLB->GetYaxis()->GetBinLowEdge(h_bin_pLB); // left-define histogram

      TString h_nVtxName;
      h_nVtxName = "h_nVtx_BCID";
      h_nVtxName += current_bcid;
      h_nVtxName += "_NTrkCut";
      h_nVtxName += nTrkCut;
      h_nVtxName += "_pLB";
      h_nVtxName += pLB;

      TH1D *h_nVtx = h_nVtxpLB->ProjectionX(h_nVtxName, h_bin_pLB, h_bin_pLB);
      if (!h_nVtx) {
        cerr << "WARNING: pLB " << pLB << " not defined in input histogram: " << h_nVtxpLB->GetName() << endl;
        continue;
      }

      //apply corrections, unfold and store results
      UnfoldRawVertexCounts(h_nVtx, p_bcid, nTrkCut, pLB, triggerType);
      delete h_nVtx;
    }
    delete h_nVtxpLB;

  } //end loop over keys of the input file

  f_in->Close();

  //Free memory from unfolding class
  delete m_unfold;
  m_unfold = 0;
}

/// Helper function to get correct mean with left-edged histograms :(
Double_t MyGetMeanError( TH1D *m_h_trueVtx, Double_t *error) {
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

void LumiVtx::UnfoldRawVertexCounts(TH1D *h_nVtx, Int_t p_bcid, Int_t nTrkCut, Int_t pLB, TString triggerType) {
  //Duration of the current pLB
  Float_t current_duration = lb_duration[pLB];//pLB_timestamps[pLB].second - pLB_timestamps[pLB].first;

  //Total number of bunch-crossings in the pLB = Delta_T * f_r * N_bunches
  Long64_t nBC = current_duration * 11245.5 * n_bcids;

  //Total number of triggerable events after prescale (i.e. 100% trigger efficiency after prescale)
  Long64_t nBC_prescaled = nBC / prescale[pLB];

  //Number of triggered events that end in the plot
  // note: if does not enter the plot, we should consider it a 0-vertex anyway (just a check now)
  Long64_t nTrig = h_nVtx->GetEntries();

  //This is the number of non-triggered events, excluding prescale effects
  Long64_t nMissedEvents = nBC_prescaled - nTrig;

  //Live-fraction accounting for dead-time
  //(note that here we're only interested in the livefraction of this pLB,
  //  it's irrelevant what other pLB do, since we get mu_vis right away)
  Double_t liveFraction(1.0);
  liveFraction = live_fraction[pLB];

  cout << "Unfolding " << h_nVtx->GetName() << endl;

  //Storing raw events
  Double_t raw_mean, raw_mean_error;
  raw_mean = MyGetMeanError(h_nVtx, &raw_mean_error);
  n_unf_raw[0][pLB] = raw_mean;
  n_unf_raw_err[0][pLB] = raw_mean_error;

#ifdef DEBUG_LUMIVTX
  cout << " Number of BC (prescaled BC)= " << nBC << "(" << nBC_prescaled << ")"
       << ", livetime: " << liveFraction << ", duration: " << current_duration << endl;
  cout << " Triggered events: " << nTrig << endl;
  cout << " (Old) Raw mean: " << h_nVtx->GetMean() << " +/- " << h_nVtx->GetMeanError() << endl;
  cout << " Raw mean: " << raw_mean << " +/- " << raw_mean_error << endl;
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
  //Note: BGRP7 is a custom alias name. This is a random trigger, typically EF_rd0_Filled_NoAlg
  if (triggerType != "BGRP7" && triggerType != "EF_rd0_Filled_NoAlg") {

    cout << " Trigger correction (0-vertex bin): " << nMissedEvents << " * " << liveFraction << "(lf)"
         << " = " << nMissedEvents*liveFraction << "(" << nMissedEvents*liveFraction / nTrig << ")" << endl;

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

  //store trigger-corrected values
  raw_mean = MyGetMeanError(h_nVtx, &raw_mean_error);
  n_unf_raw[1][pLB] = raw_mean;
  n_unf_raw_err[1][pLB] = raw_mean_error;

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
  n_unf[pLB] = m_unfold->GetMeanInteractions(&errMu);//same as h_nvtx_unfold->GetMean()
  n_unf_err[pLB] = errMu;//h_nvtx_unfold->GetMeanError();

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
  cout << "Unfolded " << h_nVtx->GetName() << "; "
       << "Unfolded mu: " << n_unf[pLB] << " +/- " << n_unf_err[pLB] << endl;
}
#endif


void LumiVtx::CorrectDataRate(TString triggerType) {
  cout << "Correcting for deadtime and prescales..." << endl;

  if (triggerType == "full") {
    #ifdef DEBUG_LUMIVTX
    cout << "Full trigger, correcting for live fraction and prescale" << endl;
    #endif
    for (vector<Int_t>::iterator ilb = lumiblocks.begin(); ilb != lumiblocks.end(); ++ilb) {
      n_vtx[*ilb] = n_vtx_raw[*ilb] / live_fraction[*ilb] * prescale[*ilb];
      n_vtx_err[*ilb] = n_vtx_raw_err[*ilb] / live_fraction[*ilb] * prescale[*ilb];

      n_evt[*ilb] = n_evt_raw[*ilb] / live_fraction[*ilb] * prescale[*ilb];
      n_trig_raw[*ilb] = lb_duration[*ilb] * 11245.5 * live_fraction[*ilb] / prescale[*ilb];
      n_trig[*ilb] = lb_duration[*ilb] * 11245.5;
      float p = n_evt_raw[*ilb] / n_trig_raw[*ilb];
      n_evt_raw_err[*ilb] = TMath::Sqrt(n_trig_raw[*ilb] * p * (1.-p));
      n_evt_err[*ilb] = n_evt_raw_err[*ilb] * live_fraction[*ilb] / prescale[*ilb];
      #ifdef DEBUG_LUMIVTX
      if ((*ilb - 1) % 50 == 0) {
        cout << "LB " << *ilb << ", n_vtx_raw = " << n_vtx_raw[*ilb] << ", n_vtx = " << n_vtx[*ilb] << ", n_trig_raw = " << n_trig_raw[*ilb] << ", n_trig = " << n_trig[*ilb] << endl;
      }
      #endif
    }
  } else if (triggerType == "fixedrate") {
    #ifdef DEBUG_LUMIVTX
    cout << "Fixed rate trigger, scaling counts to full BX rate." << endl;
    #endif
    for (vector<Int_t>::iterator ilb = lumiblocks.begin(); ilb != lumiblocks.end(); ++ilb) {
      Int_t current_nbc = 11245.5 * lb_duration[*ilb] * n_bcids;
      n_vtx[*ilb] = n_vtx_raw[*ilb] * current_nbc / n_trig_raw[*ilb];
      n_vtx_err[*ilb] = n_vtx_raw_err[*ilb] * current_nbc / n_trig_raw[*ilb];
      #ifdef DEBUG_LUMIVTX
      cout << "Set n_vtx[" << *ilb << "] = " << n_vtx[*ilb] << " = " << n_vtx_raw[*ilb] << " * " << current_nbc << " / " << n_trig_raw[*ilb] << endl;
      #endif

      n_evt[*ilb] = n_evt_raw[*ilb] * current_nbc / n_trig_raw[*ilb];
      float p = n_evt_raw[*ilb] / n_trig_raw[*ilb];
      n_evt_raw_err[*ilb] = TMath::Sqrt(n_trig_raw[*ilb] * p * (1.-p));

      n_trig[*ilb] = current_nbc;
      n_evt_err[*ilb] = n_evt_raw_err[*ilb] * current_nbc / n_trig[*ilb];
    }
  } else if (triggerType == "data" ) {
    for (vector<Int_t>::iterator ilb = lumiblocks.begin(); ilb != lumiblocks.end(); ++ilb) {
      Int_t current_nbc = 11245.5 * lb_duration[*ilb] * n_bcids;
      n_vtx[*ilb] = n_vtx_raw[*ilb] * current_nbc / n_trig_raw[*ilb];
      n_vtx_err[*ilb] = n_vtx_raw_err[*ilb] * current_nbc / n_trig_raw[*ilb];
      //n_vtx[*ilb] = n_vtx_raw[*ilb];
      //n_vtx_err[*ilb] = n_vtx_raw_err[*ilb];
      #ifdef DEBUG_LUMIVTX
      cout << "[LumiVtx] Line 698 Set n_vtx[" << *ilb << "] = " << n_vtx[*ilb] << " = " << n_vtx_raw[*ilb] << " * " << current_nbc << " / " << n_trig_raw[*ilb] << endl;
      #endif
      n_evt[*ilb] = n_evt_raw[*ilb];
      n_evt_err[*ilb] = n_evt_raw_err[*ilb];
      n_trig[*ilb] = current_nbc;
    }
  }
}

void LumiVtx::CalculateRates() {
  //Vertex rates: simply divide vertex count by LB duration.
  total_duration = 0.;
  for (vector<Int_t>::iterator ilb = lumiblocks.begin(); ilb != lumiblocks.end(); ++ilb) {
    if (lb_duration[*ilb] > 0.) {
      r_vtx[*ilb] = n_vtx[*ilb] / lb_duration[*ilb];
      cout << "[LumiVtx] INFO : r_vtx[" <<*ilb<< "] = " <<n_vtx[*ilb] << " / " << lb_duration[*ilb] << " = " << r_vtx[*ilb] << endl;
      r_vtx_err[*ilb] = n_vtx_err[*ilb] / lb_duration[*ilb];
      total_duration += lb_duration[*ilb];
    } else {
      r_vtx[*ilb] = 0.;
      r_vtx_err[*ilb] = 0.;
    }
    #ifdef DEBUG_LUMIVTX
    cout << "[LumiVtx] DEBUG : initial r_vtx[" << *ilb << "] = " << r_vtx[*ilb] << " +/- " << r_vtx_err[*ilb] << endl;
    #endif

  }
  total_r_vtx = total_n_vtx / total_duration;
  total_r_vtx_err = total_n_vtx_err / total_duration;

  //LB-level event rates: need to do the analytic poisson correction.
  //The per-lumiblock values are probably nonsense, keep that in mind :)
  for (vector<Int_t>::iterator ilb = lumiblocks.begin(); ilb != lumiblocks.end(); ++ilb) {
    if (n_evt_raw[*ilb] < n_trig_raw[*ilb]) {
      r_vtx_evt[*ilb] = -11245.5 * n_bcids * TMath::Log(1. - n_evt_raw[*ilb] / n_trig_raw[*ilb]);
      r_vtx_evt_err_up[*ilb] = -11245.5 * n_bcids * TMath::Log(1. - (n_evt_raw[*ilb] + n_evt_raw_err[*ilb]) / n_trig_raw[*ilb]) - r_vtx_evt[*ilb];
      r_vtx_evt_err_down[*ilb] = -11245.5 * n_bcids * TMath::Log(1. - (n_evt_raw[*ilb] - n_evt_raw_err[*ilb]) / n_trig_raw[*ilb]) - r_vtx_evt[*ilb];
    } else {
      r_vtx_evt[*ilb] = -1.;
      r_vtx_evt_err_up[*ilb] = 0.;
      r_vtx_evt_err_down[*ilb] = 0.;
    }
    #ifdef DEBUG_LUMIVTX
    cout << "[LumiVtx] DEBUG : r_vtx[" << *ilb << "] = " << r_vtx[*ilb] << " +/- " << r_vtx_err[*ilb] << endl;
    cout << "[LumiVtx] DEBUG : r_vtx_evt[" << *ilb << "] = " << r_vtx_evt[*ilb] << " + " << r_vtx_evt_err_up[*ilb] << " / - " << r_vtx_evt_err_down[*ilb] << endl;
    #endif
  }
  //I omit total_r_vtx_evt, because it's not really useful...

}

void LumiVtx::InitializeMaskingCorrection(TString p_energy, TString p_settings, Int_t p_ntrkcut, TH1D *h_z, TString save_name) {
  /**
    *  h_z is used to generate an expected dz distribution.
    *  Expected dz histogram is saved to include/cache.root : h_dz_expected_<save_tag>.
    */

  cout << "LumiVtx : [INFO] Initializing masking and fake corrections." << endl;

  if (p_ntrkcut == 2) {
    cout << "LumiVtx : [INFO] Skipping masking correction initialization for NTrk2" << endl;
    correct_masking = false;
    return;
  }

  TString hname = "h_pmask_dz_NTrk";
  hname += p_ntrkcut;
  TString tag = "NTrk";
  tag += p_ntrkcut;

  TString s_tag;
  if (p_energy == "7" && p_settings == "17.2-normal") {
    s_tag = "data_7TeV_17.2-normal";
  } else if (p_energy == "7" && p_settings == "17.2-VtxLumi") {
    s_tag = "data_7TeV_17.2-VtxLumi";
  } else if (p_energy == "7" && p_settings == "17.2-normal") {
    s_tag = "data_7TeV_17.2-normal";
  } else if (p_energy == "7" && p_settings == "17.2-VtxLumi") {
    s_tag = "data_7TeV_17.2-VtxLumi";
  } else {
    cerr << "[LumiVtx] ERROR : In masking correction initialization, invalid specifiers: " << endl;
    cerr << "[LumiVtx] ERROR : p_energy = " << p_energy << endl;
    cerr << "[LumiVtx] ERROR : p_settings = " << p_settings << endl;
    cerr << "[LumiVtx] ERROR : Exiting..." << endl;
    exit(1);
  }
  pmc = new PileupMaskingCorrection(s_tag, p_ntrkcut);

  pmc->GenerateDzDistribution(h_z);
  pmc->GenerateCorrection(pmc->GetExpectedDzDistribution());

  TString cache_path = "cache/";
  cache_path += save_name;
  cache_path += ".root";
  TFile *f_cache = new TFile(cache_path, "UPDATE");

  pmc->GetExpectedDzDistribution()->Write("h_dz_expected", TObject::kOverwrite);
  pmc->GetMuMap()->Write("tg_mu_obs_vs_mu_actual", TObject::kOverwrite);
  pmc->GetMuCorrection()->Write("tg_pileup_correction", TObject::kOverwrite);

  f_cache->Close();

  correct_masking = true;
}

void LumiVtx::InitializeMaskingCorrection(TString p_energy, TString p_settings, Int_t p_ntrkcut, TString save_name) {
  /**
    *  Load some time-intensive histograms from cache/<save_name>.root .
    *  dz histograms saved/loaded as h_dz_expected_<save_tag>
    */

  if (p_ntrkcut == 2) {
    cout << "LumiVtx : [INFO] Skipping masking correction initialization for NTrk2" << endl;
    correct_masking = false;
    return;
  }

  cout << "LumiVtx : [INFO] Initializing masking and fake correction, using cached masking correction." << endl;

  TString hname = "h_pmask_dz_NTrk";
  hname += p_ntrkcut;
  TString tag = "NTrk";
  tag += p_ntrkcut;
  TString s_tag;
  if (p_energy == "7" && p_settings == "17.2-normal") {
    s_tag = "data_7TeV_17.2-normal";
  } else if (p_energy == "7" && p_settings == "17.2-VtxLumi") {
    s_tag = "data_7TeV_17.2-VtxLumi";
  } else if (p_energy == "8" && p_settings == "17.2-normal") {
    s_tag = "data_8TeV_17.2-normal";
  } else if (p_energy == "8" && p_settings == "17.2-VtxLumi") {
    s_tag = "data_8TeV_17.2-VtxLumi_207216";
  } else {
    cerr << "[LumiVtx] ERROR : In masking correction initialization, invalid specifiers: " << endl;
    cerr << "[LumiVtx] ERROR : p_energy = " << p_energy << endl;
    cerr << "[LumiVtx] ERROR : p_settings = " << p_settings << endl;
    cerr << "[LumiVtx] ERROR : Exiting..." << endl;
    exit(1);
  }
  pmc = new PileupMaskingCorrection(s_tag, p_ntrkcut);

  TString cache_path = "cache/";
  cache_path += save_name;
  cache_path += ".root";
  TFile *f_cache = new TFile(cache_path, "READ");
  if (!(f_cache->IsOpen())) {
    cerr << "LumiVtx : [ERROR] Cached masking correction at " << cache_path << " was not found! Exiting." << endl;
    exit(1);
  }

  //  TH1F *h_dz_expected = (TH1F*)f_cache->Get("h_dz_expected");
  TGraphErrors *tg_mu_obs_vs_mu_actual = (TGraphErrors*)f_cache->Get("tg_mu_obs_vs_mu_actual");
  TGraphErrors *tg_pileup_correction = (TGraphErrors*)f_cache->Get("tg_pileup_correction");
  pmc->LoadCorrection(tg_mu_obs_vs_mu_actual, tg_pileup_correction);
  f_cache->Close();

  correct_masking = true;
}

void LumiVtx::InitializeFakeCorrection(TString p_energy, TString p_settings, Int_t p_nTrkCut) {
	if( p_energy == "7" ) {
		if (p_settings == "17.2-normal") {
			if (systematic_uncertainty_list["bs-45mm"]) {
				fc = new FakeCorrection("mc_7TeV_17.2_normal_pythia8_pu_bs45", p_nTrkCut);
			} else if (systematic_uncertainty_list["bs-55mm"]) {
				fc = new FakeCorrection("mc_7TeV_17.2_normal_pythia8_pu_bs55", p_nTrkCut);
			}
		} else if (p_settings == "17.2-VtxLumi") {
			fc = new FakeCorrection("mc_7TeV_17.2_VtxLumi_pythia8_pu", p_nTrkCut);
		}
	} else if( p_energy == "8" ) {
		if (p_settings == "17.2-normal") {
			fc = new FakeCorrection("mc_8TeV_17.2_normal_pythia8_pu", p_nTrkCut);
		} else if (p_settings == "17.2-VtxLumi") {
			fc = new FakeCorrection("mc_8TeV_17.2_VtxLumi_2newsets", p_nTrkCut);
		}
	}

  // Exit if FakeCorrection has not been set
	if( !fc ) {
		cerr << "[LumiVtx] ERROR : Invalid specifiers in fake correction initialization:" << endl;
		cerr << "[LumiVtx] ERROR : p_energy = " << p_energy << endl;
		cerr << "[LumiVtx] ERROR : p_settings = " << p_settings << endl;
		cerr << "[LumiVtx] ERROR : p_nTrkCut = " << p_nTrkCut << endl;
		cerr << "[LumiVtx] ERROR : Exiting..." << endl;
		exit(1);
	}

correct_fakes = true;

}

void LumiVtx::SetMuScale(Float_t p_mu_scale) {
  cout << "[LumiVtx] INFO : Setting fake correction mu scale to " << p_mu_scale << endl;
  fc->SetMuScale(p_mu_scale);
}

void LumiVtx::CorrectPileupEffects(TH2D *h_z_plb, TString p_energy, TString p_settings, Int_t p_ntrkcut, TString p_run) {
  #ifdef DEBUG_LUMIVTX
  cout << "LumiVtx::CorrectPileupEffects()" << endl;
  #endif

  /**
    *  A note on the procedure:
    *  Since fakes participate in the masking process, apply the masking correction first.
    *  Both corrections are multiplicative, so hopefully don't need to worry too much about order.
    */

  ofstream textfile;
  stringstream ss;
  ss << "/afs/hep.man.ac.uk/u/julia/vertex_lumi/mu_being_corrected_run" << p_run << "_NTrkCut" << p_ntrkcut <<".txt";
  string filename = ss.str();
  textfile.open(filename.c_str());
  textfile.precision(11);
  textfile << "pLb mu_raw mu_raw_err initial_masking_correction_factor mu_fake mu_fake_err mu_real   final_masking_correction_factor   mu_vis" << endl;
  for (vector<Int_t>::iterator ilb = lumiblocks.begin(); ilb != lumiblocks.end(); ++ilb) {

    cout << "[LumiVtx] INFO : Correcting pileup for LB " << *ilb << endl;

    // -- Vertex counting
    Float_t mu_raw = r_vtx[*ilb] / 11245.5 / n_bcids;
    Float_t mu_raw_err = r_vtx_err[*ilb] / 11245.5 / n_bcids;

    Int_t bin = h_z_plb->GetYaxis()->FindBin(*ilb);
    TString hname = "h_z_plb";
    hname += *ilb;
    TH1D *h_z = (TH1D*)h_z_plb->ProjectionX(hname, bin, bin);
    TString p_tag = "data_";
    p_tag += p_energy;
    p_tag += "TeV_";
    p_tag += p_settings;
    p_tag += "_207216";
    PileupMaskingCorrection *mc;
    Float_t initial_masking_correction_factor = 1.;
    cout << "[LumiVtx] Line 920: p_tag = " << p_tag << endl;
    if (p_ntrkcut != 2) {
      mc = new PileupMaskingCorrection(p_tag, p_ntrkcut);
      if (systematic_uncertainty_list["masking_toy_scaling"]) {
        mc->SetPmaskScale(1.02);
      }
      TString pmc_tag = "";
      pmc_tag += p_energy;
      pmc_tag += "_";
      pmc_tag += p_settings;
      pmc_tag += "_";
      pmc_tag += p_ntrkcut;
      pmc_tag += "_";
      pmc_tag += *ilb;
      cout << "[LumiVtx] Line 934: pmc_tag = " << pmc_tag << endl;
      mc->GenerateDzDistribution(h_z, pmc_tag);
      mc->GenerateCorrection(mc->GetExpectedDzDistribution());
      initial_masking_correction_factor = mc->GetCorrectionFactor(mu_raw);
      
    }

    if (initial_masking_correction_factor < 1.) {
      cout << "[LumiVtx] WARNING : initial masking correction factor = " << initial_masking_correction_factor << " < 1. Setting to 1." << endl;
      initial_masking_correction_factor = 1.;
    }

    Double_t mu_fake(0.);
    if (systematic_uncertainty_list["fake_low"]) {
      mu_fake = fc->GetFakeMuLowFromMuReconMC(mu_raw * initial_masking_correction_factor);
    } else if (systematic_uncertainty_list["fake_high"]) {
      mu_fake = fc->GetFakeMuHighFromMuReconMC(mu_raw * initial_masking_correction_factor);
    } else {
      mu_fake = fc->GetFakeMuFromMuReconMC(mu_raw * initial_masking_correction_factor);
    }
    Double_t mu_fake_err = fc->GetFakeMuUncertaintyFromMuReconMC(mu_raw * initial_masking_correction_factor);
    std::cout << "For mu_raw = " << mu_raw << " we have mu_fake = " << mu_fake << std::endl;
    if (mu_fake < 0.) {
      cout << "[LumiVtx] WARNING : mu_fake = " << mu_fake << " outside of range. Setting to 0." << endl;
      mu_fake = 0.;
    }
    if (mu_fake_err < 0.) {
      mu_fake_err = 0.;
    }
    
    Double_t mu_real = mu_raw - mu_fake;

    Double_t final_masking_correction_factor = 1.;;
    if (p_ntrkcut != 2) {
      final_masking_correction_factor = mc->GetCorrectionFactor(mu_real);
    }
    if (final_masking_correction_factor < 1.) {
      std::cout << "[LumiVtx] DEBUG : masking correction factor would be " << final_masking_correction_factor << " : setting to 1 instead" << std::endl;
      final_masking_correction_factor = 1.;
    }

    Double_t mu_vis = mu_real * final_masking_correction_factor;
    Double_t mu_vis_err = TMath::Sqrt(TMath::Power(mu_raw_err, 2) + TMath::Power(mu_fake_err, 2)) * final_masking_correction_factor;
    cout << "[LumiVtx] DEBUG : mu_vis_nvtx = " << mu_vis << " = (" << mu_raw << " - " << mu_fake << ") * " << final_masking_correction_factor << endl;

    textfile << *ilb << " " << mu_raw << " " << mu_raw_err << " " << initial_masking_correction_factor << " " << mu_fake << " " << mu_fake_err << " " << mu_real << " " << final_masking_correction_factor << " " << mu_vis << " " << mu_vis_err << endl;

    r_vtx[*ilb] = mu_vis * 11245.5 * n_bcids;
    cout << "[LumiVtx] Line 982 r_vtx[*ilb] = " << mu_vis << " * 11245.5 * " << n_bcids << " = " << r_vtx[*ilb] << endl;
    r_vtx_err[*ilb] = mu_vis_err * 11245.5 * n_bcids;
    n_vtx[*ilb] = r_vtx[*ilb] * lb_duration[*ilb];
    n_vtx_err[*ilb] = r_vtx_err[*ilb] * lb_duration[*ilb];

    mu_fake_list_nvtx[*ilb] = mu_fake;
    masking_corrections[*ilb] = final_masking_correction_factor;
    cout << "[LumiVtx] Line 983, masking_corrections["<<*ilb<<"] = " << masking_corrections[*ilb] << endl;

    // -- Event counting
    mu_raw = r_vtx_evt[*ilb] / 11245.5 / n_bcids;
    Float_t mu_raw_up = (r_vtx_evt[*ilb] + r_vtx_evt_err_up[*ilb]) / 11245.5 / n_bcids;
    Float_t mu_raw_down = (r_vtx_evt[*ilb] - r_vtx_evt_err_down[*ilb]) / 11245.5 / n_bcids;

    mu_fake = fc->GetMuFakeFromMuReconNEvt(mu_raw);
    mu_fake_err = fc->GetMuFakeUncertaintyFromMuReconNEvt(mu_raw);

    if (mu_fake < 0.) {
      cout << "[LumiVtx] WARNING : mu_fake = " << mu_fake << " < 0. Setting to 0." << endl;
      mu_fake = 0.;
      mu_fake_err = 0.;
    }

    mu_vis = mu_raw - mu_fake;

    cout << "[LumiVtx] DEBUG : mu_vis_evt = " << mu_vis << " = " << mu_raw << " - " << mu_fake << endl;
    //mu_vis = mu_raw * (1. - fake_fraction);
    Float_t mu_vis_up = TMath::Sqrt(TMath::Power(mu_raw_up, 2) + TMath::Power(mu_fake_err, 2));
    Float_t mu_vis_down = TMath::Sqrt(TMath::Power(mu_raw_down, 2) + TMath::Power(mu_fake_err, 2));

    r_vtx_evt[*ilb] = mu_vis * 11245.5 * n_bcids;
    r_vtx_evt_err_up[*ilb] = mu_vis_up * 11245.5 * n_bcids - r_vtx_evt[*ilb];
    r_vtx_evt_err_down[*ilb] = r_vtx_evt[*ilb] - mu_vis_down * 11245.5 * n_bcids;

    mu_fake_list_nevt[*ilb] = mu_fake;
    
    //delete mc;
  }
  textfile.close();
}


void LumiVtx::CorrectMasking() {
#ifdef DEBUG_LUMIVTX
  cout << "LumiVtx::CorrectPileupMasking()" << endl;
#endif

  for (vector<Int_t>::iterator ilb = lumiblocks.begin(); ilb != lumiblocks.end(); ++ilb) {
    Float_t current_mu = r_vtx[*ilb] / 11245.5;
    Float_t current_correction_factor = pmc->GetCorrectionFactor(current_mu);
    if (current_correction_factor < 1.) {
      cout << "WARNING: PileupMaskingCorrection object return a correction factor of less than 1! Manually setting correction to 1." << endl;
      current_correction_factor = 1.;
    }

    n_vtx_masked[*ilb] = n_vtx[*ilb] * (current_correction_factor - 1.);
    if (n_vtx_masked[*ilb] > 0.) {
      n_vtx[*ilb] += n_vtx_masked[*ilb];
      cout << "[LumiVtx] Line 1038 r_vtx[*ilb] = " << r_vtx[*ilb]<< " * " << current_correction_factor << " = " << r_vtx[*ilb] << endl;
      r_vtx[*ilb] = r_vtx[*ilb] * current_correction_factor;
      r_vtx_err[*ilb] = r_vtx_err[*ilb] * current_correction_factor;
    }
  }
}


void LumiVtx::ApplyMcCorrection(TString path, TString objectname) {

  cout << "Applying MC-based correction" << endl;

  TFile *f_in = new TFile(path, "READ");
  TGraphErrors *tg_mc_correction = (TGraphErrors*)f_in->Get(objectname);

  for (vector<Int_t>::iterator ilb = lumiblocks.begin(); ilb != lumiblocks.end(); ++ilb) {
    Float_t mu_uncorr = r_vtx[*ilb] / 11245.5;
    Float_t mu_corr = tg_mc_correction->Eval(mu_uncorr);
    r_vtx[*ilb] = mu_corr * 11245.5;
    //cout << "[LumiVtx] Line 1056 r_vtx[*ilb] = " << mu_corr << " * 11245.5 = " << r_vtx[*ilb] << endl; Doesn't call this
  }
}


void LumiVtx::CalculateAllLbLum() {
  //Vertex counting
  for (vector<Int_t>::iterator ilb = lumiblocks.begin(); ilb != lumiblocks.end(); ++ilb) {
    lb_inst_lum[*ilb] = r_vtx[*ilb] / sigmavis[GlobalSettings::kVtxC];
    lb_inst_lum_err[*ilb] = r_vtx_err[*ilb] / sigmavis[GlobalSettings::kVtxC];

    #ifdef DEBUG_LUMIVTX
    cout << "Set lb_inst_lum[" << *ilb << "] = " << lb_inst_lum[*ilb] << endl;
    #endif

    lb_total_lum[*ilb] = lb_inst_lum[*ilb] * lb_duration[*ilb];
    lb_total_lum_err[*ilb] = lb_inst_lum_err[*ilb] * lb_duration[*ilb];
    #ifdef DEBUG_LUMIVTX
    //cout << "[debug] LB = " << *ilb << " / Linst = " << lb_inst_lum[*ilb] << " / Linst_ofl_pref = " << ofl_lum[*ilb] << endl;
    #endif
  }

  //Event counting.
  for (vector<Int_t>::iterator ilb = lumiblocks.begin(); ilb != lumiblocks.end(); ++ilb) {
    lb_inst_lum_evt[*ilb] = r_vtx_evt[*ilb] / sigmavis[GlobalSettings::kEvtC];
    lb_inst_lum_evt_err_up[*ilb] = r_vtx_evt_err_up[*ilb] / sigmavis[GlobalSettings::kEvtC];
    lb_inst_lum_evt_err_down[*ilb] = r_vtx_evt_err_down[*ilb] / sigmavis[GlobalSettings::kEvtC];

    lb_total_lum_evt[*ilb] = lb_inst_lum_evt[*ilb] * lb_duration[*ilb];
    lb_total_lum_evt_err_up[*ilb] = lb_inst_lum_evt_err_up[*ilb] * lb_duration[*ilb];
    lb_total_lum_evt_err_down[*ilb] = lb_inst_lum_evt_err_down[*ilb] * lb_duration[*ilb];
  }

  for (vector<Int_t>::iterator ilb = lumiblocks.begin(); ilb != lumiblocks.end(); ++ilb) {
    lb_inst_lum_unf[*ilb] = lb_inst_lum_unf_err[*ilb] = 0.0;
    lb_total_lum_unf[*ilb] = lb_total_lum_unf_err[*ilb] = 0.0;
  }
  #ifdef USE_UNFOLD
  //Unfolding
  if (sigmavis.find(GlobalSettings::kUnfC) == sigmavis.end()) {
    cerr << "[LumiVtx] ERROR Requested lumi from unfolding, but sigma_vis values are not loaded.";
    return;
  }
  for (vector<Int_t>::iterator ilb = lumiblocks.begin(); ilb != lumiblocks.end(); ++ilb) {
    //LB inst lumi: \cal{L} = mu_vis * f_r * N_bunches / sigma_vis
    lb_inst_lum_unf[*ilb] = n_unf[*ilb] * 11245.5 * n_bcids / sigmavis[GlobalSettings::kUnfC];
    lb_inst_lum_unf_err[*ilb] = n_unf_err[*ilb] * 11245.5 * n_bcids / sigmavis[GlobalSettings::kUnfC];
    lb_inst_lum_unf_err[*ilb] = TMath::Sqrt(TMath::Power(lb_inst_lum_unf_err[*ilb],2) +
                                            TMath::Power(sigmavis_err[GlobalSettings::kUnfC],2));

    //LB total lumi: L = \cal{L}*Delta_T
    lb_total_lum_unf[*ilb] = lb_inst_lum_unf[*ilb] * lb_duration[*ilb];
    lb_total_lum_unf_err[*ilb] = lb_inst_lum_unf_err[*ilb] * lb_duration[*ilb];
  }
  #endif
}

void LumiVtx::CalculateAllRunLum() {
  run_trig_lum = 0;
  run_trig_lum_err = 0;
  for (vector<Int_t>::iterator ilb = lumiblocks.begin(); ilb != lumiblocks.end(); ++ilb) {
          run_trig_lum += lb_inst_lum[*ilb] * lb_duration[*ilb];
          run_trig_lum_err += TMath::Power(lb_inst_lum_err[*ilb] * lb_duration[*ilb], 2);
  }
  run_trig_lum_err = TMath::Sqrt(run_trig_lum_err);
/*  //Vertex counting
  std::cout << "CalculateAllRunLum "<<total_n_vtx << ", " << sigmavis[GlobalSettings::kVtxC] << std::endl;
  run_trig_lum = total_n_vtx / sigmavis[GlobalSettings::kVtxC];
  run_trig_lum_err = total_n_vtx_err / sigmavis[GlobalSettings::kVtxC];

  //Event counting
  run_trig_lum_evt = total_n_vtx_evt / sigmavis[GlobalSettings::kEvtC];
  run_trig_lum_evt_err_up = total_n_vtx_evt_err_up / sigmavis[GlobalSettings::kEvtC];
  run_trig_lum_evt_err_down = total_n_vtx_evt_err_down / sigmavis[GlobalSettings::kEvtC];

  //Unfolding
  run_total_lum_unf = run_total_lum_unf_err = 0.0;
  run_total_ofl_lum_unf = run_total_ofl_lum_unf_err = 0.0; //don't have enough info to get there yet..
#ifdef USE_UNFOLD
  for (vector<Int_t>::iterator ilb = lumiblocks.begin(); ilb != lumiblocks.end(); ++ilb) {
    run_total_lum_unf += lb_total_lum_unf[*ilb];
    //error is really statistical up to here. Just add in quadrature.
    run_total_lum_unf_err += lb_total_lum_unf_err[*ilb] * lb_total_lum_unf_err[*ilb];
  }
  run_total_lum_unf_err = TMath::Sqrt(run_total_lum_unf_err);
#endif
*/
}
/*
Float_t LumiVtx::PileupCorrect(Float_t n_vtx_obs) {
  return tg_pu_corr->Eval(n_vtx_obs);
}
*/

void LumiVtx::SetNBunches(Int_t p_n_bcids) {

  n_bcids = p_n_bcids;

}

void LumiVtx::SetUnfResponseMatrix(TString fileName, TString histoBaseName) {
  responseMatrixFileName = fileName;
  responseMatrixHName = histoBaseName;
  cout << "[LumiVtx] INFO Selected Response Matrix " << responseMatrixHName
       << " from file " << responseMatrixFileName << endl;
}


Double_t LumiVtx::GetLbInstLum(Int_t n, GlobalSettings::AnalysisMethods method) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NVtxRaw for invalid iterator position: " << n << endl;
    return -1;
  }

  if (method == GlobalSettings::kVtxC) {
    cout << "Line 1178, LumiVtx, requested method was VtxCounting" << endl;
    Int_t current_lb = *(lumiblocks.begin() + n);
    if (lb_inst_lum[current_lb] > 0.) {
      return lb_inst_lum[current_lb];
    } else {
      return 0.;
    }
  } else if (method == GlobalSettings::kEvtC) {
    return GetEvtLbInstLum(n);
  } else if (method == GlobalSettings::kUnfC) {
    Int_t current_lb = *(lumiblocks.begin() + n);
    if (lb_inst_lum[current_lb] > 0.) {
      return lb_inst_lum_unf[current_lb];
    } else {
      return 0.;
    }
  }

  cerr << "[LumiVtx::GetLbInstLum] ERROR Unknown method: " << method << endl;
  return 0.0;
}

Double_t LumiVtx::GetLbInstLumError(Int_t n, GlobalSettings::AnalysisMethods method) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NVtxRaw for invalid iterator position: " << n << endl;
    return -1;
  }

  if (method == GlobalSettings::kVtxC) {
    Int_t current_lb = *(lumiblocks.begin() + n);
    if (lb_inst_lum[current_lb] > 0.) {
      return lb_inst_lum_err[current_lb];
    } else {
      return 0.;
    }
  } else if (method == GlobalSettings::kEvtC) {
    //return average of up and down error (use specific function to get both)
    std::pair<Double_t, Double_t> err;
    err = GetEvtLbInstLumError(n);
    return (err.first + err.second)/2;
  } else if (method == GlobalSettings::kUnfC) {
    Int_t current_lb = *(lumiblocks.begin() + n);
    if (lb_inst_lum[current_lb] > 0.) {
      return lb_inst_lum_unf_err[current_lb];
    } else {
      return 0.;
    }
  }

  cerr << "[LumiVtx::GetLbInstLumError] ERROR Unknown method: " << method << endl;
  return 0.0;
}

Double_t LumiVtx::GetLbTotalLum(Int_t n, GlobalSettings::AnalysisMethods method) {
  if ((n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0)) {
    cout << "WARNING: Requested NVtxRaw for invalid iterator position: " << n << endl;
    return -1;
  }

  if (method == GlobalSettings::kVtxC) {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return lb_total_lum[current_lb];
  } else if (method == GlobalSettings::kEvtC) {
    return GetEvtLbTotalLum(n);
  } else if (method == GlobalSettings::kUnfC) {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return lb_total_lum_unf[current_lb];
  }

  cerr << "[LumiVtx::GetLbTotalLum] ERROR Unknown method: " << method << endl;
  return 0.0;

}

Double_t LumiVtx::GetLbTotalLumError(Int_t n, GlobalSettings::AnalysisMethods method) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NVtxRaw for invalid iterator position: " << n << endl;
    return -1;
  }

  if (method == GlobalSettings::kVtxC) {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return lb_total_lum_err[current_lb];
  } else if (method == GlobalSettings::kEvtC) {
    //return average of up and down error. Use specific function to access them separately
    std::pair<Double_t, Double_t> err;
    err = GetEvtLbTotalLumError(n);
    return (err.first + err.second)/2;
  } else if (method == GlobalSettings::kUnfC) {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return lb_total_lum_unf_err[current_lb];
  }

  cerr << "[LumiVtx::GetLbTotalLumError] ERROR Unknown method: " << method << endl;
  return 0.0;

}

Double_t LumiVtx::GetRunTotalLum(GlobalSettings::AnalysisMethods method) {
  if (method == GlobalSettings::kVtxC) {
    return run_trig_lum;
  } else if (method == GlobalSettings::kEvtC) {
    return GetEvtRunTotalLum();
  } else if (method == GlobalSettings::kUnfC) {
    return run_total_lum_unf;
  }

  cerr << "[LumiVtx::GetRunTotalLum] ERROR Unknown method: " << method << endl;
  return 0.0;

}

Double_t LumiVtx::GetRunTotalLumError(GlobalSettings::AnalysisMethods method) {
  if (method == GlobalSettings::kVtxC) {
    return run_trig_lum_err;
  } else if (method == GlobalSettings::kEvtC) {
    //return average of up and down error. Use specific function to access them separately
    std::pair<Double_t, Double_t> err;
    err = GetEvtRunTotalLumError();
    return (err.first + err.second)/2;
  } else if (method == GlobalSettings::kUnfC) {
    return run_total_lum_unf_err;
  }

  cerr << "[LumiVtx::GetRunTotalLumError] ERROR Unknown method: " << method << endl;
  return 0.0;

}

Double_t LumiVtx::GetEvtLbInstLum(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NVtxRaw for invalid iterator position: " << n << endl;
    return -1;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    if (lb_inst_lum_evt[current_lb] > 0.) {
      return lb_inst_lum_evt[current_lb];
    } else {
      return 0.;
    }
  }
}

std::pair<Double_t, Double_t> LumiVtx::GetEvtLbInstLumError(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NVtxRaw for invalid iterator position: " << n << endl;
    return std::make_pair<Double_t, Double_t>(-1.,-1.);
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    if (lb_inst_lum_evt_err_down[current_lb] > 0. && lb_inst_lum_evt_err_up[current_lb] > 0.) {
      return std::make_pair<Double_t, Double_t>(lb_inst_lum_evt_err_down[current_lb], lb_inst_lum_evt_err_up[current_lb]);
    } else {
      return std::make_pair<Double_t, Double_t>(0., 0.);
    }
  }
}

Double_t LumiVtx::GetEvtLbTotalLum(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NVtxRaw for invalid iterator position: " << n << endl;
    return -1;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return lb_total_lum_evt[current_lb];
  }
}

std::pair<Double_t, Double_t> LumiVtx::GetEvtLbTotalLumError(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NVtxRaw for invalid iterator position: " << n << endl;
    return std::make_pair<Double_t, Double_t>(-1., -1.);
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return std::make_pair<Double_t, Double_t>(lb_total_lum_evt_err_down[current_lb], lb_total_lum_evt_err_up[current_lb]);
  }
}

Double_t LumiVtx::GetEvtRunTotalLum() {
  return run_trig_lum_evt;
}

std::pair<Double_t, Double_t> LumiVtx::GetEvtRunTotalLumError() {
  return std::make_pair<Double_t, Double_t>(run_trig_lum_evt_err_down, run_trig_lum_evt_err_up);
}

Double_t LumiVtx::GetNVtxMasked(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NVtxMasked for invalid iterator position: " << n << endl;
    return -1.;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return n_vtx_masked[current_lb];
  }
}

Double_t LumiVtx::GetNVtx(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NVtxMasked for invalid iterator position: " << n << endl;
    return -1.;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return n_vtx[current_lb];
  }
}

Double_t LumiVtx::GetNEvt(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NVtxMasked for invalid iterator position: " << n << endl;
    return -1.;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return n_evt[current_lb];
  }
}

Double_t LumiVtx::GetNVtxFromNEvt(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NVtxMasked for invalid iterator position: " << n << endl;
    return -1.;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return r_vtx_evt[current_lb] * lb_duration[current_lb];
  }
}

Int_t LumiVtx::GetFirstLb() {
  return *(lumiblocks.begin());
}

Int_t LumiVtx::GetLastLb() {
  return *(lumiblocks.end());
}

Int_t LumiVtx::GetNLb() {
  return distance(lumiblocks.begin(), lumiblocks.end());
}

Int_t LumiVtx::GetLbFromDist(Int_t n) {
  //if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
  //  return -1;
  //} else {
  //  return *(lumiblocks.begin() + n);
  //}
  return lumiblocks.at(n);
}

Double_t LumiVtx::GetLbDuration(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested duration for invalid iterator position: " << n << endl;
    return -1;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return lb_duration[current_lb];
  }
}

Double_t LumiVtx::GetTotalDuration() {
  return total_duration;
}

Double_t LumiVtx::GetTotalNVtxRaw() {
  return total_n_vtx_raw;
}

Double_t LumiVtx::GetTotalNVtx() {
  return total_n_vtx;
}

Double_t LumiVtx::GetLbNVtxRaw(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NVtxRaw for invalid iterator position: " << n << endl;
    return -1;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return n_vtx_raw[current_lb];
  }
}

Double_t LumiVtx::GetLbNVtx(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NVtx for invalid iterator position: " << n << endl;
    return -1;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return n_vtx[current_lb];
  }
}

Double_t LumiVtx::GetTotalNEvtRaw() {
  return total_n_evt_raw;
}

Double_t LumiVtx::GetTotalNEvt() {
  return total_n_evt_raw * total_n_trig / total_n_trig_raw;
}

Double_t LumiVtx::GetLbNEvtRaw(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NEvtRaw for invalid iterator position: " << n << endl;
    return -1;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return n_evt_raw[current_lb];
  }
}

Double_t LumiVtx::GetLbNEvt(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NEvt for invalid iterator position: " << n << endl;
    return -1;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return n_evt_raw[current_lb] * n_trig[current_lb] / n_trig_raw[current_lb];
  }
}

Double_t LumiVtx::GetTotalNVtxEvt() {
  return total_n_vtx_evt;
}

Double_t LumiVtx::GetLbNVtxEvt(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NEvt for invalid iterator position: " << n << endl;
    return -1;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return r_vtx_evt[current_lb] * lb_duration[current_lb];
  }
}

Double_t LumiVtx::GetTotalNTotalRaw() {
  return total_n_trig_raw;
}

Double_t LumiVtx::GetTotalNTotal() {
  return total_n_trig;
}

Double_t LumiVtx::GetLbNTotalRaw(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NTotalRaw for invalid iterator position: " << n << endl;
    return -1;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return n_trig_raw[current_lb];
  }
}

Double_t LumiVtx::GetLbNTotal(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NTotal for invalid iterator position: " << n << endl;
    return -1;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return n_trig[current_lb];
  }
}

Double_t LumiVtx::GetOflLum() {
  return total_ofl_lum;
}

Double_t LumiVtx::GetOflLum(Int_t n) {
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NEvt for invalid iterator position: " << n << endl;
    return -1;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return ofl_lum[current_lb];
  }
}

Double_t LumiVtx::GetMaskingCorrection(Int_t n) {
  Double_t return_value;
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NEvt for invalid iterator position: " << n << endl;
    return_value = -1;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    return_value = masking_corrections[current_lb];
  }
  return return_value;
}

Double_t LumiVtx::GetFakeMu(Int_t n, TString p_method) {
  Double_t return_value;
  if ( (n > distance(lumiblocks.begin(), lumiblocks.end())) || (n<0) ) {
    cout << "WARNING: Requested NEvt for invalid iterator position: " << n << endl;
    return_value = -1;
  } else {
    Int_t current_lb = *(lumiblocks.begin() + n);
    if (p_method == "NVtx") {
      return_value = mu_fake_list_nvtx[current_lb];
    } else if (p_method == "NEvt") {
      return_value = mu_fake_list_nevt[current_lb];
    } else {
      cout << "LumiVtx::GetFakeCorrection: [WARNING] must specify NVtx or NEvt." << endl;
      return_value = -1;
    }
  }
  return return_value;
}

void LumiVtx::SetSystematicUncertaintyFlag(TString p_systematic_name) {

  cout << "[LumiVtx] INFO : Setting systematic uncertainty flag " << p_systematic_name << endl;
  systematic_uncertainty_list[p_systematic_name] = true;

}



void LumiVtx::SetVerbose(Int_t level) {
  verbose = level;
}

/**************************************************************************/
// Auxilliary functions
/**************************************************************************/

Double_t fitFuncSimpleGaussian(Double_t *x, Double_t *par) {
  return par[0] / (TMath::Sqrt(2. * TMath::Pi()) * par[1]) * TMath::Exp(-1. * TMath::Power(x[0] - par[2], 2) / (2 * par[1] * par[1]));
}

