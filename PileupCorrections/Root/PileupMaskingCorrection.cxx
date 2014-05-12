#define DEBUG_PMC
#define SAVE_RESIDUALS
//#define FITWITHGAUSSIAN

#ifndef PileupMaskingCorrection_cxx
#define PileupMaskingCorrection_cxx

#include "GlobalSettings/GlobalSettings.h"
#include "PileupCorrections/PileupMaskingCorrection.h"

using namespace std;

/* Fit with gaussian */
Double_t fitFunc_gaussian(Double_t *x, Double_t *par);
Double_t fitFunc_gaussian_exclude(Double_t *x, Double_t *par);

/* Fit with a dz template */
Double_t fitFunc_generic(Double_t *x, Double_t *par);
Double_t fitFunc_generic_exclude(Double_t *x, Double_t *par);
TH1D *dz_template;

Double_t meanNObs(Double_t p_mask, Double_t n_gen);
Double_t meanMuObs(Double_t x, Double_t par);
bool reject;

/***************************************************************************************/
/// - Method 1: generate new p_mask vs. dz histogram from a low pileup dz distribution
/***************************************************************************************/

PileupMaskingCorrection::PileupMaskingCorrection(TH1D *h1, TString string1) {
  /// - Constructor for generating a new differential pmask distribution from a low-pileup h_dz histogram.

  cout << "[PileupMaskingCorrection] INFO : Initializing PileupMaskingCorrection (recreate mode)" << endl;

  redo_cache = true;
  pmask_available = false;
  is_MC = kFALSE;

  h_dz = h1;
  tag = string1;
  rebin_factor = 20; 
  exclude_dz = 20.;
  finished = false;
  max_dz = 300.;
  do_ll = true;
  MaxNGenInt = 110;

  TString hname = "h_dz";
  h_dz->SetName(hname);
  h_dz->Sumw2();

  low_stats = (h_dz->Integral() > 10000 ? false : true);
  cout << "[PileupMaskingCorrection] INFO : low_stats? " << h_dz->Integral() << " > 10000 ?" << low_stats << endl;
  //cout << "[PileupMaskingCorrection] INFO : For the purpose of this test, setting it to false anyway" << endl;
  //low_stats = false;
  pmask_scale = 1.;

}

void PileupMaskingCorrection::GenerateDzDistribution(TH1D *h_z, TString tag) {

  cout << "[PileupMaskingCorrection] INFO : Generating expected delta Z distribution from input z distribution." << endl;

  h_dz_expected = new TH1D(TString("h_dz_expected_") + tag, TString("h_dz_expected_") + tag, h_z->GetXaxis()->GetNbins(), h_z->GetXaxis()->GetXmin(), h_z->GetXaxis()->GetXmax());

  for (Int_t bin=0; bin<=h_dz_expected->GetXaxis()->GetNbins(); bin++) {
    h_dz_expected->SetBinContent(bin, 0);
    h_dz_expected->SetBinError(bin, 0);
  }

  //h_dz_expected->Rebin(rebin_factor);

  Int_t trials = 1000000;
  h_z->Scale(1./h_z->Integral());

  for (int i=0; i<trials; i++) {
    Float_t z1 = h_z->GetRandom();
    Float_t z2 = h_z->GetRandom();
    h_dz_expected->Fill(z1 - z2);
  }

  dz_template = h_dz_expected;

  TString path_z = "";

  if (is_MC == 0){
    path_z = "data_h_dz_expected_"; path_z+=tag; path_z+=".root";    
  }else{
    path_z = "MC_h_dz_expected_"; path_z+=tag; path_z+=".root";    
  }
  
  TFile *f_z;
  f_z = new TFile(path_z,"RECREATE");
  h_dz_expected->Scale(1./h_dz_expected->Integral());
  h_dz_expected->Write();
  f_z->Close();
  
}

void PileupMaskingCorrection::LoadDzDistribution(TH1D *h_dz_in) {
  h_dz_expected = (TH1D*)h_dz_in->Clone();
  h_dz_expected->SetName(TString("h_dz_expected_") + tag);
  dz_template = h_dz_expected;
}

void PileupMaskingCorrection::GenerateNewPmask() {

  FitExcluded();
  FitExcludedResiduals();
  MakeFullGaussian();
  MakeDifferentialPmask();

  pmask_available = true;

}

void PileupMaskingCorrection::FitExcluded() {
  if (!h_dz_expected) {
    cerr << "ERROR: you need to set dz_template before using it to fit. Exiting..." << endl;
    exit(1);
  }
  if (!low_stats) {
    #ifdef DEBUG_PMC
    cout << "[PileupMaskingCorrection] DEBUG : FitExcluded()" << endl;
    #endif

    cout << "[PileupMaskingCorrection] INFO: rebin_factor = " << rebin_factor << ", exclude_dz = " << exclude_dz << endl;
    h_dz_rebinned = (TH1D*)h_dz->Clone();
    h_dz_rebinned->Rebin(rebin_factor); 
    TString hname = "h_dz_rebinned_";
    hname += tag;
    h_dz_rebinned->SetName(hname);

    reject = kTRUE;
    TString fname = "fit_gaussian_excluded_";
    fname += tag;
    /*//Temporarily changing back to always fitting with template
    f_dz_excluded = new TF1(fname, this, &PileupMaskingCorrection::fitFunc_generic_exclude, -1.*max_dz, max_dz, 5);
    //f_dz_excluded = new TF1(fname, fitFunc_generic_exclude, -1.*50., 50., 5);
    f_dz_excluded->SetNpx(10000);
    f_dz_excluded->SetParameter(0,h_dz_rebinned->Integral());
    f_dz_excluded->SetParLimits(0, 0., h_dz_rebinned->Integral()*100.);
    f_dz_excluded->FixParameter(1, -1.*exclude_dz);
    f_dz_excluded->FixParameter(2, exclude_dz);*/
    if ( (is_MC == 0) ){
      //cout << "[PileupMaskingCorrection] INFO: DATA" << endl;
      f_dz_excluded = new TF1(fname, this, &PileupMaskingCorrection::fitFunc_generic_exclude, -1.*max_dz, max_dz, 5);
      //f_dz_excluded = new TF1(fname, fitFunc_generic_exclude, -1.*50., 50., 5);
      f_dz_excluded->SetNpx(10000);
      f_dz_excluded->SetParameter(0,h_dz_rebinned->Integral());
      f_dz_excluded->SetParLimits(0, 0., h_dz_rebinned->Integral()*100.);
      f_dz_excluded->FixParameter(1, -1.*exclude_dz);
      f_dz_excluded->FixParameter(2, exclude_dz);
    } else {
      //cout << "[PileupMaskingCorrection] INFO: MC" << endl;
      f_dz_excluded = new TF1(fname, this, &PileupMaskingCorrection::fitFunc_gaussian, -1.*max_dz, max_dz, 3);
      f_dz_excluded->SetNpx(10000);
      f_dz_excluded->SetParameter(0,h_dz_rebinned->Integral());
      //f_dz_excluded->SetParLimits(0, -100.0*h_dz_rebinned->Integral()*100., h_dz_rebinned->Integral()*100.);
      f_dz_excluded->SetParameter(1, h_dz_rebinned->GetRMS());
      f_dz_excluded->FixParameter(2, 0.);
      f_dz_excluded->SetLineColor(kRed);
    }

    TString fit_options = "REM";
    if (do_ll) {
      fit_options += "L";
    }
    
    //Temporarily changing back to always fitting with template
    //h_dz_rebinned->Fit(f_dz_excluded,fit_options);
    if ( (is_MC == 0) ){
      //cout << "[PileupMaskingCorrection] INFO: DATA" << endl;
      h_dz_rebinned->Fit(f_dz_excluded,fit_options);
    } else {
      //cout << "[PileupMaskingCorrection] INFO: MC" << endl;
      h_dz_rebinned->Fit(f_dz_excluded,fit_options,"",-100,-20);
    }
    cout << "\t\t Excl. gaussian fit results: " << endl;
    cout << "\t\t par[0] = " << f_dz_excluded->GetParameter(0) << " +/- " << f_dz_excluded->GetParError(0) << endl;
    cout << "\t\t chi2/ndf = " << f_dz_excluded->GetChisquare() / f_dz_excluded->GetNDF() << endl;
    reject = kFALSE;
  } else {
    cout << "[PileupMaskingCorrection] INFO: else low_stats" << endl;
    h_dz_rebinned = (TH1D*)h_dz->Clone();
    h_dz_rebinned->Rebin(rebin_factor);
    TString hname = "h_dz_rebinned_";
    hname += tag;
    h_dz_rebinned->SetName(hname);

    TString fname = "fit_gaussian_excluded_";
    fname += tag;
    f_dz_excluded = new TF1(fname, this, &PileupMaskingCorrection::fitFunc_gaussian_exclude, -1.*max_dz, max_dz, 5);
    f_dz_excluded->SetNpx(10000);
    f_dz_excluded->SetParameter(0,0);
    f_dz_excluded->SetParameter(1,0.);
    f_dz_excluded->FixParameter(2, 0.);
    f_dz_excluded->FixParameter(3, -1.*exclude_dz);
    f_dz_excluded->FixParameter(4, exclude_dz);
  }
}

void PileupMaskingCorrection::FitExcludedResiduals() {
  if (!low_stats) {
    TString hname = "h_dz_rebinned_residuals_";
    hname += tag;
    h_dz_rebinned_residuals = new TH1D(hname, hname, h_dz_rebinned->GetNbinsX(), h_dz_rebinned->GetXaxis()->GetXmin(), h_dz_rebinned->GetXaxis()->GetXmax());
    for (Int_t i=1; i<=h_dz_rebinned_residuals->GetNbinsX(); i++) {
      Float_t current_dz = h_dz_rebinned_residuals->GetBinCenter(i);
      Float_t bin_error = h_dz_rebinned->GetBinError(i);

      if ((TMath::Abs(current_dz) < exclude_dz) || (bin_error <= 0.) || (TMath::Abs(current_dz) > max_dz)) {
        h_dz_rebinned_residuals->SetBinContent(i, 0.);
        continue;
      } else {
        h_dz_rebinned_residuals->SetBinContent(i, (f_dz_excluded->Eval(current_dz) - h_dz_rebinned->GetBinContent(i)) / bin_error);
      }
    }
  } else {
    TString hname = "h_dz_rebinned_residuals_";
    hname += tag;
    h_dz_rebinned_residuals = new TH1D(hname, hname, h_dz_rebinned->GetNbinsX(), h_dz_rebinned->GetXaxis()->GetXmin(), h_dz_rebinned->GetXaxis()->GetXmax());
  }
}

void PileupMaskingCorrection::MakeFullGaussian() {
	#ifdef DEBUG_PMC
  cout << "[PileupMaskingCorrection] DEBUG : MakeFullGaussian()" << endl;
	#endif

  TString fname = "fit_gaussian_full_";
  fname += tag;
  /*//Temporarily changing back to always fitting with template
  f_dz_full = new TF1(fname, this, &PileupMaskingCorrection::fitFunc_generic, h_dz->GetXaxis()->GetXmin(), h_dz->GetXaxis()->GetXmax(), 1);
  f_dz_full->SetParameter(0, f_dz_excluded->GetParameter(0) / rebin_factor);*/
  if ( (is_MC == 0) ){
    //cout << "[PileupMaskingCorrection] INFO: DATA" << endl;
    f_dz_full = new TF1(fname, this, &PileupMaskingCorrection::fitFunc_generic, h_dz->GetXaxis()->GetXmin(), h_dz->GetXaxis()->GetXmax(), 1);
    f_dz_full->SetParameter(0, f_dz_excluded->GetParameter(0) / rebin_factor);
  } else {
    //cout << "[PileupMaskingCorrection] INFO: MC" << endl;
    f_dz_full = new TF1(fname, this, &PileupMaskingCorrection::fitFunc_gaussian, h_dz->GetXaxis()->GetXmin(), h_dz->GetXaxis()->GetXmax(), 3);
    f_dz_full->SetParameter(0, f_dz_excluded->GetParameter(0) / rebin_factor);
    f_dz_full->SetParameter(1, f_dz_excluded->GetParameter(1));
    f_dz_full->SetParameter(2, f_dz_excluded->GetParameter(2));

    //f_dz_full = f_dz_excluded;
    //f_dz_full->SetName(fname);
    h_dz_random = new TH1D( h_dz_rebinned->GetName(), h_dz_rebinned->GetName(), h_dz_rebinned->GetXaxis()->GetNbins(), h_dz_rebinned->GetXaxis()->GetXmin(), h_dz_rebinned->GetXaxis()->GetXmax());
    h_dz_random->FillRandom(fname,1000000);
    h_dz_random->Sumw2();
    cout << "h_dz_random->Integral() " << h_dz_random->Integral() << endl;
  }
}

void PileupMaskingCorrection::MakeDifferentialPmask() {
	#ifdef DEBUG_PMC
  cout << "[PileupMaskingCorrection] DEBUG : MakeDifferentialPmask()" << endl;
	#endif

  if (!low_stats) {
    h_pmask_dz = (TH1D*)h_dz->Clone();
    cout << "h_pmask_dz->GetNbinsX(): " << h_pmask_dz->GetNbinsX() << endl;
    TString hname = "h_pmask_dz_";
    hname += tag;
    h_pmask_dz->SetName(hname);
    for (Int_t i=1; i<=h_dz->GetNbinsX(); i++) {
      Double_t data = h_dz->GetBinContent(i);
      Double_t data_err = h_dz->GetBinError(i);
      Double_t dz = h_dz->GetBinCenter(i);
      Double_t fit = f_dz_full->Eval(dz);
      Double_t current_pmask;
      if (fit > 0.) {
        current_pmask = (fit - data) / fit;
      } else {
        current_pmask = 0.;
      }
      h_pmask_dz->SetBinContent(i, current_pmask);
      //cout << i << " dz: " << dz << ", fit: " << fit <<", data: " << data << ", current_pmask: " << current_pmask << endl;

      //Errors?!
      Double_t err;
      if ((fit > 0) && (TMath::Abs(current_pmask) < 1.)) {
        err = TMath::Sqrt(TMath::Abs(current_pmask) * (1. - TMath::Abs(current_pmask)) / fit);
        //cout << "current_pmask = (" << fit << " - " << data << ") / " << fit << " = " << current_pmask << ", err = " << err << endl;
      } else {
        err = 0.;
      }
      //Hack: for some reason, the above sometimes sets the error to NaN and screws everything up later :(
      if (err > 0) {
        h_pmask_dz->SetBinError(i, err);
      } else {
        h_pmask_dz->SetBinError(i, 0.);
      }
    }

    // Some post-processing: if p_mask = 0, that means there just weren't enough statistics for points in that bin. So, umm, replace with the average of the nearest nonzero points...
#ifdef MC
    cout << "Ever here?" << endl;
    TH1F *h_pmask_dz_copy = (TH1F*)h_pmask_dz->Clone("h_pmask_dz_copy");
    for (Int_t bin = 1; bin <= h_pmask_dz->GetXaxis()->GetNbins(); bin++) {
      if (h_pmask_dz->GetBinContent(bin) == 0) {
        Int_t bin_minus = bin;
        Int_t bin_plus = bin;
        Float_t pmask_minus = 0.;
        Float_t pmask_plus = 0.;
        Float_t pmask_err_minus = 0.;
        Float_t pmask_err_plus = 0.;
        while (pmask_minus == 0.) {
          bin_minus--;
          if (bin_minus < 1) {
            break;
          }
          if (h_pmask_dz_copy->GetBinContent(bin_minus) != 0.) {
            pmask_minus = h_pmask_dz_copy->GetBinContent(bin_minus);
            pmask_err_minus = h_pmask_dz_copy->GetBinError(bin_minus);
          }
        }
        while (pmask_plus == 0.) {
          bin_plus++;
          if (bin_plus > h_pmask_dz->GetXaxis()->GetNbins()) {
            break;
          }
          if (h_pmask_dz_copy->GetBinContent(bin_plus) != 0.) {
            pmask_plus = h_pmask_dz_copy->GetBinContent(bin_plus);
            pmask_err_plus = h_pmask_dz_copy->GetBinError(bin_plus);
          }
        }
        if ((pmask_minus != 0.) && (pmask_plus != 0.)) {
          h_pmask_dz->SetBinContent(bin, 0.5 * (pmask_minus + pmask_plus));
          h_pmask_dz->SetBinError(bin, 0.5 * TMath::Sqrt(pmask_err_minus*pmask_err_minus + pmask_err_plus*pmask_err_plus));
        }
      }
    }
#endif
  } else {
    TString hname = "h_pmask_dz_";
    hname += tag;
    h_pmask_dz =  new TH1D(hname, hname, h_dz->GetNbinsX(), h_dz->GetXaxis()->GetXmin(), h_dz->GetXaxis()->GetXmax());
    h_pmask_dz->SetName(hname);
  }
}

/***************************************************************************************/
/// - Method 2: load cached p_mask vs. dz distribution (and everything else)
/***************************************************************************************/

PileupMaskingCorrection::PileupMaskingCorrection(TString p_tag, Int_t p_ntrkcut) {

  cout << "[PileupMaskingCorrection] INFO : Initializing PileupMaskingCorrection (load cache method)" << endl;

  redo_cache = false;
  finished = false;
  pmask_available = true;

  std::vector<TString> mc_samples;
  mc_samples.push_back("mc_7TeV_16.X_normal_pythia6_pu");
  mc_samples.push_back("mc_7TeV_16.X_hybrid_pythia6_pu");
  mc_samples.push_back("mc_7TeV_17.2_normal_pythia6_pu");
  mc_samples.push_back("mc_7TeV_17.2_VtxLumi_pythia6_pu");
  mc_samples.push_back("mc_7TeV_17.2_normal_pythia8_pu");
  mc_samples.push_back("mc_7TeV_17.2_VtxLumi_pythia8_pu");
  mc_samples.push_back("mc_7TeV_17.2_normal_pythia6_pu_noslim");
  mc_samples.push_back("mc_7TeV_17.2_normal_pythia8DL_pu");
  mc_samples.push_back("mc_7TeV_17.2_normal_pythia8DL_pu_noslim");
  mc_samples.push_back("mc_8TeV_17.2_normal_pythia8_pu");
  mc_samples.push_back("mc_8TeV_17.2_VtxLumi_pythia8_pu");
  mc_samples.push_back("mc_8TeV_17.2_VtxLumi_2newsets");
  mc_samples.push_back("mc_8TeV_17.2_VtxLumi_2newsets_10Dec");
  mc_samples.push_back("mc_8TeV_17.2_VtxLumi_mumax20");
  mc_samples.push_back("mc_8TeV_17.2_VtxLumi_mumax75");
  mc_samples.push_back("mc_7TeV_17.2_normal_pythia8_pu_bs45");
  mc_samples.push_back("mc_7TeV_17.2_VtxLumi_pythia8_pu_bs45");
  mc_samples.push_back("mc_7TeV_17.2_normal_pythia8_pu_bs55");

  std::vector<TString> data_samples;
  data_samples.push_back("data_7TeV_17.2-normal");
  data_samples.push_back("data_8TeV_17.2-normal");
  data_samples.push_back("data_7TeV_17.2-VtxLumi");
  data_samples.push_back("data_8TeV_17.2-VtxLumi");
  data_samples.push_back("data_8TeV_17.2-VtxLumi_201351");
  data_samples.push_back("data_8TeV_17.2-VtxLumi_201351_26Nov");
  data_samples.push_back("data_8TeV_17.2-VtxLumi_207216");
  data_samples.push_back("data_8TeV_17.2-VtxLumi_207219");
  data_samples.push_back("data_8TeV_17.2-VtxLumi_214984");
  data_samples.push_back("data_8TeV_17.2-VtxLumi_215021");

  if ((find(mc_samples.begin(), mc_samples.end(), p_tag) == mc_samples.end()) && (find(data_samples.begin(), data_samples.end(), p_tag) == data_samples.end())) {
    cerr << "[PileupMaskingCorrection] ERROR : Requested PMC tag " << p_tag << " not recognized. Exiting..." << endl;
    exit(1);
  }

  TString sample;
  if (p_tag == "data_7TeV_17.2-normal") {
    sample = "data_7TeV_17.2_normal";
  } else if (p_tag == "data_7TeV_17.2-VtxLumi") {
    sample = "data_7TeV_17.2_VtxLumi";
  } else if (p_tag == "data_8TeV_17.2-normal") {
    sample = "data_8TeV_17.2_normal";
  } else if (p_tag == "data_8TeV_17.2-VtxLumi") {
    sample = "data_8TeV_17.2_VtxLumi";
  } else if (p_tag == "data_8TeV_17.2-VtxLumi_201351") {
    sample = "data_8TeV_17.2_VtxLumi_201351";
  } else if (p_tag == "data_8TeV_17.2-VtxLumi_201351_26Nov") {
    sample = "data_8TeV_17.2_VtxLumi_201351_26Nov";
  } else if (p_tag == "data_8TeV_17.2-VtxLumi_207216") {
    sample = "data_8TeV_17.2_VtxLumi_207216";
  } else if (p_tag == "data_8TeV_17.2-VtxLumi_207219") {
    sample = "data_8TeV_17.2_VtxLumi_207219";  
  } else if (p_tag == "data_8TeV_17.2-VtxLumi_214984") {
    sample = "data_8TeV_17.2_VtxLumi_214984";
  } else if (p_tag == "data_8TeV_17.2-VtxLumi_215021") {
    sample = "data_8TeV_17.2_VtxLumi_215021";
  } else {
    sample = p_tag;
  }

  TString path = GlobalSettings::path_maskingCorrection;
  path += "/";
  path += sample;
  //path += "data_8TeV_17.2_VtxLumi_207216"; //trying to use July's pmask vs dz histogram for April
  path += "/pmask_cache.root";

  cout << "[PileupMaskingCorrection] INFO : Loading cache from " << path << endl;

  TFile *f_in = new TFile(path, "READ");
  TString pmask_dz_name = "h_pmask_dz_NTrk";
  pmask_dz_name += p_ntrkcut;
  h_pmask_dz = (TH1D*)f_in->Get(pmask_dz_name);
  cout << "[PileupMaskingCorrection] INFO : h_pmask_dz->GetNbinsX() " << h_pmask_dz->GetNbinsX() << endl; 
  if (h_pmask_dz) {
    h_pmask_dz->SetDirectory(0);
  } else {
    cerr << "[PileupMaskingCorrection] ERROR : Could not find histogram " << pmask_dz_name << " in file " << path << ". Exiting..." << endl;
    exit(1);
  }
  f_in->Close();

  pmask_scale = 1.;
  exclude_dz = 20.;

  /* Everything else is not really relevant if loading from cache, but set anyways to be safe */
  rebin_factor = 10; 
  max_dz = 300.;
  do_ll = false;
  low_stats = false;

}

#ifdef REMOVED_061412
PileupMaskingCorrection::PileupMaskingCorrection(TString p_source, TString p_energy, TString p_settings, TString p_ntrkcut) {

  cout << "[PileupMaskingCorrection] INFO : Initializing PileupMaskingCorrection (load cache method)" << endl;

  redo_cache = false;
  finished = false;
  pmask_available = true;

  //TString path = "/u/dryu/Luminosity/PileupMaskingCorrection/cache/";
  TString path = GlobalSettings::path_maskingCorrection;
  if (p_energy == "7" && p_source == "data") {
    path += "data_7TeV/";
  } else if (p_energy == "8" && p_source == "data") {
    path += "data_8TeV/";
  } else if (p_energy == "7" && p_source == "mc") {
    path += "mc_7TeV/";
  } else if (p_energy == "8" && p_source == "mc") {
    path += "mc_8TeV/";
  } else {
    cerr << "[PileupMaskingCorrection] ERROR : Unknown sample requested: energy = " << p_energy << " and source = " << p_source << ". Exiting..." << endl;
    exit(1);
  }

  path += p_settings;
  path += "/pmask_cache.root";

#ifdef DEBUG_PMC
  cout << "Loading cache from " << path << endl;
#endif

  TFile *f_in = new TFile(path, "READ");
  TString pmask_dz_name = "h_pmask_dz_NTrk";
  pmask_dz_name += p_ntrkcut;
  h_pmask_dz = (TH1D*)f_in->Get(pmask_dz_name);
  if (h_pmask_dz) {
    h_pmask_dz->SetDirectory(0);
  } else {
    cerr << "[PileupMaskingCorrection] ERROR : Could not find histogram " << pmask_dz_name << " in file " << path << ". Exiting..." << endl;
    exit(1);
  }
  f_in->Close();

  pmask_scale = 1.;

  /* Everything else is not really relevant if loading from cache, but set anyways to be safe */
  rebin_factor = 8;
  exclude_dz = 20.;
  max_dz = 300.;
  do_ll = false;
  low_stats = false;

}
#endif

#ifdef OLD
PileupMaskingCorrection::PileupMaskingCorrection(TString f_cache, TString pmask_name, TString string1) {
  cout << "[PileupMaskingCorrection] INFO : Initializing PileupMaskingCorrection (load cache method)" << endl;
  tag = string1;

  redo_cache = false;
  finished = false;
  pmask_available = true;

#ifdef MC
  TString key_suffix = "_mc_NTrk5";
#else
  TString key_suffix = "_BCID81_NTrkCut5";
#endif

#ifdef DEBUG_PMC
  cout << "Loading h_dz from " << f_cache << " : " << TString("h_dz") + key_suffix << endl;
#endif

  TFile *f_in = new TFile(f_cache, "READ");
  h_dz = (TH1D*)f_in->Get(TString("h_dz") + key_suffix);
  h_dz->SetDirectory(0);
  h_dz_rebinned = (TH1D*)f_in->Get(TString("h_dz_rebinned") + key_suffix);
  h_dz_rebinned->SetDirectory(0);
  f_dz_excluded = (TF1*)f_in->Get(TString("f_dz_excluded") + key_suffix);
  f_dz_full = (TF1*)f_in->Get(TString("f_dz_full") + key_suffix);
  h_pmask_dz = (TH1D*)f_in->Get(pmask_name);
  if (h_pmask_dz) {
    h_pmask_dz->SetDirectory(0);
  }
  f_in->Close();

  /* Everything else is not really relevant if loading from cache, but set anyways to be safe */
  rebin_factor = 8;
  exclude_dz = 20.;
  max_dz = 300.;
  do_ll = false;
  low_stats = false;
}
#endif

/***************************************************************************************/
/// - Generate the correction TGraphs
/***************************************************************************************/

#ifdef REMOVED_061412
void PileupMaskingCorrection::GenerateCorrection(Float_t sigma_z) {

  cout << "[PileupMaskingCorrection] INFO : Generating correction (gaussian sigma method)." << endl;

  CalculateTotalPmask(sigma_z);
  MakePuCorrTGraphs();

  finished = true;

}

#endif

void PileupMaskingCorrection::GenerateCorrection(TH1D *h_dz_input) {

  cout << "[PileupMaskingCorrection] INFO : Generating correction (input dz template method)." << endl;
  cout << "h_dz_input->GetXaxis()->GetNbins() = " << h_dz_input->GetXaxis()->GetNbins() << endl;
  cout << "h_dz_input->Integral() = " << h_dz_input->Integral() << endl;

  if (finished) {
    cout << "PileupMaskingCorrection : [WARNING] Masking correction TGraphs already exist! Attempt to generate new masking correction TGraphs ignored." << endl;
    return;
  }
  cout << "[PileupMaskingCorrection] INFO : CalculateTotalPmask " << endl;
  CalculateTotalPmask(h_dz_input);
  cout << "[PileupMaskingCorrection] INFO : MakePuCorrTGraphs() " << endl;
  MakePuCorrTGraphs();

  finished = true;

}

void PileupMaskingCorrection::LoadCorrection(TGraphErrors *tg1, TGraphErrors *tg2) {

  if (finished) {
    cout << "PileupMaskingCorrection : [WARNING] Masking correction TGraphs already exist! Attempt to load masking correction TGraphs ignored." << endl;
    return;
  }

  tg_mu_obs_vs_mu_actual = (TGraphErrors*)tg1->Clone();
  tg_pileup_correction = (TGraphErrors*)tg2->Clone();
  finished = true;

}

#ifdef REMOVED_061412
void PileupMaskingCorrection::CalculateTotalPmask(Float_t sigma_z) {
  cout << "[PileupMaskingCorrection] WARNING : Attempting to use deprecated function CalculateTotalPmask(Float_t sigma_z). Expect bugs and crashes!" << endl;
#ifdef DEBUG_PMC
  cout << "[PileupMaskingCorrection] DEBUG : CalculateTotalPmask(sigma_z = " << sigma_z << ")" << endl;
#endif
  tg_pmask_cumulative = new TGraph(h_pmask_dz->GetNbinsX());// ??
  TString tgname = "tg_pmask_cumulative_";
  tgname += tag;
  tg_pmask_cumulative->SetName(tgname);

  TString fname = "fit_gaussian_full_normalized_";
  fname += tag;
  f_dz_full_normalized = new TF1(fname, this, &PileupMaskingCorrection::fitFunc_gaussian, h_dz->GetXaxis()->GetXmin(), h_dz->GetXaxis()->GetXmax(), 3);
  f_dz_full_normalized->SetParameter(0, 1.);
  //f_dz_full_normalized->SetParameter(1, f_dz_full->GetParameter(1)); // To compare, you really should use a fixed dz for all curves.
  f_dz_full_normalized->SetParameter(1, sigma_z);
  f_dz_full_normalized->SetParameter(2, 0.);

  if (low_stats) {
    cout << "[PileupMaskingCorrection] Line 570, if low_stats" << endl;
    pmask = 0.;
    pmask_err = 0.;
  } else {
    pmask = 0;
    Double_t pmask_err2 = 0.;
    Double_t dpmask_dp0 = 0.;
    Double_t dpmask_dp1 = 0.;
    Double_t p0 = f_dz_full->GetParameter(0);
    Double_t p1 = f_dz_full->GetParameter(1);

    for (Int_t i=1; i<=h_pmask_dz->GetNbinsX(); i++) {
      Double_t dz = h_pmask_dz->GetBinCenter(i);
      if (TMath::Abs(dz) > exclude_dz) {
        continue;
      }
      Double_t bin_width =h_pmask_dz->GetBinWidth(i);
      Double_t current_pmask_dz = h_pmask_dz->GetBinContent(i);
      pmask += bin_width * f_dz_full_normalized->Eval(dz) * current_pmask_dz;

      //Errors, again, are hard! Below, they're computed from scratch, since starting from bin errors on pmask_dz ignores the fact that we're using the same gaussian in multiple places.
      Double_t data = h_dz->GetBinContent(i);
      Double_t fit = f_dz_full->Eval(dz);
      Double_t data_err;
      if ((current_pmask_dz > 0) && (current_pmask_dz < 1)) {
        data_err = TMath::Sqrt(fit * current_pmask_dz * (1.-current_pmask_dz));
      } else {
        data_err = 0.;
      }

      //pmask_err2 += TMath::Power(bin_width / p0 * data_err, 2);
      if ((current_pmask_dz > 0) && (current_pmask_dz < 1)) {
        pmask_err2 += TMath::Power(bin_width / p0 * TMath::Sqrt(data * current_pmask_dz), 2);
      } else {
        pmask_err2 += 0.;
      }
      dpmask_dp0 += data / (p0*p0) * bin_width;
      dpmask_dp1 += bin_width/TMath::Sqrt(2*TMath::Pi()) * TMath::Exp(-1.*dz*dz/(2*p1*p1)) * ((TMath::Power(dz, 2) / TMath::Power(p1, 4)) - (1. / TMath::Power(p1, 2)));

      //For in-depth analysis, store some of the partial sums.
      Int_t first_point = h_pmask_dz->FindBin(-1. * exclude_dz);
      if (i%rebin_factor == 0) {
        Int_t current_point = (i-first_point)/rebin_factor;
        tg_pmask_cumulative->SetPoint(current_point, dz, pmask);
      }
    }
    //Complete the error calculation: add in fit errors, then sqrt.
    Double_t p0_err = f_dz_full->GetParError(0);
    Double_t p1_err = f_dz_full->GetParError(1);
    pmask_err2 += TMath::Power(dpmask_dp0 * p0_err, 2) + TMath::Power(dpmask_dp1 * p1_err, 2);
    pmask_err = TMath::Sqrt(pmask_err2);
  }
}

#endif

void PileupMaskingCorrection::CalculateTotalPmask(TH1D *h_dz_input) {
  #ifdef DEBUG_PMC
  cout << "[PileupMaskingCorrection] DEBUG : Calculating p_mask using input dz template." << endl;
  #endif

  // Normalize the expected delta Z distribution
  h_dz_input = (TH1D*)h_dz_input->Clone();
  h_dz_input->Scale(1. / h_dz_input->Integral());

  // h_dz_input and h_pmask_dz have to have equal bin width. If they differ, try to fix it by rebinning. However, it's best if they are equal to begin with!!!
  if (h_dz_input->GetNbinsX() != h_pmask_dz->GetNbinsX()) {
    cout << "[PileupMaskingCorrection] WARNING : h_dz_input (" << h_dz_input->GetNbinsX() << ") and h_pmask_dz (" << h_pmask_dz->GetNbinsX() << ") have different number of bins. Trying to fix..." << endl;

    if (h_dz_input->GetNbinsX() > h_pmask_dz->GetNbinsX()) {
      if (h_dz_input->GetNbinsX() % h_pmask_dz->GetNbinsX() == 0) {
        cout << "[PileupMaskingCorrection] WARNING : Rebinning h_dz_input by " << h_dz_input->GetNbinsX() / h_pmask_dz->GetNbinsX() << endl;
        h_dz_input->Rebin(h_dz_input->GetNbinsX() / h_pmask_dz->GetNbinsX());
        cout << "[PUMC] Line 590 h_dz_input->GetNbinsX() / h_pmask_dz->GetNbinsX() = " << h_dz_input->GetNbinsX() / h_pmask_dz->GetNbinsX() << endl;
      } else {
        cerr << "Rebin failed. Exiting..." << endl;
        exit(1);
      }
    }
    if (h_dz_input->GetNbinsX() < h_pmask_dz->GetNbinsX()) {
      if (h_pmask_dz->GetNbinsX() % h_dz_input->GetNbinsX() == 0) {
        #ifdef DEBUG_PMC
        cout << "[PileupMaskingCorrection] DEBUG : Before rebinning, pmask_dz(0) = " << h_pmask_dz->GetBinContent(h_pmask_dz->FindBin(0.)) << endl;
        #endif
        cout << "[PileupMaskingCorrection] WARNING : Rebinning h_pmask_dz by " << h_pmask_dz->GetNbinsX() / h_dz_input->GetNbinsX() << endl;
        Int_t rebin_factor = h_pmask_dz->GetNbinsX() / h_dz_input->GetNbinsX();
        h_pmask_dz->Rebin(rebin_factor);
        h_pmask_dz->Scale(1. / (float)rebin_factor);
        #ifdef DEBUG_PMC
        cout << "[PileupMaskingCorrection] DEBUG : After rebinning, pmask_dz(0) = " << h_pmask_dz->GetBinContent(h_pmask_dz->FindBin(0.)) << endl;
        #endif
      } else {
        cerr << "Rebin failed. Exiting..." << endl;
        exit(1);
      }
    }
  }

  // Check if the rebinning succeeded
  if (h_dz_input->GetBinWidth(1) != h_pmask_dz->GetBinWidth(1)) {
    cerr << "[PileupMaskingCorrection] ERROR : Histogram bin widths are different! Go back to AnaVtxTree.cxx and fix it!" << endl;
    cerr << "[PileupMaskingCorrection] ERROR : h_dz_input->GetBinWidth(1) = " << h_dz_input->GetBinWidth(1) << endl;
    cerr << "[PileupMaskingCorrection] ERROR : h_pmask_dz->GetBinWidth(1) = " << h_pmask_dz->GetBinWidth(1) << endl;
    exit(1);
  }

  // Debug plot: partial integrals of p_mask(dz) * f(dz).
  tg_pmask_cumulative = new TGraph(h_pmask_dz->GetNbinsX());
  TString tgname = "tg_pmask_cumulative_";
  tgname += tag;
  tg_pmask_cumulative->SetName(tgname);

  // Compute Int[h_pmask_dz * h_dz_input, {z, -exclude_dz, exclude_dz}]
  if (low_stats) {
    cout << "[PileupMaskingCorrection] WARNING : Low stats, setting p_mask = 0." << endl;
    pmask = 0.;
    pmask_err = 0.;
  } else {
    pmask = 0;
    Double_t pmask_err2 = 0.;

    for (Int_t i=1; i<=h_pmask_dz->GetNbinsX(); i++) {
      Double_t dz = h_pmask_dz->GetBinCenter(i);
      if (TMath::Abs(dz) > exclude_dz) {
        continue;
      }
      //      Double_t bin_width =h_pmask_dz->GetBinWidth(i);
      Double_t current_pmask_dz = h_pmask_dz->GetBinContent(i);
      //Double_t current_p_dz = h_dz_input->GetBinContent(h_dz_input->FindBin(dz));
      Double_t current_p_dz = h_dz_input->GetBinContent(i);
      pmask += current_p_dz * current_pmask_dz;

      //Errors?!
      Double_t current_pmask_dz_err = h_pmask_dz->GetBinError(i);
      Double_t current_p_dz_err = h_dz_input->GetBinError(h_dz_input->FindBin(dz));
      Double_t current_err;
      if ((current_pmask_dz != 0) && (current_p_dz != 0)) {
        current_err = current_p_dz * current_pmask_dz * TMath::Sqrt(TMath::Power(current_pmask_dz_err / current_pmask_dz, 2) + TMath::Power(current_p_dz_err / current_p_dz, 2));
      } else {
        current_err = 0;
      }
      pmask_err2 += TMath::Power(current_err, 2);

      //For in-depth analysis, store some of the partial sums.
      Int_t first_point = h_pmask_dz->FindBin(-1. * exclude_dz);
      if (i%rebin_factor == 0) {
        Int_t current_point = (i-first_point)/rebin_factor;
        tg_pmask_cumulative->SetPoint(current_point, dz, pmask);
      }
    }
    //Complete the error calculation: add in fit errors, then sqrt.
    pmask_err = TMath::Sqrt(pmask_err2);
  }

  pmask = pmask * pmask_scale;
  pmask_err = pmask_err * pmask_scale;

  cout << "[PileupMaskingCorrection] INFO : p_mask = " << pmask << " +/- " << pmask_err << endl;
}

void PileupMaskingCorrection::MakePuCorrTGraphs() {
	#ifdef DEBUG_PMC
  cout << "[PileupMaskingCorrection] DEBUG : MakePuCorrTGraphs()" << endl;
	#endif

  Double_t n_points = 400;
  Double_t mumax = 0.;

  if ( (is_MC == 0) ){
    //cout << "[PileupMaskingCorrection] INFO: DATA" << endl;
    mumax = 40; 
  } else{
    //cout << "[PileupMaskingCorrection] INFO: MC" << endl;
    mumax = 71; // Maximum value of ei_actualIntPerXing for the high mu sample
  }

  cout << "[PileupMaskingCorrection] DEBUG : mumax = " << mumax << endl;

  TString tgname = "tg_mu_obs_vs_mu_actual_";
  tgname += tag;
  tg_mu_obs_vs_mu_actual = new TGraphErrors(n_points);
  tg_mu_obs_vs_mu_actual->SetName(tgname);

  float current_mu_actual = 0;

  for (int i = 1; i <= n_points; i++) {
    if (i % 100 == 0) {
      cout << "Making mu correction map: " << TMath::Nint(100. * i / n_points) << "%\r" << std::flush;
    }
      current_mu_actual = mumax * i / n_points;
    if (!low_stats) {
      tg_mu_obs_vs_mu_actual->SetPoint(i, current_mu_actual, meanMuObs(current_mu_actual, pmask));
      //      cout << "current_mu_actual = " <<current_mu_actual<< ", pmask = " << pmask << ", meanMuObs = " << meanMuObs(current_mu_actual, pmask) << endl;
      tg_mu_obs_vs_mu_actual->SetPointError(i, 0., meanMuObs(current_mu_actual, pmask + pmask_err) - meanMuObs(current_mu_actual, pmask));
    }
    else {cout << "[PileupMaskingCorrection] Line 753, else low_stats" << endl;}
  }
  #ifdef DEBUG_PMC
  cout << "[PileupMaskingCorrection] DEBUG : tg_mu_obs_vs_mu_actual " << endl;
  cout << "Number of points = " << tg_mu_obs_vs_mu_actual->GetN() << endl;
  cout << "X axis, min = " << tg_mu_obs_vs_mu_actual->GetXaxis()->GetXmin() << ", max = " << tg_mu_obs_vs_mu_actual->GetXaxis()->GetXmax() << endl;
  cout << "Y axis, min = " << tg_mu_obs_vs_mu_actual->GetYaxis()->GetXmin() << ", max = " << tg_mu_obs_vs_mu_actual->GetYaxis()->GetXmax() << endl;
   
	#endif
	
  tg_pileup_correction = new TGraphErrors(n_points);
  tgname = "tg_pileup_correction_";
  tgname += tag;
  tg_pileup_correction->SetName(tgname);


  for (int i=0; i<n_points; i++) {
    if (!low_stats) {
      Double_t current_mu_obs = tg_mu_obs_vs_mu_actual->GetY()[i];
      Double_t current_mu_actual = tg_mu_obs_vs_mu_actual->GetX()[i];
      Double_t current_mu_obs_err = tg_mu_obs_vs_mu_actual->GetEY()[i];
      if (current_mu_obs > 0) {
        tg_pileup_correction->SetPoint(i, current_mu_obs, current_mu_actual / current_mu_obs);
        tg_pileup_correction->SetPointError(i, 0., ( current_mu_obs_err * current_mu_actual / ( current_mu_obs * current_mu_obs ) ) ); 
      } else {
        tg_pileup_correction->SetPoint(i, current_mu_obs, 0.);
        tg_pileup_correction->SetPointError(i, 0., 0.);
      }
    } else {
      cout << "[PileupMaskingCorrection] Line " << __LINE__ << ", else low_stats" << endl;
      Double_t current_mu_obs = tg_mu_obs_vs_mu_actual->GetY()[i];
      tg_pileup_correction->SetPoint(i, current_mu_obs, current_mu_obs);
      tg_pileup_correction->SetPointError(i, 0., 0.);
    }
  }
  #ifdef DEBUG_PMC
  cout << "[PileupMaskingCorrection] DEBUG : tg_pileup_correction " << endl;
  cout << "Number of points = " << tg_pileup_correction->GetN() << endl;
  cout << "X axis, min = " << tg_pileup_correction->GetXaxis()->GetXmin() << ", max = " << tg_pileup_correction->GetXaxis()->GetXmax() << endl;
  cout << "Y axis, min = " << tg_pileup_correction->GetYaxis()->GetXmin() << ", max = " << tg_pileup_correction->GetYaxis()->GetXmax() << endl;
	#endif
}

void PileupMaskingCorrection::Save(TString path_rootfile, TString suffix, bool new_file) {
  if (!finished) {
    cerr << "Run the calculations before trying to save!" << endl;
    exit(1);
  }

  if (!redo_cache) {
    return;
  }

  TFile *f_out;
  if (new_file) {
    f_out = new TFile(path_rootfile, "RECREATE");
  } else {
    f_out = new TFile(path_rootfile, "UPDATE");
  }

  h_dz->Write(TString("h_dz_") + suffix);
  h_dz_rebinned->Write(TString("h_dz_rebinned_") + suffix);
  f_dz_excluded->Write(TString("f_dz_excluded_") + suffix);
  f_dz_full->Write(TString("f_dz_full_") + suffix);
  h_pmask_dz->Write(TString("h_pmask_dz_") + suffix);
  tg_mu_obs_vs_mu_actual->Write(TString("tg_mu_obs_vs_mu_actual_") + suffix);
  tg_pileup_correction->Write(TString("tg_pileup_correction_") + suffix);
  f_out->Close();

	#ifdef SAVE_RESIDUALS
  TString path_quick = path_rootfile; path_quick+="residuals.root";
  TFile *f_out_quick;
  if (new_file) {
    f_out_quick = new TFile(path_quick,"RECREATE");
  } else {
    f_out_quick = new TFile(path_quick, "UPDATE");
  }
  h_dz_rebinned_residuals->Write();
  f_out_quick->Close();
	#endif
}

void PileupMaskingCorrection::SaveDebugHistograms(TString path_rootfile, TString p_suffix, bool new_file) {

  if (!finished) {
    cerr << "Run the calculations before trying to save!" << endl;
    exit(1);
  }

  TFile *f_debug;
  if (new_file) {
    f_debug = new TFile(path_rootfile, "RECREATE");
  } else {
    f_debug = new TFile(path_rootfile, "UPDATE");
  }
  if (h_dz != 0x0) {
    h_dz->Write(h_dz->GetName() + p_suffix);
  }
  if (h_dz_rebinned != 0x0) {
    h_dz_rebinned->Write(h_dz_rebinned->GetName() + p_suffix);
  }
  if (f_dz_excluded != 0x0) {
    f_dz_excluded->Write(f_dz_excluded->GetName() + p_suffix);
  }
  if (h_pmask_dz != 0x0) {
    h_pmask_dz->Write(h_pmask_dz->GetName() + p_suffix);
  }
  if (tg_pmask_cumulative != 0x0) {
    tg_pmask_cumulative->Write(tg_pmask_cumulative->GetName() + p_suffix);
  }
  if (tg_mu_obs_vs_mu_actual != 0x0) {
    tg_mu_obs_vs_mu_actual->Write(tg_mu_obs_vs_mu_actual->GetName() + p_suffix);
  }
  if (tg_pileup_correction != 0x0) {
    tg_pileup_correction->Write(tg_pileup_correction->GetName() + p_suffix);
  }

  f_debug->Close();

}


/***************************************************************************************/
/// - Main function!
/***************************************************************************************/

Double_t PileupMaskingCorrection::GetCorrectionFactor(Double_t mu_in) {
  // Correction is not valid below the lowest point of the TGraph...
  if (mu_in < tg_pileup_correction->GetX()[0]) {
    return 1.;
  }

  return tg_pileup_correction->Eval(mu_in);
}

Double_t PileupMaskingCorrection::GetCorrectionFactorHigh(Double_t mu_in, Double_t mu_in_err) {
  if (mu_in < tg_pileup_correction->GetX()[0]) {
    return 1.;
  }

  return tg_pileup_correction->Eval(mu_in + mu_in_err);
}

Double_t PileupMaskingCorrection::GetCorrectionFactorLow(Double_t mu_in, Double_t mu_in_err) {
  if (mu_in < tg_pileup_correction->GetX()[0]) {
    return 1.;
  }

  return tg_pileup_correction->Eval(mu_in - mu_in_err);
}

/***************************************************************************************/
/// - Lastly, debugging methods
/***************************************************************************************/

#ifdef FITWITHGAUSSIAN
TH1D* PileupMaskingCorrection::GetExpectedDzDistribution() {
  cout << "GetExpectedDzDistribution: Returning h_dz_random" << endl;
  if (!h_dz_random) {
    cout << "WARNING: h_dz_random requested, but does not exist. Returning null." << endl;
    return 0;
  }
  cout << "h_dz_random->Integral() " << h_dz_random->Integral() << endl;
  return h_dz_random;
}

#else

TH1D* PileupMaskingCorrection::GetExpectedDzDistribution() {
  cout << "GetExpectedDzDistribution" << endl;
  if (!h_dz_expected) {
    cout << "WARNING: h_dz_expected requested, but does not exist. Returning null." << endl;
    return 0;
  }

  return h_dz_expected;
}

#endif

TF1* PileupMaskingCorrection::GetExcludedDzFitFunc() {
  if (!finished) {
    cerr << "Run the calculations before trying to save!" << endl;
    exit(1);
  }

  return f_dz_excluded;
}

Double_t PileupMaskingCorrection::GetGaussianSigma() {
  if (!finished) {
    cerr << "Run the calculations before trying to save!" << endl;
    exit(1);
  }

  return f_dz_full->GetParameter(1);
}


TF1* PileupMaskingCorrection::GetFullDzFitFunc() {
  if (!finished) {
    cerr << "Run the calculations before calling GetFullDzFitFunc()!" << endl;
    exit(1);
  }

  return f_dz_full;
}

TH1D* PileupMaskingCorrection::GetDifferentialPmask() {
  if (!h_pmask_dz) {
    cerr << "h_pmask_dz doesn't exist!" << endl;
    exit(1);
  }

  return h_pmask_dz;
}

Double_t PileupMaskingCorrection::GetTotalPmask() {
  if (!finished) {
    cerr << "Run the calculations before trying to save!" << endl;
    exit(1);
  }
  return pmask;
}

Double_t PileupMaskingCorrection::GetTotalPmaskError() {
  if (!finished) {
    cerr << "Run the calculations before trying to save!" << endl;
    exit(1);
  }

  return pmask_err;
}

TGraphErrors* PileupMaskingCorrection::GetMuMap() {
  if (!finished) {
    cerr << "Run the calculations before trying to save!" << endl;
    exit(1);
  }
  return tg_mu_obs_vs_mu_actual;
}

TGraphErrors* PileupMaskingCorrection::GetMuCorrection() {
  if (!finished) {
    cerr << "Run the calculations before trying to save!" << endl;
    exit(1);
  }
  return tg_pileup_correction;
}

void PileupMaskingCorrection::SetPmaskScale(float p_pmask_scale) {

  cout << "[PileupMaskingCorrection] INFO : Manually scaling p_mask by " << p_pmask_scale << endl;
  pmask_scale = p_pmask_scale;

}



/***************************************************************************************/
/// - Non-class helper functions
/***************************************************************************************/

Double_t PileupMaskingCorrection::fitFunc_gaussian_exclude(Double_t *x, Double_t *par) {
  //Gaussian, excluding a central region.
  // par[0] = normalization
  // par[1] = width
  // par[2] = center
  // par[3] = excluded region, lower edge
  // par[4] = excluded region, upper edge
  if (reject && (x[0] > par[3]) && (x[0] < par[4])) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] / (TMath::Sqrt(2. * TMath::Pi()) * par[1]) * TMath::Exp(-1. * TMath::Power(x[0] - par[2], 2) / (2. * TMath::Power(par[1], 2)));
}

Double_t PileupMaskingCorrection::fitFunc_gaussian(Double_t *x, Double_t *par) {
  // par[0] = normalization
  // par[1] = width
  // par[2] = center
  return par[0] / (TMath::Sqrt(2. * TMath::Pi()) * par[1]) * TMath::Exp(-1. * TMath::Power(x[0] - par[2], 2) / (2. * TMath::Power(par[1], 2)));
}

Double_t PileupMaskingCorrection::fitFunc_generic(Double_t *x, Double_t *par) {
  // par[0] = normalization
  return par[0] * dz_template->GetBinContent(dz_template->FindBin(x[0]));
}

Double_t PileupMaskingCorrection::fitFunc_generic_exclude(Double_t *x, Double_t *par) {
  // par[0] = normalization
  // par[1] = excluded region, lower edge
  // par[2] = excluded region, upper edge
  if (reject && (x[0] > par[1]) && (x[0] < par[2])) {
    TF1::RejectPoint();
    return 0;
  }
  return fitFunc_generic(x, par);
}

Double_t PileupMaskingCorrection::meanMuObs(Double_t x, Double_t par) {
  // par = p_mask
  // x = mu_actual
  // returns average observed mu.
  Double_t mu_obs = 0.;
  //Double_t MaxNGenInt = 110;
  if( (is_MC == 0) ){
    MaxNGenInt = 40;
  }else{
    MaxNGenInt = 110;
  }
  //if( (is_MC == 0) ){
    //cout << "[PileupMaskingCorrection] INFO: DATA" << endl;
    for (Int_t n_gen = 1; n_gen < MaxNGenInt; n_gen++) {
      Double_t poisson_factor = TMath::Exp(-1.*x) * TMath::Power(x, n_gen) / TMath::Factorial(n_gen);
      mu_obs += poisson_factor * meanNObs(par, n_gen);
    }
    //}else{
    //cout << "[PileupMaskingCorrection] INFO: MC" << endl;
    //mu_obs = meanNObs( par, x );
    //}
  return mu_obs;
}

Double_t PileupMaskingCorrection::meanNObs(Double_t p_mask, Double_t n_gen) {
  // Average number of vertices reconstructed vs. number of vertices generated, given p_mask.

  Double_t n_obs = 1.; // Start with 1, e.g. the primary vertex
  std::map<Int_t, Double_t> p_n; // Probability nth vertex observed
  p_n[1] = 1.;

  for (Int_t vtx = 2; vtx<=n_gen; vtx++) {
    p_n[vtx] = p_n[vtx-1] * (1. - p_mask * p_n[vtx-1]);
    n_obs += p_n[vtx];
  }
  return n_obs;
}

PileupMaskingCorrection::~PileupMaskingCorrection() {
  if( h_dz != 0 ) { delete h_dz; }
  if( h_dz_rebinned != 0 ) { delete h_dz_rebinned; }
  if( h_dz_expected != 0 ) { delete h_dz_expected; }
  if( f_dz_excluded != 0 ) { delete f_dz_excluded; }
  if( f_dz_full != 0 ) { delete f_dz_full; }
  if( f_dz_full_normalized != 0 ) { delete f_dz_full_normalized; }
  if( h_dz_rebinned_residuals != 0 ) { delete h_dz_rebinned_residuals; }
  if( h_pmask_dz != 0 ) { delete h_pmask_dz; }
  if( tg_pmask_cumulative != 0 ) { delete tg_pmask_cumulative; }
}

#endif
