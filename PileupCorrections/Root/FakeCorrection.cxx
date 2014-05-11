#ifndef FakeCorrection_cxx
#define FakeCorrection_cxx

#include "PileupCorrections/FakeCorrection.h"

using namespace std;

FakeCorrection::FakeCorrection(TString p_tag, Int_t p_ntrkcut) {

#ifdef INCLUDE_HACKS
  cout << "[FakeCorrection] WARNING : Including some hacks! Results should not be treated as final!!" << endl;
#endif

  /**
    *  This is the current preferred constructor.
    *  p_tag = mc_<7|8>TeV_17.2_<default|VtxLumi>
    *  p_ntrkcut = <2|5|7|10>
    */

  cout << "[FakeCorrection] INFO : Initializing fake correction with tag: " << p_tag << endl;
  cout << "[FakeCorrection] INFO : NTrkCut = " << p_ntrkcut << endl;

  std::vector<TString> known_samples;
  known_samples.push_back("mc_7TeV_17.2_normal_pythia8_pu");
  known_samples.push_back("mc_7TeV_17.2_VtxLumi_pythia8_pu");
  known_samples.push_back("mc_8TeV_17.2_normal_pythia8_pu");
  known_samples.push_back("mc_8TeV_17.2_VtxLumi_pythia8_pu");
  known_samples.push_back("mc_7TeV_17.2_normal_pythia8_pu_bs45");
  known_samples.push_back("mc_7TeV_17.2_VtxLumi_pythia8_pu_bs45");
  known_samples.push_back("mc_7TeV_17.2_normal_pythia8_pu_bs55");
  known_samples.push_back("mc_8TeV_17.2_VtxLumi_2newsets");
  known_samples.push_back("mc_8TeV_17.2_VtxLumi_2newsets_15Nov");
  known_samples.push_back("mc_8TeV_17.2_VtxLumi_mumax20");
  known_samples.push_back("mc_8TeV_17.2_VtxLumi_mumax75");

  if (find(known_samples.begin(), known_samples.end(), p_tag) == known_samples.end()) {
    cerr << "[FakeCorrection] ERROR : Requested tag " << p_tag << " not recognized. Exiting..." << endl;
    exit(1);
  }

  TString path = gs.path_fakeCorrection;
  path += "/";
  path += p_tag;
  path += "/fakerates.root";
  TFile *f_in = new TFile(path, "READ");

  TString name;
  //name = "tg_fakefraction_mu_apmc_NTrk"; name += p_ntrkcut;
  //tg_fakefraction_mu_apmc = (TGraphErrors*)f_in->Get(name);
  name = "tg_fake_fraction_MuReconMC_NTrk";
  name += p_ntrkcut;
  tg_fakefraction_MuReconMC = (TGraphErrors*)f_in->Get(name);

  name = "tg_mu_fake_MuReconMC_NTrk";
  name += p_ntrkcut;
  tg_fakemu_MuReconMC = (TGraphErrors*)f_in->Get(name);

  if (!tg_fakefraction_MuReconMC) {
    cerr << "FakeCorrection::FakeCorrection : [WARNING] Requested TGraph " << name << " was not found in input file " << path << "! Exiting..." << endl;
    exit(1);
  }

  name = "tg_nevt_mu_real_vs_mu_recon_NTrk";
  name += p_ntrkcut;
  tg_nevt_mu_real_vs_mu_recon = (TGraphErrors*)f_in->Get(name);

  if (!tg_nevt_mu_real_vs_mu_recon) {
    cerr << "[FakeCorrection] ERROR : TGraph " << name << " was not found in input file. Exiting..." << endl;
    exit(1);
  }

  mu_scale = 1.0;
  mu_scale_uncertainty = 0.0;
}

Double_t FakeCorrection::GetFakeFractionFromMuReconMC(Double_t mu_in) {
  /**
    *  This is my current "best guess" for how to apply the correction.
    *  Work with quantities after the masking correction is applied. We need to use the MC pretty heavily to get the fake correction, but masking is very different between data and MC.
    *  So, try to undo the masking correction on both quantities to get an estimate of the fake fraction.
    */

  mu_in = mu_in / mu_scale;
  Double_t val;
  if (mu_in < tg_fakefraction_MuReconMC->GetX()[0]) {
    val = 0.;
  } else {
    val = tg_fakefraction_MuReconMC->Eval(mu_in);
  }
  return val;
}

Double_t FakeCorrection::GetFakeMuFromMuReconMC(Double_t mu_in) {

  mu_in = mu_in / mu_scale;
  Double_t val;
  if (mu_in < tg_fakemu_MuReconMC->GetX()[0]) {
    val = 0.;
  } else {
    val = tg_fakemu_MuReconMC->Eval(mu_in);
  }
  return val;

}

Double_t FakeCorrection::GetFakeMuUncertaintyFromMuReconMC(Double_t mu_in) {

  mu_in = mu_in / mu_scale;
  Double_t val;
  if (mu_in < tg_fakemu_MuReconMC->GetX()[0]) {
    val = 0.;
  } else {
    TGraphErrors *tg_fakemu_fluctuate = new TGraphErrors(tg_fakemu_MuReconMC->GetN());
    for (int i = 0; i < tg_fakemu_MuReconMC->GetN(); i++) {
      tg_fakemu_fluctuate->SetPoint(i, tg_fakemu_MuReconMC->GetX()[i], tg_fakemu_MuReconMC->GetY()[i] + tg_fakemu_MuReconMC->GetEY()[i]);
      tg_fakemu_fluctuate->SetPointError(i, tg_fakemu_MuReconMC->GetEX()[i], tg_fakemu_MuReconMC->GetEY()[i]);
    }
    val = tg_fakemu_fluctuate->Eval(mu_in) - tg_fakemu_MuReconMC->Eval(mu_in);
  }
  return val;
}

Double_t FakeCorrection::GetFakeMuHighFromMuReconMC(Double_t mu_in) {

  mu_in = mu_in / (mu_scale - mu_scale_uncertainty);
  return tg_fakemu_MuReconMC->Eval(mu_in);

}

Double_t FakeCorrection::GetFakeMuLowFromMuReconMC(Double_t mu_in) {

  mu_in = mu_in / (mu_scale + mu_scale_uncertainty);
  return tg_fakemu_MuReconMC->Eval(mu_in);

}

Double_t FakeCorrection::GetMuFakeFromMuReconNEvt(Double_t mu_in) {

  mu_in = mu_in / mu_scale;
  return mu_in - tg_nevt_mu_real_vs_mu_recon->Eval(mu_in);

}

Double_t FakeCorrection::GetMuFakeUncertaintyFromMuReconNEvt(Double_t mu_in) {

  mu_in = mu_in / mu_scale;
  Double_t val;
  if (mu_in < tg_nevt_mu_real_vs_mu_recon->GetX()[0]) {
    val = 0.;
  } else {
    TGraphErrors *tg_nevt_mu_real_vs_mu_recon_fluctuate = new TGraphErrors(tg_nevt_mu_real_vs_mu_recon->GetN());
    for (int i = 0; i < tg_nevt_mu_real_vs_mu_recon->GetN(); i++) {
      tg_nevt_mu_real_vs_mu_recon_fluctuate->SetPoint(i, tg_nevt_mu_real_vs_mu_recon->GetX()[i], tg_nevt_mu_real_vs_mu_recon->GetY()[i] + tg_nevt_mu_real_vs_mu_recon->GetEY()[i]);
      tg_nevt_mu_real_vs_mu_recon_fluctuate->SetPointError(i, tg_nevt_mu_real_vs_mu_recon->GetEX()[i], tg_nevt_mu_real_vs_mu_recon->GetEY()[i]);
    }
    val = TMath::Abs(tg_nevt_mu_real_vs_mu_recon_fluctuate->Eval(mu_in) - tg_nevt_mu_real_vs_mu_recon->Eval(mu_in));
  }
  return val;

}


TGraphErrors* FakeCorrection::GetFakeCorrectionTGraph() {

  return tg_fakefraction_MuReconMC;

}

TGraph* FakeCorrection::ConvertTH1ToTGraph(TH1D *h_in) {
  Int_t n_points = 0;
  for (int bin = 1; bin <= h_in->GetXaxis()->GetNbins(); bin++) {
    if (h_in->GetBinContent(bin) > 0) {
      n_points++;
    }
  }

  TGraph *tg_out = new TGraph(n_points);
  TString old_name = (TString)(h_in->GetName());
  TString name = old_name.ReplaceAll("h_",1,"tg_",1);
  Int_t current_point = 0;
  for (int bin = 1; bin <= h_in->GetXaxis()->GetNbins(); bin++) {
    if (h_in->GetBinContent(bin) > 0) {
      tg_out->SetPoint(current_point, h_in->GetBinCenter(bin), h_in->GetBinContent(bin));
      current_point++;
    }
  }

  return tg_out;
}

void FakeCorrection::SetMuScale(float p_mu_scale) {

  cout << "[FakeCorrection] INFO : Setting mu scale to " << p_mu_scale << endl;
  mu_scale = p_mu_scale;

}

void FakeCorrection::SetMuScaleUncertainty(float p_mu_scale_uncertainty) {

  cout << "[FakeCorrection] INFO : Setting mu scale uncertainty to " << p_mu_scale_uncertainty << endl;
  mu_scale_uncertainty = p_mu_scale_uncertainty;

}


#endif
