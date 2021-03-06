#define PLOT_MC_TRUTH
#define AnaVtxTree_cxx

//Local
#include "D3PDMC/AnaVtxTree.h"
#include "D3PDMC/ATLASImport.h"
#include "D3PDMC/HistogramHelper.h"

//STL

#include <iostream>
#include <fstream>
#include <sstream>
#include <bitset>
#include <string>
#include <iomanip>

// ROOT
#include "TDataType.h"
#include "TMath.h"
#include "TAxis.h"
#include "TVector2.h"
#include "TObjString.h"
#include "TParticlePDG.h"

//Settings depending on VtxTree
//#define D3PD_PURITY_DEFINED

AnaVtxTree::AnaVtxTree() : VtxTree() {
  //Default settings
  performD0PhiFit = false;
  VtxEffCriteriaZ = 0.3; //300 um
  VtxEffCriteriaT = 5; //5 sigma
  debugLevel = 0;
  SetQualityVertexVersion("1.0.2"); // Default: Tagged vertex with >= 2 tracks
  //SetQualityVertexVersion("1.1.5"); // 
  cfgMaxDistance = 5.0;
  cfgMetric = VtxDistM_deltaZsig;
  vtxTMWThreshold = 0.7;
  vtxTMWThresholdStore = 1 - (1 - vtxTMWThreshold) / 10.;
  vtxTMTrackMatchProb = 0.7; //matching probability for tracks

  trigMetaDataTree = 0;
  triggerTool = 0;

  muFilter_min = -1; //disabled
  muFilter_max = -1; //disabled

  //nTrkCuts.push_back(2);
  nTrkCuts.push_back(3);
  nTrkCuts.push_back(4);
  nTrkCuts.push_back(5);
  //nTrkCuts.push_back(6);
  //nTrkCuts.push_back(7);
  //nTrkCuts.push_back(8);
  //nTrkCuts.push_back(10);

  dzsig_max = 0.; // 0 = no cut
  max_chi2ndf = -1.; // -1 = no cut

  skip_hs = false;
  do_timing = false;

};

Int_t AnaVtxTree::SetDebugLevel(Int_t debug) {
  std::cout << "[AnaVtxTree] Set debugLevel to " << debug << endl;
  return (debugLevel = debug);
}

void AnaVtxTree::SlaveBegin(TTree *tree) {

  VtxTree::SlaveBegin(tree);

  std::cout << "[AnaVtxTree] Init" << std::endl;

  // --- Init output text file
  outputTxt.open((outputFileName + TString(".txt")).Data());

  // --- Create histograms
  h_NTrig_NGenInt = HistogramHelper::defineHistogram("NTrig vs. NGenInt", MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "NGenInt", "NTrig", true, "h_NTrig_NGenInt");

  h_actualInt_NGenInt = HistogramHelper::define2DHistogram("actualInt vs. NGenInt", MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "actualInt","NGenInt", true, "h_actualInt_NGenInt");
  h_actualInt_mcvtxn = HistogramHelper::define2DHistogram("actualInt vs. mcvtx_n", MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "actualInt","mcvtx_n", true, "h_actualInt_mcvtx_n");

  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {

    TString hname, htitle;
    TString name_tag = cstr("_NTrk",*nTrkCut);
    //TString name_tag = "_NTrk";
    //name_tag += *nTrkCut;
    //TString title_tag = cstr(" (selection NTrk",str(*nTrkCut,")"));
    TString title_tag = " (selection NTrk";
    title_tag += *nTrkCut;
    title_tag += ")";

    hname = "h_NTrig_NVtxRecon";
    hname += name_tag;
    htitle = "Number of events vs. NVtxRecon";
    htitle += title_tag;
    h_NTrig_NVtxRecon[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNVtxRecon + 1, -0.5, MaxNVtxRecon + 0.5, "NVtxRecon", "NTrig", true, hname);

    hname = "h_truefakes_NGenInt";
    hname += name_tag;
    htitle = "Number of true fakes vs. NGenInt";
    htitle += title_tag;
    h_truefakes_NGenInt[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "NGenInt", "True Fake Vertices", true, hname);

    hname = "h_fakes_NGenInt";
    hname += name_tag;
    htitle = "Number of fakes vs. NGenInt";
    htitle += title_tag;
    h_fakes_NGenInt[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "NGenInt", "Fake Vertices", true, hname);

    hname = "h_splits_NGenInt";
    hname += name_tag;
    htitle = "Number of splits vs. NGenInt";
    htitle += title_tag;
    h_splits_NGenInt[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "NGenInt", "Split Vertices", true, hname);

    hname = "h_real_NGenInt";
    hname += name_tag;
    htitle = "Number of real vertices vs. NGenInt";
    htitle += title_tag;
    h_real_NGenInt[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "NGenInt", "Real Vertices", true, hname);

    hname = "h_all_NGenInt";
    hname += name_tag;
    htitle = "Total number of reconstructed vertices vs. NGenInt";
    htitle += title_tag;
    h_all_NGenInt[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "NGenInt", "All Vertices", true, hname);

    hname = "h_all_actualInt";
    hname += name_tag;
    htitle = "Total number of reconstructed vertices vs. ei_actualIntPerXing";
    htitle += title_tag;
    h_all_actualInt[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "ei_actualIntPerXing", "All Vertices", true, hname);

    hname = "h_all_dz_NGenInt";
    hname += name_tag;
    htitle = "Recon. vertex-pair dz vs. dz significance vs. number of generated interactions";
    htitle += title_tag;
    h_all_dz_NGenInt[*nTrkCut] = HistogramHelper::define2DHistogram(htitle, 4000, -500, 500, MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "#Deltaz (mm)", "NGenInt", true, hname);

    hname = "h_fakes_NVtxRecon";
    hname += name_tag;
    htitle = "Number of fakes vs. NVtxRecon";
    htitle += title_tag;
    h_fakes_NVtxRecon[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNVtxRecon + 1, -0.5, MaxNVtxRecon + 0.5, "NVtxRecon", "Fake Vertices", true, hname);

    hname = "h_truefakes_NVtxRecon";
    hname += name_tag;
    htitle = "Number of true fakes vs. NVtxRecon";
    htitle += title_tag;
    h_truefakes_NVtxRecon[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNVtxRecon + 1, -0.5, MaxNVtxRecon + 0.5, "NVtxRecon", "True Fake Vertices", true, hname);

    hname = "h_splits_NVtxRecon";
    hname += name_tag;
    htitle = "Number of splits vs. NVtxRecon";
    htitle += title_tag;
    h_splits_NVtxRecon[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNVtxRecon + 1, -0.5, MaxNVtxRecon + 0.5, "NVtxRecon", "Split Vertices", true, hname);

    hname = "h_real_NVtxRecon";
    hname += name_tag;
    htitle = "Number of real vertices vs. NVtxRecon";
    htitle += title_tag;
    h_real_NVtxRecon[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNVtxRecon + 1, -0.5, MaxNVtxRecon + 0.5, "NVtxRecon", "Real Vertices", true, hname);

    hname = "h_all_NVtxRecon";
    hname += name_tag;
    htitle = "Total number of reconstructed vertices vs. NVtxRecon";
    htitle += title_tag;
    h_all_NVtxRecon[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNVtxRecon + 1, -0.5, MaxNVtxRecon + 0.5, "NVtxRecon", "All Vertices", true, hname);

    hname = "h_all_dz_NVtxRecon";
    hname += name_tag;
    htitle = "Recon. vertex-pair dz vs. dz significance vs. number of generated interactions";
    htitle += title_tag;
    h_all_dz_NVtxRecon[*nTrkCut] = HistogramHelper::define2DHistogram(htitle, 4000, -500, 500, MaxNVtxRecon + 1, -0.5, MaxNVtxRecon + 0.5, "#Deltaz (mm)", "NVtxRecon", true, hname);

    hname = "h_fake_events_NGenInt";
    hname += name_tag;
    htitle = "Fake events vs. NGenInt";
    htitle += title_tag;
    h_fake_events_NGenInt[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNVtxRecon + 1, -0.5, MaxNVtxRecon + 0.5, "NVtxRecon", "Fake Events", true, hname);

    hname = "h_all_events_NGenInt";
    hname += name_tag;
    htitle = "All events vs. NGenInt";
    htitle += title_tag;
    h_all_events_NGenInt[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNVtxRecon + 1, -0.5, MaxNVtxRecon + 0.5, "NVtxRecon", "Events", true, hname);

    hname = "h_fake_events_NVtxRecon";
    hname += name_tag;
    htitle = "Fake events vs. NVtxRecon";
    htitle += title_tag;
    h_fake_events_NVtxRecon[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNVtxRecon + 1, -0.5, MaxNVtxRecon + 0.5, "NVtxRecon", "Fake Events", true, hname);

    hname = "h_all_events_NVtxRecon";
    hname += name_tag;
    htitle = "All events vs. NVtxRecon";
    htitle += title_tag;
    h_all_events_NVtxRecon[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNVtxRecon + 1, -0.5, MaxNVtxRecon + 0.5, "NVtxRecon", "Events", true, hname);

    hname = "h_event_contains_fake_NGenInt";
    hname += name_tag;
    htitle = "Events with at least one fake vertex vs. NGenInt";
    htitle += title_tag;
    h_event_contains_fake_NGenInt[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "NGenInt", "Events with #geq1 fake", true, hname);

    hname = "h_event_contains_any_NGenInt";
    hname += name_tag;
    htitle = "Events with at least one vertex vs. NGenInt";
    htitle += title_tag;
    h_event_contains_any_NGenInt[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "NGenInt", "Events", true, hname);

    hname = "h_event_contains_fake_NVtxRecon";
    hname += name_tag;
    htitle = "Events with at least one fake vertex vs. NVtxRecon";
    htitle += title_tag;
    h_event_contains_fake_NVtxRecon[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNVtxRecon + 1, -0.5, MaxNVtxRecon + 0.5, "NVtxRecon", "Events with #geq1 fake", true, hname);

    hname = "h_event_contains_any_NVtxRecon";
    hname += name_tag;
    htitle = "Events with at least one vertex vs. NVtxRecon";
    htitle += title_tag;
    h_event_contains_any_NVtxRecon[*nTrkCut] = HistogramHelper::defineHistogram(htitle, MaxNVtxRecon + 1, -0.5, MaxNVtxRecon + 0.5, "NVtxRecon", "Events", true, hname);

  }

  h_privtx_z_mu = HistogramHelper::define2DHistogram("Primary vertex Z vs. mu", 4000, -500., 500, MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "Z (mm)", "#mu", true, "h_privtx_z_mu");

  h_privtx_z_NGenInt = HistogramHelper::define2DHistogram("Primary vertex Z vs. NGenInt", 4000, -500., 500., MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "Z (mm)", "NGenInt", true, "h_privtx_z_NGenInt");

  h_privtx_z_NVtxRecon = HistogramHelper::define2DHistogram("Primary vertex Z vs. NVtxRecon", 4000, -500., 500., MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "Z (mm)", "NVtxRecon", true, "h_privtx_z_NVtxRecon");

  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {

    TString hname, htitle;
    TString name_tag = "_NTrk";
    name_tag += *nTrkCut;
    TString title_tag = " (selection NTrk";
    title_tag += *nTrkCut;
    title_tag += ")";

    hname = "h_split_dz_NGenInt";
    hname += name_tag;
    htitle = "Truth-level splits, Delta Z vs. NGenInt";
    htitle += title_tag;
    h_split_dz_NGenInt[*nTrkCut] = HistogramHelper::define2DHistogram(htitle,  4000, -500, 500, MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "#Deltaz (mm)", "NGenInt", true, hname);

    hname = "h_split_dzsig_NGenInt";
    hname += name_tag;
    htitle = "Truth-level splits, Delta Z significance vs. NGenInt";
    htitle += title_tag;
    h_split_dzsig_NGenInt[*nTrkCut] = HistogramHelper::define2DHistogram(htitle,  4000, -100., 100., MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "#frac{#Deltaz}{#sigma_{#Deltaz}}", "NGenInt", true, hname);

    hname = "h_split_dz_NVtxRecon";
    hname += name_tag;
    htitle = "Truth-level splits, Delta Z vs. NVtxRecon";
    htitle += title_tag;
    h_split_dz_NVtxRecon[*nTrkCut] = HistogramHelper::define2DHistogram(htitle,  4000, -500, 500, MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "#Deltaz (mm)", "NVtxRecon", true, hname);

    hname = "h_split_dzsig_NVtxRecon";
    hname += name_tag;
    htitle = "Truth-level splits, Delta Z significance vs. NVtxRecon";
    htitle += title_tag;
    h_split_dzsig_NVtxRecon[*nTrkCut] = HistogramHelper::define2DHistogram(htitle,  4000, -100., 100., MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "#frac{#Deltaz}{#sigma_{#Deltaz}}", "NVtxRecon", true, hname);

    hname = "h_split_dz_dzsig";
    hname += name_tag;
    htitle = "Truth-level splits, Delta Z vs. Delta Z significance";
    htitle += title_tag;
    h_split_dz_dzsig[*nTrkCut] = HistogramHelper::define2DHistogram(htitle, 4000, -500, 500, 4000, -100., 100., "#Deltaz", "#frac{#Deltaz}{#sigma_{#Deltaz}}", true, hname);

    hname = "h_real_dz_NGenInt";
    hname += name_tag;
    htitle = "Truth-level real pairs, Delta Z vs. NGenInt";
    htitle += title_tag;
    h_real_dz_NGenInt[*nTrkCut] = HistogramHelper::define2DHistogram(htitle,  1600, -500, 500, MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "#Deltaz (mm)", "NGenInt", true, hname);

    hname = "h_real_dzsig_NGenInt";
    hname += name_tag;
    htitle = "Truth-level real pairs, Delta Z significance vs. NGenInt";
    htitle += title_tag;
    h_real_dzsig_NGenInt[*nTrkCut] = HistogramHelper::define2DHistogram(htitle,  4000, -100., 100., MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "#frac{#Deltaz}{#sigma_{#Deltaz}}", "NGenInt", true, hname);

    hname = "h_real_dz_NVtxRecon";
    hname += name_tag;
    htitle = "Truth-level real pairs, Delta Z vs. NVtxRecon";
    htitle += title_tag;
    h_real_dz_NVtxRecon[*nTrkCut] = HistogramHelper::define2DHistogram(htitle,  1600, -500, 500, MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "#Deltaz (mm)", "NVtxRecon", true, hname);

    hname = "h_real_dzsig_NVtxRecon";
    hname += name_tag;
    htitle = "Truth-level real pairs, Delta Z significance vs. NVtxRecon";
    htitle += title_tag;
    h_real_dzsig_NVtxRecon[*nTrkCut] = HistogramHelper::define2DHistogram(htitle,  4000, -100., 100., MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, "#frac{#Deltaz}{#sigma_{#Deltaz}}", "NVtxRecon", true, hname);

    hname = "h_real_dz_dzsig";
    hname += name_tag;
    htitle = "Truth-level real pairs, Delta Z vs. Delta Z significance";
    htitle += title_tag;
    h_real_dz_dzsig[*nTrkCut] = HistogramHelper::define2DHistogram(htitle, 4000, -500, 500, 4000, -100., 100., "#Deltaz", "#frac{#Deltaz}{#sigma_{#Deltaz}}", true, hname);

    //------------------------------------------------------

    hname = "h_fakes_x";
    hname += name_tag;
    htitle = "Fake vertices: x";
    htitle += title_tag;
    h_fakes_x[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                          NBinsXAxis, LowXAxis, HighXAxis, "x (mm)", "", true, hname);

    hname = "h_fakes_y";
    hname += name_tag;
    htitle = "Fake vertices: y";
    htitle += title_tag;
    h_fakes_y[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                          NBinsYAxis, LowYAxis, HighYAxis, "y (mm)", "", true, hname);

    hname = "h_fakes_z";
    hname += name_tag;
    htitle = "Fake vertices: z";
    htitle += title_tag;
    h_fakes_z[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                          NBinsZAxis, LowZAxis, HighZAxis, "z (mm)", "", true, hname);

    hname = "h_fakes_cov_xx";
    hname += name_tag;
    htitle = "Fake vertices: cov_xx";
    htitle += title_tag;
    h_fakes_cov_xx[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                               100, 0., 1., "cov_{xx}", "", true, hname);

    hname = "h_fakes_cov_yy";
    hname += name_tag;
    htitle = "Fake vertices: cov_yy";
    htitle += title_tag;
    h_fakes_cov_yy[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                               100, 0., 1., "cov_{yy}", "", true, hname);

    hname = "h_fakes_cov_zz";
    hname += name_tag;
    htitle = "Fake vertices: cov_zz";
    htitle += title_tag;
    h_fakes_cov_zz[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                               200, 0., 2., "cov_{zz}", "", true, hname);

    hname = "h_fakes_cov_xy";
    hname += name_tag;
    htitle = "Fake vertices: cov_xy";
    htitle += title_tag;
    h_fakes_cov_xy[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                               100, 0., 1., "cov_{xy}", "", true, hname);

    hname = "h_fakes_cov_xz";
    hname += name_tag;
    htitle = "Fake vertices: cov_xz";
    htitle += title_tag;
    h_fakes_cov_xz[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                               100, 0., 1., "cov_{xz}", "", true, hname);

    hname = "h_fakes_cov_yz";
    hname += name_tag;
    htitle = "Fake vertices: cov_yz";
    htitle += title_tag;
    h_fakes_cov_yz[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                               200, 0., 2., "cov_{yz}", "", true, hname);

    hname = "h_fakes_chi2ndf";
    hname += name_tag;
    htitle = "Fake vertices: chi2 / ndf";
    htitle += title_tag;
    h_fakes_chi2ndf[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                100, 0., 30., "#chi^{2}/NDF", "", true, hname);

    hname = "h_fakes_chi2";
    hname += name_tag;
    htitle = "Fake vertices: chi2";
    htitle += title_tag;
    h_fakes_chi2[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                             100, 0., 100., "#chi^{2}", "", true, hname);

    hname = "h_fakes_ndf";
    hname += name_tag;
    htitle = "Fake vertices: NDF";
    htitle += title_tag;
    h_fakes_ndf[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                            100, 0., 100., "NDF", "", true, hname);

    hname = "h_fakes_ntrk";
    hname += name_tag;
    htitle = "Fake vertices: Number of tracks";
    htitle += title_tag;
    h_fakes_ntrk[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                             100, 0., 100., "NTrk", "", true, hname);

    hname = "h_fakes_sumPt";
    hname += name_tag;
    htitle = "Fake vertices: sumPt";
    htitle += title_tag;
    h_fakes_sumPt[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                              1000, 0., 500000., "#Sigma p_{T}", "", true, hname);

    //------------------------------------------------------

    hname = "h_splits_x";
    hname += name_tag;
    htitle = "Split vertices: x";
    htitle += title_tag;
    h_splits_x[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                           NBinsXAxis, LowXAxis, HighXAxis, "x (mm)", "", true, hname);

    hname = "h_splits_y";
    hname += name_tag;
    htitle = "Split vertices: y";
    htitle += title_tag;
    h_splits_y[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                           NBinsYAxis, LowYAxis, HighYAxis, "y (mm)", "", true, hname);

    hname = "h_splits_z";
    hname += name_tag;
    htitle = "Split vertices: z";
    htitle += title_tag;
    h_splits_z[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                           NBinsZAxis, LowZAxis, HighZAxis, "z (mm)", "", true, hname);

    hname = "h_splits_cov_xx";
    hname += name_tag;
    htitle = "Split vertices: cov_xx";
    htitle += title_tag;
    h_splits_cov_xx[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                100, 0., 1., "cov_{xx}", "", true, hname);

    hname = "h_splits_cov_yy";
    hname += name_tag;
    htitle = "Split vertices: cov_yy";
    htitle += title_tag;
    h_splits_cov_yy[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                100, 0., 1., "cov_{yy}", "", true, hname);

    hname = "h_splits_cov_zz";
    hname += name_tag;
    htitle = "Split vertices: cov_zz";
    htitle += title_tag;
    h_splits_cov_zz[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                200, 0., 2., "cov_{zz}", "", true, hname);

    hname = "h_splits_cov_xy";
    hname += name_tag;
    htitle = "Split vertices: cov_xy";
    htitle += title_tag;
    h_splits_cov_xy[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                100, 0., 1., "cov_{xy}", "", true, hname);

    hname = "h_splits_cov_xz";
    hname += name_tag;
    htitle = "Split vertices: cov_xz";
    htitle += title_tag;
    h_splits_cov_xz[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                100, 0., 1., "cov_{xz}", "", true, hname);

    hname = "h_splits_cov_yz";
    hname += name_tag;
    htitle = "Split vertices: cov_yz";
    htitle += title_tag;
    h_splits_cov_yz[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                200, 0., 2., "cov_{yz}", "", true, hname);

    hname = "h_splits_chi2ndf";
    hname += name_tag;
    htitle = "Split vertices: chi2 / ndf";
    htitle += title_tag;
    h_splits_chi2ndf[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                 100, 0., 30., "#chi^{2}/NDF", "", true, hname);

    hname = "h_splits_chi2";
    hname += name_tag;
    htitle = "Split vertices: chi2";
    htitle += title_tag;
    h_splits_chi2[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                              100, 0., 100., "#chi^{2}", "", true, hname);

    hname = "h_splits_ndf";
    hname += name_tag;
    htitle = "Split vertices: NDF";
    htitle += title_tag;
    h_splits_ndf[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                             100, 0., 100., "NDF", "", true, hname);

    hname = "h_splits_ntrk";
    hname += name_tag;
    htitle = "Split vertices: Number of tracks";
    htitle += title_tag;
    h_splits_ntrk[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                              100, 0., 100., "NTrk", "", true, hname);

    hname = "h_splits_sumPt";
    hname += name_tag;
    htitle = "Split vertices: sumPt";
    htitle += title_tag;
    h_splits_sumPt[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                               1000, 0., 500000., "#Sigma p_{T}", "", true, hname);


    //------------------------------------------------------

    hname = "h_real_x";
    hname += name_tag;
    htitle = "Real vertices: x";
    htitle += title_tag;
    h_real_x[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                         NBinsXAxis, LowXAxis, HighXAxis, "x (mm)", "", true, hname);

    hname = "h_real_y";
    hname += name_tag;
    htitle = "Real vertices: y";
    htitle += title_tag;
    h_real_y[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                         NBinsYAxis, LowYAxis, HighYAxis, "y (mm)", "", true, hname);

    hname = "h_real_z";
    hname += name_tag;
    htitle = "Real vertices: z";
    htitle += title_tag;
    h_real_z[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                         NBinsZAxis, LowZAxis, HighZAxis, "z (mm)", "", true, hname);

    hname = "h_real_cov_xx";
    hname += name_tag;
    htitle = "Real vertices: cov_xx";
    htitle += title_tag;
    h_real_cov_xx[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                              100, 0., 1., "cov_{xx}", "", true, hname);

    hname = "h_real_cov_yy";
    hname += name_tag;
    htitle = "Real vertices: cov_yy";
    htitle += title_tag;
    h_real_cov_yy[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                              100, 0., 1., "cov_{yy}", "", true, hname);

    hname = "h_real_cov_zz";
    hname += name_tag;
    htitle = "Real vertices: cov_zz";
    htitle += title_tag;
    h_real_cov_zz[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                              200, 0., 2., "cov_{zz}", "", true, hname);

    hname = "h_real_cov_xy";
    hname += name_tag;
    htitle = "Real vertices: cov_xy";
    htitle += title_tag;
    h_real_cov_xy[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                              100, 0., 1., "cov_{xy}", "", true, hname);

    hname = "h_real_cov_xz";
    hname += name_tag;
    htitle = "Real vertices: cov_xz";
    htitle += title_tag;
    h_real_cov_xz[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                              100, 0., 1., "cov_{xz}", "", true, hname);

    hname = "h_real_cov_yz";
    hname += name_tag;
    htitle = "Real vertices: cov_yz";
    htitle += title_tag;
    h_real_cov_yz[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                              200, 0., 2., "cov_{yz}", "", true, hname);

    hname = "h_real_chi2ndf";
    hname += name_tag;
    htitle = "Real vertices: chi2 / ndf";
    htitle += title_tag;
    h_real_chi2ndf[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                               100, 0., 30., "#chi^{2}/NDF", "", true, hname);

    hname = "h_real_chi2";
    hname += name_tag;
    htitle = "Real vertices: chi2";
    htitle += title_tag;
    h_real_chi2[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                            100, 0., 100., "#chi^{2}", "", true, hname);

    hname = "h_real_ndf";
    hname += name_tag;
    htitle = "Real vertices: NDF";
    htitle += title_tag;
    h_real_ndf[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                           100, 0., 100., "NDF", "", true, hname);

    hname = "h_real_ntrk";
    hname += name_tag;
    htitle = "Real vertices: Number of tracks";
    htitle += title_tag;
    h_real_ntrk[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                                            100, 0., 100., "NTrk", "", true, hname);

    hname = "h_real_sumPt";
    hname += name_tag;
    htitle = "Real vertices: sumPt";
    htitle += title_tag;
    h_real_sumPt[*nTrkCut] = HistogramHelper::defineHistogram(htitle,
                             1000, 0., 500000., "#Sigma p_{T}", "", true, hname);

    hname = "h_GoodVertices_actualInt";
    hname += name_tag;
    htitle = "Good Vertices vs actualIntPerXing with NTracks > ";
    htitle += title_tag;
    h_GoodVertices_actualInt[*nTrkCut] = HistogramHelper::defineProfile(htitle,
                              MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, 0., 50., "ei_actualIntPerXing", "GoodVertices", true, hname);

    hname = "h_GoodVertices_actualInt_2D";
    hname += name_tag;
    htitle = "Good Vertices vs actualIntPerXing with NTracks > ";
    htitle += title_tag;
    h_GoodVertices_actualInt_2D[*nTrkCut] = HistogramHelper::define2DHistogram(htitle,
                              MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, 50, 0., 50., "ei_actualIntPerXing", "GoodVertices", true, hname);

    hname = "h_GoodVertices_NGenInt";
    hname += name_tag;
    htitle = "Good Vertices vs NGenInt with NTracks > ";
    htitle += title_tag;
    h_GoodVertices_NGenInt[*nTrkCut] = HistogramHelper::defineProfile(htitle,
                              MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, 0., 50., "NGenInt", "GoodVertices", true, hname);

    hname = "h_GoodVertices_NGenInt_2D";
    hname += name_tag;
    htitle = "Good Vertices vs NGenInt with NTracks > ";
    htitle += title_tag;
    h_GoodVertices_NGenInt_2D[*nTrkCut] = HistogramHelper::define2DHistogram(htitle,
                              MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, 50, 0., 50., "NGenInt", "GoodVertices", true, hname);

    hname = "h_GoodVertices_mcvtxn";
    hname += name_tag;
    htitle = "Good Vertices vs mcvtx_n with NTracks > ";
    htitle += title_tag;
    h_GoodVertices_mcvtxn[*nTrkCut] = HistogramHelper::defineProfile(htitle,
                              MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, 0., 50., "mcvtx_n", "GoodVertices", true, hname);

    hname = "h_GoodVertices_mcvtxn_2D";
    hname += name_tag;
    htitle = "Good Vertices vs mcvtx_n with NTracks > ";
    htitle += title_tag;
    h_GoodVertices_mcvtxn_2D[*nTrkCut] = HistogramHelper::define2DHistogram(htitle,
                              MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, 50, 0., 50., "mcvtx_n", "GoodVertices", true, hname);

    hname = "h_NVtxRecon_actualInt";
    hname += name_tag;
    htitle = "Reconstructed Vertices vs actualIntPerXing ";
    htitle += title_tag;
    h_NVtxRecon_actualInt[*nTrkCut] = HistogramHelper::defineProfile(htitle,
                              MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, 0., 50., "ei_actualIntPerXing", "NVtxRecon", true, hname);

    hname = "h_NVtxRecon_NGenInt";
    hname += name_tag;
    htitle = "Reconstructed Vertices vs NGenInt";
    htitle += title_tag;
    h_NVtxRecon_NGenInt[*nTrkCut] = HistogramHelper::defineProfile(htitle,
                              MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, 0., 50., "NGenInt", "NVtxRecon", true, hname);

    hname = "h_NReal_actualInt";
    hname += name_tag;
    htitle = "Real Reconstructed Vertices vs actualIntPerXing ";
    htitle += title_tag;
    h_NReal_actualInt[*nTrkCut] = HistogramHelper::defineProfile(htitle,
                              MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, 0., 50., "ei_actualIntPerXing", "NReal", true, hname);

    hname = "h_NReal_NGenInt";
    hname += name_tag;
    htitle = "Real Reconstructed Vertices vs NGenInt";
    htitle += title_tag;
    h_NReal_NGenInt[*nTrkCut] = HistogramHelper::defineProfile(htitle,
                              MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, 0., 50., "NGenInt", "NReal", true, hname);

    hname = "h_NReal_NVtxRecon";
    hname += name_tag;
    htitle = "Real Reconstructed Vertices vs Reconstructed Vertices";
    htitle += title_tag;
    h_NReal_NVtxRecon[*nTrkCut] = HistogramHelper::defineProfile(htitle,
                              MaxNGenInt + 1, -0.5, MaxNGenInt + 0.5, 0., 50., "NVtxRecon", "NReal", true, hname);

  }

  h_truth_x = HistogramHelper::defineHistogram("Truth vertex X position", NBinsXAxis, LowXAxis, HighXAxis, "x (mm)", true, "h_truth_x");
  h_truth_y = HistogramHelper::defineHistogram("Truth vertex Y position", NBinsYAxis, LowYAxis, HighYAxis, "y (mm)", true, "h_truth_y");
  h_truth_z = HistogramHelper::defineHistogram("Truth vertex Z position", NBinsZAxis, LowZAxis, HighZAxis, "z (mm)", true, "h_truth_z");

  //True Matched vertex information: vertex->tracks->mcpart->mcvtx matching
  h_TM_goodUnmatched_pt = HistogramHelper::defineHistogram("P_{T} of matched tracks (prob > 0.7) with no gen particle saved.",
                                          200,0,10.,
                                          "p_{T} (GeV)", "Entries / 0.5 GeV",
                                          true, "TMgoodUnmatchedPt");
  h_TM_goodUnmatched_eta = HistogramHelper::defineHistogram("#eta of matched tracks (prob > 0.7) with no gen particle saved.",
                           60,-3.,3.,
                           "#eta", "Entries / 0.1",
                           true, "TMgoodUnmatchedEta");
  h_TM_goodUnmatched_pteta = HistogramHelper::define2DHistogram("P_{T}-#eta of matched tracks (prob > 0.7) with no gen particle saved.",
                             200,0,10., 60, -3., 3.,
                             "p_{T} (GeV)", "#eta",
                             true, "TMgoodUnmatchedPtEta");
  h_TM_highestW = HistogramHelper::defineHistogram("Highest (relative) weight from a single genVertex associated to a reconstructed vertex. If fakes, plot < 0",
                                  100, -1., 1., "Relative weight", "Entries",
                                  true, "TMhighestW");
  h_TM_2ndHighestW = HistogramHelper::defineHistogram("2nd Highest (relative) weight from a single genVertex associated to a reconstructed vertex. If fakes, plot < 0",
                                     100, -1., 1., "Relative weight", "Entries",
                                     true, "TM2ndHighestW");
  h_TM_OtherHighestW = HistogramHelper::defineHistogram("All other (>2nd) relative weight from genVertices associated to a reconstructed vertex. If fakes, plot < 0",
                                       100, -1., 1., "Relative weight", "Entries",
                                       true, "TMOtherHighestW");
  h_TM_numGenMatched = HistogramHelper::defineHistogram("Number of genVertex matched (up to totatl threshold)",
                                       15, -1.5, 13.5, "Number of genVertex matched", "Entries",
                                       true, "TMnumGenMatched");
  h_TM_numGenMatchedAll = HistogramHelper::defineHistogram("Number of genVertex matched (up to totatl threshold to Store)",
                                          15, -1.5, 13.5, "Number of ALL genVertex matched", "Entries",
                                          true, "TMnumGenMatchedAll");
  h_TM_class = HistogramHelper::defineHistogram("Truth-Match type for reconstructed vertices",
                               6, -0.5, 5.5, "Category", "Entries",
                               true, "TMclass");
  h_TM_class->GetXaxis()->SetBinLabel(1, "Match");
  h_TM_class->GetXaxis()->SetBinLabel(2, "Merge");
  h_TM_class->GetXaxis()->SetBinLabel(3, "Split");
  h_TM_class->GetXaxis()->SetBinLabel(4, "Fake");
  h_TM_class->GetXaxis()->SetBinLabel(5, "Others");
  h_TM_class->GetXaxis()->SetBinLabel(6, "ERRORS");

  h_TM_fakeRelWeight = HistogramHelper::defineHistogram("Relative fake component in reconstructed vertices",
                                       100, 0, 1., "Relative FAKE Weight", "Entries",
                                       true, "TMfakeRelWeight");

  // Old matching (distance gen-reco)
  for (int nv=0; nv<VtxZM_NMatch; nv++) {
    h_MCVrt_RecoSimu_VtxM[nv] = 0;
  }
  h_MCVrt_Num_SimuVtxGood = HistogramHelper::defineHistogram("h_MCVrt_Num_SimuVtxGood", 20, -0.5, 19.5, "Number of generated good Vertices", "Events", true, "MDnumGoodVtxSim");
  h_MCVrt_Num_RecoVtxGood = HistogramHelper::defineHistogram("h_MCVrt_Num_RecoVtxGood", 20, -0.5, 19.5, "Number of reconstructed good Vertices", "Events", true, "MDnumGoodVtxReco");
  h_MCVrt_RecoSimu_VtxM[VtxZM_Match] = HistogramHelper::defineHistogram("h_MCVrt_RecoSimu_VtxZM_Match", 20, -0.5, 19.5, "Number of generated Vertices", "Fraction of Vertices: Match", true, "MDMatch");
  h_MCVrt_RecoSimu_VtxM[VtxZM_Split] = HistogramHelper::defineHistogram("h_MCVrt_RecoSimu_VtxZM_Split", 20, -0.5, 19.5, "Number of generated Vertices", "Fraction of Vertices: Split", true, "MDSplit");
  h_MCVrt_RecoSimu_VtxM[VtxZM_Merge] = HistogramHelper::defineHistogram("h_MCVrt_RecoSimu_VtxZM_Merge", 20, -0.5, 19.5, "Number of generated Vertices", "Fraction of Vertices: Merge", true, "MDMerge");
  h_MCVrt_RecoSimu_VtxM[VtxZM_Ineff] = HistogramHelper::defineHistogram("h_MCVrt_RecoSimu_VtxZM_Ineff", 20, -0.5, 19.5, "Number of generated Vertices", "Fraction of Vertices: Ineff", true, "MDIneff");
  h_MCVrt_RecoSimu_VtxM[VtxZM_Fake] = HistogramHelper::defineHistogram("h_MCVrt_RecoSimu_VtxZM_Fake", 20, -0.5, 19.5, "Number of generated Vertices", "Fraction of Vertices: Fake", true, "MDFake");
  h_MCVrt_RecoSimu_VtxM[VtxZM_Match_2] = HistogramHelper::defineHistogram("h_MCVrt_RecoSimu_VtxZM_Match_2", 20, -0.5, 19.5, "Number of generated Vertices", "Fraction of Vertices: Match 2 Vtx", true, "MDMatch2");
  h_MCVrt_RecoSimu_VtxM[VtxZM_Others] = HistogramHelper::defineHistogram("h_MCVrt_RecoSimu_VtxZM_Others", 20, -0.5, 19.5, "Number of generated Vertices", "Fraction of Vertices: Others", true, "MDOthers");

  // Distances for reco-true matching
  h_MCVrt_RecoSimu_VtxZM_Dist_Match = HistogramHelper::defineHistogram("h_MCVrt_RecoSimu_VtxZM_Dist_Match", 100, 0, cfgMaxDistance, "Distance reco'ed-generated vertex", "Matched vertices",
                                      true, "MDdistMatch");
  h_MCVrt_RecoSimu_VtxZM_Dist_Split_gen_1 = HistogramHelper::defineHistogram("h_MCVrt_RecoSimu_VtxZM_Dist_Split_gen_1", 100, 0, cfgMaxDistance,
      "Distance nearest_reco-generated vertex", "Split vertices",
      true, "MDdistSplit1");
  h_MCVrt_RecoSimu_VtxZM_Dist_Split_gen_2 = HistogramHelper::defineHistogram("h_MCVrt_RecoSimu_VtxZM_Dist_Split_gen_2", 100, 0, cfgMaxDistance,
      "Distance second_reco-generated vertex", "Split vertices",
      true, "MDdistSplit2");
  h_MCVrt_RecoSimu_VtxZM_Dist_Split_rec = HistogramHelper::defineHistogram("h_MCVrt_RecoSimu_VtxZM_Dist_Split_rec", 100, 0, cfgMaxDistance, "Distance reco-reco vertex", "Split vertices",
                                          true, "MDdistSplitRec");
  h_MCVrt_RecoSimu_VtxZM_Dist_Merge_firstGen = HistogramHelper::defineHistogram("h_MCVrt_RecoSimu_VtxZM_Dist_Merge_firstGen", 100, 0, cfgMaxDistance,
      "Distance nearest reco'ed-generated vertex", "Merge vertices",
      true, "MDdistMerge1");
  h_MCVrt_RecoSimu_VtxZM_Dist_Merge_secondGen = HistogramHelper::defineHistogram("h_MCVrt_RecoSimu_VtxZM_Dist_Merge_secondGen", 100, 0, cfgMaxDistance,
      "Distance reco'ed-second_generated vertex", "Merge vertices",
      true, "MDdistMerge2");
  h_MCVrt_RecoSimu_VtxZM_Dist_Fake = HistogramHelper::defineHistogram("h_MCVrt_RecoSimu_VtxZM_Dist_Fake", 100, 0, cfgMaxDistance*10, "Distance nearest generated vertex", "Fake vertices",
                                     true, "MDdistFake");
  h_MCVrt_RecoSimu_VtxZM_Dist_Ineff = HistogramHelper::defineHistogram("h_MCVrt_RecoSimu_VtxZM_Dist_Ineff", 100, 0, cfgMaxDistance*10, "Distance nearest reco'ed vertex", "Ineff vertices",
                                      true, "MDdistIneff");


  //Truth matching studies

  h_TM_splitDeltaZ = HistogramHelper::defineHistogram("h_TM_splitDeltaZ", 100, -10, 10,
                                     "#Delta Z between split vertices (w.r.t. Higher #sum P_{T}^2)","Entries",
                                     true, "TM_splitDeltaZ");

  h_TM_splitDeltaZSig = HistogramHelper::defineHistogram("h_TM_splitDeltaZSig", 100, -25, 25,
                                        "#Delta Z / #sigma(#Delta Z) between split vertices (w.r.t. Higher #sum P_{T}^2)", "Entries",
                                        true, "TM_splitDeltaZSig");

  h_TM_class_hs = HistogramHelper::defineHistogram("Truth-Match type for hard-scattering interaction",
                                  6, -0.5, 5.5, "Category", "Entries",
                                  true, "TMclass_hs");
  h_TM_class_hs->GetXaxis()->SetBinLabel(1, "Events");
  h_TM_class_hs->GetXaxis()->SetBinLabel(2, "Match");
  h_TM_class_hs->GetXaxis()->SetBinLabel(3, "Merge");
  h_TM_class_hs->GetXaxis()->SetBinLabel(4, "Split");
  h_TM_class_hs->GetXaxis()->SetBinLabel(5, "ID efficiency");
  h_TM_class_hs->GetXaxis()->SetBinLabel(6, "ERRORS");

  h_TM_class_pu = HistogramHelper::defineHistogram("Truth-Match type for pile-up interactions",
                                  5, -0.5, 4.5, "Category", "Entries",
                                  true, "TMclass_pu");
  h_TM_class_pu->GetXaxis()->SetBinLabel(1, "PU interactions");
  h_TM_class_pu->GetXaxis()->SetBinLabel(2, "Match");
  h_TM_class_pu->GetXaxis()->SetBinLabel(3, "Merge");
  h_TM_class_pu->GetXaxis()->SetBinLabel(4, "Split");
  h_TM_class_pu->GetXaxis()->SetBinLabel(5, "ERRORS");

  //Init counters
  m_TotEvents = 0;
  m_TotRawEvents = 0;
  m_VtxRec = 0;
  m_VtxRecPri = 0;
  m_VtxRecSec = 0;
  m_VtxRecPU = 0;
  m_VtxSim = 0;
  m_VtxGoodSim = 0;
  m_NoVtxButTracks = 0;
  m_VtxTagRec = 0;
  m_VtxTagZSel = 0;
  m_VtxTagPurPos = 0;
  m_VtxAtLease2Trk = 0;
  m_VtxTagRecSim1 = 0;
  m_VtxOnlyTagRecSim1=0;
  m_VtxTagZSel_wrtAnyGenVtx=0;
  m_VtxIDEfficiency=0;

  m_AvgVtxRec = 0.0;
  m_AvgVtxSim = 0.0;
  m_VtxTagRecoEff = 0.0;
  m_VtxTagRecoEffErr = 0.0;
  m_VtxRecoEffPU = 0.0;
  m_VtxRecoEffPUErr = 0.0;
  m_VtxTagZSelEff = 0.0;
  m_VtxTagZSelEffErr = 0.0;
  m_VtxTagZEff = 0.0;
  m_VtxTagZEffErr = 0.0;
  m_VtxTagPurSelEff = 0.0;
  m_VtxTagPurSelEffErr = 0.0;
  m_VtxTagZSelGenAnyEff=0.0;
  m_VtxTagZSelGenAnyEffErr=0.0;
  
  for (unsigned int tveff = 0; tveff < TagVtxNumTypes; tveff++) {
    m_VtxTagThreshold[tveff] = 0;
    m_VtxTagTSelEff[tveff].SetValue(0.0,0.0);
  }

  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    m_vtx_total[*nTrkCut] = 0;
    m_vtx_fake[*nTrkCut] = 0;
    m_vtx_split[*nTrkCut] = 0;
    m_vtx_real[*nTrkCut] = 0;
    m_efficiency[*nTrkCut] = 0;
  }
  m_ngenint = 0;

  //Display settings
  cout << "[AnaVtxTree] Quality Vertex version: "
       << GetQualityVertexVersion(0) <<"." << GetQualityVertexVersion(1) <<"." << GetQualityVertexVersion(2) << std::endl;
  #ifdef PLOT_MC_TRUTH
  cout << "[AnaVtxTree] Working in MC mode: truth quantities calculated." << endl;
  cout << "[AnaVtxTree] Selection efficiency max Z distance: " << VtxEffCriteriaZ << " mm" << std::endl;
  cout << "[AnaVtxTree] Selection efficiency max T distance: " << VtxEffCriteriaT << " (sigma)" << std::endl;
  #endif
  cout << "[AnaVtxTree] Trigger filter:" << triggerName << endl;
  cout << "[AnaVtxTree] actualIntPerXing filter: " << muFilter_min << " to " << muFilter_max << endl;

  //Open output file
  outputNtp = TFile::Open((outputFileName + TString(".root")), "RECREATE");
}

Bool_t AnaVtxTree::Process(Long64_t entry) {

  //Setup trigger D3PD decoding tool
  //------------------------------
  if (trigMetaDataTree!=0 && triggerName != "" && !triggerTool) {
    triggerTool = new D3PD::TrigDecisionToolD3PD(fChain, trigMetaDataTree);
  }

  //Get current entry
  //-------------------------------
  fChain->GetTree()->GetEntry(entry);

  if (triggerTool) {
    triggerTool->GetEntry(entry);
  }
  if (m_TotRawEvents % 1000 == 0) {
    cout << "[AnaVtxTree] Processing event " << m_TotRawEvents << std::endl;
  }

  TStopwatch *timer(0);
  if (do_timing && m_TotRawEvents % 100 == 0) {
    timer = new TStopwatch();
  }

  ++m_TotRawEvents; //before GRL and trigger requirement

  // Apply trigger requirement
  //if (!trigger("L1_MBTS_1"))
  //  return false;

  // Trigger filtering
  if (triggerTool and triggerName != "") {
    if (!(triggerTool->IsPassed(triggerName.Data())) ) {
      return false;  //next event
    }
  }

  //mu-filtering
  if (muFilter_min >= 0) {
    if (ei_actualIntPerXing < muFilter_min) {
      return false;  //skip event
    }
  }
  if (muFilter_max >= 0) {
    if (ei_actualIntPerXing > muFilter_max) {
      return false;  //skip event
    }
  }

  //generator-level filters (ad-hoc)
  /*
  {
    //1. Filter for at least 2 leptons (e, mu) with pT > 10 GeV
    int passDileptonFilter=0;
    for (int i=0; i < mcpart_n; i++) {
      int pdgid = mcpart_type->at(i);
      if ( (abs(pdgid) == 11) || (abs(pdgid) == 13) )
  if (mcpart_pt->at(i) >= 10000.)
    passDileptonFilter++;
      if (passDileptonFilter >= 2)
  break; //enough
    }
    if (passDileptonFilter <2)
      return false; //skip event
  }
  */

  ++m_TotEvents;
  h_actualInt_mcvtxn->Fill(ei_actualIntPerXing, (skip_hs ? mcvtx_n - 1 : mcvtx_n));

  // --- Simulated Primary vertex
  Int_t simVtxPriIndex=0; // ALWAYS it's the first vertex
  //just x-check
  Int_t TagVtxIndex=0; //Tagged vertex, can change depending on quality requirements

  // --- Check if tagged vertex exists and has meaningful errors (dummy vertex always stored!)
  //cout << "[AnaVtxTree] Calling isGoodVertex on line 871" << endl;
  TagVtxIndex = isGoodVertex();

  Int_t NGoodRecoVertices=0;
  for (int nv=0; nv < (vxnbc_n-1); nv++) //exclude dummy
    if (isGoodVertex(nv) == nv) {
      NGoodRecoVertices++;
    }
  // --- Fill simple ones..
  m_VtxRec += (NGoodRecoVertices);

  #ifdef PLOT_MC_TRUTH
  //As first thing, run truth matching algorithms
  // New Truth-Matching
  //===================
  // - Match non beam-constrained vertices
  if (do_timing && m_TotRawEvents % 100 == 0) {
    timer->Start();
  }
  if (debugLevel > 11 or debugLevel == -16) {
    cout << "Debug non Beam-Constraint vertex matching." << endl
         << " EVENT " << ei_EventNumber << ", RUN = " << ei_RunNumber << endl;
    TmVtxNBc.SetDebug(3);
  } else {
    TmVtxNBc.SetDebug(1); //enable histograms anyway
  }
  TmVtxNBc.InitDebugHistograms("_NBC_"); //init debug histogram with a prefix
  TmVtxNBc.SetTrackMatchProbability(vtxTMTrackMatchProb);
  TmVtxNBc.SetVtxMatchWeight(vtxTMWThreshold);
  TmVtxNBc.SetVtxStoreWeight(vtxTMWThresholdStore);
  TmVtxNBc.SetGenVertexRequirementVersion(1); //in-time pile-up vertices
  //TmVtxNBc.SetGenVertexRequirementVersion(2);
  //set inputs
  TmVtxNBc.SetRecoVtxInfo(vxnbc_n,
                          vxnbc_trk_n, vxnbc_trk_weight, vxnbc_trk_index,
                          vxnbc_x, vxnbc_y, vxnbc_z,
                          vxnbc_cov_x, vxnbc_cov_y, vxnbc_cov_z);
  TmVtxNBc.SetRecoTrkInfo(trk_n,
                          trk_mcpart_probability, trk_mcpart_index,
                          trk_pt, trk_eta, trk_d0_wrtPV, trk_z0_wrtPV);
  TmVtxNBc.SetGenPartInfo(mcpart_n,
                          mcpart_mcevt_index, mcpart_mcprodvtx_index, mcpart_type,
                          mcpart_pt, mcpart_eta, mcpart_barcode, mcpart_status);
  TmVtxNBc.SetGenVtxInfo(mcvtx_n,mcvtx_mcevt_index,
                         mcvtx_x, mcvtx_y, mcvtx_z);
  TmVtxNBc.SetGenEventsInfo(mcevt_n, mcevt_pileUpType, mcevt_nparticle);

  //Sanity check
  //cout << "[AnaVtxTree] Calling isGoodGenVertex on line 959 " << endl;
  /*for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    Int_t GoodVertices = 0;
    for (int gvtx=0; gvtx < mcvtx_n; gvtx++){
      //cout << "[AnaVtxTree] Sanity check for track cut " << *nTrkCut << endl; 
      if (TmVtxNBc.isGoodGenVertex(gvtx, ei_RunNumber, ei_EventNumber, *nTrkCut)) {
        GoodVertices++;
     }
    }
    h_GoodVertices_actualInt[*nTrkCut]->Fill(ei_actualIntPerXing,GoodVertices);
    h_GoodVertices_actualInt_2D[*nTrkCut]->Fill(ei_actualIntPerXing,GoodVertices);
    h_GoodVertices_mcvtxn[*nTrkCut]->Fill(mcvtx_n,GoodVertices);
    h_GoodVertices_mcvtxn_2D[*nTrkCut]->Fill(mcvtx_n,GoodVertices);
    Int_t ngenint = 0;
    for (int ivx=0; ivx < mcvtx_n; ++ivx) {
      if (ivx == 0 && skip_hs) {
        //std::cout << "Skipping hard scatter" << std::endl;
        continue;  // In new MC, hard scatter is empty.
      }
      Int_t vtxtype = mcevt_pileUpType->at(mcvtx_mcevt_index->at(ivx));
      if ((vtxtype == 0) || (vtxtype == 1)) {
        ngenint++;
    }
  }
    h_GoodVertices_NGenInt[*nTrkCut]->Fill(ngenint,GoodVertices);
    h_GoodVertices_NGenInt_2D[*nTrkCut]->Fill(ngenint,GoodVertices);
  }*/

  //do matching
  //cout << "[AnaVtxTree] Calling MatchVertices " << endl;
  TmVtxNBc.MatchVertices(ei_RunNumber, ei_EventNumber);

  if (do_timing && m_TotRawEvents % 100 == 0) {
    cout << "[AnaVtxTree] TIMING : Matching = " << timer->RealTime() << endl;
  }

  //Now fill general truth-based properties
  //==========================
  //cout << "[AnaVtxTree] Calling isGoodGenVertex on line 943 " << endl;
  m_VtxSim += (mcvtx_n + (skip_hs ? -1 : 0));
  Int_t numSimVtxGood(0);
  for (int gvx=0; gvx < mcvtx_n; gvx++) {
    if (gvx == 0 && skip_hs) {
      continue;
    }
    if (TmVtxNBc.isGoodGenVertex(gvx, ei_RunNumber, ei_EventNumber)) {
      numSimVtxGood++;
    }
  }
  m_VtxGoodSim += numSimVtxGood;

  h_truth_z->Fill((*mcvtx_z)[simVtxPriIndex]); //Simulated vertex w/ highest Pt^2
  h_truth_x->Fill((*mcvtx_x)[simVtxPriIndex]); //Simulated vertex w/ highest Pt^2
  h_truth_y->Fill((*mcvtx_y)[simVtxPriIndex]); //Simulated vertex w/ highest Pt^2
  m_VtxRecoEffPU += Double_t(m_VtxRec) / m_VtxSim;

  // Now fill futher studies using this truth-matching
  //cout << "[AnaVtxTree] Calling isGoodGenVertex on line 963 " << endl;
  if (do_timing && m_TotRawEvents % 100 == 0) {
    timer->Start();
  }
  // Study high-pT and pile-up interactions separately
  Int_t numSimVtxGood_HS(0);
  Int_t numSimVtxGood_PU(0);
  for (int gvx=0; gvx < mcvtx_n; gvx++){ //I added this braket, there was nothing here before.
    if (TmVtxNBc.isGoodGenVertex(gvx, ei_RunNumber, ei_EventNumber)) {
      if (mcevt_pileUpType->at(mcvtx_mcevt_index->at(gvx)) == 0) {
        numSimVtxGood_HS++;
      } else {
        numSimVtxGood_PU++;
      }
    }
    }
  h_TM_class_hs->Fill(0.0, Double_t(numSimVtxGood_HS)); //total number of "good" events
  h_TM_class_pu->Fill(0.0, Double_t(numSimVtxGood_PU)); // total number of "good" pile-up interactions
  
  bool hs_used=false; //mostly for debug
  for (vector<VertexTruthMatch>::iterator itRM = TmVtxNBc.matchedVtx.begin(); itRM != TmVtxNBc.matchedVtx.end(); ++itRM) {
    if (itRM->m_type == VertexTruthMatch::VtxTM_Match) {
      if (itRM->m_matchList[0].first == simVtxPriIndex) {
        //there will be only one (if two, the latter is marked as split!)
        h_TM_class_hs->Fill(1); //matched
        if (hs_used) {
          cerr << "WARNING: Why two vertices 'Match' to hard-scattering interaction?" << endl;
        }
        hs_used = true;
        //check if tagged as PV (both for Match and Merge)
        if (itRM->m_recoVtx == TagVtxIndex) {
          h_TM_class_hs->Fill(4);
        }
      } else {
        // -> matched to a pile-up interaction
        h_TM_class_pu->Fill(1);
      }
    }
    if (itRM->m_type == VertexTruthMatch::VtxTM_Merge) {
      if (itRM->m_matchList[0].first == simVtxPriIndex) {
        h_TM_class_hs->Fill(2); //merge
        if (hs_used) {
          cerr << "WARNING: Why a 'Merge' vertex match to hard-scattering interaction after a good Match? (should be split)" << endl;
        }
        //check if tagged as PV (both for Match and Merge)
        if (itRM->m_recoVtx == TagVtxIndex) {
          h_TM_class_hs->Fill(4);
        }
      } else {
        // -> matched to a pile-up interaction
        h_TM_class_pu->Fill(2);
      }
    }
    if (itRM->m_type == VertexTruthMatch::VtxTM_Split) {
      if (itRM->m_matchList[0].first == simVtxPriIndex) {
        h_TM_class_hs->Fill(3);
      } else {
        // -> matched to a pile-up interaction
        h_TM_class_pu->Fill(3);
      }
    }
  } // end loop over matched reconstructed vertices
  if (do_timing && m_TotRawEvents % 100 == 0) {
    cout << "[AnaVtxTree] TIMING : Simone's TM histos = " << timer->RealTime() << endl;
  }
  /********************************************************************************/
  //David's stuff
  /********************************************************************************/

  if (do_timing && m_TotRawEvents % 100 == 0) {
    timer->Start();
  }

  //David's temporary filter
  // -- Number of generated interactions
  Int_t current_NGenInt = 0;
  for (int ivx=0; ivx < mcvtx_n; ++ivx) {
    if (ivx == 0 && skip_hs) {
      //std::cout << "Skipping hard scatter" << std::endl;
      continue;  // In new MC, hard scatter is empty.
    }

    //if (TmVtxNBc.isGenInteraction(ivx, ei_RunNumber, ei_EventNumber)) {
    Int_t current_type = mcevt_pileUpType->at(mcvtx_mcevt_index->at(ivx));
    if ((current_type == 0) || (current_type == 1)) {
      current_NGenInt++;
    }
  }

  // -- Total number of events
  m_ngenint += current_NGenInt;
  h_NTrig_NGenInt->Fill(current_NGenInt);
  h_actualInt_NGenInt->Fill(ei_actualIntPerXing,current_NGenInt);

  if (do_timing && m_TotRawEvents % 100 == 0) {
    cout << "[AnaVtxTree] TIMING : David 1 = " << timer->RealTime() << endl;
    timer->Start();
  }

  // -- Classify vertices as real, split, fake.
  // -- Fake: if not at least NTrkCut vertices come from a single generated interaction.
  // -- Split: non-fake, and vertex shares a dominant contribution with a previous vertex.

  std::map<Int_t, std::map<Int_t, Int_t> >  vtx_tm_result; // Map is NTrkCut: Vertex index: truth match category

  std::map<Int_t, Int_t> current_NVtxRecon;

  // Initialize vtx_tm_result to "TmUnknown" for each vertex, and also count reconstructed vertices.
  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    for (int nv=0; nv < vxnbc_n-1; nv++) {
      vtx_tm_result[*nTrkCut][nv] = TmUnknown;
      if ((*vxnbc_nTracks)[nv] >= *nTrkCut) {
        current_NVtxRecon[*nTrkCut]++;
        m_vtx_total[*nTrkCut]++;
      }
    }
    h_NTrig_NVtxRecon[*nTrkCut]->Fill(current_NVtxRecon[*nTrkCut]);
  }

  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    if (current_NVtxRecon[*nTrkCut] > 0) {
      h_all_events_NGenInt[*nTrkCut]->Fill(current_NGenInt);
      h_all_events_NVtxRecon[*nTrkCut]->Fill(current_NVtxRecon[*nTrkCut]);
    }
  }

  // Reco-truth index map
  std::map<Int_t, Int_t> reco_truth_indexmap;
  for (int nv=0; nv < vxnbc_n-1; nv++) {
    reco_truth_indexmap[nv] = -2;
    //Pre-selection: make sure the reconstructed vertex was processed in the truth matcher.
    vector<VertexTruthMatch>::iterator itRM;
    for (itRM = TmVtxNBc.matchedVtx.begin(); itRM != TmVtxNBc.matchedVtx.end(); ++itRM) {
      if (itRM->m_recoVtx == nv) {
        break;
      }
    }

    if (itRM == TmVtxNBc.matchedVtx.end()) {
      cerr << "ERROR: Good NBC vertex with no TM info at all.. needs debugging. ";
      continue; //no good match found
    }

    //Set reco_truth_indexmap with reconstructed index -> first non-fake generated vertex.
    for (vector<std::pair<int, float> >::iterator it = itRM->m_matchList.begin(); it != itRM->m_matchList.end(); ++it) {
      if ((*it).first != -1) {
        if (skip_hs && (*it).first == 0) {
          continue;
        }

        reco_truth_indexmap[nv] = (*it).first;
        break;
      }
    }
  }

  // Next, classifications

  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {

    Int_t current_NFake = 0;
    Int_t current_NSplit = 0;
    Int_t current_NReal = 0;

    for (int nv=0; nv < vxnbc_n-1; nv++) {

      if (((*vxnbc_nTracks)[nv] < *nTrkCut) || ((*vxnbc_cov_x)[nv] == 0) || ((*vxnbc_cov_y)[nv] == 0)) {
        vtx_tm_result[*nTrkCut][nv] = TmReject;
        continue;
      }

      // Also reject based on chi2/ndf
      if (max_chi2ndf > 0. && (*vxnbc_chi2)[nv] / (*vxnbc_ndof)[nv] > max_chi2ndf) {
        vtx_tm_result[*nTrkCut][nv] = TmReject;
        continue;
      }

      h_all_NGenInt[*nTrkCut]->Fill(current_NGenInt);
      h_all_actualInt[*nTrkCut]->Fill(ei_actualIntPerXing);
      h_all_NVtxRecon[*nTrkCut]->Fill(current_NVtxRecon[*nTrkCut]);
      h_NVtxRecon_NGenInt[*nTrkCut]->Fill(current_NGenInt, current_NVtxRecon[*nTrkCut]);
      h_NVtxRecon_actualInt[*nTrkCut]->Fill(ei_actualIntPerXing, current_NVtxRecon[*nTrkCut]);
      //cout << "[AnaVtxTree] INFO: current_NGenInt = " << current_NGenInt<< ", ei_actualIntPerXing = "<< ei_actualIntPerXing<< endl;

      // Get truth match for vertex nv
      vector<VertexTruthMatch>::iterator itRM;
      for (itRM = TmVtxNBc.matchedVtx.begin(); itRM != TmVtxNBc.matchedVtx.end(); ++itRM) {
        if (itRM->m_recoVtx == nv) {
          break;
        }
      }


      // -- CLASSIFICATION 1: FAKE
      bool vtx_is_fake = true;
      Int_t n_fake_tracks = 0;
      Int_t n_genvtx_matched = 0;
      Int_t dominant_genvtx_idx = -2; // If the first genint has >=10% weight, flag it for later pull calculations.
      
      vector<std::pair<int, float> >::iterator gen_vtx = (itRM->m_matchList).begin();
      
      for (vector<std::pair<int, float> >::iterator gen_vtx = (itRM->m_matchList).begin(); gen_vtx != (itRM->m_matchList).end(); ++gen_vtx) {
        if (!(&(gen_vtx))) continue;
        if ((*gen_vtx).first == 0 && skip_hs) {
          continue;
        }
        if ((*gen_vtx).first == -1) {
          n_fake_tracks = (*gen_vtx).second * (*vxnbc_nTracks)[nv];
          continue; // Skip the fake contribution
        }
        if (dominant_genvtx_idx == -2) {
          dominant_genvtx_idx = (*gen_vtx).first;
        }
        if ((*gen_vtx).second * (*vxnbc_nTracks)[nv] >= (float)*nTrkCut - 0.01) {
          vtx_is_fake = false;
        }
        n_genvtx_matched++;
      }
      if (vtx_is_fake) {
        vtx_tm_result[*nTrkCut][nv] = TmFake;
        m_vtx_fake[*nTrkCut]++;
        current_NFake++;
        h_fakes_NGenInt[*nTrkCut]->Fill(current_NGenInt);
        //cout << "m_vtx_fake[" << *nTrkCut << "] = " << m_vtx_fake[*nTrkCut] << endl;
        //cout << "current_NFake = " << current_NFake << endl;
        //cout << "current_NGenInt = " << current_NGenInt << endl;
        //cout << "current_NVtxRecon[" << *nTrkCut << "] = " << current_NVtxRecon[*nTrkCut] << endl;
        h_fakes_NVtxRecon[*nTrkCut]->Fill(current_NVtxRecon[*nTrkCut]);

        if (n_genvtx_matched == 0) {
          h_truefakes_NGenInt[*nTrkCut]->Fill(current_NGenInt);
          h_truefakes_NVtxRecon[*nTrkCut]->Fill(current_NVtxRecon[*nTrkCut]);
        }

        h_fakes_x[*nTrkCut]->Fill((*vxnbc_x)[nv]);
        h_fakes_y[*nTrkCut]->Fill((*vxnbc_y)[nv]);
        h_fakes_z[*nTrkCut]->Fill((*vxnbc_z)[nv]);
        h_fakes_cov_xx[*nTrkCut]->Fill((*vxnbc_cov_x)[nv]);
        h_fakes_cov_yy[*nTrkCut]->Fill((*vxnbc_cov_y)[nv]);
        h_fakes_cov_zz[*nTrkCut]->Fill((*vxnbc_cov_z)[nv]);
        h_fakes_cov_xy[*nTrkCut]->Fill((*vxnbc_cov_xy)[nv]);
        h_fakes_cov_xz[*nTrkCut]->Fill((*vxnbc_cov_xz)[nv]);
        h_fakes_cov_yz[*nTrkCut]->Fill((*vxnbc_cov_yz)[nv]);
        h_fakes_chi2ndf[*nTrkCut]->Fill((*vxnbc_chi2)[nv] / (*vxnbc_ndof)[nv]);
        h_fakes_chi2[*nTrkCut]->Fill((*vxnbc_chi2)[nv]);
        h_fakes_ndf[*nTrkCut]->Fill((*vxnbc_ndof)[nv]);
        h_fakes_ntrk[*nTrkCut]->Fill((*vxnbc_nTracks)[nv]);
        h_fakes_sumPt[*nTrkCut]->Fill((*vxnbc_sumPt)[nv]);

        if (current_NGenInt == 1 && *nTrkCut >= 5) {
          // We don't expect many fakes here, so dump out the information
          cout << "[AnaVtxTree] WARNING : Fake found in NGenInt = 1 bin. Information: " << endl;
          cout << "[AnaVtxTree] WARNING : \t\tTracks = " << (*vxnbc_nTracks)[nv] << endl;
          cout << "[AnaVtxTree] WARNING : \t\tmcvtx_n = " << mcvtx_n << endl;
          cout << "[AnaVtxTree] WARNING : \t\tdominant_genvtx_idx = " << dominant_genvtx_idx << endl;
          cout << "[AnaVtxTree] WARNING : \t\tMatched vertices = " << n_genvtx_matched << endl;
          for (vector<std::pair<int, float> >::iterator gen_vtx = (itRM->m_matchList).begin(); gen_vtx != (itRM->m_matchList).end(); ++gen_vtx) {
            cout << "[AnaVtxTree] WARNING : \t\t" << (*gen_vtx).first << " : weight = " << (*gen_vtx).second << " / tracks = " << (*gen_vtx).second * (*vxnbc_nTracks)[nv] << endl;
          }
          cout << "[AnaVtxTree] WARNING : \tEvent total reconstructed vertices = " << vxnbc_n-1 << endl;
          for (int i = 0; i < vxnbc_n-1; i++) {
            cout << "[AnaVtxTree] WARNING : \t\t" << i << " has nTracks = " << (*vxnbc_nTracks)[i] << endl;
          }
          cout << endl;
        }
        continue; // Vertex is fake, move on
      }

      // -- CLASSIFICATION 2: SPLIT
      for (int nv2 = 0; nv2 < nv; nv2++) {
        // Debug: make sure all the prior vertices have been classified
        if (vtx_tm_result[*nTrkCut][nv2] == TmUnknown) {
          cout << "[AnaVtxTree] WARNING : Already-processed vertex is still labelled as TmUnknown! Needs debugging!" << endl;
        }
        if ((vtx_tm_result[*nTrkCut][nv2] != TmReal) || (nv == nv2)) {
          continue;
        } else {
          if (reco_truth_indexmap[nv] == reco_truth_indexmap[nv2]) {
            vtx_tm_result[*nTrkCut][nv] = TmSplit;
            m_vtx_split[*nTrkCut]++;
            current_NSplit++;
            h_splits_NGenInt[*nTrkCut]->Fill(current_NGenInt);
            h_splits_NVtxRecon[*nTrkCut]->Fill(current_NVtxRecon[*nTrkCut]);

            h_splits_x[*nTrkCut]->Fill((*vxnbc_x)[nv]);
            h_splits_y[*nTrkCut]->Fill((*vxnbc_y)[nv]);
            h_splits_z[*nTrkCut]->Fill((*vxnbc_z)[nv]);
            h_splits_cov_xx[*nTrkCut]->Fill((*vxnbc_cov_x)[nv]);
            h_splits_cov_yy[*nTrkCut]->Fill((*vxnbc_cov_y)[nv]);
            h_splits_cov_zz[*nTrkCut]->Fill((*vxnbc_cov_z)[nv]);
            h_splits_cov_xy[*nTrkCut]->Fill((*vxnbc_cov_xy)[nv]);
            h_splits_cov_xz[*nTrkCut]->Fill((*vxnbc_cov_xz)[nv]);
            h_splits_cov_yz[*nTrkCut]->Fill((*vxnbc_cov_yz)[nv]);
            h_splits_chi2ndf[*nTrkCut]->Fill((*vxnbc_chi2)[nv] / (*vxnbc_ndof)[nv]);
            h_splits_chi2[*nTrkCut]->Fill((*vxnbc_chi2)[nv]);
            h_splits_ndf[*nTrkCut]->Fill((*vxnbc_ndof)[nv]);
            h_splits_ntrk[*nTrkCut]->Fill((*vxnbc_nTracks)[nv]);
            h_splits_sumPt[*nTrkCut]->Fill((*vxnbc_sumPt)[nv]);

            Float_t dz = (*vxnbc_z)[nv2] - (*vxnbc_z)[nv];
            Float_t err2 = (*vxnbc_cov_z)[nv2] + (*vxnbc_cov_z)[nv];
            Float_t dz_sig = dz / TMath::Sqrt(err2);
            h_split_dz_NGenInt[*nTrkCut]->Fill(dz, current_NGenInt);
            h_split_dzsig_NGenInt[*nTrkCut]->Fill(dz_sig, current_NGenInt);
            h_split_dz_NVtxRecon[*nTrkCut]->Fill(dz, current_NVtxRecon[*nTrkCut]);
            h_split_dzsig_NVtxRecon[*nTrkCut]->Fill(dz_sig, current_NVtxRecon[*nTrkCut]);
            h_split_dz_dzsig[*nTrkCut]->Fill(dz, dz_sig);

            break;
          }
        }
      }
      if (vtx_tm_result[*nTrkCut][nv] == TmSplit) {
        continue;
      }

      // -- ELSE: REAL
      vtx_tm_result[*nTrkCut][nv] = TmReal;
      m_vtx_real[*nTrkCut]++;
      current_NReal++;
      h_real_NGenInt[*nTrkCut]->Fill(current_NGenInt);
      h_real_NVtxRecon[*nTrkCut]->Fill(current_NVtxRecon[*nTrkCut]);

      h_real_x[*nTrkCut]->Fill((*vxnbc_x)[nv]);
      h_real_y[*nTrkCut]->Fill((*vxnbc_y)[nv]);
      h_real_z[*nTrkCut]->Fill((*vxnbc_z)[nv]);
      h_real_cov_xx[*nTrkCut]->Fill((*vxnbc_cov_x)[nv]);
      h_real_cov_yy[*nTrkCut]->Fill((*vxnbc_cov_y)[nv]);
      h_real_cov_zz[*nTrkCut]->Fill((*vxnbc_cov_z)[nv]);
      h_real_cov_xy[*nTrkCut]->Fill((*vxnbc_cov_xy)[nv]);
      h_real_cov_xz[*nTrkCut]->Fill((*vxnbc_cov_xz)[nv]);
      h_real_cov_yz[*nTrkCut]->Fill((*vxnbc_cov_yz)[nv]);
      h_real_chi2ndf[*nTrkCut]->Fill((*vxnbc_chi2)[nv] / (*vxnbc_ndof)[nv]);
      h_real_chi2[*nTrkCut]->Fill((*vxnbc_chi2)[nv]);
      h_real_ndf[*nTrkCut]->Fill((*vxnbc_ndof)[nv]);
      h_real_ntrk[*nTrkCut]->Fill((*vxnbc_nTracks)[nv]);
      h_real_sumPt[*nTrkCut]->Fill((*vxnbc_sumPt)[nv]);

    } //Finished clasifications
    h_NReal_NVtxRecon[*nTrkCut]->Fill(current_NVtxRecon[*nTrkCut],current_NReal);
    h_NReal_actualInt[*nTrkCut]->Fill(ei_actualIntPerXing,current_NReal);
    h_NReal_NGenInt[*nTrkCut]->Fill(current_NGenInt,current_NReal);

    if (current_NVtxRecon[*nTrkCut] > 0 && current_NReal == 0 && current_NFake > 0) {
      h_fake_events_NGenInt[*nTrkCut]->Fill(current_NGenInt);
      h_fake_events_NVtxRecon[*nTrkCut]->Fill(current_NVtxRecon[*nTrkCut]);
    }

    if (current_NVtxRecon[*nTrkCut] > 0) {
      h_event_contains_any_NGenInt[*nTrkCut]->Fill(current_NGenInt);
      h_event_contains_any_NVtxRecon[*nTrkCut]->Fill(current_NVtxRecon[*nTrkCut]);
    }

    if (current_NFake > 0) {
      h_event_contains_fake_NGenInt[*nTrkCut]->Fill(current_NGenInt);
      h_event_contains_fake_NVtxRecon[*nTrkCut]->Fill(current_NVtxRecon[*nTrkCut]);
    }

    // Loop over pairs of real vertices to look at real vtx separations
    for (int nv=0; nv < vxnbc_n-1; nv++) {
      if ((*vxnbc_nTracks)[nv] < *nTrkCut) {
        continue;
      }
      for (int nv2=0; nv2 < nv; nv2++) {
        if ((*vxnbc_nTracks)[nv2] < *nTrkCut) {
          continue;
        }

        Float_t dz = (*vxnbc_z)[nv2] - (*vxnbc_z)[nv];
        Float_t err2 = (*vxnbc_cov_z)[nv2] + (*vxnbc_cov_z)[nv];
        Float_t dz_sig = dz / TMath::Sqrt(err2);

        h_all_dz_NGenInt[*nTrkCut]->Fill(dz, current_NGenInt);
        h_all_dz_NVtxRecon[*nTrkCut]->Fill(dz, current_NVtxRecon[*nTrkCut]);

        if ((vtx_tm_result[*nTrkCut][nv] == TmReal) && (vtx_tm_result[*nTrkCut][nv2] == TmReal)) {
          h_real_dz_NGenInt[*nTrkCut]->Fill(dz, current_NGenInt);
          h_real_dzsig_NGenInt[*nTrkCut]->Fill(dz_sig, current_NGenInt);
          h_real_dz_NVtxRecon[*nTrkCut]->Fill(dz, current_NVtxRecon[*nTrkCut]);
          h_real_dzsig_NVtxRecon[*nTrkCut]->Fill(dz_sig, current_NVtxRecon[*nTrkCut]);
          h_real_dz_dzsig[*nTrkCut]->Fill(dz, dz_sig);
        }
      }
    }
  } // nTrkCut loop

  if (TagVtxIndex >= 0) {
    h_privtx_z_mu->Fill((*vxnbc_z)[TagVtxIndex], ei_actualIntPerXing);
    h_privtx_z_NGenInt->Fill((*vxnbc_z)[TagVtxIndex], current_NGenInt);
  }


  /********************************************************************************/
  //End David's stuff
  /********************************************************************************/

  #endif

  #ifdef PLOT_MC_TRUTH
  {
    //P.V. Identification efficiency
    vector<VertexTruthMatch>::iterator itRM;
    for (itRM = TmVtxNBc.matchedVtx.begin(); itRM != TmVtxNBc.matchedVtx.end(); ++itRM) {
      if (itRM->m_recoVtx == TagVtxIndex) {
        break;
      }
    }
    if (itRM != TmVtxNBc.matchedVtx.end()) {
      // Found generated vertex matched to the tagged reconstructed vertex
      if (itRM->m_matchList.size() > 0)
        if (itRM->m_matchList[0].first == simVtxPriIndex) {
          m_VtxIDEfficiency++;  //matched to the right vertex
        }
    }
  }
  #endif

  // --- Pile-up vertex variables
  //    for (vector<int>::iterator vt = vxnbc_type->begin(); vt != vxnbc_type->end(); ++vt)
  if (TagVtxIndex >= 0) {
    for (int vt = 0; vt < vxnbc_n - 1; vt++) {
      if (isGoodVertex(vt) != vt) {
        continue;
      }
      if ((*vxnbc_type)[vt] == Trk::PriVtx) {
        ++m_VtxRecPri;
      } else if ((*vxnbc_type)[vt] == Trk::SecVtx) {
        ++m_VtxRecSec;
      } else if ((*vxnbc_type)[vt] == Trk::PileUp) {
        ++m_VtxRecPU;
      }
    }
  } else { //not good vertex found
    if ((*vxnbc_nTracks)[0] > 0) {
      ++m_NoVtxButTracks;
    }
  }
  if (do_timing && m_TotRawEvents % 100 == 0) {
    cout << "[AnaVtxTree] TIMING : Simone 2 = " << timer->RealTime() << endl;
    timer->Start();
  }

  if (do_timing && m_TotRawEvents % 100 == 0) {
    cout << "[AnaVtxTree] TIMING : Final = " << timer->RealTime() << endl;
  }

  return kTRUE;

}


void AnaVtxTree::SlaveTerminate() {

/*  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    h_NTrig_NVtxRecon[*nTrkCut]
  }
*/
  cout << "[AnaVtxTree] Post-processing of results" << endl;

  cout << "[AnaVtxTree] INFO : Generated interactions: " << m_ngenint << endl;
  cout << "[AnaVtxTree] INFO : Vertex classification counts: " << endl;
  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {

    cout << "[AnaVtxTree] INFO : Fakes[" << *nTrkCut << "]\t= " << m_vtx_fake[*nTrkCut] << endl;
    cout << "[AnaVtxTree] INFO : Splits[" << *nTrkCut << "]\t= " << m_vtx_split[*nTrkCut] << endl;
    cout << "[AnaVtxTree] INFO : Real[" << *nTrkCut << "]\t= " << m_vtx_real[*nTrkCut] << endl;

  }

  cout << endl << "[AnaVtxTree] INFO : Vertex efficiencies: " << endl;
  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    cout << "[AnaVtxTree] INFO : Fakes[" << *nTrkCut << "] / NGenInt = " << m_vtx_fake[*nTrkCut] / m_ngenint << endl;
  }
  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    cout << "[AnaVtxTree] INFO : Splits[" << *nTrkCut << "] / NGenInt = " << m_vtx_split[*nTrkCut] / m_ngenint << endl;
  }
  for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    cout << "[AnaVtxTree] INFO : Real[" << *nTrkCut << "] / NGenInt = " << m_vtx_real[*nTrkCut] / m_ngenint << endl;
  }

  // --- Evaluate efficiency-related
  computeStats();

  #ifdef PLOT_MC_TRUTH
  computeMCOnlyStats();
  #endif


  // --- Print/save results
  cout << "Printing results to: " << (outputFileName + TString(".txt")).Data() << std::endl;

  ETypes::SetGlobalSignificantDigits(2);
  ETypes::SetGlobalNotation(ETypes::kFixed);
  ETypes::SetGlobalPMString(" +/- ");

  cout<<"[AnaVtxTree] RESULTS:" << endl;
  cout<<"---------- Efficiencies and related ----------"<<std::endl;
  cout<<"General:"<<endl;
  cout<<" Number of events: " << m_TotRawEvents << endl;
  cout<<" Number of events w/in GRL: " << m_TotEvents << endl;
  cout<<" Events with >1 tracks: " << m_VtxAtLease2Trk << std::endl;
  cout<<" Average number of reconstructed Vx: " << m_AvgVtxRec << endl;
  cout<<" Average number of simulated Vx: " << m_AvgVtxSim << endl;
  cout<<" Events with no Tagged Vertex but with input Tracks: " << m_NoVtxButTracks << endl;
  cout<<"Efficiency:"<<endl;
  cout<<" Reconstruction Efficiency: "<< m_VtxTagRecoEff << " +/- " << m_VtxTagRecoEffErr << endl;
  cout<<" Z-Selection Efficiency (" << VtxEffCriteriaZ << " mm): " << m_VtxTagZSelEff << " +/- " << m_VtxTagZSelEffErr << endl;
  cout<<" Z-Efficiency (Total): " << m_VtxTagZEff << " +/- " << m_VtxTagZEffErr << endl;
  cout<<" Purity-based Selection Efficiency: " <<  m_VtxTagPurSelEff << " +/- " << m_VtxTagPurSelEffErr << endl;
  cout<<" T-Efficiency:" << endl;
  cout<<"  Correct: " << m_VtxTagTSelEff[TagVtxCorrect] << endl;
  cout<<"  Split(Correct): " << m_VtxTagTSelEff[TagVtxSplitCorrect] << endl;
  cout<<"  Split(Wrong): " << m_VtxTagTSelEff[TagVtxSplitWrong] << endl;
  cout<<"  Wrong: " << m_VtxTagTSelEff[TagVtxWrong] << endl;
  cout<<"  Wrong (Fake): " << m_VtxTagTSelEff[TagVtxWrongFake] << endl;
  cout<<"Fake rates:"<<endl;
  cout<<" Single vertex fake rate: " << m_VtxFakeReco << " +/-" << m_VtxFakeRecoErr << endl;
  cout<<"Tagged Vx resolution (Reco v.s. MC)"<<endl;
  cout<<" sigma_X = " << m_VtxResX << " +/- " << m_VtxResXErr << endl;
  cout<<" sigma_Y = " << m_VtxResY << " +/- " << m_VtxResYErr << endl;
  cout<<" sigma_Z = " << m_VtxResZ << " +/- " << m_VtxResZErr << endl;
  cout<<"Pile-Up:" << endl;
  cout<<" Truth Matching Info (denominator: good generated vertices)" << endl;
  cout<<"  Matched: " << h_TM_class->GetBinContent(VertexTruthMatch::VtxTM_Match+1) << endl;
  cout<<"  Merge: " << h_TM_class->GetBinContent(VertexTruthMatch::VtxTM_Merge+1) << endl;
  cout<<"  Split: " << h_TM_class->GetBinContent(VertexTruthMatch::VtxTM_Split+1) << endl;
  cout<<"  Fake: " << h_TM_class->GetBinContent(VertexTruthMatch::VtxTM_Fake+1) << endl;
  cout<<"  Others: " << h_TM_class->GetBinContent(VertexTruthMatch::VtxTM_Others+1) << endl;
  cout<<"  ERRORS: " << h_TM_class->GetBinContent(VertexTruthMatch::VtxTM_NMatch+1) << endl;
  cout<<" Identification Efficiency: " << float(m_VtxIDEfficiency) / m_VtxTagRec << endl;
  cout<<" Z-Selection Efficiency wrt any gen Vx: " << m_VtxTagZSelGenAnyEff << " +/- " << m_VtxTagZSelGenAnyEffErr << endl;


  //Write results to a file too - waiting to implement a tee-like cout
  outputTxt << "-------- Displaying settings --------------" << std::endl;
  outputTxt << "[AnaVtxTree] Quality Vertex version: "
            << GetQualityVertexVersion(0) <<"." << GetQualityVertexVersion(1) <<"." << GetQualityVertexVersion(2) << std::endl;
  #ifdef PLOT_MC_TRUTH
  outputTxt << "[AnaVtxTree] Working in MC mode: truth quantities calculated." << endl;
  outputTxt << "[AnaVtxTree] Selection efficiency max Z distance: " << VtxEffCriteriaZ << " mm" << std::endl;
  outputTxt << "[AnaVtxTree] Selection efficiency max T distance: " << VtxEffCriteriaT << " (sigma)" << std::endl;
  outputTxt << "[AnaVtxTree] Trigger filter:" << triggerName << endl;
  outputTxt << "[AnaVtxTree] actualIntPerXing filter: " << muFilter_min << " - " << muFilter_max << endl;
  #endif
  outputTxt<<"---------- Efficiencies and related ----------"<<std::endl;
  outputTxt<<"General:"<<endl;
  outputTxt<<" Number of events: " << m_TotRawEvents << endl;
  outputTxt<<" Number of events w/in GRL: " << m_TotEvents << endl;
  outputTxt<<" Events with >1 tracks: " << m_VtxAtLease2Trk << std::endl;
  outputTxt<<" Average number of reconstructed Vx: " << m_AvgVtxRec << endl;
  outputTxt<<" Average number of simulated Vx: " << m_AvgVtxSim << endl;
  outputTxt<<" Events with no Tagged Vertex but with input Tracks: " << m_NoVtxButTracks << endl;
  outputTxt<<"Efficiency:"<<endl;
  outputTxt<<" Reconstruction Efficiency: "<< m_VtxTagRecoEff << " +/- " << m_VtxTagRecoEffErr << endl;
  outputTxt<<" Z-Selection Efficiency (" << VtxEffCriteriaZ << " mm): " << m_VtxTagZSelEff << " +/- " << m_VtxTagZSelEffErr << endl;
  outputTxt<<" Z-Efficiency (Total): " << m_VtxTagZEff << " +/- " << m_VtxTagZEffErr << endl;
  outputTxt<<" Purity-based Selection Efficiency: " <<  m_VtxTagPurSelEff << " +/- " << m_VtxTagPurSelEffErr << endl;
  outputTxt<<" T-Efficiency:" << endl;
  outputTxt<<"  Correct: " << m_VtxTagTSelEff[TagVtxCorrect] << endl;
  outputTxt<<"  Split(Correct): " << m_VtxTagTSelEff[TagVtxSplitCorrect] << endl;
  outputTxt<<"  Split(Wrong): " << m_VtxTagTSelEff[TagVtxSplitWrong] << endl;
  outputTxt<<"  Wrong: " << m_VtxTagTSelEff[TagVtxWrong] << endl;
  outputTxt<<"  Wrong (Fake): " << m_VtxTagTSelEff[TagVtxWrongFake] << endl;
  outputTxt<<"Fake rates:"<<endl;
  outputTxt<<" Single vertex fake rate: " << m_VtxFakeReco << " +/-" << m_VtxFakeRecoErr << endl;
  outputTxt<<"Tagged Vx resolution (Reco v.s. MC)"<<endl;
  outputTxt<<" sigma_X = " << m_VtxResX << " +/- " << m_VtxResXErr << endl;
  outputTxt<<" sigma_Y = " << m_VtxResY << " +/- " << m_VtxResYErr << endl;
  outputTxt<<" sigma_Z = " << m_VtxResZ << " +/- " << m_VtxResZErr << endl;
  outputTxt<<"Pile-Up:"<<endl;
  outputTxt<<" Inclusive Truth Matching Info (denominator: good generated vertices)" << endl;
  outputTxt<<"  Matched: " << h_TM_class->GetBinContent(VertexTruthMatch::VtxTM_Match+1) << endl;
  outputTxt<<"  Merge: " << h_TM_class->GetBinContent(VertexTruthMatch::VtxTM_Merge+1) << endl;
  outputTxt<<"  Split: " << h_TM_class->GetBinContent(VertexTruthMatch::VtxTM_Split+1) << endl;
  outputTxt<<"  Fake: " << h_TM_class->GetBinContent(VertexTruthMatch::VtxTM_Fake+1)  << endl;
  outputTxt<<"  Others: " << h_TM_class->GetBinContent(VertexTruthMatch::VtxTM_Others+1)  << endl;
  outputTxt<<"  ERRORS: " << h_TM_class->GetBinContent(VertexTruthMatch::VtxTM_NMatch+1) << endl;
  outputTxt<<" Identification Efficiency: " << float(m_VtxIDEfficiency) / m_VtxTagRec << endl;
  //compute efficiencies for HS, PU interactions
  double val, err;
  ETypes::EDouble eff;
  outputTxt<<" Hard-scattering Truth Matching Info (denominator: good generated vertices)" << endl;
  calculateEfficiency(val,err, int(h_TM_class_hs->GetBinContent(2)+h_TM_class_hs->GetBinContent(3)), int(h_TM_class_hs->GetBinContent(1)) );
  eff.SetValue(val, err);
  outputTxt<<"  Matched: " << eff << endl;
  calculateEfficiency(val, err, int(h_TM_class_hs->GetBinContent(3)), int(h_TM_class_hs->GetBinContent(1)));
  eff.SetValue(val, err);
  outputTxt<<"  Merge: " << eff << endl;
  calculateEfficiency(val, err, int(h_TM_class_hs->GetBinContent(4)), int(h_TM_class_hs->GetBinContent(1)));
  eff.SetValue(val, err);
  outputTxt<<"  Split: " << eff << endl;
  calculateEfficiency(val, err, int(h_TM_class_hs->GetBinContent(5)), int(h_TM_class_hs->GetBinContent(2) + h_TM_class_hs->GetBinContent(3)));
  eff.SetValue(val, err);
  outputTxt<<"  ID Efficiency: " << eff  << endl;
  outputTxt<<"  ERRORS: " << h_TM_class_hs->GetBinContent(6) / h_TM_class_hs->GetBinContent(1) << endl;

  outputTxt<<" Pile-up Truth Matching Info (denominator: good generated vertices)" << endl;
  calculateEfficiency(val,err, int(h_TM_class_pu->GetBinContent(2)+h_TM_class_pu->GetBinContent(3)), int(h_TM_class_pu->GetBinContent(1)));
  eff.SetValue(val, err);
  outputTxt<<"  Matched: " << eff << endl;
  calculateEfficiency(val, err, int(h_TM_class_pu->GetBinContent(3)), int(h_TM_class_pu->GetBinContent(1)));
  eff.SetValue(val, err);
  outputTxt<<"  Merge: " << eff << endl;
  calculateEfficiency(val, err, int(h_TM_class_pu->GetBinContent(4)), int(h_TM_class_pu->GetBinContent(1)));
  eff.SetValue(val, err);
  outputTxt<<"  Split: " << eff << endl;
  outputTxt<<"  ERRORS: " << h_TM_class_pu->GetBinContent(6) / h_TM_class_pu->GetBinContent(1) << endl;

  outputTxt<<" Average Number of Primary Vx: " << Double_t(m_VtxRecPri)/m_TotEvents << endl;
  outputTxt<<" Average Number of Secondary Vx: " << Double_t(m_VtxRecSec)/m_TotEvents << endl;
  outputTxt<<" Average Number of Pile-Up Vx: " << Double_t(m_VtxRecPU)/m_TotEvents << endl;
  outputTxt<<" Average Number of Other Vx: " << Double_t(m_VtxRec - m_VtxRecPri - m_VtxRecSec - m_VtxRecPU)/m_TotEvents << endl;
  outputTxt<<" Reconstruction Efficiency (PU): "<< m_VtxRecoEffPU << " +/- " << m_VtxRecoEffPUErr << endl;
  outputTxt<<" Z-Selection Efficiency wrt any gen Vx: " << m_VtxTagZSelGenAnyEff << " +/- " << m_VtxTagZSelGenAnyEffErr << endl;


  cout << "Saving histograms to: " << (outputFileName + TString(".root")) << std::endl;
  outputNtp->cd();
  saveHistograms();
  #ifdef PLOT_MC_TRUTH
  TmVtxNBc.WriteDebugHisto();
  #endif

  // --- Close output files
  outputNtp->Close();
  outputTxt.close();

  cout<<"[AnaVtxTree] All Done." << endl;

}

void AnaVtxTree::saveHistograms() {
  //save all histograms in the hist/ subdirectory
  TDirectory *saveDir = gDirectory;
  gDirectory->mkdir("hist"); gDirectory->cd("hist");
  for( std::vector<TH1*>::iterator hist = HistogramHelper::histQueue.begin(); hist != HistogramHelper::histQueue.end(); ++hist ) {
    (*hist)->Write();
  }
 /* for (map<string,pair<string,TH1*> >::iterator nh = storageHistQueue.begin(); nh != storageHistQueue.end(); nh++) {
    (*nh).second.second->Write();
  } */
  gDirectory = saveDir; gDirectory->cd();

  /*//create C header to retrieve Histograms from file
  ofstream outputCHdrHist;
  outputCHdrHist.open((outputFileName + TString(".h")).Data());
  //FIXME: loop over the saved vector of histogram and determine which types has to be stored, modify include accordingly
  outputCHdrHist << "// S. Pagan Griso: This file was automatically generated by AnaVtxTree class.";
  outputCHdrHist << "// P.s. The InputRootFile TFile object will be opened by left open!";
  outputCHdrHist << "#include \"TFile.h\"" << endl;
  outputCHdrHist << "#include \"TH1.h\"" << endl;
  outputCHdrHist << "#include \"TH1F.h\"" << endl;
  outputCHdrHist << "#include \"TH1D.h\"" << endl;
  outputCHdrHist << "#include \"TH2F.h\"" << endl;
  outputCHdrHist << "#include \"TH2D.h\"" << endl;
  outputCHdrHist << "" << endl;
  outputCHdrHist << "TFile *InputRootFile;" << endl;
  outputCHdrHist << "// --- Declare histograms" << endl;
  for (map<string,pair<string, TH1*> >::iterator nh = storageHistQueue.begin(); nh != storageHistQueue.end(); nh++) {
    //create 'TH1F * h_pippo;' statements
    outputCHdrHist << (*nh).second.first.c_str()
                   << " " << (*nh).first.c_str() << ";" << "  //" << (*nh).second.second->GetTitle() << endl;
  }
  outputCHdrHist << endl << endl;
  outputCHdrHist << "// --- Function to retrieve histograms" << endl;
  outputCHdrHist << "int loadHistograms(const char* InputRootFileName=\""
                 << (outputFileName + TString(".root")).Data() << "\")" << endl;
  outputCHdrHist << "{" << endl;
  outputCHdrHist << "   InputRootFile = TFile::Open(InputRootFileName);" << endl;
  cout << endl;
  for (map<string, pair<string, TH1*> >::iterator nh = storageHistQueue.begin(); nh != storageHistQueue.end(); nh++) {
    outputCHdrHist << "   " << (*nh).first.c_str() << " = (" << (*nh).second.first.c_str()
                   << ")InputRootFile->Get(\"hist/"
                   << (*nh).first.c_str() << "\");" << endl;
  }
  outputCHdrHist << endl;
  outputCHdrHist << "   return 0;" << endl;
  outputCHdrHist << "}" << endl;*/
}


void AnaVtxTree::computeStats() {
  Double_t TotEvents_denominator = m_TotEvents;

  m_AvgVtxRec = Double_t(m_VtxRec) / TotEvents_denominator;

  calculateEfficiency(m_VtxTagRecoEff,m_VtxTagRecoEffErr,m_VtxTagRec, int(TotEvents_denominator) );
  m_VtxTagRecoEff = Double_t(m_VtxTagRec) / TotEvents_denominator;
  std::cout << "In computeStats: " << m_VtxTagRecoEff << ", " << m_VtxTagRec << ", " << TotEvents_denominator << std::endl;
  m_VtxTagRecoEffErr = TMath::Sqrt(Double_t(m_VtxTagRec))/TotEvents_denominator;

  return;
}

void AnaVtxTree::computeMCOnlyStats() {
  //Compute efficiency values
  Double_t TotEvents_denominator = m_TotEvents;

  m_AvgVtxSim = Double_t(m_VtxSim) / TotEvents_denominator;
  calculateEfficiency(m_VtxRecoEffPU, m_VtxRecoEffPUErr, int(m_VtxRecoEffPU), int(TotEvents_denominator) );
  calculateEfficiency(m_VtxTagZSelEff, m_VtxTagZSelEffErr, m_VtxTagZSel, m_VtxTagRec);
  calculateEfficiency(m_VtxTagZSelGenAnyEff, m_VtxTagZSelGenAnyEffErr, m_VtxTagZSel_wrtAnyGenVtx, m_VtxTagRec);
  calculateEfficiency(m_VtxTagZEff, m_VtxTagZEffErr, m_VtxTagZSel, int(TotEvents_denominator) );
  calculateEfficiency(m_VtxTagPurSelEff, m_VtxTagPurSelEffErr, m_VtxTagPurPos, m_VtxTagRec);
  //Fake rate for single vertex: N(vtx_reco > 1)/N(vtx_reco == 1) [for vtx_sim == 1]
  if (m_VtxOnlyTagRecSim1 > 0) {
    calculateEfficiency(m_VtxFakeReco, m_VtxFakeRecoErr, m_VtxTagRecSim1 - m_VtxOnlyTagRecSim1, m_VtxOnlyTagRecSim1);
  }
  calculateEfficiency(m_VtxTagTSelEff[TagVtxCorrect].Value(), m_VtxTagTSelEff[TagVtxCorrect].Error(), m_VtxTagThreshold[TagVtxCorrect], m_VtxTagRec);
  calculateEfficiency(m_VtxTagTSelEff[TagVtxSplitCorrect].Value(), m_VtxTagTSelEff[TagVtxSplitCorrect].Error(), m_VtxTagThreshold[TagVtxSplitCorrect], m_VtxTagRec);
  calculateEfficiency(m_VtxTagTSelEff[TagVtxSplitWrong].Value(), m_VtxTagTSelEff[TagVtxSplitWrong].Error(), m_VtxTagThreshold[TagVtxSplitWrong], m_VtxTagRec);
  calculateEfficiency(m_VtxTagTSelEff[TagVtxWrong].Value(), m_VtxTagTSelEff[TagVtxWrong].Error(), m_VtxTagThreshold[TagVtxWrong], m_VtxTagRec);
  calculateEfficiency(m_VtxTagTSelEff[TagVtxWrongFake].Value(), m_VtxTagTSelEff[TagVtxWrongFake].Error(), m_VtxTagThreshold[TagVtxWrongFake], m_VtxTagRec);

  //Post-processing of Truth-Matched histograms
  //normalize to the number of good simulated vertices
  h_TM_class->Sumw2();
  for (int ib=1; ib <= h_TM_class->GetNbinsX(); ib++) {
    h_TM_class->SetBinContent(ib, h_TM_class->GetBinContent(ib) / m_VtxGoodSim);
    h_TM_class->SetBinError(ib, h_TM_class->GetBinError(ib) / m_VtxGoodSim);
  }
  h_TM_fakeRelWeight->Sumw2();
  for (int ib=1; ib <= h_TM_fakeRelWeight->GetNbinsX(); ib++) {
    h_TM_fakeRelWeight->SetBinContent(ib, h_TM_fakeRelWeight->GetBinContent(ib) / m_VtxGoodSim);
    h_TM_fakeRelWeight->SetBinError(ib, h_TM_fakeRelWeight->GetBinError(ib) / m_VtxGoodSim);
  }

}

//Conform to the data-driven method: fit gaussian core, do not use RMS
void AnaVtxTree::computeMeanRmsFit(TH1F *h_mean, TH1F *h_rms, TH2F *h_2dplot) {
  Int_t totbins = h_2dplot->GetNbinsX();
  for (int nb=1; nb <= totbins; nb++) {
    TH1D *profileY = h_2dplot->ProjectionY("_px",nb,nb);
    TString useName = profileY->GetName();
    useName += nb;
    profileY->SetName(useName);
    if (debugLevel == -11) {
      outputNtp->cd();
      profileY->Write();
    }
    if (profileY->GetEntries() > 150.) {
      profileY->Fit("gaus");
      Double_t fit_res[3];
      Double_t * fit_err;
      profileY->GetFunction("gaus")->GetParameters(fit_res);
      fit_err = profileY->GetFunction("gaus")->GetParErrors();

      h_mean->SetBinContent(nb, fit_res[1]);
      h_mean->SetBinError(nb, fit_err[1]);
      h_rms->SetBinContent(nb, fit_res[2]);
      h_rms->SetBinError(nb, fit_err[2]);
    }
    delete profileY;
  }
}

void AnaVtxTree::computeMeanRms(TH1F *h_mean, TH1F *h_rms, TH2F *h_2dplot) {
  //The only difference respect to AnaVtxTree::computeMeanRms(TH1F *h_mean, TH1F *h_rms, TH1F *h_n)
  // is that this will reject outliers outside the Axis range!
  // The good thing is that you do not need to fill h_mean and h_rms :)
  // N.B. this is subject to bin width!
  if (debugLevel == -11) {
    cout << "[AnaVtxTree::computeMeanRms] Processing Histograms ";
    if (h_mean) {
      cout << h_mean->GetName();
    } else {
      cout << "0x0";
    }
    cout << ", ";
    if (h_rms) {
      cout << h_rms->GetName();
    } else {
      cout << "0x0";
    }
    cout << endl;
  }
  for (int nb=1; nb <= h_2dplot->GetNbinsX(); nb++) {
    Double_t denom = 0.0;
    Double_t mean = 0.0;
    Double_t rms = 0.0;
    for (int inb2d=1; inb2d <= h_2dplot->GetYaxis()->GetNbins(); inb2d++) {
      denom += h_2dplot->GetBinContent(nb, inb2d);
      mean += h_2dplot->GetYaxis()->GetBinCenter(inb2d)*h_2dplot->GetBinContent(nb, inb2d);
      rms += TMath::Power(h_2dplot->GetYaxis()->GetBinCenter(inb2d),2)*h_2dplot->GetBinContent(nb, inb2d);
    }
    if (debugLevel == -11) {
      cout << "nb = " << nb << " / " << h_mean->GetNbinsX() << ": " << endl;
      cout << " Denom = " << denom << ", Pre-Mean = " << mean << ", Pre-RMS = " << rms;
    }
    if (denom == 0) {
      mean = 0.0;
      rms = 0.0;
    } else {
      mean = mean / denom;
      if (denom > 1) {
        rms = TMath::Sqrt(rms / denom - TMath::Power(mean, 2));
      } else {
        rms = TMath::Abs(mean);
      }
    }
    if (h_mean && denom > 0) {
      h_mean->SetBinContent(nb, mean);
      h_mean->SetBinError(nb, rms / TMath::Sqrt(denom));
    }
    if (h_rms) {
      h_rms->SetBinContent(nb, rms);
      if (denom > 1) {
        h_rms->SetBinError(nb, TMath::Sqrt(2*TMath::Power(rms,4)/(denom - 1)));  //Error on RMS for normal distributions
      } else {
        h_rms->SetBinError(nb, rms);  //Error on RMS
      }
    }
    if (debugLevel == -11) {
      cout << ", Mean = " << mean << ", RMS = " << rms << endl;
    }
  }
}

void AnaVtxTree::computeMeanRmsReverseAxis(TH1F *h_mean, TH1F *h_rms, TH2F *h_2dplot) {
  //The only difference respect to AnaVtxTree::computeMeanRms
  // is that this will make a plot of the average over the X axis as function of the Y axis
  if (debugLevel == -11) {
    cout << "[AnaVtxTree::computeMeanRms] Processing Histograms ";
    if (h_mean) {
      cout << h_mean->GetName();
    } else {
      cout << "0x0";
    }
    cout << ", ";
    if (h_rms) {
      cout << h_rms->GetName();
    } else {
      cout << "0x0";
    }
    cout << endl;
  }
  for (int nb=1; nb <= h_2dplot->GetNbinsY(); nb++) {
    Double_t denom = 0.0;
    Double_t mean = 0.0;
    Double_t rms = 0.0;
    for (int inb2d=1; inb2d <= h_2dplot->GetXaxis()->GetNbins(); inb2d++) {
      denom += h_2dplot->GetBinContent(inb2d, nb);
      mean += h_2dplot->GetXaxis()->GetBinCenter(inb2d)*h_2dplot->GetBinContent(inb2d, nb);
      rms += TMath::Power(h_2dplot->GetXaxis()->GetBinCenter(inb2d),2)*h_2dplot->GetBinContent(inb2d, nb);
    }
    if (debugLevel == -11) {
      cout << "nb = " << nb << " / " << h_mean->GetYaxis()->GetNbins() << ": " << endl;
      cout << " Denom = " << denom << ", Pre-Mean = " << mean << ", Pre-RMS = " << rms;
    }
    if (denom == 0) {
      mean = 0.0;
      rms = 0.0;
    } else {
      mean = mean / denom;
      if (denom > 1) {
        rms = TMath::Sqrt(rms / denom - TMath::Power(mean, 2));
      } else {
        rms = TMath::Abs(mean);
      }
    }
    if (h_mean && denom > 0) {
      h_mean->SetBinContent(nb, mean);
      h_mean->SetBinError(nb, rms / TMath::Sqrt(denom));
    }
    if (h_rms) {
      h_rms->SetBinContent(nb, rms);
      if (denom > 1) {
        h_rms->SetBinError(nb, TMath::Sqrt(2*TMath::Power(rms,4)/(denom - 1)));  //Error on RMS for normal distributions
      } else {
        h_rms->SetBinError(nb, rms);  //Error on RMS
      }
    }
    if (debugLevel == -11) {
      cout << ", Mean = " << mean << ", RMS = " << rms << endl;
    }
  }
}

void AnaVtxTree::computeMeanRms(TH1F *h_mean, TH1F *h_rms, TH1F *h_n) {
  if (debugLevel == -11) {
    cout << "[AnaVtxTree::computeMeanRms] Processing Histograms " << h_mean->GetName() << ", " << h_rms->GetName() << endl;
  }
  for (int nb=1; nb <= h_mean->GetNbinsX(); nb++) {
    Double_t denom = h_n->GetBinContent(nb);
    Double_t mean = h_mean->GetBinContent(nb);
    Double_t rms = h_rms->GetBinContent(nb);
    if (debugLevel == -11) {
      cout << "nb = " << nb << " / " << h_mean->GetNbinsX() << ": " << endl;
      cout << " Denom = " << denom << ", Pre-Mean = " << mean << ", Pre-RMS = " << rms;
    }
    if (denom == 0) {
      mean = 0.0;
      rms = 0.0;
    } else {
      mean = mean / denom;
      if (denom > 1) {
        rms = TMath::Sqrt(rms / denom - TMath::Power(mean, 2));
      } else {
        rms = mean;
      }
    }
    h_mean->SetBinContent(nb, mean);
    h_mean->SetBinError(nb, rms / TMath::Sqrt(denom));
    h_rms->SetBinContent(nb, rms);
    if (denom > 1) {
      h_rms->SetBinError(nb, TMath::Sqrt(2*TMath::Power(rms,4)/(denom - 1)));  //Error on RMS for normal distributions
    } else {
      h_rms->SetBinError(nb, rms);  //Error on RMS
    }
    if (debugLevel == -11) {
      cout << ", Mean = " << mean << ", RMS = " << rms << endl;
    }
  }
}


bool AnaVtxTree::trigger(TString triggerName) {

  return true;

}

//Calculate efficiency with (binomial) errors, handle extreme cases with Poisson
int AnaVtxTree::calculateEfficiencyAsymErrors(Double_t &eff, Double_t &effErr_p, Double_t &effErr_m, const Int_t num, const Int_t den) {
  if (num == 0) {
    if (den == 0) {
      //numerator and denominator == 0, cannot say anything
      eff = 0.0;
      effErr_p = effErr_m = 0.0;
      return -1; //return -1, so the user knows it's not meaningful
    } else {
      //numerator == 0, assume poisson denominator, 1 sigma single-side error
      eff = 0.0;
      effErr_p = 1.149 / Double_t(den); //only positive error
      effErr_m = 0.0;
      return 0;
    }
  } else if (num == den) {
    //numerator == denominator != 0
    eff = 1.0;
    effErr_p = 0.0;
    effErr_m = 1.149 / Double_t(den); //only negative error
    return 2; //only negative error set
  } else {
    //standard binomial symmetric errors.. we can do better..
    eff = Double_t(num) / Double_t(den);
    effErr_p = effErr_m = TMath::Sqrt(eff*(1-eff)/Double_t(den));
    return 0;
  }
  return -2; //will never happen
}

//Calculate efficiency with (binomial) errors
int AnaVtxTree::calculateEfficiency(Double_t &eff, Double_t &effErr, const Int_t num, const Int_t den) {
  if (den == 0) {
    eff = 0.0;
    effErr = 0.0;
    return -1; //return -1, os the user knows it's not meaningful
  } else {
    //standard binomial errors
    eff = Double_t(num) / Double_t(den);
    effErr = TMath::Sqrt(eff*(1-eff)/Double_t(den));
    return 0;
  }
  return -2; //will never happen
}

//Determine if vertexToCheck is a good vertex, if -1 check all verttices and find the first tagged good one
//qualityVertexVersion: x.y.z -- y and/or z can be omitted if 0
// x = general quality requirenents:
//    1 = standard
//    2 = DR requirement on 2-,3-track vertices
//    3 = Require at least N (default=5, otherwise=z) tracks pointing to the vertex
// y = additional vertex requirements
//    0 = none
//    1 = Accept only events where tagged vertex has N(=z) tracks
// z = Optional parameter for a given "y" mode
// examples:
// 1: standard
// 2.1.5: DR req. on 2-,3-trk vertices + tagged vertex required with at least 5 tracks
int AnaVtxTree::isGoodVertex(int vertexToCheck) {
  int taggedVertex = -1;
  //cout << "[AnaVtxTree] AnaVtxTree::isGoodVertex qualityVertexVersion = " << qualityVertexVersion << endl;
  if (qualityVertexVersion == 1) {
    //standard-definition
    Double_t minNTracks = qualityVertexParameter;
    if (minNTracks == 0) {
      minNTracks = 2;  //default value
    }
    //cout << "[AnaVtxTree] AnaVtxTree::isGoodVertex minNTracks = " << minNTracks << endl;
    if (vertexToCheck == -1) {
      //check only algorithm-tagged vertex
      //if ((*vxnbc_nTracks)[0] >= minNTracks && (*vxnbc_cov_x)[0]!=0 && (*vxnbc_cov_y)[0]!=0 && vxnbc_n > 1)
      if ((*vxnbc_nTracks)[0] >= minNTracks)
        if ((*vxnbc_cov_x)[0] != 0)
          if ((*vxnbc_cov_y)[0] != 0) {
            taggedVertex=0;
          }
    } else {
      if ((*vxnbc_nTracks)[vertexToCheck] >= minNTracks && (*vxnbc_cov_x)[vertexToCheck]!=0 && (*vxnbc_cov_y)[vertexToCheck]!=0 && vxnbc_n > 1) {
        taggedVertex=vertexToCheck;
      }
    }
  } else if (qualityVertexVersion == 2) {
    //require vertex probability
    for (int ivx=0; ivx < vxnbc_n - 1; ++ivx) {
      if (vertexToCheck > -1 && ivx != vertexToCheck) {
        continue;  //next vertex
      }
      double probVertex = TMath::Prob((*vxnbc_chi2)[ivx], (*vxnbc_ndof)[ivx]);
      //      double dz=0.0;
      //for (vector<int>::iterator nvtrk = vxnbc_trk_index->at(ivx).begin(); nvtrk != vxnbc_trk_index->at(ivx).end(); ++nvtrk) {
      //  dr += TMath::Power()
      //}
      double Dr2=-999.;
      double Dr3=-999.;
      if (vxnbc_nTracks->at(ivx) == 2) {
        //Properties for events with 2 tracks at the Tagged vertex
        //evaluate delta d0 for vertices with 2 tracks
        if (vxnbc_trk_index->at(ivx).size() == 2) {
          // (Delta r)^2 = (d0)^2 + (d'0)^2 - 2(d0)(d'0)cos(Delta phi)
          double d01 = trk_d0_wrtPV->at(vxnbc_trk_index->at(ivx).at(0));
          double d02 = trk_d0_wrtPV->at(vxnbc_trk_index->at(ivx).at(1));
          double phi1 = trk_phi_wrtPV->at(vxnbc_trk_index->at(ivx).at(0));
          double phi2 = trk_phi_wrtPV->at(vxnbc_trk_index->at(ivx).at(1));
          double dphi = TMath::Abs(TVector2::Phi_mpi_pi(phi1 - phi2));
          Dr2 = TMath::Power(d01,2) + TMath::Power(d02,2) - 2*d01*d02*TMath::Cos(dphi);
          Dr2 = TMath::Sqrt(Dr2);
        }
      } else if (vxnbc_nTracks->at(ivx) == 3) {
        if (vxnbc_trk_index->at(ivx).size() == 3) {
          //evaluate average distance
          double d01 = trk_d0_wrtPV->at(vxnbc_trk_index->at(ivx).at(0));
          double d02 = trk_d0_wrtPV->at(vxnbc_trk_index->at(ivx).at(1));
          double d03 = trk_d0_wrtPV->at(vxnbc_trk_index->at(ivx).at(2));
          double phi1 = trk_phi_wrtPV->at(vxnbc_trk_index->at(ivx).at(0));
          double phi2 = trk_phi_wrtPV->at(vxnbc_trk_index->at(ivx).at(1));
          double phi3 = trk_phi_wrtPV->at(vxnbc_trk_index->at(ivx).at(2));
          double dphi12 = TMath::Abs(TVector2::Phi_mpi_pi(phi1 - phi2));
          double dphi13 = TMath::Abs(TVector2::Phi_mpi_pi(phi1 - phi3));
          double dphi23 = TMath::Abs(TVector2::Phi_mpi_pi(phi1 - phi2));
          Dr3 = TMath::Sqrt(TMath::Power(d01,2) + TMath::Power(d02,2) - 2*d01*d02*TMath::Cos(dphi12));
          Dr3 += TMath::Sqrt(TMath::Power(d01,2) + TMath::Power(d03,2) - 2*d01*d02*TMath::Cos(dphi13));
          Dr3 += TMath::Sqrt(TMath::Power(d02,2) + TMath::Power(d03,2) - 2*d01*d02*TMath::Cos(dphi23));
          Dr3 = Dr3 / 3;
        }
      }

      if ( ((*vxnbc_nTracks)[ivx] == 2 && probVertex > 0.01 && Dr2 < 2.5) ||
           ((*vxnbc_nTracks)[ivx] == 3 && probVertex > 1E-08 && Dr3 < 2.5) ||
           ((*vxnbc_nTracks)[ivx] > 3  && probVertex > 0)) {
        taggedVertex = ivx;
        break; // Vertex with highest sumPt2 with quality requirements
      }
    }
  } else if (qualityVertexVersion == 3) {
    Double_t minNTracks = qualityVertexParameter;
    if (minNTracks == 0) {
      minNTracks = 5;  //default value
    }
    for (int ivx=0; ivx < vxnbc_n - 1; ++ivx) {
      if (vertexToCheck > -1 && ivx != vertexToCheck) {
        continue;  //next vertex
      }
      if ((*vxnbc_nTracks)[ivx] >= minNTracks && (*vxnbc_cov_x)[ivx]!=0 && (*vxnbc_cov_y)[ivx]!=0) {
        taggedVertex = ivx;
        break;
      }
    }
  }

  if (qualityVertexModifier == 1) {
    //if tagged vertex has a number of tracks != qualityVertexParameter, un-tag the vertex
    if ((*vxnbc_nTracks)[taggedVertex] != qualityVertexParameter) {
      taggedVertex = -1;
    }
  }

  return taggedVertex;
}

AnaVtxTree::TagVtxEfficiencyType AnaVtxTree::getTagVertexClass(Int_t TagVtx, Int_t simVtx, bool forceDebug) {
  if (debugLevel == -12 || forceDebug) {
    cout << "[AnaVtxTree::getTagVertexClass] Processing TagVtx = " << TagVtx << ", SimVtx = " << simVtx << endl;
  }
  Double_t T_tag = getTDistance(TagVtx, simVtx);
  Double_t T_min= 999;
  Double_t T_min_idx = -1;
  Int_t num_vtx_within_T_Threshold=0;
  //find minimum distnace vertex
  if (debugLevel == -12 || forceDebug) {
    cout << "Vertex T-distances. Raw N vertex = " << (vxnbc_n-1) << ", N tracks = " << (*vxnbc_nTracks)[TagVtx];
    cout << ". VxSim Position: " << "X = " << (*mcvtx_x)[simVtx] << ", Y = " << (*mcvtx_y)[simVtx]
         << ", Z = " << (*mcvtx_z)[simVtx] << endl;
  }
  for (int vtx=0; vtx < vxnbc_n - 1; vtx++) {
    //skip bad vertices
    if (isGoodVertex(vtx) != vtx) {
      continue;
    }
    //calculate distance
    Double_t T_dist = getTDistance(vtx, simVtx);
    if (debugLevel == -12 || forceDebug) {
      cout << "X = " << (*vxnbc_x)[vtx] << " +/- " << TMath::Sqrt((*vxnbc_cov_x)[vtx]);
      cout << ". Y = " << (*vxnbc_y)[vtx] << " +/- " << TMath::Sqrt((*vxnbc_cov_y)[vtx]);
      cout << ". Z = " << (*vxnbc_z)[vtx] << " +/- " << TMath::Sqrt((*vxnbc_cov_z)[vtx]) << endl;
      cout <<vtx << ": " << T_dist << endl;
    }
    if (T_dist < T_min) {
      T_min_idx = vtx;
      T_min = T_dist;
    }
    if (T_dist <= VtxEffCriteriaT)
      //count verticeis within threshold VtxEffCriteriaT
    {
      num_vtx_within_T_Threshold++;
    }
  }
  if (debugLevel == -12 || forceDebug) {
    cout << "T_min = " << T_min << " (idx = " << T_min_idx << ")" << endl;
  }

  //if T_min > VtxEffCriteriaT, -> TagVtxWrongFake
  if (T_min > VtxEffCriteriaT) {
    return TagVtxWrongFake;
  }
  //now T_min is within VtxEffCriteriaT
  if ( (T_min_idx == TagVtx) && (num_vtx_within_T_Threshold == 1) ) {
    return TagVtxCorrect;
  }
  if ( (T_min_idx == TagVtx) && (num_vtx_within_T_Threshold > 1) ) {
    return TagVtxSplitCorrect;
  }
  if ( (T_min_idx != TagVtx) && (T_tag <= VtxEffCriteriaT) ) {
    return TagVtxSplitWrong;
  }
  if ( (T_min_idx != TagVtx) && (T_tag > VtxEffCriteriaT) ) {
    return TagVtxWrong;
  }

  cout << "[AnaVtxTree] YOU SHOULD NOT BE HERE!!!! " << endl;
  return TagVtxNumTypes;
}

Double_t AnaVtxTree::getTDistance(Int_t VtxIdx, Int_t simVtxIdx) {
  Double_t T_distance = 0.0;
  T_distance += TMath::Power( ((*vxnbc_x)[VtxIdx] - (*mcvtx_x)[simVtxIdx]), 2) / (*vxnbc_cov_x)[VtxIdx];
  T_distance += TMath::Power( ((*vxnbc_y)[VtxIdx] - (*mcvtx_y)[simVtxIdx]), 2) / (*vxnbc_cov_y)[VtxIdx];
  T_distance += TMath::Power( ((*vxnbc_z)[VtxIdx] - (*mcvtx_z)[simVtxIdx]), 2) / (*vxnbc_cov_z)[VtxIdx];
  T_distance = TMath::Sqrt(T_distance);

  return T_distance;
}

TString AnaVtxTree::SetQualityVertexVersion(TString version) {
  if (version == "") {
    return (TString(qualityVertexVersion)+TString(".")+TString(qualityVertexModifier)+TString(".")+TString(qualityVertexParameter));
  }
  TObjArray *StrVersionArray = version.Tokenize(".");
  if (StrVersionArray->GetEntries() > 3) {
    cerr << "ERROR: qualityVertexVersion must be in the format x.y.z" << endl;
    return TString("");
  }
  TString ver="0";
  TString mod="0";
  TString par="0";
  ver = ((TObjString*)(*StrVersionArray)[0])->String();
  if (StrVersionArray->GetEntries() > 1) {
    mod = ((TObjString*)(*StrVersionArray)[1])->String();
  }
  if (StrVersionArray->GetEntries() > 2) {
    par = ((TObjString*)(*StrVersionArray)[2])->String();
  }
  qualityVertexVersion = ver.Atoi();
  qualityVertexModifier = mod.Atoi();
  qualityVertexParameter = par.Atoi();

  return version;
}

Int_t AnaVtxTree::GetQualityVertexVersion(Int_t rev) {
  switch(rev) {
    case 0:
      return qualityVertexVersion;
      break;
    case 1:
      return qualityVertexModifier;
      break;
    case 2:
      return qualityVertexParameter;
      break;
  }
  return -1;
}

void AnaVtxTree::SetTriggerFilter(TString p_triggerName, TChain *p_trigMetaDataTree) {
  triggerName = p_triggerName;
  if (p_trigMetaDataTree) {
    trigMetaDataTree = p_trigMetaDataTree;
  }
}

void AnaVtxTree::SetMuFilter(float p_muMin, float p_muMax) {
  muFilter_min = p_muMin;
  muFilter_max = p_muMax;
}

//Count number of generated tracks
//TODO: Add matching with offline tracks
//TODO: Add matching with with main interaction vertex?
Int_t AnaVtxTree::countNGenTracks(Float_t pEtaMax, Float_t pPtMin) {
  /*
  Int_t nGenTrk = 0;
  for (int nt=0; nt < mc_n; nt++)
    if ((TMath::Abs((*mc_gen_pt)[nt]) > pPtMin) &&
  (TMath::Abs((*mc_gen_eta)[nt]) < pEtaMax) )
      nGenTrk++;
  */
  // Up to now count number of tracks associated to a real track
  Int_t nGenTrk = 0;
  for (int nt=0; nt < trk_n; nt++)
    if ((TMath::Abs((*trk_pt)[nt]) > pPtMin) &&
        (TMath::Abs((*trk_eta)[nt]) < pEtaMax) &&
        ((*trk_mcpart_probability)[nt] > 0.5) ) {
      nGenTrk++;
    }

  return nGenTrk;
}

Float_t AnaVtxTree::GetVtxDistance(Int_t VxIndex1, VtxType VxType1, Int_t VxIndex2, VtxType VxType2, VtxDistanceMetric cfgMetric) {
  double x[2] = {0,0}; //positions
  double y[2] = {0,0};
  double z[2] = {0,0};
  double e2x[2] = {0,0}; //errors^2
  double e2y[2] = {0,0};
  double e2z[2] = {0,0};
  Int_t VxIndex[2] = {VxIndex1, VxIndex2};
  VtxType VxType[2] = {VxType1, VxType2};
  for (int ivx=0; ivx < 2; ivx++) {
    switch (VxType[ivx]) {
      case VtxSimulated:
        x[ivx] = mcvtx_x->at(VxIndex[ivx]);
        y[ivx] = mcvtx_y->at(VxIndex[ivx]);
        z[ivx] = mcvtx_z->at(VxIndex[ivx]);
        //leave errors at 0
        break;
      case VtxReconstructedNoBC:
        x[ivx] = vxnbc_x->at(VxIndex[ivx]);
        y[ivx] = vxnbc_y->at(VxIndex[ivx]);
        z[ivx] = vxnbc_z->at(VxIndex[ivx]);
        e2x[ivx] = vxnbc_cov_x->at(VxIndex[ivx]);
        e2y[ivx] = vxnbc_cov_y->at(VxIndex[ivx]);
        e2z[ivx] = vxnbc_cov_z->at(VxIndex[ivx]);
        break;
      default:
        std::cerr << "WARNING: Wrong type of vertex in RealTrackTree::GetVtxDistance" << std::endl;
        x[ivx] = y[ivx] = z[ivx] = 0;
    }
  }

  double dx,dy,dz;
  double dxErr, dyErr, dzErr;
  dx = x[0] - x[1];
  dy = y[0] - y[1];
  dz = z[0] - z[1];
  dxErr = TMath::Sqrt(e2x[0] + e2x[1]);
  dyErr = TMath::Sqrt(e2y[0] + e2y[1]);
  dzErr = TMath::Sqrt(e2z[0] + e2z[1]);

  if (cfgMetric == VtxDistM_deltaZ) {
    // Simple delta z
    return dz;
  }
  if (cfgMetric == VtxDistM_deltaZsig) {
    // Delta z / error
    return dz/dzErr;
  }
  if (cfgMetric == VtxDistM_3Dsig) {
    return TMath::Sqrt(dx*dx/dxErr/dxErr + dy*dy/dyErr/dyErr + dz*dz/dzErr/dzErr);
  }

  //nope
  cerr << "ERROR: Invalid cfgMetric in RealTrackTree::GetVtxDistance" << endl;
  return 0;
}

TString AnaVtxTree::MakeHistoNamesCOORD(TString t, int ia) {
  const char* AxisNames[3] = {"X", "Y", "Z"};
  TString newT = t;
  newT.ReplaceAll("COORD", AxisNames[ia]);
  return newT;
}

void AnaVtxTree::SetHS(bool p_skip_hs) {

  skip_hs = p_skip_hs;
  if (skip_hs) {
    cout << "[AnaVtxTree] INFO : Skipping hard scattering events." << endl;
  } else {
    cout << "[AnaVtxTree] INFO : Including hard scattering events." << endl;
  }

}

void AnaVtxTree::DoTiming() {

  cout << "[AnaVtxTree] INFO : Enabling timing for code performance study." << endl;
  do_timing = true;

}

void AnaVtxTree::SetMaxChi2Ndf(Float_t p_max_chi2ndf) {

  cout << "[AnaVtxTree] INFO : Filtering vertices, max chi2/ndf = " << p_max_chi2ndf << endl;
  max_chi2ndf = p_max_chi2ndf;

}


