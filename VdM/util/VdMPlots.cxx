#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <stdlib.h>


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
#include "TLine.h"

#include "GlobalSettings/GlobalSettings.h"
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
  GlobalSettings gs;
  TString vtx_method;
  Int_t p_nTrkCut;
  TString run;
  TString settings;

  // --- Scan for command line parameters
  int c;
  extern int optind;
  extern char* optarg;
  while (1) {
    int option_index = 0;
    static struct option long_options[] = { {0,0,0,0} };
    c=getopt_long(argc, argv, "r:s:m:n:",long_options,&option_index);
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
      case 'n': {
        stringstream ss;
        ss << optarg;
        try {
          ss >> p_nTrkCut;
        } catch(...) {
          cerr << "NTrkCut " << optarg << " not an integer? Exiting..." << endl;
          exit(1);
        }
        break;
      }
    }
  }

  if (!p_nTrkCut) {
    cerr << "ERROR: Specify track cut with -n. Exiting..." << endl;
    exit(1);
  }

  if (!vtx_method) {
    cerr << "ERROR: Specify vertex method with -m. Exiting..." << endl;
    exit(1);
  }

  // Default method
  TString default_fit;
  if (run == "182013") {
    default_fit = "sgcl";
  } else if (run == "201351") {
    default_fit = "spline";
  } else if (run == "215021") {
    default_fit = "sgc";
  }
  stringstream ss_folder;
  ss_folder << GlobalSettings::path_outputVdM << GlobalSettings::path_VdM_prefix << run << "/" << settings << "/" << vtx_method << "/";
  TString filename(ss_folder.str());
  filename += "vdm_results.root";

  TFile *f_in = new TFile(filename, "READ");
  if (!f_in->IsOpen()) {
    cerr << "Unable to open input file " << filename << ". Exiting..." << endl;
    exit(1);
  }

  // -- Make TCanvases showing sigma_vis, lumi_sp, Sigma_x, Sigma_y, c_x/mu_max_x, c_y/mu_max_y, chi2/ndf of all fits.

  vector<TString> fit_functions;
  fit_functions.push_back("sgcl");
  fit_functions.push_back("sgc");
  fit_functions.push_back("sg");
  fit_functions.push_back("dg");
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
  }
    else if (run == "215021") {
    bcids.push_back(1);
    bcids.push_back(2361);
    bcids.push_back(2881);
  }
  Int_t n_bcids = bcids.size();

  vector<Int_t> scans;
  if (run == "182013") {
    scans.push_back(1);
    scans.push_back(2);
  } else if (run == "201351") {
    scans.push_back(1);
    scans.push_back(2);
    scans.push_back(3);
  }
    else if (run == "215021") {
    scans.push_back(1);
  }

  vector<TString> axes;
  axes.push_back("x");
  axes.push_back("y");

  // -- Read TTree into maps
  std::map<TString, Float_t> m_sigma_vis, m_sigma_vis_err, m_lumi_sp, m_lumi_sp_err, m_lumi, m_lumi_err;
  std::map<TString, Float_t> m_Sigma_x, m_Sigma_x_err, m_Sigma_y, m_Sigma_y_err, m_mu_max_x, m_mu_max_x_err, m_mu_max_y, m_mu_max_y_err, m_c_x, m_c_x_err, m_c_y, m_c_y_err;
  std::map<TString, Float_t> m_chi2ndf_x, m_chi2ndf_y;

  TTree *t_vdm = (TTree*)f_in->Get("VdmResults");

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

    if (NTrkCut != p_nTrkCut) {
      continue;
    }

    TString tag_base = "BCID";
    tag_base += BCID;
    tag_base += "_scan";
    tag_base += Scan;

    for (vector<TString>::iterator fit_function = fit_functions.begin(); fit_function != fit_functions.end(); ++fit_function) {

      TString tag = tag_base;
      tag += "_";
      tag += *fit_function;
      Int_t vector_idx = distance(fit_functions.begin(), fit_function);

      m_sigma_vis[tag] = SigmaVis->at(vector_idx);
      m_sigma_vis_err[tag] = SigmaVisErr->at(vector_idx);
      m_lumi_sp[tag] = LumiSp->at(vector_idx);
      m_lumi_err[tag] = LumiErr->at(vector_idx);
      m_lumi[tag] = Lumi->at(vector_idx);
      m_lumi_sp_err[tag] = LumiSpErr->at(vector_idx);
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
      m_chi2ndf_x[tag] = Chi2NdfX->at(vector_idx);
      m_chi2ndf_y[tag] = Chi2NdfY->at(vector_idx);

    }
  }

  // -- sigma_vis: for each fit function and scan, make a 3-point TGraph (one point for each BCID).
  std::map<TString, TGraphErrors*> tg_sigma_vis;
  std::map<TString, std::map<Int_t, Color_t> > marker_colors;
  marker_colors["sgcl"][1] = kBlack;
  marker_colors["sgcl"][2] = kBlack;
  marker_colors["sgcl"][3] = kBlack;
  marker_colors["sgc"][1] = kViolet+10;
  marker_colors["sgc"][2] = kViolet+5;
  marker_colors["sgc"][3] = kViolet;
  marker_colors["sgc"][1] = kRed+10;
  marker_colors["sgc"][2] = kRed+5;
  marker_colors["sgc"][3] = kRed;
  marker_colors["sg"][1] = kPink+10;
  marker_colors["sg"][2] = kPink-6;
  marker_colors["sg"][3] = kPink;
  marker_colors["sg"][1] = kGreen+10;
  marker_colors["sg"][2] = kGreen-6;
  marker_colors["sg"][3] = kGreen;
  marker_colors["dg"][1] = kAzure+1;
  marker_colors["dg"][2] = kAzure-6;
  marker_colors["dg"][3] = kAzure+5;
  marker_colors["spline"][1] = kSpring;
  marker_colors["spline"][2] = kSpring-7;
  marker_colors["spline"][3] = kSpring+7;

  std::map<TString, std::map<Int_t, Int_t> > marker_styles;
  marker_styles["sgcl"][1] = 20;
  marker_styles["sgcl"][2] = 24;
  marker_styles["sgcl"][3] = 21;
  marker_styles["sgc"][1] = 33;
  marker_styles["sgc"][2] = 27;
  marker_styles["sgc"][3] = 33;
  marker_styles["sg"][1] = 23;
  marker_styles["sg"][2] = 32;
  marker_styles["sg"][3] = 23;
  marker_styles["dg"][1] = 21;
  marker_styles["dg"][2] = 25;
  marker_styles["dg"][3] = 22;
  marker_styles["spline"][1] = 22;
  marker_styles["spline"][2] = 26;
  marker_styles["spline"][3] = 23;

  std::map<TString, Float_t> marker_offsets;
  marker_offsets["sgcl"] = 0.;
  marker_offsets["sgc"] = -0.04;
  marker_offsets["sg"] = -0.08;
  marker_offsets["dg"] = 0.04;
  marker_offsets["spline"] = 0.08;

  std::map<TString, TString> fit_functions_formatted;
  fit_functions_formatted["sgcl"] = "SG+p(L)";
  fit_functions_formatted["sgc"] = "SG+p";
  fit_functions_formatted["sg"] = "SG";
  fit_functions_formatted["dg"] = "DG+p";
  fit_functions_formatted["spline"] = "Spline";

  std::map<Int_t, TString> scans_formatted;
  scans_formatted[1] = "Scan I";
  scans_formatted[2] = "Scan II";
  scans_formatted[3] = "Scan III";

  for (vector<TString>::iterator fit_function = fit_functions.begin(); fit_function != fit_functions.end(); ++fit_function) {
    for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {

      TString tag_base = "scan";
      tag_base += *scan;
      tag_base += "_";
      tag_base += *fit_function;

      TString tgname = "tg_sigma_vis_";
      tgname += tag_base;
      tg_sigma_vis[tag_base] = new TGraphErrors(bcids.size());
      tg_sigma_vis[tag_base]->SetName(tgname);
      tg_sigma_vis[tag_base]->GetHistogram()->GetYaxis()->Set(n_bcids, 0.5, n_bcids + 0.5);
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString label_tmp = "";
        label_tmp += *bcid;
        tg_sigma_vis[tag_base]->GetHistogram()->GetYaxis()->SetBinLabel(distance(bcids.begin(), bcid) + 1, label_tmp);
      }

      tg_sigma_vis[tag_base]->SetMarkerStyle(marker_styles[*fit_function][*scan]);
      tg_sigma_vis[tag_base]->SetMarkerColor(marker_colors[*fit_function][*scan]);
      tg_sigma_vis[tag_base]->SetMarkerSize(1.5);

      Int_t point = -1;

      for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {

        point++;
        TString tag = "BCID";
        tag += *bcid;
        tag += "_";
        tag += tag_base;
        tg_sigma_vis[tag_base]->SetPoint(point, m_sigma_vis[tag], point + 1. + marker_offsets[*fit_function]);
        if (*fit_function != "dg") {
          tg_sigma_vis[tag_base]->SetPointError(point, m_sigma_vis_err[tag], 0.);
        }
      }
    }
  }

  // -- Next, canvases: summary canvas, plus SGCL only canvas.

  TCanvas* c_sigma_vis_summary = new TCanvas("c_sigma_vis_summary", "c_sigma_vis_summary", 1200, 800);
  TLegend* l_sigma_vis_summary = new TLegend(0.2, 0.2, 0.4, 0.8);
  l_sigma_vis_summary->SetFillColor(0);
  l_sigma_vis_summary->SetBorderSize(0);
  bool draw_first = true;

  Float_t ymin, ymax;

  for (vector<TString>::iterator fit_function = fit_functions.begin(); fit_function != fit_functions.end(); ++fit_function) {
    for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {

      TString tag_base = "scan";
      tag_base += *scan;
      tag_base += "_";
      tag_base += *fit_function;

      TString draw_options;
      if (draw_first) {
        draw_options = "ap";
        draw_first = false;
        Float_t old_max = tg_sigma_vis[tag_base]->GetHistogram()->GetXaxis()->GetXmax();
        Float_t old_min = tg_sigma_vis[tag_base]->GetHistogram()->GetXaxis()->GetXmin();
        Float_t new_max = old_max + 0.5 * (old_max - old_min);
        Float_t new_min = old_min - 1.2 * (old_max - old_min);
        tg_sigma_vis[tag_base]->GetHistogram()->GetXaxis()->Set(100, new_min, new_max);
        tg_sigma_vis[tag_base]->GetHistogram()->GetYaxis()->Set(n_bcids, 0.5, n_bcids + 0.5);
        for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
          TString label_tmp = "";
          label_tmp += *bcid;
          tg_sigma_vis[tag_base]->GetHistogram()->GetYaxis()->SetBinLabel(distance(bcids.begin(), bcid) + 1, label_tmp);
        }
        ymin = tg_sigma_vis[tag_base]->GetHistogram()->GetMinimum();
        ymax = tg_sigma_vis[tag_base]->GetHistogram()->GetMaximum();
        tg_sigma_vis[tag_base]->GetXaxis()->SetTitle("#sigma_{vis} (mb)");
        tg_sigma_vis[tag_base]->GetYaxis()->SetTitle("BCID");
      } else {
        draw_options = "p";
      }

      tg_sigma_vis[tag_base]->Draw(draw_options);
      TString legend_entry = "";
      legend_entry += fit_functions_formatted[*fit_function];
      legend_entry += ", ";
      legend_entry += scans_formatted[*scan];
      l_sigma_vis_summary->AddEntry(tg_sigma_vis[tag_base], legend_entry, "lp");
    }
  }
  l_sigma_vis_summary->Draw();

  // -- Get lumi-weighted average of sigma_vis for sgcl
  Float_t weight = 0.;
  Float_t avg_sum = 0.;
  Float_t avg_sum2 = 0.;
  Float_t err_sum2 = 0.;
  for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {
    for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {

      TString tag = "BCID";
      tag += *bcid;
      tag += "_scan";
      tag += *scan;
      tag += "_sgcl";
      Float_t current_weight = m_lumi[tag];
      weight += current_weight;
      avg_sum += m_sigma_vis[tag] * current_weight;
      avg_sum2 += current_weight * TMath::Power(m_sigma_vis[tag], 2);
      err_sum2 += TMath::Power(current_weight * m_sigma_vis_err[tag], 2);

    }
  }
  Float_t sigma_vis_avg = avg_sum / weight;
  Float_t sigma_vis_err = TMath::Sqrt(err_sum2) / weight;
  Float_t sigma_vis_avg2 = avg_sum2 / weight;
  Float_t sigma_vis_rms = TMath::Sqrt(sigma_vis_avg2 - (sigma_vis_avg * sigma_vis_avg));
  cout << "SigmaVis results:" << endl;
  cout << "\t Avg  = " << sigma_vis_avg << endl;
  cout << "\t Err  = " << sigma_vis_err << endl;
  cout << "\t Avg2 = " << sigma_vis_avg2 << endl;
  cout << "\t RMS  = " << sigma_vis_rms << endl;
  TLine *tl_sigma_vis_avg = new TLine(sigma_vis_avg, ymin, sigma_vis_avg, ymax);
  tl_sigma_vis_avg->SetLineStyle(2);
  tl_sigma_vis_avg->SetLineColor(kBlack);
  tl_sigma_vis_avg->SetLineWidth(1);
  tl_sigma_vis_avg->Draw("same");
  TLine *tl_sigma_vis_minus = new TLine(sigma_vis_avg - sigma_vis_rms, ymin, sigma_vis_avg - sigma_vis_rms, ymax);
  tl_sigma_vis_minus->SetLineStyle(2);
  tl_sigma_vis_minus->SetLineColor(kGray);
  tl_sigma_vis_minus->SetLineWidth(1);
  tl_sigma_vis_minus->Draw("same");
  TLine *tl_sigma_vis_plus = new TLine(sigma_vis_avg + sigma_vis_rms, ymin, sigma_vis_avg + sigma_vis_rms, ymax);
  tl_sigma_vis_plus->SetLineStyle(2);
  tl_sigma_vis_plus->SetLineColor(kGray);
  tl_sigma_vis_plus->SetLineWidth(1);
  tl_sigma_vis_plus->Draw("same");

  stringstream ss_figures;
  ss_figures << ss_folder.str() << "figures/" << c_sigma_vis_summary->GetName() << "_" << vtx_method << p_nTrkCut << ".pdf";
  c_sigma_vis_summary->SaveAs(TString(ss_figures.str()));

  TCanvas* c_sigma_vis_sgcl = new TCanvas("c_sigma_vis_sgcl", "c_sigma_vis_sgcl", 1200, 800);
  TLegend* l_sigma_vis_sgcl = new TLegend(0.2, 0.4, 0.4, 0.6);
  l_sigma_vis_sgcl->SetFillColor(0);
  l_sigma_vis_sgcl->SetBorderSize(0);
  draw_first = true;
  for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {

    TString tag_base = "scan";
    tag_base += *scan;
    tag_base += "_";
    tag_base += default_fit;

    TString draw_options;
    if (draw_first) {
      draw_options = "ap";
      draw_first = false;
      Float_t old_max = tg_sigma_vis[tag_base]->GetHistogram()->GetXaxis()->GetXmax();
      Float_t old_min = tg_sigma_vis[tag_base]->GetHistogram()->GetXaxis()->GetXmin();
      Float_t new_max = old_max + 0.5 * (old_max - old_min);
      if (run  == "201351") {
        new_max = old_max + 3 * (old_max - old_min);
      }
      Float_t new_min = old_min - 1.2 * (old_max - old_min);
      tg_sigma_vis[tag_base]->GetHistogram()->GetXaxis()->Set(100, new_min, new_max);
      tg_sigma_vis[tag_base]->GetHistogram()->GetYaxis()->Set(n_bcids, 0.5, n_bcids + 0.5);
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString label_tmp = "";
        label_tmp += *bcid;
        tg_sigma_vis[tag_base]->GetHistogram()->GetYaxis()->SetBinLabel(distance(bcids.begin(), bcid) + 1, label_tmp);
      }
    } else {
      draw_options = "p";
    }

    tg_sigma_vis[tag_base]->Draw(draw_options);
    TString legend_entry = "";
    legend_entry += fit_functions_formatted[default_fit];
    legend_entry += ", ";
    legend_entry += scans_formatted[*scan];
    l_sigma_vis_sgcl->AddEntry(tg_sigma_vis[tag_base], legend_entry, "lp");
  }
  l_sigma_vis_sgcl->Draw();
  tl_sigma_vis_avg->Draw("same");
  tl_sigma_vis_minus->Draw("same");
  tl_sigma_vis_plus->Draw("same");

  ss_figures.str("");
  ss_figures << ss_folder.str() << "figures/" << c_sigma_vis_sgcl->GetName() << "_" << vtx_method << p_nTrkCut << ".pdf";

  c_sigma_vis_sgcl->SaveAs(TString(ss_figures.str()));



  // -- Chi2 / NDF plot
  std::map<TString, std::map<TString, TGraph*> > tg_chi2ndf;
  std::map<TString, TCanvas*> c_chi2ndf;
  std::map<TString, TLegend*> l_chi2ndf;

  for (vector<TString>::iterator axis = axes.begin(); axis != axes.end(); ++axis) {

    TString cname = "c_chi2ndf_";
    cname += *axis;
    cname += "_";
    cname += vtx_method;
    cname += p_nTrkCut;
    c_chi2ndf[*axis] = new TCanvas(cname, cname, 1200, 800);
    c_chi2ndf[*axis]->SetRightMargin(0.15);
    l_chi2ndf[*axis] = new TLegend(0.75* 0.85, 0.75, 0.95 * 0.85, 0.9);
    l_chi2ndf[*axis]->SetFillColor(0);
    l_chi2ndf[*axis]->SetBorderSize(1);

    bool draw_first = true;

    for (vector<TString>::iterator fit_function = fit_functions.begin(); fit_function != fit_functions.end(); ++fit_function) {

      if (*fit_function == "spline") {
        continue;
      }
      if (*fit_function == "sgc") {
        continue;
      }

      TString draw_options;
      if (draw_first) {
        draw_first = false;
        draw_options = "ap";
      } else {
        draw_options = "p";
      }

      tg_chi2ndf[*axis][*fit_function] = new TGraphErrors(bcids.size() * scans.size());
      TString tgname = "tg_chi2ndf_";
      tgname += *axis;
      tgname += "_";
      tgname += *fit_function;
      tg_chi2ndf[*axis][*fit_function]->SetName(tgname);

      Int_t current_point = 0;

      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        for (vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {
          TString tag = "BCID";
          tag += *bcid;
          tag += "_scan";
          tag += *scan;
          tag += "_";
          tag += *fit_function;
          tg_chi2ndf[*axis][*fit_function]->SetPoint(current_point, current_point + 1, (*axis == "x" ? m_chi2ndf_x[tag] : m_chi2ndf_y[tag]));
          current_point++;
        }
      }

      tg_chi2ndf[*axis][*fit_function]->SetMarkerColor(marker_colors[*fit_function][1]);
      tg_chi2ndf[*axis][*fit_function]->SetMarkerStyle(marker_styles[*fit_function][1]);
      tg_chi2ndf[*axis][*fit_function]->SetMarkerSize(1.2);
      tg_chi2ndf[*axis][*fit_function]->GetXaxis()->SetTitle("");
      tg_chi2ndf[*axis][*fit_function]->GetYaxis()->SetTitle("#chi^{2} / NDF");

      tg_chi2ndf[*axis][*fit_function]->GetHistogram()->GetXaxis()->Set(bcids.size() * scans.size(), 0.5, bcids.size() * scans.size() + 0.5);
      Int_t current_bin = 1;
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        for (vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {
          TString bin_label = "Scan ";
          bin_label += *scan;
          bin_label += " / BCID ";
          bin_label += *bcid;
          tg_chi2ndf[*axis][*fit_function]->GetHistogram()->GetXaxis()->SetBinLabel(current_bin++, bin_label);
        }
      }
      tg_chi2ndf[*axis][*fit_function]->SetMinimum(0.);
      tg_chi2ndf[*axis][*fit_function]->SetMaximum(16.);
      tg_chi2ndf[*axis][*fit_function]->Draw(draw_options);
      TString legend_entry;
      if (*fit_function != "sgcl") {
        legend_entry = fit_functions_formatted[*fit_function];
      } else {
        legend_entry = "SG+p";
      }
      l_chi2ndf[*axis]->AddEntry(tg_chi2ndf[*axis][*fit_function], legend_entry, "p");
    }
    l_chi2ndf[*axis]->Draw();
    ss_figures.str("");
    ss_figures << ss_folder.str() << "figures/" << c_chi2ndf[*axis]->GetName() << ".pdf";

    c_chi2ndf[*axis]->SaveAs(TString(ss_figures.str()));
  }

// -- Sigma_x
  std::map<TString, TGraphErrors*> tg_sigma_x;
  for (vector<TString>::iterator fit_function = fit_functions.begin(); fit_function != fit_functions.end(); ++fit_function) {
    for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {

      TString tag_base = "scan";
      tag_base += *scan;
      tag_base += "_";
      tag_base += *fit_function;

      TString tgname = "tg_sigma_x_";
      tgname += tag_base;
      tg_sigma_x[tag_base] = new TGraphErrors(bcids.size());
      tg_sigma_x[tag_base]->SetName(tgname);
      tg_sigma_x[tag_base]->GetHistogram()->GetYaxis()->Set(n_bcids, 0.5, n_bcids + 0.5);
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString label_tmp = "";
        label_tmp += *bcid;
        tg_sigma_x[tag_base]->GetHistogram()->GetYaxis()->SetBinLabel(distance(bcids.begin(), bcid) + 1, label_tmp);
      }

      tg_sigma_x[tag_base]->SetMarkerStyle(marker_styles[*fit_function][*scan]);
      tg_sigma_x[tag_base]->SetMarkerColor(marker_colors[*fit_function][*scan]);
      tg_sigma_x[tag_base]->SetMarkerSize(1.5);

      Int_t point = -1;

      for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {

        point++;
        TString tag = "BCID";
        tag += *bcid;
        tag += "_";
        tag += tag_base;
        tg_sigma_x[tag_base]->SetPoint(point, m_Sigma_x[tag], point + 1. + marker_offsets[*fit_function]);
        if (*fit_function != "dg") {
          tg_sigma_x[tag_base]->SetPointError(point, m_Sigma_x_err[tag], 0.);
        }
      }
    }
  }

  // -- sigma_x canvas
  TCanvas* c_sigma_x = new TCanvas("c_sigma_x", "c_sigma_x", 1200, 800);
  TLegend* l_sigma_x = new TLegend(0.2, 0.4, 0.4, 0.6);
  l_sigma_x->SetFillColor(0);
  l_sigma_x->SetBorderSize(0);
  draw_first = true;
  for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {

    TString tag_base = "scan";
    tag_base += *scan;
    tag_base += "_";
    tag_base += "sgcl";

    TString draw_options;
    if (draw_first) {
      draw_options = "ap";
      draw_first = false;
      Float_t old_max = tg_sigma_x[tag_base]->GetHistogram()->GetXaxis()->GetXmax();
      Float_t old_min = tg_sigma_x[tag_base]->GetHistogram()->GetXaxis()->GetXmin();
      Float_t new_max = old_max + 0.5 * (old_max - old_min);
      Float_t new_min = old_min - 1.2 * (old_max - old_min);
      tg_sigma_x[tag_base]->GetHistogram()->GetXaxis()->Set(100, new_min, new_max);
      tg_sigma_x[tag_base]->GetHistogram()->GetYaxis()->Set(n_bcids, 0.5, n_bcids + 0.5);
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString label_tmp = "";
        label_tmp += *bcid;
        tg_sigma_x[tag_base]->GetHistogram()->GetYaxis()->SetBinLabel(distance(bcids.begin(), bcid) + 1, label_tmp);
      }
    } else {
      draw_options = "p";
    }
    tg_sigma_x[tag_base]->GetXaxis()->SetTitle("#Sigma_{x} (#mum)");
    tg_sigma_x[tag_base]->GetYaxis()->SetTitle("BCID");
    tg_sigma_x[tag_base]->Draw(draw_options);
    TString legend_entry = "";
    legend_entry += fit_functions_formatted["sgcl"];
    legend_entry += ", ";
    legend_entry += scans_formatted[*scan];
    l_sigma_x->AddEntry(tg_sigma_x[tag_base], legend_entry, "lp");
  }
  l_sigma_x->Draw();

  ss_figures.str("");
  ss_figures << ss_folder.str() << "figures/" << c_sigma_x->GetName() << "_" << vtx_method << p_nTrkCut << ".pdf";

  c_sigma_x->SaveAs(TString(ss_figures.str()));


  // -- sigma_y
  std::map<TString, TGraphErrors*> tg_sigma_y;
  for (vector<TString>::iterator fit_function = fit_functions.begin(); fit_function != fit_functions.end(); ++fit_function) {
    for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {

      TString tag_base = "scan";
      tag_base += *scan;
      tag_base += "_";
      tag_base += *fit_function;

      TString tgname = "tg_sigma_y_";
      tgname += tag_base;
      tg_sigma_y[tag_base] = new TGraphErrors(bcids.size());
      tg_sigma_y[tag_base]->SetName(tgname);
      tg_sigma_y[tag_base]->GetHistogram()->GetYaxis()->Set(n_bcids, 0.5, n_bcids + 0.5);
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString label_tmp = "";
        label_tmp += *bcid;
        tg_sigma_y[tag_base]->GetHistogram()->GetYaxis()->SetBinLabel(distance(bcids.begin(), bcid) + 1, label_tmp);
      }

      tg_sigma_y[tag_base]->SetMarkerStyle(marker_styles[*fit_function][*scan]);
      tg_sigma_y[tag_base]->SetMarkerColor(marker_colors[*fit_function][*scan]);
      tg_sigma_y[tag_base]->SetMarkerSize(1.5);

      Int_t point = -1;

      for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {

        point++;
        TString tag = "BCID";
        tag += *bcid;
        tag += "_";
        tag += tag_base;
        tg_sigma_y[tag_base]->SetPoint(point, m_Sigma_y[tag], point + 1. + marker_offsets[*fit_function]);
        if (*fit_function != "dg") {
          tg_sigma_y[tag_base]->SetPointError(point, m_Sigma_y_err[tag], 0.);
        }
      }
    }
  }

  // -- sigma_y canvas
  TCanvas* c_sigma_y = new TCanvas("c_sigma_y", "c_sigma_y", 1200, 800);
  TLegend* l_sigma_y = new TLegend(0.2, 0.4, 0.4, 0.6);
  l_sigma_y->SetFillColor(0);
  l_sigma_y->SetBorderSize(0);
  draw_first = true;
  for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {

    TString tag_base = "scan";
    tag_base += *scan;
    tag_base += "_";
    tag_base += "sgcl";

    TString draw_options;
    if (draw_first) {
      draw_options = "ap";
      draw_first = false;
      Float_t old_max = tg_sigma_y[tag_base]->GetHistogram()->GetXaxis()->GetXmax();
      Float_t old_min = tg_sigma_y[tag_base]->GetHistogram()->GetXaxis()->GetXmin();
      Float_t new_max = old_max + 0.5 * (old_max - old_min);
      Float_t new_min = old_min - 1.2 * (old_max - old_min);
      tg_sigma_y[tag_base]->GetHistogram()->GetXaxis()->Set(100, new_min, new_max);
      tg_sigma_y[tag_base]->GetHistogram()->GetYaxis()->Set(n_bcids, 0.5, n_bcids + 0.5);
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString label_tmp = "";
        label_tmp += *bcid;
        tg_sigma_y[tag_base]->GetHistogram()->GetYaxis()->SetBinLabel(distance(bcids.begin(), bcid) + 1, label_tmp);
      }
    } else {
      draw_options = "p";
    }
    tg_sigma_y[tag_base]->GetXaxis()->SetTitle("#Sigma_{y} (#mum)");
    tg_sigma_y[tag_base]->GetYaxis()->SetTitle("BCID");
    tg_sigma_y[tag_base]->Draw(draw_options);
    TString legend_entry = "";
    legend_entry += fit_functions_formatted["sgcl"];
    legend_entry += ", ";
    legend_entry += scans_formatted[*scan];
    l_sigma_y->AddEntry(tg_sigma_y[tag_base], legend_entry, "lp");
  }
  l_sigma_y->Draw();

  ss_figures.str("");
  ss_figures << ss_folder.str() << "figures/" << c_sigma_y->GetName() << "_" << vtx_method << p_nTrkCut << ".pdf";

  c_sigma_y->SaveAs(TString(ss_figures.str()));


  // -- c_x
  std::map<TString, TGraphErrors*> tg_c_x;
  for (vector<TString>::iterator fit_function = fit_functions.begin(); fit_function != fit_functions.end(); ++fit_function) {
    for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {

      TString tag_base = "scan";
      tag_base += *scan;
      tag_base += "_";
      tag_base += *fit_function;

      TString tgname = "tg_c_x_";
      tgname += tag_base;
      tg_c_x[tag_base] = new TGraphErrors(bcids.size());
      tg_c_x[tag_base]->SetName(tgname);
      tg_c_x[tag_base]->GetHistogram()->GetYaxis()->Set(n_bcids, 0.5, n_bcids + 0.5);
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString label_tmp = "";
        label_tmp += *bcid;
        tg_c_x[tag_base]->GetHistogram()->GetYaxis()->SetBinLabel(distance(bcids.begin(), bcid) + 1, label_tmp);
      }

      tg_c_x[tag_base]->SetMarkerStyle(marker_styles[*fit_function][*scan]);
      tg_c_x[tag_base]->SetMarkerColor(marker_colors[*fit_function][*scan]);
      tg_c_x[tag_base]->SetMarkerSize(1.5);

      Int_t point = -1;

      for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {

        point++;
        TString tag = "BCID";
        tag += *bcid;
        tag += "_";
        tag += tag_base;
        tg_c_x[tag_base]->SetPoint(point, m_c_x[tag], point + 1. + marker_offsets[*fit_function]);
        if (*fit_function != "dg") {
          tg_c_x[tag_base]->SetPointError(point, m_c_x_err[tag], 0.);
        }
      }
    }
  }



  // -- Specific luminosity
  std::map<TString, TGraphErrors*> tg_lumi_sp;

  // -- Vertexing with different fit functions
  for (vector<TString>::iterator fit_function = fit_functions.begin(); fit_function != fit_functions.end(); ++fit_function) {
    for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {

      TString tag_base = "scan";
      tag_base += *scan;
      tag_base += "_";
      tag_base += *fit_function;

      TString tgname = "tg_lumi_sp_";
      tgname += tag_base;
      tg_lumi_sp[tag_base] = new TGraphErrors(bcids.size());
      tg_lumi_sp[tag_base]->SetName(tgname);
      tg_lumi_sp[tag_base]->GetHistogram()->GetYaxis()->Set(n_bcids, 0.5, n_bcids + 0.5);
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString label_tmp = "";
        label_tmp += *bcid;
        tg_lumi_sp[tag_base]->GetHistogram()->GetYaxis()->SetBinLabel(distance(bcids.begin(), bcid) + 1, label_tmp);
      }

      tg_lumi_sp[tag_base]->SetMarkerStyle(marker_styles[*fit_function][*scan]);
      tg_lumi_sp[tag_base]->SetMarkerColor(marker_colors[*fit_function][*scan]);
      tg_lumi_sp[tag_base]->SetMarkerSize(1.5);

      Int_t point = -1;

      for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {

        point++;
        TString tag = "BCID";
        tag += *bcid;
        tag += "_";
        tag += tag_base;
        tg_lumi_sp[tag_base]->SetPoint(point, m_lumi_sp[tag], point + 1. + marker_offsets[*fit_function]);
        if (*fit_function != "dg") {
          tg_lumi_sp[tag_base]->SetPointError(point, m_lumi_sp_err[tag], 0.);
        }
      }
    }
  }

  // -- Other detectors
  std::vector<TString> other_algs;
  other_algs.push_back("bcmhor");
  other_algs.push_back("bcmvor");
  other_algs.push_back("lucidor");

  std::map<TString, TString> other_algs_formatted;
  other_algs_formatted["bcmhor"] = "BCM_HOR";
  other_algs_formatted["bcmvor"] = "BCM_VOR";
  other_algs_formatted["lucidor"] = "LUCID_OR";

  std::map<TString, Float_t> other_algs_y_offset;
  other_algs_y_offset["bcmhor"] = 0.05;
  other_algs_y_offset["bcmvor"] = 0.10;
  other_algs_y_offset["lucidor"] = -0.05;

  marker_styles["bcmhor"][1] = 22;
  marker_styles["bcmhor"][2] = 26;
  marker_styles["bcmvor"][1] = 23;
  marker_styles["bcmvor"][2] = 32;
  marker_styles["lucidor"][1] = 21;
  marker_styles["lucidor"][2] = 25;

  //marker_colors["bcmhor"][1] = kPink - 7;
  //marker_colors["bcmhor"][2] = kPink - 7;
  //marker_colors["bcmvor"][1] = kViolet - 3;
  //marker_colors["bcmvor"][2] = kViolet - 3;
  //marker_colors["lucidor"][1] = kCyan - 8;
  //marker_colors["lucidor"][2] = kCyan - 8;
  marker_colors["bcmhor"][1] = kRed;
  marker_colors["bcmhor"][2] = kRed;
  marker_colors["bcmvor"][1] = kGreen;
  marker_colors["bcmvor"][2] = kGreen;
  marker_colors["lucidor"][1] = kBlue;
  marker_colors["lucidor"][2] = kBlue;

  for (vector<TString>::iterator alg = other_algs.begin(); alg != other_algs.end(); ++alg) {
    for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {

      TString tag_base = "scan";
      tag_base += *scan;
      tag_base += "_";
      tag_base += *alg;

      TString tgname = "tg_lumi_sp_";
      tgname += tag_base;
      tg_lumi_sp[tag_base] = new TGraphErrors(bcids.size());
      tg_lumi_sp[tag_base]->SetName(tgname);
      tg_lumi_sp[tag_base]->GetHistogram()->GetYaxis()->Set(n_bcids, 0.5, n_bcids + 0.5);
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString label_tmp = "";
        label_tmp += *bcid;
        tg_lumi_sp[tag_base]->GetHistogram()->GetYaxis()->SetBinLabel(distance(bcids.begin(), bcid) + 1, label_tmp);
      }
    }
  }
  tg_lumi_sp["scan1_bcmhor"]->SetPoint(0, 56.969, 1);
  tg_lumi_sp["scan1_bcmhor"]->SetPointError(0, 0.202, 0);
  tg_lumi_sp["scan1_bcmhor"]->SetPoint(1, 52.816, 2);
  tg_lumi_sp["scan1_bcmhor"]->SetPointError(1, 0.174, 0);
  tg_lumi_sp["scan1_bcmhor"]->SetPoint(2, 52.998, 3);
  tg_lumi_sp["scan1_bcmhor"]->SetPointError(2, 0.175, 0);
  for (int i=0; i<3; i++) {
    tg_lumi_sp["scan1_bcmhor"]->SetPoint(i, tg_lumi_sp["scan1_bcmhor"]->GetX()[i] * 10, tg_lumi_sp["scan1_bcmhor"]->GetY()[i] + other_algs_y_offset["bcmhor"]);
    tg_lumi_sp["scan1_bcmhor"]->SetPointError(i, tg_lumi_sp["scan1_bcmhor"]->GetEX()[i] * 10, tg_lumi_sp["scan1_bcmhor"]->GetEY()[i]);
  }

  tg_lumi_sp["scan2_bcmhor"]->SetPoint(0, 56.236, 1);
  tg_lumi_sp["scan2_bcmhor"]->SetPointError(0, 0.202, 0);
  tg_lumi_sp["scan2_bcmhor"]->SetPoint(1, 52.460, 2);
  tg_lumi_sp["scan2_bcmhor"]->SetPointError(1, 0.174, 0);
  tg_lumi_sp["scan2_bcmhor"]->SetPoint(2, 52.325, 3);
  tg_lumi_sp["scan2_bcmhor"]->SetPointError(2, 0.174, 0);
  for (int i=0; i<3; i++) {
    tg_lumi_sp["scan2_bcmhor"]->SetPoint(i, tg_lumi_sp["scan2_bcmhor"]->GetX()[i] * 10, tg_lumi_sp["scan2_bcmhor"]->GetY()[i] + other_algs_y_offset["bcmhor"]);
    tg_lumi_sp["scan2_bcmhor"]->SetPointError(i, tg_lumi_sp["scan2_bcmhor"]->GetEX()[i] * 10, tg_lumi_sp["scan2_bcmhor"]->GetEY()[i]);
  }

  tg_lumi_sp["scan1_bcmvor"]->SetPoint(0, 56.747, 1);
  tg_lumi_sp["scan1_bcmvor"]->SetPointError(0, 0.203, 0);
  tg_lumi_sp["scan1_bcmvor"]->SetPoint(1, 52.659, 2);
  tg_lumi_sp["scan1_bcmvor"]->SetPointError(1, 0.177, 0);
  tg_lumi_sp["scan1_bcmvor"]->SetPoint(2, 52.832, 3);
  tg_lumi_sp["scan1_bcmvor"]->SetPointError(2, 0.179, 0);
  for (int i=0; i<3; i++) {
    tg_lumi_sp["scan1_bcmvor"]->SetPoint(i, tg_lumi_sp["scan1_bcmvor"]->GetX()[i] * 10, tg_lumi_sp["scan1_bcmvor"]->GetY()[i] + other_algs_y_offset["bcmvor"]);
    tg_lumi_sp["scan1_bcmvor"]->SetPointError(i, tg_lumi_sp["scan1_bcmvor"]->GetEX()[i] * 10, tg_lumi_sp["scan1_bcmvor"]->GetEY()[i]);
  }

  tg_lumi_sp["scan2_bcmvor"]->SetPoint(0, 56.180, 1);
  tg_lumi_sp["scan2_bcmvor"]->SetPointError(0, 0.205, 0);
  tg_lumi_sp["scan2_bcmvor"]->SetPoint(1, 52.413, 2);
  tg_lumi_sp["scan2_bcmvor"]->SetPointError(1, 0.165, 0);
  tg_lumi_sp["scan2_bcmvor"]->SetPoint(2, 52.178, 3);
  tg_lumi_sp["scan2_bcmvor"]->SetPointError(2, 0.163, 0);
  for (int i=0; i<3; i++) {
    tg_lumi_sp["scan2_bcmvor"]->SetPoint(i, tg_lumi_sp["scan2_bcmvor"]->GetX()[i] * 10, tg_lumi_sp["scan2_bcmvor"]->GetY()[i] + other_algs_y_offset["bcmvor"]);
    tg_lumi_sp["scan2_bcmvor"]->SetPointError(i, tg_lumi_sp["scan2_bcmvor"]->GetEX()[i] * 10, tg_lumi_sp["scan2_bcmvor"]->GetEY()[i]);
  }

  tg_lumi_sp["scan1_lucidor"]->SetPoint(0, 5.678, 1);
  tg_lumi_sp["scan1_lucidor"]->SetPointError(0, 0.006, 0);
  tg_lumi_sp["scan1_lucidor"]->SetPoint(1, 5.265, 2);
  tg_lumi_sp["scan1_lucidor"]->SetPointError(1, 0.005, 0);
  tg_lumi_sp["scan1_lucidor"]->SetPoint(2, 5.325, 3);
  tg_lumi_sp["scan1_lucidor"]->SetPointError(2, 0.005, 0);
  for (int i=0; i<3; i++) {
    tg_lumi_sp["scan1_lucidor"]->SetPoint(i, tg_lumi_sp["scan1_lucidor"]->GetX()[i] * 100, tg_lumi_sp["scan1_lucidor"]->GetY()[i] + other_algs_y_offset["lucidor"]);
    tg_lumi_sp["scan1_lucidor"]->SetPointError(i, tg_lumi_sp["scan1_lucidor"]->GetEX()[i] * 100, tg_lumi_sp["scan1_lucidor"]->GetEY()[i]);
  }

  tg_lumi_sp["scan2_lucidor"]->SetPoint(0, 5.596, 1);
  tg_lumi_sp["scan2_lucidor"]->SetPointError(0, 0.006, 0);
  tg_lumi_sp["scan2_lucidor"]->SetPoint(1, 5.224, 2);
  tg_lumi_sp["scan2_lucidor"]->SetPointError(1, 0.005, 0);
  tg_lumi_sp["scan2_lucidor"]->SetPoint(2, 5.253, 3);
  tg_lumi_sp["scan2_lucidor"]->SetPointError(2, 0.005, 0);
  for (int i=0; i<3; i++) {
    tg_lumi_sp["scan2_lucidor"]->SetPoint(i, tg_lumi_sp["scan2_lucidor"]->GetX()[i] * 100, tg_lumi_sp["scan2_lucidor"]->GetY()[i] + other_algs_y_offset["lucidor"]);
    tg_lumi_sp["scan2_lucidor"]->SetPointError(i, tg_lumi_sp["scan2_lucidor"]->GetEX()[i] * 100, tg_lumi_sp["scan2_lucidor"]->GetEY()[i]);
  }


  // -- L_sp canvas
  TCanvas* c_lumi_sp = new TCanvas("c_lumi_sp", "c_lumi_sp", 1200, 800);
  TLegend* l_lumi_sp = new TLegend(0.2, 0.3, 0.45, 0.7);
  l_lumi_sp->SetFillColor(0);
  l_lumi_sp->SetBorderSize(1);
  draw_first = true;
  for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {

    TString tag_base = "scan";
    tag_base += *scan;
    tag_base += "_";
    tag_base += "sgcl";

    TString draw_options;
    if (draw_first) {
      draw_options = "ap";
      draw_first = false;
      Float_t old_max = tg_lumi_sp[tag_base]->GetHistogram()->GetXaxis()->GetXmax();
      Float_t old_min = tg_lumi_sp[tag_base]->GetHistogram()->GetXaxis()->GetXmin();
      Float_t new_max = old_max + 0.4 * (old_max - old_min);
      Float_t new_min = old_min - 1.5 * (old_max - old_min);
      tg_lumi_sp[tag_base]->GetHistogram()->GetXaxis()->Set(100, new_min, new_max);
      tg_lumi_sp[tag_base]->GetHistogram()->GetYaxis()->Set(n_bcids, 0.5, n_bcids + 0.5);
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString label_tmp = "";
        label_tmp += *bcid;
        tg_lumi_sp[tag_base]->GetHistogram()->GetYaxis()->SetBinLabel(distance(bcids.begin(), bcid) + 1, label_tmp);
      }
    } else {
      draw_options = "p";
    }
    tg_lumi_sp[tag_base]->GetXaxis()->SetTitle("L_{sp} (#mub^{-1}s^{-1})");
    tg_lumi_sp[tag_base]->GetYaxis()->SetTitle("BCID");
    tg_lumi_sp[tag_base]->Draw(draw_options);
    TString legend_entry = "Vertexing, ";
    legend_entry += scans_formatted[*scan];
    l_lumi_sp->AddEntry(tg_lumi_sp[tag_base], legend_entry, "p");

    for (vector<TString>::iterator alg = other_algs.begin(); alg != other_algs.end(); ++alg) {

      TString other_alg_tag = "scan";
      other_alg_tag += *scan;
      other_alg_tag += "_";
      other_alg_tag += *alg;
      tg_lumi_sp[other_alg_tag]->SetMarkerStyle(marker_styles[*alg][*scan]);
      tg_lumi_sp[other_alg_tag]->SetMarkerColor(marker_colors[*alg][*scan]);
      tg_lumi_sp[other_alg_tag]->SetMarkerSize(1.5);
      tg_lumi_sp[other_alg_tag]->Draw("p");
      legend_entry = other_algs_formatted[*alg];
      legend_entry += ", ";
      legend_entry += scans_formatted[*scan];
      l_lumi_sp->AddEntry(tg_lumi_sp[other_alg_tag], legend_entry, "p");
    }
  }
  l_lumi_sp->Draw();

  ss_figures.str("");
  ss_figures << ss_folder.str() << "figures/" << c_lumi_sp->GetName() << "_" << vtx_method << p_nTrkCut << ".pdf";

  c_lumi_sp->SaveAs(TString(ss_figures.str()));

  // -- c_x canvas
  TCanvas* c_c_x = new TCanvas("c_c_x", "c_c_x", 1200, 800);
  c_c_x->SetRightMargin(0.1);
  TLegend* l_c_x = new TLegend(0.2, 0.4, 0.4, 0.6);
  l_c_x->SetFillColor(0);
  l_c_x->SetBorderSize(0);
  draw_first = true;
  for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {

    TString tag_base = "scan";
    tag_base += *scan;
    tag_base += "_";
    tag_base += "sgcl";

    TString draw_options;
    if (draw_first) {
      draw_options = "ap";
      draw_first = false;
      Float_t old_max = tg_c_x[tag_base]->GetHistogram()->GetXaxis()->GetXmax();
      Float_t old_min = tg_c_x[tag_base]->GetHistogram()->GetXaxis()->GetXmin();
      Float_t new_max = old_max + 0.5 * (old_max - old_min);
      Float_t new_min = old_min - 1.2 * (old_max - old_min);
      tg_c_x[tag_base]->GetHistogram()->GetXaxis()->Set(100, new_min, new_max);
      tg_c_x[tag_base]->GetHistogram()->GetYaxis()->Set(n_bcids, 0.5, n_bcids + 0.5);
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString label_tmp = "";
        label_tmp += *bcid;
        tg_c_x[tag_base]->GetHistogram()->GetYaxis()->SetBinLabel(distance(bcids.begin(), bcid) + 1, label_tmp);
      }
    } else {
      draw_options = "p";
    }
    tg_c_x[tag_base]->GetXaxis()->SetTitle("c_{x}");
    tg_c_x[tag_base]->GetYaxis()->SetTitle("BCID");
    tg_c_x[tag_base]->Draw(draw_options);
    TString legend_entry = "";
    legend_entry += fit_functions_formatted["sgcl"];
    legend_entry += ", ";
    legend_entry += scans_formatted[*scan];
    l_c_x->AddEntry(tg_c_x[tag_base], legend_entry, "lp");
  }
  l_c_x->Draw();

  ss_figures.str("");
  ss_figures << ss_folder.str() << "figures/" << c_c_x->GetName() << "_" << vtx_method << p_nTrkCut << ".pdf";

  c_c_x->SaveAs(TString(ss_figures.str()));


  // -- c_y
  std::map<TString, TGraphErrors*> tg_c_y;
  for (vector<TString>::iterator fit_function = fit_functions.begin(); fit_function != fit_functions.end(); ++fit_function) {
    for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {

      TString tag_base = "scan";
      tag_base += *scan;
      tag_base += "_";
      tag_base += *fit_function;

      TString tgname = "tg_c_y_";
      tgname += tag_base;
      tg_c_y[tag_base] = new TGraphErrors(bcids.size());
      tg_c_y[tag_base]->SetName(tgname);
      tg_c_y[tag_base]->GetHistogram()->GetYaxis()->Set(n_bcids, 0.5, n_bcids + 0.5);
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString label_tmp = "";
        label_tmp += *bcid;
        tg_c_y[tag_base]->GetHistogram()->GetYaxis()->SetBinLabel(distance(bcids.begin(), bcid) + 1, label_tmp);
      }

      tg_c_y[tag_base]->SetMarkerStyle(marker_styles[*fit_function][*scan]);
      tg_c_y[tag_base]->SetMarkerColor(marker_colors[*fit_function][*scan]);
      tg_c_y[tag_base]->SetMarkerSize(1.5);

      Int_t point = -1;

      for (std::vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {

        point++;
        TString tag = "BCID";
        tag += *bcid;
        tag += "_";
        tag += tag_base;
        tg_c_y[tag_base]->SetPoint(point, m_c_y[tag], point + 1. + marker_offsets[*fit_function]);
        if (*fit_function != "dg") {
          tg_c_y[tag_base]->SetPointError(point, m_c_y_err[tag], 0.);
        }
      }
    }
  }

  // -- c_y canvas
  TCanvas* c_c_y = new TCanvas("c_c_y", "c_c_y", 1200, 800);
  c_c_y->SetRightMargin(0.1);
  TLegend* l_c_y = new TLegend(0.2, 0.4, 0.4, 0.6);
  l_c_y->SetFillColor(0);
  l_c_y->SetBorderSize(0);
  draw_first = true;
  for (std::vector<Int_t>::iterator scan = scans.begin(); scan != scans.end(); ++scan) {

    TString tag_base = "scan";
    tag_base += *scan;
    tag_base += "_";
    tag_base += "sgcl";

    TString draw_options;
    if (draw_first) {
      draw_options = "ap";
      draw_first = false;
      Float_t old_max = tg_c_y[tag_base]->GetHistogram()->GetXaxis()->GetXmax();
      Float_t old_min = tg_c_y[tag_base]->GetHistogram()->GetXaxis()->GetXmin();
      Float_t new_max = old_max + 0.5 * (old_max - old_min);
      Float_t new_min = old_min - 1.2 * (old_max - old_min);
      tg_c_y[tag_base]->GetHistogram()->GetXaxis()->Set(100, new_min, new_max);
      tg_c_y[tag_base]->GetHistogram()->GetYaxis()->Set(n_bcids, 0.5, n_bcids + 0.5);
      for (vector<Int_t>::iterator bcid = bcids.begin(); bcid != bcids.end(); ++bcid) {
        TString label_tmp = "";
        label_tmp += *bcid;
        tg_c_y[tag_base]->GetHistogram()->GetYaxis()->SetBinLabel(distance(bcids.begin(), bcid) + 1, label_tmp);
      }
    } else {
      draw_options = "p";
    }
    tg_c_y[tag_base]->GetXaxis()->SetTitle("c_{y}");
    tg_c_y[tag_base]->GetYaxis()->SetTitle("BCID");
    tg_c_y[tag_base]->Draw(draw_options);
    TString legend_entry = "";
    legend_entry += fit_functions_formatted["sgcl"];
    legend_entry += ", ";
    legend_entry += scans_formatted[*scan];
    l_c_y->AddEntry(tg_c_y[tag_base], legend_entry, "lp");
  }
  l_c_y->Draw();

  ss_figures.str("");
  ss_figures << ss_folder.str() << "figures/" << c_c_y->GetName() << "_" << vtx_method << p_nTrkCut << ".pdf";

  c_c_y->SaveAs(TString(ss_figures.str()));


}
