#include "D3PDMC/HistogramHelper.h"

// Instantiate hist queue
std::vector<TH1*> HistogramHelper::histQueue;

// 1D hist
TH1F* HistogramHelper::defineHistogram(TString name,
                                  unsigned int bins, Double_t lbin, Double_t ubin,
                                  TString xaxistitle, TString yaxistitle,
                                  bool save, TString title) {

  TH1F *histogram = new TH1F( (title==""?name:title), name, bins, lbin, ubin );

  TAxis* xaxis = histogram->GetXaxis();
  xaxis->SetTitle(xaxistitle);
  TAxis* yaxis = histogram->GetYaxis();
  yaxis->SetTitle(yaxistitle);
  //Add to the hist queue to be saved
  if (save) {
    std::string newname;
    if (title == "") { //use name as default save name
      newname = name.Data();
    } else {
      newname = title.Data();
    }

//    storageHistQueue.insert( std::pair<std::string, std::pair<std::string, TH1*> >(newname,
//                             std::pair<std::string, TH1* >(std::string("TH1F*"), dynamic_cast<TH1*>(histogram))) );
    histQueue.push_back( dynamic_cast<TH1*>(histogram) );
  }

  return histogram;
}

TH1I* HistogramHelper::defineIHistogram(TString name,
                                   unsigned int bins, Int_t lbin, Int_t ubin,
                                   TString xaxistitle, TString yaxistitle,
                                   bool save, TString save_name) {
                                   
  TH1I *histogram = new TH1I( (save_name==""?name:save_name), name, bins, lbin, ubin );

  TAxis* xaxis = histogram->GetXaxis();
  xaxis->SetTitle(xaxistitle);
  TAxis* yaxis = histogram->GetYaxis();
  yaxis->SetTitle(yaxistitle);
  //Add to the hist queue to be saved
  if (save) {
    std::string newname;
    if (save_name == "") { //use name as default save name
      newname = name.Data();
    } else {
      newname = save_name.Data();
    }

    //storageHistQueue.insert( pair<string, pair<string, TH1*> >(newname,
      //                       pair<string, TH1* >(string("TH1I*"), dynamic_cast<TH1*>(histogram))) );
      histQueue.push_back( dynamic_cast<TH1*>(histogram) );
  }

  return histogram;
}


TH2F* HistogramHelper::define2DHistogram(TString name,
                                    unsigned int xbins, Double_t lxbin, Double_t uxbin,
                                    unsigned int ybins, Double_t lybin, Double_t uybin,
                                    TString xaxistitle, TString yaxistitle,
                                    bool save, TString save_name) {
  
  TH2F *histogram = new TH2F( (save_name==""?name:save_name), name, xbins, lxbin, uxbin, ybins, lybin, uybin );

  TAxis* xaxis = histogram->GetXaxis();
  xaxis->SetTitle(xaxistitle);
  TAxis* yaxis = histogram->GetYaxis();
  yaxis->SetTitle(yaxistitle);

  if (save) {
    std::string newname;
    if (save_name == "") { //use name as default save name
      newname = name.Data();
    } else {
      newname = save_name.Data();
    }

    //storageHistQueue.insert( pair<string, pair<string, TH1*> >(newname,
      //                       pair<string, TH1* >(string("TH2F*"), dynamic_cast<TH1*>(histogram))) );
    histQueue.push_back( dynamic_cast<TH1*>(histogram) );
  }
  return histogram;
}

TH2F* HistogramHelper::define2DHistogram(TString name,
                                    unsigned int xbins, Double_t lxbin, Double_t uxbin,
                                    unsigned int ybins, Double_t lybin, Double_t uybin,
                                    TString xaxistitle, TString yaxistitle,
                                    TH1F* &h_mean, TH1F* &h_rms,
                                    bool save, TString save_name) {
  // Create 2D histogram and 1D mean/rms
  TH2F *h_main = define2DHistogram(name, xbins, lxbin, uxbin, ybins, lybin, uybin, xaxistitle, yaxistitle, save, save_name);
  TString hTitle;
  TString hSaveName;
  hTitle = name + TString(" - Mean");
  hSaveName = save_name + TString("Mean");
  h_mean = defineHistogram(hTitle, xbins, lxbin, uxbin, xaxistitle, "Mean", save, hSaveName);
  hTitle = name + TString(" - RMS");
  hSaveName = save_name + TString("Rms");
  h_rms = defineHistogram(hTitle, xbins, lxbin, uxbin, xaxistitle, "RMS", save, hSaveName);

  return h_main;
}

//3D histogram
TH3F* HistogramHelper::define3DHistogram(TString name,
                                    unsigned int xbins, Double_t lxbin, Double_t uxbin,
                                    unsigned int ybins, Double_t lybin, Double_t uybin,
                                    unsigned int zbins, Double_t lzbin, Double_t uzbin,
                                    TString xaxistitle, TString yaxistitle, TString zaxistitle,
                                    bool save, TString save_name) {

  TH3F *histogram = new TH3F( (save_name==""?name:save_name), name, xbins, lxbin, uxbin, ybins, lybin, uybin, zbins, lzbin, uzbin );

  TAxis* xaxis = histogram->GetXaxis();
  xaxis->SetTitle(xaxistitle);
  TAxis* yaxis = histogram->GetYaxis();
  yaxis->SetTitle(yaxistitle);
  TAxis* zaxis = histogram->GetZaxis();
  zaxis->SetTitle(zaxistitle);

  if (save) {
    std::string newname;
    if (save_name == "") { //use name as default save name
      newname = name.Data();
    } else {
      newname = save_name.Data();
    }

    //storageHistQueue.insert( pair<string, pair<string, TH1*> >(newname,
      //                       pair<string, TH1* >(string("TH3F*"), dynamic_cast<TH1*>(histogram))) );
      histQueue.push_back( dynamic_cast<TH1*>(histogram) );
  }
  return histogram;
}

TProfile* HistogramHelper::defineProfile(TString name,
                                    unsigned int xbins, Double_t lxbin, Double_t uxbin,
                                    Double_t lybin, Double_t uybin,
                                    TString xaxistitle, TString yaxistitle,
                                    bool save, TString save_name) {
  
  TProfile *histogram = new TProfile( (save_name==""?name:save_name), name, xbins, lxbin, uxbin, lybin, uybin );

  TAxis* xaxis = histogram->GetXaxis();
  xaxis->SetTitle(xaxistitle);
  TAxis* yaxis = histogram->GetYaxis();
  yaxis->SetTitle(yaxistitle);

  if (save) {
    std::string newname;
    if (save_name == "") { //use name as default save name
      newname = name.Data();
    } else {
      newname = save_name.Data();
    }

    //storageHistQueue.insert( pair<string, pair<string, TH1*> >(newname,
      //                       pair<string, TH1* >(string("TH2F*"), dynamic_cast<TH1*>(histogram))) );
    histQueue.push_back( dynamic_cast<TH1*>(histogram) );
  }
  return histogram;
}


