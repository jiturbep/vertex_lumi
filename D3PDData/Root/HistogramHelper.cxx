#include "D3PDData/HistogramHelper.h"

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

// 2D hist
TH2F* HistogramHelper::define2DHistogram(TString name,
                                    Int_t xbins, Double_t lxbin, Double_t uxbin,
                                    Int_t ybins, Double_t lybin, Double_t uybin,
                                    TString xaxistitle, TString yaxistitle,
                                    bool save, TString title) {

  TH2F *histogram = new TH2F( (title==""?name:title), name, xbins, lxbin, uxbin, ybins, lybin, uybin );

  TAxis* xaxis = histogram->GetXaxis();
  xaxis->SetTitle(xaxistitle);
  TAxis* yaxis = histogram->GetYaxis();
  yaxis->SetTitle(yaxistitle);
  if (save) {
    std::string newname;
    if (title == "") { //use name as default save name
      newname = name.Data();
    } else {
      newname = title.Data();
    }
//    storageHistQueue.insert( std::pair<std::string, std::pair<std::string, TH1*> >(newname,
//                             std::pair<std::string, TH1* >(std::string("TH2F*"), dynamic_cast<TH1*>(histogram))) );
    histQueue.push_back( dynamic_cast<TH1*>(histogram) );
  }
  return histogram;
}

// 3D hist
TH3F* HistogramHelper::define3DHistogram(TString name,
                                    unsigned int xbins, Double_t lxbin, Double_t uxbin,
                                    unsigned int ybins, Double_t lybin, Double_t uybin,
                                    unsigned int zbins, Double_t lzbin, Double_t uzbin,
                                    TString xaxistitle, TString yaxistitle, TString zaxistitle,
                                    bool save, TString title) {
  TH3F* histogram;
  if (title == "") {
    histogram = new TH3F(name, name, xbins, lxbin, uxbin, ybins, lybin, uybin, zbins, lzbin, uzbin);
  } else {
    histogram = new TH3F(title, name, xbins, lxbin, uxbin, ybins, lybin, uybin, zbins, lzbin, uzbin);
  }

  TAxis* xaxis = histogram->GetXaxis();
  xaxis->SetTitle(xaxistitle);
  TAxis* yaxis = histogram->GetYaxis();
  yaxis->SetTitle(yaxistitle);
  TAxis* zaxis = histogram->GetZaxis();
  zaxis->SetTitle(zaxistitle);

  if (save) {
    std::string newname;
    if (title == "") { //use name as default save name
      newname = name.Data();
    } else {
      newname = title.Data();
    }

//    storageHistQueue.insert( std::pair<std::string, std::pair<std::string, TH1*> >(newname,
//                             std::pair<std::string, TH1* >(std::string("TH3F*"), dynamic_cast<TH1*>(histogram))) );
    histQueue.push_back( dynamic_cast<TH1*>(histogram) );
  }
  return histogram;
}

