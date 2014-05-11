#ifndef HistogramHelper_h
#define HistogramHelper_h
// TODO delete hists that aren't saved to file

// STL
#include <vector>

//#include <TROOT.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

class HistogramHelper {

public:
    
    static TH1F* defineHistogram(TString name,
                                 unsigned int bins, Double_t lbin, Double_t ubin,
                                 TString xaxistitle, TString yaxistitle,
                                 bool save=true, TString title="");

    static TH2F* define2DHistogram(TString name, Int_t xbins, Double_t lxbin, Double_t uxbin,
                                   Int_t ybins, Double_t lybin, Double_t uybin,
                                   TString xaxistitle, TString yaxistitle,
                                   bool save=true, TString title="");
          
    static TH3F* define3DHistogram(TString name, unsigned
                                   int xbins, Double_t lxbin, Double_t uxbin,
                                   unsigned int ybins, Double_t lybin, Double_t uybin,
                                   unsigned int zbins, Double_t lzbin, Double_t uzbin,
                                   TString xaxistitle, TString yaxistitle, TString zaxistitle,
                                   bool save=true, TString title="");

   static std::vector<TH1*> histQueue;

};

#endif
