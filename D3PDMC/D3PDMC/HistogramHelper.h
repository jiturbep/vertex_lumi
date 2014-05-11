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
                          bool save=true, TString save_name="");
    static TH1I* defineIHistogram(TString name,
                           unsigned int bins, Int_t lbin, Int_t ubin,
                           TString xaxistitle, TString yaxistitle,
                           bool save=true, TString save_name="");
    static TH2F* define2DHistogram(TString name, unsigned
                            int xbins, Double_t lxbin, Double_t uxbin,
                            unsigned int ybins, Double_t lybin, Double_t uybin,
                            TString xaxistitle, TString yaxistitle,
                            bool save=true, TString save_name="");
    static TH2F* define2DHistogram(TString name, unsigned
                            int xbins, Double_t lxbin, Double_t uxbin,
                            unsigned int ybins, Double_t lybin, Double_t uybin,
                            TString xaxistitle, TString yaxistitle,
                            TH1F * &h_mean, TH1F * &h_rms, //created to store mean and RMS of the 2D histogram
                            bool save=true, TString save_name="");
    static TH3F* define3DHistogram(TString name,
                            unsigned int xbins, Double_t lxbin, Double_t uxbin,
                            unsigned int ybins, Double_t lybin, Double_t uybin,
                            unsigned int zbins, Double_t lzbin, Double_t uzbin,
                            TString xaxistitle, TString yaxistitle, TString zaxistitle,
                            bool save, TString save_name);

   static std::vector<TH1*> histQueue;

};

#endif
