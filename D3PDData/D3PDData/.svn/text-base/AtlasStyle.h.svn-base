#ifndef ATLAS_STYLE_H
#define ATLAS_STYLE_H

#include "TStyle.h"
#include "TROOT.h"

// the rest is taken from the ATLAS style wiki, and adjusted that it actually compiles ...
// --- Atlas style.. with some custom, lod =
//     0: ATLAS style
//     1: Add some more info
//     2: Add fitting info too
TStyle* AtlasStyle(int lod=0) {
  TStyle *atlasStyle = new TStyle("ATLAS","Atlas style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  atlasStyle->SetFrameBorderMode(icol);
  atlasStyle->SetFrameFillColor(icol);
  atlasStyle->SetCanvasBorderMode(icol);
  atlasStyle->SetCanvasColor(icol);
  atlasStyle->SetPadBorderMode(icol);
  atlasStyle->SetPadColor(icol);
  atlasStyle->SetStatColor(icol);
  //atlasStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  atlasStyle->SetPaperSize(20,26);

  // set margin sizes
  atlasStyle->SetPadTopMargin(0.05);
  atlasStyle->SetPadRightMargin(0.05);
  atlasStyle->SetPadBottomMargin(0.16);
  atlasStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  atlasStyle->SetTitleXOffset(1.4);
  atlasStyle->SetTitleYOffset(1.4);

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.05;
  atlasStyle->SetTextFont(font);

  atlasStyle->SetTextSize(tsize);
  atlasStyle->SetLabelFont(font,"x");
  atlasStyle->SetTitleFont(font,"x");
  atlasStyle->SetLabelFont(font,"y");
  atlasStyle->SetTitleFont(font,"y");
  atlasStyle->SetLabelFont(font,"z");
  atlasStyle->SetTitleFont(font,"z");

  atlasStyle->SetLabelSize(tsize,"x");
  atlasStyle->SetTitleSize(tsize,"x");
  atlasStyle->SetLabelSize(tsize,"y");
  atlasStyle->SetTitleSize(tsize,"y");
  atlasStyle->SetLabelSize(tsize,"z");
  atlasStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
  //S. Pagan Griso likes small dots for busy plots
  //atlasStyle->SetMarkerStyle(20);
  atlasStyle->SetMarkerStyle(8);
  atlasStyle->SetMarkerSize(1.2);
  atlasStyle->SetHistLineWidth(2);
  atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars
  //atlasStyle->SetErrorX(0.001);
  // get rid of error bar caps
  atlasStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  atlasStyle->SetOptTitle(0);

  // ATLAS default
  atlasStyle->SetOptStat(0);
  // A. Wildauer: change this
  if (lod == 1) {
    atlasStyle->SetOptStat(11110);
  }
  // S. Pagan Griso likes more verbosity
  if (lod == 2) {
    atlasStyle->SetOptStat(111111);
  }
  atlasStyle->SetStatW(0.3);
  atlasStyle->SetStatH(0.2);
  ///////////

  atlasStyle->SetOptFit(0);
  if (lod >=2) {
    atlasStyle->SetOptFit(1111);
  }

  // put tick marks on top and RHS of plots
  atlasStyle->SetPadTickX(1);
  atlasStyle->SetPadTickY(1);

  atlasStyle->SetPalette(1);

  return atlasStyle;

}

void SetAtlasStyle (int lod) {
  static TStyle* atlasStyle = 0;
  if( atlasStyle == 0 ) { atlasStyle = AtlasStyle(lod); }
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();
}

#endif
