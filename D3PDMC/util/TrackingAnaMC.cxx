// tope level executable
#include <iostream>
#include <fstream>
#include <sstream>
#include <getopt.h>

#include "D3PDMC/AnaVtxTree.h"
#include "TApplication.h"
#include "TFile.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TCollection.h"
#include "TObjString.h"
#include "TLegend.h"

using namespace std;

// --- Usage
void usage(const char *pname);
// --- Atlas style.. with some custom, lod =
//     0: ATLAS style
//     1: Add some more info
//     2: Add fitting info too
void SetAtlasStyle(int lod=0);

// --- global debug
int  verbose_debug;
int debugLevel;
double cfg_effCriteriaZ;
double cfg_effCriteriaT;
string outputFileName;
Long64_t maxEvents;
string cfg_qualityVertexVersion;

int main(int argc, char **argv) {
  cout << "Starting TrackingAnaMC.exe" << endl;
  // --- Create main ROOT application
  TApplication theApp("App", &argc, argv);

  // --- global variables init
  verbose_debug = 0;
  cfg_effCriteriaZ = 0.3; //300 um default Z matching for vertex selection eff calculation
  cfg_effCriteriaT = 5; //5sigma default T matching for vertex selection eff calculation
  outputFileName.clear();
  maxEvents = -1; //default = all events
  debugLevel = 0; //no debug
  cfg_qualityVertexVersion = "";
  Bool_t cfg_interactiveROOT=false;
  TString triggerFilter;
  float cfg_minMuFilter(-1);
  float cfg_maxMuFilter(-1);
  bool skip_hs(false);
  bool do_timing(false);
  Float_t max_chi2ndf(-1);

  // --- Scan for custom command-line parameters (ROOT reserves: -l,-b,-n,-q)
  int c;
  extern int optind;
  extern char* optarg;
  vector<string> inputRootFilesTxt; //input Root Files from a text file
  while (1) {
    int option_index = 0;
    static struct option long_options[] = {
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "pahDi:o:e:z:d:t:v:Im:M:c:",
                     long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {
      case 'p': {
        do_timing = true;
        break;
      }
      case 'a': {
        skip_hs = true;
        break;
      }
      case 'i': {
        string inputFileName = optarg;
        if (inputFileName[0] == '#') {
          break;  //input file commented out
        }
        ifstream inputFileStream(inputFileName.c_str());
        if (inputFileStream.fail()) {
          cerr << "ERROR: Error opening input file " << inputFileName << endl;
          return 0;
        }
        cout << "Acquiring input files from " << inputFileName << endl;
        try {
          while (!inputFileStream.eof()) {
            string currentRootInFile;
            inputFileStream >> currentRootInFile;
            if (!currentRootInFile.empty()) {
              inputRootFilesTxt.push_back(currentRootInFile);
            }
          }
        } catch ( char* errstr) {
          cerr << "ERROR: I/O Error while retrieving ROOT file list from " << inputFileName << endl;
          return 1;
        }
        break;
      }
      case 'o': {
        //set output file name
        outputFileName = optarg;
        break;
      }
      case 'e': {
        //get number of events to be processed
        stringstream maxEventsStream;
        maxEventsStream.str(optarg);
        try {
          maxEventsStream >> maxEvents;
        } catch (const char *) {
          cerr << "ERROR: Invalid number of events: " << optarg << endl;
          return 1;
        }
        break;
      }
      case 'z': {
        //get Z interval for selection efficiency evaluation
        stringstream effCriteriaStream;
        effCriteriaStream.str(optarg);
        try {
          effCriteriaStream >> cfg_effCriteriaZ;
        } catch (const char *) {
          cerr << "ERROR: Invalid Z-efficiency criteria: " << optarg << endl;
          return 1;
        }
        break;
      }
      case 't': {
        //get trigger filter
        triggerFilter = optarg;
        break;
      }
      case 'v': {
        cfg_qualityVertexVersion = optarg;
        break;
      }
      case 'm': {
        //get minimum mu filter requirement
        stringstream muFilterStream;
        muFilterStream.str(optarg);
        try {
          muFilterStream >> cfg_minMuFilter;
        } catch (const char *) {
          cerr << "ERROR: Invalid min mu filter requirement: " << optarg << endl;
          return 1;
        }
        break;
      }
      case 'M': {
        //get maximum mu filter requirement
        stringstream muFilterStream;
        muFilterStream.str(optarg);
        try {
          muFilterStream >> cfg_maxMuFilter;
        } catch (const char *) {
          cerr << "ERROR: Invalid min mu filter requirement: " << optarg << endl;
          return 1;
        }
        break;
      }
      case 'D': {
        cout << "Enabling verbose debugging." << endl;
        verbose_debug = 1;
        break;
      }
      case 'd': {
        //get debug level
        stringstream debugLevelStream;
        debugLevelStream.str(optarg);
        try {
          debugLevelStream >> debugLevel;
        } catch (const char *) {
          cerr << "ERROR: Invalid debug Level: " << optarg << endl;
          return 1;
        }
        break;
      }
      case 'I': {
        cfg_interactiveROOT=true;
        break;
      }
      case 'c': {
        stringstream ss;
        ss << optarg;
        try {
          ss >> max_chi2ndf;
        } catch (const char *) {
          cerr << "[TrackingAna] ERROR : Specified max chi2/ndf not valid: " << optarg << endl;
          return 1;
          break;
        }
        break;
      }
      default:
        usage((const char*)argv[0]);
        return 0;
        break;
    }
  } // end scan of command-line parameters

  //print warning for any unused parameter
  if (optind < argc) {
    cout << "WARNING: Discarded parameters: ";
    while (optind < argc) {
      cout << argv[optind++] << " ";
    }
    cout << endl;
  }


  // --- Set style
  SetAtlasStyle(2);

  // --- init chain
  TChain* vertexChain = new TChain("VtxTree");
  TChain* triggerChain = new TChain("VtxTreeMeta/TrigConfTree");

  // --- set input files
  //now add files from files recognized by ROOT TApplication
  TIter nextInputFile(theApp.InputFiles());
  TObjString *InputFile;
  while ( (InputFile = (TObjString*)nextInputFile()) ) {
    inputRootFilesTxt.push_back(InputFile->GetString().Data());
  }

  //Print total number of input files and file list (if verbose) and add them to the chain
  cout << "Number of input files: " << inputRootFilesTxt.size() << endl;
  for (vector<string>::iterator irf = inputRootFilesTxt.begin(); irf != inputRootFilesTxt.end(); ++irf) {
    //cout << "Input file: " << *irf << endl;
    vertexChain->Add(irf->c_str());
    triggerChain->Add(irf->c_str());
  }
  cout << "Chain created. Processing.." << endl;
  cout << "Total number of events: " << vertexChain->GetEntries() << endl;

  if (inputRootFilesTxt.size() == 0) {
    cerr << "ERROR: No ROOT input files supplied." << endl;
    usage((const char*)argv[0]);
    return 1;
  }

  // --- output file (name prefix without extension)
  // - if provided from input, check it's not with .root
  //AnaVtxTree with add .root and .txt extensions
  if (outputFileName.empty()) {
    //by default take first input file and add postfix
    outputFileName = inputRootFilesTxt[0];
    outputFileName.replace(outputFileName.end()-5, outputFileName.end(), "_results");
  } else {
    //if ends with ".root" cut it, same for ".txt"
    size_t found = outputFileName.rfind(".root");
    if (found != -1 && found == outputFileName.length()-1-5) {
      outputFileName.erase(found);
    }
    found = outputFileName.rfind(".txt");
    if (found != -1 && found == outputFileName.length()-1-4) {
      outputFileName.erase(found);
    }
  }

  // --- Create TSelector instance, setup it and initiate the Loop
  AnaVtxTree *vertexTreeSelector = new AnaVtxTree();
  vertexTreeSelector->SetEffCriteriaZ(cfg_effCriteriaZ);
  //  vertexTreeSelector->SetEffCriteriaT(cfg_effCriteriaT);
  vertexTreeSelector->SetOutputFileName(outputFileName);
  vertexTreeSelector->SetQualityVertexVersion(cfg_qualityVertexVersion);
  if (triggerFilter != "") {
    vertexTreeSelector->SetTriggerFilter(triggerFilter, triggerChain);
  }
  if (cfg_minMuFilter >=0 || cfg_maxMuFilter >=0) {
    vertexTreeSelector->SetMuFilter(cfg_minMuFilter, cfg_maxMuFilter);
  }
  std::cout << "max_chi2ndf: " << max_chi2ndf << std::endl;
  if (max_chi2ndf > 0.) {
    vertexTreeSelector->SetMaxChi2Ndf(max_chi2ndf);
  }

  if (verbose_debug) {
    vertexTreeSelector->SetDebugLevel(100);
  }
  if (debugLevel != 0) {
    vertexTreeSelector->SetDebugLevel(debugLevel);
  }
  vertexTreeSelector->SetHS(skip_hs);
  if (do_timing) {
    vertexTreeSelector->DoTiming();
  }
  cout << "maxEvents = " << maxEvents << endl;
  if (maxEvents < 0) {
    vertexChain->Process(dynamic_cast<TSelector*>(vertexTreeSelector));
  } else {
    vertexChain->Process(dynamic_cast<TSelector*>(vertexTreeSelector), "", maxEvents);
  }

  if (cfg_interactiveROOT) {
    theApp.Run();
  }

  return 0;
}

// the rest is taken from the ATLAS style wiki, and adjusted that it actually compiles ...
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
  TStyle* atlasStyle = AtlasStyle(lod);
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();
}

void usage(const char *pname) {
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << pname << " [options] [-i InputFileName | InputFile1.root [InputFile2.root [...]]]" << std::endl;
  std::cout << "\t i: Uses files listed in InputFileName as input." << std::endl;
  std::cout << "List of other options:" << std::endl;
  std::cout << "\t D  Print verbose debug(debugLevel=100)." << std::endl;
  std::cout << "\t d: Set debug level." << std::endl;
  std::cout << "\t o: Set output prefix (.root and .txt will be added)." << std::endl;
  std::cout << "\t e: Set number of events to process." << std::endl;
  std::cout << "\t z: |Z| interval for selection efficiency evaluation." << std::endl;
  std::cout << "\t v: Quality of reconstructed vertices (version.modifier.parameter)" << std::endl;
  std::cout << "\t t: Set trigger requirement" << std::endl;
  std::cout << "\t m: Set minimum mu (actualIntPerXing) filter requirement" << std::endl;
  std::cout << "\t M: Set maximum mu (actualIntPerXing) filter requirement" << std::endl;
  std::cout << "\t I  Runs interactive ROOT session, do not quit after it ends." << std::endl;
  std::cout << "\t h: Print this screen." << std::endl;
  std::cout << std::endl;
}
