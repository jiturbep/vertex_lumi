// Vertex counting for data ntuples

// Local
#include "D3PDData/AnaVtxTree.h"
#include "D3PDData/AtlasStyle.h"

// STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <string>
#include <stdio.h>

// ROOT
#include "TApplication.h"
#include "TChain.h"
#include "TObjString.h"
#include "TSystem.h"


// --- global debug
int  verbose_debug;
std::string outputFileName;
Long64_t maxEvents;
std::string asciiDumpFile;

// --- Usage
void usage(const char *pname) {
  std::cout << pname << " [-h] [options] [InputFile1.root [InputFile2.root [...]]]" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "\t e MaxEvents: Set number of events to process." << std::endl;
  std::cout << "\t o OutputPrefix: Set output prefix (.root and .txt will be added)." << std::endl;
  std::cout << "\t v: Print verbose debug." << std::endl;
  std::cout << "\t t TriggerChain: Select events passing TriggerChain." << std::endl;
  std::cout << "\t p: Enable mode for regular physics runs." << std::endl;
  std::cout << "\t c chi2Cut: Apply cut on chi2/ndf of reconstructed vertices." << std::endl;
  std::cout << "\t h: Print this screen." << std::endl;
}

int main(int argc, char **argv) {
  // --- Create main ROOT application
  TApplication theApp("App", &argc, argv);

  // --- global variables init
  verbose_debug = 0;
  outputFileName = "";
  maxEvents = -1; //default = all events
  asciiDumpFile = "";
  std::vector<std::string> inputRootFilesTxt;
  TString triggerFilter("");
  bool physics_run(false);
  Float_t max_chi2ndf(0.);

  // --- Scan for custom command-line parameters (ROOT reserves: -l,-b,-n,-q)
  int c;
  extern int optind;
  extern char* optarg;
  while (1) {
    int option_index = 0;
    static struct option long_options[] = { {0,0,0,0} };
    c=getopt_long(argc, argv, "pvhe:o:t:",long_options,&option_index);
    if (c == -1) {
      break;
    }
    switch(c) {
      case 'e': {
        std::stringstream maxEventsStream;
        maxEventsStream.str(optarg);
        try {
          maxEventsStream >> maxEvents;
        } catch (const char *) {
          std::cerr << "[TrackingAnaData] ERROR: Invalid number of events: " << optarg << std::endl;
          return 1;
          break;
        }
        break;
      }
      case 'o': {
        outputFileName = optarg;
        std::cout << "[TrackingAnaData] INFO : Setting output file name: " << outputFileName << std::endl;
        break;
      }
      case 'v': {
        std::cout << "[TrackingAnaData] INFO : Enabling verbose debugging." << std::endl;
        verbose_debug = 1;
        break;
      }
      case 'h': {
        usage((const char*)argv[0]);
        return 0;
        break;
      }
      case 't': { //Get trigger filter
        triggerFilter = optarg;
        break;
      }
      case 'p': {
        physics_run = true;
      }
      case 'c': {
        std::stringstream ss;
        ss.str(optarg);
        try {
          ss >> max_chi2ndf;
        } catch (const char *) {
          std::cerr << "[TrackingAnaData] ERROR : Invalid argument for max chi2ndf (-c): " << optarg << std::endl;
          return 1;
          break;
        }
        break;
      }
    }
  }

  //print warning for any unused parameter
  if (optind < argc) {
    std::cout << "[TrackingAnaData] WARNING : Discarded parameters: ";
    while (optind < argc) {
      std::cout << argv[optind++] << " ";
    }
    std::cout << std::endl;
  }


 
 // --- Set style
  SetAtlasStyle(2);

  // --- init chains
  TChain vertexChain("VtxTree");
  TChain triggerChain("VtxTreeMeta/TrigConfTree");

  // --- set input files
  // Create list of files in eos with a small python script
  char* runNumber = argv[3];
  std::string run = runNumber;
  gSystem->Exec( ("python /afs/cern.ch/user/j/jiturbep/vertex_lumi/D3PDData/util/CreateList.py " + run).c_str() );

  //now add files from files recognized by ROOT TApplication
  std::string line;
  std::string name = "/afs/cern.ch/user/j/jiturbep/vertex_lumi/rootfiles_fullname_";
  std::string txt = ".txt";
  std::string step1 = name + run;
  std::string fullname = step1 + txt;
  std::ifstream rootfiles(fullname.c_str());
  if(!rootfiles){
    std::cout<<"Error opening output file"<< std::endl;
    system("pause");
    return -1;
  }
  while (std::getline(rootfiles, line))
    {
      inputRootFilesTxt.push_back(line);
    }  
  //TIter nextInputFile(theApp.InputFiles());
  //TObjString *InputFile;
  //while( (InputFile = (TObjString*)nextInputFile()) ) {
  // inputRootFilesTxt.push_back(InputFile->GetString().Data());
  // }

  //Print total number of input files and file list (if verbose) and add them to the chain
  std::cout << "[TrackingAnaData] INFO : Number of input files: " << inputRootFilesTxt.size() << std::endl;
  for (std::vector<std::string>::iterator irf = inputRootFilesTxt.begin(); irf != inputRootFilesTxt.end(); ++irf) {
    //if (verbose_debug) {
      std::cout << "[TrackingAnaData] INFO : Input file: " << *irf << std::endl;
      //}
    vertexChain.Add(irf->c_str());
    triggerChain.Add(irf->c_str());
  }

  if (inputRootFilesTxt.size() == 0) {
    std::cerr << "[TrackingAnaData] ERROR: No ROOT input files supplied." << std::endl;
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
    if (found == outputFileName.length()-1-5) { outputFileName.erase(found); }
    found = outputFileName.rfind(".txt");
    if (found == outputFileName.length()-1-4) { outputFileName.erase(found); }
  }

  // --- Create TSelector instance, set it up
  AnaVtxTree *vertexTreeSelector = new AnaVtxTree();
  vertexTreeSelector->SetOutputFileName(outputFileName);
  std::cout << "[TrackingAnaData] outputFileName = " << outputFileName << std::endl;
  vertexTreeSelector->SetDumpTxtFileName(asciiDumpFile);
  if (triggerFilter != "") { vertexTreeSelector->SetTriggerFilter(triggerFilter, &triggerChain); }
  if (max_chi2ndf > 0.) { vertexTreeSelector->SetMaxChi2Ndf(max_chi2ndf); }

  // --- initiate the Loop
  if (maxEvents < 0) {
    vertexChain.Process( dynamic_cast<TSelector*>(vertexTreeSelector) );
  } else {
    std::cout << "[TrackingAnaData] INFO : Running on " << maxEvents << " events" << std::endl;
    vertexChain.Process( dynamic_cast<TSelector*>(vertexTreeSelector), "", maxEvents);
  }

  // --- clean up
  delete vertexTreeSelector;
  return 0;
}
