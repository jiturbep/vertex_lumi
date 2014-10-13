#include "GlobalSettings/GlobalSettings.h"

#include <iostream>

using namespace std;

//#define DAVID_SETTINGS
//#define SIMONE_SETTINGS
//#define MANCHESTER_SETTINGS
#define LXPLUS_SETTINGS

GlobalSettings::GlobalSettings(std::string inputConfigFile) {
  if (inputConfigFile != "") {
    //load from file
    cerr << "Not currently implemented." << endl;
    return;
  }
}

GlobalSettings::~GlobalSettings() {}

// --- Overall path/file names settings
const std::string GlobalSettings::getPrefix( const TString &runNumber ) {
  if( runNumber == "182013" ) { return path_VdM_prefix; }
  if( runNumber == "201351" ) { return path_VdM_prefix; }
  if( runNumber == "188951" ) { return path_mu_prefix; }
  return path_data_prefix;
}

//- Path prefixes at Manchester
#ifdef MANCHESTER_SETTINGS
const std::string GlobalSettings::path_VdM_prefix = "VdMScan-";
const std::string GlobalSettings::path_mu_prefix = "muScan-";
const std::string GlobalSettings::path_data_prefix = "";
#elif defined LXPLUS_SETTINGS
const std::string GlobalSettings::path_VdM_prefix = "VdMScan-";
const std::string GlobalSettings::path_mu_prefix = "muScan-";
const std::string GlobalSettings::path_data_prefix = "";
#endif

//- Raw counts from D3PDData.
#ifdef DAVID_SETTINGS
const std::string GlobalSettings::path_inputRawCount = "/eliza18/atlas/dryu/Luminosity/VertexCounts/";
const std::string GlobalSettings::path_inputRawCount_v = "v9/";
const std::string GlobalSettings::fname_inputRawCount_Tree = "InDetTrackD3PD_results_tree.root";
const std::string GlobalSettings::fname_inputRawCount_Histo = "InDetTrackD3PD_results.root";
#elif defined SIMONE_SETTINGS
const std::string GlobalSettings::path_inputRawCount = "/eliza18/atlas/spagan/VtxLumi/run/VtxUnfold/data";
const std::string GlobalSettings::path_inputRawCount_v = "";
#elif defined MANCHESTER_SETTINGS
const std::string GlobalSettings::path_inputRawCount = "/afs/hep.man.ac.uk/d/atlas-lumi/VertexCounts/";
const std::string GlobalSettings::path_inputRawCount_v = "";
const std::string GlobalSettings::fname_inputRawCount_Tree = "InDetTrackD3PD_results_tree.root";
const std::string GlobalSettings::fname_inputRawCount_Histo = "InDetTrackD3PD_results.root";
#elif defined LXPLUS_SETTINGS
const std::string GlobalSettings::path_inputRawCount = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/";
const std::string GlobalSettings::path_inputRawCount_v = "";
const std::string GlobalSettings::fname_inputRawCount_Tree = "InDetTrackD3PD_results_tree.root";
const std::string GlobalSettings::fname_inputRawCount_Histo = "InDetTrackD3PD_results.root";
#endif

// - Pileup correction storage
#ifdef DAVID_SETTINGS
const std::string GlobalSettings::path_D3PDMCResults = "/eliza18/atlas/dryu/Luminosity/VtxTruthMatchResults/";
const std::string GlobalSettings::path_fakeCorrection = "/u/dryu/Luminosity/Data/FakeCorrection/";
const std::string GlobalSettings::path_maskingCorrection = "/eliza18/atlas/dryu/Luminosity/Data/MaskingCorrection/";
const std::string GlobalSettings::path_D3PDMCResults_v = "v2/";
#elif defined MANCHESTER_SETTINGS
const std::string GlobalSettings::path_D3PDMCResults = "/afs/hep.man.ac.uk/d/atlas-lumi/VtxTruthMatchResults/";
const std::string GlobalSettings::path_fakeCorrection = "/afs/hep.man.ac.uk/d/atlas-lumi/Data/FakeCorrection/";
const std::string GlobalSettings::path_maskingCorrection = "/afs/hep.man.ac.uk/d/atlas-lumi/Data/MaskingCorrection/";
const std::string GlobalSettings::path_D3PDMCResults_v = "";
#elif defined LXPLUS_SETTINGS
const std::string GlobalSettings::path_D3PDMCResults = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VtxTruthMatchResults/";
const std::string GlobalSettings::path_fakeCorrection = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/Data/FakeCorrection/";
const std::string GlobalSettings::path_maskingCorrection = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/Data/MaskingCorrection/";
const std::string GlobalSettings::path_D3PDMCResults_v = "";
#endif

//- Unfolding settings
#ifdef MANCHESTER_SETTINGS
const std::string GlobalSettings::path_unf_ResponseMatrix = "/eliza18/atlas/spagan/VtxLumi/run/VtxUnfold/TF";
#else
const std::string GlobalSettings::path_unf_ResponseMatrix = "/eliza18/atlas/spagan/VtxLumi/run/VtxUnfold/TF";
#endif
const std::string GlobalSettings::fname_unf_ResponseMatrix = "TransferFunctionMC_2011BS.root";
const std::string GlobalSettings::other_unf_ResponseMatrixHNameBase = "h_TF";

//- Output paths
#ifdef DAVID_SETTINGS
const std::string GlobalSettings::path_outputVdM = "/u/dryu/Luminosity/Data/VdMCalibration";
#elif defined SIMONE_SETTINGS
const std::string GlobalSettings::path_outputVdM = "/eliza18/atlas/spagan/VtxLumi/run/VtxUnfold/results/vdM";
#elif defined MANCHESTER_SETTINGS
const std::string GlobalSettings::path_outputVdM = "/afs/hep.man.ac.uk/d/atlas-lumi/Data/VdMCalibration/";
#elif defined LXPLUS_SETTINGS
const std::string GlobalSettings::path_outputVdM = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/Data/VdMCalibration/";
#endif
const std::string GlobalSettings::fname_outputVdM = "vdm_results.root";
const std::string GlobalSettings::fname_outputDebugVdM = "vdm_debug.root";

// VdM external information: timestamps
#ifdef MANCHESTER_SETTINGS
const std::string GlobalSettings::path_timestamps = "/afs/hep.man.ac.uk/d/atlas-lumi/LumiVtxExtras/timestamps/";
#elif defined LXPLUS_SETTINGS
const std::string GlobalSettings::path_timestamps = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/LumiVtxExtras/timestamps/";
#else
const std::string GlobalSettings::path_timestamps = "/u/dryu/Luminosity/Data/timestamps/";
#endif

// VdM external information: deadtime corrections
#ifdef MANCHESTER_SETTINGS
const std::string GlobalSettings::path_deadtime = "/afs/hep.man.ac.uk/d/atlas-lumi/LumiVtxExtras/deadtime/";
#elif defined LXPLUS_SETTINGS
const std::string GlobalSettings::path_deadtime = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/LumiVtxExtras/deadtime/";
#else
const std::string GlobalSettings::path_deadtime = "/u/dryu/Luminosity/Data/deadtime";
#endif

// VdM external information: Official lumi ntuples from COOL info
#ifdef MANCHESTER_SETTINGS
const std::string GlobalSettings::path_lumiNtuples = "/afs/hep.man.ac.uk/d/atlas-lumi/LumiVtxExtras/CoolScanNtuple/";
#elif defined LXPLUS_SETTINGS
const std::string GlobalSettings::path_lumiNtuples = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/LumiVtxExtras/CoolScanNtuple/";
#else
const std::string GlobalSettings::path_lumiNtuples = "/u/dryu/Luminosity/Data/CoolScanNtuple";
#endif

// External information: Official lumi ntuples for physics runs
#ifdef MANCHESTER_SETTINGS
const std::string GlobalSettings::path_PhysRun_lumiNtuples = "/afs/hep.man.ac.uk/d/atlas-lumi/LumiVtxExtras/LumiNtuples/";
#elif defined LXPLUS_SETTINGS
const std::string GlobalSettings::path_PhysRun_lumiNtuples = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/LumiVtxExtras/LumiNtuples/";
#else
const std::string GlobalSettings::path_PhysRun_lumiNtuples = "/u/dryu/Luminosity/Data/LumiNtuples";
#endif

// External information: Luminosity file comparison (Stefan's format)
#ifdef MANCHESTER_SETTINGS
const std::string GlobalSettings::path_lumiTxt = "/afs/hep.man.ac.uk/d/atlas-lumi/LumiVtxExtras/LumiTxt/";
#elif defined LXPLUS_SETTINGS
const std::string GlobalSettings::path_lumiTxt = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/LumiVtxExtras/LumiTxt/";
#else
const std::string GlobalSettings::path_lumiTxt = "/u/dryu/Luminosity/Data/LumiTxt/";
#endif
const std::string GlobalSettings::fname_lumiBcmH = "vdmmayH_bkg_mu.output";
const std::string GlobalSettings::fname_lumiBcmV = "vdmmayV_bkg_mu.output";
const std::string GlobalSettings::fname_muScanAll = "allMu_muscan188951.dat";
const std::string GlobalSettings::fname_r191373 = "191373_v5.dat";

// Skip-list of pLB for various special runs.
#ifdef MANCHESTER_SETTINGS
const std::string GlobalSettings::path_plb_skiplist = "/afs/hep.man.ac.uk/d/atlas-lumi/LumiVtxExtras/plb_skiplist/";
#elif defined LXPLUS_SETTINGS
const std::string GlobalSettings::path_plb_skiplist = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/LumiVtxExtras/plb_skiplist/";
#else
const std::string GlobalSettings::path_plb_skiplist = "/u/dryu/Luminosity/Data/plb_skiplist";
#endif

//- Results of LumiVtx
#ifdef MANCHESTER_SETTINGS
const std::string GlobalSettings::path_lumiVtx = "/afs/hep.man.ac.uk/d/atlas-lumi/LumiVtxExtras/LumiVtx/";
//const std::string GlobalSettings::path_lumiVtx = "/var/tmp/julia/";
#elif defined LXPLUS_SETTINGS
const std::string GlobalSettings::path_lumiVtx = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/LumiVtxExtras/LumiVtx/";
//const std::string GlobalSettings::path_lumiVtx = "/var/tmp/julia/";
#else
const std::string GlobalSettings::path_lumiVtx = "/u/dryu/Luminosity/Data/LumiVtx/";
#endif
