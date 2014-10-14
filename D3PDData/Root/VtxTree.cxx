// The class definition in VtxTree.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("VtxTree.C")
// Root > T->Process("VtxTree.C","some options")
// Root > T->Process("VtxTree.C+")
//

#include "D3PDData/VtxTree.h"
#include <TH2.h>
#include <TStyle.h>


void VtxTree::Begin(TTree * /*tree*/) {
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

}

void VtxTree::SlaveBegin(TTree * /*tree*/) {
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

}

Bool_t VtxTree::Process(Long64_t entry) {
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either VtxTree::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.


  return kTRUE;
}

void VtxTree::SlaveTerminate() {
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}

void VtxTree::Terminate() {
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
}

Bool_t VtxTree::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void VtxTree::Init(TTree *tree) {
  // Set branch addresses and branch pointers
  if (!tree) {
    return;
  }
  fChain = tree;
  fChain->SetMakeClass(1);
  
  // Set branch addresses
  fChain->SetBranchAddress("ei_RunNumber", &ei_RunNumber, &b_ei_RunNumber);
  fChain->SetBranchAddress("ei_EventNumber", &ei_EventNumber, &b_ei_EventNumber);
  fChain->SetBranchAddress("ei_timestamp", &ei_timestamp, &b_ei_timestamp);
  fChain->SetBranchAddress("ei_timestamp_ns", &ei_timestamp_ns, &b_ei_timestamp_ns);
  fChain->SetBranchAddress("ei_lbn", &ei_lbn, &b_ei_lbn);
  fChain->SetBranchAddress("ei_bcid", &ei_bcid, &b_ei_bcid);
  fChain->SetBranchAddress("ei_actualIntPerXing", &ei_actualIntPerXing, &b_ei_actualIntPerXing);
  fChain->SetBranchAddress("ei_averageIntPerXing", &ei_averageIntPerXing, &b_ei_averageIntPerXing);
  // --- These branches aren't in the files
  //fChain->SetBranchAddress("trig_L1_TAV", &trig_L1_TAV, &b_trig_L1_TAV);
  //fChain->SetBranchAddress("trig_L2_passedPhysics", &trig_L2_passedPhysics, &b_trig_L2_passedPhysics);
  //fChain->SetBranchAddress("trig_EF_passedPhysics", &trig_EF_passedPhysics, &b_trig_EF_passedPhysics);
  //fChain->SetBranchAddress("trig_DB_SMK", &trig_DB_SMK, &b_trig_DB_SMK);
  //fChain->SetBranchAddress("trig_DB_L1PSK", &trig_DB_L1PSK, &b_trig_DB_L1PSK);
  //fChain->SetBranchAddress("trig_DB_HLTPSK", &trig_DB_HLTPSK, &b_trig_DB_HLTPSK);
  fChain->SetBranchAddress("mcvtx_n", &mcvtx_n, &b_mcvtx_n);
  fChain->SetBranchAddress("mcevt_nparticle", &mcevt_nparticle, &b_mcevt_nparticle);
  fChain->SetBranchAddress("mcvtx_mcevt_index", &mcvtx_mcevt_index, &b_mcvtx_mcevt_index);
  fChain->SetBranchAddress("mcevt_pileUpType", &mcevt_pileUpType, &b_mcevt_pileUpType);
  fChain->SetBranchAddress("vxnbc_n", &vxnbc_n, &b_vxnbc_n);
  fChain->SetBranchAddress("vxnbc_x", &vxnbc_x, &b_vxnbc_x);
  fChain->SetBranchAddress("vxnbc_y", &vxnbc_y, &b_vxnbc_y);
  fChain->SetBranchAddress("vxnbc_z", &vxnbc_z, &b_vxnbc_z);
  fChain->SetBranchAddress("vxnbc_cov_x", &vxnbc_cov_x, &b_vxnbc_cov_x);
  fChain->SetBranchAddress("vxnbc_cov_y", &vxnbc_cov_y, &b_vxnbc_cov_y);
  fChain->SetBranchAddress("vxnbc_cov_z", &vxnbc_cov_z, &b_vxnbc_cov_z);
  fChain->SetBranchAddress("vxnbc_cov_xy", &vxnbc_cov_xy, &b_vxnbc_cov_xy);
  fChain->SetBranchAddress("vxnbc_cov_xz", &vxnbc_cov_xz, &b_vxnbc_cov_xz);
  fChain->SetBranchAddress("vxnbc_cov_yz", &vxnbc_cov_yz, &b_vxnbc_cov_yz);
  fChain->SetBranchAddress("vxnbc_type", &vxnbc_type, &b_vxnbc_type);
  fChain->SetBranchAddress("vxnbc_chi2", &vxnbc_chi2, &b_vxnbc_chi2);
  fChain->SetBranchAddress("vxnbc_ndof", &vxnbc_ndof, &b_vxnbc_ndof);
  fChain->SetBranchAddress("vxnbc_px", &vxnbc_px, &b_vxnbc_px);
  fChain->SetBranchAddress("vxnbc_py", &vxnbc_py, &b_vxnbc_py);
  fChain->SetBranchAddress("vxnbc_pz", &vxnbc_pz, &b_vxnbc_pz);
  fChain->SetBranchAddress("vxnbc_E", &vxnbc_E, &b_vxnbc_E);
  fChain->SetBranchAddress("vxnbc_m", &vxnbc_m, &b_vxnbc_m);
  fChain->SetBranchAddress("vxnbc_nTracks", &vxnbc_nTracks, &b_vxnbc_nTracks);
  fChain->SetBranchAddress("vxnbc_sumPt", &vxnbc_sumPt, &b_vxnbc_sumPt);
  fChain->SetBranchAddress("vxnbc_trk_weight", &vxnbc_trk_weight, &b_vxnbc_trk_weight);
  fChain->SetBranchAddress("vxnbc_trk_n", &vxnbc_trk_n, &b_vxnbc_trk_n);
  fChain->SetBranchAddress("vxnbc_trk_index", &vxnbc_trk_index, &b_vxnbc_trk_index);
  fChain->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
  fChain->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
 
  
  // Initialise TTree
  fChain->SetBranchStatus("*", 0);
  fChain->SetBranchStatus("ei_RunNumber",1);
  fChain->SetBranchStatus("ei_EventNumber",1);
  fChain->SetBranchStatus("ei_timestamp",1);
  fChain->SetBranchStatus("ei_timestamp_ns",1);
  fChain->SetBranchStatus("ei_lbn",1);
  fChain->SetBranchStatus("ei_bcid",1);
  fChain->SetBranchStatus("ei_actualIntPerXing",1);
  fChain->SetBranchStatus("ei_averageIntPerXing",1);
  // --- These branches aren't in the files
  //fChain->SetBranchStatus("trig_L1_TAV",1);
  //fChain->SetBranchStatus("trig_L2_passedPhysics",1);
  //fChain->SetBranchStatus("trig_EF_passedPhysics",1);
  //fChain->SetBranchStatus("trig_DB_SMK",1);
  //fChain->SetBranchStatus("trig_DB_L1PSK",1);
  //fChain->SetBranchStatus("trig_DB_HLTPSK",1);
  fChain->SetBranchStatus("mcvtx_n",1);
  fChain->SetBranchStatus("mcevt_nparticle",1);
  fChain->SetBranchStatus("mcvtx_mcevt_index",1);
  fChain->SetBranchStatus("mcevt_pileUpType",1);
  fChain->SetBranchStatus("vxnbc_n",1);
  fChain->SetBranchStatus("vxnbc_x",1);
  fChain->SetBranchStatus("vxnbc_y",1);
  fChain->SetBranchStatus("vxnbc_z",1);
  fChain->SetBranchStatus("vxnbc_cov_x",1);
  fChain->SetBranchStatus("vxnbc_cov_y",1);
  fChain->SetBranchStatus("vxnbc_cov_z",1);
  fChain->SetBranchStatus("vxnbc_cov_xy",1);
  fChain->SetBranchStatus("vxnbc_cov_xz",1);
  fChain->SetBranchStatus("vxnbc_cov_yz",1);
  fChain->SetBranchStatus("vxnbc_chi2",1);
  fChain->SetBranchStatus("vxnbc_ndof",1);
  fChain->SetBranchStatus("vxnbc_px",1);
  fChain->SetBranchStatus("vxnbc_py",1);
  fChain->SetBranchStatus("vxnbc_pz",1);
  fChain->SetBranchStatus("vxnbc_E",1);
  fChain->SetBranchStatus("vxnbc_m",1);
  fChain->SetBranchStatus("vxnbc_nTracks",1);
  fChain->SetBranchStatus("vxnbc_sumPt",1);
  fChain->SetBranchStatus("vxnbc_type",1);
  fChain->SetBranchStatus("vxnbc_trk_index",1);
  fChain->SetBranchStatus("trk_pt",1);
  fChain->SetBranchStatus("trk_eta",1);
  fChain->SetBranchStatus("trk_chi2",1);
  fChain->SetBranchStatus("trk_ndof",1);
  fChain->SetBranchStatus("trk_nBLHits",1);
  fChain->SetBranchStatus("trk_nPixHits",1);
  fChain->SetBranchStatus("trk_nSCTHits",1);
  fChain->SetBranchStatus("trk_nTRTHits",1);
  fChain->SetBranchStatus("trk_nTRTHighTHits",1);
  fChain->SetBranchStatus("trk_nPixHoles",1);
  fChain->SetBranchStatus("trk_nSCTHoles",1);
  fChain->SetBranchStatus("trk_nTRTHoles",1);
  fChain->SetBranchStatus("trk_nHits",1);
}

Int_t VtxTree::GetEntry(Long64_t entry, Int_t getall ) {
  return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
}

void VtxTree::SetOption(const char *option) {
  fOption = option;
}

void VtxTree::SetObject(TObject *obj) {
  fObject = obj;
}

TList* VtxTree::GetOutputList() const {
  return fOutput;
}

