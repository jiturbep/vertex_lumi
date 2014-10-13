//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 19 15:48:31 2012 by ROOT version 5.28/00e
// from TTree InDetTrackTree/InDetTrackTree
// found on file: vertex_d3pd_pu20.root
//
// Then I chopped it up to use it more generally
//////////////////////////////////////////////////////////

#ifndef VertexTree_h
#define VertexTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <string>

using namespace std;

class VertexTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   string           fOutName;

   // Declaration of leaf types
   UInt_t          ei_RunNumber;
   UInt_t          ei_EventNumber;

   Int_t           trk_n;
   vector<float>   *trk_pt;
   vector<float>   *trk_mc_probability;
   vector<int>     *trk_mc_index;

   Int_t           vx_n;
   vector<int>     *vx_trk_n;
   vector<vector<float> > *vx_trk_weight;
   vector<vector<int> > *vx_trk_index;

   Int_t           mcevt_n;
   vector<int>     *mcevt_nparticle;
   vector<short>   *mcevt_pileUpType;

   Int_t           mcvtx_n;
   vector<int>     *mcvtx_mcevt_index;

   Int_t           mcpart_n;
   vector<float>   *mcpart_pt;
   vector<float>   *mcpart_eta;
   vector<int>     *mcpart_type;
   vector<int>     *mcpart_status;
   vector<int>     *mcpart_barcode;
   vector<int>     *mcpart_mcevt_index;


   // List of branches
   TBranch        *b_ei_RunNumber;   //!
   TBranch        *b_ei_EventNumber;   //!

   TBranch        *b_trk_n;   //!
   TBranch        *b_trk_pt;   //!
   TBranch        *b_trk_mc_probability;   //!
   TBranch        *b_trk_mc_index;   //!

   TBranch        *b_vx_n;   //!
   TBranch        *b_vx_trk_n;   //!
   TBranch        *b_vx_trk_weight;   //!
   TBranch        *b_vx_trk_index;   //!

   TBranch        *b_mcevt_n;   //!
   TBranch        *b_mcevt_nparticle;   //!
   TBranch        *b_mcevt_pileUpType;   //!

   TBranch        *b_mcvtx_n;   //!
   TBranch        *b_mcvtx_mcevt_index;   //!

   TBranch        *b_mcpart_n;   //!
   TBranch        *b_mcpart_pt;   //!
   TBranch        *b_mcpart_eta;   //!
   TBranch        *b_mcpart_type;   //!
   TBranch        *b_mcpart_status;   //!
   TBranch        *b_mcpart_barcode;   //!
   TBranch        *b_mcpart_mcevt_index;   //!

  VertexTree(TTree *tree, string outname = "TestVtxTmD3PD.root"); //require input tree
   virtual ~VertexTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef VertexTree_cxx
VertexTree::VertexTree(TTree *tree, string outname)
{
  fOutName = outname;

  Init(tree);
}

VertexTree::~VertexTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t VertexTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t VertexTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void VertexTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trk_pt = 0;
   trk_mc_probability = 0;
   trk_mc_index = 0;

   vx_trk_n = 0;
   vx_trk_weight = 0;
   vx_trk_index = 0;

   mcevt_nparticle = 0;
   mcevt_pileUpType = 0;

   mcvtx_mcevt_index = 0;

   mcpart_pt = 0;
   mcpart_eta = 0;
   mcpart_type = 0;
   mcpart_status = 0;
   mcpart_barcode = 0;
   mcpart_mcevt_index = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);  //never understood what this does...

   fChain->SetBranchAddress("ei_RunNumber", &ei_RunNumber, &b_ei_RunNumber);
   fChain->SetBranchAddress("ei_EventNumber", &ei_EventNumber, &b_ei_EventNumber);

   fChain->SetBranchAddress("trk_n", &trk_n, &b_trk_n);
   fChain->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
   fChain->SetBranchAddress("trk_mc_probability", &trk_mc_probability, &b_trk_mc_probability);
   fChain->SetBranchAddress("trk_mc_index", &trk_mc_index, &b_trk_mc_index);

   fChain->SetBranchAddress("vx_n", &vx_n, &b_vx_n);
   fChain->SetBranchAddress("vx_trk_n", &vx_trk_n, &b_vx_trk_n);
   fChain->SetBranchAddress("vx_trk_weight", &vx_trk_weight, &b_vx_trk_weight);
   fChain->SetBranchAddress("vx_trk_index", &vx_trk_index, &b_vx_trk_index);

   fChain->SetBranchAddress("mcevt_n", &mcevt_n, &b_mcevt_n);
   fChain->SetBranchAddress("mcevt_nparticle", &mcevt_nparticle, &b_mcevt_nparticle);
   fChain->SetBranchAddress("mcevt_pileUpType", &mcevt_pileUpType, &b_mcevt_pileUpType);

   fChain->SetBranchAddress("mcvtx_n", &mcvtx_n, &b_mcvtx_n);
   fChain->SetBranchAddress("mcvtx_mcevt_index", &mcvtx_mcevt_index, &b_mcvtx_mcevt_index);

   fChain->SetBranchAddress("mcpart_n", &mcpart_n, &b_mcpart_n);
   fChain->SetBranchAddress("mcpart_pt", &mcpart_pt, &b_mcpart_pt);
   fChain->SetBranchAddress("mcpart_eta", &mcpart_eta, &b_mcpart_eta);
   fChain->SetBranchAddress("mcpart_type", &mcpart_type, &b_mcpart_type);
   fChain->SetBranchAddress("mcpart_status", &mcpart_status, &b_mcpart_status);
   fChain->SetBranchAddress("mcpart_barcode", &mcpart_barcode, &b_mcpart_barcode);
   fChain->SetBranchAddress("mcpart_mcevt_index", &mcpart_mcevt_index, &b_mcpart_mcevt_index);
   Notify();
}

Bool_t VertexTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void VertexTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t VertexTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef VertexTree_cxx
