#ifndef VtxTree_h
#define VtxTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

class VtxTree : public TSelector {
  protected :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  
    // Declaration of leaf types
    std::vector<unsigned int> *trig_L1_TAV;
    std::vector<short>   *trig_L2_passedPhysics;
    std::vector<short>   *trig_EF_passedPhysics;
    UInt_t          trig_DB_SMK;
    UInt_t          trig_DB_L1PSK;
    UInt_t          trig_DB_HLTPSK;
    UInt_t          ei_RunNumber;
    UInt_t          ei_EventNumber;
    UInt_t          ei_timestamp;
    UInt_t          ei_timestamp_ns;
    UInt_t          ei_lbn;
    UInt_t          ei_bcid;
    UInt_t          ei_detmask0;
    UInt_t          ei_detmask1;
    Float_t         ei_actualIntPerXing;
    Float_t         ei_averageIntPerXing;
    Int_t           vxnbc_n;
    std::vector<float>   *vxnbc_x;
    std::vector<float>   *vxnbc_y;
    std::vector<float>   *vxnbc_z;
    std::vector<float>   *vxnbc_cov_x;
    std::vector<float>   *vxnbc_cov_y;
    std::vector<float>   *vxnbc_cov_z;
    std::vector<float>   *vxnbc_cov_xy;
    std::vector<float>   *vxnbc_cov_xz;
    std::vector<float>   *vxnbc_cov_yz;
    std::vector<int>     *vxnbc_type;
    std::vector<float>   *vxnbc_chi2;
    std::vector<int>     *vxnbc_ndof;
    std::vector<float>   *vxnbc_px;
    std::vector<float>   *vxnbc_py;
    std::vector<float>   *vxnbc_pz;
    std::vector<float>   *vxnbc_E;
    std::vector<float>   *vxnbc_m;
    std::vector<int>     *vxnbc_nTracks;
    std::vector<float>   *vxnbc_sumPt;
    std::vector<std::vector<float> > *vxnbc_trk_weight;
    std::vector<int>     *vxnbc_trk_n;
    std::vector<std::vector<int> > *vxnbc_trk_index;
    Int_t           trk_n;
    std::vector<float>   *trk_d0;
    std::vector<float>   *trk_z0;
    std::vector<float>   *trk_phi;
    std::vector<float>   *trk_theta;
    std::vector<float>   *trk_qoverp;
    std::vector<float>   *trk_pt;
    std::vector<float>   *trk_eta;
    std::vector<float>   *trk_d0_wrtPV;
    std::vector<float>   *trk_z0_wrtPV;
    std::vector<float>   *trk_phi_wrtPV;
    Int_t           mcvtx_n;
    std::vector<float>   *mcvtx_x;
    std::vector<float>   *mcvtx_y;
    std::vector<float>   *mcvtx_z;
    std::vector<int>     *mcvtx_barcode;
    std::vector<int>     *mcvtx_mcevt_index;
    std::vector<int>     *mcevt_nparticle;
    std::vector<short>   *mcevt_pileUpType;

    // List of branches
    TBranch        *b_trig_L1_TAV;   //!
    TBranch        *b_trig_L2_passedPhysics;   //!
    TBranch        *b_trig_EF_passedPhysics;   //!
    TBranch        *b_trig_DB_SMK;   //!
    TBranch        *b_trig_DB_L1PSK;   //!
    TBranch        *b_trig_DB_HLTPSK;   //!
    TBranch        *b_ei_RunNumber;   //!
    TBranch        *b_ei_EventNumber;   //!
    TBranch        *b_ei_timestamp;   //!
    TBranch        *b_ei_timestamp_ns;   //!
    TBranch        *b_ei_lbn;   //!
    TBranch        *b_ei_bcid;   //!
    TBranch        *b_ei_detmask0;   //!
    TBranch        *b_ei_detmask1;   //!
    TBranch        *b_ei_actualIntPerXing;   //!
    TBranch        *b_ei_averageIntPerXing;   //!
    TBranch        *b_vxnbc_n;   //!
    TBranch        *b_vxnbc_x;   //!
    TBranch        *b_vxnbc_y;   //!
    TBranch        *b_vxnbc_z;   //!
    TBranch        *b_vxnbc_cov_x;   //!
    TBranch        *b_vxnbc_cov_y;   //!
    TBranch        *b_vxnbc_cov_z;   //!
    TBranch        *b_vxnbc_cov_xy;   //!
    TBranch        *b_vxnbc_cov_xz;   //!
    TBranch        *b_vxnbc_cov_yz;   //!
    TBranch        *b_vxnbc_type;   //!
    TBranch        *b_vxnbc_chi2;   //!
    TBranch        *b_vxnbc_ndof;   //!
    TBranch        *b_vxnbc_px;   //!
    TBranch        *b_vxnbc_py;   //!
    TBranch        *b_vxnbc_pz;   //!
    TBranch        *b_vxnbc_E;   //!
    TBranch        *b_vxnbc_m;   //!
    TBranch        *b_vxnbc_nTracks;   //!
    TBranch        *b_vxnbc_sumPt;   //!
    TBranch        *b_vxnbc_trk_weight;   //!
    TBranch        *b_vxnbc_trk_n;   //!
    TBranch        *b_vxnbc_trk_index;   //!
    TBranch        *b_trk_n;   //!
    TBranch        *b_trk_d0;   //!
    TBranch        *b_trk_z0;   //!
    TBranch        *b_trk_phi;   //!
    TBranch        *b_trk_theta;   //!
    TBranch        *b_trk_qoverp;   //!
    TBranch        *b_trk_pt;   //!
    TBranch        *b_trk_eta;   //!
    TBranch        *b_trk_d0_wrtPV;   //!
    TBranch        *b_trk_z0_wrtPV;   //!
    TBranch        *b_trk_phi_wrtPV;   //!
  
    VtxTree(TTree * /*tree*/ =0) : fChain(0) { }
    virtual ~VtxTree() { }
    virtual Int_t   Version() const { return 2; }
    virtual void    Begin(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0);
    virtual void    SetOption(const char *option);
    virtual void    SetObject(TObject *obj);
    virtual TList  *GetOutputList() const;
    virtual void    SlaveTerminate();
    virtual void    Terminate();

    ClassDef(VtxTree,0);
};

#endif
