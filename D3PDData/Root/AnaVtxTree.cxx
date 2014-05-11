// Local
#include "D3PDData/ATLASImport.h"
#include "D3PDData/AnaVtxTree.h"
#include "D3PDData/HistogramHelper.h"
#include "GlobalSettings/GlobalSettings.h"

// STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <bitset>
#include <string>
#include <iomanip>

// ROOT
#include "TDataType.h"
#include "TMath.h"
#include "TVector2.h"
#include "TAxis.h"

////////////////////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////////////////////
AnaVtxTree::AnaVtxTree() : VtxTree()
  , m_GRLReader("../GoodRunsLists/data12_8TeV.periodAllYear_DetStatus-v54-pro13-04_DQDefects-00-00-33_PHYS_CombinedPerf_Tracking_Tracking.xml") // GRL reader tool
  , use_plbs(false) // Specify if pLBs are needed
  , split_maxdzsig(0.) // Split vertex re-merging
  , max_chi2ndf(-1.) // Apply cut on chi2/ndf of reconstructed vertices
  , current_pLB(-1) //store current pLB
  , current_pLB_run(-1) //store run for which we have loaded pLB info
  , pLBmin(100000) // Unreal value, must be set by LoadAndAddPseudoLBs at some point!
  , pLBmax(0)
  , dumpTxtFileName("") //by default no dumping on text file
  , physics_run(false)
  , qualityVertexVersion(3)
  , qualityVertexModifier(0)
  , qualityVertexParameter(5)
{};

////////////////////////////////////////////////////////////////////////////////////////
// TSelector begin
////////////////////////////////////////////////////////////////////////////////////////
void AnaVtxTree::SlaveBegin(TTree *tree) {
  //std::cout.precision(14);

  VtxTree::SlaveBegin(tree);
  flag_first = true;

  // --- Init output text file
  outputTxt.open((outputFileName + TString(".txt")).Data());

  // --- Create histograms? Can only do BCID-blind histos here...

  //Init counters
  m_TotEvents = 0;
  m_TotRawEvents = 0;
  m_triggerCounts["L1_BGRP7"] = 0;
  m_triggerCounts["L1_MBTS_2_BGRP7"] = 0;
  m_triggerCounts["L1_MBTS_2_BGRP7 && L1_BGRP7"] = 0;
  m_triggerCounts["L1_MBTS_2_BGRP7 && !L1_BGRP7"] = 0;
  m_triggerCounts["!L1_MBTS_2_BGRP7 && L1_BGRP7"] = 0;
  m_triggerCounts["!L1_MBTS_2_BGRP7 && !L1_BGRP7"] = 0;

  m_InitBcidHistograms = true;

  if (dumpTxtFileName != "") {
    std::cout << "[AnaVtxTree] INFO: Enabling dump on ASCII file: " << dumpTxtFileName << std::endl;
    //ok .. you want it.. dump on a text file
    dumpTxtFile.open(dumpTxtFileName.Data());
  }

  nTrkCuts.push_back(2);
  nTrkCuts.push_back(3); 
  nTrkCuts.push_back(4);
  nTrkCuts.push_back(5);
  nTrkCuts.push_back(6);
  nTrkCuts.push_back(7);
  nTrkCuts.push_back(8);
  nTrkCuts.push_back(10);

}

////////////////////////////////////////////////////////////////////////////////////////
// TSelector event loop
////////////////////////////////////////////////////////////////////////////////////////
Bool_t AnaVtxTree::Process(Long64_t entry) {

  // --- Setup trigger D3PD decoding tool
  if (trigMetaDataTree!=0 && triggerName != "" && !triggerTool) {
    triggerTool = new D3PD::TrigDecisionToolD3PD(fChain, trigMetaDataTree);
    std::cout << "[AnaVtxTree] INFO: Trigger tool initialized" << std::endl;
  }

  // --- Get event variables from TTree
  fChain->GetTree()->GetEntry(entry);

  // --- Print progress
  if (m_TotRawEvents % 1000 == 0) {
    std::cout << "[AnaVtxTree] INFO: Processing event " << m_TotRawEvents
              << "(Entry: " << entry << ", Tree #: " << fChain->GetTreeNumber() << ")" << std::endl;
  }
  ++m_TotRawEvents; // before GRL requirement

  // --- Don't check GRL for a VdM scan
  //if( !m_GRLReader.passedGRL( ei_RunNumber, ei_lbn ) ) {
  //  std::cout << "Failed GRL for runNumber("<<ei_RunNumber<<"), lbn("<<ei_lbn<<") " << std::endl;
  //}

  // --- Retrieve trigger tool if specified
  if (triggerTool) { triggerTool->GetEntry(entry); }

  // --- If this is the first event then set appropriate variables
  if( m_TotRawEvents == 1) { SetupLBInfo(ei_RunNumber); }

  // --- Quick rejection of events outside the scan based on LB, trigger and BCID
  if( (ei_RunNumber == 188949) && ( (ei_lbn <= 92) || (ei_lbn >= 160) ) ) { return false;}    
  if( (ei_RunNumber == 188951) && ( (ei_lbn <= 130) || (ei_lbn >= 185) ) ) {return false;}
  if (ei_RunNumber == 191373 || ei_RunNumber == 200805 || ei_RunNumber == 201351 || 
      ei_RunNumber == 215021 || ei_RunNumber == 214984 || ei_RunNumber == 207216 || 
      ei_RunNumber == 207219 || ei_RunNumber == 216399 || ei_RunNumber == 216416 || 
      ei_RunNumber == 216432 || ei_RunNumber == 207044 ) { //Added 215021, 214984 -Nov2012, 207216, 207219 -July2012 //18Nov: Deleted 206962, 206955, 206971 to treat it BCID blindly
    if (triggerTool and triggerName != "" && !(triggerTool->IsPassed(triggerName.Data())) ) { return false; } // Trigger
    if (find(bcidListCollisions.begin(), bcidListCollisions.end(), ei_bcid) == bcidListCollisions.end()) { return false; } // BCID
  }

  //Some pre-selection, based on timestamps.
  Double_t curTime = ei_timestamp + ei_timestamp_ns * TMath::Power(10.,-9);
  if ((ei_RunNumber == 188951) || (ei_RunNumber == 188949)) {
    if (timestampInScan(curTime) == 0) { return false; }
  }

  // 200805 hack trigger
  if ((ei_RunNumber == 200805) && (*trig_L1_TAV)[1] != 268435456) {
    return false;
  }

  // --- Initialise histograms
  if (m_InitBcidHistograms) {
    InitBCIDHists();
    m_InitBcidHistograms = false;
  } 
  //Find the pLB. Note that the variable current_pLB will specify the pLB if inside the scan, elsewise the normal lumiblock!
  #if defined INPUT_DATA
  //std::cout << "[AnaVtxTree] INFO: INPUT_DATA" << std::endl;
  //If input is data, decide whether to use pLBs or normal LBs.
  if (use_plbs) {
  //std::cout << "[AnaVtxTree] INFO: use_plbs" << std::endl;

    bool goodPLB=false;
    for (unsigned int idx=0; idx < pseudoLB_Timestamps.size(); idx++) {
      if ((curTime >= pseudoLB_Timestamps[idx].second.first) && (curTime <= pseudoLB_Timestamps[idx].second.second)) {
        current_pLB = pseudoLB_Timestamps[idx].first;
        if (flag_first) {
          std::cout << "[AnaVtxTree] INFO: First event: " << m_TotRawEvents << std::endl;
          flag_first = false;
        }
        goodPLB = true;
        //if (m_TotRawEvents % 1000 == 0) std::cout << "current_pLB = " << current_pLB << std::endl;
        break;
      }
    }
    // discard event if pLB info is available and event is not inside that list
    if (!goodPLB) {
      return false;
    }

  } else {
    current_pLB = ei_lbn;
   // std::cout << "ei_lbn = " << ei_lbn << std::endl;
  }

  #elif defined INPUT_MC
  
  //If input is MC, use the number of generated interactions (in-time pileup only) as the "lumiblock."
  // REVISION: How about we just use ei_actualIntPerXing? Maybe this will reduce error later on. 
  Int_t current_NGenInt = 0;

  for (int ivx=0; ivx < mcvtx_n; ++ivx) {

  if (ivx == 0 && (*mcevt_nparticle)[0] == 1 && mcevt_pileUpType->at(mcvtx_mcevt_index->at(ivx)) == 0) continue; // In new MC, hard scatter is empty.

  //if (TmVtxNBc.isGenInteraction(ivx, ei_RunNumber, ei_EventNumber)) {
  Int_t current_type = mcevt_pileUpType->at(mcvtx_mcevt_index->at(ivx));
  if ((current_type == 0) || (current_type == 1)) {
  current_NGenInt++;
    }
  }
  //current_pLB = current_NGenInt;
  current_pLB = ei_actualIntPerXing;
  #endif 
  // end INPUT_MC

  Int_t current_bcid;

  #ifdef INPUT_DATA
  
  if (physics_run) {
  current_bcid = 0;
  } else {
  current_bcid = ei_bcid;
  }
  #elif defined INPUT_MC
  current_bcid = 0;
  #endif

  /*if (m_TotRawEvents % 1000 == 0) {
    std::cout << "[AnaVtxTree] INFO: current_pLB = " << current_pLB << std::endl;
  }*/

  ++m_TotEvents;

  //Get the index of the primary vertex
  Int_t TagVtxIndex = isGoodVertex();
  
  //Set the event flags to false
  for (std::vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    m_isTightBeforeSplitCorrection[*nTrkCut] = false;
    m_isTightAfterSplitCorrection[*nTrkCut] = false;
  }
  
  //Fill primary vertex histogram(s)
  if (TagVtxIndex >= 0) {
    h_privtx_z_pLB[current_bcid]->Fill((*vxnbc_z)[TagVtxIndex], current_pLB); 
  }

  h_events_pLB[current_bcid]->Fill(current_pLB); 
  
  //Increment the total number of "triggers" map
  NTrig[current_bcid][current_pLB]++;
  for (std::vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
        NVtxBeforeSplitCorrection_pertrack[*nTrkCut] = 0; 
        NVtxBeforeSplitCorrection_pertrack_perbcid[*nTrkCut][current_bcid] = 0; 
  }
  
  //Vertex counts before applying the split correction
  for (Int_t nv = 0; nv < vxnbc_n - 1; nv++) {
    if (max_chi2ndf > 0. && (*vxnbc_chi2)[nv] / (*vxnbc_ndof)[nv] > max_chi2ndf) {
      continue;
    }

    Int_t current_ntrk = (*vxnbc_nTracks)[nv];
    for (std::vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
      if (current_ntrk >= *nTrkCut) {
        NVtxBeforeSplitCorrection[*nTrkCut][current_bcid][current_pLB]++;
        NVtxBeforeSplitCorrection_pertrack_perbcid[*nTrkCut][current_bcid]++;
        NVtxBeforeSplitCorrection_pertrack[*nTrkCut]++;
        m_isTightBeforeSplitCorrection[*nTrkCut] = true;
      }
    }
  }
  
  for (std::vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    h_nvtx_pertrack[*nTrkCut]->Fill(NVtxBeforeSplitCorrection_pertrack[*nTrkCut]);
    h_nvtx_pertrack_perbcid[*nTrkCut][current_bcid]->Fill(NVtxBeforeSplitCorrection_pertrack_perbcid[*nTrkCut][current_bcid]);
    //std::cout << "Line 218 NVtxBeforeSplitCorrection_pertrack["<<*nTrkCut<<"] = "<<NVtxBeforeSplitCorrection_pertrack[*nTrkCut]<<std::endl;
    }
  //#ifdef VERBOSE
  for (std::vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    h_nvtx_pLB[*nTrkCut][current_bcid]->Fill(NVtxBeforeSplitCorrection[*nTrkCut][current_bcid][current_pLB], current_pLB);
  }
  //#endif

  //Evil hungry vertices cannibalize nearby smaller brethren.
  std::vector<Int_t> vtx_skip_list;
  for (Int_t nv = 0; nv < vxnbc_n - 1; nv++) {
    if (find(vtx_skip_list.begin(), vtx_skip_list.end(), nv) != vtx_skip_list.end()) {
      continue;
    }

    Int_t extra_tracks = 0;
    for (Int_t nv2 = nv; nv2 < vxnbc_n - 1; nv2++) {
      if (nv == nv2) {
        continue;
      }

      if (find(vtx_skip_list.begin(), vtx_skip_list.end(), nv2) != vtx_skip_list.end()) {
        continue;
      }

      float err1 = (*vxnbc_cov_z)[nv];
      float err2 = (*vxnbc_cov_z)[nv2];
      float err = TMath::Sqrt(err1 + err2);
      float dz = (*vxnbc_z)[nv] - (*vxnbc_z)[nv2];

      Int_t ntrk11 = (*vxnbc_nTracks)[nv];
      Int_t ntrk12 = (*vxnbc_nTracks)[nv2];
      for (std::vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
        if ((ntrk11 >= *nTrkCut) && (ntrk12 >= *nTrkCut)) {
          h_vtx_dz_dzsig[*nTrkCut]->Fill(dz, dz / err);
        }
      }

      if (TMath::Abs(dz) < split_maxdzsig * err) {
        vtx_skip_list.push_back(nv2);
        extra_tracks += (*vxnbc_nTracks)[nv2];
      } else {
        #ifdef VERBOSE
        //Fill dz histograms.
        Int_t ntrk1 = (*vxnbc_nTracks)[nv];
        Int_t ntrk2 = (*vxnbc_nTracks)[nv2];
        for (std::vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
          if ((ntrk1 >= *nTrkCut) && (ntrk2 >= *nTrkCut)) {
            h_vtx_Dz_pLB[*nTrkCut][current_bcid]->Fill(dz, current_pLB);
          } else if ((ntrk1 >= *nTrkCut) && (ntrk2 < *nTrkCut)) {
            h_vtx_Dz_TightLoose_pLB[*nTrkCut][current_bcid]->Fill(dz, current_pLB);
          } else if ((ntrk1 < *nTrkCut) && (ntrk2 >= *nTrkCut)) {
            h_vtx_Dz_LooseTight_pLB[*nTrkCut][current_bcid]->Fill(dz, current_pLB);
          } else if ((ntrk1 < *nTrkCut) && (ntrk2 < *nTrkCut)) {
            h_vtx_Dz_LooseLoose_pLB[*nTrkCut][current_bcid]->Fill(dz, current_pLB);
          }
        }
        #else
        //Still need the tight-tight for pileup masking purposes.
        Int_t ntrk1 = (*vxnbc_nTracks)[nv];
        Int_t ntrk2 = (*vxnbc_nTracks)[nv2];
        for (std::vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
          if ((ntrk1 >= *nTrkCut) && (ntrk2 >= *nTrkCut)) {
            h_vtx_Dz_pLB[*nTrkCut][current_bcid]->Fill(dz, current_pLB);
          }
        }
        #endif
      }
    }
    (*vxnbc_nTracks)[nv] += extra_tracks;

    if (max_chi2ndf > 0. && (*vxnbc_chi2)[nv] / (*vxnbc_ndof)[nv] > max_chi2ndf) {
      continue;
    }

    h_vtx_nTracks_pLB[current_bcid]->Fill((*vxnbc_nTracks)[nv], current_pLB);

    for (std::vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
      if ((*vxnbc_nTracks)[nv] >= *nTrkCut) {
        NVtxAfterSplitCorrection[*nTrkCut][current_bcid][current_pLB]++;
        m_isTightAfterSplitCorrection[*nTrkCut] = true;
        h_vtx_z_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_z)[nv], current_pLB);

        #ifdef VERBOSE
        //Less important histograms
        h_vtx_x_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_x)[nv],current_pLB);
        h_vtx_y_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_y)[nv],current_pLB);

        h_vtx_xx_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_cov_x)[nv],current_pLB);
        h_vtx_yy_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_cov_y)[nv],current_pLB);
        h_vtx_zz_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_cov_z)[nv],current_pLB);
        h_vtx_xy_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_cov_xy)[nv],current_pLB);
        h_vtx_xz_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_cov_xz)[nv],current_pLB);
        h_vtx_yz_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_cov_yz)[nv],current_pLB);

        h_vtx_chi2_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_chi2)[nv],current_pLB);
        h_vtx_ndof_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_ndof)[nv],current_pLB);
        h_vtx_chi2ndof_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_chi2)[nv]/(*vxnbc_ndof)[nv],current_pLB);
        h_vtx_px_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_px)[nv],current_pLB);
        h_vtx_py_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_py)[nv],current_pLB);
        h_vtx_pz_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_pz)[nv],current_pLB);
        h_vtx_E_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_E)[nv],current_pLB);
        h_vtx_sumPt_pLB[*nTrkCut][current_bcid]->Fill((*vxnbc_sumPt)[nv],current_pLB);

        if ((*trig_L1_TAV)[1] == 268435456) {
          h_L1_BGRP7_pLB[*nTrkCut][current_bcid]->Fill(current_pLB);
        }
        if ((*trig_L1_TAV)[3] == 2147483648) {
          h_L1_3J10_pLB[*nTrkCut][current_bcid]->Fill(current_pLB);
        }

        for (int itTrk=0; itTrk < (*vxnbc_nTracks)[nv]; itTrk++) {
          int idx = (*vxnbc_trk_index)[nv][itTrk];
          h_trk_pt_pLB[*nTrkCut][current_bcid]->Fill((*trk_pt)[idx], current_pLB);
          h_trk_eta_pLB[*nTrkCut][current_bcid]->Fill((*trk_eta)[idx], current_pLB);
          h_trk_chi2ndf_pLB[*nTrkCut][current_bcid]->Fill((*trk_chi2)[idx] / (*trk_ndof)[idx], current_pLB);
          h_trk_chi2_pLB[*nTrkCut][current_bcid]->Fill((*trk_chi2)[idx], current_pLB);
          h_trk_ndf_pLB[*nTrkCut][current_bcid]->Fill((*trk_ndof)[idx], current_pLB);
          h_trk_nBLHits_pLB[*nTrkCut][current_bcid]->Fill((*trk_nBLHits)[idx], current_pLB);
          h_trk_nPixHits_pLB[*nTrkCut][current_bcid]->Fill((*trk_nPixHits)[idx], current_pLB);
          h_trk_nSCTHits_pLB[*nTrkCut][current_bcid]->Fill((*trk_nSCTHits)[idx], current_pLB);
          h_trk_nTRTHits_pLB[*nTrkCut][current_bcid]->Fill((*trk_nTRTHits)[idx], current_pLB);
          h_trk_nTRTHighTHits_pLB[*nTrkCut][current_bcid]->Fill((*trk_nTRTHighTHits)[idx], current_pLB);
          h_trk_nPixHoles_pLB[*nTrkCut][current_bcid]->Fill((*trk_nPixHoles)[idx], current_pLB);
          h_trk_nSCTHoles_pLB[*nTrkCut][current_bcid]->Fill((*trk_nSCTHoles)[idx], current_pLB);
          h_trk_nTRTHoles_pLB[*nTrkCut][current_bcid]->Fill((*trk_nTRTHoles)[idx], current_pLB);
          h_trk_nHits_pLB[*nTrkCut][current_bcid]->Fill((*trk_nHits)[idx], current_pLB);
          //h_trk_nHoles_pLB[*nTrkCut][current_bcid]->Fill((*trk_nHoles)[idx], current_pLB);
        }
        #endif
      }
    } // End NTrkCut loop
  } // End nv loop

  for (std::vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    if (m_isTightBeforeSplitCorrection[*nTrkCut]) {
      NEvtBeforeSplitCorrection[*nTrkCut][current_bcid][current_pLB]++;
    }
    if (m_isTightAfterSplitCorrection[*nTrkCut]) {
      NEvtAfterSplitCorrection[*nTrkCut][current_bcid][current_pLB]++;
    }
  }

  #ifdef VERBOSE
  h_evt_nvtx_before_remerge_pLB[current_bcid]->Fill(vxnbc_n, current_pLB);
  h_evt_nvtx_after_remerge_pLB[current_bcid]->Fill(vxnbc_n - vtx_skip_list.size(), current_pLB);
  #endif
  return kTRUE;
}

void AnaVtxTree::SlaveTerminate() {

  if (dumpTxtFileName != "") {
    //close dumping file
    dumpTxtFile.close();
  }

  // Write and save VdM tree. Put it in a different file, just in case saving seg-faults.
  std::cout << "[AnaVtxTree] INFO: Saving VdM TTree to: " << (outputFileName + TString("_tree.root")) << std::endl;
  TFile *f_outputTree = TFile::Open((outputFileName + TString("_tree.root")), "RECREATE");


  Int_t entry_pLB;
  Int_t entry_BunchId;
  Float_t entry_NVtxBeforeSplitCorrection[nTrkCuts.size()];
  Float_t entry_NVtxBeforeSplitCorrection_pertrack[nTrkCuts.size()];
  Float_t entry_NEvtBeforeSplitCorrection[nTrkCuts.size()];
  Float_t entry_NVtxAfterSplitCorrection[nTrkCuts.size()];
  Float_t entry_NEvtAfterSplitCorrection[nTrkCuts.size()];
  Float_t entry_NTrig;

  t_vdm = new TTree("t_vdm", "t_vdm");
  t_vdm->Branch( "pLB", &entry_pLB,"pLB/I");
  t_vdm->Branch( "BunchId",&entry_BunchId,"BunchId/I");
  t_vdm->Branch( "NTrig", &entry_NTrig,"NTrig/F");
  
  for (unsigned int i=0; i<nTrkCuts.size(); ++i) {
    const int trackCut( nTrkCuts.at(i) );
    t_vdm->Branch( TString::Format("NVtx%d_BSC",trackCut),
                   &(entry_NVtxBeforeSplitCorrection[i]),
                   TString::Format("NVtx%d_BSC/F",trackCut)
                 );
    t_vdm->Branch( TString::Format("NVtx%d_BSC_pertrack",trackCut),
                   &(entry_NVtxBeforeSplitCorrection_pertrack[i]),
                   TString::Format("NVtx%d_BSC/F",trackCut)
                 );
    t_vdm->Branch( TString::Format("NEvt%d_BSC",trackCut),
                   &(entry_NEvtBeforeSplitCorrection[i]),
                   TString::Format("NEvt%d_BSC/F",trackCut)
                 );
    t_vdm->Branch( TString::Format("NVtx%d_ASC",trackCut),
                   &(entry_NVtxAfterSplitCorrection[i]),
                   TString::Format("NVtx%d_ASC/F",trackCut)
                 );
    t_vdm->Branch( TString::Format("NEvt%d_ASC",trackCut),
                   &(entry_NEvtAfterSplitCorrection[i]),
                   TString::Format("NEvt%d_ASC/F",trackCut)
                 );
  }
  for (std::vector<Int_t>::iterator bcid = bcidListCollisions.begin(); bcid != bcidListCollisions.end(); ++bcid) {
    for (entry_pLB = pLBmin; entry_pLB <= pLBmax; entry_pLB++) {
      entry_BunchId = *bcid;
      for (unsigned int i=0; i<nTrkCuts.size(); i++) { 
        entry_NVtxBeforeSplitCorrection[i] = NVtxBeforeSplitCorrection[(nTrkCuts.at(i))][*bcid][entry_pLB];
        entry_NVtxBeforeSplitCorrection_pertrack[i] = NVtxBeforeSplitCorrection_pertrack[(nTrkCuts.at(i))];
        entry_NEvtBeforeSplitCorrection[i] = NEvtBeforeSplitCorrection[(nTrkCuts.at(i))][*bcid][entry_pLB];
        entry_NVtxAfterSplitCorrection[i] = NVtxAfterSplitCorrection[(nTrkCuts.at(i))][*bcid][entry_pLB];
        entry_NEvtAfterSplitCorrection[i] = NEvtAfterSplitCorrection[(nTrkCuts.at(i))][*bcid][entry_pLB];
      }
      entry_NTrig = NTrig[*bcid][entry_pLB];
      t_vdm->Fill();
    }
  }
     
  t_vdm->Write();
  f_outputTree->Close();

  // --- Print/save results
  std::cout << "[AnaVtxTree] INFO: Saving histograms to: " << (outputFileName + TString(".root")) << std::endl;
  TFile *outputNtp = TFile::Open((outputFileName + TString(".root")), "RECREATE");


  saveHistograms();

  // --- Close output files
  outputNtp->Close();
  outputTxt.close();

  std::cout<<"[AnaVtxTree] INFO: All Done." << std::endl;

}

void AnaVtxTree::saveHistograms() {
  //save all histograms in the hist/ subdirectory
  TDirectory *saveDir = gDirectory;
  gDirectory->mkdir("hist"); gDirectory->cd("hist");
  for( std::vector<TH1*>::iterator hist = HistogramHelper::histQueue.begin(); hist != HistogramHelper::histQueue.end(); ++hist ) {
    (*hist)->Write();
  }
  gDirectory = saveDir; gDirectory->cd();
}

void AnaVtxTree::LoadAndAddPseudoLB(TString fileName) {
  std::cout << "[AnaVtxTree] INFO: Opening fileName " << fileName.Data() << std::endl;
  ifstream tsFile(fileName.Data());
  if (not tsFile.is_open()) {
    std::cerr << "[AnaVtxTree] ERROR: Unable to open requested file with pseudo-LB info: " << fileName.Data() << std::endl;
    return;
  }

  Int_t pLB;
  Double_t ts_start, ts_end;

  std::vector<Int_t> pLB_list;
  std::vector<Double_t> ts_start_list;
  std::vector<Double_t> ts_end_list;

  while (tsFile.good()) {
    std::string line;
    getline(tsFile, line);
    if (line.empty()) {
      continue;
    }
    int prec=std::numeric_limits<long double>::digits10;
    std::stringstream istr;
    istr.precision(prec);
    istr << line;
    if (line.find_first_of("PseudoLB") != std::string::npos) {
      continue;  //header
    }
    istr >> pLB >> ts_start >> ts_end;
    std::pair<Double_t, Double_t> ts = std::make_pair<Double_t, Double_t>(ts_start, ts_end);
    pseudoLB_Timestamps.push_back(std::make_pair<Int_t, std::pair<Double_t, Double_t> >(pLB, ts) );
    //std::cout << "LB: " << pLB << ": " << ts_start << " - " << ts_end << std::endl;
    pLB_list.push_back(pLB);
    ts_start_list.push_back(ts_start);
    ts_end_list.push_back(ts_end);
  }

  sort(pLB_list.begin(), pLB_list.end());
  Int_t current_pLBmin = *(pLB_list.begin());
  Int_t current_pLBmax = *(pLB_list.end()-1);
  pLB_boundaries.push_back(std::make_pair<Int_t, Int_t>(current_pLBmin, current_pLBmax));
  if (current_pLBmin < pLBmin) {
    pLBmin = current_pLBmin;
  }
  if (current_pLBmax > pLBmax) {
    pLBmax = current_pLBmax;
  }
  tsFile.close();
  sort(ts_start_list.begin(), ts_start_list.end());
  sort(ts_end_list.begin(), ts_end_list.end());
  timestamp_boundaries.push_back(std::make_pair<Double_t, Double_t>(*(ts_start_list.begin()), *(ts_end_list.end()-1)));
}

int AnaVtxTree::timestampInScan(Double_t ts) {
  int scan = 0;
  for (std::vector<std::pair<Double_t, Double_t> >::iterator ts_boundary = timestamp_boundaries.begin(); ts_boundary != timestamp_boundaries.end(); ++ts_boundary) {
    if (ts >= (*ts_boundary).first && ts <= (*ts_boundary).second) {
      scan = distance(timestamp_boundaries.begin(), ts_boundary) + 1;
    }
  }
  return scan;
}

//Determine if vertexToCheck is a good vertex, if -1 check all verttices and find the first tagged good one
//qualityVertexVersion: x.y.z -- y and/or z can be omitted if 0
// x = general quality requirenents:
//    1 = standard
//    2 = DR requirement on 2-,3-track vertices
//    3 = Require at least N (default=5, otherwise=z) tracks pointing to the vertex
// y = additional vertex requirements
//    0 = none
//    1 = Accept only events where tagged vertex has N(=z) tracks
// z = Optional parameter for a given "y" mode
// examples:
// 1: standard
// 2.1.5: DR req. on 2-,3-trk vertices + tagged vertex required with at least 5 tracks
int AnaVtxTree::isGoodVertex(int vertexToCheck) {
  int taggedVertex = -1;

  if (qualityVertexVersion == 1) {
    //standard-definition
    Double_t minNTracks = qualityVertexParameter;
    if (minNTracks == 0) {
      minNTracks = 2;  //default value
    }
    if (vertexToCheck == -1) {
      //check only algorithm-tagged vertex
      if ((*vxnbc_nTracks)[0] >= minNTracks && (*vxnbc_cov_x)[0]!=0 && (*vxnbc_cov_y)[0]!=0 && vxnbc_n > 1) {
        taggedVertex=0;
      }
    } else {
      if ((*vxnbc_nTracks)[vertexToCheck] >= minNTracks && (*vxnbc_cov_x)[vertexToCheck]!=0 && (*vxnbc_cov_y)[vertexToCheck]!=0 && vxnbc_n > 1) {
        taggedVertex=vertexToCheck;
      }
    }
  } else if (qualityVertexVersion == 2) {
    //require vertex probability
    for (int ivx=0; ivx < vxnbc_n - 1; ++ivx) {
      if (vertexToCheck > -1 && ivx != vertexToCheck) {
        continue;  //next vertex
      }
      double probVertex = TMath::Prob((*vxnbc_chi2)[ivx], (*vxnbc_ndof)[ivx]);
      //      double dz=0.0;
      //for (std::vector<int>::iterator nvtrk = vxnbc_trk_index->at(ivx).begin(); nvtrk != vxnbc_trk_index->at(ivx).end(); ++nvtrk) {
      //  dr += TMath::Power()
      //}
      double Dr2=-999.;
      double Dr3=-999.;
      if (vxnbc_nTracks->at(ivx) == 2) {
        //Properties for events with 2 tracks at the Tagged vertex
        //evaluate delta d0 for vertices with 2 tracks
        if (vxnbc_trk_index->at(ivx).size() == 2) {
          // (Delta r)^2 = (d0)^2 + (d'0)^2 - 2(d0)(d'0)cos(Delta phi)
          double d01 = trk_d0_wrtPV->at(vxnbc_trk_index->at(ivx).at(0));
          double d02 = trk_d0_wrtPV->at(vxnbc_trk_index->at(ivx).at(1));
          double phi1 = trk_phi_wrtPV->at(vxnbc_trk_index->at(ivx).at(0));
          double phi2 = trk_phi_wrtPV->at(vxnbc_trk_index->at(ivx).at(1));
          double dphi = TMath::Abs(TVector2::Phi_mpi_pi(phi1 - phi2));
          Dr2 = TMath::Power(d01,2) + TMath::Power(d02,2) - 2*d01*d02*TMath::Cos(dphi);
          Dr2 = TMath::Sqrt(Dr2);
        }
      } else if (vxnbc_nTracks->at(ivx) == 3) {
        if (vxnbc_trk_index->at(ivx).size() == 3) {
          //evaluate average distance
          double d01 = trk_d0_wrtPV->at(vxnbc_trk_index->at(ivx).at(0));
          double d02 = trk_d0_wrtPV->at(vxnbc_trk_index->at(ivx).at(1));
          double d03 = trk_d0_wrtPV->at(vxnbc_trk_index->at(ivx).at(2));
          double phi1 = trk_phi_wrtPV->at(vxnbc_trk_index->at(ivx).at(0));
          double phi2 = trk_phi_wrtPV->at(vxnbc_trk_index->at(ivx).at(1));
          double phi3 = trk_phi_wrtPV->at(vxnbc_trk_index->at(ivx).at(2));
          double dphi12 = TMath::Abs(TVector2::Phi_mpi_pi(phi1 - phi2));
          double dphi13 = TMath::Abs(TVector2::Phi_mpi_pi(phi1 - phi3));
          double dphi23 = TMath::Abs(TVector2::Phi_mpi_pi(phi1 - phi2));
          Dr3 = TMath::Sqrt(TMath::Power(d01,2) + TMath::Power(d02,2) - 2*d01*d02*TMath::Cos(dphi12));
          Dr3 += TMath::Sqrt(TMath::Power(d01,2) + TMath::Power(d03,2) - 2*d01*d02*TMath::Cos(dphi13));
          Dr3 += TMath::Sqrt(TMath::Power(d02,2) + TMath::Power(d03,2) - 2*d01*d02*TMath::Cos(dphi23));
          Dr3 = Dr3 / 3;
        }
      }

      if ( ((*vxnbc_nTracks)[ivx] == 2 && probVertex > 0.01 && Dr2 < 2.5) ||
           ((*vxnbc_nTracks)[ivx] == 3 && probVertex > 1E-08 && Dr3 < 2.5) ||
           ((*vxnbc_nTracks)[ivx] > 3  && probVertex > 0)) {
        taggedVertex = ivx;
        break; // Vertex with highest sumPt2 with quality requirements
      }
    }
  } else if (qualityVertexVersion == 3) {
    Double_t minNTracks = qualityVertexParameter;
    if (minNTracks == 0) {
      minNTracks = 5;  //default value
    }
    for (int ivx=0; ivx < vxnbc_n - 1; ++ivx) {
      if (vertexToCheck > -1 && ivx != vertexToCheck) {
        continue;  //next vertex
      }
      if ((*vxnbc_nTracks)[ivx] >= minNTracks && (*vxnbc_cov_x)[ivx]!=0 && (*vxnbc_cov_y)[ivx]!=0) {
        taggedVertex = ivx;
        break;
      }
    }
  }

  if (qualityVertexModifier == 1) {
    //if tagged vertex has a number of tracks != qualityVertexParameter, un-tag the vertex
    if ((*vxnbc_nTracks)[taggedVertex] != qualityVertexParameter) {
      taggedVertex = -1;
    }
  }

  return taggedVertex;
}

void AnaVtxTree::SetDumpTxtFileName(TString name) {
  dumpTxtFileName = name;
}

void AnaVtxTree::SetTriggerFilter(TString p_triggerName, TChain *p_trigMetaDataTree) {
  triggerName = p_triggerName;
  if (p_trigMetaDataTree) {
    trigMetaDataTree = p_trigMetaDataTree;
  }
}

void AnaVtxTree::SetMaxChi2Ndf(Float_t p_max_chi2ndf) {

  std::cout << "[AnaVtxTree] INFO: Setting maximum chi2/ndf = " << p_max_chi2ndf << std::endl;
  max_chi2ndf = p_max_chi2ndf;

}

void AnaVtxTree::SetupLBInfo( const unsigned int &runNumber ) {
  // Run-by-run setup
  // - Specify BCIDs
  // - Load pLB timestamps

  // Load timestamp and pLB info
  std::cout << "[AnaVtxTree] INFO: Loading new list of pseudo LB info for run " << runNumber << std::endl;
  pseudoLB_Timestamps.clear();

  if (runNumber == 188949) { // 2011 mu-scan
    // Set pLB and BCID info
    use_plbs = true;
    bcidListCollisions.push_back(200);
    bcidListCollisions.push_back(999);
    // Load timestamps
    TString pathToPLBTimeStamps=GlobalSettings::path_timestamps+"188949/scan/";
    LoadAndAddPseudoLB((pathToPLBTimeStamps+TString("timestamps1.dat")).Data()); //add list of pseudo-LB.
    LoadAndAddPseudoLB((pathToPLBTimeStamps+TString("timestamps2.dat")).Data()); //add list of pseudo-LB.

  } else if (runNumber == 188951) { // 2011 mu-scan
    // Set pLB and BCID info
    use_plbs = true;
    bcidListCollisions.push_back(200);
    bcidListCollisions.push_back(999);
    // Load timestamps
    TString pathToPLBTimeStamps=GlobalSettings::path_timestamps+"188951/scan/";
    LoadAndAddPseudoLB((pathToPLBTimeStamps+TString("timestamps3.dat")).Data()); //add list of pseudo-LB.

  } else if (runNumber == 182013) { // 2011 VdM-scan
    // Set pLB and BCID info
    use_plbs = true;
    bcidListCollisions.push_back(81);
    bcidListCollisions.push_back(867);
    bcidListCollisions.push_back(2752);
    // Load timestamps
    TString pathToPLBTimeStamps=GlobalSettings::path_timestamps+"182013/run/";
    LoadAndAddPseudoLB((pathToPLBTimeStamps+TString("all_timestamps.dat")).Data()); //add list of pseudo-LB.

  } else if (runNumber == 191373) { // 2011/2012 something?
    // Set pLB and BCID info
    use_plbs = true;
    bcidListCollisions.push_back(1);
    // Load timestamps
    TString pathToPLBTimeStamps=GlobalSettings::path_timestamps+"191373/";
    LoadAndAddPseudoLB((pathToPLBTimeStamps+TString("191373_timestamps.dat")).Data()); //add list of pseudo-LB.

  } else if (runNumber == 200805) { // 2012 April muScan
    // Set pLB and BCID info
    use_plbs = false;
    bcidListCollisions.push_back(1);
    bcidListCollisions.push_back(149);
    bcidListCollisions.push_back(334);
    // Set lumiblocks by hand
    pLBmin = 1;
    pLBmax = 415;

  } else if (runNumber == 201351) { // 2012 April VdM scan (scan I, II & III)
    // Set pLB and BCID info
    use_plbs = true;
    bcidListCollisions.push_back(1);
    bcidListCollisions.push_back(241);
    bcidListCollisions.push_back(2881);
    bcidListCollisions.push_back(3121);
    // Load timestamps
    TString pathToPLBTimeStamps=GlobalSettings::path_timestamps+"201351/run/";
    LoadAndAddPseudoLB((pathToPLBTimeStamps+TString("timestamps.dat")).Data()); //add list of pseudo-LB.
    
  }  else if (runNumber == 207216) { // 2012 July VdM scan 4, 5, 6
    // Set pLB and BCID info
    use_plbs = true;
    bcidListCollisions.push_back(1);
    bcidListCollisions.push_back(721);
    bcidListCollisions.push_back(1821);
    // Load timestamps
    TString pathToPLBTimeStamps=GlobalSettings::path_timestamps+"207216/run/";
    LoadAndAddPseudoLB((pathToPLBTimeStamps+TString("timestamps.dat")).Data()); //add list of pseudo-LB. 
  
 }  else if (runNumber == 207219) { // 2012 July VdM scan 8
    // Set pLB and BCID info
    use_plbs = true;
    bcidListCollisions.push_back(1);
    bcidListCollisions.push_back(721);
    bcidListCollisions.push_back(1821);
    // Load timestamps
    TString pathToPLBTimeStamps=GlobalSettings::path_timestamps+"207219/run/";
    LoadAndAddPseudoLB((pathToPLBTimeStamps+TString("timestamps.dat")).Data()); //add list of pseudo-LB. 
  
 } else if (runNumber == 214984) { // 2012 November VdM scans 11 and 14
    // Set pLB and BCID info
    use_plbs = true;
    bcidListCollisions.push_back(1);
    bcidListCollisions.push_back(2361);
    bcidListCollisions.push_back(2881);
    // Load timestamps
    TString pathToPLBTimeStamps=GlobalSettings::path_timestamps+"214984/run/";
    LoadAndAddPseudoLB((pathToPLBTimeStamps+TString("timestamps.dat")).Data()); //add list of pseudo-LB.
    
 } else if (runNumber == 215021) { // 2012 November VdM scan 15
    // Set pLB and BCID info
    use_plbs = true;
    bcidListCollisions.push_back(1);
    bcidListCollisions.push_back(2361);
    bcidListCollisions.push_back(2881);
    // Load timestamps
    TString pathToPLBTimeStamps=GlobalSettings::path_timestamps+"215021/run/";
    LoadAndAddPseudoLB((pathToPLBTimeStamps+TString("timestamps.dat")).Data()); //add list of pseudo-LB.  
} else if (runNumber == 216399) { // 2012 25ns run, 97 bunches
    use_plbs = false;
    // Set lumiblocks by hand
    pLBmin = 320;
    pLBmax = 430;
    // Set pLB and BCID info
    bcidListCollisions.push_back(1);
    for (int i = 301; i<=348; i++){
      bcidListCollisions.push_back(i);
    }
    for (int j = 1195; j<=1242; j++){
      bcidListCollisions.push_back(j);
    } 
} else if (runNumber == 216416) { // 2012 25ns run, 187 bunches
    use_plbs = false;
    // Set lumiblocks by hand
    pLBmin = 120;
    pLBmax = 260;
    // Set pLB and BCID info
    for (int i = 301; i<=348 ; i++){
      bcidListCollisions.push_back(i);
    }
    for (int j = 1195; j<=1242 ; j++){
      bcidListCollisions.push_back(j);
    }
    for (int k = 2089; k<=2136 ; k++){
      bcidListCollisions.push_back(k);
    }
    for (int l = 2977; l<=3018 ; l++){
      bcidListCollisions.push_back(l);
    }
} else if (runNumber == 216432) { // 2012 25ns run, 373 bunches
    use_plbs = false;
    // Set lumiblocks by hand
    pLBmin = 500;
    pLBmax = 1100;
    // Set pLB and BCID info
    for (int i = 301; i<=348 ; i++){
      bcidListCollisions.push_back(i);
    }
    for (int j = 358; j<=405 ; j++){
      bcidListCollisions.push_back(j);
    }
    for (int k = 1195; k<=1242 ; k++){
      bcidListCollisions.push_back(k);
    }
    for (int l = 1252; l<=1299 ; l++){
      bcidListCollisions.push_back(l);
    }
    for (int m = 2089; m<=2136 ; m++){
      bcidListCollisions.push_back(m);
    }
    for (int n = 2146; n<=2193 ; n++){
      bcidListCollisions.push_back(n);
    }
    for (int o = 2977; o<=3018 ; o++){
      bcidListCollisions.push_back(o);
    }
    for (int p = 3034; p<=3075 ; p++){
      bcidListCollisions.push_back(p);
    }
} /*else if (runNumber == 206955 || runNumber == 206971 || runNumber == 207044 || runNumber == 207046 ) { // 2012 50ns run, 1368 bunches 18Nov: Deleted 206962 to treat it BCID blindly
    use_plbs = false;
    // Set lumiblocks by hand
    pLBmin = 0;
    pLBmax = 900;
    int list [1368]={66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126, 128, 130, 132, 134, 136, 146, 148, 150, 152, 154, 156, 158, 160, 162, 164, 166, 168, 170, 172, 174, 176, 178, 180, 182, 184, 186, 188, 190, 192, 194, 196, 198, 200, 202, 204, 206, 208, 210, 212, 214, 216, 226, 228, 230, 232, 234, 236, 238, 240, 242, 244, 246, 248, 250, 252, 254, 256, 258, 260, 262, 264, 266, 268, 270, 272, 274, 276, 278, 280, 282, 284, 286, 288, 290, 292, 294, 296, 306, 308, 310, 312, 314, 316, 318, 320, 322, 324, 326, 328, 330, 332, 334, 336, 338, 340, 342, 344, 346, 348, 350, 352, 354, 356, 358, 360, 362, 364, 366, 368, 370, 372, 374, 376, 413, 415, 417, 419, 421, 423, 425, 427, 429, 431, 433, 435, 437, 439, 441, 443, 445, 447, 449, 451, 453, 455, 457, 459, 461, 463, 465, 467, 469, 471, 473, 475, 477, 479, 481, 483, 493, 495, 497, 499, 501, 503, 505, 507, 509, 511, 513, 515, 517, 519, 521, 523, 525, 527, 529, 531, 533, 535, 537, 539, 541, 543, 545, 547, 549, 551, 553, 555, 557, 559, 561, 563, 573, 575, 577, 579, 581, 583, 585, 587, 589, 591, 593, 595, 597, 599, 601, 603, 605, 607, 609, 611, 613, 615, 617, 619, 621, 623, 625, 627, 629, 631, 633, 635, 637, 639, 641, 643, 653, 655, 657, 659, 661, 663, 665, 667, 669, 671, 673, 675, 677, 679, 681, 683, 685, 687, 689, 691, 693, 695, 697, 699, 701, 703, 705, 707, 709, 711, 713, 715, 717, 719, 721, 723, 773, 775, 777, 779, 781, 783, 785, 787, 789, 791, 793, 795, 797, 799, 801, 803, 805, 807, 809, 811, 813, 815, 817, 819, 821, 823, 825, 827, 829, 831, 833, 835, 837, 839, 841, 843, 853, 855, 857, 859, 861, 863, 865, 867, 869, 871, 873, 875, 877, 879, 881, 883, 885, 887, 889, 891, 893, 895, 897, 899, 901, 903, 905, 907, 909, 911, 913, 915, 917, 919, 921, 923, 960, 962, 964, 966, 968, 970, 972, 974, 976, 978, 980, 982, 984, 986, 988, 990, 992, 994, 996, 998, 1000, 1002, 1004, 1006, 1008, 1010, 1012, 1014, 1016, 1018, 1020, 1022, 1024, 1026, 1028, 1030, 1040, 1042, 1044, 1046, 1048, 1050, 1052, 1054, 1056, 1058, 1060, 1062, 1064, 1066, 1068, 1070, 1072, 1074, 1076, 1078, 1080, 1082, 1084, 1086, 1088, 1090, 1092, 1094, 1096, 1098, 1100, 1102, 1104, 1106, 1108, 1110, 1120, 1122, 1124, 1126, 1128, 1130, 1132, 1134, 1136, 1138, 1140, 1142, 1144, 1146, 1148, 1150, 1152, 1154, 1156, 1158, 1160, 1162, 1164, 1166, 1168, 1170, 1172, 1174, 1176, 1178, 1180, 1182, 1184, 1186, 1188, 1190, 1200, 1202, 1204, 1206, 1208, 1210, 1212, 1214, 1216, 1218, 1220, 1222, 1224, 1226, 1228, 1230, 1232, 1234, 1236, 1238, 1240, 1242, 1244, 1246, 1248, 1250, 1252, 1254, 1256, 1258, 1260, 1262, 1264, 1266, 1268, 1270, 1307, 1309, 1311, 1313, 1315, 1317, 1319, 1321, 1323, 1325, 1327, 1329, 1331, 1333, 1335, 1337, 1339, 1341, 1343, 1345, 1347, 1349, 1351, 1353, 1355, 1357, 1359, 1361, 1363, 1365, 1367, 1369, 1371, 1373, 1375, 1377, 1387, 1389, 1391, 1393, 1395, 1397, 1399, 1401, 1403, 1405, 1407, 1409, 1411, 1413, 1415, 1417, 1419, 1421, 1423, 1425, 1427, 1429, 1431, 1433, 1435, 1437, 1439, 1441, 1443, 1445, 1447, 1449, 1451, 1453, 1455, 1457, 1467, 1469, 1471, 1473, 1475, 1477, 1479, 1481, 1483, 1485, 1487, 1489, 1491, 1493, 1495, 1497, 1499, 1501, 1503, 1505, 1507, 1509, 1511, 1513, 1515, 1517, 1519, 1521, 1523, 1525, 1527, 1529, 1531, 1533, 1535, 1537, 1547, 1549, 1551, 1553, 1555, 1557, 1559, 1561, 1563, 1565, 1567, 1569, 1571, 1573, 1575, 1577, 1579, 1581, 1583, 1585, 1587, 1589, 1591, 1593, 1595, 1597, 1599, 1601, 1603, 1605, 1607, 1609, 1611, 1613, 1615, 1617, 1667, 1669, 1671, 1673, 1675, 1677, 1679, 1681, 1683, 1685, 1687, 1689, 1691, 1693, 1695, 1697, 1699, 1701, 1703, 1705, 1707, 1709, 1711, 1713, 1715, 1717, 1719, 1721, 1723, 1725, 1727, 1729, 1731, 1733, 1735, 1737, 1747, 1749, 1751, 1753, 1755, 1757, 1759, 1761, 1763, 1765, 1767, 1769, 1771, 1773, 1775, 1777, 1779, 1781, 1783, 1785, 1787, 1789, 1791, 1793, 1795, 1797, 1799, 1801, 1803, 1805, 1807, 1809, 1811, 1813, 1815, 1817, 1854, 1856, 1858, 1860, 1862, 1864, 1866, 1868, 1870, 1872, 1874, 1876, 1878, 1880, 1882, 1884, 1886, 1888, 1890, 1892, 1894, 1896, 1898, 1900, 1902, 1904, 1906, 1908, 1910, 1912, 1914, 1916, 1918, 1920, 1922, 1924, 1934, 1936, 1938, 1940, 1942, 1944, 1946, 1948, 1950, 1952, 1954, 1956, 1958, 1960, 1962, 1964, 1966, 1968, 1970, 1972, 1974, 1976, 1978, 1980, 1982, 1984, 1986, 1988, 1990, 1992, 1994, 1996, 1998, 2000, 2002, 2004, 2014, 2016, 2018, 2020, 2022, 2024, 2026, 2028, 2030, 2032, 2034, 2036, 2038, 2040, 2042, 2044, 2046, 2048, 2050, 2052, 2054, 2056, 2058, 2060, 2062, 2064, 2066, 2068, 2070, 2072, 2074, 2076, 2078, 2080, 2082, 2084, 2094, 2096, 2098, 2100, 2102, 2104, 2106, 2108, 2110, 2112, 2114, 2116, 2118, 2120, 2122, 2124, 2126, 2128, 2130, 2132, 2134, 2136, 2138, 2140, 2142, 2144, 2146, 2148, 2150, 2152, 2154, 2156, 2158, 2160, 2162, 2164, 2201, 2203, 2205, 2207, 2209, 2211, 2213, 2215, 2217, 2219, 2221, 2223, 2225, 2227, 2229, 2231, 2233, 2235, 2237, 2239, 2241, 2243, 2245, 2247, 2249, 2251, 2253, 2255, 2257, 2259, 2261, 2263, 2265, 2267, 2269, 2271, 2281, 2283, 2285, 2287, 2289, 2291, 2293, 2295, 2297, 2299, 2301, 2303, 2305, 2307, 2309, 2311, 2313, 2315, 2317, 2319, 2321, 2323, 2325, 2327, 2329, 2331, 2333, 2335, 2337, 2339, 2341, 2343, 2345, 2347, 2349, 2351, 2361, 2363, 2365, 2367, 2369, 2371, 2373, 2375, 2377, 2379, 2381, 2383, 2385, 2387, 2389, 2391, 2393, 2395, 2397, 2399, 2401, 2403, 2405, 2407, 2409, 2411, 2413, 2415, 2417, 2419, 2421, 2423, 2425, 2427, 2429, 2431, 2441, 2443, 2445, 2447, 2449, 2451, 2453, 2455, 2457, 2459, 2461, 2463, 2465, 2467, 2469, 2471, 2473, 2475, 2477, 2479, 2481, 2483, 2485, 2487, 2489, 2491, 2493, 2495, 2497, 2499, 2501, 2503, 2505, 2507, 2509, 2511, 2549, 2551, 2553, 2555, 2557, 2559, 2561, 2563, 2565, 2567, 2569, 2571, 2573, 2575, 2577, 2579, 2581, 2583, 2585, 2587, 2589, 2591, 2593, 2595, 2597, 2599, 2601, 2603, 2605, 2607, 2609, 2611, 2613, 2615, 2617, 2619, 2629, 2631, 2633, 2635, 2637, 2639, 2641, 2643, 2645, 2647, 2649, 2651, 2653, 2655, 2657, 2659, 2661, 2663, 2665, 2667, 2669, 2671, 2673, 2675, 2677, 2679, 2681, 2683, 2685, 2687, 2689, 2691, 2693, 2695, 2697, 2699, 2736, 2738, 2740, 2742, 2744, 2746, 2748, 2750, 2752, 2754, 2756, 2758, 2760, 2762, 2764, 2766, 2768, 2770, 2772, 2774, 2776, 2778, 2780, 2782, 2784, 2786, 2788, 2790, 2792, 2794, 2796, 2798, 2800, 2802, 2804, 2806, 2816, 2818, 2820, 2822, 2824, 2826, 2828, 2830, 2832, 2834, 2836, 2838, 2840, 2842, 2844, 2846, 2848, 2850, 2852, 2854, 2856, 2858, 2860, 2862, 2864, 2866, 2868, 2870, 2872, 2874, 2876, 2878, 2880, 2882, 2884, 2886, 2896, 2898, 2900, 2902, 2904, 2906, 2908, 2910, 2912, 2914, 2916, 2918, 2920, 2922, 2924, 2926, 2928, 2930, 2932, 2934, 2936, 2938, 2940, 2942, 2944, 2946, 2948, 2950, 2952, 2954, 2956, 2958, 2960, 2962, 2964, 2966, 2976, 2978, 2980, 2982, 2984, 2986, 2988, 2990, 2992, 2994, 2996, 2998, 3000, 3002, 3004, 3006, 3008, 3010, 3012, 3014, 3016, 3018, 3020, 3022, 3024, 3026, 3028, 3030, 3032, 3034, 3036, 3038, 3040, 3042, 3044, 3046, 3083, 3085, 3087, 3089, 3091, 3093, 3095, 3097, 3099, 3101, 3103, 3105, 3107, 3109, 3111, 3113, 3115, 3117, 3119, 3121, 3123, 3125, 3127, 3129, 3131, 3133, 3135, 3137, 3139, 3141, 3143, 3145, 3147, 3149, 3151, 3153, 3163, 3165, 3167, 3169, 3171, 3173, 3175, 3177, 3179, 3181, 3183, 3185, 3187, 3189, 3191, 3193, 3195, 3197, 3199, 3201, 3203, 3205, 3207, 3209, 3211, 3213, 3215, 3217, 3219, 3221, 3223, 3225, 3227, 3229, 3231, 3233, 3243, 3245, 3247, 3249, 3251, 3253, 3255, 3257, 3259, 3261, 3263, 3265, 3267, 3269, 3271, 3273, 3275, 3277, 3279, 3281, 3283, 3285, 3287, 3289, 3291, 3293, 3295, 3297, 3299, 3301, 3303, 3305, 3307, 3309, 3311, 3313, 3323, 3325, 3327, 3329, 3331, 3333, 3335, 3337, 3339, 3341, 3343, 3345, 3347, 3349, 3351, 3353, 3355, 3357, 3359, 3361, 3363, 3365, 3367, 3369, 3371, 3373, 3375, 3377, 3379, 3381, 3383, 3385, 3387, 3389, 3391, 3393 };
    std::cout << "Size of BCID list is" << sizeof(list)/sizeof(*list) << std::endl;
    // Set pLB and BCID info
    for (int i = 0; i<1368; i++){
      bcidListCollisions.push_back(list[i]);
    }
}*/ else { // Otherwise, we're in a physics run, so run BCID-blind.
    // Set pLB and BCID info
    use_plbs = false;
    physics_run = true;
    bcidListCollisions.push_back(0);
    // Set lumiblocks by hand
    pLBmin = 1; //changed from 1 to 0
    pLBmax = 1500; //changed from 1500 to 80 for MC Closure test
  }
  
  // Print pLB info
  nPLB = pLBmax - pLBmin + 1;
  std::cout << "/****************************************************/" << std::endl;
  std::cout << "/**            pLB Timestamp Summary               **/" << std::endl;
  std::cout << "/**     pLBmin = " << std::setw(5) << pLBmin << "                             **/" << std::endl;
  std::cout << "/**     pLBmax = " << std::setw(5) << pLBmax << "                             **/" << std::endl;
  std::cout << "/**     nPLB   = " << std::setw(5) << nPLB <<   "                             **/" << std::endl;
  std::cout << "/****************************************************/" << std::endl;
}

void AnaVtxTree::InitBCIDHists() {
  std::cout << "[AnaVtxTree] INFO: Initializing histograms" << std::endl;
  m_InitBcidHistograms = false;

TString h_bcid_name, h_bcid_title;

  // -- Maps <bcid, TH>
  for (std::vector<Int_t>::iterator ib = bcidListCollisions.begin(); ib != bcidListCollisions.end(); ++ib) {     
    
    h_privtx_z_pLB[*ib] = HistogramHelper::define2DHistogram( cstr("Primary vertex Z vs. pLB, BCID ",*ib),
                                             NBinsZAxis, LowZAxis, HighZAxis, nPLB, pLBmin, pLBmax+1,
                                             "Z (mm)", "pLB", true, cstr("PriVtxZpLB_BCID",*ib));
                                             
    h_events_pLB[*ib] = HistogramHelper::defineHistogram( cstr("NEvents per pLB, BCID ",*ib),
                                         nPLB, pLBmin, pLBmax+1,
                                         "pLB", "nEvents", true, cstr("NEvents_pLB_BCID",*ib) );

    h_vtx_nTracks_pLB[*ib] = HistogramHelper::define2DHistogram( cstr("Vertex NTrk, BCID ", *ib),
                             200., 0., 200., nPLB, pLBmin,pLBmax+1,
                             "nTracks", "pLB", true, cstr("VtxNTrk_BCID",*ib) );
                             

    #ifdef VERBOSE
    h_evt_nvtx_before_remerge_pLB[*ib] = HistogramHelper::define2DHistogram( cstr("NVtx distribution vs. pLB before remerging, BCID ",*ib),
                                         30, 0, 30, nPLB, pLBmin, pLBmax+1,
                                         "NVtx", "pLB", true, cstr("h_evt_nvtx_before_remerge_pLB_BCID",*ib) );

    h_evt_nvtx_after_remerge_pLB[*ib] = HistogramHelper::define2DHistogram(cstr("NVtx distribution vs. pLB after remerging, BCID ",*ib),
                                        30, 0, 30, nPLB, pLBmin, pLBmax+1,
                                        "NVtx", "pLB", true, cstr("h_evt_nvtx_after_remerge_pLB_BCID",*ib) );
    #endif
  }

  // -- Single histograms
  for (std::vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    h_vtx_dz_dzsig[*nTrkCut] = HistogramHelper::define2DHistogram(cstr("Vertex dz vs. dz signficance, NTrk >= ",*nTrkCut),
                               1000, -50, 50, 600, -30, 30,
                               "#Deltaz (mm)", "#Deltaz / #sigma_{#Deltaz}", true, cstr("h_vtx_dz_dzsig_NTrk",*nTrkCut) );
    h_nvtx_pertrack[*nTrkCut] = HistogramHelper::defineHistogram(cstr( "NVtx", str(", NTrk >= ",*nTrkCut) ),
                                  30, -0.5, 29.5, "NVtx", "", true, cstr( "NVtx", str("_NTrkCut_",*nTrkCut) ) );
  }


  // -- Maps <NTrkCut, <bcid, TH > >
  for (std::vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
    for (std::vector<Int_t>::iterator ib = bcidListCollisions.begin(); ib != bcidListCollisions.end(); ++ib) {
      h_vtx_z_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram( cstr( str("Vertex Z Position, BCID ",*ib), str("_NTrkCut",*nTrkCut) ),
                                   NBinsZAxis, LowZAxis, HighZAxis, nPLB, pLBmin,pLBmax+1,
                                   "Z (mm)", "pLB", true, cstr( str("VtxPosZpLB_BCID",*ib), str("_NTrkCut",*nTrkCut) ) );

      h_vtx_Dz_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr(str("Vertex Delta Z, Tight-Tight, BCID ",*ib),str(", NTrk >= ",*nTrkCut) ),
                                    NBinsDzAxis, LowDzAxis, HighDzAxis, nPLB, pLBmin,pLBmax+1,
                                    "#DeltaZ (mm)", "pLB", true, cstr( str("VtxDzTightTight_pLB_BCID",*ib), str("_NTrkCut",*nTrkCut) ));
                                    
                                    h_nvtx_pertrack_perbcid[*nTrkCut][*ib] = HistogramHelper::defineHistogram(cstr( str("NVtx, NTrk >= ",*nTrkCut), str("BCID ",*ib) ),
                                  30, -0.5, 29.5, "NVtx", "", true, cstr( str("NVtx_NTrkCut_",*nTrkCut), str("_BCID_",*ib) ) );
                                         
       h_nvtx_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("NVtx vs. pLB, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                  30, -0.5, 29.5, nPLB, pLBmin, pLBmax+1,
                                  "NVtx", "pLB", true, cstr( str("NVtxpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ) );
      //Less important histograms
      
      #ifdef VERBOSE
      h_vtx_Dz_LooseTight_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex Delta Z, Loose-Tight, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
          NBinsDzAxis, LowDzAxis, HighDzAxis, nPLB, pLBmin,pLBmax+1,
          "#DeltaZ (mm)", "pLB", true, cstr( str("VtxDzLooseTight_pLB_BCID",*ib), str("_NTrkCut",*nTrkCut) ));

      h_vtx_Dz_TightLoose_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex Delta Z, Tight-Loose, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
          NBinsDzAxis, LowDzAxis, HighDzAxis, nPLB, pLBmin,pLBmax+1,
          "#DeltaZ (mm)", "pLB", true, cstr( str("VtxDzTightLoose_pLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_vtx_Dz_LooseLoose_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex Delta Z, Loose-Loose, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
          NBinsDzAxis, LowDzAxis, HighDzAxis, nPLB, pLBmin,pLBmax+1,
          "#DeltaZ (mm)", "pLB", true, cstr( str("VtxDzLooseLoose_pLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_vtx_x_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex PosX vs. pLB, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
          NBinsXAxis, LowXAxis, HighXAxis,nPLB, pLBmin,pLBmax+1,
          "X (mm)", "pLB", true, cstr( str("VtxPosXpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_vtx_y_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex PosY vs. pLB, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
          NBinsYAxis, LowYAxis, HighYAxis,nPLB, pLBmin,pLBmax+1,
          "Y (mm)", "pLB", true, cstr( str("VtxPosYpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_vtx_xx_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex Cov XX, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                    NBinsXAxis, 0., 0.1, nPLB, pLBmin,pLBmax+1,
                                    "cov XX (mm)", "pLB", true, cstr( str("VtxCovXXpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_vtx_yy_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex Cov YY, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                    NBinsXAxis, 0., 0.1, nPLB, pLBmin,pLBmax+1,
                                    "cov YY (mm)", "pLB", true, cstr( str("VtxCovYYpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_vtx_zz_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex Cov ZZ, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                    NBinsZAxis, 0., 0.1, nPLB, pLBmin,pLBmax+1,
                                    "cov ZZ (mm)", "pLB", true, cstr( str("VtxCovZZpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_vtx_xy_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex Cov XY, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                    NBinsXAxis, 0., 0.1, nPLB, pLBmin,pLBmax+1,
                                    "cov XY (mm)", "pLB", true, cstr( str("VtxCovXYpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_vtx_xz_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex Cov XZ, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                    NBinsXAxis, 0., 0.1, nPLB, pLBmin,pLBmax+1,
                                    "cov XZ (mm)", "pLB", true, cstr( str("VtxCovXZpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_vtx_yz_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex Cov YZ, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                    NBinsXAxis, 0., 0.1, nPLB, pLBmin,pLBmax+1,
                                    "cov YZ (mm)", "pLB", true, cstr( str("VtxCovYZpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      //Vertex quality info
      h_vtx_chi2_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex Chi2, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                      200., 0., 200., nPLB, pLBmin,pLBmax+1,
                                      "#chi^{2}", "pLB", true, cstr( str("VtxChi2_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_vtx_ndof_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex Ndf, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                      200., 0., 200., nPLB, pLBmin,pLBmax+1,
                                      "NDF", "pLB", true, cstr( str("VtxNdf_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_vtx_chi2ndof_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex Chi2/Ndf, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                          200., 0., 5, nPLB, pLBmin,pLBmax+1,
                                          "#chi^{2}/NDF (mm)", "pLB", true, cstr( str("VtxChi2Ndf_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      //Energy and momentum
      h_vtx_px_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex X Momentum, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                    NBinsXAxis, LowPXAxis, HighPXAxis, nPLB, pLBmin,pLBmax+1,
                                    "p_{x} (MeV)", "pLB", true, cstr( str("VtxPXpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_vtx_py_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex Y Momentum, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                    NBinsYAxis, LowPYAxis, HighPYAxis, nPLB, pLBmin,pLBmax+1,
                                    "p_{y} (MeV)", "pLB", true, cstr( str("VtxPYpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_vtx_pz_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex Z Momentum, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                    NBinsZAxis, LowPZAxis, HighPZAxis, nPLB, pLBmin,pLBmax+1,
                                    "p_{z} (MeV)", "pLB", true, cstr( str("VtxPZpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_vtx_E_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex Energy, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                   1000, 0., 1000000, nPLB, pLBmin,pLBmax+1,
                                   "E (MeV)", "pLB", true, cstr( str("VtxEpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_vtx_sumPt_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Vertex Sum Pt, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                       1000, 0., 1000000, nPLB, pLBmin,pLBmax+1,
                                       "p_{T} (MeV)", "pLB", true, cstr( str("VtxSumPtpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));
                                  

      h_L1_BGRP7_pLB[*nTrkCut][*ib] = HistogramHelper::defineHistogram(cstr( str("Number of L1_BGRP7 triggers vs. pLB, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                      nPLB, pLBmin, pLBmax+1,
                                      "pLB", "NTrig_{L1_BGRP7}", true, cstr( str("L1_BGRP7_pLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_L1_3J10_pLB[*nTrkCut][*ib] = HistogramHelper::defineHistogram(cstr( str("Number of L1_3J10 triggers vs. pLB, BCID ",*ib), str(", NTrk >= ",*nTrkCut) ),
                                     nPLB, pLBmin, pLBmax+1,
                                     "pLB", "NTrig_{L1_3J10}", true, cstr( str("L1_3J10_pLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_trk_pt_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Track pT for tracks associated with vertices vs. pLB, BCID",*ib), str(", NTrk >= ",*nTrkCut) ),
                                    5000, 0., 100000, nPLB, pLBmin, pLBmax+1,
                                    "p_{T} (MeV)", "pLB", true, cstr( str("TrkPtpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_trk_eta_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Track eta for tracks associated with vertices vs. pLB, BCID",*ib), str(", NTrk >= ",*nTrkCut) ),
                                     100, -3., 3., nPLB, pLBmin, pLBmax+1,
                                     "#eta", "pLB", true, cstr( str("TrkEtapLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_trk_chi2ndf_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Track chi2/ndf for tracks associated with vertices vs. pLB, BCID",*ib), str(", NTrk >= ",*nTrkCut) ),
                                         100, -3., 3., nPLB, pLBmin, pLBmax+1,
                                         "#chi^{2}/NDF", "pLB", true, cstr( str("TrkChi2NdfpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_trk_chi2_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Track chi2 for tracks associated with vertices vs. pLB, BCID",*ib), str(", NTrk >= ",*nTrkCut) ),
                                      200, 0., 200., nPLB, pLBmin, pLBmax+1,
                                      "#chi^{2}", "pLB", true, cstr( str("TrkChi2pLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_trk_ndf_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Track ndf for tracks associated with vertices vs. pLB, BCID",*ib), str(", NTrk >= ",*nTrkCut) ),
                                     200, 0., 200., nPLB, pLBmin, pLBmax+1,
                                     "NDF", "pLB", true, cstr( str("TrkNdfpLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_trk_nBLHits_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Track nBLHits for tracks associated with vertices vs. pLB, BCID",*ib), str(", NTrk >= ",*nTrkCut) ),
                                         11, -0.5, 10.5, nPLB, pLBmin, pLBmax+1,
                                         "N_{hits-BL}", "pLB", true, cstr( str("TrkNBLHitspLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_trk_nPixHits_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Track nPixHits for tracks associated with vertices vs. pLB, BCID",*ib), str(", NTrk >= ",*nTrkCut) ),
                                          11, -0.5, 10.5, nPLB, pLBmin, pLBmax+1,
                                          "N_{hits-PIX}", "pLB", true, cstr( str("TrkNPIXHitspLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_trk_nSCTHits_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Track nSCTHits for tracks associated with vertices vs. pLB, BCID",*ib), str(", NTrk >= ",*nTrkCut) ),
                                          11, -0.5, 10.5,nPLB, pLBmin, pLBmax+1,
                                          "N_{hits-SCT}", "pLB", true, cstr( str("TrkNSCTHitspLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_trk_nTRTHits_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Track nTRTHits for tracks associated with vertices vs. pLB, BCID",*ib), str(", NTrk >= ",*nTrkCut) ),
                                          101, -0.5, 100.5, nPLB, pLBmin, pLBmax+1,
                                          "N_{hits-TRT}", "pLB", true,cstr( str("TrkNTRTHitspLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_trk_nTRTHighTHits_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Track nTRTHighTHits for tracks associated with vertices vs. pLB, BCID",*ib), str(", NTrk >= ",*nTrkCut) ),
                                          101, -0.5, 100.5, nPLB, pLBmin, pLBmax+1,
                                          "N_{hits-TRT-highT}", "pLB", true, cstr( str("TrkNTRTHighTHitspLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_trk_nPixHoles_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Track nPixHoles for tracks associated with vertices vs. pLB, BCID",*ib), str(", NTrk >= ",*nTrkCut) ),
                                           11, -0.5, 10.5, nPLB, pLBmin, pLBmax+1,
                                           "N_{holes-PIX}", "pLB", true, cstr( str("TrkNPixHolespLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_trk_nSCTHoles_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Track nSCTHoles for tracks associated with vertices vs. pLB, BCID",*ib), str(", NTrk >= ",*nTrkCut) ),
                                           11, -0.5, 10.5, nPLB, pLBmin, pLBmax+1,
                                           "N_{holes-SCT}", "pLB", true, cstr( str("TrkNSCTHolespLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_trk_nTRTHoles_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Track nTRTHoles for tracks associated with vertices vs. pLB, BCID",*ib), str(", NTrk >= ",*nTrkCut) ),
                                           101, -0.5, 100.5, nPLB, pLBmin, pLBmax+1,
                                           "N_{holes-TRT}", "pLB", true, cstr( str("TrkNTRTHolespLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_trk_nHits_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Track nHits for tracks associated with vertices vs. pLB, BCID",*ib), str(", NTrk >= ",*nTrkCut) ),
                                       31, -0.5, 30.5, nPLB, pLBmin, pLBmax+1,
                                       "N_{hits}", "pLB", true, cstr( str("TrkNHitspLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));

      h_trk_nHoles_pLB[*nTrkCut][*ib] = HistogramHelper::define2DHistogram(cstr( str("Track nHoles for tracks associated with vertices vs. pLB, BCID",*ib), str(", NTrk >= ",*nTrkCut) ),
                                        31, -0.5, 30.5, nPLB, pLBmin, pLBmax+1,
                                        "N_{holes}", "pLB", true, cstr( str("TrkNHolespLB_BCID",*ib), str(", NTrk >= ",*nTrkCut) ));
      #endif

      // -- Maps <NTrkCut, <bcid, <pLB, count > > >
      for (Int_t current_pLB = pLBmin; current_pLB <= pLBmax; current_pLB++) {
        NVtxBeforeSplitCorrection[*nTrkCut][*ib][current_pLB] = 0;
        NVtxBeforeSplitCorrection_pertrack_perbcid[*nTrkCut][*ib]= 0;
        NVtxBeforeSplitCorrection_pertrack[*nTrkCut] = 0;
        NEvtBeforeSplitCorrection[*nTrkCut][*ib][current_pLB] = 0;
        NVtxAfterSplitCorrection[*nTrkCut][*ib][current_pLB] = 0;
        NEvtAfterSplitCorrection[*nTrkCut][*ib][current_pLB] = 0;

        if (*nTrkCut == 2) {
          //Only once, we need to initialize the number of triggers per pLB to zero.
          NTrig[*ib][current_pLB] = 0;
        }
      } // end pLB loop
    } // end bcid loop
  } // end ntrk loop
}

