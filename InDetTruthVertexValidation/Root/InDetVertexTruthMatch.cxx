#include "InDetTruthVertexValidation/InDetVertexTruthMatch.h"

#include <iostream>
#include <map>
#include <math.h>

#include "TDatabasePDG.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;

InDetVertexTruthMatch::InDetVertexTruthMatch()
{
  //no debug
  debugTM = 0;

  //init thresholds
  vtxTMWThreshold = 0.7;
  vtxTMWThresholdStore = 1 - (1 - vtxTMWThreshold) / 10.;
  vtxTMTrackMatchProb = 0.7; //matching probability for tracks

  //Gen-vertices "good" requirement
  genVertexRequirementVersion = 1;
  
  req_ngen = 0;
  req_genpt = 400.;
  req_geneta = 2.5;
  req_ntrks = 0;

  //Init PDG Table
  m_PDG = new TDatabasePDG();

  //Init cache
  cacheValid = false;
  m_runNumber = 0;
  m_evtNumber = -1;

  //Switch off debug plots by default
//   h_TM_goodUnmatched_pt = 0;
//   h_TM_goodUnmatched_eta = 0;  
//   h_TM_goodUnmatched_pteta = 0;
//   h_TM_highestW = 0;
//   h_TM_2ndHighestW = 0;
//   h_TM_OtherHighestW = 0;
//   h_TM_numGenMatched = 0;
//   h_TM_numGenMatchedAll = 0;
//   h_TM_class = 0;
//   h_TM_fakeRelWeight = 0;  
//   h_TM_splitDeltaZ=0; ///< Delta-Z between main and split vertices
//   h_TM_splitDeltaZSig=0; ///< Delta-Z/sigma(Delta-Z) between main and split vertices  
 
}

InDetVertexTruthMatch::~InDetVertexTruthMatch()
{
  delete m_PDG;
}

void InDetVertexTruthMatch::SetTrackMatchProbability(float trkMatchProb)
{
  vtxTMTrackMatchProb = trkMatchProb;
}

void InDetVertexTruthMatch::SetVtxMatchWeight(float vtxWTh)
{
  vtxTMWThreshold = vtxWTh;
}

void InDetVertexTruthMatch::SetVtxStoreWeight(float vtxWTh)
{
  vtxTMWThresholdStore = vtxWTh;
}

void InDetVertexTruthMatch::SetGenVertexRequirementVersion(Int_t version)
{
  genVertexRequirementVersion = version;
  if(genVertexRequirementVersion == 2)
    req_ngen=2;
  else if(genVertexRequirementVersion == 3)
    req_ngen=5;
}

void InDetVertexTruthMatch::SetDebug(int dbg)
{
  debugTM = dbg;
}

void InDetVertexTruthMatch::MatchVertices(int runNumber, long evtNumber)
{
  matchedVtx.clear(); //reset container

  //create cache if needed, then set cache is valid even if run/evt number are -1 so that in this method we will use cache 
  //at end of MatchVertices cacheValid is set to false
  cacheValid = false;
  isGoodGenVertex(0, runNumber, evtNumber);
  cacheValid = true;

  //  if (debugTM >= 1)
  //    InitDebugHistograms();

  if (debugTM >= 3) {
    cout << "Reconstructed vertices" << endl;
    for (int nv=0; nv < vx_n-1; nv++) {
      cout << nv << " Ntrk= " << vx_trk_n->at(nv) << endl;
    }
    if (debugTM >= 4) {
      cout << "Total n. reco tracks: " << trk_n << endl;
      cout << " Tracks (pt, probMC, MCidx)" << endl;
      for (int it=0; it < trk_n; it++)
        cout << it << " " << trk_pt->at(it) << ", " << trk_mc_probability->at(it) << ", " << trk_mc_index->at(it) << endl;
    }
    cout << "Generated vertices (PUType = 0,1 && At least 2 particles with pT > 400, |eta| < 2.4)" << endl;
    for (int nv=0; nv < mcvtx_n; nv++) {
      if (isGoodGenVertex(nv, runNumber, evtNumber))
	cout << nv << " GenEvt= " << mcvtx_mcevt_index->at(nv) << endl;
    }
    cout << "Generated events" << endl;
    for (int ne=0; ne < mcevt_n; ne++) {
      if (mcevt_pileUpType->at(ne) > 1)
	continue;
      Int_t npacc(0); //n particles within acceptance
      for (int idxPart=0;idxPart < mcpart_n; idxPart++) {      
	if (mcpart_mcevt_index->at(idxPart) == ne && mcpart_pt->at(idxPart) > 400. && fabs(mcpart_eta->at(idxPart)) < 2.4)
	  npacc++;
      }
      cout << ne << ". PU = " << mcevt_pileUpType->at(ne) << ", NParticles = " << mcevt_nparticle->at(ne) << ", NParticlesAcc = " << npacc << endl;
    }
    cout << "Number of tracks for each reconstructed vertex" << endl;
  } //debug

  //Loop over reconstructed vertices
  for (int nv=0; nv < vx_n-1; nv++) {
    Int_t nFakes(0); //not particle matched
    map<Int_t,Int_t> nVtxG; //matched to n-th generated vertex, with ntracks (n, ntracks)
    Int_t nSec(0); //not matched to a primary vertex -- not used now
    Float_t nFakesW(0);
    map<Int_t, Float_t> nVtxGW; // matched to n-th gen vtx, with weight W (n, W)
    Int_t nSecW(0);
    Float_t nTotW(0);
    for (int nvg=0; nvg < mcvtx_n; nvg++) {
      if (isGoodGenVertex(nvg, runNumber, evtNumber)) {
	nVtxG[nvg] = 0;
	nVtxGW[nvg]  = 0.0;
      }
    }
    //cout << "Ntrk = " << vx_trk_n->at(nv) << ". C1 = " << vx_trk_weight->at(nv).size() << ". C2 = " << vx_trk_index->at(nv).size() << endl;
    for (int ntrk=0; ntrk < vx_trk_n->at(nv); ntrk++) {
      nTotW += vx_trk_weight->at(nv)[ntrk];
      //first match track to truth particle
      if (vx_trk_index->at(nv)[ntrk] < 0 || vx_trk_index->at(nv)[ntrk] >= trk_n) {
	cerr << "ERROR in link trackAtVertex -> tracks." << endl;
	cout << "TrkAtVertex pointer: " << vx_trk_index->at(nv)[ntrk] << ", NTracks: " << trk_n << endl;
	break;
      }
      if ( trk_mc_probability->at(vx_trk_index->at(nv)[ntrk]) < vtxTMTrackMatchProb) {
	nFakes++;
	nFakesW += vx_trk_weight->at(nv)[ntrk];
	continue;
      } 
      //check to what gen vertex the particle is matched
      Int_t partIdx = trk_mc_index->at(vx_trk_index->at(nv)[ntrk]);
      if (partIdx < 0 || partIdx >= mcpart_n) {
	int tmpIdx = vx_trk_index->at(nv)[ntrk];
	//cerr << "ERROR: Track with no parent particle and matched prob >= vtxTMTrackMatchProb? Idx=" << partIdx 
	//     << ". Prob = " << trk_mc_probability->at(tmpIdx) << endl;
	//cerr << "Pt, eta, phi = " << trk_pt->at(tmpIdx) << ", " << trk_eta->at(tmpIdx) << ", " << trk_phi_wrtPV->at(tmpIdx) 
	//     << ". d0, z0 = " << trk_d0_wrtPV->at(tmpIdx) << ", " << trk_z0_wrtPV->at(tmpIdx) << endl;	
	//count them as fakes :(
	nFakes++;
	nFakesW += vx_trk_weight->at(nv)[ntrk];
	continue;
      }

      // Use event Id to match vertex
      Int_t genEvtIdx = mcpart_mcevt_index->at(partIdx);
      if (mcevt_pileUpType->at(genEvtIdx) != 0 && mcevt_pileUpType->at(genEvtIdx) != 1) {
	cout << "WARNING: Fake from other Vtx type: GenEvt = " << genEvtIdx << endl;
	nFakes++;
	nFakesW += vx_trk_weight->at(nv)[ntrk];
	continue;
      }
     
      Int_t genVxIdx = -1;
      for (int nvg=0; nvg < mcvtx_n; nvg++) {
	if (mcvtx_mcevt_index->at(nvg) == genEvtIdx) {
	  genVxIdx = nvg;
	  break;
	}
      }
      if (genVxIdx < 0) {
	cerr << "ERROR: Could not find match genVtx for given genEvt with reconstructed tracks. genEvt = " << genEvtIdx 
	     << ", reco vtx = " << nv << ", track pt = " << trk_pt->at(ntrk) << endl;
	break;
      }

      nVtxG[genVxIdx]++;
      nVtxGW[genVxIdx] += vx_trk_weight->at(nv)[ntrk];
    }
    //now store track weights for each generated vertex up to a given threshold
    VertexMatchInfo tmpVtxMatch;
    tmpVtxMatch.m_totWeight = nTotW;
    map<Int_t, Float_t> tmpVtxGW = nVtxGW;
    map<Int_t,Int_t> tmpVtxG = nVtxG;
    Float_t tmpFakesW=nFakesW;
    Float_t totWStored=0.0;
    while ( totWStored < vtxTMWThresholdStore ) {
      // get maximum
      map<Int_t, Float_t>::iterator itWMax = tmpVtxGW.begin();
      map<Int_t, Float_t>::iterator itW;
      map<Int_t, Int_t>::iterator itNT = tmpVtxG.begin();
      map<Int_t, Int_t>::iterator itNTMax = tmpVtxG.begin();

      for (itW = tmpVtxGW.begin(); itW != tmpVtxGW.end(); ++itW, ++itNT) {
	if (itW->second > itWMax->second) {
	  itWMax = itW;
	  itNTMax = itNT;
	}
      }      
      if (tmpVtxGW.size() == 0 && tmpFakesW <= 0) {
	//no more vertices
	tmpVtxMatch.Add(-1, -1, 0);
      }
      // store it (if less than fake contribution
      if (itWMax->second >= tmpFakesW) {
	if (itWMax == tmpVtxGW.end()) {
	  //no maximum? all 0? Add an error entry
	  tmpVtxMatch.Add(-1, -1, 0);
	} else {
	  tmpVtxMatch.Add(itWMax->first, itWMax->second / nTotW, itNTMax->second); //store index and relative weight
	}
	totWStored += itWMax->second / nTotW;
	// remove from initial collection
	tmpVtxGW.erase(itWMax);
	tmpVtxG.erase( itNTMax );
      } else {
	//fake contribution dominant
	tmpVtxMatch.Add(-1, tmpFakesW / nTotW, nFakes);
	totWStored += tmpFakesW / nTotW;
	tmpFakesW = 0.0; //prevent double-counting :)	
      }
    }

    //finally make classification for matched, merged, fakes (at this stage split == merge)
    if (tmpVtxMatch.GetNMatched(vtxTMWThreshold) <= 0) {
      tmpVtxMatch.m_type = VertexMatchInfo::VtxTM_NMatch; //ERROR
    } else if (tmpVtxMatch.m_matchList[0].first == -1) {
      //most of the contribution is fake -> mark vertex as fake
      tmpVtxMatch.m_type = VertexMatchInfo::VtxTM_Fake;
    } else if (tmpVtxMatch.m_matchList[0].second > vtxTMWThreshold) {
      //if matched 1-1 (first weight > threshold) -> matched
      tmpVtxMatch.m_type = VertexMatchInfo::VtxTM_Match;
    } else {
      // more than one generated vertex contributing
      tmpVtxMatch.m_type = VertexMatchInfo::VtxTM_Merge;

      //OLD DEFINITION - it would be misleading in case of running on truth slimmed D3PD
//       if (tmpVtxMatch.m_matchList[1].first == -1) 
// 	tmpVtxMatch.m_type = VertexMatchInfo::VtxTM_Match; //second highest weight is from fakes -> call it matched anyway
//       else 
// 	tmpVtxMatch.m_type = VertexMatchInfo::VtxTM_Merge; //second highest weight is a true vertex -> call it merged
    }
    
    //Add to matched list
    tmpVtxMatch.m_recoVtx = nv;
    matchedVtx.push_back(tmpVtxMatch);



    if (debugTM >= 3) {

      cout << nv << "_Weights. Tot: " << nTotW
	   << ", Fakes: " << nFakesW / nTotW
	//	 << ", Secondaries: " << nSecW // not really?
	   << ", GenVtx: ";
      for (int nvg=0; nvg < mcvtx_n; nvg++) {
	if (isGoodGenVertex(nvg, runNumber, evtNumber))
	  cout << nVtxGW[nvg] / nTotW << ", ";
      }
      cout << endl;    
      if (debugTM >= 5) {
	int numMatchedThreshold=tmpVtxMatch.GetNMatched(vtxTMWThreshold);

	cout << "TM Class (before split): " << tmpVtxMatch.m_type;
	cout << ", numGenMatch: " << numMatchedThreshold;
	cout << ", numGenMatchAll: " << tmpVtxMatch.m_matchList.size();
	cout << endl;      
      }
    } //debug
    
  } //end loop over reconstructed vertices

  //now that we have all, check if there are splits
  //A Split is defined is the same generated vertex is shared *as dominant contributor* among several reconstructed vertices.
  //If that happens, the class of the vertex with highest SumPt2 is preserved (either matched or merged) the others are called Split
  for (int ng=0; ng < mcvtx_n; ng++) {

    vector<Int_t> candSplitIdx; //store index in matchedVtx of candidate split vertices
    //loop over reconstructed vertices and find the ones which share generated vertex ng as the biggest contributor
    for (Int_t idxRM  = 0; idxRM < matchedVtx.size(); idxRM++) {
      if (matchedVtx[idxRM].m_matchList.size() > 0) 
	if (matchedVtx[idxRM].m_matchList[0].first == ng) {

	  candSplitIdx.push_back(idxRM);
	}
    }

    if (candSplitIdx.size() > 1) {
      //more than one reco vertex sharing the same generated one as main contributor
      //get index of the one with max sumpt2
      float sumpt2Max=0.0;
      Int_t sumpt2MaxIdx=-1;
      for (int idxS=0; idxS < candSplitIdx.size(); idxS++) {
	    Double_t sumPt2=0.0;
	    Int_t recoIdx = matchedVtx[candSplitIdx[idxS]].m_recoVtx;

	    for (unsigned itTrk=0; itTrk < (*vx_trk_n)[recoIdx]; itTrk++)
	      sumPt2 += pow((*trk_pt)[(*vx_trk_index)[recoIdx][itTrk]],2);
	    sumPt2 = sqrt(sumPt2);
	    if (sumPt2 > sumpt2Max) {
	      sumpt2Max = sumPt2;
	      sumpt2MaxIdx = idxS;
	    }
      }
      //now set split categories for all but the highest sumpt2 vertex
      for (int idxS=0; idxS < candSplitIdx.size(); idxS++) {
	if (idxS == sumpt2MaxIdx)
	  continue;
	matchedVtx[candSplitIdx[idxS]].m_type = VertexMatchInfo::VtxTM_Split;

      }
    }
  } 
  
  //debug
  if (debugTM>=2) {
    //print matchedVtx
    cout << "Truth-Matching results" << endl;
    for (int im=0; im < matchedVtx.size(); im++) {
      cout << "Reco Vtx = " << matchedVtx[im].m_recoVtx 
	   << ". TYPE = " << matchedVtx[im].m_type 
	   << ". GenVertices: ";
      for (int img=0; img < matchedVtx[im].m_matchList.size(); img++)
	cout << "(" << matchedVtx[im].m_matchList[img].first << ", " << matchedVtx[im].m_matchList[img].second 
	     << ") ";
      cout << endl;
    }
  }

  //turn off cache override - cache will only be used when setup w/ real run and event number
  cacheValid = false;


}

int InDetVertexTruthMatch::getMotherOfSplit(int idxRecoVtx)
{
  //first getn info on requested vertex
  vector<VertexMatchInfo>::iterator itReqVtx;
  for (itReqVtx = matchedVtx.begin();
       itReqVtx != matchedVtx.end(); ++itReqVtx) {
    if (itReqVtx->m_recoVtx == idxRecoVtx)
      break;
  }
  if (itReqVtx == matchedVtx.end())
    return -1; //requested vertex not found in matched list
  if (itReqVtx->GetType() != VertexMatchInfo::VtxTM_Split)
    return -2; //requested vertex is not a split!

  //Check which is the mother-split vertex
  for (std::vector<VertexMatchInfo>::iterator itRM = matchedVtx.begin();
       itRM != matchedVtx.end(); ++itRM) {
    //check dominant contribution
    if (itRM->GetType() == VertexMatchInfo::VtxTM_Split ||
	itRM->GetType() == VertexMatchInfo::VtxTM_Fake)
      continue; //no splits, no fakes
    if (itRM->m_matchList.size() == 0)
      continue; //we need matched, this should almost never happen (do you have other vertex categories not excluded above?)
    if (itRM->m_matchList[0].first == itReqVtx->m_matchList[0].first) {
      //found it!
      return itRM->m_recoVtx;
    }
  } // loop over vertices
  return -3; //mother vertex not found
}

int InDetVertexTruthMatch::isGoodGenVertex(int vertexToCheck, int runNumber, long evtNumber)
{
  //From profiling, it seems we spend 40% of our time here... crazy.
  //Implement a caching system if a runNumber/evtNumber is provided 

  //when true we are overriding the cache
  if(cacheValid == true)
    return m_isGenVtxGood[vertexToCheck];

  if (runNumber == -1 || evtNumber == -1 || not (runNumber == m_runNumber && evtNumber == m_evtNumber)) {
    //info not cached. calculate isGenVtxGood for all the vertices
    m_isGenVtxGood.clear();
    m_runNumber = runNumber; 
    m_evtNumber = evtNumber;
    
//     Int_t minNGoodTracks=2; //min number of good (pt, eta) charged tracks at the gen vertex
//     float minGenPt=400.; //min pT of considered tracks
//     float maxGenEta=2.5; //max eta of conisdered tracks
//     if (genVertexRequirementVersion == 3) {
//       minNGoodTracks = 5;
//     }
    
    //First require it to be either a PV or in-time PileUp vertex
    for (int vtxToCheck = 0; vtxToCheck < mcvtx_n; vtxToCheck++) {
      bool isGood=false;
      Int_t genEventIdx = mcvtx_mcevt_index->at(vtxToCheck);      
      if(mcevt_pileUpType->at(genEventIdx) == 0 || 
	 mcevt_pileUpType->at(genEventIdx) == 1 ){	  //always require in time event

	if(req_ngen == 0) //dont need to check particles
	  isGood=true;
	else {

	  int nacc=0;
	  for (int idxPart=0;idxPart < mcpart_n; idxPart++) {   //count number passing acceptance

	    if (mcpart_mcevt_index->at(idxPart) != genEventIdx || mcpart_status->at(idxPart) != 1 || mcpart_barcode->at(idxPart) >= 200000) continue;

	    const TParticlePDG *pdgPart = m_PDG->GetParticle(mcpart_type->at(idxPart));
	    double charge = ( pdgPart ? pdgPart->Charge() / 3.0 : 0 );

	    if (fabs(charge) > 0.01 &&
		mcpart_pt->at(idxPart) > req_genpt && fabs(mcpart_eta->at(idxPart)) < req_geneta)
	      nacc++;

	    if( nacc >= req_ngen ) {
	      isGood=true;
	      break; //can stop looking now
	    }

	  }


	} //end else req_ngen == 0


      } //end if in time


      m_isGenVtxGood.push_back(isGood);

    } // loop over vertices
  } // refresh cache
  return m_isGenVtxGood[vertexToCheck];
}

// void InDetVertexTruthMatch::InitDebugHistograms(string prefixName)
// {
//   if (! h_TM_goodUnmatched_pt) {
//     h_TM_goodUnmatched_pt = new TH1F((string("TM")+prefixName+string("goodUnmatchedPt")).c_str(), "P_{T} of matched tracks (prob > 0.7) with no gen particle saved.",
// 				     200,0,10.);
//     h_TM_goodUnmatched_pt->GetXaxis()->SetTitle("p_{T} (GeV)");
//     h_TM_goodUnmatched_pt->GetYaxis()->SetTitle("Entries / 0.5 GeV");
//   }

//   if (!h_TM_goodUnmatched_eta) {
//     h_TM_goodUnmatched_eta = new TH1F((string("TM")+prefixName+string("goodUnmatchedEta")).c_str(), "#eta of matched tracks (prob > 0.7) with no gen particle saved.",
// 				      60,-3.,3.);
//     h_TM_goodUnmatched_eta->GetXaxis()->SetTitle("#eta");
//     h_TM_goodUnmatched_eta->GetYaxis()->SetTitle("Entries / 0.1");
//   }
  
//   if (!h_TM_goodUnmatched_pteta) {
//     h_TM_goodUnmatched_pteta = new TH2F((string("TM")+prefixName+string("goodUnmatchedPtEta")).c_str(), "P_{T}-#eta of matched tracks (prob > 0.7) with no gen particle saved.",
// 					200,0,10., 60, -3., 3.);
//     h_TM_goodUnmatched_pteta->GetXaxis()->SetTitle("p_{T} (GeV)");
//     h_TM_goodUnmatched_pteta->GetYaxis()->SetTitle("#eta");
//   }

//   if (!h_TM_highestW) {
//     h_TM_highestW = new TH1F((string("TM")+prefixName+string("highestW")).c_str(), "Highest (relative) weight from a single genVertex associated to a reconstructed vertex. If fakes plot < 0",
// 			     100, -1., 1.);
//     h_TM_highestW->GetXaxis()->SetTitle("Relative weight");
//     h_TM_highestW->GetYaxis()->SetTitle("Entries");
//   }

//   if (!h_TM_2ndHighestW) {
//     h_TM_2ndHighestW = new TH1F((string("TM")+prefixName+string("2ndHighestW")).c_str(), "2nd Highest (relative) weight from a single genVertex associated to a reconstructed vertex. If fakes plot < 0",
// 				100, -1., 1.);
//     h_TM_2ndHighestW->GetXaxis()->SetTitle("Relative weight");
//     h_TM_2ndHighestW->GetYaxis()->SetTitle("Entries");
//   }
   
//   if (!h_TM_OtherHighestW) {
//     h_TM_OtherHighestW = new TH1F((string("TM")+prefixName+string("OtherHighestW")).c_str(), "All other (>2nd) relative weight from genVertices associated to a reconstructed vertex. If fakes plot < 0",
// 				  100, -1., 1.);
//     h_TM_OtherHighestW->GetXaxis()->SetTitle("Relative weight");
//     h_TM_OtherHighestW->GetYaxis()->SetTitle("Entries");
//   }

//   if (!h_TM_numGenMatched) {
//     h_TM_numGenMatched = new TH1F((string("TM")+prefixName+string("numGenMatched")).c_str(), "Number of genVertex matched (up to totatl threshold)",
// 				  15, -1.5, 13.5);
//     h_TM_numGenMatched->GetXaxis()->SetTitle("Number of genVertex matched");
//     h_TM_numGenMatched->GetYaxis()->SetTitle("Entries");
//   }

//   if (!h_TM_numGenMatchedAll) {
//     h_TM_numGenMatchedAll = new TH1F((string("TM")+prefixName+string("numGenMatchedAll")).c_str(), "Number of genVertex matched (up to totatl threshold to Store)",
// 				     15, -1.5, 13.5);
//     h_TM_numGenMatchedAll->GetXaxis()->SetTitle("Number of ALL genVertex matched");
//     h_TM_numGenMatchedAll->GetYaxis()->SetTitle("Entries");
//   }

//   if (! h_TM_class) {
//     h_TM_class = new TH1F((string("TM")+prefixName+string("class")).c_str(), "Truth-Match type for reconstructed vertices",
// 			  6, -0.5, 5.5);
//     h_TM_class->GetXaxis()->SetTitle("Category");
//     h_TM_class->GetYaxis()->SetTitle("Entries");
//     h_TM_class->GetXaxis()->SetBinLabel(1, "Match");
//     h_TM_class->GetXaxis()->SetBinLabel(2, "Merge");
//     h_TM_class->GetXaxis()->SetBinLabel(3, "Split");
//     h_TM_class->GetXaxis()->SetBinLabel(4, "Fake");
//     h_TM_class->GetXaxis()->SetBinLabel(5, "Others");
//     h_TM_class->GetXaxis()->SetBinLabel(6, "ERRORS");
//   }

//   if (!h_TM_fakeRelWeight) {
//     h_TM_fakeRelWeight = new TH1F((string("TM")+prefixName+string("fakeRelWeight")).c_str(), "Relative fake component in reconstructed vertices",
// 				  100, 0, 1.);
//     h_TM_fakeRelWeight->GetXaxis()->SetTitle("Relative FAKE Weight");
//     h_TM_fakeRelWeight->GetYaxis()->SetTitle("Entries");
//   }

//   if (!h_TM_splitDeltaZ) {
//     h_TM_splitDeltaZ = new TH1F((string("TM")+prefixName+string("_splitDeltaZ")).c_str(), "h_TM_splitDeltaZ", 100, -10, 10);
//     h_TM_splitDeltaZ->GetXaxis()->SetTitle("#Delta Z between split vertices (w.r.t. Higher #sum P_{T}^2)");
//     h_TM_splitDeltaZ->GetYaxis()->SetTitle("Entries");
//   }

//   if (!h_TM_splitDeltaZSig) {
//     h_TM_splitDeltaZSig = new TH1F((string("TM")+prefixName+string("_splitDeltaZSig")).c_str(), "h_TM_splitDeltaZSig", 100, -25, 25);
//     h_TM_splitDeltaZSig->GetXaxis()->SetTitle("#Delta Z / #sigma(#Delta Z) between split vertices (w.r.t. Higher #sum P_{T}^2)");
//     h_TM_splitDeltaZSig->GetYaxis()->SetTitle("Entries");
//   }

// }

void InDetVertexTruthMatch::SetRecoVtxInfo(int p_vx_n, 
			      std::vector<int> *p_vx_trk_n, std::vector<std::vector<float> > *p_vx_trk_weight, std::vector<std::vector<int> > *p_vx_trk_index)//,
// 		    std::vector<float> *p_vx_x, std::vector<float> *p_vx_y, std::vector<float> *p_vx_z,
// 		    std::vector<float> *p_vx_cov_x, std::vector<float> *p_vx_cov_y, std::vector<float> *p_vx_cov_z)
{
  vx_n = p_vx_n;
//   vx_x = p_vx_x;
//   vx_y = p_vx_y;
//   vx_z = p_vx_z;
//   vx_cov_x = p_vx_cov_x;
//   vx_cov_y = p_vx_cov_y;
//   vx_cov_z = p_vx_cov_z;
  vx_trk_n = p_vx_trk_n;
  vx_trk_weight = p_vx_trk_weight;
  vx_trk_index = p_vx_trk_index;
}

///Reconstructed traks inputs (D3PD)
void InDetVertexTruthMatch::SetRecoTrkInfo(int p_trk_n,
		    vector<float> *p_trk_mc_probability, vector<int> *p_trk_mc_index,
			      vector<float> *p_trk_pt)//, vector<float> *p_trk_eta, vector<float> *p_trk_d0_wrtPV, vector<float> *p_trk_z0_wrtPV)
{
  trk_n = p_trk_n;
  trk_pt = p_trk_pt;
//   trk_eta = p_trk_eta;
//   trk_d0_wrtPV = p_trk_d0_wrtPV;
//   trk_z0_wrtPV = p_trk_z0_wrtPV;
  trk_mc_probability = p_trk_mc_probability;
  trk_mc_index = p_trk_mc_index;
}

///Generated particles inputs (D3PD)
void InDetVertexTruthMatch::SetGenPartInfo(int p_mcpart_n,
			      vector<int> *p_mcpart_mcevt_index,// vector<int> *p_mcpart_mcprodvtx_index,
			      vector<int> *p_mcpart_type,
			      vector<float> *p_mcpart_pt, vector<float> *p_mcpart_eta, vector<int> *p_mcpart_barcode, vector<int> *p_mcpart_status)
{
  mcpart_n = p_mcpart_n;
  mcpart_pt = p_mcpart_pt;
  mcpart_eta = p_mcpart_eta;
  mcpart_mcevt_index = p_mcpart_mcevt_index;
  //  mcpart_mcprodvtx_index = p_mcpart_mcprodvtx_index;
  mcpart_type = p_mcpart_type;
  mcpart_barcode = p_mcpart_barcode;
  mcpart_status = p_mcpart_status;
}

///Generated vertices inputs (D3PD)
void InDetVertexTruthMatch::SetGenVtxInfo(int p_mcvtx_n,
			     vector<int> *p_mcvtx_mcevt_index)//,
//		   vector<float> *p_mcvtx_x, vector<float> *p_mcvtx_y, vector<float> *p_mcvtx_z)
{
  mcvtx_n = p_mcvtx_n;
//   mcvtx_x = p_mcvtx_x;
//   mcvtx_y = p_mcvtx_y;
//   mcvtx_z = p_mcvtx_z;
  mcvtx_mcevt_index = p_mcvtx_mcevt_index;  
}

///Generated events inputs (D3PD)
void InDetVertexTruthMatch::SetGenEventsInfo(int p_mcevt_n, vector<short> *p_mcevt_pileUpType, vector<int> *p_mcevt_nparticle)
{
  mcevt_n = p_mcevt_n;
  mcevt_pileUpType = p_mcevt_pileUpType;
  mcevt_nparticle = p_mcevt_nparticle;
}

// void InDetVertexTruthMatch::WriteDebugHisto()
// {
//   h_TM_goodUnmatched_pt->Write();
//   h_TM_goodUnmatched_eta->Write();  
//   h_TM_goodUnmatched_pteta->Write();
//   h_TM_highestW->Write();
//   h_TM_2ndHighestW->Write();
//   h_TM_OtherHighestW->Write();
//   h_TM_numGenMatched->Write();
//   h_TM_numGenMatchedAll->Write();
//   h_TM_class->Write();
//   h_TM_fakeRelWeight->Write();  
//   h_TM_splitDeltaZ->Write(); 
//   h_TM_splitDeltaZSig->Write();
// }

VertexMatchInfo InDetVertexTruthMatch::GetRecoVtxInfo(int recoIdx)
{
  //first guess (should be the same ordering)
  if (matchedVtx[recoIdx].m_recoVtx == recoIdx)
    return matchedVtx[recoIdx];

  //loop over all values stored
  for (vector<VertexMatchInfo>::iterator iVtx = matchedVtx.begin();
       iVtx != matchedVtx.end(); ++iVtx) {
    if (iVtx->m_recoVtx == recoIdx)
      return *iVtx;
  }

  //not found. Return empty
  return VertexMatchInfo();
}

VertexMatchInfo InDetVertexTruthMatch::GetHardScatterBestMatch()
{

  //look for vertex w/ primary contribution from hard scatter
  for (vector<VertexMatchInfo>::iterator iVtx = matchedVtx.begin();
       iVtx != matchedVtx.end(); ++iVtx) {

    if(iVtx->m_type <2 && iVtx->m_matchList[0].first == 0 ) //hard scatter is first and is matched or merged (this will make sure its the best if there are multiple contributions)
      return *iVtx;

  }

  //now need to check for subdominant?

  //else no match so return empty
  return VertexMatchInfo();

}


std::vector<VertexMatchInfo> InDetVertexTruthMatch::GetHardScatterInfo()
{
  //get truth match info for all vertices with contribution from hard scatter
  std::vector<VertexMatchInfo> res;
  std::vector<double> weights;
  for (vector<VertexMatchInfo>::iterator iVtx = matchedVtx.begin();
       iVtx != matchedVtx.end(); ++iVtx) {
    for(unsigned int i=0; i< iVtx->m_matchList.size(); i++){
      if(iVtx->m_matchList[i].first == 0 ) {
	//found a vertex with truth info - now insert into proper place...
	if(res.size()==0 || iVtx->m_matchList[i].second <= weights.back()) {//either vector is empty or new addition has lowest (or equal lowest) weight
	  res.push_back(*iVtx);
	  weights.push_back(iVtx->m_matchList[i].second);
	} else if(iVtx->m_matchList[i].second >= weights.front() ) { //new addition has highest (or equal highest) weight
	  res.insert(res.begin(), *iVtx);
	  weights.insert(weights.begin(), iVtx->m_matchList[i].second);
	} else { //now need to find where it goes

	  std::vector<double>::iterator it = weights.begin();
	  std::vector<VertexMatchInfo>::iterator it2 = res.begin();
	  //shouldnt get to the end of the vector with this loop bc of first if above...
	  while( it!= weights.end() && iVtx->m_matchList[i].second < *it ) { //will always pass in once cause failed else if above
	    it++;
	    it2++; //move iterators together
	  }

	  if(it == weights.end() )
	    std::cout << "Error with HS sorting loop" << std::endl;
	  else{
	    res.insert(it2, *iVtx);
	    weights.insert( it, iVtx->m_matchList[i].second);
	  }

	}
	break; //dont need to check rest of match list since already found a HS contribution
      }
    }
  }

  //debug
  // if( res.size() > 2 ) {
  //   std::cout << "Index   Weight" << std::endl;
  //   for(unsigned int i=0; i< res.size(); i++) {
  //     std::cout << res[i].m_recoVtx << "      " << weights[i] << std::endl;
  //   }

  // }

  return res;
}
