#include "./VertexTruthMatchAlgorithm.h"
#include "InDetTruthVertexValidation/VertexMatchInfo.h"
#include "GaudiKernel/ITHistSvc.h"
#include "TH1D.h"

VertexTruthMatchAlgorithm::VertexTruthMatchAlgorithm(const std::string& name, ISvcLocator* pSvcLocator) :
  AthAlgorithm(name,pSvcLocator),
  m_VtxTmHandle("VertexTruthMatchAthAlgTool") 
{

  declareProperty( "VtxTm", m_VtxTmHandle, "public, shared IVertexTruthMatchAthAlgTool" );

}

VertexTruthMatchAlgorithm::~VertexTruthMatchAlgorithm() {}

StatusCode VertexTruthMatchAlgorithm::initialize() {
  msg(MSG::INFO) << "Initialising!" << endreq;
  //histo service
  if(service("THistSvc", m_tHistSvc).isFailure()) {
    msg(MSG::WARNING) << "Failed to get histogram service!" << endreq;
    return StatusCode::FAILURE;
  }


  // Get the VtxTm AthAlgTool
  if(m_VtxTmHandle.retrieve().isFailure()){
    msg(MSG::WARNING) << "Failed to get VertexTruthMatchAthAlgTool." << endreq;
    return StatusCode::FAILURE;
  }



  //setup the debug histos

    h_TM_goodUnmatched_pt = new TH1F("TMgoodUnmatchedPt", "P_{T} of matched tracks (prob > 0.7) with no gen particle saved.",
				     200,0,10.);
    h_TM_goodUnmatched_pt->GetXaxis()->SetTitle("p_{T} (GeV)");
    h_TM_goodUnmatched_pt->GetYaxis()->SetTitle("Entries / 0.5 GeV");



    h_TM_goodUnmatched_eta = new TH1F("TMgoodUnmatchedEta", "#eta of matched tracks (prob > 0.7) with no gen particle saved.",
				      60,-3.,3.);
    h_TM_goodUnmatched_eta->GetXaxis()->SetTitle("#eta");
    h_TM_goodUnmatched_eta->GetYaxis()->SetTitle("Entries / 0.1");

  

    h_TM_goodUnmatched_pteta = new TH2F("TMgoodUnmatchedPtEta", "P_{T}-#eta of matched tracks (prob > 0.7) with no gen particle saved.",
					200,0,10., 60, -3., 3.);
    h_TM_goodUnmatched_pteta->GetXaxis()->SetTitle("p_{T} (GeV)");
    h_TM_goodUnmatched_pteta->GetYaxis()->SetTitle("#eta");



    h_TM_highestW = new TH1F("TMhighestW", "Highest (relative) weight from a single genVertex associated to a reconstructed vertex. If fakes plot < 0",
			     100, -1., 1.);
    h_TM_highestW->GetXaxis()->SetTitle("Relative weight");
    h_TM_highestW->GetYaxis()->SetTitle("Entries");



    h_TM_2ndHighestW = new TH1F("TM2ndHighestW", "2nd Highest (relative) weight from a single genVertex associated to a reconstructed vertex. If fakes plot < 0",
				100, -1., 1.);
    h_TM_2ndHighestW->GetXaxis()->SetTitle("Relative weight");
    h_TM_2ndHighestW->GetYaxis()->SetTitle("Entries");

   

    h_TM_OtherHighestW = new TH1F("TMOtherHighestW", "All other (>2nd) relative weight from genVertices associated to a reconstructed vertex. If fakes plot < 0",
				  100, -1., 1.);
    h_TM_OtherHighestW->GetXaxis()->SetTitle("Relative weight");
    h_TM_OtherHighestW->GetYaxis()->SetTitle("Entries");



    h_TM_numGenMatched = new TH1F("TMnumGenMatched", "Number of genVertex matched (up to totatl threshold)",
				  15, -1.5, 13.5);
    h_TM_numGenMatched->GetXaxis()->SetTitle("Number of genVertex matched");
    h_TM_numGenMatched->GetYaxis()->SetTitle("Entries");



    h_TM_numGenMatchedAll = new TH1F("TMnumGenMatchedAll", "Number of genVertex matched (up to totatl threshold to Store)",
				     15, -1.5, 13.5);
    h_TM_numGenMatchedAll->GetXaxis()->SetTitle("Number of ALL genVertex matched");
    h_TM_numGenMatchedAll->GetYaxis()->SetTitle("Entries");



    h_TM_class = new TH1F("TMclass", "Truth-Match type for reconstructed vertices",
			  6, -0.5, 5.5);
    h_TM_class->GetXaxis()->SetTitle("Category");
    h_TM_class->GetYaxis()->SetTitle("Entries");
    h_TM_class->GetXaxis()->SetBinLabel(1, "Match");
    h_TM_class->GetXaxis()->SetBinLabel(2, "Merge");
    h_TM_class->GetXaxis()->SetBinLabel(3, "Split");
    h_TM_class->GetXaxis()->SetBinLabel(4, "Fake");
    h_TM_class->GetXaxis()->SetBinLabel(5, "Others");
    h_TM_class->GetXaxis()->SetBinLabel(6, "ERRORS");



    h_TM_fakeRelWeight = new TH1F("TMfakeRelWeight", "Relative fake component in reconstructed vertices",
				  100, 0, 1.);
    h_TM_fakeRelWeight->GetXaxis()->SetTitle("Relative FAKE Weight");
    h_TM_fakeRelWeight->GetYaxis()->SetTitle("Entries");



    h_TM_splitDeltaZ = new TH1F("TM_splitDeltaZ", "h_TM_splitDeltaZ", 100, -10, 10);
    h_TM_splitDeltaZ->GetXaxis()->SetTitle("#Delta Z between split vertices (w.r.t. Higher #sum P_{T}^2)");
    h_TM_splitDeltaZ->GetYaxis()->SetTitle("Entries");


    h_TM_splitDeltaZSig = new TH1F("TM_splitDeltaZSig", "h_TM_splitDeltaZSig", 100, -25, 25);
    h_TM_splitDeltaZSig->GetXaxis()->SetTitle("#Delta Z / #sigma(#Delta Z) between split vertices (w.r.t. Higher #sum P_{T}^2)");
    h_TM_splitDeltaZSig->GetYaxis()->SetTitle("Entries");




    //  m_VtxTmHandle->getMatcher()->SetDebug(1);
    //  m_VtxTmHandle->getMatcher()->InitDebugHistograms();

  //  TH1D * htype = new TH1D("h_matchType","Type of vertex match",6,-0.5,5.5);
  if (m_tHistSvc->regHist("/outfile/vtxtm/TMclass",h_TM_class).isFailure()) {
    msg(MSG::ERROR) << "Couldn't register matchType" << endmsg;
  }

  if (m_tHistSvc->regHist("/outfile/vtxtm/TMgoodUnmatchedpt",h_TM_goodUnmatched_pt).isFailure()) {
    msg(MSG::ERROR) << "Couldn't register matchType" << endmsg;
  }

  if (m_tHistSvc->regHist("/outfile/vtxtm/TMgoodUnmatchedeta",h_TM_goodUnmatched_eta).isFailure()) {
    msg(MSG::ERROR) << "Couldn't register matchType" << endmsg;
  }

  if (m_tHistSvc->regHist("/outfile/vtxtm/TMgoodUnmatchedpteta",h_TM_goodUnmatched_pteta).isFailure()) {
    msg(MSG::ERROR) << "Couldn't register matchType" << endmsg;
  }

  if (m_tHistSvc->regHist("/outfile/vtxtm/TMhighestW",h_TM_highestW).isFailure()) {
    msg(MSG::ERROR) << "Couldn't register matchType" << endmsg;
  }

  if (m_tHistSvc->regHist("/outfile/vtxtm/TM2ndHighestW",h_TM_2ndHighestW).isFailure()) {
    msg(MSG::ERROR) << "Couldn't register matchType" << endmsg;
  }

  if (m_tHistSvc->regHist("/outfile/vtxtm/TMotherHighestW",h_TM_OtherHighestW).isFailure()) {
    msg(MSG::ERROR) << "Couldn't register matchType" << endmsg;
  }

  if (m_tHistSvc->regHist("/outfile/vtxtm/TMnumGenMatched",h_TM_numGenMatched).isFailure()) {
    msg(MSG::ERROR) << "Couldn't register matchType" << endmsg;
  }

  if (m_tHistSvc->regHist("/outfile/vtxtm/TMnumGenMatchedAll",h_TM_numGenMatchedAll).isFailure()) {
    msg(MSG::ERROR) << "Couldn't register matchType" << endmsg;
  }


  return StatusCode::SUCCESS;
}

StatusCode VertexTruthMatchAlgorithm::execute() {
  msg(MSG::DEBUG) << "In execute!" << endreq;

  //Get vertex truth matching result - each VertexMatchInfo is for each reco vertex
  std::vector<VertexMatchInfo> matches = m_VtxTmHandle->getMatches();

  msg(MSG::DEBUG) << "Have " << matches.size() << " vertices with matching information." << endreq;
  //Loop over results and fill some debug histograms
  for(uint i=0; i< matches.size(); i++ ) {
    msg(MSG::DEBUG) << "Vertex " << i << " - reco index " << matches.at(i).m_recoVtx << ": is type " << matches.at(i).m_type << endreq;

    //type = matched, merged, split, etc
    h_TM_class->Fill(matches.at(i).GetType());

    //get number of vertices needed to pass threshold - default is 0.7
    //also save the highest weight, the 2nd highest weight, and other weights
    int numMatchedThreshold = matches.at(i).GetNMatched(0.7);
    if (numMatchedThreshold <= 0) {
      h_TM_numGenMatchedAll->Fill(numMatchedThreshold);
      h_TM_numGenMatched->Fill(numMatchedThreshold);
      h_TM_highestW->Fill(0);
    } else {
      h_TM_numGenMatchedAll->Fill(matches.at(i).m_matchList.size());
      h_TM_numGenMatched->Fill(numMatchedThreshold);
      if (matches.at(i).m_matchList[0].first == -1) {
	h_TM_highestW->Fill(-matches.at(i).m_matchList[0].second);
      } else {
	h_TM_highestW->Fill(matches.at(i).m_matchList[0].second);
      }
      //if (numMatchedThreshold > 1) { // -> plot only if merged or (future) split
      if (matches.at(i).m_matchList.size() > 1) { // -> plot whenever you have another significant contribution
	if (matches.at(i).m_matchList[1].first == -1) {
	  h_TM_2ndHighestW->Fill(-matches.at(i).m_matchList[1].second);
	} else {
	  h_TM_2ndHighestW->Fill(matches.at(i).m_matchList[1].second);
	}
      }
      if (matches.at(i).m_matchList.size() > 2) { // -> plot whenever you have another significant contribution
	float totWOthers=0.0;
	for (int ivg = 2; ivg < matches.at(i).m_matchList.size(); ivg++) 
	  totWOthers += matches.at(i).m_matchList[ivg].second;
	if (matches.at(i).m_matchList[2].first == -1) { //sign depends on the type of the first "other" (->3rd) vertex matched
	  h_TM_OtherHighestW->Fill(-totWOthers);
	} else {
	  h_TM_OtherHighestW->Fill(totWOthers);
	}
      }
    } // have at least one vertex matched and no weird errors

    //Fill the contribution from fakes
    int idxG=0;
    for (idxG=0; idxG < matches.at(i).m_matchList.size(); idxG++) 
      if (matches.at(i).m_matchList[idxG].first == -1) { //it's fakes
	h_TM_fakeRelWeight->Fill(matches.at(i).m_matchList[idxG].second);
	break;
      }
    if (idxG == matches.at(i).m_matchList.size())
      h_TM_fakeRelWeight->Fill(0); //fake component not found


    std::vector<VertexMatchInfo::GenVtxMatch> genvs = matches.at(i).m_matchList;
    for(uint j=0; j<genvs.size(); j++)
      msg(MSG::DEBUG) << "          Gen contribution from index " << genvs.at(j).first << " with fractional weight " << genvs.at(j).second << endreq;
  }






  /*
  //This code was studying the z separation for split vertices.  Would be able to code up something similar if you also load the vertex info in your algorithm

	//add to study histogram
	if (debugTM >= 1) {
	  double deltaz = vx_z->at(matchedVtx[candSplitIdx[sumpt2MaxIdx]].m_recoVtx);
	  deltaz = deltaz - vx_z->at(matchedVtx[candSplitIdx[idxS]].m_recoVtx);
	  double sigmaDz = pow(vx_cov_z->at(matchedVtx[candSplitIdx[sumpt2MaxIdx]].m_recoVtx), 2);
	  sigmaDz += pow(vx_cov_z->at(matchedVtx[candSplitIdx[idxS]].m_recoVtx), 2);
	  sigmaDz = sqrt(sigmaDz);
	  h_TM_splitDeltaZ->Fill(deltaz);
	  h_TM_splitDeltaZSig->Fill(deltaz / sigmaDz);
	}


  */




  return StatusCode::SUCCESS;
}

StatusCode VertexTruthMatchAlgorithm::finalize(){
  msg(MSG::INFO) << "Goodbye!" << endreq;

  return StatusCode::SUCCESS;
}
