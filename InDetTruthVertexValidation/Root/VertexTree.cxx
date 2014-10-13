#define VertexTree_cxx
#include "InDetTruthVertexValidation/VertexTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "InDetTruthVertexValidation/InDetVertexTruthMatch.h"

void VertexTree::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L VertexTree.C
//      Root > VertexTree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);
   fChain->SetBranchStatus("ei_*",1);
   fChain->SetBranchStatus("trk_*",1);
   fChain->SetBranchStatus("vx_*",1);
   fChain->SetBranchStatus("mc*",1);

   Long64_t nentries = fChain->GetEntriesFast();

   TFile fout(fOutName.c_str(),"recreate");
   InDetVertexTruthMatch m_vtxtm;


  //debug histos
  TH1F * h_TM_highestW;
  TH1F * h_TM_2ndHighestW;
  TH1F * h_TM_OtherHighestW;
  TH1F * h_TM_numGenMatched;
  TH1F * h_TM_numGenMatchedAll;
  TH1F * h_TM_class;
  TH1F * h_TM_fakeRelWeight;  

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



   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

 

      //Setuup truth matcher
      
      m_vtxtm.SetRecoVtxInfo(vx_n, vx_trk_n, vx_trk_weight, vx_trk_index);

      m_vtxtm.SetRecoTrkInfo(trk_n, trk_mc_probability, trk_mc_index, trk_pt);

      m_vtxtm.SetGenPartInfo(mcpart_n, mcpart_mcevt_index, mcpart_type,
			     mcpart_pt, mcpart_eta, mcpart_barcode, mcpart_status);

      m_vtxtm.SetGenVtxInfo(mcvtx_n,mcvtx_mcevt_index);

      m_vtxtm.SetGenEventsInfo(mcevt_n, mcevt_pileUpType, mcevt_nparticle);
      
      //run matching
      m_vtxtm.MatchVertices(ei_RunNumber, ei_EventNumber);

      //get info
      std::vector<VertexMatchInfo> matches = m_vtxtm.matchedVtx;
      cout << "Have " << matches.size() << " vertices with matching information." << endl;
      //Loop over results and fill some debug histograms
      for(unsigned int i=0; i< matches.size(); i++ ) {
	cout << "Vertex " << i << " - reco index " << matches.at(i).m_recoVtx << ": is type " << matches.at(i).m_type << endl;
	
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
	for(unsigned int j=0; j<genvs.size(); j++)
	  cout << "          Gen contribution from index " << genvs.at(j).first << " with fractional weight " << genvs.at(j).second << endl;
      }



   }

   fout.cd();
   h_TM_highestW->Write();
   h_TM_2ndHighestW->Write();
   h_TM_OtherHighestW->Write();
   h_TM_numGenMatched->Write();
   h_TM_numGenMatchedAll->Write();
   h_TM_class->Write();
   h_TM_fakeRelWeight->Write();  

   fout.Close();



}
