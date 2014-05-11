#define DEBUG
/**
  *	Script to calculate deadtime per pseudolumiblock.
  * Author: David R. Yu (dryu@lbl.gov).
  * Credit to David Berge for writing the initial python script.
  *		(this script is mostly for performance gains)
  *
  *	Input: 
  * 	- TRP trigger rates root file
  *		- TAG ntuple for event rates.
  *	Output:
  *		- Live fractions vs. pLB histogram
  *		- Saves to livefractions_<run>.root.
  *
  *	Note: segfaults in ROOT 5.32! 5.30 works. No idea why.
  */
  
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <stdlib.h>

#include "TROOT.h"
#include "TString.h"
#include "TH1D.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TLeafI.h"
#include "TLeafD.h"
#include "TMath.h"
#include "TObject.h"
#include "TChain.h"

using namespace std;

Int_t GetPlbFromTimestamp( Double_t timestamp );
Int_t GetPlbFromTimestampTrp( Double_t timestamp );
void LoadPlbTimestamps(TString path);
void LoadPlbTimestampsTrp(TString path);

std::map<Int_t, std::pair<Double_t, Double_t> > pLB_timestamps;
std::map<Int_t, std::pair<Double_t, Double_t> > pLB_timestamps_trp;
std::vector<Int_t> pLB_list;
std::vector<Int_t> pLB_list_trp;
std::map<Int_t, Float_t> lb_duration;
std::map<Int_t, Float_t> lb_duration_trp;

int main(int argc, char **argv) {
	cout << "Starting deadtime calculation" << endl;
    TString run;
    // --- Scan for command line parameters
    int c;
    extern int optind;
    extern char* optarg;
    while (1) {
		int option_index = 0;
		static struct option long_options[] = { {0,0,0,0} };
		c=getopt_long(argc, argv, "r:",long_options,&option_index);
		if (c == -1) break;
		switch(c) {
			case 'r':
				{
				run = optarg;
				break;
				}
		}
	}	
	
	std::vector<TString> ts_paths;
	TString trp_path;
	TString tag_path;
	TString D3PD_folder;
	TString out_path;
	if (run == "191373") {
	  tag_path = "/eliza18/atlas/dryu/Luminosity/ALFA/data11_7TeV.00191373.calibration_VdM.merge.TAG.x176_m841_m840/merge.root";
	  D3PD_folder = "/eliza18/atlas/dryu/Luminosity/ALFA/user.spagan.data11_7TeV.00191373.calibration_VdM.merge.VTX_MON.x176_m841.v1.0";
	  trp_path = "/eliza18/atlas/dryu/Luminosity/ALFA/deadtime/TriggerRates_ATLAS_191373.root";
	  ts_paths.push_back("/eliza18/atlas/dryu/Luminosity/ALFA/deadtime/191373_timestamps.dat");
	  out_path = "dt_191373.root";
	} else {
		cerr << "Specified run number not valid." << endl;
		exit(1);
	}
	
/***************************************************************************/
// Timestamps: take only the intersection of pLBs and TRP intervals.
/***************************************************************************/
	TFile *f_trp = new TFile(trp_path, "READ");
	TTree *t_trp = (TTree*)f_trp->Get("L1_Rate");
	Long64_t n_entries = t_trp->GetEntriesFast();

#define REDOTRPTIMESTAMPS
#ifdef REDOTRPTIMESTAMPS
	
	//Load pLB timestamps
	cout << "Loading pLB timestamps" << endl;
	for (vector<TString>::iterator it = ts_paths.begin(); it != ts_paths.end(); ++it) {
		LoadPlbTimestamps(*it);
	}
	
	//Sort pLB timestamps, and find max/min
	sort(pLB_list.begin(), pLB_list.end());
	#ifdef DEBUG
	cout << "Read in pLB timestamps: range " << *(pLB_list.begin()) << " - " << *(pLB_list.end()-1) << endl;
	#endif
	
	// Find TRP boundaries inside pLBs.
	cout << "Finding new pLB timestamps based on TRP boundaries." << endl;


	t_trp->SetBranchStatus("*", 0);
	t_trp->SetBranchStatus("TimeStamp", 1);
	Double_t current_ts;	
	t_trp->SetBranchAddress("TimeStamp", &current_ts);	
	
	//Make entire list of TRP timestamps
	vector<Double_t> trp_timestamps;
	for (Int_t i = 0; i < n_entries; i++) {
		t_trp->GetEntry(i);
		trp_timestamps.push_back(current_ts);
	}
	cout << "Sorting TRP timestamps, be patient." << endl;
	sort(trp_timestamps.begin(), trp_timestamps.end());
	
	//Keep track of whether we've found a pLB or not
	std::map<Int_t, bool> found_pLB;
	for (vector<Int_t>::iterator it = pLB_list.begin(); it != pLB_list.end(); ++it) {
		found_pLB[*it] = false;
	}
	
	// Find consecutive timestamp pairs both in one pLB.
	for (vector<Double_t>::iterator ts_end = trp_timestamps.begin()+1; ts_end != trp_timestamps.end(); ++ts_end) {
		vector<Double_t>::iterator ts_start = ts_end - 1;
		Int_t pLB_start = GetPlbFromTimestamp(*ts_start);
		Int_t pLB_end = GetPlbFromTimestamp(*ts_end);
		if (pLB_start == pLB_end) {
			if (!(found_pLB[pLB_start])) {
				found_pLB[pLB_start] = true;
				pLB_list_trp.push_back(pLB_start);
				pLB_timestamps_trp[pLB_start] = std::make_pair<Double_t, Double_t>(*ts_start, *ts_end);
			} else {
				if (*ts_start < pLB_timestamps_trp[pLB_start].first) {
					pLB_timestamps_trp[pLB_start].first = *ts_start;
				}
				if (*ts_end > pLB_timestamps_trp[pLB_start].second) {
					pLB_timestamps_trp[pLB_start].second = *ts_end;
				}
			}
		}		
	}
	
	sort(pLB_list_trp.begin(), pLB_list_trp.end());
	Int_t pLB_min = *(pLB_list_trp.begin());
	Int_t pLB_max = *(pLB_list_trp.end()-1);
	Int_t n_pLBs = pLB_max - pLB_min + 1;
	Int_t timestamp_min = TMath::FloorNint(pLB_timestamps_trp[pLB_min].first);
	Int_t timestamp_max = TMath::CeilNint(pLB_timestamps_trp[pLB_max].second);
	Int_t n_timestamps = timestamp_max - timestamp_min + 1;
	#ifdef DEBUG
	cout << "Processed timestamps based on TRP coarse timestamps:" << endl;
	cout << "\t PLB range: " << pLB_min << " - " << pLB_max << endl;
	cout << "\t Timestamp range: " << timestamp_min << " - " << timestamp_max << endl;
	#endif
	
	//Write timestamps to file
	ofstream new_ts_file;
	new_ts_file.precision(15);
	TString tsfilename = "PlbTrpTimestamps_"; tsfilename += run; tsfilename += ".dat";
	new_ts_file.open(tsfilename);

	for (vector<Int_t>::iterator it = pLB_list_trp.begin(); it != pLB_list_trp.end(); ++it) {
		new_ts_file << (*it) << " \t " << pLB_timestamps_trp[*it].first << " \t " << pLB_timestamps_trp[*it].second << endl;
		cout  << (*it) << " \t " << pLB_timestamps_trp[*it].first << " \t " << pLB_timestamps_trp[*it].second << endl;
	}
	
	new_ts_file.close();
#else
	//Load pre-existing timestamps
	TString trp_ts_path = "PlbTrpTimestamps_"; trp_ts_path += run; trp_ts_path += ".dat";
	cout << "Loading pre-existing pLB/TRP timestamps from" << trp_ts_path << endl;
	LoadPlbTimestampsTrp(trp_ts_path);
	sort(pLB_list_trp.begin(), pLB_list_trp.end());
	Int_t pLB_min = *(pLB_list_trp.begin());
	Int_t pLB_max = *(pLB_list_trp.end()-1);
	Int_t n_pLBs = pLB_max - pLB_min + 1;
	Int_t timestamp_min = TMath::FloorNint(pLB_timestamps_trp[pLB_min].first);
	Int_t timestamp_max = TMath::CeilNint(pLB_timestamps_trp[pLB_max].second);
	Int_t n_timestamps = timestamp_max - timestamp_min;
#endif
	cout << "pLB range is (" << pLB_min << ", " << pLB_max << "), n_pLBs = " << n_pLBs << endl;
	cout.precision(10);
	cout << "Timestamp range is (" << timestamp_min << ", " << timestamp_max << ")" << endl;
	
/***************************************************************************/
// Trigger rates
/***************************************************************************/
	cout << "*** Trigger Rates ***" << endl << endl;

#define REDOTRP
#ifdef REDOTRP
	//Declare histograms	
	// -- A note on the binning: the timestamp value attached to the trigger rates is the END of the interval in which triggers were counted. So, I think filling at (ts - 0.5) gives a better representation (bins correspond exactly to the intervals). 
	TH1F *h_ts_ntrig_ap = new TH1F("TsTriggersAP", ";ts; NTrigAP", n_timestamps, timestamp_min, timestamp_max);
	TH1F *h_ts_ntrig_av = new TH1F("TsTriggersAV", ";ts; NTrigAV", n_timestamps, timestamp_min, timestamp_max);
	TH1F *h_pLB_ps = new TH1F("PlbPrescale", ";pLB; PS", n_pLBs, pLB_min-0.5, pLB_max+0.5);
	
	//Fill trigger rates
	
	t_trp->SetBranchStatus("*", 0);
	t_trp->SetBranchStatus("L1_BGRP7_TAP", 1);
	t_trp->SetBranchStatus("L1_BGRP7_TAV", 1);
	t_trp->SetBranchStatus("TimeStamp", 1);
	t_trp->SetBranchStatus("L1_BGRP7_PS", 1);
	
	Double_t current_TimeStamp;
	Float_t current_L1_BGRP7_TAP;
	Float_t current_L1_BGRP7_TAV;
	Float_t current_L1_BGRP7_PS;
	Int_t current_pLB;
	
	t_trp->SetBranchAddress("L1_BGRP7_TAP", &current_L1_BGRP7_TAP);
	t_trp->SetBranchAddress("L1_BGRP7_TAV", &current_L1_BGRP7_TAV);
	t_trp->SetBranchAddress("L1_BGRP7_PS", &current_L1_BGRP7_PS);
	t_trp->SetBranchAddress("TimeStamp", &current_TimeStamp);
	
	Int_t last_pLB = -1;
	Double_t last_TimeStamp;
	for (Int_t i = 0; i < n_entries; i++) {
		if (i%10000 == 0) cout << "Entry " << i << " / " << n_entries << endl;
	
		t_trp->GetEntry(i);
		if (current_TimeStamp == last_TimeStamp) {
			#ifdef DEBUG
			cout << "Found repeated timestamp: " << current_TimeStamp << endl;
			#endif
			continue;
		}
		current_pLB = GetPlbFromTimestampTrp(current_TimeStamp);
		if (current_pLB <= 0) {
			#ifdef DEBUG
			if (i%1000 == 0) cout << "Timestamp not found: " << current_TimeStamp << endl;
			#endif
			continue;
		}
		
		h_ts_ntrig_ap->Fill(current_TimeStamp-0.5, current_L1_BGRP7_TAP);
		h_ts_ntrig_av->Fill(current_TimeStamp-0.5, current_L1_BGRP7_TAV);
		//h_pLB_ntrig_ap->Fill(current_pLB, current_L1_BGRP7_TAP);
		//h_pLB_ntrig_av->Fill(current_pLB, current_L1_BGRP7_TAV);
		
		if (current_pLB != last_pLB) h_pLB_ps->SetBinContent(current_pLB, current_L1_BGRP7_PS);
		
		last_TimeStamp = current_TimeStamp;
	}
	
	TFile *f_out1 = new TFile(out_path, "UPDATE");
	h_ts_ntrig_ap->Write("", TObject::kOverwrite);
	h_ts_ntrig_av->Write("", TObject::kOverwrite);
	h_pLB_ps->Write("", TObject::kOverwrite);
	f_out1->Close();
	
#else
	// -- Load preexisting timestamp-level trigger rates from file.
	TFile *f_trp_in = new TFile(out_path, "READ");
	TH1F *h_ts_ntrig_ap = (TH1F*)f_trp_in->Get("TsTriggersAP"); h_ts_ntrig_ap->SetDirectory(0);
	TH1F *h_ts_ntrig_av = (TH1F*)f_trp_in->Get("TsTriggersAV"); h_ts_ntrig_av->SetDirectory(0);
	TH1F *h_pLB_ps = (TH1F*)f_trp_in->Get("PlbPrescale"); h_pLB_ps->SetDirectory(0);
	f_trp_in->Close();
#endif	
	
	cout << "Done with trigger rates." << endl << endl;

/***************************************************************************/
// Event rates from TAG
/***************************************************************************/

#define D3PD
#ifdef D3PD
	cout << "*** Event Rates ***" << endl << endl;
//Declare histogram
	TH1F *h_ts_event_rate = new TH1F("TsEvents", ";ts;Events", n_timestamps, timestamp_min, timestamp_max); h_ts_event_rate->SetDirectory(0);
	
	TChain *ch = new TChain("InDetTrackTree");
	
	for (Int_t subjob = 1; subjob <= 23; subjob++) {
		TString D3PD_path = D3PD_folder; D3PD_path += "/user.spagan.005904.D3PD._000";
		if (subjob < 10) {
			D3PD_path += "0";
		}
		D3PD_path += subjob; D3PD_path += ".root";
		
		ch->Add(D3PD_path);
	}
	
	TBranch *b_ei_timestamp;   //!
	TBranch *b_ei_timestamp_ns;   //!
	TBranch *b_ei_lbn;   //!
	TBranch *b_trig_L1_TAV;
	
	std::vector<UInt_t> *current_trig_L1_TAV;
	UInt_t current_ei_timestamp;
	UInt_t current_ei_timestamp_ns;
	Double_t current_timestamp;
	UInt_t current_ei_lbn;
	
	ch->SetMakeClass(1);
	
	ch->SetBranchAddress("trig_L1_TAV", &current_trig_L1_TAV, &b_trig_L1_TAV);
	ch->SetBranchAddress("ei_timestamp", &current_ei_timestamp, &b_ei_timestamp);
	ch->SetBranchAddress("ei_timestamp_ns", &current_ei_timestamp_ns, &b_ei_timestamp_ns);
	ch->SetBranchAddress("ei_lbn", &current_ei_lbn, &b_ei_lbn);
	
	ch->SetBranchStatus("*", 0);
	ch->SetBranchStatus("trig_L1_TAV", 1);
	ch->SetBranchStatus("ei_timestamp", 1);
	ch->SetBranchStatus("ei_timestamp_ns", 1);
	ch->SetBranchStatus("ei_lbn", 1);
	

	Int_t n_entries2 = ch->GetEntries();
	cout << "TChain has " << n_entries2 << " entries." << endl;
	
	int selectItem = 60;
	Int_t total_triggers = 0;
	
	for (Long64_t i=1; i<n_entries2; i = i + 1) {
		if (i%100000 == 0) cout << "Event " << i << " / " << n_entries2 <<endl;
		
		ch->GetEntry(i);
		
		int selectWord = selectItem / 32;
		//if (!(current_trig_L1_TAV[selectWord] >> (selectItem%32) & 0x1)) continue;
		if (!((*current_trig_L1_TAV)[1] >> 28 & 0x1)) continue;
		//cout << "Found event." << endl;
		
		total_triggers++;
		
		current_timestamp = current_ei_timestamp + current_ei_timestamp_ns/1.e9;
		current_pLB = GetPlbFromTimestampTrp(current_timestamp);
		cout.precision(10);
		if ((total_triggers % 100000) == 0) cout << "On pLB " << current_pLB << " / ts = " << current_timestamp << endl;
		if (current_pLB < 0) continue;
		h_ts_event_rate->Fill(current_timestamp);
	}
	
	TFile *f_out2 = new TFile(out_path, "UPDATE");
	h_ts_event_rate->Write("", TObject::kOverwrite);
	f_out2->Close();
#endif

#ifdef REDOTAG
	cout << "*** Event Rates ***" << endl << endl;
//Declare histogram
	TH1F *h_ts_event_rate = new TH1F("TsEvents", ";ts;Events", n_timestamps, timestamp_min, timestamp_max); h_ts_event_rate->SetDirectory(0);
	
	TChain *c_InDetTrackTree = new TChain("InDetTrackTree");
	
	for (Int_t subjob = 1; subjob <= 23; subjob++) {
		TString D3PD_path = D3PD_folder; D3PD_path += "/user.spagan.005904.D3PD._000";
		if (subjob < 10) {
			D3PD_path += "0";
		}
		D3PD_path += subjob; D3PD_path += ".root";
		
		c_InDetTrackTree->Add(D3PD_path);
	}
	
	
	TFile *f_tag = new TFile(tag_path, "READ");
	
	//TTree *t_tag = (TTree*)f_tag->Get("POOLCollectionTree");
	//TTree *t_tag = (TTree*)f_tag->Get("InDetTrackTree");
	t_tag->SetBranchStatus("*", 0);
	t_tag->SetBranchStatus("L1PassedTrigMaskTAV3", 1);
	t_tag->SetBranchStatus("EventTime", 1);
	t_tag->SetBranchStatus("EventTimeNanoSec", 1);
	t_tag->SetBranchStatus("LumiBlockN", 1);
	
	UInt_t current_L1PassedTrigMaskTAV3;
	UInt_t current_EventTime;
	UInt_t current_EventTimeNanoSec;
	Double_t current_timestamp;
	UInt_t current_LumiBlockN;
	
	t_tag->SetBranchAddress("L1PassedTrigMaskTAV3", &current_L1PassedTrigMaskTAV3);
	t_tag->SetBranchAddress("EventTime", &current_EventTime);
	t_tag->SetBranchAddress("EventTimeNanoSec", &current_EventTimeNanoSec);
	t_tag->SetBranchAddress("LumiBlockN", &current_LumiBlockN);
	
	n_entries = t_tag->GetEntriesFast();
	
	int selectItem = 109;
	Int_t total_triggers = 0;
	
	for (int i=0; i<n_entries; i++) {
		if (i%100000 == 0) cout << "Event " << i << " / " << n_entries <<endl;
		t_tag->GetEntry(i);
		
		if (!(current_L1PassedTrigMaskTAV3 >> (selectItem%32) & 0x1)) continue;
		
		total_triggers++;
		
		current_timestamp = current_EventTime + current_EventTimeNanoSec/1.e9;
		current_pLB = GetPlbFromTimestampTrp(current_timestamp);
		cout.precision(10);
		if ((total_triggers % 100000) == 0) cout << "On pLB " << current_pLB << " / ts = " << current_timestamp << endl;
		if (current_pLB < 0) continue;
		h_ts_event_rate->Fill(current_timestamp);
	}
	
	TFile *f_out2 = new TFile(out_path, "UPDATE");
	h_ts_event_rate->Write("", TObject::kOverwrite);
	f_out2->Close();
	f_tag->Close();
#endif
#ifdef LOAD_TAG
	//Load pre-existing event rate histogram from output file
	TFile *f_read_event_rates = new TFile(out_path, "READ");
	TH1F *h_ts_event_rate = (TH1F*)f_read_event_rates->Get("TsEvents"); h_ts_event_rate->SetDirectory(0);
	f_read_event_rates->Close();
#endif

/***************************************************************************/
// Rebin timestamps into pLBs
/***************************************************************************/
	TH1F *h_pLB_ntrig_ap = new TH1F("PlbTriggersAP", ";pLB; NTrigAP", n_pLBs, pLB_min-0.5, pLB_max+0.5);
	TH1F *h_pLB_ntrig_av = new TH1F("PlbTriggersAV", ";pLB; NTrigAV", n_pLBs, pLB_min-0.5, pLB_max+0.5);
	TH1F *h_pLB_event_rate = new TH1F("PlbEvents", ";pLB;Events", n_pLBs, pLB_min-0.5, pLB_max+0.5);

	for (map<Int_t, std::pair<Double_t, Double_t> >::iterator it = pLB_timestamps_trp.begin(); it != pLB_timestamps_trp.end(); ++it) {
		Int_t current_pLB = (*it).first;
		Int_t ts_start = TMath::Nint((*it).second.first);
		Int_t ts_end = TMath::Nint((*it).second.second);
		
		#ifdef DEBUG
		cout << "pLB " << current_pLB << " : <" << ts_start << ", " << ts_end << ">" << endl;
		#endif

		for (Int_t current_ts = ts_start+1; current_ts < ts_end; current_ts++) {
			
			//HACK: Sometimes no triggers were recorded!! Why?
			if (h_ts_ntrig_ap->GetBinContent(h_ts_ntrig_ap->FindBin(current_ts)) < 1./1000.) continue;
			
			h_pLB_ntrig_ap->Fill(current_pLB, h_ts_ntrig_ap->GetBinContent(h_ts_ntrig_ap->FindBin(current_ts)));
			h_pLB_ntrig_av->Fill(current_pLB, h_ts_ntrig_av->GetBinContent(h_ts_ntrig_av->FindBin(current_ts)));
			h_pLB_event_rate->Fill(current_pLB, h_ts_event_rate->GetBinContent(h_ts_event_rate->FindBin(current_ts)));
			#ifdef DEBUG
			cout << "\t Events: " << h_ts_event_rate->GetBinContent(h_ts_event_rate->FindBin(current_ts)) << endl;
			cout << "\t Triggers: " << h_ts_ntrig_ap->GetBinContent(h_ts_ntrig_ap->FindBin(current_ts)) << endl;
			#endif
		}
		#ifdef DEBUG
		cout << "Totals triggers: " << h_pLB_ntrig_ap->GetBinContent(h_pLB_ntrig_ap->FindBin(current_pLB)) << endl;
		cout << "Totals events: " << h_pLB_event_rate->GetBinContent(h_pLB_event_rate->FindBin(current_pLB)) << endl;
		#endif
	}
	
	TFile *f_out3 = new TFile(out_path, "UPDATE");
	h_pLB_ntrig_ap->Write("", TObject::kOverwrite);
	h_pLB_ntrig_av->Write("", TObject::kOverwrite);
	h_pLB_event_rate->Write("", TObject::kOverwrite);
	f_out3->Close();

/***************************************************************************/
// Divide histograms to make live fraction
/***************************************************************************/
	TH1F *h_pLB_livefractions = (TH1F*)h_pLB_event_rate->Clone();
	h_pLB_livefractions->SetName("PlbLiveFractions");
	h_pLB_livefractions->SetTitle(";pLB;Live fractions;");
	
	h_pLB_livefractions->Divide(h_pLB_event_rate, h_pLB_ntrig_ap);
	
	TH1F *h_pLB_livefractions_scan;
	if (run == "182013") {
	  //Make another histogram just with scan points
	  h_pLB_livefractions_scan = (TH1F*)h_pLB_livefractions->Clone();
	  h_pLB_livefractions_scan->SetName("PlbLiveFractionsInScan");
	  
	  std::vector<std::pair<Int_t, Int_t> > scan_boundaries;
	  scan_boundaries.push_back(std::make_pair<Int_t, Int_t>(67, 115));
	  scan_boundaries.push_back(std::make_pair<Int_t, Int_t>(153, 201));
	  scan_boundaries.push_back(std::make_pair<Int_t, Int_t>(224, 272));
	  scan_boundaries.push_back(std::make_pair<Int_t, Int_t>(294, 342));
	  scan_boundaries.push_back(std::make_pair<Int_t, Int_t>(415, 463));
	  scan_boundaries.push_back(std::make_pair<Int_t, Int_t>(482, 530));
	  
	  for (Int_t i = 1; i <= h_pLB_livefractions_scan->GetNbinsX(); i++) {
	    Int_t pLB = TMath::Nint(h_pLB_livefractions_scan->GetBinCenter(i));
	    for (std::vector<std::pair<Int_t, Int_t> >::iterator it = scan_boundaries.begin(); it != scan_boundaries.end(); ++it) {
	      if ((pLB >= (*it).first) && (pLB <= (*it).second)) {
		if (((pLB - (*it).first) % 2) != 0) {
		  h_pLB_livefractions_scan->SetBinContent(i, 0);
		}
	      } else {
		h_pLB_livefractions_scan->SetBinContent(i, 0);
	      }
	    }
	  }
	}
	
	TFile *f_out4 = new TFile(out_path, "UPDATE");
	h_pLB_livefractions->Write("", TObject::kOverwrite);
	if (run == "182013") {
	  h_pLB_livefractions_scan->Write("", TObject::kOverwrite);
	}
	f_out4->Close();
	f_trp->Close();
}


void LoadPlbTimestamps(TString path) {
	ifstream f_in(path.Data());
	if (not f_in.is_open()) {
		cerr << "Unable to open timestamps text file: " << path << endl;
		exit(1);
	}
	
	Int_t current_pLB;
	Double_t ts_start, ts_end;
	TString everything_else;
	
	while (f_in.good()) {
		string line;
		getline(f_in, line);
		if (line.empty()) continue;
		stringstream istr(line);
		istr >> current_pLB >> ts_start >> ts_end >> everything_else;
		
		//HACK: No HLT trigger after pLB 700
		//if (current_pLB > 700) continue;
		
		pLB_timestamps[current_pLB] = std::make_pair<Double_t, Double_t>(ts_start, ts_end);
		
		pLB_list.push_back(current_pLB);
		
		lb_duration[current_pLB] = ts_end - ts_start;
		
		#ifdef DEBUG
		cout << "Set timestamp " << current_pLB << " : ("<< ts_start << ", " << ts_end << "), lb_duration[" << current_pLB << "] = " << lb_duration[current_pLB] << endl;
		#endif
	}
}

void LoadPlbTimestampsTrp(TString path) {
	ifstream f_in(path.Data());
	if (not f_in.is_open()) {
		cerr << "Unable to open timestamps text file: " << path << endl;
		exit(1);
	}
	
	Int_t current_pLB;
	Double_t ts_start, ts_end;
	TString everything_else;
	
	while (f_in.good()) {
		string line;
		getline(f_in, line);
		if (line.empty()) continue;
		stringstream istr(line);
		istr >> current_pLB >> ts_start >> ts_end;
		
		//HACK: No HLT trigger after pLB 700
		//if (current_pLB > 700) continue;

		pLB_timestamps_trp[current_pLB] = std::make_pair<Double_t, Double_t>(ts_start, ts_end);
		
		pLB_list_trp.push_back(current_pLB);
		
		lb_duration_trp[current_pLB] = ts_end - ts_start;
		
		#ifdef DEBUG
		cout << "Set timestamp " << current_pLB << " : ("<< ts_start << ", " << ts_end << "), lb_duration_trp[" << current_pLB << "] = " << lb_duration_trp[current_pLB] << endl;
		#endif
	}
}

Int_t GetPlbFromTimestamp( Double_t timestamp ) {
	for ( vector<Int_t>::iterator pLB = pLB_list.begin(); pLB != pLB_list.end(); ++pLB ) {
		if ((timestamp >= pLB_timestamps[*pLB].first) && (timestamp <= pLB_timestamps[*pLB].second)) return *pLB;
		
		//if ((timestamp < pLB_timestamps[*pLB].first) && (timestamp < pLB_timestamps[*pLB].second)) return -1;
	}
	//cout << "In getPLBFromTimestamp( Double_t timestamp ): failed to find timestamp " << timestamp << endl;
	return -1;
}

Int_t GetPlbFromTimestampTrp( Double_t timestamp ) {
	for ( vector<Int_t>::iterator pLB = pLB_list_trp.begin(); pLB != pLB_list_trp.end(); ++pLB ) {
		if ((timestamp > pLB_timestamps_trp[*pLB].first + 1.e-9) && (timestamp <= pLB_timestamps_trp[*pLB].second + 1.e-9)) return *pLB;
		
		//if ((timestamp < pLB_timestamps_trp[*pLB].first) && (timestamp < pLB_timestamps_trp[*pLB].second)) return -1;

	}
	//cout << "In getPLBFromTimestamp( Double_t timestamp ): failed to find timestamp " << timestamp << endl;
	return -1;
}

