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
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

using namespace std;

int main() {
	TFile *f_in = new TFile("../test2.root");
	TTree *t = (TTree*)f_in->Get("LumiVtx");
	
	std::vector<Int_t> nTrkCuts;
	nTrkCuts.push_back(2);
	nTrkCuts.push_back(5);
	nTrkCuts.push_back(7);
	nTrkCuts.push_back(10);
	
	Long64_t n_entries = t->GetEntriesFast();
	
	std::map<Int_t, std::map<Int_t, Float_t> > lumi_nevt;
	std::map<Int_t, Float_t> lumi_bcm;
	
	std::map<Int_t, TGraphAsymmErrors*> tg_lumi_nevt;
	std::map<Int_t, TGraphAsymmErrors*> tg_lumi_ratio;
	std::map<Int_t, Int_t> colorMap;
	colorMap[2] = 1;
	colorMap[5] = 4;
	colorMap[7] = 2;
	colorMap[10] = 3;
	for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
		TString name = "tg_lumi_nevt"; name += *nTrkCut;
		tg_lumi_nevt[*nTrkCut] = new TGraphAsymmErrors(n_entries);
		tg_lumi_nevt[*nTrkCut]->SetName(name);
		tg_lumi_nevt[*nTrkCut]->SetMarkerStyle(20);
		tg_lumi_nevt[*nTrkCut]->SetMarkerColor(colorMap[*nTrkCut]);
		
		name = "tg_lumi_ratio"; name += *nTrkCut;
		tg_lumi_ratio[*nTrkCut] = new TGraphAsymmErrors(n_entries);
		tg_lumi_ratio[*nTrkCut]->SetMarkerStyle(20);
		tg_lumi_ratio[*nTrkCut]->SetMarkerColor(colorMap[*nTrkCut]);
		tg_lumi_ratio[*nTrkCut]->SetName(name);
	}
	TGraph* tg_lumi_bcm = new TGraph(n_entries);
	tg_lumi_bcm->SetName("tg_lumi_bcm");
	
	
	t->SetBranchStatus("*", 1);
	for (int i=0; i<n_entries; i++) {
		t->GetEntry(i);
		Int_t current_pLB = t->GetLeaf("LB")->GetValue(0);
		
		lumi_bcm[current_pLB] = t->GetLeaf("Linst_bcm")->GetValue(0);
		tg_lumi_bcm->SetPoint(i, current_pLB, lumi_bcm[current_pLB]);

		for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
			TString branch_nevt = "Linst_NEvt"; branch_nevt += *nTrkCut;
			lumi_nevt[*nTrkCut][current_pLB] = t->GetLeaf(branch_nevt)->GetValue(0);
			tg_lumi_nevt[*nTrkCut]->SetPoint(i, current_pLB, lumi_nevt[*nTrkCut][current_pLB]);
		
			if (lumi_bcm[current_pLB] > 0) {
				cout << "i = " << i << ", ratio = " << lumi_nevt[*nTrkCut][current_pLB] / lumi_bcm[current_pLB] << endl;	
				tg_lumi_ratio[*nTrkCut]->SetPoint(i, lumi_bcm[current_pLB], (lumi_nevt[*nTrkCut][current_pLB] / lumi_bcm[current_pLB] - 1.) * 100.);
			}

			TString branch_nevt_errdown = "Linst_NEvt"; branch_nevt_errdown += *nTrkCut; branch_nevt_errdown += "_errdown";
			TString branch_nevt_errup = "Linst_NEvt"; branch_nevt_errup += *nTrkCut; branch_nevt_errup += "_errup";
			tg_lumi_nevt[*nTrkCut]->SetPointError(i, 0, 0, t->GetLeaf(branch_nevt_errdown)->GetValue(0) , t->GetLeaf(branch_nevt_errup)->GetValue(0));
			
			if (lumi_bcm[current_pLB] > 0) {
				tg_lumi_ratio[*nTrkCut]->SetPointError(i, 0, 0, t->GetLeaf(branch_nevt_errup)->GetValue(0)/lumi_bcm[current_pLB]*100., t->GetLeaf(branch_nevt_errup)->GetValue(0)/lumi_bcm[current_pLB]*100.);
			}
		}
		
	}
	
	TCanvas *c_ratio = new TCanvas("c_ratio", "c_ratio", 1000, 800);
	c_ratio->SetRightMargin(0.2);
	TLegend *l_ratio = new TLegend(0.82,0.3,1.0,0.7);
	l_ratio->SetFillColor(0);
	l_ratio->SetBorderSize(1);
	
	bool drawFirst = true;
	for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
		if (drawFirst) {
			tg_lumi_ratio[*nTrkCut]->GetHistogram()->SetMinimum(-10.);
			tg_lumi_ratio[*nTrkCut]->GetHistogram()->SetMaximum(8.);
			tg_lumi_ratio[*nTrkCut]->GetXaxis()->SetTitle("L_{bcm} (#mub^{-1}s^{-1})");
			tg_lumi_ratio[*nTrkCut]->GetYaxis()->SetTitle("L_{vtx-evt} / L_{bcm} - 1 (%)");
			tg_lumi_ratio[*nTrkCut]->Draw("ap");
			drawFirst = false;
		} else {
			tg_lumi_ratio[*nTrkCut]->Draw("p");
		}
		TString label = "NTrk >= "; label += *nTrkCut;
		l_ratio->AddEntry(tg_lumi_ratio[*nTrkCut], label, "p");
	}
	l_ratio->Draw();
	c_ratio->SaveAs("c_ratio.eps");

	TCanvas *c_ratio_log = new TCanvas("c_ratio_log", "c_ratio_log", 1000, 800);
	c_ratio_log->SetRightMargin(0.2);
	c_ratio_log->SetLogx();
	TLegend *l_ratio_log = new TLegend(0.82,0.3,1.0,0.7);
	l_ratio_log->SetFillColor(0);
	l_ratio_log->SetBorderSize(1);
	
	drawFirst = true;
	for (vector<Int_t>::iterator nTrkCut = nTrkCuts.begin(); nTrkCut != nTrkCuts.end(); ++nTrkCut) {
		if (drawFirst) {
			tg_lumi_ratio[*nTrkCut]->GetHistogram()->SetMinimum(-10.);
			tg_lumi_ratio[*nTrkCut]->GetHistogram()->SetMaximum(8.);
			tg_lumi_ratio[*nTrkCut]->GetXaxis()->SetTitle("L_{bcm} (#mub^{-1}s^{-1})");
			tg_lumi_ratio[*nTrkCut]->GetYaxis()->SetTitle("L_{vtx-evt} / L_{bcm} - 1 (%)");
			tg_lumi_ratio[*nTrkCut]->Draw("ap");
			drawFirst = false;
		} else {
			tg_lumi_ratio[*nTrkCut]->Draw("p");
		}
		TString label = "NTrk >= "; label += *nTrkCut;
		l_ratio_log->AddEntry(tg_lumi_ratio[*nTrkCut], label, "p");
	}
	l_ratio_log->Draw();
	c_ratio_log->SaveAs("c_ratio_log.eps");
}