/**
  *	Toy MC to study how whether the dz gaussian method actually gets the number of masked vertices correct.
  * Author: David Yu
  *	January 20, 2011
  **/

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <stdlib.h>

#include "TCanvas.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMath.h"

#include "atlasstyle/AtlasStyle.h"
#include "atlasstyle/AtlasLabels.h"
#include "atlasstyle/AtlasUtils.h"

#include "PileupCorrections/PileupMaskingCorrection.h"

using namespace std;

Double_t runToyMc(Float_t p_sigma_z, TF1 *p_pmask_dz, Float_t p_mu, Int_t p_events, TString p_tag, bool save_plots);
bool vertexMask(TF1 *p_pmask_dz, Float_t z1, Float_t z2);
TGraphErrors *tg_mutruth_murecon;
TGraphErrors *tg_mucorr_murecon;
TRandom *super_rand;
TH1D *h_pmask_dz;
Double_t MaskingPdf(Double_t *x, Double_t *par);

int main() {
	SetAtlasStyle();

	Float_t sigma_z = 55;
	Int_t events = 1000000;
	super_rand = new TRandom3();
	
	TF1 *f_pmask_dz = new TF1("f_pmask_dz", MaskingPdf, -500., 500., 2);
	f_pmask_dz->SetNpx(10000);
	f_pmask_dz->SetParameter(0, 3.);
	f_pmask_dz->SetParameter(1, 3.);

	TString tag = "plateau3_decay3_mu2";

	TH1F *h_pmask = new TH1F("h_pmask", "h_pmask", 500, 0.03, 0.08);
	Double_t pmask_avg = 0.;
	Double_t pmask_avg2 = 0.;

	Int_t trials = 1000;
	for (int i = 0; i < trials; i++) {
		Double_t current_pmask;
		if (i == 0) {
			current_pmask = runToyMc(sigma_z, f_pmask_dz, 1., events, tag, true);
		} else {
			current_pmask = runToyMc(sigma_z, f_pmask_dz, 1., events, tag, false);
		}
		h_pmask->Fill(current_pmask);
		pmask_avg += current_pmask;
		pmask_avg2 += TMath::Power(current_pmask, 2);
	}

	pmask_avg = pmask_avg / trials;
	pmask_avg2 = pmask_avg2 / trials;
	cout << "Average p_mask = " << pmask_avg << endl;
	Double_t pmask_rms = TMath::Sqrt(pmask_avg2 - TMath::Power(pmask_avg, 2));
	cout << "RMS = " << pmask_rms << endl;

	TCanvas *c_pmask = new TCanvas("c_pmask", "c_pmask", 1200, 800);
	TLegend *l_pmask = new TLegend(0.2, 0.6, 0.4, 0.9);
	l_pmask->SetFillColor(0);
	l_pmask->SetBorderSize(1);

	h_pmask->Draw("hist");
	TLine *pmask_truth = new TLine(0.0614055, h_pmask->GetMinimum(), 0.0614055, h_pmask->GetMaximum());
	pmask_truth->SetLineStyle(2);
	pmask_truth->SetLineColor(2);
	pmask_truth->Draw("same");

	l_pmask->AddEntry(h_pmask, "Toy MC", "l");
	l_pmask->AddEntry(pmask_truth, "Truth", "l");
	l_pmask->Draw();

	TFile *f_pmask = new TFile(TString("/u/dryu/Luminosity/Data/MaskingCorrection") + "/toy/toy_pmask.root", "RECREATE");
	h_pmask->Write();
	c_pmask->Write();
	c_pmask->SaveAs(TString("/u/dryu/Luminosity/Data/MaskingCorrection") + "/toy/c_pmask.pdf");
	f_pmask->Close();
}

Double_t runToyMc(Float_t p_sigma_z, TF1 *p_pmask_dz, Float_t p_mu, Int_t p_events, TString p_tag, bool save_plots) {

	cout << "*** Starting toy masking MC ***" << endl;
	cout << "\t mu = \t" << p_mu << endl;
	cout << "\t events = \t" << p_events << endl;
	cout << "\t tag = \t" << p_tag << endl;

	TString hname;
	
	hname = "h_dz_truth_"; hname += p_tag;
	TH1D *h_dz_truth = new TH1D(hname, hname, 4000, -500, 500);
	
	hname = "h_dz_masked_first_"; hname += p_tag;
	TH1D *h_dz_masked_first = new TH1D(hname, hname, 4000, -500, 500);
	
	hname = "h_dz_masked_notfirst_"; hname += p_tag;
	TH1D* h_dz_masked_notfirst = new TH1D(hname, hname, 4000, -500, 500);
	
	hname = "h_dz_recon_"; hname += p_tag;
	TH1D *h_dz_recon = new TH1D(hname, hname, 4000, -500, 500);

	hname = "h_z_"; hname += p_tag;
	TH1D *h_z = new TH1D(hname, hname, 4000, -500, 500);
	
	TRandom3 *rand1 = new TRandom3(0);
	TRandom3 *rand2 = new TRandom3(0);

	cout << endl << "Starting event loop" << endl;
	
	Int_t nvtx_truth = 0;
	Int_t nvtx_masked = 0;
	Int_t nvtx_recon = 0;
	
	for (Int_t i = 0; i < p_events; i++) {
		if ((i % 100000) == 0) {
			cout << "Event " << i << endl;
		}
	
		std::vector<std::pair<Float_t, bool> > z_recon; // store vertex z positions, and whether the vertex was reconstructed.
		Int_t nvtx_gen = rand1->Poisson(p_mu);
		nvtx_truth += nvtx_gen;
		
		// Generate z positions
		for (Int_t nv1 = 0; nv1 < nvtx_gen; nv1++) {
			
			Float_t current_z = rand2->Gaus(0, p_sigma_z);
			h_z->Fill(current_z);
			bool reconstructed = true;
			
			z_recon.push_back(std::make_pair<Float_t, bool>(current_z, reconstructed));
		}
		
		// Do masking
		
		for (std::vector<std::pair<Float_t, bool> >::iterator nv1 = z_recon.begin(); nv1 != z_recon.end(); nv1++) {
			
		
			for (std::vector<std::pair<Float_t, bool> >::iterator nv2 = z_recon.begin(); nv2 != nv1; ++nv2) {
				if (nv2 >= nv1) continue; // nv1 is masked by a prior vertex.
				
				h_dz_truth->Fill((*nv1).first - (*nv2).first);
				
				//If nv2 was not reconstructed, continue
				if ((*nv2).second == false) continue;
				
				//Else check if nv1 is masked
				if (vertexMask(p_pmask_dz, (*nv1).first, (*nv2).first)) {
				
					// nv1 is masked.
					// If reconstructed == true, this is the masking occurence; fill h_dz_masked_first
					if ((*nv1).second == true) {
						h_dz_masked_first->Fill((*nv1).first - (*nv2).first);
					} else {
						// Vertex has already been masked.
						h_dz_masked_notfirst->Fill((*nv1).first - (*nv2).first);
					}
					
					(*nv1).second = false;
				} else {
					//nv1 is not masked by nv2. Note that this doesn't mean nv1 is unmasked; it could be masked by a later nv2.
				}
			} // end nv2 loop
			if ((*nv1).second) {
				nvtx_recon++;
			} else {
				nvtx_masked++;
			}
		} // end nv1 loop
		
		// Fill reconstructed dz histogram
		for (std::vector<std::pair<Float_t, bool> >::iterator nv1 = z_recon.begin(); nv1 != z_recon.end(); ++nv1) {
			if ((*nv1).second == false) continue;
			
			for (std::vector<std::pair<Float_t, bool> >::iterator nv2 = z_recon.begin(); nv2 < nv1; ++nv2) {
				if (nv2 >= nv1) continue;
				if ((*nv2).second == false) continue;
				h_dz_recon->Fill((*nv1).first - (*nv2).first);
			}
		}
		
	} // end event loop
	
	// Print out summary, and set a point in the TGraph
	cout << "Summary:" << endl;
	cout << "NVtx (truth): " << nvtx_truth << endl;
	cout << "NVtx (recon): " << nvtx_recon << endl;
	cout << "NVtx (masked): " << nvtx_masked << endl;
	cout << "masked / truth = " << (float)nvtx_masked / (float)nvtx_truth << endl;

	// Pileup masking correction
	PileupMaskingCorrection *pmc = new PileupMaskingCorrection(h_dz_recon, p_tag);
	pmc->low_stats = false;
	pmc->rebin_factor = 1;
	cout << "[debug] h_z->Integral() = " << h_z->Integral() << endl;
	pmc->GenerateDzDistribution(h_z, p_tag);
	pmc->GenerateNewPmask();
	pmc->GenerateCorrection(pmc->GetExpectedDzDistribution());

	// Truth p_mask
	Float_t truth_p_mask = 0.;
	Int_t points = 100000;
	Float_t dz_min = -25;
	Float_t dz_max = 25;
	for (int i = 0; i < points; i++) {
		Float_t current_dz = dz_min + (dz_max - dz_min) * i / points;
		truth_p_mask += (dz_max - dz_min) / points * p_pmask_dz->Eval(current_dz) / (TMath::Sqrt(2 * TMath::Pi()) * TMath::Sqrt(2) * p_sigma_z) * TMath::Exp(-1. * TMath::Power(current_dz, 2) / (2 * TMath::Power(TMath::Sqrt(2) * p_sigma_z, 2)));
	}
	cout << "Truth p_mask = " << truth_p_mask << endl;
	
	// Canvases
	if (save_plots) {

		// -- DZ distributions: truth, reconstructed. 
		TString cname = "c_dz_"; cname += p_tag;
		TCanvas *c_dz = new TCanvas(cname, cname, 1200, 800);
		TLegend *l_dz = new TLegend(0.75, 0.8, 0.95, 0.95);
		l_dz->SetFillColor(0);
		l_dz->SetBorderSize(1);
		Int_t rebin_factor = 4;
		
		h_dz_truth->SetLineColor(kBlack);
		h_dz_truth->Rebin(rebin_factor);
		//h_dz_truth->Scale(1./16.);
		h_dz_truth->SetStats(0);
		h_dz_truth->GetXaxis()->SetTitle("#Deltaz (mm)");
		h_dz_truth->GetYaxis()->SetTitle("Entries / mm");
		h_dz_truth->Draw("hist");
		
		h_dz_recon->SetLineColor(kGreen);
		h_dz_recon->Rebin(rebin_factor);
		//h_dz_recon->Scale(1./16.);
		h_dz_recon->Draw("hist same");
		
		h_dz_masked_first->SetLineColor(kRed);
		h_dz_masked_first->Rebin(rebin_factor);
		//h_dz_masked_first->Scale(1./16.);
		h_dz_masked_first->Draw("hist same");	

		l_dz->AddEntry(h_dz_truth, "Expected #Delta z", "l");
		l_dz->AddEntry(h_dz_recon, "Observed #Delta z", "l");
		l_dz->AddEntry(h_dz_masked_first, "Masked vertices", "l");
		l_dz->Draw();
		
		// -- Masking PDF
		cname = "c_pmask_dz_"; cname += p_tag;
		TCanvas *c_pmask_dz = new TCanvas(cname, cname, 1200, 800);
		TH1F *frame = new TH1F("frame", "frame", 100, -50, 50);
		frame->SetMinimum(-0.1);
		frame->SetMaximum(1.1);
		frame->GetXaxis()->SetTitle("#Deltaz (mm)");
		frame->GetYaxis()->SetTitle("p_{mask}(#Deltaz)");
		frame->Draw();
		p_pmask_dz->SetLineColor(kRed);
		p_pmask_dz->Draw("same");
		TH1D *h_pmask_dz_measured = (TH1D*)pmc->GetDifferentialPmask();
		h_pmask_dz_measured->SetMarkerStyle(20);
		h_pmask_dz_measured->SetMarkerColor(kBlack);
		h_pmask_dz_measured->SetMarkerSize(1);
		h_pmask_dz_measured->Draw("p same");

		TLegend *l_pmask_dz = new TLegend(0.2, 0.7, 0.4, 0.9);
		l_pmask_dz->SetFillColor(0);
		l_pmask_dz->SetBorderSize(1);
		l_pmask_dz->AddEntry(p_pmask_dz, "Input p_{mask}(#Deltaz)", "l");
		l_pmask_dz->AddEntry(h_pmask_dz_measured, "Measured p_{mask}(#Deltaz)", "p");
		l_pmask_dz->Draw();

	/*
		pmc->GetExcludedGaussian()->SetLineColor(kCyan);
		pmc->GetExcludedGaussian()->SetLineStyle(2);
		pmc->GetExcludedGaussian()->Draw("same");
	*/

		TFile *f_out = new TFile(TString("/u/dryu/Luminosity/Data/MaskingCorrection") + TString("/toy/toy_") + p_tag + TString(".root"), "RECREATE");
		h_dz_truth->Write();
		h_dz_masked_first->Write();
		h_dz_masked_notfirst->Write();
		h_dz_recon->Write();
		c_dz->Write();
		c_dz->SaveAs(TString("/u/dryu/Luminosity/Data/MaskingCorrection") + TString("/toy/") + c_dz->GetName() + TString(".pdf"));
		c_pmask_dz->Write();
		c_pmask_dz->SaveAs(TString("/u/dryu/Luminosity/Data/MaskingCorrection") + TString("/toy/") + c_pmask_dz->GetName() + TString(".pdf"));

		pmc->Save("toy_pmc.root", p_tag, true);

		f_out->Close();
	}

	return pmc->GetTotalPmask();
}

/* Flat */
/*
bool vertexMask(Float_t dz_mask, Float_t z1, Float_t z2) {
	// Returns true if two vertices mask each other.
	if (TMath::Abs(z1 - z2) < dz_mask) {
		return true;
	} else {
		return false;
	}
}*/

/*Pyramid*/
/*
bool vertexMask(Float_t dz_mask, Float_t z1, Float_t z2) {
	// Returns true if two vertices mask each other.
	// Sawtooth masking function
	if (TMath::Abs(z1 - z2) >= dz_mask) {
		return false;
	} else {
		Float_t p_mask = -1./3 * TMath::Abs(z1 - z2) + 1.;
		Float_t r = super_rand->Rndm();
		if (r < p_mask) {
			return true;
		} else {
			return false;
		}
	}
}
*/

/*Input histogram*/
bool vertexMask(TF1 *p_pmask_dz, Float_t z1, Float_t z2) {
	// Returns true if two vertices mask each other.
	Float_t p_mask = p_pmask_dz->Eval(z1 - z2);
	Float_t r = super_rand->Rndm();
	if (r < p_mask) {
		return true;
	} else {
		return false;
	}
}

/* Input function */
Double_t MaskingPdf(Double_t *x, Double_t *par) {

	/**
	  *	Masking probability as fx of dz.
	  * par[0] = width of central plateau
	  *	par[1] = size of exponential falloff
	  */

	Double_t val;
	if (TMath::Abs(x[0]) < par[0]) {
		val = 1.;
	} else {
		Double_t y = TMath::Abs(x[0]) - par[0];
		val = TMath::Exp(-1. * y / par[1]);
	}

	return val;

}


