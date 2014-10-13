# -*- coding: utf-8 -*-
#! /usr/bin/env python
import ROOT
from ROOT import TFile, TH1F, TCanvas, TLegend

ROOT.gROOT.SetBatch(True)

sample = "8TeV_VtxLumi_mumax75"
sample = "8TeV_VtxLumi_mumax20"

if sample == "7TeV_Normal":
	rootfile = '/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VtxTruthMatchResults/mc_7TeV_17.2_normal_pythia8_pu_bs55/InDetTrackD3PD_results.root'
elif sample == "7TeV_VtxLumi":
	rootfile = '/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VtxTruthMatchResults/mc_7TeV_17.2_normal_pythia8_pu_bs55/InDetTrackD3PD_results.root'
elif sample == "7TeV_VtxLumi_narrow":
	rootfile = '/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VtxTruthMatchResults/mc_7TeV_17.2_normal_pythia8_pu_bs55/InDetTrackD3PD_results.root'
elif sample == "8TeV_Normal":
	rootfile = '/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VtxTruthMatchResults/mc_7TeV_17.2_normal_pythia8_pu_bs55/InDetTrackD3PD_results.root'
elif sample == "8TeV_VtxLumi_mumax20":
	rootfile = '/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VtxTruthMatchResults/mc_8TeV_17.2_VtxLumi_2newsets/InDetTrackD3PD_results_sample15.root'
	#rootfile = '/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VtxTruthMatchResults/mc_8TeV_17.2_VtxLumi_LowMuSample/InDetTrackD3PD_results.root'
	xupp = 30
	xlow = 0
elif sample == "8TeV_VtxLumi_mumax75":
	rootfile = '/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VtxTruthMatchResults/mc_8TeV_17.2_VtxLumi_2newsets/InDetTrackD3PD_results_sample10.root'
	#rootfile = '/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VtxTruthMatchResults/mc_8TeV_17.2_VtxLumi_HighMuSample/InDetTrackD3PD_results.root'
	xupp = 80
	xlow = 30

f_histos = TFile(rootfile,"READ")
if ("VtxLumi" in sample):
	ntrkcuts = ["2", "3", "4", "5"]
else: ntrkcuts = ["2", "5", "7","10"]

print "ntrkcuts:", ntrkcuts

colors = [ROOT.kGreen, ROOT.kRed, ROOT.kBlue, ROOT.kBlack]
hs_splitfraction = []
hs_significance = []

legend = TLegend(0.2,0.65,0.4,0.85)
legend.SetFillColor(0)

legend_2 = TLegend(0.2,0.75,0.4,0.85)
legend_2.SetFillColor(0)

print "Using sample: ", sample

for ntrkcut,color in zip(ntrkcuts,colors):
	#Split Fraction
	h_splits = f_histos.Get("hist/h_splits_NGenInt_NTrk"+ntrkcut)
	h_all = f_histos.Get("hist/h_all_NGenInt_NTrk"+ntrkcut) 
	h_ratio = h_splits.Clone()
	h_ratio.Sumw2()
	h_ratio.Divide(h_all)
	h_ratio.Scale(100)
	h_ratio.SetName("h_splitfraction_NTrk"+ntrkcut)
	h_ratio.SetMarkerColor(color)
	h_ratio.SetLineColor(color)
	h_ratio.SetStats(0)
	hs_splitfraction.append(h_ratio)
	if (ntrkcut=="2"):
		legend_2.AddEntry(h_ratio,"NTrkCut "+ntrkcut, "P")
	else: legend.AddEntry(h_ratio,"NTrkCut "+ntrkcut, "P")

	#Delta z significance
	h_real_2D = f_histos.Get("hist/h_real_dzsig_NGenInt_NTrk"+ntrkcut)
	h_split_2D = f_histos.Get("hist/h_split_dzsig_NGenInt_NTrk"+ntrkcut)
	h_real = h_real_2D.ProjectionX()
	h_split = h_split_2D.ProjectionX()
	hs_significance.append(h_real)
	hs_significance.append(h_split)

print hs_significance

for hist in hs_splitfraction:
	hist.SetTitle("")
	hist.GetXaxis().SetTitle("Generated Interactions")
	hist.GetXaxis().SetTitleOffset(0.9)
	hist.GetXaxis().SetTitleSize(0.045)
	hist.GetXaxis().SetLabelSize(0.035)
	hist.GetXaxis().SetRangeUser(xlow,xupp)
	hist.GetYaxis().SetTitle("Split Fraction (%)")
	hist.GetYaxis().SetTitleOffset(0.9)
	hist.GetYaxis().SetTitleSize(0.045)
	hist.GetYaxis().SetLabelSize(0.035)
	hist.SetMaximum(0.30)
	if (sample == "8TeV_VtxLumi_mumax75"):
		hist.SetMinimum(0)
	hist.SetMarkerStyle(8)
	hist.SetMarkerSize(0.8)
	hist.SetOption("PE")

canvas = TCanvas("canvas","canvas",1200,800)
hs_splitfraction[1].Draw("PE")
hs_splitfraction[2].Draw("PEsame")
hs_splitfraction[3].Draw("PEsame")
legend.Draw("same")
canvas.Print("SplitFractionOtherNTrks_"+sample+".pdf")

canvas_2 = TCanvas("canvas_2","canvas_2",1200,800)
hs_splitfraction[0].Draw("PE")
legend_2.Draw("same")
canvas_2.Print("SplitFractionNTrk2_"+sample+".pdf")

f_splits = TFile("SplitVerticesAnalysis_"+sample+".root","RECREATE")
for hist in hs_splitfraction:
	hist.Write()
for hist in hs_significance:
	hist.Write()
f_splits.Close()





