# -*- coding: utf-8 -*-
#! /usr/bin/env python
import ROOT
from ROOT import TFile, TH1F, TH2F, TCanvas, TLegend

ROOT.gROOT.SetBatch(True)

samples = ["MC_HighMu","MC_LowMu","VdMScanBCID1","VdMScanBCID2361","VdMScanBCID2881","PhysicsRun"]
hists = []

for sample in samples:
	print sample
	if sample == "MC_HighMu":
		rootfile = '/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VtxTruthMatchResults/mc_8TeV_17.2_VtxLumi_HighMuSample/InDetTrackD3PD_results.root'
		hist = "hist/h_privtx_z_NGenInt"
		color = ROOT.kGray+2
	if sample == "MC_LowMu":
		rootfile = '/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VtxTruthMatchResults/mc_8TeV_17.2_VtxLumi_LowMuSample/InDetTrackD3PD_results.root'
		hist = "hist/h_privtx_z_NGenInt"
		color = ROOT.kBlue
	if sample == "VdMScanBCID1":
		rootfile = '/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/VdMScan-215021/17.2-VtxLumi/InDetTrackD3PD_results_NormalD3PD_05Aug2014.root'
		hist = "hist/PriVtxZpLB_BCID1"
		color = ROOT.kBlack
	if sample == "VdMScanBCID2361":
		rootfile = '/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/VdMScan-215021/17.2-VtxLumi/InDetTrackD3PD_results_NormalD3PD_05Aug2014.root'
		hist = "hist/PriVtxZpLB_BCID2361"
		color = ROOT.kMagenta
	if sample == "VdMScanBCID2881":
		rootfile = '/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/VdMScan-215021/17.2-VtxLumi/InDetTrackD3PD_results_NormalD3PD_05Aug2014.root'
		hist = "hist/PriVtxZpLB_BCID2881"
		color = ROOT.kYellow
	if sample == "PhysicsRun":
		rootfile = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/213539/17.2-VtxLumi/InDetTrackD3PD_results.root"
		hist = "hist/PriVtxZpLB_BCID0"
		color = ROOT.kRed
	f_hist = TFile(rootfile,"READ")
	h_2d = f_hist.Get(hist)
	#if( sample == "VdMScanBCID1" or sample == "VdMScanBCID721" or sample == "VdMScanBCID1821" ):
	h_px = h_2d.ProjectionX()
	#else: h_px = h_2d.ProjectionX("px",1,10)
	h_px.SetDirectory(0)
	h_px.SetLineColor(color)
	h_px.SetMarkerColor(color)
	h_px.SetName(sample)
	h_px.Scale(1/h_px.Integral())
	hists.append(h_px)

f_zdist = TFile("TypicalZDistributions_Projection.root","RECREATE")
for hist in hists:
	hist.Write()
f_zdist.Close()

legend = ROOT.TLegend(0.65,0.65,0.85,0.85)
legend.SetFillColor(0)
legend.AddEntry(hists[1],"Low #mu MC sample", "l")
legend.AddEntry(hists[2],"VdM Scan 15, BCID 1", "l")
legend.AddEntry(hists[5],"Physics run 213539", "l")

hists[1].SetTitle("")
hists[1].GetXaxis().SetRangeUser(-200,200)
hists[1].SetStats(0)
hists[1].GetXaxis().SetLabelSize(0.03)
hists[1].GetXaxis().SetTitleSize(0.04)
hists[1].GetXaxis().SetTitleOffset(1)
canvas = ROOT.TCanvas("canvas","canvas",1000,800)
hists[1].Draw()
hists[2].Draw("same")
hists[5].Draw("same")
legend.Draw("same")
canvas.Print("TypicalZDistributions.pdf")



