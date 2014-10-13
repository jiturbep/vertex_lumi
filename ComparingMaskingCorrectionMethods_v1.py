# -*- coding: utf-8 -*-
#! /usr/bin/env python
import ROOT
import math

ROOT.TH1.AddDirectory(False)
ROOT.gROOT.SetBatch()

path = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/Data/MaskingCorrection/mc_8TeV_17.2_VtxLumi_2newsets/"

rf_Gaus = ROOT.TFile(path+"pmask_cache_2to10and1to10.root","READ")
rf_1to10PFTemp = ROOT.TFile(path+"pmask_cache_2to10and1to10FitWithTemplate.root","READ")
rf_1to10NPFTemp = ROOT.TFile(path+"pmask_cache_2to10and1to10NoPoissonFactorFitWithTemplate.root","READ")

pmask_Gaus = rf_Gaus.Get("h_pmask_dz_NTrk5")
pmask_PF_Temp = rf_1to10PFTemp.Get("h_pmask_dz_NTrk5")
pmask_NPF_Temp = rf_1to10NPFTemp.Get("h_pmask_dz_NTrk5") 

print pmask_PF_Temp.GetNbinsX(), pmask_Gaus.GetNbinsX()

tg_diff = ROOT.TGraph(pmask_PF_Temp.GetNbinsX())
for i in range(pmask_PF_Temp.GetNbinsX()+1):
	dzx = pmask_PF_Temp.GetBinCenter(i+1)
	dzy = pmask_NPF_Temp.GetBinCenter(i+1)
	if (dzx != dzy): print "dzs are different!" 
	x = pmask_PF_Temp.GetBinContent(i+1)
	y = pmask_NPF_Temp.GetBinContent(i+1)
	if (y>0):
		diff = ((x/y)-1)*100
	else: diff = 0
	tg_diff.SetPoint(i,dzx,diff)
canvas = ROOT.TCanvas("diff","diff",1200,800)
tg_diff.GetXaxis().SetRangeUser(-200,200)
tg_diff.SetMaximum(10)
tg_diff.SetMinimum(-10)
tg_diff.SetMarkerStyle(8)
tg_diff.SetMarkerSize(0.8)
tg_diff.Draw("AP")
canvas.Print("diff.pdf")






