# -*- coding: utf-8 -*-
#! /usr/bin/env python
import ROOT
from ROOT import TFile, TH1F, TCanvas, TLegend, TGraphErrors
from array import array

ROOT.gROOT.SetBatch(True)

ntrkcuts = ["3","4","5"]
markers = [24,25,26]
vdmscans = ["201351","207216","207219","214984","215021"]
colors = [ROOT.kRed,ROOT.kMagenta+3,ROOT.kMagenta+1,ROOT.kCyan+3,ROOT.kBlue+2]
fakesfile = TFile("/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/Data/FakeCorrection/mc_8TeV_17.2_VtxLumi_BothSamples/fakerates.root","READ")
tgs_pileup = []
tgs_pileup_april = []
tgs_pileup_july = []
for vdmscan,color in zip(vdmscans,colors):
  print "VdMScan: ", vdmscan
  maskingfile = TFile("/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/Data/MaskingCorrection/data_8TeV_17.2_VtxLumi_"+vdmscan+"/pmask_cache.root","READ")
  for ntrk,marker in zip(ntrkcuts,markers):
    print "Track cut:", ntrk
    tg_masking_correction = maskingfile.Get("tg_masking_correction_NTrk"+ntrk)
    npoints = tg_masking_correction.GetN()
    tg_pileup_correction = TGraphErrors(npoints)
    print "Masking correction plot has", npoints, "points"
    tg_fake_correction = fakesfile.Get("tg_fake_fraction_MuReconMC_NTrk"+ntrk)
    print "Fake correction plot has", tg_fake_correction.GetN(), "points"
    for i in range(npoints):
      murec = tg_masking_correction.GetX()[i]
      imcf = tg_masking_correction.Eval(murec)
      mufake = tg_fake_correction.Eval(murec*imcf)
      mureal = murec-mufake
      fmcf = tg_masking_correction.Eval(mureal)
      muvis = mureal*fmcf
      tg_pileup_correction.SetPoint(i,murec,muvis)
      tg_pileup_correction.SetTitle(vdmscan+", NTrkCut"+ntrk)
      tg_pileup_correction.SetName(vdmscan+"_NTrk"+ntrk)
      tg_pileup_correction.SetMarkerStyle(marker)
      tg_pileup_correction.SetMarkerColor(color)
      tg_pileup_correction.GetXaxis().SetTitle("#mu_{rec}")
      tg_pileup_correction.GetYaxis().SetTitle("#mu_{vis}")
      tg_pileup_correction.SetMarkerSize(0.8)
    tgs_pileup.append(tg_pileup_correction)
    if( vdmscan == "201351" ):
      tgs_pileup_april.append(tg_pileup_correction)
    if( vdmscan == "207216" ):
      tgs_pileup_july.append(tg_pileup_correction)
print len(tgs_pileup)

tgs_file = TFile("VertexCounting_PileupCorrectionPlots.root","recreate")
for tg in tgs_pileup:
  tg.Write()
tgs_file.Close()

legend = ROOT.TLegend(0.7,0.15,0.85,0.45)
legend.SetFillColor(0)

colors_april = [ROOT.kRed, ROOT.kRed-7, ROOT.kRed+2]
for tg,color,ntrk in zip(tgs_pileup_april,colors_april,ntrkcuts):
  tg.SetTitle("")
  tg.SetMarkerColor(color)
  tg.SetFillColor(color)
  tg.SetLineColor(0)
  tg.SetMarkerStyle(8)
  legend.AddEntry(tg,"April, NTrk"+ntrk,"PF")

colors_july = [ROOT.kBlue, ROOT.kBlue-7, ROOT.kBlue+2]
for tg,color,ntrk in zip(tgs_pileup_july,colors_july,ntrkcuts):
  tg.SetTitle("")
  tg.SetMarkerColor(color)
  tg.SetFillColor(color)
  tg.SetLineColor(0)
  tg.SetMarkerStyle(8)
  legend.AddEntry(tg,"July, NTrk"+ntrk,"PF")

canvas = ROOT.TCanvas("canvas","canvas",1200,800) 
tgs_pileup_july[0].Draw("AP")
for tg in tgs_pileup_july:
  tg.Draw("P")
for tg in tgs_pileup_april:
  tg.Draw("P")
legend.Draw("same")
canvas.Print("OverallCorrection_AprilAndJuly.pdf")



  
#  g_murecvsmu.SetMarkerStyle(8)
#  g_murecvsmu.SetMarkerSize(1)
#  g_murecvsmu.SetMarkerColor(ROOT.kRed)
#  g_murecvsmu.SetLineColor(ROOT.kRed)
#  g_murecvsmu.SetTitle("For NTrkCut "+ntrk)
#  g_murecvsmu.GetXaxis().SetTitle("#mu")
#  g_murecvsmu.GetXaxis().SetTitleOffset(0.8)
#  g_murecvsmu.GetYaxis().SetTitle("#mu_{rec}")
#  g_murecvsmu.GetYaxis().SetTitleOffset(0.8)#

#  c1 = TCanvas("c1","c1",1200,800)
#  g_murecvsmu.Draw("AP")
#  h_murecvsmu.Draw("col same")
#  c1.Print("mu_rec_vs_mu_NTrk"+ntrk+".pdf")


