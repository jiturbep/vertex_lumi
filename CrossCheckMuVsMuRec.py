# -*- coding: utf-8 -*-
#! /usr/bin/env python
import ROOT
from ROOT import TFile, TH1F, TCanvas, TLegend

ROOT.gROOT.SetBatch(True)

ntrkcuts = ["3","4","5"]

fakesfile = TFile("/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/Data/FakeCorrection/mc_8TeV_17.2_VtxLumi_BothSamples/fakerates.root","READ")
countsfile = TFile("/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VtxTruthMatchResults/mc_8TeV_17.2_VtxLumi_BothSamples/InDetTrackD3PD_results.root","READ")

h_actualInt_NGenInt = countsfile.Get("hist/h_actualInt_NGenInt")
h_NTrig_actualInt = h_actualInt_NGenInt.ProjectionX()

for ntrk in ntrkcuts:
  h_all_actualInt_NTrk = countsfile.Get("hist/h_all_actualInt_NTrk"+ntrk)
  h_murecvsmu = h_all_actualInt_NTrk.Clone()
  h_murecvsmu.Divide(h_NTrig_actualInt)
  g_murecvsmu = fakesfile.Get("tg_mu_recon_mu_NTrk"+ntrk)
  
  g_murecvsmu.SetMarkerStyle(8)
  g_murecvsmu.SetMarkerSize(1)
  g_murecvsmu.SetMarkerColor(ROOT.kRed)
  g_murecvsmu.SetLineColor(ROOT.kRed)
  g_murecvsmu.SetTitle("For NTrkCut "+ntrk)
  g_murecvsmu.GetXaxis().SetTitle("#mu")
  g_murecvsmu.GetXaxis().SetTitleOffset(0.8)
  g_murecvsmu.GetYaxis().SetTitle("#mu_{rec}")
  g_murecvsmu.GetYaxis().SetTitleOffset(0.8)

  c1 = TCanvas("c1","c1",1200,800)
  g_murecvsmu.Draw("AP")
  h_murecvsmu.Draw("col same")
  c1.Print("mu_rec_vs_mu_NTrk"+ntrk+".pdf")


