# -*- coding: utf-8 -*-
#! /usr/bin/env python
import ROOT

rfxingb = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/MC_ClosureTest/17.2-VtxLumi/InDetTrackD3PD_results_actualXing_both.root"
rfxing15 = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/MC_ClosureTest/17.2-VtxLumi/InDetTrackD3PD_results_actualXing_15.root"
rfngenb = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/MC_ClosureTest/17.2-VtxLumi/InDetTrackD3PD_results_NGenInt_both.root"
rfngen15 = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/VertexCounts/MC_ClosureTest/17.2-VtxLumi/InDetTrackD3PD_results_NGenInt_15.root"

f_rfxingb = ROOT.TFile(rfxingb,"READ")
f_rfxing15 = ROOT.TFile(rfxing15,"READ")
f_rfngenb = ROOT.TFile(rfngenb,"READ")
f_rfngen15 = ROOT.TFile(rfngen15,"READ")

hzplbxingb = f_rfxingb.Get("hist/PriVtxZpLB_BCID0")
hzplbxing15 = f_rfxing15.Get("hist/PriVtxZpLB_BCID0")
hzplbngenb = f_rfngenb.Get("hist/PriVtxZpLB_BCID0")
hzplbngen15 = f_rfngen15.Get("hist/PriVtxZpLB_BCID0")

hzxingb = hzplbxingb.ProjectionX()
hzxingb.SetName("actualIntPerXing_both")
hzxing15 = hzplbxing15.ProjectionX()
hzxing15.SetName("actualIntPerXing_15")
hzngenb = hzplbngenb.ProjectionX()
hzngenb.SetName("NGenInt_both")
hzngen15 = hzplbngen15.ProjectionX()
hzngen15.SetName("NGenInt_15")

h_ratio_b = hzxingb.Clone()
h_ratio_b.Divide(hzngenb)
h_ratio_b.SetName("Ratio_b")

h_ratio_15 = hzxing15.Clone()
h_ratio_15.Divide(hzngen15)
h_ratio_15.SetName("Ratio_15")

rf = ROOT.TFile("projections.root","RECREATE")
hzxingb.Write()
hzxing15.Write()
hzngenb.Write()
hzngen15.Write()
h_ratio_b.Write()
h_ratio_15.Write()
rf.Close()



