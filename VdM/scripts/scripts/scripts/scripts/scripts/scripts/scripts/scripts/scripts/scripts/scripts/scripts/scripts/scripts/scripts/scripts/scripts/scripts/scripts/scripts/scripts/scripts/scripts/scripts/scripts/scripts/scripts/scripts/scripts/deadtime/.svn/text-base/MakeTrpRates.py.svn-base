#! /usr/bin/python

import ROOT
import math

def GetLbStartGraph(run):
    import ROOT
    import CoolHelper    
    lbTimes = CoolHelper.getLbTimes(run)
    lbGraph = ROOT.TGraph(0)
    lbGraph.SetName("LbStarttimesGraph")
    graphCounter = 0
    for lb in lbTimes.keys():
        start = float(lbTimes[lb][0])/1e9
        lbGraph.Set(graphCounter+1)
        lbGraph.SetPoint(graphCounter,lb,start)
        graphCounter += 1
    return lbGraph

trpNtup = "TriggerRates_ATLAS_182013.root"
scanNtup = "/eliza18/atlas/dryu/GBStudies/results/ExtendedTrackingZ/MB_data7TeV_cfg-TrkMinBias_VtxStartup150_z1000/182013/2011MayScanRaw-v8.root"
run = 182013

lbStartTimeGraph = GetLbStartGraph(run)

file = ROOT.TFile(scanNtup,'read')
ROOT.gROOT.cd()
tree = file.Get('vdMScanData')

times = [] 
allTimes = [] 
br1 = tree.GetBranch('SCANDATA')
for event in range(tree.GetEntries()):
    br1.GetEntry(event)

    acq = br1.GetLeaf('AcquisitionFlag').GetValue()
    start = float(br1.GetLeaf('StartTime').GetValue())
    stop = float(br1.GetLeaf('EndTime').GetValue())
    ip = br1.GetLeaf('ScanningIP').GetValue()
    plb = br1.GetLeaf('ScanLB').GetValue()
    moving = br1.GetLeaf('MovingBeam').GetValue()
    plane = br1.GetLeaf('ScanInPlane').GetValue()
    nomSep = br1.GetLeaf('NominalSeparation').GetValue()
    allTimes.append((plb,(start/1E9,stop/1E9),(acq,ip)))
    times.append((plb,(start/1E9,stop/1E9)))
	#if acq == 1 and ip == 1:

#Get max and min pLB
pLBmin = 10000
pLBmax = 0
for plb,pair in times:
	if plb > pLBmin:
		pLBmin = plb
	if plb < pLBmax:
		pLBmax = plb
print "pLBmin = " + str(pLBmin) + " / pLBmax = " + str(pLBmax)
n_pLBs = pLBmax - pLBmin + 1

firstTime = 0
lastTime = 0
timeList = lbStartTimeGraph.GetY()
for index in range(lbStartTimeGraph.GetN()):
    now = timeList[index]
    if now>0:
        if not firstTime:
            firstTime = math.floor(now)
        lastTime = math.ceil(now)
fineRateHistAP = ROOT.TH1F("FineTriggerRateAP",";Time ( s );Rate ( Hz )",int(lastTime-firstTime),firstTime,lastTime)
fineRateHistAV = ROOT.TH1F("FineTriggerRateAV",";Time ( s );Rate ( Hz )",int(lastTime-firstTime),firstTime,lastTime)

file = ROOT.TFile(trpNtup,'read')
ROOT.gROOT.cd()
tree = file.Get('L1_Rate')
br1 = tree.GetBranch('L1_MBTS_2_BGRP7_TAP')
br2 = tree.GetBranch('L1_MBTS_2_BGRP7_TAV')
br3 = tree.GetBranch('TimeStamp')
br4 = tree.GetBranch('L1_MBTS_2_BGRP7_PS')
for event in xrange(tree.GetEntries()):
    if event%10000 == 0:
        print "Event",event,"/",tree.GetEntries()
    br1.GetEntry(event)
    br2.GetEntry(event)
    br3.GetEntry(event)
    ap = float(br1.GetLeaf('L1_MBTS_2_BGRP7_TAP').GetValue())
    av = float(br2.GetLeaf('L1_MBTS_2_BGRP7_TAV').GetValue())
    eventTime = float(br3.GetLeaf('TimeStamp').GetValue())

    histCounter = 0
    for plb,pair in times:
        histCounter += 1
        if eventTime >= pair[0] and (eventTime-1.0) <= pair[1]:
            # fill events in hist with 1 second binning, to assess fine-grained rates, and make dt correction later on
            fineRateHistAP.Fill(eventTime-0.5,ap)
            fineRateHistAV.Fill(eventTime-0.5,av)

h_ps = ROOT.TH1F("h_ps", "h_ps", n_PLBs, pLB_min - 0.5, pLB_max + 0.5)
for plb,pair in times:
	for event in xrange(tree.GetEntries()):
		br3.GetEntry(event)
		br4.GetEntry(event)
		eventTime = float(br3.GetLeaf('TimeStamp').GetValue())
		
		if eventTime < pair[0] or (eventTime-1.0) > pair[1]:
			continue
		else:
			ps = int(br4.GetLeaf('L1_MBTS_2_BGRP7_PS').GetValue())
			h_ps.Fill(plb, ps)
			break


c = ROOT.TCanvas("c","c",600,400)
c.Divide(2,1)
c.cd(1)
fineRateHistAP.Draw("hist")
c.cd(2)
fineRateHistAV.Draw("hist")
outFile = ROOT.TFile("VdMCheck.root","update")
fineRateHistAP.Write("",ROOT.TObject.kOverwrite)
fineRateHistAV.Write("",ROOT.TObject.kOverwrite)
h_ps.Write()
