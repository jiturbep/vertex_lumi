#! /usr/bin/python

import ROOT
import math
scanNtup = "/eliza18/atlas/dryu/GBStudies/results/ExtendedTrackingZ/MB_data7TeV_cfg-TrkMinBias_VtxStartup150_z1000/182013/2011MayScanRaw-v8.root"
tagNtup = [
"/eliza18/atlas/dryu/tmp/182013/data11_7TeV.00182013.calibration_VdM.merge.TAG.x123_m841_m848/data11_7TeV.00182013.calibration_VdM.merge.TAG.x123_m841_m848._0001.1", 
"/eliza18/atlas/dryu/tmp/182013/data11_7TeV.00182013.calibration_VdM.merge.TAG.x123_m841_m848/data11_7TeV.00182013.calibration_VdM.merge.TAG.x123_m841_m848._0002.1", 
"/eliza18/atlas/dryu/tmp/182013/data11_7TeV.00182013.calibration_VdM.merge.TAG.x123_m841_m848/data11_7TeV.00182013.calibration_VdM.merge.TAG.x123_m841_m848._0003.1", 
"/eliza18/atlas/dryu/tmp/182013/data11_7TeV.00182013.calibration_VdM.merge.TAG.x123_m841_m848/data11_7TeV.00182013.calibration_VdM.merge.TAG.x123_m841_m848._0004.1", 
"/eliza18/atlas/dryu/tmp/182013/data11_7TeV.00182013.calibration_VdM.merge.TAG.x123_m841_m848/data11_7TeV.00182013.calibration_VdM.merge.TAG.x123_m841_m848._0005.1", 
"/eliza18/atlas/dryu/tmp/182013/data11_7TeV.00182013.calibration_VdM.merge.TAG.x123_m841_m848/data11_7TeV.00182013.calibration_VdM.merge.TAG.x123_m841_m848._0006.1", 
"/eliza18/atlas/dryu/tmp/182013/data11_7TeV.00182013.calibration_VdM.merge.TAG.x123_m841_m848/data11_7TeV.00182013.calibration_VdM.merge.TAG.x123_m841_m848._0007.1"]
run = 182013

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

def GetPseudoLbTimes():
    file = ROOT.TFile(scanNtup,'read')
    ROOT.gROOT.cd()
    tree = file.Get('vdMScanData')

    times = [] 
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
        times.append((plb,(start/1E9,stop/1E9)))
    return times

#if acq == 1 and ip == 1:
print "Run",run
print "VdM ntuple",scanNtup
print "TAG Files",tagNtup

print "Lumi block start stop times"
lbStartTimeGraph = GetLbStartGraph(run)

print "Pseudo LB time stamps"
pLBTimes = GetPseudoLbTimes()


print "High-resolution TAG event rates"
firstTime = 0
lastTime = 0
timeList = lbStartTimeGraph.GetY()
# find begin and end timestamp for histogram
for index in range(lbStartTimeGraph.GetN()):
    now = timeList[index]
    if now>0:
        if not firstTime:
            firstTime = math.floor(now)
        lastTime = math.ceil(now)
# hist for fine resolution event rate
fineRateHist = ROOT.TH1F("FineEventRate",";Time ( s );Rate ( Hz )",int(lastTime-firstTime),firstTime,lastTime)
for fName in tagNtup:
   # open the TAG file and fill events in 1 second bins
   file = ROOT.TFile(fName,'read')
   ROOT.gROOT.cd()
   tree = file.Get('POOLCollectionTree')
   br1 = tree.GetBranch('L1PassedTrigMaskTAV3') # '3' bc below is item 109, 109/32 = 3...
   br2 = tree.GetBranch('EventTime')
   br3 = tree.GetBranch('EventTimeNanoSec')
   br4 = tree.GetBranch('BunchId')
   br5 = tree.GetBranch('LumiBlockN')
   selectItem = int(109) # 109 is L1_MBTS_2_BGRP7
   for event in xrange(tree.GetEntries()):
       if event%10000 == 0:
           print "Event",event,"/",tree.GetEntries()
       br1.GetEntry(event)
       if (int(br1.GetLeaf('L1PassedTrigMaskTAV3').GetValue()) >> (selectItem%32)) & 0x1:
           br2.GetEntry(event)
           br3.GetEntry(event)
           br4.GetEntry(event)
           br5.GetEntry(event)
           lb = int(br5.GetLeaf('LumiBlockN').GetValue())
   
           eventTime = float(br2.GetLeaf('EventTime').GetValue()) + float(br3.GetLeaf('EventTimeNanoSec').GetValue())*1e-9
           validEvent = False
           for plb,pair in pLBTimes:
               if eventTime >= pair[0] and eventTime <= pair[1]:
                   validEvent = True
           if validEvent:
               # fill events in hist with 1 second binning, to assess fine-grained rates, and make dt correction later on
               fineRateHist.Fill(eventTime)
   
fineRateHist.Draw("hist")
outFile = ROOT.TFile("VdMCheck.root","update")
fineRateHist.Write("",ROOT.TObject.kOverwrite)
