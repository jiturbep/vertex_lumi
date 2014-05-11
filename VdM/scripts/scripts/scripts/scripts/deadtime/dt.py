#! /usr/bin/python

import ROOT
import math

scanNtup = "/eliza18/atlas/dryu/GBStudies/results/ExtendedTrackingZ/MB_data7TeV_cfg-TrkMinBias_VtxStartup150_z1000/182013/2011MayScanRaw-v8.root" 
run = 182013

def readPseudoLbTimes():
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
        # if acq == 1 and ip == 1:
        times.append((plb,(start/1E9,stop/1E9)))
    return times

times = readPseudoLbTimes()

outFile = ROOT.TFile("VdMCheck.root","update")
ROOT.gROOT.cd()
fineRateHistAP = outFile.Get("FineTriggerRateAP")
fineRateHistAV = outFile.Get("FineTriggerRateAV")
fineEventRateHist = outFile.Get("FineEventRate")

fineRateHistAP.Sumw2()
fineRateHistAV.Sumw2()
fineEventRateHist.Sumw2()
# convert to kHz for convenience
fineRateHistAP.Scale(1./1000.)
fineRateHistAV.Scale(1./1000.)
fineEventRateHist.Scale(1./1000.)

c1 = ROOT.TCanvas("c1","c1",600,600)
c1.cd()
ROOT.gPad.SetTickx(1)
ROOT.gPad.SetTicky(1)
ROOT.gPad.SetGridx(1)
ROOT.gPad.SetGridy(1)
ROOT.gPad.SetLeftMargin(0.17)
ROOT.gPad.SetRightMargin(0.02)
ROOT.gPad.SetTopMargin(0.01)
fineRateHistAP.SetStats(0)
fineRateHistAP.SetLineWidth(2)
fineRateHistAP.SetMarkerStyle(20)
fineRateHistAP.SetMarkerSize(1.2)
fineRateHistAP.Draw("ep")

fineRateHistAV.SetMarkerStyle(22)
fineRateHistAV.SetMarkerSize(0.9)
fineRateHistAV.SetLineWidth(2)
fineRateHistAV.SetMarkerColor(ROOT.TColor.kRed)
fineRateHistAV.SetLineColor(fineRateHistAV.GetMarkerColor())
fineRateHistAV.Draw("ep same")

fineEventRateHist.SetMarkerStyle(20)
fineEventRateHist.SetMarkerSize(0.6)
fineEventRateHist.SetLineWidth(2)
fineEventRateHist.SetMarkerColor(ROOT.TColor.kAzure)
fineEventRateHist.SetLineColor(fineEventRateHist.GetMarkerColor())
fineEventRateHist.Draw("ep same")

l = ROOT.TLatex()
l.SetTextSize(0.05) 
l.SetNDC()
l.SetTextFont(72)
l.SetTextColor(1)
l.DrawLatex(0.36,0.93,"ATLAS");

p = ROOT.TLatex()
p.SetTextSize(0.05) 
p.SetNDC()
p.SetTextFont(42)
p.SetTextColor(1)
p.DrawLatex(0.54,0.93,"Preliminary");

#fineRateHistAP.GetXaxis().SetTimeFormat("%Ss %F1970-01-01 00:00:00");
fineRateHistAP.GetXaxis().SetTimeDisplay(1)
fineRateHistAP.GetYaxis().SetTitle("Rate ( kHz )")
fineRateHistAP.GetXaxis().SetTitle("Time (since xxx 2011 xxhxxm00s)")

fineRateHistAP.GetXaxis().SetTickLength(0.022)
fineRateHistAP.GetXaxis().SetLabelSize(0.045)
fineRateHistAP.GetYaxis().SetLabelSize(0.05)
fineRateHistAP.GetXaxis().SetTitleSize(0.045)
fineRateHistAP.GetXaxis().SetTitleOffset(1.3)
fineRateHistAP.GetYaxis().SetTitleSize(0.045)
fineRateHistAP.GetYaxis().SetTitleOffset(1.5)
fineRateHistAP.SetMaximum(fineRateHistAP.GetMaximum()*1.15)

leg2 = ROOT.TLegend(0.45,0.15,0.8,0.5)
leg2.SetFillStyle(0)
leg2.SetBorderSize(0)
leg2.AddEntry(fineRateHistAP,"TAP","pl")
leg2.AddEntry(fineRateHistAV,"TAV","pl")
leg2.AddEntry(fineEventRateHist,"Event rate","pl")
leg2.Draw("same")

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

liveHist = ROOT.TH1F("Livefractions", ";pLB;Live fraction", n_pLBs, pLBmin - 0.5, pLBmax + 0.5)
apHist = ROOT.TH1F("AfterPrescaleVdmHist",";pLB;After prescale rate ( kHz )", n_pLBs, pLBmin - 0.5, pLBmax + 0.5)
evHist = ROOT.TH1F("EventRateVdmHist",";pLB;Event rate ( kHz )", n_pLBs, pLBmin - 0.5, pLBmax + 0.5)

# now loop over pLB's and calc deadtime per pLB
histCounter = 1
for plb,pair in times:
    lowerBin = fineEventRateHist.GetXaxis().FindFixBin(pair[0])
    upperBin = fineEventRateHist.GetXaxis().FindFixBin(pair[1])
    # first get the average after-prescale rate in this pLB
    averageRate = 0.
    weight = 0.    
    for bin in range(lowerBin,upperBin+1):
        ap = fineRateHistAP.GetBinContent(bin)
        if ap > 1./1000.:
            averageRate += ap
            weight += 1.
    averageRate /= weight
    # now loop again for the deadtime, take the average rate from above in case the ap/av rates have not been published
    totalAp = 0.
    totalEv = 0.
    weight = 0.    
    for bin in range(lowerBin+1,upperBin):
        ev = fineEventRateHist.GetBinContent(bin)
        totalAp += averageRate
        totalEv += ev
        weight += 1.
    liveFraction = 0.
    if totalAp > 0.:
        liveFraction = totalEv / totalAp
    else:
        print "Error: no after-prescale rate for pLB",plb

    liveHist.SetBinContent(histCounter,liveFraction)
    apHist.SetBinContent(histCounter,averageRate)
    evHist.SetBinContent(histCounter,totalEv/weight)
    histCounter += 1

ROOT.gPad.Update()
ROOT.gPad.Modified()
outFile.cd()
liveHist.Write("",ROOT.TObject.kOverwrite)
apHist.Write("",ROOT.TObject.kOverwrite)
evHist.Write("",ROOT.TObject.kOverwrite)

