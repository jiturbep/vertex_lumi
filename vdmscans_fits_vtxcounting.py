# -*- coding: utf-8 -*-
#! /usr/bin/env python
from ROOT import *
import math
from itertools import izip
import os
import sys

gROOT.SetBatch(True)

vdmsets = ["April","July","November"]
#vdmsets = ["November"]
coordinates = ["x", "y"]
ntrks = [3,4,5]
ntrks = [5]
colors = [kRed,kBlue,kBlack]
for vdmset in vdmsets:
	fitfiles = []
	if vdmset == "April":
		vdmscans = [201351]
	elif vdmset == "July":
		vdmscans = [207216,207219]
	elif vdmset == "November":
		vdmscans = [214984,215021]
	for vdmscan in vdmscans:
		rootfile = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/Data/VdMCalibration/VdMScan-"+str(vdmscan)+"/17.2-VtxLumi/NVtx/vdm_results.root"
		for coord in coordinates:
			if vdmscan == 201351:
				scans = [1,2,3]
				BCIDs=[1,241,2881,3121]
				fitfunc = "dgc"
			if vdmscan == 207216:
				scans = [4,5,6]
				BCIDs=[1,721,1821]
				fitfunc = "dgc"
			if vdmscan == 207219:
				scans = [8]
				BCIDs=[1,721,1821]
				fitfunc = "dgc"
			if vdmscan == 214984:
				scans = [10,11,14]
				BCIDs=[1,2361,2881]
				fitfunc = "sgc"
			if vdmscan == 215021:
				scans = [15]
				BCIDs=[1,2361,2881]
				fitfunc = "sgc"
			f_scans = TFile(rootfile, "READ")
			for index,bcid in enumerate(BCIDs):
				print bcid, index
				for scan in scans:
					fits = []
					graphs = []
					graphs_res = []
					legend = TLegend(0.75,0.6,0.9,0.9)
					legend.SetFillColor(kWhite)
					legend.SetTextSize(0.018)
					topvalues = []
					bottomvalues = []
					for ntrk,color in zip(ntrks,colors):
						graph = f_scans.Get( "tg_"+coord+"_BCID"+str(bcid)+"_scan"+str(scan)+"_NTrkCut"+str(ntrk) )
						print 'Initial graph, number of points', graph.GetN()
						fit = graph.GetFunction("f_"+coord+"_"+fitfunc)
						fit.SetLineColor(color)
						fit.SetLineStyle(2)
						fits.append(fit)
						graph_new = TGraphErrors(0)
						for i in range(0,graph.GetN()) :
							x,y,xerr,yerr= Double(0),Double(0),Double(0),Double(0)
							graph.GetPoint(i,x,y)
							yerr = graph.GetErrorY(i)
							xerr = graph.GetErrorX(i)
							graph_new.SetPoint(i,x,y)
							graph_new.SetPointError(i,xerr,yerr)
						graph_new.SetMarkerStyle(8)
						#graph_new.SetMarkerSize(0.4)
						graph_new.SetMarkerColor(color)
						graph_new.SetTitle(str(coord)+"-"+str(bcid)+"-scan"+str(scan))
						graph_new.GetXaxis().SetLabelOffset(2.)
						graph_new.GetYaxis().SetTitle("\mu_{vis,sp}")
						graphs.append(graph_new)
						legend.AddEntry(graph_new,"NTrkCut "+str(ntrk),"P")
						#Residual Plot
						graph_res = TGraph(0)
						for i in range(0,graph.GetN()) :
							x,y,xerr,yerr= Double(0),Double(0),Double(0),Double(0)
							graph.GetPoint(i,x,y)
							yerr = graph.GetErrorY(i)
							res = (y - fit.Eval(x))/yerr
							graph_res.SetPoint(i,x,res)
						graph_res.SetTitle("")
						graph_res.GetXaxis().SetTitle("#Delta"+str(coord)+" (mm)")
						graph_res.GetXaxis().SetLabelSize(0.03)
						graph_res.GetYaxis().SetTitleOffset(1.2)
						graph_res.GetYaxis().SetLabelSize(0.03)
						graph_res.GetYaxis().SetTitle("data-fit/#sigma_{data}")
						graph_res.SetMarkerStyle(8)
						#graph_res.SetMarkerSize(0.4)
						graph_res.SetMarkerColor(color)
						graphs_res.append(graph_res)
						topvalues.append(graph_res.GetYaxis().GetXmax())
						bottomvalues.append(graph_res.GetYaxis().GetXmin())
					top = max(topvalues)
					bottom = min(bottomvalues)
					#Drawing fitting plots and residual plots
					FIGURE2_RATIO = 0.35
					SUBFIGURE_MARGIN = 0.15
					canvas = TCanvas( "canvas", "canvas", 1000, 1000 )
					canvas.cd()
					graphs[0].Draw("AP")
					for i in range(len(graphs)):
						graphs[i].Draw("P")
						fits[i].Draw("same")
					legend.Draw("same")
					canvas.SetBottomMargin(FIGURE2_RATIO)
					Pad = TPad( "p_test", "p_test", 0, 0, 1, 1.0-SUBFIGURE_MARGIN, 0, 0, 0) # create new pad, fullsize to have equal font-sizes in both plots                                   
					Pad.SetTopMargin(1-FIGURE2_RATIO)  # top-boundary (should be 1 - thePad->GetBottomMargin() )                                 
					Pad.SetFillStyle(0)  # needs to be transparent                                                                             
					Pad.Draw()
					Pad.cd()
					line = TLine(graph.GetXaxis().GetXmin(),0.,graph.GetXaxis().GetXmax(),0.)
					line.SetLineColor(kRed)
					line.SetLineStyle(2)
					graphs_res[0].SetMaximum(top)
					graphs_res[0].SetMinimum(bottom)
					graphs_res[0].Draw("AP")
					for i in range(len(graphs_res)):
						graphs_res[i].Draw("P")
					line.Draw("same")
					Pad.SetTicks(1,1)
					canvas.Print("FittedVdMScans/"+vdmset+"/"+str(scan)+"_"+str(bcid)+"_"+coord+".pdf")
	
	if vdmset == "April":
		scans = [1,2,3]
		BCIDs=[1,241,2881,3121]
	if vdmset == "July":
		scans = [4,5,6,8]
		BCIDs=[1,721,1821]
	if vdmset == "November":
		scans = [10,11,14,15]
		BCIDs=[1,2361,2881]
	for scan in scans:
		for bcid in BCIDs:
			for coord in coordinates:
				fitfiles.append(str(scan)+"_"+str(bcid)+"_"+coord+".pdf")
	
	latfile = open("FittedVdMScans/"+vdmset+"/"+vdmset+"_FittedScanCurves.tex", "w")
	latfile.write(r"\documentclass[8pt]{beamer}""\n")
	latfile.write(r"\usepackage{graphicx}""\n")
	latfile.write(r"\usepackage{xcolor}""\n")
	latfile.write(r"\usetheme{default}""\n")
	latfile.write(r"\begin{document}""\n")
	
	latfile.write(r"\begin{frame}""\n")
	latfile.write(r"\centering""\n")
	latfile.write(r"\large{Fitted "+vdmset+" Scan Curves}""\n")
	latfile.write(r"\end{frame}""\n")

	for f in fitfiles:
		latfile.write(r"\begin{frame}""\n")
		latfile.write(r"\center""\n")
		latfile.write(r"\includegraphics[scale=0.5]{"+f+"}""\n")
		latfile.write("\end{frame}""\n")

	latfile.write(r"\end{document}""\n")	
	latfile.close()

for vdmset in vdmsets:
	os.system("cd FittedVdMScans/"+vdmset+";pdflatex "+vdmset+"_FittedScanCurves.tex;cd -")
	os.system("scp FittedVdMScans/"+vdmset+"/"+vdmset+"_FittedScanCurves.pdf julia@pc2012.hep.manchester.ac.uk:~/WWW")
