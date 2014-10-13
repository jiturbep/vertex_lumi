# -*- coding: utf-8 -*-
#! /usr/bin/env python
import ROOT
import math 

ntrkcuts = ["3","4","5"]
colours = [ROOT.kRed, ROOT.kBlue, ROOT.kBlack]
rootfile = ROOT.TFile("muvis_ngenint.root","READ")
gs_res = []
legend = ROOT.TLegend(0.8,0.8,0.9,0.9)
legend.SetFillColor(0)
for ntrkcut,colour in zip(ntrkcuts,colours):
	print "////////////////// NTrkCut ", ntrkcut, "//////////////////"
	graph = rootfile.Get("tg_muvis_ngenint_NTrk"+ntrkcut)
	print "Number of points in the graph:", graph.GetN()
	fit = graph.GetFunction("line0")
	graph_res = ROOT.TGraph(0)
	count = 0
	for i in range(0,graph.GetN()) :
		x,y,xerr,yerr = ROOT.Double(0),ROOT.Double(0),ROOT.Double(0),ROOT.Double(0)
		graph.GetPoint(i,x,y)
		#print i,x,y
		yerr = graph.GetErrorY(i)
		if (yerr > 0):
			res = (y - fit.Eval(x))/yerr
			count +=1
			#print "Residual:", res
		else: res = 0
		graph_res.SetPoint(i, x, res)
	print "Number of points with residuals:", count
	graph_res.SetMarkerColor(colour)
	legend.AddEntry(graph_res, "NTrkCut"+ntrkcut, "P")
	gs_res.append(graph_res)
for graph_res in gs_res:
	graph_res.SetTitle("")
	graph_res.GetXaxis().SetTitle("Generated Interactions")
	graph_res.GetYaxis().SetTitle("data-fit/#sigma_{data}")
	graph_res.SetMarkerStyle(8)
	graph_res.SetMarkerSize(0.8)
	graph_res.GetXaxis().SetRangeUser(0,30)
	graph_res.SetMaximum(10)
	graph_res.SetMinimum(-10)
canvas = ROOT.TCanvas("residuals","residuals",800,600)
Pad = canvas.cd(1)
Pad.SetTicks(1,1)
gs_res[0].Draw("AP")
gs_res[1].Draw("P")
gs_res[2].Draw("P")
legend.Draw("same")
canvas.Print("FitResiduals_AllTrkCuts.pdf")
