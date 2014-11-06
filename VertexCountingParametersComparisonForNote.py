#! /usr/bin/env python
import ROOT
import math
from itertools import izip

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

scans = [1, 2, 3, 4, 5, 6, 8, 10, 11, 14, 15]
ntrkcuts = [3,4,5]
ntrkcut = 5

alg = 'vtx_counting'
bcids = []
sigma = []
sigma_err = []
lumi_sp = []
lumi_sp_err = []
capsigma_x = []
capsigma_x_err = []
capsigma_y = []
capsigma_y_err = []
mumax_x = []
mumax_x_err = []
mumax_y = []
mumax_y_err = []
chi2_x = []
chi2_y = []
for scan in scans:
	file1 = "InfoFromVdMScans/vtx_counting/results_parameters_scan"+str(scan)+"_NTrkCut"+str(ntrkcut)+".txt"
	for line in open(file1, 'r'):
		bcid_n,sigma_v,sigma_v_err,Sigma_x,Sigma_x_err,Sigma_y,Sigma_y_err,Peak_mu_x,Peak_mu_x_err,Peak_mu_y,Peak_mu_y_err,chi2x,chi2y,l_sp,l_sp_err,r_x_i,r_y_i,c_x,c_x_err,c_y,c_y_err = line.split()
		bcids.append((bcid_n))
		sigma.append(float(sigma_v))
		sigma_err.append(float(sigma_v_err))
		lumi_sp.append(float(l_sp))
		lumi_sp_err.append(float(l_sp_err))
		capsigma_x.append(float(Sigma_x))
		capsigma_x_err.append(float(Sigma_x_err))
		capsigma_y.append(float(Sigma_y))
		capsigma_y_err.append(float(Sigma_y_err))
		mumax_x.append(float(Peak_mu_x))
		mumax_x_err.append(float(Peak_mu_x_err))
		mumax_y.append(float(Peak_mu_y))
		mumax_y_err.append(float(Peak_mu_y_err))
		chi2_x.append(float(chi2x))
		chi2_y.append(float(chi2y))


print "len(bcids):", len(bcids)

############# LumiSp Plot

g_lumi_sp = ROOT.TGraphErrors(36)
g_lumi_sp.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
for i in range(len(lumi_sp)):
	g_lumi_sp.SetPoint( i, i, lumi_sp[i]/lumi_sp[0] )
	g_lumi_sp.SetPointError( i, 0.0, math.sqrt( pow(lumi_sp[i]/lumi_sp[0],2) * ( pow(lumi_sp_err[i]/lumi_sp[i],2) + pow(lumi_sp_err[0]/lumi_sp[0],2) ) ) )
g_lumi_sp.SetTitle("Comparison of L_{sp} per BCID per Scan")
g_lumi_sp.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
for i in range(len(bcids)):
	g_lumi_sp.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
g_lumi_sp.GetHistogram().LabelsOption("v","X")
g_lumi_sp.GetXaxis().SetTitle("BCIDs")
g_lumi_sp.GetXaxis().SetTitleOffset(1.3)
g_lumi_sp.GetYaxis().SetTitle("L_{sp}(vtx)/L_{sp}(1stBCID)")
g_lumi_sp.GetYaxis().SetTitleOffset(1.3)
g_lumi_sp.SetMarkerStyle(8)

ay = g_lumi_sp.GetHistogram().GetYaxis()
bottom = ay.GetBinLowEdge(1)
top = ay.GetBinUpEdge(ay.GetNbins())
bottom = 0.9
top = 1.1

l1 = ROOT.TLine(-0.5,bottom,-0.5,top)
l1.SetLineStyle(2)
l1.SetLineColor(ROOT.kRed)

l2 = ROOT.TLine(3.5,bottom,3.5,top)
l2.SetLineStyle(2)
l2.SetLineColor(ROOT.kRed)

l3 = ROOT.TLine(7.5,bottom,7.5,top)
l3.SetLineStyle(2)
l3.SetLineColor(ROOT.kRed)

l4 = ROOT.TLine(11.25,bottom,11.25,top)
l4.SetLineStyle(2)
l4.SetLineColor(ROOT.kRed)

l5 = ROOT.TLine(11.75,bottom,11.75,top)
l5.SetLineStyle(2)
l5.SetLineColor(ROOT.kBlue)

l6 = ROOT.TLine(14.5,bottom,14.5,top)
l6.SetLineStyle(2)
l6.SetLineColor(ROOT.kBlue)

l7 = ROOT.TLine(17.5,bottom,17.5,top)
l7.SetLineStyle(2)
l7.SetLineColor(ROOT.kBlue)

l8 = ROOT.TLine(20.5,bottom,20.5,top)
l8.SetLineStyle(2)
l8.SetLineColor(ROOT.kBlue)

l9 = ROOT.TLine(23.25,bottom,23.25,top)
l9.SetLineStyle(2)
l9.SetLineColor(ROOT.kBlue)

l10 = ROOT.TLine(23.75,bottom,23.75,top)
l10.SetLineStyle(2)
l10.SetLineColor(ROOT.kGreen)

l11 = ROOT.TLine(26.5,bottom,26.5,top)
l11.SetLineStyle(2)
l11.SetLineColor(ROOT.kGreen)

l12 = ROOT.TLine(29.5,bottom,29.5,top)
l12.SetLineStyle(2)
l12.SetLineColor(ROOT.kGreen)

l13 = ROOT.TLine(32.5,bottom,32.5,top)
l13.SetLineStyle(2)
l13.SetLineColor(ROOT.kGreen)

l14 = ROOT.TLine(35.25,bottom,35.25,top)
l14.SetLineStyle(2)
l14.SetLineColor(ROOT.kGreen)

c_lumi_sp = ROOT.TCanvas( "c_lumi_sp", "c_lumi_sp", 800, 600 )
#Pad = c_lumi_sp.cd(1)
#Pad.SetTicks(1,1)

g_lumi_sp.SetMaximum(1.1)
g_lumi_sp.SetMinimum(0.9)

g_lumi_sp.Draw("AP")
l1.Draw('same')
l2.Draw('same')
l3.Draw('same')
l4.Draw('same')
l5.Draw('same')
l6.Draw('same')
l7.Draw('same')
l8.Draw('same')
l9.Draw('same')
l10.Draw('same')
l11.Draw('same')
l12.Draw('same')
l13.Draw('same')
l14.Draw('same')
c_lumi_sp.Print("InfoFromVdMScans/VertexCountingLumiSpPerBCID_allScans.pdf")

############# SigmaVis Plot normalise to 1st BCID in April

g_sigma_vis = ROOT.TGraphErrors(36)
g_sigma_vis.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
for i in range(len(sigma)):
	print sigma[i]
	g_sigma_vis.SetPoint( i, i, sigma[i]/sigma[0] )
	g_sigma_vis.SetPointError( i, 0.0, math.sqrt( pow(sigma[i]/sigma[0],2) * ( pow(sigma_err[i]/sigma[i],2) + pow(sigma_err[0]/sigma[0],2) ) ) )
g_sigma_vis.SetTitle("Comparison of #sigma_{vis} per BCID per Scan")
g_sigma_vis.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
for i in range(len(bcids)):
	g_sigma_vis.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
g_sigma_vis.GetHistogram().LabelsOption("v","X")
g_sigma_vis.GetXaxis().SetTitle("BCIDs")
g_sigma_vis.GetXaxis().SetTitleOffset(1.3)
g_sigma_vis.GetYaxis().SetTitle("#sigma_{vis} normalise to 1st BCID in April")
g_sigma_vis.GetYaxis().SetTitleOffset(1.3)
g_sigma_vis.SetMarkerStyle(8)

ay = g_sigma_vis.GetHistogram().GetYaxis()
bottom = ay.GetBinLowEdge(1)
top = ay.GetBinUpEdge(ay.GetNbins())
bottom = 0.9
top = 1.1

l1 = ROOT.TLine(-0.5,bottom,-0.5,top)
l1.SetLineStyle(2)
l1.SetLineColor(ROOT.kRed)

l2 = ROOT.TLine(3.5,bottom,3.5,top)
l2.SetLineStyle(2)
l2.SetLineColor(ROOT.kRed)

l3 = ROOT.TLine(7.5,bottom,7.5,top)
l3.SetLineStyle(2)
l3.SetLineColor(ROOT.kRed)

l4 = ROOT.TLine(11.25,bottom,11.25,top)
l4.SetLineStyle(2)
l4.SetLineColor(ROOT.kRed)

l5 = ROOT.TLine(11.75,bottom,11.75,top)
l5.SetLineStyle(2)
l5.SetLineColor(ROOT.kBlue)

l6 = ROOT.TLine(14.5,bottom,14.5,top)
l6.SetLineStyle(2)
l6.SetLineColor(ROOT.kBlue)

l7 = ROOT.TLine(17.5,bottom,17.5,top)
l7.SetLineStyle(2)
l7.SetLineColor(ROOT.kBlue)

l8 = ROOT.TLine(20.5,bottom,20.5,top)
l8.SetLineStyle(2)
l8.SetLineColor(ROOT.kBlue)

l9 = ROOT.TLine(23.25,bottom,23.25,top)
l9.SetLineStyle(2)
l9.SetLineColor(ROOT.kBlue)

l10 = ROOT.TLine(23.75,bottom,23.75,top)
l10.SetLineStyle(2)
l10.SetLineColor(ROOT.kGreen)

l11 = ROOT.TLine(26.5,bottom,26.5,top)
l11.SetLineStyle(2)
l11.SetLineColor(ROOT.kGreen)

l12 = ROOT.TLine(29.5,bottom,29.5,top)
l12.SetLineStyle(2)
l12.SetLineColor(ROOT.kGreen)

l13 = ROOT.TLine(32.5,bottom,32.5,top)
l13.SetLineStyle(2)
l13.SetLineColor(ROOT.kGreen)

l14 = ROOT.TLine(35.25,bottom,35.25,top)
l14.SetLineStyle(2)
l14.SetLineColor(ROOT.kGreen)

c_sigma = ROOT.TCanvas( "c_sigma", "c_sigma", 800, 600 )
#Pad = c_sigma.cd(1)
#Pad.SetTicks(1,1)

g_sigma_vis.SetMaximum(1.1)
g_sigma_vis.SetMinimum(0.9)

g_sigma_vis.Draw("AP")
l1.Draw('same')
l2.Draw('same')
l3.Draw('same')
l4.Draw('same')
l5.Draw('same')
l6.Draw('same')
l7.Draw('same')
l8.Draw('same')
l9.Draw('same')
l10.Draw('same')
l11.Draw('same')
l12.Draw('same')
l13.Draw('same')
l14.Draw('same')
c_sigma.Print("InfoFromVdMScans/VertexCountingSigmaVisPerBCID_allScans_NormalisedToApril.pdf")

############# SigmaVis Plot normalise to 1st BCID in every scan set

g_sigma_vis_A = ROOT.TGraphErrors(12)
g_sigma_vis_A.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_sigma_vis_A.SetMarkerStyle(8)
for i in range(12):
	g_sigma_vis_A.SetPoint( i, i, sigma[i]/sigma[0] )
	g_sigma_vis_A.SetPointError( i, 0.0, math.sqrt( pow(sigma[i]/sigma[0],2) * ( pow(sigma_err[i]/sigma[i],2) + pow(sigma_err[0]/sigma[0],2) ) ) )

g_sigma_vis_J = ROOT.TGraphErrors(12)
g_sigma_vis_J.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_sigma_vis_J.SetMarkerStyle(21)
for i in range(12):
	g_sigma_vis_J.SetPoint( i, i+12, sigma[i+12]/sigma[12] )
	g_sigma_vis_J.SetPointError( i, 0.0, math.sqrt( pow(sigma[i+12]/sigma[12],2) * ( pow(sigma_err[i+12]/sigma[i+12],2) + pow(sigma_err[12]/sigma[12],2) ) ) )

g_sigma_vis_N = ROOT.TGraphErrors(12)
g_sigma_vis_N.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_sigma_vis_N.SetMarkerStyle(22)
for i in range(12):
	g_sigma_vis_N.SetPoint( i, i+24, sigma[i+24]/sigma[24] )
	g_sigma_vis_N.SetPointError( i, 0.0, math.sqrt( pow(sigma[i+24]/sigma[24],2) * ( pow(sigma_err[i+24]/sigma[i+24],2) + pow(sigma_err[24]/sigma[24],2) ) ) )

g_sigma_vis_A.SetTitle("")
for i in range(len(bcids)):
	g_sigma_vis_A.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_sigma_vis_J.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_sigma_vis_N.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
g_sigma_vis_A.GetHistogram().LabelsOption("v","X")
g_sigma_vis_J.GetHistogram().LabelsOption("v","X")
g_sigma_vis_N.GetHistogram().LabelsOption("v","X")
g_sigma_vis_A.GetXaxis().SetTitle("BCIDs")
g_sigma_vis_A.GetXaxis().SetTitleOffset(1.3)
g_sigma_vis_A.GetYaxis().SetTitle("#sigma_{vis} normalise to 1st BCID in Scan Set")
g_sigma_vis_A.GetYaxis().SetTitleOffset(1.3)

ay = g_sigma_vis_A.GetHistogram().GetYaxis()
bottom = ay.GetBinLowEdge(1)
top = ay.GetBinUpEdge(ay.GetNbins())
bottom = 0.96
top = 1.04

l1 = ROOT.TLine(-0.5,bottom,-0.5,top)
l1.SetLineStyle(2)
l1.SetLineColor(ROOT.kRed)

l2 = ROOT.TLine(3.5,bottom,3.5,top)
l2.SetLineStyle(2)
l2.SetLineColor(ROOT.kRed)

l3 = ROOT.TLine(7.5,bottom,7.5,top)
l3.SetLineStyle(2)
l3.SetLineColor(ROOT.kRed)

l4 = ROOT.TLine(11.25,bottom,11.25,top)
l4.SetLineStyle(2)
l4.SetLineColor(ROOT.kRed)

l5 = ROOT.TLine(11.75,bottom,11.75,top)
l5.SetLineStyle(2)
l5.SetLineColor(ROOT.kBlue)

l6 = ROOT.TLine(14.5,bottom,14.5,top)
l6.SetLineStyle(2)
l6.SetLineColor(ROOT.kBlue)

l7 = ROOT.TLine(17.5,bottom,17.5,top)
l7.SetLineStyle(2)
l7.SetLineColor(ROOT.kBlue)

l8 = ROOT.TLine(20.5,bottom,20.5,top)
l8.SetLineStyle(2)
l8.SetLineColor(ROOT.kBlue)

l9 = ROOT.TLine(23.25,bottom,23.25,top)
l9.SetLineStyle(2)
l9.SetLineColor(ROOT.kBlue)

l10 = ROOT.TLine(23.75,bottom,23.75,top)
l10.SetLineStyle(2)
l10.SetLineColor(ROOT.kGreen)

l11 = ROOT.TLine(26.5,bottom,26.5,top)
l11.SetLineStyle(2)
l11.SetLineColor(ROOT.kGreen)

l12 = ROOT.TLine(29.5,bottom,29.5,top)
l12.SetLineStyle(2)
l12.SetLineColor(ROOT.kGreen)

l13 = ROOT.TLine(32.5,bottom,32.5,top)
l13.SetLineStyle(2)
l13.SetLineColor(ROOT.kGreen)

l14 = ROOT.TLine(35.25,bottom,35.25,top)
l14.SetLineStyle(2)
l14.SetLineColor(ROOT.kGreen)

hline1 = ROOT.TLine(-1.5,1,36.5,1)
hline1.SetLineStyle(2)
hline1.SetLineColor(ROOT.kBlack)

c_sigma = ROOT.TCanvas( "c_sigma", "c_sigma", 800, 600 )
#Pad = c_sigma.cd(1)
#Pad.SetTicks(1,1)
g_sigma_vis_A.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_sigma_vis_A.SetMaximum(1.04)
g_sigma_vis_A.SetMinimum(0.96)

g_sigma_vis_A.Draw("AP")
g_sigma_vis_J.Draw("P")
g_sigma_vis_N.Draw("P")
l1.Draw('same')
l2.Draw('same')
l3.Draw('same')
l4.Draw('same')
l5.Draw('same')
l6.Draw('same')
l7.Draw('same')
l8.Draw('same')
l9.Draw('same')
l10.Draw('same')
l11.Draw('same')
l12.Draw('same')
l13.Draw('same')
l14.Draw('same')
hline1.Draw('same')
c_sigma.Print("InfoFromVdMScans/VertexCountingSigmaVisPerBCID_allScans_NormaliseToScanSets.pdf")

############# LumiSp Plot normalise to 1st BCID in every scan set

g_lumisp_A = ROOT.TGraphErrors(12)
g_lumisp_A.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_lumisp_A.SetMarkerStyle(8)
for i in range(12):
	g_lumisp_A.SetPoint( i, i, lumi_sp[i]/lumi_sp[0] )
	g_lumisp_A.SetPointError( i, 0.0, math.sqrt( pow(lumi_sp[i]/lumi_sp[0],2) * ( pow(lumi_sp_err[i]/lumi_sp[i],2) + pow(lumi_sp_err[0]/lumi_sp[0],2) ) ) )

g_lumisp_J = ROOT.TGraphErrors(12)
g_lumisp_J.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_lumisp_J.SetMarkerStyle(21)
for i in range(12):
	g_lumisp_J.SetPoint( i, i+12, lumi_sp[i+12]/lumi_sp[12] )
	g_lumisp_J.SetPointError( i, 0.0, math.sqrt( pow(lumi_sp[i+12]/lumi_sp[12],2) * ( pow(lumi_sp_err[i+12]/lumi_sp[i+12],2) + pow(lumi_sp_err[12]/lumi_sp[12],2) ) ) )

g_lumisp_N = ROOT.TGraphErrors(12)
g_lumisp_N.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_lumisp_N.SetMarkerStyle(22)
for i in range(12):
	g_lumisp_N.SetPoint( i, i+24, lumi_sp[i+24]/lumi_sp[24] )
	g_lumisp_N.SetPointError( i, 0.0, math.sqrt( pow(lumi_sp[i+24]/lumi_sp[24],2) * ( pow(lumi_sp_err[i+24]/lumi_sp[i+24],2) + pow(lumi_sp_err[24]/lumi_sp[24],2) ) ) )

#g_lumisp_A.SetTitle("Comparison of #sigma_{vis} per BCID per Scan")
g_lumisp_A.SetTitle("")
for i in range(len(bcids)):
	g_lumisp_A.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_lumisp_J.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_lumisp_N.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
g_lumisp_A.GetHistogram().LabelsOption("v","X")
g_lumisp_J.GetHistogram().LabelsOption("v","X")
g_lumisp_N.GetHistogram().LabelsOption("v","X")
g_lumisp_A.GetXaxis().SetTitle("BCIDs")
g_lumisp_A.GetXaxis().SetTitleOffset(1.3)
g_lumisp_A.GetYaxis().SetTitle("L_{sp} normalise to 1st BCID in Scan Set")
g_lumisp_A.GetYaxis().SetTitleOffset(1.3)

ay = g_lumisp_A.GetHistogram().GetYaxis()
bottom = ay.GetBinLowEdge(1)
top = ay.GetBinUpEdge(ay.GetNbins())
bottom = 0.85
top = 1.15

l1 = ROOT.TLine(-0.5,bottom,-0.5,top)
l1.SetLineStyle(2)
l1.SetLineColor(ROOT.kRed)

l2 = ROOT.TLine(3.5,bottom,3.5,top)
l2.SetLineStyle(2)
l2.SetLineColor(ROOT.kRed)

l3 = ROOT.TLine(7.5,bottom,7.5,top)
l3.SetLineStyle(2)
l3.SetLineColor(ROOT.kRed)

l4 = ROOT.TLine(11.25,bottom,11.25,top)
l4.SetLineStyle(2)
l4.SetLineColor(ROOT.kRed)

l5 = ROOT.TLine(11.75,bottom,11.75,top)
l5.SetLineStyle(2)
l5.SetLineColor(ROOT.kBlue)

l6 = ROOT.TLine(14.5,bottom,14.5,top)
l6.SetLineStyle(2)
l6.SetLineColor(ROOT.kBlue)

l7 = ROOT.TLine(17.5,bottom,17.5,top)
l7.SetLineStyle(2)
l7.SetLineColor(ROOT.kBlue)

l8 = ROOT.TLine(20.5,bottom,20.5,top)
l8.SetLineStyle(2)
l8.SetLineColor(ROOT.kBlue)

l9 = ROOT.TLine(23.25,bottom,23.25,top)
l9.SetLineStyle(2)
l9.SetLineColor(ROOT.kBlue)

l10 = ROOT.TLine(23.75,bottom,23.75,top)
l10.SetLineStyle(2)
l10.SetLineColor(ROOT.kGreen)

l11 = ROOT.TLine(26.5,bottom,26.5,top)
l11.SetLineStyle(2)
l11.SetLineColor(ROOT.kGreen)

l12 = ROOT.TLine(29.5,bottom,29.5,top)
l12.SetLineStyle(2)
l12.SetLineColor(ROOT.kGreen)

l13 = ROOT.TLine(32.5,bottom,32.5,top)
l13.SetLineStyle(2)
l13.SetLineColor(ROOT.kGreen)

l14 = ROOT.TLine(35.25,bottom,35.25,top)
l14.SetLineStyle(2)
l14.SetLineColor(ROOT.kGreen)

c_lumisp = ROOT.TCanvas( "c_lumisp", "c_lumisp", 800, 600 )
g_lumisp_A.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_lumisp_A.SetMaximum(1.15)
g_lumisp_A.SetMinimum(0.85)

g_lumisp_A.Draw("AP")
g_lumisp_J.Draw("P")
g_lumisp_N.Draw("P")
l1.Draw('same')
l2.Draw('same')
l3.Draw('same')
l4.Draw('same')
l5.Draw('same')
l6.Draw('same')
l7.Draw('same')
l8.Draw('same')
l9.Draw('same')
l10.Draw('same')
l11.Draw('same')
l12.Draw('same')
l13.Draw('same')
l14.Draw('same')
hline1.Draw('same')
c_lumisp.Print("InfoFromVdMScans/VertexCountingLumiSpPerBCID_allScans_NormaliseToScanSets.pdf")

############# CapSigmaX Plot normalise to 1st BCID in every scan set

g_capsigma_x_A = ROOT.TGraphErrors(12)
g_capsigma_x_A.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_capsigma_x_A.SetMarkerStyle(8)
for i in range(12):
	g_capsigma_x_A.SetPoint( i, i, capsigma_x[i]/capsigma_x[0] )
	g_capsigma_x_A.SetPointError( i, 0.0, math.sqrt( pow(capsigma_x[i]/capsigma_x[0],2) * ( pow(capsigma_x_err[i]/capsigma_x[i],2) + pow(capsigma_x_err[0]/capsigma_x[0],2) ) ) )

g_capsigma_x_J = ROOT.TGraphErrors(12)
g_capsigma_x_J.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_capsigma_x_J.SetMarkerStyle(21)
for i in range(12):
	g_capsigma_x_J.SetPoint( i, i+12, capsigma_x[i+12]/capsigma_x[12] )
	g_capsigma_x_J.SetPointError( i, 0.0, math.sqrt( pow(capsigma_x[i+12]/capsigma_x[12],2) * ( pow(capsigma_x_err[i+12]/capsigma_x[i+12],2) + pow(capsigma_x_err[12]/capsigma_x[12],2) ) ) )

g_capsigma_x_N = ROOT.TGraphErrors(12)
g_capsigma_x_N.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_capsigma_x_N.SetMarkerStyle(22)
for i in range(12):
	g_capsigma_x_N.SetPoint( i, i+24, capsigma_x[i+24]/capsigma_x[24] )
	g_capsigma_x_N.SetPointError( i, 0.0, math.sqrt( pow(capsigma_x[i+24]/capsigma_x[24],2) * ( pow(capsigma_x_err[i+24]/capsigma_x[i+24],2) + pow(capsigma_x_err[24]/capsigma_x[24],2) ) ) )

g_capsigma_x_A.SetTitle("")
for i in range(len(bcids)):
	g_capsigma_x_A.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_capsigma_x_J.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_capsigma_x_N.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
g_capsigma_x_A.GetHistogram().LabelsOption("v","X")
g_capsigma_x_J.GetHistogram().LabelsOption("v","X")
g_capsigma_x_N.GetHistogram().LabelsOption("v","X")
g_capsigma_x_A.GetXaxis().SetTitle("BCIDs")
g_capsigma_x_A.GetXaxis().SetTitleOffset(1.3)
g_capsigma_x_A.GetYaxis().SetTitle("#Sigma_{X} normalise to 1st BCID in Scan Set")
g_capsigma_x_A.GetYaxis().SetTitleOffset(1.3)

ay = g_capsigma_x_A.GetHistogram().GetYaxis()
bottom = ay.GetBinLowEdge(1)
top = ay.GetBinUpEdge(ay.GetNbins())
bottom = 0.85
top = 1.15

l1 = ROOT.TLine(-0.5,bottom,-0.5,top)
l1.SetLineStyle(2)
l1.SetLineColor(ROOT.kRed)

l2 = ROOT.TLine(3.5,bottom,3.5,top)
l2.SetLineStyle(2)
l2.SetLineColor(ROOT.kRed)

l3 = ROOT.TLine(7.5,bottom,7.5,top)
l3.SetLineStyle(2)
l3.SetLineColor(ROOT.kRed)

l4 = ROOT.TLine(11.25,bottom,11.25,top)
l4.SetLineStyle(2)
l4.SetLineColor(ROOT.kRed)

l5 = ROOT.TLine(11.75,bottom,11.75,top)
l5.SetLineStyle(2)
l5.SetLineColor(ROOT.kBlue)

l6 = ROOT.TLine(14.5,bottom,14.5,top)
l6.SetLineStyle(2)
l6.SetLineColor(ROOT.kBlue)

l7 = ROOT.TLine(17.5,bottom,17.5,top)
l7.SetLineStyle(2)
l7.SetLineColor(ROOT.kBlue)

l8 = ROOT.TLine(20.5,bottom,20.5,top)
l8.SetLineStyle(2)
l8.SetLineColor(ROOT.kBlue)

l9 = ROOT.TLine(23.25,bottom,23.25,top)
l9.SetLineStyle(2)
l9.SetLineColor(ROOT.kBlue)

l10 = ROOT.TLine(23.75,bottom,23.75,top)
l10.SetLineStyle(2)
l10.SetLineColor(ROOT.kGreen)

l11 = ROOT.TLine(26.5,bottom,26.5,top)
l11.SetLineStyle(2)
l11.SetLineColor(ROOT.kGreen)

l12 = ROOT.TLine(29.5,bottom,29.5,top)
l12.SetLineStyle(2)
l12.SetLineColor(ROOT.kGreen)

l13 = ROOT.TLine(32.5,bottom,32.5,top)
l13.SetLineStyle(2)
l13.SetLineColor(ROOT.kGreen)

l14 = ROOT.TLine(35.25,bottom,35.25,top)
l14.SetLineStyle(2)
l14.SetLineColor(ROOT.kGreen)

c_capsigmax = ROOT.TCanvas( "c_capsigmax", "c_capsigmax", 800, 600 )
g_capsigma_x_A.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_capsigma_x_A.SetMaximum(1.15)
g_capsigma_x_A.SetMinimum(0.85)

g_capsigma_x_A.Draw("AP")
g_capsigma_x_J.Draw("P")
g_capsigma_x_N.Draw("P")
l1.Draw('same')
l2.Draw('same')
l3.Draw('same')
l4.Draw('same')
l5.Draw('same')
l6.Draw('same')
l7.Draw('same')
l8.Draw('same')
l9.Draw('same')
l10.Draw('same')
l11.Draw('same')
l12.Draw('same')
l13.Draw('same')
l14.Draw('same')
hline1.Draw('same')
c_capsigmax.Print("InfoFromVdMScans/VertexCountingCapSigmaXPerBCID_allScans_NormaliseToScanSets.pdf")

############# CapSigmaY Plot normalise to 1st BCID in every scan set

g_capsigma_y_A = ROOT.TGraphErrors(12)
g_capsigma_y_A.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_capsigma_y_A.SetMarkerStyle(8)
for i in range(12):
	g_capsigma_y_A.SetPoint( i, i, capsigma_y[i]/capsigma_y[0] )
	g_capsigma_y_A.SetPointError( i, 0.0, math.sqrt( pow(capsigma_y[i]/capsigma_y[0],2) * ( pow(capsigma_y_err[i]/capsigma_y[i],2) + pow(capsigma_y_err[0]/capsigma_y[0],2) ) ) )

g_capsigma_y_J = ROOT.TGraphErrors(12)
g_capsigma_y_J.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_capsigma_y_J.SetMarkerStyle(21)
for i in range(12):
	g_capsigma_y_J.SetPoint( i, i+12, capsigma_y[i+12]/capsigma_y[12] )
	g_capsigma_y_J.SetPointError( i, 0.0, math.sqrt( pow(capsigma_y[i+12]/capsigma_y[12],2) * ( pow(capsigma_y_err[i+12]/capsigma_y[i+12],2) + pow(capsigma_y_err[12]/capsigma_y[12],2) ) ) )

g_capsigma_y_N = ROOT.TGraphErrors(12)
g_capsigma_y_N.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_capsigma_y_N.SetMarkerStyle(22)
for i in range(12):
	g_capsigma_y_N.SetPoint( i, i+24, capsigma_y[i+24]/capsigma_y[24] )
	g_capsigma_y_N.SetPointError( i, 0.0, math.sqrt( pow(capsigma_y[i+24]/capsigma_y[24],2) * ( pow(capsigma_y_err[i+24]/capsigma_y[i+24],2) + pow(capsigma_y_err[24]/capsigma_y[24],2) ) ) )

g_capsigma_y_A.SetTitle("")
for i in range(len(bcids)):
	g_capsigma_y_A.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_capsigma_y_J.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_capsigma_y_N.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
g_capsigma_y_A.GetHistogram().LabelsOption("v","X")
g_capsigma_y_J.GetHistogram().LabelsOption("v","X")
g_capsigma_y_N.GetHistogram().LabelsOption("v","X")
g_capsigma_y_A.GetXaxis().SetTitle("BCIDs")
g_capsigma_y_A.GetXaxis().SetTitleOffset(1.3)
g_capsigma_y_A.GetYaxis().SetTitle("#Sigma_{Y} normalise to 1st BCID in Scan Set")
g_capsigma_y_A.GetYaxis().SetTitleOffset(1.3)

ay = g_capsigma_y_A.GetHistogram().GetYaxis()
bottom = ay.GetBinLowEdge(1)
top = ay.GetBinUpEdge(ay.GetNbins())
bottom = 0.85
top = 1.15

l1 = ROOT.TLine(-0.5,bottom,-0.5,top)
l1.SetLineStyle(2)
l1.SetLineColor(ROOT.kRed)

l2 = ROOT.TLine(3.5,bottom,3.5,top)
l2.SetLineStyle(2)
l2.SetLineColor(ROOT.kRed)

l3 = ROOT.TLine(7.5,bottom,7.5,top)
l3.SetLineStyle(2)
l3.SetLineColor(ROOT.kRed)

l4 = ROOT.TLine(11.25,bottom,11.25,top)
l4.SetLineStyle(2)
l4.SetLineColor(ROOT.kRed)

l5 = ROOT.TLine(11.75,bottom,11.75,top)
l5.SetLineStyle(2)
l5.SetLineColor(ROOT.kBlue)

l6 = ROOT.TLine(14.5,bottom,14.5,top)
l6.SetLineStyle(2)
l6.SetLineColor(ROOT.kBlue)

l7 = ROOT.TLine(17.5,bottom,17.5,top)
l7.SetLineStyle(2)
l7.SetLineColor(ROOT.kBlue)

l8 = ROOT.TLine(20.5,bottom,20.5,top)
l8.SetLineStyle(2)
l8.SetLineColor(ROOT.kBlue)

l9 = ROOT.TLine(23.25,bottom,23.25,top)
l9.SetLineStyle(2)
l9.SetLineColor(ROOT.kBlue)

l10 = ROOT.TLine(23.75,bottom,23.75,top)
l10.SetLineStyle(2)
l10.SetLineColor(ROOT.kGreen)

l11 = ROOT.TLine(26.5,bottom,26.5,top)
l11.SetLineStyle(2)
l11.SetLineColor(ROOT.kGreen)

l12 = ROOT.TLine(29.5,bottom,29.5,top)
l12.SetLineStyle(2)
l12.SetLineColor(ROOT.kGreen)

l13 = ROOT.TLine(32.5,bottom,32.5,top)
l13.SetLineStyle(2)
l13.SetLineColor(ROOT.kGreen)

l14 = ROOT.TLine(35.25,bottom,35.25,top)
l14.SetLineStyle(2)
l14.SetLineColor(ROOT.kGreen)

c_capsigmay = ROOT.TCanvas( "c_capsigmay", "c_capsigmay", 800, 600 )
g_capsigma_y_A.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_capsigma_y_A.SetMaximum(1.15)
g_capsigma_y_A.SetMinimum(0.85)

g_capsigma_y_A.Draw("AP")
g_capsigma_y_J.Draw("P")
g_capsigma_y_N.Draw("P")
l1.Draw('same')
l2.Draw('same')
l3.Draw('same')
l4.Draw('same')
l5.Draw('same')
l6.Draw('same')
l7.Draw('same')
l8.Draw('same')
l9.Draw('same')
l10.Draw('same')
l11.Draw('same')
l12.Draw('same')
l13.Draw('same')
l14.Draw('same')
hline1.Draw('same')
c_capsigmay.Print("InfoFromVdMScans/VertexCountingCapSigmaYPerBCID_allScans_NormaliseToScanSets.pdf")

############# MuMaxX Plot normalise to 1st BCID in every scan set

g_mumax_x_A = ROOT.TGraphErrors(12)
g_mumax_x_A.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_mumax_x_A.SetMarkerStyle(8)
for i in range(12):
	g_mumax_x_A.SetPoint( i, i, mumax_x[i]/mumax_x[0] )
	g_mumax_x_A.SetPointError( i, 0.0, math.sqrt( pow(mumax_x[i]/mumax_x[0],2) * ( pow(mumax_x_err[i]/mumax_x[i],2) + pow(mumax_x_err[0]/mumax_x[0],2) ) ) )

g_mumax_x_J = ROOT.TGraphErrors(12)
g_mumax_x_J.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_mumax_x_J.SetMarkerStyle(21)
for i in range(12):
	g_mumax_x_J.SetPoint( i, i+12, mumax_x[i+12]/mumax_x[12] )
	g_mumax_x_J.SetPointError( i, 0.0, math.sqrt( pow(mumax_x[i+12]/mumax_x[12],2) * ( pow(mumax_x_err[i+12]/mumax_x[i+12],2) + pow(mumax_x_err[12]/mumax_x[12],2) ) ) )

g_mumax_x_N = ROOT.TGraphErrors(12)
g_mumax_x_N.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_mumax_x_N.SetMarkerStyle(22)
for i in range(12):
	g_mumax_x_N.SetPoint( i, i+24, mumax_x[i+24]/mumax_x[24] )
	g_mumax_x_N.SetPointError( i, 0.0, math.sqrt( pow(mumax_x[i+24]/mumax_x[24],2) * ( pow(mumax_x_err[i+24]/mumax_x[i+24],2) + pow(mumax_x_err[24]/mumax_x[24],2) ) ) )

g_mumax_x_A.SetTitle("")
for i in range(len(bcids)):
	g_mumax_x_A.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_mumax_x_J.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_mumax_x_N.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
g_mumax_x_A.GetHistogram().LabelsOption("v","X")
g_mumax_x_J.GetHistogram().LabelsOption("v","X")
g_mumax_x_N.GetHistogram().LabelsOption("v","X")
g_mumax_x_A.GetXaxis().SetTitle("BCIDs")
g_mumax_x_A.GetXaxis().SetTitleOffset(1.3)
g_mumax_x_A.GetYaxis().SetTitle("#mu_{max,x} normalise to 1st BCID in Scan Set")
g_mumax_x_A.GetYaxis().SetTitleOffset(1.3)

ay = g_mumax_x_A.GetHistogram().GetYaxis()
bottom = ay.GetBinLowEdge(1)
top = ay.GetBinUpEdge(ay.GetNbins())
bottom = 0.85
top = 1.15

l1 = ROOT.TLine(-0.5,bottom,-0.5,top)
l1.SetLineStyle(2)
l1.SetLineColor(ROOT.kRed)

l2 = ROOT.TLine(3.5,bottom,3.5,top)
l2.SetLineStyle(2)
l2.SetLineColor(ROOT.kRed)

l3 = ROOT.TLine(7.5,bottom,7.5,top)
l3.SetLineStyle(2)
l3.SetLineColor(ROOT.kRed)

l4 = ROOT.TLine(11.25,bottom,11.25,top)
l4.SetLineStyle(2)
l4.SetLineColor(ROOT.kRed)

l5 = ROOT.TLine(11.75,bottom,11.75,top)
l5.SetLineStyle(2)
l5.SetLineColor(ROOT.kBlue)

l6 = ROOT.TLine(14.5,bottom,14.5,top)
l6.SetLineStyle(2)
l6.SetLineColor(ROOT.kBlue)

l7 = ROOT.TLine(17.5,bottom,17.5,top)
l7.SetLineStyle(2)
l7.SetLineColor(ROOT.kBlue)

l8 = ROOT.TLine(20.5,bottom,20.5,top)
l8.SetLineStyle(2)
l8.SetLineColor(ROOT.kBlue)

l9 = ROOT.TLine(23.25,bottom,23.25,top)
l9.SetLineStyle(2)
l9.SetLineColor(ROOT.kBlue)

l10 = ROOT.TLine(23.75,bottom,23.75,top)
l10.SetLineStyle(2)
l10.SetLineColor(ROOT.kGreen)

l11 = ROOT.TLine(26.5,bottom,26.5,top)
l11.SetLineStyle(2)
l11.SetLineColor(ROOT.kGreen)

l12 = ROOT.TLine(29.5,bottom,29.5,top)
l12.SetLineStyle(2)
l12.SetLineColor(ROOT.kGreen)

l13 = ROOT.TLine(32.5,bottom,32.5,top)
l13.SetLineStyle(2)
l13.SetLineColor(ROOT.kGreen)

l14 = ROOT.TLine(35.25,bottom,35.25,top)
l14.SetLineStyle(2)
l14.SetLineColor(ROOT.kGreen)

c_mumax_x = ROOT.TCanvas( "c_mumax_x", "c_mumax_x", 800, 600 )
g_mumax_x_A.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_mumax_x_A.SetMaximum(1.15)
g_mumax_x_A.SetMinimum(0.85)

g_mumax_x_A.Draw("AP")
g_mumax_x_J.Draw("P")
g_mumax_x_N.Draw("P")
l1.Draw('same')
l2.Draw('same')
l3.Draw('same')
l4.Draw('same')
l5.Draw('same')
l6.Draw('same')
l7.Draw('same')
l8.Draw('same')
l9.Draw('same')
l10.Draw('same')
l11.Draw('same')
l12.Draw('same')
l13.Draw('same')
l14.Draw('same')
hline1.Draw('same')
c_mumax_x.Print("InfoFromVdMScans/VertexCountingMuMaxXPerBCID_allScans_NormaliseToScanSets.pdf")

############# MuMaxY Plot normalise to 1st BCID in every scan set

g_mumax_y_A = ROOT.TGraphErrors(12)
g_mumax_y_A.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_mumax_y_A.SetMarkerStyle(8)
for i in range(12):
	g_mumax_y_A.SetPoint( i, i, mumax_y[i]/mumax_y[0] )
	g_mumax_y_A.SetPointError( i, 0.0, math.sqrt( pow(mumax_y[i]/mumax_y[0],2) * ( pow(mumax_y_err[i]/mumax_y[i],2) + pow(mumax_y_err[0]/mumax_y[0],2) ) ) )

g_mumax_y_J = ROOT.TGraphErrors(12)
g_mumax_y_J.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_mumax_y_J.SetMarkerStyle(21)
for i in range(12):
	g_mumax_y_J.SetPoint( i, i+12, mumax_y[i+12]/mumax_y[12] )
	g_mumax_y_J.SetPointError( i, 0.0, math.sqrt( pow(mumax_y[i+12]/mumax_y[12],2) * ( pow(mumax_y_err[i+12]/mumax_y[i+12],2) + pow(mumax_y_err[12]/mumax_y[12],2) ) ) )

g_mumax_y_N = ROOT.TGraphErrors(12)
g_mumax_y_N.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_mumax_y_N.SetMarkerStyle(22)
for i in range(12):
	g_mumax_y_N.SetPoint( i, i+24, mumax_y[i+24]/mumax_y[24] )
	g_mumax_y_N.SetPointError( i, 0.0, math.sqrt( pow(mumax_y[i+24]/mumax_y[24],2) * ( pow(mumax_y_err[i+24]/mumax_y[i+24],2) + pow(mumax_y_err[24]/mumax_y[24],2) ) ) )

g_mumax_y_A.SetTitle("")
for i in range(len(bcids)):
	g_mumax_y_A.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_mumax_y_J.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_mumax_y_N.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
g_mumax_y_A.GetHistogram().LabelsOption("v","X")
g_mumax_y_J.GetHistogram().LabelsOption("v","X")
g_mumax_y_N.GetHistogram().LabelsOption("v","X")
g_mumax_y_A.GetXaxis().SetTitle("BCIDs")
g_mumax_y_A.GetXaxis().SetTitleOffset(1.3)
g_mumax_y_A.GetYaxis().SetTitle("#mu_{max,y} normalise to 1st BCID in Scan Set")
g_mumax_y_A.GetYaxis().SetTitleOffset(1.3)

ay = g_mumax_y_A.GetHistogram().GetYaxis()
bottom = ay.GetBinLowEdge(1)
top = ay.GetBinUpEdge(ay.GetNbins())
bottom = 0.85
top = 1.15

l1 = ROOT.TLine(-0.5,bottom,-0.5,top)
l1.SetLineStyle(2)
l1.SetLineColor(ROOT.kRed)

l2 = ROOT.TLine(3.5,bottom,3.5,top)
l2.SetLineStyle(2)
l2.SetLineColor(ROOT.kRed)

l3 = ROOT.TLine(7.5,bottom,7.5,top)
l3.SetLineStyle(2)
l3.SetLineColor(ROOT.kRed)

l4 = ROOT.TLine(11.25,bottom,11.25,top)
l4.SetLineStyle(2)
l4.SetLineColor(ROOT.kRed)

l5 = ROOT.TLine(11.75,bottom,11.75,top)
l5.SetLineStyle(2)
l5.SetLineColor(ROOT.kBlue)

l6 = ROOT.TLine(14.5,bottom,14.5,top)
l6.SetLineStyle(2)
l6.SetLineColor(ROOT.kBlue)

l7 = ROOT.TLine(17.5,bottom,17.5,top)
l7.SetLineStyle(2)
l7.SetLineColor(ROOT.kBlue)

l8 = ROOT.TLine(20.5,bottom,20.5,top)
l8.SetLineStyle(2)
l8.SetLineColor(ROOT.kBlue)

l9 = ROOT.TLine(23.25,bottom,23.25,top)
l9.SetLineStyle(2)
l9.SetLineColor(ROOT.kBlue)

l10 = ROOT.TLine(23.75,bottom,23.75,top)
l10.SetLineStyle(2)
l10.SetLineColor(ROOT.kGreen)

l11 = ROOT.TLine(26.5,bottom,26.5,top)
l11.SetLineStyle(2)
l11.SetLineColor(ROOT.kGreen)

l12 = ROOT.TLine(29.5,bottom,29.5,top)
l12.SetLineStyle(2)
l12.SetLineColor(ROOT.kGreen)

l13 = ROOT.TLine(32.5,bottom,32.5,top)
l13.SetLineStyle(2)
l13.SetLineColor(ROOT.kGreen)

l14 = ROOT.TLine(35.25,bottom,35.25,top)
l14.SetLineStyle(2)
l14.SetLineColor(ROOT.kGreen)

c_mumax_y = ROOT.TCanvas( "c_mumax_y", "c_mumax_y", 800, 600 )
g_mumax_y_A.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_mumax_y_A.SetMaximum(1.15)
g_mumax_y_A.SetMinimum(0.85)

g_mumax_y_A.Draw("AP")
g_mumax_y_J.Draw("P")
g_mumax_y_N.Draw("P")
l1.Draw('same')
l2.Draw('same')
l3.Draw('same')
l4.Draw('same')
l5.Draw('same')
l6.Draw('same')
l7.Draw('same')
l8.Draw('same')
l9.Draw('same')
l10.Draw('same')
l11.Draw('same')
l12.Draw('same')
l13.Draw('same')
l14.Draw('same')
hline1.Draw('same')
c_mumax_y.Print("InfoFromVdMScans/VertexCountingMuMaxYPerBCID_allScans_NormaliseToScanSets.pdf")

############# Chi2X Plot normalise to 1st BCID in every scan set

g_chi2_x_A = ROOT.TGraphErrors(12)
g_chi2_x_A.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_chi2_x_A.SetMarkerStyle(8)
for i in range(12):
	g_chi2_x_A.SetPoint( i, i, chi2_x[i] )
g_chi2_x_A.Print()

g_chi2_x_J = ROOT.TGraphErrors(12)
g_chi2_x_J.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_chi2_x_J.SetMarkerStyle(21)
for i in range(12):
	g_chi2_x_J.SetPoint( i, i+12, chi2_x[i+12] )

g_chi2_x_N = ROOT.TGraphErrors(12)
g_chi2_x_N.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_chi2_x_N.SetMarkerStyle(22)
for i in range(12):
	g_chi2_x_N.SetPoint( i, i+24, chi2_x[i+24] )

g_chi2_x_A.SetTitle("")
for i in range(len(bcids)):
	g_chi2_x_A.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_chi2_x_J.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_chi2_x_N.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
g_chi2_x_A.GetHistogram().LabelsOption("v","X")
g_chi2_x_J.GetHistogram().LabelsOption("v","X")
g_chi2_x_N.GetHistogram().LabelsOption("v","X")
g_chi2_x_A.GetXaxis().SetTitle("BCIDs")
g_chi2_x_A.GetXaxis().SetTitleOffset(1.3)
g_chi2_x_A.GetYaxis().SetTitle("#Chi^{2}/NDF X-scan ")
g_chi2_x_A.GetYaxis().SetTitleOffset(1.3)

ay = g_chi2_x_A.GetHistogram().GetYaxis()
bottom = ay.GetBinLowEdge(1)
top = ay.GetBinUpEdge(ay.GetNbins())
bottom = 0.0
top = 7.0

l1 = ROOT.TLine(-0.5,bottom,-0.5,top)
l1.SetLineStyle(2)
l1.SetLineColor(ROOT.kRed)

l2 = ROOT.TLine(3.5,bottom,3.5,top)
l2.SetLineStyle(2)
l2.SetLineColor(ROOT.kRed)

l3 = ROOT.TLine(7.5,bottom,7.5,top)
l3.SetLineStyle(2)
l3.SetLineColor(ROOT.kRed)

l4 = ROOT.TLine(11.25,bottom,11.25,top)
l4.SetLineStyle(2)
l4.SetLineColor(ROOT.kRed)

l5 = ROOT.TLine(11.75,bottom,11.75,top)
l5.SetLineStyle(2)
l5.SetLineColor(ROOT.kBlue)

l6 = ROOT.TLine(14.5,bottom,14.5,top)
l6.SetLineStyle(2)
l6.SetLineColor(ROOT.kBlue)

l7 = ROOT.TLine(17.5,bottom,17.5,top)
l7.SetLineStyle(2)
l7.SetLineColor(ROOT.kBlue)

l8 = ROOT.TLine(20.5,bottom,20.5,top)
l8.SetLineStyle(2)
l8.SetLineColor(ROOT.kBlue)

l9 = ROOT.TLine(23.25,bottom,23.25,top)
l9.SetLineStyle(2)
l9.SetLineColor(ROOT.kBlue)

l10 = ROOT.TLine(23.75,bottom,23.75,top)
l10.SetLineStyle(2)
l10.SetLineColor(ROOT.kGreen)

l11 = ROOT.TLine(26.5,bottom,26.5,top)
l11.SetLineStyle(2)
l11.SetLineColor(ROOT.kGreen)

l12 = ROOT.TLine(29.5,bottom,29.5,top)
l12.SetLineStyle(2)
l12.SetLineColor(ROOT.kGreen)

l13 = ROOT.TLine(32.5,bottom,32.5,top)
l13.SetLineStyle(2)
l13.SetLineColor(ROOT.kGreen)

l14 = ROOT.TLine(35.25,bottom,35.25,top)
l14.SetLineStyle(2)
l14.SetLineColor(ROOT.kGreen)

c_chi2_x = ROOT.TCanvas( "c_chi2_x", "c_chi2_x", 800, 600 )
g_chi2_x_A.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_chi2_x_A.SetMaximum(7.0)
g_chi2_x_A.SetMinimum(0.0)

g_chi2_x_A.Draw("AP")
g_chi2_x_J.Draw("P")
g_chi2_x_N.Draw("P")
l1.Draw('same')
l2.Draw('same')
l3.Draw('same')
l4.Draw('same')
l5.Draw('same')
l6.Draw('same')
l7.Draw('same')
l8.Draw('same')
l9.Draw('same')
l10.Draw('same')
l11.Draw('same')
l12.Draw('same')
l13.Draw('same')
l14.Draw('same')
hline1.Draw('same')
c_chi2_x.Print("InfoFromVdMScans/VertexCountingChi2XPerBCID_allScans_NormaliseToScanSets.pdf")

############# Chi2Y Plot normalise to 1st BCID in every scan set

g_chi2_y_A = ROOT.TGraphErrors(12)
g_chi2_y_A.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_chi2_y_A.SetMarkerStyle(8)
for i in range(12):
	g_chi2_y_A.SetPoint( i, i, chi2_y[i] )

g_chi2_y_J = ROOT.TGraphErrors(12)
g_chi2_y_J.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_chi2_y_J.SetMarkerStyle(21)
for i in range(12):
	g_chi2_y_J.SetPoint( i, i+12, chi2_y[i+12] )

g_chi2_y_N = ROOT.TGraphErrors(12)
g_chi2_y_N.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_chi2_y_N.SetMarkerStyle(22)
for i in range(12):
	g_chi2_y_N.SetPoint( i, i+24, chi2_y[i+24] )

g_chi2_y_A.SetTitle("")
for i in range(len(bcids)):
	g_chi2_y_A.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_chi2_y_J.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_chi2_y_N.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
g_chi2_y_A.GetHistogram().LabelsOption("v","X")
g_chi2_y_J.GetHistogram().LabelsOption("v","X")
g_chi2_y_N.GetHistogram().LabelsOption("v","X")
g_chi2_y_A.GetXaxis().SetTitle("BCIDs")
g_chi2_y_A.GetXaxis().SetTitleOffset(1.3)
g_chi2_y_A.GetYaxis().SetTitle("#Chi^{2}/NDF Y-scan ")
g_chi2_y_A.GetYaxis().SetTitleOffset(1.3)

ay = g_chi2_y_A.GetHistogram().GetYaxis()
bottom = ay.GetBinLowEdge(1)
top = ay.GetBinUpEdge(ay.GetNbins())
bottom = 0.0
top = 7.0

l1 = ROOT.TLine(-0.5,bottom,-0.5,top)
l1.SetLineStyle(2)
l1.SetLineColor(ROOT.kRed)

l2 = ROOT.TLine(3.5,bottom,3.5,top)
l2.SetLineStyle(2)
l2.SetLineColor(ROOT.kRed)

l3 = ROOT.TLine(7.5,bottom,7.5,top)
l3.SetLineStyle(2)
l3.SetLineColor(ROOT.kRed)

l4 = ROOT.TLine(11.25,bottom,11.25,top)
l4.SetLineStyle(2)
l4.SetLineColor(ROOT.kRed)

l5 = ROOT.TLine(11.75,bottom,11.75,top)
l5.SetLineStyle(2)
l5.SetLineColor(ROOT.kBlue)

l6 = ROOT.TLine(14.5,bottom,14.5,top)
l6.SetLineStyle(2)
l6.SetLineColor(ROOT.kBlue)

l7 = ROOT.TLine(17.5,bottom,17.5,top)
l7.SetLineStyle(2)
l7.SetLineColor(ROOT.kBlue)

l8 = ROOT.TLine(20.5,bottom,20.5,top)
l8.SetLineStyle(2)
l8.SetLineColor(ROOT.kBlue)

l9 = ROOT.TLine(23.25,bottom,23.25,top)
l9.SetLineStyle(2)
l9.SetLineColor(ROOT.kBlue)

l10 = ROOT.TLine(23.75,bottom,23.75,top)
l10.SetLineStyle(2)
l10.SetLineColor(ROOT.kGreen)

l11 = ROOT.TLine(26.5,bottom,26.5,top)
l11.SetLineStyle(2)
l11.SetLineColor(ROOT.kGreen)

l12 = ROOT.TLine(29.5,bottom,29.5,top)
l12.SetLineStyle(2)
l12.SetLineColor(ROOT.kGreen)

l13 = ROOT.TLine(32.5,bottom,32.5,top)
l13.SetLineStyle(2)
l13.SetLineColor(ROOT.kGreen)

l14 = ROOT.TLine(35.25,bottom,35.25,top)
l14.SetLineStyle(2)
l14.SetLineColor(ROOT.kGreen)

c_chi2_y = ROOT.TCanvas( "c_chi2_y", "c_chi2_y", 800, 600 )
g_chi2_y_A.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
g_chi2_y_A.SetMaximum(7.0)
g_chi2_y_A.SetMinimum(0.0)

g_chi2_y_A.Draw("AP")
g_chi2_y_J.Draw("P")
g_chi2_y_N.Draw("P")
l1.Draw('same')
l2.Draw('same')
l3.Draw('same')
l4.Draw('same')
l5.Draw('same')
l6.Draw('same')
l7.Draw('same')
l8.Draw('same')
l9.Draw('same')
l10.Draw('same')
l11.Draw('same')
l12.Draw('same')
l13.Draw('same')
l14.Draw('same')
hline1.Draw('same')
c_chi2_y.Print("InfoFromVdMScans/VertexCountingChi2YPerBCID_allScans_NormaliseToScanSets.pdf")

