# -*- coding: utf-8 -*-
#! /usr/bin/env python
import ROOT
import math

ROOT.TH1.AddDirectory(False)
ROOT.gROOT.SetBatch()

path = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/Data/MaskingCorrection/"
rootfile = "pmask_cache.root"
samples = ["mc_8TeV_17.2_VtxLumi_BothSamples","data_8TeV_17.2_VtxLumi_215021"]
ntrks = ["3","4","5"]
colors = [ROOT.kRed,ROOT.kBlue,ROOT.kBlack]

for sample in samples:
	print sample
	tgs_pu = []
	tgs_mu = []
	
	pu_xmax = 0
	pu_xmin = 1000
	pu_ymax = 0
	pu_ymin = 1000
	
	mu_xmax = 0
	mu_xmin = 1000
	mu_ymax = 0
	mu_ymin = 1000
	
	legend = ROOT.TLegend(0.7,0.1,0.85,0.25)
	legend.SetFillColor(ROOT.kWhite)

	tgfile = ROOT.TFile(path+sample+"/"+rootfile,"READ")
	for ntrk,color in zip(ntrks,colors):
		if (sample.find("data")==-1):
			tg_pu = tgfile.Get("tg_pileup_correction_NTrk"+ntrk)
			tg_mu = tgfile.Get("tg_mu_obs_vs_mu_actual_NTrk"+ntrk)
		else:
			bcid = "2361"
			tg_pu = tgfile.Get("tg_pileup_correction_BCID"+bcid+"_NTrkCut"+ntrk)
			tg_mu = tgfile.Get("tg_mu_obs_vs_mu_actual_BCID"+bcid+"_NTrkCut"+ntrk)
		tg_pu.SetMarkerColor(color)
		tg_mu.SetMarkerColor(color)
		print "Number of points =", tg_pu.GetN()
		print "Evaluate tg_pu at 0",tg_pu.Eval(0)," and at 0.1", tg_pu.Eval(0.1)

		legend.AddEntry(tg_pu,"NTrk "+ntrk, "P")

		if (tg_pu.GetXaxis().GetXmax()>pu_xmax): pu_xmax = tg_pu.GetXaxis().GetXmax()
		if (tg_pu.GetYaxis().GetXmax()>pu_ymax): pu_ymax = tg_pu.GetYaxis().GetXmax()
		if (tg_pu.GetXaxis().GetXmin()<pu_xmin): pu_xmin = tg_pu.GetXaxis().GetXmin()
		if (tg_pu.GetYaxis().GetXmin()<pu_ymin): pu_ymin = tg_pu.GetYaxis().GetXmin()

		if (tg_mu.GetXaxis().GetXmax()>mu_xmax): mu_xmax = tg_mu.GetXaxis().GetXmax()
		if (tg_mu.GetYaxis().GetXmax()>mu_ymax): mu_ymax = tg_mu.GetYaxis().GetXmax()
		if (tg_mu.GetXaxis().GetXmin()<mu_xmin): mu_xmin = tg_mu.GetXaxis().GetXmin()
		if (tg_mu.GetYaxis().GetXmin()<mu_ymin): mu_ymin = tg_pu.GetYaxis().GetXmin()


		tgs_pu.append(tg_pu)
		tgs_mu.append(tg_mu)

	tgs_mu[0].SetTitle("")
	tgs_mu[0].GetXaxis().SetLabelSize(0.035)
	tgs_mu[0].GetXaxis().SetTitleSize(0.035)
	tgs_mu[0].GetXaxis().SetTitle("#mu")
	tgs_mu[0].GetXaxis().SetTitleOffset(1.)
	tgs_mu[0].GetXaxis().SetRangeUser(mu_xmin,mu_xmax)
	tgs_mu[0].GetYaxis().SetLabelSize(0.035)
	tgs_mu[0].GetYaxis().SetTitleSize(0.035)
	tgs_mu[0].GetYaxis().SetTitle("#mu_{obs}")
	tgs_mu[0].GetYaxis().SetTitleOffset(1.)
	tgs_mu[0].GetYaxis().SetRangeUser(mu_ymin,mu_ymax)

	tgs_pu[0].SetTitle("")
	tgs_pu[0].GetXaxis().SetLabelSize(0.035)
	tgs_pu[0].GetXaxis().SetTitleSize(0.035)
	tgs_pu[0].GetXaxis().SetTitle("#mu_{obs}")
	tgs_pu[0].GetXaxis().SetTitleOffset(1.)
	tgs_pu[0].GetXaxis().SetRangeUser(pu_xmin,pu_xmax)
	tgs_pu[0].GetYaxis().SetLabelSize(0.035)
	tgs_pu[0].GetYaxis().SetTitleSize(0.035)
	tgs_pu[0].GetYaxis().SetTitle("#frac{#mu}{#mu_{obs}}")
	tgs_pu[0].GetYaxis().SetTitleOffset(1.)
	tgs_pu[0].GetYaxis().SetRangeUser(pu_ymin,pu_ymax)

	c_mu = ROOT.TCanvas("c_mu","c_mu",1200,800)	
	tgs_mu[0].Draw("AP")
	for tg_mu in tgs_mu:
		tg_mu.Draw("P")
	legend.Draw("same")
	c_mu.Print("tg_muobs_muactual_"+sample+".pdf")
		
	c_pu = ROOT.TCanvas("c_pu","c_pu",1200,800)	
	tgs_pu[0].Draw("AP")
	for tg_pu in tgs_pu:
		tg_pu.Draw("P")
	legend.Draw("same")
	c_pu.Print("tg_pileup_correction_"+sample+".pdf")



 
#tgs_pileup = []
#pu_xmax = 0
#pu_xmin = 1000
#pu_ymax = 0
#pu_ymin = 1000#

#hs_residuals = []
#res_ymax = 0
#res_ymin = 1000#

#hs_deltaz = []#

#legend = ROOT.TLegend(0.7,0.1,0.85,0.25)
#legend.SetFillColor(ROOT.kWhite)
#for rootfile,tag,color in zip(rootfiles,tags,colors):
#	tgfile = ROOT.TFile(path+rootfile,"READ")
#	resfile = ROOT.TFile(path+rootfile+"residuals.root","READ")
#	tg_pileup = tgfile.Get("tg_pileup_correction_NTrk5")
#	h_res = resfile.Get("h_dz_rebinned_residuals_NTrk5")
#	h_res.Sumw2()
#	hist_deltaz = tgfile.Get("h_dz_rebinned_NTrk5")
#	hist_deltaz.Sumw2() 
#	fit = tgfile.Get("f_dz_excluded_NTrk5")
#	print "tg_pileup"
#	print tg_pileup.GetXaxis().GetXmax(), tg_pileup.GetXaxis().GetXmin(), tg_pileup.GetYaxis().GetXmax(), tg_pileup.GetYaxis().GetXmin() #

#	if (tg_pileup.GetXaxis().GetXmax()>pu_xmax): pu_xmax = tg_pileup.GetXaxis().GetXmax()
#	if (tg_pileup.GetYaxis().GetXmax()>pu_ymax): pu_ymax = tg_pileup.GetYaxis().GetXmax()
#	if (tg_pileup.GetXaxis().GetXmin()<pu_xmin): pu_xmin = tg_pileup.GetXaxis().GetXmin()
#	if (tg_pileup.GetYaxis().GetXmin()<pu_ymin): pu_ymin = tg_pileup.GetYaxis().GetXmin()#

#	tg_pileup.SetMarkerColor(color)
#	tgs_pileup.append(tg_pileup)#

#	print "h_res"
#	print h_res.GetYaxis().GetXmax(), h_res.GetYaxis().GetXmin()#

#	if (h_res.GetYaxis().GetXmax()>res_ymax): res_ymax = h_res.GetYaxis().GetXmax()
#	if (h_res.GetYaxis().GetXmin()<res_ymin): res_ymin = h_res.GetYaxis().GetXmin()#

#	h_res.SetLineColor(color)
#	hs_residuals.append(h_res)#

#	legend.AddEntry(tg_pileup,tag,"P")
#	
#	canvas = ROOT.TCanvas("deltaz","deltaz",1200,800)
#	hist_deltaz.Draw()
#	fit.Draw("same")
#	canvas.Print("ComparingMaskingCorrectionMethods/DeltaZ_"+tag+".pdf")
#	hs_deltaz.append(hist_deltaz)#
#

#print "at the end "
#print pu_xmax, pu_xmin, pu_ymax, pu_ymin 
#print res_ymax, res_ymin
#canvas_pu = ROOT.TCanvas("pu_canvas", "pu_canvas", 1200, 800)
#tgs_pileup[5].GetXaxis().SetRangeUser(pu_xmin,pu_xmax)
#tgs_pileup[5].GetYaxis().SetRangeUser(pu_ymin,pu_ymax)
#tgs_pileup[5].GetXaxis().SetTitle("#mu_{obs}")
#tgs_pileup[5].GetXaxis().SetTitleOffset(0.8)
#tgs_pileup[5].GetYaxis().SetTitle("Masking Correction Factor")
#tgs_pileup[5].GetYaxis().SetTitleOffset(0.8)
#tgs_pileup[5].Draw("AP")
#for gp in tgs_pileup:
#	gp.Draw("P")
#legend.Draw("same")
#canvas_pu.Print("ComparingMaskingCorrectionMethods/tgs_pileup_correction_NTrk5.pdf")#

#canvas_res = ROOT.TCanvas("res_canvas", "res_canvas", 1200, 800)
#hs_residuals[0].GetXaxis().SetRangeUser(-100,100)
#hs_residuals[0].GetYaxis().SetRangeUser(-12, 12)
#hs_residuals[0].GetXaxis().SetTitle("delta z")
#hs_residuals[0].GetXaxis().SetTitleOffset(1)
#hs_residuals[0].GetYaxis().SetTitle("residuals")
#hs_residuals[0].GetYaxis().SetTitleOffset(1)
#hs_residuals[0].Draw("")
#for gp in hs_residuals:
#	gp.Draw("same")
#legend.Draw("same")
#canvas_res.Print("ComparingMaskingCorrectionMethods/hs_residuals_NTrk5.pdf")#

#for h,tag in zip(hs_deltaz,tags):
#	print "///////////////////////////////", tag, "////////////////////////"
#	count = 0
#	count0 = 0
#	h_ratio = ROOT.TH1D("h_ratio","h_ratio",h.GetNbinsX(),h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax())
#	for i in range(h.GetNbinsX()):
#		n = h.GetBinContent(i+1)
#		nerr = h.GetBinError(i+1)
#		d = hs_deltaz[0].GetBinContent(i+1)
#		derr = hs_deltaz[0].GetBinError(i+1)
#		if (n>0 and d>0):
#			ratio = n/d
#			ratioerr = ratio*math.sqrt(pow(nerr/n,2)+pow(derr/d,2))
#			count +=1
#		else:
#			ratio = 0
#			ratioerr = 0
#			count0+=1
#		print "dz=",h.GetBinCenter(i+1),"n=",n,"nerr=",nerr,"d",d,"derr=",derr,"ratio=",ratio,"ratioerr=",ratioerr
#		h_ratio.SetBinContent(i+1,ratio)
#		h_ratio.SetBinError(i+1,ratioerr)
#	print count, count0
#	canvas_ratio = ROOT.TCanvas("ratio","ratio",1200,800) 
#	h_ratio.SetTitle("")
#	h_ratio.GetYaxis().SetTitle("Ratio of "+tag+" to 2to10")
#	h_ratio.GetYaxis().SetTitleOffset(1)
#	h_ratio.GetXaxis().SetTitle("delta z")
#	h_ratio.GetXaxis().SetTitleOffset(1)
#	h_ratio.GetXaxis().SetRangeUser(-200,200)
#	line = ROOT.TLine(-200, 1, 200,1 )
#	line.SetLineColor(ROOT.kRed)
#	h_ratio.Draw()
#	line.Draw("same")
#	canvas_ratio.Print("ComparingMaskingCorrectionMethods/ratio_"+tag+"_TO_2to10.pdf")#

#for h,tag in zip(hs_deltaz,tags):
#	canvas_ratio = ROOT.TCanvas("ratio","ratio",1200,800)
#	hist = h.Clone() 
#	hist.Divide(h,hs_deltaz[2],1.,1.,"B")
#	hist.GetYaxis().SetTitle("Ratio of "+tag+" to 16to30")
#	hist.GetYaxis().SetTitleOffset(0.8)
#	hist.GetXaxis().SetTitle("delta z")
#	hist.GetXaxis().SetTitleOffset(0.8)
#	hist.Draw()
#	canvas_ratio.Print("ComparingMaskingCorrectionMethods/ratio_"+tag+"_TO_16to30_DivideB.pdf")


	





