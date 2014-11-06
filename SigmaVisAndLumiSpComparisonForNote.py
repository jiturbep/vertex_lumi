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

alg1 = 'vtx_counting'
algs = ['allTracks','vtxTracks','lucidEvtOr','bcmHEvtOr']
#algs = ['vtxTracks','allTracks', 'lucidEvtOr', 'lucidEvtA', 'lucidEvtC', 'bcmHEvtOr', 'bcmHEvtOrC', 'bcmVEvtOr', 'bcmVEvtOrC']

for alg2 in algs:
	print alg2
	bcids = []
	sigma1 = []
	sigma_err1 = []
	lumi_sp1 = []
	lumi_sp_err1 = []
	for scan in scans:
		file1 = "InfoFromVdMScans/vtx_counting/results_parameters_scan"+str(scan)+"_NTrkCut"+str(ntrkcut)+".txt"
		for line in open(file1, 'r'):
			bcid_n,sigma_v,sigma_v_err,Sigma_x,Sigma_x_err,Sigma_y,Sigma_y_err,Peak_mu_x,Peak_mu_x_err,Peak_mu_y,Peak_mu_y_err,chi2x,chi2y,l_sp,l_sp_err,r_x_i,r_y_i,c_x,c_x_err,c_y,c_y_err = line.split()
			bcids.append((bcid_n))
			sigma1.append(float(sigma_v))
			sigma_err1.append(float(sigma_v_err))
			lumi_sp1.append(float(l_sp))
			lumi_sp_err1.append(float(l_sp_err))

	print "len(bcids):", len(bcids)

	sigma2 = []
	sigma_err2 = []
	lumi_sp2 = []
	lumi_sp_err2 = []
	if ( alg2 == 'vtxTracks' or alg2 == 'allTracks'):
		#file2 = 'sigma_vis_and_lumi_sp_trackCounting_'+alg2+'_dgc_fit.csv'
		file2 = 'InfoFromVdMScans/sigma_vis_and_lumi_sp_trackCounting_'+alg2+'_dgc_fit.dat'
		sc = 1000.0
	else: 
		file2 = 'InfoFromVdMScans/'+alg2+'_SigmaVisLumiSpSigmas_PerBCID_AllScans.csv'
		sc = 1.0
	f = open(file2).readlines()
	firstLine = f.pop(0) #removes the first line
	for line in f:
		if ( alg2 == 'vtxTracks' or alg2 == 'allTracks'):
			scan_n, bcid_n, sigma_v, sigma_v_err, lumi_sp, lumi_sp_err, casigma_x, capsigma_y,a,b,c,d,e,f = line.split()
		else:
			scan_n, bcid_n, sigma_v, sigma_v_err, lumi_sp, lumi_sp_err, casigma_x, capsigma_y = line.split()
		if bcid_n != "all":
			sigma2.append(float(sigma_v))
			sigma_err2.append(float(sigma_v_err))
			lumi_sp2.append(float(lumi_sp)*sc)
			lumi_sp_err2.append(float(lumi_sp_err)*sc)

	if len(sigma1)!=len(sigma2): print "Different number of sigma_vis entries alg1:", len(sigma1)," and alg2:", len(sigma2) 
	if len(lumi_sp1)!=len(lumi_sp2): print "Different number of lumi_sp entries alg1:", len(lumi_sp1)," and alg2:", len(lumi_sp2) 


	############# LumiSp Plot

	g_lumi_sp = ROOT.TGraphErrors(36)
	g_lumi_sp.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
	for i in range(len(lumi_sp2)):
		g_lumi_sp.SetPoint( i, i, lumi_sp1[i]/lumi_sp2[i] )
		g_lumi_sp.SetPointError( i, 0.0, math.sqrt( pow(lumi_sp1[i]/lumi_sp2[i],2) * ( pow(lumi_sp_err1[i]/lumi_sp1[i],2) + pow(lumi_sp_err2[i]/lumi_sp2[i],2) ) ) )
	#g_lumi_sp.SetTitle("Comparison of L_{sp} per BCID per Scan")
	g_lumi_sp.SetTitle("")
	g_lumi_sp.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
	for i in range(len(bcids)):
		g_lumi_sp.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_lumi_sp.GetHistogram().LabelsOption("v","X")
	g_lumi_sp.GetXaxis().SetTitle("BCIDs")
	g_lumi_sp.GetXaxis().SetTitleOffset(1.3)
	g_lumi_sp.GetYaxis().SetTitle("L_{sp}(vtx)/L_{sp}("+alg2+")")
	g_lumi_sp.GetYaxis().SetTitleOffset(1.3)
	g_lumi_sp.SetMarkerStyle(8)

	ay = g_lumi_sp.GetHistogram().GetYaxis()
	bottom = ay.GetBinLowEdge(1)
	top = ay.GetBinUpEdge(ay.GetNbins())
	bottom = 0.97
	top = 1.03

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

	c_lumi_sp = ROOT.TCanvas( "c_lumi_sp", "c_lumi_sp", 800, 600 )
	#Pad = c_lumi_sp.cd(1)
	#Pad.SetTicks(1,1)

	g_lumi_sp.SetMaximum(1.03)
	g_lumi_sp.SetMinimum(0.97)

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
	hline1.Draw('same')
	c_lumi_sp.Print("InfoFromVdMScans/LumiSpPerBCID_allScans"+alg1+"_vs_"+alg2+".pdf")


	############# SigmaVis Plot

	g_sigma_vis = ROOT.TGraphErrors(36)
	g_sigma_vis.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
	for i in range(len(sigma1)):
		ratio1 = sigma1[i]/sigma1[30]
		ratio1_err = math.sqrt( pow(sigma1[i]/sigma1[30],2) * ( pow(sigma_err1[i]/sigma1[i],2) + pow(sigma_err1[30]/sigma1[30],2) ) )
		ratio2 = sigma2[i]/sigma2[30]
		ratio2_err = math.sqrt( pow(sigma2[i]/sigma2[30],2) * ( pow(sigma_err2[i]/sigma2[i],2) + pow(sigma_err2[30]/sigma2[30],2) ) )
		ratio = ratio1/ratio2
		ratio_err = math.sqrt( pow(ratio,2) * ( pow( ratio1_err/ratio1,2) + pow( ratio2_err/ratio2,2) ) )

		g_sigma_vis.SetPoint( i, i, ratio)
		g_sigma_vis.SetPointError( i, 0, ratio_err)

	#g_sigma_vis.SetTitle("Comparison of #sigma_{vis} per BCID per Scan")
	g_sigma_vis.SetTitle("")
	g_sigma_vis.GetHistogram().GetXaxis().Set(38,-1.5,36.5)
	for i in range(len(bcids)):
		g_sigma_vis.GetHistogram().GetXaxis().SetBinLabel(i+2,bcids[i])
	g_sigma_vis.GetHistogram().LabelsOption("v","X")
	g_sigma_vis.GetXaxis().SetTitle("BCIDs")
	g_sigma_vis.GetXaxis().SetTitleOffset(1.3)
	g_sigma_vis.GetYaxis().SetTitle("#sigma_{vis}(vtx)/#sigma_{vis}("+alg2+") Norm. to Scan 14")
	g_sigma_vis.GetYaxis().SetTitleOffset(1.3)
	g_sigma_vis.SetMarkerStyle(8)

	ay = g_sigma_vis.GetHistogram().GetYaxis()
	bottom = ay.GetBinLowEdge(1)
	top = ay.GetBinUpEdge(ay.GetNbins())
	bottom = 0.97
	top = 1.03

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

	g_sigma_vis.SetMaximum(1.03)
	g_sigma_vis.SetMinimum(0.97)

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
	c_lumi_sp.Print("InfoFromVdMScans/SigmaVisPerBCID_allScans"+alg1+"_vs_"+alg2+".pdf")

########### Averaged values
print "########### Averaged Values ###########"
scan_names = ["1","2","3","4","5","6","8","10","11","14","15"]
gs_sigma_vis = []
legend = ROOT.TLegend(0.75,0.6,0.9,0.9)
legend.SetFillColor(ROOT.kWhite)
colors = [ROOT.kRed, ROOT.kRed+2, ROOT.kBlue, ROOT.kBlue+2, ROOT.kGreen, ROOT.kGreen+2, ROOT.kCyan+1, ROOT.kCyan+3, ROOT.kMagenta+3]
for alg2,color in zip(algs,colors):
	print alg2
	sigma1 = []
	sigma_err1 = []
	lumi_sp1 = []
	lumi_sp_err1 = []
	file1 = "InfoFromVdMScans/vtx_counting_SigmaVis_PerScan_NTrkCut"+str(ntrkcut)+".txt"
	with open(file1,'r') as f:
		f.next()
		f.next()
		for line in f:
			scan,sigma_v,sigma_v_err = line.split()
			sigma1.append(float(sigma_v))
			sigma_err1.append(float(sigma_v_err))
	filelp1 = "InfoFromVdMScans/vtx_counting_LumiSp_PerScan_NTrkCut"+str(ntrkcut)+".txt"
	with open(filelp1,'r') as f:
		f.next()
		f.next()
		for line in f:
			scan,sigma_v,sigma_v_err = line.split()
			lumi_sp1.append(float(sigma_v))
			lumi_sp_err1.append(float(sigma_v_err))

	sigma2 = []
	sigma_err2 = []
	lumi_sp2 = []
	lumi_sp_err2 = []
	if ( alg2 == 'vtxTracks' or alg2 == 'allTracks'):
		#file2 = 'sigma_vis_and_lumi_sp_trackCounting_'+alg2+'_dgc_fit.csv'
		file2 = 'InfoFromVdMScans/sigma_vis_and_lumi_sp_trackCounting_'+alg2+'_dgc_fit.dat'
		sc = 1000.0
		f = open(file2).readlines()
		firstLine = f.pop(0) #removes the first line
		for line in f:
			scan_n, bcid_n, sigma_v, sigma_v_err, lumi_sp, lumi_sp_err, casigma_x, capsigma_y,a,b,c,d,e,f = line.split()
			if bcid_n == "all":
				sigma2.append(float(sigma_v))
				sigma_err2.append(float(sigma_v_err))
			if bcid_n == "1":
				lumi_sp2.append(float(lumi_sp)*sc)
				lumi_sp_err2.append(float(lumi_sp_err)*sc)
	else: 
		file2 = 'InfoFromVdMScans/'+alg2+'_SigmaVis_AllScans.csv'
		filelp2 = 'InfoFromVdMScans/'+alg2+'_LumiSp_AllScans.csv'
		sc = 1.0
		f = open(file2).readlines()
		firstLine = f.pop(0) #removes the first line
		for line in f:
			scan_n, sigma_v, sigma_v_err = line.split()
			sigma2.append(float(sigma_v))
			sigma_err2.append(float(sigma_v_err))
		flp = open(filelp2).readlines()
		firstLine = flp.pop(0) #removes the first line
		for line in flp:
			scan_n, lumi_sp, lumi_sp_err = line.split()
			lumi_sp2.append(float(lumi_sp)*sc)
			lumi_sp_err2.append(float(lumi_sp_err)*sc)

	if len(sigma1)!=len(sigma2): print "Different number of sigma_vis entries alg1:", len(sigma1)," and alg2:", len(sigma2) 
	if len(lumi_sp1)!=len(lumi_sp2): print "Different number of lumi_sp entries alg1:", len(lumi_sp1)," and alg2:", len(lumi_sp2) 

	############# SigmaVis Plot

	g_sigma_vis = ROOT.TGraphErrors(11)
	g_sigma_vis.GetHistogram().GetXaxis().Set(13,-1.5,11.5)
	for i in range(len(sigma1)):
		ratio1 = sigma1[i]/sigma1[9]
		ratio1_err = math.sqrt( pow(sigma1[i]/sigma1[9],2) * ( pow(sigma_err1[i]/sigma1[i],2) + pow(sigma_err1[9]/sigma1[9],2) ) )
		ratio2 = sigma2[i]/sigma2[9]
		ratio2_err = math.sqrt( pow(sigma2[i]/sigma2[9],2) * ( pow(sigma_err2[i]/sigma2[i],2) + pow(sigma_err2[9]/sigma2[9],2) ) )
		ratio = ratio1/ratio2
		ratio_err = math.sqrt( pow(ratio,2) * ( pow( ratio1_err/ratio1,2) + pow( ratio2_err/ratio2,2) ) )

		g_sigma_vis.SetPoint( i, i, ratio)
		g_sigma_vis.SetPointError( i, 0, ratio_err)

	g_sigma_vis.GetHistogram().GetXaxis().Set(13,-1.5,11.5)
	for i in range(len(scan_names)):
		g_sigma_vis.GetHistogram().GetXaxis().SetBinLabel(i+2,scan_names[i])
	g_sigma_vis.SetTitle("Comparison of #sigma_{vis} per Scan")
	g_sigma_vis.GetXaxis().SetTitle("Scans")
	g_sigma_vis.GetXaxis().SetTitleOffset(1.3)
	g_sigma_vis.GetYaxis().SetTitle("#sigma_{vis}(vtx)/#sigma_{vis}(alg) Norm. to Scan 14")
	g_sigma_vis.GetYaxis().SetTitleOffset(1.3)
	g_sigma_vis.SetMarkerStyle(8)
	g_sigma_vis.SetMarkerColor(color)
	g_sigma_vis.SetLineColor(color)
	gs_sigma_vis.append(g_sigma_vis)
	legend.AddEntry(g_sigma_vis, alg2, "P")

c_lumi_sp = ROOT.TCanvas( "c_sigma", "c_sigma", 800, 600 )
gs_sigma_vis[0].SetMaximum(1.03)
gs_sigma_vis[0].SetMinimum(0.96)
gs_sigma_vis[0].Draw("AP")
for i in range(len(gs_sigma_vis)):
	gs_sigma_vis[i].Draw("P")
legend.Draw("same")
c_lumi_sp.Print("InfoFromVdMScans/SigmaVis_allScans"+alg1+"_vs_otherAlgs.pdf")

########### Averaged values
print "########### Averaged Values NO RATIOS ###########"
scan_names = ["1","2","3","4","5","6","8","10","11","14","15"]
gs_sigma_vis = []
legend = ROOT.TLegend(0.7,0.15,0.85,0.45)
legend.SetFillColor(ROOT.kWhite)
colors = [ROOT.kRed, ROOT.kMagenta, ROOT.kMagenta+2, ROOT.kRed+2, ROOT.kOrange+7, ROOT.kGreen+1, ROOT.kBlue, ROOT.kBlue+2, ROOT.kCyan+1, ROOT.kCyan+2]
#colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2,ROOT.kMagenta+3]
algs_all = ['vtx_counting', 'vtxTracks','allTracks', 'lucidEvtOr', 'lucidEvtA', 'lucidEvtC', 'bcmHEvtOr', 'bcmHEvtOrC', 'bcmVEvtOr', 'bcmVEvtOrC']
#algs_all = ['vtx_counting', 'bcmHEvtOr', 'bcmHEvtOrC', 'bcmVEvtOr', 'bcmVEvtOrC']
for alg,color in zip(algs_all,colors):
	print 'Alg:', alg
	g_sigma_vis = ROOT.TGraphErrors(11)
	g_sigma_vis.GetHistogram().GetXaxis().Set(11,-0.5,10.5)
	sigma_vis = []
	sigma_vis_err = []
	for scan in scans:
		print "Scan", scan
		if alg == 'vtx_counting':
			file1 = "InfoFromVdMScans/vtx_counting/results_sigma_vis_value_scan"+str(scan)+"_NTrkCut"+str(ntrkcut)+".txt"
			for line in open(file1, 'r'):
				sigma_v,sigma_v_err = line.split()
				sigma_vis.append(float(sigma_v))
				sigma_vis_err.append(float(sigma_v_err))
				print sigma_v, sigma_v_err

		elif ( alg == 'vtxTracks' or alg == 'allTracks'):
			#file2 = 'sigma_vis_and_lumi_sp_trackCounting_'+alg+'_dgc_fit.csv'
			file2 = 'InfoFromVdMScans/sigma_vis_and_lumi_sp_trackCounting_'+alg+'_dgc_fit.dat'
			f = open(file2).readlines()
			firstLine = f.pop(0) #removes the first line
			for line in f:
				scan_n, bcid_n, sigma_v, sigma_v_err, lumi_sp, lumi_sp_err, casigma_x, capsigma_y, a, b, c, d, e, f = line.split()
				if bcid_n == "all":
					sigma_vis.append(float(sigma_v))
					sigma_vis_err.append(float(sigma_v_err))
		else: 
			file2 = 'InfoFromVdMScans/'+alg+'_SigmaVis_AllScans.csv'
			f = open(file2).readlines()
			firstLine = f.pop(0) #removes the first line
			for line in f:
				scan_n, sigma_v, sigma_v_err = line.split()
				sigma_vis.append(float(sigma_v))
				sigma_vis_err.append(float(sigma_v_err))
	for i in range(len(sigma_vis)):
		g_sigma_vis.SetPoint(i,i,sigma_vis[i]/sigma_vis[9])
		ratio_err = math.sqrt( pow(sigma_vis[i]/sigma_vis[9],2) * ( pow(sigma_vis_err[i]/sigma_vis[i],2) + pow(sigma_vis_err[9]/sigma_vis[9],2) ) )
		g_sigma_vis.SetPointError(i,0.,ratio_err)
	g_sigma_vis.GetHistogram().GetXaxis().Set(11,-0.5,10.5)
	for i in range(len(scan_names)):
		g_sigma_vis.GetHistogram().GetXaxis().SetBinLabel(i+1,scan_names[i])
	#g_sigma_vis.SetTitle("Comparison of #sigma_{vis} per Scan")
	g_sigma_vis.SetTitle("")
	g_sigma_vis.GetXaxis().SetTitle("Scans")
	g_sigma_vis.GetXaxis().SetTitleOffset(1.3)
	g_sigma_vis.GetYaxis().SetTitle("#sigma_{vis}(alg)/[#sigma_{vis}(alg) from Scan 14]")
	g_sigma_vis.GetYaxis().SetTitleOffset(1.3)
	g_sigma_vis.SetMarkerStyle(8)
	g_sigma_vis.SetMarkerColor(color)
	g_sigma_vis.SetLineColor(color)
	gs_sigma_vis.append(g_sigma_vis)
	legend.AddEntry(g_sigma_vis, alg, "P")

c_lumi_sp = ROOT.TCanvas( "c_sigma", "c_sigma", 800, 600 )
#gs_sigma_vis[0].SetMaximum(1.04)
gs_sigma_vis[0].SetMaximum(1.02) # lucid
#gs_sigma_vis[0].SetMaximum(1.01) #tracking
gs_sigma_vis[0].SetMinimum(0.94)
#gs_sigma_vis[0].SetMinimum(0.95) #tracking
gs_sigma_vis[0].Draw("AP")
for i in range(len(gs_sigma_vis)):
	gs_sigma_vis[i].Draw("P")
legend.Draw("same")
c_lumi_sp.Print("InfoFromVdMScans/SigmaVis_allScans_allAlgs.pdf")
