#! /usr/bin/env python
from ROOT import *
gROOT.SetBatch(True)

#listoffiles = ["DeltaZDistributions/Data_207216/dz_distribution.root",
			   #"DeltaZDistributions/MC/deltaz_NTrkCut3.root",
			   #"DeltaZDistributions/MC/deltaz_NTrkCut4.root",
			   #"DeltaZDistributions/MC/deltaz_NTrkCut5.root",
			   #"DeltaZDistributions/MC/deltaz_NTrkCut6.root"]
			   #,"DeltaZDistributions/MC/deltaz_NTrkCut7.root",
			   #"DeltaZDistributions/MC/deltaz_NTrkCut8.root",
			   #"DeltaZDistributions/MC/deltaz_NTrkCut10.root"]
datafiles = ["DeltaZDistributions/dz_distribution.root"]
#datafiles = ["DeltaZDistributions/Data_207216/dz_distribution.root"]
mcfiles = ["DeltaZDistributions/deltaz_NTrkCut3.root",
		   "DeltaZDistributions/deltaz_NTrkCut4.root",
		   "DeltaZDistributions/deltaz_NTrkCut5.root",
	 	   "DeltaZDistributions/deltaz_NTrkCut6.root"]
mcfiles = ["DeltaZDistributions/deltaz_NTrkCut5.root", "DeltaZDistributions/z_distribution_NTrkCut5.root"]

data_histonames = []
mc_histonames = []
data_histonames_bcidblind =[]

def CreateHistogramsPNG():
	print 'Number of mc files:', len(mcfiles)
	canvas = TCanvas('canvas','canvas',800,600)
	for filename in mcfiles:
		names = []
		rootfile = TFile(filename,'READ')
		listofhistokeys = rootfile.GetListOfKeys()
		print 'For file:',mcfiles.index(filename),', Number of histograms:',len(listofhistokeys)
		for i in range(len(listofhistokeys)):
			mc_histonames.append(listofhistokeys[i].GetName())
			names.append(listofhistokeys[i].GetName())
		mc_histograms = []
		for name in names:
			histo = rootfile.Get(name)
			#histo.Rebin(10)
			histo.GetXaxis().SetRangeUser(-100.0,100.0)
			#histo.Scale(1/histo.Integral())
			histo.GetXaxis().SetTitleOffset(1)
			histo.GetYaxis().SetTitleOffset(1)
			mc_histograms.append(histo)
			histo.Draw('ehist')
			canvas.Print("DeltaZDistributions/"+name+".pdf")

	print 'Number of data files:', len(datafiles)
	for filename in datafiles:
		rootfile = TFile(filename,'READ')
		listofhistokeys = rootfile.GetListOfKeys()
		print 'For file:',datafiles.index(filename),', Number of histograms:',len(listofhistokeys)
		for i in range(len(listofhistokeys)):
			data_histonames.append(listofhistokeys[i].GetName())
		data_histograms = []
		for name in data_histonames:
			histo = rootfile.Get(name)
			#histo.Rebin(4)
			histo.GetXaxis().SetRangeUser(-100.0,100.0)
			#histo.Scale(1/histo.Integral())
			histo.GetXaxis().SetTitleOffset(1)
			histo.GetYaxis().SetTitleOffset(1)
			data_histograms.append(histo)
			histo.Draw('ehist')
			canvas.Print("DeltaZDistributions/"+name+".pdf")
		ntrkcuts = [2,3,4,5,6,7,8,10]
		for i,ntrkcut in zip(range(0,len(data_histonames),3),ntrkcuts):
			histo = rootfile.Get(data_histonames[i])
			#print "data_histonames[i]", data_histonames[i]
			for j in range(i+1,i+3):
				histo.Add(rootfile.Get(data_histonames[j]))
				#print data_histonames[j]
			histo.Scale(1.0/3.0)
			histo.SetTitle("h_dz_NTrkCut"+str(ntrkcut))
			histo.Draw('ehist')
			canvas.Print("DeltaZDistributions/h_dz_NTrkCut"+str(ntrkcut)+".pdf")
			data_histonames_bcidblind.append("h_dz_NTrkCut"+str(ntrkcut)+".pdf")
			
def CreateLatexFile():
	print "Creating LaTeX file"
	latfile = open("DeltaZDistributions/DeltaZDistributions.tex", "w")
	latfile.write(r"\documentclass[8pt]{beamer}""\n")
	latfile.write(r"\usepackage{graphicx}""\n")
	latfile.write(r"\usepackage{xcolor}""\n")
	latfile.write(r"\usetheme{default}""\n")
	latfile.write(r"\begin{document}""\n")

	latfile.write(r"\begin{frame}""\n")
	latfile.write(r"\centering""\n")
	latfile.write(r"\large{DeltaZ distributions for data, VdM Scan 207216}""\n")
	latfile.write(r"\end{frame}""\n")


	print "len(data_histonames)", len(data_histonames)
	for ihist in range(len(data_histonames)-3):
		latfile.write(r"\begin{frame}""\n")
		latfile.write(r"\center""\n")
		latfile.write(r"\includegraphics[scale=0.5]{"+data_histonames[ihist+3]+"}""\n")
		latfile.write("\end{frame}""\n")

	latfile.write(r"\begin{frame}""\n")
	latfile.write(r"\centering""\n")
	latfile.write(r"\large{BCID-blind DeltaZ distributions for data, VdM Scan 207216}""\n")
	latfile.write(r"\end{frame}""\n")


	print "len(data_histonames_bcidblind)", len(data_histonames_bcidblind)
	for ihist in range(len(data_histonames_bcidblind)-1):
		latfile.write(r"\begin{frame}""\n")
		latfile.write(r"\center""\n")
		latfile.write(r"\includegraphics[scale=0.5]{"+data_histonames_bcidblind[ihist+1]+"}""\n")
		latfile.write("\end{frame}""\n")

	latfile.write(r"\begin{frame}""\n")
	latfile.write(r"\centering""\n")
	latfile.write(r"\large{DeltaZ distributions for MC}""\n")
	latfile.write(r"\newline""\n")
	latfile.write(r"You can see in the next plots that for low ngenint the statistics are not good and the distributions don't look like the ones we see in data. It is until ngenint 6 or 7 that the distributions start to look like data. Therefore, I used mu=5  to calculate the poisson factor used to scale the deltaz vs ngenint distributions, as oppose to 0.5 as I was doing before, when calculating the masking correction.""\n") 
	latfile.write(r"\end{frame}""\n")

	ngenints = range(2,81)*3
	ntrkcuts = [3]*len(ngenints)
	ntrkcuts += [4]*len(ngenints)
	ntrkcuts += [5]*len(ngenints)
	print "len(mc_histonames)", len(mc_histonames)
	for ngenint,ntrkcut,hist in zip(ngenints,ntrkcuts,mc_histonames):
		latfile.write(r"\begin{frame}""\n")
		latfile.write(r"\center""\n")
		latfile.write(r"\large{NTrackCut "+str(ntrkcut)+" NGenInt "+str(ngenint)+"}""\n")
		latfile.write(r"\includegraphics[scale=0.5]{"+hist+"}""\n")
		latfile.write("\end{frame}""\n")

	latfile.write(r"\end{document}""\n")
	latfile.close()

CreateHistogramsPNG()
#CreateLatexFile()

		

