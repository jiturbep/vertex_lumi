#! /usr/bin/env python

scans = [1, 2, 3, 4, 5, 6, 8, 10, 11, 14, 15]
alg = "vtx_counting"
ntrkcuts = ["3","4","5"]
for ntrkcut in ntrkcuts:
	filesigma = open("InfoFromVdMScans/"+alg+"_SigmaVis_PerScan_NTrkCut"+ntrkcut+".txt","w")
	filesigma.write(alg+"\n")
	filesigma.write("{0} {1} {2}\n".format(str().ljust(5),str("sigma_vis").ljust(15),str("sigma_vis_err").ljust(15)))
	avg_sigma_vis = []
	avg_sigma_vis_err = []
	for scan in scans:
		#sigma_vis averaged over BCIDS
		fname1 = "InfoFromVdMScans/"+alg+"/results_sigma_vis_value_scan"+str(scan)+"_NTrkCut"+ntrkcut+".txt"
		with open(fname1) as f1:
			for line in f1:
				sigma_vis,sigma_vis_err = line.split()
				avg_sigma_vis.append(float(sigma_vis))
				avg_sigma_vis_err.append(float(sigma_vis_err))
		#f3.close()
		filesigma.write("{0} {1} {2}\n".format(str(scan).ljust(5),str(sigma_vis).ljust(15),str(sigma_vis_err).ljust(15)))

	filesigmaperbcid = open("InfoFromVdMScans/"+alg+"_SigmaVis_PerBCID_NTrkCut"+ntrkcut+".txt","w")
	filesigmaperbcid.write(alg+"\n")
	filesigmaperbcid.write("{0} {1} {2} {3}\n".format(str("Scan").ljust(5), str("bcid").ljust(5),str("sigma_vis").ljust(15),str("sigma_vis_err").ljust(15)))
	for scan in scans:			
		bcids = []
		sigma_vis = []
		sigma_vis_err = []
		fname = "InfoFromVdMScans/"+alg+"/results_parameters_scan"+str(scan)+"_NTrkCut"+ntrkcut+".txt"
		for line in open(fname, "r"):
			bcid_n,sigma_v,sigma_v_err,Sigma_x,Sigma_x_err,\
			Sigma_y,Sigma_y_err,Peak_mu_x,Peak_mu_x_err,\
			Peak_mu_y,Peak_mu_y_err,chi2x,chi2y,l_sp,l_sp_err,r_x_i,r_y_i,c_x,c_x_err,c_y,c_y_err = line.split()
			bcids.append(int(bcid_n))
			sigma_vis.append(float(sigma_v))
			sigma_vis_err.append(float(sigma_v_err))
		for a,b,c in zip(bcids,sigma_vis,sigma_vis_err):
			filesigmaperbcid.write("{0} {1} {2} {3}\n".format(str(scan).ljust(5),str(a).ljust(5),str(b).ljust(15),str(c).ljust(15)))

	# Lumi sp
	filelumisp = open("InfoFromVdMScans/"+alg+"_LumiSp_PerScan_NTrkCut"+ntrkcut+".txt","w")
	filelumisp.write(alg+"\n")
	filelumisp.write("{0} {1} {2}\n".format(str().ljust(5),str("LumiSp").ljust(15),str("LumiSp_err").ljust(15)))
	avg_lumi_sp = []
	avg_lumi_sp_err = []
	for scan in scans:
		#lumi_sp averaged over BCIDS
		fname2 = "InfoFromVdMScans/"+alg+"/results_lumi_sp_value_scan"+str(scan)+"_NTrkCut"+ntrkcut+".txt"
		with open(fname2) as f2:
			for line in f2:
				lumi_sp,lumi_sp_err = line.split()
				avg_lumi_sp.append(float(lumi_sp))
				avg_lumi_sp_err.append(float(lumi_sp_err))
		#f3.close()
		filelumisp.write("{0} {1} {2}\n".format(str(scan).ljust(5),str(lumi_sp).ljust(15),str(lumi_sp_err).ljust(15)))

	filelumiperbcid = open("InfoFromVdMScans/"+alg+"_LumiSp_PerBCID_NTrkCut"+ntrkcut+".txt","w")
	filelumiperbcid.write(alg+"\n")
	filelumiperbcid.write("{0} {1} {2}\n".format(str("Scan").ljust(5),str("bcid").ljust(5),str("LumiSp").ljust(15),str("LumiSp_err").ljust(15)))
	for scan in scans:	
		bcids = []
		lumi_sp = []
		lumi_sp_err = []
		fname = "InfoFromVdMScans/"+alg+"/results_parameters_scan"+str(scan)+"_NTrkCut"+ntrkcut+".txt"
		for line in open(fname, "r"):
			bcid_n,sigma_v,sigma_v_err,Sigma_x,Sigma_x_err,\
			Sigma_y,Sigma_y_err,Peak_mu_x,Peak_mu_x_err,\
			Peak_mu_y,Peak_mu_y_err,chi2x,chi2y,l_sp,l_sp_err,r_x_i,r_y_i,c_x,c_x_err,c_y,c_y_err = line.split()
			bcids.append(int(bcid_n))
			lumi_sp.append(float(l_sp))
			lumi_sp_err.append(float(l_sp_err))
		for a,b,c in zip(bcids,lumi_sp,lumi_sp_err):
			filelumiperbcid.write("{0} {1} {2} {3}\n".format(str(scan).ljust(5),str(a).ljust(5),str(b).ljust(15),str(c).ljust(15)))
