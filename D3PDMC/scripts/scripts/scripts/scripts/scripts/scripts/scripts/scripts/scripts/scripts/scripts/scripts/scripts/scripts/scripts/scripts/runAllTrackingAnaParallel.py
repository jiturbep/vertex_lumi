#!/usr/bin/env python

__author__  = 'David Yu'
__version__ = '$Id: runAllTrackingAnaParallel.py $'
__usage__   = """%prog [options] [cmd ...]

"""

# Script to submit TrackingAna batch jobs in parallel
import os, sys, re
import glob
from array import array
from math import sqrt, floor

def SubmitJobs(p_tag, p_max_chi2ndf, p_runbatch, p_max_events):
	print "SubmitJobs called with p_tag = " + str(p_tag)

	version = "v2"

	log_pattern = re.compile("log.tgz")

	input_directory = ""
	eliza_number = ""
	extra_opt = []
	files = []

	if p_max_chi2ndf > 0.:
		extra_opt.append("-c " + str(p_max_chi2ndf))
		version = "chi2ndf_" + str(p_max_chi2ndf)
	if p_max_events > 0:
		extra_opt.append("-e " + str(p_max_events))
		version = version + "_test"



	if p_tag == "mc_8TeV_17.2_default_pythia8_pu":
		input_directory = "/eliza2/atlas/atlasdata/atlasscratchdisk/user/spagan/valid1/user.spagan.valid1.107499.singlepart_empty.recon.VTXD3PD.e603_s1469_r3560.v1.0.120419110055/*root*"
		eliza_number = "2"
		files = glob.glob(input_directory)
		print str(len(files)) + " files found."
		extra_opt.append("-a")
		files_per_job = 5
	elif p_tag == "mc_8TeV_17.2_VtxLumi_pythia8_pu":
		input_directory = "/eliza18/atlas/data/data12_8TeV/VtxNtuples_17.2.1.3/user.beate.trkReco_default.valid1.107499.singlepart_empty.recon.ESD.e603_s1469_r3560.v11.120430123934_D3PD/*root*"
		eliza_number = "18"
		files = glob.glob(input_directory)
		extra_opt.append("-a")
		files_per_job = 10
	elif p_tag == "mc_7TeV_16.X_hybrid_pythia6_pu":
		input_directory = "/eliza18/atlas/dryu/Luminosity/mc10_7TeV/user.spagan.mc10_7TeV.105000.pythia_minbias_inelastic.merge.VTXD3PD.e723_s932_s946_r2302_r2300.v6.0/*root*"
		eliza_number = "18"
		files = glob.glob(input_directory)
		files_per_job = 4
	elif p_tag == "mc_7TeV_17.2_default_pythia6_pu_noslim":
		files = glob.glob("/eliza18/atlas/spagan/VtxStudies/VtxPU/VTXD3PD/mc_Beate/MC11_7TeV.107499.singlepart_empty.pileup_pythia6_noslim_mu*.VTXD3PD.root*")
		eliza_number = "18"
		extra_opt.append("-a")
		files_per_job = 4
	elif p_tag == "mc_7TeV_17.2_default_pythia8DL_pu":
		files = glob.glob("/eliza18/atlas/spagan/VtxStudies/VtxPU/VTXD3PD/mc_Beate/MC11_7TeV.107499.singlepart_empty.pileup_Pythia8_A2M_DL_mu*.VTXD3PD.root*")
		eliza_number = "18"
		extra_opt.append("-a")
		files_per_job = 4
	elif p_tag == "mc_7TeV_17.2_default_pythia8DL_pu_noslim":
		files = glob.glob("/eliza18/atlas/spagan/VtxStudies/VtxPU/VTXD3PD/mc_Beate/MC11_7TeV.107499.singlepart_empty.pileup_Pythia8_A2M_DL_noslim_mu*.VTXD3PD.root*")
		eliza_number = "18"
		extra_opt.append("-a")
		files_per_job = 1
	elif p_tag == "mc_7TeV_17.2_default_pythia8_pu_loosetruth":
		files = glob.glob("/eliza18/atlas/spagan/VtxStudies/VtxPU/VTXD3PD/mc_Beate/LooseTruth/*mu[1-2]..LooseTruth*root*")
		eliza_number = "18"
		extra_opt.append("-a")
		files_per_job = 4
	elif p_tag == "mc_7TeV_17.2_default_pythia8_pu":
		files = glob.glob("/eliza18/atlas/spagan/VtxStudies/VtxPU/VTXD3PD/mc_Beate/MC11_7TeV.107499.singlepart_empty.pileup_Pythia8_A2M_noslim_mu*.VTXD3PD*root*")
		eliza_number = "18"
		extra_opt.append("-a")
		files_per_job = 1
	elif p_tag == "mc_7TeV_17.2_default_pythia8_pu_bs45":
		files = glob.glob("/eliza18/atlas/spagan/VtxStudies/VtxPU/VTXD3PD/mcBeate/*narrowBS*.VTXD3PD.root")
		eliza_number = "18"
		extra_opt.append("-a")
		files_per_job = 5
	elif p_tag == "mc_7TeV_17.2_default_pythia8_pu_bs55":
		files = glob.glob("/eliza18/atlas/spagan/VtxStudies/VtxPU/VTXD3PD/mcBeate/*2011BS*root*")
		eliza_number = "18"
		extra_opt.append("-a")
		files_per_job = 5
	else:
		print "ERROR: tag not recognized!"
		return 1

	n_files = len(files)

	result_directory="/eliza18/atlas/dryu/Luminosity/VtxTruthMatchResults/" + p_tag + "/" + version
	os.system("mkdir -pv " + result_directory)

	if p_runbatch == True:

		counter = 0
		job_number = 1
		counter_total = 0

		command_base = "/u/dryu/Luminosity/Code/D3PDMC/bin/TrackingAna -b -q"
		for opt in extra_opt:
			command_base = command_base + " " + opt
		command_base = command_base + " -o " + result_directory + "/InDetTrackD3PD_results_j"

		command = command_base + str(job_number)

		for filename in files:
			if log_pattern.search(filename) != None:
				continue

			counter = counter + 1
			counter_total = counter_total + 1
			command = command + " " + filename

			if counter == files_per_job or counter_total == n_files:
				counter = 0
				script_path = "tmp_" + p_tag + "_j" + str(job_number) + ".sh"
				script = open(script_path, 'w')
				script.write("#!/bin/sh\n")
				script.write(command + "\n")
				#print command
				script.close()

				submit_command = "qsub -V -l eliza" + eliza_number + "io=1 " + script_path
				submit_command
				os.system(submit_command)
				#os.system("cat " + script_path)
				#os.system("rm " + script_path)

				job_number = job_number + 1
				command = command_base + str(job_number)

				print "\n"
	else:
		command = "/u/dryu/Luminosity/Code/D3PDMC/bin/TrackingAna -b -q"
		for opt in extra_opt:
			command = command + " " + opt
		command = command + " -o " + result_directory + "/InDetTrackD3PD_results"
		for filename in files:
			command = command + " " + filename
		os.system(command)


if __name__ == '__main__':
	#
	# Save properly quoted string of original command
	#
	qargv = [ ]
	for s in sys.argv:
	    if re.search('\s|\*|\?',s):   # any white space or special characters in word so we need quoting?
	        if "'" in s:
	            qargv.append('"%s"' % re.sub('"',"'",s))
	        else:
	            qargv.append("'%s'" % re.sub("'",'"',s))
	    else:
	        qargv.append(s)
	qcmd = ' '.join(qargv)


	#
	# Argument parsing
	#
	from optparse import OptionParser
	parser = OptionParser(usage=__usage__, version=__version__)
	parser.add_option('-t', '--tag', dest='tag', default=None, help='tag: e.g. mc_8TeV_17.2_default_pythia8_pu')
	parser.add_option('-c', '--maxchi2ndf', dest='max_chi2ndf', default=-1, help='maxchi2ndf: set maximum chi2/ndf for vertices')
	parser.add_option('-a', '--batch', dest='run_batch', action='store_true', default=False, help='--batch: Run with qsub')
	parser.add_option('-e', '--maxEvents', dest='max_events', default=-1, help="maxEvents: maximum number of events to process")
	(options,args) = parser.parse_args()

	SubmitJobs(options.tag, options.max_chi2ndf, options.run_batch, options.max_events)
