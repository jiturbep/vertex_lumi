#!/usr/bin/env python

__author__  = 'David Yu'
__version__ = '$Id: RunLumiVtx.py $'
__usage__   = """%prog [options] [cmd ...]

"""

# Script to submit TrackingAna batch jobs in parallel
import os, sys, re
import glob
from array import array
from math import sqrt, floor

def SubmitJobs(p_run, p_settings, p_vtx_method, p_newFile, p_verbose, p_systematic, p_local):

	known_systematics = ["fake_low", "fake_high", "masking_toy_scaling", "bs-45mm", "bs-55mm"]
	if p_systematic != "":
		item_found = False
		for item in known_systematics:
			if p_systematic == item:
				item_found = True

		if not item_found:
			print "Unknown systematic specified: " + p_systematic
			return 1 

	output_folder = "/u/dryu/Luminosity/Data/VdMCalibration/" + p_run + "/" + p_settings
	if p_systematic != "":
		output_folder = output_folder + "_" + p_systematic
	output_folder = output_folder + "/"

	os.system("mkdir -pv " + output_folder)
	os.system("mkdir -pv " + output_folder + p_vtx_method)

	command = "runVdM -r " + p_run + " -s " + p_settings + " "
	if p_vtx_method == "NEvt":
		command = command + "-e "
	if p_newFile:
		command = command + "-n "
	if p_verbose:
		command = command + "-v "
	if p_systematic != "":
		command = command + "-u " + p_systematic

	logfile = output_folder + "log_" + p_vtx_method + ".txt"

	if p_local == True:
		print "Running VdM analysis locally:"
		print command + " >& " + logfile + "&"
		os.system(command + " >& " + logfile + "&")
	else:
		script_path = "tmp_runVdM_" + p_run + "_" + p_settings + "_" + p_vtx_method + ".sh"
		script = open(script_path, 'w')
		script.write("#!/bin/sh\n")
		script.write(command + " >& " + logfile + " \n")
		script.close()

		print "Submitting batch job:"
		os.system("cat " + script_path)
		os.system("qsub -V -l eliza18io=1 " + script_path)
		os.system("rm " + script_path)


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
	parser.add_option('-r', '--run', dest='run', type='string', default=None, help='run number')
	parser.add_option('-s', '--settings', dest='settings', type='string', default=None, help='reco settings (17.2-normal | 17.2-VtxLumi)')
	parser.add_option('-e', '--vtxMethod', dest='vtx_method', type='string', default='NVtx', help='Vertex or event counting (NVtx or NEvt)')
	parser.add_option('-n', '--newFile', dest='new_file', action='store_true', default=False, help='Clear old ROOT files')
	parser.add_option('-v', '--verbose', dest='verbose', action='store_true', default=False, help='Write out debugging histograms')
	parser.add_option('-u', '--systematic', dest='systematic', type='string', default="", help='Systematic uncertainty flag')
	parser.add_option('-l', '--local', dest='local', action='store_true', default=False, help='Run locally')
	(options,args) = parser.parse_args()

	SubmitJobs(options.run, options.settings, options.vtx_method, options.new_file, options.verbose, options.systematic, options.local)


