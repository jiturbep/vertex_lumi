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

def SubmitJobs(p_run, p_settings, p_output_prefix, p_systematic, p_ntrkcut, p_bcid):

	#Check systematic
	known_systematics = ["fake_low", "fake_high", "masking_toy_scaling", "bs-45mm", "bs-55mm"]
	if p_systematic != "":
		item_found = False
		for item in known_systematics:
			if p_systematic == item:
				item_found = True

		if not item_found:
			print "Unknown systematic specified: " + p_systematic
			return 1 

	#Make result directory
	result_directory = "/u/dryu/Luminosity/Data/LumiVtx/" + p_run + "/" + p_settings
	if p_systematic != "":
		result_directory = result_directory + "_" + p_systematic
	os.system("mkdir -pv " + result_directory)

	#Make batch submission script
	command = "GetLuminosity -r " + p_run + " -s " + p_settings
	if p_output_prefix != None:
		command = command + " -o " + p_output_prefix
	if p_systematic != "":
		command = command + " -u " + p_systematic
	if p_ntrkcut != None:
		command = command + " -n " + str(p_ntrkcut)
	if p_bcid != None:
		command = command + " -i " + str(p_bcid)

	logfile = result_directory + "/log.txt"
	command = command + " >& " + logfile

	script_path = "tmp_LumiVertex.sh"
	script = open(script_path, 'w')
	script.write("#!/bin/sh\n")
	script.write(command + "\n")
	script.close()

	print "Submitting batch job:"
	os.system("cat " + script_path)

	#Submit batch job
	os.system("qsub -V -l h_vmem=1G -l eliza18io=1 " + script_path)
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
	parser.add_option('-u', '--systematic', dest='systematic', type='string', default="", help='Systematic uncertainty flag')
	parser.add_option('-o', '--outputPrefix', dest='output_prefix', type='string', default=None, help='manual output prefix')
	parser.add_option('-n', '--nTrkCut', dest='ntrkcut', type='int', default=None, help='Single NTrkCut')
	parser.add_option('-i', '--bcid', dest='bcid', type='int', default=None, help='Single BCID')
	(options,args) = parser.parse_args()

	SubmitJobs(options.run, options.settings, options.output_prefix, options.systematic, options.ntrkcut, options.bcid)


