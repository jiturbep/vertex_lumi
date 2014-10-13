#! /usr/bin/env python

import os
import sys

vdmscans = ["April","July","November"]
for vdmscan in vdmscans:
	files = []
	for i in os.listdir("FittedVdMScans/"+vdmscan):
		files.append(i)
	command = "cd FittedVdMScans/"+vdmscan+";"
	command += "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="+vdmscan+".pdf "
	for f in files:
		command += f+" "
	print command
	#os.system(command)

