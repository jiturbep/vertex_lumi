#! /usr/bin/env python

import os
import sys

run = int(sys.argv[1])

#os.system("/afs/cern.ch/project/eos/installation/atlas/bin/eos.select ls /eos/atlas//atlascerngroupdisk/perf-lumi/NTUP_IDVTXLUMI/8TeV_PixelBeam>listofexistingD3PDs.txt")
#os.system("cp listofexistingD3PDs.txt listofexistingruns.txt;sed -i 's/data12_8TeV.00//g' listofexistingruns.txt;sed -i 's/\..[^.]*//g' listofexistingruns.txt")

runs = []
datasets = []
f=open('listofexistingD3PDs.txt','r')
for line in f:
	datasets.append(line.strip())
f.close()
g=open('listofexistingruns.txt','r')	
for line in g:
	runs.append(int(line))
g.close()

if run in runs:
	ind = runs.index(run)
	print "Run", run, ind
	command = "/afs/cern.ch/project/eos/installation/atlas/bin/eos.select ls /eos/atlas//atlascerngroupdisk/perf-lumi/NTUP_IDVTXLUMI/8TeV_PixelBeam/"+datasets[ind]+">rootfiles_"+str(ind)+".txt"
	os.system(command)
	rootfiles = []
	for line in open("rootfiles_"+str(ind)+".txt","r"):
		name = "root://eosatlas///eos/atlas//atlascerngroupdisk/perf-lumi/NTUP_IDVTXLUMI/8TeV_PixelBeam/"+datasets[ind]+"/"+str(line)
		rootfiles.append(name)
	log = open("rootfiles_fullname_"+str(run)+".txt","w")
	for name in rootfiles:
		log.write(name)
else: print "No D3PD for this run."



