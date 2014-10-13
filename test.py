#! /usr/bin/env python

import os
import sys


run = int(sys.argv[1])

os.system("/afs/cern.ch/project/eos/installation/atlas/bin/eos.select ls /eos/atlas//atlascerngroupdisk/perf-lumi/NTUP_IDVTXLUMI/8TeV_PixelBeam>listofexistingD3PDs.txt")
os.system("cp listofexistingD3PDs.txt listofexistingruns.txt;sed -i 's/data12_8TeV.00//g' listofexistingruns.txt;sed -i 's/\..[^.]*//g' listofexistingruns.txt")

runs = []
datasets = []
for line in open('listofexistingD3PDs.txt','r'):
	datasets.append(str(line))
for line in open('listofexistingruns.txt','r'):
	runs.append(int(line))
print run
print runs
ind = runs.index(run)
command = "/afs/cern.ch/project/eos/installation/atlas/bin/eos.select ls /eos/atlas//atlascerngroupdisk/perf-lumi/NTUP_IDVTXLUMI/8TeV_PixelBeam/"+datasets[ind]+" > log.txt"

print ind
print command

