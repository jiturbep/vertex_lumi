# -*- coding: utf-8 -*-
#! /usr/bin/env python
import ROOT
import math 
import os
import sys

scansets = [201351,207216,207219,214984,215021]
coords = ["x","y"]
ntrk = 5
for scanset in scansets:
  if( scanset == 201351 ):
    scans = [1,2,3]
    BCIDs = [1,241,2881,3121]
  if( scanset == 207216 ):
    scans = [4,5,6]
    BCIDs=[1,721,1821]
  if( scanset == 207219 ):
    scans = [8]
    BCIDs=[1,721,1821]
  if( scanset == 214984 ):
    scans = [10,11,14]
    BCIDs=[1,2361,2881]
  if( scanset == 215021 ):
    scans = [15]
    BCIDs=[1,2361,2881]
  rootfile = "/afs/cern.ch/work/j/jiturbep/private/atlas-lumi/Data/VdMCalibration/VdMScan-"+str(scanset)+"/17.2-VtxLumi/NVtx/vdm_results.root"
  f_scans = ROOT.TFile(rootfile, "READ")
  for scan in scans:
    for bcid in BCIDs:
      for coord in coords:
        graph = f_scans.Get( "tg_"+coord+"_BCID"+str(bcid)+"_scan"+str(scan)+"_NTrkCut"+str(ntrk) )
        for i in range(0,graph.GetN()):
          x,y,xerr,yerr= ROOT.Double(0),ROOT.Double(0),ROOT.Double(0),ROOT.Double(0)
          graph.GetPoint(i,x,y)
          yerr = graph.GetErrorY(i)
          xerr = graph.GetErrorX(i)
          if(scan==8):
            if( i==8 ):
              print scan, bcid, coord, y, yerr
          else:
            if( i ==12 ):
              print scan, bcid, coord, y, yerr
