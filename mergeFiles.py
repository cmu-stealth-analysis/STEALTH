import os
import sys
import glob
import array
import numpy as np
import ROOT
import subprocess

multiple=35
nIter=2

era='2016C'
inDir='161028_125619/0000'
eosDir='store/user/mandrews/SinglePhoton/Run%s'%(era)

procEOS = subprocess.Popen(['source /afs/cern.ch/project/eos/installation/cms/etc/setup.sh; eos ls /%s/%s/ggtree_data_*.root'%(eosDir,inDir)], stdout=subprocess.PIPE, shell=True)
inFiles = procEOS.communicate()[0].split()
print " >> Read in", len(inFiles),"files..."

sw = ROOT.TStopwatch()
sw.Start()

i = 0
while (i < nIter): 

	ggIn = ROOT.TChain("ggNtuplizer/EventTree")

	j = 0
	while (j < multiple):
		iF = j + (i*multiple)
		inName = "root://cms-xrd-global.cern.ch//%s/%s/%s"%(eosDir,inDir,inFiles[iF])
		print " >> Adding #(%02d,%03d): %s"%(j+1,iF+1,inName)
		ggIn.Add(inName)
		j += 1
		if (iF+1 == len(inFiles)):
			break

	nEntries = ggIn.GetEntries()
	outName = "ggNtuple_SinglePhoton_Run%s_%d.root"%(era,i+1)
	print " >> nEntries:",nEntries
	print " >> Writing to",outName

	outFile = ROOT.TFile(outName, "RECREATE")
	outDir = outFile.mkdir("ggNtuplizer")
	outDir.cd()
	ggOut = ggIn.CloneTree()
	outFile.Write()
	outFile.Close()

	i += 1

sw.Stop()
print "Real time: " + str(sw.RealTime() / 60.0) + " minutes"
print "CPU time:  " + str(sw.CpuTime() / 60.0) + " minutes"
