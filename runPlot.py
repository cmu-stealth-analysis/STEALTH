#!/usr/bin/env python

# Run ST plotter

import os
import sys 
import array
#import ROOT
import glob
import subprocess

#stMins = [750,1000,1250,1500]
#stMaxs = [3250,3500,3750]
stMins = [1000]
stMaxs = [3500]
#nBins = [4,5,6]
nBins = [5]
#for f in glob.glob("SinglePhoton_2016B_sT_Pho100_Jet*.root"):
#for f in glob.glob("ggNTUPLES/SinglePhoton_2016B_sT_Pho100NoPxl_JetTight30_Ht700.root"):
#for f in glob.glob("ggNTUPLES/HLTMatched/SinglePhoton_2016B_sT_Pho100MedEleVeto_Jet30TightHt700_El15Tight_Mu15Tight.root"):
for f in glob.glob("ggNTUPLES/HLTMatched/SinglePhoton_2016BCD_sT_Pho100MedEleVeto_Jet50MedHt700_El15Tight_Mu15Tight.root"):
#for f in glob.glob("ggNTUPLES/HLTMatched/SinglePhoton_2016B_sT_Pho100*.root"):
#for f in glob.glob("ggNTUPLES/SinglePhoton_2016B_sT_Pho100_JetTight30_Ht700_JetID.root"):
#for f in glob.glob("ggNTUPLES/SinglePhoton_2016B_sT_Pho100_JetTight30_Ht700.root"):
	fDir = os.path.splitext(f)[0]
	print "Processing "+str(fDir)
	fName = os.path.basename(fDir)
	for stMin in stMins:
		for stMax in stMaxs:
			for nBin in nBins:
				#f = os.path.splitext(f)[0]
				print " >> Processing St=[%f,%f] nBins=%d..." % (stMin,stMax,nBin)
				#subprocess.call("python -i plotST.py -i %s" % f, shell=True)
				subprocess.call("python plotST.py -i %s -l %f -r %f -b %d" % (fDir,stMin,stMax,nBin), shell=True)
				os.rename("sT.png","PLOTS/St_%s.png"%fName)
				#os.rename("sT.png","PLOTS/Pho100MedEleVeto_Jet30TightHt700_El15Tight_Mu15Tight/sT_%d_%d_nBins%d.png"%(stMin,stMax,nBin))
