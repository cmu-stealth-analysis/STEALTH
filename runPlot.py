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
#for f in glob.glob("SinglePhoton_2016B_sT_Pho100_Jet*.root"):
for f in glob.glob("ggNTUPLES/SinglePhoton_2016B_sT_Pho100NoPxl_JetTight30_Ht700.root"):
#for f in glob.glob("ggNTUPLES/SinglePhoton_2016B_sT_Pho100_JetTight30_Ht700_JetID.root"):
#for f in glob.glob("ggNTUPLES/SinglePhoton_2016B_sT_Pho100_JetTight30_Ht700.root"):
	for stMin in stMins:
		for stMax in stMaxs:
			f = os.path.splitext(f)[0]
			print "Processing "+str(f)
			#subprocess.call("python -i plotST.py -i %s" % f, shell=True)
			subprocess.call("python plotST.py -i %s -l %f -r %f" % (f,stMin,stMax), shell=True)
			#os.rename("sT.png","PLOTS/%s.png"%f)
			#os.rename("sT.png","PLOTS/sT_%d_%d.png"%(stMin,stMax))
