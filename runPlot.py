#!/usr/bin/env python

# Run ST plotter

import os
import sys 
import array
#import ROOT
import glob
import subprocess

#for f in glob.glob("SinglePhoton_2016B_sT_Pho100_Jet*.root"):
for f in glob.glob("ggNTUPLES/SinglePhoton_2016B_sT_Pho100_JetTight30_Ht700_JetID.root"):
	f = os.path.splitext(f)[0]
	print "Processing "+str(f)
	subprocess.call("python -i plotST.py -i %s" % f, shell=True)
	#os.rename("sT.png","PLOTS/%s.png"%f)
