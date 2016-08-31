#!/usr/bin/env python

# Run ST plotter

import os
import sys 
import array
#import ROOT
import glob
import subprocess

for f in glob.glob("SinglePhoton_2016B_sT_Pho100_Jet*.root"):
	f = os.path.splitext(f)[0]
	print "Processing "+str(f)
	subprocess.call("python plotST.py -i %s" % f, shell=True)
