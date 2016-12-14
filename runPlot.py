#!/usr/bin/env python

import os
import sys 
import glob
import subprocess

#stMins = [750,1000,1250,1500]
#stMaxs = [3250,3500,3750,4000]
stMins = [1000]
#stMins = [1250]
#stMins = [1375]
stMaxs = [3500]
#stMaxs = [3750]
#stMaxs = [4000]
#nBins = [4,5,6]
#nBins = [20]
#nBins = [20,5]
nBins = [5]

#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/DATA/JetHT/JetHT_2016*_sT_Pho30MedEleVeto_Jet30MedHt1000_El15Tight_Mu15Tight.root'
#INPATH = 'stNTUPLES/DATA/JetHT_Run2016*_selB_HT1000.root'
#INPATH = 'stNTUPLES/DATA/JetHT_Run2016*_selA_HT1000.root'
#INPATH = 'stNTUPLES/DATA/JetHT_Run2016*_selA_HT200.root'
#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/MC/*_DoubleEMEnriched_selA_HT60.root'
INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/MC/GJet_DoubleEMEnriched_selA_HT60.root'
#OUTFOLDER = 'DATA/JetHT'
OUTFOLDER = 'MC/'

print " >> Plotting ST..."

fDir = []
for f in glob.glob(INPATH):
  fDir.append(os.path.splitext(f)[0])
print " >> Input file(s):",fDir

for stMin in stMins:
  for stMax in stMaxs:
    for nBin in nBins:
			print " .. >> ST range: [ %f -> %f ) GeV" % (stMin,stMax)
			print " .. >> nBins:",nBin
			subprocess.call("python plotST.py -i %s -l %f -r %f -b %d" % (fDir,stMin,stMax,nBin), shell=True)
			#os.rename("sT.png","PLOTS/%s/sT_selB_%d_%d_nBins%d.png"%(OUTFOLDER,stMin,stMax,nBin))
			os.rename("sT.png","PLOTS/%s/sT_selA_%d_%d_nBins%d.png"%(OUTFOLDER,stMin,stMax,nBin))
