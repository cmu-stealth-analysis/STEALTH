#!/usr/bin/env python

import os
import sys 
import glob
import subprocess

#stMins = [750,1000,1250,1500]
#stMaxs = [3250,3500,3750,4000]
stMins = [1250]
#stMins = [1000]
#stMins = [500]
#stMins = [160]
#stMaxs = [3500]
stMaxs = [3750]
#stMaxs = [2960]
#nBins = [4,5,6]
#nBins = [20]
nBins = [5,20]
#nBins = [5]
#sel = 'A'
sel = 'B'
#sel = 'C'

#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/DATA/JetHT/JetHT_2016*_sT_Pho30MedEleVeto_Jet30MedHt1000_El15Tight_Mu15Tight.root'
#INPATH = 'stNTUPLES/DATA/JetHT_Run2016*_selB_HT1000.root'
#INPATH = 'stNTUPLES/DATA/JetHT_Run2016*_selA_HT1000.root'
#INPATH = 'stNTUPLES/DATA/JetHT_Run2016*_selA_HT60.root'
INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/stNTUPLES/DATA/JetHT_SepRereco_Run2016*_selB_HT1000.root'
#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/stNTUPLES/DATA/tmp/JetHT_SepRereco_Run2016*_selC_HT1000.root'
#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/stNTUPLES/MC/QCD_selB_HT1000.root'
#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/stNTUPLES/MC/GJet_selA_HT60.root'
#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/MC/GJet_DoubleEMEnriched_selA_HT60.root'
#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/MC/GJet_40toInf_DoubleEMEnriched_selA_HT60.root'
#OUTFOLDER = 'PLOTS/MC/GJet_selA'
#OUTFOLDER = 'PLOTS/MC/QCD_selB'
OUTFOLDER = 'PLOTS/DATA/JetHT_SepRereco'

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
			os.rename("sT.png","%s/sT_sel%s_%d_%d_nBins%d.png"%(OUTFOLDER,sel,stMin,stMax,nBin))
