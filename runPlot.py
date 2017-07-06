#!/usr/bin/env python

import os
import sys 
import glob
import subprocess

#stMins = [500]
#stMins = [1200]
#stMaxs = [3500]
#nBins = [30]
#nBins = [23]
#stMins = [0]
#stMins = [900]

#stMins = [100]
#stMaxs = [3100]
#nBins = [20]
stMins = [1300]
stMaxs = [3700]
nBins = [12]
stMins = [1475]
stMaxs = [3675]
nBins = [11]
stMins = [100]
stMaxs = [3700]
nBins = [18]
#stMins = [1300]
#stMaxs = [4100]
#nBins = [14]

#nBins = [18]
#sel = 'A'
sel = 'B'
#HTcuts = [800,900,1000]
HTcuts = [1000]
#iNorms = [8]
#HTcuts = [60]
#iNorms = [7]
#iNorms = [13]
#iNorms = [8]
#iNorms = [9]
iNorms = [7]
#iNorms = [1]
DR = 0.4
DR = 0.

#eosDir = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews'
eosDir = '/eos/cms/store/user/mandrews'
#INPATH = '%s/MC/stNTUPLES/GJet_selA_HT60_Pho10_DeltaR00.root'%(eosDir)
#INPATH = '%s/DATA/stNTUPLES/DoubleEG_SepRereco_Run2016*_selA_HT60_DeltaR04.root'%(eosDir)
#INPATH = '%s/DATA/stNTUPLES/JetHT_SepRereco_Run2016*_selB_HT%d_DeltaR04.root'%(eosDir,HTcut)
#INPATH = '%s/DATA/stNTUPLES/JetHT_SepRereco_Run2016*_selB_HT800_DeltaR04.root'%(eosDir)
#OUTFOLDER = 'PLOTS/MC/SUSY_selA'
#OUTFOLDER = 'PLOTS/MC/SUSY_selB'
#OUTFOLDER = 'PLOTS/MC/GJet_selA'
#OUTFOLDER = 'PLOTS/MC/QCD_selB'
#OUTFOLDER = 'PLOTS/DATA/JetHT_SepRereco'
#OUTFOLDER = 'PLOTS/DATA/DoubleEG_SepRereco'
#OUTFOLDER = 'PLOTS/DATA/DoubleEG_ReminiAOD'
OUTFOLDER = 'PLOTS/DATA/JetHT_ReminiAOD'
#OUTFOLDER = 'PLOTS/MC/DiPhoton'

print " >> Plotting ST..."

#fDir = []
#for f in glob.glob(INPATH):
#	fDir.append(os.path.splitext(f)[0])
#print " >> Input file(s):",fDir

for stMin in stMins:
  for stMax in stMaxs:
    for nBin in nBins:
			print " .. >> ST range: [ %f -> %f ) GeV" % (stMin,stMax)
			print " .. >> nBins:",nBin
			for HTcut in HTcuts:
				#INPATH = '%s/DATA/stNTUPLES/JetHT_SepRereco_Run2016*_selB_HT%d_DeltaR%02d.root'%(eosDir,HTcut,DR*10.)
				INPATH = '%s/DATA/stNTUPLES/JetHT_ReminiAOD_Run2016*_selB_HT%d_DeltaR%02d.root'%(eosDir,HTcut,DR*10.)
				#INPATH = '%s/DATA/stNTUPLES/JetHT_ReminiAOD_Run2016*_selB_HT%d_DeltaR%02d.root'%(eosDir,HTcut,DR*10.)
				#INPATH = '%s/DATA/stNTUPLES/DoubleEG_ReminiAOD_Run2016*_selA_HT%d_DeltaR%02d.root'%(eosDir,HTcut,DR*10.)
				#INPATH = '%s/MC/stNTUPLES/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_selA_HT%d_DeltaR%02d.root'%(eosDir,HTcut,DR*10.)
				#INPATH = '%s/MC/stNTUPLES/DiPhoton*_selA_HT%d_DeltaR%02d.root'%(eosDir,HTcut,DR*10.)
				#INPATH = '%s/DATA/stNTUPLES/DoubleEG_SepRereco_Run2016*_selA_HT%d_DeltaR%02d.root'%(eosDir,HTcut,DR*10.)
				#INPATH = '%s/MC/stNTUPLES/GJet_selA_HT60_Pho10_DeltaR%02d.root'%(eosDir,DR*10.)
				#INPATH = '%s/MC/stNTUPLES/sms-t7WgStealthD_sel%s_HT%d_DeltaR%02d.root'%(eosDir,sel,HTcut,DR*10.)
				#INPATH = '%s/MC/stNTUPLES/MC_QCD_selB_HT%d_DeltaR%02d.root'%(eosDir,HTcut,DR*10.)
				fDir = []
				for f in glob.glob(INPATH):
					fDir.append(os.path.splitext(f)[0])
				print " >> Input file(s):",fDir
				for iNorm in iNorms:
					subprocess.call("python plotST.py -i %s -l %f -r %f -b %d -H %d -n %d" % (fDir,stMin,stMax,nBin,HTcut,iNorm), shell=True)
					#subprocess.call("python plotST_rebin.py -i %s -l %f -r %f -b %d -H %d -n %d" % (fDir,stMin,stMax,nBin,HTcut,iNorm), shell=True)
					os.rename("sT.png","%s/sT_sel%s_HT%d_iNorm_%d_%d_%d_DeltaR%02d_nBins%d.png"%(OUTFOLDER,sel,HTcut,iNorm,stMin,stMax,DR*10.,nBin))
