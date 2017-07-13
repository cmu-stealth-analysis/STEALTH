#!/usr/bin/env python

# Script to run ST plotter over multiple scenarios

import os
import sys 
import glob
import subprocess

# ============================================
# Parameters to specify which files to process:

# Which selection
#sel = 'A'
sel = 'B'

# Which HT cuts 
#HTcuts = [800,900,1000]
HTcuts = [1000]
#HTcuts = [60]

# Which DeltaR to plot
DR = 0.4
#DR = 0.

# ============================================
# Processing parameters of plotter:

# Range over which to plot ST [GeV]
#stMins = [100]
#stMaxs = [3100]
stMins = [1300]
stMaxs = [3700]
#stMins = [1475]
#stMaxs = [3675]
#stMins = [100]
#stMaxs = [3700]
#stMins = [1300]
#stMaxs = [4100]

# Number of bins in ST
#nBins = [20]
nBins = [12]
#nBins = [11]
#nBins = [18]
#nBins = [14]
#nBins = [18]

# Which iBin index to use for ST normalization
# (0:underflow, 1:stMin, ...)
#iNorms = [8]
#iNorms = [7]
#iNorms = [13]
#iNorms = [8]
#iNorms = [9]
#iNorms = [7]
iNorms = [1]

# File output
eosDir = '/eos/cms/store/user/mandrews'
#OUTFOLDER = 'PLOTS/MC/SUSY_selA'
#OUTFOLDER = 'PLOTS/MC/SUSY_selB'
#OUTFOLDER = 'PLOTS/MC/GJet_selA'
#OUTFOLDER = 'PLOTS/MC/QCD_selB'
#OUTFOLDER = 'PLOTS/DATA/DoubleEG_ReminiAOD'
OUTFOLDER = 'PLOTS/DATA/JetHT_ReminiAOD'
#OUTFOLDER = 'PLOTS/MC/DiPhoton'
if not os.path.exists(OUTFOLDER):
    os.makedirs(OUTFOLDER)

# ============================================
# Plotting loop

print " >> Plotting ST..."

for stMin in stMins:
    for stMax in stMaxs:
        print " >> ST range: [ %f -> %f ) GeV" % (stMin,stMax)
        for nBin in nBins:
            print " >> nBins:",nBin
            for HTcut in HTcuts:
                INPATH = '%s/DATA/stNTUPLES/JetHT_ReminiAOD_Run2016*_selB_HT%d_DeltaR%02d.root'%(eosDir,HTcut,DR*10.)
                #INPATH = '%s/DATA/stNTUPLES/DoubleEG_ReminiAOD_Run2016*_selA_HT%d_DeltaR%02d.root'%(eosDir,HTcut,DR*10.)
                #INPATH = '%s/MC/stNTUPLES/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_selA_HT%d_DeltaR%02d.root'%(eosDir,HTcut,DR*10.)
                #INPATH = '%s/MC/stNTUPLES/DiPhoton*_selA_HT%d_DeltaR%02d.root'%(eosDir,HTcut,DR*10.)
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
                    #os.rename("sT.png","%s/sT_sel%s_HT%d_iNorm_%d_%d_%d_DeltaR%02d_nBins%d.png"%(OUTFOLDER,sel,HTcut,iNorm,stMin,stMax,DR*10.,nBin))
