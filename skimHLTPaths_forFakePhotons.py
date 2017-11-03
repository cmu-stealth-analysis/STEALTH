#!/usr/bin/env python

from __future__ import print_function, division

import os
import sys
# import numpy as np
import argparse
import ROOT

from tmProgressBar import tmProgressBar

# Register command line options
parser = argparse.ArgumentParser(description='Photon skimmer')
parser.add_argument('-e','--era', required=True, help='Run Era',type=str)
parser.add_argument('-n', '--maxNEvents', default=0, help="maximum number of events", type=int)
args = parser.parse_args()

# Keep time
print(" >> Running Loose Photon Skim...")
sw = ROOT.TStopwatch()
sw.Start()

runEra = args.era
listOfInputFiles = []
inputFileNamesFileObject = open('inputFiles_2016{runEra}.txt'.format(runEra=args.era), 'r')
for inputFileName in inputFileNamesFileObject:
    listOfInputFiles.append(inputFileName.strip())
inputFileNamesFileObject.close()

# Load input TTrees into TChain
ggIn = ROOT.TChain("ggNtuplizer/EventTree")
for inputFile in listOfInputFiles:
    print("Adding... " + inputFile)
    ggIn.Add(inputFile)

nEvts = ggIn.GetEntries()
print(" >> nEvts:" + str(nEvts))

# Initialize output file as empty clone
outFileStr = "/eos/cms/store/user/tmudholk/stealth/ggSKIMS/DoubleEG_Run2016{runEra}_FebReMiniAOD_SKIM_DoubleFake.root".format(runEra=args.era)
outFile = ROOT.TFile(outFileStr, "RECREATE")
outDir = outFile.mkdir("ggNtuplizer")
outDir.cd()
ggOut = ggIn.CloneTree(0) 
print(" >> Output file: " + outFileStr)

##### SKIM START #####

# Event range to process
iEvtStart = 0
iEvtEnd   = nEvts
if (args.maxNEvents > 0): iEvtEnd = args.maxNEvents
#iEvtEnd   = 10000 
print(" >> Processing entries: [{start} -> {end})".format(start=iEvtStart, end=iEvtEnd))

nAcc = 0
progressBar = tmProgressBar(nEvts)
progressBar.initializeTimer()
for jEvt in range(iEvtStart,iEvtEnd):
    # Initialize event
    if jEvt > nEvts:
        break
    treeStatus = ggIn.LoadTree(jEvt)
    if treeStatus < 0:
        break
    evtStatus = ggIn.GetEntry(jEvt)
    if evtStatus <= 0:
        continue
    if (jEvt%1000 == 0):
        progressBar.updateBar(jEvt/(iEvtEnd-iEvtStart), jEvt)

    # Photon skim by trigger path
    if (ggIn.HLTPho>>14)&1 == 0: # HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90
        continue

    # Photon skim by number of photons
    #nPhotons = 0
    #for i in range(ggIn.nPho):
    #   if (ggIn.phoEt[i] > 20.0 and ggIn.phoIDbit[i]>>0&1 != 0): # >>0:loose, >>1:medium, >>2:tight
    #       nPhotons += 1
    #if nPhotons < nPhoCut:
    #   continue 

    # Write this evt to output tree
    ggOut.Fill()
    nAcc += 1

print("")
##### SKIM END #####
outFile.Write()
outFile.Close()

sw.Stop()
print(" >> nAccepted evts: {n}/{nTot} ({percent} %)".format(n=nAcc, nTot=(iEvtEnd-iEvtStart), percent=100*nAcc/(iEvtEnd-iEvtStart)))
print(" >> Real time: {realTime} minutes".format(realTime=sw.RealTime()/60))
print(" >> CPU time: {cpuTime} minutes".format(cpuTime=sw.CpuTime()/60))
