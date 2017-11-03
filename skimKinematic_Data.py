#!/usr/bin/env python

from __future__ import print_function, division

import os
import sys
# import numpy as np
import argparse
import ROOT

from tmProgressBar import tmProgressBar

parameters = {
    "pTCutSubLeading": 25.,
    "pTCutLeading": 35.,
    "barrelEndcapTransitionEta": 1.479
}

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Create skim for double photons based on kinematic observables only.')
# inputArgumentsParser.add_argument('-e','--era', required=True, help='Run Era',type=str)
inputArgumentsParser.add_argument('--inputFromFile', action='store_true', help="Interpret inputFilePath as text file that has a list of input of files.")
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--outputFilePath', required=True, help='Prefix to output file name.',type=str)
inputArgumentsParser.add_argument('--counterStartInclusive', required=True, help="Event number from input file from which to start. The event with this index is included in the processing.", type=int)
inputArgumentsParser.add_argument('--counterEndInclusive', required=True, help="Event number from input file at which to end. The event with this index is included", type=int)
inputArguments = inputArgumentsParser.parse_args()

def passesKinematicSelection(inputTreeObject, photonIndex):
    absEta = abs(inputTreeObject.phoEta[photonIndex])
    if (absEta >= parameters["barrelEndcapTransitionEta"]): return False
    pT = inputTreeObject.phoEt[photonIndex]
    if (pT <= parameters["pTCutSubLeading"]): return False
    return True

# Keep time
print(" >> Running skim for two kinematic photons...")
sw = ROOT.TStopwatch()
sw.Start()

listOfInputFiles = []
if (inputArguments.inputFromFile):
    inputFileNamesFileObject = open(inputArguments.inputFilePath, 'r')
    for inputFileName in inputFileNamesFileObject:
        listOfInputFiles.append(inputFileName.strip())
    inputFileNamesFileObject.close()
else:
    listOfInputFiles.append(inputArguments.inputFilePath)

# Load input TTrees into TChain
ggIn = ROOT.TChain("ggNtuplizer/EventTree")

for inputFile in listOfInputFiles:
    print("Adding... " + inputFile)
    ggIn.Add(inputFile)

# Initialize output file as empty clone
outFile = ROOT.TFile(inputArguments.outputFilePath, "RECREATE")
outDir = outFile.mkdir("ggNtuplizer")
outDir.cd()
ggOut = ggIn.CloneTree(0)
print(" >> Output file: " + inputArguments.outputFilePath)

##### SKIM START #####

# Event range to process
iEvtStart = inputArguments.counterStartInclusive
iEvtEnd   = 1 + inputArguments.counterEndInclusive
print(" >> Processing entries: [{start} -> {end})".format(start=iEvtStart, end=iEvtEnd))

nAcc = 0
progressBar = tmProgressBar(iEvtEnd-iEvtStart)
progressBarUpdatePeriod = 1+((iEvtEnd-iEvtStart)//1000)
progressBar.initializeTimer()
for jEvt in range(iEvtStart,iEvtEnd):
    # Initialize event
    # if jEvt > nEvts:
    #     break
    treeStatus = ggIn.LoadTree(jEvt)
    if treeStatus < 0:
        break
    evtStatus = ggIn.GetEntry(jEvt)
    if evtStatus <= 0:
        continue
    if ((jEvt-iEvtStart) % progressBarUpdatePeriod == 0): progressBar.updateBar((jEvt-iEvtStart)/(iEvtEnd-iEvtStart), jEvt-iEvtStart)

    # Photon skim by trigger path
    # if (ggIn.HLTPho>>14)&1 == 0: # HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90
    #     continue

    # Kinematic photon skim
    nSubLeadingPhotons = 0
    nLeadingPhotons = 0
    for i in range(ggIn.nPho):
        if passesKinematicSelection(ggIn, i):
            nSubLeadingPhotons += 1
            if (ggIn.phoEt[i] > parameters["pTCutLeading"]): nLeadingPhotons += 1
    if not (nSubLeadingPhotons == 2 and nLeadingPhotons >= 1):
        continue

    # Write this evt to output tree
    ggOut.Fill()
    nAcc += 1

progressBar.terminate()
##### SKIM END #####
nEvtsOutput = ggOut.GetEntries()
if (nEvtsOutput == 0): print("WARNING: no output events detected!")
else: print("nSelectedEvents in tree: " + str(nEvtsOutput))
outFile.Write()
outFile.Close()

sw.Stop()
print(" >> nAccepted evts: {n}/{nTot} ({percent} %)".format(n=nAcc, nTot=(iEvtEnd-iEvtStart), percent=100*nAcc/(iEvtEnd-iEvtStart)))
print(" >> Real time: {realTime} minutes".format(realTime=sw.RealTime()/60))
print(" >> CPU time: {cpuTime} minutes".format(cpuTime=sw.CpuTime()/60))
