#!/usr/bin/env python

from __future__ import print_function, division

import os
import sys
import glob
import ROOT
import re
import time
import glob
import argparse

from tmProgressBar import tmProgressBar

inputArgumentsParser = argparse.ArgumentParser(description='Merge several files into a single one.')
inputArgumentsParser.add_argument('--escapedInputFilePattern', required=True, help='Escaped glob pattern to select input files.',type=str)
inputArgumentsParser.add_argument('--outputFilePath', required=True, help='Path to output file.',type=str)
inputArguments = inputArgumentsParser.parse_args()

# Keep time
sw = ROOT.TStopwatch()
sw.Start()

ggIn = ROOT.TChain("ggNtuplizer/EventTree")
ggIn.SetMaxTreeSize(100000000000) # 1 TB

listOfInputFiles = glob.glob(inputArguments.escapedInputFilePattern)

for inputFile in listOfInputFiles:
    print ("Adding... " + inputFile)
    ggIn.Add(inputFile)
nEvts = ggIn.GetEntries()
print(" >> total nEvts:" + str(nEvts))
if (nEvts == 0): sys.exit("No events found!")

outFile = ROOT.TFile(inputArguments.outputFilePath, "RECREATE")
outDir = outFile.mkdir("ggNtuplizer")
outDir.cd()
ggOut = ggIn.CloneTree(0)
print(" >> Output file: " + inputArguments.outputFilePath)

progressBar = tmProgressBar(nEvts)
progressBarUpdatePeriod = 1+(nEvts//1000)
progressBar.initializeTimer()
for jEvt in range(0, nEvts):
    treeStatus = ggIn.LoadTree(jEvt)
    if treeStatus < 0:
        break
    evtStatus = ggIn.GetEntry(jEvt)
    if evtStatus <= 0:
        continue
    if (jEvt % progressBarUpdatePeriod == 0): progressBar.updateBar(jEvt/nEvts, jEvt)
    ggOut.Fill()

progressBar.terminate()
print(" >> Merge done. Checking...")
nEvtsOutput = ggOut.GetEntries()
if (nEvtsOutput != nEvts): sys.exit("ERROR: nEvts in output file not equal to total available nEvts!")
print(" >> Writing to output file...")
outFile.Write()
outFile.Close()
print ("nEvtsOutput: " + str(nEvtsOutput))

sw.Stop()
print(" >> Real time: {realTime} minutes".format(realTime=sw.RealTime()/60))
print(" >> CPU time: {cpuTime} minutes".format(cpuTime=sw.CpuTime()/60))
