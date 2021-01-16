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
inputArgumentsParser.add_argument('--inputEscapedPattern', action='store_true', help="Interpret inputFilePath a glob pattern. WARNING: glob.glob returns list of files matching wildcard expansion IN ARBITRARY ORDER, do NOT mess around with order of events!")
inputArgumentsParser.add_argument('--inputFromFile', action='store_true', help="Interpret inputFilePath as text file that has a list of input of files.")
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--outputFilePath', required=True, help='Path to output file.',type=str)
inputArgumentsParser.add_argument('--maxNEvents', default=-1, help='Only operate over the first N events; useful for creating subsets of data.',type=int)
inputArguments = inputArgumentsParser.parse_args()

# Keep time
sw = ROOT.TStopwatch()
sw.Start()

ggIn = ROOT.TChain("ggNtuplizer/EventTree")
ggIn.SetMaxTreeSize(100000000000) # 1 TB

listOfInputFiles=[]
if(inputArguments.inputEscapedPattern):
    listOfInputFiles = glob.glob(inputArguments.inputFilePath) # WARNING!!! WARNING!!! glob.glob returns list of files matching wildcard expansion IN ARBITRARY ORDER, do NOT mess around with order of events!
elif(inputArguments.inputFromFile):
    inputFileNamesFileObject = open(inputArguments.inputFilePath, 'r')
    for inputFileName in inputFileNamesFileObject:
        listOfInputFiles.append(inputFileName.strip())
    inputFileNamesFileObject.close()
else:
    listOfInputFiles = [inputArguments.inputFilePath]

for inputFile in listOfInputFiles:
    print ("Adding... " + inputFile)
    ggIn.Add(inputFile)
nEvts = ggIn.GetEntries()
print(" >> total nEvts:" + str(nEvts))
if (nEvts == 0): sys.exit("No events found!")

outFile = ROOT.TFile.Open(inputArguments.outputFilePath, "RECREATE")
outDir = outFile.mkdir("ggNtuplizer")
outDir.cd()
ggOut = ggIn.CloneTree(0)
print(" >> Output file: " + inputArguments.outputFilePath)

nEvtsToIterateOver=nEvts
if ((inputArguments.maxNEvents > 0) and (inputArguments.maxNEvents < nEvts)): nEvtsToIterateOver = inputArguments.maxNEvents
progressBar = tmProgressBar(nEvtsToIterateOver)
progressBarUpdatePeriod = 1+(nEvtsToIterateOver//50)
progressBar.initializeTimer()
for jEvt in range(0, nEvtsToIterateOver):  # WARNING!!! WARNING!!! glob.glob returns list of files matching wildcard expansion IN ARBITRARY ORDER, do NOT mess around with order of events!
    treeStatus = ggIn.LoadTree(jEvt)
    if treeStatus < 0:
        break
    evtStatus = ggIn.GetEntry(jEvt)
    if evtStatus <= 0:
        continue
    if (jEvt % progressBarUpdatePeriod == 0): progressBar.updateBar(jEvt/nEvtsToIterateOver, jEvt)
    ggOut.Fill()

progressBar.terminate()
print(" >> Merge done. Checking...")
nEvtsOutput = ggOut.GetEntries()
if (nEvtsOutput != nEvtsToIterateOver): sys.exit("ERROR: nEvts in output file not equal to target!")
print(" >> Writing to output file...")
outFile.Write()
outFile.Close()
print ("nEvtsOutput: " + str(nEvtsOutput))

sw.Stop()
print(" >> Real time: {realTime} minutes".format(realTime=sw.RealTime()/60))
print(" >> CPU time: {cpuTime} minutes".format(cpuTime=sw.CpuTime()/60))
