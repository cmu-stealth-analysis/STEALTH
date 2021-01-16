#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, ROOT, tmROOTUtils, pdb
from tmProgressBar import tmProgressBar

ROOT.gROOT.SetBatch(ROOT.kTRUE)

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Run STEALTH selection.')
inputArgumentsParser.add_argument('--inputFromFile', action='store_true', help="Interpret inputFilePath as text file that has a list of input of files.")
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--outputFilePrefix', required=True, help='Prefix for output file.',type=str)
inputArgumentsParser.add_argument('--maxNEvents', default=-1, help='Maximum number of events to process.',type=int)
inputArgumentsParser.add_argument('--bitMaxValue', default=38, help='Maximum value of photon bit up to which to plot.', type=int)
inputArguments = inputArgumentsParser.parse_args()

def passesHLTBit(inputChainObject=None, bitIndex=None):
    if (inputChainObject is None): sys.exit("inputChainObject is None!")
    if (bitIndex < 0): sys.exit("bitIndex = {bI} is ill-defined!".format(bI = bitIndex))
    return not(((inputChainObject.HLTPho>>bitIndex)&1) == 0)

listOfInputFiles = []
if (inputArguments.inputFromFile):
    inputFileNamesFileObject = open(inputArguments.inputFilePath, 'r')
    for inputFileName in inputFileNamesFileObject:
        listOfInputFiles.append(inputFileName.strip())
    inputFileNamesFileObject.close()
else:
    listOfInputFiles.append(inputArguments.inputFilePath)
# Load input TTrees into TChain
inputTreeObject = ROOT.TChain("ggNtuplizer/EventTree")
for inputFile in listOfInputFiles:
    print("Adding... " + inputFile)
    inputTreeObject.Add(inputFile)

nEvents = inputTreeObject.GetEntries()
if (nEvents == 0): sys.exit("Number of available events is 0.")

nEventsToProcess = 0
if (inputArguments.maxNEvents > 0 and inputArguments.maxNEvents < nEvents): nEventsToProcess = inputArguments.maxNEvents
else: nEventsToProcess = nEvents

h_HLTBitsHistogram = ROOT.TH1F("h_HLTBitsHistogram", "HLT Bits;Bit Index;A.U.", 1 + inputArguments.bitMaxValue, -0.5, 0.5+inputArguments.bitMaxValue)

progressBar = tmProgressBar(nEventsToProcess)
progressBarUpdatePeriod = 1+(nEventsToProcess//50)
progressBar.initializeTimer()
# pdb.set_trace()

for eventIndex in range(0,nEventsToProcess):
    treeStatus = inputTreeObject.LoadTree(eventIndex)
    if treeStatus < 0:
        # sys.exit("Tree unreadable.")
        break
    evtStatus = inputTreeObject.GetEntry(eventIndex)
    if evtStatus <= 0:
        # sys.exit("Event in tree unreadable.")
        continue
    if (eventIndex%progressBarUpdatePeriod == 0 or eventIndex == (nEventsToProcess - 1)): progressBar.updateBar(eventIndex/nEventsToProcess, eventIndex)
    for bitIndex in range(0, 1+inputArguments.bitMaxValue):
        if (passesHLTBit(inputTreeObject, bitIndex)): h_HLTBitsHistogram.Fill(bitIndex)

progressBar.terminate()

tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[h_HLTBitsHistogram], canvasName="c_HLTBitsHistogram", outputDocumentName="temp/{prefix}_HLTPhotonBits".format(prefix=inputArguments.outputFilePrefix), customOptStat=0, enableLogY = True)
