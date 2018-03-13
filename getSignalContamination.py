#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse

from tmProgressBar import tmProgressBar
from tmGeneralUtils import prettyPrintDictionary

inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination from input MC in real data.')
inputArgumentsParser.add_argument('--inputMCPath', required=True, help='Path to input MC file.',type=str)
inputArgumentsParser.add_argument('--inputDataPath', required=True, help='Path to input data file.',type=str)
inputArgumentsParser.add_argument('--sTMin', default=1100., help='Min value of sT.',type=float)
inputArgumentsParser.add_argument('--outputDirectoryName', required=True, help='Name of output directory.',type=str)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Minimum number of jets in event.',type=str)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Least value of nJets in highest nJets bin.',type=str)
inputArgumentsParser.add_argument('--nJets', action='append', help='nJets to plot.',type=str)
inputArguments = inputArgumentsParser.parse_args()

sw = ROOT.TStopwatch()
sw.Start()

print("Analyzing data sample...")
inputDataChain = ROOT.TChain('ggNtuplizer/EventTree')
inputDataChain.Add(inputArguments.inputDataPath)
nDataEntries = inputDataChain.GetEntries()
print ("Total number of available events in data: {nDataEntries}".format(nDataEntries=nDataEntries))

nEventsInData = {}
for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    nEventsInData[nJets] = 0

progressBar = tmProgressBar(nDataEntries)
progressBarUpdatePeriod = max(1, nDataEntries//1000)
progressBar.initializeTimer()
for entryIndex in range(nDataEntries):
    entryStatus = inputDataChain.LoadTree(entryIndex)
    if entryStatus < 0: sys.exit("Tree failed to load entry at index {entryIndex}; returned status {entryStatus}".format(entryIndex=entryIndex, entryStatus=entryStatus))
    chainStatus = inputDataChain.GetEntry(entryIndex)
    if chainStatus <= 0: sys.exit("Unable to load data from chain at index {entryIndex}; GetEntry returned status {chainStatus}".format(entryIndex=entryIndex, chainStatus=chainStatus))

    if (entryIndex%progressBarUpdatePeriod == 0): progressBar.updateBar(1.0*entryIndex/nDataEntries, entryIndex)

    nStealthJets = inputDataChain.b_nJets
    nJetsBin = nStealthJets
    if (nJetsBin > inputArguments.nJetsMax): nJetsBin = inputArguments.nJetsMax
    if (nJetsBin < inputArguments.nJetsMin): sys.exit("Unexpected nJetsBin = {nJetsBin} at entry index = {entryIndex}".format(nJetsBin=nJetsBin, entryIndex=entryIndex))
    sT = inputDataChain.b_evtST
    if (sT > inputArguments.sTMin): nEventsInData[nJetsBin] += 1
progressBar.terminate()

prettyPrintDictionary(inputDict=nEventsInData, keyPrintOrder=list(range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax)))
