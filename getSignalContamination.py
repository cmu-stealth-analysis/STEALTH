#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, tmROOTUtils

from tmProgressBar import tmProgressBar
from tmGeneralUtils import prettyPrintDictionary

inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination from input MC in real data.')
inputArgumentsParser.add_argument('--inputMCPath', required=True, help='Path to input MC file.',type=str)
inputArgumentsParser.add_argument('--inputDataPath', required=True, help='Path to input data file.',type=str)
inputArgumentsParser.add_argument('--sTMin', default=1100., help='Min value of sT.',type=float)
inputArgumentsParser.add_argument('--outputPrefix', required=True, help='Prefix to output files.',type=str)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Minimum number of jets in event.',type=str)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Least value of nJets in highest nJets bin.',type=str)
inputArgumentsParser.add_argument('--nJets', action='append', help='nJets to plot.',type=str)
inputArgumentsParser.add_argument('--nGluinoMassBins', default=20, help='nBins on the gluino mass axis.',type=int) # 775 GeV --> 1775 GeV in steps of 50 GeV
inputArgumentsParser.add_argument('--minGluinoMass', default=775., help='Min gluino mass.',type=float)
inputArgumentsParser.add_argument('--maxGluinoMass', default=1775., help='Max gluino mass.',type=float)
inputArgumentsParser.add_argument('--nNeutralinoMassBins', default=133, help='nBins on the neutralino mass axis.',type=int)
inputArgumentsParser.add_argument('--minNeutralinoMass', default=93.75, help='Min neutralino mass.',type=float)
inputArgumentsParser.add_argument('--maxNeutralinoMass', default=1756.25, help='Max neutralino mass.',type=float)
inputArguments = inputArgumentsParser.parse_args()

MCPIDs = {
    "photon": 22,
    "gluino": 1000021,
    "neutralino": 1000022
}

def getGeneratedMasses(inputMCChain):
    generatedMasses = {"gluino": 0., "neutralino": 0.}
    gluinoMassIsSet = False
    neutralinoMassIsSet = False
    for genParticleIndex in range(inputMCChain.nMC):
        if (not(gluinoMassIsSet) and inputMCChain.mcPID[genParticleIndex] == MCPIDs["gluino"]):
            generatedMasses["gluino"] = inputMCChain.mcMass[genParticleIndex]
            gluinoMassIsSet = True
        if (not(neutralinoMassIsSet) and inputMCChain.mcMomPID[genParticleIndex] == MCPIDs["neutralino"]):
            generatedMasses["neutralino"] = inputMCChain.mcMomMass[genParticleIndex]
            neutralinoMassIsSet = True
        if (gluinoMassIsSet and neutralinoMassIsSet): break
    if not(gluinoMassIsSet): sys.exit("Unable to find gluino mass in an event!")
    if not(neutralinoMassIsSet): sys.exit("Unable to find neutralino mass in an event!")
    return generatedMasses

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

print("Number of events in data:")
prettyPrintDictionary(inputDict=nEventsInData, keyPrintOrder=list(range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax)))

print("Analyzing MC sample...")
inputMCChain = ROOT.TChain('ggNtuplizer/EventTree')
inputMCChain.Add(inputArguments.inputMCPath)
nMCEntries = inputMCChain.GetEntries()
print ("Total number of available events in MC: {nMCEntries}".format(nMCEntries=nMCEntries))

histograms_total_nMCEvents = {}
for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    histograms_total_nMCEvents[nJets] = ROOT.TH2F("h_total_nMCEvents_{nJets}Jets".format(nJets=nJets), "h_total_nMCEvents_{nJets}Jets".format(nJets=nJets), inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)

progressBar = tmProgressBar(nMCEntries)
progressBarUpdatePeriod = max(1, nMCEntries//1000)
progressBar.initializeTimer()
for entryIndex in range(10000):
    entryStatus = inputMCChain.LoadTree(entryIndex)
    if entryStatus < 0: sys.exit("Tree failed to load entry at index {entryIndex}; returned status {entryStatus}".format(entryIndex=entryIndex, entryStatus=entryStatus))
    chainStatus = inputMCChain.GetEntry(entryIndex)
    if chainStatus <= 0: sys.exit("Unable to load MC from chain at index {entryIndex}; GetEntry returned status {chainStatus}".format(entryIndex=entryIndex, chainStatus=chainStatus))

    if (entryIndex%progressBarUpdatePeriod == 0): progressBar.updateBar(1.0*entryIndex/nMCEntries, entryIndex)

    nStealthJets = inputMCChain.b_nJets
    nJetsBin = nStealthJets
    if (nJetsBin > inputArguments.nJetsMax): nJetsBin = inputArguments.nJetsMax
    if (nJetsBin < inputArguments.nJetsMin): sys.exit("Unexpected nJetsBin = {nJetsBin} at entry index = {entryIndex}".format(nJetsBin=nJetsBin, entryIndex=entryIndex))
    sT = inputMCChain.b_evtST
    generatedMasses = getGeneratedMasses(inputMCChain)
    if (sT > inputArguments.sTMin):
        histograms_total_nMCEvents[nJetsBin].Fill(generatedMasses["gluino"], generatedMasses["neutralino"], 1./1000)
progressBar.terminate()

for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [histograms_total_nMCEvents[nJets]], canvasName = "c_total_nMCEvents_{nJets}Jets".format(nJets=nJets), outputDocumentName = "analysis/signalContamination/{outputPrefix}_{nJets}Jets".format(outputPrefix=inputArguments.outputPrefix, nJets=nJets), customOptStat=0, customPlotOptions_firstObject="TEXTCOLZ")
