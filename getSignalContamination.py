#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, tmROOTUtils, pdb

from tmProgressBar import tmProgressBar
from tmGeneralUtils import prettyPrintDictionary

inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination from input MC in real data.')
inputArgumentsParser.add_argument('--inputMCPath', required=True, help='Path to input MC file.',type=str)
inputArgumentsParser.add_argument('--inputDataPath', required=True, help='Path to input data file.',type=str)
inputArgumentsParser.add_argument('--crossSectionsFile', default="SusyCrossSections13TevGluGlu.txt", help='Path to dat file that contains cross-sections as a function of gluino mass, to use while weighting events.',type=str)
inputArgumentsParser.add_argument('--sTMin', default=1100., help='Min value of sT.',type=float)
inputArgumentsParser.add_argument('--outputPrefix', required=True, help='Prefix to output files.',type=str)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Minimum number of jets in event.',type=str)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Least value of nJets in highest nJets bin.',type=str)
inputArgumentsParser.add_argument('--analyze_nJetsBin', action='append', default=[], help='nJets to plot.',type=int)
inputArgumentsParser.add_argument('--nGluinoMassBins', default=20, help='nBins on the gluino mass axis.',type=int) # (800 - 25) GeV --> (1750 + 25) GeV in steps of 50 GeV
inputArgumentsParser.add_argument('--minGluinoMass', default=775., help='Min gluino mass.',type=float)
inputArgumentsParser.add_argument('--maxGluinoMass', default=1775., help='Max gluino mass.',type=float)
inputArgumentsParser.add_argument('--nNeutralinoMassBins', default=133, help='nBins on the neutralino mass axis.',type=int)
inputArgumentsParser.add_argument('--minNeutralinoMass', default=93.75, help='Min neutralino mass.',type=float)
inputArgumentsParser.add_argument('--maxNeutralinoMass', default=1756.25, help='Max neutralino mass.',type=float) # (100 - 6.25) GeV --> (1750 + 6.25) GeV in steps of 12.5 GeV
inputArgumentsParser.add_argument('--totalIntegratedLuminosity', default=37760., help='Max neutralino mass.',type=float) # total integrated luminosity in 2016 = 37.76 inverse femtobarns
inputArgumentsParser.add_argument('--nGeneratedEventsPerBin', default=150000, help='Number of generated events per bin in the MC samples.',type=int)
inputArguments = inputArgumentsParser.parse_args()

nJetsBinsToAnalyze=inputArguments.analyze_nJetsBin
for nJetsBinToAnalyze in nJetsBinsToAnalyze:
    if nJetsBinToAnalyze > inputArguments.nJetsMax: sys.exit("Error: index of jets bin to analyze can be at most nJetsMax={nJetsMax}; this argument is not allowed: --analyze_nJetsBin={bin}".format(nJetsMax=inputArguments.nJetsMax, bin=nJetsBinToAnalyze))

# pdb.set_trace()

MCPIDs = {
    "photon": 22,
    "gluino": 1000021,
    "neutralino": 1000022
}

crossSectionsInputFileObject = open(inputArguments.crossSectionsFile, 'r')
crossSectionsDictionary = {}
for line in crossSectionsInputFileObject:
    crossSectionsData = line.split()
    gluinoMass = int(0.5 + float(crossSectionsData[0]))
    crossSection = float(crossSectionsData[1])
    crossSectionsDictionary[gluinoMass] = crossSection
crossSectionsInputFileObject.close()

print("Read in cross-sections as a function of gluino mass:")
prettyPrintDictionary(crossSectionsDictionary)

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

def getWeight(gluinoMass, neutralinoMass):
    try:
        weight = crossSectionsDictionary[int(0.5 + gluinoMass)]*inputArguments.totalIntegratedLuminosity/inputArguments.nGeneratedEventsPerBin
    except KeyError:
        sys.exit("Unable to find cross-section corresponding to following gluino mass: {mass}".format(mass=gluinoMass))
    return weight

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
print ("Total number of available events in MC samples: {nMCEntries}".format(nMCEntries=nMCEntries))

histograms_total_nMCEvents = {}
histograms_weighted_nMCEvents = {}
histograms_signalContamination = {}
for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    histograms_total_nMCEvents[nJets] = ROOT.TH2F("h_total_nMCEvents_{nJets}Jets".format(nJets=nJets), "h_total_nMCEvents_{nJets}Jets".format(nJets=nJets), inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
    histograms_weighted_nMCEvents[nJets] = ROOT.TH2F("h_weighted_nMCEvents_{nJets}Jets".format(nJets=nJets), "h_weighted_nMCEvents_{nJets}Jets".format(nJets=nJets), inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
    histograms_signalContamination[nJets] = ROOT.TH2F("h_signalContamination_{nJets}Jets".format(nJets=nJets), "h_signalContamination_{nJets}Jets".format(nJets=nJets), inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)

progressBar = tmProgressBar(nMCEntries)
progressBarUpdatePeriod = max(1, nMCEntries//1000)
progressBar.initializeTimer()
for entryIndex in range(nMCEntries):
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
    if (sT > inputArguments.sTMin):
        generatedMasses = getGeneratedMasses(inputMCChain)
        generated_gluinoMass = generatedMasses["gluino"]
        generated_neutralinoMass = generatedMasses["neutralino"]
        eventWeight = getWeight(generatedMasses["gluino"], generatedMasses["neutralino"])
        histograms_total_nMCEvents[nJetsBin].Fill(generated_gluinoMass, generated_neutralinoMass, 1.)
        histograms_weighted_nMCEvents[nJetsBin].Fill(generated_gluinoMass, generated_neutralinoMass, eventWeight)
        histograms_signalContamination[nJetsBin].Fill(generated_gluinoMass, generated_neutralinoMass, eventWeight/nEventsInData[nJetsBin])
progressBar.terminate()

for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    if not(nJetsBin in nJetsBinsToAnalyze): continue
    tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [histograms_total_nMCEvents[nJetsBin]], canvasName = "c_total_nMCEvents_{nJetsBin}Jets".format(nJetsBin=nJetsBin), outputDocumentName = "analysis/signalContamination/{outputPrefix}_total_nEvents_{nJetsBin}Jets".format(outputPrefix=inputArguments.outputPrefix, nJetsBin=nJetsBin), customOptStat=0, customPlotOptions_firstObject="TEXTCOLZ")
    tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [histograms_weighted_nMCEvents[nJetsBin]], canvasName = "c_weighted_nMCEvents_{nJetsBin}Jets".format(nJetsBin=nJetsBin), outputDocumentName = "analysis/signalContamination/{outputPrefix}_weighted_nEvents_{nJetsBin}Jets".format(outputPrefix=inputArguments.outputPrefix, nJetsBin=nJetsBin), customOptStat=0, customPlotOptions_firstObject="TEXTCOLZ")
    tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [histograms_signalContamination[nJetsBin]], canvasName = "c_signalContamination_{nJetsBin}Jets".format(nJetsBin=nJetsBin), outputDocumentName = "analysis/signalContamination/{outputPrefix}_signalContamination_{nJetsBin}Jets".format(outputPrefix=inputArguments.outputPrefix, nJetsBin=nJetsBin), customOptStat=0, customPlotOptions_firstObject="COLZ")
