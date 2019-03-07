#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, ROOT, argparse, array, math
import numpy as np
import ROOT, tmROOTUtils, tmStatsUtils, tmGeneralUtils, tmProgressBar

inputArgumentsParser = argparse.ArgumentParser(description='Get data event histograms and systematics.')
inputArgumentsParser.add_argument('--inputFilesList_background', required=True, help='Semicolon-separated list of paths to n-tuples containing background events, possibly containing wildcards, to pass to TChain.Add.',type=str)
inputArgumentsParser.add_argument('--inputFilesList_signal', required=True, help='Semicolon-separated list of paths to n-tuples containing signal events, possibly containing wildcards, to pass to TChain.Add.',type=str)
inputArgumentsParser.add_argument('--inputFile_MCEventHistograms', default="analysis/MCEventHistograms/MC_2018_savedObjects.root", help='Path to root file that contains weighted number of signal events per (gluino mass, neutralino mass) bin.', type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="analysis/contaminatedSamples/", help='Directory in which to store data systematics.',type=str)
inputArgumentsParser.add_argument('--outputPrefix', required=True, help='Prefix to all output file names.',type=str)
inputArgumentsParser.add_argument('--contaminationStrength', default=0., help="Contamination strength to use in making new samples.", type=float)
inputArgumentsParser.add_argument('--contaminationGluinoMassBin', default=-1, help='Gluino mass bin to use for subtracting potential signal.',type=int)
inputArgumentsParser.add_argument('--contaminationNeutralinoMassBin', default=-1, help='Neutralino mass bin to use for subtracting potential signal.',type=int)
inputArgumentsParser.add_argument('--inputFile_MCTemplate', default="plot_susyMasses_template.root", help='Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.', type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Min number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Max number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsNorm', default=2, help='Number of jets w.r.t. which to normalize the sT distributions for other jets.',type=int)
inputArguments = inputArgumentsParser.parse_args()

if (inputArguments.contaminationStrength > 0.):
    if (inputArguments.inputFilesList_signal == ""): sys.exit("ERROR: contamination strength = {s} is positive but the input signal files list is empty.".format(s=inputArguments.contaminationStrength))
    if ((inputArguments.contaminationGluinoMassBin < 0) or ((inputArguments.contaminationNeutralinoMassBin < 0))): sys.exit("ERROR: contamination strength = {s} is positive but the gluino mass bin = {gMB} or neutralino mass bin = {nMB} is negative.".format(s=inputArguments.contaminationStrength, gMB=inputArguments.contaminationGluinoMassBin, nMB=inputArguments.contaminationNeutralinoMassBin))

MCPIDs = {
    "photon": 22,
    "gluino": 1000021,
    "neutralino": 1000022
}

commonOutputString = ("{p}_contaminated_gluinoMassBin{gMB}_neutralinoMassBin{nMB}_signalStrength{s}".format(p=inputArguments.outputPrefix, gMB=inputArguments.contaminationGluinoMassBin, nMB=inputArguments.contaminationNeutralinoMassBin, s=inputArguments.contaminationStrength)).replace(".", "pt")

STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries)
STBoundaries = []
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()):
        STBoundary = float(STBoundaryString.strip())
        STBoundaries.append(STBoundary)
STBoundaries.append(14000.0) # Instead of infinity
nSTSignalBins = len(STBoundaries) - 2 # First two lines are lower and upper boundaries for the normalization bin
STNormRangeMin = STBoundaries[0]
STNormRangeMax = STBoundaries[1]
print("Using {n} signal bins for ST; norm range min: {mn}, norm range max: {mx}.".format(n = nSTSignalBins, mn=STNormRangeMin, mx=STNormRangeMax))
STRegionsAxis = ROOT.TAxis(len(STBoundaries)-1, array.array('d', STBoundaries))

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

def signalAddedFully(test_enoughEventsAdded):
    for STRegionIndex in range(1, nSTSignalBins+2):
        for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
            if not(test_enoughEventsAdded[STRegionIndex][nJetsBin]): return False
    return True

# Fetch target number of events
nSignalEventsAdded = {}
enoughEventsAdded = {}
nSignalEventsTarget = {}
signalSource = ROOT.TFile("{i}".format(i=inputArguments.inputFile_MCEventHistograms))
for STRegionIndex in range(1, nSTSignalBins+2):
    nSignalEventsAdded[STRegionIndex] = {}
    enoughEventsAdded[STRegionIndex] = {}
    nSignalEventsTarget[STRegionIndex] = {}
    for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
        nSignalEventsAdded[STRegionIndex][nJetsBin] = 0
        enoughEventsAdded[STRegionIndex][nJetsBin] = False
        sourceHist = ROOT.TH2F()
        signalSource.GetObject("h_lumiBasedYearWeightedNEvents_{n}Jets_STRegion{i}".format(n=nJetsBin, i=STRegionIndex), sourceHist)
        nSignalEventsTarget[STRegionIndex][nJetsBin] = int(0.5+sourceHist.GetBinContent(sourceHist.GetBin(inputArguments.contaminationGluinoMassBin, inputArguments.contaminationNeutralinoMassBin)))
        if (nSignalEventsTarget[STRegionIndex][nJetsBin] == 0): enoughEventsAdded[STRegionIndex][nJetsBin] = True

print("Targeting:")
tmGeneralUtils.prettyPrintDictionary(nSignalEventsTarget)

inputChain = ROOT.TChain('ggNtuplizer/EventTree')
for patternToAdd in inputArguments.inputFilesList_background.strip().split(";"):
    print("Adding pattern: {p}".format(p=patternToAdd))
    inputChain.Add(patternToAdd)
nEntries = inputChain.GetEntries()
print ("Total number of available background events: {nEntries}".format(nEntries=nEntries))

outputFile = ROOT.TFile('{oD}/{cOS}.root'.format(oD=inputArguments.outputDirectory, cOS=commonOutputString), 'RECREATE')
outputDirectory = outputFile.mkdir("ggNtuplizer")
outputDirectory.cd()
outputTree = inputChain.CloneTree(0)

progressBar = tmProgressBar.tmProgressBar(nEntries)
progressBarUpdatePeriod = max(1, nEntries//1000)
progressBar.initializeTimer()
for entryIndex in range(nEntries):
    entryStatus = inputChain.LoadTree(entryIndex)
    if entryStatus < 0: sys.exit("Tree failed to load entry at index {entryIndex}; returned status {entryStatus}".format(entryIndex=entryIndex, entryStatus=entryStatus))
    chainStatus = inputChain.GetEntry(entryIndex)
    if chainStatus <= 0: sys.exit("Unable to load data from chain at index {entryIndex}; GetEntry returned status {chainStatus}".format(entryIndex=entryIndex, chainStatus=chainStatus))
    if (entryIndex%progressBarUpdatePeriod == 0): progressBar.updateBar(1.0*entryIndex/nEntries, entryIndex)
    outputTree.Fill()
progressBar.terminate()

for patternToAdd in inputArguments.inputSignalFilesList.strip().split(";"):
    print("Adding pattern: {p}".format(p=patternToAdd))
    inputChain.Add(patternToAdd)
nTotalEntries = inputChain.GetEntries()
print ("New total number of available signal events: {nSignalEntries}".format(nSignalEntries=nSignalEntries))

generatedMCTemplate = ROOT.TFile(inputArguments.inputFile_MCTemplate, "READ")
h_MCTemplate = ROOT.TH2F()
generatedMCTemplate.GetObject("h_susyMasses_template", h_MCTemplate)
if (h_MCTemplate):
    print("Found and opened MC template.")
else:
    sys.exit("Unable to find object \"h_susyMasses_template\" in file {f}".format(f=inputArguments.MCTemplate))

nEvents_neutralinoMassUnset = 0
addedFully = False
progressBar = tmProgressBar.tmProgressBar(nSignalEntries)
progressBarUpdatePeriod = max(1, nSignalEntries//1000)
progressBar.initializeTimer()
for entryIndex in range(nEntries, nTotalEntries):
    entryStatus = inputChain.LoadTree(entryIndex)
    if entryStatus < 0: sys.exit("Tree failed to load entry at index {entryIndex}; returned status {entryStatus}".format(entryIndex=entryIndex, entryStatus=entryStatus))
    chainStatus = inputChain.GetEntry(entryIndex)
    if chainStatus <= 0: sys.exit("Unable to load data from chain at index {entryIndex}; GetEntry returned status {chainStatus}".format(entryIndex=entryIndex, chainStatus=chainStatus))

    if (entryIndex%progressBarUpdatePeriod == 0): progressBar.updateBar(1.0*entryIndex/nSignalEntries, entryIndex)

    generatedMasses = {"gluino": 0., "neutralino": 0.}
    gluinoMassIsSet = False
    gluinoMassBin = -1
    gluinoMass = -1.
    neutralinoMassIsSet = False
    neutralinoMassBin = -1
    neutralinoMass = -1.
    for genParticleIndex in range(inputChain.nMC):
        if (not(gluinoMassIsSet) and inputChain.mcPID[genParticleIndex] == MCPIDs["gluino"]):
            generatedMasses["gluino"] = inputChain.mcMass[genParticleIndex]
            gluinoMassIsSet = True
        if (not(neutralinoMassIsSet) and inputChain.mcMomPID[genParticleIndex] == MCPIDs["neutralino"]):
            generatedMasses["neutralino"] = inputChain.mcMomMass[genParticleIndex]
            neutralinoMassIsSet = True
        if (gluinoMassIsSet and neutralinoMassIsSet): break
    if gluinoMassIsSet:
        gluinoMassBin = h_MCTemplate.GetXaxis().FindFixBin(generatedMasses["gluino"])
        gluinoMass = int(0.5+h_MCTemplate.GetXaxis().GetBinCenter(gluinoMassBin))
    else:
        sys.exit("Gluino mass unset in event with index {index}".format(index=entryIndex))
    if neutralinoMassIsSet:
        neutralinoMassBin = h_MCTemplate.GetYaxis().FindFixBin(generatedMasses["neutralino"])
        neutralinoMass = int(0.5+h_MCTemplate.GetYaxis().GetBinCenter(neutralinoMassBin))
    else:
        nEvents_neutralinoMassUnset += 1
        continue

    if not((gluinoMassBin == inputArguments.contaminationGluinoMassBin) and (neutralinoMassBin == inputArguments.contaminationNeutralinoMassBin)): continue

    nStealthJets = inputChain.b_nJets
    nJetsBin = nStealthJets
    if (nJetsBin > inputArguments.nJetsMax): nJetsBin = inputArguments.nJetsMax
    if (nJetsBin < inputArguments.nJetsMin):
        print("Unexpected nJetsBin = {nJetsBin} at entry index = {entryIndex}".format(nJetsBin=nJetsBin, entryIndex=entryIndex))
        continue

    sT = inputChain.b_evtST
    if ((sT <= sTKernelEstimatorRangeMin) or (sT >= sTKernelEstimatorRangeMax)): continue
    STRegionIndex = STRegionsAxis.FindFixBin(sT)

    if not(enoughEventsAdded[STRegionIndex][nJetsBin]):
        outputTree.Fill()
        nSignalEventsAdded[STRegionIndex][nJetsBin] += 1
        if (nSignalEventsAdded[STRegionIndex][nJetsBin] == nSignalEventsTarget[STRegionIndex][nJetsBin]): enoughEventsAdded[STRegionIndex][nJetsBin] = True

    if signalAddedFully(enoughEventsAdded):
        addedFully = True
        break

progressBar.terminate()
print("nEvents with neutralino mass unset: {n}".format(n=nEvents_neutralinoMassUnset))
if not(addedFully): sys.exit("ERROR: not enough signal events available for given contamination strength.")

outputFile.Write()
outputFile.Close()

print("All done!")
