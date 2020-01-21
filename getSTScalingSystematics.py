#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, array, pdb, math
import ROOT, tmROOTUtils, tmGeneralUtils
from tmProgressBar import tmProgressBar

ROOT.gROOT.SetBatch(ROOT.kTRUE)

inputArgumentsParser = argparse.ArgumentParser(description='Get data event histograms and systematics.')
# inputArgumentsParser.add_argument('--inputFile_nEvtsObserved', required=True, help='Input text file containing the observed number of events in the control region in each (ST, nJets) bin.',type=str)
inputArgumentsParser.add_argument('--inputFile_nEvtsExpected', required=True, help='Input text file containing the observed number of events in the control region in each (ST, nJets) bin.',type=str)
inputArgumentsParser.add_argument('--inputFile_shapeSystematics', required=True, help='Input text file containing the shape systematics in each ST bin.',type=str)
inputArgumentsParser.add_argument('--inputFile_nSignalEvents', required=True, help='Input ROOT file containing a 2D histogram of the number of events in the control region in each (ST, nJets) bin in the (m_progenitor, m_neutralino) phase space.',type=str)
inputArgumentsParser.add_argument('--eventProgenitorMassBin', required=True, help='Event progenitor mass bin.', type=int)
inputArgumentsParser.add_argument('--neutralinoMassBin', required=True, help='Neutralino mass bin.', type=int)
inputArgumentsParser.add_argument('--bestFitSignalStrength', required=True, help='Signal strength to subtract.',type=float)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Min number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Max number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsNorm', default=2, help='Number of jets w.r.t. which to normalize the sT distributions for other jets.',type=int)
inputArguments = inputArgumentsParser.parse_args()

def get_systematics_dict(sourceNEvents, targetNEvents, nSTSignalBins):
    if ((len(sourceNEvents) == 0) or (len(targetNEvents) == 0)): sys.exit("ERROR: One of sourceNEvents and targetNEvents is empty. sourceNEvents={sNE}, targetNEvents={tNE}".format(sNE=sourceNEvents, tNE=targetNEvents))
    systematicsDictionary = {}
    for STRegionIndex in range(1, nSTSignalBins+2):
        sourceRatio = sourceNEvents[STRegionIndex]/sourceNEvents[1] # 1 = normalization bin
        targetRatio = targetNEvents[STRegionIndex]/targetNEvents[1] # 1 = normalization bin
        systematicsDictionary[STRegionIndex] = (sourceRatio/targetRatio)-1.0
    return systematicsDictionary

STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries)
STBoundaries = []
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()):
        STBoundary = float(STBoundaryString.strip())
        STBoundaries.append(STBoundary)
STBoundaries.append(20000.0) # Instead of infinity
nSTSignalBins = len(STBoundaries) - 2 # First two lines are lower and upper boundaries for the normalization bin
print("Using {n} signal bins in ST.".format(n=nSTSignalBins))

signalBinLabels = []
for STRegionIndex in range(1, 2 + nSTSignalBins):
    for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
        signalBinLabel = "STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin)
        signalBinLabels.append(signalBinLabel)

eventProgenitorBinIndex = inputArguments.eventProgenitorMassBin
neutralinoBinIndex = inputArguments.neutralinoMassBin

dataSTScalingSystematicsList = []
MCEventHistogramsFile = ROOT.TFile.Open(inputArguments.inputFile_nSignalEvents, "READ")
# fileContents_nEventsObserved = tmGeneralUtils.getConfigurationFromFile(inputFilePath=inputArguments.inputFile_nEvtsObserved)
fileContents_nEventsExpected = tmGeneralUtils.getConfigurationFromFile(inputFilePath=inputArguments.inputFile_nEvtsExpected)
fileContents_shapeSystematics = tmGeneralUtils.getConfigurationFromFile(inputFilePath=inputArguments.inputFile_shapeSystematics)
fractionalUncertainties_shape = {}
expectedNEvents = {}
expectedNEvents_realKernel = {}
nEventsBestFitSignalStrength = {}
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    fractionalUncertainties_shape[nJetsBin] = {}
    expectedNEvents[nJetsBin] = {}
    expectedNEvents_realKernel[nJetsBin] = {}
    nEventsBestFitSignalStrength[nJetsBin] = {}
    for STRegionIndex in range(1, 2 + nSTSignalBins):
        # Get nEvents observed and expected from the background
        signalBinLabel = "STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin)
        # observedNEvents = fileContents_nEventsObserved["observedNEvents_{l}".format(l=signalBinLabel)]
        expectedNEvents[nJetsBin][STRegionIndex] = fileContents_nEventsExpected["expectedNEvents_{l}".format(l=signalBinLabel)]

        if (nJetsBin == inputArguments.nJetsNorm):
            expectedNEvents_realKernel[nJetsBin][STRegionIndex] = fileContents_nEventsExpected["expectedNEvents_{l}".format(l=signalBinLabel)]
        else:
            expectedNEvents_realKernel[nJetsBin][STRegionIndex] = fileContents_nEventsExpected["expectedNEvents_realKernel_{l}".format(l=signalBinLabel)]

        # Get shape uncertainties
        fractionalUncertainties_shape[nJetsBin][STRegionIndex] = fileContents_shapeSystematics["fractionalUncertainty_shape_{sBL}".format(sBL=signalBinLabel)]

        # Get nEvents expected from the signal, at signal strength 1
        MCHistograms_weightedNEvents = ROOT.TH2F()
        MCEventHistogramsFile.GetObject("h_lumiBasedYearWeightedNEvents_{l}".format(l=signalBinLabel), MCHistograms_weightedNEvents)
        if (not(MCHistograms_weightedNEvents)):
            sys.exit("ERROR: Histogram MCHistograms_weightedNEvents appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
        nEventsExpectedSignal_unscaled = MCHistograms_weightedNEvents.GetBinContent(eventProgenitorBinIndex, neutralinoBinIndex)
        # Scale number of events by the best fit signal strength    
        nEventsBestFitSignalStrength[nJetsBin][STRegionIndex] = (inputArguments.bestFitSignalStrength)*nEventsExpectedSignal_unscaled

realBackground = {}
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    realBackground[nJetsBin] = {}
    for STRegionIndex in range(1, 2 + nSTSignalBins):
        # Real background = expected nEvents - best-fit signal
        realBackground[nJetsBin][STRegionIndex] = max(0., expectedNEvents_realKernel[nJetsBin][STRegionIndex] - nEventsBestFitSignalStrength[nJetsBin][STRegionIndex])

for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    source_NEvents = {}
    target_NEvents = {}
    for STRegionIndex in range(1, 2 + nSTSignalBins):
        source_NEvents[STRegionIndex] = (realBackground[nJetsBin][1])*(realBackground[inputArguments.nJetsNorm][STRegionIndex]/realBackground[inputArguments.nJetsNorm][1])
        target_NEvents[STRegionIndex] = realBackground[nJetsBin][STRegionIndex]

    fractionalUncertainties_raw_scaling = get_systematics_dict(sourceNEvents=source_NEvents, targetNEvents=target_NEvents, nSTSignalBins=nSTSignalBins)
    for STRegionIndex in range(1, 2 + nSTSignalBins):
        fractionalUncertainty_residual = max(0., abs(fractionalUncertainties_raw_scaling[STRegionIndex]) - abs(fractionalUncertainties_shape[nJetsBin][STRegionIndex]))
        signalBinLabel = "STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin)
        dataSTScalingSystematicsList.append(tuple(["float", "fractionalUncertainty_residual_scaling_{l}".format(l=signalBinLabel), fractionalUncertainty_residual]))

MCEventHistogramsFile.Close()
tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=dataSTScalingSystematicsList, outputFilePath="dataSystematics_scaling.dat")
