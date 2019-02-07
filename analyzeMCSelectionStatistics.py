#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, ROOT, tmROOTUtils, array, pdb, math
from tmProgressBar import tmProgressBar

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Merge and analyze MC selection statistics.')
inputArgumentsParser.add_argument('--inputFilesList', required=True, help='Path to file containing newline-separated list of input files.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="MCSelectionStatistics", help='Output directory.',type=str)
inputArgumentsParser.add_argument("--nGluinoMassBins", default=20, help="nBins on the gluino mass axis", type=int) # (800 - 25) GeV --> (1750 + 25) GeV in steps of 50 GeV
inputArgumentsParser.add_argument("--minGluinoMass", default=775.0, help="Min gluino mass for the 2D plots.", type=float)
inputArgumentsParser.add_argument("--maxGluinoMass", default=1775.0, help="Max gluino mass for the 2D plots.", type=float)
inputArgumentsParser.add_argument("--nNeutralinoMassBins", default=133, help="nBins on the neutralino mass axis.", type=int)
inputArgumentsParser.add_argument("--minNeutralinoMass", default=93.75, help="Min neutralino mass for the 2D plots.", type=float)
inputArgumentsParser.add_argument("--maxNeutralinoMass", default=1756.25, help="Max neutralino mass for the 2D plots.", type=float) # (100 - 6.25) GeV --> (1750 + 6.25) GeV in steps of 12.5 GeV
inputArguments = inputArgumentsParser.parse_args()

photonFailureCategories = ["eta", "pT", "hOverE", "neutralIsolation", "photonIsolation", "conversionSafeElectronVeto", "sigmaietaiataANDchargedIsolation", "sigmaietaiataANDchargedIsolationLoose"]
jetFailureCategories = ["eta", "pT", "puID", "jetID", "deltaR"]
eventFailureCategories = ["HLTPhoton", "wrongNSelectedPhotons", "incompatiblePhotonSelectionType", "lowInvariantMass", "HLTJet", "wrongNJets", "hTCut", "electronVeto", "muonVeto", "MCGenInformation"]
counterTypes = ["global", "differential"]
histograms = {}

for counterType in counterTypes:
    histograms[counterType] = {
        "photon": {},
        "jet": {},
        "event": {}
    }
    for photonFailureCategory in photonFailureCategories:
        histograms[counterType]["photon"][photonFailureCategory] = ROOT.TH2I("photonFailureCounters_MCMap_" + counterType + "_" + photonFailureCategory, counterType + " number of failing photons: " + photonFailureCategory, inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
    for jetFailureCategory in jetFailureCategories:
        histograms[counterType]["jet"][jetFailureCategory] = ROOT.TH2I("jetFailureCounters_MCMap_" + counterType + "_" + jetFailureCategory, counterType + " number of failing jets: " + jetFailureCategory, inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
    for eventFailureCategory in eventFailureCategories:
        histograms[counterType]["event"][eventFailureCategory] = ROOT.TH2I("eventFailureCounters_MCMap_" + counterType + "_" + eventFailureCategory, counterType + " number of failing events: " + eventFailureCategory, inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)

# Load files
inputFileNamesFileObject = open(inputArguments.inputFilesList, 'r')
for inputFileName in inputFileNamesFileObject:
    print("Adding histograms from file: " + inputFileName.strip())
    inputFile = ROOT.TFile.Open(inputFileName.strip(), "READ")
    if ((inputFile.IsZombie() == ROOT.kTRUE) or not(inputFile.IsOpen() == ROOT.kTRUE)):
        sys.exit("Unable to open file with name: " + inputFileName.strip())
    for counterType in counterTypes:
        for photonFailureCategory in photonFailureCategories:
            inputHistogram = ROOT.TH2I("temp", "", inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
            inputFile.GetObject("photonFailureCounters_MCMap_" + counterType + "_" + photonFailureCategory, inputHistogram)
            histograms[counterType]["photon"][photonFailureCategory].Add(inputHistogram)
        for jetFailureCategory in jetFailureCategories:
            inputHistogram = ROOT.TH2I("temp", "", inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
            inputFile.GetObject("jetFailureCounters_MCMap_begin_" + counterType + "_" + jetFailureCategory, inputHistogram)
            histograms[counterType]["jet"][jetFailureCategory].Add(inputHistogram)
        for eventFailureCategory in eventFailureCategories:
            inputHistogram = ROOT.TH2I("temp", "", inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
            inputFile.GetObject("eventFailureCounters_MCMap_" + counterType + "_" + eventFailureCategory, inputHistogram)
            histograms[counterType]["event"][eventFailureCategory].Add(inputHistogram)
    inputFile.Close()
inputFileNamesFileObject.close()

# Save outputs
for counterType in counterTypes:
    for photonFailureCategory in photonFailureCategories:
        tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [histograms[counterType]["photon"][photonFailureCategory]], canvasName = "c_photon_" + counterType + "_" + photonFailureCategory, outputDocumentName=("{oD}/MCSelectionStats_photon_" + counterType + "_" + photonFailureCategory).format(oD=inputArguments.outputDirectory), customPlotOptions_firstObject="COLZ")
    for jetFailureCategory in jetFailureCategories:
        tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [histograms[counterType]["jet"][jetFailureCategory]], canvasName = "c_jet_" + counterType + "_" + jetFailureCategory, outputDocumentName=("{oD}/MCSelectionStats_jet_" + counterType + "_" + jetFailureCategory).format(oD=inputArguments.outputDirectory), customPlotOptions_firstObject="COLZ")
    for eventFailureCategory in eventFailureCategories:
        tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [histograms[counterType]["event"][eventFailureCategory]], canvasName = "c_event_" + counterType + "_" + eventFailureCategory, outputDocumentName=("{oD}/MCSelectionStats_event_" + counterType + "_" + eventFailureCategory).format(oD=inputArguments.outputDirectory), customPlotOptions_firstObject="COLZ")

print("All done!")
