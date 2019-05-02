#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, ROOT, tmROOTUtils, array, pdb, math
from tmProgressBar import tmProgressBar

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Merge and analyze MC selection statistics.')
inputArgumentsParser.add_argument('--inputFilesList', required=True, help='Path to file containing newline-separated list of input files.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="MCSelectionStatistics", help='Output directory.',type=str)
inputArgumentsParser.add_argument('--inputFile_MCTemplate', default="plot_susyMasses_template.root", help='Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.', type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--plotPhotonFailures', action='store_true', help="Plot photon failure statistics.")
inputArgumentsParser.add_argument('--plotJetFailures', action='store_true', help="Plot jet failure statistics.")
inputArgumentsParser.add_argument('--plotEventFailures', action='store_true', help="Plot event failure statistics.")
inputArgumentsParser.add_argument('--plotAcceptanceRatios', action='store_true', help="Plot acceptance ratio statistics.")
inputArgumentsParser.add_argument("--nGluinoMassBins", default=20, help="nBins on the gluino mass axis", type=int) # (800 - 25) GeV --> (1750 + 25) GeV in steps of 50 GeV
inputArgumentsParser.add_argument("--minGluinoMass", default=775.0, help="Min gluino mass for the 2D plots.", type=float)
inputArgumentsParser.add_argument("--maxGluinoMass", default=1775.0, help="Max gluino mass for the 2D plots.", type=float)
inputArgumentsParser.add_argument("--plot_minGluinoMass", default=1000.0, help="Min gluino mass for the 2D plots (to use while plotting only).", type=float)
inputArgumentsParser.add_argument("--plot_maxGluinoMass", default=1775.0, help="Max gluino mass for the 2D plots (to use while plotting only).", type=float)
inputArgumentsParser.add_argument("--nNeutralinoMassBins", default=133, help="nBins on the neutralino mass axis.", type=int)
inputArgumentsParser.add_argument("--minNeutralinoMass", default=93.75, help="Min neutralino mass for the 2D plots.", type=float)
inputArgumentsParser.add_argument("--maxNeutralinoMass", default=1756.25, help="Max neutralino mass for the 2D plots.", type=float) # (100 - 6.25) GeV --> (1750 + 6.25) GeV in steps of 12.5 GeV
inputArguments = inputArgumentsParser.parse_args()

gluinoAxis = ROOT.TAxis(inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass)
neutralinoAxis = ROOT.TAxis(inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
generatedMCTemplate = ROOT.TFile(inputArguments.inputFile_MCTemplate)
h_MCTemplate = generatedMCTemplate.Get("h_susyMasses_template")

def plotSmoothed(inputHistogram, outputDirectory, outputFileName):
    tempGraph=ROOT.TGraph2D()
    for gluinoMassBin in range(1, 1+h_MCTemplate.GetXaxis().GetNbins()):
        for neutralinoMassBin in range(1, 1+h_MCTemplate.GetYaxis().GetNbins()):
            if not(int(0.5 + h_MCTemplate.GetBinContent(gluinoMassBin, neutralinoMassBin)) == 1): continue
            gluinoMass = h_MCTemplate.GetXaxis().GetBinCenter(gluinoMassBin)
            neutralinoMass = h_MCTemplate.GetYaxis().GetBinCenter(neutralinoMassBin)
            # print("Setting point ({gM}, {nM}): {v}".format(gM=gluinoMass, nM=neutralinoMass, v=inputHistogram.GetBinContent(inputHistogram.FindFixBin(gluinoMass, neutralinoMass))))
            tempGraph.SetPoint(tempGraph.GetN(), gluinoMass, neutralinoMass, inputHistogram.GetBinContent(inputHistogram.FindFixBin(gluinoMass, neutralinoMass)))
    tempGraph.SetNpx(160)
    tempGraph.SetNpy(266)
    outputHistogram = tempGraph.GetHistogram()
    tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [outputHistogram], canvasName = "c_" + outputFileName, outputDocumentName="{oD}/{oFN}".format(oD=outputDirectory, oFN=outputFileName), customPlotOptions_firstObject="COLZ", customXRange=[inputArguments.plot_minGluinoMass, inputArguments.plot_maxGluinoMass])

def getFileNameFormatted(fileName):
    return (fileName.strip().split("/")[-1]).split(".")[0]

photonSelectionCriteria = ["eta", "pT", "hOverE", "neutralIsolation", "photonIsolation", "conversionSafeElectronVeto", "sigmaietaiataANDchargedIsolation", "sigmaietaiataANDchargedIsolationLoose"]
jetSelectionCriteria = ["eta", "pT", "puID", "jetID", "deltaR"]
eventSelectionCriteria = ["HLTPhoton", "lowEnergyPhotons", "wrongNMediumPhotons", "lowInvariantMass", "wrongNJets", "hTCut", "electronVeto", "muonVeto", "MCGenInformation"]
counterTypes = ["global", "differential"]
histograms = {}

STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries)
nSTBoundaries = 0
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()): nSTBoundaries += 1
nSTSignalBins = nSTBoundaries - 2 + 1 # First two lines are for the normalization bin, last boundary implied at infinity
print("Using {n} signal bins for ST.".format(n = nSTSignalBins))
STRegionBoundariesFileObject.close()

# Initialize histograms
for counterType in counterTypes:
    histograms[counterType] = {
        "photon": {},
        "jet": {},
        "event": {}
    }
    for photonSelectionCriterion in photonSelectionCriteria:
        histograms[counterType]["photon"][photonSelectionCriterion] = ROOT.TH2I("photonFailureCounters_MCMap_" + counterType + "_" + photonSelectionCriterion, counterType + " number of failing photons: " + photonSelectionCriterion, inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
    for jetSelectionCriterion in jetSelectionCriteria:
        histograms[counterType]["jet"][jetSelectionCriterion] = ROOT.TH2I("jetFailureCounters_MCMap_" + counterType + "_" + jetSelectionCriterion, counterType + " number of failing jets: " + jetSelectionCriterion, inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
    for eventSelectionCriterion in eventSelectionCriteria:
        histograms[counterType]["event"][eventSelectionCriterion] = ROOT.TH2I("eventFailureCounters_MCMap_" + counterType + "_" + eventSelectionCriterion, counterType + " number of failing events: " + eventSelectionCriterion, inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)

histograms["acceptance"] = {
    "Truth": {},
    "Selection": {}
}

for STRegionIndex in range(1, 2 + nSTSignalBins):
    histograms["acceptance"]["Truth"][STRegionIndex] = ROOT.TH2I("acceptance_MCMap_eventPassesTruth_STRegion{i}".format(i=STRegionIndex), "", inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
    histograms["acceptance"]["Selection"][STRegionIndex] = ROOT.TH2I("acceptance_MCMap_eventPassesSelection_STRegion{i}".format(i=STRegionIndex), "", inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)

histograms["total"] = {
    "photon": ROOT.TH2I("photonTotalCounters_MCMap", "", inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass),
    "jet": ROOT.TH2I("jetTotalCounters_MCMap", "", inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass),
    "event": ROOT.TH2I("eventTotalCounters_MCMap", "", inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
}

# Load source files and add to histograms from each source
inputFileNamesFileObject = open(inputArguments.inputFilesList, 'r')
for inputFileName in inputFileNamesFileObject:
    formatted_fileName = getFileNameFormatted(inputFileName)
    print("Adding histograms from file: " + inputFileName.strip())
    inputFile = ROOT.TFile.Open(inputFileName.strip(), "READ")
    if ((inputFile.IsZombie() == ROOT.kTRUE) or not(inputFile.IsOpen() == ROOT.kTRUE)):
        sys.exit("Unable to open file with name: " + inputFileName.strip())
    for counterType in counterTypes:
        for photonSelectionCriterion in photonSelectionCriteria:
            inputName = ("photonFailureCounters_MCMap_" + counterType + "_" + photonSelectionCriterion)
            inputHistogram = ROOT.TH2I(inputName + "_" + formatted_fileName, "", inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
            inputFile.GetObject(inputName, inputHistogram)
            histograms[counterType]["photon"][photonSelectionCriterion].Add(inputHistogram)
        for jetSelectionCriterion in jetSelectionCriteria:
            inputName = ("jetFailureCounters_MCMap_" + counterType + "_" + jetSelectionCriterion)
            inputHistogram = ROOT.TH2I(inputName + "_" + formatted_fileName, "", inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
            inputFile.GetObject(inputName, inputHistogram)
            histograms[counterType]["jet"][jetSelectionCriterion].Add(inputHistogram)
        for eventSelectionCriterion in eventSelectionCriteria:
            inputName = ("eventFailureCounters_MCMap_" + counterType + "_" + eventSelectionCriterion)
            inputHistogram = ROOT.TH2I(inputName + "_" + formatted_fileName, "", inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
            inputFile.GetObject(inputName, inputHistogram)
            histograms[counterType]["event"][eventSelectionCriterion].Add(inputHistogram)
    for STRegionIndex in range(1, 2 + nSTSignalBins):
        for acceptanceType in ["Truth", "Selection"]:
            inputName = "acceptance_MCMap_eventPasses{t}_STRegion{i}".format(t=acceptanceType, i=STRegionIndex)
            inputHistogram = ROOT.TH2I(inputName + "_" + formatted_fileName, "", inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
            inputFile.GetObject(inputName, inputHistogram)
            histograms["acceptance"][acceptanceType][STRegionIndex].Add(inputHistogram)
    for objectType in ["photon", "jet", "event"]:
        inputName = objectType + "TotalCounters_MCMap"
        inputHistogram = ROOT.TH2I(inputName + "_" + formatted_fileName, "", inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
        inputFile.GetObject(inputName, inputHistogram)
        histograms["total"][objectType].Add(inputHistogram)
    inputFile.Close()
inputFileNamesFileObject.close()

acceptanceRatios = {}
for STRegionIndex in range(1, 2 + nSTSignalBins):
    acceptanceRatios[STRegionIndex] = tmROOTUtils.getRatioHistogram(numeratorHistogram=histograms["acceptance"]["Selection"][STRegionIndex], denominatorHistogram=histograms["acceptance"]["Truth"][STRegionIndex], name="acceptanceRatio_MCMap_STRegion{i}".format(t=acceptanceType, i=STRegionIndex))

# Save outputs
for counterType in counterTypes:
    if (inputArguments.plotPhotonFailures):
        for photonSelectionCriterion in photonSelectionCriteria:
            plotSmoothed(histograms[counterType]["photon"][photonSelectionCriterion], "{oD}".format(oD=inputArguments.outputDirectory), ("MCSelectionStats_photon_" + counterType + "_" + photonSelectionCriterion))
    if (inputArguments.plotJetFailures):
        for jetSelectionCriterion in jetSelectionCriteria:
            plotSmoothed(histograms[counterType]["jet"][jetSelectionCriterion], "{oD}".format(oD=inputArguments.outputDirectory), ("MCSelectionStats_jet_" + counterType + "_" + jetSelectionCriterion))
    if (inputArguments.plotEventFailures):
        for eventSelectionCriterion in eventSelectionCriteria:
            # plotSmoothed(histograms[counterType]["event"][eventSelectionCriterion], "{oD}".format(oD=inputArguments.outputDirectory), ("{oD}/MCSelectionStats_event_" + counterType + "_" + eventSelectionCriterion))
            acceptanceFraction = ROOT.TH2I("acceptanceFraction_MCMap_" + counterType + "_" + eventSelectionCriterion, counterType + "fraction: " + eventSelectionCriterion, inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
            for gluinoMassBin in range(1, 1+h_MCTemplate.GetXaxis().GetNbins()):
                for neutralinoMassBin in range(1, 1+h_MCTemplate.GetYaxis().GetNbins()):
                    if not(int(0.5 + h_MCTemplate.GetBinContent(gluinoMassBin, neutralinoMassBin)) == 1): continue
                    gluinoMass = h_MCTemplate.GetXaxis().GetBinCenter(gluinoMassBin)
                    neutralinoMass = h_MCTemplate.GetYaxis().GetBinCenter(neutralinoMassBin)
                    nGoodEvents = histograms["total"]["event"].GetBinContent(histograms["total"]["event"].FindFixBin(gluinoMass, neutralinoMass))
                    numerator = nGoodEvents - histograms[counterType]["event"][eventSelectionCriterion].GetBinContent(histograms[counterType]["event"][eventSelectionCriterion].FindFixBin(gluinoMass, neutralinoMass))
                    denominator = nGoodEvents
                    if (nGoodEvents > 0): acceptanceFraction.SetBinContent(acceptanceFraction.FindFixBin(gluinoMass, neutralinoMass), numerator/denominator)
            tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [acceptanceFraction], canvasName = "c_acceptanceFraction_" + counterType + "_" + eventSelectionCriterion, outputDocumentName=("{oD}/MCSelectionStats_acceptanceFraction_" + counterType + "_" + eventSelectionCriterion).format(oD=inputArguments.outputDirectory), customPlotOptions_firstObject="COLZ", customXRange=[inputArguments.plot_minGluinoMass, inputArguments.plot_maxGluinoMass])

if (inputArguments.plotAcceptanceRatios):
    for STRegionIndex in range(1, 2 + nSTSignalBins):
        plotSmoothed(acceptanceRatios[STRegionIndex], "{oD}".format(oD=inputArguments.outputDirectory), ("{oD}/MCSelectionStats_acceptanceRatio_STRegion{i}").format(i=STRegionIndex))

print("All done!")
