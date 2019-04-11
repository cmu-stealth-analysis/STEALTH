#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, pdb, math
import ROOT
import tmROOTUtils

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Merge and analyze trigger efficiency from the event selection outputs.')
inputArgumentsParser.add_argument('--inputFilesList', required=True, help='Path to file containing newline-separated list of input files.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="triggerEfficiency", help='Output directory.',type=str)
inputArgumentsParser.add_argument('--outputSuffix', required=True, help='Output prefix.',type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArguments = inputArgumentsParser.parse_args()

STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries)
nSTBoundaries = 0
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()): nSTBoundaries += 1
nSTSignalBins = nSTBoundaries - 2 + 1 # First two lines are for the normalization bin, last boundary implied at infinity
print("Using {n} signal bins for ST.".format(n = nSTSignalBins))
STRegionBoundariesFileObject.close()

# def plotSmoothed(inputHistogram, outputFileName):
#     tempGraph=ROOT.TGraph2D()
#     generatedMCTemplate = ROOT.TFile(inputArguments.inputFile_MCTemplate)
#     h_MCTemplate = generatedMCTemplate.Get("h_susyMasses_template")
#     for gluinoMassBin in range(1, 1+h_MCTemplate.GetXaxis().GetNbins()):
#         for neutralinoMassBin in range(1, 1+h_MCTemplate.GetYaxis().GetNbins()):
#             if not(int(0.5 + h_MCTemplate.GetBinContent(gluinoMassBin, neutralinoMassBin)) == 1): continue
#             gluinoMass = h_MCTemplate.GetXaxis().GetBinCenter(gluinoMassBin)
#             neutralinoMass = h_MCTemplate.GetYaxis().GetBinCenter(neutralinoMassBin)
#             print("Setting point ({gM}, {nM}): {v}".format(gM=gluinoMass, nM=neutralinoMass, v=inputHistogram.GetBinContent(inputHistogram.FindFixBin(gluinoMass, neutralinoMass))))
#             tempGraph.SetPoint(tempGraph.GetN(), gluinoMass, neutralinoMass, inputHistogram.GetBinContent(inputHistogram.FindFixBin(gluinoMass, neutralinoMass)))
#     tempGraph.SetNpx(160)
#     tempGraph.SetNpy(266)
#     outputHistogram = tempGraph.GetHistogram()
#     tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [outputHistogram], canvasName = "c_photon_" + counterType + "_" + photonFailureCategory, outputDocumentName=outputFileName, customPlotOptions_firstObject="COLZ", customXRange=[inputArguments.plot_minGluinoMass, inputArguments.plot_maxGluinoMass])

# def getFileNameFormatted(fileName):
#     return (fileName.strip().split("/")[-1]).split(".")[0]

histograms = {
    "nTriggeredEvents_cuts": ROOT.TH2I("nTriggeredEvents_cuts_output", "", 1+nSTSignalBins, 0.5, 1.5+nSTSignalBins, 5, 1.5, 6.5),
    "nTriggeredEvents_cutsANDtrigger": ROOT.TH2I("nTriggeredEvents_cutsANDtrigger_output", "", 1+nSTSignalBins, 0.5, 1.5+nSTSignalBins, 5, 1.5, 6.5)
}

print("Before adding histograms...")
tmROOTUtils.printHistogramContents(inputHistogram = histograms["nTriggeredEvents_cuts"])
tmROOTUtils.printHistogramContents(inputHistogram = histograms["nTriggeredEvents_cutsANDtrigger"])

# Load files
inputFileNamesFileObject = open(inputArguments.inputFilesList, 'r')
for inputFileName in inputFileNamesFileObject:
    # formatted_fileName = getFileNameFormatted(inputFileName)
    print("Adding histograms from file: " + inputFileName.strip())
    inputFile = ROOT.TFile.Open(inputFileName.strip(), "READ")
    if ((inputFile.IsZombie() == ROOT.kTRUE) or not(inputFile.IsOpen() == ROOT.kTRUE)):
        sys.exit("Unable to open file with name: " + inputFileName.strip())

    # inputHistogram_nWithCuts = ROOT.TH2I("nTriggeredEvents_cuts", "", 1+nSTSignalBins, 0.5, 1.5+nSTSignalBins, 5, 1.5, 6.5)
    inputHistogram_nWithCuts = ROOT.TH2I()
    inputFile.GetObject("nTriggeredEvents_cuts", inputHistogram_nWithCuts)
    print("Contents of \"nTriggeredEvents_cuts\":")
    tmROOTUtils.printHistogramContents(inputHistogram = inputHistogram_nWithCuts)
    histograms["nTriggeredEvents_cuts"].Add(inputHistogram_nWithCuts)

    # inputHistogram_nWithCutsANDtrigger = ROOT.TH2I("nTriggeredEvents_cutsANDtrigger", "", 1+nSTSignalBins, 0.5, 1.5+nSTSignalBins, 5, 1.5, 6.5)
    inputHistogram_nWithCutsANDtrigger = ROOT.TH2I()
    inputFile.GetObject("nTriggeredEvents_cutsANDtrigger", inputHistogram_nWithCutsANDtrigger)
    print("Contents of \"nTriggeredEvents_cutsANDtrigger\":")
    tmROOTUtils.printHistogramContents(inputHistogram = inputHistogram_nWithCutsANDtrigger)
    histograms["nTriggeredEvents_cutsANDtrigger"].Add(inputHistogram_nWithCutsANDtrigger)

    print("After adding: ")
    tmROOTUtils.printHistogramContents(inputHistogram = histograms["nTriggeredEvents_cuts"])
    tmROOTUtils.printHistogramContents(inputHistogram = histograms["nTriggeredEvents_cutsANDtrigger"])

    inputFile.Close()
inputFileNamesFileObject.close()

# Save outputs

tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[histograms["nTriggeredEvents_cuts"]], canvasName="c_nTriggeredEvents_cuts", outputDocumentName="{oD}/nEvents_noTriggerRequirements_{oS}".format(oD=inputArguments.outputDirectory, oS=inputArguments.outputSuffix), customPlotOptions_firstObject = "TEXTCOLZ", enableLogZ = False)
tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[histograms["nTriggeredEvents_cutsANDtrigger"]], canvasName="c_nTriggeredEvents_cutsANDtrigger", outputDocumentName="{oD}/nEvents_fullSelection_{oS}".format(oD=inputArguments.outputDirectory, oS=inputArguments.outputSuffix), customPlotOptions_firstObject = "TEXTCOLZ", enableLogZ = False)

print("All done!")
