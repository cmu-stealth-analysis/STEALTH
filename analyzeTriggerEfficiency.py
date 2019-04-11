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

histograms = {
    "nTriggeredEvents_cuts": ROOT.TH2I("nTriggeredEvents_cuts_output", "Total nEvents, no HLT selection;ST region index;nJets bin", 1+nSTSignalBins, 0.5, 1.5+nSTSignalBins, 5, 1.5, 6.5),
    "nTriggeredEvents_cutsANDtrigger": ROOT.TH2I("nTriggeredEvents_cutsANDtrigger_output", "Total nEvents, with HLT selection;ST region index;nJets bin", 1+nSTSignalBins, 0.5, 1.5+nSTSignalBins, 5, 1.5, 6.5),
    "triggerEfficiency": ROOT.TH2F("triggerEfficiency_output", "Trigger efficiency;ST region index;nJets bin", 1+nSTSignalBins, 0.5, 1.5+nSTSignalBins, 5, 1.5, 6.5)
}

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
    histograms["nTriggeredEvents_cuts"].Add(inputHistogram_nWithCuts)

    # inputHistogram_nWithCutsANDtrigger = ROOT.TH2I("nTriggeredEvents_cutsANDtrigger", "", 1+nSTSignalBins, 0.5, 1.5+nSTSignalBins, 5, 1.5, 6.5)
    inputHistogram_nWithCutsANDtrigger = ROOT.TH2I()
    inputFile.GetObject("nTriggeredEvents_cutsANDtrigger", inputHistogram_nWithCutsANDtrigger)
    histograms["nTriggeredEvents_cutsANDtrigger"].Add(inputHistogram_nWithCutsANDtrigger)

    inputFile.Close()
inputFileNamesFileObject.close()

# Calculate trigger efficiency
for xBinIndex in range(1, 1+histograms["triggerEfficiency"].GetXaxis().GetNbins()):
    xCenter = histograms["triggerEfficiency"].GetXaxis().GetBinCenter(xBinIndex)
    for yBinIndex in range(1, 1+histograms["triggerEfficiency"].GetYaxis().GetNbins()):
        yCenter = histograms["triggerEfficiency"].GetYaxis().GetBinCenter(yBinIndex)
        numerator = histograms["nTriggeredEvents_cutsANDtrigger"].GetBinContent(histograms["nTriggeredEvents_cutsANDtrigger"].FindFixBin(xCenter, yCenter))
        denominator = histograms["nTriggeredEvents_cuts"].GetBinContent(histograms["nTriggeredEvents_cuts"].FindFixBin(xCenter, yCenter))
        ratio = 0
        error = 0
        if (denominator > 0):
            ratio = numerator/denominator
            error = histograms["nTriggeredEvents_cuts"].GetBinError(histograms["nTriggeredEvents_cuts"].FindFixBin(xCenter, yCenter))/denominator
        histograms["triggerEfficiency"].SetBinContent(histograms["triggerEfficiency"].FindFixBin(xCenter, yCenter), ratio)
        histograms["triggerEfficiency"].SetBinError(histograms["triggerEfficiency"].FindFixBin(xCenter, yCenter), error)

# Save outputs

tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[histograms["nTriggeredEvents_cuts"]], canvasName="c_nTriggeredEvents_cuts", outputDocumentName="{oD}/nEvents_noTriggerRequirements_{oS}".format(oD=inputArguments.outputDirectory, oS=inputArguments.outputSuffix), customPlotOptions_firstObject = "TEXTCOLZ", enableLogZ = False)
tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[histograms["nTriggeredEvents_cutsANDtrigger"]], canvasName="c_nTriggeredEvents_cutsANDtrigger", outputDocumentName="{oD}/nEvents_fullSelection_{oS}".format(oD=inputArguments.outputDirectory, oS=inputArguments.outputSuffix), customPlotOptions_firstObject = "TEXTCOLZ", enableLogZ = False)
# tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[histograms["triggerEfficiency"]], canvasName="c_triggerEfficiency", outputDocumentName=, customTextFormat = "", customPlotOptions_firstObject = "TEXTCOLZ", enableLogZ = False)
c_triggerEfficiency = ROOT.TCanvas("c_trigEff", "c_trigEff", 1024, 768)
c_triggerEfficiency.SetBorderSize(0)
c_triggerEfficiency.SetFrameBorderMode(0)
# c_triggerEfficiency = tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[histograms["triggerEfficiency"]], canvasName="c_triggerEfficiency", outputDocumentName="".format(oD=inputArguments.outputDirectory, oS=inputArguments.outputSuffix), customTextFormat = "", customPlotOptions_firstObject = "TEXTCOLZ", enableLogZ = False)
ROOT.gStyle.SetOptStat(0)
histograms["triggerEfficiency"].Draw("COLZ")
for xBinIndex in range(1, 1+histograms["triggerEfficiency"].GetXaxis().GetNbins()):
    xCenter = histograms["triggerEfficiency"].GetXaxis().GetBinCenter(xBinIndex)
    for yBinIndex in range(1, 1+histograms["triggerEfficiency"].GetYaxis().GetNbins()):
        yCenter = histograms["triggerEfficiency"].GetYaxis().GetBinCenter(yBinIndex)
        numerator = histograms["nTriggeredEvents_cutsANDtrigger"].GetBinContent(histograms["nTriggeredEvents_cutsANDtrigger"].FindFixBin(xCenter, yCenter))
        denominator = histograms["nTriggeredEvents_cuts"].GetBinContent(histograms["nTriggeredEvents_cuts"].FindFixBin(xCenter, yCenter))
        ratio = 0
        error = 0
        if (denominator > 0):
            ratio = numerator/denominator
            error = histograms["nTriggeredEvents_cuts"].GetBinError(histograms["nTriggeredEvents_cuts"].FindFixBin(xCenter, yCenter))/denominator
        textContent = "{r:.2f} +/- {e:.2f}".format(r=ratio,e=error)
        text = ROOT.TText(xCenter, yCenter, textContent)
        text.SetTextAlign(22)
        text.Draw()
c_triggerEfficiency.SaveAs("{oD}/triggerEfficiency_{oS}.png".format(oD=inputArguments.outputDirectory, oS=inputArguments.outputSuffix))
print("All done!")
