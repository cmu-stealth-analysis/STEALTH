#!/usr/bin/env python

from __future__ import print_function, division

import argparse, os, sys, ROOT, tmROOTUtils

inputArgumentsParser = argparse.ArgumentParser(description="Script to compare two input ST distributions and datasets.")
inputArgumentsParser.add_argument("--addDatasetAndDistribution", action="append", help="Path to dataset in the format \"rootSourceFile:datasetName:distributionName:title:color\".",type=str)
inputArgumentsParser.add_argument("--title", default="Comparison of sT Distributions", help="Title of output plot.",type=str)
inputArgumentsParser.add_argument("--outputDirectory", default="analysis/STComparisons", help="Path to output directory.",type=str)
inputArgumentsParser.add_argument("--outputFile", default="STComparisons", help="Name of output file.",type=str)
inputArgumentsParser.add_argument('--sTPlotRangeMin', default=1100., help='Min value of sT to display in the plots.',type=float)
inputArgumentsParser.add_argument('--sTPlotRangeMax', default=3500., help='Max value of sT to display in the plots.',type=float)
inputArgumentsParser.add_argument('--n_sTBins', default=24, help='Number of sT bins (relevant for plotting only).',type=int)
inputArguments = inputArgumentsParser.parse_args()

colorsDict = {"red": ROOT.kRed, "blue": ROOT.kBlue, "black": ROOT.kBlack}
rooVar_sT = ROOT.RooRealVar("rooVar_sT", "rooVar_sT", inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax, "GeV")
rooVar_sT.setRange("fullRange", inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax)
STRange = ROOT.RooFit.Range(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax)
ST_normRange = ROOT.RooFit.NormRange("fullRange")
binWidth = int(0.5 + (inputArguments.sTPlotRangeMax - inputArguments.sTPlotRangeMin)/inputArguments.n_sTBins)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

inputFiles = {}
inputObjects = {"datasets": {}, "distributions": {}}
legendTitleBases = {}
objectColors = {}

# Get objects
for objectIndex in range(len(inputArguments.addDatasetAndDistribution)):
    objectInputStrings = (inputArguments.addDatasetAndDistribution[objectIndex]).split(":")
    if not(len(objectInputStrings) == 5):
        sys.exit("Five attributes needed for dataset input string. Current list: {l}".format(l = str(objectInputStrings)))
    inputFileName = objectInputStrings[0]
    inputFiles[objectIndex] = ROOT.TFile(inputFileName, "READ")
    inputDatasetName = objectInputStrings[1]
    inputObjects["datasets"][objectIndex] = ROOT.RooDataSet()
    inputFiles[objectIndex].GetObject(inputDatasetName, inputObjects["datasets"][objectIndex])
    inputDistributionName = objectInputStrings[2]
    inputObjects["distributions"][objectIndex] = ROOT.RooKeysPdf()
    inputFiles[objectIndex].GetObject(inputDistributionName, inputObjects["distributions"][objectIndex])
    legendTitleBases[objectIndex] = objectInputStrings[3]
    objectColors[objectIndex] = colorsDict[objectInputStrings[4]]

# Create plot
def customizeLegendEntry(entry, color, objectType):
    entry.ResetAttFill()
    entry.SetTextColor(color)
    entry.SetLineColor(color)
    if (objectType == "dataset"): return
    entry.ResetAttMarker()
    entry.SetMarkerColor(color)
    entry.SetMarkerStyle(8)
    entry.SetMarkerSize(0.125)
    entry.SetFillStyle(0)
    entry.SetFillColor(color)
    entry.SetLineStyle(1)
    entry.SetLineWidth(2)

def setFrameAesthetics(frame, xLabel, yLabel, title):
    frame.SetXTitle(xLabel)
    frame.SetYTitle(yLabel)
    frame.SetTitle(title)

os.system("mkdir -p {oD}".format(oD=inputArguments.outputDirectory))
outputFrame = rooVar_sT.frame(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
setFrameAesthetics(outputFrame, "#it{S}_{T} (GeV)", "Events / ({binWidth} GeV)".format(binWidth=binWidth), inputArguments.title)
outputLegend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)

for objectIndex in range(len(inputArguments.addDatasetAndDistribution)):
    inputObjects["datasets"][objectIndex].plotOn(outputFrame, ROOT.RooFit.LineColor(objectColors[objectIndex]), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.RefreshNorm())
    legendEntry_Dataset = outputLegend.AddEntry(inputObjects["datasets"][objectIndex], "Data: " + legendTitleBases[objectIndex])
    customizeLegendEntry(legendEntry_Dataset, objectColors[objectIndex], "dataset")
    # inputObjects["distributions"][objectIndex].fitTo(inputObjects["datasets"][objectIndex], STRange, ROOT.RooFit.Minos(ROOT.kTRUE), ROOT.RooFit.PrintLevel(0), ROOT.RooFit.SumW2Error(ROOT.kFALSE), ROOT.RooFit.Optimize(0))
    inputObjects["distributions"][objectIndex].plotOn(outputFrame, ROOT.RooFit.LineColor(objectColors[objectIndex]), ST_normRange)
    legendEntry_Distribution = outputLegend.AddEntry(inputObjects["distributions"][objectIndex], "PDF: " + legendTitleBases[objectIndex])
    customizeLegendEntry(legendEntry_Distribution, objectColors[objectIndex], "distribution")

# Save plots, both linear and log scale
tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [outputFrame, outputLegend], canvasName = "c_sTComparisons_linear", outputDocumentName = "{oD}/{oF}_linear".format(oD=inputArguments.outputDirectory, oF=inputArguments.outputFile))
tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [outputFrame, outputLegend], canvasName = "c_sTComparisons", outputDocumentName = "{oD}/{oF}".format(oD=inputArguments.outputDirectory, oF=inputArguments.outputFile), enableLogY = True)

# Cleanup
for objectIndex in range(len(inputArguments.addDatasetAndDistribution)):
    inputFiles[objectIndex].Close()
