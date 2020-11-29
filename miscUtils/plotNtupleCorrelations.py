#!/usr/bin/env python

from __future__ import print_function, division
import os, sys, argparse, pdb, math, json, copy
import tmGeneralUtils, tmProgressBar

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='General tool to generate a CMS-formatted comparison of various histograms; list is read in from an input JSON file whose syntax is explained in the comment immediately following the argument parser setup.')
inputArgumentsParser.add_argument('--jsonPath', required=True, help='Path to input JSON.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', required=True, help='Output directory in which to store the plots.',type=str)
inputArgumentsParser.add_argument('--printTemplate', action='store_true', help="Only print template for a skeleton JSON file and exit.")
inputArguments = inputArgumentsParser.parse_args()

if inputArguments.printTemplate:
    print("Template: ")
    print("""
{
    "targets": {
        "uniqueIDWithDuplicate": {
            "duplicate": "true",
            "semicolonSeparatedDuplicateList": "e1;e2;e3",
            "outputPath": "output_path_{d}.pdf",
            "semicolonSeparatedInputPaths": "path1{d};path2{d};path3{d}",
            "nJetsMin": "2",
            "nJetsMax": "3",
            "STMin": "1200.0",
            "title": "Title {d}",
            "xAxis": {
                "branchName": "bname_x",
                "label": "x",
                "min": "-1.0",
                "max": "1.0",
                "nBins": "50"
            },
            "yAxis": {
                "branchName": "bname_y",
                "label": "y",
                "min": "-1.0",
                "max": "1.0",
                "nBins": "50"
            }
        }
        "uniqueIDNoDuplicate": {
            "duplicate": "false",
            "outputPath": "output_path.pdf",
            "semicolonSeparatedInputPaths": "path1;path2;path3",
            "nJetsMin": "2",
            "nJetsMax": "3",
            "STMin": "1200.0",
            "title": "Title",
            "xAxis": {
                "branchName": "bname_x",
                "label": "x",
                "min": "-1.0",
                "max": "1.0",
                "nBins": "50"
            },
            "yAxis": {
                "branchName": "bname_y",
                "label": "y",
                "min": "-1.0",
                "max": "1.0",
                "nBins": "50"
            }
        }
    }
}
    """)
    sys.exit(0)

# If ROOT is imported before the input arguments parser, the default "help" message is not the right one
import ROOT
# import tdrstyle, CMS_lumi
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

def parse_target(inputData):
    outputDict = {}
    outputDict["outputPath"] = str(inputData["outputPath"])
    outputDict["inputPaths"] = str(inputData["semicolonSeparatedInputPaths"]).split(";")
    for nJetsRangeVar in ["nJetsMin", "nJetsMax"]: outputDict[nJetsRangeVar] = int(0.5 + float(str(inputData[nJetsRangeVar])))
    try:
        outputDict["STMin"] = float(str(inputData["STMin"]))
    except KeyError:
        print("STMin not set explicitly, setting it to default value 1200.0")
        outputDict["STMin"] = 1200.0
    outputDict["title"] = str(inputData["title"])
    for axisLabel in ["xAxis", "yAxis"]:
        outputDict[axisLabel] = {}
        outputDict[axisLabel]["branchName"] = str(inputData[axisLabel]["branchName"])
        outputDict[axisLabel]["label"] = str(inputData[axisLabel]["label"])
        for minMax in ["min", "max"]: outputDict[axisLabel][minMax] = float(str(inputData[axisLabel][minMax]))
        outputDict[axisLabel]["nBins"] = int(0.5 + float(str(inputData[axisLabel]["nBins"])))
    if (str(inputData["duplicate"]) == "true"):
        listToReturn = []
        for substitution in str(inputData["semicolonSeparatedDuplicateList"]).split(";"):
            outputDictCopy = copy.deepcopy(outputDict)
            outputDictCopy["outputPath"] = str(inputData["outputPath"]).format(d=substitution)
            outputDictCopy["inputPaths"] = []
            for inputPath in str(inputData["semicolonSeparatedInputPaths"]).split(";"):
                outputDictCopy["inputPaths"].append(inputPath.format(d=substitution))
            outputDictCopy["title"] = str(inputData["title"]).format(d=substitution)
            listToReturn.append(outputDictCopy)
        return listToReturn
    else:
        return [copy.deepcopy(outputDict)]

def find_and_save_correlation(correlationDetails):
    correlationName = correlationDetails["outputPath"].replace(".pdf", "")
    correlation2DHistogram = ROOT.TH2D(correlationName, correlationDetails["title"], correlationDetails["xAxis"]["nBins"], correlationDetails["xAxis"]["min"], correlationDetails["xAxis"]["max"], correlationDetails["yAxis"]["nBins"], correlationDetails["yAxis"]["min"], correlationDetails["yAxis"]["max"])
    correlation2DHistogram.GetXaxis().SetTitle(correlationDetails["xAxis"]["label"])
    correlation2DHistogram.GetYaxis().SetTitle(correlationDetails["yAxis"]["label"])
    # Load input TTrees into TChain
    inputChain = ROOT.TChain("ggNtuplizer/EventTree")

    for inputFile in correlationDetails["inputPaths"]:
        # print("Adding: {inputfile}".format(inputfile=inputFile))
        readStatus = inputChain.Add(inputFile, 0)
        if not(readStatus == 1): sys.exit("File {iF} does not exist or does not contain the correct tree.".format(iF=inputFile))

    nEvents = inputChain.GetEntries()
    print("Available nEvents: {n}".format(n=nEvents))

    progressBar = tmProgressBar.tmProgressBar(nEvents)
    progressBarUpdatePeriod = max(1, (nEvents//20))
    progressBar.initializeTimer()
    for eventIndex in range(0,nEvents):
        treeStatus = inputChain.LoadTree(eventIndex)
        if treeStatus < 0:
            # sys.exit("Tree unreadable.")
            break
        evtStatus = inputChain.GetEntry(eventIndex)
        if evtStatus <= 0:
            # sys.exit("Event in tree unreadable.")
            continue
        if (eventIndex%progressBarUpdatePeriod == 0 or eventIndex == (nEvents - 1)): progressBar.updateBar(eventIndex/nEvents, eventIndex)

        nJetsBin = inputChain.b_nJets
        if (nJetsBin > correlationDetails["nJetsMax"]): continue
        if (nJetsBin < correlationDetails["nJetsMin"]): continue
        sT = inputChain.b_evtST
        if (sT < correlationDetails["STMin"]): continue
        xValue = getattr(inputChain, correlationDetails["xAxis"]["branchName"])
        yValue = getattr(inputChain, correlationDetails["yAxis"]["branchName"])
        correlation2DHistogram.Fill(xValue, yValue)

    outputCanvas = ROOT.TCanvas("c_{n}".format(n=correlationName), "c_{n}".format(n=correlationName), 1024, 768)
    ROOT.gStyle.SetOptStat(0)
    correlation2DHistogram.Draw("COLZ")
    outputCanvas.SaveAs("{oD}/{oP}".format(oD=inputArguments.outputDirectory, oP=correlationDetails["outputPath"]))

inputFileObject = open(inputArguments.jsonPath, 'r')
inputTargets = json.load(inputFileObject)
inputFileObject.close()

for target in inputTargets["targets"]:
    print("Saving correlations for target: {t}".format(t=target))
    for correlationDetails in parse_target(inputTargets["targets"][target]):
        find_and_save_correlation(correlationDetails)
