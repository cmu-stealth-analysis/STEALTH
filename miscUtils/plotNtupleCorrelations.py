#!/usr/bin/env python

from __future__ import print_function, division
import os, sys, argparse, pdb, math, json
import tmGeneralUtils, tmProgressBar

# Register command line options
# inputArgumentsParser = argparse.ArgumentParser(description='General tool to generate a CMS-formatted comparison of various histograms; list is read in from an input JSON file whose syntax is explained in the comment immediately following the argument parser setup.')
# inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input JSON.',type=str)
# inputArgumentsParser.add_argument('--outputDirectory', required=True, help='Output directory in which to store the plots.',type=str)
# inputArguments = inputArgumentsParser.parse_args()

# If ROOT is imported before the input arguments parser, the default "help" message is not the right one
import ROOT
# import tdrstyle, CMS_lumi
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

listOfInputFiles = ["root://cmseos.fnal.gov//store/user/lpcsusystealth/selections/combined_DoublePhoton/merged_selection_data_2016_signal.root", "root://cmseos.fnal.gov//store/user/lpcsusystealth/selections/combined_DoublePhoton/merged_selection_data_2017_signal.root", "root://cmseos.fnal.gov//store/user/lpcsusystealth/selections/combined_DoublePhoton/merged_selection_data_2018_signal.root"]
maxNJets = 3
STMin = 1200.
outputPath = "~/nobackup/analysisAreas/correlations/mva_signal_data.pdf"

correlation_mva = ROOT.TH2D("MVACorr", "MVA correlation;leading;subleading", 50, -1., 1., 50, -1., 1.)

# inputFileNamesFileObject = open(inputArguments.inputFilesList, 'r')
# for inputFileName in inputFileNamesFileObject:
#     listOfInputFiles.append(inputFileName.strip())
# inputFileNamesFileObject.close()

# Load input TTrees into TChain
inputChain = ROOT.TChain("ggNtuplizer/EventTree")

for inputFile in listOfInputFiles:
    # print("Adding: {inputfile}".format(inputfile=inputFile))
    readStatus = inputChain.Add(inputFile, 0)
    if not(readStatus == 1): sys.exit("File {iF} does not exist or does not contain the correct tree.".format(iF=inputFile))

nEvents = inputChain.GetEntries()
print("Available nEvents: {n}".format(n=nEvents))

progressBar = tmProgressBar.tmProgressBar(nEvents)
progressBarUpdatePeriod = max(1, (nEvents//1000))
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
    if (nJetsBin > maxNJets): continue
    sT = inputChain.b_evtST
    if (sT < STMin): continue
    leadingPhotonMVA = inputChain.b_photonMVA_leading
    subLeadingPhotonMVA = inputChain.b_photonMVA_subLeading
    correlation_mva.Fill(leadingPhotonMVA, subLeadingPhotonMVA)

outputCanvas = ROOT.TCanvas("c_mva", "c_mva", 1024, 768)
ROOT.gStyle.SetOptStat(0)
correlation_mva.Draw("COLZ")
outputCanvas.SaveAs(outputPath)
