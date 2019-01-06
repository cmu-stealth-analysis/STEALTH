#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, ROOT, tmROOTUtils, array, pdb
from tmProgressBar import tmProgressBar

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Generate simple plot of ST distributions.')
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="STDistributionComparisons", help='Output directory.',type=str)
inputArgumentsParser.add_argument('--outputFileName', required=True, help='Name of output file.',type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--targetSTNorm', default=1150., help='Value of ST at which to normalize all histograms.',type=float)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Min nJets bin.',type=int)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Max nJets bin.',type=int)
inputArgumentsParser.add_argument('--nJetsNorm', default=2, help='Norm nJets bin.',type=int)
inputArguments = inputArgumentsParser.parse_args()

histColors = {
    2: ROOT.kBlack,
    3: ROOT.kBlue,
    4: ROOT.kRed,
    5: ROOT.kGreen,
    6: ROOT.kViolet
}

STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries)
STBoundaries = []
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()):
        STBoundary = float(STBoundaryString.strip())
        STBoundaries.append(STBoundary)
STBoundaries.append(3500.0) # Instead of infinity
n_STBins = len(STBoundaries) - 1
STRegionsAxis = ROOT.TAxis(n_STBins, array.array('d', STBoundaries))

# Load input TTrees into TChain
inputChain = ROOT.TChain("ggNtuplizer/EventTree")
inputChain.Add(inputArguments.inputFilePath)

nEvents = inputChain.GetEntries()
if (nEvents == 0): sys.exit("Number of available events is 0.")

STHistograms = {}
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    STHistograms[nJetsBin] = ROOT.TH1F("h_STDistribution_{n}Jets".format(n=nJetsBin), "ST Distribution, nJets = {n};ST (GeV);A.U.".format(n = nJetsBin), n_STBins, array.array('d', STBoundaries))
    STHistograms[nJetsBin].Sumw2()

progressBar = tmProgressBar(nEvents)
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
    if (nJetsBin > inputArguments.nJetsMax): nJetsBin = inputArguments.nJetsMax
    if (nJetsBin < inputArguments.nJetsMin): sys.exit("Unexpected nJetsBin = {nJetsBin} at entry index = {entryIndex}".format(nJetsBin=nJetsBin, entryIndex=entryIndex))

    sT = inputChain.b_evtST
    STHistograms[nJetsBin].Fill(sT)
progressBar.terminate()

outputCanvas = ROOT.TCanvas("outputCanvas", "outputCanvas", 1024, 768)
upperPad = ROOT.TPad("upperPad", "upperPad", 0.0, 0.4, 1.0, 1.0)
lowerPad = ROOT.TPad("lowerPad", "lowerPad", 0.0, 0.0, 1.0, 0.4)
upperPad.SetBottomMargin(0.025)
lowerPad.SetTopMargin(0.025)
upperPad.Draw()
lowerPad.Draw()
upperPad.cd()
outputLegend = ROOT.TLegend(0.8, 0.65, 0.9, 0.9)
ROOT.gPad.SetLogy()
ROOT.gStyle.SetOptStat(0)
STHistograms[inputArguments.nJetsNorm].SetLineColor(histColors[inputArguments.nJetsNorm])
STHistograms[inputArguments.nJetsNorm].SetTitle("Comparison of ST Distributions")
STHistograms[inputArguments.nJetsNorm].Draw("HIST E1")
STHistograms[inputArguments.nJetsNorm].GetXaxis().SetLabelOffset(999);
STHistograms[inputArguments.nJetsNorm].GetXaxis().SetLabelSize(0);
legendEntry = outputLegend.AddEntry(STHistograms[inputArguments.nJetsNorm], "nJets = {n}".format(n = inputArguments.nJetsNorm))
legendEntry.SetTextColor(histColors[inputArguments.nJetsNorm])
legendEntry.SetLineColor(histColors[inputArguments.nJetsNorm])
nTargetEntries_normBin = STHistograms[inputArguments.nJetsNorm].GetBinContent(STHistograms[inputArguments.nJetsNorm].FindFixBin(inputArguments.targetSTNorm))
ratioHistograms = {}
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    if (nJetsBin == inputArguments.nJetsNorm): continue
    nEntries_normBin = STHistograms[nJetsBin].GetBinContent(STHistograms[nJetsBin].FindFixBin(inputArguments.targetSTNorm))
    STHistograms[nJetsBin].Scale(nTargetEntries_normBin/nEntries_normBin)
    ratioHistograms[nJetsBin] = ROOT.TH1F("h_STDistributionsRatio_{n}Jets".format(n=nJetsBin), "Ratio of ST Distributions, nJets = {n};ST (GeV);ratio".format(n = nJetsBin), n_STBins, array.array('d', STBoundaries))
    ratioHistograms[nJetsBin].Divide(STHistograms[nJetsBin], STHistograms[inputArguments.nJetsNorm])
    STHistograms[nJetsBin].SetLineColor(histColors[nJetsBin])
    STHistograms[nJetsBin].Draw("HIST E1 SAME")
    legendEntry = ROOT.TLegendEntry()
    if (nJetsBin == inputArguments.nJetsMax): legendEntry = outputLegend.AddEntry(STHistograms[nJetsBin], "nJets #geq {n}".format(n = inputArguments.nJetsMax))
    else: legendEntry = outputLegend.AddEntry(STHistograms[nJetsBin], "nJets = {n}".format(n = nJetsBin))
    legendEntry.SetTextColor(histColors[nJetsBin])
    legendEntry.SetLineColor(histColors[nJetsBin])
outputLegend.Draw()
lowerPad.cd()
isFirstToBeDrawn = True
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    if (nJetsBin == inputArguments.nJetsNorm): continue
    ratioHistograms[nJetsBin].SetLineColor(histColors[nJetsBin])
    if (isFirstToBeDrawn):
        ratioHistograms[nJetsBin].SetTitle("")
        ratioHistograms[nJetsBin].Draw("E1")
        ratioHistograms[nJetsBin].GetYaxis().SetRangeUser(0.0, 3.0)
        isFirstToBeDrawn = False
    else:
        ratioHistograms[nJetsBin].Draw("E1 SAME")
lineAt1 = ROOT.TLine(STBoundaries[0], 1.0, STBoundaries[-1], 1.0)
lineAt1.SetLineColor(histColors[inputArguments.nJetsNorm])
lineAt1.Draw()
outputCanvas.SaveAs("{oD}/{oF}".format(oD=inputArguments.outputDirectory, oF=inputArguments.outputFileName))

print("All done!")
