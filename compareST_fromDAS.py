#!/usr/bin/env python
from __future__ import print_function, division

import os, sys, argparse
import numpy as np
import ROOT
from tmProgressBar import tmProgressBar
from tmROOTUtils import normalizeHistogam

parser=argparse.ArgumentParser(description = "sT comparison generator.")
parser.add_argument('--inputFilePath', action='store', help='Input file path', type=str, required=True)
parser.add_argument('--nSTBins', action='store', help='number of sT bins', type=int, default=15)
parser.add_argument('--lowEdgeFirstSTBin', action='store', help='Lower edge of first sT bin', type=float, default=900.)
parser.add_argument('--upEdgeLastSTBin', action='store', help='Upper edge of last sT bin', type=float, default=3400.)
inputArguments = parser.parse_args()

histBinningFormattingString = "({nSTBins}, {lowEdgeFirstSTBin}, {upEdgeLastSTBin})".format(nSTBins=inputArguments.nSTBins, lowEdgeFirstSTBin=inputArguments.lowEdgeFirstSTBin, upEdgeLastSTBin=inputArguments.upEdgeLastSTBin)

binColors = {2: ROOT.kBlack, 3: ROOT.kBlue, 4: ROOT.kRed, 5: ROOT.kPink, 6: ROOT.kViolet}

chain_in = ROOT.TChain('ggNtuplizer/EventTree')
chain_in.Add(inputArguments.inputFilePath)
n_entries = chain_in.GetEntries()
print('Total number of events: ' + str(n_entries))

histStack = ROOT.THStack("sTShape", "Shape of sT in various nJets bins")
legend = ROOT.TLegend(0.7, 0.6, 0.9, 0.9, "Legend:")
for nJets in range(2, 5):
    nameOfHistogram = "sT_nJetsIs" + str(nJets)
    chain_in.Draw("sT>>" + nameOfHistogram + histBinningFormattingString, "nJets == " + str(nJets))
    histogramToAdd = ROOT.gDirectory.Get(nameOfHistogram)
    normalizeHistogam(histogramToAdd)
    histogramToAdd.SetLineColor(binColors[nJets])
    histStack.Add(histogramToAdd)
    legend.AddEntry(histogramToAdd, "nJets = " + str(nJets))
chain_in.Draw("sT>>sT_nJetsIsGreaterThan6" + histBinningFormattingString, "nJets >= 6")
histogramToAdd = ROOT.gDirectory.Get("sT_nJetsIsGreaterThan6")
normalizeHistogam(histogramToAdd)
histogramToAdd.SetLineColor(binColors[6])
histStack.Add(histogramToAdd)
legend.AddEntry(histogramToAdd, "nJets >= 6")
outputCanvas = ROOT.TCanvas("STComparison", "ST Comparison", 1024, 768)
outputCanvas.cd()
histStack.Draw("nostack")
legend.Draw()
outputCanvas.SaveAs("sTComparison.png")
