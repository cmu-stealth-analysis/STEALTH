#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, array, pdb, math
import ROOT, tmROOTUtils, tdrstyle, CMS_lumi
from tmProgressBar import tmProgressBar

ROOT.gROOT.SetBatch(ROOT.kTRUE)

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Generate plot of comparison of ST distributions at various jet multiplicities.')
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="publicationPlots", help='Output directory.',type=str)
inputArgumentsParser.add_argument('--outputFileName', required=True, help='Name of output file.',type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Min nJets bin.',type=int)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Max nJets bin.',type=int)
inputArgumentsParser.add_argument('--nJetsNorm', default=2, help='Norm nJets bin.',type=int)
inputArgumentsParser.add_argument('--dataSpecialDescription', default="fake #gamma + fake #gamma", help='Special string to describe distributions.',type=str)
inputArguments = inputArgumentsParser.parse_args()

histColors = {
    2: ROOT.kBlack,
    3: ROOT.kBlue+2,
    4: ROOT.kRed+1,
    5: ROOT.kGreen+3,
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
targetSTNorm = 0.5*(STBoundaries[0] + STBoundaries[1])
print("Target norm: {tN}".format(tN=targetSTNorm))

# Load input TTrees into TChain
inputChain = ROOT.TChain("ggNtuplizer/EventTree")
inputChain.Add(inputArguments.inputFilePath)

nEvents = inputChain.GetEntries()
if (nEvents == 0): sys.exit("Number of available events is 0.")

STHistograms = {}
STHistogramsScaled = {}
ratioHistograms = {}
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    STHistograms[nJetsBin] = ROOT.TH1F("h_STDistribution_{n}Jets".format(n=nJetsBin), "", n_STBins, array.array('d', STBoundaries))
    STHistograms[nJetsBin].Sumw2()
    STHistogramsScaled[nJetsBin] = ROOT.TH1F("h_STDistribution_scaled_{n}Jets".format(n=nJetsBin), "", n_STBins, array.array('d', STBoundaries))
    STHistogramsScaled[nJetsBin].Sumw2()
    if (nJetsBin == inputArguments.nJetsNorm): continue
    ratioHistograms[nJetsBin] = ROOT.TH1F("h_ratioDistribution_{n}Jets".format(n=nJetsBin), "", n_STBins, array.array('d', STBoundaries))
    ratioHistograms[nJetsBin].Sumw2()

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

for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    tmROOTUtils.rescale1DHistogramByBinWidth(STHistograms[nJetsBin])
    normBinIndex = STHistograms[nJetsBin].GetXaxis().FindFixBin(targetSTNorm)
    normalizationFactor = STHistograms[inputArguments.nJetsNorm].GetBinContent(normBinIndex)/STHistograms[nJetsBin].GetBinContent(normBinIndex)
    for binIndex in range(1, 1+STHistograms[nJetsBin].GetXaxis().GetNbins()):
        error = STHistograms[nJetsBin].GetBinError(binIndex)
        STHistogramsScaled[nJetsBin].SetBinContent(binIndex, normalizationFactor*STHistograms[nJetsBin].GetBinContent(binIndex))
        STHistogramsScaled[nJetsBin].SetBinError(binIndex, normalizationFactor*STHistograms[nJetsBin].GetBinError(binIndex))
    if (nJetsBin == inputArguments.nJetsNorm): continue
    normBinIndex = STHistograms[nJetsBin].GetXaxis().FindFixBin(targetSTNorm)
    normalizationFactor = STHistograms[inputArguments.nJetsNorm].GetBinContent(normBinIndex)/STHistograms[nJetsBin].GetBinContent(normBinIndex)
    for binIndex in range(1, 1+STHistograms[nJetsBin].GetXaxis().GetNbins()):
        numerator = STHistograms[nJetsBin].GetBinContent(binIndex)
        numeratorError = STHistograms[nJetsBin].GetBinError(binIndex)
        denominator = STHistograms[inputArguments.nJetsNorm].GetBinContent(binIndex)
        denominatorError = STHistograms[inputArguments.nJetsNorm].GetBinError(binIndex)
        ratioHistograms[nJetsBin].SetBinContent(binIndex, 1.)
        ratioHistograms[nJetsBin].SetBinError(binIndex, 1.)
        if (denominator > 0):
            ratio = numerator/denominator
            ratioError = (1.0/denominator)*math.sqrt(math.pow(numeratorError,2) + math.pow(denominatorError*ratio,2))
            ratioHistograms[nJetsBin].SetBinContent(binIndex, normalizationFactor*ratio)
            ratioHistograms[nJetsBin].SetBinError(binIndex, normalizationFactor*ratioError)

tdrstyle.setTDRStyle()

H_ref = 600
W_ref = 800
W = W_ref
H  = H_ref
T = 0.08*H_ref
B = 0.12*H_ref
L = 0.12*W_ref
R = 0.04*W_ref

canvas = ROOT.TCanvas("c_{oFN}_{n}Jets".format(oFN=inputArguments.outputFileName, n=nJetsBin), "c_{oFN}_{n}Jets".format(oFN=inputArguments.outputFileName, n=nJetsBin), 50, 50, W, H)
canvas.SetFillColor(0)
canvas.SetBorderMode(0)
canvas.SetFrameFillStyle(0)
canvas.SetFrameBorderMode(0)
canvas.SetLeftMargin( L/W )
canvas.SetRightMargin( R/W )
canvas.SetTopMargin( T/H )
canvas.SetBottomMargin( B/H )
canvas.SetTickx(0)
canvas.SetTicky(0)
canvas.Draw()

bottomFraction = 0.25
bottomToTopRatio = bottomFraction/(1.0 - bottomFraction)
upperPad = ROOT.TPad("upperPad_{n}Jets".format(n=nJetsBin), "upperPad_{n}Jets".format(n=nJetsBin), 0., bottomFraction, 0.97, 0.97)
upperPad.SetMargin(0.12, 0.03, 0.025, 0.08) # left, right, bottom, top
lowerPad = ROOT.TPad("lowerPad_{n}Jets".format(n=nJetsBin), "lowerPad_{n}Jets".format(n=nJetsBin), 0., 0., 0.97, bottomFraction)
lowerPad.SetMargin(0.12, 0.03, 0.38, 0.03) # left, right, bottom, top
upperPad.Draw()
lowerPad.Draw()

commonTitleOffset = 0.7
# commonFillColor = ROOT.kOrange-2
# commonExpectedEventsLineColor = ROOT.kBlack
# commonExpectedEventsLineStyle = 2
commonLineWidth = 3
commonTitleSize = 0.06
commonLabelSize = 0.05

upperPad.cd()
upperPad.SetLogy()

STHistogramsScaled[inputArguments.nJetsNorm].SetLineColor(histColors[inputArguments.nJetsNorm])
STHistogramsScaled[inputArguments.nJetsNorm].SetLineWidth(commonLineWidth)
STHistogramsScaled[inputArguments.nJetsNorm].GetXaxis().SetTitleSize(commonTitleSize)
STHistogramsScaled[inputArguments.nJetsNorm].GetXaxis().SetLabelSize(commonLabelSize)
STHistogramsScaled[inputArguments.nJetsNorm].GetXaxis().SetTickLength(0)
STHistogramsScaled[inputArguments.nJetsNorm].GetXaxis().SetLabelOffset(999)
STHistogramsScaled[inputArguments.nJetsNorm].GetYaxis().SetTitle("Events/GeV")
STHistogramsScaled[inputArguments.nJetsNorm].GetYaxis().SetTitleSize(commonTitleSize)
STHistogramsScaled[inputArguments.nJetsNorm].GetYaxis().SetLabelSize(commonLabelSize)
STHistogramsScaled[inputArguments.nJetsNorm].GetYaxis().SetTitleOffset(commonTitleOffset)

for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    if (nJetsBin == inputArguments.nJetsNorm): continue
    STHistogramsScaled[nJetsBin].SetLineColor(histColors[nJetsBin])
    STHistogramsScaled[nJetsBin].SetLineWidth(commonLineWidth)

CMS_lumi.writeExtraText = False
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
CMS_lumi.lumi_13TeV = "77.8 fb^{-1}"

legend = ROOT.TLegend(0.8, 0.475, 0.95, 0.9)
legend.SetBorderSize(commonLineWidth)
legend.SetFillStyle(0)
ROOT.gStyle.SetLegendTextSize(0.05)

STHistogramsScaled[inputArguments.nJetsNorm].Draw("P0")
norm_legendEntry = legend.AddEntry(STHistogramsScaled[inputArguments.nJetsNorm], "N_{{Jets}} = {n}".format(n=inputArguments.nJetsNorm), "LPE")
norm_legendEntry.SetTextColor(histColors[inputArguments.nJetsNorm])
norm_legendEntry.SetMarkerColor(histColors[inputArguments.nJetsNorm])
STHistogramsScaled[inputArguments.nJetsNorm].GetXaxis().SetRangeUser(STBoundaries[0], STBoundaries[-1])
STHistogramsScaled[inputArguments.nJetsNorm].GetYaxis().SetRangeUser(0.0002, 11.)

if not(inputArguments.dataSpecialDescription == ""):
    latex = ROOT.TLatex()
    latex.SetTextFont(42)
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextSize(0.07)
    latex.SetTextAlign(22)
    latex.DrawLatexNDC(0.5, 0.8, inputArguments.dataSpecialDescription)

for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    if (nJetsBin == inputArguments.nJetsNorm): continue
    STHistogramsScaled[nJetsBin].Draw("AP0 SAME")
    legendText = "N_{{Jets}} = {n}".format(n=nJetsBin)
    if (nJetsBin == inputArguments.nJetsMax): legendText = "N_{{Jets}} #geq {n}".format(n=nJetsBin)
    legendEntry = legend.AddEntry(STHistogramsScaled[inputArguments.nJetsNorm], legendText, "LPE")
    legendEntry.SetTextColor(histColors[nJetsBin])
    legendEntry.SetMarkerColor(histColors[nJetsBin])
legend.Draw()
CMS_lumi.CMS_lumi(canvas, 4, 0)
upperPad.cd()
upperPad.Update()
upperPad.RedrawAxis()
frame = upperPad.GetFrame()
frame.Draw()

yTitleSize_upper = STHistogramsScaled[inputArguments.nJetsNorm].GetYaxis().GetTitleSize()
yLabelSize_upper = STHistogramsScaled[inputArguments.nJetsNorm].GetYaxis().GetLabelSize()
yTickLength_upper = STHistogramsScaled[inputArguments.nJetsNorm].GetYaxis().GetTickLength()
upperPad.Update()

lowerPad.cd()
isSet = False
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    if (nJetsBin == inputArguments.nJetsNorm): continue
    ratioHistograms[nJetsBin].SetLineColor(histColors[nJetsBin])
    ratioHistograms[nJetsBin].SetLineWidth(commonLineWidth)
    if (isSet):
        ratioHistograms[nJetsBin].Draw("AP0 SAME")
        continue
    ratioHistograms[nJetsBin].GetXaxis().SetTitle("S_{T} (GeV)")
    ratioHistograms[nJetsBin].GetXaxis().SetTitleSize(yTitleSize_upper/bottomToTopRatio)
    ratioHistograms[nJetsBin].GetXaxis().SetLabelSize(yLabelSize_upper/bottomToTopRatio)
    ratioHistograms[nJetsBin].GetXaxis().SetTickLength(yTickLength_upper)
    ratioHistograms[nJetsBin].GetXaxis().SetTitleOffset(0.86)
    ratioHistograms[nJetsBin].GetYaxis().SetTitle("#frac{N_{Jets} bin}{N_{Jets} = " + str(inputArguments.nJetsNorm) + "}")
    ratioHistograms[nJetsBin].GetYaxis().SetTitleOffset(1.4*bottomToTopRatio*commonTitleOffset)
    ratioHistograms[nJetsBin].GetYaxis().SetTitleSize(0.75*yTitleSize_upper/bottomToTopRatio)
    ratioHistograms[nJetsBin].GetYaxis().SetLabelSize(yLabelSize_upper/bottomToTopRatio)
    ratioHistograms[nJetsBin].GetYaxis().SetTickLength(yTickLength_upper)
    ratioHistograms[nJetsBin].GetYaxis().SetNdivisions(2, 0, 0)
    ratioHistograms[nJetsBin].Draw("P0")
    ratioHistograms[nJetsBin].GetXaxis().SetRangeUser(STBoundaries[0], STBoundaries[-1])
    ratioHistograms[nJetsBin].GetYaxis().SetRangeUser(-1., 3.)
    isSet = True
    
lineAt1 = ROOT.TLine(STBoundaries[0], 1., STBoundaries[-1], 1.)
lineAt1.SetLineColor(histColors[inputArguments.nJetsNorm])
lineAt1.SetLineWidth(commonLineWidth)
lineAt1.Draw()
lowerPad.cd()
lowerPad.Update()
lowerPad.RedrawAxis()
frame = lowerPad.GetFrame()
frame.Draw()

canvas.Update()
canvas.SaveAs("{oD}/{oFN}".format(oD=inputArguments.outputDirectory, oFN=inputArguments.outputFileName))

print("All done!")
