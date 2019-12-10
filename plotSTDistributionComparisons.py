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
inputArgumentsParser.add_argument('--outputFilePrefix', required=True, help='Prefix for output files.',type=str)
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
STHistogramTypes = ["total", "photon", "jet", "MET"]

STBoundaries = {}
STRegionsAxes = {}
targetSTNorms = {}
n_STBins = {}
for STHistogramType in STHistogramTypes:
    STBoundaries[STHistogramType] = []
    STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries)
    for STBoundaryString in STRegionBoundariesFileObject:
        if (STBoundaryString.strip()):
            STBoundary = float(STBoundaryString.strip())
            if (STHistogramType == "total"):
                STBoundaries[STHistogramType].append(STBoundary)
            else:
                STBoundaries[STHistogramType].append(0.5*STBoundary) # For the various components of ST, halve the bin sizes
    if (STHistogramType == "total"):
        STBoundaries[STHistogramType].append(3500.0) # Instead of infinity
    else:
        STBoundaries[STHistogramType].append(0.5*3500.0) # Instead of infinity
    n_STBins[STHistogramType] = len(STBoundaries[STHistogramType]) - 1
    STRegionsAxes[STHistogramType] = ROOT.TAxis(n_STBins[STHistogramType], array.array('d', STBoundaries[STHistogramType]))
    targetSTNorms[STHistogramType] = 0.5*(STBoundaries[STHistogramType][0] + STBoundaries[STHistogramType][1])
    print("For histogram type: {t}, target norm: {tN}".format(t=STHistogramType, tN=targetSTNorms[STHistogramType]))

# Load input TTrees into TChain
inputChain = ROOT.TChain("ggNtuplizer/EventTree")
inputChain.Add(inputArguments.inputFilePath)

nEvents = inputChain.GetEntries()
if (nEvents == 0): sys.exit("Number of available events is 0.")

STHistograms = {}
STHistogramsScaled = {}
ratioHistograms = {}
for STHistogramType in STHistogramTypes:
    STHistograms[STHistogramType] = {}
    STHistogramsScaled[STHistogramType] = {}
    ratioHistograms[STHistogramType] = {}
    for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
        STHistograms[STHistogramType][nJetsBin] = ROOT.TH1F("h_STDistribution_{t}_{n}Jets".format(t=STHistogramType, n=nJetsBin), "", n_STBins[STHistogramType], array.array('d', STBoundaries[STHistogramType]))
        STHistograms[STHistogramType][nJetsBin].Sumw2()
        STHistogramsScaled[STHistogramType][nJetsBin] = ROOT.TH1F("h_STDistribution_{t}_scaled_{n}Jets".format(t=STHistogramType, n=nJetsBin), "", n_STBins[STHistogramType], array.array('d', STBoundaries[STHistogramType]))
        STHistogramsScaled[STHistogramType][nJetsBin].Sumw2()
        if (nJetsBin == inputArguments.nJetsNorm): continue
        ratioHistograms[STHistogramType][nJetsBin] = ROOT.TH1F("h_ratioDistribution_{t}_{n}Jets".format(t=STHistogramType, n=nJetsBin), "", n_STBins[STHistogramType], array.array('d', STBoundaries[STHistogramType]))
        ratioHistograms[STHistogramType][nJetsBin].Sumw2()

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
    STHistograms["total"][nJetsBin].Fill(sT)
    sT_EM = inputChain.b_evtST_electromagnetic
    STHistograms["photon"][nJetsBin].Fill(sT_EM)
    sT_hadronic = inputChain.b_evtST_hadronic
    STHistograms["jet"][nJetsBin].Fill(sT_hadronic)
    sT_MET = inputChain.b_evtST_MET
    STHistograms["MET"][nJetsBin].Fill(sT_MET)
progressBar.terminate()

for STHistogramType in STHistogramTypes:
    for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
        tmROOTUtils.rescale1DHistogramByBinWidth(STHistograms[STHistogramType][nJetsBin])
        normBinIndex = STHistograms[STHistogramType][nJetsBin].GetXaxis().FindFixBin(targetSTNorms[STHistogramType])
        try:
            normalizationFactor = STHistograms[STHistogramType][inputArguments.nJetsNorm].GetBinContent(normBinIndex)/STHistograms[STHistogramType][nJetsBin].GetBinContent(normBinIndex)
        except ZeroDivisionError:
            continue
        for binIndex in range(1, 1+STHistograms[STHistogramType][nJetsBin].GetXaxis().GetNbins()):
            error = STHistograms[STHistogramType][nJetsBin].GetBinError(binIndex)
            STHistogramsScaled[STHistogramType][nJetsBin].SetBinContent(binIndex, normalizationFactor*STHistograms[STHistogramType][nJetsBin].GetBinContent(binIndex))
            STHistogramsScaled[STHistogramType][nJetsBin].SetBinError(binIndex, normalizationFactor*STHistograms[STHistogramType][nJetsBin].GetBinError(binIndex))
        if (nJetsBin == inputArguments.nJetsNorm): continue
        normBinIndex = STHistograms[STHistogramType][nJetsBin].GetXaxis().FindFixBin(targetSTNorms[STHistogramType])
        try:
            normalizationFactor = STHistograms[STHistogramType][inputArguments.nJetsNorm].GetBinContent(normBinIndex)/STHistograms[STHistogramType][nJetsBin].GetBinContent(normBinIndex)
        except ZeroDivisionError:
            continue
        for binIndex in range(1, 1+STHistograms[STHistogramType][nJetsBin].GetXaxis().GetNbins()):
            numerator = STHistograms[STHistogramType][nJetsBin].GetBinContent(binIndex)
            numeratorError = STHistograms[STHistogramType][nJetsBin].GetBinError(binIndex)
            denominator = STHistograms[STHistogramType][inputArguments.nJetsNorm].GetBinContent(binIndex)
            denominatorError = STHistograms[STHistogramType][inputArguments.nJetsNorm].GetBinError(binIndex)
            ratioHistograms[STHistogramType][nJetsBin].SetBinContent(binIndex, 1.)
            ratioHistograms[STHistogramType][nJetsBin].SetBinError(binIndex, 1.)
            if (denominator > 0):
                ratio = numerator/denominator
                ratioError = (1.0/denominator)*math.sqrt(math.pow(numeratorError,2) + math.pow(denominatorError*ratio,2))
                ratioHistograms[STHistogramType][nJetsBin].SetBinContent(binIndex, normalizationFactor*ratio)
                ratioHistograms[STHistogramType][nJetsBin].SetBinError(binIndex, normalizationFactor*ratioError)

tdrstyle.setTDRStyle()

for STHistogramType in STHistogramTypes:
    H_ref = 600
    W_ref = 800
    W = W_ref
    H  = H_ref
    T = 0.08*H_ref
    B = 0.12*H_ref
    L = 0.12*W_ref
    R = 0.04*W_ref

    canvas = ROOT.TCanvas("c_{oFP}_{t}_{n}Jets".format(oFP=inputArguments.outputFilePrefix, t=STHistogramType, n=nJetsBin), "c_{oFP}_{t}_{n}Jets".format(oFP=inputArguments.outputFilePrefix, t=STHistogramType, n=nJetsBin), 50, 50, W, H)
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
    upperPad = ROOT.TPad("upperPad_{t}_{n}Jets".format(t=STHistogramType, n=nJetsBin), "upperPad_{t}_{n}Jets".format(t=STHistogramType, n=nJetsBin), 0., bottomFraction, 0.97, 0.97)
    upperPad.SetMargin(0.12, 0.03, 0.025, 0.08) # left, right, bottom, top
    lowerPad = ROOT.TPad("lowerPad_{t}_{n}Jets".format(t=STHistogramType, n=nJetsBin), "lowerPad_{t}_{n}Jets".format(t=STHistogramType, n=nJetsBin), 0., 0., 0.97, bottomFraction)
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

    STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm].SetLineColor(histColors[inputArguments.nJetsNorm])
    STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm].SetLineWidth(commonLineWidth)
    STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm].GetXaxis().SetTitleSize(commonTitleSize)
    STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm].GetXaxis().SetLabelSize(commonLabelSize)
    STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm].GetXaxis().SetTickLength(0)
    STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm].GetXaxis().SetLabelOffset(999)
    STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm].GetYaxis().SetTitle("Events/GeV")
    STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm].GetYaxis().SetTitleSize(commonTitleSize)
    STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm].GetYaxis().SetLabelSize(commonLabelSize)
    STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm].GetYaxis().SetTitleOffset(commonTitleOffset)

    for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
        if (nJetsBin == inputArguments.nJetsNorm): continue
        STHistogramsScaled[STHistogramType][nJetsBin].SetLineColor(histColors[nJetsBin])
        STHistogramsScaled[STHistogramType][nJetsBin].SetLineWidth(commonLineWidth)

    CMS_lumi.writeExtraText = False
    CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    CMS_lumi.lumi_13TeV = "136.2 fb^{-1}"

    legend = ROOT.TLegend(0.8, 0.475, 0.95, 0.9)
    legend.SetBorderSize(commonLineWidth)
    legend.SetFillStyle(0)
    ROOT.gStyle.SetLegendTextSize(0.05)

    STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm].Draw("P0")
    norm_legendEntry = legend.AddEntry(STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm], "N_{{Jets}} = {n}".format(n=inputArguments.nJetsNorm), "LPE")
    norm_legendEntry.SetTextColor(histColors[inputArguments.nJetsNorm])
    norm_legendEntry.SetMarkerColor(histColors[inputArguments.nJetsNorm])
    STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm].GetXaxis().SetRangeUser(STBoundaries[STHistogramType][0], STBoundaries[STHistogramType][-1])
    STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm].GetYaxis().SetRangeUser(0.0002, 11.)

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
        STHistogramsScaled[STHistogramType][nJetsBin].Draw("AP0 SAME")
        legendText = "N_{{Jets}} = {n}".format(n=nJetsBin)
        if (nJetsBin == inputArguments.nJetsMax): legendText = "N_{{Jets}} #geq {n}".format(n=nJetsBin)
        legendEntry = legend.AddEntry(STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm], legendText, "LPE")
        legendEntry.SetTextColor(histColors[nJetsBin])
        legendEntry.SetMarkerColor(histColors[nJetsBin])
    legend.Draw()
    CMS_lumi.CMS_lumi(canvas, 4, 0)
    upperPad.cd()
    upperPad.Update()
    upperPad.RedrawAxis()
    frame = upperPad.GetFrame()
    frame.Draw()

    yTitleSize_upper = STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm].GetYaxis().GetTitleSize()
    yLabelSize_upper = STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm].GetYaxis().GetLabelSize()
    yTickLength_upper = STHistogramsScaled[STHistogramType][inputArguments.nJetsNorm].GetYaxis().GetTickLength()
    upperPad.Update()

    lowerPad.cd()
    isSet = False
    for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
        if (nJetsBin == inputArguments.nJetsNorm): continue
        ratioHistograms[STHistogramType][nJetsBin].SetLineColor(histColors[nJetsBin])
        ratioHistograms[STHistogramType][nJetsBin].SetLineWidth(commonLineWidth)
        if (isSet):
            ratioHistograms[STHistogramType][nJetsBin].Draw("AP0 SAME")
            continue
        ratioHistograms[STHistogramType][nJetsBin].GetXaxis().SetTitle("S_{T} (GeV)")
        ratioHistograms[STHistogramType][nJetsBin].GetXaxis().SetTitleSize(yTitleSize_upper/bottomToTopRatio)
        ratioHistograms[STHistogramType][nJetsBin].GetXaxis().SetLabelSize(yLabelSize_upper/bottomToTopRatio)
        ratioHistograms[STHistogramType][nJetsBin].GetXaxis().SetTickLength(yTickLength_upper)
        ratioHistograms[STHistogramType][nJetsBin].GetXaxis().SetTitleOffset(0.86)
        ratioHistograms[STHistogramType][nJetsBin].GetYaxis().SetTitle("#frac{N_{Jets} bin}{N_{Jets} = " + str(inputArguments.nJetsNorm) + "}")
        ratioHistograms[STHistogramType][nJetsBin].GetYaxis().SetTitleOffset(1.4*bottomToTopRatio*commonTitleOffset)
        ratioHistograms[STHistogramType][nJetsBin].GetYaxis().SetTitleSize(0.75*yTitleSize_upper/bottomToTopRatio)
        ratioHistograms[STHistogramType][nJetsBin].GetYaxis().SetLabelSize(yLabelSize_upper/bottomToTopRatio)
        ratioHistograms[STHistogramType][nJetsBin].GetYaxis().SetTickLength(yTickLength_upper)
        ratioHistograms[STHistogramType][nJetsBin].GetYaxis().SetNdivisions(2, 0, 0)
        ratioHistograms[STHistogramType][nJetsBin].Draw("P0")
        ratioHistograms[STHistogramType][nJetsBin].GetXaxis().SetRangeUser(STBoundaries[STHistogramType][0], STBoundaries[STHistogramType][-1])
        ratioHistograms[STHistogramType][nJetsBin].GetYaxis().SetRangeUser(-1., 3.)
        isSet = True

    lineAt1 = ROOT.TLine(STBoundaries[STHistogramType][0], 1., STBoundaries[STHistogramType][-1], 1.)
    lineAt1.SetLineColor(histColors[inputArguments.nJetsNorm])
    lineAt1.SetLineWidth(commonLineWidth)
    lineAt1.Draw()
    lowerPad.cd()
    lowerPad.Update()
    lowerPad.RedrawAxis()
    frame = lowerPad.GetFrame()
    frame.Draw()

    canvas.Update()
    canvas.SaveAs("{oD}/{oFP}_{t}.png".format(oD=inputArguments.outputDirectory, oFP=inputArguments.outputFilePrefix, t=STHistogramType))

    print("All done!")
