#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, array, pdb, math
import ROOT, tmROOTUtils, tmGeneralUtils, tdrstyle, CMS_lumi

ROOT.gROOT.SetBatch(ROOT.kTRUE)

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Generate histograms of expected and observed event distributions, based on observed data.')
inputArgumentsParser.add_argument('--path_data_observedNEvents', required=True, help='Path to file containing expected number of events in the format "int expectedNEvents_STRegionX_YJets=Z".',type=str)
inputArgumentsParser.add_argument('--path_data_expectedNEvents', required=True, help='Path to file containing observed number of events in the format "int observedNEvents_STRegionX_YJets=Z".',type=str)
inputArgumentsParser.add_argument('--path_MC_weightedNEvents', required=True, help='Path to ROOT file containing number of events expected from MC samples.',type=str)
inputArgumentsParser.add_argument('--eventProgenitor', required=True, help='Type of stealth sample. Two possible values: \"squark\" or \"gluino\".',type=str)
inputArgumentsParser.add_argument('--path_dataSystematics', default="analysis/dataSystematics/signal_dataSystematics.dat", help='Path to file containing estimated systematic due to norm events fractional uncertainty, shape, and rho.',type=str)
inputArgumentsParser.add_argument('--path_STScalingSystematics', default="analysis/dataSystematics/control_dataSystematics_scaling.dat", help='Path to file containing estimated systematic due to possible deviation from ST scaling.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', required=True, help='Output directory.',type=str)
inputArgumentsParser.add_argument('--outputFilePrefix', required=True, help='Name of output file.',type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--nJetsMin', default=4, help='Min number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Max number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsNorm', default=2, help='Number of jets w.r.t. which to normalize the sT distributions for other jets.',type=int)
inputArgumentsParser.add_argument('--plotObservedData', action='store_true', help="By default, this script does not plot observed data, only the expected number of events with error bars; except if this flag is set.")
inputArguments = inputArgumentsParser.parse_args()

# Gluino, neutralino mass bins to plot
signalBinSettings = {
    2: [],
    3: [],
    # 4: [(1000, 500, ROOT.kRed+1, 21), (1000, 950, ROOT.kBlue+2, 21), (1000, 975, ROOT.kCyan, 21), (1700, 800, ROOT.kMagenta+3, 21), (1700, 1650, ROOT.kYellow, 21), (1700, 1675, ROOT.kGreen+3, 21)],
    # 5: [(1000, 500, ROOT.kRed+1, 21), (1000, 950, ROOT.kBlue+2, 21), (1000, 975, ROOT.kCyan, 21), (1700, 800, ROOT.kMagenta+3, 21), (1700, 1650, ROOT.kYellow, 21), (1700, 1675, ROOT.kGreen+3, 21)],
    # 6: [(1000, 500, ROOT.kRed+1, 21), (1000, 950, ROOT.kBlue+2, 21), (1000, 975, ROOT.kCyan, 21), (1700, 800, ROOT.kMagenta+3, 21), (1700, 1650, ROOT.kYellow, 21), (1700, 1675, ROOT.kGreen+3, 21)]
    # 4: [(1000, 950, ROOT.kBlue+2, 21), (1700, 1675, ROOT.kGreen+3, 21)],
    # 5: [(1000, 500, ROOT.kRed+1, 11), (1000, 950, ROOT.kBlue+2, 21)],
    # 6: [(1000, 500, ROOT.kRed+1, 21), (1000, 950, ROOT.kBlue+2, 21), (1700, 800, ROOT.kMagenta+3, 21)]
    4: [(1150, 600, ROOT.kBlue+2, 21), (1950, 1900, ROOT.kRed+1, 21), (2150, 1000, ROOT.kGreen+3, 21)],
    5: [(1150, 600, ROOT.kBlue+2, 21), (1950, 1900, ROOT.kRed+1, 21), (2150, 1000, ROOT.kGreen+3, 21)],
    6: [(1150, 600, ROOT.kBlue+2, 21), (1950, 1900, ROOT.kRed+1, 21), (2150, 1000, ROOT.kGreen+3, 21)]
}

if not((inputArguments.eventProgenitor == "squark") or (inputArguments.eventProgenitor == "gluino")):
    sys.exit("ERROR: argument \"eventProgenitor\" must be one of \"squark\" or \"gluino\". Current value: {v}".format(v=inputArguments.eventProgenitor))

string_eventProgenitor = None
if (inputArguments.eventProgenitor == "gluino"):
    string_eventProgenitor = "#tilde{g}"
elif (inputArguments.eventProgenitor == "squark"):
    string_eventProgenitor = "#tilde{q}"
string_mass_eventProgenitor = "m_{" + string_eventProgenitor + "}"
string_neutralino = "#tilde{#chi}_{1}^{0}"
string_mass_neutralino = "m_{" + string_neutralino + "}"
string_singlino = "#tilde{S}"
string_mass_singlino = "m_{" + string_singlino + "}"
string_singlet = "S"
string_mass_singlet = "m_{" + string_singlet + "}"
string_gravitino = "#tilde{G}"
string_mass_gravitino = "m_{" + string_gravitino + "}"
string_photon = "#gamma"

def sqrtOfSumOfSquares(listOfNumbers):
    sumOfSquares = 0.
    for number in listOfNumbers:
        sumOfSquares += number*number
    return math.sqrt(sumOfSquares)

def getSignalBinRawText(signalBinSetting):
    meventProgenitor = signalBinSetting[0]
    mneutralino = signalBinSetting[1]
    rawText = ""
    rawText += "("
    rawText += string_mass_eventProgenitor
    rawText += ", "
    rawText += string_mass_neutralino
    rawText += ") = ("
    rawText += str(meventProgenitor)
    rawText += ", "
    rawText += str(mneutralino)
    rawText += ") GeV"
    return rawText

tdrstyle.setTDRStyle()
os.system("mkdir -p {oD}".format(oD=inputArguments.outputDirectory))

STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries)
STBoundaries = []
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()):
        STBoundary = float(STBoundaryString.strip())
        STBoundaries.append(STBoundary)
STBoundaries.append(3500.0) # Instead of infinity
n_STBins = len(STBoundaries) - 1
STRegionsAxis = ROOT.TAxis(n_STBins, array.array('d', STBoundaries))

observedEventCounters_data = tmGeneralUtils.getConfigurationFromFile(inputArguments.path_data_observedNEvents)
expectedEventCounters_data = tmGeneralUtils.getConfigurationFromFile(inputArguments.path_data_expectedNEvents)
dataSystematics = tmGeneralUtils.getConfigurationFromFile(inputArguments.path_dataSystematics)
dataScalingSystematics = tmGeneralUtils.getConfigurationFromFile(inputArguments.path_STScalingSystematics)
signalFile = ROOT.TFile.Open(inputArguments.path_MC_weightedNEvents)
if ((signalFile.IsOpen() == ROOT.kFALSE) or (signalFile.IsZombie())): sys.exit("ERROR: unable to open file with name {n}".format(n=inputArguments.path_MC_signal_weightedNEvents))
expectedNEventsPerGEVHistograms = {} # For "zero error" histograms
expectedNEventsPerGEVHistogramsCopies = {} # Copy with fill color set to white
expectedNEventsPerGEVGraphs = {}
fractionalErrorGraphs = {}
observedNEventsPerGEVGraphs = {}
signalNEventsPerGEVHistograms = {}
signalToDataRatioHistograms = {}
minSignalToExpectedFraction = {}
maxSignalToExpectedFraction = {}
ratioPlots = {}
# whiteColor = ROOT.TColor(9000, 1.0, 1.0, 1.0) # apparently SetFillColor(ROOT.kWhite) does not work (!)
for nJetsBin in range(inputArguments.nJetsMin, 1+inputArguments.nJetsMax):
    expectedNEventsPerGEVHistograms[nJetsBin] = ROOT.TH1F("h_expectedNEvents_{n}Jets".format(n=nJetsBin), "", n_STBins, array.array('d', STBoundaries))
    expectedNEventsPerGEVGraphs[nJetsBin] = ROOT.TGraphAsymmErrors(STRegionsAxis.GetNbins())
    expectedNEventsPerGEVGraphs[nJetsBin].SetName("g_expectedNEventsPerGEVGraphs_{n}Jets".format(n=nJetsBin))
    observedNEventsPerGEVGraphs[nJetsBin] = ROOT.TGraphAsymmErrors(STRegionsAxis.GetNbins())
    observedNEventsPerGEVGraphs[nJetsBin].SetName("g_observedNEventsPerGEVGraphs_{n}Jets".format(n=nJetsBin))
    observedNEventsPerGEVGraphs[nJetsBin].SetName("g_observedNEvents_{n}Jets".format(n=nJetsBin))
    observedNEventsPerGEVGraphs[nJetsBin].SetTitle("")
    fractionalErrorGraphs[nJetsBin] = ROOT.TGraphAsymmErrors(STRegionsAxis.GetNbins())
    fractionalErrorGraphs[nJetsBin].SetName("g_fractionalErrorGraphs_{n}Jets".format(n=nJetsBin))
    signalNEventsPerGEVHistograms[nJetsBin] = {}
    signalToDataRatioHistograms[nJetsBin] = {}
    minSignalToExpectedFraction[nJetsBin] = -1.0
    maxSignalToExpectedFraction[nJetsBin] = -1.0
    for signalBinIndex in range(len(signalBinSettings[nJetsBin])):
        signalNEventsPerGEVHistograms[nJetsBin][signalBinIndex] = ROOT.TH1F("h_signalNEvents_{n}Jets_index{i}".format(n=nJetsBin, i=signalBinIndex), "", n_STBins, array.array('d', STBoundaries))
        signalToDataRatioHistograms[nJetsBin][signalBinIndex] = ROOT.TH1F("h_signalToDataRatio_{n}Jets_index{i}".format(n=nJetsBin, i=signalBinIndex), "", n_STBins, array.array('d', STBoundaries))
    for STRegionIndex in range(1, 1+STRegionsAxis.GetNbins()):
        expectedNEvents = expectedEventCounters_data["expectedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
        expectedNEventsPerGEVHistograms[nJetsBin].SetBinContent(STRegionIndex, expectedNEvents)
        expectedNEventsPerGEVGraphs[nJetsBin].SetPoint(STRegionIndex-1, STRegionsAxis.GetBinCenter(STRegionIndex), expectedNEvents/STRegionsAxis.GetBinWidth(STRegionIndex))
        expectedNEventsPerGEVGraphs[nJetsBin].SetPointEXlow(STRegionIndex-1, 0.5*STRegionsAxis.GetBinWidth(STRegionIndex))
        expectedNEventsPerGEVGraphs[nJetsBin].SetPointEXhigh(STRegionIndex-1, 0.5*STRegionsAxis.GetBinWidth(STRegionIndex))
        expectedNEventsPerGEVGraphs[nJetsBin].SetPointEYlow(STRegionIndex-1, 0.)
        expectedNEventsPerGEVGraphs[nJetsBin].SetPointEYhigh(STRegionIndex-1, 0.)
        expectedNEventsPerGEVHistograms[nJetsBin].SetBinError(STRegionIndex, 0.)
        fractionalErrorGraphs[nJetsBin].SetPoint(STRegionIndex-1, STRegionsAxis.GetBinCenter(STRegionIndex), 1.)
        fractionalErrorGraphs[nJetsBin].SetPointEXlow(STRegionIndex-1, 0.5*STRegionsAxis.GetBinWidth(STRegionIndex))
        fractionalErrorGraphs[nJetsBin].SetPointEXhigh(STRegionIndex-1, 0.5*STRegionsAxis.GetBinWidth(STRegionIndex))
        fractionalErrorGraphs[nJetsBin].SetPointEYlow(STRegionIndex-1, 0.)
        fractionalErrorGraphs[nJetsBin].SetPointEYhigh(STRegionIndex-1, 0.)
        if (STRegionIndex > 1):
            expectedNEventsErrorDown_normEvents = dataSystematics["fractionalUncertaintyDown_normEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
            expectedNEventsErrorUp_normEvents = dataSystematics["fractionalUncertaintyUp_normEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
            expectedNEventsError_shape = dataSystematics["fractionalUncertainty_shape_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
            expectedNEventsError_rho = dataSystematics["fractionalUncertainty_rho_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
            expectedNEventsError_scaling = 0.
            if (nJetsBin != inputArguments.nJetsNorm): expectedNEventsError_scaling = dataScalingSystematics["fractionalUncertainty_residual_scaling_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
            expectedNEvents_netFractionalErrorDown = sqrtOfSumOfSquares([expectedNEventsErrorDown_normEvents, expectedNEventsError_shape, expectedNEventsError_rho, expectedNEventsError_scaling])
            expectedNEvents_netFractionalErrorUp = sqrtOfSumOfSquares([expectedNEventsErrorUp_normEvents, expectedNEventsError_shape, expectedNEventsError_rho, expectedNEventsError_scaling])
            expectedNEventsPerGEVGraphs[nJetsBin].SetPointEYlow(STRegionIndex-1, expectedNEvents_netFractionalErrorDown*expectedNEvents/STRegionsAxis.GetBinWidth(STRegionIndex))
            expectedNEventsPerGEVGraphs[nJetsBin].SetPointEYhigh(STRegionIndex-1, expectedNEvents_netFractionalErrorUp*expectedNEvents/STRegionsAxis.GetBinWidth(STRegionIndex))
            fractionalErrorGraphs[nJetsBin].SetPointEYlow(STRegionIndex-1, expectedNEvents_netFractionalErrorDown)
            fractionalErrorGraphs[nJetsBin].SetPointEYhigh(STRegionIndex-1, expectedNEvents_netFractionalErrorUp)
        signalNEventsHistogramSource = ROOT.TH2F()
        signalFile.GetObject("h_lumiBasedYearWeightedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin), signalNEventsHistogramSource)
        for signalBinIndex in range(len(signalBinSettings[nJetsBin])):
            signalBinSetting = signalBinSettings[nJetsBin][signalBinIndex]
            signalNEvents = signalNEventsHistogramSource.GetBinContent(signalNEventsHistogramSource.FindFixBin(signalBinSetting[0], signalBinSetting[1]))
            signalNEventsPerGEVHistograms[nJetsBin][signalBinIndex].SetBinContent(STRegionIndex, signalNEvents)
            signalNEventsPerGEVHistograms[nJetsBin][signalBinIndex].SetBinError(STRegionIndex, 0.)
            signalToDataRatioHistograms[nJetsBin][signalBinIndex].SetBinContent(STRegionIndex, signalNEvents/expectedNEvents)
            signalToDataRatioHistograms[nJetsBin][signalBinIndex].SetBinError(STRegionIndex, 0.)
            if ((minSignalToExpectedFraction[nJetsBin] < 0.0) or (signalNEvents/expectedNEvents < minSignalToExpectedFraction[nJetsBin])):
                minSignalToExpectedFraction[nJetsBin] = signalNEvents/expectedNEvents
            if ((maxSignalToExpectedFraction[nJetsBin] < 0.0) or (signalNEvents/expectedNEvents > maxSignalToExpectedFraction[nJetsBin])):
                maxSignalToExpectedFraction[nJetsBin] = signalNEvents/expectedNEvents
        if not(inputArguments.plotObservedData): continue
        observedNEvents = observedEventCounters_data["observedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
        observedNEventsPerGEVGraphs[nJetsBin].SetPoint(STRegionIndex-1, STRegionsAxis.GetBinCenter(STRegionIndex), observedNEvents/STRegionsAxis.GetBinWidth(STRegionIndex))
        # observedNEventsPerGEVGraphs[nJetsBin].SetPointEXlow(STRegionIndex-1, 0.5*STRegionsAxis.GetBinWidth(STRegionIndex))
        # observedNEventsPerGEVGraphs[nJetsBin].SetPointEXhigh(STRegionIndex-1, 0.5*STRegionsAxis.GetBinWidth(STRegionIndex))
        observedNEventsPerGEVGraphs[nJetsBin].SetPointEXlow(STRegionIndex-1, 0.)
        observedNEventsPerGEVGraphs[nJetsBin].SetPointEXhigh(STRegionIndex-1, 0.)
        poissonInterval = tmROOTUtils.getPoissonConfidenceInterval(observedNEvents=observedNEvents)
        observedNEventsPerGEVGraphs[nJetsBin].SetPointEYlow(STRegionIndex-1, (observedNEvents-poissonInterval["lower"])/STRegionsAxis.GetBinWidth(STRegionIndex))
        observedNEventsPerGEVGraphs[nJetsBin].SetPointEYhigh(STRegionIndex-1, (poissonInterval["upper"]-observedNEvents)/STRegionsAxis.GetBinWidth(STRegionIndex))

    tmROOTUtils.rescale1DHistogramByBinWidth(expectedNEventsPerGEVHistograms[nJetsBin])
    for signalBinIndex in range(len(signalBinSettings[nJetsBin])):
        tmROOTUtils.rescale1DHistogramByBinWidth(signalNEventsPerGEVHistograms[nJetsBin][signalBinIndex])

    H_ref = 600
    W_ref = 800
    W = W_ref
    H  = H_ref
    T = 0.08*H_ref
    B = 0.12*H_ref
    L = 0.12*W_ref
    R = 0.04*W_ref

    canvas = ROOT.TCanvas("c_{oFN}_{n}Jets".format(oFN=inputArguments.outputFilePrefix, n=nJetsBin), "c_{oFN}_{n}Jets".format(oFN=inputArguments.outputFilePrefix, n=nJetsBin), 50, 50, W, H)
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
    commonFillColor = ROOT.kOrange-2
    commonExpectedEventsLineColor = ROOT.kBlack
    commonExpectedEventsLineStyle = 2
    commonExpectedEventsLineWidth = 3

    upperPad.cd()
    upperPad.SetLogy()

    expectedNEventsPerGEVHistograms[nJetsBin].SetLineColor(commonExpectedEventsLineColor)
    expectedNEventsPerGEVHistograms[nJetsBin].SetLineStyle(commonExpectedEventsLineStyle)
    expectedNEventsPerGEVHistograms[nJetsBin].SetLineWidth(commonExpectedEventsLineWidth)
    expectedNEventsPerGEVHistograms[nJetsBin].SetFillColor(commonFillColor)
    expectedNEventsPerGEVHistograms[nJetsBin].SetMarkerSize(0)
    expectedNEventsPerGEVHistograms[nJetsBin].GetXaxis().SetLabelSize(0)
    expectedNEventsPerGEVHistograms[nJetsBin].GetXaxis().SetTickLength(0)
    expectedNEventsPerGEVHistograms[nJetsBin].GetXaxis().SetLabelOffset(999)
    expectedNEventsPerGEVHistograms[nJetsBin].GetYaxis().SetTitle("Events/GeV")
    expectedNEventsPerGEVHistograms[nJetsBin].GetYaxis().SetTitleOffset(commonTitleOffset)

    expectedNEventsPerGEVHistogramsCopies[nJetsBin] = expectedNEventsPerGEVHistograms[nJetsBin].Clone() # Create clone
    expectedNEventsPerGEVHistogramsCopies[nJetsBin].SetName("h_copy_expectedNEvents_{n}Jets".format(n=nJetsBin))
    expectedNEventsPerGEVHistogramsCopies[nJetsBin].SetFillColor(ROOT.kWhite)

    expectedNEventsPerGEVGraphs[nJetsBin].SetFillColor(commonFillColor)

    observedNEventsPerGEVGraphs[nJetsBin].SetLineColor(ROOT.kBlack)
    observedNEventsPerGEVGraphs[nJetsBin].SetFillColor(ROOT.kWhite)

    for signalBinIndex in range(len(signalBinSettings[nJetsBin])):
        signalBinSetting = signalBinSettings[nJetsBin][signalBinIndex]
        signalNEventsPerGEVHistograms[nJetsBin][signalBinIndex].SetLineColor(signalBinSetting[2])
        signalNEventsPerGEVHistograms[nJetsBin][signalBinIndex].SetLineStyle(5)
        signalNEventsPerGEVHistograms[nJetsBin][signalBinIndex].SetLineWidth(2)
        signalToDataRatioHistograms[nJetsBin][signalBinIndex].SetLineColor(signalBinSetting[2])
        signalToDataRatioHistograms[nJetsBin][signalBinIndex].SetLineStyle(5)
        signalToDataRatioHistograms[nJetsBin][signalBinIndex].SetLineWidth(2)

    CMS_lumi.writeExtraText = False
    CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    CMS_lumi.lumi_13TeV = "136.2 fb^{-1}"

    legend = ROOT.TLegend(0.3, 0.85, 0.95, 0.9)
    legend.SetNColumns(3)
    legend.AddEntry(None, "N_{{Jets}} = {n}".format(n=nJetsBin), "")
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    ROOT.gStyle.SetLegendTextSize(0.05)

    expectedNEventsPerGEVHistograms[nJetsBin].Draw("][") # First draw filled so that the legend entry is appropriate
    legend.AddEntry(expectedNEventsPerGEVHistograms[nJetsBin], "Predicted Background")
    expectedNEventsPerGEVHistogramsCopies[nJetsBin].Draw("][") # Next draw with white filling, overwriting previous histogram
    expectedNEventsPerGEVHistogramsCopies[nJetsBin].GetXaxis().SetRangeUser(STBoundaries[0], STBoundaries[-1])
    expectedNEventsPerGEVHistogramsCopies[nJetsBin].GetYaxis().SetRangeUser(0.0002, 11.)
    expectedNEventsPerGEVGraphs[nJetsBin].Draw("2") # For the yellow bands
    for signalBinIndex in range(len(signalBinSettings[nJetsBin])):
        signalBinSetting = signalBinSettings[nJetsBin][signalBinIndex]
        signalNEventsPerGEVHistograms[nJetsBin][signalBinIndex].Draw("A HIST SAME") # Signal distributions
        maxNSignalEvents_xpos = signalNEventsPerGEVHistograms[nJetsBin][signalBinIndex].GetBinCenter(signalNEventsPerGEVHistograms[nJetsBin][signalBinIndex].GetMaximumBin())
        maxNSignalEvents_ypos = signalNEventsPerGEVHistograms[nJetsBin][signalBinIndex].GetBinContent(signalNEventsPerGEVHistograms[nJetsBin][signalBinIndex].GetMaximumBin())
        if (signalBinSetting[3] == 11): maxNSignalEvents_xpos += (-0.4)*signalNEventsPerGEVHistograms[nJetsBin][signalBinIndex].GetBinWidth(signalNEventsPerGEVHistograms[nJetsBin][signalBinIndex].GetMaximumBin()) # For left-aligned labels
        latex = ROOT.TLatex()
        latex.SetTextFont(42)
        latex.SetTextAngle(0)
        latex.SetTextColor(signalBinSetting[2])
        latex.SetTextSize(0.045)
        latex.SetTextAlign(signalBinSetting[3])
        shiftFactor = 1.6
        # if (maxNSignalEvents_ypos < expectedNEventsPerGEVHistogramsCopies[nJetsBin].GetBinContent(signalNEventsPerGEVHistograms[nJetsBin][signalBinIndex].GetMaximumBin())): shiftFactor = 1.0/1.4
        latex.DrawLatex(maxNSignalEvents_xpos, shiftFactor*maxNSignalEvents_ypos, getSignalBinRawText(signalBinSetting))

    if (inputArguments.plotObservedData):
        observedNEventsPerGEVGraphs[nJetsBin].Draw("0PZ")
        legend.AddEntry(observedNEventsPerGEVGraphs[nJetsBin], "Data", "PE")
    expectedNEventsPerGEVHistogramsCopies[nJetsBin].Draw("][ SAME") # Have to draw again to get overlay on top of previous histograms
    legend.Draw()
    CMS_lumi.CMS_lumi(canvas, 4, 0)
    upperPad.cd()
    upperPad.Update()
    upperPad.RedrawAxis()
    frame = upperPad.GetFrame()
    frame.Draw()

    yTitleSize_upper = expectedNEventsPerGEVHistogramsCopies[nJetsBin].GetYaxis().GetTitleSize()
    yLabelSize_upper = expectedNEventsPerGEVHistogramsCopies[nJetsBin].GetYaxis().GetLabelSize()
    yTickLength_upper = expectedNEventsPerGEVHistogramsCopies[nJetsBin].GetYaxis().GetTickLength()
    upperPad.Update()

    lowerPad.cd()
    fractionalErrorGraphs[nJetsBin].GetXaxis().SetTitle("S_{T} (GeV)")
    fractionalErrorGraphs[nJetsBin].GetXaxis().SetTitleSize(yTitleSize_upper/bottomToTopRatio)
    fractionalErrorGraphs[nJetsBin].GetXaxis().SetLabelSize(yLabelSize_upper/bottomToTopRatio)
    fractionalErrorGraphs[nJetsBin].GetXaxis().SetTickLength(yTickLength_upper)
    fractionalErrorGraphs[nJetsBin].GetYaxis().SetTitle("#frac{Data}{Background}")
    fractionalErrorGraphs[nJetsBin].GetYaxis().SetTitleOffset(1.4*bottomToTopRatio*commonTitleOffset)
    fractionalErrorGraphs[nJetsBin].GetYaxis().SetTitleSize(0.75*yTitleSize_upper/bottomToTopRatio)
    fractionalErrorGraphs[nJetsBin].GetYaxis().SetLabelSize(yLabelSize_upper/bottomToTopRatio)
    fractionalErrorGraphs[nJetsBin].GetYaxis().SetTickLength(yTickLength_upper)
    fractionalErrorGraphs[nJetsBin].GetYaxis().SetNdivisions(2, 0, 0)
    fractionalErrorGraphs[nJetsBin].SetFillColor(commonFillColor)
    fractionalErrorGraphs[nJetsBin].Draw("A2")

    if (inputArguments.plotObservedData):
        ratioPlots[nJetsBin] = tmROOTUtils.getGraphOfRatioOfAsymmErrorsGraphToHistogram(numeratorGraph=observedNEventsPerGEVGraphs[nJetsBin], denominatorHistogram=expectedNEventsPerGEVHistograms[nJetsBin], outputName="g_ratioGraphs_{n}Jets".format(n=nJetsBin), outputTitle="")
        ratioPlots[nJetsBin].Draw("0PZ")
    else: # Only draw signal histogram ratios if observed data is not to be plotted; otherwise plot looks too cluttered
        for signalBinIndex in range(len(signalBinSettings[nJetsBin])):
            signalToDataRatioHistograms[nJetsBin][signalBinIndex].Draw("HIST SAME")
    lineAt1 = ROOT.TLine(STBoundaries[0], 1., STBoundaries[-1], 1.)
    lineAt1.SetLineColor(commonExpectedEventsLineColor)
    lineAt1.SetLineStyle(commonExpectedEventsLineStyle)
    lineAt1.SetLineWidth(commonExpectedEventsLineWidth)
    lineAt1.Draw()
    fractionalErrorGraphs[nJetsBin].GetXaxis().SetRangeUser(STBoundaries[0], STBoundaries[-1])
    if (inputArguments.plotObservedData):
        fractionalErrorGraphs[nJetsBin].GetYaxis().SetRangeUser(-1., 3.)
    else:
        fractionalErrorGraphs[nJetsBin].GetYaxis().SetRangeUser(0.95*minSignalToExpectedFraction[nJetsBin], 1.05*maxSignalToExpectedFraction[nJetsBin])
    lowerPad.cd()
    lowerPad.Update()
    lowerPad.RedrawAxis()
    frame = lowerPad.GetFrame()
    frame.Draw()

    canvas.Update()
    canvas.SaveAs("{oD}/{oFN}_{n}Jets.png".format(oD=inputArguments.outputDirectory, oFN=inputArguments.outputFilePrefix, n=nJetsBin))

signalFile.Close()
print("All done!")
