#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, array, pdb, math
import ROOT, tmROOTUtils, tmGeneralUtils, tdrstyle, CMS_lumi

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Generate histograms of expected and observed event distributions, based on observed data.')
inputArgumentsParser.add_argument('--path_data_observedNEvents', required=True, help='Path to file containing expected number of events in the format "int observedNEvents_STRegionX_YJets=Z".',type=str)
inputArgumentsParser.add_argument('--path_data_expectedNEvents', required=True, help='Path to file containing observed number of events in the format "float expectedNEvents_STRegionX_YJets=Z".',type=str)
inputArgumentsParser.add_argument('--path_data_adjustments', required=True, help='Path to file containing adjustments derived from MC in the format "float nominalAdjustment_STRegionX_YJets=Z".',type=str)
inputArgumentsParser.add_argument('--path_MC_weightedNEvents_gluino', required=True, help='Path to ROOT file containing number of events expected from gluino MC samples.',type=str)
inputArgumentsParser.add_argument('--path_MC_weightedNEvents_squark', required=True, help='Path to ROOT file containing number of events expected from squark MC samples.',type=str)
inputArgumentsParser.add_argument('--path_fitDiagnostics', default="/dev/null", help='Path to ROOT file containing the fit diagnostics output.',type=str)
inputArgumentsParser.add_argument('--bin_label_abbreviation', default="none", help='Bin label abbreviation to use while reading fit diagnostics.',type=str)
# inputArgumentsParser.add_argument('--eventProgenitor', required=True, help='Type of stealth sample. Two possible values: \"squark\" or \"gluino\".',type=str)
inputArgumentsParser.add_argument('--path_systematics_nominal', required=True, help='Path to file containing systematics due to norm events fractional uncertainty, shape, and rho.',type=str)
inputArgumentsParser.add_argument('--path_systematics_dataMCDiscrepancy', required=True, help='Path to file containing estimated systematics on residual MC-data discrepancy.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', required=True, help='Output directory.',type=str)
inputArgumentsParser.add_argument('--outputFilePrefix', required=True, help='Name of output file.',type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--nJetsBin', required=True, help='nJets bin for plotting.',type=int)
inputArgumentsParser.add_argument('--plotObservedData', action='store_true', help="By default, this script does not plot observed data, only the expected number of events with error bars; except if this flag is set.")
inputArgumentsParser.add_argument('--suppressSignal', action='store_true', help="By default, this script plots the signal predictions from some chosen signal bins on top of the background predictions; except if this flag is set.")
inputArgumentsParser.add_argument('--ratioMin', default=0.0, help='Max value of ratio to plot.',type=float)
inputArgumentsParser.add_argument('--ratioMax', default=2.5, help='Max value of ratio to plot.',type=float)
inputArguments = inputArgumentsParser.parse_args()

plot_signal = not(inputArguments.suppressSignal)
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
    4: [("gluino", 1100, 200, ROOT.kBlue+2, 21), ("gluino", 2000, 1900, ROOT.kRed+1, 21), ("gluino", 2000, 1000, ROOT.kGreen+3, 21)],
    5: [("gluino", 1100, 200, ROOT.kBlue+2, 21), ("gluino", 2000, 1900, ROOT.kRed+1, 21), ("gluino", 2000, 1000, ROOT.kGreen+3, 21)],
    6: [("gluino", 1100, 200, ROOT.kBlue+2, 21), ("gluino", 2000, 1900, ROOT.kRed+1, 21), ("gluino", 2000, 1000, ROOT.kGreen+3, 21)]
}
inputMCWeightedNEventsFilePaths = {
    "gluino": inputArguments.path_MC_weightedNEvents_gluino,
    "squark": inputArguments.path_MC_weightedNEvents_squark
}

# if not((inputArguments.eventProgenitor == "squark") or (inputArguments.eventProgenitor == "gluino")):
#     sys.exit("ERROR: argument \"eventProgenitor\" must be one of \"squark\" or \"gluino\". Current value: {v}".format(v=inputArguments.eventProgenitor))

def get_string_event_progenitor(event_progenitor):
    string_eventProgenitor = None
    if (event_progenitor == "gluino"):
        string_eventProgenitor = "#tilde{g}"
    elif (event_progenitor == "squark"):
        string_eventProgenitor = "#tilde{q}"
    else:
        sys.exit("ERROR: unidentified event_progenitor {e}".format(e=event_progenitor))
    return string_eventProgenitor

def get_string_mass_event_progenitor(event_progenitor):
    return ("m_{" + get_string_event_progenitor(event_progenitor) + "}")

string_neutralino = "#tilde{#chi}_{1}^{0}"
string_mass_neutralino = "m_{" + string_neutralino + "}"
string_singlino = "#tilde{S}"
string_mass_singlino = "m_{" + string_singlino + "}"
string_singlet = "S"
string_mass_singlet = "m_{" + string_singlet + "}"
string_gravitino = "#tilde{G}"
string_mass_gravitino = "m_{" + string_gravitino + "}"
string_photon = "#gamma"
nJetsBin = inputArguments.nJetsBin

def sqrtOfSumOfSquares(listOfNumbers):
    sumOfSquares = 0.
    for number in listOfNumbers:
        sumOfSquares += number*number
    return math.sqrt(sumOfSquares)

def getSignalBinRawText(signalBinSetting):
    event_progenitor = signalBinSetting[0]
    meventProgenitor = signalBinSetting[1]
    mneutralino = signalBinSetting[2]
    rawText = ""
    rawText += "("
    rawText += get_string_mass_event_progenitor(event_progenitor)
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
adjustments_data = tmGeneralUtils.getConfigurationFromFile(inputArguments.path_data_adjustments)
systematics_nominal = tmGeneralUtils.getConfigurationFromFile(inputArguments.path_systematics_nominal)
systematics_dataMCDiscrepancy = tmGeneralUtils.getConfigurationFromFile(inputArguments.path_systematics_dataMCDiscrepancy)
signalFiles = None
if plot_signal:
    signalFiles = {}
    for eventProgenitor in ["gluino", "squark"]:
        signalFiles[eventProgenitor] = ROOT.TFile.Open(inputMCWeightedNEventsFilePaths[eventProgenitor])
        if ((signalFiles[eventProgenitor].IsOpen() == ROOT.kFALSE) or (signalFiles[eventProgenitor].IsZombie())): sys.exit("ERROR: unable to open file with name {n}".format(n=inputMCWeightedNEventsFilePaths[eventProgenitor]))

fitDiagnosticsFile = None
if inputArguments.plotObservedData:
    if not(inputArguments.path_fitDiagnostics == "/dev/null"):
        fitDiagnosticsFile = ROOT.TFile.Open(inputArguments.path_fitDiagnostics)

def get_bin_label_abbreviated(STRegionIndex, nJetsBin):
    return ("{l}ST{i}J{j}".format(l=inputArguments.bin_label_abbreviation, i=STRegionIndex, j=nJetsBin))

def get_pre_fit_background(STRegionIndex, nJetsBin):
    input_histogram = ROOT.TH1F()
    fitDiagnosticsFile.GetObject("shapes_prefit/{l}/qcd".format(l=get_bin_label_abbreviated(STRegionIndex, nJetsBin)), input_histogram)
    if (not(input_histogram.GetXaxis().GetNbins() == 1)): sys.exit("ERROR: histogram shapes_prefit/{l}/qcd is in an unexpected format.".format(l=get_bin_label_abbreviated(STRegionIndex, nJetsBin)))
    return tuple([input_histogram.GetBinContent(1), input_histogram.GetBinErrorLow(1), input_histogram.GetBinErrorUp(1)])

def get_post_fit_background(STRegionIndex, nJetsBin):
    input_histogram = ROOT.TH1F()
    fitDiagnosticsFile.GetObject("shapes_fit_b/{l}/qcd".format(l=get_bin_label_abbreviated(STRegionIndex, nJetsBin)), input_histogram)
    if (not(input_histogram.GetXaxis().GetNbins() == 1)): sys.exit("ERROR: histogram shapes_fit_b/{l}/qcd is in an unexpected format.".format(l=get_bin_label_abbreviated(STRegionIndex, nJetsBin)))
    return tuple([input_histogram.GetBinContent(1), input_histogram.GetBinErrorLow(1), input_histogram.GetBinErrorUp(1)])

# whiteColor = ROOT.TColor(9000, 1.0, 1.0, 1.0) # apparently SetFillColor(ROOT.kWhite) does not work (!)
expectedNEventsPerGEVHistogram = ROOT.TH1F("h_expectedNEvents_{n}Jets".format(n=nJetsBin), "", n_STBins, array.array('d', STBoundaries))
expectedNEventsPerGEVGraph = ROOT.TGraphAsymmErrors(STRegionsAxis.GetNbins())
expectedNEventsPerGEVGraph.SetName("g_expectedNEventsPerGEVGraphs_{n}Jets".format(n=nJetsBin))
observedNEventsPerGEVGraph = ROOT.TGraphAsymmErrors(STRegionsAxis.GetNbins())
observedNEventsPerGEVGraph.SetName("g_observedNEventsPerGEVGraphs_{n}Jets".format(n=nJetsBin))
observedNEventsPerGEVGraph.SetName("g_observedNEvents_{n}Jets".format(n=nJetsBin))
observedNEventsPerGEVGraph.SetTitle("")
fractionalErrorGraph = ROOT.TGraphAsymmErrors(STRegionsAxis.GetNbins())
fractionalErrorGraph.SetName("g_fractionalErrorGraphs_{n}Jets".format(n=nJetsBin))
signalNEventsPerGEVHistograms = None
signalToDataRatioHistograms = None
minSignalToExpectedFraction = None
maxSignalToExpectedFraction = None
if plot_signal:
    signalNEventsPerGEVHistograms = {}
    signalToDataRatioHistograms = {}
    minSignalToExpectedFraction = -1.0
    maxSignalToExpectedFraction = -1.0
    for signalBinIndex in range(len(signalBinSettings[nJetsBin])):
        signalNEventsPerGEVHistograms[signalBinIndex] = ROOT.TH1F("h_signalNEvents_{n}Jets_index{i}".format(n=nJetsBin, i=signalBinIndex), "", n_STBins, array.array('d', STBoundaries))
        signalToDataRatioHistograms[signalBinIndex] = ROOT.TH1F("h_signalToDataRatio_{n}Jets_index{i}".format(n=nJetsBin, i=signalBinIndex), "", n_STBins, array.array('d', STBoundaries))
for STRegionIndex in range(1, 1+STRegionsAxis.GetNbins()):
    expectedNEvents = (expectedEventCounters_data["expectedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)])
    expectedNEventsErrorFromFitDown = None
    expectedNEventsErrorFromFitUp = None
    if (STRegionIndex > 1):
        if ((not(inputArguments.plotObservedData)) or (fitDiagnosticsFile is None)):
            expectedNEvents *= (adjustments_data["nominalAdjustment_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)])
        else:
            expectedNEvents, expectedNEventsErrorFromFitDown, expectedNEventsErrorFromFitUp = get_post_fit_background(STRegionIndex, nJetsBin)
    expectedNEventsPerGEVHistogram.SetBinContent(STRegionIndex, expectedNEvents)
    expectedNEventsPerGEVGraph.SetPoint(STRegionIndex-1, STRegionsAxis.GetBinCenter(STRegionIndex), expectedNEvents/STRegionsAxis.GetBinWidth(STRegionIndex))
    expectedNEventsPerGEVGraph.SetPointEXlow(STRegionIndex-1, 0.5*STRegionsAxis.GetBinWidth(STRegionIndex))
    expectedNEventsPerGEVGraph.SetPointEXhigh(STRegionIndex-1, 0.5*STRegionsAxis.GetBinWidth(STRegionIndex))
    expectedNEventsPerGEVGraph.SetPointEYlow(STRegionIndex-1, 0.)
    expectedNEventsPerGEVGraph.SetPointEYhigh(STRegionIndex-1, 0.)
    expectedNEventsPerGEVHistogram.SetBinError(STRegionIndex, 0.)
    fractionalErrorGraph.SetPoint(STRegionIndex-1, STRegionsAxis.GetBinCenter(STRegionIndex), 1.)
    fractionalErrorGraph.SetPointEXlow(STRegionIndex-1, 0.5*STRegionsAxis.GetBinWidth(STRegionIndex))
    fractionalErrorGraph.SetPointEXhigh(STRegionIndex-1, 0.5*STRegionsAxis.GetBinWidth(STRegionIndex))
    fractionalErrorGraph.SetPointEYlow(STRegionIndex-1, 0.)
    fractionalErrorGraph.SetPointEYhigh(STRegionIndex-1, 0.)
    if (STRegionIndex > 1):
        expectedNEventsErrorDown_normEvents = systematics_nominal["fractionalUncertaintyDown_normEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
        expectedNEventsErrorUp_normEvents = systematics_nominal["fractionalUncertaintyUp_normEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
        expectedNEventsError_shape = systematics_nominal["fractionalUncertainty_shape_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
        expectedNEventsError_rho = systematics_nominal["fractionalUncertainty_rho_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
        expectedNEventsErrorDown_adjustment_mode0 = adjustments_data["fractionalUncertaintyDown_mode0_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
        expectedNEventsErrorUp_adjustment_mode0 = adjustments_data["fractionalUncertaintyUp_mode0_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
        expectedNEventsErrorDown_adjustment_mode1 = adjustments_data["fractionalUncertaintyDown_mode1_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
        expectedNEventsErrorUp_adjustment_mode1 = adjustments_data["fractionalUncertaintyUp_mode1_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
        expectedNEventsError_DataMCDiscrepancy = 2.0*(abs((systematics_dataMCDiscrepancy["ratio_adjustment_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]) - 1.0))
        expectedNEvents_netFractionalErrorDown = sqrtOfSumOfSquares([expectedNEventsErrorDown_normEvents, expectedNEventsError_shape, expectedNEventsError_rho, expectedNEventsErrorDown_adjustment_mode0, expectedNEventsErrorDown_adjustment_mode1, expectedNEventsError_DataMCDiscrepancy])
        expectedNEvents_netFractionalErrorUp = sqrtOfSumOfSquares([expectedNEventsErrorUp_normEvents, expectedNEventsError_shape, expectedNEventsError_rho, expectedNEventsErrorUp_adjustment_mode0, expectedNEventsErrorUp_adjustment_mode1, expectedNEventsError_DataMCDiscrepancy])
        if ((expectedNEventsErrorFromFitDown is None) or (expectedNEventsErrorFromFitUp is None)):
            expectedNEventsPerGEVGraph.SetPointEYlow(STRegionIndex-1, expectedNEvents_netFractionalErrorDown*expectedNEvents/STRegionsAxis.GetBinWidth(STRegionIndex))
            expectedNEventsPerGEVGraph.SetPointEYhigh(STRegionIndex-1, expectedNEvents_netFractionalErrorUp*expectedNEvents/STRegionsAxis.GetBinWidth(STRegionIndex))
            fractionalErrorGraph.SetPointEYlow(STRegionIndex-1, expectedNEvents_netFractionalErrorDown)
            fractionalErrorGraph.SetPointEYhigh(STRegionIndex-1, expectedNEvents_netFractionalErrorUp)
        else:
            expectedNEventsPerGEVGraph.SetPointEYlow(STRegionIndex-1, expectedNEventsErrorFromFitDown/STRegionsAxis.GetBinWidth(STRegionIndex))
            expectedNEventsPerGEVGraph.SetPointEYhigh(STRegionIndex-1, expectedNEventsErrorFromFitUp/STRegionsAxis.GetBinWidth(STRegionIndex))
            fractionalErrorGraph.SetPointEYlow(STRegionIndex-1, expectedNEventsErrorFromFitDown/expectedNEvents)
            fractionalErrorGraph.SetPointEYhigh(STRegionIndex-1, expectedNEventsErrorFromFitUp/expectedNEvents)
    if plot_signal:
        signalNEventsHistogramSources = {}
        for eventProgenitor in ["gluino", "squark"]:
            signalNEventsHistogramSources[eventProgenitor] = ROOT.TH2F()
            signalFiles[eventProgenitor].GetObject("h_lumiBasedYearWeightedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin), signalNEventsHistogramSources[eventProgenitor])
        for signalBinIndex in range(len(signalBinSettings[nJetsBin])):
            signalBinSetting = signalBinSettings[nJetsBin][signalBinIndex]
            eventProgenitor = signalBinSetting[0]
            signalNEvents = signalNEventsHistogramSources[eventProgenitor].GetBinContent(signalNEventsHistogramSources[eventProgenitor].FindFixBin(signalBinSetting[1], signalBinSetting[2]))
            signalNEventsPerGEVHistograms[signalBinIndex].SetBinContent(STRegionIndex, signalNEvents)
            signalNEventsPerGEVHistograms[signalBinIndex].SetBinError(STRegionIndex, 0.)
            signalToDataRatioHistograms[signalBinIndex].SetBinContent(STRegionIndex, signalNEvents/expectedNEvents)
            signalToDataRatioHistograms[signalBinIndex].SetBinError(STRegionIndex, 0.)
            if ((minSignalToExpectedFraction < 0.0) or (signalNEvents/expectedNEvents < minSignalToExpectedFraction)):
                minSignalToExpectedFraction = signalNEvents/expectedNEvents
            if ((maxSignalToExpectedFraction < 0.0) or (signalNEvents/expectedNEvents > maxSignalToExpectedFraction)):
                maxSignalToExpectedFraction = signalNEvents/expectedNEvents
    if not(inputArguments.plotObservedData): continue
    observedNEvents = observedEventCounters_data["observedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
    observedNEventsPerGEVGraph.SetPoint(STRegionIndex-1, STRegionsAxis.GetBinCenter(STRegionIndex), observedNEvents/STRegionsAxis.GetBinWidth(STRegionIndex))
    # observedNEventsPerGEVGraph.SetPointEXlow(STRegionIndex-1, 0.5*STRegionsAxis.GetBinWidth(STRegionIndex))
    # observedNEventsPerGEVGraph.SetPointEXhigh(STRegionIndex-1, 0.5*STRegionsAxis.GetBinWidth(STRegionIndex))
    observedNEventsPerGEVGraph.SetPointEXlow(STRegionIndex-1, 0.)
    observedNEventsPerGEVGraph.SetPointEXhigh(STRegionIndex-1, 0.)
    poissonInterval = tmROOTUtils.getPoissonConfidenceInterval(observedNEvents=observedNEvents)
    observedNEventsPerGEVGraph.SetPointEYlow(STRegionIndex-1, (observedNEvents-poissonInterval["lower"])/STRegionsAxis.GetBinWidth(STRegionIndex))
    observedNEventsPerGEVGraph.SetPointEYhigh(STRegionIndex-1, (poissonInterval["upper"]-observedNEvents)/STRegionsAxis.GetBinWidth(STRegionIndex))

tmROOTUtils.rescale1DHistogramByBinWidth(expectedNEventsPerGEVHistogram)
if plot_signal:
    for signalBinIndex in range(len(signalBinSettings[nJetsBin])):
        tmROOTUtils.rescale1DHistogramByBinWidth(signalNEventsPerGEVHistograms[signalBinIndex])

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

expectedNEventsPerGEVHistogram.SetLineColor(commonExpectedEventsLineColor)
expectedNEventsPerGEVHistogram.SetLineStyle(commonExpectedEventsLineStyle)
expectedNEventsPerGEVHistogram.SetLineWidth(commonExpectedEventsLineWidth)
expectedNEventsPerGEVHistogram.SetFillColor(commonFillColor)
expectedNEventsPerGEVHistogram.SetMarkerSize(0)
expectedNEventsPerGEVHistogram.GetXaxis().SetLabelSize(0)
expectedNEventsPerGEVHistogram.GetXaxis().SetTickLength(0)
expectedNEventsPerGEVHistogram.GetXaxis().SetLabelOffset(999)
expectedNEventsPerGEVHistogram.GetYaxis().SetTitle("Events/GeV")
expectedNEventsPerGEVHistogram.GetYaxis().SetTitleOffset(commonTitleOffset)

expectedNEventsPerGEVHistogramsCopy = expectedNEventsPerGEVHistogram.Clone() # Create clone
expectedNEventsPerGEVHistogramsCopy.SetName("h_copy_expectedNEvents_{n}Jets".format(n=nJetsBin))
expectedNEventsPerGEVHistogramsCopy.SetFillColor(ROOT.kWhite)

expectedNEventsPerGEVGraph.SetFillColor(commonFillColor)

observedNEventsPerGEVGraph.SetLineColor(ROOT.kBlack)
observedNEventsPerGEVGraph.SetFillColor(ROOT.kWhite)

if plot_signal:
    for signalBinIndex in range(len(signalBinSettings[nJetsBin])):
        signalBinSetting = signalBinSettings[nJetsBin][signalBinIndex]
        signalNEventsPerGEVHistograms[signalBinIndex].SetLineColor(signalBinSetting[3])
        signalNEventsPerGEVHistograms[signalBinIndex].SetLineStyle(5)
        signalNEventsPerGEVHistograms[signalBinIndex].SetLineWidth(2)
        signalToDataRatioHistograms[signalBinIndex].SetLineColor(signalBinSetting[3])
        signalToDataRatioHistograms[signalBinIndex].SetLineStyle(5)
        signalToDataRatioHistograms[signalBinIndex].SetLineWidth(2)

# CMS_lumi.writeExtraText = False
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
CMS_lumi.lumi_13TeV = "137.2 fb^{-1}"
CMS_lumi.relPosX    = 0.15

legend = None
if inputArguments.plotObservedData:
    legend = ROOT.TLegend(0.2, 0.85, 0.95, 0.9)
else:
    legend = ROOT.TLegend(0.3, 0.85, 0.95, 0.9)
legend.SetNColumns(3)
nJetsLabel = "N_{{Jets}} = {n}".format(n=nJetsBin)
if (nJetsBin == 6): nJetsLabel = "N_{{Jets}} #geq 6"
legend.AddEntry(None, nJetsLabel, "")
legend.SetBorderSize(0)
legend.SetFillStyle(0)
ROOT.gStyle.SetLegendTextSize(0.05)

expectedNEventsPerGEVHistogram.Draw("][") # First draw filled so that the legend entry is appropriate
backgroundLabel = "Predicted Background"
if inputArguments.plotObservedData:
    backgroundLabel += " (post-fit)"
legend.AddEntry(expectedNEventsPerGEVHistogram, backgroundLabel)
expectedNEventsPerGEVHistogramsCopy.Draw("][") # Next draw with white filling, overwriting previous histogram
expectedNEventsPerGEVHistogramsCopy.GetXaxis().SetRangeUser(STBoundaries[0], STBoundaries[-1])
expectedNEventsPerGEVHistogramsCopy.GetYaxis().SetRangeUser(0.00005, 11.)
expectedNEventsPerGEVGraph.Draw("2") # For the yellow bands
if plot_signal:
    for signalBinIndex in range(len(signalBinSettings[nJetsBin])):
        signalBinSetting = signalBinSettings[nJetsBin][signalBinIndex]
        signalNEventsPerGEVHistograms[signalBinIndex].Draw("A HIST SAME") # Signal distributions
        maxNSignalEvents_xpos = signalNEventsPerGEVHistograms[signalBinIndex].GetBinCenter(signalNEventsPerGEVHistograms[signalBinIndex].GetMaximumBin())
        maxNSignalEvents_ypos = signalNEventsPerGEVHistograms[signalBinIndex].GetBinContent(signalNEventsPerGEVHistograms[signalBinIndex].GetMaximumBin())
        if (signalBinSetting[4] == 11): maxNSignalEvents_xpos += (-0.4)*signalNEventsPerGEVHistograms[signalBinIndex].GetBinWidth(signalNEventsPerGEVHistograms[signalBinIndex].GetMaximumBin()) # For left-aligned labels
        latex = ROOT.TLatex()
        latex.SetTextFont(42)
        latex.SetTextAngle(0)
        latex.SetTextColor(signalBinSetting[3])
        latex.SetTextSize(0.045)
        latex.SetTextAlign(signalBinSetting[4])
        shiftFactor = 1.6
        # if (maxNSignalEvents_ypos < expectedNEventsPerGEVHistogramsCopy.GetBinContent(signalNEventsPerGEVHistograms[signalBinIndex].GetMaximumBin())): shiftFactor = 1.0/1.4
        latex.DrawLatex(maxNSignalEvents_xpos, shiftFactor*maxNSignalEvents_ypos, getSignalBinRawText(signalBinSetting))

if (inputArguments.plotObservedData):
    observedNEventsPerGEVGraph.Draw("0PZ")
    legend.AddEntry(observedNEventsPerGEVGraph, "Data", "PE")
expectedNEventsPerGEVHistogramsCopy.Draw("][ SAME") # Have to draw again to get overlay on top of previous histograms
legend.Draw()
CMS_lumi.CMS_lumi(canvas, 4, 0)
upperPad.cd()
upperPad.Update()
upperPad.RedrawAxis()
frame = upperPad.GetFrame()
frame.Draw()

yTitleSize_upper = expectedNEventsPerGEVHistogramsCopy.GetYaxis().GetTitleSize()
yLabelSize_upper = expectedNEventsPerGEVHistogramsCopy.GetYaxis().GetLabelSize()
yTickLength_upper = expectedNEventsPerGEVHistogramsCopy.GetYaxis().GetTickLength()
upperPad.Update()

lowerPad.cd()
fractionalErrorGraph.GetXaxis().SetTitle("S_{T} (GeV)")
fractionalErrorGraph.GetXaxis().SetTitleSize(yTitleSize_upper/bottomToTopRatio)
fractionalErrorGraph.GetXaxis().SetLabelSize(yLabelSize_upper/bottomToTopRatio)
fractionalErrorGraph.GetXaxis().SetTickLength(yTickLength_upper)
if inputArguments.plotObservedData:
    fractionalErrorGraph.GetYaxis().SetTitle("#frac{Data}{Background}")
else:
    fractionalErrorGraph.GetYaxis().SetTitle("#frac{Signal}{Background}")
fractionalErrorGraph.GetYaxis().SetTitleOffset(1.4*bottomToTopRatio*commonTitleOffset)
fractionalErrorGraph.GetYaxis().SetTitleSize(0.75*yTitleSize_upper/bottomToTopRatio)
fractionalErrorGraph.GetYaxis().SetLabelSize(yLabelSize_upper/bottomToTopRatio)
fractionalErrorGraph.GetYaxis().SetTickLength(yTickLength_upper)
fractionalErrorGraph.GetYaxis().SetNdivisions(2, 0, 0)
fractionalErrorGraph.SetFillColor(commonFillColor)
fractionalErrorGraph.Draw("A2")

if (inputArguments.plotObservedData):
    ratioPlot = tmROOTUtils.getGraphOfRatioOfAsymmErrorsGraphToHistogram(numeratorGraph=observedNEventsPerGEVGraph, denominatorHistogram=expectedNEventsPerGEVHistogram, outputName="g_ratioGraphs_{n}Jets".format(n=nJetsBin), outputTitle="")
    ratioPlot.Draw("0PZ")
else: # Only draw signal histogram ratios if observed data is not to be plotted; otherwise plot looks too cluttered
    if plot_signal:
        for signalBinIndex in range(len(signalBinSettings[nJetsBin])):
            signalToDataRatioHistograms[signalBinIndex].Draw("HIST SAME")
lineAt1 = ROOT.TLine(STBoundaries[0], 1., STBoundaries[-1], 1.)
lineAt1.SetLineColor(commonExpectedEventsLineColor)
lineAt1.SetLineStyle(commonExpectedEventsLineStyle)
lineAt1.SetLineWidth(commonExpectedEventsLineWidth)
lineAt1.Draw()
fractionalErrorGraph.GetXaxis().SetRangeUser(STBoundaries[0], STBoundaries[-1])
if (inputArguments.plotObservedData):
    fractionalErrorGraph.GetYaxis().SetRangeUser(inputArguments.ratioMin, inputArguments.ratioMax)
else:
    if plot_signal:
        fractionalErrorGraph.GetYaxis().SetRangeUser(0.95*minSignalToExpectedFraction, max(2.0, 1.25*maxSignalToExpectedFraction))
    else:
        fractionalErrorGraph.GetYaxis().SetRangeUser(inputArguments.ratioMin, inputArguments.ratioMax)
lowerPad.cd()
lowerPad.Update()
lowerPad.RedrawAxis()
frame = lowerPad.GetFrame()
frame.Draw()

canvas.Update()
canvas.SaveAs("{oD}/{oFN}_{n}Jets.pdf".format(oD=inputArguments.outputDirectory, oFN=inputArguments.outputFilePrefix, n=nJetsBin))

if not(fitDiagnosticsFile is None):
    fitDiagnosticsFile.Close()

if plot_signal:
    for eventProgenitor in ["gluino", "squark"]:
        signalFiles[eventProgenitor].Close()

print("All done!")
