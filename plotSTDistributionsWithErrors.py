#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, array, pdb, math
import ROOT, tmROOTUtils, tmGeneralUtils, tmHEPDataInterface, tdrstyle, CMS_lumi

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
# inputArgumentsParser.add_argument('--path_systematics_dataMCDiscrepancy', required=True, help='Path to file containing estimated systematics on residual MC-data discrepancy.',type=str)
inputArgumentsParser.add_argument('--inputFolder_bkgCompositionUncertainties', required=True, help='Path to directory containing estimated systematics on background composition.',type=str)
inputArgumentsParser.add_argument('--signalType', required=True, choices=['signal', 'signal_loose'], help='Signal type, used while reading in some background composition systematic.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', required=True, help='Output directory.',type=str)
inputArgumentsParser.add_argument('--outputFilePrefix', required=True, help='Name of output file.',type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--nJetsBin', required=True, help='nJets bin for plotting.',type=int)
inputArgumentsParser.add_argument('--plotObservedData', action='store_true', help="By default, this script does not plot observed data, only the expected number of events with error bars; except if this flag is set.")
inputArgumentsParser.add_argument('--bkgType', default="pre", help="Sets the type of background to be plotted. If the argument \"plotObservedData\" is not passed or if \"path_fitDiagnostics\" is not set explicitly, then this argument has no effect. Otherwise, if this argument is \"pre\", then the pre-fit background and uncertainties are plotted in addition to the data; if \"post\", then post-fit background and uncertainties are plotted in addition to the data.")
inputArgumentsParser.add_argument('--suppressSignal', action='store_true', help="By default, this script plots the signal predictions from some chosen signal bins on top of the background predictions; except if this flag is set.")
inputArgumentsParser.add_argument('--ratioMin', default=0.0, help='Max value of ratio to plot.',type=float)
inputArgumentsParser.add_argument('--ratioMax', default=2.5, help='Max value of ratio to plot.',type=float)
inputArguments = inputArgumentsParser.parse_args()

plot_signal = not(inputArguments.suppressSignal)
# Gluino, neutralino mass bins to plot
signalBinSettings = {
    "c": {
        2: [],
        3: [],
        4: [],
        5: [],
        6: []
    },
    "s": {
        2: [],
        3: [],
        4: [("squark", 1150, 200, ROOT.kBlue+2, 5, "above"), ("squark", 1200, 1100, ROOT.kRed+1, 7, "above"), ("gluino", 1800, 900, ROOT.kGreen+3, 7, "above")],
        5: [("squark", 1150, 200, ROOT.kBlue+2, 7, "above"), ("squark", 1200, 1100, ROOT.kRed+1, 6, "above"), ("gluino", 1800, 900, ROOT.kGreen+3, 7, "below")],
        6: [("squark", 1150, 200, ROOT.kBlue+2, 7, "above"), ("squark", 1200, 1100, ROOT.kRed+1, 6, "above"), ("gluino", 1800, 900, ROOT.kGreen+3, 6, "below")]
    },
    "l": {
        2: [],
        3: [],
        4: [("squark", 1150, 200, ROOT.kBlue+2, 5, "below"), ("squark", 1200, 1100, ROOT.kRed+1, 7, "above"), ("gluino", 1800, 900, ROOT.kGreen+3, 7, "above")],
        5: [("squark", 1150, 200, ROOT.kBlue+2, 7, "above"), ("squark", 1200, 1100, ROOT.kRed+1, 6, "above"), ("gluino", 1800, 900, ROOT.kGreen+3, 7, "below")],
        6: [("squark", 1150, 200, ROOT.kBlue+2, 6, "above"), ("squark", 1200, 1100, ROOT.kRed+1, 4, "below"), ("gluino", 1800, 900, ROOT.kGreen+3, 7, "above")]
    }
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

def get_bkg_residual_adjustment_file_path(bkg, shift):
    return ("{i}/ratio_adjustment_all_MC_{b}_shift_{s}_{sT}.dat".format(i=inputArguments.inputFolder_bkgCompositionUncertainties, b=bkg, s=shift, sT=inputArguments.signalType))

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
# systematics_dataMCDiscrepancy = tmGeneralUtils.getConfigurationFromFile(inputArguments.path_systematics_dataMCDiscrepancy)
residual_adjustments_systematic_dict = {}
for bkg in ["Diph", "GJet", "QCD"]:
    residual_adjustments_systematic_dict[bkg] = {}
    for shift in ["up", "down"]:
        residual_adjustments_systematic_dict[bkg][shift] = tmGeneralUtils.getConfigurationFromFile(inputFilePath=get_bkg_residual_adjustment_file_path(bkg, shift))
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

def get_pre_fit_signal(STRegionIndex, nJetsBin):
    input_histogram = ROOT.TH1F()
    fitDiagnosticsFile.GetObject("shapes_prefit/{l}/stealth".format(l=get_bin_label_abbreviated(STRegionIndex, nJetsBin)), input_histogram)
    if (not(input_histogram.GetXaxis().GetNbins() == 1)): sys.exit("ERROR: histogram shapes_prefit/{l}/stealth is in an unexpected format.".format(l=get_bin_label_abbreviated(STRegionIndex, nJetsBin)))
    return tuple([input_histogram.GetBinContent(1), input_histogram.GetBinErrorLow(1), input_histogram.GetBinErrorUp(1)])

def get_post_fit_background(STRegionIndex, nJetsBin):
    input_histogram = ROOT.TH1F()
    fitDiagnosticsFile.GetObject("shapes_fit_b/{l}/qcd".format(l=get_bin_label_abbreviated(STRegionIndex, nJetsBin)), input_histogram)
    if (not(input_histogram.GetXaxis().GetNbins() == 1)): sys.exit("ERROR: histogram shapes_fit_b/{l}/qcd is in an unexpected format.".format(l=get_bin_label_abbreviated(STRegionIndex, nJetsBin)))
    return tuple([input_histogram.GetBinContent(1), input_histogram.GetBinErrorLow(1), input_histogram.GetBinErrorUp(1)])

# whiteColor = ROOT.TColor(9000, 1.0, 1.0, 1.0) # apparently SetFillColor(ROOT.kWhite) does not work (!)
expectedNEvents_preFit_raw = {}
expectedNEvents_postFit_raw = {}
observedNEvents_raw = {}
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
    for signalBinIndex in range(len(signalBinSettings[inputArguments.bin_label_abbreviation][nJetsBin])):
        signalNEventsPerGEVHistograms[signalBinIndex] = ROOT.TH1F("h_signalNEvents_{n}Jets_index{i}".format(n=nJetsBin, i=signalBinIndex), "", n_STBins, array.array('d', STBoundaries))
        signalToDataRatioHistograms[signalBinIndex] = ROOT.TH1F("h_signalToDataRatio_{n}Jets_index{i}".format(n=nJetsBin, i=signalBinIndex), "", n_STBins, array.array('d', STBoundaries))

for STRegionIndex in range(1, 1+STRegionsAxis.GetNbins()):
    expectedNEvents = (expectedEventCounters_data["expectedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)])
    expectedNEventsErrorFromFitDown = None
    expectedNEventsErrorFromFitUp = None
    if (STRegionIndex > 1):
        if ((inputArguments.plotObservedData) and not(fitDiagnosticsFile is None)):
            expectedNEvents_preFit_raw_current_bin, expectedNEvents_preFit_raw_uncDown_current_bin, expectedNEvents_preFit_raw_uncUp_current_bin = get_pre_fit_background(STRegionIndex, nJetsBin)
            expectedNEvents_preFit_raw[STRegionIndex] = tuple([expectedNEvents_preFit_raw_current_bin, expectedNEvents_preFit_raw_uncDown_current_bin, expectedNEvents_preFit_raw_uncUp_current_bin])
            expectedNEvents_postFit_raw_current_bin, expectedNEvents_postFit_raw_uncDown_current_bin, expectedNEvents_postFit_raw_uncUp_current_bin = get_post_fit_background(STRegionIndex, nJetsBin)
            expectedNEvents_postFit_raw[STRegionIndex] = tuple([expectedNEvents_postFit_raw_current_bin, expectedNEvents_postFit_raw_uncDown_current_bin, expectedNEvents_postFit_raw_uncUp_current_bin])
            if (inputArguments.bkgType == "pre"):
                expectedNEvents = expectedNEvents_preFit_raw_current_bin
                expectedNEventsErrorFromFitDown = expectedNEvents_preFit_raw_uncDown_current_bin
                expectedNEventsErrorFromFitUp = expectedNEvents_preFit_raw_uncUp_current_bin
            elif (inputArguments.bkgType == "post"):
                expectedNEvents = expectedNEvents_postFit_raw_current_bin
                expectedNEventsErrorFromFitDown = expectedNEvents_postFit_raw_uncDown_current_bin
                expectedNEventsErrorFromFitUp = expectedNEvents_postFit_raw_uncUp_current_bin
            else:
                sys.exit("ERROR: unrecognized \"bkgType\". Should be either \"pre\" or \"post\", it is currently: {a}".format(a=inputArguments.bkgType))
        else:
            expectedNEvents *= adjustments_data["nominalAdjustment_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
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
        # expectedNEventsError_DataMCDiscrepancy = 2.0*(abs((systematics_dataMCDiscrepancy["ratio_adjustment_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]) - 1.0))
        adjustment_deviationsFrom1 = []
        residual_adjustment_string = "ratio_adjustment_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)
        for bkg in ["Diph", "GJet", "QCD"]:
            for shift in ["up", "down"]:
                adjustment_ratio = residual_adjustments_systematic_dict[bkg][shift][residual_adjustment_string]
                adjustment_deviationsFrom1.append(abs(adjustment_ratio - 1.0))
        residual_systematic_bkgComposition = max(adjustment_deviationsFrom1)
        # expectedNEvents_netFractionalErrorDown = sqrtOfSumOfSquares([expectedNEventsErrorDown_normEvents, expectedNEventsError_shape, expectedNEventsError_rho, expectedNEventsErrorDown_adjustment_mode0, expectedNEventsErrorDown_adjustment_mode1, expectedNEventsError_DataMCDiscrepancy])
        expectedNEvents_netFractionalErrorDown = sqrtOfSumOfSquares([expectedNEventsErrorDown_normEvents, expectedNEventsError_shape, expectedNEventsError_rho, expectedNEventsErrorDown_adjustment_mode0, expectedNEventsErrorDown_adjustment_mode1, residual_systematic_bkgComposition])
        # expectedNEvents_netFractionalErrorUp = sqrtOfSumOfSquares([expectedNEventsErrorUp_normEvents, expectedNEventsError_shape, expectedNEventsError_rho, expectedNEventsErrorUp_adjustment_mode0, expectedNEventsErrorUp_adjustment_mode1, expectedNEventsError_DataMCDiscrepancy])
        expectedNEvents_netFractionalErrorUp = sqrtOfSumOfSquares([expectedNEventsErrorUp_normEvents, expectedNEventsError_shape, expectedNEventsError_rho, expectedNEventsErrorUp_adjustment_mode0, expectedNEventsErrorUp_adjustment_mode1, residual_systematic_bkgComposition])
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
        for signalBinIndex in range(len(signalBinSettings[inputArguments.bin_label_abbreviation][nJetsBin])):
            signalBinSetting = signalBinSettings[inputArguments.bin_label_abbreviation][nJetsBin][signalBinIndex]
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
    observedNEvents_raw[STRegionIndex] = observedNEvents
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
    for signalBinIndex in range(len(signalBinSettings[inputArguments.bin_label_abbreviation][nJetsBin])):
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
upperPad = ROOT.TPad("upperPad_{n}Jets".format(n=nJetsBin), "upperPad_{n}Jets".format(n=nJetsBin), 0., 1.01*bottomFraction, 0.97, 0.97)
upperPad.SetMargin(0.12, 0.03, 0.01, 0.08) # left, right, bottom, top
lowerPad = ROOT.TPad("lowerPad_{n}Jets".format(n=nJetsBin), "lowerPad_{n}Jets".format(n=nJetsBin), 0., 0., 0.97, 0.99*bottomFraction)
lowerPad.SetMargin(0.12, 0.03, 0.38, 0.05) # left, right, bottom, top
upperPad.Draw()
lowerPad.Draw()

commonTitleOffset = 0.9
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
    for signalBinIndex in range(len(signalBinSettings[inputArguments.bin_label_abbreviation][nJetsBin])):
        signalBinSetting = signalBinSettings[inputArguments.bin_label_abbreviation][nJetsBin][signalBinIndex]
        signalNEventsPerGEVHistograms[signalBinIndex].SetLineColor(signalBinSetting[3])
        signalNEventsPerGEVHistograms[signalBinIndex].SetLineStyle(5)
        signalNEventsPerGEVHistograms[signalBinIndex].SetLineWidth(2)
        signalToDataRatioHistograms[signalBinIndex].SetLineColor(signalBinSetting[3])
        signalToDataRatioHistograms[signalBinIndex].SetLineStyle(5)
        signalToDataRatioHistograms[signalBinIndex].SetLineWidth(2)

CMS_lumi.writeExtraText = False
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
CMS_lumi.lumi_13TeV = "138 fb^{-1}"
CMS_lumi.relPosX    = 0.15

legend = None
if inputArguments.plotObservedData:
    legend = ROOT.TLegend(0.2, 0.85, 0.95, 0.9)
else:
    legend = ROOT.TLegend(0.3, 0.85, 0.95, 0.9)
legend.SetNColumns(3)
nJetsLabel = "#it{{N}}_{{jets}} = {n}".format(n=nJetsBin)
if (nJetsBin == 6): nJetsLabel = "#it{N}_{jets} #geq 6"
legend.AddEntry(None, nJetsLabel, "")
legend.SetBorderSize(0)
legend.SetFillStyle(0)
ROOT.gStyle.SetLegendTextSize(0.05)

expectedNEventsPerGEVHistogram.Draw("][") # First draw filled so that the legend entry is appropriate
backgroundLabel = "Predicted background"
if (inputArguments.plotObservedData and not(fitDiagnosticsFile is None)):
    if (inputArguments.bkgType == "post"):
        backgroundLabel += " (post-fit)"
    else:
        backgroundLabel += " (pre-fit)"
legend.AddEntry(expectedNEventsPerGEVHistogram, backgroundLabel)
expectedNEventsPerGEVHistogramsCopy.Draw("][") # Next draw with white filling, overwriting previous histogram
expectedNEventsPerGEVHistogramsCopy.GetXaxis().SetRangeUser(STBoundaries[0], STBoundaries[-1])
expectedNEventsPerGEVHistogramsCopy.GetYaxis().SetRangeUser(0.00005, 11.)
expectedNEventsPerGEVGraph.Draw("2") # For the yellow bands
if plot_signal:
    for signalBinIndex in range(len(signalBinSettings[inputArguments.bin_label_abbreviation][nJetsBin])):
        signalBinSetting = signalBinSettings[inputArguments.bin_label_abbreviation][nJetsBin][signalBinIndex]
        signalNEventsPerGEVHistograms[signalBinIndex].Draw("A HIST SAME") # Signal distributions
        text_xpos = signalNEventsPerGEVHistograms[signalBinIndex].GetBinCenter(signalBinSetting[4])
        text_ypos = signalNEventsPerGEVHistograms[signalBinIndex].GetBinContent(signalBinSetting[4])
        # if (signalBinSetting[4] == 11): text_xpos += (-0.4)*signalNEventsPerGEVHistograms[signalBinIndex].GetBinWidth(signalNEventsPerGEVHistograms[signalBinIndex].GetMaximumBin()) # For left-aligned labels
        latex = ROOT.TLatex()
        latex.SetTextFont(42)
        latex.SetTextAngle(0)
        latex.SetTextColor(signalBinSetting[3])
        latex.SetTextSize(0.045)
        latex.SetTextAlign(22)
        if (signalBinSetting[5] == "above"):
            text_ypos = text_ypos*1.6
        elif (signalBinSetting[5] == "below"):
            text_ypos = text_ypos/2.0
        else:
            sys.exit("ERROR: signal bin setting is in an unexpected format. Expected \"above\" or \"below\", found: {s}".format(s=signalBinSetting[5]))
        latex.DrawLatex(text_xpos, text_ypos, getSignalBinRawText(signalBinSetting))

if (inputArguments.plotObservedData):
    observedNEventsPerGEVGraph.Draw("0PZ")
    legend.AddEntry(observedNEventsPerGEVGraph, "Data", "PE")
expectedNEventsPerGEVHistogramsCopy.Draw("][ SAME") # Have to draw again to get overlay on top of previous histograms
legend.Draw()
CMS_lumi.CMS_lumi(upperPad, 4, 0)
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
fractionalErrorGraph.GetXaxis().SetTitle("#it{S}_{T} (GeV)")
fractionalErrorGraph.GetXaxis().SetTitleSize(yTitleSize_upper/bottomToTopRatio)
fractionalErrorGraph.GetXaxis().SetLabelSize(yLabelSize_upper/bottomToTopRatio)
fractionalErrorGraph.GetXaxis().SetTickLength(yTickLength_upper)
if inputArguments.plotObservedData:
    fractionalErrorGraph.GetYaxis().SetTitle("#frac{obs. events}{pred. bkg.}")
else:
    fractionalErrorGraph.GetYaxis().SetTitle("#frac{signal}{pred. bkg.}")
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
        for signalBinIndex in range(len(signalBinSettings[inputArguments.bin_label_abbreviation][nJetsBin])):
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

if ((inputArguments.bkgType == "post") and
    (inputArguments.plotObservedData) and
    (not(fitDiagnosticsFile is None))):
    # Last step: tabulate observations and background predictions, and save them to a hepdata-formatted yaml file
    data_for_hepdata_yaml = {
        'ST Bin': {
            'units': 'GeV',
            'data': []
        },
        'background prediction (pre-fit)': {
            'units': None,
            'data': []
        },
        'background prediction (post-fit)': {
            'units': None,
            'data': []
        },
        'observation': {
            'units': None,
            'data': []
        }
    }
    indep_vars_for_hepdata_yaml = ['ST Bin']
    dep_vars_for_hepdata_yaml = ['background prediction (pre-fit)', 'background prediction (post-fit)', 'observation']
    out_path_for_hepdata_yaml = '{oD}/bkg_model_{fp}_{n}Jets.yaml'.format(oD=inputArguments.outputDirectory, fp=inputArguments.outputFilePrefix, n=nJetsBin)

    print("Entabulating expectations and observations...")
    output_table_file_name = "{oD}/{oFN}_{n}Jets_table.tex".format(oD=inputArguments.outputDirectory, oFN=inputArguments.outputFilePrefix, n=nJetsBin)
    output_table_file_handle = open(output_table_file_name, "w")
    output_table_file_handle.write("\\begin{tabular}{|m{0.3\\textwidth}|m{0.15\\textwidth}|m{0.15\\textwidth}|m{0.15\\textwidth}|m{0.125\\textwidth}|}\n")
    output_table_file_handle.write("  \\hline\n")
    nJetsString = "{n}".format(n=nJetsBin)
    if (nJetsBin == 6): nJetsString = "\\geq 6"
    output_table_file_handle.write("  \\multicolumn{{5}}{{|c|}}{{Predicted Background and Observations, ${s}$ Jets}} \\\\ \\hline\n".format(s=nJetsString))
    output_table_file_handle.write("  \\st{} range & \\vbox{\\hbox{\\strut background}\\hbox{\\strut prediction}\\hbox{\\strut (pre-fit)}} & \\vbox{\\hbox{\\strut example}\\hbox{\\strut signal model}\\hbox{\\strut prediction}} & \\vbox{\\hbox{\\strut background}\\hbox{\\strut prediction}\\hbox{\\strut (post-fit)}} & observation \\\\ \\hline\n")
    for STRegionIndex in range(2, 1+STRegionsAxis.GetNbins()):
        STMin = STRegionsAxis.GetBinLowEdge(STRegionIndex)
        STRegionString = None
        if (STRegionIndex == STRegionsAxis.GetNbins()):
            STRegionString = "\\st{{}} \\geq {STMin:.0f}\\gev{{}}".format(STMin=STMin)
        else:
            STMax = STRegionsAxis.GetBinUpEdge(STRegionIndex)
            STRegionString = "{STMin:.0f}\\gev{{}} \\leq \\st{{}} < {STMax:.0f}\\gev{{}}".format(STMin=STMin, STMax=STMax)
        bpre = expectedNEvents_preFit_raw[STRegionIndex][0]
        bpreerror = 0.5*(expectedNEvents_preFit_raw[STRegionIndex][1] + expectedNEvents_preFit_raw[STRegionIndex][2])
        spre, spre_error_lo, spre_error_up = get_pre_fit_signal(STRegionIndex, nJetsBin)
        spreerror = 0.5*(spre_error_lo + spre_error_up)
        bpost = expectedNEvents_postFit_raw[STRegionIndex][0]
        bposterror = 0.5*(expectedNEvents_postFit_raw[STRegionIndex][1] + expectedNEvents_postFit_raw[STRegionIndex][2])
        output_table_file_handle.write("  ${s}$ & ${bpre:.2f} \pm {bpreerror:.2f}$ & ${spre:.2f} \pm {spreerror:.2f}$ & ${bpost:.2f} \pm {bposterror:.2f}$ & {obs} \\\\ \\hline\n".format(s=STRegionString, bpre=bpre, bpreerror=bpreerror, spre=spre, spreerror=spreerror, bpost=bpost, bposterror=bposterror, obs=observedNEvents_raw[STRegionIndex]))
        if STRegionIndex == STRegionsAxis.GetNbins():
            data_for_hepdata_yaml['ST Bin']['data'].append(('> {v:.1f}'.format(v=STMin), []))
        else:
            data_for_hepdata_yaml['ST Bin']['data'].append((STMin, STMax, []))
        data_for_hepdata_yaml['background prediction (pre-fit)']['data'].append((bpre, [('total unc.', expectedNEvents_preFit_raw[STRegionIndex][2], -expectedNEvents_preFit_raw[STRegionIndex][1])]))
        data_for_hepdata_yaml['background prediction (post-fit)']['data'].append((bpost, [('total unc.', expectedNEvents_postFit_raw[STRegionIndex][2], -expectedNEvents_postFit_raw[STRegionIndex][1])]))
        data_for_hepdata_yaml['observation']['data'].append((observedNEvents_raw[STRegionIndex], []))
    output_table_file_handle.write("\\end{tabular}\n")
    output_table_file_handle.close()
    print("Table of expectations and observations written to file: {n}".format(n=output_table_file_name))
    tmHEPDataInterface.save_to_yaml(data_for_hepdata_yaml,
                                    indep_vars_for_hepdata_yaml,
                                    dep_vars_for_hepdata_yaml,
                                    out_path_for_hepdata_yaml)
    print("HEPData-formatted yaml output saved to file: {n}".format(n=out_path_for_hepdata_yaml))

if not(fitDiagnosticsFile is None):
    fitDiagnosticsFile.Close()

if plot_signal:
    for eventProgenitor in ["gluino", "squark"]:
        signalFiles[eventProgenitor].Close()

print("All done!")
