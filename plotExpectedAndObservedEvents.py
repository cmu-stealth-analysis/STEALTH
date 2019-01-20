#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, ROOT, tmROOTUtils, array, pdb, tmGeneralUtils, math

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Generate histograms of expected and observed event distributions, based on observed data.')
inputArgumentsParser.add_argument('--expectedNEventsFile', default="analysis/dataSystematics/signal_eventCounters.dat", help='Path to file containing observed number of events in the format "int observedNEvents_STRegionX_YJets=Z".',type=str)
inputArgumentsParser.add_argument('--MCSignalNEventsFile', default="analysis/MCEventHistograms/MC_2018_savedObjects.root", help='Path to ROOT file containing number of events expected from MC samples.',type=str)
inputArgumentsParser.add_argument('--gluinoMassToUse', default=1700., help='Gluino mass at which to plot ST distribution.',type=float)
inputArgumentsParser.add_argument('--neutralinoMassToUse', default=800., help='Neutralino mass at which to plot ST distribution.',type=float)
inputArgumentsParser.add_argument('--observedNEventsFile', default="analysis/dataSystematics/signal_observedEventCounters.dat", help='Path to file containing expected number of events in the format "int expectedNEvents_STRegionX_YJets=Z".',type=str)
inputArgumentsParser.add_argument('--dataSignalSystematicsFile', default="analysis/dataSystematics/signal_dataSystematics.dat", help='Path to file containing estimated systematic due to norm events fractional uncertainty, shape, and rho.',type=str)
inputArgumentsParser.add_argument('--dataControlSystematicsFile', default="analysis/dataSystematics/control_dataSystematics_sTScaling.dat", help='Path to file containing estimated systematic due to possible deviation from ST scaling.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="specialPlots", help='Output directory.',type=str)
inputArgumentsParser.add_argument('--outputFileName', required=True, help='Name of output file.',type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Min number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Max number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsNorm', default=2, help='Number of jets w.r.t. which to normalize the sT distributions for other jets.',type=int)
inputArguments = inputArgumentsParser.parse_args()
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)

def sqrtOfSumOfSquares(listOfNumbers):
    sumOfSquares = 0.
    for number in listOfNumbers:
        sumOfSquares += number*number
    return math.sqrt(sumOfSquares)

ROOT.gROOT.Macro(os.path.expanduser('loadTDR.C'))

STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries)
STBoundaries = []
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()):
        STBoundary = float(STBoundaryString.strip())
        STBoundaries.append(STBoundary)
STBoundaries.append(3500.0) # Instead of infinity
n_STBins = len(STBoundaries) - 1
STRegionsAxis = ROOT.TAxis(n_STBins, array.array('d', STBoundaries))

observedEventCounters_data = tmGeneralUtils.getConfigurationFromFile(inputArguments.observedNEventsFile)
expectedEventCounters_data = tmGeneralUtils.getConfigurationFromFile(inputArguments.expectedNEventsFile)
dataSystematics = tmGeneralUtils.getConfigurationFromFile(inputArguments.dataSignalSystematicsFile)
dataScalingSystematics = tmGeneralUtils.getConfigurationFromFile(inputArguments.dataControlSystematicsFile)
signalFile = ROOT.TFile(inputArguments.MCSignalNEventsFile)
if ((signalFile.IsOpen() == ROOT.kFALSE) or (signalFile.IsZombie())): sys.exit("ERROR: unable to open file with name {n}".format(n=inputArguments.MCSignalNEventsFile))
expectedNEventsHistograms = {}
expectedNEventsHistogramsCopies = {} # For "zero error" histograms
expectedNEventsHistogramsErrorDown = {}
expectedNEventsHistogramsErrorUp = {}
fractionalErrorsDown = {}
fractionalErrorsUp = {}
observedNEventsHistograms = {}
signalNEventsHistograms = {}
ratioPlots = {}
whiteColor = ROOT.TColor(9000, 1.0, 1.0, 1.0) # apparently SetFillColor(ROOT.kWhite) does not work (!)
for nJetsBin in range(inputArguments.nJetsMin, 1+inputArguments.nJetsMax):
    expectedNEventsHistograms[nJetsBin] = ROOT.TH1F("h_expectedNEvents_{n}Jets".format(n=nJetsBin), ";S_{T} (GeV);Events", n_STBins, array.array('d', STBoundaries))
    expectedNEventsHistogramsErrorDown[nJetsBin] = ROOT.TH1F("h_expectedNEvents_errorDown_{n}Jets".format(n=nJetsBin), ";S_{T} (GeV);Events", n_STBins, array.array('d', STBoundaries))
    expectedNEventsHistogramsErrorUp[nJetsBin] = ROOT.TH1F("h_expectedNEvents_errorUp_{n}Jets".format(n=nJetsBin), ";S_{T} (GeV);Events", n_STBins, array.array('d', STBoundaries))
    observedNEventsHistograms[nJetsBin] = ROOT.TH1F("h_observedNEvents_{n}Jets".format(n=nJetsBin), ";S_{T} (GeV);Events", n_STBins, array.array('d', STBoundaries))
    observedNEventsHistograms[nJetsBin].SetBinErrorOption(ROOT.TH1.kPoisson)
    fractionalErrorsDown[nJetsBin] = ROOT.TH1F("h_fractionalErrorsDown_{n}Jets".format(n=nJetsBin), ";S_{T} (GeV);Events", n_STBins, array.array('d', STBoundaries))
    fractionalErrorsUp[nJetsBin] = ROOT.TH1F("h_fractionalErrorsUp_{n}Jets".format(n=nJetsBin), ";S_{T} (GeV);Events", n_STBins, array.array('d', STBoundaries))
    signalNEventsHistograms[nJetsBin] = ROOT.TH1F("h_signalNEvents_{n}Jets".format(n=nJetsBin), ";S_{T} (GeV);Events", n_STBins, array.array('d', STBoundaries))
    for STRegionIndex in range(1, 1+STRegionsAxis.GetNbins()):
        STRegionMidpoint = STRegionsAxis.GetBinCenter(STRegionIndex)
        expectedNEvents = expectedEventCounters_data["expectedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
        expectedNEventsHistograms[nJetsBin].SetBinContent(STRegionIndex, expectedEventCounters_data["expectedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)])
        if (STRegionIndex != 1):
            expectedNEventsError_normEvents = dataSystematics["fractionalUncertainty_normEvents_{n}Jets".format(n=nJetsBin)]
            expectedNEventsError_shape = dataSystematics["fractionalUncertainty_Shape_STRegion{i}".format(i=STRegionIndex)]
            expectedNEventsError_rho = dataSystematics["fractionalUncertainty_rho_STRegion{i}".format(i=STRegionIndex)]
            expectedNEventsError_scaling = 0.
            if (nJetsBin != inputArguments.nJetsNorm): expectedNEventsError_scaling = max(0., dataScalingSystematics["fractionalUncertainty_sTScaling_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)] - dataSystematics["fractionalUncertainty_Shape_STRegion{i}".format(i=STRegionIndex)])
            expectedNEvents_netFractionalError = sqrtOfSumOfSquares([expectedNEventsError_normEvents, expectedNEventsError_shape, expectedNEventsError_rho, expectedNEventsError_scaling])
            expectedNEventsHistograms[nJetsBin].SetBinError(STRegionIndex, expectedNEvents_netFractionalError*expectedEventCounters_data["expectedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)])
            expectedNEventsHistogramsErrorDown[nJetsBin].SetBinContent(STRegionIndex, (1.0 - expectedNEvents_netFractionalError)*expectedEventCounters_data["expectedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)])
            expectedNEventsHistogramsErrorUp[nJetsBin].SetBinContent(STRegionIndex, (1.0 + expectedNEvents_netFractionalError)*expectedEventCounters_data["expectedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)])
            fractionalErrorsDown[nJetsBin].SetBinContent(STRegionIndex, (1.0 - expectedNEvents_netFractionalError))
            fractionalErrorsUp[nJetsBin].SetBinContent(STRegionIndex, (1.0 + expectedNEvents_netFractionalError))
        else:
            expectedNEventsHistograms[nJetsBin].SetBinError(STRegionIndex, 0.)
            expectedNEventsHistogramsErrorDown[nJetsBin].SetBinContent(STRegionIndex, expectedNEventsHistograms[nJetsBin].GetBinContent(STRegionIndex))
            expectedNEventsHistogramsErrorUp[nJetsBin].SetBinContent(STRegionIndex, expectedNEventsHistograms[nJetsBin].GetBinContent(STRegionIndex))
            fractionalErrorsDown[nJetsBin].SetBinContent(STRegionIndex, 1.)
            fractionalErrorsUp[nJetsBin].SetBinContent(STRegionIndex, 1.)
        observedNEvents = observedEventCounters_data["observedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
        for uglyHackCounter in range(0, observedNEvents):
            observedNEventsHistograms[nJetsBin].Fill(STRegionMidpoint) # Ugly hack: looking at the source code, it's not clear whether Poisson errors are returned correctly for weighted histograms
        signalNEventsHistogramSource = ROOT.TH2F()
        signalFile.GetObject("h_lumiBasedYearWeightedNEvents_{n}Jets_STRegion{i}".format(n=nJetsBin, i=STRegionIndex), signalNEventsHistogramSource)
        signalNEvents = signalNEventsHistogramSource.GetBinContent(signalNEventsHistogramSource.FindFixBin(inputArguments.gluinoMassToUse, inputArguments.neutralinoMassToUse))
        signalNEventsHistograms[nJetsBin].SetBinContent(STRegionIndex, signalNEvents)
        signalNEventsHistograms[nJetsBin].SetBinError(STRegionIndex, 0.)
        print("ST region {i}: expected: ({exp} +/- {experror}), observed: ({obs} + {obserrorup} - {obserrorlow}), signal: {sig}".format(i=STRegionIndex, exp=expectedNEventsHistograms[nJetsBin].GetBinContent(STRegionIndex), experror=expectedNEventsHistograms[nJetsBin].GetBinError(STRegionIndex), obs=observedNEventsHistograms[nJetsBin].GetBinContent(STRegionIndex), obserrorup=observedNEventsHistograms[nJetsBin].GetBinErrorUp(STRegionIndex), obserrorlow=observedNEventsHistograms[nJetsBin].GetBinErrorLow(STRegionIndex), sig=signalNEventsHistograms[nJetsBin].GetBinContent(STRegionIndex)))
    expectedNEventsHistogramsCopies[nJetsBin] = expectedNEventsHistograms[nJetsBin].Clone() # Create clone and set its errors to zero
    for STRegionIndex in range(1, 1+STRegionsAxis.GetNbins()):
        expectedNEventsHistogramsCopies[nJetsBin].SetBinError(STRegionIndex, 0.)
    canvas = ROOT.TCanvas("c_{oFN}_{n}Jets".format(oD=inputArguments.outputDirectory, oFN=inputArguments.outputFileName, n=nJetsBin), "c_{oFN}_{n}Jets".format(oD=inputArguments.outputDirectory, oFN=inputArguments.outputFileName, n=nJetsBin), 1024, 768)
    ROOT.gPad.SetLogy()
    ROOT.gStyle.SetOptStat(0)
    expectedNEventsHistogramsCopies[nJetsBin].SetLineColor(ROOT.kBlue)
    
    # expectedNEventsHistograms[nJetsBin].SetFillStyle(3144)
    # expectedNEventsHistograms[nJetsBin].SetFillColor(ROOT.kBlue)
    expectedNEventsHistogramsErrorUp[nJetsBin].SetFillStyle(3144)
    expectedNEventsHistogramsErrorUp[nJetsBin].SetFillColor(ROOT.kBlue)
    expectedNEventsHistogramsErrorUp[nJetsBin].SetLineColorAlpha(ROOT.kBlack, 0.) # ugly hacking to prevent histogram lines from being drawn...
    fractionalErrorsUp[nJetsBin].SetFillStyle(3144)
    fractionalErrorsUp[nJetsBin].SetFillColor(ROOT.kBlue)
    fractionalErrorsUp[nJetsBin].SetLineColorAlpha(ROOT.kBlack, 0.) # ugly hacking to prevent histogram lines from being drawn...
    
    expectedNEventsHistogramsErrorDown[nJetsBin].SetFillStyle(1001) # "Solid"
    expectedNEventsHistogramsErrorDown[nJetsBin].SetFillColor(9000) # Fills everything below lower template with white
    expectedNEventsHistogramsErrorDown[nJetsBin].SetLineColorAlpha(ROOT.kBlack, 0.) # ugly hacking to prevent histogram lines from being drawn...
    fractionalErrorsDown[nJetsBin].SetFillStyle(1001)
    fractionalErrorsDown[nJetsBin].SetFillColor(9000)
    fractionalErrorsDown[nJetsBin].SetLineColorAlpha(ROOT.kBlack, 0.) # ugly hacking to prevent histogram lines from being drawn...

    observedNEventsHistograms[nJetsBin].SetLineColor(ROOT.kBlack)

    signalNEventsHistograms[nJetsBin].SetLineColor(ROOT.kRed)
    signalNEventsHistograms[nJetsBin].SetLineStyle(2)

    ratioPlots[nJetsBin] = ROOT.TRatioPlot(observedNEventsHistograms[nJetsBin], expectedNEventsHistograms[nJetsBin])
    ratioPlots[nJetsBin].Draw()
    ratioPlots[nJetsBin].GetUpperPad().cd() # To change to histograms on top
    expectedNEventsHistogramsCopies[nJetsBin].Draw("A") # For the blue lines
    expectedNEventsHistogramsErrorUp[nJetsBin].Draw("SAME ][") # Region below the upper error bars is filled with blue
    expectedNEventsHistogramsErrorDown[nJetsBin].Draw("SAME ][") # Region below the lower error bars is restored to white
    observedNEventsHistograms[nJetsBin].Draw("SAME") # Observed number of events plotted with correct errors
    signalNEventsHistograms[nJetsBin].Draw("HIST SAME") # Signal distributions
    ratioPlots[nJetsBin].GetLowerPad().cd() # To change to ratio plot on bottom
    fractionalErrorsUp[nJetsBin].Draw("SAME ][") # Region below the upper error bars is filled with blue
    fractionalErrorsDown[nJetsBin].Draw("SAME ][") # Region below the lower error bars is restored to white
    ratioPlots[nJetsBin].GetLowYaxis().SetRangeUser(0., 4.)
    ROOT.gPad.Update()

    canvas.SaveAs("{oD}/{oFN}_{n}Jets.png".format(oD=inputArguments.outputDirectory, oFN=inputArguments.outputFileName, n=nJetsBin))

signalFile.Close()
print("All done!")
