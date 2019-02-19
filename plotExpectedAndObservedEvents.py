#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, array, pdb, math
import ROOT, tmROOTUtils, tmGeneralUtils, tdrstyle, CMS_lumi

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
# ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)

def sqrtOfSumOfSquares(listOfNumbers):
    sumOfSquares = 0.
    for number in listOfNumbers:
        sumOfSquares += number*number
    return math.sqrt(sumOfSquares)

tdrstyle.setTDRStyle()

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
expectedNEventsPerGEVHistograms = {} # For "zero error" histograms
expectedNEventsPerGEVHistogramsCopies = {} # Copy with fill color set to white
expectedNEventsPerGEVGraphs = {}
fractionalErrorGraphs = {}
observedNEventsPerGEVGraphs = {}
signalNEventsPerGEVHistograms = {}
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
    signalNEventsPerGEVHistograms[nJetsBin] = ROOT.TH1F("h_signalNEvents_{n}Jets".format(n=nJetsBin), "", n_STBins, array.array('d', STBoundaries))
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
            expectedNEventsError_normEvents = dataSystematics["fractionalUncertainty_normEvents_{n}Jets".format(n=nJetsBin)]
            expectedNEventsError_shape = dataSystematics["fractionalUncertainty_Shape_STRegion{i}".format(i=STRegionIndex)]
            expectedNEventsError_rho = dataSystematics["fractionalUncertainty_rho_STRegion{i}".format(i=STRegionIndex)]
            expectedNEventsError_scaling = 0.
            if (nJetsBin != inputArguments.nJetsNorm): expectedNEventsError_scaling = max(0., dataScalingSystematics["fractionalUncertainty_sTScaling_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)] - expectedNEventsError_shape)
            expectedNEvents_netFractionalError = sqrtOfSumOfSquares([expectedNEventsError_normEvents, expectedNEventsError_shape, expectedNEventsError_rho, expectedNEventsError_scaling])
            expectedNEventsPerGEVGraphs[nJetsBin].SetPointEYlow(STRegionIndex-1, expectedNEvents_netFractionalError*expectedNEvents/STRegionsAxis.GetBinWidth(STRegionIndex))
            expectedNEventsPerGEVGraphs[nJetsBin].SetPointEYhigh(STRegionIndex-1, expectedNEvents_netFractionalError*expectedNEvents/STRegionsAxis.GetBinWidth(STRegionIndex))
            fractionalErrorGraphs[nJetsBin].SetPointEYlow(STRegionIndex-1, expectedNEvents_netFractionalError)
            fractionalErrorGraphs[nJetsBin].SetPointEYhigh(STRegionIndex-1, expectedNEvents_netFractionalError)
        observedNEvents = observedEventCounters_data["observedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
        observedNEventsPerGEVGraphs[nJetsBin].SetPoint(STRegionIndex-1, STRegionsAxis.GetBinCenter(STRegionIndex), observedNEvents/STRegionsAxis.GetBinWidth(STRegionIndex))
        observedNEventsPerGEVGraphs[nJetsBin].SetPointEXlow(STRegionIndex-1, 0.5*STRegionsAxis.GetBinWidth(STRegionIndex))
        observedNEventsPerGEVGraphs[nJetsBin].SetPointEXhigh(STRegionIndex-1, 0.5*STRegionsAxis.GetBinWidth(STRegionIndex))
        poissonInterval = tmROOTUtils.getPoissonConfidenceInterval(observedNEvents=observedNEvents)
        observedNEventsPerGEVGraphs[nJetsBin].SetPointEYlow(STRegionIndex-1, (observedNEvents-poissonInterval["lower"])/STRegionsAxis.GetBinWidth(STRegionIndex))
        observedNEventsPerGEVGraphs[nJetsBin].SetPointEYhigh(STRegionIndex-1, (poissonInterval["upper"]-observedNEvents)/STRegionsAxis.GetBinWidth(STRegionIndex))
        signalNEventsHistogramSource = ROOT.TH2F()
        signalFile.GetObject("h_lumiBasedYearWeightedNEvents_{n}Jets_STRegion{i}".format(n=nJetsBin, i=STRegionIndex), signalNEventsHistogramSource)
        signalNEvents = signalNEventsHistogramSource.GetBinContent(signalNEventsHistogramSource.FindFixBin(inputArguments.gluinoMassToUse, inputArguments.neutralinoMassToUse))
        signalNEventsPerGEVHistograms[nJetsBin].SetBinContent(STRegionIndex, signalNEvents)
        signalNEventsPerGEVHistograms[nJetsBin].SetBinError(STRegionIndex, 0.)

    # print("Before rescaling: ")
    # tmROOTUtils.printHistogramContents(expectedNEventsPerGEVHistograms[nJetsBin])
    # tmROOTUtils.printHistogramContents(signalNEventsPerGEVHistograms[nJetsBin])
    
    tmROOTUtils.rescale1DHistogramByBinWidth(expectedNEventsPerGEVHistograms[nJetsBin])
    tmROOTUtils.rescale1DHistogramByBinWidth(signalNEventsPerGEVHistograms[nJetsBin])

    # print("After rescaling: ")
    # tmROOTUtils.printHistogramContents(expectedNEventsPerGEVHistograms[nJetsBin])
    # tmROOTUtils.printHistogramContents(signalNEventsPerGEVHistograms[nJetsBin])

    H_ref = 600
    W_ref = 800
    W = W_ref
    H  = H_ref
    T = 0.08*H_ref
    B = 0.12*H_ref
    L = 0.12*W_ref
    R = 0.04*W_ref

    canvas = ROOT.TCanvas("c_{oFN}_{n}Jets".format(oD=inputArguments.outputDirectory, oFN=inputArguments.outputFileName, n=nJetsBin), "c_{oFN}_{n}Jets".format(oD=inputArguments.outputDirectory, oFN=inputArguments.outputFileName, n=nJetsBin), 50, 50, W, H)
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

    upperPad = ROOT.TPad("upperPad_{n}Jets".format(n=nJetsBin), "upperPad_{n}Jets".format(n=nJetsBin), 0., 0.33, 0.97, 0.97)
    upperPad.SetMargin(0.1, 0.025, 0.02, 0.08) # left, right, bottom, top
    lowerPad = ROOT.TPad("lowerPad_{n}Jets".format(n=nJetsBin), "lowerPad_{n}Jets".format(n=nJetsBin), 0., 0., 0.97, 0.33)
    lowerPad.SetMargin(0.1, 0.025, 0.3, 0.03) # left, right, bottom, top
    upperPad.Draw()
    lowerPad.Draw()

    commonTitleOffset = 0.7
    commonFillColor = ROOT.kOrange-2
    commonExpectedEventsLineColor = ROOT.kBlack
    commonExpectedEventsLineStyle = 2
    commonExpectedEventsLineWidth = 3
    commonObservedEventsMarkerStyle = ROOT.kPlus

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

    signalNEventsPerGEVHistograms[nJetsBin].SetLineColor(ROOT.kRed)
    signalNEventsPerGEVHistograms[nJetsBin].SetLineStyle(2)

    CMS_lumi.writeExtraText = False
    CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    CMS_lumi.lumi_13TeV = "77.8 fb^{-1}"

    legend = ROOT.TLegend(0.4, 0.7, 0.95, 0.9)
    legend.SetNColumns(2)
    legend.SetBorderSize(0)

    expectedNEventsPerGEVHistograms[nJetsBin].Draw("][") # First draw filled so that the legend entry is appropriate
    expectedNEventsEntry = legend.AddEntry(expectedNEventsPerGEVHistograms[nJetsBin], "Predicted Background")
    expectedNEventsPerGEVHistogramsCopies[nJetsBin].Draw("][") # Next draw with white filling, overwriting previous histogram
    expectedNEventsPerGEVGraphs[nJetsBin].Draw("2")
    signalNEventsPerGEVHistograms[nJetsBin].Draw("A HIST SAME") # Signal distributions
    observedNEventsPerGEVGraphs[nJetsBin].Draw("0PZ")
    observedNEventsPerGEVGraphsEntry = legend.AddEntry(observedNEventsPerGEVGraphs[nJetsBin], "Data", "LPE")
    expectedNEventsPerGEVHistogramsCopies[nJetsBin].Draw("][ SAME") # Have to draw again to get overlay on top of previous histograms
    legend.Draw()
    CMS_lumi.CMS_lumi(canvas, 4, 0)
    upperPad.cd()
    upperPad.Update()
    upperPad.RedrawAxis()
    frame = upperPad.GetFrame()
    frame.Draw()

    expectedNEventsPerGEVHistogramsCopies[nJetsBin].GetXaxis().SetRangeUser(STBoundaries[0], STBoundaries[-1])
    expectedNEventsPerGEVHistogramsCopies[nJetsBin].GetYaxis().SetRangeUser(0.0002, 1.0)

    yTitleSize_upper = expectedNEventsPerGEVHistogramsCopies[nJetsBin].GetYaxis().GetTitleSize()
    yLabelSize_upper = expectedNEventsPerGEVHistogramsCopies[nJetsBin].GetYaxis().GetLabelSize()
    yTickLength_upper = expectedNEventsPerGEVHistogramsCopies[nJetsBin].GetYaxis().GetTickLength()
    upperPad.Update()

    lowerPad.cd()
    fractionalErrorGraphs[nJetsBin].GetXaxis().SetTitle("S_{T} (GeV)")
    fractionalErrorGraphs[nJetsBin].GetYaxis().SetTitle("Ratio")
    fractionalErrorGraphs[nJetsBin].GetYaxis().SetTitleOffset(0.5*commonTitleOffset)

    fractionalErrorGraphs[nJetsBin].GetXaxis().SetTitleSize(2.*yTitleSize_upper)
    fractionalErrorGraphs[nJetsBin].GetXaxis().SetLabelSize(2.*yLabelSize_upper)
    fractionalErrorGraphs[nJetsBin].GetXaxis().SetTickLength(yTickLength_upper)
    fractionalErrorGraphs[nJetsBin].GetYaxis().SetTitleSize(2.*yTitleSize_upper)
    fractionalErrorGraphs[nJetsBin].GetYaxis().SetLabelSize(2.*yLabelSize_upper)
    fractionalErrorGraphs[nJetsBin].GetYaxis().SetTickLength(yTickLength_upper)
    fractionalErrorGraphs[nJetsBin].GetYaxis().SetNdivisions(2, 0, 0)
    fractionalErrorGraphs[nJetsBin].SetFillColor(commonFillColor)
    fractionalErrorGraphs[nJetsBin].Draw("A2")

    ratioPlots[nJetsBin] = tmROOTUtils.getGraphOfRatioOfAsymmErrorsGraphToHistogram(numeratorGraph=observedNEventsPerGEVGraphs[nJetsBin], denominatorHistogram=expectedNEventsPerGEVHistograms[nJetsBin], outputName="g_ratioGraphs_{n}Jets".format(n=nJetsBin), outputTitle="")
    ratioPlots[nJetsBin].Draw("0PZ")
    lineAt1 = ROOT.TLine(STBoundaries[0], 1., STBoundaries[-1], 1.)
    lineAt1.SetLineColor(commonExpectedEventsLineColor)
    lineAt1.SetLineStyle(commonExpectedEventsLineStyle)
    lineAt1.SetLineWidth(commonExpectedEventsLineWidth)
    lineAt1.Draw()
    fractionalErrorGraphs[nJetsBin].GetXaxis().SetRangeUser(STBoundaries[0], STBoundaries[-1])
    fractionalErrorGraphs[nJetsBin].GetYaxis().SetRangeUser(-1., 3.)
    lowerPad.cd()
    lowerPad.Update()
    lowerPad.RedrawAxis()
    frame = lowerPad.GetFrame()
    frame.Draw()

    canvas.Update()
    canvas.SaveAs("{oD}/{oFN}_{n}Jets.png".format(oD=inputArguments.outputDirectory, oFN=inputArguments.outputFileName, n=nJetsBin))

signalFile.Close()
print("All done!")
