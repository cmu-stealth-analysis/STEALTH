#!/usr/bin/env python

from __future__ import print_function, division

import argparse, pdb, sys, math, array
import ROOT, tmROOTUtils, tmGeneralUtils, tdrstyle, CMS_lumi, MCTemplateReader, stealthEnv, commonFunctions

ROOT.gROOT.SetBatch(ROOT.kTRUE)

inputArgumentsParser = argparse.ArgumentParser(description='Store expected and observed limits on signal strength and cross-section.')
inputArgumentsParser.add_argument('--crossSectionsFile', required=True, help='Path to dat file that contains cross-sections as a function of eventProgenitor mass, to use while weighting events.',type=str)
inputArgumentsParser.add_argument('--MCTemplatePath', required=True, help='Path to MC template.', type=str)
inputArgumentsParser.add_argument('--eventProgenitor', required=True, help="Type of stealth sample. Two possible values: \"squark\" or \"gluino\".", type=str)
inputArgumentsParser.add_argument('--combineResultsDirectory', default="{eP}/{sER}/combineToolOutputs".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot), help='EOS path at which combine tool results can be found.',type=str)
inputArgumentsParser.add_argument('--combineOutputPrefix', default="fullChain", help='Prefix of Higgs combine results.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_rawOutput', default="limits", help='Output directory in which to store raw outputs.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_plots', default="publicationPlots", help='Output directory in which to store plots.',type=str)
inputArgumentsParser.add_argument('--outputSuffix', default="fullChain", help='Suffix to append to all results.',type=str)
inputArgumentsParser.add_argument('--maxAllowedRatio', default=10., help='Max allowed ratio for deviation between expected and observed limits.',type=float)
inputArgumentsParser.add_argument('--minNeutralinoMass', default=-1., help='Min value of the neutralino mass to plot.',type=float)
inputArgumentsParser.add_argument('--contour_signalStrength', default=1., help='Signal strength at which to obtain the contours.',type=float)
inputArgumentsParser.add_argument('--plotObserved', action='store_true', help="If this flag is set, then the observed limits are plotted in addition to the expected limits.")
inputArguments = inputArgumentsParser.parse_args()

string_gluino = "#tilde{g}"
string_mass_gluino = "m_{" + string_gluino + "}"
string_squark = "#tilde{q}"
string_mass_squark = "m_{" + string_squark + "}"
string_neutralino = "#tilde{#chi}_{1}^{0}"
string_mass_neutralino = "m_{" + string_neutralino + "}"
string_singlino = "#tilde{S}"
string_mass_singlino = "m_{" + string_singlino + "}"
string_singlet = "S"
string_mass_singlet = "m_{" + string_singlet + "}"
string_gravitino = "#tilde{G}"
string_mass_gravitino = "m_{" + string_gravitino + "}"
string_photon = "#gamma"
string_eventProgenitor = None
if (inputArguments.eventProgenitor == "gluino"):
    string_eventProgenitor = string_gluino
elif (inputArguments.eventProgenitor == "squark"):
    string_eventProgenitor = string_squark
string_mass_eventProgenitor = "m_{" + string_eventProgenitor + "}"

if (inputArguments.eventProgenitor == "gluino"):
    decayChain = "pp#rightarrow" + string_gluino + string_gluino + ", " + string_gluino + "#rightarrow" + string_squark + "q, " + string_squark + "#rightarrow" + string_neutralino + "q, " + string_neutralino + "#rightarrow" + string_photon + string_singlino + ", " + string_singlino + "#rightarrow" + string_singlet + string_gravitino + ", " + string_singlet + "#rightarrowgg"
else:
    decayChain = "pp#rightarrow" + string_squark + string_squark + ", " + string_squark + "#rightarrow" + string_neutralino + "q, " + string_neutralino + "#rightarrow" + string_photon + string_singlino + ", " + string_singlino + "#rightarrow" + string_singlet + string_gravitino + ", " + string_singlet + "#rightarrowgg"
decayChain_supplementaryInfo1 = "(" + string_mass_singlino + " = 100 GeV, " + string_mass_singlet + " = 90 GeV)"
decayChain_supplementaryInfo2 = "NLO + NLL exclusion"

def formatContours(contoursList, lineStyle, lineWidth, lineColor):
    contoursListIteratorNext = ROOT.TIter(contoursList)
    counter = 0
    while True:
        counter += 1
        contour = contoursListIteratorNext()
        if not(contour): break
        contour.SetLineStyle(lineStyle)
        contour.SetLineWidth(lineWidth)
        contour.SetLineColor(lineColor)

def drawContoursForLegend(ndc_x, delta_ndcx, ndc_y, delta_ndcy, lineWidth_middle, lineStyle_middle, lineColor_middle, lineWidth_topBottom, lineStyle_topBottom, lineColor_topBottom, middleLine_fudgeFactor):
    line_middle = ROOT.TLine()
    line_middle.SetLineWidth(lineWidth_middle)
    line_middle.SetLineStyle(lineStyle_middle)
    line_middle.SetLineColor(lineColor_middle)
    line_topBottom = ROOT.TLine()
    line_topBottom.SetLineWidth(lineWidth_topBottom)
    line_topBottom.SetLineStyle(lineStyle_topBottom)
    line_topBottom.SetLineColor(lineColor_topBottom)
    
    line_middle.DrawLineNDC(ndc_x, ndc_y+middleLine_fudgeFactor*delta_ndcy, ndc_x+delta_ndcx, ndc_y+middleLine_fudgeFactor*delta_ndcy) # middle
    line_topBottom.DrawLineNDC(ndc_x, ndc_y-0.5*delta_ndcy, ndc_x+delta_ndcx, ndc_y-0.5*delta_ndcy) # bottom
    line_topBottom.DrawLineNDC(ndc_x, ndc_y+1.5*delta_ndcy, ndc_x+delta_ndcx, ndc_y+1.5*delta_ndcy) # top

def getRGB(color):
    colorObject = ROOT.gROOT.GetColor(color)
    outputDictionary = {
        "red": colorObject.GetRed(),
        "green": colorObject.GetGreen(),
        "blue": colorObject.GetBlue()
    }
    return outputDictionary

crossSectionsInputFileObject = open(inputArguments.crossSectionsFile, 'r')
crossSectionsDictionary = {}
crossSectionsFractionalUncertaintyDictionary = {}
for line in crossSectionsInputFileObject:
    crossSectionsData = line.split()
    eventProgenitorMass = int(0.5 + float(crossSectionsData[0]))
    crossSection = float(crossSectionsData[1])
    crossSectionFractionalUncertainty = 0.01*float(crossSectionsData[2])
    crossSectionsDictionary[eventProgenitorMass] = crossSection
    crossSectionsFractionalUncertaintyDictionary[eventProgenitorMass] = crossSectionFractionalUncertainty
crossSectionsInputFileObject.close()

templateReader = MCTemplateReader.MCTemplateReader(inputArguments.MCTemplatePath)
minEventProgenitorMass = templateReader.minEventProgenitorMass
maxEventProgenitorMass = templateReader.maxEventProgenitorMass

limitsScanExpected=ROOT.TGraph2D()
limitsScanExpectedOneSigmaDown=ROOT.TGraph2D()
limitsScanExpectedOneSigmaUp=ROOT.TGraph2D()
crossSectionScanExpected = ROOT.TGraph2D()
limitsScanExpected.SetName("limitsScanExpected")
limitsScanExpectedOneSigmaDown.SetName("limitsScanExpectedOneSigmaDown")
limitsScanExpectedOneSigmaUp.SetName("limitsScanExpectedOneSigmaUp")
crossSectionScanExpected.SetName("crossSectionScanExpected")

limitsScanObserved=ROOT.TGraph2D()
limitsScanObservedOneSigmaDown=ROOT.TGraph2D()
limitsScanObservedOneSigmaUp=ROOT.TGraph2D()
crossSectionScanObserved = ROOT.TGraph2D()
limitsScanObserved.SetName("limitsScanObserved")
limitsScanObservedOneSigmaDown.SetName("limitsScanObservedOneSigmaDown")
limitsScanObservedOneSigmaUp.SetName("limitsScanObservedOneSigmaUp")
crossSectionScanObserved.SetName("crossSectionScanObserved")

signalStrengthScan = ROOT.TGraph2D()
signalStrengthScan.SetName("signalStrengthScan")

rateParamNames = []
rateParamBestFitScans = {}
for STRegionIndex in range(2, 8):
    rateParamBestFitScans[STRegionIndex] = {}
    for nJetsBin in range(4, 7):
        rateParamNames.append("rateParam_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin))
        rateParamBestFitScans[STRegionIndex][nJetsBin] = ROOT.TGraph2D()
        rateParamBestFitScans[STRegionIndex][nJetsBin].SetName("rateParamBestFitScan_STRegion{i}_{n}Jets")

expectedCrossSectionLimits = []
observedCrossSectionLimits = []

maxValue_crossSectionScanExpected = -1
minValue_crossSectionScanExpected = -1
maxValue_crossSectionScanObserved = -1
minValue_crossSectionScanObserved = -1

def passesSanityCheck(observedUpperLimits, expectedUpperLimit):
    for observedUpperLimit in observedUpperLimits:
        ratio = observedUpperLimit/expectedUpperLimit
        if ((ratio > inputArguments.maxAllowedRatio) or (ratio < (1.0/inputArguments.maxAllowedRatio))):
            return False
    return True

unavailableBins = []
anomalousBinWarnings = []
minEventProgenitorMassBin = -1
maxEventProgenitorMassBin = -1
minNeutralinoMassBin = -1
maxNeutralinoMassBin = -1
for indexPair in templateReader.nextValidBin():
    eventProgenitorMassBin = indexPair[0]
    eventProgenitorMass = (templateReader.eventProgenitorMasses)[eventProgenitorMassBin]
    neutralinoMassBin = indexPair[1]
    neutralinoMass = (templateReader.neutralinoMasses)[neutralinoMassBin]
    if (neutralinoMass < inputArguments.minNeutralinoMass): continue
    crossSection = crossSectionsDictionary[int(0.5+eventProgenitorMass)]
    print("Analyzing bin at (eventProgenitorMassBin, neutralinoMassBin) = ({gMB}, {nMB}) ==> (eventProgenitorMass, neutralinoMass) = ({gM}, {nM})".format(gMB=eventProgenitorMassBin, gM=eventProgenitorMass, nMB=neutralinoMassBin, nM=neutralinoMass))
    if ((minEventProgenitorMassBin == -1) or (eventProgenitorMassBin < minEventProgenitorMassBin)): minEventProgenitorMassBin = eventProgenitorMassBin
    if ((maxEventProgenitorMassBin == -1) or (eventProgenitorMassBin > maxEventProgenitorMassBin)): maxEventProgenitorMassBin = eventProgenitorMassBin
    if ((minNeutralinoMassBin == -1) or (neutralinoMassBin < minNeutralinoMassBin)): minNeutralinoMassBin = neutralinoMassBin
    if ((maxNeutralinoMassBin == -1) or (neutralinoMassBin > maxNeutralinoMassBin)): maxNeutralinoMassBin = neutralinoMassBin

    # # nominal
    expectedUpperLimit, expectedUpperLimitOneSigmaDown, expectedUpperLimitOneSigmaUp, observedUpperLimit = (-1, -1, -1, -1)
    try:
        expectedUpperLimit, expectedUpperLimitOneSigmaDown, expectedUpperLimitOneSigmaUp, observedUpperLimit = commonFunctions.get_expected_and_observed_limits_from_combine_output(combineOutputFilePath="{cRD}/higgsCombine_{cOP}_eventProgenitorMassBin{gMB}_neutralinoMassBin{nMB}.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=eventProgenitorMassBin, nMB=neutralinoMassBin))
    except ValueError:
        unavailableBins.append((eventProgenitorMassBin, neutralinoMassBin))
        continue

    # cross section down
    observedUpperLimitOneSigmaDown = -1
    try:
        observedUpperLimitOneSigmaDown = commonFunctions.get_observed_limit_from_combine_output(combineOutputFilePath="{cRD}/higgsCombine_{cOP}_crossSectionsDown_eventProgenitorMassBin{gMB}_neutralinoMassBin{nMB}.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=eventProgenitorMassBin, nMB=neutralinoMassBin))
    except ValueError:
        unavailableBins.append((eventProgenitorMassBin, neutralinoMassBin))
        continue

    # cross section up
    observedUpperLimitOneSigmaUp = -1
    try:
        observedUpperLimitOneSigmaUp = commonFunctions.get_observed_limit_from_combine_output(combineOutputFilePath="{cRD}/higgsCombine_{cOP}_crossSectionsUp_eventProgenitorMassBin{gMB}_neutralinoMassBin{nMB}.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=eventProgenitorMassBin, nMB=neutralinoMassBin))
    except ValueError:
        unavailableBins.append((eventProgenitorMassBin, neutralinoMassBin))
        continue

    print("Limits: Observed: ({lobsdown}, {lobs}, {lobsup}); Expected: ({lexpdown}, {lexp}, {lexpup})".format(lobsdown=observedUpperLimitOneSigmaDown, lobs=observedUpperLimit, lobsup=observedUpperLimitOneSigmaUp, lexpdown=expectedUpperLimitOneSigmaDown, lexp=expectedUpperLimit, lexpup=expectedUpperLimitOneSigmaUp))
    if (inputArguments.plotObserved and not(passesSanityCheck(observedUpperLimits=[observedUpperLimit, observedUpperLimitOneSigmaUp, observedUpperLimitOneSigmaDown], expectedUpperLimit=expectedUpperLimit))):
        anomalousBinWarnings.append("WARNING: observed limits deviate too much from expected limits at eventProgenitorMass = {gM}, neutralinoMass={nM}".format(gM=eventProgenitorMass, nM=neutralinoMass))
    limitsScanExpected.SetPoint(limitsScanExpected.GetN(), eventProgenitorMass, neutralinoMass, expectedUpperLimit)
    limitsScanExpectedOneSigmaDown.SetPoint(limitsScanExpectedOneSigmaDown.GetN(), eventProgenitorMass, neutralinoMass, expectedUpperLimitOneSigmaDown)
    limitsScanExpectedOneSigmaUp.SetPoint(limitsScanExpectedOneSigmaUp.GetN(), eventProgenitorMass, neutralinoMass, expectedUpperLimitOneSigmaUp)
    crossSectionScanExpected.SetPoint(crossSectionScanExpected.GetN(), eventProgenitorMass, neutralinoMass, expectedUpperLimit*crossSection)
    expectedCrossSectionLimits.append(((eventProgenitorMass, neutralinoMass), expectedUpperLimit*crossSection))
    if ((minValue_crossSectionScanExpected == -1) or (expectedUpperLimit*crossSection < minValue_crossSectionScanExpected)): minValue_crossSectionScanExpected = expectedUpperLimit*crossSection
    if ((maxValue_crossSectionScanExpected == -1) or (expectedUpperLimit*crossSection > maxValue_crossSectionScanExpected)): maxValue_crossSectionScanExpected = expectedUpperLimit*crossSection

    limitsScanObserved.SetPoint(limitsScanObserved.GetN(), eventProgenitorMass, neutralinoMass, observedUpperLimit)
    limitsScanObservedOneSigmaDown.SetPoint(limitsScanObservedOneSigmaDown.GetN(), eventProgenitorMass, neutralinoMass, observedUpperLimitOneSigmaDown)
    limitsScanObservedOneSigmaUp.SetPoint(limitsScanObservedOneSigmaUp.GetN(), eventProgenitorMass, neutralinoMass, observedUpperLimitOneSigmaUp)
    crossSectionScanObserved.SetPoint(crossSectionScanObserved.GetN(), eventProgenitorMass, neutralinoMass, observedUpperLimit*crossSection)
    observedCrossSectionLimits.append(((eventProgenitorMass, neutralinoMass), observedUpperLimit*crossSection))
    if ((minValue_crossSectionScanObserved == -1) or (observedUpperLimit*crossSection < minValue_crossSectionScanObserved)): minValue_crossSectionScanObserved = observedUpperLimit*crossSection
    if ((maxValue_crossSectionScanObserved == -1) or (observedUpperLimit*crossSection > maxValue_crossSectionScanObserved)): maxValue_crossSectionScanObserved = observedUpperLimit*crossSection

    if (inputArguments.plotObserved): signalStrengthScan.SetPoint(signalStrengthScan.GetN(), eventProgenitorMass, neutralinoMass, observedUpperLimit)
    else: signalStrengthScan.SetPoint(signalStrengthScan.GetN(), eventProgenitorMass, neutralinoMass, expectedUpperLimit)

    print("Now fetching best fit for rate params from multidim output...")
    try:
        rateParam_bestFits = commonFunctions.get_best_fit_rateParams_from_MultiDim_fitResult(multiDimFitResultFilePath="{cRD}/multidimfit_{cOP}_eventProgenitorMassBin{gMB}_neutralinoMassBin{nMB}.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=eventProgenitorMassBin, nMB=neutralinoMassBin), paramNames=rateParamNames)
        print("rateParam_bestFits: {rPbF}".format(rPbF=rateParam_bestFits))
        for STRegionIndex in range(2, 8):
            for nJetsBin in range(4, 7):
                rateParamBestFitScans[STRegionIndex][nJetsBin].SetPoint(rateParamBestFitScans[STRegionIndex][nJetsBin].GetN(), eventProgenitorMass, neutralinoMass, rateParam_bestFits["rateParam_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)])
    except ValueError:
        print("Unable to fetch rate params for this bin.")

for unavailableBin in unavailableBins:
    eventProgenitorMassBin = unavailableBin[0]
    eventProgenitorMass = (templateReader.eventProgenitorMasses)[eventProgenitorMassBin]
    neutralinoMassBin = unavailableBin[1]
    neutralinoMass = (templateReader.neutralinoMasses)[neutralinoMassBin]
    print("WARNING: Limits not available for (eventProgenitorMassBin, neutralinoMassBin) = ({ePMB}, {nMB}) ==> (eventProgenitorMass, neutralinoMass) = ({ePM}, {nM})".format(ePMB=eventProgenitorMassBin, ePM=eventProgenitorMass, nMB=neutralinoMassBin, nM=neutralinoMass))

for anomalousBinWarning in anomalousBinWarnings:
    print(anomalousBinWarning)

del templateReader

outputExpectedCrossSectionsFile=open("{oD}/expectedCrossSections_{s}.txt".format(oD=inputArguments.outputDirectory_rawOutput, s=inputArguments.outputSuffix), 'w')
outputExpectedCrossSectionsFile.write("{gMTitle:<19}{nMTitle:<19}{eXSTitle}\n".format(gMTitle="eventProgenitor mass", nMTitle="neutralino mass", eXSTitle="Expected limits on cross section (pb)"))
for expectedCrossSectionLimit in expectedCrossSectionLimits:
    outputExpectedCrossSectionsFile.write("{gM:<19.1f}{nM:<19.1f}{eXS:.3e}\n".format(gM=expectedCrossSectionLimit[0][0], nM=expectedCrossSectionLimit[0][1], eXS=expectedCrossSectionLimit[1]))
outputExpectedCrossSectionsFile.close()

if (inputArguments.plotObserved):
    outputObservedCrossSectionsFile=open("{oD}/observedCrossSections_{s}.txt".format(oD=inputArguments.outputDirectory_rawOutput, s=inputArguments.outputSuffix), 'w')
    outputObservedCrossSectionsFile.write("{gMTitle:<19}{nMTitle:<19}{eXSTitle}\n".format(gMTitle="eventProgenitor mass", nMTitle="neutralino mass", eXSTitle="Observed limits on cross section (pb)"))
    for observedCrossSectionLimit in observedCrossSectionLimits:
        outputObservedCrossSectionsFile.write("{gM:<19.1f}{nM:<19.1f}{oXS:.3e}\n".format(gM=observedCrossSectionLimit[0][0], nM=observedCrossSectionLimit[0][1], oXS=observedCrossSectionLimit[1]))
    outputObservedCrossSectionsFile.close()

listOf2DScans = [limitsScanExpected, limitsScanExpectedOneSigmaDown, limitsScanExpectedOneSigmaUp, crossSectionScanExpected, limitsScanObserved, limitsScanObservedOneSigmaDown, limitsScanObservedOneSigmaUp, crossSectionScanObserved, signalStrengthScan]
for STRegionIndex in range(2, 8):
    for nJetsBin in range(4, 7):
        listOf2DScans.append(rateParamBestFitScans[STRegionIndex][nJetsBin])
for scan2D in listOf2DScans:
    scan2D.SetNpx(8*(1 + maxEventProgenitorMassBin - minEventProgenitorMassBin))
    scan2D.SetNpy(2*(1 + maxNeutralinoMassBin - minNeutralinoMassBin))

histogramExpectedLimits = limitsScanExpected.GetHistogram()
histogramExpectedLimits.SetName("histogramExpectedLimits")
histogramExpectedLimitsOneSigmaDown = limitsScanExpectedOneSigmaDown.GetHistogram()
histogramExpectedLimitsOneSigmaDown.SetName("histogramExpectedLimitsOneSigmaDown")
histogramExpectedLimitsOneSigmaUp = limitsScanExpectedOneSigmaUp.GetHistogram()
histogramExpectedLimitsOneSigmaUp.SetName("histogramExpectedLimitsOneSigmaUp")
histogramCrossSectionScanExpected = crossSectionScanExpected.GetHistogram()
histogramCrossSectionScanExpected.SetName("histogramCrossSectionScanExpected")

histogramObservedLimits = limitsScanObserved.GetHistogram()
histogramObservedLimits.SetName("histogramObservedLimits")
histogramObservedLimitsOneSigmaDown = limitsScanObservedOneSigmaDown.GetHistogram()
histogramObservedLimitsOneSigmaDown.SetName("histogramObservedLimitsOneSigmaDown")
histogramObservedLimitsOneSigmaUp = limitsScanObservedOneSigmaUp.GetHistogram()
histogramObservedLimitsOneSigmaUp.SetName("histogramObservedLimitsOneSigmaUp")
histogramCrossSectionScanObserved = crossSectionScanObserved.GetHistogram()
histogramCrossSectionScanObserved.SetName("histogramCrossSectionScanObserved")

expectedLimitContours = limitsScanExpected.GetContourList(inputArguments.contour_signalStrength)
expectedLimitContours.SetName("expectedLimitContours")
expectedLimitContoursOneSigmaDown = limitsScanExpectedOneSigmaDown.GetContourList(inputArguments.contour_signalStrength)
expectedLimitContoursOneSigmaDown.SetName("expectedLimitContoursOneSigmaDown")
expectedLimitContoursOneSigmaUp = limitsScanExpectedOneSigmaUp.GetContourList(inputArguments.contour_signalStrength)
expectedLimitContoursOneSigmaUp.SetName("expectedLimitContoursOneSigmaUp")

histogramSignalStrengthScan = signalStrengthScan.GetHistogram()
histogramSignalStrengthScan.SetName("histogramSignalStrengthScan")

histogram_rateParamBestFitScans = {}
for STRegionIndex in range(2, 8):
    histogram_rateParamBestFitScans[STRegionIndex] = {}
    for nJetsBin in range(4, 7):
        histogram_rateParamBestFitScans[STRegionIndex][nJetsBin] = rateParamBestFitScans[STRegionIndex][nJetsBin].GetHistogram()
        histogram_rateParamBestFitScans[STRegionIndex][nJetsBin].SetName("best-fit rate parameter: ST region index {i}, {n} Jets".format(i=STRegionIndex, n=nJetsBin))

if (inputArguments.plotObserved):
    observedLimitContours = limitsScanObserved.GetContourList(inputArguments.contour_signalStrength)
    observedLimitContours.SetName("observedLimitContours")
    observedLimitContoursOneSigmaDown = limitsScanObservedOneSigmaDown.GetContourList(inputArguments.contour_signalStrength)
    observedLimitContoursOneSigmaDown.SetName("observedLimitContoursOneSigmaDown")
    observedLimitContoursOneSigmaUp = limitsScanObservedOneSigmaUp.GetContourList(inputArguments.contour_signalStrength)
    observedLimitContoursOneSigmaUp.SetName("observedLimitContoursOneSigmaUp")

color_expectedContours = ROOT.kRed
width_expectedContours_middle = 5
width_expectedContours_topBottom = 2
style_expectedContours_middle = 2
style_expectedContours_topBottom = 2
color_observedContours = ROOT.kBlack
width_observedContours_middle = 5
width_observedContours_topBottom = 2
style_observedContours_middle = 1
style_observedContours_topBottom = 1

formatContours(expectedLimitContours, style_expectedContours_middle, width_expectedContours_middle, color_expectedContours)
formatContours(expectedLimitContoursOneSigmaDown, style_expectedContours_topBottom, width_expectedContours_topBottom, color_expectedContours)
formatContours(expectedLimitContoursOneSigmaUp, style_expectedContours_topBottom, width_expectedContours_topBottom, color_expectedContours)
if (inputArguments.plotObserved):
    formatContours(observedLimitContours, style_observedContours_middle, width_observedContours_middle, color_observedContours)
    formatContours(observedLimitContoursOneSigmaDown, style_observedContours_topBottom, width_observedContours_topBottom, color_observedContours)
    formatContours(observedLimitContoursOneSigmaUp, style_observedContours_topBottom, width_observedContours_topBottom, color_observedContours)

# CMS_lumi.writeExtraText = False
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
CMS_lumi.lumi_13TeV = "137.2 fb^{-1}"
CMS_lumi.relPosX    = 0.15

H_ref = 600
W_ref = 800
W = W_ref
H  = H_ref
T = 0.08*H_ref
B = 0.12*H_ref
L = 0.12*W_ref
R = 0.04*W_ref

tdrstyle.setTDRStyle()
canvas = ROOT.TCanvas("c_{s}_observedLimits".format(s=inputArguments.outputSuffix), "c_{s}_observedLimits".format(s=inputArguments.outputSuffix), 50, 50, W, H)
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

ROOT.gPad.SetLogz()

# ROOT.gStyle.SetPalette(ROOT.kBird)
colorMin_RGB = getRGB(ROOT.kBlue+1)
colorMid_RGB = getRGB(ROOT.kYellow)
colorMax_RGB = getRGB(ROOT.kRed)
paletteRed = array.array('d', [colorMin_RGB["red"], colorMid_RGB["red"], colorMax_RGB["red"]])
paletteGreen = array.array('d', [colorMin_RGB["green"], colorMid_RGB["green"], colorMax_RGB["green"]])
paletteBlue = array.array('d', [colorMin_RGB["blue"], colorMid_RGB["blue"], colorMax_RGB["blue"]])
paletteStops = None
if (inputArguments.plotObserved):
    paletteStops = array.array('d', [0., math.log10(3.)/math.log10(maxValue_crossSectionScanObserved/minValue_crossSectionScanObserved), 1.]) # "mid" color is set to 5 times the min cross-section -- note that the axis is log-scaled
else:
    paletteStops = array.array('d', [0., math.log10(3.)/math.log10(maxValue_crossSectionScanExpected/minValue_crossSectionScanExpected), 1.]) # "mid" color is set to 5 times the min cross-section -- note that the axis is log-scaled

ROOT.TColor.CreateGradientColorTable(len(paletteStops), paletteStops, paletteRed, paletteGreen, paletteBlue, 999)
ROOT.gStyle.SetNumberContours(999)

ROOT.gPad.SetRightMargin(0.2)
ROOT.gPad.SetLeftMargin(0.15)
commonTitleSize = 0.046
if (inputArguments.plotObserved):
    histogramCrossSectionScanObserved.GetXaxis().SetTitle(string_mass_eventProgenitor + "(GeV)")
    histogramCrossSectionScanObserved.GetXaxis().SetTitleSize(commonTitleSize)
    histogramCrossSectionScanObserved.GetYaxis().SetTitle(string_mass_neutralino + "(GeV)")
    histogramCrossSectionScanObserved.GetYaxis().SetTitleOffset(1.)
    histogramCrossSectionScanObserved.GetYaxis().SetTitleSize(commonTitleSize)
    histogramCrossSectionScanObserved.GetZaxis().SetTitle("95% CL upper limit on cross-section (pb)")
    histogramCrossSectionScanObserved.GetZaxis().SetTitleOffset(1.)
    histogramCrossSectionScanObserved.GetZaxis().SetTitleSize(0.046)
    histogramCrossSectionScanObserved.GetXaxis().SetRangeUser(minEventProgenitorMass, maxEventProgenitorMass)
    histogramCrossSectionScanObserved.GetZaxis().SetRangeUser(minValue_crossSectionScanObserved, maxValue_crossSectionScanObserved)
    histogramCrossSectionScanObserved.Draw("colz")
else:
    histogramCrossSectionScanExpected.GetXaxis().SetTitle(string_mass_eventProgenitor + "(GeV)")
    histogramCrossSectionScanExpected.GetXaxis().SetTitleSize(commonTitleSize)
    histogramCrossSectionScanExpected.GetYaxis().SetTitle(string_mass_neutralino + "(GeV)")
    histogramCrossSectionScanExpected.GetYaxis().SetTitleOffset(1.)
    histogramCrossSectionScanExpected.GetYaxis().SetTitleSize(commonTitleSize)
    histogramCrossSectionScanExpected.GetZaxis().SetTitle("95% CL upper limit on cross-section (pb)")
    histogramCrossSectionScanExpected.GetZaxis().SetTitleOffset(1.)
    histogramCrossSectionScanExpected.GetZaxis().SetTitleSize(0.046)
    histogramCrossSectionScanExpected.GetXaxis().SetRangeUser(minEventProgenitorMass, maxEventProgenitorMass)
    histogramCrossSectionScanExpected.GetZaxis().SetRangeUser(minValue_crossSectionScanExpected, maxValue_crossSectionScanExpected)
    histogramCrossSectionScanExpected.Draw("colz")

contoursToDraw = [expectedLimitContours, expectedLimitContoursOneSigmaDown, expectedLimitContoursOneSigmaUp]
if (inputArguments.plotObserved):
    contoursToDraw.extend([observedLimitContours, observedLimitContoursOneSigmaDown, observedLimitContoursOneSigmaUp])

for contoursList in contoursToDraw:
    contoursList.Draw("SAME")

line_eventProgenitorEqualsNeutralinoMass = ROOT.TLine(minEventProgenitorMass, minEventProgenitorMass, maxEventProgenitorMass, maxEventProgenitorMass)
line_eventProgenitorEqualsNeutralinoMass.SetLineStyle(7)
line_eventProgenitorEqualsNeutralinoMass.SetLineColor(ROOT.kBlack)
line_eventProgenitorEqualsNeutralinoMass.SetLineWidth(3)
line_eventProgenitorEqualsNeutralinoMass.Draw()
ROOT.gPad.Update()

commonOffset = 0.178

latex = ROOT.TLatex()
latex.SetTextFont(42)
latex.SetTextAlign(12)
latex.SetTextColor(ROOT.kBlack)

latex.SetTextSize(0.033)
latex.DrawLatexNDC(commonOffset, 0.89, decayChain)

latex.SetTextSize(0.025)
latex.DrawLatexNDC(commonOffset, 0.845, decayChain_supplementaryInfo1)

latex.SetTextSize(0.033)
latex.DrawLatexNDC(commonOffset, 0.815, decayChain_supplementaryInfo2)

drawContoursForLegend(commonOffset, 0.03, 0.765, 0.01, width_expectedContours_middle, style_expectedContours_middle, color_expectedContours, width_expectedContours_topBottom, style_expectedContours_topBottom, color_expectedContours, 0.59)
latex.SetTextSize(0.03)
latex.SetTextColor(color_expectedContours)
latex.DrawLatexNDC(commonOffset+0.04, 0.765, "Expected #pm 1#sigma_{experiment}")

drawContoursForLegend(commonOffset, 0.03, 0.722, 0.01, width_observedContours_middle, style_observedContours_middle, color_observedContours, width_observedContours_topBottom, style_observedContours_topBottom, color_observedContours, 0.5775)
latex.SetTextSize(0.03)
latex.SetTextColor(color_observedContours)
latex.DrawLatexNDC(commonOffset+0.04, 0.722, "Observed #pm 1#sigma_{theory}")

latex.SetTextAlign(22)
latex.SetTextColor(ROOT.kBlack)
latex.SetTextSize(0.04)
latex.SetTextAngle(tmROOTUtils.getTLineAngleInDegrees(ROOT.gPad, line_eventProgenitorEqualsNeutralinoMass))
latex.DrawLatex(minEventProgenitorMass + 85., minEventProgenitorMass + 150., string_mass_eventProgenitor + " = " + string_mass_neutralino)

CMS_lumi.CMS_lumi(canvas, 4, 0)

ROOT.gPad.Update()
ROOT.gPad.RedrawAxis()
frame = ROOT.gPad.GetFrame()
frame.Draw()
canvas.Update()
if (inputArguments.plotObserved):
    canvas.SaveAs("{oD}/{s}_observedLimits.pdf".format(oD=inputArguments.outputDirectory_plots, s=inputArguments.outputSuffix))
else:
    canvas.SaveAs("{oD}/{s}_expectedLimits.pdf".format(oD=inputArguments.outputDirectory_plots, s=inputArguments.outputSuffix))

signalStrengthCanvas = ROOT.TCanvas("c_{s}_signalStrengthScan".format(s=inputArguments.outputSuffix), "c_{s}_signalStrengthScan".format(s=inputArguments.outputSuffix), 50, 50, W, H)
signalStrengthCanvas.SetFillColor(0)
signalStrengthCanvas.SetBorderMode(0)
signalStrengthCanvas.SetFrameFillStyle(0)
signalStrengthCanvas.SetFrameBorderMode(0)
signalStrengthCanvas.SetLeftMargin( L/W )
signalStrengthCanvas.SetRightMargin( R/W )
signalStrengthCanvas.SetTopMargin( T/H )
signalStrengthCanvas.SetBottomMargin( B/H )
signalStrengthCanvas.SetTickx(0)
signalStrengthCanvas.SetTicky(0)
ROOT.gPad.SetRightMargin(0.2)
ROOT.gPad.SetLeftMargin(0.15)
signalStrengthCanvas.Draw()
histogramSignalStrengthScan.GetXaxis().SetTitle(string_mass_eventProgenitor + "(GeV)")
histogramSignalStrengthScan.GetXaxis().SetTitleSize(commonTitleSize)
histogramSignalStrengthScan.GetXaxis().SetRangeUser(minEventProgenitorMass, maxEventProgenitorMass)
histogramSignalStrengthScan.GetYaxis().SetTitle(string_mass_neutralino + "(GeV)")
histogramSignalStrengthScan.GetYaxis().SetTitleOffset(1.)
histogramSignalStrengthScan.GetYaxis().SetTitleSize(commonTitleSize)
if (inputArguments.plotObserved): histogramSignalStrengthScan.GetZaxis().SetTitle("Upper limit on observed signal strength.")
else: histogramSignalStrengthScan.GetZaxis().SetTitle("Upper limit on expected signal strength.")
histogramSignalStrengthScan.GetZaxis().SetTitleOffset(1.)
histogramSignalStrengthScan.GetZaxis().SetTitleSize(0.046)
ROOT.gPad.SetLogz()
histogramSignalStrengthScan.Draw("colz")
signalStrengthCanvas.SaveAs("{oD}/{s}_signalStrength.pdf".format(oD=inputArguments.outputDirectory_plots, s=inputArguments.outputSuffix))

for STRegionIndex in range(2, 8):
    for nJetsBin in range(4, 7):
        bestFitRateParamCanvas = ROOT.TCanvas("c_{s}_bestFitRateParam_STRegion{i}_{n}Jets".format(s=inputArguments.outputSuffix, i=STRegionIndex, n=nJetsBin), "c_{s}_bestFitRateParam_STRegion{i}_{n}Jets".format(s=inputArguments.outputSuffix, i=STRegionIndex, n=nJetsBin), 50, 50, W, H)
        bestFitRateParamCanvas.SetFillColor(0)
        bestFitRateParamCanvas.SetBorderMode(0)
        bestFitRateParamCanvas.SetFrameFillStyle(0)
        bestFitRateParamCanvas.SetFrameBorderMode(0)
        bestFitRateParamCanvas.SetLeftMargin( L/W )
        bestFitRateParamCanvas.SetRightMargin( R/W )
        bestFitRateParamCanvas.SetTopMargin( T/H )
        bestFitRateParamCanvas.SetBottomMargin( B/H )
        bestFitRateParamCanvas.SetTickx(0)
        bestFitRateParamCanvas.SetTicky(0)
        ROOT.gPad.SetRightMargin(0.2)
        ROOT.gPad.SetLeftMargin(0.15)
        bestFitRateParamCanvas.Draw()
        histogram_rateParamBestFitScans[STRegionIndex][nJetsBin].GetXaxis().SetTitle(string_mass_eventProgenitor + "(GeV)")
        histogram_rateParamBestFitScans[STRegionIndex][nJetsBin].GetXaxis().SetTitleSize(commonTitleSize)
        histogram_rateParamBestFitScans[STRegionIndex][nJetsBin].GetXaxis().SetRangeUser(minEventProgenitorMass, maxEventProgenitorMass)
        histogram_rateParamBestFitScans[STRegionIndex][nJetsBin].GetYaxis().SetTitle(string_mass_neutralino + "(GeV)")
        histogram_rateParamBestFitScans[STRegionIndex][nJetsBin].GetYaxis().SetTitleOffset(1.)
        histogram_rateParamBestFitScans[STRegionIndex][nJetsBin].GetYaxis().SetTitleSize(commonTitleSize)
        histogram_rateParamBestFitScans[STRegionIndex][nJetsBin].GetZaxis().SetTitle("Best-fit value.")
        histogram_rateParamBestFitScans[STRegionIndex][nJetsBin].GetZaxis().SetTitleOffset(1.)
        histogram_rateParamBestFitScans[STRegionIndex][nJetsBin].GetZaxis().SetTitleSize(0.046)
        histogram_rateParamBestFitScans[STRegionIndex][nJetsBin].SetMaximum(10.0)
        histogram_rateParamBestFitScans[STRegionIndex][nJetsBin].SetMinimum(0.1)
        ROOT.gPad.SetLogz()
        histogram_rateParamBestFitScans[STRegionIndex][nJetsBin].Draw("colz")
        bestFitRateParamCanvas.SaveAs("{oD}/{s}_bestFitRateParam_STRegion{i}_{n}Jets.pdf".format(oD=inputArguments.outputDirectory_plots, s=inputArguments.outputSuffix, i=STRegionIndex, n=nJetsBin))

outputFileName = "{oD}/limits_{suffix}.root".format(oD=inputArguments.outputDirectory_rawOutput, suffix=inputArguments.outputSuffix)
outputFile=ROOT.TFile.Open(outputFileName, "RECREATE")
tObjectsToSave = [limitsScanExpected, limitsScanExpectedOneSigmaUp, limitsScanExpectedOneSigmaDown, crossSectionScanExpected, histogramExpectedLimits, histogramExpectedLimitsOneSigmaDown, histogramExpectedLimitsOneSigmaUp, histogramCrossSectionScanExpected, expectedLimitContours, expectedLimitContoursOneSigmaDown, expectedLimitContoursOneSigmaUp, histogramSignalStrengthScan, canvas, signalStrengthCanvas]
if (inputArguments.plotObserved):
    tObjectsToSave.extend([limitsScanObserved, limitsScanObservedOneSigmaUp, limitsScanObservedOneSigmaDown, crossSectionScanObserved, histogramObservedLimits, histogramObservedLimitsOneSigmaDown, histogramObservedLimitsOneSigmaUp, histogramCrossSectionScanObserved, observedLimitContours, observedLimitContoursOneSigmaDown, observedLimitContoursOneSigmaUp])
for tObject in tObjectsToSave:
    outputFile.WriteTObject(tObject)

outputFile.Close()
print("All done!")
