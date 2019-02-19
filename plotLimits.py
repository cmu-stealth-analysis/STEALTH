#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, tmROOTUtils, pdb, tmGeneralUtils, sys, tdrstyle, CMS_lumi

inputArgumentsParser = argparse.ArgumentParser(description='Store expected and observed limits on signal strength and cross-section.')
inputArgumentsParser.add_argument('--crossSectionsFile', default="SusyCrossSections13TevGluGlu.txt", help='Path to dat file that contains cross-sections as a function of gluino mass, to use while weighting events.',type=str)
inputArgumentsParser.add_argument('--MCTemplate', default="plot_susyMasses_template.root", help='Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.', type=str)
inputArgumentsParser.add_argument('--combineResultsDirectory', default="root://cmseos.fnal.gov//store/user/lpcsusystealth/combineToolOutputs", help='EOS path at which combine tool results can be found.',type=str)
inputArgumentsParser.add_argument('--combineOutputPrefix', default="fullChain", help='Prefix of Higgs combine results.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_rawOutput', default="limits", help='Output directory in which to store raw outputs.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_plots', default="publicationPlots", help='Output directory in which to store plots.',type=str)
inputArgumentsParser.add_argument('--outputSuffix', default="fullChain", help='Suffix to append to all results.',type=str)
inputArgumentsParser.add_argument('--minGluinoMass', default=1000., help='Minimum gluino mass on which to run.', type=float)
inputArgumentsParser.add_argument('--maxGluinoMass', default=1750., help='Max gluino mass for the 2D plots.',type=float)
inputArgumentsParser.add_argument('--maxAllowedRatio', default=10., help='Max allowed ratio for deviation between expected and observed limits.',type=float)
inputArgumentsParser.add_argument('--contour_signalStrength', default=1., help='Signal strength at which to obtain the contours.',type=float)
inputArguments = inputArgumentsParser.parse_args()

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

crossSectionsInputFileObject = open(inputArguments.crossSectionsFile, 'r')
crossSectionsDictionary = {}
crossSectionsFractionalUncertaintyDictionary = {}
for line in crossSectionsInputFileObject:
    crossSectionsData = line.split()
    gluinoMass = int(0.5 + float(crossSectionsData[0]))
    crossSection = float(crossSectionsData[1])
    crossSectionFractionalUncertainty = 0.01*float(crossSectionsData[2])
    crossSectionsDictionary[gluinoMass] = crossSection
    crossSectionsFractionalUncertaintyDictionary[gluinoMass] = crossSectionFractionalUncertainty
crossSectionsInputFileObject.close()

# print("Read in cross-sections as a function of gluino mass:")
# tmGeneralUtils.prettyPrintDictionary(crossSectionsDictionary)

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

generatedMCTemplate = ROOT.TFile(inputArguments.MCTemplate, "READ")
h_MCTemplate = generatedMCTemplate.Get("h_susyMasses_template")
for gluinoMassBin in range(1, 1+h_MCTemplate.GetXaxis().GetNbins()):
    gluinoMass = int(0.5 + h_MCTemplate.GetXaxis().GetBinCenter(gluinoMassBin))
    if (inputArguments.minGluinoMass > 0 and gluinoMass < inputArguments.minGluinoMass): continue
    for neutralinoMassBin in range(1, 1+h_MCTemplate.GetYaxis().GetNbins()):
        if not(int(0.5 + h_MCTemplate.GetBinContent(gluinoMassBin, neutralinoMassBin)) == 1): continue
        neutralinoMass = int(0.5 + h_MCTemplate.GetYaxis().GetBinCenter(neutralinoMassBin))
        crossSection = crossSectionsDictionary[gluinoMass]
        print("Analyzing bin at gluino mass bin={gMB} (gluino mass = {gM}), neutralino mass bin={nMB} (neutralino mass = {nM})".format(gMB=gluinoMassBin, gM=gluinoMass, nMB=neutralinoMassBin, nM=neutralinoMass))

        # nominal
        combineOutputFile=ROOT.TFile.Open("{cRD}/higgsCombine_{cOP}_gluinoMassBin{gMB}_neutralinoMassBin{nMB}.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin), "READ")
        if ((combineOutputFile.IsZombie() == ROOT.kTRUE) or not(combineOutputFile.IsOpen() == ROOT.kTRUE)):
            sys.exit("Error in opening file: {cRD}/higgsCombine_{cOP}_gluinoMassBin{gMB}_neutralinoMassBin{nMB}.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin))
        limitTree = ROOT.TTree()
        combineOutputFile.GetObject("limit", limitTree)
        nEntriesFound = limitTree.GetEntries()
        if not(nEntriesFound == 6): sys.exit("ERROR: limits not in proper format.")
        limitTree.GetEntry(2)
        expectedUpperLimit = limitTree.limit
        limitTree.GetEntry(1)
        expectedUpperLimitOneSigmaDown=limitTree.limit
        limitTree.GetEntry(3)
        expectedUpperLimitOneSigmaUp=limitTree.limit
        limitTree.GetEntry(5)
        observedUpperLimit = limitTree.limit
        combineOutputFile.Close()
        observedUpperLimitOneSigmaDown = 0.
        observedUpperLimitOneSigmaUp = 0.

        # cross section down
        combineOutputFile_crossSectionsDown=ROOT.TFile.Open("{cRD}/higgsCombine_{cOP}_gluinoMassBin{gMB}_neutralinoMassBin{nMB}_crossSectionsDown.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin), "READ")
        if ((combineOutputFile_crossSectionsDown.IsZombie() == ROOT.kTRUE) or not(combineOutputFile_crossSectionsDown.IsOpen() == ROOT.kTRUE)):
            sys.exit("Error in opening file: {cRD}/higgsCombine_{cOP}_gluinoMassBin{gMB}_neutralinoMassBin{nMB}.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin))
        limitTree_crossSectionsDown = ROOT.TTree()
        combineOutputFile_crossSectionsDown.GetObject("limit", limitTree_crossSectionsDown)
        nEntriesFound = limitTree_crossSectionsDown.GetEntries()
        if not(nEntriesFound == 6): sys.exit("ERROR: limits not in proper format.")
        limitTree_crossSectionsDown.GetEntry(5)
        observedUpperLimitOneSigmaUp = limitTree_crossSectionsDown.limit
        combineOutputFile_crossSectionsDown.Close()

        # cross section up
        combineOutputFile_crossSectionsUp=ROOT.TFile.Open("{cRD}/higgsCombine_{cOP}_gluinoMassBin{gMB}_neutralinoMassBin{nMB}_crossSectionsUp.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin), "READ")
        if ((combineOutputFile_crossSectionsUp.IsZombie() == ROOT.kTRUE) or not(combineOutputFile_crossSectionsUp.IsOpen() == ROOT.kTRUE)):
            sys.exit("Error in opening file: {cRD}/higgsCombine_{cOP}_gluinoMassBin{gMB}_neutralinoMassBin{nMB}.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin))
        limitTree_crossSectionsUp = ROOT.TTree()
        combineOutputFile_crossSectionsUp.GetObject("limit", limitTree_crossSectionsUp)
        nEntriesFound = limitTree_crossSectionsUp.GetEntries()
        if not(nEntriesFound == 6): sys.exit("ERROR: limits not in proper format.")
        limitTree_crossSectionsUp.GetEntry(5)
        observedUpperLimitOneSigmaDown = limitTree_crossSectionsUp.limit
        combineOutputFile_crossSectionsUp.Close()

        print("Limits: Observed: ({lobsdown}, {lobs}, {lobsup}); Expected: ({lexpdown}, {lexp}, {lexpup})".format(lobsdown=observedUpperLimitOneSigmaDown, lobs=observedUpperLimit, lobsup=observedUpperLimitOneSigmaUp, lexpdown=expectedUpperLimitOneSigmaDown, lexp=expectedUpperLimit, lexpup=expectedUpperLimitOneSigmaUp))
        observedLimitsAreSane = passesSanityCheck(observedUpperLimits=[observedUpperLimit, observedUpperLimitOneSigmaUp, observedUpperLimitOneSigmaDown], expectedUpperLimit=expectedUpperLimit)
        if not(observedLimitsAreSane): sys.exit("ERROR: observed limits deviate too much from expected limits at gluinoMass = {gM}, neutralinoMass={nM}".format(gM=gluinoMass, nM=neutralinoMass))
        limitsScanExpected.SetPoint(limitsScanExpected.GetN(), gluinoMass, neutralinoMass, expectedUpperLimit)
        limitsScanExpectedOneSigmaDown.SetPoint(limitsScanExpectedOneSigmaDown.GetN(), gluinoMass, neutralinoMass, expectedUpperLimitOneSigmaDown)
        limitsScanExpectedOneSigmaUp.SetPoint(limitsScanExpectedOneSigmaUp.GetN(), gluinoMass, neutralinoMass, expectedUpperLimitOneSigmaUp)
        crossSectionScanExpected.SetPoint(crossSectionScanExpected.GetN(), gluinoMass, neutralinoMass, expectedUpperLimit*crossSection)
        expectedCrossSectionLimits.append(((gluinoMass, neutralinoMass), expectedUpperLimit*crossSection))
        if ((minValue_crossSectionScanExpected == -1) or (expectedUpperLimit*crossSection < minValue_crossSectionScanExpected)): minValue_crossSectionScanExpected = expectedUpperLimit*crossSection
        if ((maxValue_crossSectionScanExpected == -1) or (expectedUpperLimit*crossSection > maxValue_crossSectionScanExpected)): maxValue_crossSectionScanExpected = expectedUpperLimit*crossSection


        limitsScanObserved.SetPoint(limitsScanObserved.GetN(), gluinoMass, neutralinoMass, observedUpperLimit)
        limitsScanObservedOneSigmaDown.SetPoint(limitsScanObservedOneSigmaDown.GetN(), gluinoMass, neutralinoMass, observedUpperLimitOneSigmaDown)
        limitsScanObservedOneSigmaUp.SetPoint(limitsScanObservedOneSigmaUp.GetN(), gluinoMass, neutralinoMass, observedUpperLimitOneSigmaUp)
        crossSectionScanObserved.SetPoint(crossSectionScanObserved.GetN(), gluinoMass, neutralinoMass, observedUpperLimit*crossSection)
        observedCrossSectionLimits.append(((gluinoMass, neutralinoMass), observedUpperLimit*crossSection))
        if ((minValue_crossSectionScanObserved == -1) or (observedUpperLimit*crossSection < minValue_crossSectionScanObserved)): minValue_crossSectionScanObserved = observedUpperLimit*crossSection
        if ((maxValue_crossSectionScanObserved == -1) or (observedUpperLimit*crossSection > maxValue_crossSectionScanObserved)): maxValue_crossSectionScanObserved = observedUpperLimit*crossSection

generatedMCTemplate.Close()

outputExpectedCrossSectionsFile=open("{oD}/expectedCrossSections_{s}.txt".format(oD=inputArguments.outputDirectory_rawOutput, s=inputArguments.outputSuffix), 'w')
outputExpectedCrossSectionsFile.write("{gMTitle:<19}{nMTitle:<19}{eXSTitle}\n".format(gMTitle="gluino mass", nMTitle="neutralino mass", eXSTitle="Expected limits on cross section (pb)"))
for expectedCrossSectionLimit in expectedCrossSectionLimits:
    outputExpectedCrossSectionsFile.write("{gM:<19.1f}{nM:<19.1f}{eXS:.3e}\n".format(gM=expectedCrossSectionLimit[0][0], nM=expectedCrossSectionLimit[0][1], eXS=expectedCrossSectionLimit[1]))
outputExpectedCrossSectionsFile.close()

outputObservedCrossSectionsFile=open("{oD}/observedCrossSections_{s}.txt".format(oD=inputArguments.outputDirectory_rawOutput, s=inputArguments.outputSuffix), 'w')
outputObservedCrossSectionsFile.write("{gMTitle:<19}{nMTitle:<19}{eXSTitle}\n".format(gMTitle="gluino mass", nMTitle="neutralino mass", eXSTitle="Observed limits on cross section (pb)"))
for observedCrossSectionLimit in observedCrossSectionLimits:
    outputObservedCrossSectionsFile.write("{gM:<19.1f}{nM:<19.1f}{oXS:.3e}\n".format(gM=observedCrossSectionLimit[0][0], nM=observedCrossSectionLimit[0][1], oXS=observedCrossSectionLimit[1]))
outputObservedCrossSectionsFile.close()

listOf2DScans = [limitsScanExpected, limitsScanExpectedOneSigmaDown, limitsScanExpectedOneSigmaUp, crossSectionScanExpected, limitsScanObserved, limitsScanObservedOneSigmaDown, limitsScanObservedOneSigmaUp, crossSectionScanObserved]
for scan2D in listOf2DScans:
    scan2D.SetNpx(160)
    scan2D.SetNpy(266)

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

observedLimitContours = limitsScanObserved.GetContourList(inputArguments.contour_signalStrength)
observedLimitContours.SetName("observedLimitContours")
observedLimitContoursOneSigmaDown = limitsScanObservedOneSigmaDown.GetContourList(inputArguments.contour_signalStrength)
observedLimitContoursOneSigmaDown.SetName("observedLimitContoursOneSigmaDown")
observedLimitContoursOneSigmaUp = limitsScanObservedOneSigmaUp.GetContourList(inputArguments.contour_signalStrength)
observedLimitContoursOneSigmaUp.SetName("observedLimitContoursOneSigmaUp")

formatContours(expectedLimitContours, 2, 5, ROOT.kRed)
formatContours(expectedLimitContoursOneSigmaDown, 2, 2, ROOT.kRed)
formatContours(expectedLimitContoursOneSigmaUp, 2, 2, ROOT.kRed)
formatContours(observedLimitContours, 1, 5, ROOT.kBlack)
formatContours(observedLimitContoursOneSigmaDown, 1, 2, ROOT.kBlack)
formatContours(observedLimitContoursOneSigmaUp, 1, 2, ROOT.kBlack)

CMS_lumi.writeExtraText = False
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
CMS_lumi.lumi_13TeV = "77.8 fb^{-1}"

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
histogramCrossSectionScanObserved.Draw("colz")
histogramCrossSectionScanObserved.GetXaxis().SetRangeUser(1000., 1750.)
histogramCrossSectionScanObserved.GetZaxis().SetRangeUser(0.9*minValue_crossSectionScanObserved, 1.1*maxValue_crossSectionScanObserved)
for contoursList in [expectedLimitContours, expectedLimitContoursOneSigmaDown, expectedLimitContoursOneSigmaUp, observedLimitContours, observedLimitContoursOneSigmaDown, observedLimitContoursOneSigmaUp]:
    contoursList.Draw("SAME")
CMS_lumi.CMS_lumi(canvas, 4, 0)
ROOT.gPad.Update()
ROOT.gPad.RedrawAxis()
frame = ROOT.gPad.GetFrame()
frame.Draw()
canvas.Update()
canvas.SaveAs("{oD}/{s}_observedLimits.png".format(oD=inputArguments.outputDirectory_plots, s=inputArguments.outputSuffix))

outputFileName = "{oD}/limits_{suffix}.root".format(oD=inputArguments.outputDirectory_rawOutput, suffix=inputArguments.outputSuffix)
outputFile=ROOT.TFile(outputFileName, "RECREATE")
tObjectsToSave = [limitsScanExpected, limitsScanExpectedOneSigmaUp, limitsScanExpectedOneSigmaDown, crossSectionScanExpected, histogramExpectedLimits, histogramExpectedLimitsOneSigmaDown, histogramExpectedLimitsOneSigmaUp, histogramCrossSectionScanExpected, expectedLimitContours, expectedLimitContoursOneSigmaDown, expectedLimitContoursOneSigmaUp, limitsScanObserved, limitsScanObservedOneSigmaUp, limitsScanObservedOneSigmaDown, crossSectionScanObserved, histogramObservedLimits, histogramObservedLimitsOneSigmaDown, histogramObservedLimitsOneSigmaUp, histogramCrossSectionScanObserved, observedLimitContours, observedLimitContoursOneSigmaDown, observedLimitContoursOneSigmaUp]
for tObject in tObjectsToSave:
    outputFile.WriteTObject(tObject)

outputFile.Close()
print("All done!")
