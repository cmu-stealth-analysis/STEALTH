#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, tmROOTUtils, pdb, tmGeneralUtils, sys

inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination from input MC in real data.')
inputArgumentsParser.add_argument('--crossSectionsFile', default="SusyCrossSections13TevGluGlu.txt", help='Path to dat file that contains cross-sections as a function of gluino mass, to use while weighting events.',type=str)
inputArgumentsParser.add_argument('--MCTemplate', default="plot_susyMasses_template.root", help='Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.', type=str)
inputArgumentsParser.add_argument('--combineResultsDirectory', default="root://cmseos.fnal.gov//store/user/lpcsusystealth/combineToolOutputs", help='EOS path at which combine tool results can be found.',type=str)
inputArgumentsParser.add_argument('--combineOutputPrefix', default="", help='Prefix of Higgs combine results.',type=str)
inputArgumentsParser.add_argument('--outputSuffix', default="", help='Suffix to append to all results.',type=str)
inputArgumentsParser.add_argument('--plotObservedLimits', action='store_true', help="Get observed limits in addition to expected limits.")
inputArgumentsParser.add_argument('--minGluinoMass', default=-1., help='Minimum gluino mass on which to run.', type=float)
inputArgumentsParser.add_argument('--maxGluinoMass', default=1775., help='Max gluino mass for the 2D plots.',type=float)
inputArgumentsParser.add_argument('--contour_signalStrength', default=1., help='Signal strength at which to draw the contours.',type=float)
inputArguments = inputArgumentsParser.parse_args()

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
for scan2D in [limitsScanExpected, limitsScanExpectedOneSigmaDown, limitsScanExpectedOneSigmaUp]:
    scan2D.SetTitle("Expected value of r;m_{#tilde{#it{g}}};m_{#tilde{#it{#chi_{1}^{0}}}}")
crossSectionScanExpected.SetTitle("Expected limits on cross section (pb);m_{#tilde{#it{g}}};m_{#tilde{#it{#chi_{1}^{0}}}}")

limitsScanObserved=ROOT.TGraph2D()
limitsScanObservedOneSigmaDown=ROOT.TGraph2D()
limitsScanObservedOneSigmaUp=ROOT.TGraph2D()
crossSectionScanObserved = ROOT.TGraph2D()
if inputArguments.plotObservedLimits:
    limitsScanObserved.SetName("limitsScanObserved")
    limitsScanObservedOneSigmaDown=ROOT.TGraph2D()
    limitsScanObservedOneSigmaDown.SetName("limitsScanObservedOneSigmaDown")
    limitsScanObservedOneSigmaUp=ROOT.TGraph2D()
    limitsScanObservedOneSigmaUp.SetName("limitsScanObservedOneSigmaUp")
    crossSectionScanObserved.SetName("crossSectionScanObserved")
    for scan2D in [limitsScanObserved, limitsScanObservedOneSigmaDown, limitsScanObservedOneSigmaUp]:
        scan2D.SetTitle("Observed value of r;m_{#tilde{#it{g}}};m_{#tilde{#it{#chi_{1}^{0}}}}")
    crossSectionScanObserved.SetTitle("Observed limits on cross section (pb);m_{#tilde{#it{g}}};m_{#tilde{#it{#chi_{1}^{0}}}}")

expectedCrossSectionLimits = []
observedCrossSectionLimits = []

maxValue_crossSectionScanExpected = -1
minValue_crossSectionScanExpected = -1
maxValue_crossSectionScanObserved = -1
minValue_crossSectionScanObserved = -1

generatedMCTemplate = ROOT.TFile(inputArguments.MCTemplate, "READ")
h_MCTemplate = generatedMCTemplate.Get("h_susyMasses_template")
for gluinoMassBin in range(1, 1+h_MCTemplate.GetXaxis().GetNbins()):
    gluinoMass = int(0.5 + h_MCTemplate.GetXaxis().GetBinCenter(gluinoMassBin))
    if (inputArguments.minGluinoMass > 0 and gluinoMass < inputArguments.minGluinoMass): continue
    for neutralinoMassBin in range(1, 1+h_MCTemplate.GetYaxis().GetNbins()):
        if not(int(0.5 + h_MCTemplate.GetBinContent(gluinoMassBin, neutralinoMassBin)) == 1): continue
        neutralinoMass = int(0.5 + h_MCTemplate.GetYaxis().GetBinCenter(neutralinoMassBin))
        crossSection = crossSectionsDictionary[gluinoMass]
        print("Analyzing bin at gluino mass={gM}, neutralino mass={nM}".format(gM=gluinoMass, nM=neutralinoMass))

        # nominal
        combineOutputFile=ROOT.TFile.Open("{cRD}/higgsCombine_{cOP}_gluinoMassBin{gMB}_neutralinoMassBin{nMB}.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin), "READ")
        if ((combineOutputFile.IsZombie() == ROOT.kTRUE) or not(combineOutputFile.IsOpen() == ROOT.kTRUE)):
            sys.exit("Error in opening file: {cRD}/higgsCombine_{cOP}_gluinoMassBin{gMB}_neutralinoMassBin{nMB}.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin))
        limitObject = combineOutputFile.Get("limit")
        if not(limitObject):
            combineOutputFile.Close()
            sys.exit("ERROR: nominal limits not available at gluinoMass = {gM}, neutralinoMass={nM}".format(gM=gluinoMass, nM=neutralinoMass))
        limitObject.GetEntry(2)
	expectedUpperLimit= limitObject.limit
	limitObject.GetEntry(1)
	expectedUpperLimitOneSigmaUp=limitObject.limit
	limitObject.GetEntry(3)
    	expectedUpperLimitOneSigmaDown=limitObject.limit
	limitObject.GetEntry(5)
	observedUpperLimit=limitObject.limit
        observedUpperLimitOneSigmaDown=0.
        observedUpperLimitOneSigmaUp=0.

        if inputArguments.plotObservedLimits:
            # cross section down
            combineOutputFile_crossSectionsDown=ROOT.TFile.Open("{cRD}/higgsCombine_{cOP}_gluinoMassBin{gMB}_neutralinoMassBin{nMB}_crossSectionsDown.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin), "READ")
            limitObject = combineOutputFile_crossSectionsDown.Get("limit")
            if not(limitObject):
                combineOutputFile_crossSectionsDown.Close()
                sys.exit("ERROR: limits with cross sections down not available at gluinoMass = {gM}, neutralinoMass={nM}".format(gM=gluinoMass, nM=neutralinoMass))
            limitObject.GetEntry(5)
	    observedUpperLimitOneSigmaDown=limitObject.limit

            # cross section up
            combineOutputFile_crossSectionsUp=ROOT.TFile.Open("{cRD}/higgsCombine_{cOP}_gluinoMassBin{gMB}_neutralinoMassBin{nMB}_crossSectionsUp.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin), "READ")
            limitObject = combineOutputFile_crossSectionsUp.Get("limit")
            if not(limitObject):
                combineOutputFile_crossSectionsUp.Close()
                sys.exit("ERROR: limits with cross sections up not available at gluinoMass = {gM}, neutralinoMass={nM}".format(gM=gluinoMass, nM=neutralinoMass))
            limitObject.GetEntry(5)
	    observedUpperLimitOneSigmaUp=limitObject.limit

        limitsScanExpected.SetPoint(limitsScanExpected.GetN(), gluinoMass, neutralinoMass, expectedUpperLimit)
        limitsScanExpectedOneSigmaDown.SetPoint(limitsScanExpectedOneSigmaDown.GetN(), gluinoMass, neutralinoMass, expectedUpperLimitOneSigmaDown)
        limitsScanExpectedOneSigmaUp.SetPoint(limitsScanExpectedOneSigmaUp.GetN(), gluinoMass, neutralinoMass, expectedUpperLimitOneSigmaUp)
        crossSectionScanExpected.SetPoint(crossSectionScanExpected.GetN(), gluinoMass, neutralinoMass, expectedUpperLimit*crossSection)
        expectedCrossSectionLimits.append(((gluinoMass, neutralinoMass), expectedUpperLimit*crossSection))
        if ((minValue_crossSectionScanExpected == -1) or (expectedUpperLimit*crossSection < minValue_crossSectionScanExpected)): minValue_crossSectionScanExpected = expectedUpperLimit*crossSection
        if ((maxValue_crossSectionScanExpected == -1) or (expectedUpperLimit*crossSection > maxValue_crossSectionScanExpected)): maxValue_crossSectionScanExpected = expectedUpperLimit*crossSection

        if inputArguments.plotObservedLimits:
            limitsScanObserved.SetPoint(limitsScanObserved.GetN(), gluinoMass, neutralinoMass, observedUpperLimit)
            limitsScanObservedOneSigmaDown.SetPoint(limitsScanObservedOneSigmaDown.GetN(), gluinoMass, neutralinoMass, observedUpperLimitOneSigmaDown)
            limitsScanObservedOneSigmaUp.SetPoint(limitsScanObservedOneSigmaUp.GetN(), gluinoMass, neutralinoMass, observedUpperLimitOneSigmaUp)
            crossSectionScanObserved.SetPoint(crossSectionScanObserved.GetN(), gluinoMass, neutralinoMass, observedUpperLimit*crossSection)
            observedCrossSectionLimits.append(((gluinoMass, neutralinoMass), observedUpperLimit*crossSection))
            if ((minValue_crossSectionScanObserved == -1) or (observedUpperLimit*crossSection < minValue_crossSectionScanObserved)): minValue_crossSectionScanObserved = observedUpperLimit*crossSection
            if ((maxValue_crossSectionScanObserved == -1) or (observedUpperLimit*crossSection > maxValue_crossSectionScanObserved)): maxValue_crossSectionScanObserved = observedUpperLimit*crossSection

        combineOutputFile.Close()
        if inputArguments.plotObservedLimits:    
            combineOutputFile_crossSectionsDown.Close()
            combineOutputFile_crossSectionsUp.Close()
generatedMCTemplate.Close()

outputExpectedCrossSectionsFile=open("analysis/limitPlots/expectedCrossSections_{s}.txt".format(s=inputArguments.outputSuffix), 'w')
outputExpectedCrossSectionsFile.write("{gMTitle:<19}{nMTitle:<19}{eXSTitle}\n".format(gMTitle="gluino mass", nMTitle="neutralino mass", eXSTitle="Expected limits on cross section (pb)"))
for expectedCrossSectionLimit in expectedCrossSectionLimits:
    outputExpectedCrossSectionsFile.write("{gM:<19.1f}{nM:<19.1f}{eXS:.3e}\n".format(gM=expectedCrossSectionLimit[0][0], nM=expectedCrossSectionLimit[0][1], eXS=expectedCrossSectionLimit[1]))
outputExpectedCrossSectionsFile.close()

if (inputArguments.plotObservedLimits):
    outputObservedCrossSectionsFile=open("analysis/limitPlots/observedCrossSections_{s}.txt".format(s=inputArguments.outputSuffix), 'w')
    outputObservedCrossSectionsFile.write("{gMTitle:<19}{nMTitle:<19}{eXSTitle}\n".format(gMTitle="gluino mass", nMTitle="neutralino mass", eXSTitle="Observed limits on cross section (pb)"))
    for observedCrossSectionLimit in observedCrossSectionLimits:
        outputObservedCrossSectionsFile.write("{gM:<19.1f}{nM:<19.1f}{oXS:.3e}\n".format(gM=observedCrossSectionLimit[0][0], nM=observedCrossSectionLimit[0][1], oXS=observedCrossSectionLimit[1]))
    outputObservedCrossSectionsFile.close()


for scan2D in [limitsScanExpected, limitsScanExpectedOneSigmaDown, limitsScanExpectedOneSigmaUp, crossSectionScanExpected, limitsScanObserved, limitsScanObservedOneSigmaDown, limitsScanObservedOneSigmaUp, crossSectionScanObserved]:
    scan2D.SetNpx(160)
    scan2D.SetNpy(266)

histogramExpectedLimits = limitsScanExpected.GetHistogram()
histogramExpectedLimitsOneSigmaUp = limitsScanExpectedOneSigmaUp.GetHistogram()
histogramExpectedLimitsOneSigmaDown = limitsScanExpectedOneSigmaDown.GetHistogram()
histogramCrossSectionScanExpected = crossSectionScanExpected.GetHistogram()

histogramObservedLimits = ROOT.TH2F()
histogramObservedLimitsOneSigmaDown = ROOT.TH2F()
histogramObservedLimitsOneSigmaUp = ROOT.TH2F()
histogramCrossSectionScanObserved = ROOT.TH2F()
if inputArguments.plotObservedLimits:
    histogramObservedLimits = limitsScanObserved.GetHistogram()
    histogramObservedLimitsOneSigmaDown = limitsScanObservedOneSigmaDown.GetHistogram()
    histogramObservedLimitsOneSigmaUp = limitsScanObservedOneSigmaUp.GetHistogram()
    histogramCrossSectionScanObserved = crossSectionScanObserved.GetHistogram()

# for hist2D in [histogramExpectedLimits, histogramExpectedLimitsOneSigmaDown, histogramExpectedLimitsOneSigmaUp, histogramObservedLimits, histogramObservedLimitsOneSigmaDown, histogramObservedLimitsOneSigmaUp]:
#     hist2D.GetZaxis().SetRangeUser(0.008, 110.)

# for hist2D in [histogramCrossSectionScanExpected, histogramCrossSectionScanObserved]:
#     hist2D.GetZaxis().SetRangeUser(0.0005, 5.)
histogramCrossSectionScanExpected.GetZaxis().SetRangeUser(0.9*minValue_crossSectionScanExpected, 1.1*maxValue_crossSectionScanExpected)
histogramCrossSectionScanObserved.GetZaxis().SetRangeUser(0.9*minValue_crossSectionScanObserved, 1.1*maxValue_crossSectionScanObserved)

ExpectedLimits = ROOT.TGraph()
ExpectedLimitsOneSigmaDown = ROOT.TGraph()
ExpectedLimitsOneSigmaUp = ROOT.TGraph()
ExpectedLimits.SetName("ExpectedLimits")
ExpectedLimits = limitsScanExpected.GetContourList(inputArguments.contour_signalStrength)
ExpectedLimitsContourListIteratorNext = ROOT.TIter(ExpectedLimits)
counter = 0
while True:
    counter += 1
    ExpectedLimitsContour = ExpectedLimitsContourListIteratorNext()
    if not(ExpectedLimitsContour): break
    ExpectedLimitsContour.SetName("ExpectedLimits_{c}".format(c=counter))
    ExpectedLimitsContour.SetLineStyle(2)
    ExpectedLimitsContour.SetLineWidth(5)
    ExpectedLimitsContour.SetLineColor(ROOT.kRed)
ExpectedLimitsOneSigmaDown = limitsScanExpectedOneSigmaDown.GetContourList(inputArguments.contour_signalStrength)
ExpectedLimitsOneSigmaDownContourListIteratorNext = ROOT.TIter(ExpectedLimitsOneSigmaDown)
downCounter = 0
while True:
    downCounter += 1
    ExpectedLimitsOneSigmaDownContour = ExpectedLimitsOneSigmaDownContourListIteratorNext()
    if not(ExpectedLimitsOneSigmaDownContour): break
    ExpectedLimitsOneSigmaDownContour.SetName("ExpectedLimitsOneSigmaDown_{downC}".format(downC=downCounter))
    ExpectedLimitsOneSigmaDownContour.SetLineStyle(2)
    ExpectedLimitsOneSigmaDownContour.SetLineWidth(2)
    ExpectedLimitsOneSigmaDownContour.SetLineColor(ROOT.kRed)
ExpectedLimitsOneSigmaUp = limitsScanExpectedOneSigmaUp.GetContourList(inputArguments.contour_signalStrength)
ExpectedLimitsOneSigmaUpListIteratorNext = ROOT.TIter(ExpectedLimitsOneSigmaUp)
upCounter = 0
while True:
    upCounter += 1
    ExpectedLimitsOneSigmaUpContour = ExpectedLimitsOneSigmaUpListIteratorNext()
    if not(ExpectedLimitsOneSigmaUpContour): break
    ExpectedLimitsOneSigmaUpContour.SetName("ExpectedLimitsOneSigmaDown_{upC}".format(upC=upCounter))
    ExpectedLimitsOneSigmaUpContour.SetLineStyle(2)
    ExpectedLimitsOneSigmaUpContour.SetLineWidth(2)
    ExpectedLimitsOneSigmaUpContour.SetLineColor(ROOT.kRed)

ObservedLimits = ROOT.TGraph()
ObservedLimitsOneSigmaDown = ROOT.TGraph()
ObservedLimitsOneSigmaUp = ROOT.TGraph()
if inputArguments.plotObservedLimits:
    ObservedLimits.SetName("ObservedLimits")
    ObservedLimits = limitsScanObserved.GetContourList(inputArguments.contour_signalStrength)
    ObservedLimitsContourListIteratorNext = ROOT.TIter(ObservedLimits)
    counter = 0
    while True:
        counter += 1
        ObservedLimitsContour = ObservedLimitsContourListIteratorNext()
        if not(ObservedLimitsContour): break
        ObservedLimitsContour.SetName("ObservedLimits_{c}".format(c=counter))
        ObservedLimitsContour.SetLineStyle(1)
        ObservedLimitsContour.SetLineWidth(5)
        ObservedLimitsContour.SetLineColor(ROOT.kBlack)
    ObservedLimitsOneSigmaDown = limitsScanObservedOneSigmaDown.GetContourList(inputArguments.contour_signalStrength)
    ObservedLimitsOneSigmaDownContourListIteratorNext = ROOT.TIter(ObservedLimitsOneSigmaDown)
    downCounter = 0
    while True:
        downCounter += 1
        ObservedLimitsOneSigmaDownContour = ObservedLimitsOneSigmaDownContourListIteratorNext()
        if not(ObservedLimitsOneSigmaDownContour): break
        ObservedLimitsOneSigmaDownContour.SetName("ObservedLimitsOneSigmaDown_{downC}".format(downC=downCounter))
        ObservedLimitsOneSigmaDownContour.SetLineStyle(1)
        ObservedLimitsOneSigmaDownContour.SetLineWidth(2)
        ObservedLimitsOneSigmaDownContour.SetLineColor(ROOT.kBlack)
    ObservedLimitsOneSigmaUp = limitsScanObservedOneSigmaUp.GetContourList(inputArguments.contour_signalStrength)
    ObservedLimitsOneSigmaUpListIteratorNext = ROOT.TIter(ObservedLimitsOneSigmaUp)
    upCounter = 0
    while True:
        upCounter += 1
        ObservedLimitsOneSigmaUpContour = ObservedLimitsOneSigmaUpListIteratorNext()
        if not(ObservedLimitsOneSigmaUpContour): break
        ObservedLimitsOneSigmaUpContour.SetName("ObservedLimitsOneSigmaDown_{upC}".format(upC=upCounter))
        ObservedLimitsOneSigmaUpContour.SetLineStyle(1)
        ObservedLimitsOneSigmaUpContour.SetLineWidth(2)
        ObservedLimitsOneSigmaUpContour.SetLineColor(ROOT.kBlack)

outputFileName = "analysis/limitPlots/expectedLimitPlots_{suffix}_savedObjects.root".format(suffix=inputArguments.outputSuffix)
if (inputArguments.plotObservedLimits): outputFileName = "analysis/limitPlots/fullLimitPlots_{suffix}_savedObjects.root".format(suffix=inputArguments.outputSuffix)
outputFile=ROOT.TFile(outputFileName, "RECREATE")
if (inputArguments.plotObservedLimits):
    tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[histogramCrossSectionScanObserved, ObservedLimits, ObservedLimitsOneSigmaDown, ObservedLimitsOneSigmaUp, ExpectedLimits, ExpectedLimitsOneSigmaDown, ExpectedLimitsOneSigmaUp], canvasName="c_fullLimits", outputROOTFile=outputFile, outputDocumentName="analysis/limitPlots/fullLimitPlots_{suffix}".format(suffix=inputArguments.outputSuffix), customOptStat=0, customPlotOptions_firstObject="colz", enableLogZ = True, customXRange=[inputArguments.minGluinoMass, inputArguments.maxGluinoMass])
else:
    tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[histogramCrossSectionScanExpected, ExpectedLimits, ExpectedLimitsOneSigmaDown, ExpectedLimitsOneSigmaUp], canvasName="c_expectedLimits", outputROOTFile=outputFile, outputDocumentName="analysis/limitPlots/expectedLimitPlots_{suffix}".format(suffix=inputArguments.outputSuffix), customOptStat=0, customPlotOptions_firstObject="colz", enableLogZ = True, customXRange=[inputArguments.minGluinoMass, inputArguments.maxGluinoMass])
outputFile.Close()
