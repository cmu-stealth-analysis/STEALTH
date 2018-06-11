#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, tmROOTUtils, pdb

inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination from input MC in real data.')
inputArgumentsParser.add_argument('--MCTemplate', default="plot_susyMasses_template.root", help='Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.', type=str)
inputArgumentsParser.add_argument('--combineOutputsDirectory', default="analysis/combineToolOutputs", help='Path to input MC file.',type=str)
inputArgumentsParser.add_argument('--plotObservedLimits', action='store_true', help="Get observed limits in addition to expected limits.")
inputArguments = inputArgumentsParser.parse_args()

limitsScan=ROOT.TGraph2D()
limitsScan.SetName("limitsScan")
limitsScanOneSigmaUp=ROOT.TGraph2D()
limitsScanOneSigmaUp.SetName("limitsScanOneSigmaUp")
limitsScanOneSigmaDown=ROOT.TGraph2D()
limitsScanOneSigmaDown.SetName("limitsScanOneSigmaDown")
limitsScanObserved=ROOT.TGraph2D()
limitsScanObserved.SetName("limitsScanObserved")

for scan2D in [limitsScan, limitsScanOneSigmaUp, limitsScanOneSigmaDown]:
    scan2D.SetTitle("Expected Limits;m_{#tilde{#it{g}}};m_{#tilde{#it{#chi_{1}^{0}}}}")
limitsScanObserved.SetTitle("Observed Limits;m_{#tilde{#it{g}}};m_{#tilde{#it{#chi_{1}^{0}}}}")

generatedMCTemplate = ROOT.TFile(inputArguments.MCTemplate)
h_MCTemplate = generatedMCTemplate.Get("h_susyMasses_template")
for gluinoMassBin in range(1, 1+h_MCTemplate.GetXaxis().GetNbins()):
    for neutralinoMassBin in range(1, 1+h_MCTemplate.GetYaxis().GetNbins()):
        if not(int(0.5 + h_MCTemplate.GetBinContent(gluinoMassBin, neutralinoMassBin)) == 1): continue
        gluinoMass = int(0.5 + h_MCTemplate.GetXaxis().GetBinCenter(gluinoMassBin))
        neutralinoMass = int(0.5 + h_MCTemplate.GetYaxis().GetBinCenter(neutralinoMassBin))
        print("Analyzing bin at gluino mass={gM}, neutralino mass={nM}".format(gM=gluinoMass, nM=neutralinoMass))
        combineOutputFile=ROOT.TFile("{combineOutputsDirectory}/higgsCombine_gluinoMass_{gluinoMass}_neutralinoMass_{neutralinoMass}.AsymptoticLimits.mH120.root".format(combineOutputsDirectory=inputArguments.combineOutputsDirectory, gluinoMass=gluinoMass, neutralinoMass=neutralinoMass), "READ")
        limitObject = combineOutputFile.Get("limit")
        if not(limitObject):
            print("WARNING: limits not available at gluinoMass = {gM}, neutralinoMass={nM}".format(gM=gluinoMass, nM=neutralinoMass))
            combineOutputFile.Close()
            continue
	limitObject.GetEntry(2)
	expectedUpperLimit= limitObject.limit
	limitObject.GetEntry(1)
	expectedUpperLimitOneSigmaUp=limitObject.limit
	limitObject.GetEntry(3)
    	expectedUpperLimitOneSigmaDown=limitObject.limit
	limitObject.GetEntry(5)
	observedUpperLimit=limitObject.limit
	limitsScan.SetPoint(limitsScan.GetN(), gluinoMass, neutralinoMass, expectedUpperLimit)
        limitsScanOneSigmaUp.SetPoint(limitsScanOneSigmaUp.GetN(), gluinoMass, neutralinoMass, expectedUpperLimitOneSigmaUp)
        limitsScanOneSigmaDown.SetPoint(limitsScanOneSigmaDown.GetN(), gluinoMass, neutralinoMass, expectedUpperLimitOneSigmaDown)
        limitsScanObserved.SetPoint(limitsScanObserved.GetN(), gluinoMass, neutralinoMass, observedUpperLimit)
        combineOutputFile.Close()

for scan2D in [limitsScan, limitsScanOneSigmaUp, limitsScanOneSigmaDown, limitsScanObserved]:
    scan2D.SetNpx(128)
    scan2D.SetNpy(160)

histogramExpectedLimits=limitsScan.GetHistogram()
histogramExpectedLimitsOneSigmaUp=limitsScanOneSigmaUp.GetHistogram()
histogramExpectedLimitsOneSigmaDown=limitsScanOneSigmaDown.GetHistogram()
histogramObservedLimits=limitsScanObserved.GetHistogram()

ExpectedLimits = ROOT.TGraph()
ExpectedLimits.SetName("ExpectedLimits")
ExpectedLimits = limitsScan.GetContourList(1.0)

ExpectedLimitsOneSigmaUp = limitsScanOneSigmaUp.GetContourList(1.0)
ExpectedLimitsOneSigmaUpListIteratorNext = ROOT.TIter(ExpectedLimitsOneSigmaUp)
upCounter = 0
while True:
    upCounter += 1
    ExpectedLimitsOneSigmaUpContour = ExpectedLimitsOneSigmaUpListIteratorNext()
    if not(ExpectedLimitsOneSigmaUpContour): break
    ExpectedLimitsOneSigmaUpContour.SetName("ExpectedLimitsOneSigmaDown_{upC}".format(upC=upCounter))
    ExpectedLimitsOneSigmaUpContour.SetLineStyle(5)

ExpectedLimitsOneSigmaDown = limitsScanOneSigmaDown.GetContourList(1.0)
ExpectedLimitsOneSigmaDownContourListIteratorNext = ROOT.TIter(ExpectedLimitsOneSigmaDown)
downCounter = 0
while True:
    downCounter += 1
    ExpectedLimitsOneSigmaDownContour = ExpectedLimitsOneSigmaDownContourListIteratorNext()
    if not(ExpectedLimitsOneSigmaDownContour): break
    ExpectedLimitsOneSigmaDownContour.SetName("ExpectedLimitsOneSigmaDown_{downC}".format(downC=downCounter))
    ExpectedLimitsOneSigmaDownContour.SetLineStyle(5)

ObservedLimits = ROOT.TGraph()
ObservedLimits.SetName("ObservedLimits")
ObservedLimits = limitsScanObserved.GetContourList(1.0)

outputFile=ROOT.TFile("analysis/limitPlots/expectedLimitPlots_savedObjects.root", "RECREATE")
c_expectedLimits = ROOT.TCanvas()
if (inputArguments.plotObservedLimits): c_expectedLimits = tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[histogramExpectedLimits, ObservedLimits, ExpectedLimits, ExpectedLimitsOneSigmaUp, ExpectedLimitsOneSigmaDown], canvasName="c_expectedLimits", outputROOTFile=outputFile, outputDocumentName="analysis/limitPlots/expectedLimitPlots", customOptStat=0, customPlotOptions_firstObject="colz", enableLogZ = True)
else: c_expectedLimits = tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[histogramExpectedLimits, ExpectedLimits, ExpectedLimitsOneSigmaUp, ExpectedLimitsOneSigmaDown], canvasName="c_expectedLimits", outputROOTFile=outputFile, outputDocumentName="analysis/limitPlots/expectedLimitPlots", customOptStat=0, customPlotOptions_firstObject="colz", enableLogZ = True)
outputFile.Close()
