#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, tmROOTUtils, pdb

from tmProgressBar import tmProgressBar
from tmGeneralUtils import prettyPrintDictionary

inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination from input MC in real data.')
inputArgumentsParser.add_argument('--crossSectionsFile', default="SusyCrossSections13TevGluGlu.txt", help='Path to dat file that contains cross-sections as a function of gluino mass, to use while weighting events.',type=str)
inputArgumentsParser.add_argument('--MCTemplate', default="plot_susyMasses_template.root", help='Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.', type=str)
inputArgumentsParser.add_argument('--combineOutputsDirectory', default="analysis/combineToolOutputs", help='Path to input MC file.',type=str)
inputArgumentsParser.add_argument('--plotObservedLimits', action='store_true', help="Get observed limits in addition to expected limits.")
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

print("Read in cross-sections as a function of gluino mass:")
prettyPrintDictionary(crossSectionsDictionary)

limitsScan=ROOT.TGraph2D()
limitsScan.SetName("limitsScan")
crossSectionWeightedLimitsScan=ROOT.TGraph2D()
crossSectionWeightedLimitsScan.SetName("crossSectionWeightedLimitsScan")
crossSectionWeightedExpectedLimitsScan=ROOT.TGraph2D()
crossSectionWeightedExpectedLimitsScan.SetName("crossSectionWeightedExpectedLimitsScan")

limitsScanOneSigmaUp=ROOT.TGraph2D()
limitsScanOneSigmaUp.SetName("limitsScanOneSigmaUp")
limitsScanOneSigmaDown=ROOT.TGraph2D()
limitsScanOneSigmaDown.SetName("limitsScanOneSigmaDown")
limitsScanObserved=ROOT.TGraph2D()
limitsScanObserved.SetName("limitsScanObserved")
limitsScanObservedOneSigmaUp=ROOT.TGraph2D()
limitsScanObservedOneSigmaUp.SetName("limitsScanObservedOneSigmaUp")
limitsScanObservedOneSigmaDown=ROOT.TGraph2D()
limitsScanObservedOneSigmaDown.SetName("limitsScanObservedOneSigmaDown")

generatedMCTemplate = ROOT.TFile(inputArguments.MCTemplate)
h_MCTemplate = generatedMCTemplate.Get("h_susyMasses_template")
for gluinoMassBin in range(1, 1+h_MCTemplate.GetXaxis().GetNbins()):
    for neutralinoMassBin in range(1, 1+h_MCTemplate.GetYaxis().GetNbins()):
        if not(int(0.5 + h_MCTemplate.GetBinContent(gluinoMassBin, neutralinoMassBin)) == 1): continue
        gluinoMass = int(0.5 + h_MCTemplate.GetXaxis().GetBinCenter(gluinoMassBin))
        neutralinoMass = int(0.5 + h_MCTemplate.GetYaxis().GetBinCenter(neutralinoMassBin))
        print("Analyzing bin at gluino mass={gM}, neutralino mass={nM}".format(gM=gluinoMass, nM=neutralinoMass))
        combineOutputFile=ROOT.TFile("{combineOutputsDirectory}/higgsCombine_gluinoMass_{gluinoMass}_neutralinoMass_{neutralinoMass}.AsymptoticLimits.mH120.root".format(combineOutputsDirectory=inputArguments.combineOutputsDirectory, gluinoMass=gluinoMass, neutralinoMass=neutralinoMass), "READ")
        # print("combineOutputFile contents:")
        # combineOutputFile.ls()
	limitObject = combineOutputFile.Get("limit")
        if not(limitObject):
            print("WARNING: limits not available at gluinoMass = {gM}, neutralinoMass={nM}".format(gM=gluinoMass, nM=neutralinoMass))
            combineOutputFile.Close()
            continue
	limitObject.GetEntry(2);
	expectedUpperLimit= limitObject.limit
	crossSectionWeightedExpectedUpperLimit= expectedUpperLimit * float(crossSectionsDictionary[gluinoMass])
	limitObject.GetEntry(5)
	observedUpperLimit=limitObject.limit
	crossSectionWeightedObservedUpperLimit= observedUpperLimit * float(crossSectionsDictionary[gluinoMass])
	limitObject.GetEntry(1)
	expectedUpperLimitOneSigmaUp=limitObject.limit
	limitObject.GetEntry(3)
    	expectedUpperLimitOneSigmaDown=limitObject.limit
	limitsScan.SetPoint(limitsScan.GetN(), gluinoMass, neutralinoMass, expectedUpperLimit)
        limitsScanOneSigmaUp.SetPoint(limitsScanOneSigmaUp.GetN(), gluinoMass, neutralinoMass, expectedUpperLimitOneSigmaUp)
        limitsScanOneSigmaDown.SetPoint(limitsScanOneSigmaDown.GetN(), gluinoMass, neutralinoMass, expectedUpperLimitOneSigmaDown)
        limitsScanObserved.SetPoint(limitsScanObserved.GetN(), gluinoMass, neutralinoMass, observedUpperLimit)
        # limitsScanObservedOneSigmaUp.SetPoint(limitsScanObservedOneSigmaUp.GetN(), gluinoMass, neutralinoMass, observedUpperLimit)
        # limitsScanObservedOneSigmaDown.SetPoint(limitsScanObservedOneSigmaDown.GetN(), gluinoMass, neutralinoMass, observedUpperLimit)
        crossSectionWeightedLimitsScan.SetPoint(crossSectionWeightedLimitsScan.GetN(), gluinoMass, neutralinoMass, crossSectionWeightedObservedUpperLimit)
        crossSectionWeightedExpectedLimitsScan.SetPoint(crossSectionWeightedExpectedLimitsScan.GetN(), gluinoMass, neutralinoMass, crossSectionWeightedExpectedUpperLimit)
        combineOutputFile.Close()

for scan2D in [limitsScan, limitsScanOneSigmaUp, limitsScanOneSigmaDown, limitsScanObserved, limitsScanObservedOneSigmaUp, limitsScanObservedOneSigmaDown, crossSectionWeightedLimitsScan, crossSectionWeightedExpectedLimitsScan]:
    scan2D.SetNpx(128)
    scan2D.SetNpy(160)

# limitsScan.SetNpx(128)
# limitsScan.SetNpy(160)
# limitsScanOneSigmaUp.SetNpx(128)
# limitsScanOneSigmaUp.SetNpx(160)
# limitsScanOneSigmaDown.SetNpx(128)
# limitsScanOneSigmaDown.SetNpx(160)
# limitsScanObserved.SetNpx(128)
# limitsScanObserved.SetNpy(160)
# limitsScanObservedOneSigmaUp.SetNpx(128)
# limitsScanObservedOneSigmaUp.SetNpy(160)
# limitsScanObservedOneSigmaDown.SetNpx(128)
# limitsScanObservedOneSigmaDown.SetNpy(160)
# crossSectionWeightedLimitsScan.SetNpx(128)
# crossSectionWeightedExpectedLimitsScan.SetNpy(160)

histogramExpectedLimits=limitsScan.GetHistogram()
histogramExpectedLimitsOneSigmaUp=limitsScanOneSigmaUp.GetHistogram()
histogramExpectedLimitsOneSigmaDown=limitsScanOneSigmaDown.GetHistogram()
histogramObservedLimits=limitsScanObserved.GetHistogram()
# histogramObservedLimitsOneSigmaUp=limitsScanObservedOneSigmaUp.GetHistogram()
# histogramObservedLimitsOneSigmaDown=limitsScanObservedOneSigmaDown.GetHistogram()
MassScan2D=crossSectionWeightedLimitsScan.GetHistogram()
ExpectedLimitsMassScan2D=crossSectionWeightedExpectedLimitsScan.GetHistogram()
outputCanvas=ROOT.TCanvas("outputCanvas","outputCanvas",1024, 768);
limitsScan.Draw("colz")
limitsScanOneSigmaUp.Draw("colz")
limitsScanOneSigmaDown.Draw("colz")
if (inputArguments.plotObservedLimits): limitsScanObserved.Draw("colz")
# limitsScanObservedOneSigmaUp.Draw("colz")
# limitsScanObservedOneSigmaDown.Draw("colz")
ExpectedLimits = ROOT.TGraph()
ExpectedLimits.SetName("ExpectedLimits")

ExpectedLimits = limitsScan.GetContourList(1.0);
ExpectedLimitsOneSigmaUp = limitsScanOneSigmaUp.GetContourList(1.0);
ExpectedLimitsOneSigmaDown = limitsScanOneSigmaDown.GetContourList(1.0);
ObservedLimits = limitsScanObserved.GetContourList(1.0);
# ObservedLimitsOneSigmaUp = limitsScanObservedOneSigmaUp.GetContourList(1.0);
# ObservedLimitsOneSigmaDown = limitsScanObservedOneSigmaDown.GetContourList(1.0);
limitsScan.Draw("colz")

outputFile=ROOT.TFile("analysis/limitPlots/MassScanStealth.root", "RECREATE")
ExpectedLimits.Write("ExpectedLimits")
ExpectedLimitsOneSigmaUp.Write("ExpectedLimitsOneSigmaUp")
ExpectedLimitsOneSigmaDown.Write("ExpectedLimitsOneSigmaDown")
if (inputArguments.plotObservedLimits): ObservedLimits.Write("ObservedLimits")
# ObservedLimitsOneSigmaUp.Write("ObsLimSup")
# ObservedLimitsOneSigmaDown.Write("ObsLimSdn")
MassScan2D.Write("MassScan2D")
ExpectedLimitsMassScan2D.Write("ExpectedLimitsMassScan2D")
outputFile.Close()
