#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, ROOT, argparse, array, pdb, math, array
import numpy as np
import tmROOTUtils, tmStatsUtils, tmGeneralUtils
from tmProgressBar import tmProgressBar

inputArgumentsParser = argparse.ArgumentParser(description='Get data event histograms and systematics.')
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_eventHistograms', default="analysis/dataEventHistograms/", help='Directory in which to store event histograms.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_dataSystematics', default="analysis/dataSystematics/", help='Directory in which to store data systematics.',type=str)
inputArgumentsParser.add_argument('--outputPrefix', required=True, help='Prefix to all output file names.',type=str)
inputArgumentsParser.add_argument('--ST_binWidth', default=100., help='Target bin width in sT for plotting.',type=float)
inputArgumentsParser.add_argument('--sTKernelEstimatorRangeMax', default=3500., help='Max value of sT to use in the kernel estimator.',type=float)
inputArgumentsParser.add_argument('--preNormalizationBuffer', default=200., help='Buffer to use for the kernel before the lower boundary of the normalization region.',type=float)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Min number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Max number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsNorm', default=2, help='Number of jets w.r.t. which to normalize the sT distributions for other jets.',type=int)
inputArgumentsParser.add_argument('--nToyMCs', default=1000, help='Number of toy MC samples to generate using the pdf estimators found.',type=int)
inputArgumentsParser.add_argument('--nominalRho', default=1.25, help='Value of parameter rho to be used in all the nominal adaptive Gaussian kernel estimates.',type=float)
inputArgumentsParser.add_argument('--rhoMinFactorForSystematicsEstimation', default=0.5, help='Lower value of rho to use in systematics estimation = this argument * nominalRho.',type=float)
inputArgumentsParser.add_argument('--rhoMaxFactorForSystematicsEstimation', default=2.0, help='Upper value of rho to use in systematics estimation = this argument * nominalRho.',type=float)
inputArgumentsParser.add_argument('--nRhoValuesForSystematicsEstimation', default=151, help='Number of values of rho to use in systematics estimation.',type=int)
inputArgumentsParser.add_argument('--nDatasetDivisionsForNLL', default=3, help='If this parameter is N, then for the NLL curve, the input dataset in the norm jets bin is divided into N independent datasets with the same number of events. We then find the kernel estimate combining (N-1) datasets and take its NLL with respect to the remaining 1 dataset -- there are N ways of doing this. The net NLL is the sum of each individual NLLs.',type=int)
inputArgumentsParser.add_argument('--kernelMirrorOption', default="MirrorLeft", help='Kernel mirroring option to be used in adaptive Gaussian kernel estimates',type=str)
inputArgumentsParser.add_argument('--allowHigherNJets', action='store_true', help="Allow script to look into nJets bins beyond nJets = 3. Do not use with data with the signal selections before unblinding.")
inputArgumentsParser.add_argument('--isSignal', action='store_true', help="If this flag is set, then the input file is treated as belonging to the signal region; in that case, do not compute systematics on the degree to which sT scales, but compute all other systematics. If this flag is not set, then the input file is treated as belonging to the control region; in that case, compute the systematics estimate on the degree to which sT scales, but do not compute the other systematics.")
inputArguments = inputArgumentsParser.parse_args()

kernelOptionsObjects = {"NoMirror": ROOT.RooKeysPdf.NoMirror, "MirrorLeft": ROOT.RooKeysPdf.MirrorLeft, "MirrorRight": ROOT.RooKeysPdf.MirrorRight, "MirrorBoth": ROOT.RooKeysPdf.MirrorBoth, "MirrorAsymLeft": ROOT.RooKeysPdf.MirrorAsymLeft, "MirrorAsymLeftRight": ROOT.RooKeysPdf.MirrorAsymLeftRight, "MirrorAsymRight": ROOT.RooKeysPdf.MirrorAsymRight, "MirrorLeftAsymRight": ROOT.RooKeysPdf.MirrorLeftAsymRight, "MirrorAsymBoth": ROOT.RooKeysPdf.MirrorAsymBoth}
if not(inputArguments.kernelMirrorOption in kernelOptionsObjects): sys.exit("The following element is passed as an argument for the kernel mirroring option but not in the dictionary defining the correspondence between kernel name and RooKeysPdf index: {kernelMirrorOption}".format(kernelMirrorOption=inputArguments.kernelMirrorOption))
if not(inputArguments.nJetsMax == 6): sys.exit("Only nJetsMax=6 supported temporarily. Needed to create data card template in the correct format.")

STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries)
STBoundaries = []
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()):
        STBoundary = float(STBoundaryString.strip())
        STBoundaries.append(STBoundary)
STBoundaries.append(14000.0) # Instead of infinity
nSTSignalBins = len(STBoundaries) - 2 # First two lines are lower and upper boundaries for the normalization bin
STNormRangeMin = STBoundaries[0]
STNormRangeMax = STBoundaries[1]
print("Using {n} signal bins for ST; norm range min: {mn}, norm range max: {mx}.".format(n = nSTSignalBins, mn=STNormRangeMin, mx=STNormRangeMax))
STRegionsAxis = ROOT.TAxis(len(STBoundaries)-1, array.array('d', STBoundaries))
sTKernelEstimatorRangeMin = STNormRangeMin-inputArguments.preNormalizationBuffer
sTKernelEstimatorRangeMax = inputArguments.sTKernelEstimatorRangeMax

rooVar_sT = ROOT.RooRealVar("rooVar_sT", "rooVar_sT", sTKernelEstimatorRangeMin, sTKernelEstimatorRangeMax, "GeV")
rooVar_sT.setRange("preNormalization_sTRange", sTKernelEstimatorRangeMin, STNormRangeMin)
rooVar_sT.setRange("normalization_sTRange", STNormRangeMin, STNormRangeMax)
rooVar_sT.setRange("observation_sTRange", STNormRangeMax, sTKernelEstimatorRangeMax)
for STRegionIndex in range(1, nSTSignalBins+2):
    rooVar_sT.setRange("STRange_RegionIndex{i}".format(i = STRegionIndex), STBoundaries[STRegionIndex-1], STBoundaries[STRegionIndex])
rooVar_sT.setRange("kernelEstimator_sTRange", sTKernelEstimatorRangeMin, sTKernelEstimatorRangeMax)
normalizationRange_normRange = ROOT.RooFit.NormRange("normalization_sTRange")
normalizationRange = ROOT.RooFit.Range(STNormRangeMin, STNormRangeMax)
kernelEstimatorRange = ROOT.RooFit.Range(sTKernelEstimatorRangeMin, sTKernelEstimatorRangeMax)
n_sTBins = int(0.5 + ((sTKernelEstimatorRangeMax - sTKernelEstimatorRangeMin)/inputArguments.ST_binWidth))

dataSystematicsList = []
expectedEventCountersList = []
observedEventCountersList = []

def resetSTRange():
    rooVar_sT.setRange(sTKernelEstimatorRangeMin, sTKernelEstimatorRangeMax)

def setFrameAesthetics(frame, xLabel, yLabel, title):
    frame.SetXTitle(xLabel)
    frame.SetYTitle(yLabel)
    frame.SetTitle(title)

def getExpectedNEventsFromPDFInNamedRange(nEvents_normRange, inputRooPDF, targetRangeName):
    resetSTRange()
    integralObject_targetRange = inputRooPDF.createIntegral(ROOT.RooArgSet(rooVar_sT), targetRangeName)
    integralObject_normRange = inputRooPDF.createIntegral(ROOT.RooArgSet(rooVar_sT), "normalization_sTRange")
    targetToNormRatio = integralObject_targetRange.getVal() / integralObject_normRange.getVal()
    expectedNEvents = targetToNormRatio*nEvents_normRange
    resetSTRange()
    return expectedNEvents

def getNormalizedIntegralOfPDFInNamedRange(inputRooPDF, targetRangeName):
    resetSTRange()
    integralObject = inputRooPDF.createIntegral(ROOT.RooArgSet(rooVar_sT), targetRangeName)
    integralObject_fullRange = inputRooPDF.createIntegral(ROOT.RooArgSet(rooVar_sT), "kernelEstimator_sTRange")
    normalizedIntegral = integralObject.getVal() / integralObject_fullRange.getVal()
    resetSTRange()
    return normalizedIntegral

def getKernelSystematics(sourceKernel, targetKernel):
    resetSTRange()
    systematicsDictionary = {}
    for STRegionIndex in range(1, nSTSignalBins+2):
        sourceRatio = getNormalizedIntegralOfPDFInNamedRange(sourceKernel, "STRange_RegionIndex{i}".format(i=STRegionIndex))/getNormalizedIntegralOfPDFInNamedRange(sourceKernel, "normalization_sTRange".format(i=STRegionIndex))
        targetRatio = getNormalizedIntegralOfPDFInNamedRange(targetKernel, "STRange_RegionIndex{i}".format(i=STRegionIndex))/getNormalizedIntegralOfPDFInNamedRange(targetKernel, "normalization_sTRange".format(i=STRegionIndex))
        systematicsDictionary[STRegionIndex] = (sourceRatio/targetRatio)-1.0
    resetSTRange()
    return systematicsDictionary

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

inputChain = ROOT.TChain('ggNtuplizer/EventTree')
inputChain.Add(inputArguments.inputFilePath)
nEntries = inputChain.GetEntries()
print ("Total number of available events: {nEntries}".format(nEntries=nEntries))

outputFile = ROOT.TFile('{outputDirectory}/{outputPrefix}_savedObjects.root'.format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix), 'recreate')

# Initialize TTrees
sTTrees = {}
sTArrays = {}
for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    treeName = "sTTree_{nJets}Jets".format(nJets=nJets)
    sTTrees[nJets] = ROOT.TTree(treeName, treeName)
    sTArrays[nJets] = array.array('f', [0.])
    (sTTrees[nJets]).Branch('rooVar_sT', (sTArrays[nJets]), 'rooVar_sT/F')

sTTrees_forNLLs = {}
sTArrays_forNLLs = {}
for divisionIndex in range(0, inputArguments.nDatasetDivisionsForNLL):
    sTTrees_forNLLs[divisionIndex] = ROOT.TTree("sTTree_forNLL_divisionIndex{dI}".format(dI=divisionIndex), "sTTree_forNLL_divisionIndex{dI}".format(dI=divisionIndex))
    sTArrays_forNLLs[divisionIndex] = array.array('f', [0.])
    (sTTrees_forNLLs[divisionIndex]).Branch('rooVar_sT', (sTArrays_forNLLs[divisionIndex]), 'rooVar_sT/F')

nEventsInSTRegions = {}
for STRegionIndex in range(1, nSTSignalBins+2):
    nEventsInSTRegions[STRegionIndex] = {}
    for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
        nEventsInSTRegions[STRegionIndex][nJetsBin] = 0

# prefiring weights histogram
prefiringWeightsHistogram = ROOT.TH1F("h_prefiringWeights", "Distribution of prefiring weights;prefiring weight;nEvents", 510, -0.01, 1.01)

# Fill TTrees
divisionIndex = 0
progressBar = tmProgressBar(nEntries)
progressBarUpdatePeriod = max(1, nEntries//1000)
progressBar.initializeTimer()
for entryIndex in range(nEntries):
    entryStatus = inputChain.LoadTree(entryIndex)
    if entryStatus < 0: sys.exit("Tree failed to load entry at index {entryIndex}; returned status {entryStatus}".format(entryIndex=entryIndex, entryStatus=entryStatus))
    chainStatus = inputChain.GetEntry(entryIndex)
    if chainStatus <= 0: sys.exit("Unable to load data from chain at index {entryIndex}; GetEntry returned status {chainStatus}".format(entryIndex=entryIndex, chainStatus=chainStatus))

    if (entryIndex%progressBarUpdatePeriod == 0): progressBar.updateBar(1.0*entryIndex/nEntries, entryIndex)

    prefiringWeight_fromNTuples = inputChain.b_evtScaleFactor
    prefiringWeightsHistogram.Fill(prefiringWeight_fromNTuples)

    nStealthJets = inputChain.b_nJets
    nJetsBin = nStealthJets
    if (nJetsBin > inputArguments.nJetsMax): nJetsBin = inputArguments.nJetsMax
    if (nJetsBin < inputArguments.nJetsMin): sys.exit("Unexpected nJetsBin = {nJetsBin} at entry index = {entryIndex}".format(nJetsBin=nJetsBin, entryIndex=entryIndex))

    sT = inputChain.b_evtST
    STRegionIndex = STRegionsAxis.FindFixBin(sT)
    if (STRegionIndex > 0):
        nEventsInSTRegions[STRegionIndex][nJetsBin] += 1

    if (sT >= (sTKernelEstimatorRangeMin)):
        (sTArrays[nJetsBin])[0] = sT
        (sTTrees[nJetsBin]).Fill()
        if (nJetsBin == inputArguments.nJetsNorm):
            (sTArrays_forNLLs[divisionIndex])[0] = sT
            (sTTrees_forNLLs[divisionIndex]).Fill()
            divisionIndex += 1
            if (divisionIndex == inputArguments.nDatasetDivisionsForNLL): divisionIndex = 0
progressBar.terminate()

# prefiring weight histogram for control region
if not(inputArguments.isSignal):
    prefiringWeightsCanvas = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [prefiringWeightsHistogram], canvasName = "c_prefiringWeights", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_prefiringWeights".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix), customOptStat="oume", enableLogY=True)
    statsBox = prefiringWeightsCanvas.GetPrimitive("stats")
    statsBox.SetX1NDC(0.1)
    statsBox.SetX2NDC(0.4)
    prefiringWeightsCanvas.Update()
    prefiringWeightsCanvas.SaveAs("{outputDirectory}/{outputPrefix}_prefiringWeights".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix))

# Write observed nEvents to files
# For first two nJets bins write observed nEvents in all ST bins. For higher nJets bins, if we are analyzing the signal sample, then only write observed nEvents in the normalization bin.
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    if (nJetsBin <= 3 or not(inputArguments.isSignal)):
        for STRegionIndex in range(1, nSTSignalBins+2):
            observedEventCountersList.append(tuple(["int", "observedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin), nEventsInSTRegions[STRegionIndex][nJetsBin]]))
    elif (inputArguments.isSignal):
        observedEventCountersList.append(tuple(["int", "observedNEvents_STRegion1_{n}Jets".format(n=nJetsBin), nEventsInSTRegions[1][nJetsBin]]))
tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=observedEventCountersList, outputFilePath=("{outputDirectory}/{outputPrefix}_observedEventCounters.dat".format(outputDirectory=inputArguments.outputDirectory_dataSystematics, outputPrefix=inputArguments.outputPrefix)))

# Make datasets from all sT trees
sTRooDataSets = {}
nEventsInPreNormWindows = {}
nEventsInNormWindows = {}
nEventsInObservationWindows = {}
total_nEventsInFullRange = {}
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    sTRooDataSetName = "rooDataSet_{nJets}Jets".format(nJets=nJetsBin)
    sTRooDataSets[nJetsBin] = ROOT.RooDataSet(sTRooDataSetName, sTRooDataSetName, sTTrees[nJetsBin], ROOT.RooArgSet(rooVar_sT))
    outputFile.WriteTObject(sTRooDataSets[nJetsBin])
    nEventsInPreNormWindows[nJetsBin] = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(sTRooDataSets[nJetsBin], "preNormalization_sTRange", forceNumEntries=True)
    nEventsInNormWindows[nJetsBin] = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(sTRooDataSets[nJetsBin], "normalization_sTRange", forceNumEntries=True)
    nEventsInObservationWindows[nJetsBin] = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(sTRooDataSets[nJetsBin], "observation_sTRange", forceNumEntries=True)
    total_nEventsInFullRange[nJetsBin] = nEventsInPreNormWindows[nJetsBin] + nEventsInNormWindows[nJetsBin] + nEventsInObservationWindows[nJetsBin]
    if (nJetsBin <= 3 or (nJetsBin > 3 and inputArguments.allowHigherNJets)):
        print("At nJets = {nJets}, nEventsInPreNormWindow = {preNorm}, nEventsInNormWindow = {norm}, nEventsInObservationWindow = {obs}".format(nJets = nJetsBin, preNorm = nEventsInPreNormWindows[nJetsBin], norm = nEventsInNormWindows[nJetsBin], obs = nEventsInObservationWindows[nJetsBin]))

poissonConfidenceIntervals = {}
fractionalUncertainties_nEvents_normRange = {}
fractionalUncertainties_nEvents_normRange_factors_up = {}
fractionalUncertainties_nEvents_normRange_factors_down = {}
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    poissonConfidenceIntervals[nJetsBin] = tmROOTUtils.getPoissonConfidenceInterval(observedNEvents=nEventsInNormWindows[nJetsBin])
    fractionalUncertainties_nEvents_normRange[nJetsBin] = ((poissonConfidenceIntervals[nJetsBin])["upper"] - (poissonConfidenceIntervals[nJetsBin])["lower"])/(2*nEventsInNormWindows[nJetsBin])
    fractionalUncertainties_nEvents_normRange_factors_up[nJetsBin] = ((poissonConfidenceIntervals[nJetsBin])["upper"])/(nEventsInNormWindows[nJetsBin])
    fractionalUncertainties_nEvents_normRange_factors_down[nJetsBin] = ((poissonConfidenceIntervals[nJetsBin])["lower"])/(nEventsInNormWindows[nJetsBin])
    if (inputArguments.isSignal):
        dataSystematicsList.append(tuple(["float", "fractionalUncertainty_normEvents_{n}Jets".format(n=nJetsBin), (fractionalUncertainties_nEvents_normRange[nJetsBin])]))
        print("Fractional uncertainty from Poisson errors on number of events in normalization bin at {n} jets: {a:.3f}".format(n=nJetsBin, a=fractionalUncertainties_nEvents_normRange[nJetsBin]))

rooKernel_PDF_Estimators = {
    "data": {},
    "toyMC": {},
    "rhoValues": {}
}
canvases = {
    "data": {},
    "toyMC": {},
    "rhoValues": {}
}
sTFrames = {
    "data": {},
    "toyMC": {},
    "rhoValues": {}
}

# Find estimators for norm bin
resetSTRange()
rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm] = ROOT.RooKeysPdf("normBinKernelEstimateFunction", "normBinKernelEstimateFunction", rooVar_sT, sTRooDataSets[inputArguments.nJetsNorm], kernelOptionsObjects[inputArguments.kernelMirrorOption], inputArguments.nominalRho)
rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].fitTo(sTRooDataSets[inputArguments.nJetsNorm], normalizationRange, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Optimize(0))
ratioOfIntegrals_NormNJets = getNormalizedIntegralOfPDFInNamedRange(inputRooPDF=rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm], targetRangeName="observation_sTRange")/getNormalizedIntegralOfPDFInNamedRange(inputRooPDF=rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm], targetRangeName="normalization_sTRange")
outputFile.WriteTObject(rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm])
expected_nEvents_normJetsBin_normalizationRegion = getExpectedNEventsFromPDFInNamedRange(nEvents_normRange=nEventsInNormWindows[inputArguments.nJetsNorm], inputRooPDF=rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm], targetRangeName="normalization_sTRange") # closure check, just for sanity
expected_nEvents_normJetsBin_observationRegion = getExpectedNEventsFromPDFInNamedRange(nEvents_normRange=nEventsInNormWindows[inputArguments.nJetsNorm], inputRooPDF=rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm], targetRangeName="observation_sTRange")
predictedToObservedRatio_nominalRho_observationRegion = expected_nEvents_normJetsBin_observationRegion/nEventsInObservationWindows[inputArguments.nJetsNorm]
print("Check 1 on norm factor: expected events in norm range: {nExpected}, observed: {nObserved}".format(nExpected=expected_nEvents_normJetsBin_normalizationRegion, nObserved=nEventsInNormWindows[inputArguments.nJetsNorm]))
print("Check 2 on norm factor: expected events in observation range: {nExpected}, observed: {nObserved}".format(nExpected=expected_nEvents_normJetsBin_observationRegion, nObserved=nEventsInObservationWindows[inputArguments.nJetsNorm]))
resetSTRange()
sTFrames["data"][inputArguments.nJetsNorm] = rooVar_sT.frame(sTKernelEstimatorRangeMin, sTKernelEstimatorRangeMax, n_sTBins)
sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.LineColor(ROOT.kWhite), ROOT.RooFit.RefreshNorm())
rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], ROOT.RooFit.Range("normalization_sTRange", ROOT.kFALSE), normalizationRange_normRange, ROOT.RooFit.DrawOption("F"), ROOT.RooFit.FillColor(ROOT.kYellow), ROOT.RooFit.VLines())
rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], ROOT.RooFit.LineColor(ROOT.kBlue), kernelEstimatorRange, ROOT.RooFit.Normalization(fractionalUncertainties_nEvents_normRange_factors_up[inputArguments.nJetsNorm], ROOT.RooAbsReal.Relative), ROOT.RooFit.LineStyle(ROOT.kDashed))
rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], ROOT.RooFit.LineColor(ROOT.kBlue), kernelEstimatorRange, ROOT.RooFit.Normalization(fractionalUncertainties_nEvents_normRange_factors_down[inputArguments.nJetsNorm], ROOT.RooAbsReal.Relative), ROOT.RooFit.LineStyle(ROOT.kDashed))
rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], kernelEstimatorRange, ROOT.RooFit.Normalization(1.0, ROOT.RooAbsReal.Relative))
setFrameAesthetics(sTFrames["data"][inputArguments.nJetsNorm], "#it{S}_{T} (GeV)", "Events / ({STBinWidth} GeV)".format(STBinWidth=int(0.5+inputArguments.ST_binWidth)), "Normalization bin: {nJets} Jets".format(nJets=inputArguments.nJetsNorm))
sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.RefreshNorm())
canvases["data"][inputArguments.nJetsNorm] = {}
canvases["data"][inputArguments.nJetsNorm]["linear"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["data"][inputArguments.nJetsNorm]], canvasName = "c_kernelPDF_normJetsBin_linearScale", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_kernelPDF_normJetsBin_linearScale".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix))
canvases["data"][inputArguments.nJetsNorm]["log"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["data"][inputArguments.nJetsNorm]], canvasName = "c_kernelPDF_normJetsBin", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_kernelPDF_normJetsBin".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix), enableLogY = True)
resetSTRange()

# If input sample is the control sample, then use these estimators in other nJets bins and obtain estimate of systematic on assumption that sT scales; otherwise, only obtain a prediction for number of events in all the ST regions
expected_nEventsInSTRegions = {}
for STRegionIndex in range(1, nSTSignalBins+2):
    expected_nEventsInSTRegions[STRegionIndex] = {}
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    for STRegionIndex in range(1, nSTSignalBins+2):
        expected_nEventsInSTRegions[STRegionIndex][nJetsBin] = getExpectedNEventsFromPDFInNamedRange(nEvents_normRange=nEventsInNormWindows[nJetsBin], inputRooPDF=rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm], targetRangeName=("STRange_RegionIndex{i}".format(i = STRegionIndex)))
        expectedEventCountersList.append(tuple(["float", "expectedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin), expected_nEventsInSTRegions[STRegionIndex][nJetsBin]]))
    if (nJetsBin == inputArguments.nJetsNorm): continue
    if (nJetsBin > 3 and not(inputArguments.allowHigherNJets)): continue
    sTFrames["data"][nJetsBin] = rooVar_sT.frame(sTKernelEstimatorRangeMin, sTKernelEstimatorRangeMax, n_sTBins)
    rooKernel_PDF_Estimators["data"][nJetsBin] = ROOT.RooKeysPdf("kernelEstimate_{nJetsBin}Jets".format(nJetsBin=nJetsBin), "kernelEstimate_{nJetsBin}Jets".format(nJetsBin=nJetsBin), rooVar_sT, sTRooDataSets[nJetsBin], kernelOptionsObjects[inputArguments.kernelMirrorOption], inputArguments.nominalRho)
    rooKernel_PDF_Estimators["data"][nJetsBin].fitTo(sTRooDataSets[nJetsBin], normalizationRange, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Optimize(0))
    outputFile.WriteTObject(rooKernel_PDF_Estimators["data"][nJetsBin])
    fractionalUncertaintyDict_sTScaling = getKernelSystematics(sourceKernel=rooKernel_PDF_Estimators["data"][nJetsBin], targetKernel=rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm])
    if not(inputArguments.isSignal):
        for STRegionIndex in range(1, nSTSignalBins+2):
            dataSystematicsList.append(tuple(["float", "fractionalUncertainty_sTScaling_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin), abs(fractionalUncertaintyDict_sTScaling[STRegionIndex])]))

    sTRooDataSets[nJetsBin].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.RefreshNorm(), ROOT.RooFit.LineColor(ROOT.kWhite))
    rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].fitTo(sTRooDataSets[nJetsBin], normalizationRange, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Optimize(0))
    rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.Range("normalization_sTRange", ROOT.kFALSE), normalizationRange_normRange, ROOT.RooFit.DrawOption("F"), ROOT.RooFit.FillColor(ROOT.kYellow), ROOT.RooFit.VLines())
    rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.LineColor(ROOT.kBlue), kernelEstimatorRange, normalizationRange_normRange, ROOT.RooFit.Normalization(fractionalUncertainties_nEvents_normRange_factors_up[nJetsBin], ROOT.RooAbsReal.Relative), ROOT.RooFit.LineStyle(ROOT.kDashed))
    rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.LineColor(ROOT.kBlue), kernelEstimatorRange, normalizationRange_normRange, ROOT.RooFit.Normalization(fractionalUncertainties_nEvents_normRange_factors_down[nJetsBin], ROOT.RooAbsReal.Relative), ROOT.RooFit.LineStyle(ROOT.kDashed))
    rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.LineColor(ROOT.kBlue), kernelEstimatorRange, normalizationRange_normRange, ROOT.RooFit.Normalization(1.0, ROOT.RooAbsReal.Relative))
    sTRooDataSets[nJetsBin].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.RefreshNorm())

    if (nJetsBin == inputArguments.nJetsMax): setFrameAesthetics(sTFrames["data"][nJetsBin], "#it{S}_{T} (GeV)", "Events / ({STBinWidth} GeV)".format(STBinWidth=int(0.5+inputArguments.ST_binWidth)), "#geq {nJetsBin} Jets".format(nJetsBin=nJetsBin))
    else: setFrameAesthetics(sTFrames["data"][nJetsBin], "#it{S}_{T} (GeV)", "Events / ({STBinWidth} GeV)".format(STBinWidth=int(0.5+inputArguments.ST_binWidth)), "{nJetsBin} Jets".format(nJetsBin=nJetsBin))
    canvases["data"][nJetsBin] = {}
    canvases["data"][nJetsBin]["linear"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["data"][nJetsBin]], canvasName = "c_kernelPDF_{nJetsBin}Jets_linearScale".format(nJetsBin=nJetsBin), outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_kernelPDF_{nJetsBin}Jets_linearScale".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix, nJetsBin=nJetsBin))
    canvases["data"][nJetsBin]["log"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["data"][nJetsBin]], canvasName = "c_kernelPDF_{nJetsBin}Jets".format(nJetsBin=nJetsBin), outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_kernelPDF_{nJetsBin}Jets".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix, nJetsBin=nJetsBin), enableLogY = True)
    resetSTRange()

rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].fitTo(sTRooDataSets[inputArguments.nJetsNorm], normalizationRange, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Optimize(0)) # Reset estimator normalization

if not(inputArguments.isSignal):
    tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=dataSystematicsList, outputFilePath=("{outputDirectory}/{outputPrefix}_dataSystematics_sTScaling.dat".format(outputDirectory=inputArguments.outputDirectory_dataSystematics, outputPrefix=inputArguments.outputPrefix)))
    # No need to go further if the region is the control region
    outputFile.Write()
    outputFile.Close()
    sys.exit(0)

# Generate and estimate toy MC datsets with the individual kernels
toyRooDataSets = {}
sTFrames["toyMC"]["DataAndEstimators"] = rooVar_sT.frame(sTKernelEstimatorRangeMin, sTKernelEstimatorRangeMax, n_sTBins)
sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["toyMC"]["DataAndEstimators"], ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.RefreshNorm(), ROOT.RooFit.LineColor(ROOT.kWhite))
rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["toyMC"]["DataAndEstimators"], ROOT.RooFit.Range("normalization_sTRange", ROOT.kFALSE), normalizationRange_normRange, ROOT.RooFit.DrawOption("F"), ROOT.RooFit.FillColor(ROOT.kYellow), ROOT.RooFit.VLines())
toyVsOriginalIntegralsRatioHistograms = {}
for STRegionIndex in range(1, nSTSignalBins+2):
    toyVsOriginalIntegralsRatioHistograms[STRegionIndex] = ROOT.TH1F("h_toyVsOriginalIntegralsRatioHistogram_STRegion{i}".format(i=STRegionIndex), "(Toy MC/data) - 1.0, {STMin} < ST < {STMax};ratio;Toy MC events".format(STMin = STBoundaries[STRegionIndex-1], STMax = STBoundaries[STRegionIndex]), 40, 0., 0.)
totalIntegralCheckHistogram = ROOT.TH1F("h_totalIntegralCheck", "Total integral;total integral;Toy MC events", 40, 0., 0.)
goodMCSampleIndex = 0
randomGenerator = ROOT.TRandom1()
randomGenerator.SetSeed(99) # SetSeed = 0 would set seed by using some information from a ROOT "UUID", but this ensures that results are reproducible
progressBar = tmProgressBar(inputArguments.nToyMCs)
progressBarUpdatePeriod = max(1, inputArguments.nToyMCs//1000)
progressBar.initializeTimer()
while goodMCSampleIndex < inputArguments.nToyMCs:
    resetSTRange()
    toyRooDataSets[goodMCSampleIndex] = ROOT.RooDataSet("toyMCDataSet_index{goodMCSampleIndex}".format(goodMCSampleIndex=goodMCSampleIndex), "toyMCDataSet_index{goodMCSampleIndex}".format(goodMCSampleIndex=goodMCSampleIndex), ROOT.RooArgSet(rooVar_sT))
    rooVar_sT.setRange(sTKernelEstimatorRangeMin, STNormRangeMin)
    nPreNormEventsToGenerate = nEventsInPreNormWindows[inputArguments.nJetsNorm]
    dataSet_preNormWindow = rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].generate(ROOT.RooArgSet(rooVar_sT), nPreNormEventsToGenerate)
    toyRooDataSets[goodMCSampleIndex].append(dataSet_preNormWindow)
    rooVar_sT.setRange(STNormRangeMin, STNormRangeMax)
    nNormEventsToGenerate = nEventsInNormWindows[inputArguments.nJetsNorm]
    dataSet_normWindow = rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].generate(ROOT.RooArgSet(rooVar_sT), nNormEventsToGenerate)
    toyRooDataSets[goodMCSampleIndex].append(dataSet_normWindow)
    rooVar_sT.setRange(STNormRangeMax, sTKernelEstimatorRangeMax)
    for STRegionIndex in range(2, nSTSignalBins+2): # Starts from index 2: norm bin not included
        nEventsToGenerate = randomGenerator.Poisson(nEventsInSTRegions[STRegionIndex][inputArguments.nJetsNorm])
        if (STRegionIndex == (1+nSTSignalBins)): rooVar_sT.setRange(STBoundaries[STRegionIndex-1], sTKernelEstimatorRangeMax) # For last bin estimator only goes up to range max
        else: rooVar_sT.setRange(STBoundaries[STRegionIndex-1], STBoundaries[STRegionIndex])
        dataSet_observationRegion = rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].generate(ROOT.RooArgSet(rooVar_sT), nEventsToGenerate)
        toyRooDataSets[goodMCSampleIndex].append(dataSet_observationRegion)
        resetSTRange()
    rooKernel_PDF_Estimators["toyMC"][goodMCSampleIndex] = ROOT.RooKeysPdf("toyMCKernelEstimateFunction_{index}".format(index=goodMCSampleIndex), "toyMCKernelEstimateFunction_{index}".format(index=goodMCSampleIndex), rooVar_sT, toyRooDataSets[goodMCSampleIndex], kernelOptionsObjects[inputArguments.kernelMirrorOption], inputArguments.nominalRho)
    rooKernel_PDF_Estimators["toyMC"][goodMCSampleIndex].fitTo(toyRooDataSets[goodMCSampleIndex], normalizationRange, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Optimize(0))
    rooKernel_PDF_Estimators["toyMC"][goodMCSampleIndex].plotOn(sTFrames["toyMC"]["DataAndEstimators"], kernelEstimatorRange, normalizationRange_normRange, ROOT.RooFit.Normalization(1.0, ROOT.RooAbsReal.Relative))
    fractionalUncertaintyDict_toyShape = getKernelSystematics(sourceKernel=rooKernel_PDF_Estimators["toyMC"][goodMCSampleIndex], targetKernel=rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm])
    for STRegionIndex in range(1, nSTSignalBins+2):
        toyVsOriginalIntegralsRatioHistograms[STRegionIndex].Fill(fractionalUncertaintyDict_toyShape[STRegionIndex])
    totalIntegralCheckValue = getNormalizedIntegralOfPDFInNamedRange(inputRooPDF=rooKernel_PDF_Estimators["toyMC"][goodMCSampleIndex], targetRangeName="kernelEstimator_sTRange")
    totalIntegralCheckHistogram.Fill(totalIntegralCheckValue)
    resetSTRange()
    if goodMCSampleIndex%progressBarUpdatePeriod == 0: progressBar.updateBar(1.0*goodMCSampleIndex/inputArguments.nToyMCs, goodMCSampleIndex)
    goodMCSampleIndex += 1
progressBar.terminate()
print("Plotting datasets...")
for goodMCSampleIndex in range(0, inputArguments.nToyMCs):
    toyRooDataSets[goodMCSampleIndex].plotOn(sTFrames["toyMC"]["DataAndEstimators"], ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson)) # Plotting separately so that the datapoints won't be masked by the fits
print("Done plotting datasets.")
sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["toyMC"]["DataAndEstimators"], ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["toyMC"]["DataAndEstimators"], normalizationRange_normRange, ROOT.RooFit.LineColor(ROOT.kRed), kernelEstimatorRange)
resetSTRange()

# Plot the toy MC data and estimators
setFrameAesthetics(sTFrames["toyMC"]["DataAndEstimators"], "#it{S}_{T} (GeV)", "Events / ({STBinWidth} GeV)".format(STBinWidth=int(0.5+inputArguments.ST_binWidth)), "Data and estimators: {n} Toy MCs".format(n = inputArguments.nToyMCs))
canvases["toyMC"]["DataAndEstimators"] = {}
canvases["toyMC"]["DataAndEstimators"]["linear"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["toyMC"]["DataAndEstimators"]], canvasName = "c_toyMCDataAndKernelEstimates_linearScale", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_toyMCDataAndKernelEstimates_linearScale".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix))
canvases["toyMC"]["DataAndEstimators"]["log"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["toyMC"]["DataAndEstimators"]], canvasName = "c_toyMCDataAndKernelEstimates", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_toyMCDataAndKernelEstimates".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix), enableLogY = True)
resetSTRange()

# Estimate of the uncertainty = RMS of the distribution of the predicted to observed number of events
for STRegionIndex in range(1, nSTSignalBins+2):
    dataSystematicsList.append(tuple(["float", "fractionalUncertainty_Shape_STRegion{i}".format(i=STRegionIndex), toyVsOriginalIntegralsRatioHistograms[STRegionIndex].GetRMS()]))
    # Plot the shape systematics estimate
    canvases["toyMC"]["systematics__STRegion{i}".format(i=STRegionIndex)] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [toyVsOriginalIntegralsRatioHistograms[STRegionIndex]], canvasName = "c_shapeSystematics_STRegion{i}".format(i=STRegionIndex), outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_shapeSystematics_STRegion{i}".format(i=STRegionIndex, outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix))

# Plot the integral checks, to see that the normalization is similar
canvases["toyMC"]["systematicsCheck"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [totalIntegralCheckHistogram], canvasName = "c_shapeSystematicsCheck", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_shapeSystematicsCheck".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix))

dataSets_forNLL = {}
for divisionIndex in range(0, inputArguments.nDatasetDivisionsForNLL):
    dataSets_forNLL[divisionIndex] = ROOT.RooDataSet("dataSet", sTRooDataSetName, sTTrees_forNLLs[divisionIndex], ROOT.RooArgSet(rooVar_sT))

rhoSystematicsGraph = ROOT.TGraph()
rhoSystematicsGraph.SetName("rhoSystematics")
rhoSystematicsGraph.SetTitle("(Integral, ST > {stnormhi})/(Integral, {stnormlo} < ST < {stnormhi} );rho;".format(stnormlo = STNormRangeMin, stnormhi = STNormRangeMax))
rhoNLLGraph = ROOT.TGraph()
rhoNLLGraph.SetName("rhoNLL")
rhoNLLGraph.SetTitle("NLL;rho;")
rhoMinForSystematicsEstimation = inputArguments.rhoMinFactorForSystematicsEstimation * inputArguments.nominalRho
rhoMaxForSystematicsEstimation = inputArguments.rhoMaxFactorForSystematicsEstimation * inputArguments.nominalRho
resetSTRange()
kernelSystematics_rhoLow = {}
kernelSystematics_rhoHigh = {}
progressBar = tmProgressBar(inputArguments.nRhoValuesForSystematicsEstimation)
progressBarUpdatePeriod = max(1, inputArguments.nRhoValuesForSystematicsEstimation//1000)
progressBar.initializeTimer()
for rhoCounter in range(0, inputArguments.nRhoValuesForSystematicsEstimation): # Step 1: fill toy dataset with PDF at each rho value
    resetSTRange()
    if rhoCounter%progressBarUpdatePeriod == 0: progressBar.updateBar(1.0*rhoCounter/inputArguments.nRhoValuesForSystematicsEstimation, rhoCounter)
    rhoValue = rhoMinForSystematicsEstimation + (rhoCounter/(inputArguments.nRhoValuesForSystematicsEstimation-1))*(rhoMaxForSystematicsEstimation - rhoMinForSystematicsEstimation)
    rooKernel_PDF_Estimators["rhoValues"][rhoCounter] = ROOT.RooKeysPdf("normBinKernelEstimateFunction_rhoCounter_{rhoC}".format(rhoC=rhoCounter), "normBinKernelEstimateFunction_rhoCounter_{rhoC}".format(rhoC=rhoCounter), rooVar_sT, sTRooDataSets[inputArguments.nJetsNorm], kernelOptionsObjects[inputArguments.kernelMirrorOption], rhoValue)
    ratioOfIntegrals_atThisRho = getNormalizedIntegralOfPDFInNamedRange(inputRooPDF=rooKernel_PDF_Estimators["rhoValues"][rhoCounter], targetRangeName="observation_sTRange")/getNormalizedIntegralOfPDFInNamedRange(inputRooPDF=rooKernel_PDF_Estimators["rhoValues"][rhoCounter], targetRangeName="normalization_sTRange")
    rhoSystematicsGraph.SetPoint(rhoSystematicsGraph.GetN(), rhoValue, ratioOfIntegrals_atThisRho)
    if (rhoCounter == 0): kernelSystematics_rhoLow = getKernelSystematics(sourceKernel=rooKernel_PDF_Estimators["rhoValues"][0], targetKernel=rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm])
    if (rhoCounter == (-1 + inputArguments.nRhoValuesForSystematicsEstimation)): kernelSystematics_rhoHigh = getKernelSystematics(sourceKernel=rooKernel_PDF_Estimators["rhoValues"][-1 + inputArguments.nRhoValuesForSystematicsEstimation], targetKernel=rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm])

    sumNLLs = 0.
    for testDataDivisionIndex in range(0, inputArguments.nDatasetDivisionsForNLL):
        estimationDataSet = ROOT.RooDataSet("estimationDataSet_testDataDivisionIndex{i}".format(i=testDataDivisionIndex), "estimationDataSet_testDataDivisionIndex{i}".format(i=testDataDivisionIndex), ROOT.RooArgSet(rooVar_sT))
        for estimationDataDivisionIndex in range(0, inputArguments.nDatasetDivisionsForNLL):
            if (estimationDataDivisionIndex == testDataDivisionIndex): continue
            estimationDataSet.append(dataSets_forNLL[estimationDataDivisionIndex])
        estimatedKernel = ROOT.RooKeysPdf("normBinKernelEstimateFunction_forNLL_rhoCounter_{rhoC}_testDataDivisionIndex_{i}".format(rhoC=rhoCounter, i=testDataDivisionIndex), "normBinKernelEstimateFunction_forNLL_rhoCounter_{rhoC}_testDataDivisionIndex_{i}".format(rhoC=rhoCounter, i=testDataDivisionIndex), rooVar_sT, estimationDataSet, kernelOptionsObjects[inputArguments.kernelMirrorOption], rhoValue)
        nll = (estimatedKernel).createNLL(dataSets_forNLL[testDataDivisionIndex])
        sumNLLs += nll.getVal()
    rhoNLLGraph.SetPoint(rhoNLLGraph.GetN(), rhoValue, sumNLLs)
    resetSTRange()
progressBar.terminate()
canvases["rhoValues"]["ratiosGraph"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [rhoSystematicsGraph], canvasName = "c_rhoSystematicsGraph", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_rhoSystematics".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix))
canvases["rhoValues"]["NLLGraph"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [rhoNLLGraph], canvasName = "c_rhoNLLGraph", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_rhoNLL".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix))

for STRegionIndex in range(1, nSTSignalBins+2):
    fractionalUncertainty_rho = 0.5*(abs(kernelSystematics_rhoLow[STRegionIndex]) + abs(kernelSystematics_rhoHigh[STRegionIndex]))
    dataSystematicsList.append(tuple(["float", "fractionalUncertainty_rho_STRegion{i}".format(i=STRegionIndex), fractionalUncertainty_rho]))

# Plot kernels for three values of rho
def customizeLegendEntryForLine(entry, color):
    entry.ResetAttMarker()
    entry.SetMarkerColor(color)
    entry.SetMarkerStyle(8)
    entry.SetMarkerSize(0.125)
    entry.ResetAttFill()
    entry.SetFillStyle(0)
    entry.SetFillColor(color)
    entry.SetTextColor(color)
    entry.SetLineStyle(1)
    entry.SetLineColor(color)
    entry.SetLineWidth(3)

dataAndKernelsLegend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
resetSTRange()
sTFrames["rhoValues"]["dataAndKernels"] = rooVar_sT.frame(sTKernelEstimatorRangeMin, sTKernelEstimatorRangeMax, n_sTBins)
setFrameAesthetics(sTFrames["rhoValues"]["dataAndKernels"], "#it{S}_{T} (GeV)", "Events / ({STBinWidth} GeV)".format(STBinWidth=int(0.5+inputArguments.ST_binWidth)), "Kernel estimates, nJets = {nJets}".format(nJets=inputArguments.nJetsNorm))
sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.RefreshNorm(), ROOT.RooFit.LineColor(ROOT.kWhite))
rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.Range("normalization_sTRange", ROOT.kFALSE), normalizationRange_normRange, ROOT.RooFit.DrawOption("F"), ROOT.RooFit.FillColor(ROOT.kYellow), ROOT.RooFit.VLines())
rooKernel_PDF_Estimators["rhoValues"][0].plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.LineColor(ROOT.kBlue), kernelEstimatorRange)
rhoMin_legendEntry = dataAndKernelsLegend.AddEntry(rooKernel_PDF_Estimators["rhoValues"][0], "PDF estimate: rho = {rhoMin:.2f}".format(rhoMin=rhoMinForSystematicsEstimation))
customizeLegendEntryForLine(rhoMin_legendEntry, ROOT.kBlue)
rooKernel_PDF_Estimators["rhoValues"][-1+inputArguments.nRhoValuesForSystematicsEstimation].plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.LineColor(ROOT.kRed), kernelEstimatorRange)
rhoMax_legendEntry = dataAndKernelsLegend.AddEntry(rooKernel_PDF_Estimators["rhoValues"][-1+inputArguments.nRhoValuesForSystematicsEstimation], "PDF estimate: rho = {rhoMax:.2f}".format(rhoMax=rhoMaxForSystematicsEstimation))
customizeLegendEntryForLine(rhoMax_legendEntry, ROOT.kRed)
rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.LineColor(ROOT.kBlack), normalizationRange_normRange, kernelEstimatorRange)
rho1_legendEntry = dataAndKernelsLegend.AddEntry(rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm], "PDF estimate: rho = {rhoNom:.2f}".format(rhoNom=inputArguments.nominalRho))
customizeLegendEntryForLine(rho1_legendEntry, ROOT.kBlack)
sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.RefreshNorm())
dataSet_legendEntry = dataAndKernelsLegend.AddEntry(sTRooDataSets[inputArguments.nJetsNorm], "Data")
canvases["rhoValues"]["dataAndKernels"] = {}
canvases["rhoValues"]["dataAndKernels"]["linear"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["rhoValues"]["dataAndKernels"], dataAndKernelsLegend], canvasName = "c_kernelPDFEstimate_dataAndKernels_linearScale", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_kernelPDF_rhoValues_linearScale".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix))
canvases["rhoValues"]["dataAndKernels"]["log"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["rhoValues"]["dataAndKernels"], dataAndKernelsLegend], canvasName = "c_kernelPDFEstimate_dataAndKernels", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_kernelPDF_rhoValues".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix), enableLogY = True)

for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    sTTrees[nJetsBin].Write()

outputFile.Write()
outputFile.Close()

tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=dataSystematicsList, outputFilePath=("{outputDirectory}/{outputPrefix}_dataSystematics.dat".format(outputDirectory=inputArguments.outputDirectory_dataSystematics, outputPrefix=inputArguments.outputPrefix)))
tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=expectedEventCountersList, outputFilePath=("{outputDirectory}/{outputPrefix}_eventCounters.dat".format(outputDirectory=inputArguments.outputDirectory_dataSystematics, outputPrefix=inputArguments.outputPrefix)))

print("All done!")
