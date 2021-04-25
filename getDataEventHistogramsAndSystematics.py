#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, array, pdb, math
import ROOT, tmROOTUtils, tmStatsUtils, tmGeneralUtils
from tmProgressBar import tmProgressBar

ROOT.gROOT.SetBatch(ROOT.kTRUE)

inputArgumentsParser = argparse.ArgumentParser(description='Get data event histograms and systematics.')
inputArgumentsParser.add_argument('--inputFilesList', required=True, help='Semicolon-separated list of paths, possibly containing wildcards, to pass to TChain.Add.',type=str)
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
inputArgumentsParser.add_argument('--nDatasetDivisionsForNLL', default=5, help='If this parameter is N, then for the NLL curve, the input dataset in the norm jets bin is divided into N independent datasets with the same number of events. We then find the kernel estimate combining (N-1) datasets and take its NLL with respect to the remaining 1 dataset -- there are N ways of doing this. The net NLL is the sum of each individual NLLs.',type=int)
inputArgumentsParser.add_argument('--kernelMirrorOption', default="MirrorLeft", help='Kernel mirroring option to be used in adaptive Gaussian kernel estimates',type=str)
inputArgumentsParser.add_argument('--analyzeSignalBins', action='store_true', help="If this flag is set, then the signal region data is unblinded. Specifically, the kernels and data are plotted -- and the observed event counters are stored -- for all bins rather than only the normalization bin.")
inputArguments = inputArgumentsParser.parse_args()

kernelOptionsObjects = {"NoMirror": ROOT.RooKeysPdf.NoMirror, "MirrorLeft": ROOT.RooKeysPdf.MirrorLeft, "MirrorRight": ROOT.RooKeysPdf.MirrorRight, "MirrorBoth": ROOT.RooKeysPdf.MirrorBoth, "MirrorAsymLeft": ROOT.RooKeysPdf.MirrorAsymLeft, "MirrorAsymLeftRight": ROOT.RooKeysPdf.MirrorAsymLeftRight, "MirrorAsymRight": ROOT.RooKeysPdf.MirrorAsymRight, "MirrorLeftAsymRight": ROOT.RooKeysPdf.MirrorLeftAsymRight, "MirrorAsymBoth": ROOT.RooKeysPdf.MirrorAsymBoth}
if not(inputArguments.kernelMirrorOption in kernelOptionsObjects): sys.exit("The following element is passed as an argument for the kernel mirroring option but not in the dictionary defining the correspondence between kernel name and RooKeysPdf index: {kernelMirrorOption}".format(kernelMirrorOption=inputArguments.kernelMirrorOption))
if not(inputArguments.nJetsMax == 6): sys.exit("Only nJetsMax=6 supported temporarily. Needed to create data card template in the correct format.")

STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries)
STBoundaries = []
STBoundariesToPlot = []
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()):
        STBoundary = float(STBoundaryString.strip())
        STBoundaries.append(STBoundary)
        STBoundariesToPlot.append(STBoundary)
STBoundaries.append(20000.0) # Instead of infinity
STBoundariesToPlot.append(3500)
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

def getChiSqPerDOF(sourcePDF, targetHistogram, nEvents_normRange, nFitParameters, printDebug=False):
    print("Calculating chisq:")
    nDegreesOfFreedom = -(nFitParameters)
    totalChiSq = 0.
    for STRegionIndex in range(1, nSTSignalBins+2):
        resetSTRange()
        prediction = getExpectedNEventsFromPDFInNamedRange(nEvents_normRange, sourcePDF, "STRange_RegionIndex{i}".format(i=STRegionIndex))
        if (printDebug): print("prediction: {p}".format(p=prediction))
        observation = targetHistogram.GetBinContent(STRegionIndex)
        if (printDebug): print("observation: {o}".format(o=observation))
        observationError = targetHistogram.GetBinError(STRegionIndex)
        if (observation > 0):
            nDegreesOfFreedom += 1
            totalChiSq += pow((prediction-observation)/observationError, 2)
    if (printDebug): print("totalChiSq: {tCS}, nDOF: {n}".format(tCS = totalChiSq, n=nDegreesOfFreedom))
    if (nDegreesOfFreedom <= 0): sys.exit("ERROR: check input to chisq.")
    return (totalChiSq/nDegreesOfFreedom)

def getNormalizedIntegralOfPDFInNamedRange(inputRooPDF, targetRangeName):
    resetSTRange()
    integralObject = inputRooPDF.createIntegral(ROOT.RooArgSet(rooVar_sT), targetRangeName)
    integralObject_fullRange = inputRooPDF.createIntegral(ROOT.RooArgSet(rooVar_sT), "kernelEstimator_sTRange")
    normalizedIntegral = integralObject.getVal() / integralObject_fullRange.getVal()
    resetSTRange()
    return normalizedIntegral

def getKernelSystematics(sourceKernel=None, targetKernel=None):
    if ((sourceKernel is None) or (targetKernel is None)): sys.exit("ERROR: One of sourceKernel and targetKernel is None.")
    resetSTRange()
    systematicsDictionary = {}
    for STRegionIndex in range(1, nSTSignalBins+2):
        sourceRatio = getNormalizedIntegralOfPDFInNamedRange(sourceKernel, "STRange_RegionIndex{i}".format(i=STRegionIndex))/getNormalizedIntegralOfPDFInNamedRange(sourceKernel, "normalization_sTRange".format(i=STRegionIndex))
        targetRatio = getNormalizedIntegralOfPDFInNamedRange(targetKernel, "STRange_RegionIndex{i}".format(i=STRegionIndex))/getNormalizedIntegralOfPDFInNamedRange(targetKernel, "normalization_sTRange".format(i=STRegionIndex))
        systematicsDictionary[STRegionIndex] = (sourceRatio/targetRatio)-1.0
    resetSTRange()
    return systematicsDictionary

def plotSystematicsInSTBin(systematicsDictionary=None, outputFilePath=None, outputTitlePrefix=None, sourceNEvents_numerator=None, sourceNEvents_denominator=None):
    if ((systematicsDictionary is None) or
        (outputFilePath is None) or
        (outputTitlePrefix is None) or
        (sourceNEvents_numerator is None)):
        sys.exit("ERROR: plotSystematicsInSTBin: one of systematicsDictionary, outputFilePath, outputHistName, and sourceNEvents_numerator is None.")

    outputHistName = outputTitlePrefix.replace(" ", "_").lower()
    outputTitle = outputTitlePrefix + ", #frac{f_{shape}(target)}{f_{shape}(source)} - 1.0;ST (GeV);(ratio - 1.0)"
    errorsGraph = ROOT.TGraphAsymmErrors(STRegionsAxis.GetNbins())
    errorsGraph.SetName("errors_" + outputHistName)
    errorsGraph.SetTitle(outputTitle)
    systematicsHistogram = ROOT.TH1F(outputHistName, outputTitle, len(STBoundariesToPlot)-1, array.array('d', STBoundariesToPlot))
    for STRegionIndex in range(1, nSTSignalBins+2):
        systematicsHistogram.SetBinContent(STRegionIndex, systematicsDictionary[STRegionIndex])
        systematicsHistogram.SetBinError(STRegionIndex, 0.)
        poissonInterval_numerator = tmROOTUtils.getPoissonConfidenceInterval(observedNEvents=sourceNEvents_numerator[STRegionIndex])
        netError = (poissonInterval_numerator["upper"] - poissonInterval_numerator["lower"])/(poissonInterval_numerator["upper"] + poissonInterval_numerator["lower"])
        if not(sourceNEvents_denominator is None):
            poissonInterval_denominator = tmROOTUtils.getPoissonConfidenceInterval(observedNEvents=sourceNEvents_denominator[STRegionIndex])
            netError = math.sqrt(pow((poissonInterval_numerator["upper"] - poissonInterval_numerator["lower"])/(poissonInterval_numerator["upper"] + poissonInterval_numerator["lower"]), 2) + pow((poissonInterval_denominator["upper"] - poissonInterval_denominator["lower"])/(poissonInterval_denominator["upper"] + poissonInterval_denominator["lower"]), 2))
        errorsGraph.SetPoint(STRegionIndex-1, systematicsHistogram.GetXaxis().GetBinCenter(STRegionIndex), 0.)
        errorsGraph.SetPointEXlow(STRegionIndex-1, 0.5*systematicsHistogram.GetXaxis().GetBinWidth(STRegionIndex))
        errorsGraph.SetPointEXhigh(STRegionIndex-1, 0.5*systematicsHistogram.GetXaxis().GetBinWidth(STRegionIndex))
        errorsGraph.SetPointEYlow(STRegionIndex-1, netError)
        errorsGraph.SetPointEYhigh(STRegionIndex-1, netError)
    outputCanvas = ROOT.TCanvas("systematics_" + outputHistName, "systematics_" + outputHistName, 1024, 512)
    outputCanvas.cd()
    ROOT.gPad.SetTopMargin(0.2)
    ROOT.gStyle.SetOptStat(0)
    systematicsHistogram.SetMarkerStyle(ROOT.kOpenSquare)
    systematicsHistogram.SetMarkerSize(1.)
    systematicsHistogram.Draw()
    errorsGraph.SetFillColor(ROOT.kOrange-2)
    errorsGraph.Draw("2")
    systematicsHistogram.Draw("AXIS SAME")
    systematicsHistogram.Draw("SAME")
    systematicsHistogram.Draw("P0 SAME")
    outputCanvas.Update()
    systematicsHistogram.GetYaxis().SetRangeUser(-1.45, 1.45)
    outputCanvas.Update()
    lineAt0 = ROOT.TLine(systematicsHistogram.GetXaxis().GetXmin(), 0., systematicsHistogram.GetXaxis().GetXmax(), 0.)
    lineAt0.Draw()
    outputCanvas.Update()
    outputCanvas.SaveAs(outputFilePath)

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

inputChain = ROOT.TChain('ggNtuplizer/EventTree')
for patternToAdd in inputArguments.inputFilesList.strip().split(";"):
    print("Adding pattern: {p}".format(p=patternToAdd))
    inputChain.Add(patternToAdd)
nEntries = inputChain.GetEntries()
print ("Total number of available events: {nEntries}".format(nEntries=nEntries))

outputFile = ROOT.TFile.Open('{oD}/{oP}_savedObjects.root'.format(oD=inputArguments.outputDirectory_eventHistograms, oP=inputArguments.outputPrefix), 'recreate')

# Initialize TTrees
sTTrees = {}
sTArrays = {}
sTHistograms_forChi2 = {}
for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    treeName = "sTTree_{nJets}Jets".format(nJets=nJets)
    sTTrees[nJets] = ROOT.TTree(treeName, treeName)
    sTArrays[nJets] = array.array('f', [0.])
    (sTTrees[nJets]).Branch('rooVar_sT', (sTArrays[nJets]), 'rooVar_sT/F')
    sTHistograms_forChi2[nJets] = ROOT.TH1F("h_sT_forChi2_{n}JetsBin".format(n=nJets), "sT_forChi2_{n}JetsBin".format(n=nJets), len(STBoundariesToPlot)-1, array.array('d', STBoundariesToPlot))

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
progressBarUpdatePeriod = max(1, nEntries//50)
progressBar.initializeTimer()
for entryIndex in range(nEntries):
    entryStatus = inputChain.LoadTree(entryIndex)
    if entryStatus < 0: sys.exit("Tree failed to load entry at index {entryIndex}; returned status {entryStatus}".format(entryIndex=entryIndex, entryStatus=entryStatus))
    chainStatus = inputChain.GetEntry(entryIndex)
    if chainStatus <= 0: sys.exit("Unable to load data from chain at index {entryIndex}; GetEntry returned status {chainStatus}".format(entryIndex=entryIndex, chainStatus=chainStatus))

    if (entryIndex%progressBarUpdatePeriod == 0): progressBar.updateBar(1.0*entryIndex/nEntries, entryIndex)

    prefiringWeight_fromNTuples = inputChain.b_evtPrefiringWeight
    prefiringWeightsHistogram.Fill(prefiringWeight_fromNTuples)

    nStealthJets = inputChain.b_nJetsDR
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
        sTHistograms_forChi2[nJetsBin].Fill(sT)
        if (nJetsBin == inputArguments.nJetsNorm):
            (sTArrays_forNLLs[divisionIndex])[0] = sT
            (sTTrees_forNLLs[divisionIndex]).Fill()
            divisionIndex += 1
            if (divisionIndex == inputArguments.nDatasetDivisionsForNLL): divisionIndex = 0
progressBar.terminate()

prefiringWeightsCanvas = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [prefiringWeightsHistogram], canvasName = "c_prefiringWeights", outputROOTFile = outputFile, outputDocumentName = "{oD}/{oP}_prefiringWeights".format(oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix), customOptStat="oume", enableLogY=True)
statsBox = prefiringWeightsCanvas.GetPrimitive("stats")
statsBox.SetX1NDC(0.1)
statsBox.SetX2NDC(0.4)
prefiringWeightsCanvas.Update()
prefiringWeightsCanvas.SaveAs("{oD}/{oP}_prefiringWeights.pdf".format(oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix))

# Write observed nEvents to files
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    for STRegionIndex in range(1, nSTSignalBins+2):
        writeThisBin = True
        if not(inputArguments.analyzeSignalBins):
            if ((nJetsBin != inputArguments.nJetsNorm) and (STRegionIndex > 1)):
                writeThisBin = False
        if (writeThisBin):
            observedEventCountersList.append(tuple(["int", "observedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin), nEventsInSTRegions[STRegionIndex][nJetsBin]]))
tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=observedEventCountersList, outputFilePath=("{oD}/{oP}_observedEventCounters.dat".format(oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix)))

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
    print("Fractional uncertainty from Poisson errors on number of events in normalization bin at {n} jets: {a:.3f}".format(n=nJetsBin, a=fractionalUncertainties_nEvents_normRange[nJetsBin]))
    for STRegionIndex in range(1, nSTSignalBins+2): # Same uncertainty for all ST bins
        dataSystematicsList.append(tuple(["float", "fractionalUncertaintyDown_normEvents_STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin), (fractionalUncertainties_nEvents_normRange_factors_down[nJetsBin]-1.0)]))
        dataSystematicsList.append(tuple(["float", "fractionalUncertaintyUp_normEvents_STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin), (fractionalUncertainties_nEvents_normRange_factors_up[nJetsBin]-1.0)]))

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

dataSets_forNLL = {}
for divisionIndex in range(0, inputArguments.nDatasetDivisionsForNLL):
    dataSets_forNLL[divisionIndex] = ROOT.RooDataSet("dataSet", sTRooDataSetName, sTTrees_forNLLs[divisionIndex], ROOT.RooArgSet(rooVar_sT))

# Start finding systematics wrt rho
# First find nominal value of rho
def NLLAsAFunctionOfRho(rho):
    resetSTRange()
    sumNLL = 0.
    for testDataDivisionIndex in range(0, inputArguments.nDatasetDivisionsForNLL):
        estimationDataSet = ROOT.RooDataSet("estimationDataSet_testDataDivisionIndex{i}".format(i=testDataDivisionIndex), "estimationDataSet_testDataDivisionIndex{i}".format(i=testDataDivisionIndex), ROOT.RooArgSet(rooVar_sT))
        for estimationDataDivisionIndex in range(0, inputArguments.nDatasetDivisionsForNLL):
            if (estimationDataDivisionIndex == testDataDivisionIndex): continue
            estimationDataSet.append(dataSets_forNLL[estimationDataDivisionIndex])
        estimatedKernel = ROOT.RooKeysPdf("normBinKernelEstimateFunction_forNLL_testDataDivisionIndex_{i}".format(i=testDataDivisionIndex), "normBinKernelEstimateFunction_forNLL_testDataDivisionIndex_{i}".format(i=testDataDivisionIndex), rooVar_sT, estimationDataSet, kernelOptionsObjects[inputArguments.kernelMirrorOption], rho)
        nll = (estimatedKernel).createNLL(dataSets_forNLL[testDataDivisionIndex])
        sumNLL += nll.getVal()
    return sumNLL

# rhoNominal = tmStatsUtils.getStrictlyConvexFunctionApproximateMinimum(inputFunction=NLLAsAFunctionOfRho, xRange=[0.5, 3.5], autoZeroTolerance=True)
rhoNominal = tmStatsUtils.getGlobalMinimum(inputFunction=NLLAsAFunctionOfRho, xRange=[0.5, 3.5], autoZeroTolerance=True, printDebug=True)
print("Found rhoNominal = {rN}".format(rN=rhoNominal))
# Write out rhoNominal to file
rhoNominalOutputFileObject = open("{oD}/{oP}_rhoNominal.dat".format(oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix), 'w')
rhoNominalOutputFileObject.write("{rN}\n".format(rN=rhoNominal))
rhoNominalOutputFileObject.close()

nll_at_rhoNominal = NLLAsAFunctionOfRho(rhoNominal)
resetSTRange()
kernel_at_rhoNominal = ROOT.RooKeysPdf("normBinKernelEstimate_rhoNominal", "normBinKernelEstimate_rhoNominal", rooVar_sT, sTRooDataSets[inputArguments.nJetsNorm], kernelOptionsObjects[inputArguments.kernelMirrorOption], rhoNominal)
# The one-sigma fluctuations in rho are defined by nll(rho) = nll(rhoNominal) + 0.5 on both sides of the minimum
def nllValueAboveMinimumMinusHalf(rho):
    return (NLLAsAFunctionOfRho(rho) - nll_at_rhoNominal - 0.5)
rhoOneSigmaDown = tmStatsUtils.getMonotonicFunctionApproximateZero(inputFunction=nllValueAboveMinimumMinusHalf, xRange=[0.5*rhoNominal, rhoNominal], autoZeroTolerance=True, printDebug=True)
print("Found rhoNominal - 1*sigma = {rN1SD}".format(rN1SD=rhoOneSigmaDown))
resetSTRange()
kernel_at_rhoOneSigmaDown = ROOT.RooKeysPdf("normBinKernelEstimate_rhoOneSigmaDown", "normBinKernelEstimate_rhoOneSigmaDown", rooVar_sT, sTRooDataSets[inputArguments.nJetsNorm], kernelOptionsObjects[inputArguments.kernelMirrorOption], rhoOneSigmaDown)
sourceNEvents_numerator_forKernelSystematics = {STRegionIndex: nEventsInSTRegions[STRegionIndex][inputArguments.nJetsNorm] for STRegionIndex in range(1, nSTSignalBins+2)}
kernelSystematics_rhoOneSigmaDown = getKernelSystematics(sourceKernel=kernel_at_rhoOneSigmaDown, targetKernel=kernel_at_rhoNominal)
plotSystematicsInSTBin(systematicsDictionary=kernelSystematics_rhoOneSigmaDown, outputFilePath="{oD}/{oP}_systematics_rhoOneSigmaDown.pdf".format(oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix), outputTitlePrefix="Data systematics, #rho_{nominal}-1#sigma:", sourceNEvents_numerator=sourceNEvents_numerator_forKernelSystematics, sourceNEvents_denominator=None)
rhoOneSigmaUp = tmStatsUtils.getMonotonicFunctionApproximateZero(inputFunction=nllValueAboveMinimumMinusHalf, xRange=[rhoNominal, 2.0*rhoNominal], autoZeroTolerance=True, printDebug=True)
print("Found rhoNominal + 1*sigma = {rN1SU}".format(rN1SU=rhoOneSigmaUp))
resetSTRange()
kernel_at_rhoOneSigmaUp = ROOT.RooKeysPdf("normBinKernelEstimate_rhoOneSigmaUp", "normBinKernelEstimate_rhoOneSigmaUp", rooVar_sT, sTRooDataSets[inputArguments.nJetsNorm], kernelOptionsObjects[inputArguments.kernelMirrorOption], rhoOneSigmaUp)
kernelSystematics_rhoOneSigmaUp = getKernelSystematics(sourceKernel=kernel_at_rhoOneSigmaUp, targetKernel=kernel_at_rhoNominal)
plotSystematicsInSTBin(systematicsDictionary=kernelSystematics_rhoOneSigmaUp, outputFilePath="{oD}/{oP}_systematics_rhoOneSigmaUp.pdf".format(oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix), outputTitlePrefix="Data systematics, #rho_{nominal}+1#sigma:", sourceNEvents_numerator=sourceNEvents_numerator_forKernelSystematics, sourceNEvents_denominator=None)

# Now plot the NLL
rhoNLLGraph = ROOT.TGraph()
rhoNLLGraph.SetName("rhoNLL")
rhoNLLGraph.SetTitle("NLL;#rho;")
resetSTRange()
nRhoValuesForSystematicsEstimation = 201
rhoMinForSystematicsEstimation = rhoOneSigmaDown - (rhoNominal - rhoOneSigmaDown)
rhoMaxForSystematicsEstimation = rhoOneSigmaUp + (rhoOneSigmaUp - rhoNominal)
progressBar = tmProgressBar(nRhoValuesForSystematicsEstimation)
progressBarUpdatePeriod = max(1, nRhoValuesForSystematicsEstimation//50)
progressBar.initializeTimer()
for rhoCounter in range(0, nRhoValuesForSystematicsEstimation):
    resetSTRange()
    if rhoCounter%progressBarUpdatePeriod == 0: progressBar.updateBar(1.0*rhoCounter/nRhoValuesForSystematicsEstimation, rhoCounter)
    rhoValue = rhoMinForSystematicsEstimation + (rhoCounter/(nRhoValuesForSystematicsEstimation-1))*(rhoMaxForSystematicsEstimation - rhoMinForSystematicsEstimation)
    rhoNLLGraph.SetPoint(rhoNLLGraph.GetN(), rhoValue, NLLAsAFunctionOfRho(rhoValue))
    resetSTRange()
progressBar.terminate()

canvases["rhoValues"]["NLLGraph"] = ROOT.TCanvas("c_rhoNLLGraph", "c_rhoNLLGraph", 1024, 768)
canvases["rhoValues"]["NLLGraph"].SetBorderSize(0)
canvases["rhoValues"]["NLLGraph"].SetFrameBorderMode(0)
canvases["rhoValues"]["NLLGraph"].cd()
rhoNLLGraph.Draw()

lineFactory = ROOT.TLine()
textFactory = ROOT.TLatex()
textFactory.SetTextSize(0.03)
textFactory.SetTextAlign(11)

# Horizontal lines
lineFactory.SetLineColor(ROOT.kGreen+3)
textFactory.SetTextColor(ROOT.kGreen+3)
lineFactory.DrawLine(rhoNLLGraph.GetHistogram().GetXaxis().GetXmin(), nll_at_rhoNominal + 0.5, rhoNLLGraph.GetHistogram().GetXaxis().GetXmax(), nll_at_rhoNominal + 0.5)
textFactory.DrawLatex(rhoNLLGraph.GetHistogram().GetXaxis().GetXmin() + 0.025*(rhoNLLGraph.GetHistogram().GetXaxis().GetXmax() - rhoNLLGraph.GetHistogram().GetXaxis().GetXmin()), nll_at_rhoNominal + 0.5, "min NLL + 0.5")
canvases["rhoValues"]["NLLGraph"].Update()
lineFactory.SetLineColor(ROOT.kBlack)
textFactory.SetTextColor(ROOT.kBlack)
lineFactory.DrawLine(rhoNLLGraph.GetHistogram().GetXaxis().GetXmin(), nll_at_rhoNominal, rhoNLLGraph.GetHistogram().GetXaxis().GetXmax(), nll_at_rhoNominal)
textFactory.DrawLatex(rhoNLLGraph.GetHistogram().GetXaxis().GetXmin() + 0.025*(rhoNLLGraph.GetHistogram().GetXaxis().GetXmax() - rhoNLLGraph.GetHistogram().GetXaxis().GetXmin()), nll_at_rhoNominal, "min NLL")
canvases["rhoValues"]["NLLGraph"].Update()

# Vertical lines
textFactory.SetTextAngle(90)
lineFactory.SetLineColor(ROOT.kBlack)
textFactory.SetTextColor(ROOT.kBlack)
lineFactory.DrawLine(rhoNominal, rhoNLLGraph.GetHistogram().GetYaxis().GetXmin(), rhoNominal, rhoNLLGraph.GetHistogram().GetYaxis().GetXmax())
textFactory.DrawLatex(rhoNominal, 0.5*(rhoNLLGraph.GetHistogram().GetYaxis().GetXmin() + rhoNLLGraph.GetHistogram().GetYaxis().GetXmax()), "#rho = #rho_{{nominal}} = {rhoNom:.2f}".format(rhoNom=rhoNominal))
canvases["rhoValues"]["NLLGraph"].Update()
lineFactory.SetLineColor(ROOT.kBlue)
textFactory.SetTextColor(ROOT.kBlue)
lineFactory.DrawLine(rhoOneSigmaDown, rhoNLLGraph.GetHistogram().GetYaxis().GetXmin(), rhoOneSigmaDown, rhoNLLGraph.GetHistogram().GetYaxis().GetXmax())
textFactory.DrawLatex(rhoOneSigmaDown, 0.5*(rhoNLLGraph.GetHistogram().GetYaxis().GetXmin() + rhoNLLGraph.GetHistogram().GetYaxis().GetXmax()), "#rho = #rho_{nominal} - 1#sigma")
canvases["rhoValues"]["NLLGraph"].Update()
lineFactory.SetLineColor(ROOT.kRed)
textFactory.SetTextColor(ROOT.kRed)
lineFactory.DrawLine(rhoOneSigmaUp, rhoNLLGraph.GetHistogram().GetYaxis().GetXmin(), rhoOneSigmaUp, rhoNLLGraph.GetHistogram().GetYaxis().GetXmax())
textFactory.DrawLatex(rhoOneSigmaUp, 0.5*(rhoNLLGraph.GetHistogram().GetYaxis().GetXmin() + rhoNLLGraph.GetHistogram().GetYaxis().GetXmax()), "#rho = #rho_{nominal} + 1#sigma")
canvases["rhoValues"]["NLLGraph"].Update()

canvases["rhoValues"]["NLLGraph"].SaveAs("{oD}/{oP}_rhoNLL.pdf".format(oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix))
outputFile.WriteTObject(rhoNLLGraph)
outputFile.WriteTObject(canvases["rhoValues"]["NLLGraph"])

for STRegionIndex in range(1, nSTSignalBins+2):
    fractionalUncertainty_rho = 0.5*(abs(kernelSystematics_rhoOneSigmaDown[STRegionIndex]) + abs(kernelSystematics_rhoOneSigmaUp[STRegionIndex]))
    for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax): # Same uncertainty for all nJets bins
        dataSystematicsList.append(tuple(["float", "fractionalUncertainty_rho_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin), fractionalUncertainty_rho]))

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

dataAndKernelsLegend = ROOT.TLegend(0.65, 0.7, 0.9, 0.9)
resetSTRange()
sTFrames["rhoValues"]["dataAndKernels"] = rooVar_sT.frame(sTKernelEstimatorRangeMin, sTKernelEstimatorRangeMax, n_sTBins)
setFrameAesthetics(sTFrames["rhoValues"]["dataAndKernels"], "#it{S}_{T} (GeV)", "Events / ({STBinWidth} GeV)".format(STBinWidth=int(0.5+inputArguments.ST_binWidth)), "Kernel estimates, nJets = {nJets}".format(nJets=inputArguments.nJetsNorm))
sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.RefreshNorm(), ROOT.RooFit.LineColor(ROOT.kWhite))
kernel_at_rhoNominal.fitTo(sTRooDataSets[inputArguments.nJetsNorm], normalizationRange, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Optimize(0))
kernel_at_rhoNominal.plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.Range("normalization_sTRange", ROOT.kFALSE), normalizationRange_normRange, ROOT.RooFit.DrawOption("F"), ROOT.RooFit.FillColor(ROOT.kYellow), ROOT.RooFit.VLines())
kernel_at_rhoOneSigmaDown.fitTo(sTRooDataSets[inputArguments.nJetsNorm], normalizationRange, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Optimize(0))
kernel_at_rhoOneSigmaDown.plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.LineColor(ROOT.kBlue), kernelEstimatorRange)
rhoOneSigmaDown_legendEntry = dataAndKernelsLegend.AddEntry(kernel_at_rhoOneSigmaDown, "PDF estimate: #rho = #rho_{nominal} - 1#sigma")
customizeLegendEntryForLine(rhoOneSigmaDown_legendEntry, ROOT.kBlue)
kernel_at_rhoOneSigmaUp.fitTo(sTRooDataSets[inputArguments.nJetsNorm], normalizationRange, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Optimize(0))
kernel_at_rhoOneSigmaUp.plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.LineColor(ROOT.kRed), kernelEstimatorRange)
rhoOneSigmaUp_legendEntry = dataAndKernelsLegend.AddEntry(kernel_at_rhoOneSigmaUp, "PDF estimate: #rho = #rho_{nominal} + 1#sigma")
customizeLegendEntryForLine(rhoOneSigmaUp_legendEntry, ROOT.kRed)
kernel_at_rhoNominal.fitTo(sTRooDataSets[inputArguments.nJetsNorm], normalizationRange, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Optimize(0))
kernel_at_rhoNominal.plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.LineColor(ROOT.kBlack), normalizationRange_normRange, kernelEstimatorRange)
rho1_legendEntry = dataAndKernelsLegend.AddEntry(kernel_at_rhoNominal, "PDF estimate: #rho = #rho_{{nominal}} = {rhoNom:.2f}".format(rhoNom=rhoNominal))
customizeLegendEntryForLine(rho1_legendEntry, ROOT.kBlack)
sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.RefreshNorm())
dataSet_legendEntry = dataAndKernelsLegend.AddEntry(sTRooDataSets[inputArguments.nJetsNorm], "Data")
canvases["rhoValues"]["dataAndKernels"] = {}
canvases["rhoValues"]["dataAndKernels"]["linear"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["rhoValues"]["dataAndKernels"], dataAndKernelsLegend], canvasName = "c_kernelPDFEstimate_dataAndKernels_linearScale", outputROOTFile = outputFile, outputDocumentName = "{oD}/{oP}_kernelPDF_rhoValues_linearScale".format(oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix))
canvases["rhoValues"]["dataAndKernels"]["log"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["rhoValues"]["dataAndKernels"], dataAndKernelsLegend], canvasName = "c_kernelPDFEstimate_dataAndKernels", outputROOTFile = outputFile, outputDocumentName = "{oD}/{oP}_kernelPDF_rhoValues".format(oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix), enableLogY = True)

# Find estimators for norm bin
resetSTRange()
rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm] = ROOT.RooKeysPdf("normBinKernelEstimateFunction", "normBinKernelEstimateFunction", rooVar_sT, sTRooDataSets[inputArguments.nJetsNorm], kernelOptionsObjects[inputArguments.kernelMirrorOption], rhoNominal)
# rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].fitTo(sTRooDataSets[inputArguments.nJetsNorm], normalizationRange, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Optimize(0))
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
# rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], ROOT.RooFit.Range("normalization_sTRange", ROOT.kFALSE), normalizationRange_normRange, ROOT.RooFit.DrawOption("F"), ROOT.RooFit.FillColor(ROOT.kYellow), ROOT.RooFit.VLines())
rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], ROOT.RooFit.LineColor(ROOT.kBlue), kernelEstimatorRange, ROOT.RooFit.Normalization(fractionalUncertainties_nEvents_normRange_factors_up[inputArguments.nJetsNorm], ROOT.RooAbsReal.Relative), ROOT.RooFit.LineStyle(ROOT.kDashed))
rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], ROOT.RooFit.LineColor(ROOT.kBlue), kernelEstimatorRange, ROOT.RooFit.Normalization(fractionalUncertainties_nEvents_normRange_factors_down[inputArguments.nJetsNorm], ROOT.RooAbsReal.Relative), ROOT.RooFit.LineStyle(ROOT.kDashed))
rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], kernelEstimatorRange, ROOT.RooFit.Normalization(1.0, ROOT.RooAbsReal.Relative))
setFrameAesthetics(sTFrames["data"][inputArguments.nJetsNorm], "#it{S}_{T} (GeV)", "Events / ({STBinWidth} GeV)".format(STBinWidth=int(0.5+inputArguments.ST_binWidth)), "Normalization bin: {nJets} Jets".format(nJets=inputArguments.nJetsNorm))
sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.RefreshNorm())
canvases["data"][inputArguments.nJetsNorm] = {}
canvases["data"][inputArguments.nJetsNorm]["linear"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["data"][inputArguments.nJetsNorm]], canvasName = "c_kernelPDF_normJetsBin_linearScale", outputROOTFile = outputFile, outputDocumentName = "{oD}/{oP}_kernelPDF_normJetsBin_linearScale".format(oD=inputArguments.outputDirectory_eventHistograms, oP=inputArguments.outputPrefix))
canvases["data"][inputArguments.nJetsNorm]["log"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["data"][inputArguments.nJetsNorm]], canvasName = "c_kernelPDF_normJetsBin", outputROOTFile = outputFile, outputDocumentName = "{oD}/{oP}_kernelPDF_normJetsBin".format(oD=inputArguments.outputDirectory_eventHistograms, oP=inputArguments.outputPrefix), enableLogY = True)
rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].fitTo(sTRooDataSets[inputArguments.nJetsNorm], kernelEstimatorRange, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Optimize(0)) # Reset estimator normalization
resetSTRange()

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
progressBarUpdatePeriod = max(1, inputArguments.nToyMCs//50)
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
    rooKernel_PDF_Estimators["toyMC"][goodMCSampleIndex] = ROOT.RooKeysPdf("toyMCKernelEstimateFunction_{index}".format(index=goodMCSampleIndex), "toyMCKernelEstimateFunction_{index}".format(index=goodMCSampleIndex), rooVar_sT, toyRooDataSets[goodMCSampleIndex], kernelOptionsObjects[inputArguments.kernelMirrorOption], rhoNominal)
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
canvases["toyMC"]["DataAndEstimators"]["linear"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["toyMC"]["DataAndEstimators"]], canvasName = "c_toyMCDataAndKernelEstimates_linearScale", outputROOTFile = outputFile, outputDocumentName = "{oD}/{oP}_toyMCDataAndKernelEstimates_linearScale".format(oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix))
canvases["toyMC"]["DataAndEstimators"]["log"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["toyMC"]["DataAndEstimators"]], canvasName = "c_toyMCDataAndKernelEstimates", outputROOTFile = outputFile, outputDocumentName = "{oD}/{oP}_toyMCDataAndKernelEstimates".format(oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix), enableLogY = True)
resetSTRange()

# Estimate of the uncertainty = RMS of the distribution of the predicted to observed number of events
shapeSystematics = {} # Duplicate data; for convenience in later estimating residual ST scaling uncertainty
for STRegionIndex in range(1, nSTSignalBins+2):
    fractionalUncertainty_inThisSTRegion = toyVsOriginalIntegralsRatioHistograms[STRegionIndex].GetRMS()
    shapeSystematics[STRegionIndex] = fractionalUncertainty_inThisSTRegion
    for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax): # Same uncertainty for all nJets bins
        dataSystematicsList.append(tuple(["float", "fractionalUncertainty_shape_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin), fractionalUncertainty_inThisSTRegion]))
    # Plot the shape systematics estimate
    canvases["toyMC"]["systematics_STRegion{i}".format(i=STRegionIndex)] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [toyVsOriginalIntegralsRatioHistograms[STRegionIndex]], canvasName = "c_shapeSystematics_STRegion{i}".format(i=STRegionIndex), outputROOTFile = outputFile, outputDocumentName = "{oD}/{oP}_shapeSystematics_STRegion{i}".format(i=STRegionIndex, oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix))
plotSystematicsInSTBin(systematicsDictionary=shapeSystematics, outputFilePath="{oD}/{oP}_systematics_shape.pdf".format(oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix), outputTitlePrefix="Shape systematics:".format(n=nJetsBin), sourceNEvents_numerator=sourceNEvents_numerator_forKernelSystematics, sourceNEvents_denominator=None)

# Plot the integral checks, to see that the normalization is similar
canvases["toyMC"]["systematicsCheck"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [totalIntegralCheckHistogram], canvasName = "c_shapeSystematicsCheck", outputROOTFile = outputFile, outputDocumentName = "{oD}/{oP}_shapeSystematicsCheck".format(oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix))

# If ST scaling uncertainties are requested, then use these estimators in other nJets bins and obtain estimate of systematic on assumption that sT scales; otherwise, only obtain a prediction for number of events in all the ST regions
expected_nEventsInSTRegions = {}
for STRegionIndex in range(1, nSTSignalBins+2):
    expected_nEventsInSTRegions[STRegionIndex] = {}

chi2Values = {}
chi2Values[inputArguments.nJetsNorm] = getChiSqPerDOF(sourcePDF=rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm], targetHistogram=sTHistograms_forChi2[inputArguments.nJetsNorm], nEvents_normRange=nEventsInNormWindows[inputArguments.nJetsNorm], nFitParameters=2)
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    for STRegionIndex in range(1, nSTSignalBins+2):
        expected_nEventsInSTRegions[STRegionIndex][nJetsBin] = getExpectedNEventsFromPDFInNamedRange(nEvents_normRange=nEventsInNormWindows[nJetsBin], inputRooPDF=rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm], targetRangeName=("STRange_RegionIndex{i}".format(i = STRegionIndex)))
        expectedEventCountersList.append(tuple(["float", "expectedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin), expected_nEventsInSTRegions[STRegionIndex][nJetsBin]]))
    if (nJetsBin == inputArguments.nJetsNorm): continue
    sTFrames["data"][nJetsBin] = rooVar_sT.frame(sTKernelEstimatorRangeMin, sTKernelEstimatorRangeMax, n_sTBins)
    rooKernel_PDF_Estimators["data"][nJetsBin] = ROOT.RooKeysPdf("kernelEstimate_{nJetsBin}Jets".format(nJetsBin=nJetsBin), "kernelEstimate_{nJetsBin}Jets".format(nJetsBin=nJetsBin), rooVar_sT, sTRooDataSets[nJetsBin], kernelOptionsObjects[inputArguments.kernelMirrorOption], rhoNominal)
    rooKernel_PDF_Estimators["data"][nJetsBin].fitTo(sTRooDataSets[nJetsBin], normalizationRange, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Optimize(0))
    outputFile.WriteTObject(rooKernel_PDF_Estimators["data"][nJetsBin])
    if (inputArguments.analyzeSignalBins):
        sourceNEvents_numerator_forKernelSystematics = {STRegionIndex: nEventsInSTRegions[STRegionIndex][inputArguments.nJetsNorm] for STRegionIndex in range(1, nSTSignalBins+2)}
        sourceNEvents_denominator_forKernelSystematics = {STRegionIndex: nEventsInSTRegions[STRegionIndex][nJetsBin] for STRegionIndex in range(1, nSTSignalBins+2)}
        for STRegionIndex in range(1, nSTSignalBins+2):
            kernelSmoothedNEvents = getExpectedNEventsFromPDFInNamedRange(nEvents_normRange=nEventsInNormWindows[nJetsBin], inputRooPDF=rooKernel_PDF_Estimators["data"][nJetsBin], targetRangeName=("STRange_RegionIndex{i}".format(i = STRegionIndex)))
            expectedEventCountersList.append(tuple(["float", "expectedNEvents_realKernel_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin), kernelSmoothedNEvents]))
            scalingDeviation = kernelSmoothedNEvents - expected_nEventsInSTRegions[STRegionIndex][nJetsBin]
            expectedEventCountersList.append(tuple(["float", "scalingDeviation_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin), scalingDeviation]))
        fractionalUncertaintyDict_scaling_raw = getKernelSystematics(sourceKernel=rooKernel_PDF_Estimators["data"][nJetsBin], targetKernel=rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm])
        plotSystematicsInSTBin(systematicsDictionary=fractionalUncertaintyDict_scaling_raw, outputFilePath="{oD}/{oP}_systematics_scaling_{n}Jets.pdf".format(oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix, n=nJetsBin), outputTitlePrefix="ST scaling systematics, {n} Jets:".format(n=nJetsBin), sourceNEvents_numerator=sourceNEvents_numerator_forKernelSystematics, sourceNEvents_denominator=sourceNEvents_denominator_forKernelSystematics)

    sTRooDataSets[nJetsBin].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.RefreshNorm(), ROOT.RooFit.LineColor(ROOT.kWhite))
    rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].fitTo(sTRooDataSets[nJetsBin], normalizationRange, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Optimize(0))
    rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.Range("normalization_sTRange", ROOT.kFALSE), normalizationRange_normRange, ROOT.RooFit.DrawOption("F"), ROOT.RooFit.FillColor(ROOT.kYellow), ROOT.RooFit.VLines())
    rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.LineColor(ROOT.kBlue), kernelEstimatorRange, normalizationRange_normRange, ROOT.RooFit.Normalization(fractionalUncertainties_nEvents_normRange_factors_up[nJetsBin], ROOT.RooAbsReal.Relative), ROOT.RooFit.LineStyle(ROOT.kDashed))
    rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.LineColor(ROOT.kBlue), kernelEstimatorRange, normalizationRange_normRange, ROOT.RooFit.Normalization(fractionalUncertainties_nEvents_normRange_factors_down[nJetsBin], ROOT.RooAbsReal.Relative), ROOT.RooFit.LineStyle(ROOT.kDashed))
    rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.LineColor(ROOT.kBlue), kernelEstimatorRange, normalizationRange_normRange, ROOT.RooFit.Normalization(1.0, ROOT.RooAbsReal.Relative))
    sTRooDataSets[nJetsBin].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.RefreshNorm())
    dataHist = ROOT.RooDataHist("dataSetToHist_{nJB}JetsBin".format(nJB=nJetsBin), "dataSetToHist_{nJB}JetsBin".format(nJB=nJetsBin), ROOT.RooArgSet(rooVar_sT), sTRooDataSets[nJetsBin])
    chi2Values[nJetsBin] = getChiSqPerDOF(sourcePDF=rooKernel_PDF_Estimators["data"][inputArguments.nJetsNorm], targetHistogram=sTHistograms_forChi2[nJetsBin], nEvents_normRange=nEventsInNormWindows[nJetsBin], nFitParameters=2)

    if (nJetsBin == inputArguments.nJetsMax):
        setFrameAesthetics(sTFrames["data"][nJetsBin], "#it{S}_{T} (GeV)", "Events / ({STBinWidth} GeV)".format(STBinWidth=int(0.5+inputArguments.ST_binWidth)), "#geq {nJetsBin} Jets".format(nJetsBin=nJetsBin))
    else:
        setFrameAesthetics(sTFrames["data"][nJetsBin], "#it{S}_{T} (GeV)", "Events / ({STBinWidth} GeV)".format(STBinWidth=int(0.5+inputArguments.ST_binWidth)), "{nJetsBin} Jets".format(nJetsBin=nJetsBin))

    if (inputArguments.analyzeSignalBins):
        canvases["data"][nJetsBin] = {}
        canvases["data"][nJetsBin]["linear"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["data"][nJetsBin]], canvasName = "c_kernelPDF_{nJetsBin}Jets_linearScale".format(nJetsBin=nJetsBin), outputROOTFile = outputFile, outputDocumentName = "{oD}/{oP}_kernelPDF_{nJetsBin}Jets_linearScale".format(oD=inputArguments.outputDirectory_eventHistograms, oP=inputArguments.outputPrefix, nJetsBin=nJetsBin))
        canvases["data"][nJetsBin]["log"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["data"][nJetsBin]], canvasName = "c_kernelPDF_{nJetsBin}Jets".format(nJetsBin=nJetsBin), outputROOTFile = outputFile, outputDocumentName = "{oD}/{oP}_kernelPDF_{nJetsBin}Jets".format(oD=inputArguments.outputDirectory_eventHistograms, oP=inputArguments.outputPrefix, nJetsBin=nJetsBin), enableLogY = True)
    resetSTRange()

# Print out chisq values
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    print("nJetsBin = {nJB}: chi2 = {chi2}".format(nJB=nJetsBin, chi2=chi2Values[nJetsBin]))

tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=expectedEventCountersList, outputFilePath=("{oD}/{oP}_eventCounters.dat".format(oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix)))

# Write a few other useful things to output file
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    sTTrees[nJetsBin].Write()

outputFile.Write()
outputFile.Close()

tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=dataSystematicsList, outputFilePath=("{oD}/{oP}_dataSystematics.dat".format(oD=inputArguments.outputDirectory_dataSystematics, oP=inputArguments.outputPrefix)))
print("All done!")
