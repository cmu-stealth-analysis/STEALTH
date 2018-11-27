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
# inputArgumentsParser.add_argument('--sTPlotRangeMin', default=1200., help='Min value of sT to display in the plots.',type=float)
inputArgumentsParser.add_argument('--sTPlotRangeMax', default=3500., help='Max value of sT to display in the plots.',type=float)
inputArgumentsParser.add_argument('--n_sTBins', default=24, help='Number of sT bins (relevant for plotting only).',type=int)
# inputArgumentsParser.add_argument('--sTKernelFitRangeMin', default=700., help='Min value of sT to use in the kernel fit. This should be slightly less than sTPlotRangeMin to try to get rid of boundary effects.',type=float)
inputArgumentsParser.add_argument('--sTKernelFitRangeMax', default=3500., help='Max value of sT to use in the kernel fit.',type=float)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Min number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Max number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsNorm', default=2, help='Number of jets w.r.t. which to normalize the sT distributions for other jets.',type=int)
inputArgumentsParser.add_argument('--nTargetEventsForSTThresholdOptimization', default=1., help='Number of jets w.r.t. which to normalize the sT distributions for other jets.',type=float)
inputArgumentsParser.add_argument('--nToyMCs', default=1000, help='Number of toy MC samples to generate using the pdf fits found.',type=int)
# inputArgumentsParser.add_argument('--varyNEventsInPreNormWindowInToyMCs', action='store_true', help="Vary the number of generated events in toy MC samples in the pre-normalization window in a Poisson distribution about the number of events in this window in the original sample; default is to keep it fixed.")
inputArgumentsParser.add_argument('--varyNEventsInNormWindowInToyMCs', action='store_true', help="Vary the number of generated events in the toy MC samples in the normalization window in a Poisson distribution about the number of events in this window in the original sample; default is to keep it fixed.")
inputArgumentsParser.add_argument('--varyNEventsInObservationWindowInToyMCs', action='store_true', help="Vary the number of generated events in the toy MC samples in the observation window in a Poisson distribution about the number of events in this window in the original sample; default is to keep it fixed.")
inputArgumentsParser.add_argument('--nominalRho', default=1.2, help='Value of parameter rho to be used in all the nominal adaptive Gaussian kernel estimates.',type=float)
inputArgumentsParser.add_argument('--rhoMinFactorForSystematicsEstimation', default=0.5, help='Lower value of rho to use in systematics estimation = this argument * nominalRho.',type=float)
inputArgumentsParser.add_argument('--rhoMaxFactorForSystematicsEstimation', default=2.0, help='Upper value of rho to use in systematics estimation = this argument * nominalRho.',type=float)
inputArgumentsParser.add_argument('--nRhoValuesForSystematicsEstimation', default=151, help='Number of values of rho to use in systematics estimation.',type=int)
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

rooVar_scaleFactor = ROOT.RooRealVar("rooVar_scaleFactor", "rooVar_scaleFactor", 0.0, 100.0, "")

rooVar_sT = ROOT.RooRealVar("rooVar_sT", "rooVar_sT", STNormRangeMin, inputArguments.sTKernelFitRangeMax, "GeV")
rooVar_sT.setRange("normalization_sTRange", STNormRangeMin, STNormRangeMax)
rooVar_sT.setRange("observation_sTRange", STNormRangeMax, inputArguments.sTKernelFitRangeMax)
for STRegionIndex in range(1, nSTSignalBins+2):
    rooVar_sT.setRange("STRange_RegionIndex{i}".format(i = STRegionIndex), STBoundaries[STRegionIndex-1], STBoundaries[STRegionIndex])
rooVar_sT.setRange("kernelFit_sTRange", STNormRangeMin, inputArguments.sTKernelFitRangeMax)
plotRange = ROOT.RooFit.Range(STNormRangeMin, inputArguments.sTPlotRangeMax)
kernelFitRange = ROOT.RooFit.Range(STNormRangeMin, inputArguments.sTKernelFitRangeMax)
normRange = ROOT.RooFit.Range(STNormRangeMin, STNormRangeMax)
binWidth = int(0.5 + (inputArguments.sTPlotRangeMax - STNormRangeMin)/inputArguments.n_sTBins)

dataSystematicsList = []
expectedEventCountersList = []
observedEventCountersList = []

def setFrameAesthetics(frame, xLabel, yLabel, title):
    frame.SetXTitle(xLabel)
    frame.SetYTitle(yLabel)
    frame.SetTitle(title)

def getExpectedNEventsFromPDFInNamedRange(normFactor, inputRooPDF, inputRooArgSet, targetRangeMin, targetRangeMax):
    targetRangeName = "targetRange_sTOptimization"
    rooVar_sT.setRange(targetRangeName, targetRangeMin, targetRangeMax)
    integralObject_targetRange = inputRooPDF.createIntegral(inputRooArgSet, targetRangeName)
    expectedNEvents = normFactor*integralObject_targetRange.getVal()
    rooVar_sT.setRange(STNormRangeMin, inputArguments.sTKernelFitRangeMax)
    return expectedNEvents

inputChain = ROOT.TChain('ggNtuplizer/EventTree')
inputChain.Add(inputArguments.inputFilePath)
nEntries = inputChain.GetEntries()
print ("Total number of available events: {nEntries}".format(nEntries=nEntries))

outputFile = ROOT.TFile('{outputDirectory}/{outputPrefix}_savedObjects.root'.format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix), 'recreate')

# Initialize TTrees
sTTrees = {}
sTArrays = {}
scaleFactorArrays = {}
for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    treeName = "sTTree_{nJets}Jets".format(nJets=nJets)
    sTTrees[nJets] = ROOT.TTree(treeName, treeName)
    sTArrays[nJets] = array.array('f', [0.])
    scaleFactorArrays[nJets] = array.array('f', [0.])
    (sTTrees[nJets]).Branch('rooVar_sT', (sTArrays[nJets]), 'rooVar_sT/F')
    (sTTrees[nJets]).Branch('rooVar_scaleFactor', (scaleFactorArrays[nJets]), 'rooVar_scaleFactor/F')

nEventsInSTRegions = {}
weighted_nEventsInSTRegions = {}
for STRegionIndex in range(1, nSTSignalBins+2):
    nEventsInSTRegions[STRegionIndex] = {}
    weighted_nEventsInSTRegions[STRegionIndex] = {}
    for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
        nEventsInSTRegions[STRegionIndex][nJetsBin] = 0
        weighted_nEventsInSTRegions[STRegionIndex][nJetsBin] = 0.0

# Scale factors histogram
scaleFactorsHistogram = ROOT.TH1F("h_scaleFactors", "Distribution of scale factors;scale factor;nEvents", 400, 1.0, 1.02)

# Fill TTrees
progressBar = tmProgressBar(nEntries)
progressBarUpdatePeriod = max(1, nEntries//1000)
progressBar.initializeTimer()
for entryIndex in range(nEntries):
    entryStatus = inputChain.LoadTree(entryIndex)
    if entryStatus < 0: sys.exit("Tree failed to load entry at index {entryIndex}; returned status {entryStatus}".format(entryIndex=entryIndex, entryStatus=entryStatus))
    chainStatus = inputChain.GetEntry(entryIndex)
    if chainStatus <= 0: sys.exit("Unable to load data from chain at index {entryIndex}; GetEntry returned status {chainStatus}".format(entryIndex=entryIndex, chainStatus=chainStatus))

    if (entryIndex%progressBarUpdatePeriod == 0): progressBar.updateBar(1.0*entryIndex/nEntries, entryIndex)

    scaleFactor = (1.0/inputChain.b_evtScaleFactor) # real nEvents * evtScaleFactor = observed nEvents => each event has to be weighted by 1/(evtScaleFactor)
    scaleFactorsHistogram.Fill(scaleFactor)

    nStealthJets = inputChain.b_nJets
    nJetsBin = nStealthJets
    if (nJetsBin > inputArguments.nJetsMax): nJetsBin = inputArguments.nJetsMax
    if (nJetsBin < inputArguments.nJetsMin): sys.exit("Unexpected nJetsBin = {nJetsBin} at entry index = {entryIndex}".format(nJetsBin=nJetsBin, entryIndex=entryIndex))

    sT = inputChain.b_evtST
    STRegionIndex = STRegionsAxis.FindFixBin(sT)
    if (STRegionIndex > 0):
        nEventsInSTRegions[STRegionIndex][nJetsBin] += 1
        weighted_nEventsInSTRegions[STRegionIndex][nJetsBin] += scaleFactor

    if (sT >= STNormRangeMin and sT < inputArguments.sTKernelFitRangeMax):
        (sTArrays[nJetsBin])[0] = sT
        (scaleFactorArrays[nJetsBin])[0] = scaleFactor
        (sTTrees[nJetsBin]).Fill()
progressBar.terminate()

# Scale factors histogram for control region
if not(inputArguments.isSignal): tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [scaleFactorsHistogram], canvasName = "c_scaleFactors", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_scaleFactors".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix), customOptStat="oume")

# Write observed nEvents to files
# For first two nJets bins print observed nEvents in all ST bins. For higher nJets bins, if we are analyzing the signal sample, then only print observed nEvents in the normalization bin.
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    if (nJetsBin <= 3 or not(inputArguments.isSignal)):
        for STRegionIndex in range(1, nSTSignalBins+2):
            observedEventCountersList.append(tuple(["float", "observedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin), weighted_nEventsInSTRegions[STRegionIndex][nJetsBin]]))
    elif (inputArguments.isSignal):
        observedEventCountersList.append(tuple(["float", "observedNEvents_STRegion1_{n}Jets".format(n=nJetsBin), weighted_nEventsInSTRegions[1][nJetsBin]]))
tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=observedEventCountersList, outputFilePath=("{outputDirectory}/{outputPrefix}_observedEventCounters.dat".format(outputDirectory=inputArguments.outputDirectory_dataSystematics, outputPrefix=inputArguments.outputPrefix)))

# Make datasets from all sT trees
sTRooDataSets = {}
weighted_sTRooDataSets = {}
nEventsInNormWindows = {}
nEventsInObservationWindows = {}
total_nEventsInFullRange = {}
weighted_nEventsInNormWindows = {}
weighted_nEventsInObservationWindows = {}
weighted_total_nEventsInFullRange = {}
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    sTRooDataSetName = "rooDataSet_{nJets}Jets".format(nJets=nJetsBin)
    sTRooDataSets[nJetsBin] = ROOT.RooDataSet(sTRooDataSetName, sTRooDataSetName, sTTrees[nJetsBin], ROOT.RooArgSet(rooVar_sT))
    weighted_sTRooDataSets[nJetsBin] = ROOT.RooDataSet("weighted_" + sTRooDataSetName, "weighted_" + sTRooDataSetName, ROOT.RooArgSet(rooVar_sT, rooVar_scaleFactor), ROOT.RooFit.WeightVar(rooVar_scaleFactor), ROOT.RooFit.Import(sTTrees[nJetsBin]))
    nEventsInNormWindows[nJetsBin] = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(sTRooDataSets[nJetsBin], "normalization_sTRange", forceNumEntries=True)
    nEventsInObservationWindows[nJetsBin] = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(sTRooDataSets[nJetsBin], "observation_sTRange", forceNumEntries=True)
    total_nEventsInFullRange[nJetsBin] = nEventsInNormWindows[nJetsBin] + nEventsInObservationWindows[nJetsBin]
    weighted_nEventsInNormWindows[nJetsBin] = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(weighted_sTRooDataSets[nJetsBin], "normalization_sTRange")
    weighted_nEventsInObservationWindows[nJetsBin] = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(weighted_sTRooDataSets[nJetsBin], "observation_sTRange")
    weighted_total_nEventsInFullRange[nJetsBin] = weighted_nEventsInNormWindows[nJetsBin] + weighted_nEventsInObservationWindows[nJetsBin]
    if (nJetsBin <= 3 or (nJetsBin > 3 and inputArguments.allowHigherNJets)):
        print("At nJets = {nJets}, nEventsInNormWindow = {norm}, nEventsInObservationWindow = {obs}, weighted nEventsInNormWindow = {wnorm}, weighted_nEventsInObservationWindow = {wobs}".format(nJets = nJetsBin, norm = nEventsInNormWindows[nJetsBin], obs = nEventsInObservationWindows[nJetsBin], wnorm = weighted_nEventsInNormWindows[nJetsBin], wobs = weighted_nEventsInObservationWindows[nJetsBin]))

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
    dataSystematicsList.append(tuple(["float", "fractionalUncertainty_normEvents", fractionalUncertainties_nEvents_normRange[inputArguments.nJetsNorm]]))
    print("Fractional uncertainty from Poisson errors on number of events in normalization bin: {a:.3f}".format(a=fractionalUncertainties_nEvents_normRange[inputArguments.nJetsNorm]))

rooKernel_PDF_Fits = {
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

# Find fits for norm bin
rooVar_sT.setRange(STNormRangeMin, inputArguments.sTKernelFitRangeMax)
rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm] = ROOT.RooKeysPdf("normBinKernelEstimateFunction", "normBinKernelEstimateFunction", rooVar_sT, weighted_sTRooDataSets[inputArguments.nJetsNorm], kernelOptionsObjects[inputArguments.kernelMirrorOption], inputArguments.nominalRho)
integralObject_normalizationRange = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].createIntegral(ROOT.RooArgSet(rooVar_sT), "normalization_sTRange")
backgroundFit_normRange_integralValue = integralObject_normalizationRange.getVal()
normFactor_checks = weighted_nEventsInNormWindows[inputArguments.nJetsNorm]/backgroundFit_normRange_integralValue
expected_nEvents_normJetsBin_observationRegion = getExpectedNEventsFromPDFInNamedRange(normFactor=normFactor_checks, inputRooPDF=rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], inputRooArgSet=ROOT.RooArgSet(rooVar_sT), targetRangeMin=STNormRangeMax, targetRangeMax=inputArguments.sTKernelFitRangeMax)
observed_nEvents_normJetsBin_observationRegion = weighted_nEventsInObservationWindows[inputArguments.nJetsNorm]
predictedToObservedRatio_nominalRho_observationRegion = expected_nEvents_normJetsBin_observationRegion/observed_nEvents_normJetsBin_observationRegion
print("Check 1 on norm factor: expected events in norm range: {nExpected}, observed: {nObserved}".format(nExpected=getExpectedNEventsFromPDFInNamedRange(normFactor=normFactor_checks, inputRooPDF=rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], inputRooArgSet=ROOT.RooArgSet(rooVar_sT), targetRangeMin=STNormRangeMin, targetRangeMax=STNormRangeMax), nObserved=nEventsInNormWindows[inputArguments.nJetsNorm]))
print("Check 2 on norm factor: expected events in observation range: {nExpected}, observed: {nObserved}".format(nExpected=expected_nEvents_normJetsBin_observationRegion, nObserved=observed_nEvents_normJetsBin_observationRegion))
rooVar_sT.setRange(STNormRangeMin, inputArguments.sTKernelFitRangeMax)
sTFrames["data"][inputArguments.nJetsNorm] = rooVar_sT.frame(STNormRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
weighted_sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], ROOT.RooFit.RefreshNorm())
rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], ROOT.RooFit.LineColor(ROOT.kBlue), plotRange, ROOT.RooFit.Normalization(fractionalUncertainties_nEvents_normRange_factors_up[inputArguments.nJetsNorm], ROOT.RooAbsReal.Relative), ROOT.RooFit.LineStyle(ROOT.kDashed))
rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], ROOT.RooFit.LineColor(ROOT.kBlue), plotRange, ROOT.RooFit.Normalization(fractionalUncertainties_nEvents_normRange_factors_down[inputArguments.nJetsNorm], ROOT.RooAbsReal.Relative), ROOT.RooFit.LineStyle(ROOT.kDashed))
rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], plotRange, ROOT.RooFit.Normalization(1.0, ROOT.RooAbsReal.Relative))
setFrameAesthetics(sTFrames["data"][inputArguments.nJetsNorm], "#it{S}_{T} (GeV)", "Events / ({binWidth} GeV)".format(binWidth=binWidth), "Normalization bin: {nJets} Jets".format(nJets=inputArguments.nJetsNorm))
canvases["data"][inputArguments.nJetsNorm] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["data"][inputArguments.nJetsNorm]], canvasName = "c_kernelPDF_normJetsBin", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_kernelPDF_normJetsBin".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix), enableLogY = True)

# If input sample is the control sample, then use these fits in other nJets bins and obtain estimate of systematic on assumption that sT scales; otherwise, only obtain a prediction for number of events in all the ST regions
rooVar_nEventsInNormBin = {}
# rooKernel_extendedPDF_Fits = {}
integralObject_normJets_normRange = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].createIntegral(ROOT.RooArgSet(rooVar_sT), "normalization_sTRange")
integralObject_normJets_observationRange = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].createIntegral(ROOT.RooArgSet(rooVar_sT), "observation_sTRange")
integralsRatio_normJets = 1.0*integralObject_normJets_observationRange.getVal() / integralObject_normJets_normRange.getVal()

integralObjects_normJets_STRegions = {}
integralRatios_normJets_STRegions = {}
for STRegionIndex in range(1, nSTSignalBins+2):
    integralObjects_normJets_STRegions[STRegionIndex] = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].createIntegral(ROOT.RooArgSet(rooVar_sT), "STRange_RegionIndex{i}".format(i = STRegionIndex))
    integralRatios_normJets_STRegions[STRegionIndex] = integralObjects_normJets_STRegions[STRegionIndex].getVal() / integralObject_normJets_normRange.getVal()

expected_nEventsInSTRegions = {}
for STRegionIndex in range(1, nSTSignalBins+2):
    expected_nEventsInSTRegions[STRegionIndex] = {}

# if inputArguments.isSignal: # Calculate optimal sT threshold
#     normFactor_optimization = nEventsInNormWindows[nJetsBin]/backgroundFit_normRange_integralValue
#     print("Check 1 on norm factor for optimization: expected events in norm range: {nExpected}, observed: {nObserved}".format(nExpected=getExpectedNEventsFromPDFInNamedRange(normFactor=normFactor_optimization, inputRooPDF=rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], inputRooArgSet=ROOT.RooArgSet(rooVar_sT), targetRangeMin=STNormRangeMin, targetRangeMax=STNormRangeMax), nObserved=nEventsInNormWindows[nJetsBin]))
#     # print("Check 2 on norm factor for optimization: expected events in observation range: {nExpected}, observed: {nObserved}".format(nExpected=getExpectedNEventsFromPDFInNamedRange(normFactor=normFactor_optimization, inputRooPDF=rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], inputRooArgSet=ROOT.RooArgSet(rooVar_sT), targetRangeMin=STNormRangeMax, targetRangeMax=inputArguments.sTKernelFitRangeMax), nObserved=nEventsInObservationWindows[nJetsBin])) # Check 2 is disabled because we don't want to unblind ourselves yet
#     optimal_sTThreshold = tmStatsUtils.getMonotonicFunctionApproximateZero(inputFunction=(lambda sTThreshold: (-1*inputArguments.nTargetEventsForSTThresholdOptimization)+getExpectedNEventsFromPDFInNamedRange(normFactor=normFactor_optimization, inputRooPDF=rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], inputRooArgSet=ROOT.RooArgSet(rooVar_sT), targetRangeMin=sTThreshold, targetRangeMax=inputArguments.sTKernelFitRangeMax)), xRange=[1.1*STNormRangeMax, 0.99*inputArguments.sTKernelFitRangeMax], autoZeroTolerance=True, printDebug=True)
#     optimalThresholdOutputFile = open("{outputDirectory}/{outputPrefix}_sTOptimalThreshold.dat".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix), 'w')
#     optimalThresholdOutputFile.write("Optimal sT Threshold: {thr}\n".format(thr=optimal_sTThreshold))
#     optimalThresholdOutputFile.close()

for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    for STRegionIndex in range(1, nSTSignalBins+2):
        expected_nEventsInSTRegions[STRegionIndex][nJetsBin] = integralRatios_normJets_STRegions[STRegionIndex]*weighted_nEventsInNormWindows[nJetsBin]
        expectedEventCountersList.append(tuple(["float", "expectedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin), expected_nEventsInSTRegions[STRegionIndex][nJetsBin]]))
    if (nJetsBin == inputArguments.nJetsNorm): continue
    if (nJetsBin > 3 and not(inputArguments.allowHigherNJets)): continue
    sTFrames["data"][nJetsBin] = rooVar_sT.frame(STNormRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
    rooKernel_PDF_Fits["data"][nJetsBin] = ROOT.RooKeysPdf("kernelEstimate_{nJetsBin}Jets".format(nJetsBin=nJetsBin), "kernelEstimate_{nJetsBin}Jets".format(nJetsBin=nJetsBin), rooVar_sT, weighted_sTRooDataSets[nJetsBin], kernelOptionsObjects[inputArguments.kernelMirrorOption], inputArguments.nominalRho)
    integralObject_observationRange = rooKernel_PDF_Fits["data"][nJetsBin].createIntegral(ROOT.RooArgSet(rooVar_sT), "observation_sTRange")
    integralObject_normRange = rooKernel_PDF_Fits["data"][nJetsBin].createIntegral(ROOT.RooArgSet(rooVar_sT), "normalization_sTRange")
    integralsRatio = 1.0 * integralObject_observationRange.getVal()/integralObject_normRange.getVal()
    fraction_predictedToActual = integralsRatio/integralsRatio_normJets
    fractionalUncertainty_sTScaling = abs(fraction_predictedToActual - 1.0)

    weighted_sTRooDataSets[nJetsBin].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.RefreshNorm())
    rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.LineColor(ROOT.kBlue), plotRange, ROOT.RooFit.NormRange("normalization_sTRange"), ROOT.RooFit.Normalization(fractionalUncertainties_nEvents_normRange_factors_up[nJetsBin], ROOT.RooAbsReal.Relative), ROOT.RooFit.LineStyle(ROOT.kDashed))
    rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.LineColor(ROOT.kBlue), plotRange, ROOT.RooFit.NormRange("normalization_sTRange"), ROOT.RooFit.Normalization(fractionalUncertainties_nEvents_normRange_factors_down[nJetsBin], ROOT.RooAbsReal.Relative), ROOT.RooFit.LineStyle(ROOT.kDashed))
    rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.LineColor(ROOT.kBlue), plotRange, ROOT.RooFit.NormRange("normalization_sTRange"), ROOT.RooFit.Normalization(1.0, ROOT.RooAbsReal.Relative))
    # rooVar_nEventsInNormBin[nJetsBin] = ROOT.RooRealVar("rooVar_nEventsInNormBin_{nJetsBin}Jets".format(nJetsBin=nJetsBin), "rooVar_nEventsInNormBin_{nJetsBin}Jets".format(nJetsBin=nJetsBin), 100, 0, 10000)
    # rooKernel_extendedPDF_Fits[nJetsBin] = ROOT.RooExtendPdf("extendedKernelPDF_{nJetsBin}Jets".format(nJetsBin=nJetsBin), "extendedKernelPDF_{nJetsBin}Jets".format(nJetsBin=nJetsBin), rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], rooVar_nEventsInNormBin[nJetsBin], "kernelFit_sTRange")
    # rooKernel_extendedPDF_Fits[nJetsBin].fitTo(weighted_sTRooDataSets[nJetsBin], normRange, ROOT.RooFit.Minos(ROOT.kTRUE), ROOT.RooFit.PrintLevel(0))
    # rooVar_nEventsInNormBin[nJetsBin].Print()
    # rooKernel_extendedPDF_Fits[nJetsBin].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.LineColor(ROOT.kBlue), plotRange, ROOT.RooFit.Normalization(fractionalUncertainties_nEvents_normRange_factors_up[nJetsBin], ROOT.RooAbsReal.Relative), ROOT.RooFit.LineStyle(ROOT.kDashed))
    # rooKernel_extendedPDF_Fits[nJetsBin].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.LineColor(ROOT.kBlue), plotRange, ROOT.RooFit.Normalization(fractionalUncertainties_nEvents_normRange_factors_down[nJetsBin], ROOT.RooAbsReal.Relative), ROOT.RooFit.LineStyle(ROOT.kDashed))
    # rooKernel_extendedPDF_Fits[nJetsBin].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.LineColor(ROOT.kBlue), plotRange, ROOT.RooFit.Normalization(1.0, ROOT.RooAbsReal.Relative))

    if not(inputArguments.isSignal): dataSystematicsList.append(tuple(["float", "fractionalUncertainty_sTScaling_{nJetsBin}Jets".format(nJetsBin=nJetsBin), fractionalUncertainty_sTScaling]))

    if (nJetsBin == inputArguments.nJetsMax): setFrameAesthetics(sTFrames["data"][nJetsBin], "#it{S}_{T} (GeV)", "Events / ({binWidth} GeV)".format(binWidth=binWidth), "#geq {nJetsBin} Jets".format(nJetsBin=nJetsBin))
    else: setFrameAesthetics(sTFrames["data"][nJetsBin], "#it{S}_{T} (GeV)", "Events / ({binWidth} GeV)".format(binWidth=binWidth), "{nJetsBin} Jets".format(nJetsBin=nJetsBin))
    canvases["data"][nJetsBin] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["data"][nJetsBin]], canvasName = "c_kernelPDF_{nJetsBin}Jets".format(nJetsBin=nJetsBin), outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_kernelPDF_{nJetsBin}Jets".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix, nJetsBin=nJetsBin), enableLogY = True)

if not(inputArguments.isSignal):
    tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=dataSystematicsList, outputFilePath=("{outputDirectory}/{outputPrefix}_dataSystematics_sTScaling.dat".format(outputDirectory=inputArguments.outputDirectory_dataSystematics, outputPrefix=inputArguments.outputPrefix)))
    # No need to go further if the region is the control region
    outputFile.Write()
    outputFile.Close()
    sys.exit(0)

# Generate and fit toy MC datsets with the fit kernels
toyRooDataSets = {}
sTFrames["toyMC"]["DataAndFits"] = rooVar_sT.frame(STNormRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
predictedToObservedNEventsHistogram = ROOT.TH1F("h_predictedToObservedNEvents", "predicted/observed;fraction;Toy MC events", 40, 0., 0.)
totalIntegralCheckHistogram = ROOT.TH1F("h_totalIntegralCheck", "Total integral;total integral;Toy MC events", 40, 0., 0.)
goodMCSampleIndex = 0
randomGenerator = ROOT.TRandom1()
randomGenerator.SetSeed(0) # Sets seed by using some information from a ROOT "UUID"
progressBar = tmProgressBar(inputArguments.nToyMCs)
progressBarUpdatePeriod = max(1, inputArguments.nToyMCs//1000)
progressBar.initializeTimer()
while goodMCSampleIndex < inputArguments.nToyMCs:
    nEventsToGenerate = total_nEventsInFullRange[inputArguments.nJetsNorm]
    toyRooDataSets[goodMCSampleIndex] = ROOT.RooDataSet("toyMCDataSet_index{goodMCSampleIndex}".format(goodMCSampleIndex=goodMCSampleIndex), "toyMCDataSet_index{goodMCSampleIndex}".format(goodMCSampleIndex=goodMCSampleIndex), ROOT.RooArgSet(rooVar_sT))
    rooVar_sT.setRange(STNormRangeMin, STNormRangeMax)
    nNormEventsToGenerate = nEventsInNormWindows[inputArguments.nJetsNorm]
    if (inputArguments.varyNEventsInNormWindowInToyMCs): nNormEventsToGenerate = randomGenerator.Poisson(nNormEventsToGenerate)
    dataSet_normWindow = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].generate(ROOT.RooArgSet(rooVar_sT), nNormEventsToGenerate)
    toyRooDataSets[goodMCSampleIndex].append(dataSet_normWindow)
    rooVar_sT.setRange(STNormRangeMax, inputArguments.sTKernelFitRangeMax)
    nObsEventsToGenerate = nEventsToGenerate - nNormEventsToGenerate
    if (inputArguments.varyNEventsInObservationWindowInToyMCs): nObsEventsToGenerate = randomGenerator.Poisson(nObsEventsToGenerate)
    dataSet_observationWindow = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].generate(ROOT.RooArgSet(rooVar_sT), nObsEventsToGenerate)
    toyRooDataSets[goodMCSampleIndex].append(dataSet_observationWindow)
    rooVar_sT.setRange(STNormRangeMin, inputArguments.sTKernelFitRangeMax)
    nToyEventsInNormWindow = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(toyRooDataSets[goodMCSampleIndex], "normalization_sTRange", forceNumEntries=True)
    nToyEventsInObservationWindow = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(toyRooDataSets[goodMCSampleIndex], "observation_sTRange", forceNumEntries=True)
    throwAwayEvent = False
    if (not(inputArguments.varyNEventsInNormWindowInToyMCs) and not(nToyEventsInNormWindow == nNormEventsToGenerate)):
        if (abs(nToyEventsInNormWindow - nNormEventsToGenerate) == 1):
            throwAwayEvent = True
            print("WARNING: incorrect data generation: nToyEventsInNormWindow = {n1}, nEventsInNormWindows[inputArguments.nJetsNorm] = {n2}".format(n1=nToyEventsInNormWindow, n2=nEventsInNormWindows[inputArguments.nJetsNorm]))
        else:
            sys.exit("ERROR: Wildly incorrect data generation: nToyEventsInNormWindow = {n1}, nEventsInNormWindows[inputArguments.nJetsNorm] = {n2}".format(n1=nToyEventsInNormWindow, n2=nEventsInNormWindows[inputArguments.nJetsNorm]))
    if (not(inputArguments.varyNEventsInObservationWindowInToyMCs) and not(nToyEventsInObservationWindow == nObsEventsToGenerate)):
        if (abs(nToyEventsInObservationWindow - nObsEventsToGenerate) == 1):
            throwAwayEvent = True
            print("WARNING: incorrect data generation: nToyEventsInObservationWindow = {n1}, nEventsInObservationWindows[inputArguments.nJetsNorm] = {n2}".format(n1=nToyEventsInObservationWindow, n2=nEventsInObservationWindows[inputArguments.nJetsNorm]))
        else:
            sys.exit("ERROR: Wildly incorrect data generation: nToyEventsInObservationWindow = {n1}, nEventsInObservationWindows[inputArguments.nJetsNorm] = {n2}".format(n1=nToyEventsInObservationWindow, n2=nEventsInObservationWindows[inputArguments.nJetsNorm]))
    if throwAwayEvent: continue
    toyRooDataSets[goodMCSampleIndex].plotOn(sTFrames["toyMC"]["DataAndFits"])
    rooKernel_PDF_Fits["toyMC"][goodMCSampleIndex] = ROOT.RooKeysPdf("toyMCKernelEstimateFunction_{index}".format(index=goodMCSampleIndex), "toyMCKernelEstimateFunction_{index}".format(index=goodMCSampleIndex), rooVar_sT, toyRooDataSets[goodMCSampleIndex], kernelOptionsObjects[inputArguments.kernelMirrorOption], inputArguments.nominalRho)
    rooKernel_PDF_Fits["toyMC"][goodMCSampleIndex].plotOn(sTFrames["toyMC"]["DataAndFits"], plotRange, ROOT.RooFit.NormRange("kernelFit_sTRange"), ROOT.RooFit.Normalization(1.0, ROOT.RooAbsReal.Relative))
    integralObject_observationRange = rooKernel_PDF_Fits["toyMC"][goodMCSampleIndex].createIntegral(ROOT.RooArgSet(rooVar_sT), "observation_sTRange")
    integralObject_normRange = rooKernel_PDF_Fits["toyMC"][goodMCSampleIndex].createIntegral(ROOT.RooArgSet(rooVar_sT), "normalization_sTRange")
    integralsRatio = 1.0*integralObject_observationRange.getVal() / integralObject_normRange.getVal()
    predicted_nEvents = integralsRatio*nToyEventsInNormWindow
    observed_nEvents = nToyEventsInObservationWindow
    ratio_predictedToObservedNEvents = 1.0*predicted_nEvents/observed_nEvents
    predictedToObservedNEventsHistogram.Fill(ratio_predictedToObservedNEvents)
    totalIntegralCheckObject = rooKernel_PDF_Fits["toyMC"][goodMCSampleIndex].createIntegral(ROOT.RooArgSet(rooVar_sT), "kernelFit_sTRange")
    totalIntegralCheckHistogram.Fill(totalIntegralCheckObject.getVal())
    if goodMCSampleIndex%progressBarUpdatePeriod == 0: progressBar.updateBar(1.0*goodMCSampleIndex/inputArguments.nToyMCs, goodMCSampleIndex)
    goodMCSampleIndex += 1
progressBar.terminate()
sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["toyMC"]["DataAndFits"], ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.RefreshNorm())
rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].plotOn(sTFrames["toyMC"]["DataAndFits"], ROOT.RooFit.LineColor(ROOT.kRed), plotRange, ROOT.RooFit.NormRange("kernelFit_sTRange"))

# Plot the toy MC data and fits
setFrameAesthetics(sTFrames["toyMC"]["DataAndFits"], "#it{S}_{T} (GeV)", "Events / ({binWidth} GeV)".format(binWidth=binWidth), "Data and fits: {n} Toy MCs".format(n = inputArguments.nToyMCs))
canvases["toyMC"]["DataAndFits"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["toyMC"]["DataAndFits"]], canvasName = "c_toyMCDataAndKernelEstimates", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_toyMCDataAndKernelEstimates".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix), enableLogY = True)

# Estimate of the uncertainty = RMS of the distribution of the predicted to observed number of events
fractionalUncertainty_Shape = predictedToObservedNEventsHistogram.GetRMS()
dataSystematicsList.append(tuple(["float", "fractionalUncertainty_Shape", fractionalUncertainty_Shape]))

# Plot the shape systematics estimate
canvases["toyMC"]["systematics"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [predictedToObservedNEventsHistogram], canvasName = "c_shapeSystematics", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_shapeSystematics".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix))

# Plot the integral checks, to see that the normalization is similar
canvases["toyMC"]["systematicsCheck"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [totalIntegralCheckHistogram], canvasName = "c_shapeSystematicsCheck", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_shapeSystematicsCheck".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix))

rhoSystematicsGraph = ROOT.TGraph()
rhoSystematicsGraph.SetName("rhoSystematics")
rhoSystematicsGraph.SetTitle("predicted/observed;rho;")
rhoChi2Graph = ROOT.TGraph()
rhoChi2Graph.SetName("rhoChi2")
rhoChi2Graph.SetTitle("Chi2;rho;")
progressBar = tmProgressBar(inputArguments.nRhoValuesForSystematicsEstimation)
progressBarUpdatePeriod = max(1, inputArguments.nRhoValuesForSystematicsEstimation//1000)
fractionalUncertainty_rhoHigh = 0
deviationAtRhoHigh_rhoLow = 0
rhoMinForSystematicsEstimation = inputArguments.rhoMinFactorForSystematicsEstimation * inputArguments.nominalRho
rhoMaxForSystematicsEstimation = inputArguments.rhoMaxFactorForSystematicsEstimation * inputArguments.nominalRho
histogrammedNormJetsDataSet = (weighted_sTRooDataSets[inputArguments.nJetsNorm]).binnedClone("h_" + (weighted_sTRooDataSets[inputArguments.nJetsNorm]).GetName(), "Binned " + (weighted_sTRooDataSets[inputArguments.nJetsNorm]).GetTitle())
progressBar.initializeTimer()
for rhoCounter in range(0, inputArguments.nRhoValuesForSystematicsEstimation):
    if rhoCounter%progressBarUpdatePeriod == 0: progressBar.updateBar(1.0*rhoCounter/inputArguments.nRhoValuesForSystematicsEstimation, rhoCounter)
    rhoValue = rhoMinForSystematicsEstimation + (rhoCounter/(inputArguments.nRhoValuesForSystematicsEstimation-1))*(rhoMaxForSystematicsEstimation - rhoMinForSystematicsEstimation)
    rooKernel_PDF_Fits["rhoValues"][rhoCounter] = ROOT.RooKeysPdf("normBinKernelEstimateFunction_rhoCounter_{rhoC}".format(rhoC=rhoCounter), "normBinKernelEstimateFunction_rhoCounter_{rhoC}".format(rhoC=rhoCounter), rooVar_sT, weighted_sTRooDataSets[inputArguments.nJetsNorm], kernelOptionsObjects[inputArguments.kernelMirrorOption], rhoValue)
    chi2 = (rooKernel_PDF_Fits["rhoValues"][rhoCounter]).createChi2(histogrammedNormJetsDataSet, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    rhoChi2Graph.SetPoint(rhoSystematicsGraph.GetN(), rhoValue, chi2.getVal())
    integralObject_rhoValues_observationRange = rooKernel_PDF_Fits["rhoValues"][rhoCounter].createIntegral(ROOT.RooArgSet(rooVar_sT), "observation_sTRange")
    integralObject_rhoValues_normalizationRange = rooKernel_PDF_Fits["rhoValues"][rhoCounter].createIntegral(ROOT.RooArgSet(rooVar_sT), "normalization_sTRange")
    integralsRatio_rhoValues_normJets = integralObject_rhoValues_observationRange.getVal()/integralObject_rhoValues_normalizationRange.getVal()
    predicted_nEvents_rhoValues = integralsRatio_rhoValues_normJets*weighted_nEventsInNormWindows[inputArguments.nJetsNorm]
    predictedToObserved_atThisRho = predicted_nEvents_rhoValues/weighted_nEventsInObservationWindows[inputArguments.nJetsNorm]
    rhoSystematicsGraph.SetPoint(rhoSystematicsGraph.GetN(), rhoValue, predictedToObserved_atThisRho)
    if (rhoCounter == 0): fractionalUncertainty_rhoLow = abs((predictedToObserved_atThisRho/predictedToObservedRatio_nominalRho_observationRegion) - 1)
    if (rhoCounter == (-1 + inputArguments.nRhoValuesForSystematicsEstimation)): fractionalUncertainty_rhoHigh = abs((predictedToObserved_atThisRho/predictedToObservedRatio_nominalRho_observationRegion) - 1)
progressBar.terminate()
canvases["rhoValues"]["ratiosGraph"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [rhoSystematicsGraph], canvasName = "c_rhoSystematicsGraph", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_rhoSystematics".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix))
canvases["rhoValues"]["chi2Graph"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [rhoChi2Graph], canvasName = "c_rhoChi2Graph", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_rhoChi2".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix))

fractionalUncertainty_rho = 0.5*(fractionalUncertainty_rhoLow + fractionalUncertainty_rhoHigh)
dataSystematicsList.append(tuple(["float", "fractionalUncertainty_rho", fractionalUncertainty_rho]))

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
rooVar_sT.setRange(STNormRangeMin, inputArguments.sTKernelFitRangeMax)
sTFrames["rhoValues"]["dataAndKernels"] = rooVar_sT.frame(STNormRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
setFrameAesthetics(sTFrames["rhoValues"]["dataAndKernels"], "#it{S}_{T} (GeV)", "Events / ({binWidth} GeV)".format(binWidth=binWidth), "Kernel estimates, nJets = {nJets}".format(nJets=inputArguments.nJetsNorm))
weighted_sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.RefreshNorm())
dataSet_legendEntry = dataAndKernelsLegend.AddEntry(weighted_sTRooDataSets[inputArguments.nJetsNorm], "Data")
rooKernel_PDF_Fits["rhoValues"][0].plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.LineColor(ROOT.kBlue), plotRange)
rhoMin_legendEntry = dataAndKernelsLegend.AddEntry(rooKernel_PDF_Fits["rhoValues"][0], "PDF fit: rho = {rhoMin:.1f}".format(rhoMin=rhoMinForSystematicsEstimation))
customizeLegendEntryForLine(rhoMin_legendEntry, ROOT.kBlue)
rooKernel_PDF_Fits["rhoValues"][-1+inputArguments.nRhoValuesForSystematicsEstimation].plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.LineColor(ROOT.kRed), plotRange)
rhoMax_legendEntry = dataAndKernelsLegend.AddEntry(rooKernel_PDF_Fits["rhoValues"][-1+inputArguments.nRhoValuesForSystematicsEstimation], "PDF fit: rho = {rhoMax:.1f}".format(rhoMax=rhoMaxForSystematicsEstimation))
customizeLegendEntryForLine(rhoMax_legendEntry, ROOT.kRed)
rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.LineColor(ROOT.kBlack), plotRange)
rho1_legendEntry = dataAndKernelsLegend.AddEntry(rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], "PDF fit: rho = {rhoNom:.1f}".format(rhoNom=inputArguments.nominalRho))
customizeLegendEntryForLine(rho1_legendEntry, ROOT.kBlack)
canvases["rhoValues"]["dataAndKernels"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["rhoValues"]["dataAndKernels"], dataAndKernelsLegend], canvasName = "c_kernelPDFEstimate_dataAndKernels", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_kernelPDF_rhoValues".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix), enableLogY = False)

for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    sTTrees[nJetsBin].Write()

outputFile.Write()
outputFile.Close()

tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=dataSystematicsList, outputFilePath=("{outputDirectory}/{outputPrefix}_dataSystematics.dat".format(outputDirectory=inputArguments.outputDirectory_dataSystematics, outputPrefix=inputArguments.outputPrefix)))
tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=expectedEventCountersList, outputFilePath=("{outputDirectory}/{outputPrefix}_eventCounters.dat".format(outputDirectory=inputArguments.outputDirectory_dataSystematics, outputPrefix=inputArguments.outputPrefix)))

print("All done!")
