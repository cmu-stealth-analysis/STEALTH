#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, ROOT, argparse, array, pdb, math
import numpy as np
import tmROOTUtils, tmStatsUtils, tmGeneralUtils
from tmProgressBar import tmProgressBar

inputArgumentsParser = argparse.ArgumentParser(description='Get data event histograms and systematics.')
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_eventHistograms', default="analysis/dataEventHistograms/", help='Directory in which to store event histograms.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_dataSystematics', default="analysis/dataSystematics/", help='Directory in which to store data systematics.',type=str)
inputArgumentsParser.add_argument('--outputPrefix', required=True, help='Prefix to all output file names.',type=str)
inputArgumentsParser.add_argument('--sTPlotRangeMin', default=1200., help='Min value of sT to display in the plots.',type=float)
inputArgumentsParser.add_argument('--sTPlotRangeMax', default=3500., help='Max value of sT to display in the plots.',type=float)
inputArgumentsParser.add_argument('--n_sTBins', default=23, help='Number of sT bins (relevant for plotting only).',type=int)
inputArgumentsParser.add_argument('--sTKernelFitRangeMin', default=700., help='Min value of sT to use in the kernel fit. This should be slightly less than sTPlotRangeMin to try to get rid of boundary effects.',type=float)
inputArgumentsParser.add_argument('--sTKernelFitRangeMax', default=3500., help='Max value of sT to use in the kernel fit.',type=float)
inputArgumentsParser.add_argument('--sTNormRangeMin', default=1200., help='Min value of sT for normalization. For all sT distributions in nJets bins except the normalization bin, this value is the min of the range in which to scale the kernel estimator.',type=float)
inputArgumentsParser.add_argument('--sTNormRangeMax', default=1300., help='Max value of sT for normalization. For all sT distributions in nJets bins except the normalization bin, this value is the max of the range in which to scale the kernel estimator.',type=float)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Min number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Max number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsNorm', default=2, help='Number of jets w.r.t. which to normalize the sT distributions for other jets.',type=int)
inputArgumentsParser.add_argument('--nTargetEventsForSTThresholdOptimization', default=1., help='Number of jets w.r.t. which to normalize the sT distributions for other jets.',type=float)
inputArgumentsParser.add_argument('--nToyMCs', default=1000, help='Number of toy MC samples to generate using the pdf fits found.',type=int)
inputArgumentsParser.add_argument('--varyNEventsInPreNormWindowInToyMCs', action='store_true', help="Vary the number of generated events in toy MC samples in the pre-normalization window in a Poisson distribution about the number of events in this window in the original sample; default is to keep it fixed.")
inputArgumentsParser.add_argument('--varyNEventsInNormWindowInToyMCs', action='store_true', help="Vary the number of generated events in the toy MC samples in the normalization window in a Poisson distribution about the number of events in this window in the original sample; default is to keep it fixed.")
inputArgumentsParser.add_argument('--varyNEventsInObservationWindowInToyMCs', action='store_true', help="Vary the number of generated events in the toy MC samples in the observation window in a Poisson distribution about the number of events in this window in the original sample; default is to keep it fixed.")
inputArgumentsParser.add_argument('--rho', default=1., help='Value of parameter rho to be used in adaptive Gaussian kernel estimates.',type=float)
inputArgumentsParser.add_argument('--rhoMinForSystematicsEstimation', default=0.5, help='Lower value of rho to use in systematics estimation.',type=float)
inputArgumentsParser.add_argument('--rhoMaxForSystematicsEstimation', default=2.0, help='Upper value of rho to use in systematics estimation.',type=float)
inputArgumentsParser.add_argument('--nRhoValuesForSystematicsEstimation', default=151, help='Number of values of rho to use in systematics estimation.',type=int)
inputArgumentsParser.add_argument('--kernelMirrorOption', default="MirrorLeft", help='Kernel mirroring option to be used in adaptive Gaussian kernel estimates',type=str)
inputArgumentsParser.add_argument('--allowHigherNJets', action='store_true', help="Allow script to look into nJets bins beyond nJets = 3. Do not use with data with the signal selections before unblinding.")
inputArgumentsParser.add_argument('--isSignal', action='store_true', help="If this flag is set, then the input file is treated as belonging to the signal region; in that case, do not compute systematics on the degree to which sT scales, but compute all other systematics. If this flag is not set, then the input file is treated as belonging to the control region; in that case, compute the systematics estimate on the degree to which sT scales, but do not compute the other systematics.")
inputArgumentsParser.add_argument('--sTStartMainRegion', default=2500., help="Value of sT separating main signal region from subordinate region.", type=float)
inputArguments = inputArgumentsParser.parse_args()
if (inputArguments.sTNormRangeMin < inputArguments.sTKernelFitRangeMin or inputArguments.sTNormRangeMax > inputArguments.sTKernelFitRangeMax):
    sys.exit("Normalization interval: ({nmin}, {nmax}) seems incompatible with kernel fitting range: ({smin, smax})".format(nmin=inputArguments.sTNormRangeMin, nmax=inputArguments.sTNormRangeMax, smin=inputArguments.sTKernelFitRangeMin, smax=inputArguments.sTKernelFitRangeMax))
kernelOptionsObjects = {"NoMirror": ROOT.RooKeysPdf.NoMirror, "MirrorLeft": ROOT.RooKeysPdf.MirrorLeft, "MirrorRight": ROOT.RooKeysPdf.MirrorRight, "MirrorBoth": ROOT.RooKeysPdf.MirrorBoth, "MirrorAsymLeft": ROOT.RooKeysPdf.MirrorAsymLeft, "MirrorAsymLeftRight": ROOT.RooKeysPdf.MirrorAsymLeftRight, "MirrorAsymRight": ROOT.RooKeysPdf.MirrorAsymRight, "MirrorLeftAsymRight": ROOT.RooKeysPdf.MirrorLeftAsymRight, "MirrorAsymBoth": ROOT.RooKeysPdf.MirrorAsymBoth}
if not(inputArguments.kernelMirrorOption in kernelOptionsObjects): sys.exit("The following element is passed as an argument for the kernel mirroring option but not in the dictionary defining the correspondence between kernel name and RooKeysPdf index: {kernelMirrorOption}".format(kernelMirrorOption=inputArguments.kernelMirrorOption))
# if (inputArguments.generateDataCardTemplate and ((inputArguments.sTScalingFractionalUncertainty_4Jets < 0. or inputArguments.sTScalingFractionalUncertainty_5Jets < 0.) or inputArguments.sTScalingFractionalUncertainty_geq6Jets < 0.)): sys.exit("If a data card template needs to be generated then the sT scaling uncertainties in all nJets bins must be passed as an argument.")

if not(inputArguments.nJetsMax == 6): sys.exit("Only nJetsMax=6 supported temporarily. Needed to create data template in the correct format.")

rooVar_sT = ROOT.RooRealVar("rooVar_sT", "rooVar_sT", inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax, "GeV")
rooVar_sT.setRange("preNormalization_sTRange", inputArguments.sTKernelFitRangeMin, inputArguments.sTNormRangeMin)
rooVar_sT.setRange("normalization_sTRange", inputArguments.sTNormRangeMin, inputArguments.sTNormRangeMax)
rooVar_sT.setRange("observation_sTRange", inputArguments.sTNormRangeMax, inputArguments.sTKernelFitRangeMax)
rooVar_sT.setRange("subordinateSignal_sTRange", inputArguments.sTNormRangeMax, inputArguments.sTStartMainRegion)
rooVar_sT.setRange("mainSignal_sTRange", inputArguments.sTStartMainRegion, inputArguments.sTKernelFitRangeMax)
rooVar_sT.setRange("kernelFit_sTRange", inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
plotRange = ROOT.RooFit.Range(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax)
kernelFitRange = ROOT.RooFit.Range(inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
normRange = ROOT.RooFit.Range(inputArguments.sTNormRangeMin, inputArguments.sTNormRangeMax)
binWidth = int(0.5 + (inputArguments.sTPlotRangeMax - inputArguments.sTPlotRangeMin)/inputArguments.n_sTBins)

# plotRangeString = "plotRange_{plotMin:2.1f}_{plotMax:2.1f}_".format(plotMin=inputArguments.sTPlotRangeMin, plotMax=inputArguments.sTPlotRangeMax)
# sTBinsString = "{n}sTBins_".format(n=inputArguments.n_sTBins)
# kernelFitRangeString = "kernelFitRange_{kMin:2.1f}_{kMax:2.1f}_".format(kMin=inputArguments.sTKernelFitRangeMin, kMax=inputArguments.sTKernelFitRangeMax)
# normRangeString = "normRange_{normMin:2.1f}_{normMax:2.1f}_".format(normMin=inputArguments.sTNormRangeMin, normMax=inputArguments.sTNormRangeMax)
# nJetsMaxString = ""
# if not(inputArguments.nJetsMax == 6): nJetsMaxString = "{n}JetsMax_".format(n=inputArguments.nJetsMax)
# nJetsNormString = ""
# if not(inputArguments.nJetsNorm == 3): nJetsNormString = "{n}JetsNorm_".format(n=inputArguments.nJetsNorm)
# nToyMCsString = ""
# if not(inputArguments.nToyMCs == 1000): nToyMCsString = "{n}ToyMCs_".format(n=inputArguments.nToyMCs)
# preNormEventsString = ""
# if (inputArguments.varyNEventsInPreNormWindowInToyMCs): preNormEventsString = "variableNEventsPreNorm_"
# normEventsString = ""
# if (inputArguments.varyNEventsInNormWindowInToyMCs): normEventsString = "variableNEventsNorm_"
# obsEventsString = ""
# if (inputArguments.varyNEventsInObservationWindowInToyMCs): obsEventsString = "variableNEventsObs_"
# rhoString = ""
# if not(inputArguments.rho == 1.): rhoString = "rho{rho:2.1f}_".format(rho=inputArguments.rho)
# kernelMirrorOptionString = ""
# if not(inputArguments.kernelMirrorOption == "MirrorLeft"): kernelMirrorOptionString = "option{opt}_".format(opt=inputArguments.kernelMirrorOption)
# sTStartMainRegionString = ""
# if not(inputArguments.sTStartMainRegion == 2500.): sTStartMainRegionString = "sTStartMainRegion_{ststart:.1f}_".format(ststart=inputArguments.sTStartMainRegion)
# concatenatedString = ("{plotRangeString}{sTBinsString}{kernelFitRangeString}{normRangeString}{nJetsMaxString}{nJetsNormString}{nToyMCsString}{preNormEventsString}{normEventsString}{obsEventsString}{rhoString}{kernelMirrorOptionString}{sTStartMainRegionString}".format(**locals())).rstrip('_')
# outputDirectoryName = ("{outputDirectoryPrefix}_{concatenatedString}".format(outputDirectoryPrefix=inputArguments.outputDirectoryPrefix, concatenatedString=concatenatedString)).replace('.', 'pt')
# if (len(outputDirectoryName) > 255): sys.exit("Length of directory name should be no more than 255 characters long. Current name: {currentName}".format(currentName=outputDirectoryName))
# outputDirectoryPath = "analysis/{outputDirectoryName}".format(outputDirectoryName=outputDirectoryName)

# if not(os.path.isdir(outputDirectoryPath)): os.system("mkdir -p {outputDirectoryPath}".format(outputDirectoryPath=outputDirectoryPath))

dataSystematicsList = []
eventCountersList = []

def setFrameAesthetics(frame, xLabel, yLabel, title):
    frame.SetXTitle(xLabel)
    frame.SetYTitle(yLabel)
    frame.SetTitle(title)

def getExpectedNEventsFromPDFInNamedRange(normFactor=None, inputRooPDF=None, inputRooArgSet=None, targetRangeMin=None, targetRangeMax=None):
    targetRangeName = "targetRange_sTOptimization"
    rooVar_sT.setRange(targetRangeName, targetRangeMin, targetRangeMax)
    integralObject_targetRange = inputRooPDF.createIntegral(inputRooArgSet, targetRangeName)
    expectedNEvents = normFactor*integralObject_targetRange.getVal()
    rooVar_sT.setRange(inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
    return expectedNEvents

sw = ROOT.TStopwatch()
sw.Start()

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

nEventsInSubordinateSignalWindows = {}
nEventsInMainSignalWindows = {}
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    nEventsInSubordinateSignalWindows[nJetsBin] = 0
    nEventsInMainSignalWindows[nJetsBin] = 0
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

    nStealthJets = inputChain.b_nJets
    nJetsBin = nStealthJets
    if (nJetsBin > inputArguments.nJetsMax): nJetsBin = inputArguments.nJetsMax
    if (nJetsBin < inputArguments.nJetsMin): sys.exit("Unexpected nJetsBin = {nJetsBin} at entry index = {entryIndex}".format(nJetsBin=nJetsBin, entryIndex=entryIndex))
    sT = inputChain.b_evtST

    if (sT > inputArguments.sTNormRangeMax and sT < inputArguments.sTStartMainRegion): nEventsInSubordinateSignalWindows[nJetsBin] += 1
    if (sT > inputArguments.sTStartMainRegion): nEventsInMainSignalWindows[nJetsBin] += 1

    if (sT > inputArguments.sTKernelFitRangeMin and sT < inputArguments.sTKernelFitRangeMax):
        (sTArrays[nJetsBin])[0] = sT
        (sTTrees[nJetsBin]).Fill()
progressBar.terminate()

# Make datasets from all sT trees
sTRooDataSets = {}
nEventsInPreNormWindows = {}
nEventsInNormWindows = {}
nEventsInObservationWindows = {}
total_nEventsInFullRange = {}
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    sTRooDataSetName = "rooDataSet_{nJets}Jets".format(nJets=nJetsBin)
    sTRooDataSets[nJetsBin] = ROOT.RooDataSet(sTRooDataSetName, sTRooDataSetName, sTTrees[nJetsBin], ROOT.RooArgSet(rooVar_sT))
    nEventsInPreNormWindows[nJetsBin] = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(sTRooDataSets[nJetsBin], "preNormalization_sTRange")
    nEventsInNormWindows[nJetsBin] = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(sTRooDataSets[nJetsBin], "normalization_sTRange")
    nEventsInObservationWindows[nJetsBin] = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(sTRooDataSets[nJetsBin], "observation_sTRange")
    total_nEventsInFullRange[nJetsBin] = nEventsInPreNormWindows[nJetsBin] + nEventsInNormWindows[nJetsBin] + nEventsInObservationWindows[nJetsBin]
    if (nJetsBin <= 3 or (nJetsBin > 3 and inputArguments.allowHigherNJets)): print("At nJets = {nJets}, nEventsInPreNormWindow = {preNorm}, nEventsInNormWindow = {norm}, nEventsInObservationWindow = {obs}, nEventsInSurbordinateSignalWindow={sub}, nEventsInMainSignalWindow={mainsig}".format(nJets = nJetsBin, preNorm = nEventsInPreNormWindows[nJetsBin], norm = nEventsInNormWindows[nJetsBin], obs = nEventsInObservationWindows[nJetsBin], sub = nEventsInSubordinateSignalWindows[nJetsBin], mainsig=nEventsInMainSignalWindows[nJetsBin]))

poissonConfidenceInterval_normEvents = tmROOTUtils.getPoissonConfidenceInterval(observedNEvents=nEventsInNormWindows[inputArguments.nJetsNorm])
fractionalUncertainty_normEvents = (poissonConfidenceInterval_normEvents["upper"] - poissonConfidenceInterval_normEvents["lower"])/(2*nEventsInNormWindows[inputArguments.nJetsNorm])
fractionalUncertainty_normEvents_normalizationFactor_Up = (poissonConfidenceInterval_normEvents["upper"])/(nEventsInNormWindows[inputArguments.nJetsNorm])
fractionalUncertainty_normEvents_normalizationFactor_Down = (poissonConfidenceInterval_normEvents["lower"])/(nEventsInNormWindows[inputArguments.nJetsNorm])
if (inputArguments.isSignal):
    dataSystematicsList.append(tuple(["float", "fractionalUncertainty_normEvents", fractionalUncertainty_normEvents]))
    print("Fractional uncertainties from Poisson errors on number of events in normalization bin: {a:.3f}".format(a=fractionalUncertainty_normEvents))

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
rooVar_sT.setRange(inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm] = ROOT.RooKeysPdf("normBinKernelEstimateFunction", "normBinKernelEstimateFunction", rooVar_sT, sTRooDataSets[inputArguments.nJetsNorm], kernelOptionsObjects[inputArguments.kernelMirrorOption], inputArguments.rho)
integralObject_normalizationRange = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].createIntegral(ROOT.RooArgSet(rooVar_sT), "normalization_sTRange")
backgroundFit_normRange_integralValue = integralObject_normalizationRange.getVal()
normFactor_checks = nEventsInNormWindows[inputArguments.nJetsNorm]/backgroundFit_normRange_integralValue
expected_nEvents_normJetsBin_observationRegion = getExpectedNEventsFromPDFInNamedRange(normFactor=normFactor_checks, inputRooPDF=rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], inputRooArgSet=ROOT.RooArgSet(rooVar_sT), targetRangeMin=inputArguments.sTNormRangeMax, targetRangeMax=inputArguments.sTKernelFitRangeMax)
observed_nEvents_normJetsBin_observationRegion = nEventsInObservationWindows[inputArguments.nJetsNorm]
predictedToObservedRatio_nominalRho_observationRegion = expected_nEvents_normJetsBin_observationRegion/observed_nEvents_normJetsBin_observationRegion
print("Check 1 on norm factor: expected events in norm range: {nExpected}, observed: {nObserved}".format(nExpected=getExpectedNEventsFromPDFInNamedRange(normFactor=normFactor_checks, inputRooPDF=rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], inputRooArgSet=ROOT.RooArgSet(rooVar_sT), targetRangeMin=inputArguments.sTNormRangeMin, targetRangeMax=inputArguments.sTNormRangeMax), nObserved=nEventsInNormWindows[inputArguments.nJetsNorm]))
print("Check 2 on norm factor: expected events in observation range: {nExpected}, observed: {nObserved}".format(nExpected=expected_nEvents_normJetsBin_observationRegion, nObserved=observed_nEvents_normJetsBin_observationRegion))
rooVar_sT.setRange(inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
sTFrames["data"][inputArguments.nJetsNorm] = rooVar_sT.frame(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm])
rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], plotRange)
setFrameAesthetics(sTFrames["data"][inputArguments.nJetsNorm], "#it{S}_{T} (GeV)", "Events / ({binWidth} GeV)".format(binWidth=binWidth), "Normalization bin: {nJets} Jets".format(nJets=inputArguments.nJetsNorm))
canvases["data"][inputArguments.nJetsNorm] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["data"][inputArguments.nJetsNorm]], canvasName = "c_kernelPDF_normJetsBin", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_kernelPDF_normJetsBin".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix))

# If input sample is the control sample, then use these fits in other nJets bins and obtain estimate of systematic on assumption that sT scales; otherwise, only obtain a prediction for number of events in subordinate and main signal region
rooVar_nEventsInNormBin = {}
rooKernel_extendedPDF_Fits = {}
integralObject_normJets_normRange = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].createIntegral(ROOT.RooArgSet(rooVar_sT), "normalization_sTRange")
integralObject_normJets_observationRange = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].createIntegral(ROOT.RooArgSet(rooVar_sT), "observation_sTRange")
integralsRatio_normJets = 1.0*integralObject_normJets_observationRange.getVal() / integralObject_normJets_normRange.getVal()

integralObject_normJets_subordinateSignalRange = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].createIntegral(ROOT.RooArgSet(rooVar_sT), "subordinateSignal_sTRange")
integralObject_normJets_mainSignalRange = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].createIntegral(ROOT.RooArgSet(rooVar_sT), "mainSignal_sTRange")
integralsRatio_normJets_subordinateSignalRange = integralObject_normJets_subordinateSignalRange.getVal() / integralObject_normJets_normRange.getVal()
integralsRatio_normJets_mainSignalRange = integralObject_normJets_mainSignalRange.getVal() / integralObject_normJets_normRange.getVal()

expected_nEventsInSubordinateSignalRegions = {}
expected_nEventsInMainSignalRegions = {}

if inputArguments.isSignal: # Calculate optimal sT threshold
    normFactor_optimization = nEventsInNormWindows[nJetsBin]/backgroundFit_normRange_integralValue
    print("Check 1 on norm factor for optimization: expected events in norm range: {nExpected}, observed: {nObserved}".format(nExpected=getExpectedNEventsFromPDFInNamedRange(normFactor=normFactor_optimization, inputRooPDF=rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], inputRooArgSet=ROOT.RooArgSet(rooVar_sT), targetRangeMin=inputArguments.sTNormRangeMin, targetRangeMax=inputArguments.sTNormRangeMax), nObserved=nEventsInNormWindows[nJetsBin]))
    # print("Check 2 on norm factor for optimization: expected events in observation range: {nExpected}, observed: {nObserved}".format(nExpected=getExpectedNEventsFromPDFInNamedRange(normFactor=normFactor_optimization, inputRooPDF=rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], inputRooArgSet=ROOT.RooArgSet(rooVar_sT), targetRangeMin=inputArguments.sTNormRangeMax, targetRangeMax=inputArguments.sTKernelFitRangeMax), nObserved=nEventsInObservationWindows[nJetsBin])) # Check 2 is disabled because we don't want to unblind ourselves yet
    optimal_sTThreshold = tmStatsUtils.getMonotonicFunctionApproximateZero(inputFunction=(lambda sTThreshold: (-1*inputArguments.nTargetEventsForSTThresholdOptimization)+getExpectedNEventsFromPDFInNamedRange(normFactor=normFactor_optimization, inputRooPDF=rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], inputRooArgSet=ROOT.RooArgSet(rooVar_sT), targetRangeMin=sTThreshold, targetRangeMax=inputArguments.sTKernelFitRangeMax)), xRange=[1.1*inputArguments.sTNormRangeMax, 0.99*inputArguments.sTKernelFitRangeMax], autoZeroTolerance=True, printDebug=True)
    optimalThresholdOutputFile = open("{outputDirectory}/{outputPrefix}_sTOptimalThreshold.dat".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix), 'w')
    optimalThresholdOutputFile.write("Optimal sT Threshold: {thr}\n".format(thr=optimal_sTThreshold))
    optimalThresholdOutputFile.close()

for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    expected_nEventsInSubordinateSignalRegions[nJetsBin] = integralsRatio_normJets_subordinateSignalRange*nEventsInNormWindows[nJetsBin]
    eventCountersList.append(tuple(["float", "expectedNEvents_subordinateRegion_{nJetsBin}Jets".format(nJetsBin=nJetsBin), (expected_nEventsInSubordinateSignalRegions[nJetsBin])]))
    expected_nEventsInMainSignalRegions[nJetsBin] = integralsRatio_normJets_mainSignalRange*nEventsInNormWindows[nJetsBin]
    eventCountersList.append(tuple(["float", "expectedNEvents_mainRegion_{nJetsBin}Jets".format(nJetsBin=nJetsBin), (expected_nEventsInMainSignalRegions[nJetsBin])]))
    if (nJetsBin == inputArguments.nJetsNorm): continue
    if (nJetsBin > 3 and not(inputArguments.allowHigherNJets)): continue
    sTFrames["data"][nJetsBin] = rooVar_sT.frame(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
    rooKernel_PDF_Fits["data"][nJetsBin] = ROOT.RooKeysPdf("kernelEstimate_{nJetsBin}Jets".format(nJetsBin=nJetsBin), "kernelEstimate_{nJetsBin}Jets".format(nJetsBin=nJetsBin), rooVar_sT, sTRooDataSets[nJetsBin], kernelOptionsObjects[inputArguments.kernelMirrorOption], inputArguments.rho)
    integralObject_observationRange = rooKernel_PDF_Fits["data"][nJetsBin].createIntegral(ROOT.RooArgSet(rooVar_sT), "observation_sTRange")
    integralObject_normRange = rooKernel_PDF_Fits["data"][nJetsBin].createIntegral(ROOT.RooArgSet(rooVar_sT), "normalization_sTRange")
    integralsRatio = 1.0 * integralObject_observationRange.getVal()/integralObject_normRange.getVal()
    sTRooDataSets[nJetsBin].plotOn(sTFrames["data"][nJetsBin])
    rooVar_nEventsInNormBin[nJetsBin] = ROOT.RooRealVar("rooVar_nEventsInNormBin_{nJetsBin}Jets".format(nJetsBin=nJetsBin), "rooVar_nEventsInNormBin_{nJetsBin}Jets".format(nJetsBin=nJetsBin), 100, 0, 10000)
    rooKernel_extendedPDF_Fits[nJetsBin] = ROOT.RooExtendPdf("extendedKernelPDF_{nJetsBin}Jets".format(nJetsBin=nJetsBin), "extendedKernelPDF_{nJetsBin}Jets".format(nJetsBin=nJetsBin), rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], rooVar_nEventsInNormBin[nJetsBin], "kernelFit_sTRange")
    rooKernel_extendedPDF_Fits[nJetsBin].fitTo(sTRooDataSets[nJetsBin], normRange, ROOT.RooFit.Minos(ROOT.kTRUE), ROOT.RooFit.PrintLevel(0))
    rooVar_nEventsInNormBin[nJetsBin].Print()
    rooKernel_extendedPDF_Fits[nJetsBin].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.LineColor(ROOT.kBlue), plotRange, ROOT.RooFit.Normalization(fractionalUncertainty_normEvents_normalizationFactor_Up, ROOT.RooAbsReal.Relative), ROOT.RooFit.LineStyle(ROOT.kDashed))
    rooKernel_extendedPDF_Fits[nJetsBin].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.LineColor(ROOT.kBlue), plotRange, ROOT.RooFit.Normalization(fractionalUncertainty_normEvents_normalizationFactor_Down, ROOT.RooAbsReal.Relative), ROOT.RooFit.LineStyle(ROOT.kDashed))
    rooKernel_extendedPDF_Fits[nJetsBin].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.LineColor(ROOT.kBlue), plotRange)
    fraction_predictedToActual = integralsRatio/integralsRatio_normJets
    fractionalUncertainty_sTScaling = abs(fraction_predictedToActual - 1.0)

    if not(inputArguments.isSignal): dataSystematicsList.append(tuple(["float", "fractionalUncertainty_sTScaling_{nJetsBin}Jets".format(nJetsBin=nJetsBin), fractionalUncertainty_sTScaling]))
    
    # scalingSystematicsOutputFile.write("{nJetsBin}    {fraction}\n".format(nJetsBin=nJetsBin, fraction=fraction_predictedToActual))
    if (nJetsBin == inputArguments.nJetsMax): setFrameAesthetics(sTFrames["data"][nJetsBin], "#it{S}_{T} (GeV)", "Events / ({binWidth} GeV)".format(binWidth=binWidth), "#geq {nJetsBin} Jets".format(nJetsBin=nJetsBin))
    else: setFrameAesthetics(sTFrames["data"][nJetsBin], "#it{S}_{T} (GeV)", "Events / ({binWidth} GeV)".format(binWidth=binWidth), "{nJetsBin} Jets".format(nJetsBin=nJetsBin))
    canvases["data"][nJetsBin] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["data"][nJetsBin]], canvasName = "c_kernelPDF_{nJetsBin}Jets".format(nJetsBin=nJetsBin), outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_kernelPDF_{nJetsBin}Jets".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix, nJetsBin=nJetsBin))

if not(inputArguments.isSignal):
    tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=dataSystematicsList, outputFilePath=("{outputDirectory}/{outputPrefix}_dataSystematics_sTScaling.dat".format(outputDirectory=inputArguments.outputDirectory_dataSystematics, outputPrefix=inputArguments.outputPrefix)))
    # No need to go further if the region is the control region
    outputFile.Write()
    outputFile.Close()
    sys.exit(0)

# Generate and fit toy MC datsets with the fit kernels
toyRooDataSets = {}
sTFrames["toyMC"]["DataAndFits"] = rooVar_sT.frame(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
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
    rooVar_sT.setRange(inputArguments.sTKernelFitRangeMin, inputArguments.sTNormRangeMin)
    nPreNormEventsToGenerate = nEventsInPreNormWindows[inputArguments.nJetsNorm]
    if (inputArguments.varyNEventsInPreNormWindowInToyMCs): nPreNormEventsToGenerate = randomGenerator.Poisson(nEventsInPreNormWindows[inputArguments.nJetsNorm])
    dataSet_preNormWindow = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].generate(ROOT.RooArgSet(rooVar_sT), nPreNormEventsToGenerate)
    toyRooDataSets[goodMCSampleIndex].append(dataSet_preNormWindow)
    rooVar_sT.setRange(inputArguments.sTNormRangeMin, inputArguments.sTNormRangeMax)
    nNormEventsToGenerate = nEventsInNormWindows[inputArguments.nJetsNorm]
    if (inputArguments.varyNEventsInNormWindowInToyMCs): nNormEventsToGenerate = randomGenerator.Poisson(nEventsInNormWindows[inputArguments.nJetsNorm])
    dataSet_normWindow = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].generate(ROOT.RooArgSet(rooVar_sT), nNormEventsToGenerate)
    toyRooDataSets[goodMCSampleIndex].append(dataSet_normWindow)
    rooVar_sT.setRange(inputArguments.sTNormRangeMax, inputArguments.sTKernelFitRangeMax)
    nObsEventsToGenerate = nEventsInObservationWindows[inputArguments.nJetsNorm]
    if (inputArguments.varyNEventsInObservationWindowInToyMCs): nObsEventsToGenerate = nEventsInObservationWindows[inputArguments.nJetsNorm]
    dataSet_observationWindow = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].generate(ROOT.RooArgSet(rooVar_sT), nObsEventsToGenerate)
    toyRooDataSets[goodMCSampleIndex].append(dataSet_observationWindow)
    rooVar_sT.setRange(inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
    nToyEventsInNormWindow = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(toyRooDataSets[goodMCSampleIndex], "normalization_sTRange")
    nToyEventsInObservationWindow = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(toyRooDataSets[goodMCSampleIndex], "observation_sTRange")
    throwAwayEvent = False
    if (not(inputArguments.varyNEventsInNormWindowInToyMCs) and not(nToyEventsInNormWindow == nEventsInNormWindows[inputArguments.nJetsNorm])):
        if (abs(nToyEventsInNormWindow - nEventsInNormWindows[inputArguments.nJetsNorm]) == 1):
            throwAwayEvent = True
            print("WARNING: incorrect data generation: nToyEventsInNormWindow = {n1}, nEventsInNormWindows[inputArguments.nJetsNorm] = {n2}".format(n1=nToyEventsInNormWindow, n2=nEventsInNormWindows[inputArguments.nJetsNorm]))
        else:
            sys.exit("ERROR: Wildly incorrect data generation: nToyEventsInNormWindow = {n1}, nEventsInNormWindows[inputArguments.nJetsNorm] = {n2}".format(n1=nToyEventsInNormWindow, n2=nEventsInNormWindows[inputArguments.nJetsNorm]))
    if (not(inputArguments.varyNEventsInObservationWindowInToyMCs) and not(nToyEventsInObservationWindow == nEventsInObservationWindows[inputArguments.nJetsNorm])):
        if (abs(nToyEventsInObservationWindow - nEventsInObservationWindows[inputArguments.nJetsNorm]) == 1):
            throwAwayEvent = True
            print("WARNING: incorrect data generation: nToyEventsInObservationWindow = {n1}, nEventsInObservationWindows[inputArguments.nJetsNorm] = {n2}".format(n1=nToyEventsInObservationWindow, n2=nEventsInObservationWindows[inputArguments.nJetsNorm]))
        else:
            sys.exit("ERROR: Wildly incorrect data generation: nToyEventsInObservationWindow = {n1}, nEventsInObservationWindows[inputArguments.nJetsNorm] = {n2}".format(n1=nToyEventsInObservationWindow, n2=nEventsInObservationWindows[inputArguments.nJetsNorm]))
    if throwAwayEvent: continue
    toyRooDataSets[goodMCSampleIndex].plotOn(sTFrames["toyMC"]["DataAndFits"])
    rooKernel_PDF_Fits["toyMC"][goodMCSampleIndex] = ROOT.RooKeysPdf("toyMCKernelEstimateFunction_{index}".format(index=goodMCSampleIndex), "toyMCKernelEstimateFunction_{index}".format(index=goodMCSampleIndex), rooVar_sT, toyRooDataSets[goodMCSampleIndex], kernelOptionsObjects[inputArguments.kernelMirrorOption], inputArguments.rho)
    rooKernel_PDF_Fits["toyMC"][goodMCSampleIndex].plotOn(sTFrames["toyMC"]["DataAndFits"])
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
rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].plotOn(sTFrames["toyMC"]["DataAndFits"], ROOT.RooFit.LineColor(ROOT.kRed), plotRange)

# Plot the toy MC data and fits
setFrameAesthetics(sTFrames["toyMC"]["DataAndFits"], "#it{S}_{T} (GeV)", "Events / ({binWidth} GeV)".format(binWidth=binWidth), "Data and fits: {n} Toy MCs".format(n = inputArguments.nToyMCs))
canvases["toyMC"]["DataAndFits"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["toyMC"]["DataAndFits"]], canvasName = "c_toyMCDataAndKernelEstimates", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_toyMCDataAndKernelEstimates".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix))

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
progressBar = tmProgressBar(inputArguments.nRhoValuesForSystematicsEstimation)
progressBarUpdatePeriod = max(1, inputArguments.nRhoValuesForSystematicsEstimation//1000)
fractionalUncertainty_rhoHigh = 0
deviationAtRhoHigh_rhoLow = 0
progressBar.initializeTimer()
for rhoCounter in range(0, inputArguments.nRhoValuesForSystematicsEstimation):
    if rhoCounter%progressBarUpdatePeriod == 0: progressBar.updateBar(1.0*rhoCounter/inputArguments.nRhoValuesForSystematicsEstimation, rhoCounter)
    rhoValue = inputArguments.rhoMinForSystematicsEstimation + (rhoCounter/(inputArguments.nRhoValuesForSystematicsEstimation-1))*(inputArguments.rhoMaxForSystematicsEstimation - inputArguments.rhoMinForSystematicsEstimation)
    rooKernel_PDF_Fits["rhoValues"][rhoCounter] = ROOT.RooKeysPdf("normBinKernelEstimateFunction_rhoCounter_{rhoC}".format(rhoC=rhoCounter), "normBinKernelEstimateFunction_rhoCounter_{rhoC}".format(rhoC=rhoCounter), rooVar_sT, sTRooDataSets[inputArguments.nJetsNorm], kernelOptionsObjects[inputArguments.kernelMirrorOption], rhoValue)
    integralObject_rhoValues_observationRange = rooKernel_PDF_Fits["rhoValues"][rhoCounter].createIntegral(ROOT.RooArgSet(rooVar_sT), "observation_sTRange")
    integralObject_rhoValues_normalizationRange = rooKernel_PDF_Fits["rhoValues"][rhoCounter].createIntegral(ROOT.RooArgSet(rooVar_sT), "normalization_sTRange")
    integralsRatio_rhoValues_normJets = integralObject_rhoValues_observationRange.getVal()/integralObject_rhoValues_normalizationRange.getVal()
    predicted_nEvents_rhoValues = integralsRatio_rhoValues_normJets*nEventsInNormWindows[inputArguments.nJetsNorm]
    predictedToObserved_atThisRho = predicted_nEvents_rhoValues/nEventsInObservationWindows[inputArguments.nJetsNorm]
    rhoSystematicsGraph.SetPoint(rhoSystematicsGraph.GetN(), rhoValue, predictedToObserved_atThisRho)
    if (rhoCounter == 0): fractionalUncertainty_rhoLow = abs((predictedToObserved_atThisRho/predictedToObservedRatio_nominalRho_observationRegion) - 1)
    if (rhoCounter == (-1 + inputArguments.nRhoValuesForSystematicsEstimation)): fractionalUncertainty_rhoHigh = abs((predictedToObserved_atThisRho/predictedToObservedRatio_nominalRho_observationRegion) - 1)
progressBar.terminate()
canvases["rhoValues"]["ratiosGraph"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [rhoSystematicsGraph], canvasName = "c_rhoSystematicsGraph", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_rhoSystematics".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix))

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
rooVar_sT.setRange(inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
sTFrames["rhoValues"]["dataAndKernels"] = rooVar_sT.frame(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
setFrameAesthetics(sTFrames["rhoValues"]["dataAndKernels"], "#it{S}_{T} (GeV)", "Events / ({binWidth} GeV)".format(binWidth=binWidth), "Kernel estimates, nJets = {nJets}".format(nJets=inputArguments.nJetsNorm))
sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["rhoValues"]["dataAndKernels"])
dataSet_legendEntry = dataAndKernelsLegend.AddEntry(sTRooDataSets[inputArguments.nJetsNorm], "Data")
rooKernel_PDF_Fits["rhoValues"][0].plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.LineColor(ROOT.kBlue), plotRange)
rhoMin_legendEntry = dataAndKernelsLegend.AddEntry(rooKernel_PDF_Fits["rhoValues"][0], "PDF fit: rho = {rhoMin:.1f}".format(rhoMin=inputArguments.rhoMinForSystematicsEstimation))
customizeLegendEntryForLine(rhoMin_legendEntry, ROOT.kBlue)
rooKernel_PDF_Fits["rhoValues"][-1+inputArguments.nRhoValuesForSystematicsEstimation].plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.LineColor(ROOT.kRed), plotRange)
rhoMax_legendEntry = dataAndKernelsLegend.AddEntry(rooKernel_PDF_Fits["rhoValues"][-1+inputArguments.nRhoValuesForSystematicsEstimation], "PDF fit: rho = {rhoMax:.1f}".format(rhoMax=inputArguments.rhoMaxForSystematicsEstimation))
customizeLegendEntryForLine(rhoMax_legendEntry, ROOT.kRed)
rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].plotOn(sTFrames["rhoValues"]["dataAndKernels"], ROOT.RooFit.LineColor(ROOT.kBlack), plotRange)
rho1_legendEntry = dataAndKernelsLegend.AddEntry(rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], "PDF fit: rho = 1.0")
customizeLegendEntryForLine(rho1_legendEntry, ROOT.kBlack)
canvases["rhoValues"]["dataAndKernels"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["rhoValues"]["dataAndKernels"], dataAndKernelsLegend], canvasName = "c_kernelPDFEstimate_dataAndKernels", outputROOTFile = outputFile, outputDocumentName = "{outputDirectory}/{outputPrefix}_kernelPDF_rhoValues".format(outputDirectory=inputArguments.outputDirectory_eventHistograms, outputPrefix=inputArguments.outputPrefix))

for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    sTTrees[nJetsBin].Write()

outputFile.Write()
outputFile.Close()

tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=dataSystematicsList, outputFilePath=("{outputDirectory}/{outputPrefix}_dataSystematics.dat".format(outputDirectory=inputArguments.outputDirectory_dataSystematics, outputPrefix=inputArguments.outputPrefix)))
tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=eventCountersList, outputFilePath=("{outputDirectory}/{outputPrefix}_eventCounters.dat".format(outputDirectory=inputArguments.outputDirectory_dataSystematics, outputPrefix=inputArguments.outputPrefix)))

sw.Stop()
print ('Real time: ' + str(sw.RealTime() / 60.0) + ' minutes')
print ('CPU time:  ' + str(sw.CpuTime() / 60.0) + ' minutes')
