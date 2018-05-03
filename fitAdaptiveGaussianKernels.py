#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, ROOT, argparse, array, pdb, math
import numpy as np
import tmROOTUtils, tmStatsUtils
from tmProgressBar import tmProgressBar

inputArgumentsParser = argparse.ArgumentParser(description='Run STEALTH selection.')
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--outputDirectoryPrefix', required=True, help='String to prefix to output directory name.',type=str)
inputArgumentsParser.add_argument('--sTPlotRangeMin', default=1000., help='Min value of sT to display in the plots.',type=float)
inputArgumentsParser.add_argument('--sTPlotRangeMax', default=2500., help='Max value of sT to display in the plots.',type=float)
inputArgumentsParser.add_argument('--n_sTBins', default=15, help='Number of sT bins (relevant for plotting only).',type=int)
inputArgumentsParser.add_argument('--sTKernelFitRangeMin', default=850., help='Min value of sT to use in the kernel fit. This should be slightly less than sTPlotRangeMin to try to get rid of boundary effects.',type=float)
inputArgumentsParser.add_argument('--sTKernelFitRangeMax', default=2500., help='Max value of sT to use in the kernel fit.',type=float)
inputArgumentsParser.add_argument('--sTNormRangeMin', default=1000., help='Min value of sT for normalization. For all sT distributions in nJets bins except the normalization bin, this value is the min of the range in which to scale the kernel estimator.',type=float)
inputArgumentsParser.add_argument('--sTNormRangeMax', default=1100., help='Max value of sT for normalization. For all sT distributions in nJets bins except the normalization bin, this value is the max of the range in which to scale the kernel estimator.',type=float)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Min number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsMax', default=4, help='Max number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsNorm', default=2, help='Number of jets w.r.t. which to normalize the sT distributions for other jets.',type=int)
inputArgumentsParser.add_argument('--nTargetEventsForSTThresholdOptimization', default=1., help='Number of jets w.r.t. which to normalize the sT distributions for other jets.',type=float)
inputArgumentsParser.add_argument('--nToyMCs', default=1000, help='Number of toy MC samples to generate using the pdf fits found.',type=int)
inputArgumentsParser.add_argument('--varyNEventsInPreNormWindowInToyMCs', action='store_true', help="Vary the number of generated events in toy MC samples in the pre-normalization window in a Poisson distribution about the number of events in this window in the original sample; default is to keep it fixed.")
inputArgumentsParser.add_argument('--varyNEventsInNormWindowInToyMCs', action='store_true', help="Vary the number of generated events in the toy MC samples in the normalization window in a Poisson distribution about the number of events in this window in the original sample; default is to keep it fixed.")
inputArgumentsParser.add_argument('--varyNEventsInObservationWindowInToyMCs', action='store_true', help="Vary the number of generated events in the toy MC samples in the observation window in a Poisson distribution about the number of events in this window in the original sample; default is to keep it fixed.")
inputArgumentsParser.add_argument('--rho', default=1., help='Value of parameter rho to be used in adaptive Gaussian kernel estimates.',type=float)
inputArgumentsParser.add_argument('--kernelMirrorOption', default="MirrorLeft", help='Kernel mirroring option to be used in adaptive Gaussian kernel estimates',type=str)
inputArguments = inputArgumentsParser.parse_args()
if (inputArguments.sTNormRangeMin < inputArguments.sTKernelFitRangeMin or inputArguments.sTNormRangeMax > inputArguments.sTKernelFitRangeMax):
    sys.exit("Normalization interval: ({nmin}, {nmax}) seems incompatible with kernel fitting range: ({smin, smax})".format(nmin=inputArguments.sTNormRangeMin, nmax=inputArguments.sTNormRangeMax, smin=inputArguments.sTKernelFitRangeMin, smax=inputArguments.sTKernelFitRangeMax))
kernelOptionsObjects = {"NoMirror": ROOT.RooKeysPdf.NoMirror, "MirrorLeft": ROOT.RooKeysPdf.MirrorLeft, "MirrorRight": ROOT.RooKeysPdf.MirrorRight, "MirrorBoth": ROOT.RooKeysPdf.MirrorBoth, "MirrorAsymLeft": ROOT.RooKeysPdf.MirrorAsymLeft, "MirrorAsymLeftRight": ROOT.RooKeysPdf.MirrorAsymLeftRight, "MirrorAsymRight": ROOT.RooKeysPdf.MirrorAsymRight, "MirrorLeftAsymRight": ROOT.RooKeysPdf.MirrorLeftAsymRight, "MirrorAsymBoth": ROOT.RooKeysPdf.MirrorAsymBoth}
if not(inputArguments.kernelMirrorOption in kernelOptionsObjects): sys.exit("The following element is passed as an argument for the kernel mirroring option but not in the dictionary defining the correspondence between kernel name and RooKeysPdf index: {kernelMirrorOption}".format(kernelMirrorOption=inputArguments.kernelMirrorOption))

rooVar_sT = ROOT.RooRealVar("rooVar_sT", "rooVar_sT", inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax, "GeV")
rooVar_sT.setRange("preNormalization_sTRange", inputArguments.sTKernelFitRangeMin, inputArguments.sTNormRangeMin)
rooVar_sT.setRange("normalization_sTRange", inputArguments.sTNormRangeMin, inputArguments.sTNormRangeMax)
rooVar_sT.setRange("observation_sTRange", inputArguments.sTNormRangeMax, inputArguments.sTKernelFitRangeMax)
rooVar_sT.setRange("kernelFit_sTRange", inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
plotRange = ROOT.RooFit.Range(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax)
kernelFitRange = ROOT.RooFit.Range(inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
normRange = ROOT.RooFit.Range(inputArguments.sTNormRangeMin, inputArguments.sTNormRangeMax)
binWidth = int(0.5 + (inputArguments.sTPlotRangeMax - inputArguments.sTPlotRangeMin)/inputArguments.n_sTBins)

plotRangeString = "plotRange_{plotMin:2.1f}_{plotMax:2.1f}_".format(plotMin=inputArguments.sTPlotRangeMin, plotMax=inputArguments.sTPlotRangeMax)
sTBinsString = "{n}sTBins_".format(n=inputArguments.n_sTBins)
kernelFitRangeString = "kernelFitRange_{kMin:2.1f}_{kMax:2.1f}_".format(kMin=inputArguments.sTKernelFitRangeMin, kMax=inputArguments.sTKernelFitRangeMax)
normRangeString = "normRange_{normMin:2.1f}_{normMax:2.1f}_".format(normMin=inputArguments.sTNormRangeMin, normMax=inputArguments.sTNormRangeMax)
nJetsMaxString = ""
if not(inputArguments.nJetsMax == 6): nJetsMaxString = "{n}JetsMax_".format(n=inputArguments.nJetsMax)
nJetsNormString = ""
if not(inputArguments.nJetsNorm == 3): nJetsNormString = "{n}JetsNorm_".format(n=inputArguments.nJetsNorm)
nToyMCsString = ""
if not(inputArguments.nToyMCs == 1000): nToyMCsString = "{n}ToyMCs_".format(n=inputArguments.nToyMCs)
preNormEventsString = ""
if (inputArguments.varyNEventsInPreNormWindowInToyMCs): preNormEventsString = "variableNEventsPreNorm_"
normEventsString = ""
if (inputArguments.varyNEventsInNormWindowInToyMCs): normEventsString = "variableNEventsNorm_"
obsEventsString = ""
if (inputArguments.varyNEventsInObservationWindowInToyMCs): obsEventsString = "variableNEventsObs_"
rhoString = ""
if not(inputArguments.rho == 1.): rhoString = "rho{rho:2.1f}_".format(rho=inputArguments.rho)
kernelMirrorOptionString = ""
if not(inputArguments.kernelMirrorOption == "MirrorLeft"): kernelMirrorOptionString = "option{opt}_".format(opt=inputArguments.kernelMirrorOption)
concatenatedString = ("{plotRangeString}{sTBinsString}{kernelFitRangeString}{normRangeString}{nJetsMaxString}{nJetsNormString}{nToyMCsString}{preNormEventsString}{normEventsString}{obsEventsString}{rhoString}{kernelMirrorOptionString}".format(**locals())).rstrip('_')
outputDirectoryName = ("{outputDirectoryPrefix}_{concatenatedString}".format(outputDirectoryPrefix=inputArguments.outputDirectoryPrefix, concatenatedString=concatenatedString)).replace('.', 'pt')
if (len(outputDirectoryName) > 255): sys.exit("Length of directory name should be no more than 255 characters long. Current name: {currentName}".format(currentName=outputDirectoryName))
outputDirectoryPath = "analysis/{outputDirectoryName}".format(outputDirectoryName=outputDirectoryName)

if not(os.path.isdir(outputDirectoryPath)): os.system("mkdir -p {outputDirectoryPath}".format(outputDirectoryPath=outputDirectoryPath))

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

outputFile = ROOT.TFile('{outputDirectoryPath}/savedObjects.root'.format(outputDirectoryPath=outputDirectoryPath), 'recreate')

# Initialize TTrees
sTTrees = {}
sTArrays = {}
for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    treeName = "sTTree_{nJets}Jets".format(nJets=nJets)
    sTTrees[nJets] = ROOT.TTree(treeName, treeName)
    sTArrays[nJets] = array.array('f', [0.])
    (sTTrees[nJets]).Branch('rooVar_sT', (sTArrays[nJets]), 'rooVar_sT/F')

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
for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    sTRooDataSetName = "rooDataSet_{nJets}Jets".format(nJets=nJets)
    sTRooDataSets[nJets] = ROOT.RooDataSet(sTRooDataSetName, sTRooDataSetName, sTTrees[nJets], ROOT.RooArgSet(rooVar_sT))
    nEventsInPreNormWindows[nJets] = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(sTRooDataSets[nJets], "preNormalization_sTRange")
    nEventsInNormWindows[nJets] = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(sTRooDataSets[nJets], "normalization_sTRange")
    nEventsInObservationWindows[nJets] = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(sTRooDataSets[nJets], "observation_sTRange")
    total_nEventsInFullRange[nJets] = nEventsInPreNormWindows[nJets] + nEventsInNormWindows[nJets] + nEventsInObservationWindows[nJets]
    print("At nJets = {nJets}, nEventsInPreNormWindow = {preNorm}, nEventsInNormWindow = {norm}, nEventsInObservationWindow = {obs}".format(nJets = nJets, preNorm = nEventsInPreNormWindows[nJets], norm = nEventsInNormWindows[nJets], obs = nEventsInObservationWindows[nJets]))

rooKernel_PDF_Fits = {
    "data": {},
    "toyMC": {}
}
canvases = {
    "data": {},
    "toyMC": {}
}
sTFrames = {
    "data": {},
    "toyMC": {}
}

# Find fits for norm bin
rooVar_sT.setRange(inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm] = ROOT.RooKeysPdf("normBinKernelEstimateFunction", "normBinKernelEstimateFunction", rooVar_sT, sTRooDataSets[inputArguments.nJetsNorm], kernelOptionsObjects[inputArguments.kernelMirrorOption], inputArguments.rho)
integralObject_normalizationRange = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].createIntegral(ROOT.RooArgSet(rooVar_sT), "normalization_sTRange")
normFactor = nEventsInNormWindows[inputArguments.nJetsNorm]/integralObject_normalizationRange.getVal()
print("Check 1 on norm factor: expected events in norm range: {nExpected}, observed: {nObserved}".format(nExpected=getExpectedNEventsFromPDFInNamedRange(normFactor=normFactor, inputRooPDF=rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], inputRooArgSet=ROOT.RooArgSet(rooVar_sT), targetRangeMin=inputArguments.sTNormRangeMin, targetRangeMax=inputArguments.sTNormRangeMax), nObserved=nEventsInNormWindows[inputArguments.nJetsNorm]))
print("Check 2 on norm factor: expected events in observation range: {nExpected}, observed: {nObserved}".format(nExpected=getExpectedNEventsFromPDFInNamedRange(normFactor=normFactor, inputRooPDF=rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], inputRooArgSet=ROOT.RooArgSet(rooVar_sT), targetRangeMin=inputArguments.sTNormRangeMax, targetRangeMax=inputArguments.sTKernelFitRangeMax), nObserved=nEventsInObservationWindows[inputArguments.nJetsNorm]))
optimal_sTThreshold = tmStatsUtils.getMonotonicFunctionApproximateZero(inputFunction=(lambda sTThreshold: (-1*inputArguments.nTargetEventsForSTThresholdOptimization)+getExpectedNEventsFromPDFInNamedRange(normFactor=normFactor, inputRooPDF=rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], inputRooArgSet=ROOT.RooArgSet(rooVar_sT), targetRangeMin=sTThreshold, targetRangeMax=inputArguments.sTKernelFitRangeMax)), xRange=[1.1*inputArguments.sTNormRangeMax, 0.99*inputArguments.sTKernelFitRangeMax], autoZeroTolerance=True, printDebug=True)
optimalThresholdOutputFile = open("{outputDirectoryPath}/sTOptimalThreshold.dat".format(outputDirectoryPath=outputDirectoryPath), 'w')
optimalThresholdOutputFile.write("Optimal sT Threshold: {thr}\n".format(thr=optimal_sTThreshold))
optimalThresholdOutputFile.close()
rooVar_sT.setRange(inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
sTFrames["data"][inputArguments.nJetsNorm] = rooVar_sT.frame(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm])
rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], plotRange)
setFrameAesthetics(sTFrames["data"][inputArguments.nJetsNorm], "#it{S}_{T} (GeV)", "Events / ({binWidth} GeV)".format(binWidth=binWidth), "Normalization bin: {nJets} Jets".format(nJets=inputArguments.nJetsNorm))
canvases["data"][inputArguments.nJetsNorm] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["data"][inputArguments.nJetsNorm]], canvasName = "c_kernelPDF_normJetsBin", outputROOTFile = outputFile, outputDocumentName = "{outputDirectoryPath}/kernelPDF_normJetsBin".format(outputDirectoryPath=outputDirectoryPath))

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
    if (not(inputArguments.varyNEventsInNormWindowInToyMCs) and not(nToyEventsInNormWindow == nEventsInNormWindows[inputArguments.nJetsNorm])): sys.exit("Error: check nGeneratedEvents") # leaving in for now
    if (not(inputArguments.varyNEventsInObservationWindowInToyMCs) and not(nToyEventsInObservationWindow == nEventsInObservationWindows[inputArguments.nJetsNorm])): sys.exit("Error: check nGeneratedEvents") # leaving in for now
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
canvases["toyMC"]["DataAndFits"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["toyMC"]["DataAndFits"]], canvasName = "c_toyMCDataAndKernelEstimates", outputROOTFile = outputFile, outputDocumentName = "{outputDirectoryPath}/toyMCDataAndKernelEstimates".format(outputDirectoryPath=outputDirectoryPath))

# Plot the shape systematics estimate
canvases["toyMC"]["systematics"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [predictedToObservedNEventsHistogram], canvasName = "c_shapeSystematics", outputROOTFile = outputFile, outputDocumentName = "{outputDirectoryPath}/shapeSystematics".format(outputDirectoryPath=outputDirectoryPath))

# Plot the integral checks, to see that the normalization is similar
canvases["toyMC"]["systematicsCheck"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [totalIntegralCheckHistogram], canvasName = "c_shapeSystematicsCheck", outputROOTFile = outputFile, outputDocumentName = "{outputDirectoryPath}/shapeSystematicsCheck".format(outputDirectoryPath=outputDirectoryPath))

# Finally use these fits in other nJets bins and obtain estimate of systematic on assumption that sT scales
rooVar_nEventsInNormBin = {}
rooKernel_extendedPDF_Fits = {}
scalingSystematicsOutputFile = open("{outputDirectoryPath}/sTScalingSystematics.dat".format(outputDirectoryPath=outputDirectoryPath), 'w')
integralObject_normJets_observationRange = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].createIntegral(ROOT.RooArgSet(rooVar_sT), "observation_sTRange")
integralObject_normJets_normRange = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].createIntegral(ROOT.RooArgSet(rooVar_sT), "normalization_sTRange")
integralsRatio_normJets = 1.0*integralObject_normJets_observationRange.getVal() / integralObject_normJets_normRange.getVal()
scalingSystematicsOutputFile.write("{nJetsBin}    {fraction}\n".format(nJetsBin=inputArguments.nJetsNorm, fraction=integralsRatio_normJets))
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    if (nJetsBin == inputArguments.nJetsNorm): continue
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
    rooKernel_extendedPDF_Fits[nJetsBin].plotOn(sTFrames["data"][nJetsBin], ROOT.RooFit.LineColor(ROOT.kBlue), plotRange)
    # nEventsInNormWindow = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(sTRooDataSets[nJetsBin], "normalization_sTRange")
    # predicted_nEventsInObservationWindow = integralsRatio_normJets*nEventsInNormWindow
    # observed_nEventsInObservationWindow = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(sTRooDataSets[nJetsBin], "observation_sTRange")
    # fraction_predictedToActual = 1.0*predicted_nEventsInObservationWindow / observed_nEventsInObservationWindow
    # fraction_predictedToActual = integralsRatio/integralsRatio_normJets
    fraction_predictedToActual = integralsRatio
    scalingSystematicsOutputFile.write("{nJetsBin}    {fraction}\n".format(nJetsBin=nJetsBin, fraction=fraction_predictedToActual))
    if (nJetsBin == inputArguments.nJetsMax): setFrameAesthetics(sTFrames["data"][nJetsBin], "#it{S}_{T} (GeV)", "Events / ({binWidth} GeV)".format(binWidth=binWidth), "#geq {nJetsBin} Jets".format(nJetsBin=nJetsBin))
    else: setFrameAesthetics(sTFrames["data"][nJetsBin], "#it{S}_{T} (GeV)", "Events / ({binWidth} GeV)".format(binWidth=binWidth), "{nJetsBin} Jets".format(nJetsBin=nJetsBin))
    canvases["data"][nJetsBin] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["data"][nJetsBin]], canvasName = "c_kernelPDF_{nJetsBin}Jets".format(nJetsBin=nJetsBin), outputROOTFile = outputFile, outputDocumentName = "{outputDirectoryPath}/kernelPDF_{nJetsBin}Jets".format(outputDirectoryPath=outputDirectoryPath, nJetsBin=nJetsBin))
scalingSystematicsOutputFile.close()

for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    sTTrees[nJetsBin].Write()

outputFile.Write()
outputFile.Close()

sw.Stop()
print ('Real time: ' + str(sw.RealTime() / 60.0) + ' minutes')
print ('CPU time:  ' + str(sw.CpuTime() / 60.0) + ' minutes')
