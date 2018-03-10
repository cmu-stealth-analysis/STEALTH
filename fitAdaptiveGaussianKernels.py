#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, ROOT, argparse, array, pdb, math
import numpy as np
import tmROOTUtils
from tmProgressBar import tmProgressBar

inputArgumentsParser = argparse.ArgumentParser(description='Run STEALTH selection.')
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--outputDirectoryPrefix', required=True, help='String to prefix to output directory name.',type=str)
inputArgumentsParser.add_argument('--sTPlotRangeMin', default=800., help='Min value of sT to display in the plots.',type=float)
inputArgumentsParser.add_argument('--sTPlotRangeMax', default=2500., help='Max value of sT to display in the plots.',type=float)
inputArgumentsParser.add_argument('--n_sTBins', default=15, help='Number of sT bins (relevant for plotting only).',type=int)
inputArgumentsParser.add_argument('--sTKernelFitRangeMin', default=750., help='Min value of sT to use in the kernel fit. This should be slightly less than sTPlotRangeMin to try to get rid of boundary effects.',type=float)
inputArgumentsParser.add_argument('--sTKernelFitRangeMax', default=2500., help='Max value of sT to use in the kernel fit.',type=float)
inputArgumentsParser.add_argument('--sTNormRangeMin', default=800., help='Min value of sT for normalization. For all sT distributions in nJets bins except the normalization bin, this value is the min of the range in which to scale the kernel estimator.',type=float)
inputArgumentsParser.add_argument('--sTNormRangeMax', default=900., help='Max value of sT for normalization. For all sT distributions in nJets bins except the normalization bin, this value is the max of the range in which to scale the kernel estimator.',type=float)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Min number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Max number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsNorm', default=3, help='Number of jets w.r.t. which to normalize the sT distributions for other jets.',type=int)
inputArgumentsParser.add_argument('--nToyMCs', default=1000, help='Number of toy MC samples to generate using the pdf fits found.',type=int)
inputArgumentsParser.add_argument('--fixTotalNEventsInToyMCs', action='store_true', help="Keep the number of generated events in the toy MC samples equal to the number of events in the original sample. If this argument is not passed, the default behavior is to vary the number of events generated in the toy MCs following a Poisson distribution about the expected mean.")
inputArgumentsParser.add_argument('--fixNEventsInNormWindowInToyMCs', action='store_true', help="Keep the number of generated events in the toy MC samples in the normalization window equal to the number of events in the original sample in the normalization window. If this argument is not passed, the default behavior is to accept all generated MC samples regardless of the number of events in the normalization window; in that case the systematics check plots the ratio of integrals times the number of events in the normalization window. Use option \"--forceNoScaling\" to turn this additional scaling off.")
inputArgumentsParser.add_argument('--forceNoScaling', action='store_true', help="Force no scaling of the systematics by the number of events in the normalization window; see option \"--fixNEventsInNormWindowInToyMCs\".")
inputArgumentsParser.add_argument('--rho', default=1., help='Value of parameter rho to be used in adaptive Gaussian kernel estimates.',type=float)
inputArgumentsParser.add_argument('--kernelMirrorOption', default="MirrorLeftAsymRight", help='Kernel mirroring option to be used in adaptive Gaussian kernel estimates',type=str)
inputArguments = inputArgumentsParser.parse_args()
if (inputArguments.sTNormRangeMin < inputArguments.sTKernelFitRangeMin or inputArguments.sTNormRangeMax > inputArguments.sTKernelFitRangeMax):
    sys.exit("Normalization interval: ({nmin}, {nmax}) seems incompatible with kernel fitting range: ({smin, smax})".format(nmin=inputArguments.sTNormRangeMin, nmax=inputArguments.sTNormRangeMax, smin=inputArguments.sTKernelFitRangeMin, smax=inputArguments.sTKernelFitRangeMax))
if (inputArguments.forceNoScaling and inputArguments.fixNEventsInNormWindowInToyMCs): sys.exit("Option \"--fixNEventsInNormWindowInToyMCs\" implies that the number of events in the normalization window is to be fixed, but option \"forceNoScaling\" implies otherwise. Please check arguments.")
kernelOptionsObjects = {"NoMirror": ROOT.RooKeysPdf.NoMirror, "MirrorLeft": ROOT.RooKeysPdf.MirrorLeft, "MirrorRight": ROOT.RooKeysPdf.MirrorRight, "MirrorBoth": ROOT.RooKeysPdf.MirrorBoth, "MirrorAsymLeft": ROOT.RooKeysPdf.MirrorAsymLeft, "MirrorAsymLeftRight": ROOT.RooKeysPdf.MirrorAsymLeftRight, "MirrorAsymRight": ROOT.RooKeysPdf.MirrorAsymRight, "MirrorLeftAsymRight": ROOT.RooKeysPdf.MirrorLeftAsymRight, "MirrorAsymBoth": ROOT.RooKeysPdf.MirrorAsymBoth}
if not(inputArguments.kernelMirrorOption in kernelOptionsObjects): sys.exit("The following element is passed as an argument for the kernel mirroring option but not in the dictionary defining the correspondence between kernel name and RooKeysPdf index: {kernelMirrorOption}".format(kernelMirrorOption=inputArguments.kernelMirrorOption))

rooVar_sT = ROOT.RooRealVar("rooVar_sT", "rooVar_sT", inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax, "GeV")
rooVar_sT.setRange("normalization_sTRange", inputArguments.sTNormRangeMin, inputArguments.sTNormRangeMax)
rooVar_sT.setRange("observation_sTRange", inputArguments.sTNormRangeMax, inputArguments.sTKernelFitRangeMax)
rooVar_sT.setRange("kernelFit_sTRange", inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
plotRange = ROOT.RooFit.Range(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax)
kernelFitRange = ROOT.RooFit.Range(inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
normRange = ROOT.RooFit.Range(inputArguments.sTNormRangeMin, inputArguments.sTNormRangeMax)

totalEventsString = "variableTotalNEventsInToyMCs"
if (inputArguments.fixTotalNEventsInToyMCs): totalEventsString = "fixedTotalNEventsInToyMCs"
normEventsString = "variableNEventsInNormWindowInToyMCs"
if (inputArguments.fixNEventsInNormWindowInToyMCs): totalEventsString = "fixedNEventsInNormWindowInToyMCs"
outputDirectoryName = ("{outputDirectoryPrefix}_normMin_{normMin:2.1f}_normMax_{normMax:2.1f}_{nJetsMax}JetsMax_{nJetsNorm}JetsNorm_{nToyMCs}ToyMCs_{totalEventsString}_{normEventsString}_rho_{rho:2.1f}_{kernelMirrorOption}".format(outputDirectoryPrefix=inputArguments.outputDirectoryPrefix, normMin=inputArguments.sTNormRangeMin, normMax=inputArguments.sTNormRangeMax, nJetsMax=inputArguments.nJetsMax, nJetsNorm=inputArguments.nJetsNorm, nToyMCs=inputArguments.nToyMCs, totalEventsString=totalEventsString, normEventsString=normEventsString, rho=inputArguments.rho, kernelMirrorOption=inputArguments.kernelMirrorOption)).replace('.', 'pt')
if (len(outputDirectoryName) > 255): sys.exit("Length of directory name should be no more than 255 characters long. Current name: {currentName}".format(currentName=outputDirectoryName))
outputDirectoryPath = "analysis/{outputDirectoryName}".format(outputDirectoryName=outputDirectoryName)

if not(os.path.isdir(outputDirectoryPath)): os.system("mkdir -p {outputDirectoryPath}".format(outputDirectoryPath=outputDirectoryPath))

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
nEventsInNormWindows = {}
for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    sTRooDataSetName = "rooDataSet_{nJets}Jets".format(nJets=nJets)
    sTRooDataSets[nJets] = ROOT.RooDataSet(sTRooDataSetName, sTRooDataSetName, sTTrees[nJets], ROOT.RooArgSet(rooVar_sT))
    nEventsInNormWindows[nJets] = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(sTRooDataSets[nJets], "normalization_sTRange")
    print("At nJets = {nJets}, nEventsInNormWindow = {nEventsInWindow}".format(nJets = nJets, nEventsInWindow = nEventsInNormWindows[nJets]))

nEventsInNormJetsBin = sTTrees[inputArguments.nJetsNorm].GetEntries()
print("nEventsInNormJetsBin = {nEventsInNormJetsBin}".format(nEventsInNormJetsBin = nEventsInNormJetsBin))

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
rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm] = ROOT.RooKeysPdf("normBinKernelEstimateFunction", "normBinKernelEstimateFunction", rooVar_sT, sTRooDataSets[inputArguments.nJetsNorm], kernelOptionsObjects[inputArguments.kernelMirrorOption], inputArguments.rho)
sTFrames["data"][inputArguments.nJetsNorm] = rooVar_sT.frame(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm])
rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][inputArguments.nJetsNorm], plotRange)
canvases["data"][inputArguments.nJetsNorm] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["data"][inputArguments.nJetsNorm]], canvasName = "c_kernelPDF_normJetsBin", outputROOTFile = outputFile, outputDocumentName = "{outputDirectoryPath}/kernelPDF_normJetsBin".format(outputDirectoryPath=outputDirectoryPath))

# Generate and fit toy MC datsets with the fit kernels
toyRooDataSets = {}
sTFrames["toyMC"]["DataAndFits"] = rooVar_sT.frame(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
kernelIntegralRatiosHistogram = ROOT.TH1F("h_kernelIntegralRatios", "h_kernelIntegralRatios", 40, 0., 0.)
totalIntegralCheckHistogram = ROOT.TH1F("h_totalIntegralCheck", "h_totalIntegralCheck", 40, 0., 0.)
nToyMCEventsHistogram = ROOT.TH1F("h_nToyMCEvents", "h_nToyMCEvents", 1 + int(0.5 + nEventsInNormJetsBin + 5*math.sqrt(nEventsInNormJetsBin)) - int(0.5 + nEventsInNormJetsBin - 5*math.sqrt(nEventsInNormJetsBin)), -0.5 + int(0.5 + nEventsInNormJetsBin - 5*math.sqrt(nEventsInNormJetsBin)), 0.5 + int(0.5 + nEventsInNormJetsBin + 5*math.sqrt(nEventsInNormJetsBin))) # Expected mean: nEventsInNormJetsBin, expected sigma = sqrt(nEventsInNormJetsBin)
goodMCSampleIndex = 0
randomGenerator = ROOT.TRandom1()
randomGenerator.SetSeed(0) # Sets seed by using some information from a ROOT "UUID"
progressBar = tmProgressBar(inputArguments.nToyMCs)
progressBarUpdatePeriod = max(1, inputArguments.nToyMCs//1000)
progressBar.initializeTimer()
while goodMCSampleIndex < inputArguments.nToyMCs:
    nEventsToGenerate = randomGenerator.Poisson(nEventsInNormJetsBin)
    if (inputArguments.fixTotalNEventsInToyMCs): nEventsToGenerate = nEventsInNormJetsBin
    toyRooDataSets[goodMCSampleIndex] = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].generate(ROOT.RooArgSet(rooVar_sT), nEventsToGenerate)
    nToyEventsInNormWindow = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(toyRooDataSets[goodMCSampleIndex], "normalization_sTRange")
    if (not(nToyEventsInNormWindow == nEventsInNormWindows[inputArguments.nJetsNorm]) and inputArguments.fixNEventsInNormWindowInToyMCs): continue
    nToyMCEventsHistogram.Fill(1.0*nEventsToGenerate)
    toyRooDataSets[goodMCSampleIndex].plotOn(sTFrames["toyMC"]["DataAndFits"])
    rooKernel_PDF_Fits["toyMC"][goodMCSampleIndex] = ROOT.RooKeysPdf("toyMCKernelEstimateFunction_{index}".format(index=goodMCSampleIndex), "toyMCKernelEstimateFunction_{index}".format(index=goodMCSampleIndex), rooVar_sT, toyRooDataSets[goodMCSampleIndex], kernelOptionsObjects[inputArguments.kernelMirrorOption], inputArguments.rho)
    rooKernel_PDF_Fits["toyMC"][goodMCSampleIndex].plotOn(sTFrames["toyMC"]["DataAndFits"])
    integralObject_observationRange = rooKernel_PDF_Fits["toyMC"][goodMCSampleIndex].createIntegral(ROOT.RooArgSet(rooVar_sT), "observation_sTRange")
    integralObject_normRange = rooKernel_PDF_Fits["toyMC"][goodMCSampleIndex].createIntegral(ROOT.RooArgSet(rooVar_sT), "normalization_sTRange")
    integralsRatio = integralObject_observationRange.getVal() / integralObject_normRange.getVal()
    if (not(inputArguments.fixNEventsInNormWindowInToyMCs) and not(inputArguments.forceNoScaling)):
        integralsRatio = integralsRatio*nToyEventsInNormWindow
    kernelIntegralRatiosHistogram.Fill(integralsRatio)
    totalIntegralCheckObject = rooKernel_PDF_Fits["toyMC"][goodMCSampleIndex].createIntegral(ROOT.RooArgSet(rooVar_sT), "kernelFit_sTRange")
    totalIntegralCheckHistogram.Fill(totalIntegralCheckObject.getVal())
    if goodMCSampleIndex%progressBarUpdatePeriod == 0: progressBar.updateBar(1.0*goodMCSampleIndex/inputArguments.nToyMCs, goodMCSampleIndex)
    goodMCSampleIndex += 1
progressBar.terminate()
sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames["toyMC"]["DataAndFits"], ROOT.RooFit.LineColor(ROOT.kRed))
rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].plotOn(sTFrames["toyMC"]["DataAndFits"], ROOT.RooFit.LineColor(ROOT.kRed), plotRange)

# Plot the toy MC data and fits
canvases["toyMC"]["DataAndFits"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["toyMC"]["DataAndFits"]], canvasName = "c_toyMCDataAndKernelEstimates", outputROOTFile = outputFile, outputDocumentName = "{outputDirectoryPath}/toyMCDataAndKernelEstimates".format(outputDirectoryPath=outputDirectoryPath))

# Plot the shape systematics estimate
canvases["toyMC"]["systematics"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [kernelIntegralRatiosHistogram], canvasName = "c_shapeSystematics", outputROOTFile = outputFile, outputDocumentName = "{outputDirectoryPath}/shapeSystematics".format(outputDirectoryPath=outputDirectoryPath))

# Plot the integral checks, to see that the normalization is similar
canvases["toyMC"]["systematicsCheck"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [totalIntegralCheckHistogram], canvasName = "c_shapeSystematicsCheck", outputROOTFile = outputFile, outputDocumentName = "{outputDirectoryPath}/shapeSystematicsCheck".format(outputDirectoryPath=outputDirectoryPath))

# Plot a histogram of the generated number of events
canvases["toyMC"]["nToyMCEvents"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [nToyMCEventsHistogram], canvasName = "c_nToyMCEvents", outputROOTFile = outputFile, outputDocumentName = "{outputDirectoryPath}/nToyMCEvents".format(outputDirectoryPath=outputDirectoryPath))
        
# Finally use these fits in other nJets bins and obtain estimate of systematic on assumption that sT scales
rooVar_nEventsInNormBin = {}
rooKernel_extendedPDF_Fits = {}
scalingSystematicsOutputFile = open("{outputDirectoryPath}/sTScalingSystematics.dat".format(outputDirectoryPath=outputDirectoryPath), 'w')
for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    if (nJets == inputArguments.nJetsNorm): continue
    sTFrames["data"][nJets] = rooVar_sT.frame(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
    sTRooDataSets[nJets].plotOn(sTFrames["data"][nJets])
    rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].plotOn(sTFrames["data"][nJets], ROOT.RooFit.LineColor(ROOT.kRed), plotRange)
    rooVar_nEventsInNormBin[nJets] = ROOT.RooRealVar("rooVar_nEventsInNormBin_{nJets}Jets".format(nJets=nJets), "rooVar_nEventsInNormBin_{nJets}Jets".format(nJets=nJets), 100, 0, 10000)
    rooKernel_extendedPDF_Fits[nJets] = ROOT.RooExtendPdf("extendedKernelPDF_{nJets}Jets".format(nJets=nJets), "extendedKernelPDF_{nJets}Jets".format(nJets=nJets), rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm], rooVar_nEventsInNormBin[nJets], "kernelFit_sTRange")
    rooKernel_extendedPDF_Fits[nJets].fitTo(sTRooDataSets[nJets], normRange, ROOT.RooFit.Minos(ROOT.kTRUE), ROOT.RooFit.PrintLevel(0))
    rooVar_nEventsInNormBin[nJets].Print()
    rooKernel_extendedPDF_Fits[nJets].plotOn(sTFrames["data"][nJets], ROOT.RooFit.LineColor(ROOT.kBlue), plotRange)
    integralObject_observationRange = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].createIntegral(ROOT.RooArgSet(rooVar_sT), "observation_sTRange")
    integralObject_normRange = rooKernel_PDF_Fits["data"][inputArguments.nJetsNorm].createIntegral(ROOT.RooArgSet(rooVar_sT), "normalization_sTRange")
    predicted_nEventsInObservationWindow = integralObject_observationRange.getVal() / integralObject_normRange.getVal()
    nEventsInNormWindow = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(sTRooDataSets[nJets], "normalization_sTRange")
    if (not(inputArguments.forceNoScaling)):
        predicted_nEventsInObservationWindow = predicted_nEventsInObservationWindow*nEventsInNormWindow
    nEventsInObservationWindow = tmROOTUtils.getNEventsInNamedRangeInRooDataSet(sTRooDataSets[nJets], "observation_sTRange")
    fraction_predictedToActual = predicted_nEventsInObservationWindow / nEventsInObservationWindow
    scalingSystematicsOutputFile.write("{nJets}    {fraction}\n".format(nJets=nJets, fraction=fraction_predictedToActual))
    canvases["data"][nJets] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [sTFrames["data"][nJets]], canvasName = "c_kernelPDF_{nJets}Jets".format(nJets=nJets), outputROOTFile = outputFile, outputDocumentName = "{outputDirectoryPath}/kernelPDF_{nJets}Jets".format(outputDirectoryPath=outputDirectoryPath, nJets=nJets))
scalingSystematicsOutputFile.close()

for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    sTTrees[nJets].Write()

outputFile.Write()
outputFile.Close()

sw.Stop()
print ('Real time: ' + str(sw.RealTime() / 60.0) + ' minutes')
print ('CPU time:  ' + str(sw.CpuTime() / 60.0) + ' minutes')
