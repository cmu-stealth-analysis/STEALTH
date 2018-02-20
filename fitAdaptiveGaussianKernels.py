#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, ROOT, argparse, array, pdb
import numpy as np
import tmROOTUtils
from tmProgressBar import tmProgressBar

inputArgumentsParser = argparse.ArgumentParser(description='Run STEALTH selection.')
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
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
inputArgumentsParser.add_argument('--outputFilesString', required=True, help='String to include in all output file names.',type=str)
inputArgumentsParser.add_argument('--enableRho', action='append', help='Value of the adaptive Gaussian fit parameter rho to be enabled; repeat argument multiple times for multiple values.', type=float)
inputArgumentsParser.add_argument('--enableKernel', action='append', help='Type of kernel used by adaptive Gaussian fit to be enabled; repeat argument multiple times for multiple values.', type=str)
inputArguments = inputArgumentsParser.parse_args()
if (inputArguments.sTNormRangeMin < inputArguments.sTKernelFitRangeMin or inputArguments.sTNormRangeMax > inputArguments.sTKernelFitRangeMax):
    print ("Normalization interval: ({nmin}, {nmax}) seems incompatible with kernel fitting range: ({smin, smax})".format(nmin=inputArguments.sTNormRangeMin, nmax=inputArguments.sTNormRangeMax, smin=inputArguments.sTKernelFitRangeMin, smax=inputArguments.sTKernelFitRangeMax))
enabledRhos = inputArguments.enableRho
enabledKernels = inputArguments.enableKernel
if (inputArguments.enableRho is None): enabledRhos = [1.0]
if (inputArguments.enableKernel is None): enabledKernels = ["MirrorLeftAsymRight"]
kernels = {"NoMirror": ROOT.RooKeysPdf.NoMirror, "MirrorLeft": ROOT.RooKeysPdf.MirrorLeft, "MirrorRight": ROOT.RooKeysPdf.MirrorRight, "MirrorBoth": ROOT.RooKeysPdf.MirrorBoth, "MirrorAsymLeft": ROOT.RooKeysPdf.MirrorAsymLeft, "MirrorAsymLeftRight": ROOT.RooKeysPdf.MirrorAsymLeftRight, "MirrorAsymRight": ROOT.RooKeysPdf.MirrorAsymRight, "MirrorLeftAsymRight": ROOT.RooKeysPdf.MirrorLeftAsymRight, "MirrorAsymBoth": ROOT.RooKeysPdf.MirrorAsymBoth}
for enabledKernel in enabledKernels:
    if not(enabledKernel in kernels.keys()): sys.exit("The following element is present in list of enabled kernels but not in the dictionary defining the correspondence between kernel name and RooKeysPdf index: {enabledKernel}".format(enabledKernel=enabledKernel))

rooVar_sT = ROOT.RooRealVar("rooVar_sT", "rooVar_sT", inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax, "GeV")
rooVar_sT.setRange("restricted_sTRange", inputArguments.sTNormRangeMin, inputArguments.sTNormRangeMax)
rooVar_sT.setRange("full_sTRange", inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
rooVar_sT.setRange("observation_sTRange", inputArguments.sTNormRangeMax, inputArguments.sTKernelFitRangeMax)
rooVar_nEventsInNormRegion = ROOT.RooRealVar("rooVar_nEventsInNormRegion", "rooVar_nEventsInNormRegion", 20, 0, 10000.)
plotRange = ROOT.RooFit.Range(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax)
kernelFitRange = ROOT.RooFit.Range(inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
normRange = ROOT.RooFit.Range(inputArguments.sTNormRangeMin, inputArguments.sTNormRangeMax)

sw = ROOT.TStopwatch()
sw.Start()

inputChain = ROOT.TChain('ggNtuplizer/EventTree')
inputChain.Add(inputArguments.inputFilePath)
nEntries = inputChain.GetEntries()
print ("Total number of available events: {nEntries}".format(nEntries=nEntries))

outputFile = ROOT.TFile('analysis/sTDistributions_{outputFilesString}_{nJetsNorm}JetsNorm_{nJetsMax}JetsMax_{n_sTBins}Bins.root'.format(outputFilesString=inputArguments.outputFilesString, nJetsNorm=inputArguments.nJetsNorm, nJetsMax=inputArguments.nJetsMax, n_sTBins=inputArguments.n_sTBins), 'recreate')

# Initialize TTrees
sTTrees = {}
sTArrays = {}
for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    treeName = "sTTree_{nJets}Jets".format(nJets=nJets)
    sTTrees[nJets] = ROOT.TTree(treeName, treeName)
    sTArrays[nJets] = array.array('f', [0.])
    (sTTrees[nJets]).Branch('rooVar_sT', (sTArrays[nJets]), 'rooVar_sT/F')

# Fill trees
progressBar = tmProgressBar(nEntries)
progressBarUpdatePeriod = nEntries//1000
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
for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    sTRooDataSetName = "rooDataSet_{nJets}Jets".format(nJets=nJets)
    sTRooDataSets[nJets] = ROOT.RooDataSet(sTRooDataSetName, sTRooDataSetName, sTTrees[nJets], ROOT.RooArgSet(rooVar_sT))

nEventsInNormJetsBin = sTTrees[inputArguments.nJetsNorm].GetEntries()
# First find the kernel fits in background nJets bin
rooKernel_PDF_Fits = {}
canvases = {}
sTFrames = {}
toyFits = {}
toy_sTFrames = {}
for kernelType in enabledKernels:
    print("Performing fits for kernel type: {kernelType}".format(kernelType=kernelType))
    rooKernel_PDF_Fits[kernelType] = {}
    canvases[kernelType] = {}
    sTFrames[kernelType] = {}
    toyFits[kernelType] = {}
    toy_sTFrames[kernelType] = {}
    for rho in enabledRhos:
        print("Performing fits for rho: {rho}".format(rho=rho))
        rhoStr = ("rho_{rho:2.1f}".format(rho=rho)).replace('.', 'pt')
        canvases[kernelType][rhoStr] = {}
        sTFrames[kernelType][rhoStr] = {}
        # First find fits for norm bin
        functionName = "rooKernelFunction_{kernelType}_{rhoStr}".format(kernelType=kernelType, rhoStr=rhoStr)
        rooKernel_PDF_Fits[kernelType][rhoStr] = ROOT.RooKeysPdf(functionName, functionName, rooVar_sT, sTRooDataSets[inputArguments.nJetsNorm], kernels[kernelType], rho)
        sTFrames[kernelType][rhoStr][inputArguments.nJetsNorm] = rooVar_sT.frame(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
        sTRooDataSets[inputArguments.nJetsNorm].plotOn(sTFrames[kernelType][rhoStr][inputArguments.nJetsNorm], plotRange)
        rooKernel_PDF_Fits[kernelType][rhoStr].plotOn(sTFrames[kernelType][rhoStr][inputArguments.nJetsNorm], plotRange)
        canvasName = "c_sTUnbinnedFit_{kernelType}_{rhoStr}_norm".format(kernelType=kernelType, rhoStr=rhoStr)
        outputFileName = "analysis/plot_sT_UnbinnedFit_{outputFilesString}_{nJetsNorm}JetsNorm_{nJetsMax}JetsMax_{n_sTBins}Bins_{kernelType}_{rhoStr}_norm".format(outputFilesString=inputArguments.outputFilesString, nJetsNorm=inputArguments.nJetsNorm, nJetsMax=inputArguments.nJetsMax, n_sTBins=inputArguments.n_sTBins, kernelType=kernelType, rhoStr=rhoStr)
        plotList = [sTFrames[kernelType][rhoStr][inputArguments.nJetsNorm]]
        canvases[kernelType][rhoStr][inputArguments.nJetsNorm] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = plotList, canvasName = canvasName, outputROOTFile = outputFile, outputDocumentName = outputFileName)
        # Next use these fits in other nJets bins
        for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
            if (nJets == inputArguments.nJetsNorm): continue
            sTFrames[kernelType][rhoStr][nJets] = rooVar_sT.frame(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
            sTRooDataSets[nJets].plotOn(sTFrames[kernelType][rhoStr][nJets], plotRange)
            rooKernel_PDF_Fits[kernelType][rhoStr].plotOn(sTFrames[kernelType][rhoStr][nJets], ROOT.RooFit.LineColor(ROOT.kRed), plotRange)
            extendedPDFName = "extendedPDF_{kernelType}_{rhoStr}_{nJets}Jets".format(kernelType=kernelType, rhoStr=rhoStr, nJets=nJets)
            extendedPDF = ROOT.RooExtendPdf(extendedPDFName, extendedPDFName, rooKernel_PDF_Fits[kernelType][rhoStr], rooVar_nEventsInNormRegion, "restricted_sTRange")
            extendedPDF.fitTo(sTRooDataSets[nJets], normRange, ROOT.RooFit.Minos(ROOT.kTRUE), ROOT.RooFit.PrintLevel(0))
            rooVar_nEventsInNormRegion.Print()
            rooVar_sT.Print()
            extendedPDF.plotOn(sTFrames[kernelType][rhoStr][nJets], ROOT.RooFit.LineColor(ROOT.kBlue), plotRange)
            canvasName = "c_sTUnbinnedFit_{kernelType}_{rhoStr}_{nJets}Jets".format(kernelType=kernelType, rhoStr=rhoStr, nJets=nJets)
            outputFileName = "analysis/plot_sT_UnbinnedFit_{outputFilesString}_{nJetsNorm}JetsNorm_{nJetsMax}JetsMax_{n_sTBins}Bins_{kernelType}_{rhoStr}_{nJets}Jets".format(outputFilesString=inputArguments.outputFilesString, nJetsNorm=inputArguments.nJetsNorm, nJetsMax=inputArguments.nJetsMax, n_sTBins=inputArguments.n_sTBins, kernelType=kernelType, rhoStr=rhoStr, nJets=nJets)
            plotList = [sTFrames[kernelType][rhoStr][nJets]]
            canvases[kernelType][rhoStr][nJets] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = plotList, canvasName = canvasName, outputROOTFile = outputFile, outputDocumentName = outputFileName)
        # Finally, generate and fit new datsets with the fit kernels
        toyRooDataSets = {}
        toyFits[kernelType][rhoStr] = {}
        toy_sTFrames[kernelType][rhoStr] = {}
        toy_sTFrames[kernelType][rhoStr]["DataAndFits"] = rooVar_sT.frame(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
        integralValuesHistogramName = "h_integralValues_{kernelType}_{rhoStr}".format(kernelType=kernelType, rhoStr=rhoStr)
        integralValuesHistogram = ROOT.TH1F(integralValuesHistogramName, integralValuesHistogramName, 40, 0., 0.)
        totalIntegralCheckHistogramName = "h_totalIntegralCheck_{kernelType}_{rhoStr}".format(kernelType=kernelType, rhoStr=rhoStr)
        totalIntegralCheckHistogram = ROOT.TH1F(totalIntegralCheckHistogramName, totalIntegralCheckHistogramName, 40, 0., 0.)
        for counter in range(0, inputArguments.nToyMCs):
            toyRooDataSets[counter] = rooKernel_PDF_Fits[kernelType][rhoStr].generate(ROOT.RooArgSet(rooVar_sT), nEventsInNormJetsBin)
            toyRooDataSets[counter].plotOn(toy_sTFrames[kernelType][rhoStr]["DataAndFits"], plotRange)
            nEntries = toyRooDataSets[counter].numEntries()
            toyFitName = "toyFit_{kernelType}_{rhoStr}_{counter}".format(kernelType=kernelType, rhoStr=rhoStr, counter=counter)
            toyFits[kernelType][rhoStr][counter] = ROOT.RooKeysPdf(toyFitName, toyFitName, rooVar_sT, toyRooDataSets[counter], kernels[kernelType], rho)
            toyFits[kernelType][rhoStr][counter].plotOn(toy_sTFrames[kernelType][rhoStr]["DataAndFits"], plotRange)
            toyExtendedPDFName = "toyExtendedPDF_{kernelType}_{rhoStr}_{counter}".format(kernelType=kernelType, rhoStr=rhoStr, counter=counter)
            toyExtendedPDF = ROOT.RooExtendPdf(toyExtendedPDFName, toyExtendedPDFName, toyFits[kernelType][rhoStr][counter], rooVar_nEventsInNormRegion, "restricted_sTRange")
            toyExtendedPDF.fitTo(toyRooDataSets[counter], normRange, ROOT.RooFit.Minos(ROOT.kTRUE), ROOT.RooFit.PrintLevel(0))
            integralObject = toyExtendedPDF.createIntegral(ROOT.RooArgSet(rooVar_sT), "observation_sTRange")
            integralValuesHistogram.Fill(integralObject.getVal())
            totalIntegralCheckObject = toyExtendedPDF.createIntegral(ROOT.RooArgSet(rooVar_sT), "full_sTRange")
            totalIntegralCheckHistogram.Fill(totalIntegralCheckObject.getVal())
        rooKernel_PDF_Fits[kernelType][rhoStr].plotOn(toy_sTFrames[kernelType][rhoStr]["DataAndFits"], ROOT.RooFit.LineColor(ROOT.kRed), plotRange)
        sTRooDataSets[inputArguments.nJetsNorm].plotOn(toy_sTFrames[kernelType][rhoStr]["DataAndFits"], ROOT.RooFit.LineColor(ROOT.kRed), plotRange)

        # First plot the toy MC data and fits
        canvasName = "c_sT_toyData_{kernelType}_{rhoStr}".format(kernelType=kernelType, rhoStr=rhoStr)
        outputFileName = "analysis/plot_sT_MCToys_{outputFilesString}_{nJetsNorm}JetsNorm_{nJetsMax}JetsMax_{n_sTBins}Bins_{kernelType}_{rhoStr}".format(outputFilesString=inputArguments.outputFilesString, nJetsNorm=inputArguments.nJetsNorm, nJetsMax=inputArguments.nJetsMax, n_sTBins=inputArguments.n_sTBins, kernelType=kernelType, rhoStr=rhoStr)
        plotList = [toy_sTFrames[kernelType][rhoStr]["DataAndFits"]]
        canvases[kernelType][rhoStr]["DataAndFits"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = plotList, canvasName = canvasName, outputROOTFile = outputFile, outputDocumentName = outputFileName)

        # Next plot the systematics estimate
        canvasName = "c_systematics_{kernelType}_{rhoStr}".format(kernelType=kernelType, rhoStr=rhoStr)
        outputFileName = "analysis/plot_systematics_{outputFilesString}_{nJetsNorm}JetsNorm_{nJetsMax}JetsMax_{n_sTBins}Bins_{kernelType}_{rhoStr}".format(outputFilesString=inputArguments.outputFilesString, nJetsNorm=inputArguments.nJetsNorm, nJetsMax=inputArguments.nJetsMax, n_sTBins=inputArguments.n_sTBins, kernelType=kernelType, rhoStr=rhoStr)
        plotList = [integralValuesHistogram]
        canvases[kernelType][rhoStr]["systematics"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = plotList, canvasName = canvasName, outputROOTFile = outputFile, outputDocumentName = outputFileName)

        # Finally plot the integral checks, to see that the normalization is similar
        canvasName = "c_systematicsCheck_{kernelType}_{rhoStr}".format(kernelType=kernelType, rhoStr=rhoStr)
        outputFileName = "analysis/plot_systematicsCheck_{outputFilesString}_{nJetsNorm}JetsNorm_{nJetsMax}JetsMax_{n_sTBins}Bins_{kernelType}_{rhoStr}".format(outputFilesString=inputArguments.outputFilesString, nJetsNorm=inputArguments.nJetsNorm, nJetsMax=inputArguments.nJetsMax, n_sTBins=inputArguments.n_sTBins, kernelType=kernelType, rhoStr=rhoStr)
        plotList = [totalIntegralCheckHistogram]
        canvases[kernelType][rhoStr]["systematicsCheck"] = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = plotList, canvasName = canvasName, outputROOTFile = outputFile, outputDocumentName = outputFileName)

for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    sTTrees[nJets].Write()

outputFile.Write()
outputFile.Close()

sw.Stop()
print ('Real time: ' + str(sw.RealTime() / 60.0) + ' minutes')
print ('CPU time:  ' + str(sw.CpuTime() / 60.0) + ' minutes')
