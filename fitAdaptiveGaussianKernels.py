#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, ROOT, argparse, array
import numpy as np
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
inputArgumentsParser.add_argument('--outputFilesString', required=True, help='String to include in all output file names.',type=str)
inputArguments = inputArgumentsParser.parse_args()
if (inputArguments.sTNormRangeMin < inputArguments.sTKernelFitRangeMin or inputArguments.sTNormRangeMax > inputArguments.sTKernelFitRangeMax):
    print ("Normalization interval: ({nmin}, {nmax}) seems incompatible with kernel fitting range: ({smin, smax})".format(nmin=inputArguments.sTNormRangeMin, nmax=inputArguments.sTNormRangeMax, smin=inputArguments.sTKernelFitRangeMin, smax=inputArguments.sTKernelFitRangeMax))

rooVar_sT = ROOT.RooRealVar("rooVar_sT", "rooVar_sT", inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax, "GeV")
rooVar_sT.setRange("restricted_sTRange", inputArguments.sTNormRangeMin, inputArguments.sTNormRangeMax)
rooVar_sT.setRange("full_sTRange", inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
plotRange = ROOT.RooFit.Range(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax)
kernelFitRange = ROOT.RooFit.Range(inputArguments.sTKernelFitRangeMin, inputArguments.sTKernelFitRangeMax)
normRange = ROOT.RooFit.Range(inputArguments.sTNormRangeMin, inputArguments.sTNormRangeMax)

kernels = {"NoMirror": ROOT.RooKeysPdf.NoMirror, "MirrorLeft": ROOT.RooKeysPdf.MirrorLeft, "MirrorRight": ROOT.RooKeysPdf.MirrorRight, "MirrorBoth": ROOT.RooKeysPdf.MirrorBoth, "MirrorAsymLeft": ROOT.RooKeysPdf.MirrorAsymLeft, "MirrorAsymLeftRight": ROOT.RooKeysPdf.MirrorAsymLeftRight, "MirrorAsymRight": ROOT.RooKeysPdf.MirrorAsymRight, "MirrorLeftAsymRight": ROOT.RooKeysPdf.MirrorLeftAsymRight, "MirrorAsymBoth": ROOT.RooKeysPdf.MirrorAsymBoth}
enabledKernels = ["MirrorLeftAsymRight"]
enabledRhos = [0.5, 0.8, 1.0, 1.2]

for enabledKernel in enabledKernels:
    if not(enabledKernel in kernels.keys()): sys.exit("The following element is present in list of enabled kernels but not in the dictionary defining the correspondence between kernel name and RooKeysPdf index: {enabledKernel}".format(enabledKernel=enabledKernel))

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

# First find the kernel fits in background nJets bin
normTree = sTTrees[inputArguments.nJetsNorm]
dataSetName = "rooDataSet_{nJetsNorm}Jets".format(nJetsNorm=inputArguments.nJetsNorm)
normNJetsRooDataSet = ROOT.RooDataSet(dataSetName, dataSetName, normTree, ROOT.RooArgSet(rooVar_sT))
rooKernel_PDF_Fits = {}
print("Performing fits for background nJets bin")
for kernelType in enabledKernels:
    print("Defining fit for kernel type: {kernelType}".format(kernelType=kernelType))
    rooKernel_PDF_Fits[kernelType] = {}
    for rho in enabledRhos:
        rhoStr = ("rho_{rho:2.1f}".format(rho=rho)).replace('.', 'pt')
        print("Defining fit for rho: {rhoStr}".format(rhoStr=rhoStr))
        functionName = "rooKernelFunction_{kernelType}_{rhoStr}".format(kernelType=kernelType, rhoStr=rhoStr)
        (rooKernel_PDF_Fits[kernelType])[rhoStr] = ROOT.RooKeysPdf(functionName, functionName, rooVar_sT, normNJetsRooDataSet, kernels[kernelType], rho)
        sTFrame = rooVar_sT.frame(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
        normNJetsRooDataSet.plotOn(sTFrame, plotRange)
        (rooKernel_PDF_Fits[kernelType][rhoStr]).plotOn(sTFrame, plotRange)
        # normNJetsRooDataSet.plotOn(sTFrame)
        # (rooKernel_PDF_Fits[kernelType][rhoStr]).plotOn(sTFrame)
        c_sTUnbinnedFitName = "c_sTUnbinnedFit_{kernelType}_{rhoStr}_norm".format(kernelType=kernelType, rhoStr=rhoStr)
        c_sTUnbinnedFit = ROOT.TCanvas(c_sTUnbinnedFitName, c_sTUnbinnedFitName, 1024, 768)
        c_sTUnbinnedFit.SetBorderSize(0)
        c_sTUnbinnedFit.SetFrameBorderMode(0)
        sTFrame.Draw()
        c_sTUnbinnedFit.SaveAs("analysis/plot_sT_UnbinnedFit_{outputFilesString}_{nJetsNorm}JetsNorm_{nJetsMax}JetsMax_{n_sTBins}Bins_{kernelType}_{rhoStr}_norm.png".format(outputFilesString=inputArguments.outputFilesString, nJetsNorm=inputArguments.nJetsNorm, nJetsMax=inputArguments.nJetsMax, n_sTBins=inputArguments.n_sTBins, kernelType=kernelType, rhoStr=rhoStr))
        c_sTUnbinnedFit.Write()

print("Performing fits and scaling in non-background nJet bins")
rooVar_scales = {}
rooVar_nEventsInNormRegion = ROOT.RooRealVar("rooVar_nEventsInNormRegion", "rooVar_nEventsInNormRegion", 20, 0, 10000.)
for nJets in range(inputArguments.nJetsMin, inputArguments.nJetsMax + 1):
    if (nJets == inputArguments.nJetsNorm): continue
    sTRooDataSetName = "rooDataSet_{nJets}Jets".format(nJets=nJets)
    sTRooDataSet = ROOT.RooDataSet(sTRooDataSetName, sTRooDataSetName, sTTrees[nJets], ROOT.RooArgSet(rooVar_sT))
    for kernelType in enabledKernels:
        print("Defining fit for kernel type: {kernelType}".format(kernelType=kernelType))
        for rho in enabledRhos:
            rhoStr = ("rho_{rho:2.1f}".format(rho=rho)).replace('.', 'pt')
            print("Defining fit for rho: {rhoStr}".format(rhoStr=rhoStr))
            sTFrame = rooVar_sT.frame(inputArguments.sTPlotRangeMin, inputArguments.sTPlotRangeMax, inputArguments.n_sTBins)
            sTRooDataSet.plotOn(sTFrame, plotRange)
            (rooKernel_PDF_Fits[kernelType][rhoStr]).plotOn(sTFrame, ROOT.RooFit.LineColor(ROOT.kBlue), plotRange)
            # sTRooDataSet.plotOn(sTFrame)
            # (rooKernel_PDF_Fits[kernelType][rhoStr]).plotOn(sTFrame, ROOT.RooFit.LineColor(ROOT.kBlue))
            extendedPDFName = "extendedPDF_{kernelType}_{rhoStr}_{nJets}Jets".format(kernelType=kernelType, rhoStr=rhoStr, nJets=nJets)
            extendedPDF = ROOT.RooExtendPdf(extendedPDFName, extendedPDFName, rooKernel_PDF_Fits[kernelType][rhoStr], rooVar_nEventsInNormRegion, "restricted_sTRange")
            extendedPDF.fitTo(sTRooDataSet, normRange, ROOT.RooFit.Minos(ROOT.kTRUE))
            rooVar_nEventsInNormRegion.Print()
            rooVar_sT.Print()
            extendedPDF.plotOn(sTFrame, ROOT.RooFit.LineColor(ROOT.kRed), plotRange)
            # extendedPDF.plotOn(sTFrame, ROOT.RooFit.LineColor(ROOT.kRed))
            c_sTUnbinnedFitName = "c_sTUnbinnedFit_{kernelType}_{rhoStr}_{nJets}Jets".format(kernelType=kernelType, rhoStr=rhoStr, nJets=nJets)
            c_sTUnbinnedFit = ROOT.TCanvas(c_sTUnbinnedFitName, c_sTUnbinnedFitName, 1024, 768)
            c_sTUnbinnedFit.SetBorderSize(0)
            c_sTUnbinnedFit.SetFrameBorderMode(0)
            sTFrame.Draw()
            c_sTUnbinnedFit.SaveAs("analysis/plot_sT_UnbinnedFit_{outputFilesString}_{nJetsNorm}JetsNorm_{nJetsMax}JetsMax_{n_sTBins}Bins_{kernelType}_{rhoStr}_{nJets}Jets.png".format(outputFilesString=inputArguments.outputFilesString, nJetsNorm=inputArguments.nJetsNorm, nJetsMax=inputArguments.nJetsMax, n_sTBins=inputArguments.n_sTBins, kernelType=kernelType, rhoStr=rhoStr, nJets=nJets))
            c_sTUnbinnedFit.Write()

for nJets in range(inputArguments.nJetsMin, inputArguments.nJetsMax + 1):
    (sTTrees[nJets]).Write()

outputFile.Write()
outputFile.Close()

sw.Stop()
print ('Real time: ' + str(sw.RealTime() / 60.0) + ' minutes')
print ('CPU time:  ' + str(sw.CpuTime() / 60.0) + ' minutes')
