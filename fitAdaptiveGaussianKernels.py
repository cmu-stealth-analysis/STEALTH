#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, ROOT, argparse, array
import numpy as np
from tmProgressBar import tmProgressBar

inputArgumentsParser = argparse.ArgumentParser(description='Run STEALTH selection.')
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--n_sTBins', default=8, help='Number of sT bins.',type=int)
inputArgumentsParser.add_argument('--sTMin', default=800., help='Min value of sT to plot. This is also used as the min of the range in the fit performed by the kernel estimator in the nJets bin to which to normalize.',type=float)
inputArgumentsParser.add_argument('--sTMax', default=1600., help='Max value of sT to plot. This is also used as the max of the range in the fit performed by the kernel estimator in the nJets bin to which to normalize.',type=float)
inputArgumentsParser.add_argument('--sTNormRangeMin', default=800., help='Min value of sT for normalization. For all sT distributions in nJets bins except the normalization bin, this value is the min of the range in which to scale the kernel estimator.',type=float)
inputArgumentsParser.add_argument('--sTNormRangeMax', default=900., help='Max value of sT for normalization. For all sT distributions in nJets bins except the normalization bin, this value is the max of the range in which to scale the kernel estimator.',type=float)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Min number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Max number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsNorm', default=3, help='Number of jets w.r.t. which to normalize the sT distributions for other jets.',type=int)
inputArgumentsParser.add_argument('--outputFilesSuffix', required=True, help='Prefix for output files.',type=str)
inputArguments = inputArgumentsParser.parse_args()
if (inputArguments.sTNormRangeMin < inputArguments.sTMin or inputArguments.sTNormRangeMax > inputArguments.sTMax):
    print ("Normalization interval: ({nmin}, {nmax}) is incompatible with sT range: ({smin, smax})".format(nmin=inputArguments.sTNormRangeMin, nmax=inputArguments.sTNormRangeMax, smin=inputArguments.sTMin, smax=inputArguments.sTMax))

rooVar_sT = ROOT.RooRealVar("rooVar_sT", "rooVar_sT", inputArguments.sTMin, inputArguments.sTMax)

kernels = {"NoMirror": ROOT.RooKeysPdf.NoMirror, "MirrorLeft": ROOT.RooKeysPdf.MirrorLeft, "MirrorRight": ROOT.RooKeysPdf.MirrorRight, "MirrorBoth": ROOT.RooKeysPdf.MirrorBoth, "MirrorAsymLeft": ROOT.RooKeysPdf.MirrorAsymLeft, "MirrorAsymLeftRight": ROOT.RooKeysPdf.MirrorAsymLeftRight, "MirrorAsymRight": ROOT.RooKeysPdf.MirrorAsymRight, "MirrorLeftAsymRight": ROOT.RooKeysPdf.MirrorLeftAsymRight, "MirrorAsymBoth": ROOT.RooKeysPdf.MirrorAsymBoth}
enabledKernels = ["MirrorLeftAsymRight"]
enabledRhos = [0.5, 0.8, 1.0, 1.2]

for enabledKernel in enabledKernels:
    if not(enabledKernel in kernels.keys()): sys.exit("The following element is present in list of enabled kernels but not in the dictionary defining the correspondence between kernel name and RooKeysPdf index: " + enabledKernel)

sw = ROOT.TStopwatch()
sw.Start()

inputChain = ROOT.TChain('ggNtuplizer/EventTree')
inputChain.Add(inputArguments.inputFilePath)
n_entries = inputChain.GetEntries()
print ('Total number of available events: ' + str(n_entries))

outputFile = ROOT.TFile('analysis/sTDistributions_{outputFilesSuffix}_{nJetsNorm}JetsNorm_{nJetsMax}JetsMax_{n_sTBins}Bins.root'.format(outputFilesSuffix=inputArguments.outputFilesSuffix, nJetsNorm=inputArguments.nJetsNorm, nJetsMax=inputArguments.nJetsMax, n_sTBins=inputArguments.n_sTBins), 'recreate')

# Initialize TTrees
sTTrees = {}
restricted_sTTrees = {}
sTArrays = {}
restricted_sTArrays = {}
for nJets in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    treeName = "sTTree_{nJets}Jets".format(nJets=nJets)
    sTTrees[nJets] = ROOT.TTree(treeName, treeName)
    treeName = "restricted_" + treeName
    restricted_sTTrees[nJets] = ROOT.TTree(treeName, treeName)
    sTArrays[nJets] = array.array('f', [0.])
    restricted_sTArrays[nJets] = array.array('f', [0.])
    (sTTrees[nJets]).Branch('rooVar_sT', (sTArrays[nJets]), 'rooVar_sT/F')
    (restricted_sTTrees[nJets]).Branch('rooVar_sT', (restricted_sTArrays[nJets]), 'rooVar_sT/F')

# Fill histograms
progressBar = tmProgressBar(n_entries)
progressBarUpdatePeriod = n_entries//1000
progressBar.initializeTimer()
for entryIndex in range(n_entries):
    entryStatus = inputChain.LoadTree(entryIndex)
    if entryStatus < 0: sys.exit("Tree failed to load entry at index {entryIndex}; returned status {entryStatus}".format(entryIndex=entryIndex, entryStatus=entryStatus))
    chainStatus = inputChain.GetEntry(entryIndex)
    if chainStatus <= 0: sys.exit("Unable to load data from chain at index {entryIndex}; GetEntry returned status {chainStatus}".format(entryIndex=entryIndex, chainStatus=chainStatus))

    if (entryIndex%progressBarUpdatePeriod == 0): progressBar.updateBar(1.0*entryIndex/n_entries, entryIndex)

    nStealthJets = inputChain.b_nJets
    nJetsBin = nStealthJets
    if (nJetsBin > inputArguments.nJetsMax): nJetsBin = inputArguments.nJetsMax
    if (nJetsBin < inputArguments.nJetsMin):
        print("nJetsBin = {nJetsBin} at entry index = {entryIndex}".format(nJetsBin=nJetsBin, entryIndex=entryIndex))
        continue

    sT = inputChain.b_evtST
    if (sT > inputArguments.sTMin and sT < inputArguments.sTMax):
        (sTArrays[nJetsBin])[0] = sT
        (sTTrees[nJetsBin]).Fill()
        if (sT > inputArguments.sTNormRangeMin and sT < inputArguments.sTNormRangeMax):
            (restricted_sTArrays[nJetsBin])[0] = sT
            (restricted_sTTrees[nJetsBin]).Fill()

progressBar.terminate()

# First find the kernel fits in background nJets bin
normTree = sTTrees[inputArguments.nJetsNorm]
dataSetName = "rooDataSet_" + str(inputArguments.nJetsNorm) + "Jets"
normNJetsRooDataSet = ROOT.RooDataSet(dataSetName, dataSetName, normTree, ROOT.RooArgSet(rooVar_sT))
rooKernel_PDF_Fits = {}
print("Performing fits for background nJets bin")
for kernelType in enabledKernels:
    print("Defining fit for kernel type: " + kernelType)
    for rho in enabledRhos:
        rhoStr = ("rho_{rho:2.1f}".format(rho=rho)).replace('.', 'pt')
        print("Defining fit for rho: " + rhoStr)
        functionName = "rooKernelFunction_{kernelType}_{rhoStr}".format(kernelType=kernelType, rhoStr=rhoStr)
        rooKernel_PDF_Fits[functionName] = ROOT.RooKeysPdf(functionName, functionName, rooVar_sT, normNJetsRooDataSet, kernels[kernelType], rho)
        sTFrame = rooVar_sT.frame(inputArguments.sTMin, inputArguments.sTMax, inputArguments.n_sTBins)
        normNJetsRooDataSet.plotOn(sTFrame)
        rooKernel_PDF_Fits[functionName].plotOn(sTFrame)
        c_sTUnbinnedFit = ROOT.TCanvas('c_sTUnbinnedFit_' + kernelType + "_" + rhoStr + "_norm", 'c_sTUnbinnedFit_' + kernelType + "_" + rhoStr + "_norm", 1024, 768)
        c_sTUnbinnedFit.SetBorderSize(0)
        c_sTUnbinnedFit.SetFrameBorderMode(0)
        sTFrame.Draw()
        c_sTUnbinnedFit.SaveAs("analysis/plot_sT_UnbinnedFit_{outputFilesSuffix}_{nJetsNorm}JetsNorm_{nJetsMax}JetsMax_{n_sTBins}Bins_{kernelType}_{rhoStr}_norm.png".format(outputFilesSuffix=inputArguments.outputFilesSuffix, nJetsNorm=inputArguments.nJetsNorm, nJetsMax=inputArguments.nJetsMax, n_sTBins=inputArguments.n_sTBins, kernelType=kernelType, rhoStr=rhoStr))
        c_sTUnbinnedFit.Write()

# print("Performing fits and scaling in non-background nJet bins")
# rooVar_scales = {}
# for nJets in range(inputArguments.nJetsMin, inputArguments.nJetsMax + 1):
#     if (nJets == inputArguments.nJetsNorm): continue
#     restricted_sTTree = restricted_sTTrees[nJets]
#     restricted_sTRooDataSet = ROOT.RooDataSet("restrictedInputData_" + str(i) + "Jets", "restrictedInputData_" + str(i) + "Jets", restricted_sTTree, ROOT.RooArgSet(rooVar_sT))
#     for kernelType in enabledKernels:
#         print("Defining fit for kernel type: " + kernelType)
#         for rho in enabledRhos:
#             rhoStr = ("rho_{rho:2.1f}".format(rho=rho)).replace('.', 'pt')
#             print("Defining fit for rho: " + rhoStr)
#             restricted_sTFrame = rooVar_sT.frame(inputArguments.sTNormRangeMin, inputArguments.sTNormRangeMax, 1)
#             restricted_sTRooDataSet.plotOn(restricted_sTFrame)
#             unscaledFunctionName = "rooKernelFunction_" + kernelType + "_" + rhoStr
#             clonedFunctionName = unscaledFunctionName + "_cloned_{i}Jets".format(i=i)
#             clonedFunction = (rooKernel_PDF_Fits[unscaledFunctionName]).clone(clonedFunctionName)
#             clonedFunction_asTF = clonedFunction.asTF(ROOT.RooArgList(rooVar_sT))
#             # clonedFunctionName = clonedFunction_asTF.GetName()
#             # lambdaExpressionForScaledFunction = "[&](double *sT, double *scale){return scale[0]*" + clonedFunctionName + "(sT)};"
#             lambdaExpressionForScaledFunction = "[&](double *x, double *p){ return p[0]*clonedFunction_asTF(x); }"
#             print("Lambda expression:")
#             print(lambdaExpressionForScaledFunction)
#             print("Checking clonedFunction_asTF:")
#             print(str(clonedFunction_asTF(900.)))
#             scaledTF1 = ROOT.TF1("scaledTF1_{kernelType}_{rhoStr}_{i}Jets".format(kernelType=kernelType, rhoStr=rhoStr, i=i), lambdaExpressionForScaledFunction, inputArguments.sTNormRangeMin, inputArguments.sTNormRangeMax, 1)
#             canvasName = 'c_sTUnbinnedFit_{kernelType}_{rhoStr}'.format(kernelType=kernelType, rhoStr=rhoStr)
#             c_sTUnbinnedFit = ROOT.TCanvas(canvasName, canvasName, 1024, 768)
#             c_sTUnbinnedFit.SetBorderSize(0)
#             c_sTUnbinnedFit.SetFrameBorderMode(0)
#             restricted_sTFrame.Draw()
#             c_sTUnbinnedFit.SaveAs("analysis/plot_sT_UnbinnedFit_{outputFilesSuffix}_{nJetsNorm}JetsNorm_{nJetsMax}JetsMax_{n_sTBins}_Bins_{kernelType}_{rhoStr}_{nJets}Jets.png".format(outputFilesSuffix=inputArguments.outputFilesSuffix, kernelType=kernelType, rhoStr=rhoStr, nJets=nJets))
#             c_sTUnbinnedFit.Write()

for nJets in range(inputArguments.nJetsMin, inputArguments.nJetsMax + 1):
    (sTTrees[nJets]).Write()
    (restricted_sTTrees[nJets]).Write()

outputFile.Write()
outputFile.Close()

sw.Stop()
print ('Real time: ' + str(sw.RealTime() / 60.0) + ' minutes')
print ('CPU time:  ' + str(sw.CpuTime() / 60.0) + ' minutes')
