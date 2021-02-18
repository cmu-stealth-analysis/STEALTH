#!/usr/bin/env python

from __future__ import print_function, division
import os, sys, argparse, pdb, math, json, subprocess, array
import tmProgressBar
import stealthEnv, ROOT

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

evtSTEM_minAllowed = -1.0

STMin = 700.
STMax = 3500.
nSTBins = 28

colors = {
    2: ROOT.kBlack,
    3: ROOT.kBlue+2,
    4: ROOT.kRed+1,
    5: ROOT.kGreen+3,
    6: ROOT.kViolet,
}

# selection = "singlemedium"
# identifier = "data"
# sourceFile  = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_tightenedLooseSignal_singlePhotonTrigger_lowerSTThreshold/merged_selection_{i}_singlephoton_2017_control_{s}.root".format(i=identifier, s=selection)
# getMCWeights = False
# outputDirectory = "~/nobackup/analysisAreas/STDistributions_singlephoton"
# STBoundariesSourceFile = "STRegionBoundaries_forNormOptimization.dat"

selection = "control"
identifier = "data_control"
sourceFile  = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_tightenedLooseSignal_singlePhotonTrigger_lowerSTThreshold/merged_selection_data_2017_{s}.root".format(s=selection)
getMCWeights = False
outputDirectory = "~/nobackup/analysisAreas/STDistributions_doublephoton"
STBoundariesSourceFile = "STRegionBoundaries_forNormOptimization_wider.dat"

if not(os.path.isdir(outputDirectory)): subprocess.check_call("mkdir -p {oD}".format(oD=outputDirectory), shell=True, executable="/bin/bash")

STRegionBoundariesFileObject = open(STBoundariesSourceFile)
STBoundaries = []
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()):
        STBoundary = float(STBoundaryString.strip())
        STBoundaries.append(STBoundary)
# STBoundaries.append(3500)
nSTSignalBins = len(STBoundaries) - 2 # First two lines are lower and upper boundaries for the normalization bin
n_STBins = len(STBoundaries) - 1

print("Getting ST datasets for source: {sF}".format(sF=sourceFile))

inputFile = ROOT.TFile.Open(sourceFile)

inputTree = ROOT.TTree()
inputFile.GetObject("ggNtuplizer/EventTree", inputTree)
nEntries = inputTree.GetEntries()
print("Available nEvts: {n}".format(n=nEntries))

outputFile = ROOT.TFile.Open("{oD}/distributions_{s}_{i}.root".format(oD=outputDirectory, i=identifier, s=selection), "RECREATE")
ROOT.RooAbsData.setDefaultStorageType(ROOT.RooAbsData.Tree)

rooSTVar = ROOT.RooRealVar("rooVar_sT", "rooVar_sT", STMin, STMax, "GeV")
rooWeightVar = ROOT.RooRealVar("rooVar_weight", "rooVar_weight", 0., 100000., "")
dataSets = {}
STDistributions = {}
for nJetsBin in range(2, 7):
    STDistributions[nJetsBin] = ROOT.TH1F("h_ST_{n}JetsBin".format(n=nJetsBin), "ST distribution: {n} Jets;ST".format(n=nJetsBin), n_STBins, array.array('d', STBoundaries))
    STDistributions[nJetsBin].Sumw2()
    dataSets[nJetsBin] = ROOT.RooDataSet("dataSet_{n}Jets".format(n=nJetsBin), "dataSet_{n}Jets".format(n=nJetsBin), ROOT.RooArgSet(rooSTVar, rooWeightVar), "rooVar_weight")

progressBar = tmProgressBar.tmProgressBar(nEntries)
progressBarUpdatePeriod = max(1, nEntries//50)
progressBar.initializeTimer()
for eventIndex in range(0, nEntries):
    evtStatus = inputTree.GetEntry(eventIndex)
    if (evtStatus <= 0):
        continue
    if (eventIndex % progressBarUpdatePeriod == 0): progressBar.updateBar(eventIndex/nEntries, eventIndex)
    ST = inputTree.b_evtST
    nJetsDR = inputTree.b_nJetsDR
    if ((ST < STMin) or (ST > STMax)): continue
    nJetsBin = min(nJetsDR, 6)
    if (nJetsBin < 2): continue

    eventWeight = 1.0
    if getMCWeights: eventWeight = inputTree.b_MCCustomWeight

    evtSTEM = inputTree.b_evtST_electromagnetic
    if ((evtSTEM_minAllowed > 0.) and (evtSTEM <= evtSTEM_minAllowed)): continue

    STBinIndex = STDistributions[nJetsBin].FindFixBin(ST)
    STBinWidth = STDistributions[nJetsBin].GetXaxis().GetBinUpEdge(STBinIndex) - STDistributions[nJetsBin].GetXaxis().GetBinLowEdge(STBinIndex)
    eventWeight_histograms = eventWeight/STBinWidth

    rooSTVar.setVal(ST)
    dataSets[nJetsBin].add(ROOT.RooArgSet(rooSTVar), eventWeight)
    STDistributions[nJetsBin].Fill(ST, eventWeight_histograms)

print()
for nJetsBin in range(2, 7):
    outputFile.WriteTObject(STDistributions[nJetsBin])
    outputFile.WriteTObject(dataSets[nJetsBin])

inputFile.Close()
outputFile.Close()

print("Done!")
