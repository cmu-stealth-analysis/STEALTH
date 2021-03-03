#!/usr/bin/env python

from __future__ import print_function, division
import os, sys, argparse, pdb, math, json, subprocess, array
import tmProgressBar
import stealthEnv, ROOT

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

STMin = 700.
STMax = 3500.
STBoundariesSourceFile = "STRegionBoundaries_normOptimization.dat"

colors = {
    2: ROOT.kBlack,
    3: ROOT.kBlue+2,
    4: ROOT.kRed+1,
    5: ROOT.kGreen+3,
    6: ROOT.kViolet,
}

# selection = "singlemedium"
# identifier = "data"
# year = "2017"
# yearPattern = "{y}".format(y=year)
# if (year == "all"): yearPattern = "*"
# sourceFilePattern  = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_singlePhotonTrigger_lowerSTThreshold/merged_selection_{i}_singlephoton_{yP}_control_{s}.root".format(i=identifier, s=selection, yP=yearPattern)
# getMCWeights = True
# if (identifier[0:4] == "data"): getMCWeights = False
# outputDirectory = stealthEnv.analysisRoot + "/STDistributions_singlephoton"
# evtSTEM_minAllowed = 200.

selection = "signal"
identifier = "MC_GJet"
year = "all"
yearPattern = "{y}".format(y=year)
if (year == "all"): yearPattern = "*"
sourceFilePattern  = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_lowerSTThreshold/merged_selection_{i}_{yP}_{s}.root".format(i=identifier, s=selection, yP=yearPattern)
if (identifier == "MC_GJet"):
    if (year == "all"):
        sourceFilePattern  = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_lowerSTThreshold/merged_selection_MC_GJet*_{s}.root".format(i=identifier, s=selection, yP=yearPattern)
    else:
        sys.exit("ERROR: Unrecognized (year, identifier) combo: ({y}, {i})".format(y=year, i=identifier))
getMCWeights = True
if (identifier[0:4] == "data"): getMCWeights = False
outputDirectory = stealthEnv.analysisRoot + "/STDistributions_doublephoton"
evtSTEM_minAllowed = -1.0

if not(os.path.isdir(outputDirectory)): subprocess.check_call("mkdir -p {oD}".format(oD=outputDirectory), shell=True, executable="/bin/bash")

STRegionBoundariesFileObject = open(STBoundariesSourceFile)
STBoundaries = []
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()):
        STBoundary = float(STBoundaryString.strip())
        STBoundaries.append(STBoundary)
STBoundaries.append(3500)
nSTSignalBins = len(STBoundaries) - 2 # First two lines are lower and upper boundaries for the normalization bin
n_STBins = len(STBoundaries) - 1

print("Getting ST datasets for source: {sF}".format(sF=sourceFilePattern))

inputChain = ROOT.TChain("ggNtuplizer/EventTree")
inputChain.SetMaxTreeSize(100000000000) # 1 TB
inputChain.Add(sourceFilePattern)

nEntries = inputChain.GetEntries()
print("Available nEvts: {n}".format(n=nEntries))

outputFile = ROOT.TFile.Open("{oD}/distributions_{y}_{s}_{i}.root".format(oD=outputDirectory, y=year, i=identifier, s=selection), "RECREATE")
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
    if (eventIndex % progressBarUpdatePeriod == 0): progressBar.updateBar(eventIndex/nEntries, eventIndex)
    treeStatus = inputChain.LoadTree(eventIndex)
    if (treeStatus < 0):
        break
    evtStatus = inputChain.GetEntry(eventIndex)
    if (evtStatus <= 0):
        continue
    ST = inputChain.b_evtST
    nJetsDR = inputChain.b_nJetsDR
    if ((ST < STMin) or (ST > STMax)): continue
    nJetsBin = min(nJetsDR, 6)
    if (nJetsBin < 2): continue

    eventWeight = 1.0
    if getMCWeights: eventWeight = inputChain.b_MCCustomWeight

    evtSTEM = inputChain.b_evtST_electromagnetic
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

outputFile.Close()

print("Done!")