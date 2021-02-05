#!/usr/bin/env python

from __future__ import print_function, division
import os, sys, argparse, pdb, math, json, subprocess, array
import tmProgressBar
import stealthEnv, ROOT

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

normalizeInFirstBin = True
outputDirectory = "~/nobackup/analysisAreas/GJetQCDSTMakeup_higherNorm"

selection = "singleloose"
blinded = False
source_QCD  = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_tightenedLooseSignal/merged_selection_MC_QCD_singlephoton_2017_control_{s}.root".format(s=selection)
source_GJet = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_tightenedLooseSignal/merged_selection_MC_GJet17_singlephoton_2017_control_{s}.root".format(s=selection)
source_data = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_tightenedLooseSignal/merged_selection_data_singlephoton_2016_control_{s}.root".format(s=selection)

# selection = "control"
# blinded = False
# source_QCD  = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_tightenedLooseSignal/merged_selection_MC_QCD_2017_{s}.root".format(s=selection)
# source_GJet = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_tightenedLooseSignal/merged_selection_MC_GJet17_2017_{s}.root".format(s=selection)
# source_data = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_tightenedLooseSignal/merged_selection_data_2017_{s}.root".format(s=selection)

if not(os.path.isdir(outputDirectory)): subprocess.check_call("mkdir -p {oD}".format(oD=outputDirectory), shell=True, executable="/bin/bash")

STRegionBoundariesFileObject = open("STRegionBoundaries_forSinglePhoton.dat")
STBoundaries = []
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()):
        STBoundary = float(STBoundaryString.strip())
        STBoundaries.append(STBoundary)
STBoundaries.append(3500)
nSTSignalBins = len(STBoundaries) - 2 # First two lines are lower and upper boundaries for the normalization bin
n_STBins = len(STBoundaries) - 1

def get_ST_distributions(source, histPrefix, histTitle, getMCWeights):
    inputFile = ROOT.TFile.Open(source)
    inputTree = ROOT.TTree()
    inputFile.GetObject("ggNtuplizer/EventTree", inputTree)
    nEntries = inputTree.GetEntries()
    print("Available nEvts: {n}".format(n=nEntries))
    STDistributions = {}
    for nJetsBin in range(2, 7):
        STDistributions[nJetsBin] = ROOT.TH1F(histPrefix + "_{n}JetsBin".format(n=nJetsBin), histTitle + ": {n} Jets;ST".format(n=nJetsBin), n_STBins, array.array('d', STBoundaries))
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
        if (ST < STBoundaries[0]): continue
        nJetsBin = min(nJetsDR, 6)
        if (nJetsBin < 2): continue

        STBinIndex = STDistributions[nJetsBin].FindFixBin(ST)
        STBinWidth = STDistributions[nJetsBin].GetXaxis().GetBinUpEdge(STBinIndex) - STDistributions[nJetsBin].GetXaxis().GetBinLowEdge(STBinIndex)
        eventWeight = 1.0/STBinWidth
        if getMCWeights: eventWeight = inputTree.b_MCCustomWeight/STBinWidth

        STDistributions[nJetsBin].Fill(ST, eventWeight)
    print()
    inputFile.Close()
    return STDistributions

print("Getting ST distributions for QCD background...")
STDistributions_QCD = get_ST_distributions(source_QCD, "ST_QCD", "", True)
print("Getting ST distributions for GJet background...")
STDistributions_GJet = get_ST_distributions(source_GJet, "ST_GJet", "", True)
STDistributions_data = None
if not(blinded):
    print("Getting ST distributions for data...")
    STDistributions_data = get_ST_distributions(source_data, "ST_data", "", False)

if (normalizeInFirstBin and not(blinded)):
    for nJetsBin in range(2, 7):
        nEvents_normBin_data = STDistributions_data[nJetsBin].GetBinContent(1)
        nEvents_normBin_QCD  = STDistributions_QCD[nJetsBin].GetBinContent(1)
        nEvents_normBin_GJet = STDistributions_GJet[nJetsBin].GetBinContent(1)
        normFactor = (nEvents_normBin_data)/(nEvents_normBin_QCD + nEvents_normBin_GJet)
        STDistributions_QCD[nJetsBin].Scale(normFactor)
        STDistributions_GJet[nJetsBin].Scale(normFactor)

for nJetsBin in range(2, 7):
    outputCanvas = ROOT.TCanvas("o_{n}Jets".format(n=nJetsBin), "o_{n}Jets".format(n=nJetsBin), 1024, 768)
    outputStack = ROOT.THStack()
    STDistributions_QCD[nJetsBin].SetFillColor(ROOT.kBlue)
    outputStack.Add(STDistributions_QCD[nJetsBin])
    STDistributions_GJet[nJetsBin].SetFillColor(ROOT.kRed)
    outputStack.Add(STDistributions_GJet[nJetsBin])
    outputStack.Draw("hist")
    if not(blinded): STDistributions_data[nJetsBin].Draw("SAME")
    titleText = ROOT.TText()
    titleText.SetTextFont(42)
    titleText.SetTextAlign(21)
    titleText.DrawTextNDC(0.5, 0.95, "ST distributions, {n} Jets".format(n=nJetsBin))
    outputStack.GetXaxis().SetTitle("ST")
    outputStack.GetYaxis().SetTitle("Events/GeV")
    ROOT.gPad.SetLogy()
    outputCanvas.Update()
    outputCanvas.SaveAs("{oD}/comparison_{s}_{n}Jets.pdf".format(oD=outputDirectory, s=selection, n=nJetsBin))
