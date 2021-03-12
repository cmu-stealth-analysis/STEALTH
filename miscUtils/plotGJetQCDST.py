#!/usr/bin/env python

from __future__ import print_function, division
import os, sys, argparse, pdb, math, json, subprocess, array
import tmProgressBar
import stealthEnv, ROOT

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

evtSTEM_minAllowed = 200.0

# normalizeInFirstBin = True
normalizeSTInFirstBin = False

selection = "singlemedium"
blinded = False
outputDirectory = "~/nobackup/analysisAreas/GJetQCDMakeup_singlePhoton"
source_QCD  = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_singlePhotonTrigger_lowerSTThreshold/merged_selection_MC_QCD_singlephoton_2017_control_{s}.root".format(s=selection)
source_GJet = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_singlePhotonTrigger_lowerSTThreshold/merged_selection_MC_GJet17_singlephoton_2017_control_{s}.root".format(s=selection)
source_data = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_singlePhotonTrigger_lowerSTThreshold/merged_selection_data_singlephoton_2017_control_{s}.root".format(s=selection)

# selection = "signal_loose"
# blinded = True
# outputDirectory = "~/nobackup/analysisAreas/GJetQCDMakeup_doublePhoton"
# source_QCD  = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_lowerSTThreshold/merged_selection_MC_QCD*_{s}.root".format(s=selection)
# source_GJet = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_lowerSTThreshold/merged_selection_MC_GJet*_{s}.root".format(s=selection)
# source_data = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_lowerSTThreshold/merged_selection_data*_{s}.root".format(s=selection)

if not(os.path.isdir(outputDirectory)): subprocess.check_call("mkdir -p {oD}".format(oD=outputDirectory), shell=True, executable="/bin/bash")

STRegionBoundariesFileObject = open("STRegionBoundaries.dat")
STBoundaries = []
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()):
        STBoundary = float(STBoundaryString.strip())
        STBoundaries.append(STBoundary)
STBoundaries.append(3500)
nSTSignalBins = len(STBoundaries) - 2 # First two lines are lower and upper boundaries for the normalization bin
n_STBins = len(STBoundaries) - 1

def get_distributions(source, histPrefix, histTitle, getMCWeights):
    inputChain = ROOT.TChain("ggNtuplizer/EventTree")
    inputChain.SetMaxTreeSize(100000000000) # 1 TB
    inputChain.Add(source)

    nEntries = inputChain.GetEntries()
    print("Available nEvts: {n}".format(n=nEntries))
    distributions = {"ST": {}, "PT_leadingJet": {}}
    for nJetsBin in range(2, 7):
        distributions["ST"][nJetsBin] = ROOT.TH1F(histPrefix + "_ST_{n}JetsBin".format(n=nJetsBin), histTitle + ": {n} Jets;ST;Events/GeV".format(n=nJetsBin), n_STBins, array.array('d', STBoundaries))
        distributions["PT_leadingJet"][nJetsBin] = ROOT.TH1F(histPrefix + "_PT_leadingJet_{n}JetsBin".format(n=nJetsBin), histTitle + ": {n} Jets;PT;Events/GeV".format(n=nJetsBin), 5, 50, 1550)
    progressBar = tmProgressBar.tmProgressBar(nEntries)
    progressBarUpdatePeriod = max(1, nEntries//50)
    progressBar.initializeTimer()
    for eventIndex in range(0, nEntries):
        treeStatus = inputChain.LoadTree(eventIndex)
        if (treeStatus < 0):
            break
        evtStatus = inputChain.GetEntry(eventIndex)
        if (evtStatus <= 0):
            continue
        if (eventIndex % progressBarUpdatePeriod == 0): progressBar.updateBar(eventIndex/nEntries, eventIndex)
        ST = inputChain.b_evtST
        nJetsDR = inputChain.b_nJetsDR
        if (ST < STBoundaries[0]): continue
        nJetsBin = min(nJetsDR, 6)
        if (nJetsBin < 2): continue
        # PT_leadingJet = inputChain.b_jetPT_leading
        PT_leadingJet = 200.0 # temporary

        STBinIndex = distributions["ST"][nJetsBin].FindFixBin(ST)
        STBinWidth = distributions["ST"][nJetsBin].GetXaxis().GetBinUpEdge(STBinIndex) - distributions["ST"][nJetsBin].GetXaxis().GetBinLowEdge(STBinIndex)
        eventWeight = 1.0/STBinWidth
        if getMCWeights: eventWeight = inputChain.b_MCCustomWeight/STBinWidth

        evtSTEM = inputChain.b_evtST_electromagnetic

        if ((evtSTEM_minAllowed > 0.) and (evtSTEM <= evtSTEM_minAllowed)): continue

        distributions["ST"][nJetsBin].Fill(ST, eventWeight)

        PTBinIndex = distributions["PT_leadingJet"][nJetsBin].FindFixBin(PT_leadingJet)
        PTBinWidth = (distributions["PT_leadingJet"][nJetsBin].GetXaxis().GetBinUpEdge(PTBinIndex) - distributions["PT_leadingJet"][nJetsBin].GetXaxis().GetBinLowEdge(PTBinIndex))
        eventWeight = 1.0/PTBinWidth
        if getMCWeights: eventWeight = inputChain.b_MCCustomWeight/PTBinWidth
        distributions["PT_leadingJet"][nJetsBin].Fill(PT_leadingJet, eventWeight)
    print()
    return distributions

print("Getting distributions for QCD background...")
print("source: {s}".format(s=source_QCD))
distributions_QCD = get_distributions(source_QCD, "QCD", "", True)
print("Getting distributions for GJet background...")
distributions_GJet = get_distributions(source_GJet, "GJet", "", True)
print("source: {s}".format(s=source_GJet))
distributions_data = None
if not(blinded):
    print("Getting distributions for data...")
    print("source: {s}".format(s=source_data))
    distributions_data = get_distributions(source_data, "data", "", False)

if (normalizeSTInFirstBin and not(blinded)):
    for nJetsBin in range(2, 7):
        nEvents_normBin_data = distributions_data["ST"][nJetsBin].GetBinContent(1)
        nEvents_normBin_QCD  = distributions_QCD["ST"][nJetsBin].GetBinContent(1)
        nEvents_normBin_GJet = distributions_GJet["ST"][nJetsBin].GetBinContent(1)
        try:
            normFactor = (nEvents_normBin_data)/(nEvents_normBin_QCD + nEvents_normBin_GJet)
            distributions_QCD["ST"][nJetsBin].Scale(normFactor)
            distributions_GJet["ST"][nJetsBin].Scale(normFactor)
        except ZeroDivisionError:
            print("WARNING: normalizeSTInFirstBin is set to True, but there are no MC events in the first bin. Not doing any normalization.")

for distributionType in ["ST", "PT_leadingJet"]:
    for nJetsBin in range(2, 7):
        outputCanvas = ROOT.TCanvas("o_{dT}_{n}Jets".format(dT=distributionType, n=nJetsBin), "o_{dT}_{n}Jets".format(dT=distributionType, n=nJetsBin), 1024, 768)
        outputStack = ROOT.THStack()
        distributions_QCD[distributionType][nJetsBin].SetFillColor(ROOT.kBlue)
        outputStack.Add(distributions_QCD[distributionType][nJetsBin])
        distributions_GJet[distributionType][nJetsBin].SetFillColor(ROOT.kRed)
        outputStack.Add(distributions_GJet[distributionType][nJetsBin])
        outputStack.Draw("hist")
        if not(blinded): distributions_data[distributionType][nJetsBin].Draw("SAME")
        titleText = ROOT.TText()
        titleText.SetTextFont(42)
        titleText.SetTextAlign(21)
        titleText.DrawTextNDC(0.5, 0.95, "{dT} distributions, {n} Jets".format(dT=distributionType, n=nJetsBin))
        outputStack.GetXaxis().SetTitle(distributionType)
        outputStack.GetYaxis().SetTitle("Events/GeV")
        ROOT.gPad.SetLogy()
        outputCanvas.Update()
        outputCanvas.SaveAs("{oD}/comparison_{dT}_{s}_{n}Jets.pdf".format(oD=outputDirectory, dT=distributionType, s=selection, n=nJetsBin))
