#!/usr/bin/env python

from __future__ import print_function, division
import os, sys, argparse, pdb, math, json, subprocess, array
import tmProgressBar
import stealthEnv, ROOT

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

optional_identifier_string = "_newBkgMC"

# selection = "singleloose"
# blinded = False
# outputDirectory = "~/nobackup/analysisAreas/GJetQCDMakeup_singlePhoton"
# sources_QCD  = [
#     stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_MC_QCD18_singlephoton_2018_control_{s}.root".format(s=selection, ois=optional_identifier_string),
#     stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_MC_QCD17_singlephoton_2017_control_{s}.root".format(s=selection, ois=optional_identifier_string),
#     stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_MC_QCD16_singlephoton_2016_control_{s}.root".format(s=selection, ois=optional_identifier_string)
# ]
# sources_GJet = [
#     stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_MC_GJet18_singlephoton_2018_control_{s}.root".format(s=selection, ois=optional_identifier_string),
#     stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_MC_GJet17_singlephoton_2017_control_{s}.root".format(s=selection, ois=optional_identifier_string),
#     stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_MC_GJet16_singlephoton_2016_control_{s}.root".format(s=selection, ois=optional_identifier_string)
# ]
# sources_data = [
#     stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_data_singlephoton_2018_control_{s}.root".format(s=selection, ois=optional_identifier_string),
#     stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_data_singlephoton_2017_control_{s}.root".format(s=selection, ois=optional_identifier_string),
#     stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_data_singlephoton_2016_control_{s}.root".format(s=selection, ois=optional_identifier_string)
# ]
# normalizeSTInFirstBin = False
# evtSTEM_minAllowed = 200.0

selection = "signal_loose"
blinded = True
outputDirectory = "~/nobackup/analysisAreas/BackgroundMC_STDistributions"
normalizeSTInFirstBin = False
evtSTEM_minAllowed = -1.0

sources_diphoton  = [
    stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_MC_DiPhotonJets_2016_{s}.root".format(s=selection, ois=optional_identifier_string),
    stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_MC_DiPhotonJets_2017_{s}.root".format(s=selection, ois=optional_identifier_string),
    stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_MC_DiPhotonJets_2018_{s}.root".format(s=selection, ois=optional_identifier_string),
]
sources_GJet = [
    stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_MC_GJetHT16_2016_{s}.root".format(s=selection, ois=optional_identifier_string),
    stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_MC_GJetHT17_2017_{s}.root".format(s=selection, ois=optional_identifier_string),
    stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_MC_GJetHT18_2018_{s}.root".format(s=selection, ois=optional_identifier_string)
]
sources_QCD = [
    stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_MC_HighHTQCD16_2016_{s}.root".format(s=selection, ois=optional_identifier_string),
    stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_MC_HighHTQCD17_2017_{s}.root".format(s=selection, ois=optional_identifier_string),
    stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_MC_HighHTQCD18_2018_{s}.root".format(s=selection, ois=optional_identifier_string)
]
sources_data = [
    stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_data_2016_{s}.root".format(s=selection, ois=optional_identifier_string),
    stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_data_2017_{s}.root".format(s=selection, ois=optional_identifier_string),
    stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton{ois}/merged_selection_data_2018_{s}.root".format(s=selection, ois=optional_identifier_string)
]

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

def get_distributions(sources, histPrefix, histTitle, getMCWeights):
    inputChain = ROOT.TChain("ggNtuplizer/EventTree")
    inputChain.SetMaxTreeSize(100000000000) # 1 TB
    for source in sources:
        inputChain.Add(source)

    nEntries = inputChain.GetEntries()
    print("Available nEvts: {n}".format(n=nEntries))
    distributions = {"ST": {}}
    for nJetsBin in range(2, 7):
        distributions["ST"][nJetsBin] = ROOT.TH1D(histPrefix + "_ST_{n}JetsBin".format(n=nJetsBin), histTitle + ": {n} Jets;ST;Events/GeV".format(n=nJetsBin), n_STBins, array.array('d', STBoundaries))
        distributions["ST"][nJetsBin].Sumw2()
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
        # PT_leadingJet = 200.0 # temporary

        STBinIndex = distributions["ST"][nJetsBin].FindFixBin(ST)
        STBinWidth = distributions["ST"][nJetsBin].GetXaxis().GetBinUpEdge(STBinIndex) - distributions["ST"][nJetsBin].GetXaxis().GetBinLowEdge(STBinIndex)
        eventWeight = 1.0/STBinWidth
        if getMCWeights:
            eventWeight *= inputChain.b_MCXSecWeight
            eventWeight *= inputChain.genWeight
            eventWeight *= inputChain.b_evtPrefiringWeight
            eventWeight *= inputChain.b_evtphotonMCScaleFactor
            eventWeight *= inputChain.b_PUWeightNoSelection

        evtSTEM = inputChain.b_evtST_electromagnetic

        if ((evtSTEM_minAllowed > 0.) and (evtSTEM <= evtSTEM_minAllowed)): continue

        distributions["ST"][nJetsBin].Fill(ST, eventWeight)
    progressBar.terminate()
    return distributions

print("Getting distributions for diphoton background...")
print("sources: {s}".format(s=sources_diphoton))
distributions_diphoton = get_distributions(sources_diphoton, "diphoton", "Diphoton", True)
print("Getting distributions for GJet background...")
print("sources: {s}".format(s=sources_GJet))
distributions_GJet = get_distributions(sources_GJet, "GJet", "GJet", True)
print("Getting distributions for QCD background...")
print("sources: {s}".format(s=sources_QCD))
distributions_QCD = get_distributions(sources_QCD, "QCD", "QCD", True)

distributions_data = None
if not(blinded):
    print("Getting distributions for data...")
    print("sources: {s}".format(s=sources_data))
    distributions_data = get_distributions(sources_data, "data", "", False)

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

for distributionType in ["ST"# , "PT_leadingJet"
]:
    for nJetsBin in range(2, 7):
        outputCanvas = ROOT.TCanvas("o_{dT}_{n}Jets".format(dT=distributionType, n=nJetsBin), "o_{dT}_{n}Jets".format(dT=distributionType, n=nJetsBin), 1440, 1200)
        ROOT.gStyle.SetOptStat(0)
        legend = ROOT.TLegend(0.6, 0.6, 0.9, 0.9, "")
        distributions_diphoton[distributionType][nJetsBin].SetLineColor(ROOT.kGreen+2)
        distributions_diphoton[distributionType][nJetsBin].SetTitle("{dT} distributions, {n} Jets".format(dT=distributionType, n=nJetsBin))
        distributions_diphoton[distributionType][nJetsBin].Draw("HIST E0 P")
        distributions_diphoton[distributionType][nJetsBin].GetXaxis().SetTitle(distributionType)
        distributions_diphoton[distributionType][nJetsBin].GetYaxis().SetRangeUser(0.00001, 10)
        sum_distributions = distributions_diphoton[distributionType][nJetsBin].Clone()
        sum_distributions.SetName("sum_{dt}_{n}JetsBin".format(dt=distributionType, n=nJetsBin))
        legendEntry = legend.AddEntry(distributions_diphoton[distributionType][nJetsBin], "Diphoton MC")
        legendEntry.SetTextColor(ROOT.kGreen+2)
        legendEntry.SetLineColor(ROOT.kGreen+2)
        distributions_GJet[distributionType][nJetsBin].SetLineColor(ROOT.kRed+1)
        distributions_GJet[distributionType][nJetsBin].Draw("HIST E0 P SAME")
        legendEntry = legend.AddEntry(distributions_GJet[distributionType][nJetsBin], "GJet MC")
        legendEntry.SetTextColor(ROOT.kRed+1)
        legendEntry.SetLineColor(ROOT.kRed+1)
        sum_distributions.Add(distributions_GJet[distributionType][nJetsBin])
        distributions_QCD[distributionType][nJetsBin].SetLineColor(ROOT.kBlue+1)
        distributions_QCD[distributionType][nJetsBin].Draw("HIST E0 P SAME")
        legendEntry = legend.AddEntry(distributions_QCD[distributionType][nJetsBin], "QCD MC")
        legendEntry.SetTextColor(ROOT.kBlue+1)
        legendEntry.SetLineColor(ROOT.kBlue+1)
        sum_distributions.Add(distributions_QCD[distributionType][nJetsBin])
        sum_distributions.SetLineColor(ROOT.kMagenta+2)
        sum_distributions.Draw("HIST SAME")
        legendEntry = legend.AddEntry(distributions_QCD[distributionType][nJetsBin], "Diphoton+GJet+QCD MC")
        legendEntry.SetTextColor(ROOT.kMagenta+2)
        legendEntry.SetLineColor(ROOT.kMagenta+2)
        if not(blinded): distributions_data[distributionType][nJetsBin].Draw("SAME")
        # titleText = ROOT.TText()
        # titleText.SetTextFont(42)
        # titleText.SetTextAlign(21)
        # titleText.DrawTextNDC(0.5, 0.95, )
        # outputStack.GetXaxis().SetTitle(distributionType)
        # outputStack.GetYaxis().SetTitle("Events/GeV")
        ROOT.gPad.SetLogy()
        ROOT.gPad.Update()
        # outputStack.GetYaxis().SetRangeUser(0.00001, 10.)
        ROOT.gPad.Update()
        legend.Draw()
        ROOT.gPad.Update()
        outputCanvas.SaveAs("{oD}/distributions_{dT}_{s}_{n}Jets.pdf".format(oD=outputDirectory, dT=distributionType, s=selection, n=nJetsBin))
