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

STNormMax = 1612.5

colors = {
    2: ROOT.kBlack,
    3: ROOT.kBlue+2,
    4: ROOT.kRed+1,
    5: ROOT.kGreen+3,
    6: ROOT.kViolet,
}

outputDirectory = "~/nobackup/analysisAreas/normBinOptimization_control"

# selection = "singlefake"
# prefix = "MC_GJet"
# sourceFile  = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_tightenedLooseSignal_singlePhotonTrigger_lowerSTThreshold/merged_selection_MC_GJet17_singlephoton_2017_control_{s}.root".format(s=selection)
# getMCWeights = True

selection = "control"
prefix = "data"
sourceFile  = stealthEnv.EOSPrefix + "/store/user/lpcsusystealth/selections/combined_DoublePhoton_tightenedLooseSignal_singlePhotonTrigger_lowerSTThreshold/merged_selection_data_2017_{s}.root".format(s=selection)
getMCWeights = False

if not(os.path.isdir(outputDirectory)): subprocess.check_call("mkdir -p {oD}".format(oD=outputDirectory), shell=True, executable="/bin/bash")

STRegionBoundariesFileObject = open("STRegionBoundaries_forNormOptimization_wider.dat")
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

outputFile = ROOT.TFile.Open("{oD}/{s}_{p}_outputs.root".format(oD=outputDirectory, p=prefix, s=selection), "RECREATE")
ROOT.RooAbsData.setDefaultStorageType(ROOT.RooAbsData.Tree)

rooSTVar = ROOT.RooRealVar("rooVar_sT", "rooVar_sT", STMin, STMax, "GeV")
rooWeightVar = ROOT.RooRealVar("rooVar_weight", "rooVar_weight", 0., 100000., "")
dataSets = {}
STDistributions = {}
for nJetsBin in range(2, 7):
    STDistributions[nJetsBin] = ROOT.TH1F("h_ST_{n}JetsBin".format(n=nJetsBin), "ST distribution: {n} Jets;ST".format(n=nJetsBin), n_STBins, array.array('d', STBoundaries))

for nJetsBin in range(2, 7):
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
inputFile.Close()

# kernels = {}
# shape_cumulatives = {}
# rooSTVar.setBins(nSTBins)
# STFrame = rooSTVar.frame(ROOT.RooFit.Bins(nSTBins), ROOT.RooFit.Title("kernels, cumulative".format(n=nJetsBin)))

# outputCanvas = ROOT.TCanvas("c_cumulativeShapes", "c_cumulativeShapes", 1024, 768)
# for nJetsBin in range(2, 7):
#     dataHist = ROOT.RooDataHist("dataHist_{n}Jets".format(n=nJetsBin), "data, {n} Jets".format(n=nJetsBin), ROOT.RooArgSet(rooSTVar), dataSets[nJetsBin])
#     kernels[nJetsBin] = ROOT.RooKeysPdf("kernelEstimate_{n}Jets".format(n=nJetsBin), "kernelEstimate_{n}Jets".format(n=nJetsBin), rooSTVar, dataSets[nJetsBin], ROOT.RooKeysPdf.MirrorLeft, 1.8)
#     shape_cumulatives[nJetsBin] = kernels[nJetsBin].createCdf(ROOT.RooArgSet(rooSTVar))
#     shape_cumulatives[nJetsBin].plotOn(STFrame, ROOT.RooFit.LineColor(colors[nJetsBin]))
# STFrame.Draw()
# outputCanvas.SaveAs("{oD}/STCumulativeShapes_{p}_{s}.pdf".format(oD=outputDirectory, p=prefix, s=selection))

normBinIndex = 5
outputCanvas = ROOT.TCanvas("c_chisqPerNDF", "c_chisqPerNDF", 1024, 768)
chisqPerNDFGraph = ROOT.TGraph()
while True:
    STNorm = STDistributions[2].GetXaxis().GetBinCenter(normBinIndex)
    if (STNorm > STNormMax): break
    print ("Trying ST norm: {n}".format(n=STNorm))
    histogramCopy_2Jets = ROOT.TH1F("h_ST_2JetsBin_normBin{i}".format(i=normBinIndex), "ST distribution: {2} Jets;ST", n_STBins-normBinIndex, array.array('d', STBoundaries[normBinIndex:]))
    histogramCopy_2Jets.Sumw2()
    for binCounter in range(1, 1 + histogramCopy_2Jets.GetXaxis().GetNbins()):
        binCenter = histogramCopy_2Jets.GetXaxis().GetBinCenter(binCounter)
        original_binIndex = STDistributions[2].FindFixBin(binCenter)
        histogramCopy_2Jets.SetBinContent(binCounter, STDistributions[2].GetBinContent(original_binIndex))
        histogramCopy_2Jets.SetBinError(binCounter, STDistributions[2].GetBinError(original_binIndex))
    chisqPerNDF_avg = 0.
    for nJetsBin in range(3, 7):
        histogramCopy = ROOT.TH1F("h_ST_{n}JetsBin_normBin{i}".format(n=nJetsBin, i=normBinIndex), "ST distribution: {2} Jets;ST", n_STBins-normBinIndex, array.array('d', STBoundaries[normBinIndex:]))
        histogramCopy.Sumw2()
        for binCounter in range(1, 1 + histogramCopy.GetXaxis().GetNbins()):
            binCenter = histogramCopy.GetXaxis().GetBinCenter(binCounter)
            original_binIndex = STDistributions[nJetsBin].FindFixBin(binCenter)
            histogramCopy.SetBinContent(binCounter, STDistributions[nJetsBin].GetBinContent(original_binIndex))
            histogramCopy.SetBinError(binCounter, STDistributions[nJetsBin].GetBinError(original_binIndex))
        # comparisonType = "UU"
        # if getMCWeights: comparisonType = "WW"
        comparisonType = "WW"
        chi2 = histogramCopy.Chi2Test(histogramCopy_2Jets, "{cT} P CHI2/NDF".format(cT=comparisonType))
        chisqPerNDF_avg += chi2/4 # four nJets bins
    chisqPerNDFGraph.SetPoint(chisqPerNDFGraph.GetN(), STNorm, chisqPerNDF_avg)
    normBinIndex += 1
chisqPerNDFGraph.SetMarkerStyle(ROOT.kCircle)
chisqPerNDFGraph.SetMarkerSize(0.5)
chisqPerNDFGraph.Draw("AP")
chisqPerNDFGraph.GetXaxis().SetTitle("ST norm")
chisqPerNDFGraph.GetYaxis().SetTitle("#chi^{2}/ndf")
outputCanvas.Update()
outputCanvas.SaveAs("{oD}/chiSqPerNDF_{p}_{s}.pdf".format(oD=outputDirectory, p=prefix, s=selection))

STNormTarget = 1300.0
outputCanvas = ROOT.TCanvas("c_STDistributions", "c_STDistributions", 1024, 768)
legend = ROOT.TLegend(0.7, 0.6, 0.9, 0.9)
STDistributions[2].SetLineColor(colors[2])
STDistributions[2].Draw()
legendEntry = legend.AddEntry(STDistributions[2], "nJets = {nJetsBin}".format(nJetsBin=2))
legendEntry.SetMarkerColor(colors[2])
legendEntry.SetLineColor(colors[2])
legendEntry.SetTextColor(colors[2])
for nJetsBin in range(3, 7):
    STDistributions[nJetsBin].Scale(STDistributions[2].GetBinContent(STDistributions[2].FindFixBin(STNormTarget))/STDistributions[nJetsBin].GetBinContent(STDistributions[nJetsBin].FindFixBin(STNormTarget)))
    STDistributions[nJetsBin].SetLineColor(colors[nJetsBin])
    STDistributions[nJetsBin].Draw("SAME")
    legendText = "nJets = {nJetsBin}".format(nJetsBin=nJetsBin)
    if (nJetsBin == 6): legendText = "nJets #geq 6"
    legendEntry = legend.AddEntry(STDistributions[nJetsBin], legendText)
    legendEntry.SetMarkerColor(colors[nJetsBin])
    legendEntry.SetLineColor(colors[nJetsBin])
    legendEntry.SetTextColor(colors[nJetsBin])
legend.Draw()
ROOT.gPad.SetLogy()
outputCanvas.Update()
outputCanvas.SaveAs("{oD}/STDistributions_{p}_{s}.pdf".format(oD=outputDirectory, p=prefix, s=selection))
outputFile.Close()
