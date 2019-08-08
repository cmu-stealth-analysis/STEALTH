#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, os, sys, argparse
inputArgumentsParser = argparse.ArgumentParser(description='Extract a few statistics histograms and save them to output image.')
inputArgumentsParser.add_argument('--inputFilePath', default="mergedStatistics/statistics_MC_2017.root", help='Path to file containing merged statistics.',type=str)
inputArgumentsParser.add_argument('--outputFolder', default="mergedStatistics", help='Path to folder in which to store output plots.',type=str)
inputArguments = inputArgumentsParser.parse_args()

selectionsList = ["signal", "control_mediumfake", "control_fakefake"]
colorsDict = {
    "signal": ROOT.kBlack,
    "control_mediumfake": ROOT.kBlue,
    "control_fakefake": ROOT.kRed
}

MCRegionsList = ["MC_bulk", "MC_lowNeutralinoMass", "MC_gluinoNeutralinoDegenerate"]

inputFile = ROOT.TFile.Open(inputArguments.inputFilePath, "READ")
if ((inputFile.IsZombie() == ROOT.kTRUE) or not(inputFile.IsOpen() == ROOT.kTRUE)): sys.exit("ERROR: Unable to open file \"{f}\"".format(f=inputArguments.inputFilePath))

def saveHistograms(outputFolder, prefix, suffix, additionalFormatting):
    inputHistograms = {}
    runningMaxValue = -1.
    for selection in selectionsList:
        inputHistograms[selection] = ROOT.TH1F()
        inputFile.GetObject(prefix+selection+suffix, inputHistograms[selection])
        if (inputHistograms[selection]):
            outputCanvas = ROOT.TCanvas("oC", "oC", 1024, 768)
            inputHistograms[selection].Draw()
            ROOT.gPad.SetLogy()
            ROOT.gPad.Update()
            outputCanvas.SaveAs(outputFolder + "/" + (prefix+selection+suffix) + ".png")
        else:
            print("ERROR: histogram named \"{n}\" not found in file \"{f}\"".format(n=prefix+selection+suffix, f=inputArguments.inputFilePath))
    
    # outputCanvas = ROOT.TCanvas("oC", "oC", 1024, 768)
    # legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    # inputHistograms["signal"].Draw()
    # inputHistograms["signal"].SetLineColor(colorsDict["signal"])
    # legendEntry = legend.AddEntry(inputHistograms["signal"], "signal")
    # legendEntry.SetLineColor(colorsDict["signal"])
    # legendEntry.SetTextColor(colorsDict["signal"])
    # ROOT.gStyle.SetOptStat(0)
    # ROOT.gPad.SetLogy()
    # ROOT.gPad.Update()
    # # inputHistograms["control_mediumfake"].Draw("same")
    # # inputHistograms["control_mediumfake"].SetLineColor(colorsDict["control_mediumfake"])
    # # legendEntry = legend.AddEntry(inputHistograms["control_mediumfake"], "control_mediumfake")
    # # legendEntry.SetLineColor(colorsDict["control_mediumfake"])
    # # legendEntry.SetTextColor(colorsDict["control_mediumfake"])
    # # inputHistograms["control_fakefake"].Draw("same")
    # # inputHistograms["control_mediumfake"].SetLineColor(colorsDict["control_fakefake"])
    # # legendEntry = legend.AddEntry(inputHistograms["control_fakefake"], "control_fakefake")
    # # legendEntry.SetLineColor(colorsDict["control_fakefake"])
    # # legendEntry.SetTextColor(colorsDict["control_fakefake"])
    # # ROOT.gPad.Update()
    # legend.Draw("same")
    # ROOT.gPad.Update()
    # outputCanvas.SaveAs(outputFolder + "/" + (prefix+suffix).replace("__", "_") + ".png")

for MCRegion in MCRegionsList:
    saveHistograms(outputFolder=inputArguments.outputFolder, prefix="deltaR_nearestGenJet_", suffix="_truePhotons_"+MCRegion, additionalFormatting="")    
    saveHistograms(outputFolder=inputArguments.outputFolder, prefix="deltaR_nearestSingletMomGenJet_", suffix="_truePhotons_"+MCRegion, additionalFormatting="")
    saveHistograms(outputFolder=inputArguments.outputFolder, prefix="deltaR_nearestGluinoMomGenJet_", suffix="_truePhotons_"+MCRegion, additionalFormatting="")
    saveHistograms(outputFolder=inputArguments.outputFolder, prefix="partonMomID_", suffix="_all_genJets_"+MCRegion, additionalFormatting="")
    saveHistograms(outputFolder=inputArguments.outputFolder, prefix="partonID_", suffix="_all_genJets_"+MCRegion, additionalFormatting="")
    saveHistograms(outputFolder=inputArguments.outputFolder, prefix="MC_nGluinoMomGenJets_", suffix="_selectedEvents_"+MCRegion, additionalFormatting="")
    saveHistograms(outputFolder=inputArguments.outputFolder, prefix="MC_nGenJets_", suffix="_selectedEvents_"+MCRegion, additionalFormatting="")

inputFile.Close()
print("All done!")
