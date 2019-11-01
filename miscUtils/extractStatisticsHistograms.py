#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, os, sys, argparse, stealthEnv
inputArgumentsParser = argparse.ArgumentParser(description='Extract a few statistics histograms and save them to output image.')
inputArgumentsParser.add_argument('--inputFilePath', default="{eP}/{sER}/statistics/combined_DoublePhoton/merged_statistics_MC_stealth_t5_2017.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot), help='Path to file containing merged statistics.',type=str)
inputArgumentsParser.add_argument('--outputFolder', default="mergedStatistics", help='Path to folder in which to store output plots.',type=str)
inputArguments = inputArgumentsParser.parse_args()

# selectionsList = ["signal", "control_fakefake"]
# colorsDict = {
#     "signal": ROOT.kBlack,
#     "signal_loose": ROOT.kBlue,
#     "control_fakefake": ROOT.kRed
# }

# MCRegionsList = ["bulk_closeToContours"]

ROOT.gROOT.SetBatch(ROOT.kTRUE)
inputFile = ROOT.TFile.Open(inputArguments.inputFilePath, "READ")
if ((inputFile.IsZombie() == ROOT.kTRUE) or not(inputFile.IsOpen() == ROOT.kTRUE)): sys.exit("ERROR: Unable to open file \"{f}\"".format(f=inputArguments.inputFilePath))
print("Opened file: {iFP}".format(iFP=inputArguments.inputFilePath))

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
    # # inputHistograms["signal_loose"].Draw("same")
    # # inputHistograms["signal_loose"].SetLineColor(colorsDict["signal_loose"])
    # # legendEntry = legend.AddEntry(inputHistograms["signal_loose"], "signal_loose")
    # # legendEntry.SetLineColor(colorsDict["signal_loose"])
    # # legendEntry.SetTextColor(colorsDict["signal_loose"])
    # # inputHistograms["control_fakefake"].Draw("same")
    # # inputHistograms["signal_loose"].SetLineColor(colorsDict["control_fakefake"])
    # # legendEntry = legend.AddEntry(inputHistograms["control_fakefake"], "control_fakefake")
    # # legendEntry.SetLineColor(colorsDict["control_fakefake"])
    # # legendEntry.SetTextColor(colorsDict["control_fakefake"])
    # # ROOT.gPad.Update()
    # legend.Draw("same")
    # ROOT.gPad.Update()
    # outputCanvas.SaveAs(outputFolder + "/" + (prefix+suffix).replace("__", "_") + ".png")

# for MCRegion in MCRegionsList:
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="deltaR_nearestGenJet_", suffix="_truePhotons_"+MCRegion, additionalFormatting="")    
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="deltaR_nearestSingletMomGenJet_", suffix="_truePhotons_"+MCRegion, additionalFormatting="")
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="deltaR_nearestGluinoMomGenJet_", suffix="_truePhotons_"+MCRegion, additionalFormatting="")
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="partonMomID_", suffix="_all_genJets_"+MCRegion, additionalFormatting="")
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="partonID_", suffix="_all_genJets_"+MCRegion, additionalFormatting="")
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="MC_nGluinoMomGenJets_", suffix="_selectedEvents_"+MCRegion, additionalFormatting="")
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="MC_nGenJets_", suffix="_selectedEvents_"+MCRegion, additionalFormatting="")

for nJetsBin in range(2, 7):
    efficiency_signal_histogramName = "IDEfficiency_{nJB}Jets_signal_MC_bulk_closeToContours".format(nJB=nJetsBin)
    h_efficiency_signal = ROOT.TEfficiency()
    inputFile.GetObject(efficiency_signal_histogramName, h_efficiency_signal)
    if not(h_efficiency_signal): sys.exit("ERROR: Unable to open histogram with name: {hN}".format(hN=efficiency_signal_histogramName))
    efficiency_control_histogramName = "IDEfficiency_{nJB}Jets_control_fakefake_MC_bulk_closeToContours".format(nJB=nJetsBin)
    h_efficiency_control = ROOT.TEfficiency()
    inputFile.GetObject(efficiency_control_histogramName, h_efficiency_control)
    if not(h_efficiency_control): sys.exit("ERROR: Unable to open histogram with name: {hN}".format(hN=efficiency_control_histogramName))
    outputCanvas = ROOT.TCanvas("oC", "oC", 1024, 768)
    legend = ROOT.TLegend(0.6, 0.85, 0.9, 0.9)
    legend.SetNColumns(2)

    h_efficiency_signal.SetLineColor(ROOT.kBlack)
    h_efficiency_signal.Draw()
    ROOT.gPad.Update()
    signal_graphObject = h_efficiency_signal.GetPaintedGraph()
    signal_graphObject.SetMinimum(0.)
    signal_graphObject.SetMaximum(0.21)
    ROOT.gPad.Update()
    legendEntry = legend.AddEntry(h_efficiency_signal, "signal")
    legendEntry.SetLineColor(ROOT.kBlack)
    legendEntry.SetTextColor(ROOT.kBlack)

    h_efficiency_control.SetLineColor(ROOT.kRed)
    h_efficiency_control.Draw("SAME")
    ROOT.gPad.Update()
    control_graphObject = h_efficiency_control.GetPaintedGraph()
    control_graphObject.SetMinimum(0.)
    control_graphObject.SetMaximum(0.21)
    ROOT.gPad.Update()
    legendEntry = legend.AddEntry(h_efficiency_control, "control")
    legendEntry.SetLineColor(ROOT.kRed)
    legendEntry.SetTextColor(ROOT.kRed)

    legend.Draw("SAME")
    ROOT.gPad.Update()
    outputCanvas.SaveAs("{oF}/efficiencies_{n}Jets.png".format(oF=inputArguments.outputFolder, n=nJetsBin))

inputFile.Close()
print("All done!")
