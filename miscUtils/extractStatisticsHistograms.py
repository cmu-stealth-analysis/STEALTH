#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, os, sys, argparse, stealthEnv, array, math
inputArgumentsParser = argparse.ArgumentParser(description='Extract a few statistics histograms and save them to output image.')
inputArgumentsParser.add_argument('--inputFilePath', default="{eP}/{sER}/statistics/combined_DoublePhoton/merged_statistics_MC_stealth_t5_2017.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot), help='Path to file containing merged statistics.',type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
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

STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries)
STBoundaries = []
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()):
        STBoundary = float(STBoundaryString.strip())
        STBoundaries.append(STBoundary)
STBoundaries.append(3500.0) # Instead of infinity
n_STBins = len(STBoundaries) - 1
STRegionsAxis = ROOT.TAxis(n_STBins, array.array('d', STBoundaries))

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
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="deltaR_nearestEventProgenitorMomGenJet_", suffix="_truePhotons_"+MCRegion, additionalFormatting="")
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="partonMomID_", suffix="_all_genJets_"+MCRegion, additionalFormatting="")
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="partonID_", suffix="_all_genJets_"+MCRegion, additionalFormatting="")
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="MC_nEventProgenitorMomGenJets_", suffix="_selectedEvents_"+MCRegion, additionalFormatting="")
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="MC_nGenJets_", suffix="_selectedEvents_"+MCRegion, additionalFormatting="")

for nJetsBin in range(2, 7):
    efficiency_signal_histogramName = "IDEfficiency_{nJB}Jets_signal_MC_bulk_closeToContours".format(nJB=nJetsBin)
    h_efficiency_signal = ROOT.TEfficiency()
    inputFile.GetObject(efficiency_signal_histogramName, h_efficiency_signal)
    if not(h_efficiency_signal): sys.exit("ERROR: Unable to open histogram with name: {hN}".format(hN=efficiency_signal_histogramName))
    h_efficiencyClone_signal = h_efficiency_signal.Clone("IDEfficiencyClone_{nJB}Jets_signal_MC_bulk_closeToContours".format(nJB=nJetsBin))
    efficiency_control_histogramName = "IDEfficiency_{nJB}Jets_control_fakefake_MC_bulk_closeToContours".format(nJB=nJetsBin)
    h_efficiency_control = ROOT.TEfficiency()
    inputFile.GetObject(efficiency_control_histogramName, h_efficiency_control)
    if not(h_efficiency_control): sys.exit("ERROR: Unable to open histogram with name: {hN}".format(hN=efficiency_control_histogramName))
    h_efficiencyClone_control = h_efficiency_control.Clone("IDEfficiencyClone_{nJB}Jets_control_fakefake_MC_bulk_closeToContours".format(nJB=nJetsBin))

    outputCanvas = ROOT.TCanvas("oC", "oC", 1024, 768)
    ROOT.gStyle.SetOptStat(0)
    outputCanvas.Divide(1, 2)
    outputCanvas.cd(1)
    legend = ROOT.TLegend(0.6, 0.85, 0.9, 0.9)
    legend.SetNColumns(2)

    h_efficiency_signal.SetLineColor(ROOT.kBlack)
    h_efficiency_signal.Draw()
    ROOT.gPad.Update()
    signal_graphObject = h_efficiency_signal.GetPaintedGraph()
    signal_graphObject.GetXaxis().SetRangeUser(STBoundaries[0], STBoundaries[-1])
    signal_graphObject.GetXaxis().SetTitle("")
    signal_graphObject.SetMinimum(0.)
    signal_graphObject.SetMaximum(0.3)
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

    outputCanvas.cd(2)
    h_efficiencyRatio = ROOT.TH1F("efficiencyRatio_{nJB}".format(nJB=nJetsBin), "", n_STBins, array.array('d', STBoundaries))
    for binIndex in range(1, 1 + h_efficiencyRatio.GetXaxis().GetNbins()):
        binCenter = h_efficiencyRatio.GetXaxis().GetBinCenter(binIndex)
        signalEfficiency = h_efficiencyClone_signal.GetEfficiency(binIndex)
        controlEfficiency = h_efficiencyClone_control.GetEfficiency(binIndex)
        if ((signalEfficiency > 0.) and (controlEfficiency > 0.)):
            signalEfficiencyFractionalError = 0.5*(h_efficiencyClone_signal.GetEfficiencyErrorLow(binIndex) + h_efficiencyClone_signal.GetEfficiencyErrorUp(binIndex))/signalEfficiency
            controlEfficiencyFractionalError = 0.5*(h_efficiencyClone_control.GetEfficiencyErrorLow(binIndex) + h_efficiencyClone_control.GetEfficiencyErrorUp(binIndex))/controlEfficiency
            h_efficiencyRatio.SetBinContent(binIndex, (signalEfficiency/controlEfficiency))
            h_efficiencyRatio.SetBinError(binIndex, (signalEfficiency/controlEfficiency)*math.sqrt(pow(signalEfficiencyFractionalError, 2) + pow(controlEfficiencyFractionalError, 2)))
    h_efficiencyRatio.SetLineColor(ROOT.kBlack)
    h_efficiencyRatio.Draw()
    h_efficiencyRatio.GetXaxis().SetTitle("ST")
    h_efficiencyRatio.GetYaxis().SetRangeUser(0., 5.)
    ROOT.gPad.Update()
    linePlotter = ROOT.TLine()
    linePlotter.SetLineColor(ROOT.kBlack)
    linePlotter.SetLineStyle(ROOT.kDashed)
    linePlotter.DrawLine(h_efficiencyRatio.GetXaxis().GetXmin(), 1., h_efficiencyRatio.GetXaxis().GetXmax(), 1.)
    ROOT.gPad.Update()
    outputCanvas.SaveAs("{oF}/efficiencies_{n}Jets.png".format(oF=inputArguments.outputFolder, n=nJetsBin))

inputFile.Close()
print("All done!")
