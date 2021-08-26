#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, os, sys, argparse, stealthEnv, array, math, pdb, subprocess
inputArgumentsParser = argparse.ArgumentParser(description='Extract a few statistics histograms and save them to output image.')
inputArgumentsParser.add_argument('--inputFilePath', default="{eP}/{sER}/statistics/combined_DoublePhoton/merged_statistics_MC_stealth_t5_2017.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot), help='Path to file containing merged statistics.',type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--outputFolder', default="{aR}/IDEfficiencies", help='Path to folder in which to store output plots.',type=str)
# inputArgumentsParser.add_argument('--restrictToMCBulk', action='store_true', help="ID efficiency plots are sourced as the \"bulk\" MC ones.")
inputArguments = inputArgumentsParser.parse_args()

# selectionsList = ["signal", "control_fakefake"]
# colorsDict = {
#     "signal": ROOT.kBlack,
#     "signal_loose": ROOT.kBlue,
#     "control_fakefake": ROOT.kRed
# }

# MCRegionsList = ["bulk_closeToContours"]

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

if not(os.path.isdir(inputArguments.outputFolder)): subprocess.check_call("mkdir -p {oD}".format(oD=inputArguments.outputFolder), shell=True, executable="/bin/bash")

STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries)
STBoundaries = []
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()):
        STBoundary = float(STBoundaryString.strip())
        STBoundaries.append(STBoundary)
STBoundaries.append(3500)
nSTSignalBins = len(STBoundaries) - 2 # First two lines are lower and upper boundaries for the normalization bin
n_STBins = len(STBoundaries) - 1
STRegionsAxis = ROOT.TAxis(len(STBoundaries)-1, array.array('d', STBoundaries))

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
            outputCanvas.SaveAs(outputFolder + "/" + (prefix+selection+suffix) + ".pdf")
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
    # outputCanvas.SaveAs(outputFolder + "/" + (prefix+suffix).replace("__", "_") + ".pdf")

# for MCRegion in MCRegionsList:
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="deltaR_nearestGenJet_", suffix="_truePhotons_"+MCRegion, additionalFormatting="")    
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="deltaR_nearestSingletMomGenJet_", suffix="_truePhotons_"+MCRegion, additionalFormatting="")
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="deltaR_nearestEventProgenitorMomGenJet_", suffix="_truePhotons_"+MCRegion, additionalFormatting="")
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="partonMomID_", suffix="_all_genJets_"+MCRegion, additionalFormatting="")
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="partonID_", suffix="_all_genJets_"+MCRegion, additionalFormatting="")
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="MC_nEventProgenitorMomGenJets_", suffix="_selectedEvents_"+MCRegion, additionalFormatting="")
#     saveHistograms(outputFolder=inputArguments.outputFolder, prefix="MC_nGenJets_", suffix="_selectedEvents_"+MCRegion, additionalFormatting="")

def plotAndSaveIDEfficiencies(efficiencies, outputFolder):
    print("Plotting efficiencies...")
    clonesForRatios = {# "control": {},
                       "signal": {},
                       "signal_loose": {}}
    for nJetsBin in range(2, 7):
        h_efficiency_signal = efficiencies["signal"][nJetsBin]
        # h_efficiencyClone_signal = h_efficiency_signal.Clone("IDEfficiencyClone_{nJB}Jets_signal".format(nJB=nJetsBin))
        # clonesForRatios["signal"][nJetsBin] = h_efficiency_signal.Clone("IDEfficiencyCloneForRatio_{nJB}Jets_signal".format(nJB=nJetsBin))
        h_efficiency_signal_loose = efficiencies["signal_loose"][nJetsBin]
        # h_efficiencyClone_signal_loose = h_efficiency_signal_loose.Clone("IDEfficiencyClone_{nJB}Jets_signal_loose".format(nJB=nJetsBin))
        # clonesForRatios["signal_loose"][nJetsBin] = h_efficiency_signal_loose.Clone("IDEfficiencyCloneForRatio_{nJB}Jets_signal_loose".format(nJB=nJetsBin))
        # h_efficiency_control = efficiencies["control"][nJetsBin]
        # h_efficiencyClone_control = h_efficiency_control.Clone("IDEfficiencyClone_{nJB}Jets_control_fakefake".format(nJB=nJetsBin))
        # clonesForRatios["control"][nJetsBin] = h_efficiency_control.Clone("IDEfficiencyCloneForRatio_{nJB}Jets_control".format(nJB=nJetsBin))

        outputCanvas = ROOT.TCanvas("oC_{nJB}".format(nJB=nJetsBin), "oC_{nJB}".format(nJB=nJetsBin), 1024, 768)
        ROOT.gStyle.SetOptStat(0)
        # outputCanvas.Divide(1, 2)
        # outputCanvas.cd(1)
        legend = ROOT.TLegend(0.4, 0.85, 0.9, 0.9)
        # legend.SetNColumns(3)
        legend.SetNColumns(2)

        h_efficiency_signal.SetLineColor(ROOT.kBlack)
        h_efficiency_signal.Draw()
        ROOT.gPad.Update()
        signal_graphObject = h_efficiency_signal.GetPaintedGraph()
        signal_graphObject.GetXaxis().SetRangeUser(STBoundaries[0], STBoundaries[-1])
        signal_graphObject.GetXaxis().SetTitle("")
        signal_graphObject.SetMinimum(0.)
        signal_graphObject.SetMaximum(0.5)
        ROOT.gPad.Update()
        legendEntry = legend.AddEntry(h_efficiency_signal, "signal")
        legendEntry.SetLineColor(ROOT.kBlack)
        legendEntry.SetTextColor(ROOT.kBlack)

        h_efficiency_signal_loose.SetLineColor(ROOT.kBlue)
        h_efficiency_signal_loose.Draw("SAME")
        ROOT.gPad.Update()
        signal_loose_graphObject = h_efficiency_signal_loose.GetPaintedGraph()
        signal_loose_graphObject.SetMinimum(0.)
        # signal_loose_graphObject.SetMaximum(0.21)
        ROOT.gPad.Update()
        legendEntry = legend.AddEntry(h_efficiency_signal_loose, "signal_loose")
        legendEntry.SetLineColor(ROOT.kBlue)
        legendEntry.SetTextColor(ROOT.kBlue)

        # h_efficiency_control.SetLineColor(ROOT.kRed)
        # h_efficiency_control.Draw("SAME")
        # ROOT.gPad.Update()
        # control_graphObject = h_efficiency_control.GetPaintedGraph()
        # control_graphObject.SetMinimum(0.)
        # # control_graphObject.SetMaximum(0.21)
        # ROOT.gPad.Update()
        # legendEntry = legend.AddEntry(h_efficiency_control, "control")
        # legendEntry.SetLineColor(ROOT.kRed)
        # legendEntry.SetTextColor(ROOT.kRed)

        legend.Draw("SAME")
        ROOT.gPad.Update()

        # outputCanvas.cd(2)

        # legendRatio = ROOT.TLegend(0.4, 0.85, 0.9, 0.9)
        # legendRatio.SetNColumns(2)

        # h_efficiencyRatio_signalLooseOverSignal = ROOT.TH1F("efficiencyRatio_signalLooseOverSignal_{nJB}".format(nJB=nJetsBin), "", n_STBins, array.array('d', STBoundaries))
        # for binIndex in range(1, 1 + h_efficiencyRatio_signalLooseOverSignal.GetXaxis().GetNbins()):
        #     binCenter = h_efficiencyRatio_signalLooseOverSignal.GetXaxis().GetBinCenter(binIndex)
        #     signalEfficiency = h_efficiencyClone_signal.GetEfficiency(binIndex)
        #     signalLooseEfficiency = h_efficiencyClone_signal_loose.GetEfficiency(binIndex)
        #     if ((signalEfficiency > 0.) and (signalLooseEfficiency > 0.)):
        #         signalEfficiencyFractionalError = 0.5*(h_efficiencyClone_signal.GetEfficiencyErrorLow(binIndex) + h_efficiencyClone_signal.GetEfficiencyErrorUp(binIndex))/signalEfficiency
        #         signalLooseEfficiencyFractionalError = 0.5*(h_efficiencyClone_signal_loose.GetEfficiencyErrorLow(binIndex) + h_efficiencyClone_signal_loose.GetEfficiencyErrorUp(binIndex))/signalLooseEfficiency
        #         h_efficiencyRatio_signalLooseOverSignal.SetBinContent(binIndex, (signalLooseEfficiency/signalEfficiency))
        #         h_efficiencyRatio_signalLooseOverSignal.SetBinError(binIndex, (signalLooseEfficiency/signalEfficiency)*math.sqrt(pow(signalEfficiencyFractionalError, 2) + pow(signalLooseEfficiencyFractionalError, 2)))
        # h_efficiencyRatio_signalLooseOverSignal.SetLineColor(ROOT.kBlue)
        # h_efficiencyRatio_signalLooseOverSignal.Draw()
        # h_efficiencyRatio_signalLooseOverSignal.GetXaxis().SetTitle("ST")
        # h_efficiencyRatio_signalLooseOverSignal.GetYaxis().SetRangeUser(-1., 15.)
        # legendRatioEntry = legendRatio.AddEntry(h_efficiencyRatio_signalLooseOverSignal, "loose signal/signal")
        # legendRatioEntry.SetLineColor(ROOT.kBlue)
        # legendRatioEntry.SetTextColor(ROOT.kBlue)
        # ROOT.gPad.Update()

        # h_efficiencyRatio_controlOverSignal = ROOT.TH1F("efficiencyRatio_controlOverSignal_{nJB}".format(nJB=nJetsBin), "", n_STBins, array.array('d', STBoundaries))
        # for binIndex in range(1, 1 + h_efficiencyRatio_controlOverSignal.GetXaxis().GetNbins()):
        #     binCenter = h_efficiencyRatio_controlOverSignal.GetXaxis().GetBinCenter(binIndex)
        #     signalEfficiency = h_efficiencyClone_signal.GetEfficiency(binIndex)
        #     controlEfficiency = h_efficiencyClone_control.GetEfficiency(binIndex)
        #     if ((signalEfficiency > 0.) and (controlEfficiency > 0.)):
        #         signalEfficiencyFractionalError = 0.5*(h_efficiencyClone_signal.GetEfficiencyErrorLow(binIndex) + h_efficiencyClone_signal.GetEfficiencyErrorUp(binIndex))/signalEfficiency
        #         controlEfficiencyFractionalError = 0.5*(h_efficiencyClone_control.GetEfficiencyErrorLow(binIndex) + h_efficiencyClone_control.GetEfficiencyErrorUp(binIndex))/controlEfficiency
        #         h_efficiencyRatio_controlOverSignal.SetBinContent(binIndex, (controlEfficiency/signalEfficiency))
        #         h_efficiencyRatio_controlOverSignal.SetBinError(binIndex, (controlEfficiency/signalEfficiency)*math.sqrt(pow(signalEfficiencyFractionalError, 2) + pow(controlEfficiencyFractionalError, 2)))
        # h_efficiencyRatio_controlOverSignal.SetLineColor(ROOT.kRed)
        # h_efficiencyRatio_controlOverSignal.Draw("SAME")
        # legendRatioEntry = legendRatio.AddEntry(h_efficiencyRatio_controlOverSignal, "control/signal")
        # legendRatioEntry.SetLineColor(ROOT.kRed)
        # legendRatioEntry.SetTextColor(ROOT.kRed)
        # ROOT.gPad.Update()

        # linePlotter = ROOT.TLine()
        # linePlotter.SetLineColor(ROOT.kBlack)
        # linePlotter.SetLineStyle(ROOT.kDashed)
        # linePlotter.DrawLine(h_efficiencyRatio_signalLooseOverSignal.GetXaxis().GetXmin(), 1., h_efficiencyRatio_signalLooseOverSignal.GetXaxis().GetXmax(), 1.)
        # ROOT.gPad.Update()

        # legendRatio.Draw("SAME")
        # ROOT.gPad.Update()

        outputCanvas.SaveAs("{oF}/efficiencies_{nJB}Jets.pdf".format(oF=outputFolder, nJB=nJetsBin))

    # print("Now plotting ratios of efficiency shapes between nJets bins:")

    # # First fill a dictionary with the shape values
    # efficiencyShapes = {"control": {},
    #                     "signal": {},
    #                     "signal_loose": {}}
    # for nJetsBin in range(2, 7):
    #     for selection in ["signal", "signal_loose", "control"]:
    #         STNorm = 0.5*(STRegionsAxis.GetBinLowEdge(1) + STRegionsAxis.GetBinUpEdge(1))
    #         efficiency_norm = clonesForRatios[selection][nJetsBin].GetEfficiency(clonesForRatios[selection][nJetsBin].FindFixBin(STNorm))
    #         efficiencyError_norm = 0.5*(clonesForRatios[selection][nJetsBin].GetEfficiencyErrorLow(clonesForRatios[selection][nJetsBin].FindFixBin(STNorm)) + clonesForRatios[selection][nJetsBin].GetEfficiencyErrorUp(clonesForRatios[selection][nJetsBin].FindFixBin(STNorm)))
    #         efficiencyShapes[selection][nJetsBin] = {}
    #         for STBinIndex in range(2, 1 + STRegionsAxis.GetNbins()):
    #             STBinMidpoint = 0.5*(STRegionsAxis.GetBinLowEdge(STBinIndex) + STRegionsAxis.GetBinUpEdge(STBinIndex))
    #             efficiency_STBin = clonesForRatios[selection][nJetsBin].GetEfficiency(clonesForRatios[selection][nJetsBin].FindFixBin(STBinMidpoint))
    #             efficiencyError_STBin = 0.5*(clonesForRatios[selection][nJetsBin].GetEfficiencyErrorLow(clonesForRatios[selection][nJetsBin].FindFixBin(STBinMidpoint)) + clonesForRatios[selection][nJetsBin].GetEfficiencyErrorUp(clonesForRatios[selection][nJetsBin].FindFixBin(STBinMidpoint)))
    #             efficiencyShapes[selection][nJetsBin][STBinIndex] = {"value": efficiency_STBin/efficiency_norm,
    #                                                                  "error": efficiency_STBin/efficiency_norm*math.sqrt(pow(efficiencyError_STBin/efficiency_STBin, 2) + pow(efficiencyError_norm/efficiency_norm, 2))}

    # # Next plot the ratios
    # for nJetsBin in range(3, 7):
    #     for selection in ["signal", "signal_loose", "control"]:            
    #         outputCanvas = ROOT.TCanvas("oCRatio_{nJB}_{s}".format(nJB=nJetsBin, s=selection), "oCRatio_{nJB}_{s}".format(nJB=nJetsBin, s=selection), 1024, 768)
    #         outputHistogram = ROOT.TH1F("efficiencyShapeRatio_{s}_{nJB}Jets".format(s=selection, nJB=nJetsBin), "{s} selection: efficiency shape ratio, {nJB} Jets/2 Jets;ST;ratio".format(s=selection, nJB=nJetsBin), len(STBoundaries)-1, array.array('d', STBoundaries))
    #         outputHistogram.SetBinContent(1, 1.)
    #         outputHistogram.SetBinError(1, 0.)
    #         for STBinIndex in range(2, 1 + STRegionsAxis.GetNbins()):
    #             # pdb.set_trace()
    #             outputHistogram.SetBinContent(STBinIndex, efficiencyShapes[selection][nJetsBin][STBinIndex]["value"]/efficiencyShapes[selection][2][STBinIndex]["value"])
    #             outputHistogram.SetBinError(STBinIndex, (efficiencyShapes[selection][nJetsBin][STBinIndex]["value"]/efficiencyShapes[selection][2][STBinIndex]["value"])*math.sqrt(pow(efficiencyShapes[selection][nJetsBin][STBinIndex]["error"]/efficiencyShapes[selection][nJetsBin][STBinIndex]["value"], 2) + pow(efficiencyShapes[selection][2][STBinIndex]["error"]/efficiencyShapes[selection][2][STBinIndex]["value"], 2)))
    #         outputHistogram.Draw()
    #         linePlotter = ROOT.TLine()
    #         linePlotter.SetLineColor(ROOT.kBlack)
    #         linePlotter.SetLineStyle(ROOT.kDashed)
    #         linePlotter.DrawLine(STRegionsAxis.GetXmin(), 1., STRegionsAxis.GetXmax(), 1.)
    #         ROOT.gPad.Update()
    #         outputCanvas.SaveAs("{oF}/efficiencyShapeRatio_{s}_{nJB}Jets.pdf".format(oF=outputFolder, s=selection, nJB=nJetsBin))

def plotAndSave2DMiscHistograms(outputFilePath_chIso_neutIso, outputFilePath_chIso_phoIso, outputFilePath_neutIso_phoIso):
    outputCanvas_chIso_neutIso = ROOT.TCanvas("oC1", "oC1", 1024, 768)
    h_chIso_neutIso = ROOT.TH2F()
    h_chIso_neutIso.SetName("chIso_neutIso")
    inputFile.GetObject("chIso_neutIso_fake", h_chIso_neutIso)
    h_chIso_neutIso.SetName("chIso_neutIso")
    h_chIso_neutIso.Draw("colz")
    ROOT.gPad.SetLogz()
    h_chIso_neutIso.GetXaxis().SetRangeUser(0., 10.)
    h_chIso_neutIso.GetYaxis().SetRangeUser(0., 10.)
    outputCanvas_chIso_neutIso.Update()
    outputCanvas_chIso_neutIso.SaveAs(outputFilePath_chIso_neutIso)

    outputCanvas_chIso_phoIso = ROOT.TCanvas("oC2", "oC2", 1024, 768)
    h_chIso_phoIso = ROOT.TH2F()
    h_chIso_phoIso.SetName("chIso_phoIso")
    inputFile.GetObject("chIso_phoIso_fake", h_chIso_phoIso)
    h_chIso_phoIso.SetName("chIso_phoIso")
    h_chIso_phoIso.Draw("colz")
    ROOT.gPad.SetLogz()
    h_chIso_phoIso.GetXaxis().SetRangeUser(0., 10.)
    h_chIso_phoIso.GetYaxis().SetRangeUser(0., 10.)
    outputCanvas_chIso_phoIso.Update()
    outputCanvas_chIso_phoIso.SaveAs(outputFilePath_chIso_phoIso)

    outputCanvas_neutIso_phoIso = ROOT.TCanvas("oC3", "oC3", 1024, 768)
    h_neutIso_phoIso = ROOT.TH2F()
    h_neutIso_phoIso.SetName("neutIso_phoIso")
    inputFile.GetObject("neutIso_phoIso_fake", h_neutIso_phoIso)
    h_neutIso_phoIso.SetName("neutIso_phoIso")
    h_neutIso_phoIso.Draw("colz")
    ROOT.gPad.SetLogz()
    h_neutIso_phoIso.GetXaxis().SetRangeUser(0., 10.)
    h_neutIso_phoIso.GetYaxis().SetRangeUser(0., 10.)
    outputCanvas_neutIso_phoIso.Update()
    outputCanvas_neutIso_phoIso.SaveAs(outputFilePath_neutIso_phoIso)

def saveJetFakeProbabilities(name_efficiencyToExtract, outputFilePath, rangeMax):
    outputCanvas = ROOT.TCanvas("c_{n}".format(n=name_efficiencyToExtract), "c_{n}".format(n=name_efficiencyToExtract), 1024, 768)
    efficiencyToExtract = ROOT.TEfficiency()
    efficiencyToExtract.SetName(name_efficiencyToExtract)
    inputFile.GetObject(name_efficiencyToExtract, efficiencyToExtract)
    efficiencyToExtract.Draw()
    outputCanvas.Update()
    efficiencyToExtract.GetPaintedGraph().SetMaximum(rangeMax)
    outputCanvas.Update()
    outputCanvas.SaveAs(outputFilePath)

def saveJetFakeProbabilityRatios(name_efficiencyToExtract_numerator, name_efficiencyToExtract_denominator, outputFilePath, rangeMax=-1):
    outputCanvas = ROOT.TCanvas("c_ratio_{n}To{d}".format(n=name_efficiencyToExtract_numerator, d=name_efficiencyToExtract_denominator), "c_ratio_{n}To{d}".format(n=name_efficiencyToExtract_numerator, d=name_efficiencyToExtract_denominator), 1024, 768)

    efficiencyToExtract_numerator = ROOT.TEfficiency()
    efficiencyToExtract_numerator.SetName(name_efficiencyToExtract_numerator)
    inputFile.GetObject(name_efficiencyToExtract_numerator, efficiencyToExtract_numerator)

    efficiencyToExtract_denominator = ROOT.TEfficiency()
    efficiencyToExtract_denominator.SetName(name_efficiencyToExtract_denominator)
    inputFile.GetObject(name_efficiencyToExtract_denominator, efficiencyToExtract_denominator)

    efficiencyToExtract_denominator_totalHistCopy = efficiencyToExtract_denominator.GetCopyTotalHisto()
    ratioAxis = efficiencyToExtract_denominator_totalHistCopy.GetXaxis()
    ratioGraph = ROOT.TGraphErrors()
    for binIndex in range(1, 1 + ratioAxis.GetNbins()):
        numerator = efficiencyToExtract_numerator.GetEfficiency(binIndex)
        denominator = efficiencyToExtract_denominator.GetEfficiency(binIndex)
        if ((numerator > 0) and (denominator > 0)):
            ratio = numerator/denominator
            numeratorError = 0.5*(efficiencyToExtract_numerator.GetEfficiencyErrorLow(binIndex) + efficiencyToExtract_numerator.GetEfficiencyErrorUp(binIndex))
            denominatorError = 0.5*(efficiencyToExtract_denominator.GetEfficiencyErrorLow(binIndex) + efficiencyToExtract_denominator.GetEfficiencyErrorUp(binIndex))
            ratioError = ratio*math.sqrt(pow(numeratorError/numerator, 2) + pow(denominatorError/denominator, 2))
            ratioGraphBinIndex = ratioGraph.GetN()
            ratioGraph.SetPoint(ratioGraphBinIndex, ratioAxis.GetBinCenter(binIndex), ratio)
            ratioGraph.SetPointError(ratioGraphBinIndex, 0.5*(ratioAxis.GetBinUpEdge(binIndex) - ratioAxis.GetBinLowEdge(binIndex)), ratioError)
    ratioGraph.Draw()
    ratioGraph.GetXaxis().SetTitle(ratioAxis.GetTitle())
    ratioGraph.GetYaxis().SetTitle("ratio")
    outputCanvas.SaveAs(outputFilePath)

inputIDEfficiencies = {# "control": {},
                       "signal": {},
                       "signal_loose": {}}

for nJetsBin in range(2, 7):
    signalSource = "IDEfficiency_{nJB}Jets_signal".format(nJB=nJetsBin)
    # if (inputArguments.restrictToMCBulk): signalSource += "_MC_bulk_closeToContours"
    signalSource += "_MC_bulk_closeToContours"
    inputIDEfficiencies["signal"][nJetsBin] = ROOT.TEfficiency()
    inputIDEfficiencies["signal"][nJetsBin].SetName(signalSource)
    inputFile.GetObject(signalSource, inputIDEfficiencies["signal"][nJetsBin])
    if not(inputIDEfficiencies["signal"][nJetsBin]): sys.exit("ERROR: Unable to open histogram with name: {hN}".format(hN=signalSource))
    inputIDEfficiencies["signal"][nJetsBin].SetName(signalSource)

    signalLooseSource = "IDEfficiency_{nJB}Jets_signal_loose".format(nJB=nJetsBin)
    # if (inputArguments.restrictToMCBulk): signalLooseSource += "_MC_bulk_closeToContours"
    signalLooseSource += "_MC_bulk_closeToContours"
    inputIDEfficiencies["signal_loose"][nJetsBin] = ROOT.TEfficiency()
    inputIDEfficiencies["signal_loose"][nJetsBin].SetName(signalLooseSource)
    inputFile.GetObject(signalLooseSource, inputIDEfficiencies["signal_loose"][nJetsBin])
    if not(inputIDEfficiencies["signal_loose"][nJetsBin]): sys.exit("ERROR: Unable to open histogram with name: {hN}".format(hN=signalLooseSource))
    inputIDEfficiencies["signal_loose"][nJetsBin].SetName(signalLooseSource)

    # controlSource = "IDEfficiency_{nJB}Jets_control_fakefake".format(nJB=nJetsBin)
    # # if (inputArguments.restrictToMCBulk): controlSource += "_MC_bulk_closeToContours"
    # controlSource += "_MC_bulk_closeToContours"
    # inputIDEfficiencies["control"][nJetsBin] = ROOT.TEfficiency()
    # inputIDEfficiencies["control"][nJetsBin].SetName(controlSource)
    # inputFile.GetObject(controlSource, inputIDEfficiencies["control"][nJetsBin])
    # if not(inputIDEfficiencies["control"][nJetsBin]): sys.exit("ERROR: Unable to open histogram with name: {hN}".format(hN=controlSource))
    # inputIDEfficiencies["control"][nJetsBin].SetName(controlSource)

plotAndSaveIDEfficiencies(efficiencies=inputIDEfficiencies, outputFolder=inputArguments.outputFolder)

# plotAndSave2DMiscHistograms(outputFilePath_chIso_neutIso="{oF}/chIso_neutIso.pdf".format(oF=inputArguments.outputFolder), outputFilePath_chIso_phoIso="{oF}/chIso_phoIso.pdf".format(oF=inputArguments.outputFolder), outputFilePath_neutIso_phoIso="{oF}/neutIso_phoIso.pdf".format(oF=inputArguments.outputFolder))

# saveJetFakeProbabilities(name_efficiencyToExtract="recoEfficiency_fake_cleanJet", outputFilePath="{oF}/recoEfficiency_fake_cleanJet.pdf".format(oF=inputArguments.outputFolder), rangeMax=0.025)
# saveJetFakeProbabilities(name_efficiencyToExtract="recoEfficiency_fake_jetOverlappingPhoton", outputFilePath="{oF}/recoEfficiency_fake_jetOverlappingPhoton.pdf".format(oF=inputArguments.outputFolder), rangeMax=0.003)
# saveJetFakeProbabilities(name_efficiencyToExtract="recoEfficiency_loose_cleanJet", outputFilePath="{oF}/recoEfficiency_loose_cleanJet.pdf".format(oF=inputArguments.outputFolder), rangeMax=0.004)
# saveJetFakeProbabilities(name_efficiencyToExtract="recoEfficiency_loose_jetOverlappingPhoton", outputFilePath="{oF}/recoEfficiency_loose_jetOverlappingPhoton.pdf".format(oF=inputArguments.outputFolder), rangeMax=0.0002)
# saveJetFakeProbabilities(name_efficiencyToExtract="recoEfficiency_nominal_cleanJet", outputFilePath="{oF}/recoEfficiency_nominal_cleanJet.pdf".format(oF=inputArguments.outputFolder), rangeMax=0.004)
# saveJetFakeProbabilities(name_efficiencyToExtract="recoEfficiency_nominal_jetOverlappingPhoton", outputFilePath="{oF}/recoEfficiency_nominal_jetOverlappingPhoton.pdf".format(oF=inputArguments.outputFolder), rangeMax=0.0002)

# saveJetFakeProbabilityRatios("recoEfficiency_fake_cleanJet", "recoEfficiency_nominal_cleanJet", "{oF}/recoEfficiencyRatio_fakeToNominal_cleanJet.pdf".format(oF=inputArguments.outputFolder))
# saveJetFakeProbabilityRatios("recoEfficiency_loose_cleanJet", "recoEfficiency_nominal_cleanJet", "{oF}/recoEfficiencyRatio_looseToNominal_cleanJet.pdf".format(oF=inputArguments.outputFolder))
# saveJetFakeProbabilityRatios("recoEfficiency_fake_jetOverlappingPhoton", "recoEfficiency_nominal_jetOverlappingPhoton", "{oF}/recoEfficiencyRatio_fakeToNominal_jetOverlappingPhoton.pdf".format(oF=inputArguments.outputFolder))
# saveJetFakeProbabilityRatios("recoEfficiency_loose_jetOverlappingPhoton", "recoEfficiency_nominal_jetOverlappingPhoton", "{oF}/recoEfficiencyRatio_looseToNominal_jetOverlappingPhoton.pdf".format(oF=inputArguments.outputFolder))

inputFile.Close()
print("All done!")
