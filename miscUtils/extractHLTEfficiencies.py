#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, os, sys, argparse, stealthEnv
inputArgumentsParser = argparse.ArgumentParser(description='Extract HLT efficiencies histograms and save them to output image.')
inputArgumentsParser.add_argument('--outputFolder', default="~/nobackup/analysisAreas/HLTEfficiencies", help='Path to folder in which to store output files.',type=str)
inputArgumentsParser.add_argument('--outputPrefix', default="HLTEfficiencies", help='Path to folder in which to store output files.',type=str)
inputArguments = inputArgumentsParser.parse_args()

ROOT.gROOT.SetBatch(ROOT.kTRUE)

sources = {
    2016: {
        "clean": "{eP}/{sER}/statistics/combined_DoublePhoton_hgg/merged_statistics_MC_hgg_noJetSelection_2017.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot)
    },
    2017: {
        "clean": "{eP}/{sER}/statistics/combined_DoublePhoton_hgg/merged_statistics_MC_hgg_noJetSelection_2017.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot)
    },
    2018: {
        "clean": "{eP}/{sER}/statistics/combined_DoublePhoton_hgg/merged_statistics_MC_hgg_noJetSelection_2017.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot)
    }
}

targets = {
    "signal": "hltEfficiency1D_leadingPhoton_signal",
    "signal_loose": "hltEfficiency1D_leadingPhoton_signal_loose",
    "control": "hltEfficiency1D_leadingPhoton_control_fakefake"
}

for selection, efficiencyName in targets.items():
    print("Extracting selection: {s}, efficiencyName: {eN}".format(s=selection, eN=efficiencyName))
    for year, sourceTypePathDict in sources.items():
        print("Fetching efficiencies for year: {y}".format(y=year))
        # outputFile = ROOT.TFile.Open(inputArguments.outputFolder + "/" + inputArguments.outputPrefix + "_" + selection + ".root", "RECREATE")
        outputFileName = "{oF}/{oP}_{s}_{y}.root".format(oF=inputArguments.outputFolder, oP=inputArguments.outputPrefix, s=selection, y=year)
        outputFile = ROOT.TFile.Open(outputFileName, "RECREATE")
        if ((outputFile.IsZombie() == ROOT.kTRUE) or not(outputFile.IsOpen() == ROOT.kTRUE)): sys.exit("ERROR: Unable to open file \"{oFN}\"".format(oFN=outputFileName))
        for sourceType, path in sourceTypePathDict.items():
            print("Fetching efficiency of type: \"{sT}\" from path: {p}".format(sT=sourceType, p=path))
            inputFile = ROOT.TFile.Open(path, "READ")
            if ((inputFile.IsZombie() == ROOT.kTRUE) or not(inputFile.IsOpen() == ROOT.kTRUE)): sys.exit("ERROR: Unable to open file at path \"{p}\"".format(p=path))
            efficiencyToFetch = ROOT.TEfficiency()
            efficiency_label = (efficiencyName.replace("hltEfficiency1D_leadingPhoton_", "")).replace("_" + selection, "") + "_" + sourceType
            efficiencyToFetch.SetName(efficiency_label)
            inputFile.GetObject(efficiencyName, efficiencyToFetch)
            efficiencyToFetch.SetName("hltEfficiency_{s}".format(s=sourceType))
            c = ROOT.TCanvas("output_" + efficiency_label + "_" + efficiencyName + "_" + str(year), "output_" + efficiency_label + "_" + efficiencyName + "_" + str(year), 1024, 768)
            efficiencyToFetch.Draw()
            c.SaveAs("{oF}/{oP}_{l}_{y}.pdf".format(oF=inputArguments.outputFolder, oP=inputArguments.outputPrefix, l=efficiency_label, y=year))
            outputFile.WriteTObject(efficiencyToFetch)
            inputFile.Close()
        outputFile.Close()

print("All done!")
