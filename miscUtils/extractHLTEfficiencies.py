#!/usr/bin/env python

from __future__ import print_function, division

import subprocess, os, sys, argparse
import ROOT
import stealthEnv

inputArgumentsParser = argparse.ArgumentParser(description='Extract HLT efficiencies histograms and save them to output image.')
inputArgumentsParser.add_argument('--outputFolderPlots', default="~/nobackup/analysisAreas/HLTEfficiencies", help='Path to folder in which to store output plots.',type=str)
inputArgumentsParser.add_argument('--outputFolderEOSPath', default="{sER}/HLTEfficiencies".format(sER=stealthEnv.stealthEOSRoot), help='Path to folder in which to store output ROOT files.',type=str)
inputArgumentsParser.add_argument('--inputPrefix', default="", help='Input prefix used in paths to statistics files.',type=str)
inputArgumentsParser.add_argument('--outputPrefix', default="HLTEfficiencies", help='Prefix for output files.',type=str)
inputArguments = inputArgumentsParser.parse_args()

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

if not(os.path.isdir("{oF}".format(oF=inputArguments.outputFolderPlots))): subprocess.check_call("mkdir -p {oF}".format(oF=inputArguments.outputFolderPlots), shell=True, executable="/bin/bash")
if not(os.path.isdir("{sA}/HLTEfficiencies".format(sA=stealthEnv.scratchArea))): subprocess.check_call("mkdir -p {sA}/HLTEfficiencies".format(sA=stealthEnv.scratchArea), shell=True, executable="/bin/bash")
subprocess.check_call("eos {eP} mkdir -p {oFEP}".format(eP=stealthEnv.EOSPrefix, oFEP=inputArguments.outputFolderEOSPath), shell=True, executable="/bin/bash")

sources = {
    2016: {
        "clean": "{eP}/{sER}/statistics/combined_DoublePhoton_{iP}/merged_statistics_MC_hgg_noJetSelection_2016.root".format(eP=stealthEnv.EOSPrefix, iP=inputArguments.inputPrefix, sER=stealthEnv.stealthEOSRoot)
    },
    2017: {
        "clean": "{eP}/{sER}/statistics/combined_DoublePhoton_{iP}/merged_statistics_MC_hgg_noJetSelection_2017.root".format(eP=stealthEnv.EOSPrefix, iP=inputArguments.inputPrefix, sER=stealthEnv.stealthEOSRoot)
    },
    2018: {
        "clean": "{eP}/{sER}/statistics/combined_DoublePhoton_{iP}/merged_statistics_MC_hgg_noJetSelection_2018.root".format(eP=stealthEnv.EOSPrefix, iP=inputArguments.inputPrefix, sER=stealthEnv.stealthEOSRoot)
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
        outputFileName = "{oP}_{s}_{y}.root".format(oP=inputArguments.outputPrefix, s=selection, y=year)
        outputFilePath = "{sA}/HLTEfficiencies/{oFN}".format(sA=stealthEnv.scratchArea, oFN=outputFileName)
        outputFile = ROOT.TFile.Open(outputFilePath, "RECREATE")
        if ((outputFile.IsZombie() == ROOT.kTRUE) or not(outputFile.IsOpen() == ROOT.kTRUE)): sys.exit("ERROR: Unable to open file \"{oFP}\"".format(oFP=outputFilePath))
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
            c.SaveAs("{oF}/{oP}_{l}_{y}.pdf".format(oF=inputArguments.outputFolderPlots, oP=inputArguments.outputPrefix, l=efficiency_label, y=year))
            outputFile.WriteTObject(efficiencyToFetch)
            inputFile.Close()
        outputFile.Close()
        subprocess.check_call("xrdcp --nopbar --silent --force --path --streams 15 {oFP} {eP}/{oFEP}/{oFN} && rm {oFP}".format(oFP=outputFilePath, eP=stealthEnv.EOSPrefix, oFEP=inputArguments.outputFolderEOSPath, oFN=outputFileName), shell=True, executable="/bin/bash")
print("All done!")
