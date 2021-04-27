#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, os, sys, argparse, stealthEnv, subprocess
inputArgumentsParser = argparse.ArgumentParser(description='Extract HLT efficiencies histograms and save them to output image.')
inputArgumentsParser.add_argument('--outputFolder', default="/uscms/home/tmudholk/nobackup/analysisAreas/PUWeights", help='Path to folder in which to store output files.',type=str)
inputArguments = inputArgumentsParser.parse_args()

ROOT.gROOT.SetBatch(ROOT.kTRUE)
if not(os.path.isdir(inputArguments.outputFolder)): subprocess.check_call("mkdir -p {oF}".format(oF=inputArguments.outputFolder), shell=True, executable="/bin/bash")

sourceFolder_MC = "{eP}/{sER}/analysisEOSAreas/analysis/".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot)
sourceFolder_data = "{sR}/getMCSystematics/data/".format(sR=stealthEnv.stealthRoot)

# targets = {
#     "signal": "hltEfficiency1D_leadingPhoton_signal",
#     "signal_loose": "hltEfficiency1D_leadingPhoton_signal_loose",
#     "control": "hltEfficiency1D_leadingPhoton_control_fakefake"
# }

# for selection, efficiencyName in targets.items():
#     print("Extracting selection: {s}, efficiencyName: {eN}".format(s=selection, eN=efficiencyName))
#     for year, sourceTypePathDict in sources.items():
#         print("Fetching efficiencies for year: {y}".format(y=year))
#         # outputFile = ROOT.TFile.Open(inputArguments.outputFolder + "/" + inputArguments.outputPrefix + "_" + selection + ".root", "RECREATE")
#         outputFileName = "{oF}/{oP}_{s}_{y}.root".format(oF=inputArguments.outputFolder, oP=inputArguments.outputPrefix, s=selection, y=year)
#         outputFile = ROOT.TFile.Open(outputFileName, "RECREATE")
#         if ((outputFile.IsZombie() == ROOT.kTRUE) or not(outputFile.IsOpen() == ROOT.kTRUE)): sys.exit("ERROR: Unable to open file \"{oFN}\"".format(oFN=outputFileName))
#         for sourceType, path in sourceTypePathDict.items():
#             print("Fetching efficiency of type: \"{sT}\" from path: {p}".format(sT=sourceType, p=path))
#             inputFile = ROOT.TFile.Open(path, "READ")
#             if ((inputFile.IsZombie() == ROOT.kTRUE) or not(inputFile.IsOpen() == ROOT.kTRUE)): sys.exit("ERROR: Unable to open file at path \"{p}\"".format(p=path))
#             efficiencyToFetch = ROOT.TEfficiency()
#             efficiency_label = (efficiencyName.replace("hltEfficiency1D_leadingPhoton_", "")).replace("_" + selection, "") + "_" + sourceType
#             efficiencyToFetch.SetName(efficiency_label)
#             inputFile.GetObject(efficiencyName, efficiencyToFetch)
#             efficiencyToFetch.SetName("hltEfficiency_{s}".format(s=sourceType))
#             c = ROOT.TCanvas("output_" + efficiency_label + "_" + efficiencyName + "_" + str(year), "output_" + efficiency_label + "_" + efficiencyName + "_" + str(year), 1024, 768)
#             efficiencyToFetch.Draw()
#             c.SaveAs("{oF}/{oP}_{l}_{y}.pdf".format(oF=inputArguments.outputFolder, oP=inputArguments.outputPrefix, l=efficiency_label, y=year))
#             outputFile.WriteTObject(efficiencyToFetch)
#             inputFile.Close()
#         outputFile.Close()

selection_names = {
    "signal": "signal",
    "signal_loose": "loose signal",
    "control": "diphoton control"
}

for year in ["2016", "2017", "2018"]:
    # First data histograms
    source_file_path = "{sF}/dataPU_{y}.root".format(sF=sourceFolder_data, y=year)
    inputFileObject = ROOT.TFile.Open(source_file_path, "READ")
    if ((inputFileObject.IsZombie() == ROOT.kTRUE) or not(inputFileObject.IsOpen() == ROOT.kTRUE)):
        sys.exit("ERROR: Unable to open file {f}".format(f=source_file_path))
    puDataHistogram = ROOT.TH1D()
    inputFileObject.GetObject("pileup", puDataHistogram)
    puDataHistogram.GetXaxis().SetTitle("PU")
    puDataHistogram.GetYaxis().SetTitle("Events/100")
    puDataHistogram.SetTitle(("Pileup distribution, data, {y}").format(y=year))
    outputObjectName = "pileup_data_{y}".format(y=year)
    c = ROOT.TCanvas("output_" + outputObjectName, "output_" + outputObjectName, 1024, 768)
    puDataHistogram.Draw()
    c.SaveAs("{oF}/{n}.pdf".format(oF=inputArguments.outputFolder, n=outputObjectName))
    inputFileObject.Close()
    # Next 
    for production_type in ["gluino", "squark"]:
        for signal_selection in ["signal", "signal_loose"]:
            source_file_path = "{sF}/PUWeights_{y}_{p}_{s}.root".format(sF=sourceFolder_MC, y=year, p=production_type, s=signal_selection)
            print("Extracting plots from file: {s_f_p}".format(s_f_p=source_file_path))
            inputFileObject = ROOT.TFile.Open(source_file_path, "READ")
            if ((inputFileObject.IsZombie() == ROOT.kTRUE) or not(inputFileObject.IsOpen() == ROOT.kTRUE)):
                sys.exit("ERROR: Unable to open file {f}".format(f=source_file_path))
            puWeightsHistogram = ROOT.TH1D()
            inputFileObject.GetObject("pileupWeights", puWeightsHistogram)
            puWeightsHistogram.GetXaxis().SetTitle("PU")
            puWeightsHistogram.GetYaxis().SetTitle("weight")
            puWeightsHistogram.SetTitle(("PU weights, {y}, {n} selection, di{p} production").format(y=year, n=selection_names[signal_selection], p=production_type))
            outputObjectName = "pileup_weights_vs_pu_{y}_{p}_{s}".format(y=year, p=production_type, s=signal_selection)
            c = ROOT.TCanvas("output_" + outputObjectName, "output_" + outputObjectName, 1024, 768)
            puWeightsHistogram.Draw()
            c.SaveAs("{oF}/{n}.pdf".format(oF=inputArguments.outputFolder, n=outputObjectName))
            inputFileObject.Close()

print("All done!")
