#!/usr/bin/env python

from __future__ import print_function, division

import os

os.system("rm -vf *.pyc")

import sys, signal, argparse, subprocess
import stealthEnv

import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

inputArgumentsParser = argparse.ArgumentParser(description='Run script to get relative MC norms.')
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArguments = inputArgumentsParser.parse_args()

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)
analysisOutputDirectory = "{aR}/analysis{oI}".format(aR=stealthEnv.analysisRoot, oI=optional_identifier)
histogramsSourceDirectory = "{sER}/analysisEOSAreas/analysis{oI}/MCNorms".format(sER=stealthEnv.stealthEOSRoot, oI=optional_identifier)
output_folder = "{aOD}/MCNorms".format(aOD=analysisOutputDirectory)
if not(os.path.isdir(output_folder)): subprocess.check_call("mkdir -p {oF}".format(oF=output_folder), shell=True, executable="/bin/bash")

def signal_handler(sig, frame):
    sys.exit("Terminated by user.")
signal.signal(signal.SIGINT, signal_handler)

def get_weighted_sum_events(input_th1, range_min, range_max):
    weighted_sum_nevents = 0
    for bin_counter in range(input_th1.GetXaxis().FindFixBin(range_min), 1+input_th1.GetXaxis().FindFixBin(range_max)):
        weighted_sum_nevents += input_th1.GetBinContent(bin_counter)
    return weighted_sum_nevents

processes_BKG = ["DiPhotonJets", "GJetHT", "HighHTQCD"]
sources = {}
for selection in ["pureQCD", "singlephoton"]:
    sources[selection] = {}
    for process in (processes_BKG + ["data"]):
        sources[selection][process] = "{eP}/{i}/histograms_{p}_{s}.root".format(eP=stealthEnv.EOSPrefix, i=histogramsSourceDirectory, p=process, s=selection)

colors = {
    "data": ROOT.kBlack,
    "DiPhotonJets": ROOT.kRed,
    "GJetHT": ROOT.kBlue,
    "HighHTQCD": ROOT.kGreen
}

source_file_objects = {
    "pureQCD": {},
    "singlephoton": {}
}

integrals = {}
for selection in ["pureQCD", "singlephoton"]:
    integrals[selection] = {}
    for process in (processes_BKG + ["data"]):
        integrals[selection][process] = {}

namePrefixes_histogramsToGet = {
    "pureQCD": "pT_leadingJet",
    "singlephoton": "pT_leadingPhoton"
}
titlePrefixes = {
    "pureQCD": "pT of  leading jet",
    "singlephoton": "pT of leading photon"
}
normalization_ranges = {
    "pureQCD": (200.5, 999.5),
    "singlephoton": (300.5, 699.5)
}

for selection in ["pureQCD", "singlephoton"]:
    for process in (processes_BKG + ["data"]):
        source_file_objects[selection][process] = ROOT.TFile.Open(sources[selection][process], "READ")
        if (((source_file_objects[selection][process]).IsZombie() == ROOT.kTRUE) or (not((source_file_objects[selection][process].IsOpen()) == ROOT.kTRUE))):
            sys.exit("ERROR: Unable to open file {f}".format(f=sources[selection][process]))

    for nJetsBin in range(2, 7):
        input_histograms = {}
        output_canvas = ROOT.TCanvas("{s}_{pr}_{n}JetsBin".format(s=selection, pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), "{s}_{pr}_{n}JetsBin".format(s=selection, pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), 1200, 1024)
        output_stack = ROOT.THStack("{s}_{pr}_{n}JetsBin".format(s=selection, pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), "{tp}, {n} jets bin;pT;events/bin".format(tp=titlePrefixes[selection], n=nJetsBin))
        ROOT.gPad.SetLogy()
        for process in (processes_BKG):
            input_histograms[process] = ROOT.TH1D()
            (source_file_objects[selection][process]).GetObject("{pr}_{n}JetsBin".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), input_histograms[process])
            if input_histograms[process]:
                input_histograms[process].SetLineColor(colors[process])
                input_histograms[process].SetFillColorAlpha(colors[process], 0.75)
                output_stack.Add(input_histograms[process])
                integrals[selection][process][nJetsBin] = get_weighted_sum_events(input_histograms[process], normalization_ranges[selection][0], normalization_ranges[selection][1])
                print("integrals[\"{s}\"][\"{p}\"][{n}]: {i}".format(s=selection, p=process, n=nJetsBin, i=integrals[selection][process][nJetsBin]))
            else:
                sys.exit("ERROR: unable to find histogram named \"{pr}_{n}JetsBin\" in input file for process {p}".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin, p=process))
        # output_stack.Draw("nostack")
        output_stack.Draw()
        input_histograms["data"] = ROOT.TH1D()
        (source_file_objects[selection]["data"]).GetObject("{pr}_{n}JetsBin".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), input_histograms["data"])
        if input_histograms["data"]:
            input_histograms["data"].SetLineColor(colors["data"])
            input_histograms["data"].Draw("SAME")
            integrals[selection]["data"][nJetsBin] = get_weighted_sum_events(input_histograms["data"], normalization_ranges[selection][0], normalization_ranges[selection][1])
            print("integrals[\"{s}\"][\"data\"][{n}]: {i}".format(s=selection, n=nJetsBin, i=integrals[selection]["data"][nJetsBin]))
        else:
            sys.exit("ERROR: unable to find histogram named \"{pr}_{n}JetsBin\" in input file for data.".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin))
        ROOT.gPad.Update()
        output_canvas.SaveAs("{o}/{s}_{pr}_{n}JetsBin_preKCorrection.pdf".format(s=selection, pr=namePrefixes_histogramsToGet[selection], o=output_folder, n=nJetsBin))

    for process in (processes_BKG + ["data"]):
        source_file_objects[selection][process].Close()

print("All done!")
