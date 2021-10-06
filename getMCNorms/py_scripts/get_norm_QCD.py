#!/usr/bin/env python

from __future__ import print_function, division

import os

os.system("rm -vf *.pyc")

import sys, signal, argparse, subprocess
import stealthEnv

import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

inputArgumentsParser = argparse.ArgumentParser(description='Run analysis chain.')
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

processes_BKG = ["DiPhotonJets", "GJetHT", "HighHTQCD"]
sources = {
    "data": "{eP}/{i}/histograms_data_pureQCD.root".format(eP=stealthEnv.EOSPrefix, i=histogramsSourceDirectory),
    "DiPhotonJets": "{eP}/{i}/histograms_DiPhotonJets_pureQCD.root".format(eP=stealthEnv.EOSPrefix, i=histogramsSourceDirectory),
    "GJetHT": "{eP}/{i}/histograms_GJetHT_pureQCD.root".format(eP=stealthEnv.EOSPrefix, i=histogramsSourceDirectory),
    "HighHTQCD": "{eP}/{i}/histograms_HighHTQCD_pureQCD.root".format(eP=stealthEnv.EOSPrefix, i=histogramsSourceDirectory)
}
colors = {
    "data": ROOT.kBlack,
    "DiPhotonJets": ROOT.kRed,
    "GJetHT": ROOT.kBlue,
    "HighHTQCD": ROOT.kGreen
}

source_file_objects = {}
for process in (processes_BKG + ["data"]):
    source_file_objects[process] = ROOT.TFile.Open(sources[process], "READ")
    if (((source_file_objects[process]).IsZombie() == ROOT.kTRUE) or (not((source_file_objects[process].IsOpen()) == ROOT.kTRUE))):
        sys.exit("ERROR: Unable to open file {f}".format(f=sources[process]))

for nJetsBin in range(2, 7):
    input_histograms = {}
    output_canvas = ROOT.TCanvas("pT_leadingJet_{n}JetsBin".format(n=nJetsBin), "pT_leadingJet_{n}JetsBin".format(n=nJetsBin), 1200, 1024)
    output_stack = ROOT.THStack("pT_leadingJet_{n}JetsBin".format(n=nJetsBin), "pT of leading jet, {n} jets bin;pT;events/bin".format(n=nJetsBin))
    ROOT.gPad.SetLogy()
    for process in (processes_BKG):
        input_histograms[process] = ROOT.TH1D()
        (source_file_objects[process]).GetObject("pT_leadingJet_{n}JetsBin".format(n=nJetsBin), input_histograms[process])
        if input_histograms[process]:
            input_histograms[process].SetLineColor(colors[process])
            input_histograms[process].SetFillColorAlpha(colors[process], 0.75)
            output_stack.Add(input_histograms[process])
        else:
            sys.exit("ERROR: unable to find histogram named \"pT_leadingJet_{n}JetsBin\" in input file for process {p}".format(n=nJetsBin, p=process))
    output_stack.Draw("nostack")
    input_histograms["data"] = ROOT.TH1D()
    (source_file_objects["data"]).GetObject("pT_leadingJet_{n}JetsBin".format(n=nJetsBin), input_histograms["data"])
    if input_histograms["data"]:
        input_histograms["data"].SetLineColor(colors["data"])
        input_histograms["data"].Draw("SAME")
    else:
        sys.exit("ERROR: unable to find histogram named \"pT_leadingJet_{n}JetsBin\" in input file for data.".format(n=nJetsBin))
    ROOT.gPad.Update()
    output_canvas.SaveAs("{o}/leadingJetPT_pureQCD_{n}JetsBin.pdf".format(o=output_folder, n=nJetsBin))

for process in (processes_BKG + ["data"]):
    source_file_objects[process].Close()

print("All done!")
