#!/usr/bin/env python

from __future__ import print_function, division

import os

os.system("rm -vf *.pyc")

import sys, signal, argparse, subprocess
import numpy
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
# selections = ["pureQCD", "singlephoton", "diphoton"]
selections = ["pureQCD", "singlephoton"]
sources = {}
for selection in selections:
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
    "singlephoton": {},
    "diphoton": {}
}

integrals = {}
for selection in selections:
    integrals[selection] = {}
    for process in (processes_BKG + ["data"]):
        integrals[selection][process] = {}

namePrefixes_histogramsToGet = {
    "pureQCD": "pT_leadingJet",
    "singlephoton": "pT_leadingPhoton",
    "diphoton": "evtST"
}
titlePrefixes = {
    "pureQCD": "pT of  leading jet",
    "singlephoton": "pT of leading photon",
    "diphoton": "Event ST"
}
normalization_ranges = {
    "pureQCD": (300.5, 699.5),
    "singlephoton": (300.5, 699.5),
    "diphoton": (1000.5, 1299.5)
}
xLabels = {
    "pureQCD": "jet pT",
    "singlephoton": "photon pT",
    "diphoton": "ST"
}

# Step 1: Get coefficients for equation
for selection in selections:
    for process in (processes_BKG + ["data"]):
        source_file_objects[selection][process] = ROOT.TFile.Open(sources[selection][process], "READ")
        if (((source_file_objects[selection][process]).IsZombie() == ROOT.kTRUE) or (not((source_file_objects[selection][process].IsOpen()) == ROOT.kTRUE))):
            sys.exit("ERROR: Unable to open file {f}".format(f=sources[selection][process]))

    for nJetsBin in range(2, 7):
        input_histograms = {}
        output_canvas = ROOT.TCanvas("{s}_{pr}_{n}JetsBin".format(s=selection, pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), "{s}_{pr}_{n}JetsBin".format(s=selection, pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), 1200, 1024)
        output_stack = ROOT.THStack("{s}_{pr}_{n}JetsBin".format(s=selection, pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), "{tp}, {n} jets bin;{xl};events/bin".format(tp=titlePrefixes[selection], xl=xLabels[selection], n=nJetsBin))
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
        output_stack.Draw("HIST")
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

# Step 2: solve the equation
# We have, in each nJets bin, for each sample:
# K_QCD*QCD + K_GJet*GJet + K_diphoton*diphoton = data
# so that:
# /QCD_sel1    GJet_sel1    diphoton_sel1\  /K_QCD      \    /data_sel1\
# |                                      |  |           |    |         |
# |QCD_sel2    GJet_sel2    diphoton_sel2|  |K_GJet     | =  |data_sel2|
# |                                      |  |           |    |         |
# \QCD_sel3    GJet_sel3    diphoton_sel3/  \K_diphoton /    \data_sel3/
#
# where sel1 = pure QCD, sel2 = single photon, sel3 = diphoton.
#
# Then:
# /K_QCD      \   /QCD_sel1    GJet_sel1    diphoton_sel1\ -1  /data_sel1\
# |           |   |                                      |     |         |
# |K_GJet     | = |QCD_sel2    GJet_sel2    diphoton_sel2|     |data_sel2|
# |           |   |                                      |     |         |
# \K_diphoton /   \QCD_sel3    GJet_sel3    diphoton_sel3/     \data_sel3/

# OR just solve the 2D equation for the QCD and GJet norm factors
# /K_QCD \     /QCD_sel1    GJet_sel1\-1          /data_sel1\
# |      |  =  |                     |      X     |         |
# \K_GJet/     \QCD_sel2    GJet_sel2/            \data_sel2/

K_fit = {
    "HighHTQCD": {},
    "GJetHT": {}
}
for nJetsBin in range(2, 7):
    # MCObservationsMatrix = numpy.zeros((3, 3))
    # MCObservationsMatrix[0] = [integrals["pureQCD"]["HighHTQCD"][nJetsBin], integrals["pureQCD"]["GJetHT"][nJetsBin], integrals["pureQCD"]["DiPhotonJets"][nJetsBin]]
    # MCObservationsMatrix[1] = [integrals["singlephoton"]["HighHTQCD"][nJetsBin], integrals["singlephoton"]["GJetHT"][nJetsBin], integrals["singlephoton"]["DiPhotonJets"][nJetsBin]]
    # MCObservationsMatrix[2] = [integrals["diphoton"]["HighHTQCD"][nJetsBin], integrals["diphoton"]["GJetHT"][nJetsBin], integrals["diphoton"]["DiPhotonJets"][nJetsBin]]
    MCObservationsMatrix = numpy.zeros((2, 2))
    MCObservationsMatrix[0] = [integrals["pureQCD"]["HighHTQCD"][nJetsBin], integrals["pureQCD"]["GJetHT"][nJetsBin]]
    MCObservationsMatrix[1] = [integrals["singlephoton"]["HighHTQCD"][nJetsBin], integrals["singlephoton"]["GJetHT"][nJetsBin]]

    # dataColumn = numpy.zeros((3, 1))
    # dataColumn[0] = [integrals["pureQCD"]["data"][nJetsBin]]
    # dataColumn[1] = [integrals["singlephoton"]["data"][nJetsBin]]
    # dataColumn[2] = [integrals["diphoton"]["data"][nJetsBin]]
    dataColumn = numpy.zeros((2, 1))
    dataColumn[0] = [integrals["pureQCD"]["data"][nJetsBin]]
    dataColumn[1] = [integrals["singlephoton"]["data"][nJetsBin]]

    print("-"*100)
    print("nJets bin: {n}".format(n=nJetsBin))
    print("Equation to solve: Ax = y, where A:")
    print(MCObservationsMatrix)
    print("y:")
    print(dataColumn)
    Kvalues_nparray = numpy.matmul(numpy.linalg.inv(MCObservationsMatrix), dataColumn)
    print("Solution found... x:")
    print(Kvalues_nparray)
    Kvalues = Kvalues_nparray.tolist()
    if not(len(Kvalues) == 2): sys.exit("ERROR: Kvalues_nparray = {k} in unexpected format.".format(k=str(Kvalues_nparray)))
    K_fit["HighHTQCD"][nJetsBin] = float(Kvalues[0][0])
    K_fit["GJetHT"][nJetsBin] = float(Kvalues[1][0])

# Step 3: Plot scaled histograms
for selection in selections:
    for process in (processes_BKG + ["data"]):
        source_file_objects[selection][process] = ROOT.TFile.Open(sources[selection][process], "READ")
        if (((source_file_objects[selection][process]).IsZombie() == ROOT.kTRUE) or (not((source_file_objects[selection][process].IsOpen()) == ROOT.kTRUE))):
            sys.exit("ERROR: Unable to open file {f}".format(f=sources[selection][process]))

    for nJetsBin in range(2, 7):
        input_histograms_unscaled = {}
        histograms_scaled = {}
        output_canvas = ROOT.TCanvas("{s}_{pr}_{n}JetsBin_postScaleFix".format(s=selection, pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), "{s}_{pr}_{n}JetsBin_postScaleFix".format(s=selection, pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), 1200, 1024)
        output_stack = ROOT.THStack("{s}_{pr}_{n}JetsBin_postScaleFix".format(s=selection, pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), "{tp}, scaled, {n} jets bin;{xl};events/bin".format(tp=titlePrefixes[selection], xl=xLabels[selection], n=nJetsBin))
        ROOT.gPad.SetLogy()
        for process in (processes_BKG):
            if not(process in K_fit): continue
            input_histograms_unscaled[process] = ROOT.TH1D()
            (source_file_objects[selection][process]).GetObject("{pr}_{n}JetsBin".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), input_histograms_unscaled[process])
            if input_histograms_unscaled[process]:
                input_histograms_unscaled[process].SetLineColor(colors[process])
                input_histograms_unscaled[process].SetFillColorAlpha(colors[process], 0.75)
                # output_stack.Add(input_histograms_unscaled[process])
                histograms_scaled[process] = input_histograms_unscaled[process].Clone()
                histograms_scaled[process].SetName(input_histograms_unscaled[process].GetName() + "_scaled")
                print("Scaling by K_fit: {k}".format(k=K_fit[process][nJetsBin]))
                histograms_scaled[process].Scale(K_fit[process][nJetsBin])
                output_stack.Add(histograms_scaled[process])
            else:
                sys.exit("ERROR: unable to find histogram named \"{pr}_{n}JetsBin\" in input file for process {p}".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin, p=process))
        # output_stack.Draw("nostack")
        output_stack.Draw("HIST")
        input_histograms_unscaled["data"] = ROOT.TH1D()
        (source_file_objects[selection]["data"]).GetObject("{pr}_{n}JetsBin".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), input_histograms_unscaled["data"])
        if input_histograms_unscaled["data"]:
            input_histograms_unscaled["data"].SetLineColor(colors["data"])
            input_histograms_unscaled["data"].Draw("SAME")
        else:
            sys.exit("ERROR: unable to find histogram named \"{pr}_{n}JetsBin\" in input file for data.".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin))
        ROOT.gPad.Update()
        output_canvas.SaveAs("{o}/{s}_{pr}_{n}JetsBin_postKCorrection.pdf".format(s=selection, pr=namePrefixes_histogramsToGet[selection], o=output_folder, n=nJetsBin))

    for process in (processes_BKG + ["data"]):
        source_file_objects[selection][process].Close()

print("All done!")
