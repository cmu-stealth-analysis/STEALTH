#!/usr/bin/env python

from __future__ import print_function, division

import os

os.system("rm -vf *.pyc")

import sys, signal, argparse, subprocess
import numpy
import stealthEnv
import tmGeneralUtils

import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

inputArgumentsParser = argparse.ArgumentParser(description='Run script to get relative MC norms.')
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArgumentsParser.add_argument('--runUnblinded', action='store_true', help="If this flag is set, then the signal region ST distributions are unblinded.")
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
for selection in (selections + ["diphoton"]):
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
    "singlephoton": "pT_leadingPhoton"
}
titlePrefixes = {
    "pureQCD": "pT of  leading jet",
    "singlephoton": "pT of leading photon",
}
normalization_ranges = {
    "pureQCD": (300.5, 699.5),
    "singlephoton": (300.5, 699.5),
}
xLabels = {
    "pureQCD": "jet pT",
    "singlephoton": "photon pT",
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
        input_histograms["data"] = ROOT.TH1D()
        (source_file_objects[selection]["data"]).GetObject("{pr}_{n}JetsBin".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), input_histograms["data"])
        if input_histograms["data"]:
            input_histograms["data"].SetLineColor(colors["data"])
            input_histograms["data"].Draw()
            ROOT.gPad.Update()
            input_histograms["data"].GetXaxis().SetRangeUser(normalization_ranges[selection][0], normalization_ranges[selection][1])
            integrals[selection]["data"][nJetsBin] = get_weighted_sum_events(input_histograms["data"], normalization_ranges[selection][0], normalization_ranges[selection][1])
            print("integrals[\"{s}\"][\"data\"][{n}]: {i}".format(s=selection, n=nJetsBin, i=integrals[selection]["data"][nJetsBin]))
        else:
            sys.exit("ERROR: unable to find histogram named \"{pr}_{n}JetsBin\" in input file for data.".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin))
        ROOT.gPad.Update()
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
        output_stack.Draw("HIST SAME")
        ROOT.gPad.Update()
        output_stack.GetXaxis().SetRangeUser(normalization_ranges[selection][0], normalization_ranges[selection][1])
        ROOT.gPad.Update()
        input_histograms["data"].Draw("SAME")
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
        input_histograms_unscaled["data"] = ROOT.TH1D()
        (source_file_objects[selection]["data"]).GetObject("{pr}_{n}JetsBin".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), input_histograms_unscaled["data"])
        if input_histograms_unscaled["data"]:
            input_histograms_unscaled["data"].SetLineColor(colors["data"])
            input_histograms_unscaled["data"].Draw()
            ROOT.gPad.Update()
            input_histograms_unscaled["data"].GetXaxis().SetRangeUser(normalization_ranges[selection][0], normalization_ranges[selection][1])
        else:
            sys.exit("ERROR: unable to find histogram named \"{pr}_{n}JetsBin\" in input file for data.".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin))
        ROOT.gPad.Update()
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
        ROOT.gPad.Update()
        output_stack.GetXaxis().SetRangeUser(normalization_ranges[selection][0], normalization_ranges[selection][1])
        ROOT.gPad.Update()
        input_histograms_unscaled["data"].Draw("SAME")
        ROOT.gPad.Update()
        output_canvas.SaveAs("{o}/{s}_{pr}_{n}JetsBin_postKCorrection.pdf".format(s=selection, pr=namePrefixes_histogramsToGet[selection], o=output_folder, n=nJetsBin))

    for process in (processes_BKG + ["data"]):
        source_file_objects[selection][process].Close()

# Step 4: diphoton plots, pre-normalization
K_normSTBin = {}
diphoton_plots_to_extract = ["diphoton_invMass_zeroJets", "diphoton_nJets_in_normST"]
diphoton_plots_to_extract_source_names = {
    "diphoton_invMass_zeroJets": "invMass_zeroJets",
    "diphoton_nJets_in_normST": "nJets_in_normST"
}
diphoton_plots_to_extract_source_titles = {
    "diphoton_invMass_zeroJets": "diphoton invariant mass (2 tight #gamma);m;nEvents/bin",
    "diphoton_nJets_in_normST": "nJets distribution, 1200.0 GeV < ST < 1300.0 GeV;nJets bin;nEvents"
}
diphoton_plots_to_extract_xranges = {
    "diphoton_invMass_zeroJets": (140., 200.),
    "diphoton_nJets_in_normST": (0., 7.)
}
diphoton_plots_to_extract_yranges = {
    "diphoton_invMass_zeroJets": (0., 9000.),
    "diphoton_nJets_in_normST": (0., 400.)
}
diphoton_plots_to_extract_logScale = {
    "diphoton_invMass_zeroJets": False,
    "diphoton_nJets_in_normST": False
}
for process in (processes_BKG + ["data"]):
    source_file_objects["diphoton"][process] = ROOT.TFile.Open(sources["diphoton"][process], "READ")
    if (((source_file_objects["diphoton"][process]).IsZombie() == ROOT.kTRUE) or (not((source_file_objects["diphoton"][process].IsOpen()) == ROOT.kTRUE))):
        sys.exit("ERROR: Unable to open file {f}".format(f=sources["diphoton"][process]))

for diphoton_plot_to_extract in diphoton_plots_to_extract:
    output_canvas = ROOT.TCanvas(diphoton_plot_to_extract, diphoton_plot_to_extract, 1200, 1024)
    output_stack = ROOT.THStack(diphoton_plot_to_extract, diphoton_plots_to_extract_source_titles[diphoton_plot_to_extract])
    if diphoton_plots_to_extract_logScale[diphoton_plot_to_extract]:
        ROOT.gPad.SetLogy()
    else:
        ROOT.gPad.SetLogy(0)
    input_histograms_raw = {}
    input_histograms_scaled = {}
    input_histograms_raw["data"] = ROOT.TH1D()
    (source_file_objects["diphoton"]["data"]).GetObject(diphoton_plots_to_extract_source_names[diphoton_plot_to_extract], input_histograms_raw["data"])
    if input_histograms_raw["data"]:
        input_histograms_raw["data"].SetLineColor(colors["data"])
        input_histograms_raw["data"].Draw()
        ROOT.gPad.Update()
        input_histograms_raw["data"].GetXaxis().SetRangeUser(diphoton_plots_to_extract_xranges[diphoton_plot_to_extract][0], diphoton_plots_to_extract_xranges[diphoton_plot_to_extract][1])
        input_histograms_raw["data"].GetYaxis().SetRangeUser(diphoton_plots_to_extract_yranges[diphoton_plot_to_extract][0], diphoton_plots_to_extract_yranges[diphoton_plot_to_extract][1])
    else:
        sys.exit("ERROR: unable to find histogram named \"{n}\" in input file for data.".format(n=diphoton_plots_to_extract_source_names[diphoton_plot_to_extract]))
    ROOT.gPad.Update()
    for process in (processes_BKG):
        input_histograms_raw[process] = ROOT.TH1D()
        (source_file_objects["diphoton"][process]).GetObject(diphoton_plots_to_extract_source_names[diphoton_plot_to_extract], input_histograms_raw[process])
        if input_histograms_raw[process]:
            input_histograms_raw[process].SetLineColor(colors[process])
            input_histograms_raw[process].SetFillColorAlpha(colors[process], 0.75)
            input_histograms_scaled[process] = input_histograms_raw[process].Clone()
            input_histograms_scaled[process].SetName((input_histograms_raw[process]).GetName() + "_scaled")
            if (diphoton_plot_to_extract == "diphoton_invMass_zeroJets"):
                if (process in K_fit):
                    scale = K_fit[process][2]
                    input_histograms_scaled[process].Scale(scale)
            elif (diphoton_plot_to_extract == "diphoton_nJets_in_normST"):
                if (process in K_fit):
                    for nJetsBin in range(2, 7):
                        scale = K_fit[process][nJetsBin]
                        bin_content_unscaled = (input_histograms_raw[process]).GetBinContent((input_histograms_raw[process]).FindFixBin(1.0*nJetsBin))
                        input_histograms_scaled[process].SetBinContent((input_histograms_scaled[process]).FindFixBin(1.0*nJetsBin), scale*bin_content_unscaled)
                        bin_error_unscaled = (input_histograms_raw[process]).GetBinError((input_histograms_raw[process]).FindFixBin(1.0*nJetsBin))
                        input_histograms_scaled[process].SetBinError((input_histograms_scaled[process]).FindFixBin(1.0*nJetsBin), scale*bin_error_unscaled)
            else:
                sys.exit("ERROR: diphoton plot {p} not supported.".format(p=diphoton_plot_to_extract))
            output_stack.Add(input_histograms_scaled[process])
        else:
            sys.exit("ERROR: unable to find histogram named \"{n}\" in input file for process {p}.".format(n=diphoton_plots_to_extract_source_names[diphoton_plot_to_extract], p=process))
    if (diphoton_plot_to_extract == "diphoton_nJets_in_normST"):
        for nJetsBin in range(2, 7):
            sum_norms = 0.
            for process in (processes_BKG):
                sum_norms += (input_histograms_scaled[process]).GetBinContent((input_histograms_scaled[process]).FindFixBin(1.0*nJetsBin))
            data_norm = (input_histograms_raw["data"]).GetBinContent((input_histograms_raw["data"]).FindFixBin(1.0*nJetsBin))
            print("In {n} jets bin, sum_norms: {sn}, data_norm: {dn}".format(n=nJetsBin, sn=sum_norms, dn=data_norm))
            K_normSTBin[nJetsBin] = data_norm/sum_norms
    output_stack.Draw("HIST SAME")
    ROOT.gPad.Update()
    output_stack.GetXaxis().SetRangeUser(diphoton_plots_to_extract_xranges[diphoton_plot_to_extract][0], diphoton_plots_to_extract_xranges[diphoton_plot_to_extract][1])
    ROOT.gPad.Update()
    input_histograms_raw["data"].Draw("SAME")
    ROOT.gPad.Update()
    output_canvas.SaveAs("{o}/{p}_uncorrected.pdf".format(o=output_folder, p=diphoton_plot_to_extract))

for process in (processes_BKG + ["data"]):
    source_file_objects["diphoton"][process].Close()

# Step 5: plot ST distributions
plots_to_extract = ["ST_2JetsBin", "ST_3JetsBin", "ST_4JetsBin", "ST_5JetsBin", "ST_6JetsBin"]
plots_to_extract_source_names = {
    "ST_2JetsBin": "ST_2JetsBin",
    "ST_3JetsBin": "ST_3JetsBin",
    "ST_4JetsBin": "ST_4JetsBin",
    "ST_5JetsBin": "ST_5JetsBin",
    "ST_6JetsBin": "ST_6JetsBin"
}
plots_to_extract_source_nJetsBins = {
    "ST_2JetsBin": 2,
    "ST_3JetsBin": 3,
    "ST_4JetsBin": 4,
    "ST_5JetsBin": 5,
    "ST_6JetsBin": 6
}
plots_to_extract_source_titles = {
    "ST_2JetsBin": "ST distribution, 2 Jets;ST;nEvents/bin",
    "ST_3JetsBin": "ST distribution, 3 Jets;ST;nEvents/bin",
    "ST_4JetsBin": "ST distribution, 4 Jets;ST;nEvents/bin",
    "ST_5JetsBin": "ST distribution, 5 Jets;ST;nEvents/bin",
    "ST_6JetsBin": "ST distribution, 6 Jets;ST;nEvents/bin",
}
plots_to_extract_yranges = {
    "ST_2JetsBin": (0.001, 5.),
    "ST_3JetsBin": (0.001, 5.),
    "ST_4JetsBin": (0.001, 5.),
    "ST_5JetsBin": (0.001, 5.),
    "ST_6JetsBin": (0.001, 5.)
}
plots_to_extract_logScale = {
    "ST_2JetsBin": True,
    "ST_3JetsBin": True,
    "ST_4JetsBin": True,
    "ST_5JetsBin": True,
    "ST_6JetsBin": True
}
ST_distribution_is_blinded = {
    "ST_2JetsBin": False,
    "ST_3JetsBin": False,
    "ST_4JetsBin": not(inputArguments.runUnblinded),
    "ST_5JetsBin": not(inputArguments.runUnblinded),
    "ST_6JetsBin": not(inputArguments.runUnblinded)
}
for process in (processes_BKG + ["data"]):
    source_file_objects["diphoton"][process] = ROOT.TFile.Open(sources["diphoton"][process], "READ")
    if (((source_file_objects["diphoton"][process]).IsZombie() == ROOT.kTRUE) or (not((source_file_objects["diphoton"][process].IsOpen()) == ROOT.kTRUE))):
        sys.exit("ERROR: Unable to open file {f}".format(f=sources["diphoton"][process]))

for plot_to_extract in plots_to_extract:
    output_canvas = ROOT.TCanvas(plot_to_extract, plot_to_extract, 1200, 1024)
    output_stack = ROOT.THStack(plot_to_extract, plots_to_extract_source_titles[plot_to_extract])
    if plots_to_extract_logScale[plot_to_extract]:
        ROOT.gPad.SetLogy()
    else:
        ROOT.gPad.SetLogy(0)
    input_histograms_raw = {}
    input_histograms_scaled = {}
    if ST_distribution_is_blinded[plot_to_extract]:
        print("ST distribution is blinded. Not plotting data.")
    else:
        input_histograms_raw["data"] = ROOT.TH1D()
        (source_file_objects["diphoton"]["data"]).GetObject(plots_to_extract_source_names[plot_to_extract], input_histograms_raw["data"])
        if input_histograms_raw["data"]:
            input_histograms_raw["data"].SetLineColor(colors["data"])
            input_histograms_raw["data"].Draw()
            ROOT.gPad.Update()
            input_histograms_raw["data"].GetYaxis().SetRangeUser(plots_to_extract_yranges[plot_to_extract][0], plots_to_extract_yranges[plot_to_extract][1])
        else:
            sys.exit("ERROR: unable to find histogram named \"{n}\" in input file for data.".format(n=plots_to_extract_source_names[plot_to_extract]))
        ROOT.gPad.Update()
    histograms_sum = None
    for process in (processes_BKG):
        input_histograms_raw[process] = ROOT.TH1D()
        (source_file_objects["diphoton"][process]).GetObject(plots_to_extract_source_names[plot_to_extract], input_histograms_raw[process])
        if input_histograms_raw[process]:
            input_histograms_raw[process].SetLineColor(colors[process])
            input_histograms_raw[process].SetFillColorAlpha(colors[process], 0.75)
            input_histograms_scaled[process] = input_histograms_raw[process].Clone()
            input_histograms_scaled[process].SetName((input_histograms_raw[process]).GetName() + "_scaled")
            if (process in K_fit):
                K_scale = K_fit[process][plots_to_extract_source_nJetsBins[plot_to_extract]]
                input_histograms_scaled[process].Scale(K_scale)
            if (histograms_sum is None):
                histograms_sum = (input_histograms_scaled[process]).Clone()
                histograms_sum.SetName((input_histograms_scaled[process]).GetName() + "_sum")
            else:
                histograms_sum.Add(input_histograms_scaled[process])
        else:
            sys.exit("ERROR: unable to find histogram named \"{n}\" in input file for process {p}.".format(n=plots_to_extract_source_names[plot_to_extract], p=process))
    # normalizations
    integral_histograms_sum = None
    integral_data = None
    if ST_distribution_is_blinded[plot_to_extract]:
        print("ST distribution is blinded. Not normalizing.")
    else:
        integral_histograms_sum = histograms_sum.Integral(2, histograms_sum.GetXaxis().GetNbins(), "width")
        integral_data = (input_histograms_raw["data"]).Integral(2, (input_histograms_raw["data"]).GetXaxis().GetNbins(), "width")
    for process in (processes_BKG):
        if not(ST_distribution_is_blinded[plot_to_extract]):
            input_histograms_scaled[process].Scale(integral_data/integral_histograms_sum)
        output_stack.Add(input_histograms_scaled[process])
    output_stack_draw_options = "HIST"
    if not(ST_distribution_is_blinded[plot_to_extract]):
        output_stack_draw_options += " SAME"
    output_stack.Draw(output_stack_draw_options)
    ROOT.gPad.Update()
    if not(ST_distribution_is_blinded[plot_to_extract]):
        input_histograms_raw["data"].Draw("SAME")
        ROOT.gPad.Update()
    output_canvas.SaveAs("{o}/{p}_scaled.pdf".format(o=output_folder, p=plot_to_extract))

for process in (processes_BKG + ["data"]):
    source_file_objects["diphoton"][process].Close()

# Print out normalizations to a LaTeX-formatted table, and store them in a dictionary
norms_to_save = []
norms_tex_file_handle = open("{o}/norm_values.tex".format(o=output_folder), 'w')
norms_tex_file_handle.write("\\begin{tabular}{|p{0.25\\textwidth}|p{0.175\\textwidth}|p{0.175\\textwidth}|}\n")
norms_tex_file_handle.write("  \\hline\n")
norms_tex_file_handle.write("  nJets bin & HighHTQCD & GJetHT \\\\ \\hline\n")
for nJetsBin in range(2, 7):
    nJetsString = "{n}".format(n=nJetsBin)
    if (nJetsBin == 6): nJetsString = "$\\geq$ 6"
    norms_tex_file_handle.write("  {n} & {nQCD:.2f} & {nGJet:.2f} \\\\ \\hline\n".format(n=nJetsString, nQCD=K_fit["HighHTQCD"][nJetsBin], nGJet=K_fit["GJetHT"][nJetsBin], nDiph=1.0, nOverall=K_normSTBin[nJetsBin]))
    for process in processes_BKG:
        norm_value = K_normSTBin[nJetsBin]
        if (process in K_fit):
            norm_value *= K_fit[process][nJetsBin]
        # norms_to_save.append(tuple(["float", "norm_values_{p}_{n}JetsBin".format(p=process, n=nJetsBin), norm_value]))
        norms_to_save.append(tuple(["float", "norm_values_{p}_{n}JetsBin".format(p=process, n=nJetsBin), 1.0]))
norms_tex_file_handle.write("\\end{tabular}\n")
norms_tex_file_handle.close()
tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=norms_to_save, outputFilePath=("{o}/norm_values_nominal.dat".format(o=output_folder)))
print("K_fit:")
tmGeneralUtils.prettyPrintDictionary(K_fit)
print("K_normSTBin:")
tmGeneralUtils.prettyPrintDictionary(K_normSTBin)

print("All done!")
