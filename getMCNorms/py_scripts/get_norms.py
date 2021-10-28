#!/usr/bin/env python

from __future__ import print_function, division

import os

os.system("rm -vf *.pyc")

import sys, signal, argparse, subprocess, math
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

N_JETS_NORM = 2
FRACTIONAL_TOLERANCE_FOR_CHECKS=0.001

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
    # "pureQCD": "ST_fineBinned",
    # "singlephoton": "ST_fineBinned"
}
titlePrefixes = {
    "pureQCD": "pT of leading jet",
    "singlephoton": "pT of leading photon"
    # "pureQCD": "ST",
    # "singlephoton": "ST"
}
normalization_ranges = {
    "pureQCD": (300.5, 699.5),
    "singlephoton": (300.5, 699.5)
    # "pureQCD": (1275.0, 3250.0),
    # "singlephoton": (1275.0, 3250.0)
}
xRanges = {
    "pureQCD": (300.0, 700.0),
    "singlephoton": (300.0, 700.0)
    # "pureQCD": (1200.0, 3500.0),
    # "singlephoton": (1200.0, 3500.0)
}
xLabels = {
    "pureQCD": "jet pT",
    "singlephoton": "photon pT"
    # "pureQCD": "ST",
    # "singlephoton": "ST"
}
yRanges = {
    "pureQCD": (1.0, 10000.0),
    "singlephoton": (1.0, 500.0)
    # "pureQCD": (0.001, 100.),
    # "singlephoton": (0.05, 500.0),
}

# Step 1: Get coefficients for equation
print("-"*200)
print("Building coefficients for two-variable equation...")
for selection in selections:
    for process in (processes_BKG + ["data"]):
        source_file_objects[selection][process] = ROOT.TFile.Open(sources[selection][process], "READ")
        if (((source_file_objects[selection][process]).IsZombie() == ROOT.kTRUE) or (not((source_file_objects[selection][process].IsOpen()) == ROOT.kTRUE))):
            sys.exit("ERROR: Unable to open file {f}".format(f=sources[selection][process]))

    for nJetsBin in range(2, 7):
        input_histograms = {}
        output_canvas = ROOT.TCanvas("{s}_{pr}_{n}JetsBin".format(s=selection, pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), "{s}_{pr}_{n}JetsBin".format(s=selection, pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), 1200, 1024)
        ROOT.gStyle.SetOptStat(0)
        output_stack = ROOT.THStack("{s}_{pr}_{n}JetsBin".format(s=selection, pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), "{tp}, {n} jets bin;{xl};events/bin".format(tp=titlePrefixes[selection], xl=xLabels[selection], n=nJetsBin))
        ROOT.gPad.SetLogy()
        input_histograms["data"] = ROOT.TH1D()
        (source_file_objects[selection]["data"]).GetObject("{pr}_{n}JetsBin".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), input_histograms["data"])
        if input_histograms["data"]:
            input_histograms["data"].SetLineColor(colors["data"])
            input_histograms["data"].Draw()
            input_histograms["data"].GetYaxis().SetRangeUser(yRanges[selection][0], yRanges[selection][1])
            ROOT.gPad.Update()
            input_histograms["data"].GetXaxis().SetRangeUser(xRanges[selection][0], xRanges[selection][1])
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
        output_stack.GetXaxis().SetRangeUser(xRanges[selection][0], xRanges[selection][1])
        ROOT.gPad.Update()
        input_histograms["data"].Draw("SAME")
        ROOT.gPad.Update()
        output_canvas.SaveAs("{o}/{s}_{pr}_{n}JetsBin_preKCorrection.pdf".format(s=selection, pr=namePrefixes_histogramsToGet[selection], o=output_folder, n=nJetsBin))

    for process in (processes_BKG + ["data"]):
        source_file_objects[selection][process].Close()
print("-"*200)

# Step 2: solve the equation
print("-"*200)
print("Solving two-variable equation for normalizations...")
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

# OR just sub diphoton scale = 1.0 and solve the 2D equation for the QCD and GJet norm factors
# /K_QCD \     /QCD_sel1    GJet_sel1\-1          /data_sel1 - diphoton_sel1\
# |      |  =  |                     |      X     |                         |
# \K_GJet/     \QCD_sel2    GJet_sel2/            \data_sel2 - diphoton_sel2/

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
    dataColumn[0] = [integrals["pureQCD"]["data"][nJetsBin] - integrals["pureQCD"]["DiPhotonJets"][nJetsBin]
    ]
    dataColumn[1] = [integrals["singlephoton"]["data"][nJetsBin] - integrals["singlephoton"]["DiPhotonJets"][nJetsBin]
    ]

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
print("-"*200)

# Step 3: Plot scaled histograms
print("-"*200)
print("Plotting post-correction histograms...")
for selection in selections:
    for process in (processes_BKG + ["data"]):
        source_file_objects[selection][process] = ROOT.TFile.Open(sources[selection][process], "READ")
        if (((source_file_objects[selection][process]).IsZombie() == ROOT.kTRUE) or (not((source_file_objects[selection][process].IsOpen()) == ROOT.kTRUE))):
            sys.exit("ERROR: Unable to open file {f}".format(f=sources[selection][process]))

    for nJetsBin in range(2, 7):
        input_histograms_unscaled = {}
        histograms_scaled = {}
        output_canvas = ROOT.TCanvas("{s}_{pr}_{n}JetsBin_postScaleFix".format(s=selection, pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), "{s}_{pr}_{n}JetsBin_postScaleFix".format(s=selection, pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), 1200, 1024)
        ROOT.gStyle.SetOptStat(0)
        output_stack = ROOT.THStack("{s}_{pr}_{n}JetsBin_postScaleFix".format(s=selection, pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), "{tp}, scaled, {n} jets bin;{xl};events/bin".format(tp=titlePrefixes[selection], xl=xLabels[selection], n=nJetsBin))
        ROOT.gPad.SetLogy()
        input_histograms_unscaled["data"] = ROOT.TH1D()
        (source_file_objects[selection]["data"]).GetObject("{pr}_{n}JetsBin".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), input_histograms_unscaled["data"])
        if input_histograms_unscaled["data"]:
            input_histograms_unscaled["data"].SetLineColor(colors["data"])
            input_histograms_unscaled["data"].Draw()
            input_histograms_unscaled["data"].GetYaxis().SetRangeUser(yRanges[selection][0], yRanges[selection][1])
            ROOT.gPad.Update()
            input_histograms_unscaled["data"].GetXaxis().SetRangeUser(xRanges[selection][0], xRanges[selection][1])
        else:
            sys.exit("ERROR: unable to find histogram named \"{pr}_{n}JetsBin\" in input file for data.".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin))
        ROOT.gPad.Update()
        for process in (processes_BKG):
            # if not(process in K_fit): continue
            input_histograms_unscaled[process] = ROOT.TH1D()
            (source_file_objects[selection][process]).GetObject("{pr}_{n}JetsBin".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin), input_histograms_unscaled[process])
            if input_histograms_unscaled[process]:
                input_histograms_unscaled[process].SetLineColor(colors[process])
                input_histograms_unscaled[process].SetFillColorAlpha(colors[process], 0.75)
                # output_stack.Add(input_histograms_unscaled[process])
                histograms_scaled[process] = input_histograms_unscaled[process].Clone()
                histograms_scaled[process].SetName(input_histograms_unscaled[process].GetName() + "_scaled")
                if process in K_fit:
                    print("Scaling by K_fit: {k}".format(k=K_fit[process][nJetsBin]))
                    histograms_scaled[process].Scale(K_fit[process][nJetsBin])
                output_stack.Add(histograms_scaled[process])
            else:
                sys.exit("ERROR: unable to find histogram named \"{pr}_{n}JetsBin\" in input file for process {p}".format(pr=namePrefixes_histogramsToGet[selection], n=nJetsBin, p=process))
        # output_stack.Draw("nostack")
        output_stack.Draw("HIST SAME")
        ROOT.gPad.Update()
        output_stack.GetXaxis().SetRangeUser(normalization_ranges[selection][0], normalization_ranges[selection][1])
        ROOT.gPad.Update()
        input_histograms_unscaled["data"].Draw("SAME")
        ROOT.gPad.Update()
        output_canvas.SaveAs("{o}/{s}_{pr}_{n}JetsBin_postKCorrection.pdf".format(s=selection, pr=namePrefixes_histogramsToGet[selection], o=output_folder, n=nJetsBin))

    for process in (processes_BKG + ["data"]):
        source_file_objects[selection][process].Close()
print("-"*200)

# Step 4: diphoton plots, pre-normalization
print("-"*200)
print("Making pre-normalization diphoton plots...")
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
    ROOT.gStyle.SetOptStat(0)
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
print("-"*200)

# Step 5: plot ST distributions for the GJet and QCD selections
print("-"*200)
print("Plotting ST distributions for GJet and QCD selections...")
plots_to_extract = ["ST_fineBinned_2JetsBin", "ST_fineBinned_3JetsBin", "ST_fineBinned_4JetsBin", "ST_fineBinned_5JetsBin", "ST_fineBinned_6JetsBin"]
plots_to_extract_source_names = {
    "ST_fineBinned_2JetsBin": "ST_fineBinned_2JetsBin",
    "ST_fineBinned_3JetsBin": "ST_fineBinned_3JetsBin",
    "ST_fineBinned_4JetsBin": "ST_fineBinned_4JetsBin",
    "ST_fineBinned_5JetsBin": "ST_fineBinned_5JetsBin",
    "ST_fineBinned_6JetsBin": "ST_fineBinned_6JetsBin"
}
plots_to_extract_source_nJetsBins = {
    "ST_fineBinned_2JetsBin": 2,
    "ST_fineBinned_3JetsBin": 3,
    "ST_fineBinned_4JetsBin": 4,
    "ST_fineBinned_5JetsBin": 5,
    "ST_fineBinned_6JetsBin": 6
}
plots_to_extract_source_titles = {
    "ST_fineBinned_2JetsBin": "ST distribution (fine-binned), 2 Jets;ST;nEvents/bin",
    "ST_fineBinned_3JetsBin": "ST distribution (fine-binned), 3 Jets;ST;nEvents/bin",
    "ST_fineBinned_4JetsBin": "ST distribution (fine-binned), 4 Jets;ST;nEvents/bin",
    "ST_fineBinned_5JetsBin": "ST distribution (fine-binned), 5 Jets;ST;nEvents/bin",
    "ST_fineBinned_6JetsBin": "ST distribution (fine-binned), 6 Jets;ST;nEvents/bin",
}
plots_to_extract_yranges = {
    "pureQCD": {
        "ST_fineBinned_2JetsBin": (0.001, 100.),
        "ST_fineBinned_3JetsBin": (0.001, 100.),
        "ST_fineBinned_4JetsBin": (0.001, 100.),
        "ST_fineBinned_5JetsBin": (0.001, 100.),
        "ST_fineBinned_6JetsBin": (0.001, 100.)
    },
    "singlephoton": {
        "ST_fineBinned_2JetsBin": (0.05, 500.),
        "ST_fineBinned_3JetsBin": (0.05, 500.),
        "ST_fineBinned_4JetsBin": (0.05, 500.),
        "ST_fineBinned_5JetsBin": (0.05, 500.),
        "ST_fineBinned_6JetsBin": (0.05, 500.)
    }
}
plots_to_extract_logScale = {
    "ST_fineBinned_2JetsBin": True,
    "ST_fineBinned_3JetsBin": True,
    "ST_fineBinned_4JetsBin": True,
    "ST_fineBinned_5JetsBin": True,
    "ST_fineBinned_6JetsBin": True
}
for selection in selections:
    for process in (processes_BKG + ["data"]):
        source_file_objects[selection][process] = ROOT.TFile.Open(sources[selection][process], "READ")
        if (((source_file_objects[selection][process]).IsZombie() == ROOT.kTRUE) or (not((source_file_objects[selection][process].IsOpen()) == ROOT.kTRUE))):
            sys.exit("ERROR: Unable to open file {f}".format(f=sources[selection][process]))

    ratio_values_and_errors = {}
    for plot_to_extract in plots_to_extract:
        nJetsBin = plots_to_extract_source_nJetsBins[plot_to_extract]
        output_canvas = ROOT.TCanvas(plot_to_extract, plot_to_extract, 1200, 600)
        ROOT.gStyle.SetOptStat(0)
        output_stack = ROOT.THStack(plot_to_extract, plots_to_extract_source_titles[plot_to_extract])
        if plots_to_extract_logScale[plot_to_extract]:
            ROOT.gPad.SetLogy()
        else:
            ROOT.gPad.SetLogy(0)
        input_histograms_raw = {}
        input_histograms_scaled = {}
        input_histograms_raw["data"] = ROOT.TH1D()
        (source_file_objects[selection]["data"]).GetObject(plots_to_extract_source_names[plot_to_extract], input_histograms_raw["data"])
        if input_histograms_raw["data"]:
            input_histograms_raw["data"].SetLineColor(colors["data"])
            input_histograms_raw["data"].Draw()
            ROOT.gPad.Update()
            input_histograms_raw["data"].GetYaxis().SetRangeUser(plots_to_extract_yranges[selection][plot_to_extract][0], plots_to_extract_yranges[selection][plot_to_extract][1])
        else:
            sys.exit("ERROR: unable to find histogram named \"{n}\" in input file for data.".format(n=plots_to_extract_source_names[plot_to_extract]))
        ROOT.gPad.Update()
        histograms_sum = None
        for process in (processes_BKG):
            input_histograms_raw[process] = ROOT.TH1D()
            (source_file_objects[selection][process]).GetObject(plots_to_extract_source_names[plot_to_extract], input_histograms_raw[process])
            if input_histograms_raw[process]:
                input_histograms_raw[process].SetLineColor(colors[process])
                input_histograms_raw[process].SetFillColorAlpha(colors[process], 0.75)
                input_histograms_scaled[process] = input_histograms_raw[process].Clone()
                input_histograms_scaled[process].SetName((input_histograms_raw[process]).GetName() + "_scaled")
                if (process in K_fit):
                    K_scale = K_fit[process][nJetsBin]
                    input_histograms_scaled[process].Scale(K_scale)
                if (histograms_sum is None):
                    histograms_sum = (input_histograms_scaled[process]).Clone()
                    histograms_sum.SetName((input_histograms_scaled[process]).GetName() + "_sum")
                else:
                    histograms_sum.Add(input_histograms_scaled[process])
            else:
                sys.exit("ERROR: unable to find histogram named \"{n}\" in input file for process {p}.".format(n=plots_to_extract_source_names[plot_to_extract], p=process))
        # normalizations
        integral_histograms_sum = histograms_sum.Integral(2, histograms_sum.GetXaxis().GetNbins(), "width")
        integral_data = (input_histograms_raw["data"]).Integral(2, (input_histograms_raw["data"]).GetXaxis().GetNbins(), "width")
        histograms_sum.Scale(integral_data/integral_histograms_sum)
        for process in (processes_BKG):
            print("Scaling all bkg histograms by factor: {f:.2f}".format(f=integral_data/integral_histograms_sum))
            input_histograms_scaled[process].Scale(integral_data/integral_histograms_sum)
            output_stack.Add(input_histograms_scaled[process])
        # Get ratio and ratio errors
        ratio_values_and_errors[nJetsBin] = []
        for STBinIndex in range(1, 1 + histograms_sum.GetXaxis().GetNbins()):
            # a = MC_diphoton (scaled)
            # b = MC_GJet (scaled)
            # c = MC_QCD (scaled)
            # s = a+b+c ( = sum of scaled MCs)
            # d = data
            # We have:
            # r = s/d
            # delta_r = r*sqrt((delta_s/s)^2 + (delta_d/d)^2)
            # Here, delta_s^2 = delta_a^2 + delta_b^2 + delta_c^2
            # so that
            # delta_r = r*sqrt(((delta_a^2 + delta_b^2 + delta_c^2)/s^2) + (delta_d^2/d^2))
            # Turns out this is all probably handled correctly by the Add method above...
            MC_sum_yields = histograms_sum.GetBinContent(STBinIndex)
            MC_sum_yields_error = histograms_sum.GetBinError(STBinIndex)
            data_observed = (input_histograms_raw["data"]).GetBinContent(STBinIndex)
            data_observed_error = (input_histograms_raw["data"]).GetBinError(STBinIndex)
            ratio = (data_observed/MC_sum_yields)
            ratio_error = ratio*math.sqrt(pow((1.0*data_observed_error)/data_observed, 2) + pow((1.0*MC_sum_yields_error)/MC_sum_yields, 2))
            bin_center = histograms_sum.GetXaxis().GetBinCenter(STBinIndex)
            bin_width = histograms_sum.GetXaxis().GetBinUpEdge(STBinIndex) - histograms_sum.GetXaxis().GetBinLowEdge(STBinIndex)
            ratio_values_and_errors[nJetsBin].append((STBinIndex, bin_center, ratio, bin_width/math.sqrt(12), ratio_error))
            # print("data: {d:.2f} +/- {dd:.2f}, sum: {s:.2f} +/- {ds:.2f}, diphoton: {di:.2f} +/- {ddi:.2f}, gjet: {gj:.2f} +/- {dgj:.2f}, qcd: {q:.2f} +/- {dq:.2f}; ratio_values_and_errors[nJetsBin] element: {e}".format(d=data_observed, dd=data_observed_error, s=MC_sum_yields, ds=MC_sum_yields_error, di=input_histograms_scaled["DiPhotonJets"].GetBinContent(STBinIndex), ddi=input_histograms_scaled["DiPhotonJets"].GetBinError(STBinIndex), gj=input_histograms_scaled["GJetHT"].GetBinContent(STBinIndex), dgj=input_histograms_scaled["GJetHT"].GetBinError(STBinIndex), q=input_histograms_scaled["HighHTQCD"].GetBinContent(STBinIndex), dq=input_histograms_scaled["HighHTQCD"].GetBinError(STBinIndex), e=str(ratio_values_and_errors[nJetsBin][-1])))
        output_stack.Draw("HIST SAME")
        ROOT.gPad.Update()
        input_histograms_raw["data"].Draw("SAME")
        ROOT.gPad.Update()
        output_canvas.SaveAs("{o}/{s}_{p}_postKCorrection.pdf".format(o=output_folder, s=selection, p=plot_to_extract))
        output_canvas_ratio = ROOT.TCanvas(plot_to_extract + "_ratio", plot_to_extract + "_ratio", 1200, 300)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gPad.SetLogy(0)
        data_mc_ratio_tgraph = ROOT.TGraphErrors()
        data_mc_ratio_tgraph.SetName("ratio_data_mc_{s}_{p}".format(s=selection, p=plot_to_extract))
        data_mc_ratio_tgraph.SetTitle("sum_MC / data")
        for bin_index, STVal, ratio, delta_STVal, delta_ratio in ratio_values_and_errors[nJetsBin]:
            ratioGraphBinIndex = data_mc_ratio_tgraph.GetN()
            data_mc_ratio_tgraph.SetPoint(ratioGraphBinIndex, STVal, ratio)
            data_mc_ratio_tgraph.SetPointError(ratioGraphBinIndex, delta_STVal, delta_ratio)
        data_mc_ratio_tgraph.GetXaxis().SetTitle(histograms_sum.GetXaxis().GetTitle())
        data_mc_ratio_tgraph.GetXaxis().SetLimits(histograms_sum.GetXaxis().GetXmin(), histograms_sum.GetXaxis().GetXmax())
        data_mc_ratio_tgraph.GetYaxis().SetTitle("ratio")
        data_mc_ratio_tgraph.GetHistogram().SetMinimum(-0.5)
        data_mc_ratio_tgraph.GetHistogram().SetMaximum(3.5)
        data_mc_ratio_tgraph.Draw("AP0")
        ROOT.gPad.Update()
        lineAt1 = ROOT.TLine(histograms_sum.GetXaxis().GetXmin(), 1., histograms_sum.GetXaxis().GetXmax(), 1.)
        lineAt1.SetLineColor(ROOT.kBlack)
        lineAt1.SetLineStyle(ROOT.kDashed)
        lineAt1.Draw()
        fit_function_string = "[0] + [1]*((x/{c:.4f}) - 1.0)".format(c=histograms_sum.GetXaxis().GetBinCenter(1))
        # print("Fitting function: {f}".format(f=fit_function_string))
        linear_fit_tf1 = ROOT.TF1("fit_linear_" + data_mc_ratio_tgraph.GetName(), fit_function_string, histograms_sum.GetXaxis().GetBinLowEdge(1), histograms_sum.GetXaxis().GetBinUpEdge(histograms_sum.GetXaxis().GetNbins()))
        linear_fit_tf1.SetParName(0, "const_{s}_{p}".format(s=selection, p=plot_to_extract))
        linear_fit_tf1.SetParameter(0, 1.0)
        linear_fit_tf1.SetParLimits(0, 0.0, 5.0)
        linear_fit_tf1.SetParName(1, "slope_{s}_{p}".format(s=selection, p=plot_to_extract))
        linear_fit_tf1.SetParameter(1, 0.0)
        linear_fit_tf1.SetParLimits(1, -5.0, 5.0)
        linear_fit_result = data_mc_ratio_tgraph.Fit(linear_fit_tf1, "QS0+")
        if not(linear_fit_result.Status() == 0):
            print("Warning: fit failed with fit options \"QS0+\". Now trying options \"QEX0S0+\"...")
            linear_fit_result = data_mc_ratio_tgraph.Fit(linear_fit_tf1, "QEX0S0+")
            if not(linear_fit_result.Status() == 0): sys.exit("ERROR: Unable to find fit for selection: {s}, plot_to_extract: {p}".format(s=selection, p=plot_to_extract))
        best_fit_const = linear_fit_result.Value(0)
        best_fit_const_error = linear_fit_result.ParError(0)
        best_fit_slope = linear_fit_result.Value(1)
        best_fit_slope_error = linear_fit_result.ParError(1)
        legend_data_mc_ratio = ROOT.TLegend(0.1, 0.7, 0.9, 0.9);
        legend_data_mc_ratio.SetFillStyle(0);
        linear_fit_tf1.SetLineColor(ROOT.kBlue)
        linear_fit_tf1.Draw("LSAME")
        ROOT.gPad.Update()
        legend_entry = legend_data_mc_ratio.AddEntry(linear_fit_tf1, "nominal fit: ({c:.3f} #pm {dc:.3f}) + ({s:.3f} #pm {ds:.3f})*((ST/{n:.2f}) - 1.0)".format(c=best_fit_const, dc=best_fit_const_error, s=best_fit_slope, ds=best_fit_slope_error, n=histograms_sum.GetXaxis().GetBinCenter(1)));
        legend_entry.SetLineColor(ROOT.kBlue)
        legend_entry.SetMarkerColor(ROOT.kBlue)
        legend_entry.SetTextColor(ROOT.kBlue)
        legend_data_mc_ratio.Draw()
        ROOT.gPad.Update()
        output_canvas_ratio.SaveAs("{o}/{s}_{p}_dataMCRatio_postKCorrection.pdf".format(o=output_folder, s=selection, p=plot_to_extract))

    for plot_to_extract in plots_to_extract:
        nJetsBin = plots_to_extract_source_nJetsBins[plot_to_extract]
        if (nJetsBin < 4): continue
        nJetsBinTitle = "{n} Jets".format(n=nJetsBin)
        if (nJetsBin == 6): nJetsBinTitle = "#geq 6 Jets"

        # Just to get x-axis info...
        input_data_histogram = ROOT.TH1D()
        (source_file_objects[selection]["data"]).GetObject(plots_to_extract_source_names[plot_to_extract], input_data_histogram)
        if not(input_data_histogram):
            sys.exit("ERROR: unable to find histogram named \"{n}\" in input file for data.".format(n=plots_to_extract_source_names[plot_to_extract]))

        output_canvas_mismodeling_ratio = ROOT.TCanvas("mismodeling_ratio_" + plot_to_extract, "mismodeling_ratio_" + plot_to_extract, 1200, 800)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gPad.SetLogy(0)
        mc_mismodeling_ratio_tgraph = ROOT.TGraphErrors()
        mc_mismodeling_ratio_tgraph.SetName("ratio_data_mc_{s}_{p}".format(s=selection, p=plot_to_extract))
        mc_mismodeling_ratio_tgraph.SetTitle("(sum_MC / data) ({nt}) / ((sum_MC / data) ({nn} Jets))".format(nt=nJetsBinTitle, nn=N_JETS_NORM))
        if not(len(ratio_values_and_errors[nJetsBin]) == len(ratio_values_and_errors[N_JETS_NORM])):
            sys.exit("ERROR: at nJetsBin={n}, len(ratio_values_and_errors) = {len1}; while in the norm bin, len(ratio_values_and_errors) = {len2}".format(n=nJetsBin, len1=len(ratio_values_and_errors[nJetsBin]), len2=len(ratio_values_and_errors[N_JETS_NORM])))
        for ratio_index in range(len(ratio_values_and_errors[nJetsBin])):
            bin_index, STVal, ratio, delta_STVal, delta_ratio = ratio_values_and_errors[nJetsBin][ratio_index]
            bin_index_norm, STVal_norm, ratio_norm, delta_STVal_norm, delta_ratio_norm = ratio_values_and_errors[N_JETS_NORM][ratio_index]
            if ((not(bin_index == bin_index_norm)) or
                (math.fabs(1.0 - (STVal/STVal_norm)) > FRACTIONAL_TOLERANCE_FOR_CHECKS) or
                (math.fabs(1.0 - (delta_STVal/delta_STVal_norm)) > FRACTIONAL_TOLERANCE_FOR_CHECKS)):
                sys.exit("Error: incompatible ratio values and errors at nJetsBin={n}, ratio_index: {ri}; ratio_values_and_errors[nJetsBin][ratio_index]: {r1}, ratio_values_and_errors[N_JETS_NORM][ratio_index]: {r2}".format(n=nJetsBin, ri=ratio_index, r1=str(ratio_values_and_errors[nJetsBin][ratio_index]), r2=str(ratio_values_and_errors[N_JETS_NORM][ratio_index])))
            mismodeling_ratio_tgraph_index = mc_mismodeling_ratio_tgraph.GetN()
            mismodeling_ratio = ratio/ratio_norm
            mismodeling_ratio_error = mismodeling_ratio*math.sqrt(pow(delta_ratio/ratio, 2) + pow(delta_ratio_norm/ratio_norm, 2))
            mc_mismodeling_ratio_tgraph.SetPoint(mismodeling_ratio_tgraph_index, STVal, mismodeling_ratio)
            mc_mismodeling_ratio_tgraph.SetPointError(mismodeling_ratio_tgraph_index, delta_STVal, mismodeling_ratio_error)
        mc_mismodeling_ratio_tgraph.GetXaxis().SetTitle(input_data_histogram.GetXaxis().GetTitle())
        mc_mismodeling_ratio_tgraph.GetXaxis().SetLimits(input_data_histogram.GetXaxis().GetXmin(), input_data_histogram.GetXaxis().GetXmax())
        mc_mismodeling_ratio_tgraph.GetYaxis().SetTitle("ratio of ratios")
        mc_mismodeling_ratio_tgraph.GetHistogram().SetMinimum(-0.5)
        mc_mismodeling_ratio_tgraph.GetHistogram().SetMaximum(3.5)
        mc_mismodeling_ratio_tgraph.Draw("AP0")
        ROOT.gPad.Update()
        # lineAt1 = ROOT.TLine(input_data_histogram.GetXaxis().GetXmin(), 1., input_data_histogram.GetXaxis().GetXmax(), 1.)
        # lineAt1.SetLineColor(ROOT.kBlack)
        # lineAt1.SetLineStyle(ROOT.kDashed)
        # lineAt1.Draw()
        # draw_horizontal_line_at_yval(yval=1.0, xmin=input_data_histogram.GetXaxis().GetXmin(), xmax=input_data_histogram.GetXaxis().GetXmax(), line_color=ROOT.kBlack, line_style=ROOT.kDashed)
        # fit_function_string = "[0] + [1]*((x/{c:.4f}) - 1.0)".format(c=input_data_histogram.GetXaxis().GetBinCenter(1))
        fit_function_string = "[0]"
        # print("Fitting function: {f}".format(f=fit_function_string))
        const_fit_tf1 = ROOT.TF1("fit_const_" + mc_mismodeling_ratio_tgraph.GetName(), fit_function_string, input_data_histogram.GetXaxis().GetBinLowEdge(1), input_data_histogram.GetXaxis().GetBinUpEdge(input_data_histogram.GetXaxis().GetNbins()))
        const_fit_tf1.SetParName(0, "const_{s}_{p}".format(s=selection, p=plot_to_extract))
        const_fit_tf1.SetParameter(0, 1.0)
        const_fit_tf1.SetParLimits(0, 0.0, 5.0)
        # const_fit_tf1.SetParName(1, "slope_{s}_{p}".format(s=selection, p=plot_to_extract))
        # const_fit_tf1.SetParameter(1, 0.0)
        # const_fit_tf1.SetParLimits(1, -5.0, 5.0)
        const_fit_result = mc_mismodeling_ratio_tgraph.Fit(const_fit_tf1, "QS0+")
        if not(const_fit_result.Status() == 0):
            print("Warning: fit failed with fit options \"QS0+\". Now trying options \"QEX0S0+\"...")
            const_fit_result = mc_mismodeling_ratio_tgraph.Fit(const_fit_tf1, "QEX0S0+")
            if not(const_fit_result.Status() == 0): sys.exit("ERROR: Unable to find fit for selection: {s}, plot_to_extract: {p}".format(s=selection, p=plot_to_extract))
        best_fit_const = const_fit_result.Value(0)
        best_fit_const_error = const_fit_result.ParError(0)
        # best_fit_slope = const_fit_result.Value(1)
        # best_fit_slope_error = const_fit_result.ParError(1)
        legend_mismodeling_ratio = ROOT.TLegend(0.1, 0.7, 0.9, 0.9)
        legend_mismodeling_ratio.SetFillStyle(0)
        const_fit_tf1.SetLineStyle(ROOT.kSolid)
        const_fit_tf1.SetLineColor(ROOT.kBlue)
        const_fit_tf1.Draw("LSAME")
        const_fit_tf1_clone = const_fit_tf1.Clone()
        const_fit_tf1_clone.SetName(const_fit_tf1.GetName() + "_clone")
        const_fit_tf1_clone.SetLineStyle(ROOT.kDashed)
        const_fit_tf1_clone.SetParameter(0, best_fit_const + 0.1)
        const_fit_tf1_clone.DrawCopy("LSAME")
        const_fit_tf1_clone.SetParameter(0, best_fit_const - 0.1)
        const_fit_tf1_clone.DrawCopy("LSAME")
        ROOT.gPad.Update()
        # legend_entry = legend_mismodeling_ratio.AddEntry(const_fit_tf1, "nominal fit: ({c:.3f} #pm {dc:.3f}) + ({s:.3f} #pm {ds:.3f})*((ST/{n:.2f}) - 1.0)".format(c=best_fit_const, dc=best_fit_const_error, s=best_fit_slope, ds=best_fit_slope_error, n=input_data_histogram.GetXaxis().GetBinCenter(1)))
        legend_entry = legend_mismodeling_ratio.AddEntry(const_fit_tf1, "nominal fit: ({c:.3f} #pm {dc:.3f})".format(c=best_fit_const, dc=best_fit_const_error))
        legend_entry.SetLineStyle(ROOT.kSolid)
        legend_entry.SetLineColor(ROOT.kBlue)
        legend_entry.SetMarkerColor(ROOT.kBlue)
        legend_entry.SetTextColor(ROOT.kBlue)
        legend_entry = legend_mismodeling_ratio.AddEntry(const_fit_tf1_clone, "nominal fit #pm 0.1")
        legend_entry.SetLineStyle(ROOT.kDashed)
        legend_entry.SetLineColor(ROOT.kBlue)
        legend_entry.SetMarkerColor(ROOT.kBlue)
        legend_entry.SetTextColor(ROOT.kBlue)
        legend_mismodeling_ratio.Draw()
        ROOT.gPad.Update()
        output_canvas_mismodeling_ratio.SaveAs("{o}/{s}_{n}JetsBin_mismodeling_ratio.pdf".format(o=output_folder, s=selection, n=nJetsBin))

    for process in (processes_BKG + ["data"]):
        source_file_objects[selection][process].Close()
print("-"*200)

# Step 6: plot ST distributions for the diphoton selections
print("-"*200)
print("Plotting ST distributions for diphoton selections...")
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
    ROOT.gStyle.SetOptStat(0)
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
print("-"*200)

# Print out normalizations to a LaTeX-formatted table, and store them in a dictionary
print("-"*200)
print("Storing normalizations to TeX-formatted table...")
norms_to_save = []
norms_tex_file_handle_file_path = "{o}/norm_values.tex".format(o=output_folder)
norms_tex_file_handle = open(norms_tex_file_handle_file_path, 'w')
norms_tex_file_handle.write("\\begin{tabular}{|p{0.25\\textwidth}|p{0.175\\textwidth}|p{0.175\\textwidth}|}\n")
norms_tex_file_handle.write("  \\hline\n")
norms_tex_file_handle.write("  nJets bin & HighHTQCD & GJetHT \\\\ \\hline\n")
for nJetsBin in range(2, 7):
    nJetsString = "{n}".format(n=nJetsBin)
    if (nJetsBin == 6): nJetsString = "$\\geq$ 6"
    norms_tex_file_handle.write("  {n} & {nQCD:.2f} & {nGJet:.2f} \\\\ \\hline\n".format(n=nJetsString, nQCD=K_fit["HighHTQCD"][nJetsBin], nGJet=K_fit["GJetHT"][nJetsBin], nDiph=1.0, nOverall=K_normSTBin[nJetsBin]))
    for process in processes_BKG:
        value_to_save = 1.0
        norm_value = K_normSTBin[nJetsBin]
        if (process in K_fit):
            norm_value *= K_fit[process][nJetsBin]
            value_to_save = K_fit[process][nJetsBin]
        # norms_to_save.append(tuple(["float", "norm_values_{p}_{n}JetsBin".format(p=process, n=nJetsBin), norm_value]))
        norms_to_save.append(tuple(["float", "norm_values_{p}_{n}JetsBin".format(p=process, n=nJetsBin), value_to_save]))
norms_tex_file_handle.write("\\end{tabular}\n")
norms_tex_file_handle.close()
print("Wrote normalization values to file: {o}".format(o=norms_tex_file_handle_file_path))
tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=norms_to_save, outputFilePath=("{o}/norm_values_nominal.dat".format(o=output_folder)))

print("-"*200)
print("K_fit:")
tmGeneralUtils.prettyPrintDictionary(K_fit)
print("K_normSTBin:")
tmGeneralUtils.prettyPrintDictionary(K_normSTBin)

print("-"*200)
print("All done!")
