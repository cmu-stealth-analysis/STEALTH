#!/usr/bin/env python

from __future__ import print_function, division
import os, sys, pdb, math, subprocess
import stealthEnv, ROOT

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

TOLERANCE=0.0001

output_folder = "/uscms/home/tmudholk/nobackup/analysisAreas/lowSTHistogramComparisons"
if (not(os.path.isdir(output_folder))): subprocess.check_call("mkdir -p {o}".format(o=output_folder), shell=True, executable="/bin/bash")

input_folder = "/uscms/home/tmudholk/nobackup/analysisAreas/lowSTHistograms"
dataset_names = ["EMEnrichedGJet", "Diphoton"]
dataset_denominator = "EMEnrichedGJet"
selections = ["signal", "signal_loose"]

histogram_names = {}
histogram_names["dist2D"] = "dist2D"
histogram_names["dist_nJets"] = "dist_nJets"

def get_ratio_with_error(numerator, numeratorErrorDown, numeratorErrorUp, denominator, denominatorErrorDown, denominatorErrorUp):
    if ((numerator <= TOLERANCE) or (denominator <= TOLERANCE)):
        return (0., 0., False)
    ratio = numerator/denominator
    fractionalError_numerator = (0.5*(numeratorErrorUp + numeratorErrorDown))/numerator
    fractionalError_denominator = (0.5*(denominatorErrorUp + denominatorErrorDown))/denominator
    fractionalError_ratio = math.sqrt(math.pow(fractionalError_numerator, 2) + math.pow(fractionalError_denominator, 2))
    return (ratio, ratio*fractionalError_ratio, True)

def save_2D_distribution(dist2D, title, name_to_save_as):
    output_canvas = ROOT.TCanvas("c_" + dist2D.GetName(), "c_" + dist2D.GetName(), 1200, 1200)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gPad.SetLogz()
    dist2D.SetTitle(title)
    # ROOT.gStyle.SetPaintTextFormat(".2f")
    ROOT.gPad.Update()
    ROOT.gPad.SetBottomMargin(0.1)
    ROOT.gPad.SetLeftMargin(0.1)
    ROOT.gPad.SetTopMargin(0.1)
    ROOT.gPad.SetRightMargin(0.175)
    ROOT.gPad.Update()
    dist2D.Draw("COLZ TEXT66")
    ROOT.gPad.Update()
    output_canvas.SaveAs("{o}/{n}.pdf".format(o=output_folder, n=name_to_save_as))

def save_th2D_ratio(th2d_numerator, th2d_denominator, title, name_to_save_as):
    hist_ratio = th2d_numerator.Clone()
    hist_ratio.SetName("ratio_" + th2d_numerator.GetName() + "_" + th2d_denominator.GetName())
    for xBinCounter in range(1, 1 + hist_ratio.GetXaxis().GetNbins()):
        for yBinCounter in range(1, 1 + hist_ratio.GetYaxis().GetNbins()):
            numerator, numeratorErrorDown, numeratorErrorUp = th2d_numerator.GetBinContent(xBinCounter, yBinCounter), th2d_numerator.GetBinErrorLow(xBinCounter, yBinCounter), th2d_numerator.GetBinErrorUp(xBinCounter, yBinCounter)
            denominator, denominatorErrorDown, denominatorErrorUp = th2d_denominator.GetBinContent(xBinCounter, yBinCounter), th2d_denominator.GetBinErrorLow(xBinCounter, yBinCounter), th2d_denominator.GetBinErrorUp(xBinCounter, yBinCounter)
            ratio, ratioError, isMeaningful = get_ratio_with_error(numerator, numeratorErrorDown, numeratorErrorUp, denominator, denominatorErrorDown, denominatorErrorUp)
            if isMeaningful:
                hist_ratio.SetBinContent(xBinCounter, yBinCounter, ratio)
                hist_ratio.SetBinError(xBinCounter, yBinCounter, ratioError)
            else:
                hist_ratio.SetBinContent(xBinCounter, yBinCounter, 0.)
                hist_ratio.SetBinError(xBinCounter, yBinCounter, 0.)

    output_canvas = ROOT.TCanvas("c_ratio_" + th2d_numerator.GetName() + "_" + th2d_denominator.GetName(), "c_ratio_" + th2d_numerator.GetName() + "_" + th2d_denominator.GetName(), 1200, 1200)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gPad.SetLogz()
    hist_ratio.SetTitle(title)
    hist_ratio.GetZaxis().SetTitle("")
    # ROOT.gStyle.SetPaintTextFormat(".3f")
    ROOT.gPad.Update()
    ROOT.gPad.SetBottomMargin(0.1)
    ROOT.gPad.SetLeftMargin(0.1)
    ROOT.gPad.SetTopMargin(0.1)
    ROOT.gPad.SetRightMargin(0.175)
    ROOT.gPad.Update()
    hist_ratio.Draw("COLZ TEXT66")
    ROOT.gPad.Update()
    output_canvas.SaveAs("{o}/{n}.pdf".format(o=output_folder, n=name_to_save_as))

def save_th1ds_with_ratio(th1d_numerator, numerator_name, th1d_denominator, denominator_name, title_overall, title_ratio, name_to_save_as):
    hist_ratio = th1d_numerator.Clone()
    output_canvas = ROOT.TCanvas("c_ratio_" + th1d_numerator.GetName() + "_" + th1d_denominator.GetName(), "c_ratio_" + th1d_numerator.GetName() + "_" + th1d_denominator.GetName(), 900, 1600)
    output_canvas.SetBottomMargin(0.1)
    output_canvas.SetLeftMargin(0.1)
    output_canvas.SetTopMargin(0.1)
    output_canvas.SetRightMargin(0.1)
    ROOT.gStyle.SetOptStat(0)
    output_canvas.Divide(1, 2, 0., 0.)
    output_canvas.cd(1)
    ROOT.gPad.SetLogy()
    ROOT.gPad.Update()
    legend = ROOT.TLegend(0.3, 0.9, 0.7, 1.0)
    legend.SetNColumns(2)
    th1d_numerator.SetTitle("")
    th1d_numerator.GetXaxis().SetTitle("")
    th1d_numerator.SetLineColor(ROOT.kBlue)
    th1d_numerator.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.Update()
    th1d_numerator.Draw("HIST PE0")
    ROOT.gPad.Update()
    legend_entry_numerator = legend.AddEntry(th1d_numerator, numerator_name)
    legend_entry_numerator.SetLineColor(ROOT.kBlue)
    legend_entry_numerator.SetTextColor(ROOT.kBlue)
    legend_entry_numerator.SetMarkerColor(ROOT.kBlue)
    th1d_denominator.SetLineColor(ROOT.kRed)
    th1d_denominator.SetMarkerColor(ROOT.kRed)
    th1d_denominator.Draw("HIST PE0 SAME")
    ROOT.gPad.Update()
    legend_entry_denominator = legend.AddEntry(th1d_denominator, denominator_name)
    legend_entry_denominator.SetLineColor(ROOT.kRed)
    legend_entry_denominator.SetTextColor(ROOT.kRed)
    legend_entry_denominator.SetMarkerColor(ROOT.kRed)
    ROOT.gPad.Update()
    legend.Draw()
    ROOT.gPad.Update()
    output_canvas.cd(2)
    ROOT.gPad.Update()
    hist_ratio.SetName("ratio_" + th1d_numerator.GetName() + "_" + th1d_denominator.GetName())
    for xBinCounter in range(1, 1 + hist_ratio.GetXaxis().GetNbins()):
        numerator, numeratorErrorDown, numeratorErrorUp = th1d_numerator.GetBinContent(xBinCounter), th1d_numerator.GetBinErrorLow(xBinCounter), th1d_numerator.GetBinErrorUp(xBinCounter)
        denominator, denominatorErrorDown, denominatorErrorUp = th1d_denominator.GetBinContent(xBinCounter), th1d_denominator.GetBinErrorLow(xBinCounter), th1d_denominator.GetBinErrorUp(xBinCounter)
        ratio, ratioError, isMeaningful = get_ratio_with_error(numerator, numeratorErrorDown, numeratorErrorUp, denominator, denominatorErrorDown, denominatorErrorUp)
        if isMeaningful:
            hist_ratio.SetBinContent(xBinCounter, ratio)
            hist_ratio.SetBinError(xBinCounter, ratioError)
        else:
            hist_ratio.SetBinContent(xBinCounter, 0.)
            hist_ratio.SetBinError(xBinCounter, 0.)

    ROOT.gStyle.SetOptStat(0)
    ROOT.gPad.SetLogy()
    hist_ratio.SetTitle("")
    hist_ratio.GetYaxis().SetTitle(title_ratio)
    # ROOT.gStyle.SetPaintTextFormat(".3f")
    ROOT.gPad.Update()
    # ROOT.gPad.SetBottomMargin(0.1)
    # ROOT.gPad.SetLeftMargin(0.1)
    # ROOT.gPad.SetTopMargin(0.1)
    # ROOT.gPad.SetRightMargin(0.175)
    # ROOT.gPad.Update()
    hist_ratio.Draw("")
    ROOT.gPad.Update()
    output_canvas.SaveAs("{o}/{n}.pdf".format(o=output_folder, n=name_to_save_as))

for selection in selections:
    for year_string in ["2016", "2017", "2018", "all"]:
        dists2D = {}
        distsNJets = {}
        for dataset_name in dataset_names:
            input_root_path = "{i}/{d}_{y}_{s}.root".format(i=input_folder, d=dataset_name, y=year_string, s=selection)
            inputFileHandle = ROOT.TFile.Open(input_root_path, "READ")
            if ((inputFileHandle.IsZombie() == ROOT.kTRUE) or not(inputFileHandle.IsOpen() == ROOT.kTRUE)): sys.exit("ERROR: Unable to open file at path \"{p}\"".format(p=input_root_path))
            dists2D[dataset_name] = ROOT.TH2D()
            inputFileHandle.GetObject(histogram_names["dist2D"], dists2D[dataset_name])
            save_2D_distribution(dists2D[dataset_name], "{n}, {y}, {s}".format(n=dataset_name, y=year_string, s=selection), "dist2D_{d}_{y}_{s}".format(d=dataset_name, y=year_string, s=selection))
            distsNJets[dataset_name] = ROOT.TH1D()
            inputFileHandle.GetObject(histogram_names["dist_nJets"], distsNJets[dataset_name])
            inputFileHandle.Close()
        for dataset_name in dataset_names:
            if (dataset_name == dataset_denominator): continue
            save_th2D_ratio(dists2D[dataset_name], dists2D[dataset_denominator], "{n} / {d}, {y}, {s}".format(n=dataset_name, d=dataset_denominator, y=year_string, s=selection), "ratio_2D_{n}_to_{d}_{y}_{s}".format(n=dataset_name, d=dataset_denominator, y=year_string, s=selection))
            save_th1ds_with_ratio(distsNJets[dataset_name], dataset_name, distsNJets[dataset_denominator], dataset_denominator, "NJets distributions with ratio: {n} / {d}, {y}, {s}".format(n=dataset_name, d=dataset_denominator, y=year_string, s=selection), "{n} / {d}".format(n=dataset_name, d=dataset_denominator), "dist_nJetsWithRatio_{n}_to_{d}_{y}_{s}".format(n=dataset_name, d=dataset_denominator, y=year_string, s=selection))
