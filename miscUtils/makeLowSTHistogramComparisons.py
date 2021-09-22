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
dataset_names = ["DiPhotonJets", "EMEnrichedGJetPt", "HighHTQCD"]
dataset_colors = {
    "DiPhotonJets": ROOT.kBlack,
    "EMEnrichedGJetPt": ROOT.kBlue,
    "HighHTQCD": ROOT.kRed
}
dataset_denominator = "DiPhotonJets"
selections = ["signal", "signal_loose"]
years_for_nJets_ratios = ["all"]

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

def save_th1ds_with_ratio(numerators_th1_details, denominator_th1_details, title_overall, title_ratio, name_to_save_as):
    th1d_denominator = denominator_th1_details["th1d"]
    output_canvas = ROOT.TCanvas("c_{n}".format(n=name_to_save_as), "c_{n}".format(n=name_to_save_as), 1600, 1600)
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
    legend.SetNColumns(1+len(numerators_th1_details))
    th1d_denominator.SetTitle("")
    th1d_denominator.GetXaxis().SetTitle("")
    th1d_denominator.SetLineColor(denominator_th1_details["color"])
    th1d_denominator.SetMarkerColor(denominator_th1_details["color"])
    ROOT.gPad.Update()
    th1d_denominator.Draw("HIST E0 P")
    ROOT.gPad.Update()
    legend_entry_denominator = legend.AddEntry(th1d_denominator, denominator_th1_details["name"])
    legend_entry_denominator.SetLineColor(denominator_th1_details["color"])
    legend_entry_denominator.SetTextColor(denominator_th1_details["color"])
    legend_entry_denominator.SetMarkerColor(denominator_th1_details["color"])
    for numerator_th1_details in numerators_th1_details:
        th1d_numerator = numerator_th1_details["th1d"]
        th1d_numerator.SetLineColor(numerator_th1_details["color"])
        th1d_numerator.SetMarkerColor(numerator_th1_details["color"])
        th1d_numerator.Draw("HIST E0 P SAME")
        ROOT.gPad.Update()
        legend_entry_numerator = legend.AddEntry(th1d_numerator, numerator_th1_details["name"])
        legend_entry_numerator.SetLineColor(numerator_th1_details["color"])
        legend_entry_numerator.SetTextColor(numerator_th1_details["color"])
        legend_entry_numerator.SetMarkerColor(numerator_th1_details["color"])
        ROOT.gPad.Update()
    legend.Draw()
    ROOT.gPad.Update()
    output_canvas.cd(2)
    ROOT.gPad.Update()
    ROOT.gPad.SetLogy(0)
    ROOT.gPad.Update()
    is_first_ratio = True
    hist_ratios = {}
    for numerator_th1_details in numerators_th1_details:
        th1d_numerator = numerator_th1_details["th1d"]
        hist_ratios[numerator_th1_details["name"]] = ROOT.TH1D("ratio_" + th1d_numerator.GetName() + "_" + th1d_denominator.GetName(), "", th1d_denominator.GetXaxis().GetNbins(), th1d_denominator.GetXaxis().GetBinLowEdge(1), th1d_denominator.GetXaxis().GetBinUpEdge(th1d_denominator.GetXaxis().GetNbins()))
        hist_ratios[numerator_th1_details["name"]].GetYaxis().SetRangeUser(-0.5, 3.5)
        ROOT.gPad.Update()
        for xBinCounter in range(1, 1 + hist_ratios[numerator_th1_details["name"]].GetXaxis().GetNbins()):
            numerator, numeratorErrorDown, numeratorErrorUp = th1d_numerator.GetBinContent(xBinCounter), th1d_numerator.GetBinErrorLow(xBinCounter), th1d_numerator.GetBinErrorUp(xBinCounter)
            denominator, denominatorErrorDown, denominatorErrorUp = th1d_denominator.GetBinContent(xBinCounter), th1d_denominator.GetBinErrorLow(xBinCounter), th1d_denominator.GetBinErrorUp(xBinCounter)
            ratio, ratioError, isMeaningful = get_ratio_with_error(numerator, numeratorErrorDown, numeratorErrorUp, denominator, denominatorErrorDown, denominatorErrorUp)
            if isMeaningful:
                hist_ratios[numerator_th1_details["name"]].SetBinContent(xBinCounter, ratio)
                hist_ratios[numerator_th1_details["name"]].SetBinError(xBinCounter, ratioError)
            else:
                hist_ratios[numerator_th1_details["name"]].SetBinContent(xBinCounter, 0.)
                hist_ratios[numerator_th1_details["name"]].SetBinError(xBinCounter, 0.)
        hist_ratios[numerator_th1_details["name"]].SetLineColor(numerator_th1_details["color"])
        hist_ratios[numerator_th1_details["name"]].SetMarkerColor(numerator_th1_details["color"])
        ROOT.gPad.Update()
        if is_first_ratio:
            is_first_ratio = False
            hist_ratios[numerator_th1_details["name"]].SetTitle("")
            hist_ratios[numerator_th1_details["name"]].GetXaxis().SetTitle("nJets")
            hist_ratios[numerator_th1_details["name"]].GetYaxis().SetTitle(title_ratio)
            hist_ratios[numerator_th1_details["name"]].Draw("HIST E0 P")
        else:
            hist_ratios[numerator_th1_details["name"]].Draw("HIST E0 P SAME")
        ROOT.gPad.Update()
    ROOT.gPad.RedrawAxis()
    lineObject = ROOT.TLine()
    lineObject.SetLineColor(ROOT.kBlack)
    lineObject.SetLineStyle(ROOT.kDashed)
    lineObject.DrawLine(th1d_denominator.GetXaxis().GetBinLowEdge(1), 1.0, th1d_denominator.GetXaxis().GetBinUpEdge(th1d_denominator.GetXaxis().GetNbins()), 1.0)
    ROOT.gPad.Update()
    output_canvas.SaveAs("{o}/{n}.pdf".format(o=output_folder, n=name_to_save_as))

for selection in selections:
    for year_string in ["2016", "2017", "2018", "all"]:
        denominator_th1_details = None
        numerators_th1_details = []
        dists2D = {}
        # distsNJets = {}
        for dataset_name in dataset_names:
            input_root_path = "{i}/{d}_{y}_{s}.root".format(i=input_folder, d=dataset_name, y=year_string, s=selection)
            inputFileHandle = ROOT.TFile.Open(input_root_path, "READ")
            if ((inputFileHandle.IsZombie() == ROOT.kTRUE) or not(inputFileHandle.IsOpen() == ROOT.kTRUE)): sys.exit("ERROR: Unable to open file at path \"{p}\"".format(p=input_root_path))
            dists2D[dataset_name] = ROOT.TH2D()
            inputFileHandle.GetObject(histogram_names["dist2D"], dists2D[dataset_name])
            save_2D_distribution(dists2D[dataset_name], "{n}, {y}, {s}".format(n=dataset_name, y=year_string, s=selection), "dist2D_{d}_{y}_{s}".format(d=dataset_name, y=year_string, s=selection))
            th1_details = {
                "th1d": ROOT.TH1D(),
                "color": dataset_colors[dataset_name],
                "name": dataset_name
            }
            inputFileHandle.GetObject(histogram_names["dist_nJets"], th1_details["th1d"])
            if (dataset_name == dataset_denominator):
                denominator_th1_details = th1_details
            else:
                numerators_th1_details.append(th1_details)
                save_th2D_ratio(dists2D[dataset_name], dists2D[dataset_denominator], "{n} / {d}, {y}, {s}".format(n=dataset_name, d=dataset_denominator, y=year_string, s=selection), "ratio_2D_{n}_to_{d}_{y}_{s}".format(n=dataset_name, d=dataset_denominator, y=year_string, s=selection))
            inputFileHandle.Close()
        if not(year_string in years_for_nJets_ratios): continue
        save_th1ds_with_ratio(numerators_th1_details, denominator_th1_details, "NJets distributions, {y}, {s}".format(y=year_string, s=selection), "Ratios w.r.t. {d}".format(d=dataset_denominator), "dist_nJetsWithRatio_{y}_{s}".format(y=year_string, s=selection))
