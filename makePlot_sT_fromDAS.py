#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, ROOT, argparse
import numpy as np
from tmProgressBar import tmProgressBar

inputArgumentsParser = argparse.ArgumentParser(description='Run STEALTH selection.')
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--outputFilesSuffix', required=True, help='Prefix for output files.',type=str)
inputArguments = inputArgumentsParser.parse_args()

def make_ratio_graph(g_name, h_num, h_den):
    gae = ROOT.TGraphAsymmErrors()
    gae.SetName(g_name)
    h_rat = h_num.Clone()
    h_rat.Divide(h_den)
    for i in range(1, h_rat.GetNbinsX() + 1):
        n_rat = 0
        r_high = 0

        # tail = (1 - cl) / 2; for 95% CL, tail = (1 - 0.95) / 2 = 0.025
        tail = 0.16
        if h_num.GetBinError(i) == 0.0 or h_den.GetBinError(i) == 0.0:
            continue

        n_num = pow(h_num.GetBinContent(i) / h_num.GetBinError(i), 2)
        n_den = pow(h_den.GetBinContent(i) / h_den.GetBinError(i), 2)
        q_low = ROOT.Math.fdistribution_quantile_c(1 - tail, n_num * 2, (n_den + 1) * 2)
        r_low = q_low * n_num / (n_den + 1)
        if n_den > 0:
            n_rat = n_num / n_den
            q_high = ROOT.Math.fdistribution_quantile_c(tail, (n_num + 1) * 2, n_den * 2)
            r_high = q_high * (n_num + 1) / n_den
            gae.SetPoint(i-1,
                h_rat.GetBinCenter(i),
                h_rat.GetBinContent(i)
            )
            gae.SetPointError(i-1,
                h_rat.GetBinWidth(i) / 2,
                h_rat.GetBinWidth(i) / 2,
                n_rat - r_low,
                r_high - n_rat
            )
        # print ("i = " + str(i) + ": point = (" + str(h_rat.GetBinCenter(i)) + "," + str(h_rat.GetBinContent(i)) + "), point errors = (" + str(h_rat.GetBinWidth(i) / 2) + "," + str(h_rat.GetBinWidth(i) / 2) + "," + str(n_rat - r_low) + "," + str(r_high - n_rat) + ")")
    gae.GetXaxis().SetRangeUser(
        h_rat.GetBinLowEdge(1),
        h_rat.GetBinLowEdge(h_rat.GetNbinsX()) +
        h_rat.GetBinWidth(h_rat.GetNbinsX()))
    return gae

sw = ROOT.TStopwatch()
sw.Start()

chain_in = ROOT.TChain('ggNtuplizer/EventTree')
chain_in.Add(inputArguments.inputFilePath)
n_entries = chain_in.GetEntries()
print ('Total number of events: ' + str(n_entries))

# Initialize histograms
n_st_bins = 5
st_min = 900.0
st_norm = 1100.0
st_max = 3000.0
n_jets_min = 2
n_jets_max = 6
file_out = ROOT.TFile('analysis/hSTs_{outputFilesSuffix}.root'.format(outputFilesSuffix=inputArguments.outputFilesSuffix), 'recreate')
histograms = {}
for i in range(n_jets_min, n_jets_max + 1):
    hist_name = 'st_' + str(i) + 'Jets'
    histograms[hist_name] = ROOT.TH1F(
        'h_' + hist_name, ';S_{T} (GeV);AU', n_st_bins, st_min, st_max
        )
    histograms[hist_name].Sumw2()

# Fill histograms
progressBar = tmProgressBar(n_entries)
progressBarUpdatePeriod = n_entries//1000
progressBar.initializeTimer()
for j_entry in range(n_entries):
    i_entry = chain_in.LoadTree(j_entry)
    if i_entry < 0:
        break
    nb = chain_in.GetEntry(j_entry)
    if nb <= 0:
        continue

    if (j_entry%progressBarUpdatePeriod == 0): progressBar.updateBar(1.0*j_entry/n_entries, j_entry)

    n_stealth_jets = chain_in.b_nJets
    st = chain_in.b_evtST
    if n_stealth_jets >= n_jets_min and n_stealth_jets <= n_jets_max:
        histograms['st_' + str(n_stealth_jets) + 'Jets'].Fill(st)
    elif n_stealth_jets >= n_jets_min and n_stealth_jets > n_jets_max:
        histograms['st_' + str(n_jets_max) + 'Jets'].Fill(st)

print("\n") # proceed to next line after progress bar
# Scale histograms
for i in range(n_jets_min, n_jets_max + 1):
    hist_name = 'st_' + str(i) + 'Jets'
    norm_bin = histograms[hist_name].GetXaxis().FindBin(st_norm)
    max_bin = histograms[hist_name].GetNbinsX()
    norm = histograms[hist_name].Integral(norm_bin, max_bin)
    histograms[hist_name].Scale(1.0 / norm)

# Plot histograms
c_st = ROOT.TCanvas('c_st', '', 580, 620)
c_st.SetBorderSize(0)
c_st.SetFrameBorderMode(0)
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetOptStat(0)
p_top = ROOT.TPad('p_top', '', 0.005, 0.27, 0.995, 0.995)
p_bottom = ROOT.TPad('p_bottom', '', 0.005, 0.005, 0.995, 0.27)
p_top.Draw()
p_bottom.Draw()
p_top.SetMargin(12.0e-02, 3.0e-02, 5.0e-03, 2.0e-02)
p_bottom.SetMargin(12.0e-02, 3.0e-02, 29.0e-02, 4.0e-02)
p_top.SetTicky()
p_bottom.SetTicky()
p_bottom.SetGridy()
p_top.cd()
ROOT.gPad.SetLogy()

#Draw histograms
line_colors = [
    ROOT.kRed+1, ROOT.kGreen-3, ROOT.kAzure-2, ROOT.kOrange-3, ROOT.kGray
    ]
hist_name_0 = 'st_' + str(n_jets_min) + 'Jets'
histograms[hist_name_0].GetYaxis().SetRangeUser(2.0e-4, 9.0)
histograms[hist_name_0].GetXaxis().SetTitleOffset(1.1)
histograms[hist_name_0].SetLineColor(line_colors[0])
histograms[hist_name_0].Draw('E')
for i in range(n_jets_min + 1, n_jets_max + 1):
    hist_name = 'st_' + str(i) + 'Jets'
    histograms[hist_name].SetLineColor(line_colors[i - n_jets_min])
    histograms[hist_name].Draw('SAME E')
    c_st.Update()

#Draw legend and labels
legend = ROOT.TLegend(0.75, 0.65, 0.86, 0.92)
for i in range(n_jets_min, n_jets_max):
    hist_name = 'st_' + str(i) + 'Jets'
    legend.AddEntry(histograms[hist_name], str(i) + ' jets')
legend.AddEntry(histograms['st_' + str(n_jets_max) + 'Jets'],
                '#geq ' + str(n_jets_max) + ' jets')
legend.SetBorderSize(0)
legend.Draw()

# Draw ratios
p_bottom.cd()
f_unity = ROOT.TF1('f_unity', '[0]', st_min, st_max)
f_unity.SetTitle("")
f_unity.SetParameter(0, 1.0)
f_unity.GetXaxis().SetLabelSize(0.1)
f_unity.GetXaxis().SetTickLength(0.1)
f_unity.GetXaxis().SetTitle('S_{T} (GeV)')
f_unity.GetXaxis().SetTitleOffset(1.14)
f_unity.GetXaxis().SetTitleSize(0.12)
f_unity.GetYaxis().SetLabelSize(0.1)
f_unity.GetYaxis().SetNdivisions(305)
f_unity.GetYaxis().SetRangeUser(0.25, 4.0)
f_unity.GetYaxis().SetTickLength(0.04)
f_unity.GetYaxis().SetTitle('n_{j} / n_{j} = ' + str(n_jets_min))
f_unity.GetYaxis().SetTitleSize(0.11)
f_unity.GetYaxis().SetTitleOffset(0.5)
f_unity.SetLineStyle(4)
f_unity.Draw()
graphsToPlot = {}
for i in range(n_jets_min + 1, n_jets_max + 1):
    hist_name = 'st_' + str(i) + 'Jets'
    graphsToPlot[i] = make_ratio_graph('gae_' + str(i) + 'Jets', histograms[hist_name], histograms[hist_name_0])
    graphsToPlot[i].SetMarkerStyle(20)
    graphsToPlot[i].SetMarkerColor(line_colors[i - n_jets_min])
    graphsToPlot[i].SetLineColor(line_colors[i - n_jets_min])
    graphsToPlot[i].SetMarkerSize(0.9)
    graphsToPlot[i].Draw('P SAME')
    c_st.Update()
c_st.SaveAs("analysis/plot_st_{outputFilesSuffix}.png".format(outputFilesSuffix=inputArguments.outputFilesSuffix))
c_st.Write()
file_out.Write()
file_out.Close()

sw.Stop()
print ('Real time: ' + str(sw.RealTime() / 60.0) + ' minutes')
print ('CPU time:  ' + str(sw.CpuTime() / 60.0) + ' minutes')
