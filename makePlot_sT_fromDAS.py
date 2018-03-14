#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, ROOT, argparse
import numpy as np
from tmProgressBar import tmProgressBar

inputArgumentsParser = argparse.ArgumentParser(description='Make basic histograms of sT and nJets.')
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--nSTBins_sTHistograms', default=10, help='Number of sT bins in sT histograms.',type=int)
inputArgumentsParser.add_argument('--sTMin_sTHistograms', default=500., help='Min value of sT.',type=float)
inputArgumentsParser.add_argument('--sTMax_sTHistograms', default=2500., help='Max value of sT.',type=float)
inputArgumentsParser.add_argument('--sTNormRangeMin_sTHistograms', default=1000., help='Min value of sT for normalization.',type=float)
inputArgumentsParser.add_argument('--sTNormRangeMax_sTHistograms', default=2500., help='Max value of sT for normalization.',type=float)
inputArgumentsParser.add_argument('--nJetsMin_sTHistograms', default=2, help='Min number of nJets for sT histograms.',type=int)
inputArgumentsParser.add_argument('--nJetsMax_sTHistograms', default=6, help='Least value of nJets in highest nJets bin.',type=int)
inputArgumentsParser.add_argument('--nJetsNorm_sTHistograms', default=2, help='Number of jets w.r.t. which to normalize the sT distributions for other jets.',type=int)
inputArgumentsParser.add_argument('--sTMin_nJetsHistograms', default=1100., help='Min value of sT to include in nJets histograms.',type=float)
inputArgumentsParser.add_argument('--nJetsMin_nJetsHistograms', default=2, help='Min number of jets in nJets histogram.',type=int)
inputArgumentsParser.add_argument('--nJetsMax_nJetsHistograms', default=10, help='Max number of jets in nJets histogram.',type=int)
inputArgumentsParser.add_argument('--outputFilesSuffix', required=True, help='Prefix for output files.',type=str)
inputArguments = inputArgumentsParser.parse_args()

lineColors = {
    2: ROOT.kBlack,
    3: ROOT.kRed+1,
    4: ROOT.kGreen-3,
    5: ROOT.kAzure-2,
    6: ROOT.kOrange-3
}

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
            # print ("i = " + str(i) + ": point = (" + str(h_rat.GetBinCenter(i)) + "," + str(h_rat.GetBinContent(i)) + "), point errors = (" + str(h_rat.GetBinWidth(i) / 2) + "," + str(h_rat.GetBinWidth(i) / 2) + "," + str(n_rat
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
file_out = ROOT.TFile('analysis/basicHistograms/hSTs_{outputFilesSuffix}.root'.format(outputFilesSuffix=inputArguments.outputFilesSuffix), 'recreate')
histograms = {
    "sT": {},
    "nJets": None
}
for nJetsBin in range(inputArguments.nJetsMin_sTHistograms, inputArguments.nJetsMax_sTHistograms + 1):
    hist_name = 'st_' + str(nJetsBin) + 'Jets'
    histograms["sT"][nJetsBin] = ROOT.TH1F('h_' + hist_name, ';S_{T} (GeV);AU', inputArguments.nSTBins_sTHistograms, inputArguments.sTMin_sTHistograms, inputArguments.sTMax_sTHistograms)
    histograms["sT"][nJetsBin].Sumw2()
    hist_name = 'nJets_' + str(nJetsBin)
    histograms["nJets"] = ROOT.TH1F('h_' + hist_name, ';nJets;nEvents', 1 + inputArguments.nJetsMax_nJetsHistograms - inputArguments.nJetsMin_nJetsHistograms, -0.5 + inputArguments.nJetsMin_nJetsHistograms, 0.5 + inputArguments.nJetsMax_nJetsHistograms)
    histograms["nJets"].Sumw2()

# Fill histograms
progressBar = tmProgressBar(n_entries)
progressBarUpdatePeriod = max(1, n_entries//1000)
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
    histograms["nJets"].Fill(n_stealth_jets)
    nJetsBin = n_stealth_jets
    if (nJetsBin > inputArguments.nJetsMax_sTHistograms): nJetsBin = inputArguments.nJetsMax_sTHistograms
    st = chain_in.b_evtST
    histograms["sT"][nJetsBin].Fill(st)

progressBar.terminate()
# Scale histograms
for nJetsBin in range(inputArguments.nJetsMin_sTHistograms, inputArguments.nJetsMax_sTHistograms + 1):
    norm_bin_start = histograms["sT"][nJetsBin].GetXaxis().FindBin(inputArguments.sTNormRangeMin_sTHistograms)
    norm_bin_end = histograms["sT"][nJetsBin].GetXaxis().FindBin(inputArguments.sTNormRangeMax_sTHistograms)
    normValue = histograms["sT"][nJetsBin].Integral(norm_bin_start, norm_bin_end)
    histograms["sT"][nJetsBin].Scale(1.0 / normValue)

# Plot histograms
c_st = ROOT.TCanvas('c_st', 'c_st', 1024, 768)
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

# histograms["sT"][inputArguments.nJetsMin_sTHistograms].GetYaxis().SetRangeUser(2.0e-4, 9.0)
histograms["sT"][inputArguments.nJetsMin_sTHistograms].GetXaxis().SetTitleOffset(1.1)
histograms["sT"][inputArguments.nJetsMin_sTHistograms].SetLineColor(lineColors[inputArguments.nJetsMin_sTHistograms])
histograms["sT"][inputArguments.nJetsMin_sTHistograms].Draw('E')
for nJetsBin in range(inputArguments.nJetsMin_sTHistograms + 1, inputArguments.nJetsMax_sTHistograms + 1):
    histograms["sT"][nJetsBin].SetLineColor(lineColors[nJetsBin])
    histograms["sT"][nJetsBin].Draw('SAME E')
    c_st.Update()

#Draw legend and labels
legend = ROOT.TLegend(0.75, 0.65, 0.86, 0.92)
for nJetsBin in range(inputArguments.nJetsMin_sTHistograms, inputArguments.nJetsMax_sTHistograms):
    legend.AddEntry(histograms["sT"][nJetsBin], str(nJetsBin) + ' jets')
legend.AddEntry(histograms["sT"][inputArguments.nJetsMax_sTHistograms], '#geq ' + str(inputArguments.nJetsMax_sTHistograms) + ' jets')
legend.SetBorderSize(0)
legend.Draw()
c_st.Update()

# Draw ratios
p_bottom.cd()
f_unity = ROOT.TF1('f_unity', '[0]', inputArguments.sTMin_sTHistograms, inputArguments.sTMax_sTHistograms)
f_unity.SetTitle("")
f_unity.SetParameter(0, 1.0)
f_unity.GetXaxis().SetLabelSize(0.1)
f_unity.GetXaxis().SetTickLength(0.1)
f_unity.GetXaxis().SetTitle('S_{T} (GeV)')
f_unity.GetXaxis().SetTitleOffset(1.14)
f_unity.GetXaxis().SetTitleSize(0.12)
f_unity.GetYaxis().SetLabelSize(0.1)
f_unity.GetYaxis().SetNdivisions(603)
f_unity.GetYaxis().SetRangeUser(0.0, 2.0)
f_unity.GetYaxis().SetTickLength(0.04)
f_unity.GetYaxis().SetTitle('n_{j} / n_{j} = ' + str(inputArguments.nJetsNorm_sTHistograms))
f_unity.GetYaxis().SetTitleSize(0.11)
f_unity.GetYaxis().SetTitleOffset(0.2)
f_unity.SetLineStyle(4)
f_unity.SetLineColor(lineColors[inputArguments.nJetsNorm_sTHistograms])
f_unity.Draw()
graphsToPlot = {}
for nJetsBin in range(inputArguments.nJetsMin_sTHistograms, 1 + inputArguments.nJetsMax_sTHistograms):
    if (nJetsBin == inputArguments.nJetsNorm_sTHistograms): continue
    graphsToPlot[nJetsBin] = make_ratio_graph('gae_' + str(nJetsBin) + 'Jets', histograms["sT"][nJetsBin], histograms["sT"][inputArguments.nJetsNorm_sTHistograms])
    graphsToPlot[nJetsBin].SetMarkerStyle(20)
    graphsToPlot[nJetsBin].SetMarkerColor(lineColors[nJetsBin])
    graphsToPlot[nJetsBin].SetLineColor(lineColors[nJetsBin])
    graphsToPlot[nJetsBin].SetMarkerSize(0.9)
    graphsToPlot[nJetsBin].Draw('P SAME')
    c_st.Update()
c_st.Print("analysis/basicHistograms/histograms_sT_{outputFilesSuffix}.png".format(outputFilesSuffix=inputArguments.outputFilesSuffix))
c_st.Write()

c_nJets = c_st = ROOT.TCanvas('c_nJets', 'c_nJets', 1024, 768)
c_nJets.SetBorderSize(0)
c_nJets.SetFrameBorderMode(0)
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gPad.SetLogy()
histograms["nJets"].Draw('')
c_nJets.Print("analysis/basicHistograms/histograms_nJets_{outputFilesSuffix}.png".format(outputFilesSuffix=inputArguments.outputFilesSuffix))

file_out.Write()
file_out.Close()

sw.Stop()
print ('Real time: ' + str(sw.RealTime() / 60.0) + ' minutes')
print ('CPU time:  ' + str(sw.CpuTime() / 60.0) + ' minutes')
