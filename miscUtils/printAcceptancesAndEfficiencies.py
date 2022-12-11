#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, subprocess, argparse

import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

INPUT_FILE_PATH = 'root://cmseos.fnal.gov//store/user/lpcsusystealth/statistics/combined_DoublePhoton/merged_statistics_MC_stealth_t5_2017.root'
input_file_handle = ROOT.TFile.Open(INPUT_FILE_PATH, 'READ')
if ((input_file_handle.IsZombie() == ROOT.kTRUE) or not(input_file_handle.IsOpen() == ROOT.kTRUE)): sys.exit('ERROR: Unable to open file {f}'.format(f=INPUT_FILE_PATH))

for mass_bin_descr in ["eventProgenitor1950_neutralino1850", "eventProgenitor2100_neutralino1000", "eventProgenitor1850_neutralino200"]:
    acceptances = ROOT.TEfficiency()
    label = "Acceptance_6Jets_all_MC_" + mass_bin_descr
    acceptances.SetName(label)
    input_file_handle.GetObject(label, acceptances)
    acc = acceptances.GetEfficiency(acceptances.FindFixBin(3000.))
    acc_error_up = acceptances.GetEfficiencyErrorUp(acceptances.FindFixBin(3000.))
    acc_error_lo = acceptances.GetEfficiencyErrorLow(acceptances.FindFixBin(3000.))
    efficiencies = ROOT.TEfficiency()
    label = "IDEfficiency_6Jets_all_MC_" + mass_bin_descr
    efficiencies.SetName(label)
    input_file_handle.GetObject(label, efficiencies)
    eff = efficiencies.GetEfficiency(efficiencies.FindFixBin(3000.))
    eff_error_up = efficiencies.GetEfficiencyErrorUp(efficiencies.FindFixBin(3000.))
    eff_error_lo = efficiencies.GetEfficiencyErrorLow(efficiencies.FindFixBin(3000.))
    print("{b}: acceptance = {acc:.3f} (+{acc_error_up:.3f}/-{acc_error_lo:.3f}), efficiency = {eff:.3f} (+{eff_error_up:.3f}/-{eff_error_lo:.3f})".format(b=mass_bin_descr.replace("eventProgenitor", "gluino"), acc=acc, acc_error_up=acc_error_up, acc_error_lo=acc_error_lo, eff=eff, eff_error_up=eff_error_up, eff_error_lo=eff_error_lo))

input_file_handle.Close()
print("All done!")
