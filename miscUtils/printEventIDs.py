#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, subprocess, argparse

inputArgumentsParser = argparse.ArgumentParser(description='Print event IDs for all events in a given n-tuple.')
inputArgumentsParser.add_argument('--inputPath_with_prefix', required=True, action='append', help='Path to sample.', type=str)
inputArguments = inputArgumentsParser.parse_args()

import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

inputChain = ROOT.TChain("ggNtuplizer/EventTree")
inputChain.SetMaxTreeSize(100000000000) # 1 TB
for inputPath in inputArguments.inputPath_with_prefix:
    inputChain.Add(inputPath)
nEntries = inputChain.GetEntries()

for eventIndex in range(0, nEntries):
    treeStatus = inputChain.LoadTree(eventIndex)
    if (treeStatus < 0):
        break
    evtStatus = inputChain.GetEntry(eventIndex)
    if (evtStatus <= 0):
        continue
    # print("{r}:{l}:{e}-{r}:{l}:{e}".format(r=inputChain.run, l=inputChain.lumis, e=inputChain.event))
    index_leading = inputChain.b_photonIndex_leading
    index_subLeading = inputChain.b_photonIndex_subLeading
    print("{r} {l} {e} {ptl:.3f} {etal:.3f} {phil:.3f} {ptsl:.3f} {etasl:.3f} {phisl:.3f}".format(r=inputChain.run, l=inputChain.lumis, e=inputChain.event,
                                                                                                  ptl=inputChain.phoEt[index_leading], etal=inputChain.phoEta[index_leading], phil=inputChain.phoPhi[index_leading],
                                                                                                  ptsl=inputChain.phoEt[index_subLeading], etasl=inputChain.phoEta[index_subLeading], phisl=inputChain.phoPhi[index_subLeading]))
# ./miscUtils/printEventIDs.py --inputPath_with_prefix ${EOSPREFIX}/store/user/lpcsusystealth/selections/combined_DoublePhoton/merged_selection_MC_DiPhotonJets_2017_signal.root > ~/nobackup/cmssw/minimal_analyzer/eventLists/events_DiPhotonJets17_signal.txt
# ./miscUtils/printEventIDs.py --inputPath_with_prefix ${EOSPREFIX}/store/user/lpcsusystealth/selections/combined_DoublePhoton/merged_selection_MC_HighHTQCD17_2017_signal.root > ~/nobackup/cmssw/minimal_analyzer/eventLists/events_HighHTQCD17_signal.txt
# ./miscUtils/printEventIDs.py --inputPath_with_prefix ${EOSPREFIX}/store/user/lpcsusystealth/selections/combined_DoublePhoton/merged_selection_MC_GJetHT17_2017_signal.root > ~/nobackup/cmssw/minimal_analyzer/eventLists/events_GJetHT17_signal.txt
