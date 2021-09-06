#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, subprocess, argparse
import tmProgressBar
import stealthEnv

inputArgumentsParser = argparse.ArgumentParser(description='Print weighted number of events from MC sample with weights.')
inputArgumentsParser.add_argument('--inputPath_with_prefix', required=True, action='append', help='Path to MC sample.', type=str)
inputArgumentsParser.add_argument('--addMCWeights', action='store_true', help="If this flag is set, then MC event weights are read in as well.")
inputArguments = inputArgumentsParser.parse_args()

import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

inputChain = ROOT.TChain("ggNtuplizer/EventTree")
inputChain.SetMaxTreeSize(100000000000) # 1 TB
for inputPath in inputArguments.inputPath_with_prefix:
    inputChain.Add(inputPath)
nEntries = inputChain.GetEntries()
print("Available nEvts: {n}".format(n=nEntries))

nEventsDict = {}
for nJetsBin in range(2, 7):
    nEventsDict[nJetsBin] = {
        "weighted": {
            "norm": 0.0,
            "obs": 0.0
        },
        "unweighted": {
            "norm": 0,
            "obs": 0
        }
    }

progressBar = tmProgressBar.tmProgressBar(nEntries)
progressBarUpdatePeriod = max(1, nEntries//50)
progressBar.initializeTimer()
for eventIndex in range(0, nEntries):
    treeStatus = inputChain.LoadTree(eventIndex)
    if (treeStatus < 0):
        break
    evtStatus = inputChain.GetEntry(eventIndex)
    if (evtStatus <= 0):
        continue
    if (eventIndex % progressBarUpdatePeriod == 0): progressBar.updateBar(eventIndex/nEntries, eventIndex)
    ST = inputChain.b_evtST
    ST_region = None
    if ((ST >= 1200.0) and (ST < 1300.0)):
        ST_region = "norm"
    elif (ST >= 1300.0):
        ST_region = "obs"
    else:
        continue
    nJetsDR = inputChain.b_nJetsDR
    nJetsBin = min(nJetsDR, 6)
    if (nJetsBin < 2): continue
    eventWeight = 1.0
    if (inputArguments.addMCWeights):
        eventWeight = inputChain.b_MCCustomWeight
    nEventsDict[nJetsBin]["weighted"][ST_region] += eventWeight
    nEventsDict[nJetsBin]["unweighted"][ST_region] += 1

for nJetsBin in range(2, 7):
    print("In {n} jets bin, nEvents_norm (unweighted): {nnormu}, nEvents_obs (unweighted) = {nobsu}, nEvents_norm (weighted): {nnormw}, nEvents_obs (weighted) = {nobsw}".format(n=nJetsBin, nnormu=nEventsDict[nJetsBin]["unweighted"]["norm"], nobsu=nEventsDict[nJetsBin]["unweighted"]["obs"], nnormw=nEventsDict[nJetsBin]["weighted"]["norm"], nobsw=nEventsDict[nJetsBin]["weighted"]["obs"]))

print("All done!")
