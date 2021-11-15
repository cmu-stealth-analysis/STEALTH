#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, subprocess, argparse

inputArgumentsParser = argparse.ArgumentParser(description='Print weighted number of events from MC sample with weights.')
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
    print("{r},{e}".format(r=inputChain.run, e=inputChain.event))
