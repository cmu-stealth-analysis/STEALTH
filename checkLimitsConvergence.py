#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, sys

inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination from input MC in real data.')
inputArgumentsParser.add_argument('--inputROOTFile', required=True, help='Name of input ROOT file containing observed and expected limits.',type=str)
inputArgumentsParser.add_argument('--maxAllowedRatio', default=10., help='Max allowed ratio for deviation between expected and observed limits.',type=float)
inputArguments = inputArgumentsParser.parse_args()

def passesSanityChecks(observedUpperLimit, expectedUpperLimit):
    ratio = observedUpperLimit/expectedUpperLimit
    if ((ratio > inputArguments.maxAllowedRatio) or (ratio < (1.0/inputArguments.maxAllowedRatio))): return False
    return True

inputFile=ROOT.TFile.Open("{fName}".format(fName=inputArguments.inputROOTFile), "READ")
if ((inputFile.IsZombie() == ROOT.kTRUE) or not(inputFile.IsOpen() == ROOT.kTRUE)):
    sys.exit("Error in opening file: {fName}".format(fName=inputArguments.inputROOTFile))
limitTree = ROOT.TTree()
inputFile.GetObject("limit", limitTree)
nEntriesFound = limitTree.GetEntries()
if not(nEntriesFound == 6): sys.exit("ERROR: limits not in proper format.")
limitTree.GetEntry(2)
expectedUpperLimit = limitTree.limit
limitTree.GetEntry(1)
expectedUpperLimitOneSigmaDown=limitTree.limit
limitTree.GetEntry(3)
expectedUpperLimitOneSigmaUp=limitTree.limit
limitTree.GetEntry(5)
observedUpperLimit = limitTree.limit
inputFile.Close()

if (passesSanityChecks(observedUpperLimit, expectedUpperLimit)): print("true")
else: print("false")
