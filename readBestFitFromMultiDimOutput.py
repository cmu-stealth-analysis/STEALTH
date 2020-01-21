#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, sys

inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination from input MC in real data.')
inputArgumentsParser.add_argument('--inputROOTFile', required=True, help='Name of input ROOT file containing observed and expected limits.',type=str)
inputArguments = inputArgumentsParser.parse_args()

inputFile=ROOT.TFile.Open("{fName}".format(fName=inputArguments.inputROOTFile), "READ")
if ((inputFile.IsZombie() == ROOT.kTRUE) or not(inputFile.IsOpen() == ROOT.kTRUE)):
    sys.exit("Error in opening file: {fName}".format(fName=inputArguments.inputROOTFile))
limitTree = ROOT.TTree()
inputFile.GetObject("limit", limitTree)
nEntriesFound = limitTree.GetEntries()
if not(nEntriesFound == 1): sys.exit("Error: multidim fit output not in expected format.")
limitTree.GetEntry(0)
print("{lTr:.6f}".format(lTr=limitTree.r))
inputFile.Close()
