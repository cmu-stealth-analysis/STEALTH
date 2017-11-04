#!/usr/bin/env python

from __future__ import print_function, division

import glob, argparse, ROOT

inputArgumentsParser = argparse.ArgumentParser(description='Print total number of events in a file pattern to a file.')
inputArgumentsParser.add_argument('--escapedInputFilePattern', required=True, help='Escaped glob pattern to select input files IN ARBITRARY ORDER.',type=str)  # WARNING!!! WARNING!!! glob.glob returns list of files matching wildcard expansion IN ARBITRARY ORDER, do NOT mess around with order of events!
inputArgumentsParser.add_argument('--outputFilePath', required=True, help='Path to output file.',type=str)
inputArguments = inputArgumentsParser.parse_args()

ggIn = ROOT.TChain("ggNtuplizer/EventTree")

listOfInputFiles = glob.glob(inputArguments.escapedInputFilePattern) # WARNING!!! WARNING!!! glob.glob returns list of files matching wildcard expansion IN ARBITRARY ORDER, do NOT mess around with order of events!

for inputFile in listOfInputFiles:
    print ("Adding... " + inputFile)
    ggIn.Add(inputFile)
nEvts = ggIn.GetEntries()
print(" >> total nEvts:" + str(nEvts))
if (nEvts == 0): sys.exit("No events found!")
