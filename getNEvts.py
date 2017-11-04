#!/usr/bin/env python

from __future__ import print_function, division

import glob, argparse, ROOT, sys

inputArgumentsParser = argparse.ArgumentParser(description='Print total number of events in a file pattern to a file.')
inputArgumentsParser.add_argument('--escapedInputFilePattern', required=True, help='Escaped glob pattern to select input files IN ARBITRARY ORDER.',type=str)  # WARNING!!! WARNING!!! glob.glob returns list of files matching wildcard expansion IN ARBITRARY ORDER, do NOT mess around with order of events!
inputArgumentsParser.add_argument('--noGlob', action='store_true', help='If set, the glob expansion is disabled and the input file pattern is parsed as a single file.')  # WARNING!!! WARNING!!! glob.glob returns list of files matching wildcard expansion IN ARBITRARY ORDER, do NOT mess around with order of events!
inputArgumentsParser.add_argument('--outputFilePath', required=True, help='Path to output file.',type=str)
inputArguments = inputArgumentsParser.parse_args()

ggIn = ROOT.TChain("ggNtuplizer/EventTree")

if (inputArguments.noGlob):
    listOfInputFiles = [inputArguments.escapedInputFilePattern]
else:
    listOfInputFiles = glob.glob(inputArguments.escapedInputFilePattern) # WARNING!!! WARNING!!! glob.glob returns list of files matching wildcard expansion IN ARBITRARY ORDER, do NOT mess around with order of events!

for inputFile in listOfInputFiles:
    print ("Adding... " + inputFile)
    ggIn.Add(inputFile)
nEvts = ggIn.GetEntries()

if (nEvts == 0): sys.exit("No events found!")
outputFile = open(inputArguments.outputFilePath, 'w')
outputFile.write("Total nEvts in " + inputArguments.escapedInputFilePattern + ": " + str(nEvts))
outputFile.close()
