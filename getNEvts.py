#!/usr/bin/env python

from __future__ import print_function, division

import argparse, ROOT, sys

inputArgumentsParser = argparse.ArgumentParser(description='Print total number of events in a file pattern to a file.')
inputArgumentsParser.add_argument('--inputFilesList', required=True, help="Path to file containing list of input files.", type=str)
inputArguments = inputArgumentsParser.parse_args()

listOfInputFiles = []
inputFileNamesFileObject = open(inputArguments.inputFilesList, 'r')
for inputFileName in inputFileNamesFileObject:
    listOfInputFiles.append(inputFileName.strip())
inputFileNamesFileObject.close()

# Load input TTrees into TChain
ggIn = ROOT.TChain("ggNtuplizer/EventTree")

for inputFile in listOfInputFiles:
    # print("Adding: {inputfile}".format(inputfile=inputFile))
    readStatus = ggIn.Add(inputFile, 0)
    if not(readStatus == 1): sys.exit("File {iF} does not exist or does not contain the correct tree.".format(iF=inputFile))

nEvts = ggIn.GetEntries()
print("{n}".format(n=nEvts))
