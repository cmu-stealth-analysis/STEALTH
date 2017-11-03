#!/usr/bin/env python

from __future__ import print_function, division

import os
import sys
# import numpy as np
import argparse
import ROOT

from tmProgressBar import tmProgressBar

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Submit jobs for photon skim based on kinematic observables only.')
# inputArgumentsParser.add_argument('-e','--era', required=True, help='Run Era',type=str)
inputArgumentsParser.add_argument('--inputFromFile', action='store_true', help="Interpret inputFilePath as text file that has a list of input of files.")
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--outputFolder', required=True, help='Prefix to output file name.',type=str)
inputArgumentsParser.add_argument('--outputFilePrefix', required=True, help='Prefix to output file name.',type=str)
inputArgumentsParser.add_argument('--logDirectory', required=True, help='Prefix to output file name.',type=str)
inputArgumentsParser.add_argument('--nEvtsPerOutputFile', default=(10**6), help="number of events per output file", type=int)
inputArgumentsParser.add_argument('--queue', default="1nh", help="queue to submit to for each set of events", type=str)
inputArgumentsParser.add_argument('--isDryRun', action='store_true', help="Do not submit the actual jobs; instead print what would have been submitted.")
inputArguments = inputArgumentsParser.parse_args()

print(" >> Submitting jobs for running skim for two kinematic photons...")

listOfInputFiles = []
if (inputArguments.inputFromFile):
    inputFileNamesFileObject = open(inputArguments.inputFilePath, 'r')
    for inputFileName in inputFileNamesFileObject:
        listOfInputFiles.append(inputFileName.strip())
    inputFileNamesFileObject.close()
else:
    listOfInputFiles.append(inputArguments.inputFilePath)
    
# Load input TTrees into TChain
ggIn = ROOT.TChain("ggNtuplizer/EventTree")
for inputFile in listOfInputFiles:
    print("Adding: " + inputFile)
    ggIn.Add(inputFile)

nEvts = ggIn.GetEntries()
print(" >> total nEvts:" + str(nEvts))

if not(nEvts > 0): sys.exit("Found 0 events!")

startCounter = 0
endCounter = 0
outputIndex = 1
while endCounter < nEvts:
    endCounter = startCounter + inputArguments.nEvtsPerOutputFile - 1
    isLastIteration = (endCounter >= nEvts)
    if isLastIteration: endCounter = (nEvts - 1)
    outputFilePath = inputArguments.outputFolder + "/" + inputArguments.outputFilePrefix + "_" + str(outputIndex) + ".root"
    logFilePath = inputArguments.logDirectory + "/log_" + inputArguments.outputFilePrefix + "_" + str(outputIndex) + ".log"
    commandToCall = "bsub -q {queue} submitJobs_skimKinematic_Data_Helper.sh {inputFilePath} {outputFilePath} {startCounter} {endCounter} {logFilePath}".format(queue=inputArguments.queue, inputFilePath=inputArguments.inputFilePath, outputFilePath=outputFilePath, startCounter=startCounter, endCounter=endCounter, logFilePath=logFilePath)
    if inputArguments.inputFromFile: commandToCall += " inputFromFile"
    print ("Generated command: " + commandToCall)
    if (inputArguments.isDryRun):
        print("Not submitting due to dryRun flag.")
    else:
        os.system(commandToCall)
        print ("Submitted.")
    if isLastIteration: break
    startCounter = 1+endCounter
    if (startCounter >= nEvts): break
    outputIndex += 1
