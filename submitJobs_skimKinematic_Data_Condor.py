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
# inputArgumentsParser.add_argument('--inputFromFile', action='store_true', help="Interpret inputFilePath as text file that has a list of input of files.")
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--workingDirectory', required=True, help='Path to working directory.',type=str)
inputArgumentsParser.add_argument('--outputFilePrefix', required=True, help='Prefix to output file name.',type=str)
inputArgumentsParser.add_argument('--nEvtsPerOutputFile', default=(10**6), help="Number of events per output file.", type=int)
inputArgumentsParser.add_argument('--isDryRun', action='store_true', help="Do not submit the actual jobs; instead, only print the shell command that would have been called.")
inputArguments = inputArgumentsParser.parse_args()

print(" >> Submitting jobs for running skim for two kinematic photons...")

# listOfInputFiles = []
# if (inputArguments.inputFromFile):
#     inputFileNamesFileObject = open(inputArguments.inputFilePath, 'r')
#     for inputFileName in inputFileNamesFileObject:
#         listOfInputFiles.append(inputFileName.strip())
#     inputFileNamesFileObject.close()
# else:
#     listOfInputFiles.append(inputArguments.inputFilePath)
    
# Load input TTrees into TChain
ggIn = ROOT.TChain("ggNtuplizer/EventTree")
# for inputFile in listOfInputFiles:
#     print("Adding: " + inputFile)
#     ggIn.Add(inputFile)

print("Adding: " + inputArguments.inputFilePath)
ggIn.Add(inputArguments.inputFilePath)

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
    jdlFileName = "skimKinematic_begin_{startCounter}_end_{endCounter}.jdl".format(startCounter=startCounter, endCounter=endCounter)
    outputFileName = inputArguments.outputFilePrefix + "_begin_{startCounter}_end_{endCounter}".format(startCounter=startCounter, endCounter=endCounter) + ".root"
    jdlCreationCommand = "./createJDL.sh {workingDirectory}/{jdlFileName} {inputFilePath} {outputFileName} {startCounter} {endCounter}".format(workingDirectory=inputArguments.workingDirectory, jdlFileName=jdlFileName, inputFilePath=inputArguments.inputFilePath, outputFileName=outputFileName, startCounter=startCounter, endCounter=endCounter)
    commandToCall = "cd {workingDirectory} && condor_submit {jdlFileName} && cd -".format(workingDirectory=inputArguments.workingDirectory, jdlFileName=jdlFileName)
    print ("Creating JDL: " + jdlCreationCommand)
    os.system(jdlCreationCommand)
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
