#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, ROOT

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Submit jobs for final event selection.')
inputArgumentsParser.add_argument('--eventSelectHelperScriptName', default="submitJobs_selectEvents_Helper_Condor.sh", help='Path to helper script for event selection.', type=str)
inputArgumentsParser.add_argument('--inputFromFile', action='store_true', help="Interpret inputFilePath as text file that has a list of input of files.")
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--workingDirectory', default='/uscms/home/tmudholk/private/stealth/STEALTH/condor_working_directory', help='Path to working directory.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', default='finalSelection', help='Output directory name.',type=str)
inputArgumentsParser.add_argument('--outputFilePrefix', default='DoubleEG_FebReminiAOD_finalSelection', help='Prefix to output file name.',type=str)
inputArgumentsParser.add_argument('--nEvtsPerOutputFile', default=(10**6), help="Number of events per output file.", type=int)
inputArgumentsParser.add_argument('--photonSelectionType', default="fake", help='Takes value fake for fake photon selection and medium for selection based on medium ID.',type=str)
inputArgumentsParser.add_argument('--year', default=-1, help='Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger, and the photon ID cuts which are based on year-dependent recommendations. Default year: -1, which is used for MC and means the trigger is disabled and the 2017 photon ID recommendations are implemented.', type=int) # Bit 14 for 2016 data: HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90
inputArgumentsParser.add_argument('--JECUncertainty', default=0, help='Apply a uniform upward or downward jet energy uncertainty correction to jet pt. Default: 0, i.e. do not apply any other correction. +/-1 are allowed as well, shifting all jet pt up or down respectively by the relevant jet energy correction.', type=int)
inputArgumentsParser.add_argument('--isDryRun', action='store_true', help="Do not submit the actual jobs: instead, only print the shell command that would have been called.")
inputArguments = inputArgumentsParser.parse_args()

print(" >> Submitting jobs for running event selection...")

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

currentWorkingDirectory = os.getcwd()
copyCommand = "cp -u {eventSelectHelperScriptName} {workingDirectory}/. && cd {currentWorkingDirectory}".format(eventSelectHelperScriptName=inputArguments.eventSelectHelperScriptName, workingDirectory=inputArguments.workingDirectory, currentWorkingDirectory=currentWorkingDirectory)
os.system(copyCommand)

startCounter = 0
endCounter = 0
outputIndex = 1
while endCounter < nEvts:
    endCounter = startCounter + inputArguments.nEvtsPerOutputFile - 1
    isLastIteration = (endCounter >= nEvts)
    if isLastIteration: endCounter = (nEvts - 1)
    outputFileName = inputArguments.outputFilePrefix + "_begin_{startCounter}_end_{endCounter}".format(startCounter=startCounter, endCounter=endCounter) + ".root"
    jdlPrefix = "selectEvents_outputFile_{outputFilePrefix}_begin_{startCounter}_end_{endCounter}".format(outputFilePrefix=inputArguments.outputFilePrefix, startCounter=startCounter, endCounter=endCounter)
    jdlFileName = jdlPrefix + ".jdl"
    jdlLogPrefix = jdlPrefix
    jdlCreationCommand = "./createJDL_eventSelection.sh {workingDirectory}/{jdlFileName} {inputFilePath} {outputFileName} {startCounter} {endCounter} {selectionType} {year} {JECUncertainty} {outputDirectory} {jdlLogPrefix}".format(workingDirectory=inputArguments.workingDirectory, jdlFileName=jdlFileName, inputFilePath=inputArguments.inputFilePath, outputFileName=outputFileName, startCounter=startCounter, endCounter=endCounter, selectionType=inputArguments.photonSelectionType, year=inputArguments.year, JECUncertainty=inputArguments.JECUncertainty, outputDirectory=inputArguments.outputDirectory, jdlLogPrefix=jdlLogPrefix)
    if (inputArguments.inputFromFile): jdlCreationCommand += " inputFromFile"
    else: jdlCreationCommand += " DisableInputFromFile"
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
