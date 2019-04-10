#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, ROOT, tmJDLInterface

# Make sure that at most one instance is running at a time
if (os.path.isfile("submitEventSelectionJobs.lock")):
    sys.exit("ERROR: only one instance of event selector can run at a time!")
else:
    os.system("touch submitEventSelectionJobs.lock")

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Submit jobs for final event selection.')
inputArgumentsParser.add_argument('--inputFilesList', required=True, help="Path to file containing list of input files.", type=str)
inputArgumentsParser.add_argument('--nEventsInInputFilesList', default=-1, help='Do not explicitly calculate the number of events in input files list; take it without verification to be this value.',type=int)
inputArgumentsParser.add_argument('--isMC', required=True, help="Takes values \"true\" or \"false\" indicating whether or not input file is a MC sample -- if so, disable HLT photon trigger and enable additional MC selection.", type=str)
inputArgumentsParser.add_argument('--photonSelectionType', required=True, help="Photon selection type: can be any one of: \"fake\", \"medium\", \"mediumfake\", or \"singlemedium\".", type=str)
inputArgumentsParser.add_argument('--year', required=True, help="Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger (unless sample is MC), and the photon ID cuts which are based on year-dependent recommendations.", type=str)
inputArgumentsParser.add_argument('--outputFilePrefix', required=True, help='Prefix to output file name.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', required=True, help='Output directory name.',type=str)
inputArgumentsParser.add_argument('--isDryRun', action='store_true', help="Do not submit the actual jobs: instead, only print the shell command that would have been called.")
inputArguments = inputArgumentsParser.parse_args()

nEvtsPerOutputFile = 0
if (inputArguments.isMC == "true"): nEvtsPerOutputFile = (10**6)
elif (inputArguments.isMC == "false"): nEvtsPerOutputFile = (10**7)
else:
    os.system("rm -f submitEventSelectionJobs.lock")
    sys.exit("argument \"isMC\" can be one of \"true\" or \"false\".")

print(" >> Submitting jobs for running event selection...")

nEvts = inputArguments.nEventsInInputFilesList

if (nEvts < 0):
    listOfInputFiles = []
    inputFileNamesFileObject = open(inputArguments.inputFilesList, 'r')
    for inputFileName in inputFileNamesFileObject:
        listOfInputFiles.append(inputFileName.strip())
    inputFileNamesFileObject.close()

    # Load input TTrees into TChain
    ggIn = ROOT.TChain("ggNtuplizer/EventTree")

    for inputFile in listOfInputFiles:
        # print("Adding: " + inputFile)
        ggIn.Add(inputFile)

    nEvts = ggIn.GetEntries()

print(" >> total nEvts:" + str(nEvts))

if not(nEvts > 0):
    os.system("rm -f submitEventSelectionJobs.lock")
    sys.exit("Found 0 events!")

currentWorkingDirectory = os.getcwd()
# Make sure the tarballs to transfer are up to date
updateCommand = "cd ~/private && ./update_tmUtilsTarball.sh && cd {cWD} && ./update_eventSelectionTarball.sh && cd {cWD}".format(cWD=currentWorkingDirectory)
os.system(updateCommand)
# Copy event selection helper script into the working directory
copyCommand = "cp -u eventSelectionHelper.sh condor_working_directory/."
os.system(copyCommand)

filesToTransfer = ["/uscms/home/tmudholk/private/tmUtils.tar.gz", "/uscms/home/tmudholk/private/extract_tmUtilsTarball.sh", "/uscms/home/tmudholk/private/stealth/STEALTH/eventSelection.tar.gz", "/uscms/home/tmudholk/private/stealth/STEALTH/extract_eventSelectionTarball.sh", "/uscms/home/tmudholk/private/stealth/STEALTH/STRegionBoundaries.dat", "/uscms/home/tmudholk/private/stealth/STEALTH/{iFL}".format(iFL=inputArguments.inputFilesList)]

startCounter = 0
endCounter = 0
outputIndex = 1
while endCounter < nEvts:
    endCounter = startCounter + nEvtsPerOutputFile - 1
    isLastIteration = (endCounter >= nEvts)
    if isLastIteration: endCounter = (nEvts - 1)
    processIdentifier = "selectionJob_{prefix}_begin_{sC}_end_{eC}".format(prefix=inputArguments.outputFilePrefix, sC=startCounter, eC=endCounter)
    outputFileName = "{prefix}_begin_{sC}_end_{eC}.root".format(prefix=inputArguments.outputFilePrefix, sC=startCounter, eC=endCounter)
    jdlInterface = tmJDLInterface.tmJDLInterface(processName=processIdentifier, scriptPath="eventSelectionHelper.sh", outputDirectoryRelativePath="condor_working_directory")
    jdlInterface.addFilesToTransferFromList(filesToTransfer)
    # Arguments for script:
    # Note: it seems simpler and certainly more readable to just include the "=" signs with the argument names, but I'm not sure whether that is allowed by JDL.
    jdlInterface.addScriptArgument("{iFL}".format(iFL=inputArguments.inputFilesList)) # Argument 1: inputFilesList
    jdlInterface.addScriptArgument("{oFN}".format(oFN=outputFileName)) # Argument 2: outputFilePath
    jdlInterface.addScriptArgument("{isMC}".format(isMC=inputArguments.isMC)) # Argument 3: isMC
    jdlInterface.addScriptArgument("{sC}".format(sC=startCounter)) # Argument 4: counterStartInclusive
    jdlInterface.addScriptArgument("{eC}".format(eC=endCounter)) # Argument 5: counterEndInclusive
    jdlInterface.addScriptArgument("{pST}".format(pST=inputArguments.photonSelectionType)) # Argument 6: photonSelectionType
    jdlInterface.addScriptArgument("{year}".format(year=inputArguments.year)) # Argument 7: year

    # Other arguments:
    jdlInterface.addScriptArgument("{oD}".format(oD=inputArguments.outputDirectory)) # Argument 8: output folder path

    # Write JDL
    jdlInterface.writeToFile()

    submissionCommand = "cd condor_working_directory && condor_submit {pI}.jdl && cd ..".format(pI=processIdentifier)
    print ("Generated command: {sC}".format(sC=submissionCommand))
    if (inputArguments.isDryRun):
        print("Not submitting due to dryRun flag.")
    else:
        os.system(submissionCommand)
        print ("Submitted.")
    if isLastIteration: break
    startCounter = 1+endCounter
    if (startCounter >= nEvts): break
    outputIndex += 1

os.system("rm -f submitEventSelectionJobs.lock")
