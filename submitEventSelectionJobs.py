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
inputArgumentsParser.add_argument('--runOnlyDataOrMC', default="all", help="Takes values \"data\" or \"MC\" if selection is to be run only on data or only on MC samples. For MC selections, disable HLT photon trigger and enable additional MC selection. Default is \"all\", which means both selectionsare run.", type=str)
inputArgumentsParser.add_argument('--year', default="all", help="Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger (unless sample is MC), and the photon ID cuts which are based on year-dependent recommendations.", type=str)
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_selections', default="/store/user/lpcsusystealth/selections/DoublePhoton", help='Output directory name in which to store event selections.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_statistics', default="/store/user/lpcsusystealth/statistics/DoublePhoton", help='Output directory name in which to store statistics histograms.',type=str)
inputArgumentsParser.add_argument('--enable_cache', action='store_true', help="Read in number of input events as previously cached values.")
inputArgumentsParser.add_argument('--isProductionRun', action='store_true', help="By default, this script does not submit the actual jobs and instead only prints the shell command that would have been called. Passing this switch will execute the commands.")
inputArgumentsParser.add_argument('--preserveLogs', action='store_true', help="By default, this script moves all event selection logs to the archives. This switch will keep the logs where they are (but they may be overwritten).")
inputArgumentsParser.add_argument('--preserveInputFileLists', action='store_true', help="By default, this script regenerates the input file lists to be fed to the statistics and selection merging scripts. This switch will preserve the input file lists.")
inputArguments = inputArgumentsParser.parse_args()

def execute_in_env(command):
    env_setup_command = "cd /uscms/home/tmudholk/private/stealth/STEALTH && source setupEnv.sh"
    os.system("{e_s_c} && set -x && {c} && set +x".format(e_s_c = env_setup_command, c=command))

def get_nEvts_from_file(fileName):
    nEvts_array = []
    nEventsFileObject = open(cached_nEvents_list, 'r')
    for line in nEventsFileObject:
        nEvts_array.append(line.strip())
    nEventsFileObject.close()
    if (not(len(nEvts_array) == 1)):
        print("cache file contains no lines or more than one line. Array of values in file: {s}".format(s = str(nEvts_array)))
        return (-1)
    elif (not((nEvts_array[0]).isdigit())):
        print("cache file contains nEvents in an invalid format. Found string: {s}".format(s = str(nEvts_array[0])))
        return (-1)
    else:
        return int(nEvts_array[0])

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)

if not(inputArguments.preserveLogs):
    os.system("set -x && mkdir -p ~/nobackup/archives/logs && mv condor_working_directory/* ~/nobackup/archives/logs/. && set +x")

selectionTypesToRun = []
if (inputArguments.runOnlyDataOrMC == "data"):
    selectionTypesToRun.append("data")
elif (inputArguments.runOnlyDataOrMC == "MC"):
    selectionTypesToRun.append("MC")
elif (inputArguments.runOnlyDataOrMC == "all"):
    selectionTypesToRun.append("data")
    selectionTypesToRun.append("MC")
else:
    sys.exit("ERROR: invalid value for argument \"runOnlyDataOrMC\": {v}".format(v=inputArguments.runOnlyDataOrMC))

yearsToRun = []
if (inputArguments.year == "2016"):
    yearsToRun.append(2016)
elif (inputArguments.year == "2017"):
    yearsToRun.append(2017)
elif (inputArguments.year == "all"):
    yearsToRun.append(2016)
    yearsToRun.append(2017)
else:
    sys.exit("ERROR: invalid value for argument \"year\": {v}".format(v=inputArguments.year))

fileLists = {
    "data": {
        2016: "fileLists/inputFileList_data_DoubleEG_2016_ntuplized80X.txt",
        2017: "fileLists/inputFileList_data_DoubleEG_2017_ntuplized949cand2.txt"
    },
    "MC": {
        2016: "fileLists/inputFileList_MC_2018Production_ntuplized80X.txt",
        2017: "fileLists/inputFileList_MC_2018Production_ntuplized949cand2.txt"
    }
}

cached_nEvents_lists = {
    "data": {
        2016: "cached/cached_nEvents_inputFileList_data_DoubleEG_2016_ntuplized80X.txt",
        2017: "cached/cached_nEvents_inputFileList_data_DoubleEG_2017_ntuplized949cand2.txt"
    },
    "MC": {
        2016: "cached/cached_nEvents_inputFileList_MC_2018Production_ntuplized80X.txt",
        2017: "cached/cached_nEvents_inputFileList_MC_2018Production_ntuplized949cand2.txt"
    }
}

execute_in_env("eos root://cmseos.fnal.gov mkdir -p {oD}{oI}".format(oD=inputArguments.outputDirectory_selections, oI=optional_identifier))
execute_in_env("eos root://cmseos.fnal.gov mkdir -p {oD}{oI}".format(oD=inputArguments.outputDirectory_statistics, oI=optional_identifier))

currentWorkingDirectory = os.getcwd()
# Make sure the tarballs to transfer are up to date
updateCommand = "cd ~/private && ./update_tmUtilsTarball.sh && cd {cWD} && ./update_eventSelectionTarball.sh && cd {cWD}".format(cWD=currentWorkingDirectory)
os.system(updateCommand)
# Copy event selection helper script into the working directory
copyCommand = "cp -u eventSelectionHelper.sh condor_working_directory/."
os.system(copyCommand)

if not(inputArguments.preserveInputFileLists):
    os.system("rm fileLists/inputFileList_selections_*.txt && rm fileLists/inputFileList_statistics_*.txt")

for selectionType in selectionTypesToRun:
    isMC="none"
    if (selectionType == "MC"):
        isMC = "true" # the string "true", not the boolean -- because this is the format expected by the event selection script
    elif (selectionType == "data"):
        isMC = "false"
    for year in yearsToRun:
        inputFilesList = fileLists[selectionType][year]
        cached_nEvents_list = cached_nEvents_lists[selectionType][year]
        print("Submitting jobs for year={y}, selection type={t}".format(y=year, t=selectionType))
        nEvtsPerOutputFile = 0
        if (selectionType == "MC"):
            nEvtsPerOutputFile = (5*(10**5))
        elif (selectionType == "data"):
            nEvtsPerOutputFile = (5*(10**6))
        nEvts = 0
        cache_needs_regeneration = True
        if (inputArguments.enable_cache):
            if not(os.path.isfile(cached_nEvents_list)):
                print("Unable to find file containing cached number of events: tried searching for \"{f}\"".format(f=cached_nEvents_list))
            else:
                nEvts = get_nEvts_from_file(cached_nEvents_list)
                if (nEvts < 0):
                    cache_needs_regeneration = True
                else:
                    cache_needs_regeneration = False
        if (cache_needs_regeneration):
            print("Caching number of events for selectionType = {sT}, year = {y}".format(sT=selectionType, y=year))
            execute_in_env("./getNEvts.py --inputFilesList {iFL} | tee {o}".format(iFL=inputFilesList, o=cached_nEvents_list))
            nEvts = get_nEvts_from_file("{o}".format(o=cached_nEvents_list))

        print("Total available nEvts:" + str(nEvts))

        if not(nEvts > 0):
            os.system("rm -f submitEventSelectionJobs.lock")
            sys.exit("Found 0 events!")

        filesToTransfer = ["/uscms/home/tmudholk/private/tmUtils.tar.gz", "/uscms/home/tmudholk/private/extract_tmUtilsTarball.sh", "/uscms/home/tmudholk/private/stealth/STEALTH/eventSelection.tar.gz", "/uscms/home/tmudholk/private/stealth/STEALTH/extract_eventSelectionTarball.sh", "/uscms/home/tmudholk/private/stealth/STEALTH/{iFL}".format(iFL=inputFilesList)]
        formatted_iFL = (inputFilesList.split("/"))[-1]

        startCounter = 0
        endCounter = 0
        while endCounter < nEvts:
            endCounter = startCounter + nEvtsPerOutputFile - 1
            isLastIteration = (endCounter >= nEvts)
            if isLastIteration: endCounter = (nEvts - 1)
            processIdentifier = "selectionJob_{sT}_{y}_begin_{sC}_end_{eC}".format(sT=selectionType, y=year, sC=startCounter, eC=endCounter)
            jdlInterface = tmJDLInterface.tmJDLInterface(processName=processIdentifier, scriptPath="eventSelectionHelper.sh", outputDirectoryRelativePath="condor_working_directory")
            jdlInterface.addFilesToTransferFromList(filesToTransfer)
            # Arguments for script:
            # Note: it seems simpler and certainly more readable to just include the "=" signs with the argument names, but I'm not sure whether that is allowed by JDL.
            jdlInterface.addScriptArgument("{iFL}".format(iFL=formatted_iFL)) # Argument 1: inputFilesList
            jdlInterface.addScriptArgument("{iMC}".format(iMC=isMC)) # Argument 2: isMC
            jdlInterface.addScriptArgument("{sC}".format(sC=startCounter)) # Argument 3: counterStartInclusive
            jdlInterface.addScriptArgument("{eC}".format(eC=endCounter)) # Argument 4: counterEndInclusive
            jdlInterface.addScriptArgument("{y}".format(y=year)) # Argument 5: year

            # Other arguments:
            jdlInterface.addScriptArgument("{oD}{oI}".format(oD=inputArguments.outputDirectory_selections, oI=optional_identifier)) # Argument 6: selections output folder path
            jdlInterface.addScriptArgument("{oD}{oI}".format(oD=inputArguments.outputDirectory_statistics, oI=optional_identifier)) # Argument 7: statistics output folder path

            # Write JDL
            jdlInterface.writeToFile()

            submissionCommand = "cd condor_working_directory && condor_submit {pI}.jdl && cd ..".format(pI=processIdentifier)
            print ("Generated command: {sC}".format(sC=submissionCommand))
            if (inputArguments.isProductionRun):
                os.system(submissionCommand)
                print ("Submitted.")
            else:
                print("Not submitting because productionRun flag was not set.")
            if not(inputArguments.preserveInputFileLists):
                os.system("echo \"root://cmseos.fnal.gov/{oD}{oI}/selection_{t}_{y}_signal_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}_{y}{oI}_signal.txt".format(oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, y=year, b=startCounter, e=endCounter))
                os.system("echo \"root://cmseos.fnal.gov/{oD}{oI}/selection_{t}_{y}_control_fakefake_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}_{y}{oI}_control_fakefake.txt".format(oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, y=year, b=startCounter, e=endCounter))
                os.system("echo \"root://cmseos.fnal.gov/{oD}{oI}/selection_{t}_{y}_control_mediumfake_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}_{y}{oI}_control_mediumfake.txt".format(oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, y=year, b=startCounter, e=endCounter))
                os.system("echo \"root://cmseos.fnal.gov/{oD}{oI}/statistics_{t}_{y}_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_statistics_{t}_{y}{oI}.txt".format(oD=inputArguments.outputDirectory_statistics, oI=optional_identifier, t=selectionType, y=year, b=startCounter, e=endCounter))
            if isLastIteration: break
            startCounter = 1+endCounter
            if (startCounter >= nEvts): break

os.system("rm -f submitEventSelectionJobs.lock")
