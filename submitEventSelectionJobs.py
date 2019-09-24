#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, re, ROOT, tmJDLInterface

# Environment variables
stealthRoot = os.getenv("STEALTH_ROOT")
stealthEOSRoot = os.getenv("STEALTH_EOS_ROOT")
stealthArchives = os.getenv("STEALTH_ARCHIVES")
EOSPrefix = os.getenv("EOSPREFIX")
tmUtilsParent = os.getenv("TM_UTILS_PARENT")
hostname = os.getenv("HOSTNAME")
x509Proxy = os.getenv("X509_USER_PROXY")
habitat = ""
if ("lxplus" in hostname):
    habitat = "lxplus"
elif ("fnal" in hostname):
    habitat = "fnal"
else:
    sys.exit("ERROR: Unrecognized hostname: {h}, seems to be neither lxplus nor fnal.".format(h=hostname))

print("Environment variables:")
print("stealthRoot={sR}".format(sR=stealthRoot))
print("stealthEOSRoot={sER}".format(sER=stealthEOSRoot))
print("stealthArchives={sA}".format(sA=stealthArchives))
print("EOSPrefix={eP}".format(eP=EOSPrefix))
print("tmUtilsParent={tUP}".format(tUP=tmUtilsParent))
print("hostname={hN}".format(hN=hostname))
print("x509Proxy={xP}".format(xP=x509Proxy))

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
inputArgumentsParser.add_argument('--outputDirectory_selections', default="{sER}/selections/DoublePhoton".format(sER=stealthEOSRoot), help='Output directory name in which to store event selections.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_statistics', default="{sER}/statistics/DoublePhoton".format(sER=stealthEOSRoot), help='Output directory name in which to store statistics histograms.',type=str)
inputArgumentsParser.add_argument('--enable_cache', action='store_true', help="Read in number of input events as previously cached values.")
inputArgumentsParser.add_argument('--disableJetSelection', action='store_true', help="Disable jet selection.")
inputArgumentsParser.add_argument('--isProductionRun', action='store_true', help="By default, this script does not submit the actual jobs and instead only prints the shell command that would have been called. Passing this switch will execute the commands.")
inputArgumentsParser.add_argument('--preserveLogs', action='store_true', help="By default, this script moves all event selection logs to the archives. This switch will keep the logs where they are (but they may be overwritten).")
inputArgumentsParser.add_argument('--preserveInputFileLists', action='store_true', help="By default, this script regenerates the input file lists to be fed to the statistics and selection merging scripts. This switch will preserve the input file lists.")
inputArguments = inputArgumentsParser.parse_args()

def execute_in_env(commandToRun, printDebug=False):
    env_setup_command = "bash -c \"cd {sR} && source setupEnv.sh".format(sR=stealthRoot)
    runInEnv = "{e_s_c} && set -x && {c} && set +x\"".format(e_s_c=env_setup_command, c=commandToRun)
    if (printDebug):
        print("About to execute command:")
        print("{c}".format(c=runInEnv))
    os.system(runInEnv)

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
    os.system("set -x && mkdir -p {sA}/logs && rsync --quiet --progress -a condor_working_directory/ {sA}/logs/ && rm -rf condor_working_directory/* && set +x".format(sA=stealthArchives))

selectionTypesToRun = []
if (inputArguments.runOnlyDataOrMC == "data"):
    selectionTypesToRun.append("data")
elif (inputArguments.runOnlyDataOrMC == "MC"):
    # selectionTypesToRun.append("MC_stealth_t6")
    selectionTypesToRun.append("MC_stealth_t5")
    # selectionTypesToRun.append("MC_hgg")
elif (inputArguments.runOnlyDataOrMC == "all"):
    selectionTypesToRun.append("data")
    # selectionTypesToRun.append("MC_stealth_t6")
    selectionTypesToRun.append("MC_stealth_t5")
    # selectionTypesToRun.append("MC_hgg")
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
        2017: "fileLists/inputFileList_data_DoubleEG_2017_producedAug2019.txt"
    },
    "MC_stealth_t6": {
        2017: "fileLists/inputFileList_MC_Fall17_stealth_t6Wg.txt"
    },
    "MC_stealth_t5": {
        2017: "fileLists/inputFileList_MC_Fall17_stealth_t5Wg.txt"
    },
    "MC_hgg": {
        2017: "fileLists/inputFileList_MC_Fall17_hgg.txt"
    }
}

cached_nEvents_lists = {
    "data": {
        2016: "cached/cached_nEvents_inputFileList_data_DoubleEG_2016_ntuplized80X.txt",
        2017: "cached/cached_nEvents_inputFileList_data_DoubleEG_2017_jul2019.txt"
    },
    "MC_stealth_t6": {
        2017: "cached/cached_nEvents_inputFileList_MC_Fall17_stealth_t6Wg.txt"
    },
    "MC_stealth_t5": {
        2017: "cached/cached_nEvents_inputFileList_MC_Fall17_stealth_t5Wg.txt"
    },
    "MC_hgg": {
        2017: "cached/cached_nEvents_inputFileList_MC_Fall17_hgg.txt"
    }
}

execute_in_env("eos {eP} mkdir -p {oD}{oI}".format(eP=EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier), printDebug=True)
execute_in_env("eos {eP} mkdir -p {oD}{oI}".format(eP=EOSPrefix, oD=inputArguments.outputDirectory_statistics, oI=optional_identifier), printDebug=True)

# Make sure the tarballs to transfer are up to date
updateCommand = "cd {tUP} && ./update_tmUtilsTarball.sh && cd {sR} && ./update_eventSelectionTarball.sh && cd {sR}".format(tUP=tmUtilsParent, sR=stealthRoot)
os.system(updateCommand)
# Copy event selection helper script into the working directory
copyCommand = "cp -u eventSelectionHelper.sh condor_working_directory/."
os.system(copyCommand)

if not(inputArguments.preserveInputFileLists):
    os.system("rm fileLists/inputFileList_selections_*.txt && rm fileLists/inputFileList_statistics_*.txt")

for selectionType in selectionTypesToRun:
    isMC="none"
    if (re.match(r"MC_", selectionType)):
        isMC = "true" # the string "true", not the boolean -- because this is the format expected by the event selection script
    elif (selectionType == "data"):
        isMC = "false"
    disableJetSelectionString = "none"
    if (inputArguments.disableJetSelection): disableJetSelectionString = "true"
    else: disableJetSelectionString = "false"
    for year in yearsToRun:
        inputFilesList = fileLists[selectionType][year]
        cached_nEvents_list = cached_nEvents_lists[selectionType][year]
        print("Submitting jobs for year={y}, selection type={t}".format(y=year, t=selectionType))
        nEvtsPerOutputFile = 0
        if (re.match(r"MC_", selectionType)):
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

        filesToTransfer = ["{xP}".format(xP=x509Proxy), "{tUP}/tmUtils.tar.gz".format(tUP=tmUtilsParent), "{tUP}/extract_tmUtilsTarball.sh".format(tUP=tmUtilsParent), "{sR}/eventSelection.tar.gz".format(sR=stealthRoot), "{sR}/extract_eventSelectionTarball.sh".format(sR=stealthRoot), "{sR}/{iFL}".format(sR=stealthRoot, iFL=inputFilesList), "{sR}/STRegionBoundaries.dat".format(sR=stealthRoot)]
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
            jdlInterface.addScriptArgument("{dJS}".format(dJS=disableJetSelectionString)) # Argument 3: disableJetSelection
            jdlInterface.addScriptArgument("{sC}".format(sC=startCounter)) # Argument 4: counterStartInclusive
            jdlInterface.addScriptArgument("{eC}".format(eC=endCounter)) # Argument 5: counterEndInclusive
            jdlInterface.addScriptArgument("{y}".format(y=year)) # Argument 6: year

            # Other arguments:
            jdlInterface.addScriptArgument("{oD}{oI}".format(oD=inputArguments.outputDirectory_selections, oI=optional_identifier)) # Argument 7: selections output folder path
            jdlInterface.addScriptArgument("{oD}{oI}".format(oD=inputArguments.outputDirectory_statistics, oI=optional_identifier)) # Argument 8: statistics output folder path
            jdlInterface.addScriptArgument("{sT}".format(sT=selectionType)) # Argument 9: selection type

            if (habitat == "lxplus"):
                jdlInterface.setFlavor("tomorrow")
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
                os.system("echo \"{eP}/{oD}{oI}/selection_{t}_{y}_signal_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}_{y}{oI}_signal.txt".format(eP=EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, y=year, b=startCounter, e=endCounter))
                os.system("echo \"{eP}/{oD}{oI}/selection_{t}_{y}_control_fakefake_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}_{y}{oI}_control_fakefake.txt".format(eP=EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, y=year, b=startCounter, e=endCounter))
                os.system("echo \"{eP}/{oD}{oI}/selection_{t}_{y}_control_mediumfake_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}_{y}{oI}_control_mediumfake.txt".format(eP=EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, y=year, b=startCounter, e=endCounter))
                os.system("echo \"{eP}/{oD}{oI}/statistics_{t}_{y}_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_statistics_{t}_{y}{oI}.txt".format(eP=EOSPrefix, oD=inputArguments.outputDirectory_statistics, oI=optional_identifier, t=selectionType, y=year, b=startCounter, e=endCounter))
            if isLastIteration: break
            startCounter = 1+endCounter
            if (startCounter >= nEvts): break

os.system("rm -f submitEventSelectionJobs.lock")
