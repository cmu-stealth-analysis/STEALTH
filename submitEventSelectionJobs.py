#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, re, ROOT, tmJDLInterface, stealthEnv, commonFunctions

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Submit jobs for final event selection.')
inputArgumentsParser.add_argument('--runOnlyDataOrMC', default="all", help="Takes values \"data\" or \"MC\" if selection is to be run only on data or only on MC samples. For MC selections, disable HLT photon trigger and enable additional MC selection. Default is \"all\", which means both selections are run.", type=str)
inputArgumentsParser.add_argument('--year', default="all", help="Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger (unless sample is MC), and the photon ID cuts which are based on year-dependent recommendations.", type=str)
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_selections', default="{sER}/selections/DoublePhoton".format(sER=stealthEnv.stealthEOSRoot), help='Output directory name in which to store event selections.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_statistics', default="{sER}/statistics/DoublePhoton".format(sER=stealthEnv.stealthEOSRoot), help='Output directory name in which to store statistics histograms.',type=str)
# inputArgumentsParser.add_argument('--enable_cache', action='store_true', help="Read in number of input events as previously cached values.")
inputArgumentsParser.add_argument('--disableJetSelection', action='store_true', help="Disable jet selection.")
inputArgumentsParser.add_argument('--isProductionRun', action='store_true', help="By default, this script does not submit the actual jobs and instead only prints the shell command that would have been called. Passing this switch will execute the commands.")
inputArgumentsParser.add_argument('--preserveLogs', action='store_true', help="By default, this script moves all event selection logs to the archives. This switch will keep the logs where they are (but they may be overwritten).")
inputArgumentsParser.add_argument('--preserveInputFileLists', action='store_true', help="By default, this script regenerates the input file lists to be fed to the statistics and selection merging scripts. This switch will preserve the input file lists.")
inputArguments = inputArgumentsParser.parse_args()

def execute_in_env(commandToRun, printDebug=False):
    env_setup_command = "bash -c \"cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot)
    runInEnv = "{e_s_c} && set -x && {c} && set +x\"".format(e_s_c=env_setup_command, c=commandToRun)
    if (printDebug):
        print("About to execute command:")
        print("{c}".format(c=runInEnv))
    os.system(runInEnv)

# def read_nEvts_from_file(fileName):
#     nEvts_array = []
#     nEventsFileObject = open(fileName, 'r')
#     for line in nEventsFileObject:
#         nEvts_array.append(line.strip())
#     nEventsFileObject.close()
#     if (not(len(nEvts_array) == 1)):
#         print("cache file contains no lines or more than one line. Array of values in file: {s}".format(s = str(nEvts_array)))
#         return (-1)
#     elif (not((nEvts_array[0]).isdigit())):
#         print("cache file contains nEvents in an invalid format. Found string: {s}".format(s = str(nEvts_array[0])))
#         return (-1)
#     else:
#         return int(nEvts_array[0])

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)

if not(inputArguments.preserveLogs):
    os.system("set -x && mkdir -p {sA}/logs && rsync --quiet --progress -a {cWAR}/selection{oI}/ {sA}/logs/ && rm -rf {cWAR}/selection{oI}/* && set +x".format(cWAR=stealthEnv.condorWorkAreaRoot, sA=stealthEnv.stealthArchives, oI=optional_identifier))

os.system("mkdir -p {cWAR}/selection{oI}".format(cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier))

selectionTypesToRun = []
if (inputArguments.runOnlyDataOrMC == "data"):
    selectionTypesToRun.append("data")
elif (inputArguments.runOnlyDataOrMC == "MC"):
    selectionTypesToRun.append("MC_stealth_t5")
elif (inputArguments.runOnlyDataOrMC == "all"):
    selectionTypesToRun.append("data")
    selectionTypesToRun.append("MC_stealth_t5")
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
    # "MC_hgg": {
    #     2017: "fileLists/inputFileList_MC_Fall17_hgg.txt"
    # },
    # "MC_stealth_t6": {
    #     2017: "fileLists/inputFileList_MC_Fall17_stealth_t6Wg.txt"
    # },
    "MC_stealth_t5": {
        2016: "fileLists/inputFileList_MC_Fall17_stealth_t5Wg.txt",
        2017: "fileLists/inputFileList_MC_Fall17_stealth_t5Wg.txt"
    },
    "data": {
        2016: "fileLists/inputFileList_data_DoubleEG_2016_ntuplizedOct2019.txt",
        2017: "fileLists/inputFileList_data_DoubleEG_2017_ntuplizedOct2019.txt"
    }
}

target_nFilesPerJob = {
    "MC_stealth_t5": {
        2016: 25,
        2017: 25,
        # 2018: 25
    },
    "data": {
        2016: 150,
        2017: 150,
        # 2018: 200
    }
}

# cached_nEvents_lists = {
#     "data": {
#         2016: "cached/cached_nEvents_inputFileList_data_DoubleEG_2016_ntuplizedOct2019.txt",
#         2017: "cached/cached_nEvents_inputFileList_data_DoubleEG_2017_ntuplizedOct2019.txt"
#     },
#     # "MC_stealth_t6": {
#     #     2017: "cached/cached_nEvents_inputFileList_MC_Fall17_stealth_t6Wg.txt"
#     # },
#     "MC_stealth_t5": {
#         2016: "cached/cached_nEvents_inputFileList_MC_Fall17_stealth_t5Wg.txt",
#         2017: "cached/cached_nEvents_inputFileList_MC_Fall17_stealth_t5Wg.txt"
#     }# ,
#     # "MC_hgg": {
#     #     2017: "cached/cached_nEvents_inputFileList_MC_Fall17_hgg.txt"
#     # }
# }

execute_in_env("eos {eP} mkdir -p {oD}{oI}".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier), printDebug=True)
execute_in_env("eos {eP} mkdir -p {oD}{oI}".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_statistics, oI=optional_identifier), printDebug=True)

# Make sure the tarballs to transfer are up to date
updateCommand = "cd {tUP} && ./update_tmUtilsTarball.sh && cd {sR} && ./update_eventSelectionTarball.sh && cd {sR}".format(tUP=stealthEnv.tmUtilsParent, sR=stealthEnv.stealthRoot)
os.system(updateCommand)
# Copy event selection helper script into the working directory
copyCommand = "cd {sR} && cp -u eventSelectionHelper.sh {cWAR}/selection{oI}/.".format(sR=stealthEnv.stealthRoot, cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier)
os.system(copyCommand)

disableJetSelectionString = "none"
if (inputArguments.disableJetSelection): disableJetSelectionString = "true"
else: disableJetSelectionString = "false"

for selectionType in selectionTypesToRun:
    isMC="none"
    if (re.match(r"MC_", selectionType)):
        isMC = "true" # the string "true", not the boolean -- because this is the format expected by the event selection script
    elif (selectionType == "data"):
        isMC = "false"
    for year in yearsToRun:
        if not(inputArguments.preserveInputFileLists):
            os.system("cd {sR} && rm fileLists/inputFileList_selections_{t}_{y}{oI}_*.txt && rm fileLists/inputFileList_statistics_{t}_{y}{oI}.txt".format(oI=optional_identifier, t=selectionType, y=year, sR=stealthEnv.stealthRoot))
        inputPathsFile = fileLists[selectionType][year]
        # cached_nEvents_list = cached_nEvents_lists[selectionType][year]
        # nEvtsPerOutputFile = 0
        # if (re.match(r"MC_", selectionType)):
        #     nEvtsPerOutputFile = (5*(10**5))
        # elif (selectionType == "data"):
        #     nEvtsPerOutputFile = (5*(10**6))
        nFilesPerJob = target_nFilesPerJob[selectionType][year]
        print("Submitting jobs for year={y}, selection type={t}".format(y=year, t=selectionType))
        # nEvts = 0
        # cache_needs_regeneration = True
        # if (inputArguments.enable_cache):
        #     if (os.path.isfile(cached_nEvents_list)):
        #         nEvts = read_nEvts_from_file(fileName=cached_nEvents_list)
        #         if (nEvts < 0):
        #             print("Unable to fetch number of events from cache. Regenerating...")
        #             cache_needs_regeneration = True
        #         else:
        #             cache_needs_regeneration = False
        #     else:
        #         print("Unable to find file containing cached number of events: tried searching for \"{f}\"".format(f=cached_nEvents_list))
        # if (cache_needs_regeneration):
        #     print("Caching number of events for selectionType = {sT}, year = {y}".format(sT=selectionType, y=year))
        #     nEvts = commonFunctions.get_nEvts_from_fileList(inputPathsFile=inputPathsFile, printDebug=True)
        #     outputFile = open("{f}".format(f=cached_nEvents_list), 'w')
        #     outputFile.write("{n}\n".format(n=nEvts))
        #     outputFile.close()

        total_nLines = commonFunctions.get_number_of_lines_in_file(inputFilePath=inputPathsFile)
        print("Total available nLines: {n}".format(n=total_nLines))

        if not(total_nLines > 0):
            os.system("rm -f submitEventSelectionJobs.lock")
            sys.exit("ERROR: Found 0 lines in input path {p}.".format(p=inputPathsFile))

        filesToTransfer = ["{xP}".format(xP=stealthEnv.x509Proxy), "{tUP}/tmUtils.tar.gz".format(tUP=stealthEnv.tmUtilsParent), "{tUP}/extract_tmUtilsTarball.sh".format(tUP=stealthEnv.tmUtilsParent), "{sR}/eventSelection.tar.gz".format(sR=stealthEnv.stealthRoot), "{sR}/extract_eventSelectionTarball.sh".format(sR=stealthEnv.stealthRoot), "{sR}/{iPF}".format(sR=stealthEnv.stealthRoot, iPF=inputPathsFile), "{sR}/STRegionBoundaries.dat".format(sR=stealthEnv.stealthRoot)]
        formatted_iPF = (inputPathsFile.split("/"))[-1]

        startLine = 1
        endLine = 0
        while endLine < total_nLines:
            endLine = startLine + nFilesPerJob - 1
            isLastIteration = (endLine >= total_nLines)
            if isLastIteration:
                endLine = total_nLines
            processIdentifier = "selectionJob_{sT}_{y}_begin_{sL}_end_{eL}".format(sT=selectionType, y=year, sL=startLine, eL=endLine)
            jdlInterface = tmJDLInterface.tmJDLInterface(processName=processIdentifier, scriptPath="eventSelectionHelper.sh", outputDirectoryRelativePath="{cWAR}/selection{oI}".format(cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier)) # works even if "outputDirectoryRelativePath" is an absolute path
            jdlInterface.addFilesToTransferFromList(filesToTransfer)
            # Arguments for script:
            # Note: it seems simpler and certainly more readable to just include the "=" signs with the argument names, but I'm not sure whether that is allowed by JDL.
            jdlInterface.addScriptArgument("{iPF}".format(iPF=formatted_iPF)) # Argument 1: inputPathsFile
            jdlInterface.addScriptArgument("{iMC}".format(iMC=isMC)) # Argument 2: isMC
            jdlInterface.addScriptArgument("{dJS}".format(dJS=disableJetSelectionString)) # Argument 3: disableJetSelection
            jdlInterface.addScriptArgument("{sL}".format(sL=startLine)) # Argument 4: lineNumberStartInclusive
            jdlInterface.addScriptArgument("{eL}".format(eL=endLine)) # Argument 5: lineNumberEndInclusive
            jdlInterface.addScriptArgument("{y}".format(y=year)) # Argument 6: year

            # Other arguments:
            jdlInterface.addScriptArgument("{eP}".format(eP=stealthEnv.EOSPrefix)) # Argument 7: EOS prefix
            jdlInterface.addScriptArgument("{oD}{oI}".format(oD=inputArguments.outputDirectory_selections, oI=optional_identifier)) # Argument 8: selections output folder path
            jdlInterface.addScriptArgument("{oD}{oI}".format(oD=inputArguments.outputDirectory_statistics, oI=optional_identifier)) # Argument 9: statistics output folder path
            jdlInterface.addScriptArgument("{sT}".format(sT=selectionType)) # Argument 10: selection type

            if (stealthEnv.habitat == "lxplus"):
                jdlInterface.setFlavor("tomorrow")
            # Write JDL
            jdlInterface.writeToFile()

            submissionCommand = "cd {cWAR}/selection{oI}/ && condor_submit {pI}.jdl && cd {sR}".format(pI=processIdentifier, cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier, sR=stealthEnv.stealthRoot)
            print ("Generated command: {sC}".format(sC=submissionCommand))
            if (inputArguments.isProductionRun):
                os.system(submissionCommand)
                print ("Submitted.")
            else:
                print("Not submitting because isProductionRun flag was not set.")
            if not(inputArguments.preserveInputFileLists):
                os.system("echo \"{eP}/{oD}{oI}/selection_{t}_{y}_signal_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}_{y}{oI}_signal.txt".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, y=year, b=startLine, e=endLine))
                os.system("echo \"{eP}/{oD}{oI}/selection_{t}_{y}_control_fakefake_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}_{y}{oI}_control_fakefake.txt".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, y=year, b=startLine, e=endLine))
                os.system("echo \"{eP}/{oD}{oI}/selection_{t}_{y}_control_mediumfake_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}_{y}{oI}_control_mediumfake.txt".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, y=year, b=startLine, e=endLine))
                os.system("echo \"{eP}/{oD}{oI}/statistics_{t}_{y}_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_statistics_{t}_{y}{oI}.txt".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_statistics, oI=optional_identifier, t=selectionType, y=year, b=startLine, e=endLine))
            if isLastIteration: break
            startLine = 1+endLine
            if (startLine > total_nLines): break

os.system("rm -f submitEventSelectionJobs.lock")
