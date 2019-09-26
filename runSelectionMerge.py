#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, re, subprocess, time

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
if (os.path.isfile("runSelectionMerge.lock")):
    sys.exit("ERROR: only one instance of event merger can run at a time!")
else:
    os.system("touch runSelectionMerge.lock")

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Run event selection merging scripts.')
inputArgumentsParser.add_argument('--runOnlyDataOrMC', default="all", help="Takes values \"data\" or \"MC\" if selection is to be run only on data or only on MC samples. For MC selections, disable HLT photon trigger and enable additional MC selection. Default is \"all\", which means both selections are run.", type=str)
inputArgumentsParser.add_argument('--year', default="all", help="Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger (unless sample is MC), and the photon ID cuts which are based on year-dependent recommendations.", type=str)
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArguments = inputArgumentsParser.parse_args()

def execute_in_env(commandToRun, printDebug=False):
    env_setup_command = "bash -c \"cd {sR} && source setupEnv.sh".format(sR=stealthRoot)
    runInEnv = "{e_s_c} && set -x && {c} && set +x\"".format(e_s_c=env_setup_command, c=commandToRun)
    if (printDebug):
        print("About to execute command:")
        print("{c}".format(c=runInEnv))
    os.system(runInEnv)

def spawnMerge(inputFilesList, outputFolder, outputFileName):
    tmp_outputFileName = outputFileName
    if (outputFileName[-5:] == ".root"):
        tmp_outputFileName = outputFileName[:-5]
    logFileName = "mergeLog_{oFN}.txt".format(oFN=tmp_outputFileName)
    shellCommand = "eventSelection/bin/mergeEventSelections inputFilesList={iFL} outputFolder={oF} outputFileName={oFN} > mergeLogs/{lFN} 2>&1".format(iFL=inputFilesList, oF=outputFolder, oFN=outputFileName, lFN=logFileName)
    print("About to spawn: {c}".format(c=shellCommand))
    process = subprocess.Popen("{c}".format(c=shellCommand), shell=True)
    return (logFileName, process)

def monitor(processes):
    current_runningProcessesLogsList = [key for key in processes.keys()]
    print("current_runningProcessesLogsList: {rPLL}".format(rPLL=current_runningProcessesLogsList))
    while True:
        time.sleep(10)
        for i in range(0, 10): print("\n")
        if (len(current_runningProcessesLogsList) == 0): break
        next_runningProcessesLogsList = []
        for runningProcessLogsCounter in range(0, len(current_runningProcessesLogsList)):
            logFileName = current_runningProcessesLogsList[runningProcessLogsCounter]
            # print("Checking: {lFN}".format(lFN=logFileName))
            if ((processes[logFileName] is None) or not((processes[logFileName].poll() is None))):
                print("Done: mergeLogs/{lFN}".format(lFN=logFileName))
                os.system("tail -20 mergeLogs/{lFN}".format(lFN=logFileName))
            else: # process hasn't terminated
                print("Output of mergeLogs/{lFN}:".format(lFN=logFileName))
                os.system("tail -2 mergeLogs/{lFN}".format(lFN=logFileName))
                for i in range(0, 3): print("\n")
                next_runningProcessesLogsList.append(logFileName)
        current_runningProcessesLogsList = [logFileName for logFileName in next_runningProcessesLogsList]

def killAll(processes):
    runningProcessesLogsList = [key for key in processes.keys()]
    for log in runningProcessesLogsList:
        (processes[log]).kill()

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)

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

execute_in_env("eos {eP} mkdir -p {sER}/combined_DoublePhoton{oI}".format(eP=EOSPrefix, sER=stealthEOSRoot, oI=optional_identifier), printDebug=True)

processes = {}
for selectionType in selectionTypesToRun:
    for year in yearsToRun:
        for selectionRegion in ["signal", "control_fakefake", "control_mediumfake"]:
            inputFilesList = "fileLists/inputFileList_selections_{t}_{y}{oI}_{r}.txt".format(oI=optional_identifier, t=selectionType, y=year, r=selectionRegion)
            print("Spawning merge job for year={y}, selection type={t}, selection region={r}".format(y=year, t=selectionType, r=selectionRegion))
            processTuple=spawnMerge(inputFilesList=inputFilesList, outputFolder="{eP}/{sER}/selections/combined_DoublePhoton{oI}".format(eP=EOSPrefix, sER=stealthEOSRoot, oI=optional_identifier), outputFileName="merged_selection_{t}_{y}_{r}.root".format(t=selectionType, y=year, r=selectionRegion))
            if (processTuple[0] in processes.keys()):
                killAll(processes)
                sys.exit("ERROR: found duplicate: {k}".format(k=processTuple[0]))
            processes[processTuple[0]] = processTuple[1]

monitor(processes)

os.system("rm -f runSelectionMerge.lock")
