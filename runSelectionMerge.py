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
inputArgumentsParser.add_argument('--selectionsToRun', default="data,MC", help="Comma-separated list of selections to run. Allowed: \"data\", \"data_singlemedium\", \"data_jetHT\", \"MC\", \"MC_EMEnrichedQCD\", or \"MC_QCD\". For MC selections, disable HLT photon trigger and enable additional MC selection. Default is \"data,MC\".", type=str)
inputArgumentsParser.add_argument('--year', default="all", help="Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger (unless sample is MC), and the photon ID cuts which are based on year-dependent recommendations.", type=str)
inputArgumentsParser.add_argument('--disableJetSelection', action='store_true', help="Disable jet selection.")
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArguments = inputArgumentsParser.parse_args()

def execute_in_env(commandToRun, printDebug=False):
    env_setup_command = "bash -c \"cd {sR} && source setupEnv.sh".format(sR=stealthRoot)
    runInEnv = "{e_s_c} && set -x && {c} && set +x\"".format(e_s_c=env_setup_command, c=commandToRun)
    if (printDebug):
        print("About to execute command:")
        print("{c}".format(c=runInEnv))
    os.system(runInEnv)

def spawnMerge(scriptPath, inputFilesList, outputFolder, outputFileName, extraArguments=""):
    tmp_outputFileName = outputFileName
    if (outputFileName[-5:] == ".root"):
        tmp_outputFileName = outputFileName[:-5]
    logFileName = "mergeLog_{oFN}.txt".format(oFN=tmp_outputFileName)
    shellCommand = ""
    shellCommand += "{sP} inputFilesList={iFL} outputFolder={oF} outputFileName={oFN}".format(sP=scriptPath, iFL=inputFilesList, oF=outputFolder, oFN=outputFileName)
    if not(extraArguments == ""):
        shellCommand += (" " + extraArguments)
    shellCommand += " > mergeLogs/{lFN} 2>&1".format(lFN=logFileName)
    print("About to spawn: {c}".format(c=shellCommand))
    process = subprocess.Popen("{c}".format(c=shellCommand), shell=True)
    return (logFileName, process)

def monitor(processes):
    current_runningProcessesLogsList = [key for key in processes.keys()]
    print("current_runningProcessesLogsList: {rPLL}".format(rPLL=current_runningProcessesLogsList))
    while True:
        time.sleep(10)
        for i in range(0, 5): print("\n")
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
                for i in range(0, 2): print("\n")
                next_runningProcessesLogsList.append(logFileName)
        current_runningProcessesLogsList = [logFileName for logFileName in next_runningProcessesLogsList]

def killAll(processes):
    runningProcessesLogsList = [key for key in processes.keys()]
    for log in runningProcessesLogsList:
        (processes[log]).kill()

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)

allJetsString = ""
if (inputArguments.disableJetSelection):
    allJetsString = "_allJets"

selectionTypesToRun = []
for inputSelectionToRun in (inputArguments.selectionsToRun.split(",")):
    if (inputSelectionToRun == "data"):
        selectionTypesToRun.append("data")
    elif (inputSelectionToRun == "data_singlemedium"):
        selectionTypesToRun.append("data_singlemedium")
    elif (inputSelectionToRun == "data_jetHT"):
        selectionTypesToRun.append("data_jetHT")
    elif (inputSelectionToRun == "MC"):
        selectionTypesToRun.append("MC_stealth_t5")
        selectionTypesToRun.append("MC_stealth_t6")
    elif (inputSelectionToRun == "MC_EMEnrichedQCD"):
        selectionTypesToRun.append("MC_EMEnrichedQCD")
    elif (inputSelectionToRun == "MC_QCD"):
        selectionTypesToRun.append("MC_QCD")
    else:
        sys.exit("ERROR: invalid value for argument \"selectionsToRun\": {v}".format(v=inputSelectionToRun))

yearsToRun = []
if (inputArguments.year == "2016"):
    yearsToRun.append(2016)
elif (inputArguments.year == "2017"):
    yearsToRun.append(2017)
elif (inputArguments.year == "2018"):
    yearsToRun.append(2018)
elif (inputArguments.year == "all"):
    yearsToRun.append(2016)
    yearsToRun.append(2017)
    yearsToRun.append(2018)
else:
    sys.exit("ERROR: invalid value for argument \"year\": {v}".format(v=inputArguments.year))

execute_in_env("eos {eP} mkdir -p {sER}/selections/combined_DoublePhoton{oI}".format(eP=EOSPrefix, sER=stealthEOSRoot, oI=optional_identifier), printDebug=True)
execute_in_env("eos {eP} mkdir -p {sER}/statistics/combined_DoublePhoton{oI}".format(eP=EOSPrefix, sER=stealthEOSRoot, oI=optional_identifier), printDebug=True)

processes = {}
for selectionType in selectionTypesToRun:
    isMC = True
    isMCString = "true"
    if (("data" in selectionType) or ((selectionType == "MC_EMEnrichedQCD") or (selectionType == "MC_QCD"))):
        isMC = False
        isMCString = "false"
    for year in yearsToRun:
        if ((selectionType == "MC_QCD") or (selectionType == "MC_EMEnrichedQCD")):
            if (year != 2017): continue # The only reason we need these is to calculate ID efficiencies
        if not(selectionType == "data_singlemedium"):
            inputFilesList_statistics = "fileLists/inputFileList_statistics_{t}{aJS}_{y}{oI}.txt".format(oI=optional_identifier, t=selectionType, aJS=allJetsString, y=year)
            print("Spawning statistics merge job for year={y}, selection type={t}".format(t=selectionType, y=year))
            processTuple_statistics = spawnMerge(scriptPath="eventSelection/bin/mergeStatisticsHistograms", inputFilesList=inputFilesList_statistics, outputFolder="{eP}/{sER}/statistics/combined_DoublePhoton{oI}".format(eP=EOSPrefix, sER=stealthEOSRoot, oI=optional_identifier), outputFileName="merged_statistics_{t}{aJS}_{y}.root".format(t=selectionType, aJS=allJetsString, y=year), extraArguments="isMC={iMC}".format(iMC=isMCString))
            if (processTuple_statistics[0] in processes.keys()):
                killAll(processes)
                sys.exit("ERROR: found duplicate: {k}".format(k=processTuple_statistics[0]))
            processes[processTuple_statistics[0]] = processTuple_statistics[1]
        for selectionRegion in ["signal", "signal_loose", "control_fakefake", "control_singlemedium"]:
            selectionRegionString = "{sR}".format(sR=selectionRegion)
            if (selectionRegion == "control_fakefake"): selectionRegionString = "control"
            run_merge = True
            if (selectionRegion == "control_singlemedium"):
                run_merge = False
                if ((selectionType == "data_singlemedium") and not(isMC)): run_merge = True
            else:
                if (selectionType == "data_singlemedium"): run_merge = False
            if (selectionType == "data_jetHT"): run_merge = False
            if not(run_merge): continue
            inputFilesList_eventMerge = "fileLists/inputFileList_selections_{t}{aJS}_{y}{oI}_{r}.txt".format(oI=optional_identifier, t=selectionType, aJS=allJetsString, y=year, r=selectionRegion)
            print("Spawning merge job for year={y}, selection type={t}, selection region={r}".format(y=year, t=selectionType, r=selectionRegion))
            processTuple_eventMerge = spawnMerge(scriptPath="eventSelection/bin/mergeEventSelections", inputFilesList=inputFilesList_eventMerge, outputFolder="{eP}/{sER}/selections/combined_DoublePhoton{oI}".format(eP=EOSPrefix, sER=stealthEOSRoot, oI=optional_identifier), outputFileName="merged_selection_{t}{aJS}_{y}_{sRS}.root".format(t=selectionType, aJS=allJetsString, y=year, sRS=selectionRegionString))
            if (processTuple_eventMerge[0] in processes.keys()):
                killAll(processes)
                sys.exit("ERROR: found duplicate: {k}".format(k=processTuple_eventMerge[0]))
            processes[processTuple_eventMerge[0]] = processTuple_eventMerge[1]

monitor(processes)

os.system("rm -f runSelectionMerge.lock")
