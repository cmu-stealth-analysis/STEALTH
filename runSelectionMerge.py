#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, signal
import tmMultiProcessLauncher # from tmPyUtils
import stealthEnv # from this folder

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Run event selection merging scripts.')
inputArgumentsParser.add_argument('--selectionsToRun', default="data,MC", help="Comma-separated list of selections to run. Allowed: \"data\", \"data_singlemedium\", \"data_jetHT\", \"MC\", \"MC_EMEnrichedQCD\", \"MC_GJet\", \"MC_QCD\", or \"MC_hgg\". For MC selections, disable HLT photon trigger and enable additional MC selection. Default is \"data,MC\".", type=str)
inputArgumentsParser.add_argument('--year', default="all", help="Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger (unless sample is MC), and the photon ID cuts which are based on year-dependent recommendations.", type=str)
inputArgumentsParser.add_argument('--disableJetSelection', action='store_true', help="Disable jet selection.")
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArguments = inputArgumentsParser.parse_args()

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)
mergeLogsDirectory = "{aR}/analysis{oi}/mergeLogs".format(aR=stealthEnv.analysisRoot, oi=optional_identifier)

multiProcessLauncher = None
def checkAndEstablishLock(): # Make sure that at most one instance is running at a time
    global multiProcessLauncher
    if (os.path.isfile("{mLD}/runSelectionMerge.lock".format(mLD=mergeLogsDirectory))):
        sys.exit("ERROR: only one instance of event merger can run at a time! If you're sure this is not a problem, remove this file: {mLD}/runSelectionMerge.lock".format(mLD=mergeLogsDirectory))
    else:
        os.system("mkdir -p {mLD} && touch {mLD}/runSelectionMerge.lock".format(mLD=mergeLogsDirectory))
        multiProcessLauncher = tmMultiProcessLauncher.tmMultiProcessLauncher(logOutputFolder=mergeLogsDirectory, printDebug=True)
checkAndEstablishLock()

def removeLock():
    global multiProcessLauncher
    multiProcessLauncher.killAll()
    os.system("rm -f {mLD}/runSelectionMerge.lock".format(mLD=mergeLogsDirectory))

def signal_handler(sig, frame):
    removeLock()
    sys.exit("Terminated by user.")
signal.signal(signal.SIGINT, signal_handler)

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
    elif (inputSelectionToRun == "MC_GJet"):
        selectionTypesToRun.append("MC_GJet")
    elif (inputSelectionToRun == "MC_QCD"):
        selectionTypesToRun.append("MC_QCD")
    elif (inputSelectionToRun == "MC_hgg"):
        selectionTypesToRun.append("MC_hgg")
    else:
        removeLock()
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
    removeLock()
    sys.exit("ERROR: invalid value for argument \"year\": {v}".format(v=inputArguments.year))

stealthEnv.execute_in_env("eos {eP} mkdir -p {sER}/selections/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), functionToCallIfCommandExitsWithError=removeLock)
stealthEnv.execute_in_env("eos {eP} mkdir -p {sER}/statistics/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), functionToCallIfCommandExitsWithError=removeLock)

for selectionType in selectionTypesToRun:
    isMC = True
    isMCString = "true"
    if (("data" in selectionType)
        or (selectionType == "MC_EMEnrichedQCD")
        or (selectionType == "MC_GJet")
        or (selectionType == "MC_QCD")
        or (selectionType == "MC_hgg")
    ):
        isMC = False
        isMCString = "false"
    for year in yearsToRun:
        mergeStatistics = True
        if ((selectionType == "MC_QCD") or (selectionType == "MC_EMEnrichedQCD") or (selectionType == "MC_hgg")):
            if (year != 2017): # The only reason we need these is to calculate ID efficiencies
                mergeStatistics = False
        if (selectionType == "MC_GJet"):
            if (year != 2016): # The only reason we need these is to calculate ID efficiencies
                mergeStatistics = False
        if (selectionType == "data_singlemedium"):
            mergeStatistics = False
        if mergeStatistics:
            inputFilesList_statistics = "fileLists/inputFileList_statistics_{t}{aJS}_{y}{oI}.txt".format(oI=optional_identifier, t=selectionType, aJS=allJetsString, y=year)
            mergeStatisticsCommand = "eventSelection/bin/mergeStatisticsHistograms inputFilesList={iFL} outputFolder={oF} outputFileName={oFN} isMC={iMCS}".format(iFL=inputFilesList_statistics, oF="{eP}/{sER}/statistics/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), oFN="merged_statistics_{t}{aJS}_{y}.root".format(t=selectionType, aJS=allJetsString, y=year), iMCS=isMCString)
            multiProcessLauncher.spawn(shellCommands=mergeStatisticsCommand, logFileName="mergeLog_statistics_{t}{aJS}_{y}.log".format(t=selectionType, aJS=allJetsString, y=year), printDebug=True)
        for selectionRegion in ["signal", "signal_loose", "control_fakefake", "control_singlemedium"]:
            selectionRegionString = "{sR}".format(sR=selectionRegion)
            if (selectionRegion == "control_fakefake"): selectionRegionString = "control"
            mergeSelection = True
            if (selectionRegion == "control_singlemedium"):
                mergeSelection = False
                if ((selectionType == "data_singlemedium") and not(isMC)): mergeSelection = True
            else:
                if (selectionType == "data_singlemedium"): mergeSelection = False
            if ((selectionType == "data_jetHT") or (selectionType == "MC_QCD") or (selectionType == "MC_hgg") or (selectionType == "MC_EMEnrichedQCD")): mergeSelection = False
            if ((selectionType == "MC_GJet") and (year != 2016)): mergeSelection = False
            if not(mergeSelection): continue
            inputFilesList_selection = "fileLists/inputFileList_selections_{t}{aJS}_{y}{oI}_{r}.txt".format(oI=optional_identifier, t=selectionType, aJS=allJetsString, y=year, r=selectionRegion)
            mergeSelectionCommand = "eventSelection/bin/mergeEventSelections inputFilesList={iFL} outputFolder={oF} outputFileName={oFN}".format(iFL=inputFilesList_selection, oF="{eP}/{sER}/selections/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), oFN="merged_selection_{t}{aJS}_{y}_{sRS}.root".format(t=selectionType, aJS=allJetsString, y=year, sRS=selectionRegionString))
            multiProcessLauncher.spawn(shellCommands=mergeSelectionCommand, logFileName="mergeLog_selection_{t}{aJS}_{y}_{sRS}.log".format(t=selectionType, aJS=allJetsString, y=year, sRS=selectionRegionString), printDebug=True)
multiProcessLauncher.monitorToCompletion()
removeLock()
