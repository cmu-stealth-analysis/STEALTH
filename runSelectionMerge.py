#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, signal, re
import tmMultiProcessLauncher # from tmPyUtils
import stealthEnv # from this folder

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Run event selection merging scripts.')
inputArgumentsParser.add_argument('--selectionsToRun', default="data", help="Comma-separated list of selections to run. Allowed: \"data\", \"data_singlephoton\", \"data_jetHT\", \"MC\", \"MC_DiPhotonJets\", \"MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)(16|17|18)(|_singlephoton)\", or \"MC_hgg\".", type=str)
inputArgumentsParser.add_argument('--year', default="all", help="Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger (unless sample is MC), and the photon ID cuts which are based on year-dependent recommendations.", type=str)
inputArgumentsParser.add_argument('--disablePhotonSelection', action='store_true', help="Disable photon selection.")
inputArgumentsParser.add_argument('--disableJetSelection', action='store_true', help="Disable jet selection.")
inputArgumentsParser.add_argument('--invertElectronVeto', action='store_true', help="Invert electron veto.")
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArguments = inputArgumentsParser.parse_args()

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)
mergeLogsDirectory = "{aR}/analysis{oi}/mergeLogs".format(aR=stealthEnv.analysisRoot, oi=optional_identifier)

n_subsamples = {
    "MC_EMEnrichedGJetPt16": 3,
    "MC_EMEnrichedGJetPt17": 3,
    "MC_EMEnrichedGJetPt18": 3,
    "MC_HighHTQCD16": 7,
    "MC_HighHTQCD17": 8,
    "MC_HighHTQCD18": 8,
    "MC_GJetHT16": 5,
    "MC_GJetHT17": 5,
    "MC_GJetHT18": 5
}
year_dependent_mc_labels = list(n_subsamples.keys())
for year_dependent_mc_label in year_dependent_mc_labels:
    n_subsamples["{k}_singlephoton".format(k=year_dependent_mc_label)] = n_subsamples[year_dependent_mc_label]

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
    if not(multiProcessLauncher is None): multiProcessLauncher.killAll()
    os.system("rm -f {mLD}/runSelectionMerge.lock".format(mLD=mergeLogsDirectory))

def signal_handler(sig, frame):
    removeLock()
    multiProcessLauncher.killAll()
    sys.exit("Terminated by user.")
signal.signal(signal.SIGINT, signal_handler)

photonSelectionString = ""
if (inputArguments.disablePhotonSelection):
    photonSelectionString = "_noPhotonSelection"

jetSelectionString = ""
if (inputArguments.disableJetSelection):
    jetSelectionString = "_noJetSelection"

electronVetoString = ""
if (inputArguments.invertElectronVeto):
    electronVetoString = "_invertElectronVeto"

overallIdentificationString = "{pSS}{jSS}{eVS}".format(pSS=photonSelectionString, jSS=jetSelectionString, eVS=electronVetoString)

selectionTypesToRun = []
selectionTypesToRun_Step2 = []
for inputSelectionToRun in (inputArguments.selectionsToRun.split(",")):
    if (inputSelectionToRun == "data"):
        selectionTypesToRun.append("data")
    elif (inputSelectionToRun == "data_singlephoton"):
        selectionTypesToRun.append("data_singlephoton")
    elif (inputSelectionToRun == "data_jetHT"):
        selectionTypesToRun.append("data_jetHT")
    elif (inputSelectionToRun == "MC"):
        selectionTypesToRun.append("MC_stealth_t5")
        selectionTypesToRun.append("MC_stealth_t6")
    elif (inputSelectionToRun == "MC_hgg"):
        selectionTypesToRun.append("MC_hgg")
    elif (inputSelectionToRun == "MC_DiPhotonJets"):
        selectionTypesToRun.append("MC_DiPhotonJets")
    elif (inputSelectionToRun == "MC_DiPhotonJets_singlephoton"):
        selectionTypesToRun.append("MC_DiPhotonJets_singlephoton")
    else:
        MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)(|_singlephoton)$", inputSelectionToRun)
        if MCBKGMatch:
            full_match = MCBKGMatch.group(0)
            dataset_id = MCBKGMatch.group(1)
            year_last_two_digits_str = MCBKGMatch.group(2)
            singlephoton_match = MCBKGMatch.group(3)
            for index_subsample in range(n_subsamples[full_match], 0, -1):
                selectionTypesToRun.append("{m}_{i}".format(m=full_match, i=index_subsample))
            if (singlephoton_match == "_singlephoton"):
                for selectionRegion in ["control_singlemedium", "control_singleloose", "control_singlefake"]:
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_{m}{oIS}_20{y2}{oI}_{r}.txt".format(m=full_match, oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion, y2=year_last_two_digits_str)
                    os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
                os.system("rm -f fileLists/inputFileList_step2Merge_statistics_{m}{oIS}_20{y2}{oI}.txt && touch fileLists/inputFileList_step2Merge_statistics_{m}{oIS}_20{y2}{oI}.txt".format(m=full_match, oI=optional_identifier, oIS=overallIdentificationString, y2=year_last_two_digits_str))
            else:
                for selectionRegion in ["signal", "signal_loose", "control_fakefake"]:
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_{m}{oIS}_20{y2}{oI}_{r}.txt".format(m=full_match, oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion, y2=year_last_two_digits_str)
                    os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
                os.system("rm -f fileLists/inputFileList_step2Merge_statistics_{m}{oIS}_20{y2}{oI}.txt && touch fileLists/inputFileList_step2Merge_statistics_{m}{oIS}_20{y2}{oI}.txt".format(m=full_match, oI=optional_identifier, oIS=overallIdentificationString, y2=year_last_two_digits_str))
            selectionTypesToRun_Step2.append(full_match)
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

stealthEnv.execute_in_env(commandToRun="eos {eP} mkdir -p {sER}/selections/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), functionToCallIfCommandExitsWithError=removeLock)
stealthEnv.execute_in_env(commandToRun="eos {eP} mkdir -p {sER}/statistics/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), functionToCallIfCommandExitsWithError=removeLock)

# Step 1 merge: sufficient for everything except year-dependent MC datasets
filesToCleanup = []
for selectionType in selectionTypesToRun:
    isMC = ((selectionType == "MC_stealth_t5") or (selectionType == "MC_stealth_t6"))
    isMCString = None
    if isMC: isMCString = "true"
    else: isMCString = "false"
    for year in yearsToRun:
        mergeStatistics = True
        MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)(|_singlephoton)_([0-9]*)$", selectionType)
        if MCBKGMatch:
            year_last_two_digits_str = MCBKGMatch.group(2)
            year_from_match = 2000 + int(0.5 + float(year_last_two_digits_str))
            if not(year_from_match == year): mergeStatistics = False
        if ((selectionType == "data_singlephoton") or
            (bool(re.match(r"^MC_DiPhotonJets_singlephoton$", selectionType))) or
            (bool(re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)_singlephoton_([0-9]*)$", selectionType)))):
            mergeStatistics = False
        if mergeStatistics:
            inputFilesList_statistics = "fileLists/inputFileList_statistics_{t}{oIS}_{y}{oI}.txt".format(oI=optional_identifier, t=selectionType, oIS=overallIdentificationString, y=year)
            outputFilePath = "merged_statistics_{t}{oIS}_{y}.root".format(t=selectionType, oIS=overallIdentificationString, y=year)
            outputFolder = "{eP}/{sER}/statistics/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier)
            mergeStatisticsCommand = "eventSelection/bin/mergeStatisticsHistograms inputFilesList={iFL} outputFolder={oF} outputFileName={oFP} isMC={iMCS}".format(iFL=inputFilesList_statistics, oF=outputFolder, oFP=outputFilePath, iMCS=isMCString)
            MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)_([0-9]*)$", selectionType)
            if MCBKGMatch:
                full_match = MCBKGMatch.group(0)
                dataset_id = MCBKGMatch.group(1)
                year_last_two_digits_str = MCBKGMatch.group(2)
                mergeStep2FilePath = "fileLists/inputFileList_step2Merge_statistics_MC_{did}{y2}{oIS}_20{y2}{oI}.txt".format(did=dataset_id, oI=optional_identifier, oIS=overallIdentificationString, y2=year_last_two_digits_str)
                os.system("echo {oF}/{oFP} >> {mS2FP}".format(oF=outputFolder, oFP=outputFilePath, mS2FP=mergeStep2FilePath))
                filesToCleanup.append("{sER}/statistics/combined_DoublePhoton{oI}/{oFP}".format(sER=stealthEnv.stealthEOSRoot, oI=optional_identifier, oFP=outputFilePath))
            multiProcessLauncher.spawn(shellCommands=mergeStatisticsCommand, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="mergeLog_statistics_{t}{oIS}_{y}.log".format(t=selectionType, oIS=overallIdentificationString, y=year), printDebug=True)
        for selectionRegion in ["signal", "signal_loose", "control_fakefake", "control_singlemedium", "control_singleloose", "control_singlefake"]:
            selectionRegionString = selectionRegion
            if (selectionRegion == "control_fakefake"): selectionRegionString = "control"
            mergeSelection = True
            if ((selectionRegion == "control_singlemedium") or (selectionRegion == "control_singleloose") or (selectionRegion == "control_singlefake")):
                mergeSelection = False
                if ((selectionType == "data_singlephoton") or
                    (bool(re.match(r"^MC_DiPhotonJets_singlephoton$", selectionType))) or
                    (bool(re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)_singlephoton_([0-9]*)$", selectionType)))): mergeSelection = True
            else:
                if ((selectionType == "data_singlephoton") or
                    (bool(re.match(r"^MC_DiPhotonJets_singlephoton$", selectionType))) or
                    (bool(re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)_singlephoton_([0-9]*)$", selectionType)))): mergeSelection = False
            if ((selectionType == "data_jetHT") or (selectionType == "MC_hgg")): mergeSelection = False
            MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)(|_singlephoton)_([0-9]*)$", selectionType)
            if MCBKGMatch:
                year_last_two_digits_str = MCBKGMatch.group(2)
                year_from_match = 2000 + int(0.5 + float(year_last_two_digits_str))
                if not(year_from_match == year): mergeSelection = False
            if not(mergeSelection): continue
            inputFilesList_selection = "fileLists/inputFileList_selections_{t}{oIS}_{y}{oI}_{r}.txt".format(oI=optional_identifier, t=selectionType, oIS=overallIdentificationString, y=year, r=selectionRegion)
            outputFolder = "{eP}/{sER}/selections/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier)
            outputFilePath = "merged_selection_{t}{oIS}_{y}_{sRS}.root".format(t=selectionType, oIS=overallIdentificationString, y=year, sRS=selectionRegionString)
            mergeSelectionCommand = "eventSelection/bin/mergeEventSelections inputFilesList={iFL} outputFolder={oF} outputFileName={oFP}".format(iFL=inputFilesList_selection, oF=outputFolder, oFP=outputFilePath)
            MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)(|_singlephoton)_([0-9]*)$", selectionType)
            if MCBKGMatch:
                full_match = MCBKGMatch.group(0)
                dataset_id = MCBKGMatch.group(1)
                year_last_two_digits_str = MCBKGMatch.group(2)
                singlephoton_match = MCBKGMatch.group(3)
                mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_{did}{y2}{spm}{oIS}_20{y2}{oI}_{r}.txt".format(oI=optional_identifier, did=dataset_id, spm=singlephoton_match, oIS=overallIdentificationString, r=selectionRegion, y2=year_last_two_digits_str)
                os.system("echo {oF}/{oFP} >> {mS2FP}".format(oF=outputFolder, oFP=outputFilePath, mS2FP=mergeStep2FilePath))
                filesToCleanup.append("{sER}/selections/combined_DoublePhoton{oI}/{oFP}".format(sER=stealthEnv.stealthEOSRoot, oI=optional_identifier, oFP=outputFilePath))
            multiProcessLauncher.spawn(shellCommands=mergeSelectionCommand, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="mergeLog_selection_{t}{oIS}_{y}_{sRS}.log".format(t=selectionType, oIS=overallIdentificationString, y=year, sRS=selectionRegionString), printDebug=True)
multiProcessLauncher.monitorToCompletion()

# Step 2 merge: for year-dependent datasets
monitoringNeeded = False
for selectionType in selectionTypesToRun_Step2:
    isMC = ((selectionType == "MC_stealth_t5") or (selectionType == "MC_stealth_t6"))
    isMCString = None
    if isMC: isMCString = "true"
    else: isMCString = "false"
    for year in yearsToRun:
        mergeStatistics = True
        MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)$", selectionType)
        if MCBKGMatch:
            year_last_two_digits_str = MCBKGMatch.group(2)
            year_from_match = 2000 + int(0.5 + float(year_last_two_digits_str))
            if not(year_from_match == year): mergeStatistics = False
        if (bool(re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)_singlephoton$", selectionType))):
            mergeStatistics = False
        if mergeStatistics:
            inputFilesList_statistics = "fileLists/inputFileList_step2Merge_statistics_{t}{oIS}_{y}{oI}.txt".format(t=selectionType, y=year, oI=optional_identifier, oIS=overallIdentificationString)
            outputFilePath = "merged_statistics_{t}{oIS}_{y}.root".format(t=selectionType, oIS=overallIdentificationString, y=year)
            mergeStatisticsCommand = "eventSelection/bin/mergeStatisticsHistograms inputFilesList={iFL} outputFolder={oF} outputFileName={oFP} isMC={iMCS}".format(iFL=inputFilesList_statistics, oF="{eP}/{sER}/statistics/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), oFP=outputFilePath, iMCS=isMCString)
            multiProcessLauncher.spawn(shellCommands=mergeStatisticsCommand, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="mergeLog_statistics_{t}{oIS}_{y}.log".format(t=selectionType, oIS=overallIdentificationString, y=year), printDebug=True)
            monitoringNeeded = True
        for selectionRegion in ["signal", "signal_loose", "control_fakefake", "control_singlemedium", "control_singleloose", "control_singlefake"]:
            selectionRegionString = "{sR}".format(sR=selectionRegion)
            if (selectionRegion == "control_fakefake"): selectionRegionString = "control"
            mergeSelection = True
            if ((selectionRegion == "control_singlemedium") or (selectionRegion == "control_singleloose") or (selectionRegion == "control_singlefake")):
                mergeSelection = False
                if (bool(re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)_singlephoton$", selectionType))):
                    mergeSelection = True
            else:
                if (bool(re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)_singlephoton$", selectionType))):
                    mergeSelection = False
            MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)(|_singlephoton)$", selectionType)
            if MCBKGMatch:
                year_last_two_digits_str = MCBKGMatch.group(2)
                year_from_match = 2000 + int(0.5 + float(year_last_two_digits_str))
                if not(year_from_match == year): mergeSelection = False
            if not(mergeSelection): continue
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_{t}{oIS}_{y}{oI}_{r}.txt".format(t=selectionType, y=year, oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
            outputFolder = "{eP}/{sER}/selections/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier)
            outputFilePath = "merged_selection_{t}{oIS}_{y}_{sRS}.root".format(t=selectionType, oIS=overallIdentificationString, y=year, sRS=selectionRegionString)
            mergeSelectionCommand = "eventSelection/bin/mergeEventSelections inputFilesList={mS2FP} outputFolder={oF} outputFileName={oFP}".format(mS2FP=mergeStep2FilePath, oF=outputFolder, oFP=outputFilePath)
            multiProcessLauncher.spawn(shellCommands=mergeSelectionCommand, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="mergeLog_selection_{t}{oIS}_{y}_{sRS}.log".format(t=selectionType, oIS=overallIdentificationString, y=year, sRS=selectionRegionString), printDebug=True)
            monitoringNeeded = True
if monitoringNeeded: multiProcessLauncher.monitorToCompletion()

for fileToCleanup in filesToCleanup:
    stealthEnv.execute_in_env(commandToRun="eos {eP} rm {fTC}".format(eP=stealthEnv.EOSPrefix, fTC=fileToCleanup), functionToCallIfCommandExitsWithError=removeLock)

removeLock()
