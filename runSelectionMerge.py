#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, signal, re
import tmMultiProcessLauncher # from tmPyUtils
import stealthEnv # from this folder

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Run event selection merging scripts.')
inputArgumentsParser.add_argument('--selectionsToRun', default="data,MC,MC_GJet", help="Comma-separated list of selections to run. Allowed: \"data\", \"data_singlephoton\", \"data_jetHT\", \"MC\", \"MC_EMEnrichedQCD\", \"MC_GJet\", \"MC_GJet_singlephoton\", \"MC_QCD\", \"MC_QCD_singlephoton\", or \"MC_hgg\". For MC selections, disable HLT photon trigger and enable additional MC selection. Default is \"data,MC,MC_GJet\".", type=str)
inputArgumentsParser.add_argument('--year', default="all", help="Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger (unless sample is MC), and the photon ID cuts which are based on year-dependent recommendations.", type=str)
inputArgumentsParser.add_argument('--disableJetSelection', action='store_true', help="Disable jet selection.")
inputArgumentsParser.add_argument('--invertElectronVeto', action='store_true', help="Invert electron veto.")
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
    if not(multiProcessLauncher is None): multiProcessLauncher.killAll()
    os.system("rm -f {mLD}/runSelectionMerge.lock".format(mLD=mergeLogsDirectory))

def signal_handler(sig, frame):
    removeLock()
    multiProcessLauncher.killAll()
    sys.exit("Terminated by user.")
signal.signal(signal.SIGINT, signal_handler)

allJetsString = ""
if (inputArguments.disableJetSelection):
    allJetsString = "_allJets"

electronVetoString = ""
if (inputArguments.invertElectronVeto):
    electronVetoString = "_invertElectronVeto"

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
    elif (inputSelectionToRun == "MC_EMEnrichedQCD"):
        selectionTypesToRun.append("MC_EMEnrichedQCD")
    elif (inputSelectionToRun == "MC_GJet"):
        # selectionTypesToRun.append("MC_GJet")
        selectionTypesToRun.append("MC_GJet5")
        selectionTypesToRun.append("MC_GJet4")
        selectionTypesToRun.append("MC_GJet3")
        selectionTypesToRun.append("MC_GJet2")
        selectionTypesToRun.append("MC_GJet1")
        for selectionRegion in ["signal", "signal_loose", "control_fakefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_GJet{aJS}{eVS}_2016{oI}_{r}.txt".format(oI=optional_identifier, aJS=allJetsString, eVS=electronVetoString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        selectionTypesToRun_Step2.append("MC_GJet")
    elif (inputSelectionToRun == "MC_GJet_singlephoton"):
        # selectionTypesToRun.append("MC_GJet_singlephoton")
        selectionTypesToRun.append("MC_GJet_singlephoton5")
        selectionTypesToRun.append("MC_GJet_singlephoton4")
        selectionTypesToRun.append("MC_GJet_singlephoton3")
        selectionTypesToRun.append("MC_GJet_singlephoton2")
        selectionTypesToRun.append("MC_GJet_singlephoton1")
        for selectionRegion in ["control_singlemedium", "control_singleloose", "control_singlefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_GJet_singlephoton{aJS}{eVS}_2016{oI}_{r}.txt".format(oI=optional_identifier, aJS=allJetsString, eVS=electronVetoString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        selectionTypesToRun_Step2.append("MC_GJet_singlephoton")
    elif (inputSelectionToRun == "MC_QCD"):
        # selectionTypesToRun.append("MC_QCD")
        selectionTypesToRun.append("MC_QCD6")
        selectionTypesToRun.append("MC_QCD5")
        selectionTypesToRun.append("MC_QCD4")
        selectionTypesToRun.append("MC_QCD3")
        selectionTypesToRun.append("MC_QCD2")
        selectionTypesToRun.append("MC_QCD1")
        for selectionRegion in ["signal", "signal_loose", "control_fakefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_QCD{aJS}{eVS}_2017{oI}_{r}.txt".format(oI=optional_identifier, aJS=allJetsString, eVS=electronVetoString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        selectionTypesToRun_Step2.append("MC_QCD")
    elif (inputSelectionToRun == "MC_QCD_singlephoton"):
        # selectionTypesToRun.append("MC_QCD_singlephoton")
        selectionTypesToRun.append("MC_QCD_singlephoton6")
        selectionTypesToRun.append("MC_QCD_singlephoton5")
        selectionTypesToRun.append("MC_QCD_singlephoton4")
        selectionTypesToRun.append("MC_QCD_singlephoton3")
        selectionTypesToRun.append("MC_QCD_singlephoton2")
        selectionTypesToRun.append("MC_QCD_singlephoton1")
        for selectionRegion in ["control_singlemedium", "control_singleloose", "control_singlefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_QCD_singlephoton{aJS}{eVS}_2017{oI}_{r}.txt".format(oI=optional_identifier, aJS=allJetsString, eVS=electronVetoString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        selectionTypesToRun_Step2.append("MC_QCD_singlephoton")
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

effectiveLuminosities = {
    # DAS query for MC_GJet: dataset dataset=/GJets_DR-0p4_HT-*_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_qcut19_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
    "MC_GJet1": 0.05745,
    "MC_GJet_singlephoton1": 0.05745,
    "MC_GJet2": 0.1865,
    "MC_GJet_singlephoton2": 0.1865,
    "MC_GJet3": 0.849,
    "MC_GJet_singlephoton3": 0.849,
    "MC_GJet4": 7.588,
    "MC_GJet_singlephoton4": 7.588,
    "MC_GJet5": 22.59,
    "MC_GJet_singlephoton5": 22.59,
    # DAS query for MC_QCD: dataset dataset=/QCD_HT*_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v*/MINIAODSIM
    "MC_QCD1": 0.003089,
    "MC_QCD_singlephoton1": 0.003089,
    "MC_QCD2": 0.03316,
    "MC_QCD_singlephoton2": 0.03316,
    "MC_QCD3": 0.1566,
    "MC_QCD_singlephoton3": 0.1566,
    "MC_QCD4": 0.9067,
    "MC_QCD_singlephoton4": 0.9067,
    "MC_QCD5": 9.867,
    "MC_QCD_singlephoton5": 9.867,
    "MC_QCD6": 47.5,
    "MC_QCD_singlephoton6": 47.5
}

def getMCWeight(selectionType):
    # actual weight = crossSection * integrated_lumi_run2 / nGeneratedEvents
    # However, nGeneratedEvents is not directly available from XSDB
    # The "effective lumi" is available, though
    # So we use nGeneratedEvents = cross_section * effective_lumi
    # Therefore, weight = integrated_lumi_run2 / effective_lumi
    return ((35.9182 + 41.5273 + 59.7360)/(effectiveLuminosities[selectionType]))

stealthEnv.execute_in_env(commandToRun="eos {eP} mkdir -p {sER}/selections/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), functionToCallIfCommandExitsWithError=removeLock)
stealthEnv.execute_in_env(commandToRun="eos {eP} mkdir -p {sER}/statistics/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), functionToCallIfCommandExitsWithError=removeLock)

# Step 1 merge: sufficient for everything except datasets with different event weights
filesToCleanup = []
for selectionType in selectionTypesToRun:
    isMC = True
    isMCString = "true"
    if (("data" in selectionType)
        or (selectionType == "MC_EMEnrichedQCD")
        or (bool(re.match(r"^MC_GJet[0-9]*$", selectionType)))
        or (bool(re.match(r"^MC_GJet_singlephoton[0-9]*$", selectionType)))
        or (bool(re.match(r"^MC_QCD[0-9]*$", selectionType)))
        or (bool(re.match(r"^MC_QCD_singlephoton[0-9]*$", selectionType)))
        or (selectionType == "MC_hgg")
    ):
        isMC = False
        isMCString = "false"
    for year in yearsToRun:
        mergeStatistics = True
        if ((bool(re.match(r"^MC_QCD[0-9]*$", selectionType))) or (selectionType == "MC_EMEnrichedQCD") or (selectionType == "MC_hgg")):
            if (year != 2017): # The only reason we need these is to calculate ID efficiencies
                mergeStatistics = False
        if (bool(re.match(r"^MC_GJet[0-9]*$", selectionType))):
            if (year != 2016): # The only reason we need these is to calculate scaling systematics
                mergeStatistics = False
        if ((selectionType == "data_singlephoton") or
            (bool(re.match(r"^MC_GJet_singlephoton[0-9]*$", selectionType))) or
            (bool(re.match(r"^MC_QCD_singlephoton[0-9]*$", selectionType)))):
            mergeStatistics = False
        if mergeStatistics:
            inputFilesList_statistics = "fileLists/inputFileList_statistics_{t}{aJS}{eVS}_{y}{oI}.txt".format(oI=optional_identifier, t=selectionType, aJS=allJetsString, eVS=electronVetoString, y=year)
            mergeStatisticsCommand = "eventSelection/bin/mergeStatisticsHistograms inputFilesList={iFL} outputFolder={oF} outputFileName={oFN} isMC={iMCS}".format(iFL=inputFilesList_statistics, oF="{eP}/{sER}/statistics/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), oFN="merged_statistics_{t}{aJS}{eVS}_{y}.root".format(t=selectionType, aJS=allJetsString, eVS=electronVetoString, y=year), iMCS=isMCString)
            multiProcessLauncher.spawn(shellCommands=mergeStatisticsCommand, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="mergeLog_statistics_{t}{aJS}{eVS}_{y}.log".format(t=selectionType, aJS=allJetsString, eVS=electronVetoString, y=year), printDebug=True)
        for selectionRegion in ["signal", "signal_loose", "control_fakefake", "control_singlemedium", "control_singleloose", "control_singlefake"]:
            selectionRegionString = selectionRegion
            if (selectionRegion == "control_fakefake"): selectionRegionString = "control"
            mergeSelection = True
            if ((selectionRegion == "control_singlemedium") or (selectionRegion == "control_singleloose") or (selectionRegion == "control_singlefake")):
                mergeSelection = False
                if (((selectionType == "data_singlephoton") or
                     (bool(re.match(r"^MC_GJet_singlephoton[0-9]*$", selectionType))) or
                     (bool(re.match(r"^MC_QCD_singlephoton[0-9]*$", selectionType)))) and not(isMC)): mergeSelection = True
            else:
                if ((selectionType == "data_singlephoton") or (bool(re.match(r"^MC_GJet_singlephoton[0-9]*$", selectionType))) or (bool(re.match(r"^MC_QCD_singlephoton[0-9]*$", selectionType)))): mergeSelection = False
            if ((selectionType == "data_jetHT") or (selectionType == "MC_hgg") or (selectionType == "MC_EMEnrichedQCD")): mergeSelection = False
            if (((bool(re.match(r"^MC_GJet[0-9]*$", selectionType))) or (bool(re.match(r"^MC_GJet_singlephoton[0-9]*$", selectionType)))) and (year != 2016)): mergeSelection = False # The only reason we need these is to calculate scaling systematics
            if (((bool(re.match(r"^MC_QCD[0-9]*$", selectionType))) or (bool(re.match(r"^MC_QCD_singlephoton[0-9]*$", selectionType))) or (selectionType == "MC_EMEnrichedQCD")) and (year != 2017)): mergeSelection = False # The only reason we need these is to calculate ID efficiencies
            if not(mergeSelection): continue
            inputFilesList_selection = "fileLists/inputFileList_selections_{t}{aJS}{eVS}_{y}{oI}_{r}.txt".format(oI=optional_identifier, t=selectionType, aJS=allJetsString, eVS=electronVetoString, y=year, r=selectionRegion)
            outputFolder = "{eP}/{sER}/selections/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier)
            outputFilePath = "merged_selection_{t}{aJS}{eVS}_{y}_{sRS}.root".format(t=selectionType, aJS=allJetsString, eVS=electronVetoString, y=year, sRS=selectionRegionString)
            mergeSelectionCommand = "eventSelection/bin/mergeEventSelections inputFilesList={iFL} outputFolder={oF} outputFileName={oFP}".format(iFL=inputFilesList_selection, oF=outputFolder, oFP=outputFilePath)
            if ((bool(re.match(r"^MC_GJet[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_GJet_singlephoton[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_QCD[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_QCD_singlephoton[0-9]*$", selectionType)))):
                mergeSelectionCommand += " addWeightBranch={w:.9f}".format(w=getMCWeight(selectionType))
                mergeStep2FilePath = ""
                if (bool(re.match(r"^MC_GJet[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_GJet{aJS}{eVS}_2016{oI}_{r}.txt".format(oI=optional_identifier, aJS=allJetsString, eVS=electronVetoString, r=selectionRegion)
                elif (bool(re.match(r"^MC_GJet_singlephoton[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_GJet_singlephoton{aJS}{eVS}_2016{oI}_{r}.txt".format(oI=optional_identifier, aJS=allJetsString, eVS=electronVetoString, r=selectionRegion)
                elif (bool(re.match(r"^MC_QCD[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_QCD{aJS}{eVS}_2017{oI}_{r}.txt".format(oI=optional_identifier, aJS=allJetsString, eVS=electronVetoString, r=selectionRegion)
                elif (bool(re.match(r"^MC_QCD_singlephoton[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_QCD_singlephoton{aJS}{eVS}_2017{oI}_{r}.txt".format(oI=optional_identifier, aJS=allJetsString, eVS=electronVetoString, r=selectionRegion)
                os.system("echo {oF}/{oFP} >> {mS2FP}".format(oF=outputFolder, oFP=outputFilePath, mS2FP=mergeStep2FilePath))
                filesToCleanup.append("{sER}/selections/combined_DoublePhoton{oI}/{oFP}".format(sER=stealthEnv.stealthEOSRoot, oI=optional_identifier, oFP=outputFilePath))
            multiProcessLauncher.spawn(shellCommands=mergeSelectionCommand, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="mergeLog_selection_{t}{aJS}{eVS}_{y}_{sRS}.log".format(t=selectionType, aJS=allJetsString, eVS=electronVetoString, y=year, sRS=selectionRegionString), printDebug=True)
multiProcessLauncher.monitorToCompletion()

# Step 2 merge: for datasets with different event weights
monitoringNeeded = False
for selectionType in selectionTypesToRun_Step2:
    for year in yearsToRun:
        for selectionRegion in ["signal", "signal_loose", "control_fakefake", "control_singlemedium", "control_singleloose", "control_singlefake"]:
            selectionRegionString = "{sR}".format(sR=selectionRegion)
            if (selectionRegion == "control_fakefake"): selectionRegionString = "control"
            mergeSelection = True
            if ((selectionRegion == "control_singlemedium") or (selectionRegion == "control_singleloose") or (selectionRegion == "control_singlefake")):
                mergeSelection = False
                if ((selectionType == "MC_GJet_singlephoton") or
                    (selectionType == "MC_QCD_singlephoton")):
                    mergeSelection = True
            else:
                if ((selectionType == "MC_GJet_singlephoton") or
                    (selectionType == "MC_QCD_singlephoton")):
                    mergeSelection = False
            if (((selectionType == "MC_GJet") or (selectionType == "MC_GJet_singlephoton")) and (year != 2016)): mergeSelection = False
            if (((selectionType == "MC_QCD") or (selectionType == "MC_QCD_singlephoton")) and (year != 2017)): mergeSelection = False
            if not(mergeSelection): continue
            # mergeStep2FilePath = "fileLists/inputFileList_step2Merge_{t}{aJS}{eVS}_2016{oI}_{r}.txt".format(t=selectionType, oI=optional_identifier, aJS=allJetsString, eVS=electronVetoString, r=selectionRegion)
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_{t}{aJS}{eVS}_{y}{oI}_{r}.txt".format(t=selectionType, y=year, oI=optional_identifier, aJS=allJetsString, eVS=electronVetoString, r=selectionRegion)
            outputFolder = "{eP}/{sER}/selections/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier)
            outputFilePath = "merged_selection_{t}{aJS}{eVS}_{y}_{sRS}.root".format(t=selectionType, aJS=allJetsString, eVS=electronVetoString, y=year, sRS=selectionRegionString)
            mergeSelectionCommand = "eventSelection/bin/mergeEventSelections inputFilesList={mS2FP} outputFolder={oF} outputFileName={oFP}".format(mS2FP=mergeStep2FilePath, oF=outputFolder, oFP=outputFilePath)
            multiProcessLauncher.spawn(shellCommands=mergeSelectionCommand, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="mergeLog_selection_{t}{aJS}{eVS}_{y}_{sRS}.log".format(t=selectionType, aJS=allJetsString, eVS=electronVetoString, y=year, sRS=selectionRegionString), printDebug=True)
            monitoringNeeded = True
if monitoringNeeded: multiProcessLauncher.monitorToCompletion()

for fileToCleanup in filesToCleanup:
    stealthEnv.execute_in_env(commandToRun="eos {eP} rm {fTC}".format(eP=stealthEnv.EOSPrefix, fTC=fileToCleanup), functionToCallIfCommandExitsWithError=removeLock)

removeLock()
