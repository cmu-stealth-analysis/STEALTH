#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, signal, re
import tmMultiProcessLauncher # from tmPyUtils
import stealthEnv # from this folder

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Run event selection merging scripts.')
inputArgumentsParser.add_argument('--selectionsToRun', default="data", help="Comma-separated list of selections to run. Allowed: \"data\", \"data_singlephoton\", \"data_jetHT\", \"MC\", \"MC_EMEnrichedQCD\", \"MC_GJet16\", \"MC_GJet17\", \"MC_GJet18\", \"MC_GJet16_singlephoton\", \"MC_GJet17_singlephoton\", \"MC_GJet18_singlephoton\", \"MC_QCD16\", \"MC_QCD17\", \"MC_QCD18\", \"MC_QCD16_singlephoton\", \"MC_QCD17_singlephoton\", \"MC_QCD18_singlephoton\", \"MC_DiPhotonJets\", \"MC_EMEnrichedGJetPt16\", \"MC_EMEnrichedGJetPt17\", \"MC_EMEnrichedGJetPt18\", \"MC_HighHTQCD16\", \"MC_HighHTQCD17\", \"MC_HighHTQCD18\", or \"MC_hgg\". For MC selections, disable HLT photon trigger and enable additional MC selection. Default is \"data\".", type=str)
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
    "MC_HighHTQCD18": 8
}

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
    elif (inputSelectionToRun == "MC_EMEnrichedQCD"):
        selectionTypesToRun.append("MC_EMEnrichedQCD3")
        selectionTypesToRun.append("MC_EMEnrichedQCD2")
        selectionTypesToRun.append("MC_EMEnrichedQCD1")
        for selectionRegion in ["signal", "signal_loose", "control_fakefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_EMEnrichedQCD{oIS}_2017{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        os.system("rm -f fileLists/inputFileList_step2Merge_statistics_MC_EMEnrichedQCD{oIS}_2017{oI}.txt && touch fileLists/inputFileList_step2Merge_statistics_MC_EMEnrichedQCD{oIS}_2017{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString))
        selectionTypesToRun_Step2.append("MC_EMEnrichedQCD")
    elif (inputSelectionToRun == "MC_GJet16"):
        selectionTypesToRun.append("MC_GJet16_5")
        selectionTypesToRun.append("MC_GJet16_4")
        selectionTypesToRun.append("MC_GJet16_3")
        selectionTypesToRun.append("MC_GJet16_2")
        selectionTypesToRun.append("MC_GJet16_1")
        for selectionRegion in ["signal", "signal_loose", "control_fakefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_GJet16{oIS}_2016{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        os.system("rm -f fileLists/inputFileList_step2Merge_statistics_MC_GJet16{oIS}_2016{oI}.txt && touch fileLists/inputFileList_step2Merge_statistics_MC_GJet16{oIS}_2016{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString))
        selectionTypesToRun_Step2.append("MC_GJet16")
    elif (inputSelectionToRun == "MC_GJet16_singlephoton"):
        selectionTypesToRun.append("MC_GJet16_singlephoton5")
        selectionTypesToRun.append("MC_GJet16_singlephoton4")
        selectionTypesToRun.append("MC_GJet16_singlephoton3")
        selectionTypesToRun.append("MC_GJet16_singlephoton2")
        selectionTypesToRun.append("MC_GJet16_singlephoton1")
        for selectionRegion in ["control_singlemedium", "control_singleloose", "control_singlefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_GJet16_singlephoton{oIS}_2016{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        selectionTypesToRun_Step2.append("MC_GJet16_singlephoton")
    elif (inputSelectionToRun == "MC_GJet17"):
        selectionTypesToRun.append("MC_GJet17_4")
        selectionTypesToRun.append("MC_GJet17_3")
        selectionTypesToRun.append("MC_GJet17_2")
        selectionTypesToRun.append("MC_GJet17_1")
        for selectionRegion in ["signal", "signal_loose", "control_fakefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_GJet17{oIS}_2017{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        os.system("rm -f fileLists/inputFileList_step2Merge_statistics_MC_GJet17{oIS}_2017{oI}.txt && touch fileLists/inputFileList_step2Merge_statistics_MC_GJet17{oIS}_2017{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString))
        selectionTypesToRun_Step2.append("MC_GJet17")
    elif (inputSelectionToRun == "MC_GJet17_singlephoton"):
        selectionTypesToRun.append("MC_GJet17_singlephoton4")
        selectionTypesToRun.append("MC_GJet17_singlephoton3")
        selectionTypesToRun.append("MC_GJet17_singlephoton2")
        selectionTypesToRun.append("MC_GJet17_singlephoton1")
        for selectionRegion in ["control_singlemedium", "control_singleloose", "control_singlefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_GJet17_singlephoton{oIS}_2017{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        selectionTypesToRun_Step2.append("MC_GJet17_singlephoton")
    elif (inputSelectionToRun == "MC_GJet18"):
        selectionTypesToRun.append("MC_GJet18_4")
        selectionTypesToRun.append("MC_GJet18_3")
        selectionTypesToRun.append("MC_GJet18_2")
        selectionTypesToRun.append("MC_GJet18_1")
        for selectionRegion in ["signal", "signal_loose", "control_fakefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_GJet18{oIS}_2018{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        os.system("rm -f fileLists/inputFileList_step2Merge_statistics_MC_GJet18{oIS}_2018{oI}.txt && touch fileLists/inputFileList_step2Merge_statistics_MC_GJet18{oIS}_2018{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString))
        selectionTypesToRun_Step2.append("MC_GJet18")
    elif (inputSelectionToRun == "MC_GJet18_singlephoton"):
        selectionTypesToRun.append("MC_GJet18_singlephoton4")
        selectionTypesToRun.append("MC_GJet18_singlephoton3")
        selectionTypesToRun.append("MC_GJet18_singlephoton2")
        selectionTypesToRun.append("MC_GJet18_singlephoton1")
        for selectionRegion in ["control_singlemedium", "control_singleloose", "control_singlefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_GJet18_singlephoton{oIS}_2018{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        selectionTypesToRun_Step2.append("MC_GJet18_singlephoton")
    elif (inputSelectionToRun == "MC_QCD16"):
        selectionTypesToRun.append("MC_QCD16_6")
        selectionTypesToRun.append("MC_QCD16_5")
        selectionTypesToRun.append("MC_QCD16_4")
        selectionTypesToRun.append("MC_QCD16_3")
        selectionTypesToRun.append("MC_QCD16_2")
        selectionTypesToRun.append("MC_QCD16_1")
        for selectionRegion in ["signal", "signal_loose", "control_fakefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_QCD16{oIS}_2016{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        os.system("rm -f fileLists/inputFileList_step2Merge_statistics_MC_QCD16{oIS}_2016{oI}.txt && touch fileLists/inputFileList_step2Merge_statistics_MC_QCD16{oIS}_2016{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString))
        selectionTypesToRun_Step2.append("MC_QCD16")
    elif (inputSelectionToRun == "MC_QCD16_singlephoton"):
        selectionTypesToRun.append("MC_QCD16_singlephoton6")
        selectionTypesToRun.append("MC_QCD16_singlephoton5")
        selectionTypesToRun.append("MC_QCD16_singlephoton4")
        selectionTypesToRun.append("MC_QCD16_singlephoton3")
        selectionTypesToRun.append("MC_QCD16_singlephoton2")
        selectionTypesToRun.append("MC_QCD16_singlephoton1")
        for selectionRegion in ["control_singlemedium", "control_singleloose", "control_singlefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_QCD16_singlephoton{oIS}_2016{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        selectionTypesToRun_Step2.append("MC_QCD16_singlephoton")
    elif (inputSelectionToRun == "MC_QCD17"):
        selectionTypesToRun.append("MC_QCD17_6")
        selectionTypesToRun.append("MC_QCD17_5")
        selectionTypesToRun.append("MC_QCD17_4")
        selectionTypesToRun.append("MC_QCD17_3")
        selectionTypesToRun.append("MC_QCD17_2")
        selectionTypesToRun.append("MC_QCD17_1")
        for selectionRegion in ["signal", "signal_loose", "control_fakefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_QCD17{oIS}_2017{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        os.system("rm -f fileLists/inputFileList_step2Merge_statistics_MC_QCD17{oIS}_2017{oI}.txt && touch fileLists/inputFileList_step2Merge_statistics_MC_QCD17{oIS}_2017{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString))
        selectionTypesToRun_Step2.append("MC_QCD17")
    elif (inputSelectionToRun == "MC_QCD17_singlephoton"):
        selectionTypesToRun.append("MC_QCD17_singlephoton6")
        selectionTypesToRun.append("MC_QCD17_singlephoton5")
        selectionTypesToRun.append("MC_QCD17_singlephoton4")
        selectionTypesToRun.append("MC_QCD17_singlephoton3")
        selectionTypesToRun.append("MC_QCD17_singlephoton2")
        selectionTypesToRun.append("MC_QCD17_singlephoton1")
        for selectionRegion in ["control_singlemedium", "control_singleloose", "control_singlefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_QCD17_singlephoton{oIS}_2017{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        selectionTypesToRun_Step2.append("MC_QCD17_singlephoton")
    elif (inputSelectionToRun == "MC_QCD18"):
        selectionTypesToRun.append("MC_QCD18_6")
        selectionTypesToRun.append("MC_QCD18_5")
        selectionTypesToRun.append("MC_QCD18_4")
        selectionTypesToRun.append("MC_QCD18_3")
        selectionTypesToRun.append("MC_QCD18_2")
        selectionTypesToRun.append("MC_QCD18_1")
        for selectionRegion in ["signal", "signal_loose", "control_fakefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_QCD18{oIS}_2018{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        os.system("rm -f fileLists/inputFileList_step2Merge_statistics_MC_QCD18{oIS}_2018{oI}.txt && touch fileLists/inputFileList_step2Merge_statistics_MC_QCD18{oIS}_2018{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString))
        selectionTypesToRun_Step2.append("MC_QCD18")
    elif (inputSelectionToRun == "MC_QCD18_singlephoton"):
        selectionTypesToRun.append("MC_QCD18_singlephoton6")
        selectionTypesToRun.append("MC_QCD18_singlephoton5")
        selectionTypesToRun.append("MC_QCD18_singlephoton4")
        selectionTypesToRun.append("MC_QCD18_singlephoton3")
        selectionTypesToRun.append("MC_QCD18_singlephoton2")
        selectionTypesToRun.append("MC_QCD18_singlephoton1")
        for selectionRegion in ["control_singlemedium", "control_singleloose", "control_singlefake"]:
            mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_QCD18_singlephoton{oIS}_2018{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
            os.system("rm -f {mS2FP} && touch {mS2FP}".format(mS2FP=mergeStep2FilePath))
        selectionTypesToRun_Step2.append("MC_QCD18_singlephoton")
    elif (inputSelectionToRun == "MC_hgg"):
        selectionTypesToRun.append("MC_hgg")
    elif (inputSelectionToRun == "MC_DiPhotonJets"):
        selectionTypesToRun.append("MC_DiPhotonJets")
    else:
        MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD)([0-9]*)$", inputSelectionToRun)
        if MCBKGMatch:
            full_match = MCBKGMatch.group(0)
            dataset_id = MCBKGMatch.group(1)
            year_last_two_digits_str = MCBKGMatch.group(2)
            for index_subsample in range(n_subsamples[full_match], 0, -1):
                selectionTypesToRun.append("{m}_{i}".format(m=full_match, i=index_subsample))
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

# EDIT: The following logic was incorrect, because the effective lumi on XSDB is only for 1 million events (!)
# what a waste of time... commenting out

# effectiveLuminosities = {
#     # DAS query for MC_GJet16: dataset dataset=/GJets_DR-0p4_HT-*_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_qcut19_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
#     "MC_GJet16_1": 0.05745,
#     "MC_GJet16_singlephoton1": 0.05745,
#     "MC_GJet16_2": 0.1865,
#     "MC_GJet16_singlephoton2": 0.1865,
#     "MC_GJet16_3": 0.849,
#     "MC_GJet16_singlephoton3": 0.849,
#     "MC_GJet16_4": 7.588,
#     "MC_GJet16_singlephoton4": 7.588,
#     "MC_GJet16_5": 22.59,
#     "MC_GJet16_singlephoton5": 22.59,

#     # DAS query for MC_GJet17: dataset dataset=/GJets_DR-0p4_HT*/RunIIFall17*/MINIAODSIM
#     "MC_GJet17_1": 0.1986,
#     "MC_GJet17_singlephoton1": 0.1986,
#     "MC_GJet17_2": 0.8849,
#     "MC_GJet17_singlephoton2": 0.8849,
#     "MC_GJet17_3": 7.986,
#     "MC_GJet17_singlephoton3": 7.986,
#     "MC_GJet17_4": 24.35,
#     "MC_GJet17_singlephoton4": 24.35,

#     # DAS query for MC_GJet18: dataset dataset=/GJets_DR-0p4_HT-*_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM
#     "MC_GJet18_1": 0.1986,
#     "MC_GJet18_singlephoton1": 0.1986,
#     "MC_GJet18_2": 0.8849,
#     "MC_GJet18_singlephoton2": 0.8849,
#     "MC_GJet18_3": 7.986,
#     "MC_GJet18_singlephoton3": 7.986,
#     "MC_GJet18_4": 24.35,
#     "MC_GJet18_singlephoton4": 24.35,

#     # DAS query for MC_QCD16: dataset dataset=/QCD_HT*_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM
#     "MC_QCD16_1": 0.002878,
#     "MC_QCD16_singlephoton1": 0.002878,
#     "MC_QCD16_2": 0.03119,
#     "MC_QCD16_singlephoton2": 0.03119,
#     "MC_QCD16_3": 0.1464,
#     "MC_QCD16_singlephoton3": 0.1464,
#     "MC_QCD16_4": 0.8285,
#     "MC_QCD16_singlephoton4": 0.8285,
#     "MC_QCD16_5": 8.335,
#     "MC_QCD16_singlephoton5": 8.335,
#     "MC_QCD16_6": 39.6,
#     "MC_QCD16_singlephoton6": 39.6,

#     # DAS query for MC_QCD17: dataset dataset=/QCD_HT*_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v*/MINIAODSIM
#     "MC_QCD17_1": 0.003089,
#     "MC_QCD17_singlephoton1": 0.003089,
#     "MC_QCD17_2": 0.03316,
#     "MC_QCD17_singlephoton2": 0.03316,
#     "MC_QCD17_3": 0.1566,
#     "MC_QCD17_singlephoton3": 0.1566,
#     "MC_QCD17_4": 0.9067,
#     "MC_QCD17_singlephoton4": 0.9067,
#     "MC_QCD17_5": 9.867,
#     "MC_QCD17_singlephoton5": 9.867,
#     "MC_QCD17_6": 47.5,
#     "MC_QCD17_singlephoton6": 47.5,

#     # DAS query for MC_QCD18: dataset dataset=/QCD_HT*_TuneCP5_13TeV-madgraph-pythia8/RunIISpring18MiniAOD-100X_upgrade2018_realistic_v10-v1/MINIAODSIM
#     # slightly different from MC_QCD17 for some reason
#     "MC_QCD18_1": 0.00308,
#     "MC_QCD18_singlephoton1": 0.00308,
#     "MC_QCD18_2": 0.033,
#     "MC_QCD18_singlephoton2": 0.033,
#     "MC_QCD18_3": 0.1572,
#     "MC_QCD18_singlephoton3": 0.1572,
#     "MC_QCD18_4": 0.9015,
#     "MC_QCD18_singlephoton4": 0.9015,
#     "MC_QCD18_5": 9.848,
#     "MC_QCD18_singlephoton5": 9.848,
#     "MC_QCD18_6": 47.61,
#     "MC_QCD18_singlephoton6": 47.61,

#     # DAS query for MC_EMEnrichedQCD: dataset dataset=/QCD_Pt-*_DoubleEMEnriched_MGG*/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_*/MINIAODSIM
#     "MC_EMEnrichedQCD1": 0.0404,
#     "MC_EMEnrichedQCD2": 0.004121,
#     "MC_EMEnrichedQCD3": 0.008517
# }

# # integratedLuminosities = {
# #     2016: 35.9182,
# #     2017: 41.5273,
# #     2018: 59.7360
# # }

# lumi_values_raw_json = None
# integratedLuminosities = {}
# with open("xSecLumiInfo/lumi_run2.json", 'r') as lumi_json_file_handle:
#     lumi_values_raw_json = json.load(lumi_json_file_handle)
#     for year in [2016, 2017, 2018]:
#         integratedLuminosities[year] = lumi_values_raw_json[str(year)]

# def getMCWeight(selectionType, year):
#     # actual weight = crossSection * integrated_lumi / nGeneratedEvents
#     # However, nGeneratedEvents is not directly available from XSDB
#     # The "effective lumi" is available, though
#     # So we use nGeneratedEvents = cross_section * effective_lumi
#     # Therefore, weight = integrated_lumi / effective_lumi
#     return ((integratedLuminosities[year])/(effectiveLuminosities[selectionType]))

stealthEnv.execute_in_env(commandToRun="eos {eP} mkdir -p {sER}/selections/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), functionToCallIfCommandExitsWithError=removeLock)
stealthEnv.execute_in_env(commandToRun="eos {eP} mkdir -p {sER}/statistics/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), functionToCallIfCommandExitsWithError=removeLock)

# Step 1 merge: sufficient for everything except datasets with different event weights
filesToCleanup = []
for selectionType in selectionTypesToRun:
    isMC = ((selectionType == "MC_stealth_t5") or (selectionType == "MC_stealth_t6"))
    isMCString = None
    if isMC: isMCString = "true"
    else: isMCString = "false"
    # if (("data" in selectionType)
    #     or (bool(re.match(r"^MC_EMEnrichedQCD[0-9]*$", selectionType)))
    #     or (selectionType == "MC_hgg")
    #     or (bool(re.match(r"^MC_GJet16_singlephoton[0-9]*$", selectionType)))
    #     or (bool(re.match(r"^MC_GJet17_singlephoton[0-9]*$", selectionType)))
    #     or (bool(re.match(r"^MC_GJet18_singlephoton[0-9]*$", selectionType)))
    #     or (bool(re.match(r"^MC_QCD16_singlephoton[0-9]*$", selectionType)))
    #     or (bool(re.match(r"^MC_QCD17_singlephoton[0-9]*$", selectionType)))
    #     or (bool(re.match(r"^MC_QCD18_singlephoton[0-9]*$", selectionType)))):
    #     isMC = False
    #     isMCString = "false"
    for year in yearsToRun:
        mergeStatistics = True
        if (bool(re.match(r"^MC_EMEnrichedQCD[0-9]*$", selectionType))):
            if (year != 2017): # The only reason we need these is to calculate ID efficiencies
                mergeStatistics = False
        if (bool(re.match(r"^MC_GJet16_[0-9]*$", selectionType))):
            if (year != 2016):
                mergeStatistics = False
        if (bool(re.match(r"^MC_GJet17_[0-9]*$", selectionType))):
            if (year != 2017):
                mergeStatistics = False
        if (bool(re.match(r"^MC_GJet18_[0-9]*$", selectionType))):
            if (year != 2018):
                mergeStatistics = False
        if (bool(re.match(r"^MC_QCD16_[0-9]*$", selectionType))):
            if (year != 2016):
                mergeStatistics = False
        if (bool(re.match(r"^MC_QCD17_[0-9]*$", selectionType))):
            if (year != 2017):
                mergeStatistics = False
        if (bool(re.match(r"^MC_QCD18_[0-9]*$", selectionType))):
            if (year != 2018):
                mergeStatistics = False
        MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD)([0-9]*)_([0-9]*)$", selectionType)
        if MCBKGMatch:
            year_last_two_digits_str = MCBKGMatch.group(2)
            year_from_match = 2000 + int(0.5 + float(year_last_two_digits_str))
            if not(year_from_match == year): mergeStatistics = False
        if ((selectionType == "data_singlephoton") or
            (bool(re.match(r"^MC_GJet16_singlephoton[0-9]*$", selectionType))) or
            (bool(re.match(r"^MC_GJet17_singlephoton[0-9]*$", selectionType))) or
            (bool(re.match(r"^MC_GJet18_singlephoton[0-9]*$", selectionType))) or
            (bool(re.match(r"^MC_QCD16_singlephoton[0-9]*$", selectionType))) or
            (bool(re.match(r"^MC_QCD17_singlephoton[0-9]*$", selectionType))) or
            (bool(re.match(r"^MC_QCD18_singlephoton[0-9]*$", selectionType)))):
            mergeStatistics = False
        if mergeStatistics:
            inputFilesList_statistics = "fileLists/inputFileList_statistics_{t}{oIS}_{y}{oI}.txt".format(oI=optional_identifier, t=selectionType, oIS=overallIdentificationString, y=year)
            outputFilePath = "merged_statistics_{t}{oIS}_{y}.root".format(t=selectionType, oIS=overallIdentificationString, y=year)
            outputFolder = "{eP}/{sER}/statistics/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier)
            mergeStatisticsCommand = "eventSelection/bin/mergeStatisticsHistograms inputFilesList={iFL} outputFolder={oF} outputFileName={oFP} isMC={iMCS}".format(iFL=inputFilesList_statistics, oF=outputFolder, oFP=outputFilePath, iMCS=isMCString)
            if (bool(re.match(r"^MC_EMEnrichedQCD[0-9]*$", selectionType)) or
                (bool(re.match(r"^MC_GJet16_[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_GJet17_[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_GJet18_[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_GJet16_singlephoton[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_GJet17_singlephoton[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_GJet18_singlephoton[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_QCD16_[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_QCD17_[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_QCD18_[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_QCD16_singlephoton[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_QCD17_singlephoton[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_QCD18_singlephoton[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD)([0-9]*)_([0-9]*)$", selectionType)))
            ):
                mergeStep2FilePath = ""
                if (bool(re.match(r"^MC_EMEnrichedQCD[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_statistics_MC_EMEnrichedQCD{oIS}_2017{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString)
                if (bool(re.match(r"^MC_GJet16_[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_statistics_MC_GJet16{oIS}_2016{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString)
                elif (bool(re.match(r"^MC_GJet16_singlephoton[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_statistics_MC_GJet16_singlephoton{oIS}_2016{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString)
                elif (bool(re.match(r"^MC_GJet17_[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_statistics_MC_GJet17{oIS}_2017{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString)
                elif (bool(re.match(r"^MC_GJet17_singlephoton[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_statistics_MC_GJet17_singlephoton{oIS}_2017{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString)
                elif (bool(re.match(r"^MC_GJet18_[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_statistics_MC_GJet18{oIS}_2018{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString)
                elif (bool(re.match(r"^MC_GJet18_singlephoton[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_statistics_MC_GJet18_singlephoton{oIS}_2018{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString)
                elif (bool(re.match(r"^MC_QCD16_[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_statistics_MC_QCD16{oIS}_2016{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString)
                elif (bool(re.match(r"^MC_QCD16_singlephoton[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_statistics_MC_QCD16_singlephoton{oIS}_2016{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString)
                elif (bool(re.match(r"^MC_QCD17_[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_statistics_MC_QCD17{oIS}_2017{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString)
                elif (bool(re.match(r"^MC_QCD17_singlephoton[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_statistics_MC_QCD17_singlephoton{oIS}_2017{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString)
                elif (bool(re.match(r"^MC_QCD18_[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_statistics_MC_QCD18{oIS}_2018{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString)
                elif (bool(re.match(r"^MC_QCD18_singlephoton[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_statistics_MC_QCD18_singlephoton{oIS}_2018{oI}.txt".format(oI=optional_identifier, oIS=overallIdentificationString)
                else:
                    MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD)([0-9]*)_([0-9]*)$", selectionType)
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
                     (bool(re.match(r"^MC_GJet16_singlephoton[0-9]*$", selectionType))) or
                     (bool(re.match(r"^MC_GJet17_singlephoton[0-9]*$", selectionType))) or
                     (bool(re.match(r"^MC_GJet18_singlephoton[0-9]*$", selectionType))) or
                     (bool(re.match(r"^MC_QCD16_singlephoton[0-9]*$", selectionType))) or
                     (bool(re.match(r"^MC_QCD17_singlephoton[0-9]*$", selectionType))) or
                     (bool(re.match(r"^MC_QCD18_singlephoton[0-9]*$", selectionType)))): mergeSelection = True
            else:
                if ((selectionType == "data_singlephoton") or (bool(re.match(r"^MC_GJet16_singlephoton[0-9]*$", selectionType))) or (bool(re.match(r"^MC_GJet17_singlephoton[0-9]*$", selectionType))) or (bool(re.match(r"^MC_GJet18_singlephoton[0-9]*$", selectionType))) or (bool(re.match(r"^MC_QCD16_singlephoton[0-9]*$", selectionType))) or (bool(re.match(r"^MC_QCD17_singlephoton[0-9]*$", selectionType))) or (bool(re.match(r"^MC_QCD18_singlephoton[0-9]*$", selectionType)))): mergeSelection = False
            if ((selectionType == "data_jetHT") or (selectionType == "MC_hgg")): mergeSelection = False
            if (((bool(re.match(r"^MC_GJet16_[0-9]*$", selectionType))) or (bool(re.match(r"^MC_GJet16_singlephoton[0-9]*$", selectionType)))) and (year != 2016)): mergeSelection = False
            if (((bool(re.match(r"^MC_GJet17_[0-9]*$", selectionType))) or (bool(re.match(r"^MC_GJet17_singlephoton[0-9]*$", selectionType)))) and (year != 2017)): mergeSelection = False
            if (((bool(re.match(r"^MC_GJet18_[0-9]*$", selectionType))) or (bool(re.match(r"^MC_GJet18_singlephoton[0-9]*$", selectionType)))) and (year != 2018)): mergeSelection = False
            if (((bool(re.match(r"^MC_QCD16_[0-9]*$", selectionType))) or (bool(re.match(r"^MC_QCD16_singlephoton[0-9]*$", selectionType)))) and (year != 2016)): mergeSelection = False
            if (((bool(re.match(r"^MC_QCD17_[0-9]*$", selectionType))) or (bool(re.match(r"^MC_QCD17_singlephoton[0-9]*$", selectionType)))) and (year != 2017)): mergeSelection = False
            if (((bool(re.match(r"^MC_QCD18_[0-9]*$", selectionType))) or (bool(re.match(r"^MC_QCD18_singlephoton[0-9]*$", selectionType)))) and (year != 2018)): mergeSelection = False
            if ((bool(re.match(r"^MC_EMEnrichedQCD[0-9]*$", selectionType))) and (year != 2017)): mergeSelection = False
            MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD)([0-9]*)_([0-9]*)$", selectionType)
            if MCBKGMatch:
                year_last_two_digits_str = MCBKGMatch.group(2)
                year_from_match = 2000 + int(0.5 + float(year_last_two_digits_str))
                if not(year_from_match == year): mergeSelection = False
            if not(mergeSelection): continue
            inputFilesList_selection = "fileLists/inputFileList_selections_{t}{oIS}_{y}{oI}_{r}.txt".format(oI=optional_identifier, t=selectionType, oIS=overallIdentificationString, y=year, r=selectionRegion)
            outputFolder = "{eP}/{sER}/selections/combined_DoublePhoton{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier)
            outputFilePath = "merged_selection_{t}{oIS}_{y}_{sRS}.root".format(t=selectionType, oIS=overallIdentificationString, y=year, sRS=selectionRegionString)
            mergeSelectionCommand = "eventSelection/bin/mergeEventSelections inputFilesList={iFL} outputFolder={oF} outputFileName={oFP}".format(iFL=inputFilesList_selection, oF=outputFolder, oFP=outputFilePath)
            if ((bool(re.match(r"^MC_EMEnrichedQCD[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_GJet16_[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_GJet16_singlephoton[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_GJet17_[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_GJet17_singlephoton[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_GJet18_[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_GJet18_singlephoton[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_QCD16_[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_QCD16_singlephoton[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_QCD17_[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_QCD17_singlephoton[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_QCD18_[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_QCD18_singlephoton[0-9]*$", selectionType))) or
                (bool(re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD)([0-9]*)_([0-9]*)$", selectionType)))
            ):
                # mergeSelectionCommand += " addWeightBranch={w:.9f}".format(w=getMCWeight(selectionType, year))
                mergeStep2FilePath = ""
                if (bool(re.match(r"^MC_EMEnrichedQCD[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_EMEnrichedQCD{oIS}_2017{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
                if (bool(re.match(r"^MC_GJet16_[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_GJet16{oIS}_2016{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
                elif (bool(re.match(r"^MC_GJet16_singlephoton[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_GJet16_singlephoton{oIS}_2016{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
                elif (bool(re.match(r"^MC_GJet17_[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_GJet17{oIS}_2017{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
                elif (bool(re.match(r"^MC_GJet17_singlephoton[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_GJet17_singlephoton{oIS}_2017{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
                elif (bool(re.match(r"^MC_GJet18_[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_GJet18{oIS}_2018{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
                elif (bool(re.match(r"^MC_GJet18_singlephoton[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_GJet18_singlephoton{oIS}_2018{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
                elif (bool(re.match(r"^MC_QCD16_[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_QCD16{oIS}_2016{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
                elif (bool(re.match(r"^MC_QCD16_singlephoton[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_QCD16_singlephoton{oIS}_2016{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
                elif (bool(re.match(r"^MC_QCD17_[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_QCD17{oIS}_2017{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
                elif (bool(re.match(r"^MC_QCD17_singlephoton[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_QCD17_singlephoton{oIS}_2017{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
                elif (bool(re.match(r"^MC_QCD18_[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_QCD18{oIS}_2018{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
                elif (bool(re.match(r"^MC_QCD18_singlephoton[0-9]*$", selectionType))):
                    mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_QCD18_singlephoton{oIS}_2018{oI}_{r}.txt".format(oI=optional_identifier, oIS=overallIdentificationString, r=selectionRegion)
                else:
                    MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD)([0-9]*)_([0-9]*)$", selectionType)
                    if MCBKGMatch:
                        full_match = MCBKGMatch.group(0)
                        dataset_id = MCBKGMatch.group(1)
                        year_last_two_digits_str = MCBKGMatch.group(2)
                        mergeStep2FilePath = "fileLists/inputFileList_step2Merge_MC_{did}{y2}{oIS}_20{y2}{oI}_{r}.txt".format(oI=optional_identifier, did=dataset_id, oIS=overallIdentificationString, r=selectionRegion, y2=year_last_two_digits_str)
                os.system("echo {oF}/{oFP} >> {mS2FP}".format(oF=outputFolder, oFP=outputFilePath, mS2FP=mergeStep2FilePath))
                filesToCleanup.append("{sER}/selections/combined_DoublePhoton{oI}/{oFP}".format(sER=stealthEnv.stealthEOSRoot, oI=optional_identifier, oFP=outputFilePath))
            multiProcessLauncher.spawn(shellCommands=mergeSelectionCommand, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="mergeLog_selection_{t}{oIS}_{y}_{sRS}.log".format(t=selectionType, oIS=overallIdentificationString, y=year, sRS=selectionRegionString), printDebug=True)
multiProcessLauncher.monitorToCompletion()

# Step 2 merge: for datasets with different event weights
monitoringNeeded = False
for selectionType in selectionTypesToRun_Step2:
    isMC = ((selectionType == "MC_stealth_t5") or (selectionType == "MC_stealth_t6"))
    isMCString = None
    if isMC: isMCString = "true"
    else: isMCString = "false"
    for year in yearsToRun:
        mergeStatistics = True
        if (selectionType == "MC_EMEnrichedQCD"):
            if (year != 2017):
                mergeStatistics = False
        if (selectionType == "MC_GJet16"):
            if (year != 2016):
                mergeStatistics = False
        if (selectionType == "MC_GJet17"):
            if (year != 2017):
                mergeStatistics = False
        if (selectionType == "MC_GJet18"):
            if (year != 2018):
                mergeStatistics = False
        if (selectionType == "MC_QCD16"):
            if (year != 2016):
                mergeStatistics = False
        if (selectionType == "MC_QCD17"):
            if (year != 2017):
                mergeStatistics = False
        if (selectionType == "MC_QCD18"):
            if (year != 2018):
                mergeStatistics = False
        MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD)([0-9]*)$", selectionType)
        if MCBKGMatch:
            year_last_two_digits_str = MCBKGMatch.group(2)
            year_from_match = 2000 + int(0.5 + float(year_last_two_digits_str))
            if not(year_from_match == year): mergeStatistics = False
        if ((selectionType == "data_singlephoton") or
            (selectionType == "MC_GJet16_singlephoton") or
            (selectionType == "MC_GJet17_singlephoton") or
            (selectionType == "MC_GJet18_singlephoton") or
            (selectionType == "MC_QCD16_singlephoton") or
            (selectionType == "MC_QCD17_singlephoton") or
            (selectionType == "MC_QCD18_singlephoton")):
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
                if ((selectionType == "MC_GJet16_singlephoton") or
                    (selectionType == "MC_GJet17_singlephoton") or
                    (selectionType == "MC_GJet18_singlephoton") or
                    (selectionType == "MC_QCD16_singlephoton") or
                    (selectionType == "MC_QCD17_singlephoton") or
                    (selectionType == "MC_QCD18_singlephoton")):
                    mergeSelection = True
            else:
                if ((selectionType == "MC_GJet16_singlephoton") or
                    (selectionType == "MC_GJet17_singlephoton") or
                    (selectionType == "MC_GJet18_singlephoton") or
                    (selectionType == "MC_QCD16_singlephoton") or
                    (selectionType == "MC_QCD17_singlephoton") or
                    (selectionType == "MC_QCD18_singlephoton")):
                    mergeSelection = False
            if ((selectionType == "MC_EMEnrichedQCD") and (year != 2017)): mergeSelection = False
            if (((selectionType == "MC_GJet16") or (selectionType == "MC_GJet16_singlephoton")) and (year != 2016)): mergeSelection = False
            if (((selectionType == "MC_GJet17") or (selectionType == "MC_GJet17_singlephoton")) and (year != 2017)): mergeSelection = False
            if (((selectionType == "MC_GJet18") or (selectionType == "MC_GJet18_singlephoton")) and (year != 2018)): mergeSelection = False
            if (((selectionType == "MC_QCD16") or (selectionType == "MC_QCD16_singlephoton")) and (year != 2016)): mergeSelection = False
            if (((selectionType == "MC_QCD17") or (selectionType == "MC_QCD17_singlephoton")) and (year != 2017)): mergeSelection = False
            if (((selectionType == "MC_QCD18") or (selectionType == "MC_QCD18_singlephoton")) and (year != 2018)): mergeSelection = False
            MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD)([0-9]*)$", selectionType)
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
