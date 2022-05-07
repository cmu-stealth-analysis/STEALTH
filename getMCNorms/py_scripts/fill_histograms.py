#!/usr/bin/env python

from __future__ import print_function, division

import os

os.system("rm -vf *.pyc")

import sys, signal, argparse, subprocess
import tmMultiProcessLauncher # from tmPyUtils
import stealthEnv # from this folder

inputArgumentsParser = argparse.ArgumentParser(description='Run analysis chain.')
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArgumentsParser.add_argument('--histCategory', required=True, choices=["diphoton", "singlephoton", "pureQCD"], help="Category of histograms to create.",type=str)
inputArgumentsParser.add_argument('--isDryRun', action='store_true', help="Only print the commands to run, do not actually run them.")
inputArguments = inputArgumentsParser.parse_args()

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)
analysisOutputDirectory = "{aR}/analysis{oI}".format(aR=stealthEnv.analysisRoot, oI=optional_identifier)
outputDirectoryEOS = "{sER}/analysisEOSAreas/analysis{oI}/MCNorms".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier)
analysisLogsDirectory = "{aOD}/analysisLogs".format(aOD=analysisOutputDirectory)

subprocess.check_call("mkdir -p {d}".format(d=analysisLogsDirectory), shell=True, executable="/bin/bash")
subprocess.check_call("eos {eP} mkdir -p {o}".format(eP=stealthEnv.EOSPrefix, o=outputDirectoryEOS), shell=True, executable="/bin/bash")
subprocess.check_call("cd getMCNorms && make && cd ..", shell=True, executable="/bin/bash")

def signal_handler(sig, frame):
    sys.exit("Terminated by user.")
signal.signal(signal.SIGINT, signal_handler)

multiProcessLauncher = tmMultiProcessLauncher.tmMultiProcessLauncher(logOutputFolder=analysisLogsDirectory, printDebug=True)

processes_BKG = [# "DiPhotonJets", 
                 "DiPhotonJetsBox", "GJetHT", "HighHTQCD"]
is_year_dependent = {
    # "DiPhotonJets": False,
    "DiPhotonJetsBox": False,
    "GJetHT": True,
    "HighHTQCD": True
}
n_subsamples = {
    "MC_HighHTQCD16": 7,
    "MC_HighHTQCD17": 8,
    "MC_HighHTQCD18": 8,
    "MC_GJetHT16": 5,
    "MC_GJetHT17": 4,
    "MC_GJetHT18": 4
}

sources = {}

# diphoton, to estimate DiPhotonJetsBox contributions
sources["diphoton"] = {}
sources["diphoton"]["data"] = []
for year in [2016, 2017, 2018]:
    (sources["diphoton"]["data"]).append("fileLists/inputFileList_selections_data_noJetSelection_{y4}{oI}_signal.txt".format(y4=year, oI=optional_identifier))
for process_BKG in processes_BKG:
    sources["diphoton"][process_BKG] = []
    for year in [2016, 2017, 2018]:
        year_last_two_digits_str = ""
        if is_year_dependent[process_BKG]:
            year_last_two_digits_str = str(year-2000)
            for index_subsample in range(1, 1+n_subsamples["MC_" + process_BKG + year_last_two_digits_str]):
                (sources["diphoton"][process_BKG]).append("fileLists/inputFileList_selections_MC_{p}{y2}_{i}_noJetSelection_{y4}{oI}_signal.txt".format(p=process_BKG, y2=year_last_two_digits_str, i=index_subsample, y4=year, oI=optional_identifier))
        else:
            (sources["diphoton"][process_BKG]).append("fileLists/inputFileList_selections_MC_{p}_noJetSelection_{y4}{oI}_signal.txt".format(p=process_BKG, y4=year, oI=optional_identifier))

# singlephoton, to estimate GJet contribution
sources["singlephoton"] = {}
sources["singlephoton"]["data"] = []
for year in [2016, 2017, 2018]:
    # (sources["singlephoton"]["data"]).append("{i}/merged_selection_data_singlephoton_{y4}_control_singlemedium.root".format(i=source_directory_data, y4=year))
    (sources["singlephoton"]["data"]).append("fileLists/inputFileList_selections_data_singlephoton_{y4}{oI}_control_singlemedium.txt".format(y4=year, oI=optional_identifier))
for process_BKG in processes_BKG:
    sources["singlephoton"][process_BKG] = []
    for year in [2016, 2017, 2018]:
        year_last_two_digits_str = ""
        if is_year_dependent[process_BKG]:
            year_last_two_digits_str = str(year-2000)
            for index_subsample in range(1, 1+n_subsamples["MC_" + process_BKG + year_last_two_digits_str]):
                (sources["singlephoton"][process_BKG]).append("fileLists/inputFileList_selections_MC_{p}{y2}_singlephoton_{i}_{y4}{oI}_control_singlemedium.txt".format(p=process_BKG, y2=year_last_two_digits_str, i=index_subsample, y4=year, oI=optional_identifier))
        else:
            (sources["singlephoton"][process_BKG]).append("fileLists/inputFileList_selections_MC_{p}_singlephoton_{y4}{oI}_control_singlemedium.txt".format(p=process_BKG, y4=year, oI=optional_identifier))

# QCD, to estimate pure QCD contribution
sources["pureQCD"] = {}
sources["pureQCD"]["data"] = []
for year in [2016, 2017, 2018]:
    # (sources["pureQCD"]["data"]).append("{i}/merged_selection_data_pureQCD_{y4}_control_singlemedium.root".format(i=source_directory_data, y4=year))
    (sources["pureQCD"]["data"]).append("fileLists/inputFileList_selections_data_noPhotonSelection_{y4}{oI}_unified.txt".format(y4=year, oI=optional_identifier))
for process_BKG in processes_BKG:
    sources["pureQCD"][process_BKG] = []
    for year in [2016, 2017, 2018]:
        year_last_two_digits_str = ""
        if is_year_dependent[process_BKG]:
            year_last_two_digits_str = str(year-2000)
            for index_subsample in range(1, 1+n_subsamples["MC_" + process_BKG + year_last_two_digits_str]):
                (sources["pureQCD"][process_BKG]).append("fileLists/inputFileList_selections_MC_{p}{y2}_{i}_noPhotonSelection_{y4}{oI}_unified.txt".format(p=process_BKG, y2=year_last_two_digits_str, i=index_subsample, y4=year, oI=optional_identifier))
        else:
            (sources["pureQCD"][process_BKG]).append("fileLists/inputFileList_selections_MC_{p}_noPhotonSelection_{y4}{oI}_unified.txt".format(p=process_BKG, y4=year, oI=optional_identifier))

for process in (["data"] + processes_BKG):
    command_to_run = "./getMCNorms/bin/{c}".format(c=inputArguments.histCategory)
    command_to_run += " inputPathsFiles={i}".format(i=(",".join(sources[inputArguments.histCategory][process])))
    command_to_run += " outputFolder={eP}/{o}".format(eP=stealthEnv.EOSPrefix, o=outputDirectoryEOS)
    command_to_run += " outputFileName=histograms_{p}_{c}.root".format(p=process, c=inputArguments.histCategory)
    if (process == "data"):
        command_to_run += " addMCWeights=false"
    else:
        command_to_run += " addMCWeights=true"
    if (inputArguments.isDryRun):
        print("Not spawning due to dry run flag: {c}".format(c=command_to_run))
    else:
        multiProcessLauncher.spawn(shellCommands=command_to_run, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step_MCNorms_fill_histograms_{p}_{c}.log".format(p=process, c=inputArguments.histCategory), printDebug=True)
if not(inputArguments.isDryRun): multiProcessLauncher.monitorToCompletion()

# ./getMCNorms/py_scripts/fill_histograms.py --histCategory "diphoton" && ./getMCNorms/py_scripts/fill_histograms.py --histCategory "singlephoton" && ./getMCNorms/py_scripts/fill_histograms.py --histCategory "pureQCD"
