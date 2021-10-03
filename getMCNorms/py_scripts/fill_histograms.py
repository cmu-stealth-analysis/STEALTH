#!/usr/bin/env python

from __future__ import print_function, division

import os

os.system("rm -vf *.pyc")

import sys, signal, argparse, subprocess
import tmMultiProcessLauncher # from tmPyUtils
import stealthEnv # from this folder

inputArgumentsParser = argparse.ArgumentParser(description='Run analysis chain.')
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArgumentsParser.add_argument('--selectionSuffix', default="", help='If set, the input n-tuples are read with this suffix.',type=str)
inputArgumentsParser.add_argument('--histCategory', required=True, choices=["singlephoton"], help="Category of histograms to create.",type=str)
inputArgumentsParser.add_argument('--isDryRun', action='store_true', help="Only print the commands to run, do not actually run them.")
inputArguments = inputArgumentsParser.parse_args()

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)
selection_suffix = ""
if (inputArguments.selectionSuffix != ""): selection_suffix = "_{sS}".format(sS=inputArguments.selectionSuffix)
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

processes_BKG = ["DiPhotonJets", "GJetHT", "HighHTQCD"]
is_year_dependent = {
    "DiPhotonJets": False,
    "GJetHT": True,
    "HighHTQCD": True
}

source_directory_data = "{eP}/{sER}/selections/combined_DoublePhoton{sS}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix)
sources = {}
for hist_category in ["singlephoton"]:
    sources[hist_category] = {}
    sources[hist_category]["data"] = []
    for year in [2016, 2017, 2018]:
        (sources[hist_category]["data"]).append("{i}/merged_selection_data_singlephoton_{y4}_control_singlemedium.root".format(i=source_directory_data, y4=year))
    for process_BKG in processes_BKG:
        sources[hist_category][process_BKG] = []
        for year in [2016, 2017, 2018]:
            year_last_two_digits_str = ""
            if is_year_dependent[process_BKG]: year_last_two_digits_str = str(year-2000)
            (sources[hist_category][process_BKG]).append("{i}/merged_selection_MC_{p}{y2}_singlephoton_{y4}_control_singlemedium.root".format(i=source_directory_data, p=process_BKG, y2=year_last_two_digits_str, y4=year))

for process in (["data"] + processes_BKG):
    command_to_run = "./getMCNorms/bin/{c}".format(c=inputArguments.histCategory)
    command_to_run += " inputPaths={i}".format(i=(",".join(sources[inputArguments.histCategory][process])))
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
