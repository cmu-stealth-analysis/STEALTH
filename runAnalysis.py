#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, signal, argparse, re, subprocess, time, stealthEnv

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination for the control region.')
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArgumentsParser.add_argument('--chain', default="all", help="Chain to run: can be \"data\", \"MC\", or \"combine\". Default: \"all\".",type=str)
inputArguments = inputArgumentsParser.parse_args()

if not((inputArguments.chain == "data") or
       (inputArguments.chain == "MC") or
       (inputArguments.chain == "combine") or
       (inputArguments.chain == "all")):
    inputArgumentsParser.print_help()
    sys.exit("Unexpected value for argument \"chain\": {a}. See help message above.".format(a=inputArguments.chain))

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)

# L_total = L_2016 + L_2017 => deltaL_total/L_total = (deltaL_2016/L_2016)*(L_2016/L_total) + (deltaL_2017/L_2017)*(L_2017/L_total)
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-17-004/index.html, the 2017 uncertainty is 2.3 percent
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-17-001/index.html, the 2016 uncertainty is 2.5 percent
# L_total = L_2016 + L_2017 => deltaL_total/L_total = (deltaL_2016/L_2016)*(L_2016/L_total) + (deltaL_2017/L_2017)*(L_2017/L_total)
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-17-004/index.html, the 2017 uncertainty is 2.3 percent
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-17-001/index.html, the 2016 uncertainty is 2.5 percent
lumi_uncertainty = 0.024
integrated_lumi_2017_string = "41900.0"
hltefficiency_pattern_leading = "hltEfficiency_leadingPhoton_control_fakefake"
hltefficiency_pattern_subLeading = "hltEfficiency_subLeadingPhoton_control_fakefake"

def checkAndEstablishLock(optional_identifier):
    if (os.path.isfile("analysis{oI}_running.lock".format(oI=optional_identifier))):
        sys.exit("ERROR: analysis{oI}_running.lock already exists!".format(oI=optional_identifier))
    os.system("touch analysis{oI}_running.lock".format(oI=optional_identifier))

def removeLock(optional_identifier):
    os.system("rm -f analysis{oI}_running.lock".format(oI=optional_identifier))

def signal_handler(sig, frame):
    removeLock(optional_identifier)
    sys.exit("Terminated by user.")
signal.signal(signal.SIGINT, signal_handler)

def execute_in_env(commandToRun, optional_identifier, printDebug=False):
    env_setup_command = "bash -c \"cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot)
    runInEnv = "{e_s_c} && set -x && {c} && set +x\"".format(e_s_c=env_setup_command, c=commandToRun)
    if (printDebug):
        print("About to execute command:")
        print("{c}".format(c=runInEnv))
    returnCode = os.system(runInEnv)
    if (returnCode != 0):
        removeLock(optional_identifier)
        sys.exit("ERROR: command \"{c}\" returned status {rC}".format(c=commandToRun, rC=returnCode))

def run_data_step(outputDirectory, inputFilesList, outputPrefix, isSignal, optional_identifier):
    for outputSubdirectory in ["dataEventHistograms", "dataSystematics"]:
        os.system("mkdir -p {oD}/{oS}".format(oD=outputDirectory, oS=outputSubdirectory))
    command = ("./getDataEventHistogramsAndSystematics.py --inputFilesList {iFL} --outputDirectory_eventHistograms {oD}/dataEventHistograms/ --outputDirectory_dataSystematics {oD}/dataSystematics/ --outputPrefix {oP}".format(iFL=inputFilesList, oD=outputDirectory, oP=outputPrefix))
    if (isSignal): command += " --isSignal"
    execute_in_env(command, optional_identifier)

def run_MC_chain(outputDirectory, dataPrefix, outputPrefix, inputMCPathMain, integratedLuminosityMainString, HLTEfficiencySources, MCTemplatePath, allowUnrestrictedSignalContamination, optional_identifier):
    for outputSubdirectory in ["MCEventHistograms", "MCSystematics", "signalContamination"]:
        os.system("mkdir -p {oD}/{oS}".format(oD=outputDirectory, oS=outputSubdirectory))
    command_update = ("cd getMCSystematics && make && cd ..")
    execute_in_env(command_update, optional_identifier)
    command_getHists = ("./getMCSystematics/bin/getEventHistograms inputMCPathMain={iMCPM} integratedLuminosityMain={iLM} outputDirectory={oD}/MCEventHistograms/ outputPrefix={oP} HLTEfficiencySources={HES} MCTemplatePath={MTP}".format(iMCPM=inputMCPathMain, iLM=integratedLuminosityMainString, oD=outputDirectory, oP=outputPrefix, HES=HLTEfficiencySources, MTP=MCTemplatePath))
    execute_in_env(command_getHists, optional_identifier)
    unrestrictedSignalContamination = "false" # the string, not the bool
    if allowUnrestrictedSignalContamination: unrestrictedSignalContamination = "true"
    command_getSystematics = ("./getMCSystematics/bin/getMCUncertainties inputPath={oD}/MCEventHistograms/{oP}_savedObjects.root MCTemplatePath={MTP} inputNEventsFile={oD}/dataSystematics/{dP}_observedEventCounters.dat outputDirectory={oD}/MCSystematics/ outputDirectory_signalContamination={oD}/signalContamination/ outputPrefix={oP} unrestrictedSignalContamination={uSC}".format(oD=outputDirectory, oP=outputPrefix, MTP=MCTemplatePath, dP=dataPrefix, uSC=unrestrictedSignalContamination))
    execute_in_env(command_getSystematics, optional_identifier)

def run_combine_chain(outputDirectory, outputPrefix, MCTemplatePath, luminosity_uncertainty, optional_identifier):
    os.system("mkdir -p {oD}/dataCards".format(oD=outputDirectory))
    command_createCards = ("./createDataCards.py --outputPrefix {oP} --outputDirectory {oD}/dataCards/ --MCTemplatePath {MTP} --inputFile_MCEventHistograms {oD}/MCEventHistograms/{oP}_savedObjects.root --inputFile_MCUncertainties {oD}/MCSystematics/{oP}_MCUncertainties_savedObjects.root --inputFile_dataSystematics {oD}/dataSystematics/signal_dataSystematics.dat --inputFile_dataSystematics_sTScaling {oD}/dataSystematics/control_dataSystematics_sTScaling.dat --inputFile_dataSystematics_expectedEventCounters {oD}/dataSystematics/signal_eventCounters.dat --inputFile_dataSystematics_observedEventCounters {oD}/dataSystematics/signal_observedEventCounters.dat --luminosityUncertainty {lU}".format(oP=outputPrefix, oD=outputDirectory, MTP=MCTemplatePath, lU=luminosity_uncertainty))
    execute_in_env(command_createCards, optional_identifier)
    command_updateEOSDirectory = ("eos root://cmseos.fnal.gov ls {sER}/combineToolOutputs && eos root://cmseos.fnal.gov rm -r {sER}/combineToolOutputs && eos root://cmseos.fnal.gov mkdir -p {sER}/combineToolOutputs".format(sER=stealthEnv.stealthEOSRoot))
    execute_in_env(command_updateEOSDirectory)
    command_submitCombineJobs = ("./submitCombineToolJobs.py --dataCardsPrefix {oP} --outputDirectory {eP}/{sER}/combineToolOutputs/".format(oP=outputPrefix, eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot))
    execute_in_env(command_submitCombineJobs)

checkAndEstablishLock(optional_identifier)

if (inputArguments.chain == "all"):
    runSequence = ["data", "MC", "combine"]
else:
    runSequence = [inputArguments.chain]

for step in runSequence:
    if (step == "data"):
        run_data_step(outputDirectory="analysis{oI}".format(oI=optional_identifier), inputFilesList="{eP}/{sER}/selections/combined_DoublePhoton{oI}/merged_selection_data_2017_control_fakefake.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), outputPrefix="control", isSignal=False, optional_identifier=optional_identifier)
        run_data_step(outputDirectory="analysis{oI}".format(oI=optional_identifier), inputFilesList="{eP}/{sER}/selections/combined_DoublePhoton{oI}/merged_selection_data_2017_signal.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), outputPrefix="signal", isSignal=True, optional_identifier=optional_identifier)
    elif (step == "MC"):
        run_MC_chain(outputDirectory="analysis{oI}".format(oI=optional_identifier), dataPrefix="signal", outputPrefix="MC_stealth_t5_2017", inputMCPathMain="{eP}/{sER}/selections/combined_DoublePhoton{oI}/merged_selection_MC_stealth_t5_2017_signal.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), integratedLuminosityMainString="{iL}".format(iL=integrated_lumi_2017_string), HLTEfficiencySources="{eP}/{sER}/statistics/combined_DoublePhoton{oI}/merged_statistics_MC_stealth_t5_2017.root:{hltEPL}:{hltEPsL}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier, hltEPL=hltefficiency_pattern_leading, hltEPsL=hltefficiency_pattern_subLeading), MCTemplatePath="{eP}/{sER}/MCGeneratedMasses/MCGeneratedMasses/MCGeneratedMasses_mc_Fall17_stealth_t5Wg_savedObjects.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot), allowUnrestrictedSignalContamination=False, optional_identifier=optional_identifier)
    elif (step == "combine"):
        run_combine_chain(outputDirectory="analysis{oI}".format(oI=optional_identifier), outputPrefix="fullChain", MCTemplatePath="{eP}/{sER}/MCGeneratedMasses/MCGeneratedMasses/MCGeneratedMasses_mc_Fall17_stealth_t5Wg_savedObjects.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot), luminosity_uncertainty=lumi_uncertainty, optional_identifier=optional_identifier)

removeLock(optional_identifier)
