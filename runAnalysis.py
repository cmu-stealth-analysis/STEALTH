#!/usr/bin/env python

from __future__ import print_function, division

import os

os.system("rm -vf *.pyc")

import sys, signal, argparse, re, subprocess, time, stealthEnv

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination for the control region.')
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArgumentsParser.add_argument('--selectionSuffix', default="", help='If set, the input n-tuples are read with this suffix.',type=str)
inputArgumentsParser.add_argument('--chain', default="all", help="Chain to run: can be \"data\", \"MC\", \"combine\", \"signalContamination\", \"ancillaryPlots\", or \"limits\". Default: \"all\", which runs everything except \"limits\" (because that needs the condor jobs for the combine tool to finish).",type=str)
inputArgumentsParser.add_argument('--year', default="all", help="Year of data-taking. Can be \"2016\", \"2017\", or (default) \"all\".", type=str)
inputArgumentsParser.add_argument('--runUnblinded', action='store_true', help="If this flag is set, then the signal region data is unblinded.")
inputArguments = inputArgumentsParser.parse_args()

if not((inputArguments.chain == "data") or
       (inputArguments.chain == "MC") or
       (inputArguments.chain == "combine") or
       (inputArguments.chain == "signalContamination") or
       (inputArguments.chain == "ancillaryPlots") or
       (inputArguments.chain == "limits") or
       (inputArguments.chain == "all")):
    inputArgumentsParser.print_help()
    sys.exit("Unexpected value for argument \"chain\": {a}. See help message above.".format(a=inputArguments.chain))

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)

selection_suffix = ""
if (inputArguments.selectionSuffix != ""): optional_identifier = "_{sS}".format(sS=inputArguments.selectionSuffix)

# L_total = L_2016 + L_2017 => deltaL_total/L_total = (deltaL_2016/L_2016)*(L_2016/L_total) + (deltaL_2017/L_2017)*(L_2017/L_total)
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-17-004/index.html, the 2017 uncertainty is 2.3 percent
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-17-001/index.html, the 2016 uncertainty is 2.5 percent
# L_total = L_2016 + L_2017 => deltaL_total/L_total = (deltaL_2016/L_2016)*(L_2016/L_total) + (deltaL_2017/L_2017)*(L_2017/L_total)
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-17-004/index.html, the 2017 uncertainty is 2.3 percent
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-17-001/index.html, the 2016 uncertainty is 2.5 percent
lumi_uncertainty = 0.024
integrated_lumi_strings = {
    "2016": "35920.0",
    "2017": "41900.0"
}

def checkAndEstablishLock(optional_identifier):
    if (os.path.isfile("analysis{oI}_running.lock".format(oI=optional_identifier))):
        sys.exit("ERROR: analysis{oI}_running.lock already exists!".format(oI=optional_identifier))
    os.system("touch analysis{oI}_running.lock".format(oI=optional_identifier))

def removeLock(optional_identifier):
    os.system("rm -v -f analysis{oI}_running.lock".format(oI=optional_identifier))

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

def run_data_step(outputDirectory, inputFilesList, outputPrefix, getSystematics, getSTScalingSystematics, analyzeSignalBins, optional_identifier):
    for outputSubdirectory in ["dataEventHistograms", "dataSystematics"]:
        os.system("mkdir -p {oD}/{oS}".format(oD=outputDirectory, oS=outputSubdirectory))
    command = ("./getDataEventHistogramsAndSystematics.py --inputFilesList {iFL} --outputDirectory_eventHistograms {oD}/dataEventHistograms/ --outputDirectory_dataSystematics {oD}/dataSystematics/ --outputPrefix {oP}".format(iFL=inputFilesList, oD=outputDirectory, oP=outputPrefix))
    if (getSystematics): command += " --getSystematics"
    if (getSTScalingSystematics): command += " --getSTScalingSystematics"
    if (analyzeSignalBins): command += " --analyzeSignalBins"
    execute_in_env(command, optional_identifier)

def run_MC_chain(outputDirectory, dataPrefix, outputPrefix, inputMCPathMain, integratedLuminosityMainString, HLTEfficiencySourceMain, inputMCPathsAux, integratedLuminositiesAux, HLTEfficiencySourcesAux, MCTemplatePath, getSignalContaminationOutsideSidebands, optional_identifier):
    for outputSubdirectory in ["MCEventHistograms", "MCSystematics", "signalContamination"]:
        os.system("mkdir -p {oD}/{oS}".format(oD=outputDirectory, oS=outputSubdirectory))
    command_update = ("cd getMCSystematics && make && cd ..")
    execute_in_env(command_update, optional_identifier)
    command_getHists = ("./getMCSystematics/bin/getEventHistograms inputMCPathMain={iMCPM} integratedLuminosityMain={iLM} HLTEfficiencySourceMain={HESM} outputDirectory={oD}/MCEventHistograms/ outputPrefix={oP} MCTemplatePath={MTP}".format(iMCPM=inputMCPathMain, iLM=integratedLuminosityMainString, oD=outputDirectory, oP=outputPrefix, HESM=HLTEfficiencySourceMain, MTP=MCTemplatePath))
    if (len(inputMCPathsAux) != 0):
        command_getHists += " inputMCPathsAux="
        for inputMCPathAux in inputMCPathsAux:
            command_getHists += (inputMCPathAux + ";")
        command_getHists = command_getHists[:-1] # to remove the last ";"
        command_getHists += " integratedLuminositiesAux="
        for integratedLuminosityAux in integratedLuminositiesAux:
            command_getHists += (integratedLuminosityAux + ";")
        command_getHists = command_getHists[:-1] # to remove the last ";"
        command_getHists += " HLTEfficiencySourcesAux="
        for HLTEfficiencySourceAux in HLTEfficiencySourcesAux:
            command_getHists += (HLTEfficiencySourceAux + ";")
        command_getHists = command_getHists[:-1] # to remove the last ";"
    execute_in_env(command_getHists, optional_identifier)
    signalContaminationOutsideSidebandsString = "false" # the string, not the bool
    if getSignalContaminationOutsideSidebands:
        signalContaminationOutsideSidebandsString = "true"
    command_getSystematics = ("./getMCSystematics/bin/getMCUncertainties inputPath={oD}/MCEventHistograms/{oP}_savedObjects.root MCTemplatePath={MTP} inputNEventsFile={oD}/dataSystematics/{dP}_observedEventCounters.dat outputDirectory={oD}/MCSystematics/ outputDirectory_signalContamination={oD}/signalContamination/ outputPrefix={oP} getSignalContaminationOutsideSidebands={sCOSS}".format(oD=outputDirectory, oP=outputPrefix, MTP=MCTemplatePath, dP=dataPrefix, sCOSS=signalContaminationOutsideSidebandsString))
    execute_in_env(command_getSystematics, optional_identifier)

def run_combine_chain(outputDirectory, combineResultsDirectory, prefix_MCChainStep, outputPrefix, MCTemplatePath, luminosity_uncertainty, runUnblinded, optional_identifier):
    os.system("mkdir -p {oD}/dataCards".format(oD=outputDirectory))
    command_createCards = ("./createDataCards.py --outputPrefix {oP} --outputDirectory {oD}/dataCards/ --MCTemplatePath {MTP} --inputFile_MCEventHistograms {oD}/MCEventHistograms/{pMCS}_savedObjects.root --inputFile_MCUncertainties {oD}/MCSystematics/{pMCS}_MCUncertainties_savedObjects.root --inputFile_dataSystematics {oD}/dataSystematics/signal_dataSystematics.dat --inputFile_dataSystematics_sTScaling {oD}/dataSystematics/control_dataSystematics_sTScaling.dat --inputFile_dataSystematics_expectedEventCounters {oD}/dataSystematics/signal_eventCounters.dat --inputFile_dataSystematics_observedEventCounters {oD}/dataSystematics/signal_observedEventCounters.dat --luminosityUncertainty {lU}".format(oP=outputPrefix, pMCS=prefix_MCChainStep, oD=outputDirectory, MTP=MCTemplatePath, lU=luminosity_uncertainty))
    if (runUnblinded): command_createCards += " --runUnblinded"
    execute_in_env(command_createCards, optional_identifier)
    os.system("eos {eP} ls {sER}/combineToolOutputs/{cRD} && eos {eP} rm -r {sER}/combineToolOutputs/{cRD}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, cRD=combineResultsDirectory))
    command_updateEOSDirectory = ("eos {eP} mkdir -p {sER}/combineToolOutputs/{cRD}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, cRD=combineResultsDirectory))
    execute_in_env(command_updateEOSDirectory, optional_identifier)
    command_submitCombineJobs = ("./submitCombineToolJobs.py --dataCardsDirectory {oD}/dataCards/ --dataCardsPrefix {oP} --outputDirectory {eP}/{sER}/combineToolOutputs/{cRD}/ --MCTemplatePath {MTP}".format(oD=outputDirectory, oP=outputPrefix, eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, MTP=MCTemplatePath, cRD=combineResultsDirectory, oI=inputArguments.optionalIdentifier))
    if (inputArguments.optionalIdentifier != ""): command_submitCombineJobs += " --optionalIdentifier {oI}".format(oI=inputArguments.optionalIdentifier) # Just "inputArguments.optionalIdentifier", without the underscore
    execute_in_env(command_submitCombineJobs, optional_identifier)

def produce_ancillary_plots(outputDirectory, controlDataPath, outputFileName_STComparisons, path_data_expectedNEvents_control, path_data_observedNEvents_control, path_MC_weightedNEvents_control, path_dataSystematics_control, path_STScalingSystematics, path_data_expectedNEvents_signal, path_data_observedNEvents_signal, path_MC_weightedNEvents_signal, path_dataSystematics_signal, runUnblinded, optional_identifier):
    for outputSubdirectory in ["dataEventHistograms", "dataSystematics", "publicationPlots"]:
        os.system("mkdir -p {oD}/{oS}".format(oD=outputDirectory, oS=outputSubdirectory))

    # For signal:
    command_signalSTDistributions_dataAndSignal = "./plotSTDistributionsWithErrors.py --path_data_expectedNEvents {pDENE} --path_data_observedNEvents {pDONE} --path_MC_weightedNEvents {pMCWNE} --path_dataSystematics {pDS} --path_STScalingSystematics {pSTSS} --outputDirectory {oD}/publicationPlots/ --outputFilePrefix {oFP}".format(pDENE=path_data_expectedNEvents_signal, pDONE=path_data_observedNEvents_signal, pMCWNE=path_MC_weightedNEvents_signal, pDS=path_dataSystematics_signal, pSTSS=path_STScalingSystematics, oD=outputDirectory, oFP="STDistributions_signal")
    if runUnblinded:
        command_signalSTDistributions_dataAndSignal += " --plotObservedData"
    execute_in_env(command_signalSTDistributions_dataAndSignal, optional_identifier)

    # For control:
    command_controlDistributions = "./plotSTDistributionComparisons.py --inputFilePath {cDP} --outputDirectory {oD}/publicationPlots --outputFileName {oFN}".format(cDP=controlDataPath, oD=outputDirectory, oFN=outputFileName_STComparisons)
    execute_in_env(command_controlDistributions, optional_identifier)
    command_controlSystematics = ("./getDataEventHistogramsAndSystematics.py --inputFilesList {cDP} --outputDirectory_eventHistograms {oD}/dataEventHistograms/ --outputDirectory_dataSystematics {oD}/dataSystematics/ --outputPrefix control --getSystematics --analyzeSignalBins".format(cDP=controlDataPath, oD=outputDirectory))
    execute_in_env(command_controlSystematics, optional_identifier) # To create control region systematics with the same method as signal systematics; required as input for the next step
    command_controlSTDistributions_dataAndSignal = "./plotSTDistributionsWithErrors.py --path_data_expectedNEvents {pDENE} --path_data_observedNEvents {pDONE} --path_MC_weightedNEvents {pMCWNE} --path_dataSystematics {pDS} --path_STScalingSystematics {pSTSS} --outputDirectory {oD}/publicationPlots/ --outputFilePrefix {oFP} --plotObservedData".format(pDENE=path_data_expectedNEvents_control, pDONE=path_data_observedNEvents_control, pMCWNE=path_MC_weightedNEvents_control, pDS=path_dataSystematics_control, pSTSS=path_STScalingSystematics, oD=outputDirectory, oFP="STDistributions_control")
    execute_in_env(command_controlSTDistributions_dataAndSignal, optional_identifier)

def get_signal_contamination(outputDirectory, dataPrefix, outputPrefix, inputMCPathMain, integratedLuminosityMainString, HLTEfficiencySourceMain, inputMCPathsAux, integratedLuminositiesAux, HLTEfficiencySourcesAux, MCTemplatePath, optional_identifier):
    for outputSubdirectory in ["MCEventHistograms", "MCSystematics", "signalContamination"]:
        os.system("mkdir -p {oD}/{oS}".format(oD=outputDirectory, oS=outputSubdirectory))
    command_update = ("cd getMCSystematics && make && cd ..")
    execute_in_env(command_update, optional_identifier)
    command_getHists = ("./getMCSystematics/bin/getEventHistograms inputMCPathMain={iMCPM} integratedLuminosityMain={iLM} HLTEfficiencySourceMain={HESM} outputDirectory={oD}/MCEventHistograms/ outputPrefix={oP} MCTemplatePath={MTP}".format(iMCPM=inputMCPathMain, iLM=integratedLuminosityMainString, HESM=HLTEfficiencySourceMain, oD=outputDirectory, oP=outputPrefix, MTP=MCTemplatePath))
    if (len(inputMCPathsAux) != 0):
        command_getHists += " inputMCPathsAux="
        for inputMCPathAux in inputMCPathsAux:
            command_getHists += (inputMCPathAux + ";")
        command_getHists = command_getHists[:-1] # to remove the last ";"
        command_getHists += " integratedLuminositiesAux="
        for integratedLuminosityAux in integratedLuminositiesAux:
            command_getHists += (integratedLuminosityAux + ";")
        command_getHists = command_getHists[:-1] # to remove the last ";"
        command_getHists += " HLTEfficiencySourcesAux="
        for HLTEfficiencySourceAux in HLTEfficiencySourcesAux:
            command_getHists += (HLTEfficiencySourceAux + ";")
        command_getHists = command_getHists[:-1] # to remove the last ";"
    execute_in_env(command_getHists, optional_identifier)
    command_getSignalContamination = ("./getMCSystematics/bin/getMCUncertainties inputPath={oD}/MCEventHistograms/{oP}_savedObjects.root MCTemplatePath={MTP} inputNEventsFile={oD}/dataSystematics/{dP}_observedEventCounters.dat outputDirectory={oD}/MCSystematics/ outputDirectory_signalContamination={oD}/signalContamination/ outputPrefix={oP} getSignalContaminationOutsideSidebands=true".format(oD=outputDirectory, oP=outputPrefix, MTP=MCTemplatePath, dP=dataPrefix))
    execute_in_env(command_getSignalContamination, optional_identifier)

def plot_limits(outputDirectory, combineResultsDirectory, MCTemplatePath, runUnblinded, optional_identifier):
    os.system("mkdir -p {oD}/publicationPlots".format(oD=outputDirectory))
    command_plotLimits = "condor_q && ./plotLimits.py --MCTemplatePath {MTP} --combineResultsDirectory {eP}/{sER}/combineToolOutputs/{cRD} --combineOutputPrefix fullChain --outputDirectory_rawOutput limits --outputDirectory_plots {oD}/publicationPlots --outputSuffix fullChain".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, MTP=MCTemplatePath, cRD=combineResultsDirectory, oD=outputDirectory)
    if (runUnblinded): command_plotLimits += " --plotObserved"
    execute_in_env(command_plotLimits, optional_identifier)

checkAndEstablishLock(optional_identifier)

if (inputArguments.chain == "all"):
    runSequence = ["data", "MC", "combine", "signalContamination", "ancillaryPlots"] # "limits" should be run manually at the end once all the combine jobs have finished running
else:
    runSequence = [inputArguments.chain]

yearPattern = inputArguments.year
if (yearPattern == "all"):
    yearPattern = "*"

for step in runSequence:
    if (step == "data"):
        run_data_step(outputDirectory="analysis{oI}".format(oI=optional_identifier), inputFilesList="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_data_{yP}_control_fakefake.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, yP=yearPattern), outputPrefix="control", getSystematics=False, getSTScalingSystematics=True, analyzeSignalBins=True, optional_identifier=optional_identifier)
        run_data_step(outputDirectory="analysis{oI}".format(oI=optional_identifier), inputFilesList="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_data_{yP}_signal.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, yP=yearPattern), outputPrefix="signal", getSystematics=True, getSTScalingSystematics=False, analyzeSignalBins=inputArguments.runUnblinded, optional_identifier=optional_identifier)
    elif (step == "MC"):
        hltefficiency_pattern_leading = "hltEfficiency_leadingPhoton_signal"
        hltefficiency_pattern_subLeading = "hltEfficiency_subLeadingPhoton_signal"
        MCPathMain = ""
        lumiMain = ""
        HLTEfficiencySourceMain = ""
        MCPathsAux = []
        lumisAux = []
        HLTEfficiencySourcesAux = []
        if (inputArguments.year == "all"):
            MCPathMain = "{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_t5_2017_signal.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix)
            lumiMain = integrated_lumi_strings["2017"]
            HLTEfficiencySourceMain = "{eP}/{sER}/statistics/combined_DoublePhoton{sS}/merged_statistics_MC_stealth_t5_2017.root\|{hltEPL}\|{hltEPsL}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, hltEPL=hltefficiency_pattern_leading, hltEPsL=hltefficiency_pattern_subLeading)
            MCPathsAux = ["{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_t5_2016_signal.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix)]
            lumisAux = [integrated_lumi_strings["2016"]]
            HLTEfficiencySourcesAux = ["{eP}/{sER}/statistics/combined_DoublePhoton{sS}/merged_statistics_MC_stealth_t5_2016.root\|{hltEPL}\|{hltEPsL}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, hltEPL=hltefficiency_pattern_leading, hltEPsL=hltefficiency_pattern_subLeading)]
        else:
            MCPathMain = "{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_t5_{y}_signal.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, y=inputArguments.year)
            lumiMain = integrated_lumi_strings[inputArguments.year]
            HLTEfficiencySourceMain = "{eP}/{sER}/statistics/combined_DoublePhoton{sS}/merged_statistics_MC_stealth_t5_{y}.root\|{hltEPL}\|{hltEPsL}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, hltEPL=hltefficiency_pattern_leading, hltEPsL=hltefficiency_pattern_subLeading, y=inputArguments.year)
        run_MC_chain(outputDirectory="analysis{oI}".format(oI=optional_identifier), dataPrefix="signal", outputPrefix="MC_stealth_t5_{y}_signal".format(y=inputArguments.year), inputMCPathMain=MCPathMain, integratedLuminosityMainString=lumiMain, HLTEfficiencySourceMain=HLTEfficiencySourceMain, inputMCPathsAux=MCPathsAux, integratedLuminositiesAux=lumisAux, HLTEfficiencySourcesAux=HLTEfficiencySourcesAux, MCTemplatePath="{eP}/{sER}/MCGeneratedMasses/MCGeneratedMasses/MCGeneratedMasses_mc_Fall17_stealth_t5Wg_savedObjects.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot), getSignalContaminationOutsideSidebands=False, optional_identifier=optional_identifier)
    elif (step == "combine"):
        run_combine_chain(outputDirectory="analysis{oI}".format(oI=optional_identifier), combineResultsDirectory="{cRD}".format(cRD="combineResults{oI}".format(oI=optional_identifier)), prefix_MCChainStep="MC_stealth_t5_{y}_signal".format(y = inputArguments.year), outputPrefix="fullChain", MCTemplatePath="{eP}/{sER}/MCGeneratedMasses/MCGeneratedMasses/MCGeneratedMasses_mc_Fall17_stealth_t5Wg_savedObjects.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot), luminosity_uncertainty=lumi_uncertainty, runUnblinded=inputArguments.runUnblinded, optional_identifier=optional_identifier)
    elif (step == "ancillaryPlots"):
        produce_ancillary_plots(outputDirectory="analysis{oI}".format(oI=optional_identifier), controlDataPath="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_data_{yP}_control_fakefake.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, yP=yearPattern), outputFileName_STComparisons="control_STComparison.png", path_data_expectedNEvents_control="analysis{oI}/dataSystematics/control_eventCounters.dat".format(oI=optional_identifier), path_data_observedNEvents_control="analysis{oI}/dataSystematics/control_observedEventCounters.dat".format(oI=optional_identifier), path_MC_weightedNEvents_control="analysis{oI}/MCEventHistograms/MC_stealth_t5_{y}_control_savedObjects.root".format(oI=optional_identifier, y=inputArguments.year), path_dataSystematics_control="analysis{oI}/dataSystematics/control_dataSystematics.dat".format(oI=optional_identifier), path_STScalingSystematics="analysis{oI}/dataSystematics/control_dataSystematics_sTScaling.dat".format(oI=optional_identifier), path_data_expectedNEvents_signal="analysis{oI}/dataSystematics/signal_eventCounters.dat".format(oI=optional_identifier), path_data_observedNEvents_signal="analysis{oI}/dataSystematics/signal_observedEventCounters.dat".format(oI=optional_identifier), path_MC_weightedNEvents_signal="analysis{oI}/MCEventHistograms/MC_stealth_t5_{y}_signal_savedObjects.root".format(oI=optional_identifier, y=inputArguments.year), path_dataSystematics_signal="analysis{oI}/dataSystematics/signal_dataSystematics.dat".format(oI=optional_identifier), runUnblinded=inputArguments.runUnblinded, optional_identifier=optional_identifier)
    elif (step == "signalContamination"):
        hltefficiency_pattern_leading = "hltEfficiency_leadingPhoton_control_fakefake"
        hltefficiency_pattern_subLeading = "hltEfficiency_subLeadingPhoton_control_fakefake"
        MCPathMain = ""
        lumiMain = ""
        HLTEfficiencySourceMain = ""
        MCPathsAux = []
        lumisAux = []
        HLTEfficiencySourcesAux = []
        if (inputArguments.year == "all"):
            MCPathMain = "{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_t5_2017_control_fakefake.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix)
            lumiMain = integrated_lumi_strings["2017"]
            HLTEfficiencySourceMain = "{eP}/{sER}/statistics/combined_DoublePhoton{sS}/merged_statistics_MC_stealth_t5_2017.root\|{hltEPL}\|{hltEPsL}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, hltEPL=hltefficiency_pattern_leading, hltEPsL=hltefficiency_pattern_subLeading)
            MCPathsAux = ["{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_t5_2016_control_fakefake.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix)]
            lumisAux = [integrated_lumi_strings["2016"]]
            HLTEfficiencySourcesAux = ["{eP}/{sER}/statistics/combined_DoublePhoton{sS}/merged_statistics_MC_stealth_t5_2016.root\|{hltEPL}\|{hltEPsL}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, hltEPL=hltefficiency_pattern_leading, hltEPsL=hltefficiency_pattern_subLeading)]
        else:
            MCPathMain = "{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_t5_{y}_control_fakefake.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, y=inputArguments.year)
            lumiMain = integrated_lumi_strings[inputArguments.year]
            HLTEfficiencySourceMain = "{eP}/{sER}/statistics/combined_DoublePhoton{sS}/merged_statistics_MC_stealth_t5_{y}.root\|{hltEPL}\|{hltEPsL}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, hltEPL=hltefficiency_pattern_leading, hltEPsL=hltefficiency_pattern_subLeading, y=inputArguments.year)
        get_signal_contamination(outputDirectory="analysis{oI}".format(oI=optional_identifier), dataPrefix="control", outputPrefix="MC_stealth_t5_{y}_control".format(y=inputArguments.year), inputMCPathMain=MCPathMain, integratedLuminosityMainString=lumiMain, HLTEfficiencySourceMain=HLTEfficiencySourceMain, inputMCPathsAux=MCPathsAux, integratedLuminositiesAux=lumisAux, HLTEfficiencySourcesAux=HLTEfficiencySourcesAux, MCTemplatePath="{eP}/{sER}/MCGeneratedMasses/MCGeneratedMasses/MCGeneratedMasses_mc_Fall17_stealth_t5Wg_savedObjects.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot), optional_identifier=optional_identifier)
    elif (step == "limits"):
        # plot_limits(outputDirectory, combineResultsDirectory, MCTemplatePath, runUnblinded, optional_identifier)
        plot_limits(outputDirectory="analysis{oI}".format(oI=optional_identifier), combineResultsDirectory="{cRD}".format(cRD="combineResults{oI}".format(oI=optional_identifier)), MCTemplatePath="{eP}/{sER}/MCGeneratedMasses/MCGeneratedMasses/MCGeneratedMasses_mc_Fall17_stealth_t5Wg_savedObjects.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot), runUnblinded=inputArguments.runUnblinded, optional_identifier=optional_identifier)
    else:
        removeLock(optional_identifier)
        sys.exit("ERROR: Unrecognized step: {s}".format(s=step))

removeLock(optional_identifier)
