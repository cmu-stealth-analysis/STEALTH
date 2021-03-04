#!/usr/bin/env python

from __future__ import print_function, division

import os

os.system("rm -vf *.pyc")

import sys, signal, argparse
import tmMultiProcessLauncher # from tmPyUtils
import stealthEnv # from this folder

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Run analysis chain.')
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArgumentsParser.add_argument('--selectionSuffix', default="", help='If set, the input n-tuples are read with this suffix.',type=str)
inputArgumentsParser.add_argument('--chain', default="all", help="Chain to run: can be \"data\", \"MC\", \"combine\", \"ancillaryPlots\", or \"limits\". Default: \"all\", which runs everything except \"limits\" (because that needs the condor jobs for the combine tool to finish).",type=str)
inputArgumentsParser.add_argument('--year', default="all", help="Year of data-taking. Can be \"2016\", \"2017\", \"2018\", or (default) \"all\".", type=str)
inputArgumentsParser.add_argument('--runUnblinded', action='store_true', help="If this flag is set, then the signal region data is unblinded.")
inputArgumentsParser.add_argument('--noLooseSignal', action='store_true', help="Do not add loose photons in a different signal bin. Run on a single signal type. By default data cards are created with signal, loose signal, and control selections.")
inputArgumentsParser.add_argument('--isDryRun', action='store_true', help="Only print the commands to run, do not actually run them.")
inputArguments = inputArgumentsParser.parse_args()

if not((inputArguments.chain == "data") or
       (inputArguments.chain == "MC") or
       (inputArguments.chain == "combine") or
       (inputArguments.chain == "ancillaryPlots") or
       (inputArguments.chain == "limits") or
       (inputArguments.chain == "all")):
    inputArgumentsParser.print_help()
    sys.exit("Unexpected value for argument \"chain\": {a}. See help message above.".format(a=inputArguments.chain))

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)
analysisOutputDirectory = "{aR}/analysis{oI}".format(aR=stealthEnv.analysisRoot, oI=optional_identifier)
analysisEOSOutputDirectory = "{sER}/analysisEOSAreas/analysis{oI}".format(sER=stealthEnv.stealthEOSRoot, oI=optional_identifier)
combineResultsEOSOutputDirectory = "{sER}/combineToolOutputs/combineResults{oI}".format(sER=stealthEnv.stealthEOSRoot, oI=optional_identifier)
analysisLogsDirectory = "{aOD}/analysisLogs".format(aOD=analysisOutputDirectory)
slopeAdjustmentsDirectory = "{aR}/normBinOptimization_doublephoton".format(aR=stealthEnv.analysisRoot)

selection_suffix = ""
if (inputArguments.selectionSuffix != ""): selection_suffix = "_{sS}".format(sS=inputArguments.selectionSuffix)

multiProcessLauncher = None
def removeLock():
    global multiProcessLauncher
    if not(multiProcessLauncher is None): multiProcessLauncher.killAll()
    os.system("rm -f {aLD}/analysisRunning.lock".format(aLD=analysisLogsDirectory))

def checkAndEstablishLock(): # Make sure that at most one instance is running at a time
    global multiProcessLauncher
    if (os.path.isfile("{aLD}/analysisRunning.lock".format(aLD=analysisLogsDirectory))):
        sys.exit("ERROR: only one instance of analysis running script can run at a time! If you're sure this is not a problem, remove this file: {aLD}/analysisRunning.lock".format(aLD=analysisLogsDirectory))
    else:
        command_createAnalysisParentDirectory = "mkdir -p {aOD}".format(aOD=analysisOutputDirectory)
        stealthEnv.execute_in_env(commandToRun=command_createAnalysisParentDirectory, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)
        for outputSubdirectory in ["dataEventHistograms", "dataSystematics", "MCEventHistograms", "MCSystematics", "signalContamination", "publicationPlots", "limits", "analysisLogs"]:
            command_createAnalysisSubdirectory = "mkdir -p {aOD}/{oS}".format(aOD=analysisOutputDirectory, oS=outputSubdirectory)
            stealthEnv.execute_in_env(commandToRun=command_createAnalysisSubdirectory, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)
        command_createEOSAnalysisArea = ("eos {eP} mkdir -p {aEOD}".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory))
        stealthEnv.execute_in_env(commandToRun=command_createEOSAnalysisArea, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)
        os.system("touch {aLD}/analysisRunning.lock".format(aLD=analysisLogsDirectory))
        multiProcessLauncher = tmMultiProcessLauncher.tmMultiProcessLauncher(logOutputFolder=analysisLogsDirectory, printDebug=True)
checkAndEstablishLock()

def signal_handler(sig, frame):
    removeLock()
    sys.exit("Terminated by user.")
signal.signal(signal.SIGINT, signal_handler)

# L_total = L_2016 + L_2017 + L_2018 => deltaL_total/L_total = (deltaL_2016/L_2016)*(L_2016/L_total) + (deltaL_2017/L_2017)*(L_2017/L_total) + (deltaL_2018/L_2018)*(L_2018/L_total)
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-17-001/index.html, the 2016 uncertainty is 2.5 percent
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-17-004/index.html, the 2017 uncertainty is 2.3 percent
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-18-002/index.html, the 2018 uncertainty is 2.5 percent

lumi_uncertainty = 0.025
integrated_lumi_strings = {
    "2016": "35918.2",
    "2017": "41527.3",
    "2018": "59736.0"
}

if (inputArguments.chain == "all"):
    runSequence = ["data", "MC", "combine", "ancillaryPlots"] # "limits" should be run manually at the end once all the combine jobs have finished running
else:
    runSequence = [inputArguments.chain]

yearPattern = inputArguments.year
if (yearPattern == "all"):
    yearPattern = "201*"

list_signalTypes = ["signal"]
if not(inputArguments.noLooseSignal):
    list_signalTypes.append("signal_loose")

eventProgenitors = ["gluino", "squark"]
crossSectionsForProgenitor = {
    "gluino": "SusyCrossSections13TevGluGlu.txt",
    "squark": "SusyCrossSections13TevSquarkSquark.txt"
}
tDesignationsForProgenitor = {
    "gluino": "t5",
    "squark": "t6"
}
MCTemplatesForProgenitor = {
    "gluino": "{eP}/{sER}/MCGeneratedMasses/MCGeneratedMasses/MCGeneratedMasses_stealth_t5Wg_savedObjects.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot),
    "squark": "{eP}/{sER}/MCGeneratedMasses/MCGeneratedMasses/MCGeneratedMasses_stealth_t6Wg_savedObjects.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot)
}
minNeutralinoMassToPlot = {
    "gluino": -1.,
    "squark": 195.
}

def transfer_file_to_EOS_area(sourceFile, targetDirectory):
    sourceFileName = (sourceFile.split("/"))[-1]
    command = "xrdcp --verbose --force --path --streams 15 {sF} {tD}/{sFN}".format(sF=sourceFile, tD=targetDirectory, sFN=sourceFileName)
    stealthEnv.execute_in_env(commandToRun=command, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)

def get_commands_data_chain(inputFilesList, outputPrefix, analyzeSignalBins):
    command_getEventHistogramsAndSystematics = ("./getDataEventHistogramsAndSystematics.py --inputFilesList {iFL} --outputDirectory_eventHistograms {aOD}/dataEventHistograms/ --outputDirectory_dataSystematics {aOD}/dataSystematics/ --outputPrefix {oP}".format(iFL=inputFilesList, aOD=analysisOutputDirectory, oP=outputPrefix))
    if (analyzeSignalBins): command_getEventHistogramsAndSystematics += " --analyzeSignalBins"
    return [command_getEventHistogramsAndSystematics]

def get_commands_MC_chain(eventProgenitor, dataPrefix, outputPrefix, inputMCPathMain, inputDataPUSourceMain, PUWeightsOutputPathMain, inputHLTEfficienciesPathMain, integratedLuminosityMainString, inputMCPathsAux, inputDataPUSourcesAux, PUWeightsOutputPathsAux, inputHLTEfficienciesPathsAux, integratedLuminositiesAux, getSignalContaminationOutsideSidebands):
    commands_MC_chain = []
    command_getPUWeightsMain = ("./getMCSystematics/bin/makePUWeights inputDataPath={iDPUSM} inputMCPath={iMCPM} outputFolder={eP}/{aEOD} outputFileName={PUWOPM}".format(iDPUSM=inputDataPUSourceMain, iMCPM=inputMCPathMain, eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory, PUWOPM=PUWeightsOutputPathMain))
    commands_MC_chain.append(command_getPUWeightsMain)
    for auxIndex in range(len(inputDataPUSourcesAux)):
        inputDataPUSource = inputDataPUSourcesAux[auxIndex]
        inputMCPathAux = inputMCPathsAux[auxIndex]
        PUWeightsOutputPathAux = PUWeightsOutputPathsAux[auxIndex]
        command_getPUWeightsAux = ("./getMCSystematics/bin/makePUWeights inputDataPath={iDPUS} inputMCPath={iMCPA} outputFolder={eP}/{aEOD} outputFileName={PUWOPA}".format(iDPUS=inputDataPUSource, iMCPA=inputMCPathAux, eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory, PUWOPA=PUWeightsOutputPathAux))
        commands_MC_chain.append(command_getPUWeightsAux)

    command_getHists = ("./getMCSystematics/bin/getEventHistograms eventProgenitor={eP} crossSectionsFilePath={cSFP} inputMCPathMain={iMCPM} inputHLTEfficienciesPathMain={iHLTEPM} integratedLuminosityMain={iLM} outputDirectory={aOD}/MCEventHistograms/ outputPrefix={oP} MCTemplatePath={MTP} inputPUWeightsPathMain={eP2}/{aEOD}/{iPUWPM}".format(eP=eventProgenitor, eP2=stealthEnv.EOSPrefix, cSFP=crossSectionsForProgenitor[eventProgenitor], iMCPM=inputMCPathMain, iHLTEPM=inputHLTEfficienciesPathMain, iLM=integratedLuminosityMainString, aOD=analysisOutputDirectory, aEOD=analysisEOSOutputDirectory, oP=outputPrefix, MTP=MCTemplatesForProgenitor[eventProgenitor], iPUWPM=PUWeightsOutputPathMain))
    if (len(inputMCPathsAux) != 0):
        command_getHists += " inputMCPathsAux="
        for inputMCPathAux in inputMCPathsAux:
            command_getHists += (inputMCPathAux + "\;")
        command_getHists = command_getHists[:-2] # to remove the last "\;"
        command_getHists += " inputHLTEfficienciesPathsAux="
        for inputHLTEfficienciesPathAux in inputHLTEfficienciesPathsAux:
            command_getHists += (inputHLTEfficienciesPathAux + "\;")
        command_getHists = command_getHists[:-2] # to remove the last "\;"
        command_getHists += " integratedLuminositiesAux="
        for integratedLuminosityAux in integratedLuminositiesAux:
            command_getHists += (integratedLuminosityAux + "\;")
        command_getHists = command_getHists[:-2] # to remove the last "\;"
        command_getHists += " inputPUWeightsPathsAux="
        for inputPUWeightPathAux in PUWeightsOutputPathsAux:
            command_getHists += ("{eP}/{aEOD}/".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory) + inputPUWeightPathAux + "\;")
        command_getHists = command_getHists[:-2] # to remove the last "\;"
    commands_MC_chain.append(command_getHists)

    signalContaminationOutsideSidebandsString = "false" # this is a string, not a bool
    if getSignalContaminationOutsideSidebands:
        signalContaminationOutsideSidebandsString = "true"
    command_getSystematics = ("./getMCSystematics/bin/getMCUncertainties inputPath={aOD}/MCEventHistograms/{oP}_savedObjects.root MCTemplatePath={MTP} inputNEventsFile={aOD}/dataSystematics/{dP}_observedEventCounters.dat inputDataUncertaintiesFile={aOD}/dataSystematics/{dP}_dataSystematics.dat outputDirectory={aOD}/MCSystematics/ outputDirectory_signalContamination={aOD}/signalContamination/ outputPrefix={oP} getSignalContaminationOutsideSidebands={sCOSS}".format(aOD=analysisOutputDirectory, oP=outputPrefix, MTP=MCTemplatesForProgenitor[eventProgenitor], dP=dataPrefix, sCOSS=signalContaminationOutsideSidebandsString))
    commands_MC_chain.append(command_getSystematics)
    return commands_MC_chain

def run_combine_chain(eventProgenitor, path_dataSystematics_signal, path_dataSystematics_signal_loose, path_dataSystematics_control, path_dataSystematics_scaling_signal, path_dataSystematics_scaling_signal_loose, path_dataSystematics_scaling_control, path_dataSystematics_scalingQuality, path_slopeAdjustment_signal, path_slopeAdjustment_signal_loose, path_slopeAdjustment_control, path_dataObservedEventCounters_signal, path_dataObservedEventCounters_signal_loose, path_dataObservedEventCounters_control, path_dataExpectedEventCounters_signal, path_dataExpectedEventCounters_signal_loose, path_dataExpectedEventCounters_control, MCPrefix_signal, MCPrefix_signal_loose, MCPrefix_control, outputPrefix):
    command_submitCombineJobs = ("./submitCombineToolJobs.py --dataCardsPrefix {oP} --outputDirectory {eP}/{cREOSOD}/ --eventProgenitor {eP2} --MCTemplatePath {MTP} --minNeutralinoMass {mNM} --crossSectionsFileName {cSFP} --path_dataSystematics_signal {pDSS} --path_dataSystematics_signal_loose {pDSSL} --path_dataSystematics_control {pDSC} --path_dataSystematics_scaling_signal {pDSSS} --path_dataSystematics_scaling_signal_loose {pDSSSL} --path_dataSystematics_scaling_control {pDSSC} --path_dataSystematics_scalingQuality {pDSSQ} --path_slopeAdjustment_signal {pSAS} --path_slopeAdjustment_signal_loose {pSASL} --path_slopeAdjustment_control {pSAC} --path_dataObservedEventCounters_signal {pDOECS} --path_dataObservedEventCounters_signal_loose {pDOECSL} --path_dataObservedEventCounters_control {pDOECC} --path_dataExpectedEventCounters_signal {pDEECS} --path_dataExpectedEventCounters_signal_loose {pDEECSL} --path_dataExpectedEventCounters_control {pDEECC} --MCHistogramsSignal {EAA}/{pS}_savedObjects.root --MCHistogramsSignalLoose {EAA}/{pSL}_savedObjects.root --MCHistogramsControl {EAA}/{pC}_savedObjects.root --MCUncertaintiesSignal {EAA}/{pS}_MCUncertainties_savedObjects.root --MCUncertaintiesSignalLoose {EAA}/{pSL}_MCUncertainties_savedObjects.root --MCUncertaintiesControl {EAA}/{pC}_MCUncertainties_savedObjects.root --luminosityUncertainty {lU} --EOSAnalysisArea {EAA}".format(oP=outputPrefix, eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, MTP=MCTemplatesForProgenitor[eventProgenitor], mNM=minNeutralinoMassToPlot[eventProgenitor], eP2=eventProgenitor, cREOSOD=combineResultsEOSOutputDirectory, cSFP=crossSectionsForProgenitor[eventProgenitor], pDSS=path_dataSystematics_signal, pDSSL=path_dataSystematics_signal_loose, pDSC=path_dataSystematics_control, pDSSS=path_dataSystematics_scaling_signal, pDSSSL=path_dataSystematics_scaling_signal_loose, pDSSC=path_dataSystematics_scaling_control, pDSSQ=path_dataSystematics_scalingQuality, pSAS=path_slopeAdjustment_signal, pSASL=path_slopeAdjustment_signal_loose, pSAC=path_slopeAdjustment_control, pDOECS=path_dataObservedEventCounters_signal, pDOECSL=path_dataObservedEventCounters_signal_loose, pDOECC=path_dataObservedEventCounters_control, pDEECS=path_dataExpectedEventCounters_signal, pDEECSL=path_dataExpectedEventCounters_signal_loose, pDEECC=path_dataExpectedEventCounters_control, pS=MCPrefix_signal, pSL=MCPrefix_signal_loose, pC=MCPrefix_control, lU=lumi_uncertainty, EAA="{eP}/{aEOD}".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory)))
    if (inputArguments.optionalIdentifier != ""): command_submitCombineJobs += " --optionalIdentifier {oI}".format(oI=inputArguments.optionalIdentifier) # Just "inputArguments.optionalIdentifier", without the underscore
    if (inputArguments.runUnblinded): command_submitCombineJobs += " --runUnblinded"
    if not(inputArguments.noLooseSignal): command_submitCombineJobs += " --addLooseSignal"
    if (inputArguments.isDryRun): command_submitCombineJobs += " --isDryRun"
    stealthEnv.execute_in_env(commandToRun=command_submitCombineJobs, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)

def produce_STComparisons(dataPath, outputFilePrefix_STComparisons, analyzeSignalBins, useWeights):
    command_produceComparisons = "./plotSTDistributionComparisons.py --inputFilePath {dP} --outputDirectory_plots {aOD}/publicationPlots --outputDirectory_dataSystematics {aOD}/dataSystematics --outputFilePrefix {oFP}".format(dP=dataPath, aOD=analysisOutputDirectory, oFP=outputFilePrefix_STComparisons)
    if not(analyzeSignalBins): command_produceComparisons += " --nJetsMaxPlot 3"
    if (useWeights): command_produceComparisons += " --weighted"
    stealthEnv.execute_in_env(commandToRun=command_produceComparisons, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)

def produce_ancillary_plots_control(eventProgenitor, path_data_expectedNEvents, path_data_observedNEvents, path_MC_weightedNEvents, path_dataSystematics):
    command_controlSTDistributions_dataAndSignal = "./plotSTDistributionsWithErrors.py --eventProgenitor {eP} --path_data_expectedNEvents {pDENE} --path_data_observedNEvents {pDONE} --path_MC_weightedNEvents {pMCWNE} --path_dataSystematics {pDS} --outputDirectory {aOD}/publicationPlots/ --outputFilePrefix {oFP} --plotObservedData".format(eP=eventProgenitor, pDENE=path_data_expectedNEvents, pDONE=path_data_observedNEvents, pMCWNE=path_MC_weightedNEvents, pDS=path_dataSystematics, aOD=analysisOutputDirectory, oFP="STDistributions_{eP}_control".format(eP=eventProgenitor))
    stealthEnv.execute_in_env(commandToRun=command_controlSTDistributions_dataAndSignal, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)

def produce_ancillary_plots_signal(eventProgenitor, signalType, path_data_expectedNEvents, path_data_observedNEvents, path_MC_weightedNEvents, path_dataSystematics):
    command_signalSTDistributions_dataAndSignal = "./plotSTDistributionsWithErrors.py --eventProgenitor {eP} --path_data_expectedNEvents {pDENE} --path_data_observedNEvents {pDONE} --path_MC_weightedNEvents {pMCWNE} --path_dataSystematics {pDS} --outputDirectory {aOD}/publicationPlots/ --outputFilePrefix {oFP}".format(eP=eventProgenitor, pDENE=path_data_expectedNEvents, pDONE=path_data_observedNEvents, pMCWNE=path_MC_weightedNEvents, pDS=path_dataSystematics, aOD=analysisOutputDirectory, oFP="STDistributions_{eP}_{sT}".format(eP=eventProgenitor, sT=signalType))
    if inputArguments.runUnblinded: command_signalSTDistributions_dataAndSignal += " --plotObservedData"
    stealthEnv.execute_in_env(commandToRun=command_signalSTDistributions_dataAndSignal, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)

def plot_limits(eventProgenitor):
    selectionsString = "signal"
    if not(inputArguments.noLooseSignal): selectionsString += ",signal_loose"
    command_plotLimits = "condor_q && ./plotLimits.py --crossSectionsFile {cSFP} --MCTemplatePath {MTP} --eventProgenitor {eP2} --combineResultsDirectory {eP}/{cREOSOD} --combineOutputPrefix {eP2} --outputDirectory_rawOutput {aOD}/limits --outputDirectory_plots {aOD}/publicationPlots --outputSuffix {eP2} --minNeutralinoMass {mNM} --selectionsList {sS}".format(eP=stealthEnv.EOSPrefix, MTP=MCTemplatesForProgenitor[eventProgenitor], cREOSOD=combineResultsEOSOutputDirectory, eP2=eventProgenitor, cSFP=crossSectionsForProgenitor[eventProgenitor], aOD=analysisOutputDirectory, mNM=minNeutralinoMassToPlot[eventProgenitor], sS=selectionsString)
    if (inputArguments.runUnblinded): command_plotLimits += " --plotObserved"
    stealthEnv.execute_in_env(commandToRun=command_plotLimits, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)

for step in runSequence:
    if (step == "data"):
        for signalType in (list_signalTypes + ["control"]):
            produce_STComparisons(dataPath="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_GJet17_2017_{sT}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, sT=signalType), outputFilePrefix_STComparisons="{sT}_GJet17_STComparisons".format(sT=signalType), analyzeSignalBins=True, useWeights=True)
            produce_STComparisons(dataPath="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_QCD_2017_{sT}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, sT=signalType), outputFilePrefix_STComparisons="{sT}_QCD_STComparisons".format(sT=signalType), analyzeSignalBins=True, useWeights=True)
        shellCommands_control = get_commands_data_chain(inputFilesList="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_data_{yP}_control.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, yP=yearPattern), outputPrefix="control", analyzeSignalBins=True)
        if (inputArguments.isDryRun): print("Not spawning due to dry run flag: {sC_c}".format(sC_c=shellCommands_control))
        else: multiProcessLauncher.spawn(shellCommands=shellCommands_control, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step1_data_control.log", printDebug=True)
        for signalType in list_signalTypes:
            shellCommands_signal = get_commands_data_chain(inputFilesList="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_data_{yP}_{sT}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, yP=yearPattern, sT=signalType), outputPrefix="{sT}".format(sT=signalType), analyzeSignalBins=inputArguments.runUnblinded)
            if (inputArguments.isDryRun): print("Not spawning due to dry run flag: {sC_s}".format(sC_s=shellCommands_signal))
            else: multiProcessLauncher.spawn(shellCommands=shellCommands_signal, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step1_data_{sT}.log".format(sT=signalType), printDebug=True)
        if not(inputArguments.isDryRun): multiProcessLauncher.monitorToCompletion()
    elif (step == "MC"):
        command_update = ("cd getMCSystematics && make && cd ..")
        stealthEnv.execute_in_env(commandToRun=command_update, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)
        for eventProgenitor in eventProgenitors:
            crossSectionsPath = crossSectionsForProgenitor[eventProgenitor]
            for signalType in (list_signalTypes + ["control"]):
                MCPathMain = ""
                dataPUSourceMain = ""
                HLTEfficienciesPathMain = ""
                PUWeightsOutputPathMain = ""
                lumiMain = ""
                MCPathsAux = []
                dataPUSourcesAux = []
                PUWeightsOutputPathsAux = []
                HLTEfficienciesPathsAux = []
                lumisAux = []
                if (inputArguments.year == "all"):
                    MCPathMain = "{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_{tD}_2017_{sT}.root".format(sT=signalType, eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, tD=tDesignationsForProgenitor[eventProgenitor])
                    dataPUSourceMain = "getMCSystematics/data/dataPU_2017.root"
                    HLTEfficienciesPathMain = "{eP}/{sER}/HLTEfficiencies/HLTEfficiencies_{sT}_2017.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sT=signalType)
                    PUWeightsOutputPathMain = "PUWeights_2017_{eP}_{sT}.root".format(eP=eventProgenitor, sT=signalType)
                    lumiMain = integrated_lumi_strings["2017"]
                    MCPathsAux = ["{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_{tD}_2016_{sT}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, sT=signalType, tD=tDesignationsForProgenitor[eventProgenitor]), "{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_{tD}_2018_{sT}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, sT=signalType, tD=tDesignationsForProgenitor[eventProgenitor])]
                    dataPUSourcesAux = ["getMCSystematics/data/dataPU_2016.root", "getMCSystematics/data/dataPU_2018.root"]
                    PUWeightsOutputPathsAux = ["PUWeights_2016_{eP}_{sT}.root".format(eP=eventProgenitor, sT=signalType), "PUWeights_2018_{eP}_{sT}.root".format(eP=eventProgenitor, sT=signalType)]
                    HLTEfficienciesPathsAux = ["{eP}/{sER}/HLTEfficiencies/HLTEfficiencies_{sT}_2016.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sT=signalType), "{eP}/{sER}/HLTEfficiencies/HLTEfficiencies_{sT}_2018.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sT=signalType)]
                    lumisAux = [integrated_lumi_strings["2016"], integrated_lumi_strings["2018"]]
                else:
                    MCPathMain = "{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_{tD}_{y}_{sT}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, y=inputArguments.year, sT=signalType, tD=tDesignationsForProgenitor[eventProgenitor])
                    dataPUSourceMain = "getMCSystematics/data/dataPU_{y}.root".format(y=inputArguments.year)
                    HLTEfficienciesPathMain = "{eP}/{sER}/HLTEfficiencies/HLTEfficiencies_{sT}_{y}.root".format(sT=signalType, y=inputArguments.year)
                    PUWeightsOutputPathMain = "PUWeights_{y}_{eP}_{sT}.root".format(y=inputArguments.year, eP=eventProgenitor, sT=signalType)
                    lumiMain = integrated_lumi_strings[inputArguments.year]
                get_signalContamination_outside_sidebands = False
                if (signalType == "control"): get_signalContamination_outside_sidebands = True
                shellCommands_MC = get_commands_MC_chain(eventProgenitor=eventProgenitor, dataPrefix=signalType, outputPrefix="MC_stealth_{eP}_{y}_{sT}".format(eP=eventProgenitor, y=inputArguments.year, sT=signalType), inputMCPathMain=MCPathMain, inputDataPUSourceMain=dataPUSourceMain, PUWeightsOutputPathMain=PUWeightsOutputPathMain, inputHLTEfficienciesPathMain=HLTEfficienciesPathMain, integratedLuminosityMainString=lumiMain, inputMCPathsAux=MCPathsAux, inputDataPUSourcesAux=dataPUSourcesAux, PUWeightsOutputPathsAux=PUWeightsOutputPathsAux, inputHLTEfficienciesPathsAux=HLTEfficienciesPathsAux, integratedLuminositiesAux=lumisAux, getSignalContaminationOutsideSidebands=get_signalContamination_outside_sidebands)
                if (inputArguments.isDryRun): print("Not spawning due to dry run flag: {sC_MC}".format(sC_MC=shellCommands_MC))
                else: multiProcessLauncher.spawn(shellCommands=shellCommands_MC, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step2_MC_{eP}_{sT}.log".format(eP=eventProgenitor, sT=signalType), printDebug=True)
        if not(inputArguments.isDryRun): multiProcessLauncher.monitorToCompletion()
        for eventProgenitor in eventProgenitors:
            for signalType in (list_signalTypes + ["control"]):
                transfer_file_to_EOS_area(sourceFile="{aOD}/MCEventHistograms/{oP}_savedObjects.root".format(aOD=analysisOutputDirectory, oP="MC_stealth_{eP}_{y}_{sT}".format(eP=eventProgenitor, y=inputArguments.year, sT=signalType)), targetDirectory="{eP}/{sER}/analysisEOSAreas/analysis{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier))
                transfer_file_to_EOS_area(sourceFile="{aOD}/MCSystematics/{oP}_MCUncertainties_savedObjects.root".format(aOD=analysisOutputDirectory, oP="MC_stealth_{eP}_{y}_{sT}".format(eP=eventProgenitor, y=inputArguments.year, sT=signalType)), targetDirectory="{eP}/{sER}/analysisEOSAreas/analysis{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier))
    elif (step == "combine"):
        # Make sure the CMSSW source tarball is the latest version
        print("Updating and uploading CMSSW source tarball...")
        updateCommand = "cd {sCB}/.. && ./uploadTarball.sh && cd {sR}".format(sCB=stealthEnv.stealthCMSSWBase, sR=stealthEnv.stealthRoot)
        os.system(updateCommand)
        command_cleanCombineResultsDirectory = "eos {eP} ls {sER}/combineToolOutputs/combineResults{oI}; eos {eP} rm -r {sER}/combineToolOutputs/combineResults{oI} || echo \"combine results output directory does not exist.\"".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier)
        stealthEnv.execute_in_env(commandToRun=command_cleanCombineResultsDirectory, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)
        command_createEOSDirectory = ("eos {eP} mkdir -p {sER}/combineToolOutputs/combineResults{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier))
        stealthEnv.execute_in_env(commandToRun=command_createEOSDirectory, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)
        for eventProgenitor in eventProgenitors:
            crossSectionsPath = crossSectionsForProgenitor[eventProgenitor]
            run_combine_chain(eventProgenitor=eventProgenitor, path_dataSystematics_signal="{aOD}/dataSystematics/signal_dataSystematics.dat".format(aOD=analysisOutputDirectory), path_dataSystematics_signal_loose="{aOD}/dataSystematics/signal_loose_dataSystematics.dat".format(aOD=analysisOutputDirectory), path_dataSystematics_control="{aOD}/dataSystematics/control_dataSystematics.dat".format(aOD=analysisOutputDirectory), path_dataSystematics_scaling_signal="{aOD}/dataSystematics/signal_GJet17_STComparisons_scalingSystematics.dat".format(aOD=analysisOutputDirectory), path_dataSystematics_scaling_signal_loose="{aOD}/dataSystematics/signal_loose_GJet17_STComparisons_scalingSystematics.dat".format(aOD=analysisOutputDirectory), path_dataSystematics_scaling_control="{aOD}/dataSystematics/control_GJet17_STComparisons_scalingSystematics.dat".format(aOD=analysisOutputDirectory), path_dataSystematics_scalingQuality="{aOD}/../scalingQualitySystematics/scalingQualitySystematics_combined.dat".format(aOD=analysisOutputDirectory), path_slopeAdjustment_signal="{sAD}/fitParameters_all_MC_GJet_signal.dat".format(sAD=slopeAdjustmentsDirectory), path_slopeAdjustment_signal_loose="{sAD}/fitParameters_all_MC_GJet_signal_loose.dat".format(sAD=slopeAdjustmentsDirectory), path_slopeAdjustment_control="{sAD}/fitParameters_all_MC_GJet_control.dat".format(sAD=slopeAdjustmentsDirectory), path_dataObservedEventCounters_signal="{aOD}/dataSystematics/signal_observedEventCounters.dat".format(aOD=analysisOutputDirectory), path_dataObservedEventCounters_signal_loose="{aOD}/dataSystematics/signal_loose_observedEventCounters.dat".format(aOD=analysisOutputDirectory), path_dataObservedEventCounters_control="{aOD}/dataSystematics/control_observedEventCounters.dat".format(aOD=analysisOutputDirectory), path_dataExpectedEventCounters_signal="{aOD}/dataSystematics/signal_eventCounters.dat".format(aOD=analysisOutputDirectory), path_dataExpectedEventCounters_signal_loose="{aOD}/dataSystematics/signal_loose_eventCounters.dat".format(aOD=analysisOutputDirectory), path_dataExpectedEventCounters_control="{aOD}/dataSystematics/control_eventCounters.dat".format(aOD=analysisOutputDirectory), MCPrefix_signal="MC_stealth_{eP}_{y}_signal".format(eP=eventProgenitor, y=inputArguments.year), MCPrefix_signal_loose="MC_stealth_{eP}_{y}_signal_loose".format(eP=eventProgenitor, y=inputArguments.year), MCPrefix_control="MC_stealth_{eP}_{y}_control".format(eP=eventProgenitor, y=inputArguments.year), outputPrefix="{eP}".format(eP=eventProgenitor))
    elif (step == "ancillaryPlots"):
        produce_STComparisons(dataPath="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_data_{yP}_control.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, yP=yearPattern), outputFilePrefix_STComparisons="control_STComparisons", analyzeSignalBins=True, useWeights=False)
        for signalType in list_signalTypes:
            produce_STComparisons(dataPath="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_data_{yP}_{sT}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, yP=yearPattern, sT=signalType), outputFilePrefix_STComparisons="{sT}_STComparisons".format(sT=signalType), analyzeSignalBins=inputArguments.runUnblinded, useWeights=False)
        for eventProgenitor in eventProgenitors:
            crossSectionsPath = crossSectionsForProgenitor[eventProgenitor]
            produce_ancillary_plots_control(eventProgenitor=eventProgenitor, path_data_expectedNEvents="{aOD}/dataSystematics/control_eventCounters.dat".format(aOD=analysisOutputDirectory), path_data_observedNEvents="{aOD}/dataSystematics/control_observedEventCounters.dat".format(aOD=analysisOutputDirectory), path_MC_weightedNEvents="{aOD}/MCEventHistograms/MC_stealth_{eP}_{y}_control_savedObjects.root".format(aOD=analysisOutputDirectory, eP=eventProgenitor, y=inputArguments.year), path_dataSystematics="{aOD}/dataSystematics/control_dataSystematics.dat".format(aOD=analysisOutputDirectory))
            for signalType in list_signalTypes:
                produce_ancillary_plots_signal(eventProgenitor=eventProgenitor, signalType=signalType, path_data_expectedNEvents="{aOD}/dataSystematics/{sT}_eventCounters.dat".format(aOD=analysisOutputDirectory, sT=signalType), path_data_observedNEvents="{aOD}/dataSystematics/{sT}_observedEventCounters.dat".format(aOD=analysisOutputDirectory, sT=signalType), path_MC_weightedNEvents="{aOD}/MCEventHistograms/MC_stealth_{eP}_{y}_{sT}_savedObjects.root".format(aOD=analysisOutputDirectory, eP=eventProgenitor, y=inputArguments.year, sT=signalType), path_dataSystematics="{aOD}/dataSystematics/{sT}_dataSystematics.dat".format(aOD=analysisOutputDirectory, sT=signalType))
    elif (step == "limits"):
        for eventProgenitor in eventProgenitors:
            plot_limits(eventProgenitor)
    else:
        removeLock()
        sys.exit("ERROR: Unrecognized step: {s}".format(s=step))

removeLock()
