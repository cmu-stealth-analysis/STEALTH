#!/usr/bin/env python

from __future__ import print_function, division

import os

os.system("rm -vf *.pyc")

import sys, signal, argparse, json, subprocess
import tmMultiProcessLauncher, tmGeneralUtils # from tmPyUtils
import stealthEnv # from this folder

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Run analysis chain.')
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArgumentsParser.add_argument('--selectionSuffix', default="", help='If set, the input n-tuples are read with this suffix.',type=str)
inputArgumentsParser.add_argument('--chain', default="all", help="Chain to run: can be \"data\", \"BKGMC\", \"MC\", \"combine\", \"ancillaryPlots\", \"observations\", \"limits\", or \"statisticsChecks\". Default: \"all\", which runs everything except \"observations\", \"limits\", and \"statisticsChecks\" (because those require the condor jobs running the combine tool to finish).",type=str)
inputArgumentsParser.add_argument('--year', default="all", help="Year of data-taking. Can be \"2016\", \"2017\", \"2018\", or (default) \"all\".", type=str)
inputArgumentsParser.add_argument('--runUnblinded', action='store_true', help="If this flag is set, then the signal region data is unblinded.")
inputArgumentsParser.add_argument('--noLooseSignal', action='store_true', help="Do not add loose photons in a different signal bin. Run on a single signal type. By default data cards are created with signal, loose signal, and control selections.")
inputArgumentsParser.add_argument('--isDryRun', action='store_true', help="Only print the commands to run, do not actually run them.")
inputArguments = inputArgumentsParser.parse_args()

if not((inputArguments.chain == "data") or
       (inputArguments.chain == "BKGMC") or
       (inputArguments.chain == "MC") or
       (inputArguments.chain == "combine") or
       (inputArguments.chain == "ancillaryPlots") or
       (inputArguments.chain == "observations") or
       (inputArguments.chain == "limits") or
       (inputArguments.chain == "statisticsChecks") or
       (inputArguments.chain == "all")):
    inputArgumentsParser.print_help()
    sys.exit("Unexpected value for argument \"chain\": {a}. See help message above.".format(a=inputArguments.chain))

if (inputArguments.year != "all"): sys.exit("ERROR: Currently the BKGMC step only works if run on all years.")

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)
analysisOutputDirectory = "{aR}/analysis{oI}".format(aR=stealthEnv.analysisRoot, oI=optional_identifier)
analysisEOSOutputDirectory = "{sER}/analysisEOSAreas/analysis{oI}".format(sER=stealthEnv.stealthEOSRoot, oI=optional_identifier)
HLTEfficienciesSource = "{sER}/HLTEfficiencies{oI}".format(sER=stealthEnv.stealthEOSRoot, oI=optional_identifier)
combineResultsEOSOutputDirectory = "{sER}/combineToolOutputs/combineResults{oI}".format(sER=stealthEnv.stealthEOSRoot, oI=optional_identifier)
analysisLogsDirectory = "{aOD}/analysisLogs".format(aOD=analysisOutputDirectory)

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
        for outputSubdirectory in ["dataEventHistograms", "dataSystematics", "fits_doublephoton", "fits_doublephoton_decoupled", "fits_singlephoton", "MCEventHistograms", "MCSystematics", "signalContamination", "publicationPlots", "limits", "statisticsChecks", "analysisLogs"]:
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

def check_execution_statuses_manually(executionStatusFilesList):
    for executionStatusFilePath in executionStatusFilesList:
        with open(executionStatusFilePath, 'r') as executionStatusFileHandle:
            executionStatusString = executionStatusFileHandle.read()
            executionSuccessful = (executionStatusString == "0\n")
            if not(executionSuccessful):
                sys.exit("ERROR: status for file {fname} set to unsuccessful".format(fname=executionStatusFilePath))

lumi_uncertainty = 0.018
lumi_values_raw_json = None
with open("xSecLumiInfo/lumi_run2.json", 'r') as lumi_json_file_handle:
    lumi_values_raw_json = json.load(lumi_json_file_handle)
integrated_lumi_strings = {}
for year_string in ["2016", "2017", "2018"]:
    integrated_lumi_strings[year_string] = str(lumi_values_raw_json[year_string])

if (inputArguments.chain == "all"):
    runSequence = ["data", "BKGMC", "MC", "combine", "ancillaryPlots"] # "observations", "limits", and "statisticsChecks" should be run manually at the end once all the combine jobs have finished running
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
MCNominalAdjustmentFilePaths = {
    "signal": "{aOD}/fits_doublephoton/adjustments_all_MC_Bkg_signal.dat".format(aOD=analysisOutputDirectory),
    "signal_loose": "{aOD}/fits_doublephoton/adjustments_all_MC_Bkg_signal_loose.dat".format(aOD=analysisOutputDirectory),
    "control": "{aOD}/fits_doublephoton/adjustments_all_MC_Bkg_control.dat".format(aOD=analysisOutputDirectory)
}
# DataMCRatioAdjustmentsFilePaths = {
#     "QCD": {
#         "signal": "{aOD}/fits_doublephoton/ratio_adjustment_all_MC_QCD_Bkg_signal.dat".format(aOD=analysisOutputDirectory),
#         "signal_loose": "{aOD}/fits_doublephoton/ratio_adjustment_all_MC_QCD_Bkg_signal_loose.dat".format(aOD=analysisOutputDirectory),
#         "control": "{aOD}/fits_doublephoton/ratio_adjustment_all_MC_QCD_Bkg_control.dat".format(aOD=analysisOutputDirectory)
#     },
#     "diphoton": {
#         "signal": "{aOD}/fits_doublephoton/ratio_adjustment_all_MC_Diph_Bkg_signal.dat".format(aOD=analysisOutputDirectory),
#         "signal_loose": "{aOD}/fits_doublephoton/ratio_adjustment_all_MC_Diph_Bkg_signal_loose.dat".format(aOD=analysisOutputDirectory),
#         "control": "{aOD}/fits_doublephoton/ratio_adjustment_all_MC_Diph_Bkg_control.dat".format(aOD=analysisOutputDirectory)
#     }
# }
binLabelAbbreviations = {
    "signal": "s",
    "signal_loose": "l",
    "control": "c"
}
ratioRangesToPlot = {
    "pre": {
        "signal": {
            "min": 0.0,
            "max": 2.5
        },
        "signal_loose": {
            "min": 0.0,
            "max": 6.5
        }
    },
    "post": {
        "signal": {
            "min": 0.25,
            "max": 1.75
        },
        "signal_loose": {
            "min": 0.0,
            "max": 2.5
        }
    }
}
minNeutralinoMassToPlot = {
    "gluino": 120.,
    "squark": 195.
}
eventProgenitorMassOffsets = {
    "gluino": 225.,
    "squark": 275.
}
minMassDifferences = {
    "gluino": 100.0,
    "squark": 100.0
}

manualAdjustmentRatios = {
    "combined": {
        "signal_loose": (-0.5, 8.5)
    }
}

def transfer_file_to_EOS_area(sourceFile, targetDirectory):
    sourceFileName = (sourceFile.split("/"))[-1]
    command = "xrdcp --verbose --force --path --streams 15 {sF} {tD}/{sFN}".format(sF=sourceFile, tD=targetDirectory, sFN=sourceFileName)
    stealthEnv.execute_in_env(commandToRun=command, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)

def get_commands_data_chain(inputFilesList, outputPrefix, analyzeSignalBins):
    command_getEventHistogramsAndSystematics = ("./getDataEventHistogramsAndSystematics.py --inputFilesList {iFL} --outputDirectory_eventHistograms {aOD}/dataEventHistograms/ --outputDirectory_dataSystematics {aOD}/dataSystematics/ --outputPrefix {oP}".format(iFL=inputFilesList, aOD=analysisOutputDirectory, oP=outputPrefix))
    if (analyzeSignalBins): command_getEventHistogramsAndSystematics += " --analyzeSignalBins"
    return [command_getEventHistogramsAndSystematics]

def read_rho_nominal_from_file(rhoNominalFilePath=None):
    linesInFile = []
    with open(rhoNominalFilePath) as rhoNominalFileObject:
        for line in rhoNominalFileObject:
            linesInFile.append(line)
    if not(len(linesInFile) == 1):
        removeLock()
        sys.exit("ERROR: This file is in an unexpected format: {fp}".format(fp=rhoNominalFilePath))
    rho_nominal = float((linesInFile[0]).strip())
    return rho_nominal

def get_commands_doublephoton_BKGMC_chain(sourceData_BKGMC, readParametersExplicitlyFromSource, adjustmentPlots_min, adjustmentPlots_max, sourceData_data, identifier, outputFolder, selectionString, rhoNominal, disableStrictChecks):
    commands_doublephoton_BKGMC = []
    command_BKGMC_doublephoton = "./fitScripts/bin/runFits sourceData={sD} outputFolder={oF} selection={sS} getJECShiftedDistributions=true identifier={i} yearString=all STBoundariesSourceFile={sR}/STRegionBoundaries.dat PDF_nSTBins=25 rhoNominal={rN} adjustmentPlots_min={amin} adjustmentPlots_max={amax} minAllowedEMST=-1.0".format(sD=sourceData_BKGMC, oF=outputFolder, sS=selectionString, i=identifier, rN=rhoNominal, sR=stealthEnv.stealthRoot, amin=adjustmentPlots_min, amax=adjustmentPlots_max)
    if not(readParametersExplicitlyFromSource is None):
        command_BKGMC_doublephoton += " readParametersFromFiles={ps},{sR}/STRegionBoundaries.dat plotConcise=true".format(ps=readParametersExplicitlyFromSource, sR=stealthEnv.stealthRoot)
    if disableStrictChecks:
        command_BKGMC_doublephoton += " disableStrictChecks=true"
    commands_doublephoton_BKGMC.append(command_BKGMC_doublephoton)
    if not(sourceData_data is None):
        command_data_doublephoton = "./fitScripts/bin/runFits sourceData={sD} outputFolder={oF} selection={sS} getJECShiftedDistributions=false identifier=data yearString=all STBoundariesSourceFile={sR}/STRegionBoundaries.dat PDF_nSTBins=25 rhoNominal={rN} adjustmentPlots_min={amin} adjustmentPlots_max={amax} minAllowedEMST=-1.0 readParametersFromFiles={oF}/binned_fitParameters_all_MC_Bkg_{sS}.dat,{sR}/STRegionBoundaries.dat plotConcise=true".format(sD=sourceData_data, oF=outputFolder, sS=selectionString, rN=rhoNominal, amin=adjustmentPlots_min, amax=adjustmentPlots_max, sR=stealthEnv.stealthRoot)
        if disableStrictChecks:
            command_data_doublephoton += " disableStrictChecks=true"
        commands_doublephoton_BKGMC.append(command_data_doublephoton)
    return commands_doublephoton_BKGMC

def get_commands_singlephoton_BKGMC_chain(sourceData_BKGMC, adjustmentPlots_min, adjustmentPlots_max, sourceData_data, outputFolder, selectionString, yearString, rhoNominal, disableStrictChecks):
    commands_singlephoton_BKGMC = []
    command_BKGMC_singlephoton = "./fitScripts/bin/runFits sourceData={sD} outputFolder={oF} selection={sS} getJECShiftedDistributions=false identifier=MC_Bkg yearString={yS} STBoundariesSourceFile={sR}/STRegionBoundariesFineBinned.dat PDF_nSTBins=50 rhoNominal={rN} adjustmentPlots_min={amin} adjustmentPlots_max={amax} minAllowedEMST=200.0 plotConcise=true".format(sD=sourceData_BKGMC, oF=outputFolder, sS=selectionString, yS=yearString, rN=rhoNominal, amin=adjustmentPlots_min, amax=adjustmentPlots_max, sR=stealthEnv.stealthRoot)
    if disableStrictChecks:
        command_BKGMC_singlephoton += " disableStrictChecks=true"
    commands_singlephoton_BKGMC.append(command_BKGMC_singlephoton)
    command_data_singlephoton = "./fitScripts/bin/runFits sourceData={sD} outputFolder={oF} selection={sS} getJECShiftedDistributions=false identifier=data yearString={yS} STBoundariesSourceFile={sR}/STRegionBoundariesFineBinned.dat PDF_nSTBins=50 rhoNominal={rN} adjustmentPlots_min={amin} adjustmentPlots_max={amax} minAllowedEMST=200.0 readParametersFromFiles={oF}/binned_fitParameters_{yS}_MC_Bkg_{sS}.dat,{sR}/STRegionBoundaries.dat plotConcise=true".format(sD=sourceData_data, oF=outputFolder, sS=selectionString, yS=yearString, rN=rhoNominal, amin=adjustmentPlots_min, amax=adjustmentPlots_max, sR=stealthEnv.stealthRoot)
    if disableStrictChecks:
        command_data_singlephoton += " disableStrictChecks=true"
    commands_singlephoton_BKGMC.append(command_data_singlephoton)
    return commands_singlephoton_BKGMC

def get_commands_singlephoton_modulated_BKGMC_chain(sourceData_BKGMC, readParametersExplicitlyFromSource, identifier, adjustmentPlots_min, adjustmentPlots_max, outputFolder, selectionString, yearString, rhoNominal, disableStrictChecks):
    commands_singlephoton_BKGMC = []
    command_BKGMC_singlephoton = "./fitScripts/bin/runFits sourceData={sD} outputFolder={oF} selection={sS} getJECShiftedDistributions=false identifier={i} yearString={yS} STBoundariesSourceFile={sR}/STRegionBoundariesFineBinned.dat PDF_nSTBins=50 rhoNominal={rN} adjustmentPlots_min={amin} adjustmentPlots_max={amax} minAllowedEMST=200.0 plotConcise=true".format(sD=sourceData_BKGMC, oF=outputFolder, sS=selectionString, i=identifier, yS=yearString, rN=rhoNominal, amin=adjustmentPlots_min, amax=adjustmentPlots_max, sR=stealthEnv.stealthRoot)
    if disableStrictChecks:
        command_BKGMC_singlephoton += " disableStrictChecks=true"
    command_BKGMC_singlephoton += " readParametersFromFiles={ps},{sR}/STRegionBoundaries.dat plotConcise=true".format(ps=readParametersExplicitlyFromSource, sR=stealthEnv.stealthRoot)
    commands_singlephoton_BKGMC.append(command_BKGMC_singlephoton)
    return commands_singlephoton_BKGMC

def get_commands_MC_chain(eventProgenitor, dataPrefix, outputPrefix, inputMCPathMain, # inputDataPUSourceMain, PUWeightsOutputPathMain,
                          inputHLTEfficienciesPathMain, integratedLuminosityMainString, inputMCPathsAux, # inputDataPUSourcesAux, PUWeightsOutputPathsAux,
                          inputHLTEfficienciesPathsAux, integratedLuminositiesAux, getSignalContaminationOutsideSidebands):
    commands_MC_chain = []
    # command_getPUWeightsMain = ("./getPUWeights/bin/makePUWeights inputDataPath={iDPUSM} inputMCPaths={iMCPM} outputFolder={eP}/{aEOD} outputFileName={PUWOPM} addMCXSecWeight=false".format(iDPUSM=inputDataPUSourceMain, iMCPM=inputMCPathMain, eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory, PUWOPM=PUWeightsOutputPathMain))
    # commands_MC_chain.append(command_getPUWeightsMain)
    # for auxIndex in range(len(inputDataPUSourcesAux)):
    #     inputDataPUSource = inputDataPUSourcesAux[auxIndex]
    #     inputMCPathAux = inputMCPathsAux[auxIndex]
    #     PUWeightsOutputPathAux = PUWeightsOutputPathsAux[auxIndex]
    #     command_getPUWeightsAux = ("./getPUWeights/bin/makePUWeights inputDataPath={iDPUS} inputMCPaths={iMCPA} outputFolder={eP}/{aEOD} outputFileName={PUWOPA} addMCXSecWeight=false".format(iDPUS=inputDataPUSource, iMCPA=inputMCPathAux, eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory, PUWOPA=PUWeightsOutputPathAux))
    #     commands_MC_chain.append(command_getPUWeightsAux)

    # command_getHists = ("./getMCSystematics/bin/getEventHistograms eventProgenitor={eP} crossSectionsFilePath={cSFP} inputMCPathMain={iMCPM} inputHLTEfficienciesPathMain={iHLTEPM} integratedLuminosityMain={iLM} outputDirectory={aOD}/MCEventHistograms/ outputPrefix={oP} MCTemplatePath={MTP} inputPUWeightsPathMain={eP2}/{aEOD}/{iPUWPM}".format(eP=eventProgenitor, eP2=stealthEnv.EOSPrefix, cSFP=crossSectionsForProgenitor[eventProgenitor], iMCPM=inputMCPathMain, iHLTEPM=inputHLTEfficienciesPathMain, iLM=integratedLuminosityMainString, aOD=analysisOutputDirectory, aEOD=analysisEOSOutputDirectory, oP=outputPrefix, MTP=MCTemplatesForProgenitor[eventProgenitor], iPUWPM=PUWeightsOutputPathMain))
    command_getHists = ("./getMCSystematics/bin/getEventHistograms eventProgenitor={eP} crossSectionsFilePath={cSFP} inputMCPathMain={iMCPM} inputHLTEfficienciesPathMain={iHLTEPM} integratedLuminosityMain={iLM} outputDirectory={aOD}/MCEventHistograms/ outputPrefix={oP} MCTemplatePath={MTP}".format(eP=eventProgenitor, eP2=stealthEnv.EOSPrefix, cSFP=crossSectionsForProgenitor[eventProgenitor], iMCPM=inputMCPathMain, iHLTEPM=inputHLTEfficienciesPathMain, iLM=integratedLuminosityMainString, aOD=analysisOutputDirectory, aEOD=analysisEOSOutputDirectory, oP=outputPrefix, MTP=MCTemplatesForProgenitor[eventProgenitor]))
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
        # command_getHists += " inputPUWeightsPathsAux="
        # for inputPUWeightPathAux in PUWeightsOutputPathsAux:
        #     command_getHists += ("{eP}/{aEOD}/".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory) + inputPUWeightPathAux + "\;")
        # command_getHists = command_getHists[:-2] # to remove the last "\;"
    commands_MC_chain.append(command_getHists)

    signalContaminationOutsideSidebandsString = "false" # this is a string, not a bool
    if getSignalContaminationOutsideSidebands:
        signalContaminationOutsideSidebandsString = "true"
    command_getSystematics = ("./getMCSystematics/bin/getMCUncertainties inputPath={aOD}/MCEventHistograms/{oP}_savedObjects.root MCTemplatePath={MTP} inputNEventsFile={aOD}/dataSystematics/{dP}_observedEventCounters.dat inputDataUncertaintiesFile={aOD}/dataSystematics/{dP}_dataSystematics.dat outputDirectory={aOD}/MCSystematics/ outputDirectory_signalContamination={aOD}/signalContamination/ outputPrefix={oP} getSignalContaminationOutsideSidebands={sCOSS}".format(aOD=analysisOutputDirectory, oP=outputPrefix, MTP=MCTemplatesForProgenitor[eventProgenitor], dP=dataPrefix, sCOSS=signalContaminationOutsideSidebandsString))
    commands_MC_chain.append(command_getSystematics)
    return commands_MC_chain

def run_combine_chain(eventProgenitor, path_dataSystematics_signal, path_dataSystematics_signal_loose, # path_dataSystematics_control,
                      path_MCShapeAdjustment_signal, path_MCShapeAdjustment_signal_loose, # path_MCShapeAdjustment_control,
                      paths_bkgCompositionSystematic, path_dataObservedEventCounters_signal, path_dataObservedEventCounters_signal_loose, # path_dataObservedEventCounters_control,
                      path_dataExpectedEventCounters_signal, path_dataExpectedEventCounters_signal_loose, # path_dataExpectedEventCounters_control,
                      MCPrefix_signal, MCPrefix_signal_loose, # MCPrefix_control,
                      outputPrefix):
    # command_submitCombineJobs = ("./submitCombineToolJobs.py --dataCardsPrefix {oP} --outputDirectory {eP}/{cREOSOD}/ --eventProgenitor {eP2} --MCTemplatePath {MTP} --minNeutralinoMass {mNM} --eventProgenitorMassOffset {ePMO} --minMassDifference {mMD} --crossSectionsFileName {cSFP} --path_dataSystematics_signal {pDSS} --path_dataSystematics_signal_loose {pDSSL} --path_dataSystematics_control {pDSC} --path_MCShapeAdjustment_signal {pSAS} --path_MCShapeAdjustment_signal_loose {pSASL} --path_MCShapeAdjustment_control {pSAC} --paths_bkgCompositionSystematic {pBCS} --path_dataObservedEventCounters_signal {pDOECS} --path_dataObservedEventCounters_signal_loose {pDOECSL} --path_dataObservedEventCounters_control {pDOECC} --path_dataExpectedEventCounters_signal {pDEECS} --path_dataExpectedEventCounters_signal_loose {pDEECSL} --path_dataExpectedEventCounters_control {pDEECC} --MCHistogramsSignal {EAA}/{pS}_savedObjects.root --MCHistogramsSignalLoose {EAA}/{pSL}_savedObjects.root --MCHistogramsControl {EAA}/{pC}_savedObjects.root --MCUncertaintiesSignal {EAA}/{pS}_MCUncertainties_savedObjects.root --MCUncertaintiesSignalLoose {EAA}/{pSL}_MCUncertainties_savedObjects.root --MCUncertaintiesControl {EAA}/{pC}_MCUncertainties_savedObjects.root --luminosityUncertainty {lU} --EOSAnalysisArea {EAA}".format(oP=outputPrefix, eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, MTP=MCTemplatesForProgenitor[eventProgenitor], mNM=minNeutralinoMassToPlot[eventProgenitor], ePMO=eventProgenitorMassOffsets[eventProgenitor], mMD=minMassDifferences[eventProgenitor], eP2=eventProgenitor, cREOSOD=combineResultsEOSOutputDirectory, cSFP=crossSectionsForProgenitor[eventProgenitor], pDSS=path_dataSystematics_signal, pDSSL=path_dataSystematics_signal_loose, pDSC=path_dataSystematics_control, pSAS=path_MCShapeAdjustment_signal, pSASL=path_MCShapeAdjustment_signal_loose, pSAC=path_MCShapeAdjustment_control, pBCS=paths_bkgCompositionSystematic, pDOECS=path_dataObservedEventCounters_signal, pDOECSL=path_dataObservedEventCounters_signal_loose, pDOECC=path_dataObservedEventCounters_control, pDEECS=path_dataExpectedEventCounters_signal, pDEECSL=path_dataExpectedEventCounters_signal_loose, pDEECC=path_dataExpectedEventCounters_control, pS=MCPrefix_signal, pSL=MCPrefix_signal_loose, pC=MCPrefix_control, lU=lumi_uncertainty, EAA="{eP}/{aEOD}".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory)))
    command_submitCombineJobs = ("./submitCombineToolJobs.py --dataCardsPrefix {oP} --outputDirectory {eP}/{cREOSOD}/ --eventProgenitor {eP2} --MCTemplatePath {MTP} --minNeutralinoMass {mNM} --eventProgenitorMassOffset {ePMO} --minMassDifference {mMD} --crossSectionsFileName {cSFP} --path_dataSystematics_signal {pDSS} --path_dataSystematics_signal_loose {pDSSL} --path_MCShapeAdjustment_signal {pSAS} --path_MCShapeAdjustment_signal_loose {pSASL} --paths_bkgCompositionSystematic {pBCS} --path_dataObservedEventCounters_signal {pDOECS} --path_dataObservedEventCounters_signal_loose {pDOECSL} --path_dataExpectedEventCounters_signal {pDEECS} --path_dataExpectedEventCounters_signal_loose {pDEECSL} --MCHistogramsSignal {EAA}/{pS}_savedObjects.root --MCHistogramsSignalLoose {EAA}/{pSL}_savedObjects.root --MCUncertaintiesSignal {EAA}/{pS}_MCUncertainties_savedObjects.root --MCUncertaintiesSignalLoose {EAA}/{pSL}_MCUncertainties_savedObjects.root --luminosityUncertainty {lU} --EOSAnalysisArea {EAA}".format(oP=outputPrefix, eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, MTP=MCTemplatesForProgenitor[eventProgenitor], mNM=minNeutralinoMassToPlot[eventProgenitor], ePMO=eventProgenitorMassOffsets[eventProgenitor], mMD=minMassDifferences[eventProgenitor], eP2=eventProgenitor, cREOSOD=combineResultsEOSOutputDirectory, cSFP=crossSectionsForProgenitor[eventProgenitor], pDSS=path_dataSystematics_signal, pDSSL=path_dataSystematics_signal_loose, pSAS=path_MCShapeAdjustment_signal, pSASL=path_MCShapeAdjustment_signal_loose, pBCS=paths_bkgCompositionSystematic, pDOECS=path_dataObservedEventCounters_signal, pDOECSL=path_dataObservedEventCounters_signal_loose, pDEECS=path_dataExpectedEventCounters_signal, pDEECSL=path_dataExpectedEventCounters_signal_loose, pS=MCPrefix_signal, pSL=MCPrefix_signal_loose, lU=lumi_uncertainty, EAA="{eP}/{aEOD}".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory)))
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

# def get_commands_ancillary_plots_control(path_data_expectedNEvents, path_data_observedNEvents, path_data_adjustments, path_MC_weightedNEvents_gluino, path_MC_weightedNEvents_squark, path_systematics_nominal, path_systematics_dataMCDiscrepancy, nJetsBin):
#     commands_ancillary_plots_control = []
#     command_controlSTDistributions_dataAndSignal = "./plotSTDistributionsWithErrors.py --path_data_expectedNEvents {pDENE} --path_data_observedNEvents {pDONE} --path_data_adjustments {pDA} --path_MC_weightedNEvents_gluino {pMCWNEg} --path_MC_weightedNEvents_squark {pMCWNEs} --bin_label_abbreviation c --path_systematics_nominal {pSN} --path_systematics_dataMCDiscrepancy {pSDMCD} --outputDirectory {aOD}/publicationPlots/ --outputFilePrefix {oFP} --nJetsBin {n} --plotObservedData --suppressSignal".format(pDENE=path_data_expectedNEvents, pDONE=path_data_observedNEvents, pDA=path_data_adjustments, pMCWNEg=path_MC_weightedNEvents_gluino, pMCWNEs=path_MC_weightedNEvents_squark, pSN=path_systematics_nominal, pSDMCD=path_systematics_dataMCDiscrepancy, aOD=analysisOutputDirectory, oFP="STDistributions_control", n=nJetsBin)
#     # stealthEnv.execute_in_env(commandToRun=command_controlSTDistributions_dataAndSignal, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)
#     commands_ancillary_plots_control.append(command_controlSTDistributions_dataAndSignal)
#     return commands_ancillary_plots_control

def get_commands_ancillary_plots_signal(signalType, path_data_expectedNEvents, path_data_observedNEvents, path_data_adjustments, path_MC_weightedNEvents_gluino, path_MC_weightedNEvents_squark, path_systematics_nominal, inputFolder_bkgCompositionUncertainties, nJetsBin, bkgType, run_unblinded):
    commands_ancillary_plots_signal = []
    command_signalSTDistributions_dataAndSignal = "./plotSTDistributionsWithErrors.py --path_data_expectedNEvents {pDENE} --path_data_observedNEvents {pDONE} --path_data_adjustments {pDA} --path_MC_weightedNEvents_gluino {pMCWNEg} --path_MC_weightedNEvents_squark {pMCWNEs} --bin_label_abbreviation {bla} --path_systematics_nominal {pSN} --inputFolder_bkgCompositionUncertainties {iFBCU} --signalType {sT} --outputDirectory {aOD}/publicationPlots/ --nJetsBin {n}".format(pDENE=path_data_expectedNEvents, pDONE=path_data_observedNEvents, pDA=path_data_adjustments, pMCWNEg=path_MC_weightedNEvents_gluino, pMCWNEs=path_MC_weightedNEvents_squark, bla=binLabelAbbreviations[signalType], pSN=path_systematics_nominal, iFBCU=inputFolder_bkgCompositionUncertainties, sT=signalType, aOD=analysisOutputDirectory, n=nJetsBin)
    ratioRangeToPlot = None
    if ((bkgType is None) or (bkgType == "pre")):
        ratioRangeToPlot = ratioRangesToPlot["pre"][signalType]
    elif (bkgType == "post"):
        ratioRangeToPlot = ratioRangesToPlot["post"][signalType]
    else:
        removeLock()
        sys.exit("ERROR: unrecognized bkgType: {t}".format(t=bkgType))
    command_signalSTDistributions_dataAndSignal += " --ratioMin {rmin:.3f} --ratioMax {rmax:.3f}".format(rmin=ratioRangeToPlot["min"], rmax=ratioRangeToPlot["max"])
    if run_unblinded:
        command_signalSTDistributions_dataAndSignal += " --plotObservedData --path_fitDiagnostics {p} --bkgType {t} --outputFilePrefix {oFP}".format(p="{eP}/{aEOD}/fitDiagnostics.root".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory), t=bkgType, oFP="STDistributions_{t}Fit_{sT}".format(t=bkgType, sT=signalType))
    else:
        command_signalSTDistributions_dataAndSignal += " --outputFilePrefix {oFP}".format(oFP="STDistributions_blinded_{sT}".format(sT=signalType))
    # stealthEnv.execute_in_env(commandToRun=command_signalSTDistributions_dataAndSignal, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)
    commands_ancillary_plots_signal.append(command_signalSTDistributions_dataAndSignal)
    return commands_ancillary_plots_signal

def get_commands_plot_limits(eventProgenitor):
    commands_plot_limits = []
    selectionsString = "signal"
    if not(inputArguments.noLooseSignal): selectionsString += ",signal_loose"
    command_plotLimits = "./plotLimits.py --crossSectionsFile {cSFP} --MCTemplatePath {MTP} --eventProgenitor {eP2} --combineResultsDirectory {eP}/{cREOSOD} --combineOutputPrefix {eP2} --outputDirectory_rawOutput {aOD}/limits --outputDirectory_plots {aOD}/publicationPlots --outputSuffix {eP2} --minNeutralinoMass {mNM} --eventProgenitorMassOffset {ePMO} --minMassDifference {mMD} --selectionsList {sS} --signalContaminationSource_signal {sCSS} --signalContaminationSource_signal_loose {sCSSL} --signalContaminationMonitor_source_folder_eos {sCMSF}".format(eP=stealthEnv.EOSPrefix, MTP=MCTemplatesForProgenitor[eventProgenitor], cREOSOD=combineResultsEOSOutputDirectory, eP2=eventProgenitor, cSFP=crossSectionsForProgenitor[eventProgenitor], aOD=analysisOutputDirectory, mNM=minNeutralinoMassToPlot[eventProgenitor], ePMO=eventProgenitorMassOffsets[eventProgenitor], mMD=minMassDifferences[eventProgenitor], sS=selectionsString, sCSS="{eP}/{aEOD}/MC_stealth_{eP2}_all_signal_MCUncertainties_savedObjects.root".format(eP=stealthEnv.EOSPrefix, eP2=eventProgenitor, aEOD=analysisEOSOutputDirectory), sCSSL="{eP}/{aEOD}/MC_stealth_{eP2}_all_signal_loose_MCUncertainties_savedObjects.root".format(eP=stealthEnv.EOSPrefix, eP2=eventProgenitor, aEOD=analysisEOSOutputDirectory), sCMSF="{aEOD}/signalContaminationMonitor".format(aEOD=analysisEOSOutputDirectory))
    commands_plot_limits.append(command_plotLimits)
    if inputArguments.runUnblinded:
        command_plotLimits += " --plotObserved"
        commands_plot_limits.append(command_plotLimits)
    return commands_plot_limits

for step in runSequence:
    if (step == "data"):
        shellCommands_control = get_commands_data_chain(inputFilesList="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_data_{yP}_control.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, yP=yearPattern), outputPrefix="control", analyzeSignalBins=True)
        if (inputArguments.isDryRun): print("Not spawning due to dry run flag: {sC_c}".format(sC_c=shellCommands_control))
        else: multiProcessLauncher.spawn(shellCommands=shellCommands_control, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step_data_control.log", printDebug=True)
        for signalType in list_signalTypes:
            shellCommands_signal = get_commands_data_chain(inputFilesList="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_data_{yP}_{sT}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, yP=yearPattern, sT=signalType), outputPrefix="{sT}".format(sT=signalType), analyzeSignalBins=inputArguments.runUnblinded)
            if (inputArguments.isDryRun): print("Not spawning due to dry run flag: {sC_s}".format(sC_s=shellCommands_signal))
            else: multiProcessLauncher.spawn(shellCommands=shellCommands_signal, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step_data_{sT}.log".format(sT=signalType), printDebug=True)
        if not(inputArguments.isDryRun): multiProcessLauncher.monitorToCompletion()
    elif (step == "BKGMC"):
        # command_update = ("cd getPUWeights && make && cd ../fitScripts && make && cd ..")
        # stealthEnv.execute_in_env(commandToRun=command_update, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)

        # Initialize some paths
        inputBKGMCPaths = {}
        inputBKGMCPaths_singlephoton = {}
        # PUReweightingFileNames = {}
        for yearString in ["16", "17", "18"]:
            inputBKGMCPaths[yearString] = {}
            inputBKGMCPaths_singlephoton[yearString] = {}
            # PUReweightingFileNames[yearString] = {}
            for signalType in (list_signalTypes + ["control"]):
                background_names = ["DiPhotonJets", "GJetHT", "HighHTQCD"]
                background_name_contains_y2 = {
                    "DiPhotonJets": False,
                    "GJetHT": True,
                    "HighHTQCD": True
                }
                selection_string_singlephoton = None
                if (signalType == "signal"):
                    selection_string_singlephoton = "singlemedium"
                elif (signalType == "signal_loose"):
                    selection_string_singlephoton = "singleloose"
                elif (signalType == "control"):
                    selection_string_singlephoton = "singlefake"

                inputBKGMCPaths[yearString][signalType] = {}
                inputBKGMCPaths_singlephoton[yearString][signalType] = {}
                # PUReweightingFileNames[yearString][signalType] = {}
                for background_name in background_names:
                    year_int = 2000 + int(0.5 + float(yearString))
                    yearStringInPath = ""
                    if background_name_contains_y2[background_name]: yearStringInPath = yearString
                    inputBKGMCPaths[yearString][signalType][background_name] = "{eP}/store/user/lpcsusystealth/selections/combined_DoublePhoton{s_s}/merged_selection_MC_{b}{ys}_{y}_{sT}.root".format(eP=stealthEnv.EOSPrefix, s_s=selection_suffix, b=background_name, ys=yearStringInPath, y=year_int, sT=signalType)
                    inputBKGMCPaths_singlephoton[yearString][signalType][background_name] = "{eP}/store/user/lpcsusystealth/selections/combined_DoublePhoton{s_s}/merged_selection_MC_{b}{ys}_singlephoton_{y}_control_{s_s_s}.root".format(eP=stealthEnv.EOSPrefix, s_s=selection_suffix, b=background_name, ys=yearStringInPath, y=year_int, s_s_s=selection_string_singlephoton)
                    # PUReweightingFileNames[yearString][signalType][background_name] = "PUWeights_bkg_{b}_{y}_{s}.root".format(b=background_name, y=yearString, s=signalType)
                    # command_getPUReweightingHistograms = "./getPUWeights/bin/makePUWeights"
                    # command_getPUReweightingHistograms += " inputDataPath={i}".format(i=("getPUWeights/data/dataPU_20{y}.root".format(y=yearString)))
                    # command_getPUReweightingHistograms += " inputMCPaths={i}".format(i=inputBKGMCPaths_singlephoton[yearString][signalType][background_name])
                    # command_getPUReweightingHistograms += " outputFolder={eP}/{aEOD}".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory)
                    # command_getPUReweightingHistograms += " outputFileName={o}".format(o=PUReweightingFileNames[yearString][signalType][background_name])
                    # command_getPUReweightingHistograms += " addMCXSecWeight=true"
                    # if (inputArguments.isDryRun):
                    #     print("Not spawning due to dry run flag: {c_gPURH}".format(c_gPURH=command_getPUReweightingHistograms))
                    # else:
                    #     multiProcessLauncher.spawn(shellCommands=[command_getPUReweightingHistograms], optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step_BKGMC_PUReweighting_{b}_{y}_{s}.log".format(b=background_name, y=yearString, s=signalType), printDebug=True)
            # if not(inputArguments.isDryRun): multiProcessLauncher.monitorToCompletion()
        
        # Step 0: Create dat file containing MC norms
        command_getMCNorms = "./getMCNorms/py_scripts/get_norms.py"
        if (inputArguments.optionalIdentifier != ""): command_getMCNorms += " --optionalIdentifier {o}".format(o=inputArguments.optionalIdentifier)
        if (inputArguments.runUnblinded): command_getMCNorms += " --runUnblinded"
        stealthEnv.execute_in_env(commandToRun=command_getMCNorms, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)
        norm_values_cfg = tmGeneralUtils.getConfigurationFromFile("{aOD}/MCNorms/norm_values_nominal.dat".format(aOD=analysisOutputDirectory))
        nominal_norm_value_strings = {}
        nominal_norm_value_strings_singlephoton = {}
        for background_name in ["DiPhotonJets", "GJetHT", "HighHTQCD"]:
            nominal_norm_value_strings[background_name] = ""
            nominal_norm_value_strings_singlephoton[background_name] = ""
            for nJetsBin in range(2, 7):
                nominal_norm_value_strings[background_name] += ("1.0#")
                nominal_norm_value_strings_singlephoton[background_name] += (str(norm_values_cfg["norm_values_{p}_{n}JetsBin".format(p=background_name, n=nJetsBin)]) + "#")
            nominal_norm_value_strings[background_name] = (nominal_norm_value_strings[background_name])[:-1] # To remove the last #
            nominal_norm_value_strings_singlephoton[background_name] = (nominal_norm_value_strings_singlephoton[background_name])[:-1] # To remove the last #
        sourceData_BKGMC_dict = {}
        sourceData_BKGMC_singlephoton_dict = {}
        for signalType in (list_signalTypes + ["control"]):
            selection_string_singlephoton = None
            if (signalType == "signal"): selection_string_singlephoton = "singlemedium"
            elif (signalType == "signal_loose"): selection_string_singlephoton = "singleloose"
            elif (signalType == "control"): selection_string_singlephoton = "singlefake"
            sourceData_BKGMC_dict[signalType] = {}
            sourceData_BKGMC_singlephoton_dict[signalType] = {}
            for yearString in ["16", "17", "18"]:
                sourceData_BKGMC_dict[signalType]["bkgDiph{ys}".format(ys=yearString)] = inputBKGMCPaths[yearString][signalType]["DiPhotonJets"]
                # sourceData_BKGMC_dict[signalType]["PUDiph{ys}".format(ys=yearString)] = "{eP}/{aEOD}/".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory) + PUReweightingFileNames[yearString][signalType]["DiPhotonJets"]
                sourceData_BKGMC_dict[signalType]["bkgGJet{ys}".format(ys=yearString)] = inputBKGMCPaths[yearString][signalType]["GJetHT"]
                # sourceData_BKGMC_dict[signalType]["PUGJet{ys}".format(ys=yearString)] = "{eP}/{aEOD}/".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory) + PUReweightingFileNames[yearString][signalType]["GJetHT"]
                sourceData_BKGMC_dict[signalType]["bkgQCD{ys}".format(ys=yearString)] = inputBKGMCPaths[yearString][signalType]["HighHTQCD"]
                # sourceData_BKGMC_dict[signalType]["PUQCD{ys}".format(ys=yearString)] = "{eP}/{aEOD}/".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory) + PUReweightingFileNames[yearString][signalType]["HighHTQCD"]
                sourceData_BKGMC_singlephoton_dict[signalType]["bkgDiph{ys}".format(ys=yearString)] = inputBKGMCPaths_singlephoton[yearString][signalType]["DiPhotonJets"]
                # sourceData_BKGMC_singlephoton_dict[signalType]["PUDiph{ys}".format(ys=yearString)] = "{eP}/{aEOD}/".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory) + PUReweightingFileNames[yearString][signalType]["DiPhotonJets"]
                sourceData_BKGMC_singlephoton_dict[signalType]["bkgGJet{ys}".format(ys=yearString)] = inputBKGMCPaths_singlephoton[yearString][signalType]["GJetHT"]
                # sourceData_BKGMC_singlephoton_dict[signalType]["PUGJet{ys}".format(ys=yearString)] = "{eP}/{aEOD}/".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory) + PUReweightingFileNames[yearString][signalType]["GJetHT"]
                sourceData_BKGMC_singlephoton_dict[signalType]["bkgQCD{ys}".format(ys=yearString)] = inputBKGMCPaths_singlephoton[yearString][signalType]["HighHTQCD"]
                # sourceData_BKGMC_singlephoton_dict[signalType]["PUQCD{ys}".format(ys=yearString)] = "{eP}/{aEOD}/".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory) + PUReweightingFileNames[yearString][signalType]["HighHTQCD"]
            sourceData_BKGMC_dict[signalType]["wgtDiph"] = nominal_norm_value_strings["DiPhotonJets"]
            sourceData_BKGMC_dict[signalType]["wgtGJet"] = nominal_norm_value_strings["GJetHT"]
            sourceData_BKGMC_dict[signalType]["wgtQCD"] = nominal_norm_value_strings["HighHTQCD"]
            sourceData_BKGMC_singlephoton_dict[signalType]["wgtDiph"] = nominal_norm_value_strings_singlephoton["DiPhotonJets"]
            sourceData_BKGMC_singlephoton_dict[signalType]["wgtGJet"] = nominal_norm_value_strings_singlephoton["GJetHT"]
            sourceData_BKGMC_singlephoton_dict[signalType]["wgtQCD"] = nominal_norm_value_strings_singlephoton["HighHTQCD"]

        # Step 1: Run over full background MC double photon selections
        # for signalType in (list_signalTypes + ["control"]):
        for signalType in (list_signalTypes):
            compare_data_to_MC_prediction = (signalType == "control")
            sourceData_data = None
            if compare_data_to_MC_prediction:
                sourceData_data = "{eP}/store/user/lpcsusystealth/selections/combined_DoublePhoton{s_s}/merged_selection_data_2016_{sT}.root,{eP}/store/user/lpcsusystealth/selections/combined_DoublePhoton{s_s}/merged_selection_data_2017_{sT}.root,{eP}/store/user/lpcsusystealth/selections/combined_DoublePhoton{s_s}/merged_selection_data_2018_{sT}.root".format(eP=stealthEnv.EOSPrefix, s_s=selection_suffix, sT=signalType)
            rho_nominal = read_rho_nominal_from_file(rhoNominalFilePath="{aOD}/dataSystematics/{sT}_rhoNominal.dat".format(aOD=analysisOutputDirectory, sT=signalType))
            sourceData_BKGMC = "{bkgDiph16}!true!{wgtDiph},{bkgGJet16}!true!{wgtGJet},{bkgQCD16}!true!{wgtQCD},{bkgDiph17}!true!{wgtDiph},{bkgGJet17}!true!{wgtGJet},{bkgQCD17}!true!{wgtQCD},{bkgDiph18}!true!{wgtDiph},{bkgGJet18}!true!{wgtGJet},{bkgQCD18}!true!{wgtQCD}".format(**(sourceData_BKGMC_dict[signalType]))
            adjustmentPlots_min = -0.5
            adjustmentPlots_max = 5.5
            if signalType in manualAdjustmentRatios["combined"]:
                adjustmentPlots_min = manualAdjustmentRatios["combined"][signalType][0]
                adjustmentPlots_max = manualAdjustmentRatios["combined"][signalType][1]
            shellCommands_BKGMC_doublephoton = get_commands_doublephoton_BKGMC_chain(sourceData_BKGMC=sourceData_BKGMC, readParametersExplicitlyFromSource=None, adjustmentPlots_min=adjustmentPlots_min, adjustmentPlots_max=adjustmentPlots_max, sourceData_data=sourceData_data, identifier="MC_Bkg", outputFolder="{aOD}/fits_doublephoton".format(aOD=analysisOutputDirectory), selectionString=signalType, rhoNominal=rho_nominal, disableStrictChecks=False)
            if (inputArguments.isDryRun): print("Not spawning due to dry run flag: {sC_BKGMC_d}".format(sC_BKGMC_d=shellCommands_BKGMC_doublephoton))
            else: multiProcessLauncher.spawn(shellCommands=shellCommands_BKGMC_doublephoton, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step_BKGMC_doublephoton_{sT}.log".format(sT=signalType), printDebug=True)
        if not(inputArguments.isDryRun): multiProcessLauncher.monitorToCompletion()

        # Step 2: Run over the six different plausible background compositions
        # for signalType in (list_signalTypes + ["control"]):
        for signalType in (list_signalTypes):
            compare_data_to_MC_prediction = (signalType == "control")
            rho_nominal = read_rho_nominal_from_file(rhoNominalFilePath="{aOD}/dataSystematics/{sT}_rhoNominal.dat".format(aOD=analysisOutputDirectory, sT=signalType))
            outputFolder="{aOD}/fits_doublephoton".format(aOD=analysisOutputDirectory)

            executionStatusFilesToMonitor = []
            for bkg_to_modulate in ["Diph", "GJet", "QCD"]:
                for modulation in ["up", "down"]:
                    source_data_string = ""
                    for year_string_to_add in ["16", "17", "18"]:
                        for bkg_to_add in ["Diph", "GJet", "QCD"]:
                            weight_string = None
                            if (bkg_to_add == bkg_to_modulate):
                                if (modulation == "up"):
                                    weight_string = "2.0"
                                elif (modulation == "down"):
                                    weight_string = "0.5"
                            source_data_string += "{bkg" + bkg_to_add + year_string_to_add + "}!true!{wgt" + bkg_to_add + "}"
                            if not(weight_string is None): source_data_string += "!" + weight_string
                            source_data_string += ","
                    source_data_string = source_data_string[:-1] # To remove the last comma
                    adjustmentPlots_min = -0.5
                    adjustmentPlots_max = 5.5
                    shellCommands_modulated_bkg = get_commands_doublephoton_BKGMC_chain(sourceData_BKGMC=(source_data_string.format(**(sourceData_BKGMC_dict[signalType]))), readParametersExplicitlyFromSource="{oF}/binned_fitParameters_all_MC_Bkg_{s}.dat".format(oF=outputFolder, s=signalType), adjustmentPlots_min=adjustmentPlots_min, adjustmentPlots_max=adjustmentPlots_max, sourceData_data=None, identifier="MC_{b}_shift_{ud}".format(b=bkg_to_modulate, ud=modulation), outputFolder=outputFolder, selectionString=signalType, rhoNominal=rho_nominal, disableStrictChecks=False)
                    if (inputArguments.isDryRun):
                        print("Not running due to dry run flag: {c}".format(c=shellCommands_modulated_bkg))
                    else:
                        multiProcessLauncher.spawn(shellCommands=shellCommands_modulated_bkg, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step_BKGMC_doublephoton_{b}_shift_{ud}_{sT}.log".format(b=bkg_to_modulate, ud=modulation, sT=signalType), printDebug=True)
                        executionStatusFilesToMonitor.append("{oF}/execution_status_all_{i}_{s}.txt".format(oF=outputFolder, i="MC_{b}_shift_{ud}".format(b=bkg_to_modulate, ud=modulation) ,s=signalType))
                        # Command above finishes running and then crashes for no apparent reason
                        # for command in shellCommands_modulated_bkg:
                        #     subprocess.check_call(command, executable="/bin/bash", shell=True)
            if not(inputArguments.isDryRun): multiProcessLauncher.monitorToCompletion(killAllOnOneFailure=False)
            check_execution_statuses_manually(executionStatusFilesToMonitor)

        # Step 3: Do single photon selections
        executionStatusFilesToMonitor = []
        for signalType in (list_signalTypes):
            selection_string_singlephoton = None
            if (signalType == "signal"): selection_string_singlephoton = "singlemedium"
            elif (signalType == "signal_loose"): selection_string_singlephoton = "singleloose"
            elif (signalType == "control"): selection_string_singlephoton = "singlefake"
            # rho_nominal = read_rho_nominal_from_file(rhoNominalFilePath="{aOD}/dataSystematics/{sT}_rhoNominal.dat".format(aOD=analysisOutputDirectory, sT=signalType))
            sourceData_BKGMC_singlephoton="{bkgDiph16}!true!{wgtDiph},{bkgGJet16}!true!{wgtGJet},{bkgQCD16}!true!{wgtQCD},{bkgDiph17}!true!{wgtDiph},{bkgGJet17}!true!{wgtGJet},{bkgQCD17}!true!{wgtQCD},{bkgDiph18}!true!{wgtDiph},{bkgGJet18}!true!{wgtGJet},{bkgQCD18}!true!{wgtQCD}".format(**(sourceData_BKGMC_singlephoton_dict[signalType]))
            sourceData_data_singlephoton="{eP}/store/user/lpcsusystealth/selections/combined_DoublePhoton{s_s}/merged_selection_data_singlephoton_2016_control_{s_s_s}.root,{eP}/store/user/lpcsusystealth/selections/combined_DoublePhoton{s_s}/merged_selection_data_singlephoton_2017_control_{s_s_s}.root,{eP}/store/user/lpcsusystealth/selections/combined_DoublePhoton{s_s}/merged_selection_data_singlephoton_2018_control_{s_s_s}.root".format(eP=stealthEnv.EOSPrefix, s_s=selection_suffix, s_s_s=selection_string_singlephoton)
            adjustmentPlots_min = -0.5
            adjustmentPlots_max = 5.5
            shellCommands_BKGMC_singlephoton = get_commands_singlephoton_BKGMC_chain(sourceData_BKGMC=sourceData_BKGMC_singlephoton, sourceData_data=sourceData_data_singlephoton, adjustmentPlots_min=adjustmentPlots_min, adjustmentPlots_max=adjustmentPlots_max, outputFolder="{aOD}/fits_singlephoton".format(aOD=analysisOutputDirectory), selectionString=selection_string_singlephoton, yearString="all", rhoNominal=1.2, disableStrictChecks=False)
            if (inputArguments.isDryRun):
                print("Not running due to dry run flag: {sC_BKGMC_s}".format(sC_BKGMC_s=shellCommands_BKGMC_singlephoton))
            else:
                multiProcessLauncher.spawn(shellCommands=shellCommands_BKGMC_singlephoton, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step_BKGMC_singlephoton_{sT}.log".format(sT=signalType), printDebug=True)
                executionStatusFilesToMonitor.append("{aOD}/fits_singlephoton/execution_status_all_MC_Bkg_{s_s_s}.txt".format(aOD=analysisOutputDirectory, s_s_s=selection_string_singlephoton))
                executionStatusFilesToMonitor.append("{aOD}/fits_singlephoton/execution_status_all_data_{s_s_s}.txt".format(aOD=analysisOutputDirectory, s_s_s=selection_string_singlephoton))
                # for command in shellCommands_BKGMC_singlephoton:
                #     subprocess.check_call(command, executable="/bin/bash", shell=True)
        if not(inputArguments.isDryRun): multiProcessLauncher.monitorToCompletion(killAllOnOneFailure=False)
        check_execution_statuses_manually(executionStatusFilesToMonitor)

        # Step 4: Single photon selections over six plausible background compositions
        for signalType in (list_signalTypes):
            selection_string_singlephoton = None
            if (signalType == "signal"): selection_string_singlephoton = "singlemedium"
            elif (signalType == "signal_loose"): selection_string_singlephoton = "singleloose"
            elif (signalType == "control"): selection_string_singlephoton = "singlefake"
            # rho_nominal = read_rho_nominal_from_file(rhoNominalFilePath="{aOD}/dataSystematics/{sT}_rhoNominal.dat".format(aOD=analysisOutputDirectory, sT=signalType))

            executionStatusFilesToMonitor = []
            for bkg_to_modulate in ["Diph", "GJet", "QCD"]:
                for modulation in ["up", "down"]:
                    sourceData_BKGMC_singlephoton_modulated = ""
                    for year_string_to_add in ["16", "17", "18"]:
                        for bkg_to_add in ["Diph", "GJet", "QCD"]:
                            weight_string = None
                            if (bkg_to_add == bkg_to_modulate):
                                if (modulation == "up"):
                                    weight_string = "2.0"
                                elif (modulation == "down"):
                                    weight_string = "0.5"
                            sourceData_BKGMC_singlephoton_modulated += "{bkg" + bkg_to_add + year_string_to_add + "}!true!{wgt" + bkg_to_add + "}"
                            if not(weight_string is None): sourceData_BKGMC_singlephoton_modulated += "!" + weight_string
                            sourceData_BKGMC_singlephoton_modulated += ","
                    sourceData_BKGMC_singlephoton_modulated = sourceData_BKGMC_singlephoton_modulated[:-1] # To remove the last comma
                    adjustmentPlots_min = -0.5
                    adjustmentPlots_max = 5.5
                    shellCommands_BKGMC_modulated_singlephoton = get_commands_singlephoton_modulated_BKGMC_chain(sourceData_BKGMC=sourceData_BKGMC_singlephoton_modulated.format(**(sourceData_BKGMC_singlephoton_dict[signalType])), readParametersExplicitlyFromSource="{oF}/binned_fitParameters_all_MC_Bkg_{s_s_s}.dat".format(oF="{aOD}/fits_singlephoton".format(aOD=analysisOutputDirectory), s_s_s=selection_string_singlephoton), identifier="MC_{b}_shift_{ud}".format(b=bkg_to_modulate, ud=modulation), adjustmentPlots_min=adjustmentPlots_min, adjustmentPlots_max=adjustmentPlots_max, outputFolder="{aOD}/fits_singlephoton".format(aOD=analysisOutputDirectory), selectionString=selection_string_singlephoton, yearString="all", rhoNominal=1.2, disableStrictChecks=False)
                    if (inputArguments.isDryRun):
                        print("Not spawning due to dry run flag: {sC_BKGMC_s}".format(sC_BKGMC_s=shellCommands_BKGMC_modulated_singlephoton))
                    else:
                        multiProcessLauncher.spawn(shellCommands=shellCommands_BKGMC_modulated_singlephoton, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step_BKGMC_singlephoton_{b}_shift_{ud}_{sT}.log".format(b=bkg_to_modulate, ud=modulation, sT=signalType), printDebug=True)
                        executionStatusFilesToMonitor.append("{aOD}/fits_singlephoton/execution_status_all_{i}_{s_s_s}.txt".format(aOD=analysisOutputDirectory, i="MC_{b}_shift_{ud}".format(b=bkg_to_modulate, ud=modulation), s_s_s=selection_string_singlephoton))
                        # for command in shellCommands_BKGMC_modulated_singlephoton:
                        #     subprocess.check_call(command, executable="/bin/bash", shell=True)
            if not(inputArguments.isDryRun): multiProcessLauncher.monitorToCompletion(killAllOnOneFailure=False)
            check_execution_statuses_manually(executionStatusFilesToMonitor)

        # Step 5: Find adjustments for three backgrounds separately
        for signalType in (list_signalTypes):
            # compare_data_to_MC_prediction = (signalType == "control")
            compare_data_to_MC_prediction = False
            rho_nominal = read_rho_nominal_from_file(rhoNominalFilePath="{aOD}/dataSystematics/{sT}_rhoNominal.dat".format(aOD=analysisOutputDirectory, sT=signalType))
            outputFolder="{aOD}/fits_doublephoton_decoupled".format(aOD=analysisOutputDirectory)

            executionStatusFilesToMonitor = []
            for bkg in ["Diph", "GJet", "QCD"]:
                source_data_string = ""
                for year_string_to_add in ["16", "17", "18"]:
                    source_data_string += "{bkg" + bkg + year_string_to_add + "}!true!{wgt" + bkg_to_add + "},"
                source_data_string = source_data_string[:-1] # To remove the last comma
                adjustmentPlots_min = -0.5
                adjustmentPlots_max = 5.5
                shellCommands_decoupled_bkg = get_commands_doublephoton_BKGMC_chain(sourceData_BKGMC=(source_data_string.format(**(sourceData_BKGMC_dict[signalType]))), readParametersExplicitlyFromSource=None, adjustmentPlots_min=adjustmentPlots_min, adjustmentPlots_max=adjustmentPlots_max, sourceData_data=None, identifier="MC_{b}".format(b=bkg), outputFolder=outputFolder, selectionString=signalType, rhoNominal=rho_nominal, disableStrictChecks=True)
                if (inputArguments.isDryRun):
                    print("Not spawning due to dry run flag: {c}".format(c=shellCommands_decoupled_bkg))
                else:
                    multiProcessLauncher.spawn(shellCommands=shellCommands_decoupled_bkg, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step_BKGMC_doublephoton_decoupled_{b}_{sT}.log".format(b=bkg, sT=signalType), printDebug=True)
                    executionStatusFilesToMonitor.append("{oF}/execution_status_all_{i}_{s}.txt".format(oF=outputFolder, i="MC_{b}".format(b=bkg), s=signalType))
                    # for command in shellCommands_decoupled_bkg:
                    #     subprocess.check_call(command, executable="/bin/bash", shell=True)
            if not(inputArguments.isDryRun): multiProcessLauncher.monitorToCompletion(killAllOnOneFailure=False)
            check_execution_statuses_manually(executionStatusFilesToMonitor)

    elif (step == "MC"):
        # command_update = ("cd getPUWeights && make && cd ../getMCSystematics && make && cd ..")
        # stealthEnv.execute_in_env(commandToRun=command_update, isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)
        for eventProgenitor in eventProgenitors:
            crossSectionsPath = crossSectionsForProgenitor[eventProgenitor]
            for signalType in (list_signalTypes + ["control"]):
                MCPathMain = ""
                # dataPUSourceMain = ""
                HLTEfficienciesPathMain = ""
                # PUWeightsOutputPathMain = ""
                lumiMain = ""
                MCPathsAux = []
                # dataPUSourcesAux = []
                # PUWeightsOutputPathsAux = []
                HLTEfficienciesPathsAux = []
                lumisAux = []
                if (inputArguments.year == "all"):
                    MCPathMain = "{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_{tD}_2017_{sT}.root".format(sT=signalType, eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, tD=tDesignationsForProgenitor[eventProgenitor])
                    # dataPUSourceMain = "getPUWeights/data/dataPU_2017.root"
                    HLTEfficienciesPathMain = "{eP}/{HES}/HLTEfficiencies_{sT}_2017.root".format(eP=stealthEnv.EOSPrefix, HES=HLTEfficienciesSource, sT=signalType)
                    # PUWeightsOutputPathMain = "PUWeights_2017_{eP}_{sT}.root".format(eP=eventProgenitor, sT=signalType)
                    lumiMain = integrated_lumi_strings["2017"]
                    MCPathsAux = ["{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_{tD}_2016_{sT}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, sT=signalType, tD=tDesignationsForProgenitor[eventProgenitor]), "{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_{tD}_2018_{sT}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, sT=signalType, tD=tDesignationsForProgenitor[eventProgenitor])]
                    # dataPUSourcesAux = ["getPUWeights/data/dataPU_2016.root", "getPUWeights/data/dataPU_2018.root"]
                    # PUWeightsOutputPathsAux = ["PUWeights_2016_{eP}_{sT}.root".format(eP=eventProgenitor, sT=signalType), "PUWeights_2018_{eP}_{sT}.root".format(eP=eventProgenitor, sT=signalType)]
                    HLTEfficienciesPathsAux = ["{eP}/{HES}/HLTEfficiencies_{sT}_2016.root".format(eP=stealthEnv.EOSPrefix, HES=HLTEfficienciesSource, sT=signalType), "{eP}/{HES}/HLTEfficiencies_{sT}_2018.root".format(eP=stealthEnv.EOSPrefix, HES=HLTEfficienciesSource, sT=signalType)]
                    lumisAux = [integrated_lumi_strings["2016"], integrated_lumi_strings["2018"]]
                else:
                    MCPathMain = "{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_{tD}_{y}_{sT}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, y=inputArguments.year, sT=signalType, tD=tDesignationsForProgenitor[eventProgenitor])
                    # dataPUSourceMain = "getPUWeights/data/dataPU_{y}.root".format(y=inputArguments.year)
                    HLTEfficienciesPathMain = "{eP}/{HES}/HLTEfficiencies_{sT}_{y}.root".format(sT=signalType, y=inputArguments.year)
                    # PUWeightsOutputPathMain = "PUWeights_{y}_{eP}_{sT}.root".format(y=inputArguments.year, eP=eventProgenitor, sT=signalType)
                    lumiMain = integrated_lumi_strings[inputArguments.year]
                get_signalContamination_outside_sidebands = False
                if (signalType == "control"): get_signalContamination_outside_sidebands = True
                shellCommands_MC = get_commands_MC_chain(eventProgenitor=eventProgenitor, dataPrefix=signalType, outputPrefix="MC_stealth_{eP}_{y}_{sT}".format(eP=eventProgenitor, y=inputArguments.year, sT=signalType), inputMCPathMain=MCPathMain, # inputDataPUSourceMain=dataPUSourceMain, PUWeightsOutputPathMain=PUWeightsOutputPathMain, 
                                                         inputHLTEfficienciesPathMain=HLTEfficienciesPathMain, integratedLuminosityMainString=lumiMain, inputMCPathsAux=MCPathsAux, # inputDataPUSourcesAux=dataPUSourcesAux, PUWeightsOutputPathsAux=PUWeightsOutputPathsAux, 
                                                         inputHLTEfficienciesPathsAux=HLTEfficienciesPathsAux, integratedLuminositiesAux=lumisAux, getSignalContaminationOutsideSidebands=get_signalContamination_outside_sidebands)
                if (inputArguments.isDryRun): print("Not spawning due to dry run flag: {sC_MC}".format(sC_MC=shellCommands_MC))
                else: multiProcessLauncher.spawn(shellCommands=shellCommands_MC, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step_MC_{eP}_{sT}.log".format(eP=eventProgenitor, sT=signalType), printDebug=True)
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
            paths_bkgCompositionSystematic = ""
            for bkg_to_modulate in ["Diph", "GJet", "QCD"]:
                for modulation in ["up", "down"]:
                    paths_bkgCompositionSystematic += "{p},".format(p="{aOD}/fits_doublephoton/ratio_adjustment_all_MC_{b}_shift_{s}_signal.dat".format(aOD=analysisOutputDirectory, b=bkg_to_modulate, s=modulation))
                    if not(inputArguments.noLooseSignal): paths_bkgCompositionSystematic += "{p},".format(p="{aOD}/fits_doublephoton/ratio_adjustment_all_MC_{b}_shift_{s}_signal_loose.dat".format(aOD=analysisOutputDirectory, b=bkg_to_modulate, s=modulation))
            paths_bkgCompositionSystematic = paths_bkgCompositionSystematic[:-1] # To remove trailing comma
            # run_combine_chain(eventProgenitor=eventProgenitor, path_dataSystematics_signal="{aOD}/dataSystematics/signal_dataSystematics.dat".format(aOD=analysisOutputDirectory), path_dataSystematics_signal_loose="{aOD}/dataSystematics/signal_loose_dataSystematics.dat".format(aOD=analysisOutputDirectory), path_dataSystematics_control="{aOD}/dataSystematics/control_dataSystematics.dat".format(aOD=analysisOutputDirectory), path_MCShapeAdjustment_signal=MCNominalAdjustmentFilePaths["signal"], path_MCShapeAdjustment_signal_loose=MCNominalAdjustmentFilePaths["signal_loose"], path_MCShapeAdjustment_control=MCNominalAdjustmentFilePaths["control"], paths_bkgCompositionSystematic=paths_bkgCompositionSystematic, path_dataObservedEventCounters_signal="{aOD}/dataSystematics/signal_observedEventCounters.dat".format(aOD=analysisOutputDirectory), path_dataObservedEventCounters_signal_loose="{aOD}/dataSystematics/signal_loose_observedEventCounters.dat".format(aOD=analysisOutputDirectory), path_dataObservedEventCounters_control="{aOD}/dataSystematics/control_observedEventCounters.dat".format(aOD=analysisOutputDirectory), path_dataExpectedEventCounters_signal="{aOD}/dataSystematics/signal_eventCounters.dat".format(aOD=analysisOutputDirectory), path_dataExpectedEventCounters_signal_loose="{aOD}/dataSystematics/signal_loose_eventCounters.dat".format(aOD=analysisOutputDirectory), path_dataExpectedEventCounters_control="{aOD}/dataSystematics/control_eventCounters.dat".format(aOD=analysisOutputDirectory), MCPrefix_signal="MC_stealth_{eP}_{y}_signal".format(eP=eventProgenitor, y=inputArguments.year), MCPrefix_signal_loose="MC_stealth_{eP}_{y}_signal_loose".format(eP=eventProgenitor, y=inputArguments.year), MCPrefix_control="MC_stealth_{eP}_{y}_control".format(eP=eventProgenitor, y=inputArguments.year), outputPrefix="{eP}".format(eP=eventProgenitor))
            run_combine_chain(eventProgenitor=eventProgenitor, path_dataSystematics_signal="{aOD}/dataSystematics/signal_dataSystematics.dat".format(aOD=analysisOutputDirectory), path_dataSystematics_signal_loose="{aOD}/dataSystematics/signal_loose_dataSystematics.dat".format(aOD=analysisOutputDirectory), path_MCShapeAdjustment_signal=MCNominalAdjustmentFilePaths["signal"], path_MCShapeAdjustment_signal_loose=MCNominalAdjustmentFilePaths["signal_loose"], paths_bkgCompositionSystematic=paths_bkgCompositionSystematic, path_dataObservedEventCounters_signal="{aOD}/dataSystematics/signal_observedEventCounters.dat".format(aOD=analysisOutputDirectory), path_dataObservedEventCounters_signal_loose="{aOD}/dataSystematics/signal_loose_observedEventCounters.dat".format(aOD=analysisOutputDirectory), path_dataExpectedEventCounters_signal="{aOD}/dataSystematics/signal_eventCounters.dat".format(aOD=analysisOutputDirectory), path_dataExpectedEventCounters_signal_loose="{aOD}/dataSystematics/signal_loose_eventCounters.dat".format(aOD=analysisOutputDirectory), MCPrefix_signal="MC_stealth_{eP}_{y}_signal".format(eP=eventProgenitor, y=inputArguments.year), MCPrefix_signal_loose="MC_stealth_{eP}_{y}_signal_loose".format(eP=eventProgenitor, y=inputArguments.year), outputPrefix="{eP}".format(eP=eventProgenitor))
    elif (step == "ancillaryPlots"):
        # produce_STComparisons(dataPath="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_data_{yP}_control.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, yP=yearPattern), outputFilePrefix_STComparisons="control_STComparisons", analyzeSignalBins=True, useWeights=False)
        for signalType in list_signalTypes:
            produce_STComparisons(dataPath="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_data_{yP}_{sT}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, yP=yearPattern, sT=signalType), outputFilePrefix_STComparisons="{sT}_STComparisons".format(sT=signalType), analyzeSignalBins=inputArguments.runUnblinded, useWeights=False)
        for nJetsBin in range(4, 7):
            # commands_ancillary_plots_control = get_commands_ancillary_plots_control(path_data_expectedNEvents="{aOD}/dataSystematics/control_eventCounters.dat".format(aOD=analysisOutputDirectory), path_data_observedNEvents="{aOD}/dataSystematics/control_observedEventCounters.dat".format(aOD=analysisOutputDirectory), path_data_adjustments=MCNominalAdjustmentFilePaths["control"], path_MC_weightedNEvents_gluino="{aOD}/MCEventHistograms/MC_stealth_gluino_{y}_control_savedObjects.root".format(aOD=analysisOutputDirectory, y=inputArguments.year), path_MC_weightedNEvents_squark="{aOD}/MCEventHistograms/MC_stealth_squark_{y}_control_savedObjects.root".format(aOD=analysisOutputDirectory, y=inputArguments.year), path_systematics_nominal="{aOD}/dataSystematics/control_dataSystematics.dat".format(aOD=analysisOutputDirectory), path_systematics_dataMCDiscrepancy="{aOD}/fits_singlephoton/{DMRAF}".format(aOD=analysisOutputDirectory, DMRAF=DataMCRatioAdjustmentsFilePaths["control"]), nJetsBin=nJetsBin)
            # if inputArguments.isDryRun:
            #     print("Not spawning due to dry run flag: {c}".format(c=commands_ancillary_plots_control))
            # else:
            #     multiProcessLauncher.spawn(shellCommands=commands_ancillary_plots_control, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step_ancillary_control_{n}Jets.log".format(n=nJetsBin), printDebug=True)
            for signalType in list_signalTypes:
                commands_ancillary_plots_signal = get_commands_ancillary_plots_signal(signalType=signalType, path_data_expectedNEvents="{aOD}/dataSystematics/{sT}_eventCounters.dat".format(aOD=analysisOutputDirectory, sT=signalType), path_data_observedNEvents="{aOD}/dataSystematics/{sT}_observedEventCounters.dat".format(aOD=analysisOutputDirectory, sT=signalType), path_data_adjustments=MCNominalAdjustmentFilePaths[signalType], path_MC_weightedNEvents_gluino="{aOD}/MCEventHistograms/MC_stealth_gluino_{y}_{sT}_savedObjects.root".format(aOD=analysisOutputDirectory, y=inputArguments.year, sT=signalType), path_MC_weightedNEvents_squark="{aOD}/MCEventHistograms/MC_stealth_squark_{y}_{sT}_savedObjects.root".format(aOD=analysisOutputDirectory, y=inputArguments.year, sT=signalType), path_systematics_nominal="{aOD}/dataSystematics/{sT}_dataSystematics.dat".format(aOD=analysisOutputDirectory, sT=signalType), inputFolder_bkgCompositionUncertainties="{aOD}/fits_doublephoton".format(aOD=analysisOutputDirectory), nJetsBin=nJetsBin, bkgType=None, run_unblinded=False)
                if inputArguments.isDryRun:
                    print("Not spawning due to dry run flag: {c}".format(c=commands_ancillary_plots_signal))
                else:
                    multiProcessLauncher.spawn(shellCommands=commands_ancillary_plots_signal, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step_ancillary_{sT}_{n}Jets.log".format(sT=signalType, n=nJetsBin), printDebug=True)
        if not(inputArguments.isDryRun): multiProcessLauncher.monitorToCompletion()
    elif (step == "observations"):
        if (inputArguments.runUnblinded):
            print("Producing fit diagnostics file...")
            path_data_card_template_folder = "{eP}/{aEOD}/dataCards/combinedFit".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory)
            path_data_card_template_file = "gluino_dataCard_eventProgenitorMassBin21_neutralinoMassBin73.txt".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory)
            output_folder_with_eos_prefix = "{eP}/{aEOD}".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory)
            stealthEnv.execute_in_env(commandToRun="{f}/{s} --outputFolder {o} --datacardTemplateParentFolderWithPrefix {tfolder} --datacardTemplateFileName {tfile} --identifier data".format(f=stealthEnv.stealthRoot, s="runStatisticsChecksUnblinded.py", o="{aOD}/statisticsChecks".format(aOD=analysisOutputDirectory) ,tfolder=path_data_card_template_folder, tfile=path_data_card_template_file), isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)
            transfer_file_to_EOS_area(sourceFile="{aOD}/statisticsChecks/data/fitDiagnostics.root".format(aOD=analysisOutputDirectory), targetDirectory="{eP}/{sER}/analysisEOSAreas/analysis{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier))
            for nJetsBin in range(4, 7):
                for signalType in list_signalTypes:
                    commands_ancillary_plots_signal_preFit = get_commands_ancillary_plots_signal(signalType=signalType, path_data_expectedNEvents="{aOD}/dataSystematics/{sT}_eventCounters.dat".format(aOD=analysisOutputDirectory, sT=signalType), path_data_observedNEvents="{aOD}/dataSystematics/{sT}_observedEventCounters.dat".format(aOD=analysisOutputDirectory, sT=signalType), path_data_adjustments=MCNominalAdjustmentFilePaths[signalType], path_MC_weightedNEvents_gluino="{aOD}/MCEventHistograms/MC_stealth_gluino_{y}_{sT}_savedObjects.root".format(aOD=analysisOutputDirectory, y=inputArguments.year, sT=signalType), path_MC_weightedNEvents_squark="{aOD}/MCEventHistograms/MC_stealth_squark_{y}_{sT}_savedObjects.root".format(aOD=analysisOutputDirectory, y=inputArguments.year, sT=signalType), path_systematics_nominal="{aOD}/dataSystematics/{sT}_dataSystematics.dat".format(aOD=analysisOutputDirectory, sT=signalType), inputFolder_bkgCompositionUncertainties="{aOD}/fits_doublephoton".format(aOD=analysisOutputDirectory), nJetsBin=nJetsBin, bkgType="pre", run_unblinded=True)
                    if inputArguments.isDryRun:
                        print("Not spawning due to dry run flag: {c}".format(c=commands_ancillary_plots_signal_preFit))
                    else:
                        multiProcessLauncher.spawn(shellCommands=commands_ancillary_plots_signal_preFit, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step_ancillary_with_observations_preFit_{sT}_{n}Jets.log".format(sT=signalType, n=nJetsBin), printDebug=True)
                    commands_ancillary_plots_signal_postFit = get_commands_ancillary_plots_signal(signalType=signalType, path_data_expectedNEvents="{aOD}/dataSystematics/{sT}_eventCounters.dat".format(aOD=analysisOutputDirectory, sT=signalType), path_data_observedNEvents="{aOD}/dataSystematics/{sT}_observedEventCounters.dat".format(aOD=analysisOutputDirectory, sT=signalType), path_data_adjustments=MCNominalAdjustmentFilePaths[signalType], path_MC_weightedNEvents_gluino="{aOD}/MCEventHistograms/MC_stealth_gluino_{y}_{sT}_savedObjects.root".format(aOD=analysisOutputDirectory, y=inputArguments.year, sT=signalType), path_MC_weightedNEvents_squark="{aOD}/MCEventHistograms/MC_stealth_squark_{y}_{sT}_savedObjects.root".format(aOD=analysisOutputDirectory, y=inputArguments.year, sT=signalType), path_systematics_nominal="{aOD}/dataSystematics/{sT}_dataSystematics.dat".format(aOD=analysisOutputDirectory, sT=signalType), inputFolder_bkgCompositionUncertainties="{aOD}/fits_doublephoton".format(aOD=analysisOutputDirectory), nJetsBin=nJetsBin, bkgType="post", run_unblinded=True)
                    if inputArguments.isDryRun:
                        print("Not spawning due to dry run flag: {c}".format(c=commands_ancillary_plots_signal_postFit))
                    else:
                        multiProcessLauncher.spawn(shellCommands=commands_ancillary_plots_signal_postFit, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step_ancillary_with_observations_postFit_{sT}_{n}Jets.log".format(sT=signalType, n=nJetsBin), printDebug=True)
            if not(inputArguments.isDryRun): multiProcessLauncher.monitorToCompletion()
        else:
            print("Doing nothing, runUnblinded flag is not set.")
    elif (step == "limits"):
        stealthEnv.execute_in_env(commandToRun="condor_q", isDryRun=inputArguments.isDryRun, functionToCallIfCommandExitsWithError=removeLock)
        for eventProgenitor in eventProgenitors:
            commands_plot_limits = get_commands_plot_limits(eventProgenitor)
            if inputArguments.isDryRun:
                print("Not spawning due to dry run flag: {c}".format(c=commands_plot_limits))
            else:
                multiProcessLauncher.spawn(shellCommands=commands_plot_limits, optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot), logFileName="step_limits_{eP}.log".format(eP=eventProgenitor), printDebug=True)
        if not(inputArguments.isDryRun): multiProcessLauncher.monitorToCompletion()
    elif (step == "statisticsChecks"):
        if (inputArguments.runUnblinded):
            datacards_to_check_dict = {
                "gluinoA": "gluino_dataCard_eventProgenitorMassBin23_neutralinoMassBin73.txt", # Gluino mass 2100 GeV, neutralino mass 1000 GeV
                "gluinoB": "gluino_dataCard_eventProgenitorMassBin20_neutralinoMassBin141.txt", # Gluino mass 1950 GeV, neutralino mass 1850 GeV
                "gluinoC": "gluino_dataCard_eventProgenitorMassBin18_neutralinoMassBin9.txt", # Gluino mass 1850 GeV, neutralino mass 200 GeV
                "squarkA": "squark_dataCard_eventProgenitorMassBin21_neutralinoMassBin57.txt", # Squark mass 1850 GeV, neutralino mass 800 GeV
                "squarkB": "squark_dataCard_eventProgenitorMassBin19_neutralinoMassBin125.txt", # Squark mass 1750 GeV, neutralino mass 1650 GeV
                "squarkC": "squark_dataCard_eventProgenitorMassBin14_neutralinoMassBin9.txt" # Squark mass 1500 GeV, neutralino mass 200 GeV
            }
            for datacard_id in datacards_to_check_dict:
                commands_run_statistics_checks = []
                commands_run_statistics_checks.append("./runStatisticsChecks.py --outputFolder {oF} --datacardTemplateParentFolderWithPrefix {dTPFWP} --datacardTemplateFileName {dTFN} --identifier {ident}".format(oF="{aOD}/statisticsChecks".format(aOD=analysisOutputDirectory), dTPFWP="{eP}/{aEOD}/dataCards/combinedFit".format(eP=stealthEnv.EOSPrefix, aEOD=analysisEOSOutputDirectory), dTFN=datacards_to_check_dict[datacard_id], ident=datacard_id))
                if inputArguments.isDryRun:
                    print("Not spawning due to dry run flag: {c}".format(c=commands_run_statistics_checks))
                else:
                    multiProcessLauncher.spawn(shellCommands=commands_run_statistics_checks,
                                               optionalEnvSetup="cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot),
                                               logFileName="step_statisticsChecks_{ident}.log".format(ident=datacard_id), printDebug=True)
            if not(inputArguments.isDryRun): multiProcessLauncher.monitorToCompletion()
        else:
            print("Doing nothing, runUnblinded flag is not set.")
    else:
        removeLock()
        sys.exit("ERROR: Unrecognized step: {s}".format(s=step))

removeLock()
