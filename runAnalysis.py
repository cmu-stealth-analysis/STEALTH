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
inputArgumentsParser.add_argument('--year', default="all", help="Year of data-taking. Can be \"2016\", \"2017\", \"2018\", or (default) \"all\".", type=str)
inputArgumentsParser.add_argument('--runUnblinded', action='store_true', help="If this flag is set, then the signal region data is unblinded.")
inputArgumentsParser.add_argument('--addLooseSignal', action='store_true', help="Add loose photons in a different signal bin. Run on a single signal type. By default data cards are created with only medium photons.")
inputArgumentsParser.add_argument('--isDryRun', action='store_true', help="Only print the commands to run, do not actually run them.")
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
if (inputArguments.selectionSuffix != ""): selection_suffix = "_{sS}".format(sS=inputArguments.selectionSuffix)

# L_total = L_2016 + L_2017 + L_2018 => deltaL_total/L_total = (deltaL_2016/L_2016)*(L_2016/L_total) + (deltaL_2017/L_2017)*(L_2017/L_total) + (deltaL_2018/L_2018)*(L_2018/L_total)
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-17-001/index.html, the 2016 uncertainty is 2.5 percent
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-17-004/index.html, the 2017 uncertainty is 2.3 percent
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-18-002/index.html, the 2018 uncertainty is 2.5 percent

lumi_uncertainty = 0.025
integrated_lumi_strings = {
    "2016": "35524.0",
    "2017": "41857.4",
    "2018": "58872.7"
}

if (inputArguments.chain == "all"):
    runSequence = ["data", "MC", "combine", "signalContamination", "ancillaryPlots"] # "limits" should be run manually at the end once all the combine jobs have finished running
else:
    runSequence = [inputArguments.chain]

yearPattern = inputArguments.year
if (yearPattern == "all"):
    yearPattern = "*"

list_signalTypes = ["signal"]
abbreviated_signalTypes = {
    "signal": "s",
    "signal_loose": "l"
}
if (inputArguments.addLooseSignal):
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

def checkAndEstablishLock(optional_identifier):
    if (os.path.isfile("{aR}/analysis{oI}_running.lock".format(aR=stealthEnv.analysisRoot, oI=optional_identifier))):
        sys.exit("ERROR: {aR}/analysis{oI}_running.lock already exists!".format(aR=stealthEnv.analysisRoot, oI=optional_identifier))
    os.system("touch {aR}/analysis{oI}_running.lock".format(aR=stealthEnv.analysisRoot, oI=optional_identifier))

def removeLock(optional_identifier):
    os.system("rm -v -f {aR}/analysis{oI}_running.lock".format(aR=stealthEnv.analysisRoot, oI=optional_identifier))

def signal_handler(sig, frame):
    removeLock(optional_identifier)
    sys.exit("Terminated by user.")
signal.signal(signal.SIGINT, signal_handler)

def execute_in_env(commandToRun, optional_identifier, dryRun):
    env_setup_command = "bash -c \"cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot)
    runInEnv = "{e_s_c} && set -x && {c} && set +x\"".format(e_s_c=env_setup_command, c=commandToRun)
    if (dryRun):
        print("Dry-run, not executing:")
        print("{c}".format(c=runInEnv))
    else:
        print("Executing:")
        print("{c}".format(c=runInEnv))
        returnCode = os.system(runInEnv)
        if (returnCode != 0):
            removeLock(optional_identifier)
            sys.exit("ERROR: command \"{c}\" returned status {rC}".format(c=commandToRun, rC=returnCode))

def transfer_file_to_EOS_area(sourceFile, targetDirectory, optional_identifier, dryRun):
    sourceFileName = (sourceFile.split("/"))[-1]
    command = "xrdcp --verbose --force --path --streams 15 {sF} {tD}/{sFN}".format(sF=sourceFile, tD=targetDirectory, sFN=sourceFileName)
    execute_in_env(command, optional_identifier, dryRun)

def run_data_step(outputDirectory, inputFilesList, outputPrefix, analyzeSignalBins, optional_identifier, isDryRun):
    for outputSubdirectory in ["dataEventHistograms", "dataSystematics"]:
        os.system("mkdir -p {oD}/{oS}".format(oD=outputDirectory, oS=outputSubdirectory))
    command = ("./getDataEventHistogramsAndSystematics.py --inputFilesList {iFL} --outputDirectory_eventHistograms {oD}/dataEventHistograms/ --outputDirectory_dataSystematics {oD}/dataSystematics/ --outputPrefix {oP}".format(iFL=inputFilesList, oD=outputDirectory, oP=outputPrefix))
    if (analyzeSignalBins): command += " --analyzeSignalBins"
    execute_in_env(command, optional_identifier, isDryRun)

def run_MC_chain(outputDirectory, eventProgenitor, crossSectionsFilePath, dataPrefix, outputPrefix, inputMCPathMain, integratedLuminosityMainString, inputMCPathsAux, integratedLuminositiesAux, MCTemplatePath, getSignalContaminationOutsideSidebands, optional_identifier, isDryRun):
    for outputSubdirectory in ["MCEventHistograms", "MCSystematics", "signalContamination"]:
        os.system("mkdir -p {oD}/{oS}".format(oD=outputDirectory, oS=outputSubdirectory))
    command_update = ("cd getMCSystematics && make && cd ..")
    execute_in_env(command_update, optional_identifier, isDryRun)
    command_getHists = ("./getMCSystematics/bin/getEventHistograms eventProgenitor={eP} crossSectionsFilePath={cSFP} inputMCPathMain={iMCPM} integratedLuminosityMain={iLM} outputDirectory={oD}/MCEventHistograms/ outputPrefix={oP} MCTemplatePath={MTP}".format(eP=eventProgenitor, cSFP=crossSectionsFilePath, iMCPM=inputMCPathMain, iLM=integratedLuminosityMainString, oD=outputDirectory, oP=outputPrefix, MTP=MCTemplatePath))
    if (len(inputMCPathsAux) != 0):
        command_getHists += " inputMCPathsAux="
        for inputMCPathAux in inputMCPathsAux:
            command_getHists += (inputMCPathAux + "\;")
        command_getHists = command_getHists[:-2] # to remove the last "\;"
        command_getHists += " integratedLuminositiesAux="
        for integratedLuminosityAux in integratedLuminositiesAux:
            command_getHists += (integratedLuminosityAux + "\;")
        command_getHists = command_getHists[:-2] # to remove the last "\;"
    execute_in_env(command_getHists, optional_identifier, isDryRun)
    signalContaminationOutsideSidebandsString = "false" # the string, not the bool
    if getSignalContaminationOutsideSidebands:
        signalContaminationOutsideSidebandsString = "true"
    command_getSystematics = ("./getMCSystematics/bin/getMCUncertainties inputPath={oD}/MCEventHistograms/{oP}_savedObjects.root MCTemplatePath={MTP} inputNEventsFile={oD}/dataSystematics/{dP}_observedEventCounters.dat inputDataUncertaintiesFile={oD}/dataSystematics/{dP}_dataSystematics.dat inputDataSTScalingUncertaintiesFile={oD}/dataSystematics/control_dataSystematics_scaling.dat outputDirectory={oD}/MCSystematics/ outputDirectory_signalContamination={oD}/signalContamination/ outputPrefix={oP} getSignalContaminationOutsideSidebands={sCOSS}".format(oD=outputDirectory, oP=outputPrefix, MTP=MCTemplatePath, dP=dataPrefix, sCOSS=signalContaminationOutsideSidebandsString))
    execute_in_env(command_getSystematics, optional_identifier, isDryRun)
    transfer_file_to_EOS_area("{oD}/MCEventHistograms/{oP}_savedObjects.root".format(oD=outputDirectory, oP=outputPrefix), "{eP}/{sER}/analysisEOSAreas/analysis{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), optional_identifier, isDryRun)
    transfer_file_to_EOS_area("{oD}/MCSystematics/{oP}_MCUncertainties_savedObjects.root".format(oD=outputDirectory, oP=outputPrefix), "{eP}/{sER}/analysisEOSAreas/analysis{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), optional_identifier, isDryRun)

def run_combine_chain(outputDirectory, eventProgenitor, minNeutralinoMass, crossSectionsFilePath, combineResultsDirectory, path_dataSystematics_signal, path_dataSystematics_signal_loose, path_dataSystematics_control, path_dataObservedEventCounters_signal, path_dataObservedEventCounters_signal_loose, path_dataObservedEventCounters_control, path_dataExpectedEventCounters_signal, path_dataExpectedEventCounters_signal_loose, path_dataExpectedEventCounters_control, MCPrefix_signal, MCPrefix_signal_loose, MCPrefix_control, outputPrefix, MCTemplatePath, luminosity_uncertainty, addLooseSignal, runUnblinded, EOSAnalysisArea, optional_identifier, isDryRun):
    command_submitCombineJobs = ("./submitCombineToolJobs.py --dataCardsPrefix {oP} --outputDirectory {eP}/{sER}/combineToolOutputs/{cRD}/ --eventProgenitor {eP2} --MCTemplatePath {MTP} --minNeutralinoMass {mNM} --crossSectionsFileName {cSFP} --path_dataSystematics_signal {pDSS} --path_dataSystematics_signal_loose {pDSSL} --path_dataSystematics_control {pDSC} --path_dataObservedEventCounters_signal {pDOECS} --path_dataObservedEventCounters_signal_loose {pDOECSL} --path_dataObservedEventCounters_control {pDOECC} --path_dataExpectedEventCounters_signal {pDEECS} --path_dataExpectedEventCounters_signal_loose {pDEECSL} --path_dataExpectedEventCounters_control {pDEECC} --MCHistogramsSignal {EAA}/{pS}_savedObjects.root --MCHistogramsSignalLoose {EAA}/{pSL}_savedObjects.root --MCHistogramsControl {EAA}/{pC}_savedObjects.root --MCUncertaintiesSignal {EAA}/{pS}_MCUncertainties_savedObjects.root --MCUncertaintiesSignalLoose {EAA}/{pSL}_MCUncertainties_savedObjects.root --MCUncertaintiesControl {EAA}/{pC}_MCUncertainties_savedObjects.root --luminosityUncertainty {lU} --EOSAnalysisArea {EAA}".format(oD=outputDirectory, oP=outputPrefix, eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, MTP=MCTemplatePath, mNM=minNeutralinoMass, eP2=eventProgenitor, cRD=combineResultsDirectory, cSFP=crossSectionsFilePath, pDSS=path_dataSystematics_signal, pDSSL=path_dataSystematics_signal_loose, pDSC=path_dataSystematics_control, pDOECS=path_dataObservedEventCounters_signal, pDOECSL=path_dataObservedEventCounters_signal_loose, pDOECC=path_dataObservedEventCounters_control, pDEECS=path_dataExpectedEventCounters_signal, pDEECSL=path_dataExpectedEventCounters_signal_loose, pDEECC=path_dataExpectedEventCounters_control, pS=MCPrefix_signal, pSL=MCPrefix_signal_loose, pC=MCPrefix_control, lU=luminosity_uncertainty, EAA=EOSAnalysisArea))
    if (inputArguments.optionalIdentifier != ""): command_submitCombineJobs += " --optionalIdentifier {oI}".format(oI=inputArguments.optionalIdentifier) # Just "inputArguments.optionalIdentifier", without the underscore
    if (runUnblinded): command_submitCombineJobs += " --runUnblinded"
    if (addLooseSignal): command_submitCombineJobs += " --addLooseSignal"
    if (isDryRun): command_submitCombineJobs += " --isDryRun"
    execute_in_env(command_submitCombineJobs, optional_identifier, isDryRun)

def produce_STComparisons(outputDirectory, controlDataPath, outputFilePrefix_STComparisons, optional_identifier, isDryRun):
    for outputSubdirectory in ["dataEventHistograms", "dataSystematics", "publicationPlots"]:
        os.system("mkdir -p {oD}/{oS}".format(oD=outputDirectory, oS=outputSubdirectory))
    command_controlDistributions = "./plotSTDistributionComparisons.py --inputFilePath {cDP} --outputDirectory {oD}/publicationPlots --outputFilePrefix {oFP}".format(cDP=controlDataPath, oD=outputDirectory, oFP=outputFilePrefix_STComparisons)
    execute_in_env(command_controlDistributions, optional_identifier, isDryRun)

def produce_ancillary_plots_control(outputDirectory, eventProgenitor, path_data_expectedNEvents, path_data_observedNEvents, path_MC_weightedNEvents, path_dataSystematics, path_STScalingSystematics, optional_identifier, isDryRun):
    for outputSubdirectory in ["dataEventHistograms", "dataSystematics", "publicationPlots"]:
        os.system("mkdir -p {oD}/{oS}".format(oD=outputDirectory, oS=outputSubdirectory))
    command_controlSTDistributions_dataAndSignal = "./plotSTDistributionsWithErrors.py --eventProgenitor {eP} --path_data_expectedNEvents {pDENE} --path_data_observedNEvents {pDONE} --path_MC_weightedNEvents {pMCWNE} --path_dataSystematics {pDS} --path_STScalingSystematics {pSTSS} --outputDirectory {oD}/publicationPlots/ --outputFilePrefix {oFP} --plotObservedData".format(eP=eventProgenitor, pDENE=path_data_expectedNEvents, pDONE=path_data_observedNEvents, pMCWNE=path_MC_weightedNEvents, pDS=path_dataSystematics, pSTSS=path_STScalingSystematics, oD=outputDirectory, oFP="STDistributions_{eP}_control".format(eP=eventProgenitor))
    execute_in_env(command_controlSTDistributions_dataAndSignal, optional_identifier, isDryRun)

def produce_ancillary_plots_signal(outputDirectory, eventProgenitor, signalType, path_data_expectedNEvents, path_data_observedNEvents, path_MC_weightedNEvents, path_dataSystematics, path_STScalingSystematics, runUnblinded, optional_identifier, isDryRun):
    for outputSubdirectory in ["dataEventHistograms", "dataSystematics", "publicationPlots"]:
        os.system("mkdir -p {oD}/{oS}".format(oD=outputDirectory, oS=outputSubdirectory))
    command_signalSTDistributions_dataAndSignal = "./plotSTDistributionsWithErrors.py --eventProgenitor {eP} --path_data_expectedNEvents {pDENE} --path_data_observedNEvents {pDONE} --path_MC_weightedNEvents {pMCWNE} --path_dataSystematics {pDS} --path_STScalingSystematics {pSTSS} --outputDirectory {oD}/publicationPlots/ --outputFilePrefix {oFP}".format(eP=eventProgenitor, pDENE=path_data_expectedNEvents, pDONE=path_data_observedNEvents, pMCWNE=path_MC_weightedNEvents, pDS=path_dataSystematics, pSTSS=path_STScalingSystematics, oD=outputDirectory, oFP="STDistributions_{eP}_{sT}".format(eP=eventProgenitor, sT=signalType))
    if runUnblinded:
        command_signalSTDistributions_dataAndSignal += " --plotObservedData"
    execute_in_env(command_signalSTDistributions_dataAndSignal, optional_identifier, isDryRun)

def get_signal_contamination(outputDirectory, crossSectionsFilePath, eventProgenitor, dataPrefix, outputPrefix, inputMCPathMain, integratedLuminosityMainString, inputMCPathsAux, integratedLuminositiesAux, MCTemplatePath, optional_identifier, isDryRun):
    for outputSubdirectory in ["MCEventHistograms", "MCSystematics", "signalContamination"]:
        os.system("mkdir -p {oD}/{oS}".format(oD=outputDirectory, oS=outputSubdirectory))
    command_update = ("cd getMCSystematics && make && cd ..")
    execute_in_env(command_update, optional_identifier, isDryRun)
    command_getHists = ("./getMCSystematics/bin/getEventHistograms eventProgenitor={eP} crossSectionsFilePath={cSFP} inputMCPathMain={iMCPM} integratedLuminosityMain={iLM} outputDirectory={oD}/MCEventHistograms/ outputPrefix={oP} MCTemplatePath={MTP}".format(eP=eventProgenitor, cSFP=crossSectionsFilePath, iMCPM=inputMCPathMain, iLM=integratedLuminosityMainString, oD=outputDirectory, oP=outputPrefix, MTP=MCTemplatePath))
    if (len(inputMCPathsAux) != 0):
        command_getHists += " inputMCPathsAux="
        for inputMCPathAux in inputMCPathsAux:
            command_getHists += (inputMCPathAux + "\;")
        command_getHists = command_getHists[:-2] # to remove the last "\;"
        command_getHists += " integratedLuminositiesAux="
        for integratedLuminosityAux in integratedLuminositiesAux:
            command_getHists += (integratedLuminosityAux + "\;")
        command_getHists = command_getHists[:-2] # to remove the last "\;"
    execute_in_env(command_getHists, optional_identifier, isDryRun)
    command_getSignalContamination = ("./getMCSystematics/bin/getMCUncertainties inputPath={oD}/MCEventHistograms/{oP}_savedObjects.root MCTemplatePath={MTP} inputNEventsFile={oD}/dataSystematics/{dP}_observedEventCounters.dat inputDataUncertaintiesFile={oD}/dataSystematics/{dP}_dataSystematics.dat inputDataSTScalingUncertaintiesFile={oD}/dataSystematics/control_dataSystematics_scaling.dat outputDirectory={oD}/MCSystematics/ outputDirectory_signalContamination={oD}/signalContamination/ outputPrefix={oP} getSignalContaminationOutsideSidebands=true".format(oD=outputDirectory, oP=outputPrefix, MTP=MCTemplatePath, dP=dataPrefix))
    execute_in_env(command_getSignalContamination, optional_identifier, isDryRun)

def plot_limits(outputDirectory, crossSectionsFilePath, eventProgenitor, combineResultsDirectory, MCTemplatePath, minNeutralinoMass, runUnblinded, optional_identifier, isDryRun):
    os.system("mkdir -p {oD}/publicationPlots && mkdir -p {oD}/limits".format(oD=outputDirectory))
    command_plotLimits = "condor_q && ./plotLimits.py --crossSectionsFile {cSFP} --MCTemplatePath {MTP} --eventProgenitor {eP2} --combineResultsDirectory {eP}/{sER}/combineToolOutputs/{cRD} --combineOutputPrefix {eP2} --outputDirectory_rawOutput {oD}/limits --outputDirectory_plots {oD}/publicationPlots --outputSuffix {eP2} --minNeutralinoMass {mNM}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, MTP=MCTemplatePath, cRD=combineResultsDirectory, eP2=eventProgenitor, cSFP=crossSectionsFilePath, oD=outputDirectory, mNM=minNeutralinoMass)
    if (runUnblinded): command_plotLimits += " --plotObserved"
    execute_in_env(command_plotLimits, optional_identifier, isDryRun)

checkAndEstablishLock(optional_identifier)
command_createEOSAnalysisAreas = ("eos {eP} mkdir -p {sER}/analysisEOSAreas/analysis{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier))
execute_in_env(command_createEOSAnalysisAreas, optional_identifier, inputArguments.isDryRun)

for step in runSequence:
    if (step == "data"):
        run_data_step(outputDirectory="{aR}/analysis{oI}".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), inputFilesList="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_data_{yP}_control_fakefake.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, yP=yearPattern), outputPrefix="control", analyzeSignalBins=True, optional_identifier=optional_identifier, isDryRun=inputArguments.isDryRun)
        for signalType in list_signalTypes:
            run_data_step(outputDirectory="{aR}/analysis{oI}".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), inputFilesList="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_data_{yP}_{sT}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, yP=yearPattern, sT=signalType), outputPrefix="{sT}".format(sT=signalType), analyzeSignalBins=inputArguments.runUnblinded, optional_identifier=optional_identifier, isDryRun=inputArguments.isDryRun)
    elif (step == "MC"):
        for eventProgenitor in eventProgenitors:
            crossSectionsPath = crossSectionsForProgenitor[eventProgenitor]
            for signalType in list_signalTypes:
                MCPathMain = ""
                lumiMain = ""
                MCPathsAux = []
                lumisAux = []
                if (inputArguments.year == "all"):
                    MCPathMain = "{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_{tD}_2017_{sT}.root".format(sT=signalType, eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, tD=tDesignationsForProgenitor[eventProgenitor])
                    lumiMain = integrated_lumi_strings["2017"]
                    MCPathsAux = ["{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_{tD}_2016_{sT}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, sT=signalType, tD=tDesignationsForProgenitor[eventProgenitor]), "{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_{tD}_2018_{sT}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, sT=signalType, tD=tDesignationsForProgenitor[eventProgenitor])]
                    lumisAux = [integrated_lumi_strings["2016"], integrated_lumi_strings["2018"]]
                else:
                    MCPathMain = "{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_{tD}_{y}_{sT}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, y=inputArguments.year, sT=signalType, tD=tDesignationsForProgenitor[eventProgenitor])
                    lumiMain = integrated_lumi_strings[inputArguments.year]
                run_MC_chain(outputDirectory="{aR}/analysis{oI}".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), eventProgenitor=eventProgenitor, crossSectionsFilePath=crossSectionsPath, dataPrefix="{sT}".format(sT=signalType), outputPrefix="MC_stealth_{eP}_{y}_{sT}".format(eP=eventProgenitor, y=inputArguments.year, sT=signalType), inputMCPathMain=MCPathMain, integratedLuminosityMainString=lumiMain, inputMCPathsAux=MCPathsAux, integratedLuminositiesAux=lumisAux, MCTemplatePath=MCTemplatesForProgenitor[eventProgenitor], getSignalContaminationOutsideSidebands=False, optional_identifier=optional_identifier, isDryRun=inputArguments.isDryRun)
    elif (step == "combine"):
        # Make sure the CMSSW source tarball is the latest version
        print("Updating and uploading CMSSW source tarball...")
        updateCommand = "cd {sCB}/.. && ./uploadTarball.sh && cd {sR}".format(sCB=stealthEnv.stealthCMSSWBase, sR=stealthEnv.stealthRoot)
        os.system(updateCommand)
        command_cleanCombineResultsDirectory = "eos {eP} ls {sER}/combineToolOutputs/{cRD} && eos {eP} rm -r {sER}/combineToolOutputs/combineResults{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, cRD="combineResults{oI}".format(oI=optional_identifier), oI=optional_identifier)
        execute_in_env(command_cleanCombineResultsDirectory, optional_identifier, inputArguments.isDryRun)
        command_createEOSDirectory = ("eos {eP} mkdir -p {sER}/combineToolOutputs/combineResults{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier))
        execute_in_env(command_createEOSDirectory, optional_identifier, inputArguments.isDryRun)
        for eventProgenitor in eventProgenitors:
            crossSectionsPath = crossSectionsForProgenitor[eventProgenitor]
            run_combine_chain(outputDirectory="{aR}/analysis{oI}".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), eventProgenitor=eventProgenitor, minNeutralinoMass=minNeutralinoMassToPlot[eventProgenitor], crossSectionsFilePath=crossSectionsPath, combineResultsDirectory="combineResults{oI}".format(oI=optional_identifier), path_dataSystematics_signal="{aR}/analysis{oI}/dataSystematics/signal_dataSystematics.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), path_dataSystematics_signal_loose="{aR}/analysis{oI}/dataSystematics/signal_loose_dataSystematics.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), path_dataSystematics_control="{aR}/analysis{oI}/dataSystematics/control_dataSystematics.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), path_dataObservedEventCounters_signal="{aR}/analysis{oI}/dataSystematics/signal_observedEventCounters.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), path_dataObservedEventCounters_signal_loose="{aR}/analysis{oI}/dataSystematics/signal_loose_observedEventCounters.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), path_dataObservedEventCounters_control="{aR}/analysis{oI}/dataSystematics/control_observedEventCounters.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), path_dataExpectedEventCounters_signal="{aR}/analysis{oI}/dataSystematics/signal_eventCounters.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), path_dataExpectedEventCounters_signal_loose="{aR}/analysis{oI}/dataSystematics/signal_loose_eventCounters.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), path_dataExpectedEventCounters_control="{aR}/analysis{oI}/dataSystematics/control_eventCounters.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), MCPrefix_signal="MC_stealth_{eP}_{y}_signal".format(eP=eventProgenitor, y=inputArguments.year), MCPrefix_signal_loose="MC_stealth_{eP}_{y}_signal_loose".format(eP=eventProgenitor, y=inputArguments.year), MCPrefix_control="MC_stealth_{eP}_{y}_control".format(eP=eventProgenitor, y=inputArguments.year), outputPrefix="{eP}".format(eP=eventProgenitor), MCTemplatePath=MCTemplatesForProgenitor[eventProgenitor], luminosity_uncertainty=lumi_uncertainty, addLooseSignal=inputArguments.addLooseSignal, runUnblinded=inputArguments.runUnblinded, EOSAnalysisArea="{eP}/{sER}/analysisEOSAreas/analysis{oI}".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, oI=optional_identifier), optional_identifier=optional_identifier, isDryRun=inputArguments.isDryRun)
    elif (step == "ancillaryPlots"):
        produce_STComparisons(outputDirectory="{aR}/analysis{oI}".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), controlDataPath="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_data_{yP}_control_fakefake.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, yP=yearPattern), outputFilePrefix_STComparisons="control_STComparisons", optional_identifier=optional_identifier, isDryRun=inputArguments.isDryRun)
        produce_STComparisons(outputDirectory="{aR}/analysis{oI}".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), controlDataPath="{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_data_singlemedium_{yP}_control_singlemedium.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, yP=yearPattern), outputFilePrefix_STComparisons="control_singlemedium_STComparisons", optional_identifier=optional_identifier, isDryRun=inputArguments.isDryRun)
        for eventProgenitor in eventProgenitors:
            crossSectionsPath = crossSectionsForProgenitor[eventProgenitor]
            produce_ancillary_plots_control(outputDirectory="{aR}/analysis{oI}".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), eventProgenitor=eventProgenitor, path_data_expectedNEvents="{aR}/analysis{oI}/dataSystematics/control_eventCounters.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), path_data_observedNEvents="{aR}/analysis{oI}/dataSystematics/control_observedEventCounters.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), path_MC_weightedNEvents="{aR}/analysis{oI}/MCEventHistograms/MC_stealth_{eP}_{y}_control_savedObjects.root".format(aR=stealthEnv.analysisRoot, oI=optional_identifier, eP=eventProgenitor, y=inputArguments.year), path_dataSystematics="{aR}/analysis{oI}/dataSystematics/control_dataSystematics.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), path_STScalingSystematics="{aR}/analysis{oI}/dataSystematics/control_dataSystematics_scaling.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), optional_identifier=optional_identifier, isDryRun=inputArguments.isDryRun)
            for signalType in list_signalTypes:
                produce_ancillary_plots_signal(outputDirectory="{aR}/analysis{oI}".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), eventProgenitor=eventProgenitor, signalType=signalType, path_data_expectedNEvents="{aR}/analysis{oI}/dataSystematics/{sT}_eventCounters.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier, sT=signalType), path_data_observedNEvents="{aR}/analysis{oI}/dataSystematics/{sT}_observedEventCounters.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier, sT=signalType), path_MC_weightedNEvents="{aR}/analysis{oI}/MCEventHistograms/MC_stealth_{eP}_{y}_{sT}_savedObjects.root".format(aR=stealthEnv.analysisRoot, oI=optional_identifier, eP=eventProgenitor, y=inputArguments.year, sT=signalType), path_dataSystematics="{aR}/analysis{oI}/dataSystematics/{sT}_dataSystematics.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier, sT=signalType), path_STScalingSystematics="{aR}/analysis{oI}/dataSystematics/control_dataSystematics_scaling.dat".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), runUnblinded=inputArguments.runUnblinded, optional_identifier=optional_identifier, isDryRun=inputArguments.isDryRun)
    elif (step == "signalContamination"):
        for eventProgenitor in eventProgenitors:
            crossSectionsPath = crossSectionsForProgenitor[eventProgenitor]
            MCPathMain = ""
            lumiMain = ""
            MCPathsAux = []
            lumisAux = []
            if (inputArguments.year == "all"):
                MCPathMain = "{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_{tD}_2017_control_fakefake.root".format(eP=stealthEnv.EOSPrefix, tD=tDesignationsForProgenitor[eventProgenitor], sER=stealthEnv.stealthEOSRoot, sS=selection_suffix)
                lumiMain = integrated_lumi_strings["2017"]
                MCPathsAux = ["{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_{tD}_2016_control_fakefake.root".format(eP=stealthEnv.EOSPrefix, tD=tDesignationsForProgenitor[eventProgenitor], sER=stealthEnv.stealthEOSRoot, sS=selection_suffix), "{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_{tD}_2018_control_fakefake.root".format(eP=stealthEnv.EOSPrefix, tD=tDesignationsForProgenitor[eventProgenitor], sER=stealthEnv.stealthEOSRoot, sS=selection_suffix)]
                lumisAux = [integrated_lumi_strings["2016"], integrated_lumi_strings["2018"]]
            else:
                MCPathMain = "{eP}/{sER}/selections/combined_DoublePhoton{sS}/merged_selection_MC_stealth_{tD}_{y}_control_fakefake.root".format(eP=stealthEnv.EOSPrefix, tD=tDesignationsForProgenitor[eventProgenitor], sER=stealthEnv.stealthEOSRoot, sS=selection_suffix, y=inputArguments.year)
                lumiMain = integrated_lumi_strings[inputArguments.year]
            get_signal_contamination(outputDirectory="{aR}/analysis{oI}".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), crossSectionsFilePath=crossSectionsPath, eventProgenitor=eventProgenitor, dataPrefix="control", outputPrefix="MC_stealth_{eP}_{y}_control".format(eP=eventProgenitor, y=inputArguments.year), inputMCPathMain=MCPathMain, integratedLuminosityMainString=lumiMain, inputMCPathsAux=MCPathsAux, integratedLuminositiesAux=lumisAux, MCTemplatePath=MCTemplatesForProgenitor[eventProgenitor], optional_identifier=optional_identifier, isDryRun=inputArguments.isDryRun)
    elif (step == "limits"):
        for eventProgenitor in eventProgenitors:
            crossSectionsPath = crossSectionsForProgenitor[eventProgenitor]
            plot_limits(outputDirectory="{aR}/analysis{oI}".format(aR=stealthEnv.analysisRoot, oI=optional_identifier), crossSectionsFilePath=crossSectionsPath, eventProgenitor=eventProgenitor, combineResultsDirectory="combineResults{oI}".format(oI=optional_identifier), MCTemplatePath=MCTemplatesForProgenitor[eventProgenitor], minNeutralinoMass=minNeutralinoMassToPlot[eventProgenitor], runUnblinded=inputArguments.runUnblinded, optional_identifier=optional_identifier, isDryRun=inputArguments.isDryRun)
    else:
        removeLock(optional_identifier)
        sys.exit("ERROR: Unrecognized step: {s}".format(s=step))

removeLock(optional_identifier)
