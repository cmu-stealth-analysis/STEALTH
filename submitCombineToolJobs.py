#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, os, sys, argparse, tmJDLInterface, MCTemplateReader, stealthEnv

inputArgumentsParser = argparse.ArgumentParser(description='Run Higgs combine tool on already generated datacards.')
inputArgumentsParser.add_argument('--dataCardsPrefix', default="", help='Data cards prefix.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="root://cmseos.fnal.gov//store/user/lpcsusystealth/combineToolOutputs", help='EOS path on which to store combine tool outputs.',type=str)
inputArgumentsParser.add_argument('--eventProgenitor', required=True, help="Type of stealth sample. Two possible values: \"squark\" or \"gluino\".", type=str)
inputArgumentsParser.add_argument('--MCTemplatePath', default="plot_susyMasses_template.root", help='Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.', type=str)
inputArgumentsParser.add_argument('--minNeutralinoMass', default=-1., help='Min value of the neutralino mass to plot.',type=float)
inputArgumentsParser.add_argument('--eventProgenitorMassOffset', default=-1., help='Min value of the event progenitor mass to plot is obtained by adding this offset to the template.',type=float)
inputArgumentsParser.add_argument('--minMassDifference', default=-1., help='Min difference between the masses of the event progenitor and neutralino.',type=float)
inputArgumentsParser.add_argument('--crossSectionsFileName', required=True, help='Path to dat file that contains cross-sections as a function of eventProgenitor mass, to use while weighting events.',type=str)
inputArgumentsParser.add_argument('--path_dataSystematics_signal', required=True, help='Path to root file with systematics for the signal region.', type=str)
inputArgumentsParser.add_argument('--path_dataSystematics_signal_loose', required=True, help='Path to root file with systematics for the loose signal region.', type=str)
# inputArgumentsParser.add_argument('--path_dataSystematics_control', required=True, help='Path to root file with systematics for the control region.', type=str)
inputArgumentsParser.add_argument('--path_MCShapeAdjustment_signal', required=True, help='Path to dat file containing MC shape adjustments for the signal selection.', type=str)
inputArgumentsParser.add_argument('--path_MCShapeAdjustment_signal_loose', required=True, help='Path to dat file containing MC shape adjustments for the loose signal selection.', type=str)
# inputArgumentsParser.add_argument('--path_MCShapeAdjustment_control', required=True, help='Path to dat file containing MC shape adjustment for the control selection.', type=str)
# inputArgumentsParser.add_argument('--path_DataMCRatioAdjustment_QCD_signal', required=True, help='Path to dat file containing data/QCD MC ratio adjustments for the signal selection.', type=str)
# inputArgumentsParser.add_argument('--path_DataMCRatioAdjustment_QCD_signal_loose', required=True, help='Path to dat file containing data/QCD MC ratio adjustments for the loose signal selection.', type=str)
# inputArgumentsParser.add_argument('--path_DataMCRatioAdjustment_diphoton_signal', required=True, help='Path to dat file containing data/diphoton MC ratio adjustments for the signal selection.', type=str)
# inputArgumentsParser.add_argument('--path_DataMCRatioAdjustment_diphoton_signal_loose', required=True, help='Path to dat file containing data/diphoton MC ratio adjustments for the loose signal selection.', type=str)
inputArgumentsParser.add_argument('--paths_bkgCompositionSystematic', required=True, help='Comma-separated list of paths to dat files containing diphoton MC ratio residuals for the background composition systematic.', type=str)
inputArgumentsParser.add_argument('--path_dataObservedEventCounters_signal', required=True, help='Path to root file with observed event counters for the signal region.', type=str)
inputArgumentsParser.add_argument('--path_dataObservedEventCounters_signal_loose', required=True, help='Path to root file with observed event counters for the loose signal region.', type=str)
# inputArgumentsParser.add_argument('--path_dataObservedEventCounters_control', required=True, help='Path to root file with observed event counters for the control region.', type=str)
inputArgumentsParser.add_argument('--path_dataExpectedEventCounters_signal', required=True, help='Path to root file with expected event counters for the signal region.', type=str)
inputArgumentsParser.add_argument('--path_dataExpectedEventCounters_signal_loose', required=True, help='Path to root file with expected event counters for the loose signal region.', type=str)
# inputArgumentsParser.add_argument('--path_dataExpectedEventCounters_control', required=True, help='Path to root file with expected event counters for the control region.', type=str)
inputArgumentsParser.add_argument('--MCHistogramsSignal', required=True, help='Path to root file with MC event histograms for the signal region.', type=str)
inputArgumentsParser.add_argument('--MCHistogramsSignalLoose', required=True, help='Path to root file with MC event histograms for the loose signal region.', type=str)
# inputArgumentsParser.add_argument('--MCHistogramsControl', required=True, help='Path to root file with MC event histograms for the control region.', type=str)
inputArgumentsParser.add_argument('--MCUncertaintiesSignal', required=True, help='Path to root file with MC uncertainties for the signal region.', type=str)
inputArgumentsParser.add_argument('--MCUncertaintiesSignalLoose', required=True, help='Path to root file with MC uncertainties for the loose signal region.', type=str)
# inputArgumentsParser.add_argument('--MCUncertaintiesControl', required=True, help='Path to root file with MC uncertainties for the control region.', type=str)
inputArgumentsParser.add_argument('--luminosityUncertainty', required=True, help='Uncertainty on the luminosity.', type=float)
inputArgumentsParser.add_argument('--EOSAnalysisArea', required=True, help='Path to EOS analysis area to use for saving intermediate outputs.', type=str)
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArgumentsParser.add_argument('--runUnblinded', action='store_true', help="Pass runUnblinded flag while creating datacards.")
inputArgumentsParser.add_argument('--addLooseSignal', action='store_true', help="Add loose photons in a different signal bin. By default data cards are created with only medium photons.")
inputArgumentsParser.add_argument('--isDryRun', action='store_true', help="Do not submit the actual jobs: instead, only print the shell command that would have been called.")
inputArguments = inputArgumentsParser.parse_args()

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)

# Copy combine tool helper script into the working directory
copyCommand = "mkdir -p {cWAR}/combine{oI} && cd {sR} && cp -u combineToolHelper.sh {cWAR}/combine{oI}/.".format(sR=stealthEnv.stealthRoot, cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier)
os.system(copyCommand)

# Make sure the tarballs to transfer are up to date
updateCommand = "cd {tUP} && ./update_tmUtilsTarball.sh && cd {sR}".format(tUP=stealthEnv.tmUtilsParent, sR=stealthEnv.stealthRoot)
os.system(updateCommand)

templateReader = MCTemplateReader.MCTemplateReader(inputArguments.MCTemplatePath)
for indexPair in templateReader.nextValidBin():
    eventProgenitorMassBin = indexPair[0]
    eventProgenitorMass = (templateReader.eventProgenitorMasses)[eventProgenitorMassBin]
    neutralinoMassBin = indexPair[1]
    neutralinoMass = (templateReader.neutralinoMasses)[neutralinoMassBin]
    if (eventProgenitorMass < (templateReader.minEventProgenitorMass + inputArguments.eventProgenitorMassOffset)): continue
    if (neutralinoMass < inputArguments.minNeutralinoMass): continue
    if ((eventProgenitorMass - neutralinoMass) < inputArguments.minMassDifference): continue
    print("eventProgenitor mass: {gM}, neutralino mass: {nM}".format(gM=eventProgenitorMass, nM=neutralinoMass))
    x509ProxyPath = stealthEnv.x509Proxy
    tmUtilsTarballPath = "{tUP}/tmUtils.tar.gz".format(tUP=stealthEnv.tmUtilsParent)
    tmUtilsExtractionScriptPath = "{tUP}/extract_tmUtilsTarball.sh".format(tUP=stealthEnv.tmUtilsParent)
    remoteEnvSetupScriptPath = "{sR}/setup_environment_remote.sh".format(sR=stealthEnv.stealthRoot)
    MCTemplateReaderPath = "{sR}/MCTemplateReader.py".format(sR=stealthEnv.stealthRoot)
    crossSectionsFilePath = "{sR}/{cSFN}".format(sR=stealthEnv.stealthRoot, cSFN=inputArguments.crossSectionsFileName)
    STRegionBoundariesFilePath = "{sR}/STRegionBoundaries.dat".format(sR=stealthEnv.stealthRoot)
    dataSystematicsPath_signal = inputArguments.path_dataSystematics_signal
    dataSystematicsPath_signal_loose = inputArguments.path_dataSystematics_signal_loose
    # dataSystematicsPath_control = inputArguments.path_dataSystematics_control
    MCShapeAdjustmentPath_signal = inputArguments.path_MCShapeAdjustment_signal
    MCShapeAdjustmentPath_signal_loose = inputArguments.path_MCShapeAdjustment_signal_loose
    # MCShapeAdjustmentPath_control = inputArguments.path_MCShapeAdjustment_control
    # DataMCRatioAdjustmentPath_QCD_signal = inputArguments.path_DataMCRatioAdjustment_QCD_signal
    # DataMCRatioAdjustmentPath_QCD_signal_loose = inputArguments.path_DataMCRatioAdjustment_QCD_signal_loose
    # DataMCRatioAdjustmentPath_diphoton_signal = inputArguments.path_DataMCRatioAdjustment_diphoton_signal
    # DataMCRatioAdjustmentPath_diphoton_signal_loose = inputArguments.path_DataMCRatioAdjustment_diphoton_signal_loose
    bkgCompositionSystematicFilePaths = (inputArguments.paths_bkgCompositionSystematic).split(",")
    dataObservedEventCountersPath_signal = inputArguments.path_dataObservedEventCounters_signal
    dataObservedEventCountersPath_signal_loose = inputArguments.path_dataObservedEventCounters_signal_loose
    # dataObservedEventCountersPath_control = inputArguments.path_dataObservedEventCounters_control
    dataExpectedEventCountersPath_signal = inputArguments.path_dataExpectedEventCounters_signal
    dataExpectedEventCountersPath_signal_loose = inputArguments.path_dataExpectedEventCounters_signal_loose
    # dataExpectedEventCountersPath_control = inputArguments.path_dataExpectedEventCounters_control
    createDataCardScriptPath = "{sR}/createDataCard.py".format(sR=stealthEnv.stealthRoot)
    commonPyFunctionsFilePath = "{sR}/commonFunctions.py".format(sR=stealthEnv.stealthRoot)
    limitsConvergenceCheckScriptPath = "{sR}/checkLimitsConvergence.py".format(sR=stealthEnv.stealthRoot)

    # filesToTransfer = [x509ProxyPath, tmUtilsTarballPath, tmUtilsExtractionScriptPath, remoteEnvSetupScriptPath, MCTemplateReaderPath, crossSectionsFilePath, STRegionBoundariesFilePath, dataSystematicsPath_signal, dataSystematicsPath_control, MCShapeAdjustmentPath_signal, MCShapeAdjustmentPath_signal_loose, MCShapeAdjustmentPath_control, dataObservedEventCountersPath_signal, dataObservedEventCountersPath_control, dataExpectedEventCountersPath_signal, dataExpectedEventCountersPath_control, createDataCardScriptPath, commonPyFunctionsFilePath, limitsConvergenceCheckScriptPath]
    filesToTransfer = [x509ProxyPath, tmUtilsTarballPath, tmUtilsExtractionScriptPath, remoteEnvSetupScriptPath, MCTemplateReaderPath, crossSectionsFilePath, STRegionBoundariesFilePath, dataSystematicsPath_signal, MCShapeAdjustmentPath_signal, MCShapeAdjustmentPath_signal_loose, dataObservedEventCountersPath_signal, dataExpectedEventCountersPath_signal, createDataCardScriptPath, commonPyFunctionsFilePath, limitsConvergenceCheckScriptPath]
    filesToTransfer.extend(bkgCompositionSystematicFilePaths)
    if (inputArguments.addLooseSignal): filesToTransfer.extend([dataSystematicsPath_signal_loose, dataObservedEventCountersPath_signal_loose, dataExpectedEventCountersPath_signal_loose])
    processIdentifier = "combineJob_{prefix}_eventProgenitorMassBin{gMB}_neutralinoMassBin{nMB}".format(prefix=inputArguments.dataCardsPrefix, gMB=eventProgenitorMassBin, nMB=neutralinoMassBin)
    jdlInterface = tmJDLInterface.tmJDLInterface(processName=processIdentifier, scriptPath="combineToolHelper.sh", outputDirectoryRelativePath="{cWAR}/combine{oI}".format(cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier))  # works even if "outputDirectoryRelativePath" is an absolute path
    jdlInterface.addFilesToTransferFromList(filesToTransfer)
    # Arguments for script:
    jdlInterface.addScriptArgument("{oD}".format(oD=inputArguments.outputDirectory)) # Argument 1: output directory
    jdlInterface.addScriptArgument("{dCP}".format(dCP=inputArguments.dataCardsPrefix)) # Argument 2: data cards prefix
    jdlInterface.addScriptArgument("{gMB}".format(gMB=eventProgenitorMassBin)) # Argument 3: eventProgenitor mass bin index
    jdlInterface.addScriptArgument("{nMB}".format(nMB=neutralinoMassBin)) # Argument 4: neutralino mass bin index
    jdlInterface.addScriptArgument("{cSFN}".format(cSFN=inputArguments.crossSectionsFileName)) # Argument 5: cross-sections file name
    jdlInterface.addScriptArgument("{MCTP}".format(MCTP=inputArguments.MCTemplatePath)) # Argument 6: MC template path
    jdlInterface.addScriptArgument("{MCHS}".format(MCHS=inputArguments.MCHistogramsSignal)) # Argument 7: path to MC event histograms, signal
    jdlInterface.addScriptArgument("{MCHSl}".format(MCHSl=inputArguments.MCHistogramsSignalLoose)) # Argument 8: path to MC event histograms, loose signal
    # jdlInterface.addScriptArgument("{MCHC}".format(MCHC=inputArguments.MCHistogramsControl)) # Argument 9: path to MC event histograms, control
    jdlInterface.addScriptArgument("{MCUS}".format(MCUS=inputArguments.MCUncertaintiesSignal)) # Argument 9: path to MC uncertainties, signal
    jdlInterface.addScriptArgument("{MCUSl}".format(MCUSl=inputArguments.MCUncertaintiesSignalLoose)) # Argument 10: path to MC uncertainties, loose signal
    # jdlInterface.addScriptArgument("{MCUC}".format(MCUC=inputArguments.MCUncertaintiesControl)) # Argument 11: path to MC uncertainties, control
    jdlInterface.addScriptArgument("{LU:.4f}".format(LU=inputArguments.luminosityUncertainty)) # Argument 11: lumi uncertainty
    jdlInterface.addScriptArgument("{EAA}".format(EAA=inputArguments.EOSAnalysisArea)) # Argument 12: EOS analysis area path
    runUnblindedString = "false"
    if (inputArguments.runUnblinded): runUnblindedString = "true"
    jdlInterface.addScriptArgument("{aLSS}".format(aLSS=runUnblindedString)) # Argument 13: run unblinded switch
    addLooseSignalString = "false"
    if (inputArguments.addLooseSignal): addLooseSignalString = "true"
    jdlInterface.addScriptArgument("{aLSS}".format(aLSS=addLooseSignalString)) # Argument 14: add loose signal switch
    # Write JDL
    jdlInterface.writeToFile()
    submissionCommand = "cd {cWAR}/combine{oI}/ && condor_submit {pI}.jdl && cd {sR}".format(pI=processIdentifier, cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier, sR=stealthEnv.stealthRoot)
    print ("Generated command: {sC}".format(sC=submissionCommand))
    if (inputArguments.isDryRun):
        print("Not submitting due to dryRun flag.")
    else:
        os.system(submissionCommand)
        print ("Submitted.")
