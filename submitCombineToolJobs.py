#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, os, sys, argparse, tmJDLInterface, MCTemplateReader, stealthEnv

inputArgumentsParser = argparse.ArgumentParser(description='Run Higgs combine tool on already generated datacards.')
inputArgumentsParser.add_argument('--dataCardsPrefix', default="", help='Data cards prefix.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="root://cmseos.fnal.gov//store/user/lpcsusystealth/combineToolOutputs", help='EOS path on which to store combine tool outputs.',type=str)
inputArgumentsParser.add_argument('--eventProgenitor', required=True, help="Type of stealth sample. Two possible values: \"squark\" or \"gluino\".", type=str)
inputArgumentsParser.add_argument('--MCTemplatePath', default="plot_susyMasses_template.root", help='Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.', type=str)
inputArgumentsParser.add_argument('--crossSectionsFileName', required=True, help='Path to dat file that contains cross-sections as a function of eventProgenitor mass, to use while weighting events.',type=str)
inputArgumentsParser.add_argument('--path_dataSystematics_signal', required=True, help='Path to root file with systematics for the signal region.', type=str)
inputArgumentsParser.add_argument('--path_dataSystematics_signal_loose', required=True, help='Path to root file with systematics for the loose signal region.', type=str)
inputArgumentsParser.add_argument('--path_dataSystematics_control', required=True, help='Path to root file with systematics for the control region.', type=str)
inputArgumentsParser.add_argument('--path_dataObservedEventCounters_signal', required=True, help='Path to root file with observed event counters for the signal region.', type=str)
inputArgumentsParser.add_argument('--path_dataObservedEventCounters_signal_loose', required=True, help='Path to root file with observed event counters for the loose signal region.', type=str)
inputArgumentsParser.add_argument('--path_dataObservedEventCounters_control', required=True, help='Path to root file with observed event counters for the control region.', type=str)
inputArgumentsParser.add_argument('--path_dataExpectedEventCounters_signal', required=True, help='Path to root file with expected event counters for the signal region.', type=str)
inputArgumentsParser.add_argument('--path_dataExpectedEventCounters_signal_loose', required=True, help='Path to root file with expected event counters for the loose signal region.', type=str)
inputArgumentsParser.add_argument('--path_dataExpectedEventCounters_control', required=True, help='Path to root file with expected event counters for the control region.', type=str)
inputArgumentsParser.add_argument('--MCHistogramsSignal', required=True, help='Path to root file with MC event histograms for the signal region.', type=str)
inputArgumentsParser.add_argument('--MCHistogramsSignalLoose', required=True, help='Path to root file with MC event histograms for the loose signal region.', type=str)
inputArgumentsParser.add_argument('--MCHistogramsControl', required=True, help='Path to root file with MC event histograms for the control region.', type=str)
inputArgumentsParser.add_argument('--MCUncertaintiesSignal', required=True, help='Path to root file with MC uncertainties for the signal region.', type=str)
inputArgumentsParser.add_argument('--MCUncertaintiesSignalLoose', required=True, help='Path to root file with MC uncertainties for the loose signal region.', type=str)
inputArgumentsParser.add_argument('--MCUncertaintiesControl', required=True, help='Path to root file with MC uncertainties for the control region.', type=str)
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

# The following are semi-educated guesses based on a test running of the combine tool. They are simply to try and reach convergence without fiddling with rMax, they should not change the tool output.
def get_rmax_lowNeutralinoMass(gluinoMass):
    rMax = 20.
    if (gluinoMass < 1125.0):
        rMax = 20.0
    elif (gluinoMass < 1375.0):
        rMax = 50.0
    elif (gluinoMass < 1625.0):
        rMax = 100.0
    elif (gluinoMass < 1825.0):
        rMax = 500.0
    elif (gluinoMass < 2075.0):
        rMax = 1000.0
    elif (gluinoMass < 2125.0):
        rMax = 5000.0
    elif (gluinoMass < 2375.0):
        rMax = 10000.0
    return rMax

def get_rmax_bulk(gluinoMass):
    rMax = 20.
    if (gluinoMass < 1125.0):
        rMax = 0.1
    elif (gluinoMass < 1475.0):
        rMax = 0.5
    elif (gluinoMass < 1625.0):
        rMax = 1.0
    elif (gluinoMass < 1825.0):
        rMax = 5.0
    elif (gluinoMass < 1975.0):
        rMax = 10.0
    elif (gluinoMass < 2125.0):
        rMax = 20.0
    elif (gluinoMass < 2375.0):
        rMax = 100.0
    return rMax

def get_rmax(eventProgenitorMass, neutralinoMass, eventProgenitor):
    rMax = 20.
    if (neutralinoMass < 118.75):
        rMax = get_rmax_lowNeutralinoMass(eventProgenitorMass)
    elif (neutralinoMass < 137.5):
        rMax = 0.075*get_rmax_lowNeutralinoMass(eventProgenitorMass) # 0.075 is an empirical(ly inspired) guess
    else:
        rMax = get_rmax_bulk(eventProgenitorMass)
    if (eventProgenitor == "gluino"): return rMax
    return (10*rMax)

templateReader = MCTemplateReader.MCTemplateReader(inputArguments.MCTemplatePath)
for indexPair in templateReader.nextValidBin():
    eventProgenitorMassBin = indexPair[0]
    eventProgenitorMass = (templateReader.eventProgenitorMasses)[eventProgenitorMassBin]
    neutralinoMassBin = indexPair[1]
    neutralinoMass = (templateReader.neutralinoMasses)[neutralinoMassBin]
    initial_rMax = 10.0*get_rmax(eventProgenitorMass, neutralinoMass, inputArguments.eventProgenitor) # with safety margin
    print("eventProgenitor mass: {gM}, neutralino mass: {nM}".format(gM=eventProgenitorMass, nM=neutralinoMass))
    # if ((inputArguments.dataCardsDirectory)[0] == "/"): # inputArguments.dataCardsDirectory is likely an absolute path, not a relative path
    #     dataCardPathsPrefix = "{dCD}/{dCP}_dataCard_eventProgenitorMassBin{gMB}_neutralinoMassBin{nMB}".format(dCD=inputArguments.dataCardsDirectory, dCP=inputArguments.dataCardsPrefix, gMB=eventProgenitorMassBin, nMB=neutralinoMassBin)
    # else:
    #     dataCardPathsPrefix = "{sR}/{dCD}/{dCP}_dataCard_eventProgenitorMassBin{gMB}_neutralinoMassBin{nMB}".format(sR=stealthEnv.stealthRoot, dCD=inputArguments.dataCardsDirectory, dCP=inputArguments.dataCardsPrefix, gMB=eventProgenitorMassBin, nMB=neutralinoMassBin)
    x509ProxyPath = stealthEnv.x509Proxy
    tmUtilsTarballPath = "{tUP}/tmUtils.tar.gz".format(tUP=stealthEnv.tmUtilsParent)
    tmUtilsExtractionScriptPath = "{tUP}/extract_tmUtilsTarball.sh".format(tUP=stealthEnv.tmUtilsParent)
    remoteEnvSetupScriptPath = "{sR}/setup_environment_remote.sh".format(sR=stealthEnv.stealthRoot)
    MCTemplateReaderPath = "{sR}/MCTemplateReader.py".format(sR=stealthEnv.stealthRoot)
    crossSectionsFilePath = "{sR}/{cSFN}".format(sR=stealthEnv.stealthRoot, cSFN=inputArguments.crossSectionsFileName)
    STRegionBoundariesFilePath = "{sR}/STRegionBoundaries.dat".format(sR=stealthEnv.stealthRoot)
    dataSystematicsPath_signal = inputArguments.path_dataSystematics_signal
    dataSystematicsPath_signal_loose = inputArguments.path_dataSystematics_signal_loose
    dataSystematicsPath_control = inputArguments.path_dataSystematics_control
    dataObservedEventCountersPath_signal = inputArguments.path_dataObservedEventCounters_signal
    dataObservedEventCountersPath_signal_loose = inputArguments.path_dataObservedEventCounters_signal_loose
    dataObservedEventCountersPath_control = inputArguments.path_dataObservedEventCounters_control
    dataExpectedEventCountersPath_signal = inputArguments.path_dataExpectedEventCounters_signal
    dataExpectedEventCountersPath_signal_loose = inputArguments.path_dataExpectedEventCounters_signal_loose
    dataExpectedEventCountersPath_control = inputArguments.path_dataExpectedEventCounters_control
    createDataCardScriptPath = "{sR}/createDataCard.py".format(sR=stealthEnv.stealthRoot)
    commonPyFunctionsFilePath = "{sR}/commonFunctions.py".format(sR=stealthEnv.stealthRoot)
    readBestFitScriptPath = "{sR}/readBestFitFromMultiDimOutput.py".format(sR=stealthEnv.stealthRoot)
    scalingSystematicsScriptPath = "{sR}/getSTScalingSystematics.py".format(sR=stealthEnv.stealthRoot)
    limitsConvergenceCheckScriptPath = "{sR}/checkLimitsConvergence.py".format(sR=stealthEnv.stealthRoot)

    filesToTransfer = [x509ProxyPath, tmUtilsTarballPath, tmUtilsExtractionScriptPath, remoteEnvSetupScriptPath, MCTemplateReaderPath, crossSectionsFilePath, STRegionBoundariesFilePath, dataSystematicsPath_signal, dataSystematicsPath_control, dataSystematicsPath_signal, dataObservedEventCountersPath_signal, dataObservedEventCountersPath_control, dataExpectedEventCountersPath_signal, dataExpectedEventCountersPath_control, createDataCardScriptPath, commonPyFunctionsFilePath, readBestFitScriptPath, scalingSystematicsScriptPath, limitsConvergenceCheckScriptPath]
    if (inputArguments.addLooseSignal): filesToTransfer.extend([dataSystematicsPath_signal_loose, dataObservedEventCountersPath_signal_loose, dataExpectedEventCountersPath_signal_loose])
    processIdentifier = "combineJob_{prefix}_eventProgenitorMassBin{gMB}_neutralinoMassBin{nMB}".format(prefix=inputArguments.dataCardsPrefix, gMB=eventProgenitorMassBin, nMB=neutralinoMassBin)
    jdlInterface = tmJDLInterface.tmJDLInterface(processName=processIdentifier, scriptPath="combineToolHelper.sh", outputDirectoryRelativePath="{cWAR}/combine{oI}".format(cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier))  # works even if "outputDirectoryRelativePath" is an absolute path
    jdlInterface.addFilesToTransferFromList(filesToTransfer)
    # Arguments for script:
    jdlInterface.addScriptArgument("{oD}".format(oD=inputArguments.outputDirectory)) # Argument 1: output directory
    jdlInterface.addScriptArgument("{dCP}".format(dCP=inputArguments.dataCardsPrefix)) # Argument 2: data cards prefix
    jdlInterface.addScriptArgument("{gMB}".format(gMB=eventProgenitorMassBin)) # Argument 3: eventProgenitor mass bin index
    jdlInterface.addScriptArgument("{nMB}".format(nMB=neutralinoMassBin)) # Argument 4: neutralino mass bin index
    jdlInterface.addScriptArgument("{irM:.1f}".format(irM=initial_rMax)) # Argument 5: initial rMax
    jdlInterface.addScriptArgument("{cSFN}".format(cSFN=inputArguments.crossSectionsFileName)) # Argument 6: cross-sections file name
    jdlInterface.addScriptArgument("{MCTP}".format(MCTP=inputArguments.MCTemplatePath)) # Argument 7: MC template path
    jdlInterface.addScriptArgument("{MCHS}".format(MCHS=inputArguments.MCHistogramsSignal)) # Argument 8: path to MC event histograms, signal
    jdlInterface.addScriptArgument("{MCHSl}".format(MCHSl=inputArguments.MCHistogramsSignalLoose)) # Argument 9: path to MC event histograms, loose signal
    jdlInterface.addScriptArgument("{MCHC}".format(MCHC=inputArguments.MCHistogramsControl)) # Argument 10: path to MC event histograms, control
    jdlInterface.addScriptArgument("{MCUS}".format(MCUS=inputArguments.MCUncertaintiesSignal)) # Argument 11: path to MC uncertainties, signal
    jdlInterface.addScriptArgument("{MCUSl}".format(MCUSl=inputArguments.MCUncertaintiesSignalLoose)) # Argument 12: path to MC uncertainties, loose signal
    jdlInterface.addScriptArgument("{MCUC}".format(MCUC=inputArguments.MCUncertaintiesControl)) # Argument 13: path to MC uncertainties, control
    jdlInterface.addScriptArgument("{LU:.4f}".format(LU=inputArguments.luminosityUncertainty)) # Argument 14: lumi uncertainty
    jdlInterface.addScriptArgument("{EAA}".format(EAA=inputArguments.EOSAnalysisArea)) # Argument 15: EOS analysis area path
    runUnblindedString = "false"
    if (inputArguments.runUnblinded): runUnblindedString = "true"
    jdlInterface.addScriptArgument("{aLSS}".format(aLSS=runUnblindedString)) # Argument 16: run unblinded switch
    addLooseSignalString = "false"
    if (inputArguments.addLooseSignal): addLooseSignalString = "true"
    jdlInterface.addScriptArgument("{aLSS}".format(aLSS=addLooseSignalString)) # Argument 17: add loose signal switch
    # Write JDL
    jdlInterface.writeToFile()
    submissionCommand = "cd {cWAR}/combine{oI}/ && condor_submit {pI}.jdl && cd {sR}".format(pI=processIdentifier, cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier, sR=stealthEnv.stealthRoot)
    print ("Generated command: {sC}".format(sC=submissionCommand))
    if (inputArguments.isDryRun):
        print("Not submitting due to dryRun flag.")
    else:
        os.system(submissionCommand)
        print ("Submitted.")
