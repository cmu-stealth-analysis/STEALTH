#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, pdb, tmGeneralUtils, tmCombineDataCardInterface, os, sys, MCTemplateReader

inputArgumentsParser = argparse.ArgumentParser(description='Create data card from MC and data systematics and nEvents data.')
inputArgumentsParser.add_argument('--outputPrefix', required=True, help='Prefix to output files.', type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="analysis/dataCards/", help='Output directory.', type=str)
inputArgumentsParser.add_argument('--eventProgenitorMassBin', required=True, help='Event progenitor mass bin.', type=int)
inputArgumentsParser.add_argument('--neutralinoMassBin', required=True, help='Neutralino mass bin.', type=int)
inputArgumentsParser.add_argument('--crossSectionsFile', required=True, help='Path to dat file that contains cross-sections as a function of eventProgenitor mass, from which to get the fractional uncertainties.',type=str)
inputArgumentsParser.add_argument('--crossSectionsScale', required=True, help='Cross sections shift scale. Accepted values: +/- 1 or 0. Signal is scaled by the factor cross_section*(1+this_scale*cross_section_uncertainty)', type=int)
inputArgumentsParser.add_argument('--MCTemplatePath', default="plot_susyMasses_template.root", help='Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.', type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--inputFile_MCEventHistograms_signal', required=True, help='Input MC event histograms, signal.', type=str)
inputArgumentsParser.add_argument('--inputFile_MCEventHistograms_signal_loose', required=True, help='Input MC event histograms, loose signal.', type=str)
inputArgumentsParser.add_argument('--inputFile_MCEventHistograms_control', required=True, help='Input MC event histograms, control.', type=str)
inputArgumentsParser.add_argument('--inputFile_MCUncertainties_signal', required=True, help='Input MC uncertainties, signal.', type=str)
inputArgumentsParser.add_argument('--inputFile_MCUncertainties_signal_loose', required=True, help='Input MC uncertainties, loose signal.', type=str)
inputArgumentsParser.add_argument('--inputFile_MCUncertainties_control', required=True, help='Input MC uncertainties, control.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics_signal', required=True, help='Input file containing fractional uncertainties from signal data.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics_signal_loose', required=True, help='Input file containing fractional uncertainties from loose signal data.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics_control', required=True, help='Input file containing fractional uncertainties from control data.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics_expectedEventCounters_signal', required=True, help='Input file containing expected number of events from signal data.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics_expectedEventCounters_signal_loose', required=True, help='Input file containing expected number of events from loose signal data.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics_expectedEventCounters_control', required=True, help='Input file containing expected number of events from control data.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics_observedEventCounters_signal', required=True, help='Input file containing observed number of events from signal data.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics_observedEventCounters_signal_loose', required=True, help='Input file containing observed number of events from loose signal data.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics_observedEventCounters_control', required=True, help='Input file containing observed number of events from control data.', type=str)
inputArgumentsParser.add_argument('--luminosityUncertainty', required=True, help='Uncertainty on the luminosity.', type=float)
inputArgumentsParser.add_argument('--runUnblinded', action='store_true', help="If this flag is set, then the signal region data is unblinded. Specifically, the entry for the observed number of events is filled from the data, rather than from the expectation values.")
inputArgumentsParser.add_argument('--regionsToUse', required=True, help="Comma-separated list of regions to run on.", type=str)
inputArguments = inputArgumentsParser.parse_args()

SYSTEMATIC_SIGNIFICANCE_THRESHOLD = 0.0005
LOGNORMAL_REGULARIZE_THRESHOLD = 0.01 # To "regularize" the lognormal error; the error in the column is supposed to be (1 plus deltaX/X), and if this is 0, the combine tool complains

if not((inputArguments.crossSectionsScale == 0) or (abs(inputArguments.crossSectionsScale) == 1)): sys.exit("ERROR: argument crossSectionsScale can only take values 0 or +/- 1. Currently, crossSectionsScale={cSS}".format(cSS=inputArguments.crossSectionsScale))

def get_dict_nEvents(inputPath, localSignalLabels, inputPrefix):
    fileContents_nEvents = tmGeneralUtils.getConfigurationFromFile(inputFilePath=inputPath)
    outputDict = {}
    for signalBinLabel in localSignalLabels:
        outputDict[signalBinLabel] = fileContents_nEvents["{iP}_{l}".format(iP=inputPrefix, l=signalBinLabel)]
    return outputDict

def get_dict_expectedNEvents_stealth(stealthNEventsHistograms, eventProgenitorMassBin, neutralinoMassBin, localSignalLabels, scaleFactor):
    outputDict = {}
    for signalBinLabel in localSignalLabels:
        outputDict[signalBinLabel] = scaleFactor*((stealthNEventsHistograms[signalBinLabel]).GetBinContent(eventProgenitorMassBin, neutralinoMassBin))
        if (outputDict[signalBinLabel] < 0.001): outputDict[signalBinLabel] = 0.001 # To avoid possible pathologies in the combine algorithm
    return outputDict

def get_symmetric_data_systematics_from_file(localSignalLabels, dataSystematicLabels, sourceFile):
    sourceSystematics = tmGeneralUtils.getConfigurationFromFile(sourceFile)
    outputDict = {}
    for dataSystematicLabel in dataSystematicLabels:
        outputDict[dataSystematicLabel] = {}
        for signalBinLabel in localSignalLabels:
            outputDict[dataSystematicLabel][signalBinLabel] = 1.0 + sourceSystematics["fractionalUncertainty_{dSL}_{sBL}".format(dSL=dataSystematicLabel, sBL=signalBinLabel)]
    return outputDict

def get_asymmetric_data_systematics_from_file(localSignalLabels, dataSystematicLabels, sourceFile):
    sourceSystematics = tmGeneralUtils.getConfigurationFromFile(sourceFile)
    outputDict = {}
    for dataSystematicLabel in dataSystematicLabels:
        outputDict[dataSystematicLabel] = {}
        for signalBinLabel in localSignalLabels:
            outputDict[dataSystematicLabel][signalBinLabel] = {}
            for upDownLabel in ["Down", "Up"]:
                outputDict[dataSystematicLabel][signalBinLabel][upDownLabel] = 1.0 + sourceSystematics["fractionalUncertainty{uDL}_{dSL}_{sBL}".format(uDL=upDownLabel, dSL=dataSystematicLabel, sBL=signalBinLabel)]
    return outputDict

def build_data_systematic_with_check(list_signalTypes, dict_localToGlobalBinLabels, dict_localSignalLabelsToUse, dict_sources_dataSystematics):
    outputDict = {}
    isSignificant = False
    for signalType in list_signalTypes:
        if (not(signalType in dict_localSignalLabelsToUse)): continue
        if (len(dict_localSignalLabelsToUse[signalType]) == 0): continue
        sourceDict_dataSystematics = dict_sources_dataSystematics[signalType]
        for localSignalBinLabel in dict_localSignalLabelsToUse[signalType]:
            globalSignalBinLabel = dict_localToGlobalBinLabels[signalType][localSignalBinLabel]
            outputDict[globalSignalBinLabel] = {}
            if isinstance(sourceDict_dataSystematics[localSignalBinLabel], dict):
                outputDict[globalSignalBinLabel]["qcd"] = {}
                outputDict[globalSignalBinLabel]["scalingDeviation"] = {}
                for upDownLabel in ["Down", "Up"]:
                    try:
                        for backgroundProcessLabel in ["qcd", "scalingDeviation"]:
                            outputDict[globalSignalBinLabel][backgroundProcessLabel][upDownLabel] = sourceDict_dataSystematics[localSignalBinLabel][upDownLabel]
                            if (abs(outputDict[globalSignalBinLabel][backgroundProcessLabel][upDownLabel] - 1.0) > SYSTEMATIC_SIGNIFICANCE_THRESHOLD): isSignificant = True
                            if (abs(outputDict[globalSignalBinLabel][backgroundProcessLabel][upDownLabel]) < LOGNORMAL_REGULARIZE_THRESHOLD): outputDict[globalSignalBinLabel][backgroundProcessLabel][upDownLabel] = LOGNORMAL_REGULARIZE_THRESHOLD
                    except KeyError:
                        sys.exit("ERROR: sourceDict_dataSystematics[localSignalBinLabel] is a dict but does not have elements named \"Up\" or \"Down\". dict contents: {c}".format(c=sourceDict_dataSystematics[localSignalBinLabel]))
            else:
                for backgroundProcessLabel in ["qcd", "scalingDeviation"]:
                    outputDict[globalSignalBinLabel][backgroundProcessLabel] = sourceDict_dataSystematics[localSignalBinLabel]
                    if (abs(outputDict[globalSignalBinLabel][backgroundProcessLabel] - 1.0) > SYSTEMATIC_SIGNIFICANCE_THRESHOLD): isSignificant = True
                    if (abs(outputDict[globalSignalBinLabel][backgroundProcessLabel]) < LOGNORMAL_REGULARIZE_THRESHOLD): outputDict[globalSignalBinLabel][backgroundProcessLabel] = LOGNORMAL_REGULARIZE_THRESHOLD
    return (isSignificant, outputDict)

def build_MC_constant_systematic(list_signalTypes, dict_localToGlobalBinLabels, constantFractionalUncertainty):
    outputDict = {}
    for signalType in list_signalTypes:
        for localSignalBinLabel in (dict_localToGlobalBinLabels[signalType].keys()):
            globalSignalBinLabel = dict_localToGlobalBinLabels[signalType][localSignalBinLabel]
            outputDict[globalSignalBinLabel] = {}
            outputDict[globalSignalBinLabel]["stealth"] = 1.0 + constantFractionalUncertainty
    return outputDict

def get_MC_systematic_from_histogram(localSignalLabels, inputHistograms, eventProgenitorMassBin, neutralinoMassBin):
    outputDict = {}
    for signalBinLabel in localSignalLabels:
        outputDict[signalBinLabel] = 1.0 + (inputHistograms[signalBinLabel]).GetBinContent(eventProgenitorMassBin, neutralinoMassBin)
    return outputDict

def get_asymmetric_MC_systematic_from_histogram(localSignalLabels, inputHistograms, eventProgenitorMassBin, neutralinoMassBin):
    outputDict = {}
    for signalBinLabel in localSignalLabels:
        outputDict[signalBinLabel] = {}
        for UpDownShift in ["Down", "Up"]:
            outputDict[signalBinLabel][UpDownShift] = 1.0 + (inputHistograms[signalBinLabel][UpDownShift]).GetBinContent(eventProgenitorMassBin, neutralinoMassBin)
    return outputDict

def build_MC_systematic_with_check(list_signalTypes, dict_localToGlobalBinLabels, dict_localSignalLabelsToUse, dict_sources_dataSystematics):
    outputDict = {}
    isSignificant = False
    for signalType in list_signalTypes:
        if (not(signalType in dict_localSignalLabelsToUse)): continue
        if (len(dict_localSignalLabelsToUse[signalType]) == 0): continue
        sourceDict_MCSystematics = dict_sources_dataSystematics[signalType]
        for localSignalBinLabel in dict_localSignalLabelsToUse[signalType]:
            globalSignalBinLabel = dict_localToGlobalBinLabels[signalType][localSignalBinLabel]
            outputDict[globalSignalBinLabel] = {}
            if isinstance(sourceDict_MCSystematics[localSignalBinLabel], dict):
                outputDict[globalSignalBinLabel]["stealth"] = {}
                for upDownLabel in ["Down", "Up"]:
                    try:
                        outputDict[globalSignalBinLabel]["stealth"][upDownLabel] = sourceDict_MCSystematics[localSignalBinLabel][upDownLabel]
                        if (abs(outputDict[globalSignalBinLabel]["stealth"][upDownLabel] - 1.0) > SYSTEMATIC_SIGNIFICANCE_THRESHOLD): isSignificant = True
                        if (abs(outputDict[globalSignalBinLabel]["stealth"][upDownLabel]) < LOGNORMAL_REGULARIZE_THRESHOLD): outputDict[globalSignalBinLabel]["stealth"][upDownLabel] = LOGNORMAL_REGULARIZE_THRESHOLD
                    except KeyError:
                        sys.exit("ERROR: sourceDict_MCSystematics[localSignalBinLabel] is a dict but does not have elements named \"Up\" or \"Down\". dict contents: {c}".format(c=sourceDict_MCSystematics[localSignalBinLabel]))
            else:
                outputDict[globalSignalBinLabel]["stealth"] = sourceDict_MCSystematics[localSignalBinLabel]
                if (abs(outputDict[globalSignalBinLabel]["stealth"] - 1.0) > SYSTEMATIC_SIGNIFICANCE_THRESHOLD): isSignificant = True
                if (abs(outputDict[globalSignalBinLabel]["stealth"]) < LOGNORMAL_REGULARIZE_THRESHOLD): outputDict[globalSignalBinLabel]["stealth"] = LOGNORMAL_REGULARIZE_THRESHOLD
    return (isSignificant, outputDict)

def createDataCard(outputPath,
                   signalBinLabels,
                   observedNEvents,
                   expectedNEvents_qcd,
                   expectedRate_scalingDeviation,
                   systematics_data,
                   systematics_data_labels,
                   systematics_data_types,
                   expectedNEvents_stealth,
                   systematics_MC,
                   systematics_MC_labels,
                   systematics_MC_types,
                   rateParamLabels,
                   rateParamProperties
):
    expectedNEvents = {}
    for signalBinLabel in signalBinLabels:
        expectedNEvents[signalBinLabel] = {}
        expectedNEvents[signalBinLabel]["qcd"] = expectedNEvents_qcd[signalBinLabel]
        expectedNEvents[signalBinLabel]["scalingDeviation"] = expectedRate_scalingDeviation[signalBinLabel]
        expectedNEvents[signalBinLabel]["stealth"] = expectedNEvents_stealth[signalBinLabel]

    systematics = {}
    systematicsLabels = []
    systematicsTypes = {}
    for dataSystematicLabel in systematics_data_labels:
        systematics[dataSystematicLabel] = systematics_data[dataSystematicLabel]
        systematicsLabels.append(dataSystematicLabel)
        systematicsTypes[dataSystematicLabel] = systematics_data_types[dataSystematicLabel]
    for MCSystematicLabel in systematics_MC_labels:
        systematics[MCSystematicLabel] = systematics_MC[MCSystematicLabel]
        systematicsLabels.append(MCSystematicLabel)
        systematicsTypes[MCSystematicLabel] = systematics_MC_types[MCSystematicLabel]

    combineInterface = tmCombineDataCardInterface.tmCombineDataCardInterface(list_signalBinLabels=signalBinLabels, list_backgroundProcessLabels=["qcd", "scalingDeviation"], list_signalProcessLabels=["stealth"], list_systematicsLabels=systematicsLabels, list_rateParamLabels=rateParamLabels, dict_rateParamProperties=rateParamProperties, dict_observedNEvents=observedNEvents, dict_expectedNEvents=expectedNEvents, dict_systematicsTypes=systematicsTypes, dict_systematics=systematics)
    combineInterface.writeToFile(outputFilePath=outputPath)

crossSectionsInputFileObject = open(inputArguments.crossSectionsFile, 'r')
crossSectionsDictionary = {}
crossSectionsFractionalUncertaintyDictionary = {}
for line in crossSectionsInputFileObject:
    crossSectionsData = line.split()
    eventProgenitorMassInt = int(0.5 + float(crossSectionsData[0]))
    crossSection = float(crossSectionsData[1])
    crossSectionFractionalUncertainty = 0.01*float(crossSectionsData[2])
    crossSectionsDictionary[eventProgenitorMassInt] = crossSection
    crossSectionsFractionalUncertaintyDictionary[eventProgenitorMassInt] = crossSectionFractionalUncertainty
crossSectionsInputFileObject.close()

STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries, 'r')
nSTBoundaries = 0
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()): nSTBoundaries += 1
nSTSignalBins = nSTBoundaries - 2 + 1 # First two lines are for the normalization bin, last boundary implied at infinity
print("Using {n} signal bins for ST.".format(n = nSTSignalBins))
STRegionBoundariesFileObject.close()

# signalTypesToUse = ["signal", "control"]
signalTypesToUse = []
for signalType in ((inputArguments.regionsToUse).strip()).split(","):
    if not(signalType in ["signal", "signal_loose", "control"]): sys.exit("ERROR: Unrecognized region to use: {r}".format(r=signalType))
    signalTypesToUse.append(signalType)

abbreviated_signalTypes = {
    "signal": "s",
    "signal_loose": "l",
    "control": "c"
}

inputDataFilePaths = {
    "signal": {
        "observations": inputArguments.inputFile_dataSystematics_observedEventCounters_signal,
        "expectations": inputArguments.inputFile_dataSystematics_expectedEventCounters_signal
    },
    "signal_loose": {
        "observations": inputArguments.inputFile_dataSystematics_observedEventCounters_signal_loose,
        "expectations": inputArguments.inputFile_dataSystematics_expectedEventCounters_signal_loose
    },
    "control": {
        "observations": inputArguments.inputFile_dataSystematics_observedEventCounters_control,
        "expectations": inputArguments.inputFile_dataSystematics_expectedEventCounters_control
    }
}

inputDataSystematicsFilePaths = {
    "signal": inputArguments.inputFile_dataSystematics_signal,
    "signal_loose": inputArguments.inputFile_dataSystematics_signal_loose,
    "control": inputArguments.inputFile_dataSystematics_control
    # "scaling": inputArguments.inputFile_dataSystematics_sTScaling
}

inputMCFilePaths = {
    "signal": {
        "eventHistograms": inputArguments.inputFile_MCEventHistograms_signal,
        "uncertainties": inputArguments.inputFile_MCUncertainties_signal
    },
    "signal_loose": {
        "eventHistograms": inputArguments.inputFile_MCEventHistograms_signal_loose,
        "uncertainties": inputArguments.inputFile_MCUncertainties_signal_loose
    },
    "control": {
        "eventHistograms": inputArguments.inputFile_MCEventHistograms_control,
        "uncertainties": inputArguments.inputFile_MCUncertainties_control
    }
}

localSignalBinLabels = []
globalSignalBinLabels = []
dict_localToGlobalBinLabels = {}
for signalType in signalTypesToUse:
    dict_localToGlobalBinLabels[signalType] = {}

for STRegionIndex in range(2, 2 + nSTSignalBins):
    for nJetsBin in range(4, 7):
        localSignalBinLabel = "STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin)
        localSignalBinLabels.append(localSignalBinLabel)
        for signalType in signalTypesToUse:
            dict_localToGlobalBinLabels[signalType][localSignalBinLabel] = (abbreviated_signalTypes[signalType] + "ST{r}J{n}".format(r=STRegionIndex, n=nJetsBin))
            globalSignalBinLabels.append(dict_localToGlobalBinLabels[signalType][localSignalBinLabel])

expectedNEvents_qcd = {}
expectedRate_scalingDeviation = {}
for signalType in signalTypesToUse:
    expectedNEventsLocal_qcd = get_dict_nEvents(inputPath=inputDataFilePaths[signalType]["expectations"], localSignalLabels=localSignalBinLabels, inputPrefix="expectedNEvents")
    for localLabel in localSignalBinLabels:
        globalLabel = dict_localToGlobalBinLabels[signalType][localLabel]
        expectedNEvents_qcd[globalLabel] = expectedNEventsLocal_qcd[localLabel]
        if (signalType == "control"): expectedRate_scalingDeviation[globalLabel] = 1

for signalType in signalTypesToUse:
    if (signalType == "control"): continue
    expectedNEventsLocal_qcd = get_dict_nEvents(inputPath=inputDataFilePaths[signalType]["expectations"], localSignalLabels=localSignalBinLabels, inputPrefix="expectedNEvents")
    expectedNEventsLocal_qcd_control = get_dict_nEvents(inputPath=inputDataFilePaths["control"]["expectations"], localSignalLabels=localSignalBinLabels, inputPrefix="expectedNEvents")
    for localLabel in localSignalBinLabels:
        globalLabel = dict_localToGlobalBinLabels[signalType][localLabel]
        expectedRate_scalingDeviation[globalLabel] = expectedNEventsLocal_qcd[localLabel]/expectedNEventsLocal_qcd_control[localLabel]

observedNEvents = {}
initialDeviations_STScaling = get_dict_nEvents(inputPath=inputDataFilePaths["control"]["expectations"], localSignalLabels=localSignalBinLabels, inputPrefix="scalingDeviation")
if (inputArguments.runUnblinded):
    for signalType in signalTypesToUse:
        observedNEventsLocal = get_dict_nEvents(inputPath=inputDataFilePaths[signalType]["observations"], localSignalLabels=localSignalBinLabels, inputPrefix="observedNEvents")
        for localLabel in localSignalBinLabels:
            globalLabel = dict_localToGlobalBinLabels[signalType][localLabel]
            observedNEvents[globalLabel] = observedNEventsLocal[localLabel]
else: # force Asimov dataset
    for signalType in signalTypesToUse:
        for localLabel in localSignalBinLabels:
            globalLabel = dict_localToGlobalBinLabels[signalType][localLabel]
            observedNEvents[globalLabel] = expectedNEvents_qcd[globalLabel] + expectedRate_scalingDeviation[globalLabel]*initialDeviations_STScaling[localLabel]

# rate params for ST scaling
rateParamLabels = []
rateParamProperties = {}
for localSignalBinLabel in localSignalBinLabels:
    rateParamLabels.append("rateParam_{lSBL}".format(lSBL=localSignalBinLabel))
    list_globalLabels_rateParams = []
    minPhysicallyPermissibleDeviation = -1.0
    for signalType in signalTypesToUse:
        globalSignalBinLabel = dict_localToGlobalBinLabels[signalType][localSignalBinLabel]
        list_globalLabels_rateParams.append(globalSignalBinLabel)
        physicallyPermissibleDeviation = expectedNEvents_qcd[globalSignalBinLabel]/expectedRate_scalingDeviation[globalSignalBinLabel]
        if ((minPhysicallyPermissibleDeviation < 0.0) or (physicallyPermissibleDeviation < minPhysicallyPermissibleDeviation)): minPhysicallyPermissibleDeviation = physicallyPermissibleDeviation
    rateParamProperties["rateParam_{lSBL}".format(lSBL=localSignalBinLabel)] = (list_globalLabels_rateParams, ["scalingDeviation"], initialDeviations_STScaling[localSignalBinLabel], -1.0*minPhysicallyPermissibleDeviation, 9.0*minPhysicallyPermissibleDeviation)
    if ((initialDeviations_STScaling[localSignalBinLabel] < -1.0*minPhysicallyPermissibleDeviation) or (initialDeviations_STScaling[localSignalBinLabel] > 9.0*minPhysicallyPermissibleDeviation)): sys.exit("Sanity check problem: initialDeviations_STScaling[localSignalBinLabel] = {a}, minPhysicallyPermissibleDeviation = {b}".format(a=initialDeviations_STScaling[localSignalBinLabel], b=minPhysicallyPermissibleDeviation))

# Data systematics
systematics_data = {}
systematics_data_labels = []
systematics_data_types = {}

for signalType in signalTypesToUse:
    sources_symmetricDataSystematicsLocal = get_symmetric_data_systematics_from_file(localSignalLabels=localSignalBinLabels, dataSystematicLabels=["shape", "rho"], sourceFile=inputDataSystematicsFilePaths[signalType])
    sources_asymmetricDataSystematicsLocal = get_asymmetric_data_systematics_from_file(localSignalLabels=localSignalBinLabels, dataSystematicLabels=["normEvents"], sourceFile=inputDataSystematicsFilePaths[signalType])
    sources_dataSystematics = {}
    for dataSystematic in ["shape", "rho"]:
        sources_dataSystematics[dataSystematic] = {signalType: sources_symmetricDataSystematicsLocal[dataSystematic]}
    for dataSystematic in ["normEvents"]:
        sources_dataSystematics[dataSystematic] = {signalType: sources_asymmetricDataSystematicsLocal[dataSystematic]}
    for nJetsBin in range(4, 7):
        localLabelsToUse = {signalType: []}
        for STRegionIndex in range(2, 2+nSTSignalBins):
            signalBinLabel = "STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin)
            localLabelsToUse[signalType].append(signalBinLabel)
        sources_normSystematics = {signalType: sources_asymmetricDataSystematicsLocal["normEvents"]}
        tmp = build_data_systematic_with_check(list_signalTypes=[signalType], dict_localToGlobalBinLabels=dict_localToGlobalBinLabels, dict_localSignalLabelsToUse=localLabelsToUse, dict_sources_dataSystematics=sources_dataSystematics["normEvents"])
        if (tmp[0]):
            systematicsLabel = "normEvents_{n}Jets_{sT}".format(n=nJetsBin, sT=signalType)
            systematics_data_labels.append(systematicsLabel)
            systematics_data_types[systematicsLabel] = "lnN"
            systematics_data[systematicsLabel] = tmp[1]
    # shape and rho systematics are correlated across all nJets bins
    for dataSystematic in ["shape", "rho"]:
        for STRegionIndex in range(2, 2+nSTSignalBins):
            localLabelsToUse = {signalType: []}
            for nJetsBin in range(4, 7):
                signalBinLabel = "STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin)
                localLabelsToUse[signalType].append(signalBinLabel)
            tmp = build_data_systematic_with_check(list_signalTypes=[signalType], dict_localToGlobalBinLabels=dict_localToGlobalBinLabels, dict_localSignalLabelsToUse=localLabelsToUse, dict_sources_dataSystematics=sources_dataSystematics[dataSystematic])
            if (tmp[0]):
                systematicsLabel = "{dS}_STRegion{r}_{sT}".format(dS=dataSystematic, r=STRegionIndex, sT=signalType)
                systematics_data_labels.append(systematicsLabel)
                systematics_data_types[systematicsLabel] = "lnN"
                systematics_data[systematicsLabel] = tmp[1]

MCHistograms_weightedNEvents = {}
MCHistograms_MCStatUncertainties = {}
MCHistograms_JECUncertainties = {}
MCHistograms_UnclusteredMETUncertainties = {}
MCHistograms_JERMETUncertainties = {}
MCHistograms_prefiringWeightsUncertainties = {}
MCHistograms_HLTUncertainties = {}
MCHistograms_photonScaleFactorUncertainties = {}
for signalType in signalTypesToUse:
    MCHistograms_weightedNEvents[signalType] = {}
    MCHistograms_MCStatUncertainties[signalType] = {}
    MCHistograms_JECUncertainties[signalType] = {}
    MCHistograms_UnclusteredMETUncertainties[signalType] = {}
    MCHistograms_JERMETUncertainties[signalType] = {}
    MCHistograms_prefiringWeightsUncertainties[signalType] = {}
    MCHistograms_HLTUncertainties[signalType] = {}
    MCHistograms_photonScaleFactorUncertainties[signalType] = {}

MCEventHistograms = {}
MCUncertainties = {}
expectedNEvents_stealth = None
for signalType in signalTypesToUse:
    MCEventHistograms[signalType] = ROOT.TFile.Open(inputMCFilePaths[signalType]["eventHistograms"], "READ")
    MCUncertainties[signalType] = ROOT.TFile.Open(inputMCFilePaths[signalType]["uncertainties"], "READ")
    for signalBinLabel in localSignalBinLabels:
        MCHistograms_weightedNEvents[signalType][signalBinLabel] = ROOT.TH2F()
        MCEventHistograms[signalType].GetObject("h_lumiBasedYearWeightedNEvents_{sBL}".format(sBL=signalBinLabel), MCHistograms_weightedNEvents[signalType][signalBinLabel])
        if (not(MCHistograms_weightedNEvents[signalType][signalBinLabel])):
            sys.exit("ERROR: Histogram MCHistograms_weightedNEvents[signalType][signalBinLabel] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
        MCHistograms_MCStatUncertainties[signalType][signalBinLabel] = {}
        MCHistograms_JECUncertainties[signalType][signalBinLabel] = {}
        MCHistograms_UnclusteredMETUncertainties[signalType][signalBinLabel] = {}
        MCHistograms_JERMETUncertainties[signalType][signalBinLabel] = {}
        MCHistograms_prefiringWeightsUncertainties[signalType][signalBinLabel] = {}
        MCHistograms_HLTUncertainties[signalType][signalBinLabel] = {}
        MCHistograms_photonScaleFactorUncertainties[signalType][signalBinLabel] = {}
        for UpDownShift in ["Down", "Up"]:
            MCHistograms_MCStatUncertainties[signalType][signalBinLabel][UpDownShift] = ROOT.TH2F()
            MCUncertainties[signalType].GetObject("h_MCStatisticsFractionalError{UDS}_{sBL}".format(UDS=UpDownShift, sBL=signalBinLabel), MCHistograms_MCStatUncertainties[signalType][signalBinLabel][UpDownShift])
            if (not(MCHistograms_MCStatUncertainties[signalType][signalBinLabel][UpDownShift])):
                sys.exit("ERROR: Histogram MCHistograms_MCStatUncertainties[signalType][signalBinLabel][UpDownShift] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
            MCHistograms_JECUncertainties[signalType][signalBinLabel][UpDownShift] = ROOT.TH2F()
            MCUncertainties[signalType].GetObject("h_JECUncertainty{UDS}_{sBL}".format(UDS=UpDownShift, sBL=signalBinLabel), MCHistograms_JECUncertainties[signalType][signalBinLabel][UpDownShift])
            if (not(MCHistograms_JECUncertainties[signalType][signalBinLabel][UpDownShift])):
                sys.exit("ERROR: Histogram MCHistograms_JECUncertainties[signalType][signalBinLabel][UpDownShift] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
            MCHistograms_UnclusteredMETUncertainties[signalType][signalBinLabel][UpDownShift] = ROOT.TH2F()
            MCUncertainties[signalType].GetObject("h_UnclusteredMETUncertainty{UDS}_{sBL}".format(UDS=UpDownShift, sBL=signalBinLabel), MCHistograms_UnclusteredMETUncertainties[signalType][signalBinLabel][UpDownShift])
            if (not(MCHistograms_UnclusteredMETUncertainties[signalType][signalBinLabel][UpDownShift])):
                sys.exit("ERROR: Histogram MCHistograms_UnclusteredMETUncertainties[signalType][signalBinLabel][UpDownShift] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
            MCHistograms_JERMETUncertainties[signalType][signalBinLabel][UpDownShift] = ROOT.TH2F()
            MCUncertainties[signalType].GetObject("h_JERMETUncertainty{UDS}_{sBL}".format(UDS=UpDownShift, sBL=signalBinLabel), MCHistograms_JERMETUncertainties[signalType][signalBinLabel][UpDownShift])
            if (not(MCHistograms_JERMETUncertainties[signalType][signalBinLabel][UpDownShift])):
                sys.exit("ERROR: Histogram MCHistograms_JERMETUncertainties[signalType][signalBinLabel][UpDownShift] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
            MCHistograms_prefiringWeightsUncertainties[signalType][signalBinLabel][UpDownShift] = ROOT.TH2F()
            MCUncertainties[signalType].GetObject("h_prefiringWeightsUncertainty{UDS}_{sBL}".format(UDS=UpDownShift, sBL=signalBinLabel), MCHistograms_prefiringWeightsUncertainties[signalType][signalBinLabel][UpDownShift])
            if (not(MCHistograms_prefiringWeightsUncertainties[signalType][signalBinLabel][UpDownShift])):
                sys.exit("ERROR: Histogram MCHistograms_prefiringWeightsUncertainties[signalType][signalBinLabel][UpDownShift] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
            MCHistograms_HLTUncertainties[signalType][signalBinLabel][UpDownShift] = ROOT.TH2F()
            MCUncertainties[signalType].GetObject("h_HLTUncertainty{UDS}_{sBL}".format(UDS=UpDownShift, sBL=signalBinLabel), MCHistograms_HLTUncertainties[signalType][signalBinLabel][UpDownShift])
            if (not(MCHistograms_HLTUncertainties[signalType][signalBinLabel][UpDownShift])):
                sys.exit("ERROR: Histogram MCHistograms_HLTUncertainties[signalType][signalBinLabel][UpDownShift] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
            MCHistograms_photonScaleFactorUncertainties[signalType][signalBinLabel][UpDownShift] = ROOT.TH2F()
            MCUncertainties[signalType].GetObject("h_photonMCScaleFactorUncertainty{UDS}_{sBL}".format(UDS=UpDownShift, sBL=signalBinLabel), MCHistograms_photonScaleFactorUncertainties[signalType][signalBinLabel][UpDownShift])
            if (not(MCHistograms_photonScaleFactorUncertainties[signalType][signalBinLabel][UpDownShift])):
                sys.exit("ERROR: Histogram MCHistograms_photonScaleFactorUncertainties[signalType][signalBinLabel][UpDownShift] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))

templateReader = MCTemplateReader.MCTemplateReader(inputArguments.MCTemplatePath)
# for indexPair in templateReader.nextValidBin():
# eventProgenitorBinIndex = indexPair[0]
eventProgenitorBinIndex = inputArguments.eventProgenitorMassBin
eventProgenitorMass = (templateReader.eventProgenitorMasses)[eventProgenitorBinIndex]
# neutralinoBinIndex = indexPair[1]
neutralinoBinIndex = inputArguments.neutralinoMassBin
neutralinoMass = (templateReader.neutralinoMasses)[neutralinoBinIndex]

# This extra scale factor is the correction to the nominal signal (calculated assuming the nominal cross-section).
# To move the cross-section down(up), we use cross_section_corrected = cross_section_nominal*{1 -(+) fractional_uncertainty}
crossSectionCorrectionScale = 1.0 + (inputArguments.crossSectionsScale)*crossSectionsFractionalUncertaintyDictionary[int(0.5+eventProgenitorMass)]

# MC systematics
systematics_MC = {}
systematics_MC_labels = []
systematics_MC_types = {}
# lumi systematic is fully correlated across all bins
systematics_MC_labels.append("lumi")
systematics_MC_types["lumi"] = "lnN"
systematics_MC["lumi"] = build_MC_constant_systematic(list_signalTypes=signalTypesToUse, dict_localToGlobalBinLabels=dict_localToGlobalBinLabels, constantFractionalUncertainty=inputArguments.luminosityUncertainty)

MCSystematicsSource_MCStatistics = {}
MCSystematicsSource_JEC = {}
MCSystematicsSource_Unclstrd = {}
MCSystematicsSource_JER = {}
MCSystematicsSource_pref = {}
MCSystematicsSource_HLT = {}
MCSystematicsSource_phoSF = {}
for signalType in signalTypesToUse:
    MCSystematicsSource_MCStatistics[signalType] = get_asymmetric_MC_systematic_from_histogram(localSignalLabels=localSignalBinLabels, inputHistograms=MCHistograms_MCStatUncertainties[signalType], eventProgenitorMassBin=eventProgenitorBinIndex, neutralinoMassBin=neutralinoBinIndex)
    MCSystematicsSource_JEC[signalType] = get_asymmetric_MC_systematic_from_histogram(localSignalLabels=localSignalBinLabels, inputHistograms=MCHistograms_JECUncertainties[signalType], eventProgenitorMassBin=eventProgenitorBinIndex, neutralinoMassBin=neutralinoBinIndex)
    MCSystematicsSource_Unclstrd[signalType] = get_asymmetric_MC_systematic_from_histogram(localSignalLabels=localSignalBinLabels, inputHistograms=MCHistograms_UnclusteredMETUncertainties[signalType], eventProgenitorMassBin=eventProgenitorBinIndex, neutralinoMassBin=neutralinoBinIndex)
    MCSystematicsSource_JER[signalType] = get_asymmetric_MC_systematic_from_histogram(localSignalLabels=localSignalBinLabels, inputHistograms=MCHistograms_JERMETUncertainties[signalType], eventProgenitorMassBin=eventProgenitorBinIndex, neutralinoMassBin=neutralinoBinIndex)
    MCSystematicsSource_pref[signalType] = get_asymmetric_MC_systematic_from_histogram(localSignalLabels=localSignalBinLabels, inputHistograms=MCHistograms_prefiringWeightsUncertainties[signalType], eventProgenitorMassBin=eventProgenitorBinIndex, neutralinoMassBin=neutralinoBinIndex)
    MCSystematicsSource_HLT[signalType] = get_asymmetric_MC_systematic_from_histogram(localSignalLabels=localSignalBinLabels, inputHistograms=MCHistograms_HLTUncertainties[signalType], eventProgenitorMassBin=eventProgenitorBinIndex, neutralinoMassBin=neutralinoBinIndex)
    MCSystematicsSource_phoSF[signalType] = get_asymmetric_MC_systematic_from_histogram(localSignalLabels=localSignalBinLabels, inputHistograms=MCHistograms_photonScaleFactorUncertainties[signalType], eventProgenitorMassBin=eventProgenitorBinIndex, neutralinoMassBin=neutralinoBinIndex)

for signalType in signalTypesToUse:
    # MCStatistics and HLT systematics in each bin are uncorrelated
    for signalBinLabel in localSignalBinLabels:
        localLabelsToUse = {signalType: [signalBinLabel]}
        tmp = build_MC_systematic_with_check(list_signalTypes=signalTypesToUse, dict_localToGlobalBinLabels=dict_localToGlobalBinLabels, dict_localSignalLabelsToUse=localLabelsToUse, dict_sources_dataSystematics=MCSystematicsSource_MCStatistics)
        if (tmp[0]):
            systematicLabel = "MCStatistics_{sBL}_{sT}".format(sBL=signalBinLabel, sT=signalType)
            systematics_MC_labels.append(systematicLabel)
            systematics_MC_types[systematicLabel] = "lnN"
            systematics_MC[systematicLabel] = tmp[1]
        tmp = build_MC_systematic_with_check(list_signalTypes=signalTypesToUse, dict_localToGlobalBinLabels=dict_localToGlobalBinLabels, dict_localSignalLabelsToUse=localLabelsToUse, dict_sources_dataSystematics=MCSystematicsSource_HLT)
        if (tmp[0]):
            systematicLabel = "HLT_{sBL}_{sT}".format(sBL=signalBinLabel, sT=signalType)
            systematics_MC_labels.append(systematicLabel)
            systematics_MC_types[systematicLabel] = "lnN"
            systematics_MC[systematicLabel] = tmp[1]

# All other MC uncertainties are correlated across all bins
localLabelsToUse = {signalType: localSignalBinLabels}
tmp = build_MC_systematic_with_check(list_signalTypes=signalTypesToUse, dict_localToGlobalBinLabels=dict_localToGlobalBinLabels, dict_localSignalLabelsToUse=localLabelsToUse, dict_sources_dataSystematics=MCSystematicsSource_JEC)
if (tmp[0]):
    systematics_MC_labels.append("JEC")
    systematics_MC_types["JEC"] = "lnN"
    systematics_MC["JEC"] = tmp[1]
tmp = build_MC_systematic_with_check(list_signalTypes=signalTypesToUse, dict_localToGlobalBinLabels=dict_localToGlobalBinLabels, dict_localSignalLabelsToUse=localLabelsToUse, dict_sources_dataSystematics=MCSystematicsSource_Unclstrd)
if (tmp[0]):
    systematics_MC_labels.append("Unclstrd")
    systematics_MC_types["Unclstrd"] = "lnN"
    systematics_MC["Unclstrd"] = tmp[1]
tmp = build_MC_systematic_with_check(list_signalTypes=signalTypesToUse, dict_localToGlobalBinLabels=dict_localToGlobalBinLabels, dict_localSignalLabelsToUse=localLabelsToUse, dict_sources_dataSystematics=MCSystematicsSource_JER)
if (tmp[0]):
    systematics_MC_labels.append("JER")
    systematics_MC_types["JER"] = "lnN"
    systematics_MC["JER"] = tmp[1]
tmp = build_MC_systematic_with_check(list_signalTypes=signalTypesToUse, dict_localToGlobalBinLabels=dict_localToGlobalBinLabels, dict_localSignalLabelsToUse=localLabelsToUse, dict_sources_dataSystematics=MCSystematicsSource_pref)
if (tmp[0]):
    systematics_MC_labels.append("pref")
    systematics_MC_types["pref"] = "lnN"
    systematics_MC["pref"] = tmp[1]
tmp = build_MC_systematic_with_check(list_signalTypes=signalTypesToUse, dict_localToGlobalBinLabels=dict_localToGlobalBinLabels, dict_localSignalLabelsToUse=localLabelsToUse, dict_sources_dataSystematics=MCSystematicsSource_phoSF)
if (tmp[0]):
    systematics_MC_labels.append("phoSF")
    systematics_MC_types["phoSF"] = "lnN"
    systematics_MC["phoSF"] = tmp[1]

print("Creating data card for eventProgenitor mass = {gM}, neutralino mass = {nM}".format(gM=eventProgenitorMass, nM=neutralinoMass))
expectedNEvents_stealth = {}
for signalType in signalTypesToUse:
    expectedNEventsLocal_stealth = get_dict_expectedNEvents_stealth(stealthNEventsHistograms=MCHistograms_weightedNEvents[signalType], eventProgenitorMassBin=eventProgenitorBinIndex, neutralinoMassBin=neutralinoBinIndex, localSignalLabels=localSignalBinLabels, scaleFactor=crossSectionCorrectionScale)
    for localLabel in localSignalBinLabels:
        globalLabel = dict_localToGlobalBinLabels[signalType][localLabel]
        expectedNEvents_stealth[globalLabel] = expectedNEventsLocal_stealth[localLabel]
createDataCard(outputPath="{oD}/{oP}_dataCard_eventProgenitorMassBin{gBI}_neutralinoMassBin{nBI}.txt".format(oD=inputArguments.outputDirectory, oP=inputArguments.outputPrefix, gBI=eventProgenitorBinIndex, nBI=neutralinoBinIndex),
               signalBinLabels=globalSignalBinLabels,
               observedNEvents=observedNEvents,
               expectedNEvents_qcd=expectedNEvents_qcd,
               expectedRate_scalingDeviation=expectedRate_scalingDeviation,
               systematics_data=systematics_data,
               systematics_data_labels=systematics_data_labels,
               systematics_data_types=systematics_data_types,
               expectedNEvents_stealth=expectedNEvents_stealth,
               systematics_MC=systematics_MC,
               systematics_MC_labels=systematics_MC_labels,
               systematics_MC_types=systematics_MC_types,
               rateParamLabels=rateParamLabels,
               rateParamProperties=rateParamProperties)

for signalType in signalTypesToUse:
    MCEventHistograms[signalType].Close()
    MCUncertainties[signalType].Close()
