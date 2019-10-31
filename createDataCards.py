#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, pdb, tmGeneralUtils, tmCombineDataCardInterface, os, sys, MCTemplateReader

inputArgumentsParser = argparse.ArgumentParser(description='Create data cards from MC and data systematics and nEvents data.')
inputArgumentsParser.add_argument('--outputPrefix', required=True, help='Prefix to output files.', type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="analysis/dataCards/", help='Output directory.', type=str)
inputArgumentsParser.add_argument('--crossSectionsFile', default="SusyCrossSections13TevGluGlu.txt", help='Path to dat file that contains cross-sections as a function of gluino mass, from which to get the fractional uncertainties.',type=str)
inputArgumentsParser.add_argument('--MCTemplatePath', default="plot_susyMasses_template.root", help='Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.', type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--inputFile_MCEventHistograms', required=True, help='Input MC event histograms.', type=str)
inputArgumentsParser.add_argument('--inputFile_MCUncertainties', required=True, help='Input MC uncertainties.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics', required=True, help='Input file containing fractional uncertainties from signal data.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics_sTScaling', required=True, help='Input file containing sT scaling systematics from control data.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics_expectedEventCounters', required=True, help='Input file containing expected number of events from signal data.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics_observedEventCounters', required=True, help='Input file containing observed number of events from signal data.', type=str)
inputArgumentsParser.add_argument('--luminosityUncertainty', required=True, help='Uncertainty on the luminosity.', type=float)
inputArgumentsParser.add_argument('--runUnblinded', action='store_true', help="If this flag is set, then the signal region data is unblinded. Specifically, the entry for the observed number of events is filled from the data, rather than from the expectation values.")
inputArguments = inputArgumentsParser.parse_args()

SYSTEMATIC_SIGNIFICANCE_THRESHOLD = 0.0005
LOGNORMAL_REGULARIZE_THRESHOLD = 0.01 # To "regularize" the lognormal error; the error in the column is supposed to be (1 plus deltaX/X), and if this is 0, the combine tool complains

def get_dict_nEvents(inputPath, signalLabels, nEventsPrefix):
    fileContents_nEvents = tmGeneralUtils.getConfigurationFromFile(inputFilePath=inputPath)
    outputDict = {}
    for signalBinLabel in signalLabels:
        outputDict[signalBinLabel] = fileContents_nEvents["{nEP}NEvents_{l}".format(nEP=nEventsPrefix, l=signalBinLabel)]
    return outputDict

def get_dict_expectedNEvents_stealth(stealthNEventsHistograms, gluinoMassBin, neutralinoMassBin, signalLabels, scaleFactor):
    outputDict = {}
    for signalBinLabel in signalLabels:
        outputDict[signalBinLabel] = scaleFactor*((stealthNEventsHistograms[signalBinLabel]).GetBinContent(gluinoMassBin, neutralinoMassBin))
    return outputDict

def get_symmetric_data_systematics_from_file(signalLabels, dataSystematicLabels, sourceFile):
    sourceSystematics = tmGeneralUtils.getConfigurationFromFile(sourceFile)
    outputDict = {}
    for dataSystematicLabel in dataSystematicLabels:
        outputDict[dataSystematicLabel] = {}
        for signalBinLabel in signalLabels:
            outputDict[dataSystematicLabel][signalBinLabel] = 1.0 + sourceSystematics["fractionalUncertainty_{dSL}_{sBL}".format(dSL=dataSystematicLabel, sBL=signalBinLabel)]
    return outputDict

def get_asymmetric_data_systematics_from_file(signalLabels, dataSystematicLabels, sourceFile):
    sourceSystematics = tmGeneralUtils.getConfigurationFromFile(sourceFile)
    outputDict = {}
    for dataSystematicLabel in dataSystematicLabels:
        outputDict[dataSystematicLabel] = {}
        for signalBinLabel in signalLabels:
            outputDict[dataSystematicLabel][signalBinLabel] = {}
            for upDownLabel in ["Down", "Up"]:
                outputDict[dataSystematicLabel][signalBinLabel][upDownLabel] = 1.0 + sourceSystematics["fractionalUncertainty{uDL}_{dSL}_{sBL}".format(uDL=upDownLabel, dSL=dataSystematicLabel, sBL=signalBinLabel)]
    return outputDict

def build_data_systematic_with_check(signalLabelsToUse, sourceDict_dataSystematics):
    outputDict = {}
    isSignificant = False
    for signalBinLabel in signalLabelsToUse:
        outputDict[signalBinLabel] = {}
        if isinstance(sourceDict_dataSystematics[signalBinLabel], dict):
            outputDict[signalBinLabel]["qcd"] = {}
            for upDownLabel in ["Down", "Up"]:
                try:
                    outputDict[signalBinLabel]["qcd"][upDownLabel] = sourceDict_dataSystematics[signalBinLabel][upDownLabel]
                    if (abs(outputDict[signalBinLabel]["qcd"][upDownLabel] - 1.0) > SYSTEMATIC_SIGNIFICANCE_THRESHOLD): isSignificant = True
                    if (abs(outputDict[signalBinLabel]["qcd"][upDownLabel]) < LOGNORMAL_REGULARIZE_THRESHOLD): outputDict[signalBinLabel]["qcd"][upDownLabel] = LOGNORMAL_REGULARIZE_THRESHOLD
                except KeyError:
                    sys.exit("ERROR: sourceDict_dataSystematics[signalBinLabel] is a dict but does not have elements named \"Up\" or \"Down\". dict contents: {c}".format(c=sourceDict_dataSystematics[signalBinLabel]))
        else:
            outputDict[signalBinLabel]["qcd"] = sourceDict_dataSystematics[signalBinLabel]
            if (abs(outputDict[signalBinLabel]["qcd"] - 1.0) > SYSTEMATIC_SIGNIFICANCE_THRESHOLD): isSignificant = True
            if (abs(outputDict[signalBinLabel]["qcd"]) < LOGNORMAL_REGULARIZE_THRESHOLD): outputDict[signalBinLabel]["qcd"] = LOGNORMAL_REGULARIZE_THRESHOLD
    return (isSignificant, outputDict)

def build_data_systematic_as_residual_with_check(signalLabelsToUse, sourceDict_dataSystematics, sourceDict_toTakeResidualWithRespectTo):
    outputDict = {}
    isSignificant = False
    for signalBinLabel in signalLabelsToUse:
        if ((isinstance(sourceDict_dataSystematics[signalBinLabel], dict)) or (isinstance(sourceDict_toTakeResidualWithRespectTo[signalBinLabel], dict))):
            sys.exit("ERROR: Asymmetric data systematics not yet implemented for residual data uncertainties.")
        outputDict[signalBinLabel] = {}
        outputDict[signalBinLabel]["qcd"] = 1.0 + max(0.0, sourceDict_dataSystematics[signalBinLabel] - sourceDict_toTakeResidualWithRespectTo[signalBinLabel])
        if (abs(outputDict[signalBinLabel]["qcd"] - 1.0) > SYSTEMATIC_SIGNIFICANCE_THRESHOLD): isSignificant = True
        if (abs(outputDict[signalBinLabel]["qcd"]) < LOGNORMAL_REGULARIZE_THRESHOLD): outputDict[signalBinLabel]["qcd"] = LOGNORMAL_REGULARIZE_THRESHOLD
    return (isSignificant, outputDict)

def build_MC_constant_systematic(signalLabels, constantFractionalUncertainty):
    outputDict = {}
    for signalBinLabel in signalLabels:
        outputDict[signalBinLabel] = {}
        outputDict[signalBinLabel]["stealth"] = 1.0 + constantFractionalUncertainty
    return outputDict

def get_MC_systematic_from_histogram(signalLabels, inputHistograms, gluinoMassBin, neutralinoMassBin):
    outputDict = {}
    for signalBinLabel in signalLabels:
        outputDict[signalBinLabel] = 1.0 + (inputHistograms[signalBinLabel]).GetBinContent(gluinoMassBin, neutralinoMassBin)
    return outputDict

def get_asymmetric_MC_systematic_from_histogram(signalLabels, inputHistograms, gluinoMassBin, neutralinoMassBin):
    outputDict = {}
    for signalBinLabel in signalLabels:
        outputDict[signalBinLabel] = {}
        for UpDownShift in ["Down", "Up"]:
            outputDict[signalBinLabel][UpDownShift] = 1.0 + (inputHistograms[signalBinLabel][UpDownShift]).GetBinContent(gluinoMassBin, neutralinoMassBin)
    return outputDict

def build_MC_systematic_with_check(signalLabelsToUse, sourceDict_MCSystematics):
    outputDict = {}
    isSignificant = False
    for signalBinLabel in signalLabelsToUse:
        outputDict[signalBinLabel] = {}
        if isinstance(sourceDict_MCSystematics[signalBinLabel], dict):
            outputDict[signalBinLabel]["stealth"] = {}
            for upDownLabel in ["Down", "Up"]:
                try:
                    outputDict[signalBinLabel]["stealth"][upDownLabel] = sourceDict_MCSystematics[signalBinLabel][upDownLabel]
                    if (abs(outputDict[signalBinLabel]["stealth"][upDownLabel] - 1.0) > SYSTEMATIC_SIGNIFICANCE_THRESHOLD): isSignificant = True
                    if (abs(outputDict[signalBinLabel]["stealth"][upDownLabel]) < LOGNORMAL_REGULARIZE_THRESHOLD): outputDict[signalBinLabel]["stealth"][upDownLabel] = LOGNORMAL_REGULARIZE_THRESHOLD
                except KeyError:
                    sys.exit("ERROR: sourceDict_MCSystematics[signalBinLabel] is a dict but does not have elements named \"Up\" or \"Down\". dict contents: {c}".format(c=sourceDict_MCSystematics[signalBinLabel]))
        else:
            outputDict[signalBinLabel]["stealth"] = sourceDict_MCSystematics[signalBinLabel]
            if (abs(outputDict[signalBinLabel]["stealth"] - 1.0) > SYSTEMATIC_SIGNIFICANCE_THRESHOLD): isSignificant = True
            if (abs(outputDict[signalBinLabel]["stealth"]) < LOGNORMAL_REGULARIZE_THRESHOLD): outputDict[signalBinLabel]["stealth"] = LOGNORMAL_REGULARIZE_THRESHOLD
    return (isSignificant, outputDict)

def createDataCard(outputPath,
                   signalBinLabels,
                   observedNEvents,
                   expectedNEvents_qcd,
                   systematics_data,
                   systematics_data_labels,
                   systematics_data_types,
                   expectedNEvents_stealth,
                   systematics_MC,
                   systematics_MC_labels,
                   systematics_MC_types):
    expectedNEvents = {}
    for signalBinLabel in signalBinLabels:
        expectedNEvents[signalBinLabel] = {}
        expectedNEvents[signalBinLabel]["qcd"] = expectedNEvents_qcd[signalBinLabel]
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

    combineInterface = tmCombineDataCardInterface.tmCombineDataCardInterface(list_signalBinLabels=signalBinLabels, list_backgroundProcessLabels=["qcd"], list_signalProcessLabels=["stealth"], list_systematicsLabels=systematicsLabels, dict_observedNEvents=observedNEvents, dict_expectedNEvents=expectedNEvents, dict_systematicsTypes=systematicsTypes, dict_systematics=systematics)
    combineInterface.writeToFile(outputFilePath=outputPath)

crossSectionsInputFileObject = open(inputArguments.crossSectionsFile, 'r')
crossSectionsDictionary = {}
crossSectionsFractionalUncertaintyDictionary = {}
for line in crossSectionsInputFileObject:
    crossSectionsData = line.split()
    gluinoMassInt = int(0.5 + float(crossSectionsData[0]))
    crossSection = float(crossSectionsData[1])
    crossSectionFractionalUncertainty = 0.01*float(crossSectionsData[2])
    crossSectionsDictionary[gluinoMassInt] = crossSection
    crossSectionsFractionalUncertaintyDictionary[gluinoMassInt] = crossSectionFractionalUncertainty
crossSectionsInputFileObject.close()

STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries)
nSTBoundaries = 0
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()): nSTBoundaries += 1
nSTSignalBins = nSTBoundaries - 2 + 1 # First two lines are for the normalization bin, last boundary implied at infinity
print("Using {n} signal bins for ST.".format(n = nSTSignalBins))
STRegionBoundariesFileObject.close()

signalBinLabels = []
for STRegionIndex in range(2, 2 + nSTSignalBins):
    for nJetsBin in range(4, 7):
        signalBinLabel = "STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin)
        signalBinLabels.append(signalBinLabel)

observedNEvents = {}
if (inputArguments.runUnblinded):
    observedNEvents = get_dict_nEvents(inputPath=inputArguments.inputFile_dataSystematics_observedEventCounters, signalLabels=signalBinLabels, nEventsPrefix="observed")
else:
    observedNEvents = get_dict_nEvents(inputPath=inputArguments.inputFile_dataSystematics_expectedEventCounters, signalLabels=signalBinLabels, nEventsPrefix="expected")
expectedNEvents_qcd = get_dict_nEvents(inputPath=inputArguments.inputFile_dataSystematics_expectedEventCounters, signalLabels=signalBinLabels, nEventsPrefix="expected")

# Data systematics
systematics_data = {}
systematics_data_labels = []
systematics_data_types = {}

sources_symmetricDataSystematics = get_symmetric_data_systematics_from_file(signalLabels=signalBinLabels, dataSystematicLabels=["shape", "rho"], sourceFile=inputArguments.inputFile_dataSystematics)
sources_asymmetricDataSystematics = get_asymmetric_data_systematics_from_file(signalLabels=signalBinLabels, dataSystematicLabels=["normEvents"], sourceFile=inputArguments.inputFile_dataSystematics)
sources_dataSystematics_scaling = get_symmetric_data_systematics_from_file(signalLabels=signalBinLabels, dataSystematicLabels=["scaling"], sourceFile=inputArguments.inputFile_dataSystematics_sTScaling)
# normEvents systematic is correlated across all ST bins
for nJetsBin in range(4, 7):
    labelsToUse = []
    for STRegionIndex in range(2, 2+nSTSignalBins):
        signalBinLabel = "STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin)
        labelsToUse.append(signalBinLabel)
    tmp = build_data_systematic_with_check(signalLabelsToUse=labelsToUse, sourceDict_dataSystematics=sources_asymmetricDataSystematics["normEvents"])
    if (tmp[0]):
        systematics_data_labels.append("normEvents_{n}Jets".format(n=nJetsBin))
        systematics_data_types["normEvents_{n}Jets".format(n=nJetsBin)] = "lnN"
        systematics_data["normEvents_{n}Jets".format(n=nJetsBin)] = tmp[1]
# shape and rho systematics are correlated across all nJets bins
for dataSystematic in ["shape", "rho"]:
    for STRegionIndex in range(2, 2+nSTSignalBins):
        labelsToUse = []
        for nJetsBin in range(4, 7):
            signalBinLabel = "STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin)
            labelsToUse.append(signalBinLabel)
        tmp = build_data_systematic_with_check(signalLabelsToUse=labelsToUse, sourceDict_dataSystematics=sources_symmetricDataSystematics[dataSystematic])
        if (tmp[0]):
            systematics_data_labels.append("{dS}_STRegion{r}".format(dS=dataSystematic, r=STRegionIndex))
            systematics_data_types["{dS}_STRegion{r}".format(dS=dataSystematic, r=STRegionIndex)] = "lnN"
            systematics_data["{dS}_STRegion{r}".format(dS=dataSystematic, r=STRegionIndex)] = tmp[1]
# scaling systematic is uncorrelated across all bins
for signalBinLabel in signalBinLabels:
    tmp = build_data_systematic_as_residual_with_check(signalLabelsToUse=[signalBinLabel], sourceDict_dataSystematics=sources_dataSystematics_scaling["scaling"], sourceDict_toTakeResidualWithRespectTo=sources_symmetricDataSystematics["shape"])
    if (tmp[0]):
        systematics_data_labels.append("scaling_{l}".format(l=signalBinLabel))
        systematics_data_types["scaling_{l}".format(l=signalBinLabel)] = "lnN"
        systematics_data["scaling_{l}".format(l=signalBinLabel)] = tmp[1]

MCHistograms_weightedNEvents = {}
MCHistograms_MCStatUncertainties = {}
MCHistograms_JECUncertainties = {}
MCHistograms_UnclusteredMETUncertainties = {}
MCHistograms_JERMETUncertainties = {}
MCHistograms_prefiringWeightsUncertainties = {}
MCHistograms_photonScaleFactorUncertainties = {}

MCEventHistograms = ROOT.TFile.Open(inputArguments.inputFile_MCEventHistograms)
MCUncertainties = ROOT.TFile.Open(inputArguments.inputFile_MCUncertainties)
expectedNEvents_stealth = None
for signalBinLabel in signalBinLabels:
    MCHistograms_weightedNEvents[signalBinLabel] = ROOT.TH2F()
    MCEventHistograms.GetObject("h_lumiBasedYearWeightedNEvents_{sBL}".format(sBL=signalBinLabel), MCHistograms_weightedNEvents[signalBinLabel])
    if (not(MCHistograms_weightedNEvents[signalBinLabel])):
        sys.exit("ERROR: Histogram MCHistograms_weightedNEvents[signalBinLabel] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))

    MCHistograms_MCStatUncertainties[signalBinLabel] = {}
    MCHistograms_JECUncertainties[signalBinLabel] = {}
    MCHistograms_UnclusteredMETUncertainties[signalBinLabel] = {}
    MCHistograms_JERMETUncertainties[signalBinLabel] = {}
    MCHistograms_prefiringWeightsUncertainties[signalBinLabel] = {}
    MCHistograms_photonScaleFactorUncertainties[signalBinLabel] = {}
    for UpDownShift in ["Down", "Up"]:
        MCHistograms_MCStatUncertainties[signalBinLabel][UpDownShift] = ROOT.TH2F()
        MCUncertainties.GetObject("h_MCStatisticsFractionalError{UDS}_{sBL}".format(UDS=UpDownShift, sBL=signalBinLabel), MCHistograms_MCStatUncertainties[signalBinLabel][UpDownShift])
        if (not(MCHistograms_MCStatUncertainties[signalBinLabel][UpDownShift])):
            sys.exit("ERROR: Histogram MCHistograms_MCStatUncertainties[signalBinLabel][UpDownShift] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
        MCHistograms_JECUncertainties[signalBinLabel][UpDownShift] = ROOT.TH2F()
        MCUncertainties.GetObject("h_JECUncertainty{UDS}_{sBL}".format(UDS=UpDownShift, sBL=signalBinLabel), MCHistograms_JECUncertainties[signalBinLabel][UpDownShift])
        if (not(MCHistograms_JECUncertainties[signalBinLabel][UpDownShift])):
            sys.exit("ERROR: Histogram MCHistograms_JECUncertainties[signalBinLabel][UpDownShift] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
        MCHistograms_UnclusteredMETUncertainties[signalBinLabel][UpDownShift] = ROOT.TH2F()
        MCUncertainties.GetObject("h_UnclusteredMETUncertainty{UDS}_{sBL}".format(UDS=UpDownShift, sBL=signalBinLabel), MCHistograms_UnclusteredMETUncertainties[signalBinLabel][UpDownShift])
        if (not(MCHistograms_UnclusteredMETUncertainties[signalBinLabel][UpDownShift])):
            sys.exit("ERROR: Histogram MCHistograms_UnclusteredMETUncertainties[signalBinLabel][UpDownShift] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
        MCHistograms_JERMETUncertainties[signalBinLabel][UpDownShift] = ROOT.TH2F()
        MCUncertainties.GetObject("h_JERMETUncertainty{UDS}_{sBL}".format(UDS=UpDownShift, sBL=signalBinLabel), MCHistograms_JERMETUncertainties[signalBinLabel][UpDownShift])
        if (not(MCHistograms_JERMETUncertainties[signalBinLabel][UpDownShift])):
            sys.exit("ERROR: Histogram MCHistograms_JERMETUncertainties[signalBinLabel][UpDownShift] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
        MCHistograms_prefiringWeightsUncertainties[signalBinLabel][UpDownShift] = ROOT.TH2F()
        MCUncertainties.GetObject("h_prefiringWeightsUncertainty{UDS}_{sBL}".format(UDS=UpDownShift, sBL=signalBinLabel), MCHistograms_prefiringWeightsUncertainties[signalBinLabel][UpDownShift])
        if (not(MCHistograms_prefiringWeightsUncertainties[signalBinLabel][UpDownShift])):
            sys.exit("ERROR: Histogram MCHistograms_prefiringWeightsUncertainties[signalBinLabel][UpDownShift] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
        MCHistograms_photonScaleFactorUncertainties[signalBinLabel][UpDownShift] = ROOT.TH2F()
        MCUncertainties.GetObject("h_photonMCScaleFactorUncertainty{UDS}_{sBL}".format(UDS=UpDownShift, sBL=signalBinLabel), MCHistograms_photonScaleFactorUncertainties[signalBinLabel][UpDownShift])
        if (not(MCHistograms_photonScaleFactorUncertainties[signalBinLabel][UpDownShift])):
            sys.exit("ERROR: Histogram MCHistograms_photonScaleFactorUncertainties[signalBinLabel][UpDownShift] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))

templateReader = MCTemplateReader.MCTemplateReader(inputArguments.MCTemplatePath)
for indexPair in templateReader.nextValidBin():
    gluinoBinIndex = indexPair[0]
    gluinoMass = (templateReader.gluinoMasses)[gluinoBinIndex]
    neutralinoBinIndex = indexPair[1]
    neutralinoMass = (templateReader.neutralinoMasses)[neutralinoBinIndex]
    crossSectionFractionalUncertaintyScaleFactor = 1.0 + crossSectionsFractionalUncertaintyDictionary[gluinoMass]

    # MC systematics
    systematics_MC = {}
    systematics_MC_labels = []
    systematics_MC_types = {}
    # lumi systematic is fully correlated across all bins
    systematics_MC_labels.append("lumi")
    systematics_MC_types["lumi"] = "lnN"
    systematics_MC["lumi"] = build_MC_constant_systematic(signalLabels=signalBinLabels, constantFractionalUncertainty=inputArguments.luminosityUncertainty)
    # MCStatistics systematic in each bin is uncorrelated
    MCSystematicsSource_MCStatistics=get_asymmetric_MC_systematic_from_histogram(signalLabels=signalBinLabels, inputHistograms=MCHistograms_MCStatUncertainties, gluinoMassBin=gluinoBinIndex, neutralinoMassBin=neutralinoBinIndex)
    for signalBinLabel in signalBinLabels:
        tmp = build_MC_systematic_with_check(signalLabelsToUse=[signalBinLabel], sourceDict_MCSystematics=MCSystematicsSource_MCStatistics)
        if (tmp[0]):
            systematics_MC_labels.append("MCStatistics_{sBL}".format(sBL=signalBinLabel))
            systematics_MC_types["MCStatistics_{sBL}".format(sBL=signalBinLabel)] = "lnN"
            systematics_MC["MCStatistics_{sBL}".format(sBL=signalBinLabel)] = tmp[1]
    # All other uncertainties are correlated across all bins
    tmp = build_MC_systematic_with_check(signalLabelsToUse=signalBinLabels, sourceDict_MCSystematics=get_asymmetric_MC_systematic_from_histogram(signalLabels=signalBinLabels, inputHistograms=MCHistograms_JECUncertainties, gluinoMassBin=gluinoBinIndex, neutralinoMassBin=neutralinoBinIndex))
    if (tmp[0]):
        systematics_MC_labels.append("JEC")
        systematics_MC_types["JEC"] = "lnN"
        systematics_MC["JEC"] = tmp[1]
    tmp = build_MC_systematic_with_check(signalLabelsToUse=signalBinLabels, sourceDict_MCSystematics=get_asymmetric_MC_systematic_from_histogram(signalLabels=signalBinLabels, inputHistograms=MCHistograms_UnclusteredMETUncertainties, gluinoMassBin=gluinoBinIndex, neutralinoMassBin=neutralinoBinIndex))
    if (tmp[0]):
        systematics_MC_labels.append("Unclstrd")
        systematics_MC_types["Unclstrd"] = "lnN"
        systematics_MC["Unclstrd"] = tmp[1]
    tmp = build_MC_systematic_with_check(signalLabelsToUse=signalBinLabels, sourceDict_MCSystematics=get_asymmetric_MC_systematic_from_histogram(signalLabels=signalBinLabels, inputHistograms=MCHistograms_JERMETUncertainties, gluinoMassBin=gluinoBinIndex, neutralinoMassBin=neutralinoBinIndex))
    if (tmp[0]):
        systematics_MC_labels.append("JER")
        systematics_MC_types["JER"] = "lnN"
        systematics_MC["JER"] = tmp[1]
    tmp = build_MC_systematic_with_check(signalLabelsToUse=signalBinLabels, sourceDict_MCSystematics=get_asymmetric_MC_systematic_from_histogram(signalLabels=signalBinLabels, inputHistograms=MCHistograms_prefiringWeightsUncertainties, gluinoMassBin=gluinoBinIndex, neutralinoMassBin=neutralinoBinIndex))
    if (tmp[0]):
        systematics_MC_labels.append("pref")
        systematics_MC_types["pref"] = "lnN"
        systematics_MC["pref"] = tmp[1]
    tmp = build_MC_systematic_with_check(signalLabelsToUse=signalBinLabels, sourceDict_MCSystematics=get_asymmetric_MC_systematic_from_histogram(signalLabels=signalBinLabels, inputHistograms=MCHistograms_photonScaleFactorUncertainties, gluinoMassBin=gluinoBinIndex, neutralinoMassBin=neutralinoBinIndex))
    if (tmp[0]):
        systematics_MC_labels.append("phoSF")
        systematics_MC_types["phoSF"] = "lnN"
        systematics_MC["phoSF"] = tmp[1]

    print("Creating data cards for gluino mass = {gM}, neutralino mass = {nM}".format(gM=gluinoMass, nM=neutralinoMass))
    expectedNEvents_stealth = get_dict_expectedNEvents_stealth(stealthNEventsHistograms=MCHistograms_weightedNEvents, gluinoMassBin=gluinoBinIndex, neutralinoMassBin=neutralinoBinIndex, signalLabels=signalBinLabels, scaleFactor=1.0)
    createDataCard(outputPath="{oD}/{oP}_dataCard_gluinoMassBin{gBI}_neutralinoMassBin{nBI}.txt".format(oD=inputArguments.outputDirectory, oP=inputArguments.outputPrefix, gBI=gluinoBinIndex, nBI=neutralinoBinIndex),
                   signalBinLabels=signalBinLabels,
                   observedNEvents=observedNEvents,
                   expectedNEvents_qcd=expectedNEvents_qcd,
                   systematics_data=systematics_data,
                   systematics_data_labels=systematics_data_labels,
                   systematics_data_types=systematics_data_types,
                   expectedNEvents_stealth=expectedNEvents_stealth,
                   systematics_MC=systematics_MC,
                   systematics_MC_labels=systematics_MC_labels,
                   systematics_MC_types=systematics_MC_types)
    expectedNEvents_stealth = get_dict_expectedNEvents_stealth(stealthNEventsHistograms=MCHistograms_weightedNEvents, gluinoMassBin=gluinoBinIndex, neutralinoMassBin=neutralinoBinIndex, signalLabels=signalBinLabels, scaleFactor=1.0/crossSectionFractionalUncertaintyScaleFactor)
    createDataCard(outputPath="{oD}/{oP}_dataCard_gluinoMassBin{gBI}_neutralinoMassBin{nBI}_crossSectionsDown.txt".format(oD=inputArguments.outputDirectory, oP=inputArguments.outputPrefix, gBI=gluinoBinIndex, nBI=neutralinoBinIndex),
                   signalBinLabels=signalBinLabels,
                   observedNEvents=observedNEvents,
                   expectedNEvents_qcd=expectedNEvents_qcd,
                   systematics_data=systematics_data,
                   systematics_data_labels=systematics_data_labels,
                   systematics_data_types=systematics_data_types,
                   expectedNEvents_stealth=expectedNEvents_stealth,
                   systematics_MC=systematics_MC,
                   systematics_MC_labels=systematics_MC_labels,
                   systematics_MC_types=systematics_MC_types)
    expectedNEvents_stealth = get_dict_expectedNEvents_stealth(stealthNEventsHistograms=MCHistograms_weightedNEvents, gluinoMassBin=gluinoBinIndex, neutralinoMassBin=neutralinoBinIndex, signalLabels=signalBinLabels, scaleFactor=crossSectionFractionalUncertaintyScaleFactor)
    createDataCard(outputPath="{oD}/{oP}_dataCard_gluinoMassBin{gBI}_neutralinoMassBin{nBI}_crossSectionsUp.txt".format(oD=inputArguments.outputDirectory, oP=inputArguments.outputPrefix, gBI=gluinoBinIndex, nBI=neutralinoBinIndex),
                   signalBinLabels=signalBinLabels,
                   observedNEvents=observedNEvents,
                   expectedNEvents_qcd=expectedNEvents_qcd,
                   systematics_data=systematics_data,
                   systematics_data_labels=systematics_data_labels,
                   systematics_data_types=systematics_data_types,
                   expectedNEvents_stealth=expectedNEvents_stealth,
                   systematics_MC=systematics_MC,
                   systematics_MC_labels=systematics_MC_labels,
                   systematics_MC_types=systematics_MC_types)
MCEventHistograms.Close()
MCUncertainties.Close()
