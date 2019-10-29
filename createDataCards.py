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

def build_data_systematic(signalLabelsToUse, sourceDict_dataSystematics):
    outputDict = {}
    for signalBinLabel in signalLabelsToUse:
        outputDict[signalBinLabel] = {}
        if isinstance(sourceDict_dataSystematics[signalBinLabel], dict):
            outputDict[signalBinLabel]["qcd"] = {}
            for upDownLabel in ["Down", "Up"]:
                try:
                    outputDict[signalBinLabel]["qcd"][upDownLabel] = sourceDict_dataSystematics[signalBinLabel][upDownLabel]
                except KeyError:
                    sys.exit("ERROR: sourceDict_dataSystematics[signalBinLabel] is a dict but does not have elements named \"Up\" or \"Down\". dict contents: {c}".format(c=sourceDict_dataSystematics[signalBinLabel]))
        else:
            outputDict[signalBinLabel]["qcd"] = sourceDict_dataSystematics[signalBinLabel]
    return outputDict

def build_data_systematic_as_residual(signalLabelsToUse, sourceDict_dataSystematics, sourceDict_toTakeResidualWithRespectTo):
    outputDict = {}
    for signalBinLabel in signalLabelsToUse:
        if ((isinstance(sourceDict_dataSystematics[signalBinLabel], dict)) or (isinstance(sourceDict_toTakeResidualWithRespectTo[signalBinLabel], dict))):
            sys.exit("ERROR: Asymmetric data systematics not yet implemented for residual data uncertainties.")
        outputDict[signalBinLabel] = {}
        outputDict[signalBinLabel]["qcd"] = 1.0 + max(0.0, sourceDict_dataSystematics[signalBinLabel] - sourceDict_toTakeResidualWithRespectTo[signalBinLabel])
    return outputDict

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

def build_MC_systematic(signalLabels, sourceDict_MCSystematics):
    outputDict = {}
    for signalBinLabel in signalLabels:
        outputDict[signalBinLabel] = {}
        outputDict[signalBinLabel]["stealth"] = sourceDict_MCSystematics[signalBinLabel]
    return outputDict

def createDataCard(runUnblinded,
                   outputPath,
                   nJetsMin, nJetsMax, STRegionMin, STRegionMax,
                   signalBinLabels,
                   inputPath_observedNEvents,
                   inputPath_expectedNEvents,
                   inputPath_dataSystematics,
                   inputPath_dataSystematics_STScaling,
                   gluinoMassBin, neutralinoMassBin,
                   inputHistograms_stealthNEvents,
                   crossSectionsScaleFactor,
                   lumiUncertainty,
                   inputHistograms_MCUncertainty_MCStats,
                   inputHistograms_MCUncertainty_JEC,
                   inputHistograms_MCUncertainty_Unclstrd,
                   inputHistograms_MCUncertainty_JER,
                   inputHistograms_MCUncertainty_pref,
                   inputHistograms_MCUncertainty_phoSF):
    observedNEvents = {}
    if (runUnblinded):
        observedNEvents = get_dict_nEvents(inputPath=inputPath_observedNEvents, signalLabels=signalBinLabels, nEventsPrefix="observed")
    else:
        observedNEvents = get_dict_nEvents(inputPath=inputPath_expectedNEvents, signalLabels=signalBinLabels, nEventsPrefix="expected")
    backgroundProcessLabels = ["qcd"]
    signalProcessLabels = ["stealth"]
    expectedNEvents_qcd = get_dict_nEvents(inputPath=inputPath_expectedNEvents, signalLabels=signalBinLabels, nEventsPrefix="expected")
    expectedNEvents_stealth = get_dict_expectedNEvents_stealth(stealthNEventsHistograms=inputHistograms_stealthNEvents, gluinoMassBin=gluinoMassBin, neutralinoMassBin=neutralinoMassBin, signalLabels=signalBinLabels, scaleFactor=crossSectionsScaleFactor)
    expectedNEvents = {}
    for signalBinLabel in signalBinLabels:
        expectedNEvents[signalBinLabel] = {}
        expectedNEvents[signalBinLabel]["qcd"] = expectedNEvents_qcd[signalBinLabel]
        expectedNEvents[signalBinLabel]["stealth"] = expectedNEvents_stealth[signalBinLabel]
    # systematicsLabels = ["normEvents", "shape", "rho", "scaling", "lumi", "MCStats", "JEC", "Unclstrd", "JER", "pref", "phoSF"]
    systematics = {}
    systematicsLabels = []
    systematicsTypes = {}

    # Data systematics
    sources_symmetricDataSystematics = get_symmetric_data_systematics_from_file(signalLabels=signalBinLabels, dataSystematicLabels=["shape", "rho"], sourceFile=inputPath_dataSystematics)
    sources_asymmetricDataSystematics = get_asymmetric_data_systematics_from_file(signalLabels=signalBinLabels, dataSystematicLabels=["normEvents"], sourceFile=inputPath_dataSystematics)
    sources_dataSystematics_scaling = get_symmetric_data_systematics_from_file(signalLabels=signalBinLabels, dataSystematicLabels=["scaling"], sourceFile=inputPath_dataSystematics_STScaling)
    # normEvents systematic is correlated across all ST bins
    for nJetsBin in range(nJetsMin, 1+nJetsMax):
        labelsToUse = []
        for STRegionIndex in range(STRegionMin, 1+STRegionMax):
            signalBinLabel = "STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin)
            labelsToUse.append(signalBinLabel)
        systematicsLabels.append("normEvents_{n}Jets".format(n=nJetsBin))
        systematics["normEvents_{n}Jets".format(n=nJetsBin)] = build_data_systematic(signalLabelsToUse=labelsToUse, sourceDict_dataSystematics=sources_asymmetricDataSystematics["normEvents"])
    # shape and rho systematics are correlated across all nJets bins
    for dataSystematic in ["shape", "rho"]:
        for STRegionIndex in range(STRegionMin, 1+STRegionMax):
            labelsToUse = []
            for nJetsBin in range(nJetsMin, 1+nJetsMax):
                signalBinLabel = "STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin)
                labelsToUse.append(signalBinLabel)
            systematicsLabels.append("{dS}_STRegion{r}".format(dS=dataSystematic, r=STRegionIndex))
            systematics["{dS}_STRegion{r}".format(dS=dataSystematic, r=STRegionIndex)] = build_data_systematic(signalLabelsToUse=labelsToUse, sourceDict_dataSystematics=sources_symmetricDataSystematics[dataSystematic])
    # scaling systematic is uncorrelated across all bins
    for signalBinLabel in signalBinLabels:
        systematicsLabels.append("scaling_{l}".format(l=signalBinLabel))
        systematics["scaling_{l}".format(l=signalBinLabel)] = build_data_systematic_as_residual(signalLabelsToUse=[signalBinLabel], sourceDict_dataSystematics=sources_dataSystematics_scaling["scaling"], sourceDict_toTakeResidualWithRespectTo=sources_symmetricDataSystematics["shape"])

    # MC systematics
    systematicsLabels.append("lumi")
    systematics["lumi"] = build_MC_constant_systematic(signalLabels=signalBinLabels, constantFractionalUncertainty=lumiUncertainty)
    systematicsLabels.append("MCStats")
    systematics["MCStats"] = build_MC_systematic(signalLabels=signalBinLabels, sourceDict_MCSystematics=get_MC_systematic_from_histogram(signalLabels=signalBinLabels, inputHistograms=inputHistograms_MCUncertainty_MCStats, gluinoMassBin=gluinoMassBin, neutralinoMassBin=neutralinoMassBin))
    systematicsLabels.append("JEC")
    systematics["JEC"] = build_MC_systematic(signalLabels=signalBinLabels, sourceDict_MCSystematics=get_MC_systematic_from_histogram(signalLabels=signalBinLabels, inputHistograms=inputHistograms_MCUncertainty_JEC, gluinoMassBin=gluinoMassBin, neutralinoMassBin=neutralinoMassBin))
    systematicsLabels.append("Unclstrd")
    systematics["Unclstrd"] = build_MC_systematic(signalLabels=signalBinLabels, sourceDict_MCSystematics=get_MC_systematic_from_histogram(signalLabels=signalBinLabels, inputHistograms=inputHistograms_MCUncertainty_Unclstrd, gluinoMassBin=gluinoMassBin, neutralinoMassBin=neutralinoMassBin))
    systematicsLabels.append("JER")
    systematics["JER"] = build_MC_systematic(signalLabels=signalBinLabels, sourceDict_MCSystematics=get_MC_systematic_from_histogram(signalLabels=signalBinLabels, inputHistograms=inputHistograms_MCUncertainty_JER, gluinoMassBin=gluinoMassBin, neutralinoMassBin=neutralinoMassBin))
    systematicsLabels.append("pref")
    systematics["pref"] = build_MC_systematic(signalLabels=signalBinLabels, sourceDict_MCSystematics=get_MC_systematic_from_histogram(signalLabels=signalBinLabels, inputHistograms=inputHistograms_MCUncertainty_pref, gluinoMassBin=gluinoMassBin, neutralinoMassBin=neutralinoMassBin))
    systematicsLabels.append("phoSF")
    systematics["phoSF"] = build_MC_systematic(signalLabels=signalBinLabels, sourceDict_MCSystematics=get_MC_systematic_from_histogram(signalLabels=signalBinLabels, inputHistograms=inputHistograms_MCUncertainty_phoSF, gluinoMassBin=gluinoMassBin, neutralinoMassBin=neutralinoMassBin))

    for systematic in systematicsLabels:
        systematicsTypes[systematic] = "lnN"

    combineInterface = tmCombineDataCardInterface.tmCombineDataCardInterface(list_signalBinLabels=signalBinLabels, list_backgroundProcessLabels=backgroundProcessLabels, list_signalProcessLabels=signalProcessLabels, list_systematicsLabels=systematicsLabels, dict_observedNEvents=observedNEvents, dict_expectedNEvents=expectedNEvents, dict_systematicsTypes=systematicsTypes, dict_systematics=systematics)
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

histograms_weightedNEvents = {}
histograms_MCStatUncertainties = {}
histograms_JECUncertainties = {}
histograms_UnclusteredMETUncertainties = {}
histograms_JERMETUncertainties = {}
histograms_prefiringWeightsUncertainties = {}
histograms_photonScaleFactorUncertainties = {}

MCEventHistograms = ROOT.TFile.Open(inputArguments.inputFile_MCEventHistograms)
MCUncertainties = ROOT.TFile.Open(inputArguments.inputFile_MCUncertainties)
signalBinLabels = []
for STRegionIndex in range(2, 2 + nSTSignalBins):
    for nJetsBin in range(4, 7):
        signalBinLabel = "STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin)
        signalBinLabels.append(signalBinLabel)
        histograms_weightedNEvents[signalBinLabel] = ROOT.TH2F()
        MCEventHistograms.GetObject("h_lumiBasedYearWeightedNEvents_{n}Jets_STRegion{r}".format(n=nJetsBin, r=STRegionIndex), histograms_weightedNEvents[signalBinLabel])
        if (not(histograms_weightedNEvents[signalBinLabel])):
            sys.exit("ERROR: Histogram histograms_weightedNEvents[signalBinLabel] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
        histograms_MCStatUncertainties[signalBinLabel] = ROOT.TH2F()
        MCUncertainties.GetObject("h_MCStatisticsFractionalError_{n}Jets_STRegion{r}".format(n=nJetsBin, r=STRegionIndex), histograms_MCStatUncertainties[signalBinLabel])
        if (not(histograms_MCStatUncertainties[signalBinLabel])):
            sys.exit("ERROR: Histogram histograms_MCStatUncertainties[signalBinLabel] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
        histograms_JECUncertainties[signalBinLabel] = ROOT.TH2F()
        MCUncertainties.GetObject("h_JECUncertainty_{n}Jets_STRegion{r}".format(n=nJetsBin, r=STRegionIndex), histograms_JECUncertainties[signalBinLabel])
        if (not(histograms_JECUncertainties[signalBinLabel])):
            sys.exit("ERROR: Histogram histograms_JECUncertainties[signalBinLabel] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
        histograms_UnclusteredMETUncertainties[signalBinLabel] = ROOT.TH2F()
        MCUncertainties.GetObject("h_UnclusteredMETUncertainty_{n}Jets_STRegion{r}".format(n=nJetsBin, r=STRegionIndex), histograms_UnclusteredMETUncertainties[signalBinLabel])
        if (not(histograms_UnclusteredMETUncertainties[signalBinLabel])):
            sys.exit("ERROR: Histogram histograms_UnclusteredMETUncertainties[signalBinLabel] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
        histograms_JERMETUncertainties[signalBinLabel] = ROOT.TH2F()
        MCUncertainties.GetObject("h_JERMETUncertainty_{n}Jets_STRegion{r}".format(n=nJetsBin, r=STRegionIndex), histograms_JERMETUncertainties[signalBinLabel])
        if (not(histograms_JERMETUncertainties[signalBinLabel])):
            sys.exit("ERROR: Histogram histograms_JERMETUncertainties[signalBinLabel] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
        histograms_prefiringWeightsUncertainties[signalBinLabel] = ROOT.TH2F()
        MCUncertainties.GetObject("h_prefiringWeightsUncertainty_{n}Jets_STRegion{r}".format(n=nJetsBin, r=STRegionIndex), histograms_prefiringWeightsUncertainties[signalBinLabel])
        if (not(histograms_prefiringWeightsUncertainties[signalBinLabel])):
            sys.exit("ERROR: Histogram histograms_prefiringWeightsUncertainties[signalBinLabel] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))
        histograms_photonScaleFactorUncertainties[signalBinLabel] = ROOT.TH2F()
        MCUncertainties.GetObject("h_photonMCScaleFactorUncertainty_{n}Jets_STRegion{r}".format(n=nJetsBin, r=STRegionIndex), histograms_photonScaleFactorUncertainties[signalBinLabel])
        if (not(histograms_photonScaleFactorUncertainties[signalBinLabel])):
            sys.exit("ERROR: Histogram histograms_photonScaleFactorUncertainties[signalBinLabel] appears to be a nullptr at STRegionIndex={r}, nJets={n}".format(r=STRegionIndex, n=nJetsBin))

templateReader = MCTemplateReader.MCTemplateReader(inputArguments.MCTemplatePath)
for indexPair in templateReader.nextValidBin():
    gluinoBinIndex = indexPair[0]
    gluinoMass = (templateReader.gluinoMasses)[gluinoBinIndex]
    neutralinoBinIndex = indexPair[1]
    neutralinoMass = (templateReader.neutralinoMasses)[neutralinoBinIndex]
    crossSectionFractionalUncertaintyScaleFactor = 1.0 + crossSectionsFractionalUncertaintyDictionary[gluinoMass]
    print("Creating data cards for gluino mass = {gM}, neutralino mass = {nM}".format(gM=gluinoMass, nM=neutralinoMass))
    createDataCard(runUnblinded=inputArguments.runUnblinded,
                   outputPath="{oD}/{oP}_dataCard_gluinoMassBin{gBI}_neutralinoMassBin{nBI}.txt".format(oD=inputArguments.outputDirectory, oP=inputArguments.outputPrefix, gBI=gluinoBinIndex, nBI=neutralinoBinIndex),
                   nJetsMin=4, nJetsMax=6, STRegionMin=2, STRegionMax=1+nSTSignalBins,
                   signalBinLabels=signalBinLabels,
                   inputPath_observedNEvents=inputArguments.inputFile_dataSystematics_observedEventCounters,
                   inputPath_expectedNEvents=inputArguments.inputFile_dataSystematics_expectedEventCounters,
                   inputPath_dataSystematics=inputArguments.inputFile_dataSystematics,
                   inputPath_dataSystematics_STScaling=inputArguments.inputFile_dataSystematics_sTScaling,
                   gluinoMassBin=gluinoBinIndex, neutralinoMassBin=neutralinoBinIndex,
                   inputHistograms_stealthNEvents=histograms_weightedNEvents,
                   crossSectionsScaleFactor=1.0,
                   lumiUncertainty=inputArguments.luminosityUncertainty,
                   inputHistograms_MCUncertainty_MCStats=histograms_MCStatUncertainties,
                   inputHistograms_MCUncertainty_JEC=histograms_JECUncertainties,
                   inputHistograms_MCUncertainty_Unclstrd=histograms_UnclusteredMETUncertainties,
                   inputHistograms_MCUncertainty_JER=histograms_JERMETUncertainties,
                   inputHistograms_MCUncertainty_pref=histograms_prefiringWeightsUncertainties,
                   inputHistograms_MCUncertainty_phoSF=histograms_photonScaleFactorUncertainties)
    createDataCard(runUnblinded=inputArguments.runUnblinded,
                   outputPath="{oD}/{oP}_dataCard_gluinoMassBin{gBI}_neutralinoMassBin{nBI}_crossSectionsDown.txt".format(oD=inputArguments.outputDirectory, oP=inputArguments.outputPrefix, gBI=gluinoBinIndex, nBI=neutralinoBinIndex),
                   nJetsMin=4, nJetsMax=6, STRegionMin=2, STRegionMax=1+nSTSignalBins,
                   signalBinLabels=signalBinLabels,
                   inputPath_observedNEvents=inputArguments.inputFile_dataSystematics_observedEventCounters,
                   inputPath_expectedNEvents=inputArguments.inputFile_dataSystematics_expectedEventCounters,
                   inputPath_dataSystematics=inputArguments.inputFile_dataSystematics,
                   inputPath_dataSystematics_STScaling=inputArguments.inputFile_dataSystematics_sTScaling,
                   gluinoMassBin=gluinoBinIndex, neutralinoMassBin=neutralinoBinIndex,
                   inputHistograms_stealthNEvents=histograms_weightedNEvents,
                   crossSectionsScaleFactor=1.0/crossSectionFractionalUncertaintyScaleFactor,
                   lumiUncertainty=inputArguments.luminosityUncertainty,
                   inputHistograms_MCUncertainty_MCStats=histograms_MCStatUncertainties,
                   inputHistograms_MCUncertainty_JEC=histograms_JECUncertainties,
                   inputHistograms_MCUncertainty_Unclstrd=histograms_UnclusteredMETUncertainties,
                   inputHistograms_MCUncertainty_JER=histograms_JERMETUncertainties,
                   inputHistograms_MCUncertainty_pref=histograms_prefiringWeightsUncertainties,
                   inputHistograms_MCUncertainty_phoSF=histograms_photonScaleFactorUncertainties)
    createDataCard(runUnblinded=inputArguments.runUnblinded,
                   outputPath="{oD}/{oP}_dataCard_gluinoMassBin{gBI}_neutralinoMassBin{nBI}_crossSectionsUp.txt".format(oD=inputArguments.outputDirectory, oP=inputArguments.outputPrefix, gBI=gluinoBinIndex, nBI=neutralinoBinIndex),
                   nJetsMin=4, nJetsMax=6, STRegionMin=2, STRegionMax=1+nSTSignalBins,
                   signalBinLabels=signalBinLabels,
                   inputPath_observedNEvents=inputArguments.inputFile_dataSystematics_observedEventCounters,
                   inputPath_expectedNEvents=inputArguments.inputFile_dataSystematics_expectedEventCounters,
                   inputPath_dataSystematics=inputArguments.inputFile_dataSystematics,
                   inputPath_dataSystematics_STScaling=inputArguments.inputFile_dataSystematics_sTScaling,
                   gluinoMassBin=gluinoBinIndex, neutralinoMassBin=neutralinoBinIndex,
                   inputHistograms_stealthNEvents=histograms_weightedNEvents,
                   crossSectionsScaleFactor=crossSectionFractionalUncertaintyScaleFactor,
                   lumiUncertainty=inputArguments.luminosityUncertainty,
                   inputHistograms_MCUncertainty_MCStats=histograms_MCStatUncertainties,
                   inputHistograms_MCUncertainty_JEC=histograms_JECUncertainties,
                   inputHistograms_MCUncertainty_Unclstrd=histograms_UnclusteredMETUncertainties,
                   inputHistograms_MCUncertainty_JER=histograms_JERMETUncertainties,
                   inputHistograms_MCUncertainty_pref=histograms_prefiringWeightsUncertainties,
                   inputHistograms_MCUncertainty_phoSF=histograms_photonScaleFactorUncertainties)
MCEventHistograms.Close()
MCUncertainties.Close()
