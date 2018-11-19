#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, tmROOTUtils, pdb, tmGeneralUtils, os, sys

inputArgumentsParser = argparse.ArgumentParser(description='Create data cards from MC and data systematics and nEvents data.')
inputArgumentsParser.add_argument('--outputPrefix', required=True, help='Prefix to output files.', type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="analysis/dataCards/", help='Output directory.', type=str)
inputArgumentsParser.add_argument('--inputFile_MCTemplate', default="plot_susyMasses_template.root", help='Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.', type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--inputFile_MCEventHistograms', required=True, help='Input MC event histograms.', type=str)
inputArgumentsParser.add_argument('--inputFile_MCUncertainties', required=True, help='Input MC uncertainties.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics', required=True, help='Input file containing fractional uncertainties from signal data.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics_sTScaling', required=True, help='Input file containing sT scaling systematics from control data.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics_eventCounters', required=True, help='Input file containing expected number of events from signal data.', type=str)
inputArgumentsParser.add_argument('--luminosityUncertainty', required=True, help='Uncertainty on the luminosity.', type=float)
inputArgumentsParser.add_argument('--defaultValue_JECUncertainty', default=0.10, help='Default value for the JEC uncertainty -- this value is used if the number of events in the MC sample passing the selection is 0.', type=float)
inputArgumentsParser.add_argument('--defaultValue_statUncertainty', default=1.00, help='Default value for the MC statistical uncertainty -- this value is used if the number of events in the MC sample passing the selection is 0.', type=float)
inputArguments = inputArgumentsParser.parse_args()

def alignFixedWidthFloatLeft(width, precision, number):
    if not(isinstance(number, float)): sys.exit("alignFixedWidthFloatLeft called with non-float object: {o}".format(o=number))
    formatter = ""
    if (width == 0):
        formatter = "{{n:.{p}f}}".format(p=precision)
    else:
        formatter = "{{n:<{w}.{p}f}}".format(w=width, p=precision)
    returnString =formatter.format(n=number)
    # if (number == 0): return returnString
    # returnStringValue = float(returnString)
    # fractionalError = (returnStringValue - number)/number
    # if (fractionalError > 0.01):
    #     sys.exit("ERROR: the number {n} is not accurately translated into a floating point representation: {rep}".format(n=number, rep=returnString))
    return returnString

def alignFixedWidthStringLeft(width, inputString):
    if not(isinstance(inputString, str)): sys.exit("alignFixedWidthStringLeft called with non-string object: {o}".format(o=inputString))
    formatter = ""
    formatter = "{{s:<{w}}}".format(w=width)
    returnString = formatter.format(s=inputString)
    return returnString

def createDataCard(outputDirectory, outputFileName, lookupTable, nSTSignalBins):
    dataCardTemplate = open('{outputDirectory}/{outputFileName}.txt'.format(outputDirectory=outputDirectory, outputFileName=outputFileName), 'w')
    dataCardTemplate.write("# Auto-generated by the script \"createDataCards.py\"\n")
    dataCardTemplate.write("imax {nC}  number of channels\n".format(nC=nSTSignalBins*3)) # 3 jets bins in range (4, 7)
    dataCardTemplate.write("jmax 1   number of backgrounds\n")
    dataCardTemplate.write("kmax 7   number of nuisance parameters (sources of systematic uncertainties)\n")
    dataCardTemplate.write("------------\n")

    binDescriptions = alignFixedWidthStringLeft(15, "bin") # Example:                       "bin            STReg2_4Jets    STReg3_4Jets    STReg4_4Jets    STReg2_5Jets    STReg3_5Jets    STReg4_5Jets\n"
    binObservations = alignFixedWidthStringLeft(15, "observation") # Example:               "observation    ndata_r2_4J     ndata_r3_4J     ndata_r4_4J     ndata_r2_5J     ndata_r3_5J     ndata_r4_5J\n"
    for STRegionIndex in range(2, 2 + nSTSignalBins):
        for nJetsBin in range(4, 7):
            binDescriptions += alignFixedWidthStringLeft(20, ("STReg{i}_{n}Jets").format(i=STRegionIndex, n=nJetsBin))
            binObservations += alignFixedWidthFloatLeft(20, 3, lookupTable["ndata_r{i}_{n}J".format(i=STRegionIndex, n=nJetsBin)]) # temporary, while data is unblinded -- useful for expected limit plots
    dataCardTemplate.write(binDescriptions.rstrip() + "\n")
    dataCardTemplate.write(binObservations.rstrip() + "\n")
    dataCardTemplate.write("------------\n")
    binTitles = alignFixedWidthStringLeft(20, "bin") # Example:          "bin                   STReg2_4Jets    STReg2_4Jets    STReg3_4Jets    STReg3_4Jets    STReg4_4Jets    STReg4_4Jets    STReg2_5Jets    STReg2_5Jets    STReg3_5Jets    STReg3_5Jets    STReg4_5Jets    STReg4_5Jets\n"
    processLabels = alignFixedWidthStringLeft(20, "process") # Example:  "process               t7Wg            qcd             t7Wg            qcd             t7Wg            qcd             t7Wg            qcd             t7Wg            qcd             t7Wg            qcd\n"
    processIndices = alignFixedWidthStringLeft(20, "process") # Example: "process               0               1               0               1               0               1               0               1               0               1               0               1\n"
    processRates = alignFixedWidthStringLeft(20, "rate") # Example:      "rate                  nmc_r2_4J       ndata_r2_4J     nmc_r3_4J       ndata_r3_4J     nmc_r4_4J       ndata_r4_4J     nmc_r2_5J       ndata_r2_5J     nmc_r3_5J       ndata_r3_5J     nmc_r4_5J       ndata_r4_5J\n"

    for nJetsBin in range(4, 7):
        for STRegionIndex in range(2, 2 + nSTSignalBins):
            for processIndex in range(2):
                binTitles += alignFixedWidthStringLeft(17, ("STReg{i}_{n}Jets").format(i=STRegionIndex, n=nJetsBin))
                processIndices += alignFixedWidthStringLeft(17, "{i}".format(i=processIndex))
                if (processIndex == 0): # MC
                    processLabels += alignFixedWidthStringLeft(17, "t7Wg")
                    processRates += alignFixedWidthFloatLeft(17, 3, lookupTable["nmc_r{i}_{n}J".format(i=STRegionIndex, n=nJetsBin)])
                else: # data
                    processLabels += alignFixedWidthStringLeft(17, "qcd")
                    processRates += alignFixedWidthFloatLeft(17, 3, lookupTable["ndata_r{i}_{n}J".format(i=STRegionIndex, n=nJetsBin)])
    dataCardTemplate.write(binTitles.rstrip() + "\n")
    dataCardTemplate.write(processLabels.rstrip() + "\n")
    dataCardTemplate.write(processIndices.rstrip() + "\n")
    dataCardTemplate.write(processRates.rstrip() + "\n")
    dataCardTemplate.write("------------\n")

    # Data uncertainties
    normUncertainties = alignFixedWidthStringLeft(13, "normEvents") + alignFixedWidthStringLeft(7, "lnN")
    shapeUncertainties = alignFixedWidthStringLeft(13, "shape") + alignFixedWidthStringLeft(7, "lnN")
    scaleUncertainties = alignFixedWidthStringLeft(13, "scaling") + alignFixedWidthStringLeft(7, "lnN")
    rhoUncertainties = alignFixedWidthStringLeft(13, "rho") + alignFixedWidthStringLeft(7, "lnN")

    # MC uncertainties
    jetEUncertainties = alignFixedWidthStringLeft(13, "jetE") + alignFixedWidthStringLeft(7, "lnN")
    lumiUncertainties = alignFixedWidthStringLeft(13, "lumi") + alignFixedWidthStringLeft(7, "lnN")
    MCStatsUncertainties = alignFixedWidthStringLeft(13, "MCStats") + alignFixedWidthStringLeft(7, "lnN")
    for nJetsBin in range(4, 7):
        for STRegionIndex in range(2, 2 + nSTSignalBins):
            for processIndex in range(2):
                if (processIndex == 0): # MC
                    normUncertainties += alignFixedWidthStringLeft(17, "-")
                    shapeUncertainties += alignFixedWidthStringLeft(17, "-")
                    scaleUncertainties += alignFixedWidthStringLeft(17, "-")
                    rhoUncertainties += alignFixedWidthStringLeft(17, "-")

                    jetEUncertainties += alignFixedWidthFloatLeft(17, 3, lookupTable["jec_r{i}_{n}J".format(i=STRegionIndex, n=nJetsBin)])
                    lumiUncertainties += alignFixedWidthFloatLeft(17, 3, lookupTable["lumiUnc"])
                    MCStatsUncertainties += alignFixedWidthFloatLeft(17, 3, lookupTable["stat_r{i}_{n}J".format(i=STRegionIndex, n=nJetsBin)])
                else: # data
                    normUncertainties += alignFixedWidthFloatLeft(17, 3, lookupTable["normUnc"])
                    shapeUncertainties += alignFixedWidthFloatLeft(17, 3, lookupTable["shapeUnc"])
                    scaleUncertainties += alignFixedWidthFloatLeft(17, 3, lookupTable["scaleUnc{n}J".format(n=nJetsBin)])
                    rhoUncertainties += alignFixedWidthFloatLeft(17, 3, lookupTable["rhoUnc"])

                    jetEUncertainties += alignFixedWidthStringLeft(17, "-")
                    lumiUncertainties += alignFixedWidthStringLeft(17, "-")
                    MCStatsUncertainties += alignFixedWidthStringLeft(17, "-")
    dataCardTemplate.write(normUncertainties.rstrip() + "\n")
    dataCardTemplate.write(shapeUncertainties.rstrip() + "\n")
    dataCardTemplate.write(scaleUncertainties.rstrip() + "\n")
    dataCardTemplate.write(rhoUncertainties.rstrip() + "\n")
    dataCardTemplate.write(jetEUncertainties.rstrip() + "\n")
    dataCardTemplate.write(lumiUncertainties.rstrip() + "\n")
    dataCardTemplate.write(MCStatsUncertainties.rstrip() + "\n")
    dataCardTemplate.write("\n")
    dataCardTemplate.close()

STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries)
nSTBoundaries = 0
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()): nSTBoundaries += 1
nSTSignalBins = nSTBoundaries - 2 + 1 # First two lines are for the normalization bin, last boundary implied at infinity
print("Using {n} signal bins for ST.".format(n = nSTSignalBins))

dataSystematics = tmGeneralUtils.getConfigurationFromFile(inputArguments.inputFile_dataSystematics)
dataSystematics_sTScaling = tmGeneralUtils.getConfigurationFromFile(inputArguments.inputFile_dataSystematics_sTScaling)
eventCounters_data = tmGeneralUtils.getConfigurationFromFile(inputArguments.inputFile_dataSystematics_eventCounters)

lookupTable = {}
for nJetsBin in range(4, 7):
    for STRegionIndex in range(2, 2 + nSTSignalBins): # region index 1 is for norm bin
        lookupTable["ndata_r{i}_{n}J".format(i=STRegionIndex, n=nJetsBin)] = eventCounters_data["expectedNEvents_STRegion{i}_{n}Jets".format(i=STRegionIndex, n=nJetsBin)]
    lookupTable["scaleUnc{n}J".format(n=nJetsBin)] = 1.0 + dataSystematics_sTScaling["fractionalUncertainty_sTScaling_{n}Jets".format(n=nJetsBin)]
lookupTable["normUnc"] = 1.0 + dataSystematics["fractionalUncertainty_normEvents"]
lookupTable["shapeUnc"] = 1.0 + dataSystematics["fractionalUncertainty_Shape"]
lookupTable["rhoUnc"] = 1.0 + dataSystematics["fractionalUncertainty_rho"]
lookupTable["lumiUnc"] = 1.0 + inputArguments.luminosityUncertainty

# lookupTable = {}
# for lookupItem in lookupTable_unformatted.keys():
#     nCharacters_lookupKey = len(lookupItem)
#     lookupTable[lookupItem] = getFixedWidthString(2 + nCharacters_lookupKey, 3, lookupTable_unformatted[lookupItem]) # nCharacters_lookupKey + 2 for braces; precision defaults to 3 everywhere

histograms_weightedNEvents = {}
histograms_JECUncertainties = {}
histograms_MCStatUncertainties = {}

for STRegionIndex in range(2, 2 + nSTSignalBins):
    histograms_weightedNEvents[STRegionIndex] = {}
    histograms_JECUncertainties[STRegionIndex] = {}
    histograms_MCStatUncertainties[STRegionIndex] = {}

MCEventHistograms = ROOT.TFile(inputArguments.inputFile_MCEventHistograms)
MCUncertainties = ROOT.TFile(inputArguments.inputFile_MCUncertainties)
for nJetsBin in range(4, 7):
    for STRegionIndex in range(2, 2 + nSTSignalBins):
        # histograms_weightedNEvents[STRegionIndex][nJetsBin] = MCEventHistograms.Get("h_weighted_nMCEvents_JECNominal_{n}Jets_STRegion{r}".format(n=nJetsBin, r=STRegionIndex))
        histograms_weightedNEvents[STRegionIndex][nJetsBin] = MCEventHistograms.Get("h_lumiBasedYearWeighted_nMCEvents_{n}Jets_STRegion{r}".format(n=nJetsBin, r=STRegionIndex))
        histograms_JECUncertainties[STRegionIndex][nJetsBin] = MCUncertainties.Get("h_JECUncertainty_{n}Jets_STRegion{r}".format(n=nJetsBin, r=STRegionIndex))
        histograms_MCStatUncertainties[STRegionIndex][nJetsBin] = MCUncertainties.Get("h_MCStatisticsFractionalError_{n}Jets_STRegion{r}".format(n=nJetsBin, r=STRegionIndex))

generatedMCTemplate = ROOT.TFile(inputArguments.inputFile_MCTemplate)
h_MCTemplate = generatedMCTemplate.Get("h_susyMasses_template")
for gluinoMassBin in range(1, 1+h_MCTemplate.GetXaxis().GetNbins()):
    for neutralinoMassBin in range(1, 1+h_MCTemplate.GetYaxis().GetNbins()):
        if not(int(0.5 + h_MCTemplate.GetBinContent(gluinoMassBin, neutralinoMassBin)) == 1): continue
        gluinoMass = h_MCTemplate.GetXaxis().GetBinCenter(gluinoMassBin)
        neutralinoMass = h_MCTemplate.GetYaxis().GetBinCenter(neutralinoMassBin)
        print("Creating data card for gluino mass = {gM}, neutralino mass = {nM}".format(gM=gluinoMass, nM=neutralinoMass))
        tempLookupTable = {}
        for nJetsBin in range(4, 7):
            for STRegionIndex in range(2, 2 + nSTSignalBins):
                weightedNEvents = (histograms_weightedNEvents[STRegionIndex][nJetsBin]).GetBinContent((histograms_weightedNEvents[STRegionIndex][nJetsBin]).FindFixBin(gluinoMass, neutralinoMass))
                tempLookupTable["nmc_r{i}_{n}J".format(i=STRegionIndex, n=nJetsBin)] = weightedNEvents
                jecUncertainty = inputArguments.defaultValue_JECUncertainty
                statUncertainty = inputArguments.defaultValue_statUncertainty
                if (weightedNEvents > 0):
                    jecUncertainty = (histograms_JECUncertainties[STRegionIndex][nJetsBin]).GetBinContent((histograms_JECUncertainties[STRegionIndex][nJetsBin]).FindFixBin(gluinoMass, neutralinoMass))
                    statUncertainty = (histograms_MCStatUncertainties[STRegionIndex][nJetsBin]).GetBinContent((histograms_MCStatUncertainties[STRegionIndex][nJetsBin]).FindFixBin(gluinoMass, neutralinoMass))
                    if (jecUncertainty < 0.001): jecUncertainty = 0.001 # see below
                    # JEC uncertainty is small in two cases:
                    # (1) Not enough statistics in a given bin, so statistics drives this uncertainty down to 0, in which case the statistical uncertainty will be high, so it's fine to set jecUncertainty to 0.1%
                    # (2) High statistics but genuinely low JEC uncertainty in a bin --> in this case as well JEC will not be the dominant uncertainty, so 0.1% is OK
                tempLookupTable["jec_r{i}_{n}J".format(i=STRegionIndex, n=nJetsBin)] = 1.0 + jecUncertainty
                tempLookupTable["stat_r{i}_{n}J".format(i=STRegionIndex, n=nJetsBin)] = 1.0 + statUncertainty
        for lookupItem in tempLookupTable.keys():
            lookupTable[lookupItem] = tempLookupTable[lookupItem]
        createDataCard(inputArguments.outputDirectory, "{outputPrefix}_dataCard_gluinoMassBin{gMB}_neutralinoMassBin{nMB}".format(outputPrefix=inputArguments.outputPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin), lookupTable, nSTSignalBins)
        for lookupItem in tempLookupTable.keys():
            lookupTable.pop(lookupItem)
        tempLookupTable.clear()
