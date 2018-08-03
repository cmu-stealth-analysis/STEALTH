#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, tmROOTUtils, pdb, tmGeneralUtils, os, sys

inputArgumentsParser = argparse.ArgumentParser(description='Create data cards from MC and data systematics and nEvents data.')
inputArgumentsParser.add_argument('--outputPrefix', required=True, help='Prefix to output files.', type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="analysis/dataCards/", help='Output directory.', type=str)
inputArgumentsParser.add_argument('--inputFile_MCTemplate', default="plot_susyMasses_template.root", help='Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.', type=str)
inputArgumentsParser.add_argument('--inputFile_MCEventHistograms', required=True, help='Input MC event histograms.', type=str)
inputArgumentsParser.add_argument('--inputFile_MCUncertainties', required=True, help='Input MC uncertainties.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics', required=True, help='Input file containing fractional uncertainties from signal data.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics_sTScaling', required=True, help='Input file containing sT scaling systematics from control data.', type=str)
inputArgumentsParser.add_argument('--inputFile_dataSystematics_eventCounters', required=True, help='Input file containing expected number of events from signal data.', type=str)
inputArgumentsParser.add_argument('--luminosityUncertainty', required=True, help='Uncertainty on the luminosity.', type=float)
inputArgumentsParser.add_argument('--defaultValue_JECUncertainty', default=0.10, help='Default value for the JEC uncertainty -- this value is used if the number of events in the MC sample passing the selection is 0.', type=float)
inputArgumentsParser.add_argument('--defaultValue_statUncertainty', default=1.00, help='Default value for the MC statistical uncertainty -- this value is used if the number of events in the MC sample passing the selection is 0.', type=float)
inputArguments = inputArgumentsParser.parse_args()

def getFixedWidthString(width, precision, number):
    formatter = ""
    if (width == 0):
        formatter = "{{n:.{p}f}}".format(p=precision)
    else:
        formatter = "{{n:<{w}.{p}f}}".format(w=width, p=precision)
    returnString =formatter.format(n=number)
    returnStringValue = float(returnString)
    fractionalError = (returnStringValue - number)/(1 + number) # In most cases, we're going to add 1.0 to "number" anyway -- this stops the code from complaining when, for example, 1.0045 is  translated to 1.005
    if (fractionalError > 0.01):
        sys.exit("ERROR: the number {n} is not accurately translated into a floating point representation: {rep}".format(n=number, rep=returnString))
    return returnString

def createDataCard(outputDirectory, outputFileName, lookupTable):
    dataCardTemplate = open('{outputDirectory}/{outputFileName}.txt'.format(outputDirectory=outputDirectory, outputFileName=outputFileName), 'w')
    dataCardTemplate.write("# Auto-generated by the script \"createDataCards.py\"\n")
    dataCardTemplate.write("imax 6  number of channels\n")
    dataCardTemplate.write("jmax 1  number of backgrounds\n")
    dataCardTemplate.write("kmax 7  number of nuisance parameters (sources of systematic uncertainties)\n")
    dataCardTemplate.write("------------\n")
    dataCardTemplate.write("bin            sub4Jets     main4Jets    sub5Jets     main5Jets    sub6Jets     main6Jets\n")
    dataCardTemplate.write("observation    {sub4Jets}   {main4Jets}  {sub5Jets}   {main5Jets}  {sub6Jets}   {main6Jets}\n".format(**lookupTable)) # temporary, while data is unblinded -- useful for expected limit plots
    dataCardTemplate.write("------------\n")
    dataCardTemplate.write("bin                  sub4Jets      sub4Jets      main4Jets     main4Jets      sub5Jets      sub5Jets      main5Jets     main5Jets      sub6Jets      sub6Jets      main6Jets     main6Jets\n")
    dataCardTemplate.write("process              t7Wg          qcd           t7Wg          qcd            t7Wg          qcd           t7Wg          qcd            t7Wg          qcd           t7Wg          qcd\n")
    dataCardTemplate.write("process              0             1             0             1              0             1             0             1              0             1             0             1\n")
    dataCardTemplate.write("rate                 {mcsub_4J}    {sub4Jets}    {mcmain_4J}   {main4Jets}    {mcsub_5J}    {sub5Jets}    {mcmain_5J}   {main5Jets}    {mcsub_6J}    {sub6Jets}    {mcmain_6J}   {main6Jets}\n".format(**lookupTable))
    dataCardTemplate.write("------------\n")
    dataCardTemplate.write("normEvents   lnN     -             {normUnc}     -             {normUnc}      -             {normUnc}     -             {normUnc}      -             {normUnc}     -             {normUnc}\n".format(**lookupTable))
    dataCardTemplate.write("shape        lnN     -             {shapeUnc}    -             {shapeUnc}     -             {shapeUnc}    -             {shapeUnc}     -             {shapeUnc}    -             {shapeUnc}\n".format(**lookupTable))
    dataCardTemplate.write("scaling      lnN     -             {scaleUnc4J}  -             {scaleUnc4J}   -             {scaleUnc5J}  -             {scaleUnc5J}   -             {scaleUnc6J}  -             {scaleUnc6J}\n".format(**lookupTable))
    dataCardTemplate.write("rho          lnN     -             {rhoUnc}      -             {rhoUnc}       -             {rhoUnc}      -             {rhoUnc}       -             {rhoUnc}      -             {rhoUnc}\n".format(**lookupTable))
    dataCardTemplate.write("jetE         lnN     {jec_sub4J}   -             {jec_main4J}  -              {jec_sub5J}   -             {jec_main5J}  -              {jec_sub6J}   -             {jec_main6J}  -\n".format(**lookupTable))
    dataCardTemplate.write("lumi         lnN     {lumiUnc}     -             {lumiUnc}     -              {lumiUnc}     -             {lumiUnc}     -              {lumiUnc}     -             {lumiUnc}     -\n".format(**lookupTable))
    dataCardTemplate.write("MCStats      lnN     {stat_sub4J}  -             {stat_main4J} -              {stat_sub5J}  -             {stat_main5J} -              {stat_sub6J}  -             {stat_main6J} -\n".format(**lookupTable))
    dataCardTemplate.close()

dataSystematics = tmGeneralUtils.getConfigurationFromFile(inputArguments.inputFile_dataSystematics)
dataSystematics_sTScaling = tmGeneralUtils.getConfigurationFromFile(inputArguments.inputFile_dataSystematics_sTScaling)
eventCounters_data = tmGeneralUtils.getConfigurationFromFile(inputArguments.inputFile_dataSystematics_eventCounters)

lookupTable_unformatted = {}
for nJetsBin in range(4, 7):
    lookupTable_unformatted["sub{n}Jets".format(n=nJetsBin)] = eventCounters_data["expectedNEvents_subordinateRegion_{n}Jets".format(n=nJetsBin)]
    lookupTable_unformatted["main{n}Jets".format(n=nJetsBin)] = eventCounters_data["expectedNEvents_mainRegion_{n}Jets".format(n=nJetsBin)]
    lookupTable_unformatted["scaleUnc{n}J".format(n=nJetsBin)] = 1.0 + dataSystematics_sTScaling["fractionalUncertainty_sTScaling_{n}Jets".format(n=nJetsBin)]
lookupTable_unformatted["normUnc"] = 1.0 + dataSystematics["fractionalUncertainty_normEvents"]
lookupTable_unformatted["shapeUnc"] = 1.0 + dataSystematics["fractionalUncertainty_Shape"]
lookupTable_unformatted["rhoUnc"] = 1.0 + dataSystematics["fractionalUncertainty_rho"]
lookupTable_unformatted["lumiUnc"] = 1.0 + inputArguments.luminosityUncertainty

lookupTable = {}
for lookupItem in lookupTable_unformatted.keys():
    lookupItemLength = len(lookupItem)
    lookupTable[lookupItem] = getFixedWidthString(2 + lookupItemLength, 3, lookupTable_unformatted[lookupItem]) # lookupItemLength + 2 for braces; precision defaults to 3 everywhere

generatedMCTemplate = ROOT.TFile(inputArguments.inputFile_MCTemplate)
h_MCTemplate = generatedMCTemplate.Get("h_susyMasses_template")

histograms_weightedNEvents = {
    "main": {},
    "sub": {}
}
histograms_JECUncertainties = {
    "main": {},
    "sub": {}
}
histograms_MCStatUncertainties = {
    "main": {},
    "sub": {}
}
MCEventHistograms = ROOT.TFile(inputArguments.inputFile_MCEventHistograms)
MCUncertainties = ROOT.TFile(inputArguments.inputFile_MCUncertainties)
for nJetsBin in range(4, 7):
    for region in ["sub", "main"]:
        histograms_weightedNEvents[region][nJetsBin] = MCEventHistograms.Get("h_weighted_nMCEvents_JECNominal_{n}Jets_{r}".format(n=nJetsBin, r=region))
        histograms_JECUncertainties[region][nJetsBin] = MCUncertainties.Get("h_JECUncertainty_{n}Jets_{r}".format(n=nJetsBin, r=region))
        histograms_MCStatUncertainties[region][nJetsBin] = MCUncertainties.Get("h_MCStatisticsFractionalError_{n}Jets_{r}".format(n=nJetsBin, r=region))
    
for gluinoMassBin in range(1, 1+h_MCTemplate.GetXaxis().GetNbins()):
    for neutralinoMassBin in range(1, 1+h_MCTemplate.GetYaxis().GetNbins()):
        if not(int(0.5 + h_MCTemplate.GetBinContent(gluinoMassBin, neutralinoMassBin)) == 1): continue
        gluinoMass = h_MCTemplate.GetXaxis().GetBinCenter(gluinoMassBin)
        neutralinoMass = h_MCTemplate.GetYaxis().GetBinCenter(neutralinoMassBin)
        print("Creating data card for gluino mass = {gM}, neutralino mass={nM}".format(gM=gluinoMass, nM=neutralinoMass))
        tempLookupTable_unformatted = {}
        for nJetsBin in range(4, 7):
            for region in ["sub", "main"]:
                weightedNEvents = (histograms_weightedNEvents[region][nJetsBin]).GetBinContent((histograms_weightedNEvents[region][nJetsBin]).FindFixBin(gluinoMass, neutralinoMass))
                tempLookupTable_unformatted["mc{r}_{n}J".format(r=region, n=nJetsBin)] = weightedNEvents
                jecUncertainty = inputArguments.defaultValue_JECUncertainty
                statUncertainty = inputArguments.defaultValue_statUncertainty
                if (weightedNEvents > 0):
                    statUncertainty = (histograms_MCStatUncertainties[region][nJetsBin]).GetBinContent((histograms_MCStatUncertainties[region][nJetsBin]).FindFixBin(gluinoMass, neutralinoMass))
                    jecUncertainty = (histograms_JECUncertainties[region][nJetsBin]).GetBinContent((histograms_JECUncertainties[region][nJetsBin]).FindFixBin(gluinoMass, neutralinoMass))
                    if (jecUncertainty < 0.001): jecUncertainty = 0.001 # see below
                    # JEC uncertainty is small in two cases:
                    # (1) Not enough statistics in a given bin, in which case the statistical uncertainty will be high, so it's fine to set jecUncertainty to 0.001
                    # (2) High statistics but genuinely low JEC uncertainty in a bin, in which case a 0.1% uncertainty is defendable
                tempLookupTable_unformatted["jec_{r}{n}J".format(r=region, n=nJetsBin)] = 1.0 + jecUncertainty
                tempLookupTable_unformatted["stat_{r}{n}J".format(r=region, n=nJetsBin)] = 1.0 + statUncertainty
        for lookupItem in tempLookupTable_unformatted.keys():
            lookupItemLength = len(lookupItem)
            lookupTable[lookupItem] = getFixedWidthString(2 + lookupItemLength, 3, tempLookupTable_unformatted[lookupItem]) # lookupItemLength + 2 for braces; precision defaults to 3 everywhere
        createDataCard(inputArguments.outputDirectory, "{outputPrefix}_dataCard_gluinoMassBin{gMB}_neutralinoMassBin{nMB}".format(outputPrefix=inputArguments.outputPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin), lookupTable)
        for lookupItem in tempLookupTable_unformatted.keys():
            lookupTable.pop(lookupItem)
        tempLookupTable_unformatted.clear()