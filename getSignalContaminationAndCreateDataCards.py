#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, tmROOTUtils, pdb

from tmProgressBar import tmProgressBar
from tmGeneralUtils import prettyPrintDictionary

inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination from input MC in real data.')
inputArgumentsParser.add_argument('--inputMCPath', required=True, help='Path to input MC file.',type=str)
inputArgumentsParser.add_argument('--maxMCEvents', default=0, help='Set a custom maximum number of MC events.',type=int)
inputArgumentsParser.add_argument('--inputDataPath', required=True, help='Path to input data file.',type=str)
inputArgumentsParser.add_argument('--crossSectionsFile', default="SusyCrossSections13TevGluGlu.txt", help='Path to dat file that contains cross-sections as a function of gluino mass, to use while weighting events.',type=str)
inputArgumentsParser.add_argument('--MCTemplate', default="plot_susyMasses_template.root", help='Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.', type=str)
inputArgumentsParser.add_argument('--dataCardTemplate', required=True, help='Path to Higgs Combine Tool datacard template.', type=str)
inputArgumentsParser.add_argument('--sTMin_normWindow', default=1000., help='Min value of sT.',type=float)
inputArgumentsParser.add_argument('--sTMax_normWindow', default=1100., help='Max value of sT.',type=float)
inputArgumentsParser.add_argument('--sTStartMainRegion', default=2500., help='Lowest value of sT in main observation bin.',type=float)
inputArgumentsParser.add_argument('--outputPrefix', required=True, help='Prefix to output files.',type=str)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Minimum number of jets in event.',type=int)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Least value of nJets in highest nJets bin.',type=int)
inputArgumentsParser.add_argument('--analyze_nJetsBin', action='append', default=[], help='nJets to plot.',type=int)
inputArgumentsParser.add_argument('--nGluinoMassBins', default=20, help='nBins on the gluino mass axis.',type=int) # (800 - 25) GeV --> (1750 + 25) GeV in steps of 50 GeV
inputArgumentsParser.add_argument('--minGluinoMass', default=775., help='Min gluino mass for the 2D plots.',type=float)
inputArgumentsParser.add_argument('--maxGluinoMass', default=1775., help='Max gluino mass for the 2D plots.',type=float)
inputArgumentsParser.add_argument('--nNeutralinoMassBins', default=133, help='nBins on the neutralino mass axis.',type=int)
inputArgumentsParser.add_argument('--minNeutralinoMass', default=93.75, help='Min neutralino mass for the 2D plots.',type=float)
inputArgumentsParser.add_argument('--maxNeutralinoMass', default=1756.25, help='Max neutralino mass for the 2D plots.',type=float) # (100 - 6.25) GeV --> (1750 + 6.25) GeV in steps of 12.5 GeV
inputArgumentsParser.add_argument('--totalIntegratedLuminosity', default=37760., help='Total integrated luminosity for the total data-taking period.',type=float) # total integrated luminosities for 2016 + 2017 = 37.76 (2016) + 46.02 (2017) fb^{-1} = 83780 pb^{-1}; default = 2016 only
inputArgumentsParser.add_argument('--nGeneratedEventsPerBin', default=150000, help='Number of generated events per bin in the MC samples.',type=int)
inputArguments = inputArgumentsParser.parse_args()

if not(inputArguments.nJetsMax == 6): sys.exit("Only nJetsMax=6 supported temporarily. Needed to fill in data template in the correct format.")
nJetsBinsToAnalyze=inputArguments.analyze_nJetsBin
for nJetsBinToAnalyze in nJetsBinsToAnalyze:
    if nJetsBinToAnalyze > inputArguments.nJetsMax: sys.exit("Error: index of jets bin to analyze can be at most nJetsMax={nJetsMax}; this argument is not allowed: --analyze_nJetsBin={bin}".format(nJetsMax=inputArguments.nJetsMax, bin=nJetsBinToAnalyze))

# pdb.set_trace()

MCPIDs = {
    "photon": 22,
    "gluino": 1000021,
    "neutralino": 1000022
}

crossSectionsInputFileObject = open(inputArguments.crossSectionsFile, 'r')
crossSectionsDictionary = {}
crossSectionsFractionalUncertaintyDictionary = {}
for line in crossSectionsInputFileObject:
    crossSectionsData = line.split()
    gluinoMass = int(0.5 + float(crossSectionsData[0]))
    crossSection = float(crossSectionsData[1])
    crossSectionFractionalUncertainty = 0.01*float(crossSectionsData[2])
    crossSectionsDictionary[gluinoMass] = crossSection
    crossSectionsFractionalUncertaintyDictionary[gluinoMass] = crossSectionFractionalUncertainty
crossSectionsInputFileObject.close()

print("Read in cross-sections as a function of gluino mass:")
prettyPrintDictionary(crossSectionsDictionary)

def getGeneratedMasses(inputMCChain):
    generatedMasses = {"gluino": 0., "neutralino": 0.}
    gluinoMassIsSet = False
    neutralinoMassIsSet = False
    for genParticleIndex in range(inputMCChain.nMC):
        if (not(gluinoMassIsSet) and inputMCChain.mcPID[genParticleIndex] == MCPIDs["gluino"]):
            generatedMasses["gluino"] = inputMCChain.mcMass[genParticleIndex]
            gluinoMassIsSet = True
        if (not(neutralinoMassIsSet) and inputMCChain.mcMomPID[genParticleIndex] == MCPIDs["neutralino"]):
            generatedMasses["neutralino"] = inputMCChain.mcMomMass[genParticleIndex]
            neutralinoMassIsSet = True
        if (gluinoMassIsSet and neutralinoMassIsSet): break
    if not(gluinoMassIsSet): sys.exit("Unable to find gluino mass in an event!")
    if not(neutralinoMassIsSet): sys.exit("Unable to find neutralino mass in an event!")
    return generatedMasses

def getWeight(gluinoMass, neutralinoMass):
    try:
        weight = crossSectionsDictionary[int(0.5 + gluinoMass)]*inputArguments.totalIntegratedLuminosity/inputArguments.nGeneratedEventsPerBin
    except KeyError:
        sys.exit("Unable to find cross-section corresponding to following gluino mass: {mass}".format(mass=gluinoMass))
    return weight

def getHistogramTitle(histogramType, nJetsBin, zone):
    histogramTypeString = ""
    if (histogramType == "total"): histogramTypeString = "Total MC Events"
    elif (histogramType == "weighted"): histogramTypeString = "Weighted MC Events"
    elif (histogramType == "signalContamination"): histogramTypeString = "Signal contamination"
    else: sys.exit("Unknown histogram type: {histType}".format(histType=histogramType))

    nJetsString = ""
    if (nJetsBin < inputArguments.nJetsMax): nJetsString = "{nJetsBin} Jets".format(nJetsBin = nJetsBin)
    elif (nJetsBin == inputArguments.nJetsMax): nJetsString = "#geq {nJetsBin} Jets".format(nJetsBin = nJetsBin)
    else: sys.exit("Unknown nJets bin: {nJetsBin}".format(nJetsBin = nJetsBin))

    sTRangeString = ""
    if (zone == "norm"): sTRangeString = "{sTNormMin} < #it{{S}}_T < {sTNormMax}".format(sTNormMin = inputArguments.sTMin_normWindow, sTNormMax = inputArguments.sTMax_normWindow)
    elif (zone == "obs"): sTRangeString = "#it{{S}}_T > {sTNormMax}".format(sTNormMax = inputArguments.sTMax_normWindow)
    elif (zone == "sub"): sTRangeString = "{sTNormMax} < #it{{S}}_T < {sTStartMainRegion}".format(sTNormMax = inputArguments.sTMax_normWindow, sTStartMainRegion = inputArguments.sTStartMainRegion)
    elif (zone == "main"): sTRangeString = "#it{{S}}_T > {sTStartMainRegion}".format(sTStartMainRegion = inputArguments.sTStartMainRegion)
    else: sys.exit("Unknown zone: {zone}".format(zone=zone))

    axesLabelsString = ";m_{#tilde{#it{g}}};m_{#tilde{#it{#chi_{1}^{0}}}}"

    title = "{typeString}, {nJetsString}, {sTRangeString}{axesLabelsString}".format(typeString = histogramTypeString, nJetsString = nJetsString, sTRangeString = sTRangeString, axesLabelsString = axesLabelsString)
    return title

def createDataCard(templateFileName, outputFileName, nEvents_subordinateRegions, nEvents_mainRegions, fractionalError_subordinateRegions, fractionalError_mainRegions):
    replacementTuplesList = []
    for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
        if (nJetsBin == 2 or nJetsBin == 3): continue
        formattedString_MC_SUBORD = "{nSubRegion:<11.3f}".format(nSubRegion=nEvents_subordinateRegions[nJetsBin])
        replacementTuplesList.append(("MC_SUBORD_{nJetsBin}".format(nJetsBin=nJetsBin), formattedString_MC_SUBORD))
        formattedString_MC_MN_REG = "{nMainRegion:<11.3f}".format(nMainRegion=nEvents_mainRegions[nJetsBin])
        replacementTuplesList.append(("MC_MN_REG_{nJetsBin}".format(nJetsBin=nJetsBin), formattedString_MC_MN_REG))
        formattedString_SU_S = "{fESubRegion:5.3f}".format(fESubRegion=(1+fractionalError_subordinateRegions[nJetsBin]))
        replacementTuplesList.append(("SU_S{nJetsBin}".format(nJetsBin=nJetsBin), formattedString_SU_S))
        formattedString_SU_M = "{fEMainRegion:5.3f}".format(fEMainRegion=(1+fractionalError_mainRegions[nJetsBin]))
        replacementTuplesList.append(("SU_M{nJetsBin}".format(nJetsBin=nJetsBin), formattedString_SU_M))
    
    templateFile = open(templateFileName, 'r')
    outputFile = open(outputFileName, 'w')
    for line in templateFile:
        runningString = line.strip()
        nextString = runningString
        for replacementTuple in replacementTuplesList:
            nextString = runningString.replace(replacementTuple[0], replacementTuple[1], 1)
            runningString = nextString
        outputFile.write(runningString + "\n")
    outputFile.close()
    templateFile.close()

sw = ROOT.TStopwatch()
sw.Start()

print("Analyzing data sample...")
inputDataChain = ROOT.TChain('ggNtuplizer/EventTree')
inputDataChain.Add(inputArguments.inputDataPath)
nDataEntries = inputDataChain.GetEntries()
print ("Total number of available events in data: {nDataEntries}".format(nDataEntries=nDataEntries))

nEventsInData = {
    "norm": {},
    "obs": {},
    "sub": {},
    "main": {}
}
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    nEventsInData["norm"][nJetsBin] = 0
    nEventsInData["obs"][nJetsBin] = 0
    nEventsInData["sub"][nJetsBin] = 0
    nEventsInData["main"][nJetsBin] = 0

progressBar = tmProgressBar(nDataEntries)
progressBarUpdatePeriod = max(1, nDataEntries//1000)
progressBar.initializeTimer()
for entryIndex in range(nDataEntries):
    entryStatus = inputDataChain.LoadTree(entryIndex)
    if entryStatus < 0: sys.exit("Tree failed to load entry at index {entryIndex}; returned status {entryStatus}".format(entryIndex=entryIndex, entryStatus=entryStatus))
    chainStatus = inputDataChain.GetEntry(entryIndex)
    if chainStatus <= 0: sys.exit("Unable to load data from chain at index {entryIndex}; GetEntry returned status {chainStatus}".format(entryIndex=entryIndex, chainStatus=chainStatus))

    if (entryIndex%progressBarUpdatePeriod == 0): progressBar.updateBar(1.0*entryIndex/nDataEntries, entryIndex)

    nStealthJets = inputDataChain.b_nJets
    nJetsBin = nStealthJets
    if (nJetsBin > inputArguments.nJetsMax): nJetsBin = inputArguments.nJetsMax
    if (nJetsBin < inputArguments.nJetsMin): sys.exit("Unexpected nJetsBin = {nJetsBin} at entry index = {entryIndex}".format(nJetsBin=nJetsBin, entryIndex=entryIndex))
    sT = inputDataChain.b_evtST
    if (sT > inputArguments.sTMax_normWindow and sT <= inputArguments.sTStartMainRegion): nEventsInData["sub"][nJetsBin] += 1
    if (sT > inputArguments.sTStartMainRegion): nEventsInData["main"][nJetsBin] += 1

    if (sT > inputArguments.sTMax_normWindow): nEventsInData["obs"][nJetsBin] += 1
    elif (sT > inputArguments.sTMin_normWindow): nEventsInData["norm"][nJetsBin] += 1
progressBar.terminate()

print("Number of events in data in norm sT range:")
prettyPrintDictionary(inputDict=nEventsInData["norm"], keyPrintOrder=list(range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax)))
# print("Number of events in data in observation sT range:")
# prettyPrintDictionary(inputDict=nEventsInData["obs"], keyPrintOrder=list(range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax)))

print("Analyzing MC sample...")
inputMCChain = ROOT.TChain('ggNtuplizer/EventTree')
inputMCChain.Add(inputArguments.inputMCPath)
nMCEntries = inputMCChain.GetEntries()
print ("Total number of available events in MC samples: {nMCEntries}".format(nMCEntries=nMCEntries))

if not(inputArguments.maxMCEvents == 0):
    nMCEntries = inputArguments.maxMCEvents
    print("Limiting loop over MC entries to {n} MC events.".format(n = nMCEntries))

histograms_total_nMCEvents = {
    "norm": {},
    "obs": {},
    "sub": {},
    "main": {}
}
histograms_weighted_nMCEvents = {
    "norm": {},
    "obs": {},
    "sub": {},
    "main": {}
}
histograms_signalContamination = {
    "norm": {},
    "obs": {},
    "sub": {},
    "main": {}
}

for zone in ["norm", "obs", "sub", "main"]:
    for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
        histograms_total_nMCEvents[zone][nJetsBin] = ROOT.TH2F("h_total_nMCEvents_{nJetsBin}Jets_{zone}".format(nJetsBin=nJetsBin, zone=zone), getHistogramTitle("total", nJetsBin, zone), inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
        histograms_weighted_nMCEvents[zone][nJetsBin] = ROOT.TH2F("h_weighted_nMCEvents_{nJetsBin}Jets_{zone}".format(nJetsBin=nJetsBin, zone=zone), getHistogramTitle("weighted", nJetsBin, zone), inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
        histograms_weighted_nMCEvents[zone][nJetsBin].Sumw2()
        histograms_signalContamination[zone][nJetsBin] = ROOT.TH2F("h_signalContamination_{nJetsBin}Jets_{zone}".format(nJetsBin=nJetsBin, zone=zone), getHistogramTitle("signalContamination", nJetsBin, zone), inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)

progressBar = tmProgressBar(nMCEntries)
progressBarUpdatePeriod = max(1, nMCEntries//1000)
progressBar.initializeTimer()
for entryIndex in range(nMCEntries):
    entryStatus = inputMCChain.LoadTree(entryIndex)
    if entryStatus < 0: sys.exit("Tree failed to load entry at index {entryIndex}; returned status {entryStatus}".format(entryIndex=entryIndex, entryStatus=entryStatus))
    chainStatus = inputMCChain.GetEntry(entryIndex)
    if chainStatus <= 0: sys.exit("Unable to load MC from chain at index {entryIndex}; GetEntry returned status {chainStatus}".format(entryIndex=entryIndex, chainStatus=chainStatus))

    if (entryIndex%progressBarUpdatePeriod == 0): progressBar.updateBar(1.0*entryIndex/nMCEntries, entryIndex)

    nStealthJets = inputMCChain.b_nJets
    nJetsBin = nStealthJets
    if (nJetsBin > inputArguments.nJetsMax): nJetsBin = inputArguments.nJetsMax
    if (nJetsBin < inputArguments.nJetsMin): sys.exit("Unexpected nJetsBin = {nJetsBin} at entry index = {entryIndex}".format(nJetsBin=nJetsBin, entryIndex=entryIndex))
    sT = inputMCChain.b_evtST
    zonesToFill = []

    if (sT > inputArguments.sTMax_normWindow and sT <= inputArguments.sTStartMainRegion): zonesToFill.append("sub")
    if (sT > inputArguments.sTStartMainRegion): zonesToFill.append("main")

    if sT > inputArguments.sTMax_normWindow: zonesToFill.append("obs")
    elif sT > inputArguments.sTMin_normWindow: zonesToFill.append("norm")

    if (len(zonesToFill) == 0): continue

    generatedMasses = getGeneratedMasses(inputMCChain)
    generated_gluinoMass = generatedMasses["gluino"]
    generated_neutralinoMass = generatedMasses["neutralino"]
    eventWeight = getWeight(generatedMasses["gluino"], generatedMasses["neutralino"])
    for zone in zonesToFill:
        histograms_total_nMCEvents[zone][nJetsBin].Fill(generated_gluinoMass, generated_neutralinoMass, 1.)
        histograms_weighted_nMCEvents[zone][nJetsBin].Fill(generated_gluinoMass, generated_neutralinoMass, eventWeight)
        histograms_signalContamination[zone][nJetsBin].Fill(generated_gluinoMass, generated_neutralinoMass, eventWeight/nEventsInData[zone][nJetsBin])
progressBar.terminate()

generatedMCTemplate = ROOT.TFile(inputArguments.MCTemplate)
h_MCTemplate = generatedMCTemplate.Get("h_susyMasses_template")
for gluinoMassBin in range(1, 1+h_MCTemplate.GetXaxis().GetNbins()):
    for neutralinoMassBin in range(1, 1+h_MCTemplate.GetYaxis().GetNbins()):
        if not(int(0.5 + h_MCTemplate.GetBinContent(gluinoMassBin, neutralinoMassBin)) == 1): continue
        gluinoMass = int(0.5 + h_MCTemplate.GetXaxis().GetBinCenter(gluinoMassBin))
        neutralinoMass = int(0.5 + h_MCTemplate.GetYaxis().GetBinCenter(neutralinoMassBin))
        MCFractionalError_subordinateRegions = {}
        nMCEvents_subordinateRegions = {}
        MCFractionalError_mainRegions = {}
        nMCEvents_mainRegions = {}
        for nJetsBin in range(inputArguments.nJetsMin, 1+inputArguments.nJetsMax):
            if (nJetsBin == 2 or nJetsBin == 3): continue

            MCEventsAndErrors_subordinateRegion = tmROOTUtils.get2DHistogramContentAndErrorAtCoordinates(inputTH2=histograms_total_nMCEvents["sub"][nJetsBin], xValue=1.0*gluinoMass, yValue=1.0*neutralinoMass)
            if (MCEventsAndErrors_subordinateRegion["content"] == 0):
                print("WARNING:  at gluino mass = {gM}, neutralino mass={nM}, nJetsBin = {nJetsBin}, total number of MC events is 0!".format(gM=gluinoMass, nM=neutralinoMass, nJetsBin=nJetsBin))
                MCFractionalError_subordinateRegions[nJetsBin] = 0.5
                nMCEvents_subordinateRegions[nJetsBin] = 0.
            else:
                MCFractionalError_subordinateRegions[nJetsBin] = MCEventsAndErrors_subordinateRegion["error"]/MCEventsAndErrors_subordinateRegion["content"]
                weightedMCEventsAndErrors_subordinateRegion = tmROOTUtils.get2DHistogramContentAndErrorAtCoordinates(inputTH2=histograms_weighted_nMCEvents["sub"][nJetsBin], xValue=1.0*gluinoMass, yValue=1.0*neutralinoMass)
                nMCEvents_subordinateRegions[nJetsBin] = weightedMCEventsAndErrors_subordinateRegion["content"]

            MCEventsAndErrors_mainRegion = tmROOTUtils.get2DHistogramContentAndErrorAtCoordinates(inputTH2=histograms_total_nMCEvents["main"][nJetsBin], xValue=1.0*gluinoMass, yValue=1.0*neutralinoMass)
            if (MCEventsAndErrors_mainRegion["content"] == 0):
                print("WARNING:  at gluino mass = {gM}, neutralino mass={nM}, nJetsBin = {nJetsBin}, total number of MC events is 0!".format(gM=gluinoMass, nM=neutralinoMass, nJetsBin=nJetsBin))
                MCFractionalError_mainRegions[nJetsBin] = 0.5
                nMCEvents_mainRegions[nJetsBin] = 0.
            else:
                MCFractionalError_mainRegions[nJetsBin] = MCEventsAndErrors_mainRegion["error"]/MCEventsAndErrors_mainRegion["content"]
                weightedMCEventsAndErrors_mainRegion = tmROOTUtils.get2DHistogramContentAndErrorAtCoordinates(inputTH2=histograms_weighted_nMCEvents["main"][nJetsBin], xValue=1.0*gluinoMass, yValue=1.0*neutralinoMass)
                nMCEvents_mainRegions[nJetsBin] = weightedMCEventsAndErrors_mainRegion["content"]

        createDataCard(inputArguments.dataCardTemplate, "analysis/dataCards/dataCard_gluinoMass_{gM}_neutralinoMass_{nM}.txt".format(gM=gluinoMass, nM=neutralinoMass), nMCEvents_subordinateRegions, nMCEvents_mainRegions, MCFractionalError_subordinateRegions, MCFractionalError_mainRegions)

outputFile = ROOT.TFile("analysis/signalContamination/{prefix}_savedObjects.root".format(prefix=inputArguments.outputPrefix), "RECREATE")
for zone in ["norm", "obs", "sub", "main"]:
    for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
        tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [histograms_total_nMCEvents[zone][nJetsBin]], canvasName = "c_total_nMCEvents_{nJetsBin}Jets_{zone}".format(nJetsBin=nJetsBin, zone=zone), outputROOTFile=outputFile, outputDocumentName = "analysis/signalContamination/{outputPrefix}_total_nEvents_{nJetsBin}Jets_{zone}".format(outputPrefix=inputArguments.outputPrefix, nJetsBin=nJetsBin, zone=zone), customOptStat=0, customTextFormat=".0f", customPlotOptions_firstObject="TEXTCOLZ", enableLogZ = True)
        tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [histograms_weighted_nMCEvents[zone][nJetsBin]], canvasName = "c_weighted_nMCEvents_{nJetsBin}Jets_{zone}".format(nJetsBin=nJetsBin, zone=zone), outputROOTFile=outputFile, outputDocumentName = "analysis/signalContamination/{outputPrefix}_weighted_nEvents_{nJetsBin}Jets_{zone}".format(outputPrefix=inputArguments.outputPrefix, nJetsBin=nJetsBin, zone=zone), customOptStat=0, customTextFormat=".0f", customPlotOptions_firstObject="TEXTCOLZ", enableLogZ = True)
        if (zone == "obs" and nJetsBin == inputArguments.nJetsMax): tmROOTUtils.extractTH2Contents(histograms_weighted_nMCEvents[zone][nJetsBin], ("analysis/signalContamination/{outputPrefix}_weighted_nEvents_{nJetsBin}Jets_{zone}_extracted.txt").format(outputPrefix=inputArguments.outputPrefix, nJetsBin=nJetsBin, zone=zone), quantityName = "Weighted number of MC events", includeOverflow = True, formatSpecifiers = ["%.1f", "%.1f", "%.3f"])
        if (not(zone == "norm") and not(nJetsBin in nJetsBinsToAnalyze)): continue
        tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [histograms_signalContamination[zone][nJetsBin]], canvasName = "c_signalContamination_{nJetsBin}Jets_{zone}".format(nJetsBin=nJetsBin, zone=zone), outputROOTFile=outputFile, outputDocumentName = "analysis/signalContamination/{outputPrefix}_signalContamination_{nJetsBin}Jets_{zone}".format(outputPrefix=inputArguments.outputPrefix, nJetsBin=nJetsBin, zone=zone), customOptStat=0, customTextFormat=".3f", customPlotOptions_firstObject="COLZ", enableLogZ = True)
outputFile.Close()
