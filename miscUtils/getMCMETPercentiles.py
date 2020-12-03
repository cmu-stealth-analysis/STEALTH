#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, tmROOTUtils, pdb, numpy

from tmProgressBar import tmProgressBar
from tmGeneralUtils import prettyPrintDictionary

inputArgumentsParser = argparse.ArgumentParser(description='Get MET percentiles from input MC.')
inputArgumentsParser.add_argument('--inputMCPath', required=True, help='Path to input MC file.',type=str)
inputArgumentsParser.add_argument('--sTMin', default=1100., help='Min value of sT.',type=float)
inputArgumentsParser.add_argument('--outputPrefix', required=True, help='Prefix to output files.',type=str)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Minimum number of jets in event.',type=str)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Least value of nJets in highest nJets bin.',type=str)
inputArgumentsParser.add_argument('--analyze_nJetsBin', action='append', default=[], help='nJets to plot.',type=int)
inputArgumentsParser.add_argument('--nGluinoMassBins', default=20, help='nBins on the gluino mass axis.',type=int) # (800 - 25) GeV --> (1750 + 25) GeV in steps of 50 GeV
inputArgumentsParser.add_argument('--minGluinoMass', default=775., help='Min gluino mass.',type=float)
inputArgumentsParser.add_argument('--maxGluinoMass', default=1775., help='Max gluino mass.',type=float)
inputArgumentsParser.add_argument('--nNeutralinoMassBins', default=133, help='nBins on the neutralino mass axis.',type=int)
inputArgumentsParser.add_argument('--minNeutralinoMass', default=93.75, help='Min neutralino mass.',type=float)
inputArgumentsParser.add_argument('--maxNeutralinoMass', default=1756.25, help='Max neutralino mass.',type=float) # (100 - 6.25) GeV --> (1750 + 6.25) GeV in steps of 12.5 GeV
inputArgumentsParser.add_argument('--addPercentile', action='append', help='MET percentile to generate.',type=float)
inputArgumentsParser.add_argument('--addParticularMETDistribution', action='append', help='Get 1-D MET distribution particularly in this (gluino mass, neutralino mass) bin. Gluino and neutralino masses should be comma-separated.',type=float)
inputArguments = inputArgumentsParser.parse_args()

nJetsBinsToAnalyze=inputArguments.analyze_nJetsBin
for nJetsBinToAnalyze in nJetsBinsToAnalyze:
    if nJetsBinToAnalyze > inputArguments.nJetsMax: sys.exit("Error: index of jets bin to analyze can be at most nJetsMax={nJetsMax}; this argument is not allowed: --analyze_nJetsBin={bin}".format(nJetsMax=inputArguments.nJetsMax, bin=nJetsBinToAnalyze))

percentilesToGenerate = inputArguments.addPercentile

MCPIDs = {
    "photon": 22,
    "gluino": 1000021,
    "neutralino": 1000022
}

gluinoAxis = ROOT.TAxis(inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass)
neutralinoAxis = ROOT.TAxis(inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
particularIDsList = []
particularMETDistributionStrings = inputArguments.addParticularMETDistributions
for particularMETDistributionString in particularMETDistributionStrings:
    particularMETDistributionString_split = particularMETDistributionString.split(',')
    if not(len(particularMETDistributionString_split) == 2): sys.exit("Error: particular MET distribution must specify exactly two values separated by a comma: gluino and neutralino mass.")
    particularGluinoID = int(0.5 + gluinoAxis.GetBinCenter(gluinoAxis.FindBin(float(particularMETDistributionString_split[0]))))
    particularNeutralinoID = int(0.5 + neutralinoAxis.GetBinCenter(neutralinoAxis.FindBin(float(particularMETDistributionString_split[1]))))
    particularIDsList.append(tuple([particularGluinoID, particularNeutralinoID]))

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

def getHistogramTitle(percentile, nJetsBin):
    percentileString = "MC MET {pc:2.2f} percentile".format(pc=percentile)
    
    nJetsString = ""
    if (nJetsBin < inputArguments.nJetsMax): nJetsString = "{nJetsBin} Jets".format(nJetsBin = nJetsBin)
    elif (nJetsBin == inputArguments.nJetsMax): nJetsString = "#geq {nJetsBin} Jets".format(nJetsBin = nJetsBin)
    else: sys.exit("Unknown nJets bin: {nJetsBin}".format(nJetsBin = nJetsBin))

    axesLabelsString = ";m_{#tilde{#it{g}}};m_{#tilde{#it{#chi_{1}^{0}}}}"
    title = "{percentileString}, {nJetsString}{axesLabelsString}".format(percentileString = percentileString, nJetsString = nJetsString, axesLabelsString = axesLabelsString)
    return title

sw = ROOT.TStopwatch()
sw.Start()

METArraysDictionaries = {}
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    METArraysDictionaries[nJetsBin] = {}

print("Analyzing MC sample...")
inputMCChain = ROOT.TChain('ggNtuplizer/EventTree')
inputMCChain.Add(inputArguments.inputMCPath)
nMCEntries = inputMCChain.GetEntries()
print ("Total number of available events in MC samples: {nMCEntries}".format(nMCEntries=nMCEntries))

particularHistograms = {}
for particularIDsCounter in range(len(particularIDsList)):
    particularIDs = particularIDsList[particularIDsCounter]
    gluinoMass = particularIDs[0]
    neutralinoMass = particularIDs[1]
    histogramTitle = "MET distribution, m_{#tilde{#it{g}}} = " + str(gluinoMass) + ", m_{#tilde{#it{#chi_{1}^{0}}}} = " + str(neutralinoMass) + ";MET(GeV);nEvents"
    particularHistograms[particularIDsCounter] = ROOT.TH1F("h_METDistribution_{counter}".format(counter = particularIDsCounter), histogramTitle, 80, 0., 0.)

progressBar = tmProgressBar(nMCEntries)
progressBarUpdatePeriod = max(1, nMCEntries//1000)
progressBar.initializeTimer()
for entryIndex in range(nMCEntries):
    entryStatus = inputMCChain.LoadTree(entryIndex)
    if entryStatus < 0: sys.exit("Tree failed to load entry at index {entryIndex}; returned status {entryStatus}".format(entryIndex=entryIndex, entryStatus=entryStatus))
    chainStatus = inputMCChain.GetEntry(entryIndex)
    if chainStatus <= 0: sys.exit("Unable to load MC from chain at index {entryIndex}; GetEntry returned status {chainStatus}".format(entryIndex=entryIndex, chainStatus=chainStatus))

    if (entryIndex%progressBarUpdatePeriod == 0): progressBar.updateBar(1.0*entryIndex/nMCEntries, entryIndex)

    nStealthJets = inputMCChain.b_nJetsAll
    nJetsBin = nStealthJets
    if (nJetsBin > inputArguments.nJetsMax): nJetsBin = inputArguments.nJetsMax
    if (nJetsBin < inputArguments.nJetsMin): sys.exit("Unexpected nJetsBin = {nJetsBin} at entry index = {entryIndex}".format(nJetsBin=nJetsBin, entryIndex=entryIndex))
    if (not(nJetsBin in nJetsBinsToAnalyze)): continue
    METArraysDictionary = METArraysDictionaries[nJetsBin]
    sT = inputMCChain.b_evtST
    if sT < inputArguments.sTMin: continue
    generatedMasses = getGeneratedMasses(inputMCChain)
    generated_gluinoMass = generatedMasses["gluino"]
    gluinoMassID = int(0.5 + gluinoAxis.GetBinCenter(gluinoAxis.FindBin(generated_gluinoMass)))
    if not(gluinoMassID in METArraysDictionary):
        METArraysDictionary[gluinoMassID] = {}
    generated_neutralinoMass = generatedMasses["neutralino"]
    neutralinoMassID = int(0.5 + neutralinoAxis.GetBinCenter(neutralinoAxis.FindBin(generated_neutralinoMass)))
    if not(neutralinoMassID in METArraysDictionary[gluinoMassID]):
        METArraysDictionary[gluinoMassID][neutralinoMassID] = []
    MET = inputMCChain.pfMET
    METArraysDictionary[gluinoMassID][neutralinoMassID].append(MET)
progressBar.terminate()

for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    if (not(nJetsBin in nJetsBinsToAnalyze)): continue
    histograms_metPercentile = {}
    for percentileCounter in range(len(percentilesToGenerate)):
        histograms_metPercentile[percentileCounter] = ROOT.TH2F(("h_metPercentile_{pc:2.2f}_{nJetsBin}Jets".format(pc=percentilesToGenerate[percentileCounter], nJetsBin=nJetsBin)).replace('.', 'pt'), getHistogramTitle(percentilesToGenerate[percentileCounter], nJetsBin), inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
    METArraysDictionary = METArraysDictionaries[nJetsBin]
    for gluinoMassID in METArraysDictionary.keys():
        for neutralinoMassID in METArraysDictionary[gluinoMassID].keys():
            MET_numpyArray = numpy.array(METArraysDictionary[gluinoMassID][neutralinoMassID])
            percentileValues = numpy.percentile(MET_numpyArray, percentilesToGenerate, overwrite_input=True)
            for percentileCounter in range(len(percentilesToGenerate)):
                histograms_metPercentile[percentileCounter].SetBinContent(histograms_metPercentile[percentileCounter].FindBin(gluinoMassID, neutralinoMassID), percentileValues[percentileCounter])
            METArraysDictionary[gluinoMassID].pop(neutralinoMassID)
        METArraysDictionary.pop(gluinoMassID)
    METArraysDictionaries.pop(nJetsBin)
    for percentileCounter in range(len(percentilesToGenerate)):
        tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [histograms_metPercentile[percentileCounter]], canvasName = ("c_metPercentile_{pc:2.2f}_{nJetsBin}Jets".format(pc=percentilesToGenerate[percentileCounter], nJetsBin=nJetsBin)).replace('.', 'pt'), outputDocumentName = ("analysis/MCMETPercentiles/{outputPrefix}_METPercentiles_{pc:2.2f}_percentile_{nJetsBin}Jets".format(outputPrefix = inputArguments.outputPrefix, pc=percentilesToGenerate[percentileCounter], nJetsBin=nJetsBin)).replace('.', 'pt'), customOptStat=0, customPlotOptions_firstObject="COLZ", enableLogZ = True)
