#!/usr/bin/env python

import ROOT, tmROOTUtils, sys, argparse, os
from tmProgressBar import tmProgressBar

ROOT.gROOT.SetBatch(ROOT.kTRUE)

inputArgumentsParser = argparse.ArgumentParser(description='Plot generated gluino and neutralino masses from a set of generated MC files.')
inputArgumentsParser.add_argument('--inputFilesList', required=True, help="Path to file containing list of input files.", type=str)
inputArgumentsParser.add_argument('--outputFolder', default="MCGeneratedMasses", help='Output folder.',type=str)
inputArgumentsParser.add_argument('--outputPrefix', default="", required=True, help='Output prefix.',type=str)
inputArgumentsParser.add_argument('--prescale', default=1, help='Prescale on number of events. Default: no prescale',type=int)
inputArgumentsParser.add_argument('--nGluinoMassBins', default=28, help='nBins on the gluino mass axis.',type=int) # (1000 - 25) GeV --> (2350 + 25) GeV in steps of 50 GeV
inputArgumentsParser.add_argument('--minGluinoMass', default=975., help='Min gluino mass.',type=float)
inputArgumentsParser.add_argument('--maxGluinoMass', default=2375., help='Max gluino mass.',type=float)
inputArgumentsParser.add_argument('--nNeutralinoMassBins', default=181, help='nBins on the neutralino mass axis.',type=int) # (100 - 6.25) GeV --> (2350 + 6.25) GeV in steps of 12.5 GeV
inputArgumentsParser.add_argument('--minNeutralinoMass', default=93.75, help='Min neutralino mass.',type=float)
inputArgumentsParser.add_argument('--maxNeutralinoMass', default=2356.25, help='Max neutralino mass.',type=float)
inputArguments = inputArgumentsParser.parse_args()

os.system("mkdir -p {oF}".format(oF=inputArguments.outputFolder))

MCPIDs = {
    "photon": 22,
    "gluino": 1000021,
    "neutralino": 1000022
}

prescaleInverse = 1./inputArguments.prescale
prescaleEnabled = not(inputArguments.prescale == 1)

h_gluinoMass = ROOT.TH1F("h_gluinoMass", ";m_{#tilde{#it{g}}};Total nEvents", inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass)
h_neutralinoMass = ROOT.TH1F("h_neutralinoMass", ";m_{#tilde{#it{#chi_{1}^{0}}}};Total nEvents", inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)
h_masses = ROOT.TH2F("h_masses", "Total nEvents;m_{#tilde{#it{g}}};m_{#tilde{#it{#chi_{1}^{0}}}}", inputArguments.nGluinoMassBins, inputArguments.minGluinoMass, inputArguments.maxGluinoMass, inputArguments.nNeutralinoMassBins, inputArguments.minNeutralinoMass, inputArguments.maxNeutralinoMass)

def fillHistogramsFromFile(fileName, histGluinoMass, histNeutralinoMass, hist2D, randomNumberGenerator):
    global nEvents_neutralinoMassUnset, nEvents_analyzed, minGluinoMassFound, maxGluinoMassFound, minNeutralinoMassFound, maxNeutralinoMassFound
    inputMCChain = ROOT.TChain('ggNtuplizer/EventTree')
    inputMCChain.Add(fileName)
    nMCEntries = inputMCChain.GetEntries()
    print ("Number of events in file {f}: {nMCEntries}".format(f=fileName, nMCEntries=nMCEntries))
    for entryIndex in range(nMCEntries):
        if prescaleEnabled:
            randomNumber = randomNumberGenerator.Uniform()
            if randomNumber > prescaleInverse: continue
        nEvents_analyzed += 1
        entryStatus = inputMCChain.LoadTree(entryIndex)
        if entryStatus < 0: sys.exit("Tree failed to load entry at index {entryIndex}; returned status {entryStatus}".format(entryIndex=entryIndex, entryStatus=entryStatus))
        chainStatus = inputMCChain.GetEntry(entryIndex)
        if chainStatus <= 0: sys.exit("Unable to load MC from chain at index {entryIndex}; GetEntry returned status {chainStatus}".format(entryIndex=entryIndex, chainStatus=chainStatus))
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
        if gluinoMassIsSet:
            h_gluinoMass.Fill(generatedMasses["gluino"])
        else:
            sys.exit("Gluino mass unset in event with index {index}".format(index=entryIndex))
        if neutralinoMassIsSet:
            h_neutralinoMass.Fill(generatedMasses["neutralino"])
            h_masses.Fill(generatedMasses["gluino"], generatedMasses["neutralino"])
        else:
            nEvents_neutralinoMassUnset += 1 # sys.exit("Unable to find neutralino mass in an event!")
        if ((minGluinoMassFound < 0) or (generatedMasses["gluino"] < minGluinoMassFound)):
            minGluinoMassFound = generatedMasses["gluino"]
        if ((maxGluinoMassFound < 0) or (generatedMasses["gluino"] > maxGluinoMassFound)):
            maxGluinoMassFound = generatedMasses["gluino"]
        if ((minNeutralinoMassFound < 0) or (generatedMasses["neutralino"] < minNeutralinoMassFound)):
            minNeutralinoMassFound = generatedMasses["neutralino"]
        if ((maxNeutralinoMassFound < 0) or (generatedMasses["neutralino"] > maxNeutralinoMassFound)):
            maxNeutralinoMassFound = generatedMasses["neutralino"]
        if ((generatedMasses["gluino"] < inputArguments.minGluinoMass) or
            (generatedMasses["gluino"] > inputArguments.maxGluinoMass) or
            (generatedMasses["neutralino"] > inputArguments.maxNeutralinoMass)):
            sys.exit("ERROR: Unexpected (gluino, neutralino) masses: {g,n}".format(g=generatedMasses["gluino"], n=generatedMasses["neutralino"]))

randomNumberGenerator = ROOT.TRandom2()
randomNumberGenerator.SetSeed(100)

listOfInputFiles = []
inputFileNamesFileObject = open(inputArguments.inputFilesList, 'r')
for inputFileName in inputFileNamesFileObject:
    listOfInputFiles.append(inputFileName.strip())
inputFileNamesFileObject.close()

nEvents_neutralinoMassUnset = 0
nEvents_analyzed = 0
minGluinoMassFound = -1
maxGluinoMassFound = -1
minNeutralinoMassFound = -1
maxNeutralinoMassFound = -1

for fname in listOfInputFiles:
    fillHistogramsFromFile(fname, h_gluinoMass, h_neutralinoMass, h_masses, randomNumberGenerator)

print(("Min gluino mass found: {minGMF}, Max gluino mass found: {maxGMF}, Min neutralino mass found: {minNMF}, Max neutralino mass found: {maxNMF}").format(minGMF=minGluinoMassFound, maxGMF=maxGluinoMassFound, minNMF=minNeutralinoMassFound, maxNMF=maxNeutralinoMassFound))

outputFile = ROOT.TFile.Open("{oF}/MCGeneratedMasses_{p}_savedObjects.root".format(oF=inputArguments.outputFolder, p=inputArguments.outputPrefix), "RECREATE")
c_gluinoMass = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [h_gluinoMass], canvasName = "c_gluinoMass", outputROOTFile=outputFile, outputDocumentName="{oF}/gluinoMassDistribution_{p}".format(oF=inputArguments.outputFolder, p=inputArguments.outputPrefix))
c_neutralinoMass = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [h_neutralinoMass], canvasName = "c_neutralinoMass", outputROOTFile=outputFile, outputDocumentName="{oF}/neutralinoMassDistribution_{p}".format(oF=inputArguments.outputFolder, p=inputArguments.outputPrefix))
c_masses = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [h_masses], canvasName = "c_masses", outputROOTFile=outputFile, outputDocumentName="{oF}/generatedMassesDistribution_{p}".format(oF=inputArguments.outputFolder, p=inputArguments.outputPrefix), customOptStat=0, customPlotOptions_firstObject="TEXTCOLZ")
outputFile.Close()

print("nEvents with neutralino mass unset: {n}, total nEvents analyzed = {tot}".format(n=nEvents_neutralinoMassUnset, tot=nEvents_analyzed))
