#!/usr/bin/env python

import ROOT, tmROOTUtils, sys, argparse
from tmProgressBar import tmProgressBar

inputArgumentsParser = argparse.ArgumentParser(description='Plot generated gluino and neutralino masses from a set of generated MC files.')
inputArgumentsParser.add_argument('--inputDataPath', required=True, help='List of files in the format of a string to pass to chain.Add().',type=str) # e.g. "root://cmseos.fnal.gov//store/user/lpcsusystealth/stealth2018Ntuples/job_SMS-T7WgStealth/SMS-T7WgStealth_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_job_SMS-T7WgStealth/180525_031712/0000/ggtree_mc_*.root"
inputArgumentsParser.add_argument('--outputPrefix', required=True, help='Output prefix.', type=str)
inputArgumentsParser.add_argument('--prescale', default=1, help='Prescale on number of events. Default: no prescale',type=int)
inputArgumentsParser.add_argument('--nGluinoMassBins', default=20, help='nBins on the gluino mass axis.',type=int)
inputArgumentsParser.add_argument('--minGluinoMass', default=775., help='Min gluino mass.',type=float)
inputArgumentsParser.add_argument('--maxGluinoMass', default=1775., help='Max gluino mass.',type=float)
inputArgumentsParser.add_argument('--nNeutralinoMassBins', default=133, help='nBins on the neutralino mass axis.',type=int)
inputArgumentsParser.add_argument('--minNeutralinoMass', default=93.75, help='Min neutralino mass.',type=float)
inputArgumentsParser.add_argument('--maxNeutralinoMass', default=1756.25, help='Max neutralino mass.',type=float)
inputArguments = inputArgumentsParser.parse_args()

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

inputMCChain = ROOT.TChain('ggNtuplizer/EventTree')
inputMCChain.Add(inputArguments.inputDataPath)
nMCEntries = inputMCChain.GetEntries()
print ("Total number of available events in MC: {nMCEntries}".format(nMCEntries=nMCEntries))

randomNumberGenerator = ROOT.TRandom2()
randomNumberGenerator.SetSeed(100)

nEvents_neutralinoMassUnset = 0
nEvents_analyzed = 0
progressBar = tmProgressBar(nMCEntries)
progressBarUpdatePeriod = max(1, nMCEntries//1000)
progressBar.initializeTimer()
for entryIndex in range(nMCEntries):
    if (entryIndex%progressBarUpdatePeriod == 0): progressBar.updateBar(1.0*entryIndex/nMCEntries, entryIndex)
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
progressBar.terminate()

outputFile = ROOT.TFile("analysis/MCMassDistributions/{prefix}_savedObjects.root".format(prefix=inputArguments.outputPrefix), "RECREATE")
c_gluinoMass = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [h_gluinoMass], canvasName = "c_gluinoMass", outputROOTFile=outputFile, outputDocumentName="analysis/MCMassDistributions/{prefix}gluinoMassDistribution".format(prefix=inputArguments.outputPrefix))
c_neutralinoMass = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [h_neutralinoMass], canvasName = "c_neutralinoMass", outputROOTFile=outputFile, outputDocumentName="analysis/MCMassDistributions/{prefix}neutralinoMassDistribution".format(prefix=inputArguments.outputPrefix))
c_masses = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [h_masses], canvasName = "c_masses", outputROOTFile=outputFile, outputDocumentName="analysis/MCMassDistributions/{prefix}generatedMassesDistribution".format(prefix=inputArguments.outputPrefix), customOptStat=0, customPlotOptions_firstObject="TEXTCOLZ")
outputFile.Close()

print("nEvents with neutralino mass unset: {n}, total nEvents analyzed = {tot}".format(n=nEvents_neutralinoMassUnset, tot=nEvents_analyzed))
