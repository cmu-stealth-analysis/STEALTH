#!/usr/bin/env python

import ROOT, tmROOTUtils, sys, argparse, os

ROOT.gROOT.SetBatch(ROOT.kTRUE)

inputArgumentsParser = argparse.ArgumentParser(description='Plot generated gluino and neutralino masses from a set of generated MC files.')
inputArgumentsParser.add_argument('--inputFilesList', required=True, help="Path to file containing list of input files.", type=str)
inputArgumentsParser.add_argument('--outputFolder', default="MCGeneratedMasses", help='Output folder.',type=str)
inputArgumentsParser.add_argument('--outputPrefix', default="", required=True, help='Output prefix.',type=str)
inputArgumentsParser.add_argument('--prescale', default=1, help='Prescale on number of events. Default: no prescale',type=int)
inputArgumentsParser.add_argument('--MCType', required=True, help="Type of MC input. Currently supported: disquark production (\"squark\") and digluino production(\"gluino\").", type=str)
inputArguments = inputArgumentsParser.parse_args()

os.system("mkdir -p {oF}".format(oF=inputArguments.outputFolder))

MCPIDs = {
    "photon": 22,
    "gluino": 1000021,
    "neutralino": 1000022,
    "sdown": 1000001,
    "sup": 1000002,
    "sstrange": 1000003,
    "scharm": 1000004,
    "sbottom": 1000005,
    "stop": 1000006
}

nGluinoMassBins = -1
minGluinoMass = -1.
maxGluinoMass = -1.
nSquarkMassBins = -1
minSquarkMass = -1.
maxSquarkMass = -1.
nNeutralinoMassBins = -1
minNeutralinoMass = -1.
maxNeutralinoMass = -1.
if (inputArguments.MCType == "gluino"):
    nGluinoMassBins = 28 # (1000 - 25) GeV --> (2350 + 25) GeV in steps of 50 GeV
    minGluinoMass = 975.
    maxGluinoMass = 2375.
    nNeutralinoMassBins = 181 # (100 - 6.25) GeV --> (2350 + 6.25) GeV in steps of 12.5 GeV
    minNeutralinoMass = 93.75
    maxNeutralinoMass = 2356.25
elif (inputArguments.MCType == "squark"):
    nSquarkMassBins = 25 # (850 - 25) GeV --> (2050 + 25) GeV in steps of 50 GeV
    minSquarkMass = 825.
    maxSquarkMass = 2075.
    nNeutralinoMassBins = 157 # (100 - 6.25) GeV --> (2050 + 6.25) GeV in steps of 12.5 GeV
    minNeutralinoMass = 93.75
    maxNeutralinoMass = 2056.25
else:
    sys.exit("ERROR: inputArguments.MCType has to be one of \"gluino\" and \"squark\". Currently, it is: {MT}".format(MT=inputArguments.MCType))

# For both of these, we check whether the input pid is in a range. I'm not sure whether checking for equality works because the int may not be typecast correctly.
def is_neutralino_pid(input_pid):
    return ((abs(input_pid) > MCPIDs["neutralino"] - 0.5) and (abs(input_pid) < MCPIDs["neutralino"] + 0.5))

def is_gluino_pid(input_pid):
    return ((abs(input_pid) > MCPIDs["gluino"] - 0.5) and (abs(input_pid) < MCPIDs["gluino"] + 0.5))

def is_squark_pid(input_pid):
    return ((abs(input_pid) > MCPIDs["sdown"] - 0.5) and (abs(input_pid) < MCPIDs["stop"] + 0.5))

prescaleInverse = 1./inputArguments.prescale
prescaleEnabled = not(inputArguments.prescale == 1)

h_gluinoMass = None
h_squarkMass = None
h_masses = None
h_neutralinoMass = ROOT.TH1F("h_neutralinoMass", ";m_{#tilde{#it{#chi_{1}^{0}}}};Total nEvents", nNeutralinoMassBins, minNeutralinoMass, maxNeutralinoMass)

if (inputArguments.MCType == "gluino"):
    h_gluinoMass = ROOT.TH1F("h_gluinoMass", ";m_{#tilde{#it{g}}};Total nEvents", nGluinoMassBins, minGluinoMass, maxGluinoMass)
    h_masses = ROOT.TH2F("h_masses", "Total nEvents;m_{#tilde{#it{g}}};m_{#tilde{#it{#chi_{1}^{0}}}}", nGluinoMassBins, minGluinoMass, maxGluinoMass, nNeutralinoMassBins, minNeutralinoMass, maxNeutralinoMass)
elif (inputArguments.MCType == "squark"):
    h_squarkMass = ROOT.TH1F("h_squarkMass", ";m_{#tilde{#it{q}}};Total nEvents", nSquarkMassBins, minSquarkMass, maxSquarkMass)
    h_masses = ROOT.TH2F("h_masses", "Total nEvents;m_{#tilde{#it{q}}};m_{#tilde{#it{#chi_{1}^{0}}}}", nSquarkMassBins, minSquarkMass, maxSquarkMass, nNeutralinoMassBins, minNeutralinoMass, maxNeutralinoMass)

def fillHistogramsFromFile(fileName, histGluinoMass, histSquarkMass, histNeutralinoMass, hist2D, randomNumberGenerator, MCType):
    global nEvents_gluinoMassUnset, nEvents_squarkMassUnset, nEvents_neutralinoMassUnset, nEvents_analyzed, minGluinoMassFound, maxGluinoMassFound, minSquarkMassFound, maxSquarkMassFound, minNeutralinoMassFound, maxNeutralinoMassFound
    inputMCChain = ROOT.TChain('ggNtuplizer/EventTree')
    inputMCChain.Add(fileName)
    nMCEntries = inputMCChain.GetEntries()
    print ("Number of events in file {f}: {nMCEntries}".format(f=fileName, nMCEntries=nMCEntries))
    inputIsGluino = (MCType == "gluino")
    inputIsSquark = (MCType == "squark")
    for entryIndex in range(nMCEntries):
        if prescaleEnabled:
            randomNumber = randomNumberGenerator.Uniform()
            if randomNumber > prescaleInverse: continue
        nEvents_analyzed += 1
        entryStatus = inputMCChain.LoadTree(entryIndex)
        if entryStatus < 0: sys.exit("Tree failed to load entry at index {entryIndex}; returned status {entryStatus}".format(entryIndex=entryIndex, entryStatus=entryStatus))
        chainStatus = inputMCChain.GetEntry(entryIndex)
        if chainStatus <= 0: sys.exit("Unable to load MC from chain at index {entryIndex}; GetEntry returned status {chainStatus}".format(entryIndex=entryIndex, chainStatus=chainStatus))
        generatedMasses = {"gluino": 0., "squark": 0., "neutralino": 0.}
        gluinoMassIsSet = False
        squarkMassIsSet = False
        neutralinoMassIsSet = False
        for genParticleIndex in range(inputMCChain.nMC):
            # if (abs(inputMCChain.mcPID[genParticleIndex]) > 1000000): print("Found possible SUSY particle, PID: {pid}".format(pid=inputMCChain.mcPID[genParticleIndex]))
            if inputIsGluino:
                if (not(gluinoMassIsSet) and is_gluino_pid(inputMCChain.mcPID[genParticleIndex])):
                    generatedMasses["gluino"] = inputMCChain.mcMass[genParticleIndex]
                    gluinoMassIsSet = True
            elif inputIsSquark:
                if (not(squarkMassIsSet) and is_squark_pid(inputMCChain.mcPID[genParticleIndex])):
                    generatedMasses["squark"] = inputMCChain.mcMass[genParticleIndex]
                    squarkMassIsSet = True
            if (not(neutralinoMassIsSet) and is_neutralino_pid(inputMCChain.mcMomPID[genParticleIndex])):
                generatedMasses["neutralino"] = inputMCChain.mcMomMass[genParticleIndex]
                neutralinoMassIsSet = True
            if ((gluinoMassIsSet or squarkMassIsSet) and neutralinoMassIsSet): break
        if inputIsGluino:
            if gluinoMassIsSet:
                h_gluinoMass.Fill(generatedMasses["gluino"])
            else:
                nEvents_gluinoMassUnset += 1
                #     sys.exit("Gluino mass unset in event with index {index}".format(index=entryIndex))
        elif inputIsSquark:
            if squarkMassIsSet:
                h_squarkMass.Fill(generatedMasses["squark"])
            else:
                nEvents_squarkMassUnset += 1
            # else: # Commenting out because for some strange reason the squark mass is unset in most events --
            #     # in fact it seems to be unset in all events with charginos. This is not a problem for the moment.
            #     # sys.exit("Squark mass unset in event with index {index}".format(index=entryIndex))
            #     print("Squark mass unset in event with index {index}".format(index=entryIndex))
        if neutralinoMassIsSet:
            h_neutralinoMass.Fill(generatedMasses["neutralino"])
        else:
            nEvents_neutralinoMassUnset += 1 # sys.exit("Unable to find neutralino mass in an event!")
        if inputIsGluino:
            if (gluinoMassIsSet and neutralinoMassIsSet):
                h_masses.Fill(generatedMasses["gluino"], generatedMasses["neutralino"])
        elif inputIsSquark:
            if (squarkMassIsSet and neutralinoMassIsSet):
                h_masses.Fill(generatedMasses["squark"], generatedMasses["neutralino"])
        if (inputIsGluino and gluinoMassIsSet):
            if ((minGluinoMassFound < 0) or (generatedMasses["gluino"] < minGluinoMassFound)):
                minGluinoMassFound = generatedMasses["gluino"]
            if ((maxGluinoMassFound < 0) or (generatedMasses["gluino"] > maxGluinoMassFound)):
                maxGluinoMassFound = generatedMasses["gluino"]
        elif (inputIsSquark and squarkMassIsSet):
            if ((minSquarkMassFound < 0) or (generatedMasses["squark"] < minSquarkMassFound)):
                minSquarkMassFound = generatedMasses["squark"]
            if ((maxSquarkMassFound < 0) or (generatedMasses["squark"] > maxSquarkMassFound)):
                maxSquarkMassFound = generatedMasses["squark"]
        if (neutralinoMassIsSet and ((minNeutralinoMassFound < 0) or (generatedMasses["neutralino"] < minNeutralinoMassFound))):
            minNeutralinoMassFound = generatedMasses["neutralino"]
        if (neutralinoMassIsSet and ((maxNeutralinoMassFound < 0) or (generatedMasses["neutralino"] > maxNeutralinoMassFound))):
            maxNeutralinoMassFound = generatedMasses["neutralino"]
        if (inputIsGluino and gluinoMassIsSet):
            if ((generatedMasses["gluino"] < minGluinoMass) or
                (generatedMasses["gluino"] > maxGluinoMass) or
                (generatedMasses["neutralino"] > maxNeutralinoMass)):
                sys.exit("ERROR: Unexpected (gluino, neutralino) masses: {g},{n}".format(g=generatedMasses["gluino"], n=generatedMasses["neutralino"]))
        elif (inputIsSquark and squarkMassIsSet):
            if ((generatedMasses["squark"] < minSquarkMass) or
                (generatedMasses["squark"] > maxSquarkMass) or
                (generatedMasses["neutralino"] > maxNeutralinoMass)):
                sys.exit("ERROR: Unexpected (squark, neutralino) masses: {q},{n}".format(q=generatedMasses["squark"], n=generatedMasses["neutralino"]))

randomNumberGenerator = ROOT.TRandom2()
randomNumberGenerator.SetSeed(100)

listOfInputFiles = []
inputFileNamesFileObject = open(inputArguments.inputFilesList, 'r')
for inputFileName in inputFileNamesFileObject:
    listOfInputFiles.append(inputFileName.strip())
inputFileNamesFileObject.close()

nEvents_gluinoMassUnset = 0
nEvents_squarkMassUnset = 0
nEvents_neutralinoMassUnset = 0
nEvents_analyzed = 0
minGluinoMassFound = -1
maxGluinoMassFound = -1
minSquarkMassFound = -1
maxSquarkMassFound = -1
minNeutralinoMassFound = -1
maxNeutralinoMassFound = -1

for fname in listOfInputFiles:
    fillHistogramsFromFile(fname, h_gluinoMass, h_squarkMass, h_neutralinoMass, h_masses, randomNumberGenerator, inputArguments.MCType)

if (inputArguments.MCType == "gluino"):
    print(("Min gluino mass found: {minGMF}, Max gluino mass found: {maxGMF}, Min neutralino mass found: {minNMF}, Max neutralino mass found: {maxNMF}").format(minGMF=minGluinoMassFound, maxGMF=maxGluinoMassFound, minNMF=minNeutralinoMassFound, maxNMF=maxNeutralinoMassFound))
elif (inputArguments.MCType == "squark"):
    print(("Min squark mass found: {minGMF}, Max squark mass found: {maxGMF}, Min neutralino mass found: {minNMF}, Max neutralino mass found: {maxNMF}").format(minGMF=minSquarkMassFound, maxGMF=maxSquarkMassFound, minNMF=minNeutralinoMassFound, maxNMF=maxNeutralinoMassFound))

outputFile = ROOT.TFile.Open("{oF}/MCGeneratedMasses_{p}_savedObjects.root".format(oF=inputArguments.outputFolder, p=inputArguments.outputPrefix), "RECREATE")
if (inputArguments.MCType == "gluino"):
    c_gluinoMass = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [h_gluinoMass], canvasName = "c_gluinoMass", outputROOTFile=outputFile, outputDocumentName="{oF}/gluinoMassDistribution_{p}".format(oF=inputArguments.outputFolder, p=inputArguments.outputPrefix))
elif (inputArguments.MCType == "squark"):
    c_squarkMass = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [h_squarkMass], canvasName = "c_squarkMass", outputROOTFile=outputFile, outputDocumentName="{oF}/squarkMassDistribution_{p}".format(oF=inputArguments.outputFolder, p=inputArguments.outputPrefix))
c_neutralinoMass = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [h_neutralinoMass], canvasName = "c_neutralinoMass", outputROOTFile=outputFile, outputDocumentName="{oF}/neutralinoMassDistribution_{p}".format(oF=inputArguments.outputFolder, p=inputArguments.outputPrefix))
c_masses = tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [h_masses], canvasName = "c_masses", outputROOTFile=outputFile, outputDocumentName="{oF}/generatedMassesDistribution_{p}".format(oF=inputArguments.outputFolder, p=inputArguments.outputPrefix), customOptStat=0, customPlotOptions_firstObject="TEXTCOLZ")
outputFile.Close()

print("nEvents with neutralino mass unset: {n}, nEvents with gluino mass unset: {ng}, nEvents with squark mass unset: {ns}, total nEvents analyzed = {tot}".format(n=nEvents_neutralinoMassUnset, ng=nEvents_gluinoMassUnset, ns=nEvents_squarkMassUnset, tot=nEvents_analyzed))
