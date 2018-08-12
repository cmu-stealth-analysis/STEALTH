#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, ROOT
import numpy as np

from tmProgressBar import tmProgressBar
from tmGeneralUtils import prettyPrintDictionary

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Run STEALTH selection.')
inputArgumentsParser.add_argument('--inputFromFile', action='store_true', help="Interpret inputFilePath as text file that has a list of input of files.")
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--outputFilePath', required=True, help='Path to output file.',type=str)
inputArgumentsParser.add_argument('--counterStartInclusive', required=True, help="Event number from input file from which to start. The event with this index is included in the processing.", type=int)
inputArgumentsParser.add_argument('--counterEndInclusive', required=True, help="Event number from input file at which to end. The event with this index is included", type=int)
inputArgumentsParser.add_argument('--photonSelectionType', default="fake", help='Takes value fake for fake photon selection and medium for selection based on medium ID.',type=str)
inputArgumentsParser.add_argument('--year', default=-1, help='Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger, and the photon ID cuts which are based on year-dependent recommendations. Default year: -1, which is used for MC and means the trigger is disabled and the 2017 photon ID recommendations are implemented.', type=int) # Bit 14 for 2016 data: HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90
inputArgumentsParser.add_argument('--JECUncertainty', default=0, help='Apply a uniform upward or downward jet energy uncertainty correction to jet pt. Default: 0, i.e. do not apply any other correction. +/-1 are allowed as well, shifting all jet pt up or down respectively by the relevant jet energy correction.', type=int)
inputArguments = inputArgumentsParser.parse_args()

parameters = {
    "pTCutSubLeading": 25.,
    "pTCutLeading": 35.,
    "photonEtaCut": 1.442,
    "nSubLeadingPhotons": 2,
    "nLeadingPhotons": 1,
    "jetEtaCut": 2.4,
    "jetpTCut": 30.,
    "jetPUIDThreshold": 0.61,
    "jetSelectionID": 6,
    "minDeltaRCut": 0.4,
    "nJetsCut": 2,
    "HTCut": 60.,
    "electronPtCut": 15.,
    "electronEtaCut": 2.5,
    "electronDzCut": 0.1,
    "electronPFPUIsoCut": 0.1,
    "nElectronsCut": 0,
    "muonPtCut": 15.,
    "muonPFPUIsoCut": 0.12,
    "nMuonsCut": 0,
    "METThreshold": 15.,
    "region1UpperBoundEA": 1.0,
    "region2UpperBoundEA": 1.479,
    "PIDs": {
        "photon": 22,
        "gluino": 1000021,
        "neutralino": 1000022
    },
    "nPhotonsWithNeutralinoMom": 2
}

if (inputArguments.year == 2016):
    parameters["towerHOverECut"] = 0.0396
    parameters["sigmaietaietaRange"] = [0.01022, 0.015]
    parameters["chargedIsolationRange"] = [0.441, 15.]
    parameters["neutralIsolationCutCoefficients"] = [2.725, 0.0148, 0.000017]
    parameters["photonIsolationCutCoefficients"] = [2.571, 0.0047]
    parameters["region1EAValues"] = {
        "chargedHadrons": 0.036,
        "neutralHadrons": 0.0597,
        "photons": 0.121
    }
    parameters["region2EAValues"] = {
        "chargedHadrons": 0.0377,
        "neutralHadrons": 0.0807,
        "photons": 0.1107
    }
    parameters["HLTPhotonBit"] = 14
elif (inputArguments.year == 2017 or inputArguments.year == -1):
    parameters["towerHOverECut"] = 0.035
    parameters["sigmaietaietaRange"] = [0.0103, 0.015]
    parameters["chargedIsolationRange"] = [1.416, 15.]
    parameters["neutralIsolationCutCoefficients"] = [2.491, 0.0126, 0.000026]
    parameters["photonIsolationCutCoefficients"] = [2.952, 0.004]
    parameters["region1EAValues"] = {
        "chargedHadrons": 0.0385,
        "neutralHadrons": 0.0636,
        "photons": 0.124
    }
    parameters["region2EAValues"] = {
        "chargedHadrons": 0.0468,
        "neutralHadrons": 0.1103,
        "photons": 0.1093
    }
    parameters["HLTPhotonBit"] = 36
else:
    sys.exit("Only years 2016, 2017, and -1 (for MC) supported at the moment.")

photonFailureCategories = ["eta", "pT", "hOverE", "neutralIsolation", "photonIsolation", "conversionSafeElectronVeto", "sigmaietaiataXORchargedIsolation", "mediumIDCut"]
globalPhotonChecksFailDictionary = {photonFailureCategory: 0 for photonFailureCategory in photonFailureCategories}
differentialPhotonChecksFailDictionary = {photonFailureCategory: 0 for photonFailureCategory in photonFailureCategories}

jetFailureCategories = ["eta", "pT", "PFLooseID", "puID", "jetID"]
globalJetChecksFailDictionary = {jetFailureCategory: 0 for jetFailureCategory in jetFailureCategories}
differentialJetChecksFailDictionary = {jetFailureCategory: 0 for jetFailureCategory in jetFailureCategories}

eventFailureCategories = ["HLTPhoton", "wrongNMediumOrFakePhotons", "wrongNPhotons", "HLTJet", "wrongNJets", "hTCut", "electronVeto", "muonVeto", "MCGenInformation"]
globalEventChecksFailDictionary = {eventFailureCategory: 0 for eventFailureCategory in eventFailureCategories}
differentialEventChecksFailDictionary = {eventFailureCategory: 0 for eventFailureCategory in eventFailureCategories}

counterNames = ["failingPhotons", "failingJets", "failingEvents", "acceptedEvents"]
counters = {counterName: 0 for counterName in counterNames}

evtST = 0
nJetsDR = 0

if (not(inputArguments.JECUncertainty == 0) and not(abs(inputArguments.JECUncertainty) == 1)): sys.exit("The only values currently supported for JECUncertainty are 0 and +/-1.")

def getEffectiveArea(objectType, absEta):
    if (absEta <= parameters["region1UpperBoundEA"]):
        return parameters["region1EAValues"][objectType]
    elif (absEta <= parameters["region2UpperBoundEA"]):
        return parameters["region2EAValues"][objectType]
    else:
        return 0. # Potentially dangerous

def getRhoCorrectedPFIsolation(objectType, absEta, nominalPFIso, eventRho):
    return max(nominalPFIso - eventRho*getEffectiveArea(objectType, absEta), 0.)

def passesMediumPhotonSelection(inputTreeObject, photonIndex, eventRho):
    passesSelection = True

    absEta = abs(inputTreeObject.phoEta[photonIndex])
    if (absEta >= parameters["photonEtaCut"]):
        globalPhotonChecksFailDictionary["eta"] += 1
        if passesSelection:
            differentialPhotonChecksFailDictionary["eta"] += 1
            passesSelection = False

    pT = inputTreeObject.phoEt[photonIndex]
    if pT <= parameters["pTCutSubLeading"]:
        globalPhotonChecksFailDictionary["pT"] += 1
        if passesSelection:
            differentialPhotonChecksFailDictionary["pT"] += 1
            passesSelection = False

    passesMediumID = ((inputTreeObject.phoIDbit[photonIndex]>>1&1) == 1)
    if not(passesMediumID):
        globalPhotonChecksFailDictionary["mediumIDCut"] += 1
        if passesSelection:
            differentialPhotonChecksFailDictionary["mediumIDCut"] += 1
            passesSelection = False

    passesElectronVeto = (inputTreeObject.phoEleVeto[photonIndex] == 1)
    if not(passesElectronVeto):
        globalPhotonChecksFailDictionary["conversionSafeElectronVeto"] += 1
        if passesSelection:
            differentialPhotonChecksFailDictionary["conversionSafeElectronVeto"] += 1
            passesSelection = False

    if not(passesSelection): counters["failingPhotons"] += 1
    return passesSelection

def passesFakePhotonSelection(inputTreeObject, photonIndex, eventRho):
    passesSelection = True

    absEta = abs(inputTreeObject.phoEta[photonIndex])
    if (absEta >= parameters["photonEtaCut"]):
        globalPhotonChecksFailDictionary["eta"] += 1
        if passesSelection:
            differentialPhotonChecksFailDictionary["eta"] += 1
            passesSelection = False

    pT = inputTreeObject.phoEt[photonIndex]
    if pT <= parameters["pTCutSubLeading"]:
        globalPhotonChecksFailDictionary["pT"] += 1
        if passesSelection:
            differentialPhotonChecksFailDictionary["pT"] += 1
            passesSelection = False

    towerHOverE = inputTreeObject.phoHoverE[photonIndex]
    if towerHOverE >= parameters["towerHOverECut"]:
        globalPhotonChecksFailDictionary["hOverE"] += 1
        if passesSelection:
            differentialPhotonChecksFailDictionary["hOverE"] += 1
            passesSelection = False

    PFNeutralIsolation = getRhoCorrectedPFIsolation("neutralHadrons", absEta, inputTreeObject.phoPFNeuIso[photonIndex], eventRho)
    neutralIsolationCutCoefficients = parameters["neutralIsolationCutCoefficients"]
    if (PFNeutralIsolation >= neutralIsolationCutCoefficients[0] + pT*neutralIsolationCutCoefficients[1] + pT*pT*neutralIsolationCutCoefficients[2]) :
        globalPhotonChecksFailDictionary["neutralIsolation"] += 1
        if passesSelection:
            differentialPhotonChecksFailDictionary["neutralIsolation"] += 1
            passesSelection = False

    PFPhotonIsolation = getRhoCorrectedPFIsolation("photons", absEta, inputTreeObject.phoPFPhoIso[photonIndex], eventRho)
    photonIsolationCutCoefficients = parameters["photonIsolationCutCoefficients"]
    if (PFPhotonIsolation >= photonIsolationCutCoefficients[0] + pT*photonIsolationCutCoefficients[1]):
        globalPhotonChecksFailDictionary["photonIsolation"] += 1
        if passesSelection:
            differentialPhotonChecksFailDictionary["photonIsolation"] += 1
            passesSelection = False

    passesElectronVeto = (inputTreeObject.phoEleVeto[photonIndex] == 1)
    if not(passesElectronVeto):
        globalPhotonChecksFailDictionary["conversionSafeElectronVeto"] += 1
        if passesSelection:
            differentialPhotonChecksFailDictionary["conversionSafeElectronVeto"] += 1
            passesSelection = False

    sigmaietaieta = inputTreeObject.phoSigmaIEtaIEtaFull5x5[photonIndex]
    sigmaietaietaRange = parameters["sigmaietaietaRange"]
    within_sigmaietaietaRange = (sigmaietaieta > sigmaietaietaRange[0] and sigmaietaieta < sigmaietaietaRange[1])

    chargedIsolation = getRhoCorrectedPFIsolation("chargedHadrons", absEta, inputTreeObject.phoPFChIso[photonIndex], eventRho)
    chargedIsolationRange = parameters["chargedIsolationRange"]
    within_chargedIsolationRange = (chargedIsolation > chargedIsolationRange[0] and chargedIsolation < chargedIsolationRange[1])

    # if not(within_sigmaietaietaRange^within_chargedIsolationRange):# ^ operator = xor
    if not(within_sigmaietaietaRange or within_chargedIsolationRange):
        globalPhotonChecksFailDictionary["sigmaietaiataXORchargedIsolation"] += 1
        if passesSelection:
            differentialPhotonChecksFailDictionary["sigmaietaiataXORchargedIsolation"] += 1
            passesSelection = False

    if not(passesSelection): counters["failingPhotons"] += 1
    return passesSelection

def passesExtraMCSelection(inputTreeObject):
    nPhotonsWithNeutralinoMom = 0
    for generatedParticleIndex in range(inputTreeObject.nMC):
        if ((inputTreeObject.mcPID[generatedParticleIndex] == parameters["PIDs"]["photon"]) and (inputTreeObject.mcMomPID[generatedParticleIndex] == parameters["PIDs"]["neutralino"])): nPhotonsWithNeutralinoMom += 1
    return (nPhotonsWithNeutralinoMom == parameters["nPhotonsWithNeutralinoMom"])

def passesJetSelection(inputTreeObject, jetIndex):
    passesSelection = True

    absEta = abs(inputTreeObject.jetEta[jetIndex])
    if (absEta >= parameters["jetEtaCut"]):
        globalJetChecksFailDictionary["eta"] += 1
        if passesSelection:
            differentialJetChecksFailDictionary["eta"] += 1
            passesSelection = False

    pT = inputTreeObject.jetPt[jetIndex]
    if ((pT + inputArguments.JECUncertainty*inputTreeObject.jetJECUnc[jetIndex]*inputTreeObject.jetPt[jetIndex]) <= parameters["jetpTCut"]):
        globalJetChecksFailDictionary["pT"] += 1
        if passesSelection:
            differentialJetChecksFailDictionary["pT"] += 1
            passesSelection = False

    PFLooseID = inputTreeObject.jetPFLooseId[jetIndex]
    if not(PFLooseID):
        globalJetChecksFailDictionary["PFLooseID"] += 1
        if passesSelection:
            differentialJetChecksFailDictionary["PFLooseID"] += 1
            passesSelection = False

    puID = inputTreeObject.jetPUID[jetIndex]
    if (puID <= parameters["jetPUIDThreshold"]):
        globalJetChecksFailDictionary["puID"] += 1
        if passesSelection:
            differentialJetChecksFailDictionary["puID"] += 1
            passesSelection = False

    jetID = inputTreeObject.jetID[jetIndex]
    if not(jetID == 6):
        globalJetChecksFailDictionary["jetID"] += 1
        if passesSelection:
            differentialJetChecksFailDictionary["jetID"] += 1
            passesSelection = False

    if not(passesSelection): counters["failingJets"] += 1
    return passesSelection

def countTightElectrons(inputTreeObject):
    nElectrons = 0
    for electronIndex in range(inputTreeObject.nEle):
        if not(inputTreeObject.elePt[electronIndex] > parameters["electronPtCut"]
                and inputTreeObject.eleIDbit[electronIndex]>>3&1 == 1 # >>0:veto, >>1:loose, >>2:medium, >>3:tight
                and abs(inputTreeObject.eleEta[electronIndex]) < parameters["electronEtaCut"]
                and abs(inputTreeObject.eleDz[electronIndex]) < parameters["electronDzCut"]
                and inputTreeObject.elePFPUIso[electronIndex] < parameters["electronPFPUIsoCut"]):
            continue
        nElectrons += 1
    return nElectrons

def countTightMuons(inputTreeObject):
    nMuons = 0
    for muonIndex in range(inputTreeObject.nMu):
        if not (inputTreeObject.muPt[muonIndex] > parameters["muonPtCut"]
                and inputTreeObject.muPFPUIso[muonIndex] < parameters["muonPFPUIsoCut"]
                and inputTreeObject.muIDbit[muonIndex]>>2&1 == 1 # >>0:loose, >>1:med, >>2:tight, >>3:soft, >>4:highpT
        ):
            continue
        nMuons += 1
    return nMuons

def eventPassesSelection(inputTreeObject):
    global evtST, nJetsDR
    passesSelection = True

    isMCSelection = (inputArguments.photonSelectionType == "mediumMC" or inputArguments.photonSelectionType == "fakeMC" or inputArguments.photonSelectionType == "mediumfakeMC")

    if (not(isMCSelection) and not((parameters["HLTPhotonBit"]) < 0)): # Apply HLT photon cut for data, not for MC, and only if the HLT bit is passed explicitly overwriting its default value of -1
        if (((inputTreeObject.HLTPho>>(parameters["HLTPhotonBit"]))&1) == 0):
            globalEventChecksFailDictionary["HLTPhoton"] += 1
            if passesSelection:
                differentialEventChecksFailDictionary["HLTPhoton"] += 1
                passesSelection = False

    photonPassingSelectionIndices = [] # for DeltaR check: keep a list of photon indices passing photon selection

    eventRho = inputTreeObject.rho

    nSubLeadingPhotons = 0
    nLeadingPhotons = 0
    required_nMediumPhotons = 0
    required_nFakePhotons = 0
    if (inputArguments.photonSelectionType == "medium" or inputArguments.photonSelectionType == "mediumMC"):
        required_nMediumPhotons = 2
        required_nFakePhotons = 0
    elif (inputArguments.photonSelectionType == "fake" or inputArguments.photonSelectionType == "fakeMC"):
        required_nMediumPhotons = 0
        required_nFakePhotons = 2
    elif (inputArguments.photonSelectionType == "mediumfake" or inputArguments.photonSelectionType == "mediumfakeMC"):
        required_nMediumPhotons = 1
        required_nFakePhotons = 1
    else:
        sys.exit("Unrecognized selection type: {type}. Supported: \"medium\", \"fake\", \"mediumfake\", \"mediumMC\", \"fakeMC\", \"mediumfakeMC\".")

    nMediumPhotons = 0
    nFakePhotons = 0
    for photonIndex in range(inputTreeObject.nPho):
        isMedium = False
        isFake = False
        if (passesMediumPhotonSelection(inputTreeObject, photonIndex, eventRho)):
            isMedium = True
            nMediumPhotons += 1
        elif (passesFakePhotonSelection(inputTreeObject, photonIndex, eventRho)):
            isFake = True
            nFakePhotons += 1
        if not(isMedium or isFake): continue
        if inputTreeObject.phoEt[photonIndex] > parameters["pTCutLeading"]: nLeadingPhotons += 1
        nSubLeadingPhotons += 1
        evtST += inputTreeObject.phoEt[photonIndex]
        photonPassingSelectionIndices.append(photonIndex)

    if not(nMediumPhotons == required_nMediumPhotons and nFakePhotons == required_nFakePhotons):
        globalEventChecksFailDictionary["wrongNMediumOrFakePhotons"] += 1
        if passesSelection:
            differentialEventChecksFailDictionary["wrongNMediumOrFakePhotons"] += 1
            passesSelection = False

    if not(nSubLeadingPhotons == parameters["nSubLeadingPhotons"] and nLeadingPhotons >= parameters["nLeadingPhotons"]):
        globalEventChecksFailDictionary["wrongNPhotons"] += 1
        if passesSelection:
            differentialEventChecksFailDictionary["wrongNPhotons"] += 1
            passesSelection = False

    if (isMCSelection):
        if not(passesExtraMCSelection(inputTreeObject)):
            globalEventChecksFailDictionary["MCGenInformation"] += 1
            if passesSelection:
                differentialEventChecksFailDictionary["MCGenInformation"] += 1
                passesSelection = False

    # if (inputTreeObject.HLTJet>>33&1 == 0 and inputTreeObject.HLTJet>>18&1 == 0): # HLT_PFHT900,PFJet450, resp.
    if (False):
        globalEventChecksFailDictionary["HLTJet"] += 1
        if passesSelection:
            differentialEventChecksFailDictionary["HLTJet"] += 1
            passesSelection = False

    nJets = 0
    evtHT = 0
    for jetIndex in range(inputTreeObject.nJet):
        if not(passesJetSelection(inputTreeObject, jetIndex)): continue
        nJets += 1
        evtHT += (inputTreeObject.jetPt[jetIndex] + inputArguments.JECUncertainty*inputTreeObject.jetJECUnc[jetIndex]*inputTreeObject.jetPt[jetIndex]) # Add jet pT to HT (even though not sure if it's photon)
        # DeltaR check: ensure this jet is well-separated from any of the good photons
        # To avoid double-counting, only add jet pT to ST if we're sure its not a photon 
        minDeltaRij = 100.
        for photonIndex in photonPassingSelectionIndices: # loop over "good" photon indices
            deltaR = np.hypot(inputTreeObject.phoEta[photonIndex]-inputTreeObject.jetEta[jetIndex],inputTreeObject.phoPhi[photonIndex]-inputTreeObject.jetPhi[jetIndex]) #DeltaR(pho[photonIndex],jet[jetIndex])
            if deltaR < minDeltaRij:
                minDeltaRij = deltaR
        if minDeltaRij < parameters["minDeltaRCut"]:
            continue
        nJetsDR += 1 # nJets passing the DeltaR check
        evtST += (inputTreeObject.jetPt[jetIndex] + inputArguments.JECUncertainty*inputTreeObject.jetJECUnc[jetIndex]*inputTreeObject.jetPt[jetIndex])

    if (nJetsDR < parameters["nJetsCut"]):
        globalEventChecksFailDictionary["wrongNJets"] += 1
        if passesSelection:
            differentialEventChecksFailDictionary["wrongNJets"] += 1
            passesSelection = False

    if (evtHT < parameters["HTCut"]):
        globalEventChecksFailDictionary["hTCut"] += 1
        if passesSelection:
            differentialEventChecksFailDictionary["hTCut"] += 1
            passesSelection = False

    if (countTightElectrons(inputTreeObject) != parameters["nElectronsCut"]):
        globalEventChecksFailDictionary["electronVeto"] += 1
        if passesSelection:
            differentialEventChecksFailDictionary["electronVeto"] += 1
            passesSelection = False

    if (countTightMuons(inputTreeObject) != parameters["nMuonsCut"]):
        globalEventChecksFailDictionary["muonVeto"] += 1
        if passesSelection:
            differentialEventChecksFailDictionary["muonVeto"] += 1
            passesSelection = False

    if inputTreeObject.pfMET > parameters["METThreshold"]:
        evtST += inputTreeObject.pfMET

    if not(passesSelection): counters["failingEvents"] += 1
    return passesSelection

## MAIN ##
def main():

    global evtST, nJetsDR

    # Keep time
    sw = ROOT.TStopwatch()
    sw.Start()

    # For DeltaR cut
    # minDeltaRij = 100.
    nJetsTot = 0

    listOfInputFiles = []
    if (inputArguments.inputFromFile):
        inputFileNamesFileObject = open(inputArguments.inputFilePath, 'r')
        for inputFileName in inputFileNamesFileObject:
            listOfInputFiles.append(inputFileName.strip())
        inputFileNamesFileObject.close()
    else:
        listOfInputFiles.append(inputArguments.inputFilePath)

    print("Started adding input files into chain from this list of input files: {lIF}".format(lIF=str(listOfInputFiles)))
    # Load input TTrees into TChain
    inputTreeObject = ROOT.TChain("ggNtuplizer/EventTree")
    for inputFile in listOfInputFiles:
        # print("Adding... " + inputFile)
        inputTreeObject.Add(inputFile)
    print("Finished adding input files to chain.")

    # Initialize output file as empty clone
    outFile = ROOT.TFile(inputArguments.outputFilePath, "RECREATE")
    outDir = outFile.mkdir("ggNtuplizer")
    outDir.cd()
    outputTreeObject = inputTreeObject.CloneTree(0)
    print(" >> Output file: " + inputArguments.outputFilePath)

    # Initialize output branches 
    nJets_  = np.zeros(1, dtype=int)
    evtST_  = np.zeros(1, dtype=float)
    b_nJets = outputTreeObject.Branch("b_nJets", nJets_, "b_nJets/I")
    b_evtST = outputTreeObject.Branch("b_evtST", evtST_, "b_evtST/D")

    ##### EVENT SELECTION START #####

    nEvts = inputTreeObject.GetEntries()
    print(" >> nEvts: " + str(nEvts))

    # Event range to process
    eventIndexStart = inputArguments.counterStartInclusive
    eventIndexEnd   = 1 + inputArguments.counterEndInclusive
    # if inputArguments.maxNEvts > 0: eventIndexEnd = inputArguments.maxNEvts
    print(" >> Processing entries: [{start} -> {end})".format(start=eventIndexStart, end=eventIndexEnd))

    progressBar = tmProgressBar(eventIndexEnd - eventIndexStart)
    progressBarUpdatePeriod = 1+((eventIndexEnd-eventIndexStart)//1000)
    progressBar.initializeTimer()
    for eventIndex in range(eventIndexStart,eventIndexEnd):
        # Initialize event
        if eventIndex > nEvts:
            sys.exit("Event counter falls outside event range. Counter: " + str(eventIndex) + ", available number of events: " + str(nEvts))
            # break
        treeStatus = inputTreeObject.LoadTree(eventIndex)
        if treeStatus < 0:
            # sys.exit("Tree unreadable.")
            break
        evtStatus = inputTreeObject.GetEntry(eventIndex)
        if evtStatus <= 0:
            # sys.exit("Event in tree unreadable.")
            continue

        if ((eventIndex-eventIndexStart)%progressBarUpdatePeriod == 0 or eventIndex == (eventIndexEnd-1)): progressBar.updateBar((eventIndex-eventIndexStart)/(eventIndexEnd-eventIndexStart), eventIndex-eventIndexStart)

        evtST = 0
        nJetsDR = 0
        passesSelection = eventPassesSelection(inputTreeObject)
        if not(passesSelection): continue
        # Write this evt to output tree
        evtST_[0] = evtST
        nJets_[0] = nJetsDR
        outputTreeObject.Fill()
        counters["acceptedEvents"] += 1

    ##### EVENT SELECTION END #####
    progressBar.terminate()
    outFile.Write()
    outFile.Close()

    sw.Stop()
    print("_"*200)
    print("_"*200)
    print("Parameters:")
    prettyPrintDictionary(parameters)
    print("Arguments passed: {args}".format(args=inputArguments))
    print("_"*200)
    print("-"*200)
    print(" Counters: ")
    prettyPrintDictionary(counters)
    globalPhotonPercentagesDictionary = {key: 100*globalPhotonChecksFailDictionary[key]/counters["failingPhotons"] for key in globalPhotonChecksFailDictionary.keys()}
    globalJetPercentagesDictionary = {key: 100*globalJetChecksFailDictionary[key]/counters["failingJets"] for key in globalJetChecksFailDictionary.keys()}
    globalEventPercentagesDictionary = {key: 100*globalEventChecksFailDictionary[key]/counters["failingEvents"] for key in globalEventChecksFailDictionary}
    differentialPhotonPercentagesDictionary = {key: 100*differentialPhotonChecksFailDictionary[key]/counters["failingPhotons"] for key in differentialPhotonChecksFailDictionary.keys()}
    differentialJetPercentagesDictionary = {key: 100*differentialJetChecksFailDictionary[key]/counters["failingJets"] for key in differentialJetChecksFailDictionary.keys()}
    differentialEventPercentagesDictionary = {key: 100*differentialEventChecksFailDictionary[key]/counters["failingEvents"] for key in differentialEventChecksFailDictionary}
    print("_"*200)
    print("-"*200)
    print(" >> Global photon failure percentages:")
    prettyPrintDictionary(globalPhotonPercentagesDictionary, "10.7f", photonFailureCategories)
    print(" >> Global jet failure percentages:")
    prettyPrintDictionary(globalJetPercentagesDictionary, "10.7f", jetFailureCategories)
    print(" >> Global event failure percentages:")
    prettyPrintDictionary(globalEventPercentagesDictionary, "10.7f", eventFailureCategories)
    print("-"*200)
    print("_"*200)
    print("")
    print(" >> Differential photon failure percentages:")
    prettyPrintDictionary(differentialPhotonPercentagesDictionary, "10.7f", photonFailureCategories)
    print(" >> Differential jet failure percentages:")
    prettyPrintDictionary(differentialJetPercentagesDictionary, "10.7f", jetFailureCategories)
    print(" >> Differential event failure percentages:")
    prettyPrintDictionary(differentialEventPercentagesDictionary, "10.7f", eventFailureCategories)
    print("-"*200)
    print("_"*200)
    print("")
    print("-"*200)
    print("_"*200)
    print("")
    print(" >> Differential photon failure counters:")
    prettyPrintDictionary(inputDict=differentialPhotonChecksFailDictionary, keyPrintOrder=photonFailureCategories)
    print(" >> Differential jet failure counters:")
    prettyPrintDictionary(inputDict=differentialJetChecksFailDictionary, keyPrintOrder=jetFailureCategories)
    print(" >> Differential event failure counters:")
    prettyPrintDictionary(inputDict=differentialEventChecksFailDictionary, keyPrintOrder=eventFailureCategories)
    print("-"*200)
    print("_"*200)
    print("")
    print(" >> Real time: {realTime} minutes".format(realTime=sw.RealTime()/60))
    print(" >> CPU time: {cpuTime} minutes".format(cpuTime=sw.CpuTime()/60))
    print(" >> Accepted Event Statistics: {n}/{nTot} ({percent} %)".format(n=counters["acceptedEvents"], nTot=(eventIndexEnd-eventIndexStart), percent=100*counters["acceptedEvents"]/(eventIndexEnd-eventIndexStart)))

#_____ Call main() ______#
if __name__ == '__main__':
    main()
