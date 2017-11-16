#!/usr/bin/env python

from __future__ import print_function, division

import os
import sys
import numpy as np
import argparse
import ROOT

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
inputArguments = inputArgumentsParser.parse_args()

parameters = {
    "pTCutSubLeading": 25.,
    "pTCutLeading": 35.,
    "barrelEndcapTransitionEta": 1.479,
    "towerHOverECut": 0.0396,
    "neutralIsolationCutCoefficients": [2.725, 0.0148, 0.000017],
    "photonIsolationCutCoefficients": [2.571, 0.0047],
    "R9Cut": 1.0,
    "sigmaietaietaRange": [0.01022, 0.015],
    "chargedIsolationRange": [0.441, 15.],
    "nSubLeadingPhotons": 2,
    "nLeadingPhotons": 1,
    "jetpTCut": 30.,
    "jetPUIDThreshold": 0.62,
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
    "region1EAValues": {
        "chargedHadrons": 0.036,
        "neutralHadrons": 0.0597,
        "photons": 0.121
    },
    "region2UpperBoundEA": 1.479,
    "region2EAValues": {
        "chargedHadrons": 0.0377,
        "neutralHadrons": 0.0807,
        "photons": 0.1107
    }
}

photonFailureCategories = ["eta", "pT", "hOverE", "neutralIsolation", "photonIsolation", "pixelSeedVeto", "R9", "sigmaietaiataXORchargedIsolation", "mediumIDCut", "conversionSafeElectronVeto"]
globalPhotonChecksFailDictionary = {photonFailureCategory: 0 for photonFailureCategory in photonFailureCategories}
differentialPhotonChecksFailDictionary = {photonFailureCategory: 0 for photonFailureCategory in photonFailureCategories}

jetFailureCategories = ["eta", "pT", "PFLooseID", "puID", "jetID"]
globalJetChecksFailDictionary = {jetFailureCategory: 0 for jetFailureCategory in jetFailureCategories}
differentialJetChecksFailDictionary = {jetFailureCategory: 0 for jetFailureCategory in jetFailureCategories}

eventFailureCategories = ["wrongNPhotons", "HLTJet", "wrongNJets", "hTCut", "electronVeto", "muonVeto"]
globalEventChecksFailDictionary = {eventFailureCategory: 0 for eventFailureCategory in eventFailureCategories}
differentialEventChecksFailDictionary = {eventFailureCategory: 0 for eventFailureCategory in eventFailureCategories}

counterNames = ["failingPhotons", "failingJets", "failingEvents", "acceptedEvents"]
counters = {counterName: 0 for counterName in counterNames}

evtST = 0
nJetsDR = 0

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
    if (absEta >= parameters["barrelEndcapTransitionEta"]):
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
    if (absEta >= parameters["barrelEndcapTransitionEta"]):
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

    hasPixelSeed = (inputTreeObject.phohasPixelSeed[photonIndex] > 0)
    if (hasPixelSeed):
        globalPhotonChecksFailDictionary["pixelSeedVeto"] += 1
        if passesSelection:
            differentialPhotonChecksFailDictionary["pixelSeedVeto"] += 1
            passesSelection = False

    R9 = inputTreeObject.phoR9[photonIndex]
    if (R9 >= parameters["R9Cut"]):
        globalPhotonChecksFailDictionary["R9"] += 1
        if passesSelection:
            differentialPhotonChecksFailDictionary["R9"] += 1
            passesSelection = False

    sigmaietaieta = inputTreeObject.phoSigmaIEtaIEtaFull5x5[photonIndex]
    sigmaietaietaRange = parameters["sigmaietaietaRange"]
    within_sigmaietaietaRange = (sigmaietaieta > sigmaietaietaRange[0] and sigmaietaieta < sigmaietaietaRange[1])

    chargedIsolation = getRhoCorrectedPFIsolation("chargedHadrons", absEta, inputTreeObject.phoPFChIso[photonIndex], eventRho)
    chargedIsolationRange = parameters["chargedIsolationRange"]
    within_chargedIsolationRange = (chargedIsolation > chargedIsolationRange[0] and chargedIsolation < chargedIsolationRange[1])

    if not(within_sigmaietaietaRange^within_chargedIsolationRange):# ^ operator = xor
        globalPhotonChecksFailDictionary["sigmaietaiataXORchargedIsolation"] += 1
        if passesSelection:
            differentialPhotonChecksFailDictionary["sigmaietaiataXORchargedIsolation"] += 1
            passesSelection = False

    if not(passesSelection): counters["failingPhotons"] += 1
    return passesSelection

if (inputArguments.photonSelectionType == "fake"): passesPhotonSelection = passesFakePhotonSelection
elif (inputArguments.photonSelectionType == "medium"): passesPhotonSelection = passesMediumPhotonSelection
else: sys.exit("Undefined photon selection type: " + inputArguments.photonSelectionType + ". Accepted types: \"fake\" or \"medium\"")

def passesJetSelection(inputTreeObject, jetIndex):

    passesSelection = True

    absEta = abs(inputTreeObject.jetEta[jetIndex])
    if (absEta >= parameters["barrelEndcapTransitionEta"]):
        globalJetChecksFailDictionary["eta"] += 1
        if passesSelection:
            differentialJetChecksFailDictionary["eta"] += 1
            passesSelection = False

    pT = inputTreeObject.jetPt[jetIndex]
    if (pT <= parameters["jetpTCut"]):
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

    # if inputTreeObject.HLTPho>>14&1 == 0: # HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90
    #     # we shouldn't enter here because we're already dealing with skims, so throw an error
    #     sys.exit("Found event failing the HLT Photon bit 14 check -- this should never happen for already skimmed inputs!")

    photonPassingSelectionIndices = [] # for DeltaR check: keep a list of photon indices passing photon selection
    nSubLeadingPhotons = 0
    nLeadingPhotons = 0

    eventRho = inputTreeObject.rho
    # print("here1: passesSelection = " + str(passesSelection))
    # print("number of photons: " + str(inputTreeObject.nPho))
    for photonIndex in range(inputTreeObject.nPho):
        if not(passesPhotonSelection(inputTreeObject, photonIndex, eventRho)): continue
        if inputTreeObject.phoEt[photonIndex] > parameters["pTCutLeading"]: nLeadingPhotons += 1
        nSubLeadingPhotons += 1
        evtST += inputTreeObject.phoEt[photonIndex]
        photonPassingSelectionIndices.append(photonIndex)
    # print("here2: passesSelection = " + str(passesSelection))
    if not(nSubLeadingPhotons == parameters["nSubLeadingPhotons"] and nLeadingPhotons >= parameters["nLeadingPhotons"]):
        globalEventChecksFailDictionary["wrongNPhotons"] += 1
        if passesSelection:
            differentialEventChecksFailDictionary["wrongNPhotons"] += 1
            passesSelection = False

    # print("here3: passesSelection = " + str(passesSelection))
    # if (inputTreeObject.HLTJet>>33&1 == 0 and inputTreeObject.HLTJet>>18&1 == 0): # HLT_PFHT900,PFJet450, resp.
    if (False):
        globalEventChecksFailDictionary["HLTJet"] += 1
        if passesSelection:
            differentialEventChecksFailDictionary["HLTJet"] += 1
            passesSelection = False

    # print("here4: passesSelection = " + str(passesSelection))
    nJets = 0
    evtHT = 0
    for jetIndex in range(inputTreeObject.nJet):
        if not(passesJetSelection(inputTreeObject, jetIndex)): continue
        nJets += 1
        evtHT += inputTreeObject.jetPt[jetIndex] # Add jet pT to HT (even though not sure if it's photon)
        # DeltaR check: ensure this jet is well-separated from any of the good photons
        # To avoid double-counting, only add jet pT to ST if we're sure its not a photon 
        minDeltaRij = 100.
        for photonIndex in photonPassingSelectionIndices: # loop over "good" photon indices
            dR = np.hypot(inputTreeObject.phoEta[photonIndex]-inputTreeObject.jetEta[jetIndex],inputTreeObject.phoPhi[photonIndex]-inputTreeObject.jetPhi[jetIndex]) #DeltaR(pho[photonIndex],jet[jetIndex])
            if dR < minDeltaRij:
                minDeltaRij = dR
        if minDeltaRij < parameters["minDeltaRCut"]:
            continue
        nJetsDR += 1 # nJets passing the DeltaR check
        evtST += inputTreeObject.jetPt[jetIndex]

    # print("here5: passesSelection = " + str(passesSelection))
    if (nJetsDR < parameters["nJetsCut"]):
        globalEventChecksFailDictionary["wrongNJets"] += 1
        if passesSelection:
            differentialEventChecksFailDictionary["wrongNJets"] += 1
            passesSelection = False

    # print("here6: passesSelection = " + str(passesSelection))
    if (evtHT < parameters["HTCut"]):
        globalEventChecksFailDictionary["hTCut"] += 1
        if passesSelection:
            differentialEventChecksFailDictionary["hTCut"] += 1
            passesSelection = False

    # print("here7: passesSelection = " + str(passesSelection))
    if (countTightElectrons(inputTreeObject) != parameters["nElectronsCut"]):
        globalEventChecksFailDictionary["electronVeto"] += 1
        if passesSelection:
            differentialEventChecksFailDictionary["electronVeto"] += 1
            passesSelection = False

    # print("here8: passesSelection = " + str(passesSelection))
    if (countTightMuons(inputTreeObject) != parameters["nMuonsCut"]):
        globalEventChecksFailDictionary["muonVeto"] += 1
        if passesSelection:
            differentialEventChecksFailDictionary["muonVeto"] += 1
            passesSelection = False

    # print("here9: passesSelection = " + str(passesSelection))
    if inputTreeObject.pfMET > parameters["METThreshold"]:
        evtST += inputTreeObject.pfMET

    # print("here10: passesSelection = " + str(passesSelection))
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

    # Load input TTrees into TChain
    inputTreeObject = ROOT.TChain("ggNtuplizer/EventTree")
    for inputFile in listOfInputFiles:
        print("Adding... " + inputFile)
        inputTreeObject.Add(inputFile)

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
    print(" >> Real time: {realTime} minutes".format(realTime=sw.RealTime()/60))
    print(" >> CPU time: {cpuTime} minutes".format(cpuTime=sw.CpuTime()/60))
    print(" >> Accepted Event Statistics: {n}/{nTot} ({percent} %)".format(n=counters["acceptedEvents"], nTot=(eventIndexEnd-eventIndexStart), percent=100*counters["acceptedEvents"]/(eventIndexEnd-eventIndexStart)))

#_____ Call main() ______#
if __name__ == '__main__':
    main()
