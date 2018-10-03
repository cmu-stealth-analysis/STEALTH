#!/usr/bin/env python

import os
import sys
import numpy as np
import ROOT
import argparse

from tmGeneralUtils import prettyPrintDictionary

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

def find_deltaR(v1, vList):
    min_deltaR = -1.0
    for i in vList:
        deltaR = v1.DeltaR(i)
        if ((min_deltaR < 0) or (deltaR < min_deltaR)):
            min_deltaR = deltaR
    return min_deltaR


def find_m(vList):
    v1 = ROOT.TLorentzVector()
    for i in vList:
        v1 += i
    return v1.M()


def find_weight(chain, nEvtsDict, xSecDict):
    m = 0.0
    for i in range(chain.nMC):
        if chain.mcPID[i] == 1000021:    # MC particle no for gluino
            m = float(chain.mcMass[i])
            break

    n_events = 1
    for i, j in nEvtsDict.items():
        if m < i:
            n_events = j
            break

    return xSecDict[m] / n_events

eventFailureCategories = ["HLTPhoton", "wrongNSelectedPhotons", "incompatiblePhotonSelectionType", "lowInvariantMass", "HLTJet", "wrongNJets", "hTCut", "electronVeto", "muonVeto", "MCGenInformation"]
eventChecksFailDictionary = {eventFailureCategory: 0 for eventFailureCategory in eventFailureCategories}

sw = ROOT.TStopwatch()
sw.Start()

listOfInputFiles = []
if (inputArguments.inputFromFile):
    inputFileNamesFileObject = open(inputArguments.inputFilePath, 'r')
    for inputFileName in inputFileNamesFileObject:
        listOfInputFiles.append(inputFileName.strip())
    inputFileNamesFileObject.close()
else:
    listOfInputFiles.append(inputArguments.inputFilePath)

chain_in = ROOT.TChain('ggNtuplizer/EventTree')
# chain_in.Add('test_*.root')
#chain_in.Add('root://cmseos.fnal.gov://store/user/weinberg/test_*.root')
#chain_in.SetBranchStatus('tau*', 0)
for inputFile in listOfInputFiles:
    # print("Adding... " + inputFile)
    chain_in.Add(inputFile)
n_entries = chain_in.GetEntries()
print('Total entries: ' + str(n_entries))

# file_out = ROOT.TFile('skim_test.root', 'recreate')
file_out = ROOT.TFile(inputArguments.outputFilePath, "RECREATE")
dir_out = file_out.mkdir('ggNtuplizer')
dir_out.cd()
tree_out = chain_in.CloneTree(0)
#tree_out = chain_in.CloneTree()    # Copy entire tree

n_tightPhotons = np.zeros(1, dtype=int)
n_fakes = np.zeros(1, dtype=int)
n_tightJets = np.zeros(1, dtype=int)
st = np.zeros(1, dtype=float)
weight = np.zeros(1, dtype=float)
mDiphoton = np.zeros(1, dtype=float)
deltaR_jets = ROOT.std.vector('float')()

b_n_tightPhotons = tree_out.Branch('n_tightPhotons', n_tightPhotons, 'n_tightPhotons/I')
b_n_fakes = tree_out.Branch('n_fakes', n_fakes, 'n_fakes/I')
b_n_tightJets = tree_out.Branch('n_tightJets', n_tightJets, 'n_tightJets/I')
b_st = tree_out.Branch('st', st, 'st/D')
b_weight = tree_out.Branch('weight', weight, 'weight/D')
b_mDiphoton = tree_out.Branch('mDiphoton', mDiphoton, 'mDiphoton/D')
b_deltaR_jets = tree_out.Branch('deltaR_jets', deltaR_jets)

# For MC with single cross section
# crossSection = 1.0

# For SUSY MC scan
# n_eventsDict = {975: 150000, 9999: 20000}
# crossSectionDict = {}
# with open('./SusyCrossSections13TevGluGlu.txt') as f_crossSection:
#     for l in f_crossSection:
#         line = l.split()
#         crossSectionDict[float(line[0])] = float(line[1])

for j_entry in range(inputArguments.counterStartInclusive, 1 + inputArguments.counterEndInclusive):
# for j_entry in range(100000):
    i_entry = chain_in.LoadTree(j_entry)
    if i_entry < 0:
        break
    nb = chain_in.GetEntry(j_entry)
    if nb <= 0:
        continue
    if j_entry % 100000 == 0:
        print('Processing entry ' + str(j_entry)
              + ' (' + str(round(100.0 * j_entry / n_entries, 2)) + '%)')

    st[0] = 0.0

    # n_genPhotons = 0
    # for i in range(chain_in.nMC):
    #     if (chain_in.mcPID[i] == 22
    #         and chain_in.mcMomPID[i] == 1000023
    #         and chain_in.mcStatusFlag[i] == 7):
    #         n_genPhotons += 1

    n_tightPhotons[0] = 0
    n_fakes[0] = 0
    # v_tightPhotonList = []
    # v_fakeList = []
    # v_loosePhotonList = []
    v_selectedPhotonsList = []
    for i in range(chain_in.nPho):
        if (chain_in.phoEt[i] > 25.0
            and abs(chain_in.phoEta[i]) < 1.442
            and chain_in.phoEleVeto[i] == 1):
            ea_ch = 0.0385
            ea_neu = 0.0636
            ea_pho = 0.1240
            if abs(chain_in.phoEta[i]) > 1.0:
                ea_ch = 0.0468
                ea_neu = 0.1103
                ea_pho = 0.1093
            phoPFChIso_corr = max(chain_in.phoPFChIso[i] - chain_in.rho * ea_ch, 0.0)
            phoPFNeuIso_corr = max(chain_in.phoPFNeuIso[i] - chain_in.rho * ea_neu, 0.0)
            phoPFPhoIso_corr = max(chain_in.phoPFPhoIso[i] - chain_in.rho * ea_pho, 0.0)

            if ((((chain_in.phoIDbit[i]) >> 1) & 1) == 1):   # Medium photon ID
                st[0] += chain_in.phoEt[i]
                n_tightPhotons[0] += 1
                v_photon = ROOT.TLorentzVector()
                v_photon.SetPtEtaPhiE(chain_in.phoEt[i], chain_in.phoEta[i],
                                      chain_in.phoPhi[i], chain_in.phoE[i])
                # v_tightPhotonList.append(v_photon)
                # v_loosePhotonList.append(v_photon)
                v_selectedPhotonsList.append(v_photon)
            elif (chain_in.phoHoverE[i] < 0.035
                  and phoPFNeuIso_corr < 2.491 + 0.0126 * chain_in.phoEt[i] + 0.000026 * chain_in.phoEt[i]**2
                  and phoPFPhoIso_corr < 2.952 + 0.0040 * chain_in.phoEt[i]):
                within_sigmaietaietaRange = ((chain_in.phoSigmaIEtaIEtaFull5x5[i] > 0.0103) and (chain_in.phoSigmaIEtaIEtaFull5x5[i] < 0.02)) # Max shape cut
                within_chargedIsolationRange = ((phoPFChIso_corr > 1.416) and (phoPFChIso_corr < 6.0)) # Max charged iso cut
                if (within_sigmaietaietaRange or within_chargedIsolationRange):
                    st[0] += chain_in.phoEt[i]
                    n_fakes[0] += 1
                    v_fake = ROOT.TLorentzVector()
                    v_fake.SetPtEtaPhiE(chain_in.phoEt[i], chain_in.phoEta[i],
                                        chain_in.phoPhi[i], chain_in.phoE[i])
                    # v_fakeList.append(v_fake)
                    # v_loosePhotonList.append(v_fake)
                    v_selectedPhotonsList.append(v_fake)
    if len(v_selectedPhotonsList) >= 2:
        mDiphoton[0] = find_m([v_selectedPhotonsList[0], v_selectedPhotonsList[1]])

    deltaR_jets.clear()
    n_tightJets[0] = 0
    for i in range(chain_in.nJet):
        if (chain_in.jetPt[i] > 30.0
            and abs(chain_in.jetEta[i]) < 2.5
            and chain_in.jetID[i] == 6    # Tight jet ID
            and chain_in.jetPUID[i] > 0.61):    #  Medium jet PU ID
            v_jet = ROOT.TLorentzVector()
            v_jet.SetPtEtaPhiE(chain_in.jetPt[i], chain_in.jetEta[i],
                               chain_in.jetPhi[i], chain_in.jetEn[i])
            # deltaR = min(find_deltaR(v_jet, v_tightPhotonList),
            #              find_deltaR(v_jet, v_fakeList))
            deltaR = find_deltaR(v_jet, v_selectedPhotonsList)
            if (deltaR > 0.): deltaR_jets.push_back(deltaR)
            if ((deltaR > 0.4) or (deltaR < 0.)):
                st[0] += chain_in.jetPt[i]
                n_tightJets[0] += 1

    st[0] += chain_in.pfMET

    # weight[0] = crossSection / n_entries
    # weight[0] = find_weight(chain_in, n_eventsDict, crossSectionDict)

    # passes_trigger = (chain_in.HLTPho >> 7 & 1 == 1)    # HLT_Photon175_v*
    # passes_trigger = (chain_in.HLTPho >> 16 & 1 == 1
    #                   or chain_in.HLTPho >> 17 & 1 == 1
    #                   or chain_in.HLTPho >> 37 & 1 == 1
    #                   or chain_in.HLTPho >> 38 & 1 == 1)
    passes_trigger = ((chain_in.HLTPho >> 16 & 1 == 1) or (chain_in.HLTPho >> 37 & 1 == 1))
    #passes_trigger = True

    passesPhotonSelection = False
    if (inputArguments.photonSelectionType == "fake"):
        passesPhotonSelection = (n_tightPhotons[0] == 0 and n_fakes[0] == 2)
    elif (inputArguments.photonSelectionType == "mediumfake"):
        passesPhotonSelection = (n_tightPhotons[0] == 1 and n_fakes[0] == 1)
    elif (inputArguments.photonSelectionType == "medium"):
        passesPhotonSelection = (n_tightPhotons[0] == 2 and n_fakes[0] == 0)
    else:
        sys.exit("Unknown selection type: {t}".format(t=inputArguments.photonSelectionType))

    eventPassesSelectionSoFar = True
    if (eventPassesSelectionSoFar and not(passes_trigger)):
        eventPassesSelectionSoFar = False
        eventChecksFailDictionary["HLTPhoton"] += 1
    if (eventPassesSelectionSoFar and not((n_tightPhotons[0] + n_fakes[0] == 2) and v_selectedPhotonsList[0].Et() > 35.0)):
        eventPassesSelectionSoFar = False
        eventChecksFailDictionary["wrongNSelectedPhotons"] += 1
    if (eventPassesSelectionSoFar and not(passesPhotonSelection)):
        eventPassesSelectionSoFar = False
        eventChecksFailDictionary["incompatiblePhotonSelectionType"] += 1
    if ((n_tightPhotons[0] + n_fakes[0] == 2) and (eventPassesSelectionSoFar and not(mDiphoton[0] > 60.0))):
        eventPassesSelectionSoFar = False
        eventChecksFailDictionary["lowInvariantMass"] += 1
    if (eventPassesSelectionSoFar and not(n_tightJets[0] >= 2)):
        eventPassesSelectionSoFar = False
        eventChecksFailDictionary["wrongNJets"] += 1

    if (eventPassesSelectionSoFar): tree_out.Fill()

file_out.Write()
file_out.Close()

sw.Stop()
print('Real time: ' + str(round(sw.RealTime() / 60.0, 2)) + ' minutes')
print('CPU time:  ' + str(round(sw.CpuTime() / 60.0, 2)) + ' minutes')

print("Failure statistics:")
print("_"*200)
prettyPrintDictionary(inputDict=eventChecksFailDictionary, keyPrintOrder=eventFailureCategories)
