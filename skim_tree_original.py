import os
import sys
import numpy as np
import ROOT


def find_deltaR(v1, vList):
    min_deltaR = 999.9
    for i in vList:
        deltaR = v1.DeltaR(i)
        if deltaR < min_deltaR:
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


sw = ROOT.TStopwatch()
sw.Start()

chain_in = ROOT.TChain('ggNtuplizer/EventTree')
chain_in.Add('test_*.root')
#chain_in.Add('root://cmseos.fnal.gov://store/user/weinberg/test_*.root')
#chain_in.SetBranchStatus('tau*', 0)
n_entries = chain_in.GetEntries()
print('Total entries: ' + str(n_entries))

file_out = ROOT.TFile('skim_test.root', 'recreate')
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

for j_entry in range(n_entries):
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
    v_tightPhotonList = []
    v_fakeList = []
    v_loosePhotonList = []
    for i in range(chain_in.nPho):
        if (chain_in.phoEt[i] > 25.0
            and abs(chain_in.phoEta[i]) < 1.444
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

            if chain_in.phoIDbit[i] >> 1 & 1 == 1:   # Medium photon ID
                st[0] += chain_in.phoEt[i]
                n_tightPhotons[0] += 1
                v_photon = ROOT.TLorentzVector()
                v_photon.SetPtEtaPhiE(chain_in.phoEt[i], chain_in.phoEta[i],
                                      chain_in.phoPhi[i], chain_in.phoE[i])
                v_tightPhotonList.append(v_photon)
                v_loosePhotonList.append(v_photon)
            elif (chain_in.phoHoverE[i] < 0.035
                  and chain_in.phoSigmaIEtaIEtaFull5x5[i] < 0.02    # Max shape cut
                  and phoPFChIso_corr < 6.0    # Max charged iso cut
                  and phoPFNeuIso_corr < 2.491 + 0.0126 * chain_in.phoEt[i] + 0.000026 * chain_in.phoEt[i]**2
                  and phoPFPhoIso_corr < 2.952 + 0.0040 * chain_in.phoEt[i]):
                st[0] += chain_in.phoEt[i]
                n_fakes[0] += 1
                v_fake = ROOT.TLorentzVector()
                v_fake.SetPtEtaPhiE(chain_in.phoEt[i], chain_in.phoEta[i],
                                    chain_in.phoPhi[i], chain_in.phoE[i])
                v_fakeList.append(v_fake)
                v_loosePhotonList.append(v_fake)
    if len(v_tightPhotonList) >= 2:
        mDiphoton[0] = find_m([v_tightPhotonList[0], v_tightPhotonList[1]])

    deltaR_jets.clear()
    n_tightJets[0] = 0
    for i in range(chain_in.nJet):
        if (chain_in.jetPt[i] > 30.0
            and abs(chain_in.jetEta[i]) < 2.5
            and chain_in.jetID[i] >> 1 & 1 == 1    # Tight jet ID
            and chain_in.jetPUID[i] > 0.61):    #  Medium jet PU ID
            v_jet = ROOT.TLorentzVector()
            v_jet.SetPtEtaPhiE(chain_in.jetPt[i], chain_in.jetEta[i],
                               chain_in.jetPhi[i], chain_in.jetEn[i])
            deltaR = min(find_deltaR(v_jet, v_tightPhotonList),
                         find_deltaR(v_jet, v_fakeList))
            deltaR_jets.push_back(deltaR)
            if deltaR > 0.4:
                st[0] += chain_in.jetPt[i]
                n_tightJets[0] += 1

    st[0] += chain_in.pfMET

    # weight[0] = crossSection / n_entries
    # weight[0] = find_weight(chain_in, n_eventsDict, crossSectionDict)

    # passes_trigger = (chain_in.HLTPho >> 7 & 1 == 1)    # HLT_Photon175_v*
    passes_trigger = (chain_in.HLTPho >> 16 & 1 == 1
                      or chain_in.HLTPho >> 17 & 1 == 1
                      or chain_in.HLTPho >> 37 & 1 == 1
                      or chain_in.HLTPho >> 38 & 1 == 1)
    #passes_trigger = True

    if (n_tightPhotons <= 1    # Not in signal region
        and v_loosePhotonList[0].Et() > 35.0    # Has high-ET photon
        and mDiphoton[0] > 60.0    # Diphoton cut fully efficient
        and passes_trigger):    # Passes diphoton trigger
        tree_out.Fill()

file_out.Write()
file_out.Close()

sw.Stop()
print('Real time: ' + str(round(sw.RealTime() / 60.0, 2)) + ' minutes')
print('CPU time:  ' + str(round(sw.CpuTime() / 60.0, 2)) + ' minutes')
