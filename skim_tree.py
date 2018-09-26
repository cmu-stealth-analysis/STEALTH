import os
import sys
import numpy as np
import ROOT

isMC = True

# region = "onemediumonefake"
# region = "doublefake"
region = "doublemedium"

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

nFailingEvents_HLTTrigger = 0
nFailingEvents_incompatibleSelectionType = 0
nFailingEvents_wrongNSelectedPhotons = 0
nFailingEvents_lowInvariantMass = 0
nFailingEvents_wrongNJets = 0
nFailingEvents_hTCut = 0
nFailingEvents_electronVeto = 0
nFailingEvents_muonVeto = 0
nFailingEvents_MCGenInfo = 0
nSelectedEvents = 0
nTotalEvents = 0

nFailingJets_eta = 0
nFailingJets_pT = 0
nFailingJets_PFLooseID = 0
nFailingJets_puID = 0
nFailingJets_jetID = 0
nFailingJets_deltaR = 0
nSelectedJets = 0
nFailingJets = 0
nTotalJets = 0

nPhotonsFailingpT = 0
nPhotonsFailingEta = 0
nPhotonsFailingElectronVeto = 0
nPhotonsFailingMediumID = 0
nPhotonsFailingHOverE = 0
nPhotonsFailingSigIEIEORChIso = 0
nPhotonsFailingNeuIso = 0
nPhotonsFailingPhoIso = 0
nSelectedPhotons = 0
nFailingPhotons = 0
nTotalPhotons = 0

sw = ROOT.TStopwatch()
sw.Start()

chain_in = ROOT.TChain('ggNtuplizer/EventTree')
# chain_in.Add('test_*.root')
if (isMC):
    chain_in.Add('root://cmseos.fnal.gov//store/user/lpcsusystealth/stealth2018Ntuples/job_SMS-T7WgStealth/SMS-T7WgStealth_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_job_SMS-T7WgStealth/180525_031712/0000/ggtree_mc_96.root')
    chain_in.Add('root://cmseos.fnal.gov//store/user/lpcsusystealth/stealth2018Ntuples/job_SMS-T7WgStealth/SMS-T7WgStealth_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_job_SMS-T7WgStealth/180525_031712/0000/ggtree_mc_97.root')
    chain_in.Add('root://cmseos.fnal.gov//store/user/lpcsusystealth/stealth2018Ntuples/job_SMS-T7WgStealth/SMS-T7WgStealth_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_job_SMS-T7WgStealth/180525_031712/0000/ggtree_mc_98.root')
    chain_in.Add('root://cmseos.fnal.gov//store/user/lpcsusystealth/stealth2018Ntuples/job_SMS-T7WgStealth/SMS-T7WgStealth_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_job_SMS-T7WgStealth/180525_031712/0000/ggtree_mc_99.root')
else:
    chain_in.Add('root://cmseos.fnal.gov//store/user/lpcsusystealth/stealth2018Ntuples/job_DoubleEG_newHLT_v2_Run2017B/DoubleEG/crab_job_DoubleEG_newHLT_v2_Run2017B/180712_154152/0000/ggtree_data_1.root')
    chain_in.Add('root://cmseos.fnal.gov//store/user/lpcsusystealth/stealth2018Ntuples/job_DoubleEG_newHLT_v2_Run2017B/DoubleEG/crab_job_DoubleEG_newHLT_v2_Run2017B/180712_154152/0000/ggtree_data_10.root')
    chain_in.Add('root://cmseos.fnal.gov//store/user/lpcsusystealth/stealth2018Ntuples/job_DoubleEG_newHLT_v2_Run2017B/DoubleEG/crab_job_DoubleEG_newHLT_v2_Run2017B/180712_154152/0000/ggtree_data_100.root')
    chain_in.Add('root://cmseos.fnal.gov//store/user/lpcsusystealth/stealth2018Ntuples/job_DoubleEG_newHLT_v2_Run2017B/DoubleEG/crab_job_DoubleEG_newHLT_v2_Run2017B/180712_154152/0000/ggtree_data_101.root')
    chain_in.Add('root://cmseos.fnal.gov//store/user/lpcsusystealth/stealth2018Ntuples/job_DoubleEG_newHLT_v2_Run2017B/DoubleEG/crab_job_DoubleEG_newHLT_v2_Run2017B/180712_154152/0000/ggtree_data_102.root')
    chain_in.Add('root://cmseos.fnal.gov//store/user/lpcsusystealth/stealth2018Ntuples/job_DoubleEG_newHLT_v2_Run2017B/DoubleEG/crab_job_DoubleEG_newHLT_v2_Run2017B/180712_154152/0000/ggtree_data_103.root')
    chain_in.Add('root://cmseos.fnal.gov//store/user/lpcsusystealth/stealth2018Ntuples/job_DoubleEG_newHLT_v2_Run2017B/DoubleEG/crab_job_DoubleEG_newHLT_v2_Run2017B/180712_154152/0000/ggtree_data_104.root')
    chain_in.Add('root://cmseos.fnal.gov//store/user/lpcsusystealth/stealth2018Ntuples/job_DoubleEG_newHLT_v2_Run2017B/DoubleEG/crab_job_DoubleEG_newHLT_v2_Run2017B/180712_154152/0000/ggtree_data_105.root')
    chain_in.Add('root://cmseos.fnal.gov//store/user/lpcsusystealth/stealth2018Ntuples/job_DoubleEG_newHLT_v2_Run2017B/DoubleEG/crab_job_DoubleEG_newHLT_v2_Run2017B/180712_154152/0000/ggtree_data_106.root')
    chain_in.Add('root://cmseos.fnal.gov//store/user/lpcsusystealth/stealth2018Ntuples/job_DoubleEG_newHLT_v2_Run2017B/DoubleEG/crab_job_DoubleEG_newHLT_v2_Run2017B/180712_154152/0000/ggtree_data_107.root')
    #chain_in.Add('root://cmseos.fnal.gov://store/user/weinberg/test_*.root')
    #chain_in.SetBranchStatus('tau*', 0)
n_entries = chain_in.GetEntries()
print('Total entries: ' + str(n_entries))

outputFileName = "/uscms_data/d3/tmudholk/temp/skim_fromMarc_{selType}.root".format(selType = region)
if (isMC): outputFileName = "/uscms_data/d3/tmudholk/temp/skim_fromMarc_{selType}_isMC.root".format(selType = region)

file_out = ROOT.TFile(outputFileName, 'recreate')
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
# for j_entry in range(700000, 800001):
    nTotalEvents += 1
    i_entry = chain_in.LoadTree(j_entry)
    if i_entry < 0:
        break
    nb = chain_in.GetEntry(j_entry)
    if nb <= 0:
        continue
    if j_entry % 10000 == 0:
        print('Processing entry ' + str(j_entry)
              + ' (' + str(round(100.0 * j_entry / n_entries, 2)) + '%)')

    st[0] = 0.0

    passesMCSelection = True
    if (isMC):
        n_genPhotons = 0
        for i in range(chain_in.nMC):
            if (chain_in.mcPID[i] == 22
                and chain_in.mcMomPID[i] == 1000022
                and chain_in.mcStatusFlag[i] == 7):
                n_genPhotons += 1
        if not (n_genPhotons == 2):
            nFailingEvents_MCGenInfo += 1
            passesMCSelection = False

    n_tightPhotons[0] = 0
    n_fakes[0] = 0
    v_allSelectedPhotonsList = []
    nSelectedPhotonsPassingSubLeadingPtCut = 0
    nSelectedPhotonsPassingLeadingPtCut = 0
    for i in range(chain_in.nPho):
        nTotalPhotons += 1
        phopassessubleadingpT = (chain_in.phoEt[i] > 25.0)
        if not(phopassessubleadingpT): nPhotonsFailingpT += 1
        phopassesleadingpT = (chain_in.phoEt[i] > 35.0)
        phopasseseta = (abs(chain_in.phoEta[i]) < 1.442)
        if not(phopasseseta): nPhotonsFailingEta += 1
        phopasseselectronveto = (chain_in.phoEleVeto[i] == 1)
        if not(phopasseselectronveto): nPhotonsFailingElectronVeto += 1
        phopassescommoncuts = ((phopassessubleadingpT and phopasseseta) and phopasseselectronveto)
        ea_ch = 0.0385
        ea_neu = 0.0636
        ea_pho = 0.1240
        if ((abs(chain_in.phoEta[i]) > 1.0) and (abs(chain_in.phoEta[i]) <= 1.479)):
            ea_ch = 0.0468
            ea_neu = 0.1103
            ea_pho = 0.1093
        elif (abs(chain_in.phoEta[i]) > 1.479): # Compatibility
            ea_ch = 0.
            ea_neu = 0.
            ea_pho = 0.
        phoPFChIso_corr = max(chain_in.phoPFChIso[i] - chain_in.rho * ea_ch, 0.0)
        phoPFNeuIso_corr = max(chain_in.phoPFNeuIso[i] - chain_in.rho * ea_neu, 0.0)
        phoPFPhoIso_corr = max(chain_in.phoPFPhoIso[i] - chain_in.rho * ea_pho, 0.0)
        phopassesmedID = ((((chain_in.phoIDbit[i])>>1)&1) == 1)   # Medium photon ID
        if (phopassesmedID):
            if (phopassescommoncuts):
                nSelectedPhotonsPassingSubLeadingPtCut += 1
                if (phopassesleadingpT): nSelectedPhotonsPassingLeadingPtCut += 1
                nSelectedPhotons += 1
                st[0] += chain_in.phoEt[i]
                n_tightPhotons[0] += 1
                v_photon = ROOT.TLorentzVector()
                v_photon.SetPtEtaPhiE(chain_in.phoEt[i], chain_in.phoEta[i],
                                      chain_in.phoPhi[i], chain_in.phoE[i])
                # v_tightPhotonList.append(v_photon)
                # v_loosePhotonList.append(v_photon)
                v_allSelectedPhotonsList.append(v_photon)
            else:
                nFailingPhotons += 1
        else:
            nPhotonsFailingMediumID += 1
            phopasseshovere = (chain_in.phoHoverE[i] < 0.035)
            if not(phopasseshovere): nPhotonsFailingHOverE += 1
            phopassessigieie = ((chain_in.phoSigmaIEtaIEtaFull5x5[i] >= 0.0103) and (chain_in.phoSigmaIEtaIEtaFull5x5[i] < 0.02)) # Max shape cut
            phopasseschiso = ((phoPFChIso_corr >= 1.416) and (phoPFChIso_corr < 6.0))     # Max charged iso cut
            phopassessigieieorchiso = (phopassessigieie or phopasseschiso)
            if not(phopassessigieieorchiso): nPhotonsFailingSigIEIEORChIso += 1
            phopassesneutiso = (phoPFNeuIso_corr < (2.491 + 0.0126 * chain_in.phoEt[i] + 0.000026 * chain_in.phoEt[i]**2))
            if not(phopassesneutiso): nPhotonsFailingNeuIso += 1
            phopassesphoiso = (phoPFPhoIso_corr < (2.952 + 0.0040 * chain_in.phoEt[i]))
            if not(phopassesphoiso): nPhotonsFailingPhoIso += 1
            if (phopasseshovere
                and phopassessigieieorchiso
                and phopassesneutiso
                and phopassesphoiso
                and phopassescommoncuts):
                nSelectedPhotonsPassingSubLeadingPtCut += 1
                if (phopassesleadingpT): nSelectedPhotonsPassingLeadingPtCut += 1
                nSelectedPhotons += 1
                st[0] += chain_in.phoEt[i]
                n_fakes[0] += 1
                v_fake = ROOT.TLorentzVector()
                v_fake.SetPtEtaPhiE(chain_in.phoEt[i], chain_in.phoEta[i],
                                    chain_in.phoPhi[i], chain_in.phoE[i])
                # v_fakeList.append(v_fake)
                # v_loosePhotonList.append(v_fake)
                v_allSelectedPhotonsList.append(v_fake)
            else:
                nFailingPhotons += 1
    # if len(v_tightPhotonList) >= 2:
    if ((nSelectedPhotonsPassingSubLeadingPtCut == 2) and (nSelectedPhotonsPassingLeadingPtCut >= 1)):
        mDiphoton[0] = find_m([v_allSelectedPhotonsList[0], v_allSelectedPhotonsList[1]])

    deltaR_jets.clear()
    n_tightJets[0] = 0
    for i in range(chain_in.nJet):
        jetPassesQuality = False
        nTotalJets += 1
        jetpassespT = (chain_in.jetPt[i] > 30.0)
        if not(jetpassespT): nFailingJets_pT += 1
        jetpasseseta = (abs(chain_in.jetEta[i]) < 2.5)
        if not(jetpasseseta): nFailingJets_eta += 1
        jetpassesjetid = (chain_in.jetID[i] == 6)     # Tight jet ID
        if not(jetpassesjetid): nFailingJets_jetID += 1
        jetpassespflooseid = (chain_in.jetPFLooseId[i])  # PF loose ID
        if not(jetpassespflooseid): nFailingJets_PFLooseID += 1
        jetpassespuid = (chain_in.jetPUID[i] > 0.61)    #  Medium jet PU ID
        if not(jetpassespuid): nFailingJets_puID += 1
        if (jetpassespT
            and jetpasseseta
            and jetpassesjetid
            and jetpassespflooseid
            and jetpassespuid):
            jetPassesQuality = True
            nSelectedJets += 1
        else:
            nFailingJets += 1
        v_jet = ROOT.TLorentzVector()
        v_jet.SetPtEtaPhiE(chain_in.jetPt[i], chain_in.jetEta[i],
                           chain_in.jetPhi[i], chain_in.jetEn[i])
        # deltaR = min(find_deltaR(v_jet, v_tightPhotonList),
        #              find_deltaR(v_jet, v_fakeList))
        deltaR = find_deltaR(v_jet, v_allSelectedPhotonsList)
        if (not(deltaR > 0.4)):
            nFailingJets_deltaR += 1
        if (jetPassesQuality): deltaR_jets.push_back(deltaR)
        if (jetPassesQuality and (deltaR > 0.4)):
            st[0] += chain_in.jetPt[i]
            n_tightJets[0] += 1

    st[0] += chain_in.pfMET

    # weight[0] = crossSection / n_entries
    # weight[0] = find_weight(chain_in, n_eventsDict, crossSectionDict)

    passes_trigger = True
    if not(isMC): passes_trigger = (chain_in.HLTPho >> 38 & 1 == 1)
    # passes_trigger = (chain_in.HLTPho >> 16 & 1 == 1
    #                   or chain_in.HLTPho >> 17 & 1 == 1
    #                   or chain_in.HLTPho >> 37 & 1 == 1
    #                   or chain_in.HLTPho >> 38 & 1 == 1)
    #passes_trigger = True
    
    if not(passes_trigger): nFailingEvents_HLTTrigger += 1
    hasCorrectNMediumAndFakePhotons = False
    if (region == "doublemedium"):
        hasCorrectNMediumAndFakePhotons = ((n_tightPhotons[0] == 2) and (n_fakes[0] == 0))
    elif (region == "onemediumonefake"):
        hasCorrectNMediumAndFakePhotons = ((n_tightPhotons[0] == 1) and (n_fakes[0] == 1))
    elif (region == "doublefake"):
        hasCorrectNMediumAndFakePhotons = ((n_tightPhotons[0] == 0) and (n_fakes[0] == 2))

    if (not(hasCorrectNMediumAndFakePhotons)): nFailingEvents_incompatibleSelectionType += 1
    if (not((nSelectedPhotonsPassingSubLeadingPtCut == 2) and (nSelectedPhotonsPassingLeadingPtCut >= 1))): nFailingEvents_wrongNSelectedPhotons += 1
    else:
        if not(mDiphoton[0] > 60.0): nFailingEvents_lowInvariantMass += 1

    if not(n_tightJets[0] >= 2): nFailingEvents_wrongNJets += 1
    if (passesMCSelection
        and hasCorrectNMediumAndFakePhotons
        and (nSelectedPhotonsPassingLeadingPtCut >= 1)       # Has high-ET photon
        and (mDiphoton[0] > 60.0)                    # Diphoton cut fully efficient
        and passes_trigger                         # Passes diphoton trigger
        and n_tightJets[0] >= 2):                  # At least 2 tight jets
        tree_out.Fill()
        nSelectedEvents += 1

file_out.Write()
file_out.Close()

sw.Stop()
print('Real time: ' + str(round(sw.RealTime() / 60.0, 2)) + ' minutes')
print('CPU time:  ' + str(round(sw.CpuTime() / 60.0, 2)) + ' minutes')

print("nFailingEvents_HLTTrigger = {n}".format(n = nFailingEvents_HLTTrigger))
print("nFailingEvents_wrongNSelectedPhotons = {n}".format(n = nFailingEvents_wrongNSelectedPhotons))
print("nFailingEvents_incompatibleSelectionType = {n}".format(n = nFailingEvents_incompatibleSelectionType))
print("nFailingEvents_lowInvariantMass = {n}".format(n = nFailingEvents_lowInvariantMass))
print("nFailingEvents_wrongNJets = {n}".format(n = nFailingEvents_wrongNJets))
print("nFailingEvents_MCGenInfo = {n}".format(n = nFailingEvents_MCGenInfo))
print("nSelectedEvents: {n1}/{ntot}".format(n1 = nSelectedEvents, ntot = nTotalEvents))
print("")

print("jets info: (failed/passed/total) = ({nf}/{np}/{ntot})".format(nf = nFailingJets, np = nSelectedJets, ntot = nTotalJets))
print("nFailingJets_eta = {n}".format(n = nFailingJets_eta))
print("nFailingJets_pT = {n}".format(n = nFailingJets_pT))
print("nFailingJets_PFLooseID = {n}".format(n = nFailingJets_PFLooseID))
print("nFailingJets_puID = {n}".format(n = nFailingJets_puID))
print("nFailingJets_jetID = {n}".format(n = nFailingJets_jetID))
print("nFailingJets_deltaR = {n}".format(n = nFailingJets_deltaR))
print("")

print("photons info: (failed/passed/total) = ({nf}/{np}/{ntot})".format(nf = nFailingPhotons, np = nSelectedPhotons, ntot = nTotalPhotons))
print("nPhotonsFailingpT = {n}".format(n = nPhotonsFailingpT))
print("nPhotonsFailingEta = {n}".format(n = nPhotonsFailingEta))
print("nPhotonsFailingElectronVeto = {n}".format(n = nPhotonsFailingElectronVeto))
print("nPhotonsFailingMediumID = {n}".format(n = nPhotonsFailingMediumID))
print("nPhotonsFailingHOverE = {n}".format(n = nPhotonsFailingHOverE))
print("nPhotonsFailingSigIEIEORChIso = {n}".format(n = nPhotonsFailingSigIEIEORChIso))
print("nPhotonsFailingNeuIso = {n}".format(n = nPhotonsFailingNeuIso))
print("nPhotonsFailingPhoIso = {n}".format(n = nPhotonsFailingPhoIso))
