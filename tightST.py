import os
import sys
import glob
import array
import numpy as np
import ROOT


sw = ROOT.TStopwatch()
sw.Start()

ggIn = ROOT.TChain("ggNtuplizer/EventTree")
ggIn.Add("ggNTUPLES/SinglePhoton_2016B_sT_Pho100_JetLoose30_Ht700.root")

ggIn.SetBranchStatus("b_*", 0)

#outFile = ROOT.TFile("SinglePhoton_2016B_sT_Pho100_JetLoose30_Ht700.root", "RECREATE")
#outFile = ROOT.TFile("ggNTUPLES/HLTMatched/SinglePhoton_2016B_sT_Pho100Med_Jet30TightHt700_El15Tight_Mu15Tight.root", "RECREATE")
outFile = ROOT.TFile("ggNTUPLES/HLTMatched/SinglePhoton_2016B_sT_Pho100MedEleVeto_Jet50MedHt700_El15Tight_Mu15Tight.root", "RECREATE")
#outFile = ROOT.TFile("DoubleEG_2015D_evtSt.root", "RECREATE")
#outFile = ROOT.TFile("SingleElectron_evtSt.root", "RECREATE")
outDir = outFile.mkdir("ggNtuplizer")
outDir.cd()
ggOut = ggIn.CloneTree(0)


# evtSt = array.array('f', [0])
nJets_  = np.zeros(1, dtype=int)
evtSt_  = np.zeros(1, dtype=float)
StDiff_ = np.zeros(1, dtype=float)
b_nJets = ggOut.Branch("b_nJet", nJets_, "b_nJet/I")
b_evtSt = ggOut.Branch("b_evtST", evtSt_, "b_evtST/D")
b_StDiff = ggOut.Branch("b_STDiff", StDiff_, "b_STDiff/D")

#ggIn.SetBranchStatus("tau*", 0)
nEntries = ggIn.GetEntries()

count = 0
print "nEntries: " + str(nEntries) 
#for jEvt in range(500000):
for jEvt in range(nEntries):

    # Initialize event
    iEvt = ggIn.LoadTree(jEvt)
    if iEvt < 0:
        break
    nb = ggIn.GetEntry(jEvt)
    if nb <= 0:
        continue
    if jEvt % 10000 == 0:
        print "Processing entry " + str(jEvt)

    evtSt = 0.

    # Photon selection
    nPhotons = 0
    for i in range(ggIn.nPho):
        if (ggIn.phoEt[i] > 100.0 and 
            abs(ggIn.phoEta[i]) < 1.479 and #isEB
            #ggIn.phohasPixelSeed[i] == False and
            ggIn.phoEleVeto[i] == True and
            ggIn.phoIDbit[i]>>1&1 == 1): #medium photonID
            nPhotons += 1
            evtSt += ggIn.phoEt[i]
    if nPhotons != 1:
       continue 
    if (ggIn.HLTPho>>27)&1 == 0 and (ggIn.HLTPho>>28)&1 == 0:
      continue

    # Jet selection
    nJets = 0
    sumHt = 0
    for i in range(ggIn.nJet):
        if (ggIn.jetPt[i] > 50.0 and
            #ggIn.jetNHF[i] < 0.90 and
            #ggIn.jetNEF[i] < 0.90 and
            ggIn.jetNHF[i] < 0.99 and
            ggIn.jetNEF[i] < 0.99 and
            ggIn.jetCHF[i] > 0. and
            ggIn.jetCEF[i] < 0.99 and
            ggIn.jetNCH[i] > 0. and
            abs(ggIn.jetEta[i]) < 2.4 and
			#ggIn.jetPFLooseId[i] == 1):
			ggIn.jetPFLooseId[i] == 1 and
			ggIn.jetPUidFullDiscriminant[i] > 0.62):
            nJets += 1
            evtSt += ggIn.jetPt[i]
            sumHt += ggIn.jetPt[i]
    if nJets < 2 or sumHt < 700:
        continue

    # Electron veto
    nEle = 0
    for i in range(ggIn.nEle):
        if (ggIn.elePt[i] > 15.0 and 
            abs(ggIn.eleEta[i]) < 2.5 and
            ggIn.eleIDbit[i]>>3&1 == 1 and # tight electronID
            #ggIn.eleIDbit[i]>>0&1 == 1 and # veto electronID
            abs(ggIn.eleDz[i]) < 0.1 and
            ggIn.elePFPUIso[i] < 0.1):
            nEle += 1
    if nEle != 0:
       continue 

    # Muon veto
    nMu = 0
    for i in range(ggIn.nMu):
        if (ggIn.muPt[i] > 15.0 and 
            ggIn.muIsTightID[i] == 1 and
            #ggIn.muIsLooseID[i] == 1 and
            ggIn.muPFPUIso[i] < 0.12):
            nMu += 1
    if nMu != 0:
        continue

    # MET selection
    if ggIn.pfMET > 15.:
        evtSt += ggIn.pfMET

    # Fill output tree
    evtSt_[0] = evtSt
    nJets_[0] = nJets
    StDiff_[0] = evtSt - ggIn.pfMETsumEt
    ggOut.Fill()
    count += 1
    #if count % 10 == 0:
    #    print " >> S_T = " + str(evtSt_[0])


outFile.Write()
outFile.Close()

sw.Stop()
print "Total events accepted:",count
print "Real time: " + str(sw.RealTime() / 60.0) + " minutes"
print "CPU time:  " + str(sw.CpuTime() / 60.0) + " minutes"
