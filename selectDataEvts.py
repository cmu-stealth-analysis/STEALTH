import os
import sys
import numpy as np
import ROOT

# Keep time
print " >> Running STEALTH 2016 Data Selection..."
sw = ROOT.TStopwatch()
sw.Start()

# Load input TTrees into TChain
#ggInStr = "/eos/uscms/store/user/mba2012/JetHT/JetHT_2016G_Pho30Loose_*.root"
ggInStr = "~/eos/cms/store/user/mandrews/JetHT/JetHT_2016*_Pho30Loose*.root"
ggIn = ROOT.TChain("ggNtuplizer/EventTree")
ggIn.Add(ggInStr)
nEvts = ggIn.GetEntries()
print " >> Input file(s):",ggInStr
print " >> nEvts:",nEvts

# Initialize output file as empty clone
outFileStr = "stNTUPLES/DATA/JetHTCtoH_selA_2to3jt.root"
outFile = ROOT.TFile(outFileStr, "RECREATE")
outDir = outFile.mkdir("ggNtuplizer")
outDir.cd()
ggOut = ggIn.CloneTree(0) 
print " >> Output file:",outFileStr

# Initialize output branches 
nJets_  = np.zeros(1, dtype=int)
evtST_  = np.zeros(1, dtype=float)
b_nJets = ggOut.Branch("b_nJets", nJets_, "b_nJets/I")
b_evtST = ggOut.Branch("b_evtST", evtST_, "b_evtST/D")

##### EVENT SELECTION START #####
nAcc = 0
iEvtStart = 0
iEvtEnd   = nEvts
#iEvtEnd   = 20000
print " >> Processing entries: [",iEvtStart,"->",iEvtEnd,")"
for jEvt in range(iEvtStart,iEvtEnd):

	# Initialize event
	if jEvt > nEvts:
		break
	treeStatus = ggIn.LoadTree(jEvt)
	if treeStatus < 0:
		break
	evtStatus = ggIn.GetEntry(jEvt)
	if evtStatus <= 0:
		continue
	if jEvt % 10000 == 0:
		print " .. Processing entry",jEvt

	evtST = 0.

	# Photon selection
	nPhotons = 0
	for i in range(ggIn.nPho):
		if (ggIn.phoEt[0] < 35.):
			continue
		if (ggIn.phoEt[i] > 30.0
				and ggIn.phoIDbit[i]>>1&1 != 0 # >>0:loose, >>1:medium, >>2:tight
				and ggIn.phoEleVeto[i] == True
				and abs(ggIn.phoEta[i]) < 1.479 # isEB
				):
			nPhotons += 1
			evtST += ggIn.phoEt[i]
	#if nPhotons != 1:
	if nPhotons != 2:
		continue 
	#if (ggIn.HLTPho>>27)&1 == 0 and (ggIn.HLTPho>>28)&1 == 0:
	#  continue

	# Jet selection
	nJets = 0
	evtHT = 0
	for i in range(ggIn.nJet):
		if (ggIn.jetPt[i] > 30.0
				and abs(ggIn.jetEta[i]) < 2.4
				and ggIn.jetPFLooseId[i] != 0
				and ggIn.jetPUidFullDiscriminant[i] > 0.62
				# Medium or Tight:
				and	ggIn.jetCHF[i] > 0.
				and	ggIn.jetCEF[i] < 0.99
				and	ggIn.jetNCH[i] > 0.
				and ggIn.jetNHF[i] < 0.99
				and ggIn.jetNEF[i] < 0.99
				# Tight only:
				#and ggIn.jetNHF[i] < 0.90
				#and ggIn.jetNEF[i] < 0.90
				):
			nJets += 1
			evtST += ggIn.jetPt[i]
			evtHT += ggIn.jetPt[i]
	#if nJets < 2 or evtHT < 1000:
	if nJets < 2 or evtHT < 1000 or nJets > 3:
		continue

	# Electron veto
	nEle = 0
	for i in range(ggIn.nEle):
		if (ggIn.elePt[i] > 15.0
				and ggIn.eleIDbit[i]>>3&1 == 1 # >>0:veto, >>1:loose, >>2:medium, >>3:tight
				and abs(ggIn.eleEta[i]) < 2.5
				and abs(ggIn.eleDz[i]) < 0.1
				and ggIn.elePFPUIso[i] < 0.1
				):
			nEle += 1
	if nEle != 0:
		continue 

	# Muon veto
	nMu = 0
	for i in range(ggIn.nMu):
		if (ggIn.muPt[i] > 15.0
				and ggIn.muIsTightID[i] == 1
				#and ggIn.muIsLooseID[i] == 1
				and ggIn.muPFPUIso[i] < 0.12
				):
			nMu += 1
	if nMu != 0:
		continue

	# MET selection
	if ggIn.pfMET > 15.:
		evtST += ggIn.pfMET

	# Write this evt to output tree
	evtST_[0] = evtST
	nJets_[0] = nJets
	ggOut.Fill()
	nAcc += 1

##### EVENT SELECTION END #####
outFile.Write()
outFile.Close()

sw.Stop()
print " >> nAccepted evts:",nAcc,"/",iEvtEnd-iEvtStart,"(",100.*nAcc/(iEvtEnd-iEvtStart),"% )"
print " >> Real time:",sw.RealTime()/60.,"minutes"
print " >> CPU time: ",sw.CpuTime() /60.,"minutes"
