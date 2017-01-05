import os
import sys
import numpy as np
import argparse
import ROOT

# Register command line options
parser = argparse.ArgumentParser(description='Run STEALTH selection.')
parser.add_argument('-s','--sel', default='B', help='Selection scheme: A or B.',type=str)
parser.add_argument('-e','--era', default='', help='Run Era.',type=str)
parser.add_argument('-H','--ht',  default=60., help='HT cut.',type=float)
args = parser.parse_args()

# Initialize selection params
nPhoCut_  = 1
PhoEtLead = 20. 
PhoEtAll  = 20. 
doPhoTrg  = False
if args.sel == 'A':
	nPhoCut_ = 2
	PhoEtLead = 35.
	PhoEtAll  = 25. 
	doPhoTrg = True
if args.sel == 'C':
	nPhoCut_ = 0
HTcut_   = args.ht
runEra   = args.era
massCut_ = 90.
print " >> Running STEALTH 2016 Data Selection:",args.sel
print " >> HT cut:",HTcut_
print " >> Era: 2016%s"%runEra

## MAIN ##
def main():

	# Keep time
	sw = ROOT.TStopwatch()
	sw.Start()

	# Load input TTrees into TChain
	#ggInStr = "~/eos/cms/store/user/mandrews/DATA/JetHT/JetHT_2016%s_Pho20Loose_*.root"%runEra
	#ggInStr = "~/eos/cms/store/user/mandrews/DATA/JetHT/JetHT_2016%s_SKIM_1Pho20Loose_*.root"%runEra
	#ggInStr = "~/eos/cms/store/user/mandrews/MC/MC_QCD_St*_SKIM_1Pho20Loose_*.root"
	#ggInStr = "~/eos/cms/store/user/mandrews/MC/GJet_*_DoubleEMEnriched_*_SKIM_2Pho25Loose_*.root"
	ggInStr = "~/eos/cms/store/user/mandrews/DATA/JetHT_SepRereco/JetHT_Run2016%s_SepRereco_HLTPFHT200250900_Merge.root"%runEra
	ggIn = ROOT.TChain("ggNtuplizer/EventTree")
	ggIn.Add(ggInStr)
	nEvts = ggIn.GetEntries()
	print " >> Input file(s):",ggInStr
	print " >> nEvts:",nEvts

	# Initialize output file as empty clone
	#outFileStr = "~/eos/cms/store/user/mandrews/stNTUPLES/DATA/JetHT_SepRereco_Run2016%s_sel%s_HT%d.root"%(runEra,args.sel,HTcut_)
	#outFileStr = "~/eos/cms/store/user/mandrews/stNTUPLES/MC/QCD_sel%s_HT%d.root"%(args.sel,HTcut_)
	#outFileStr = "~/eos/cms/store/user/mandrews/stNTUPLES/MC/GJet_sel%s_HT%d.root"%(args.sel,HTcut_)
	outFileStr = "test.root"
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

	# For invariant mass cut
	#vPho0 = ROOT.TLorentzVector()
	#vPho1 = ROOT.TLorentzVector()

	##### EVENT SELECTION START #####
	nAcc = 0
	iEvtStart = 0
	iEvtEnd   = nEvts
	#iEvtEnd   = 50000
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
		if doPhoTrg and (ggIn.HLTPho>>14)&1 == 0: # HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90
		  continue
		nPhotons = 0
		for i in range(ggIn.nPho):
			if (ggIn.phoEt[0] < PhoEtLead):
				continue
			if (ggIn.phoEt[i] > PhoEtAll
					and ggIn.phoIDbit[i]>>1&1 == 1 # >>0:loose, >>1:medium, >>2:tight
					and ggIn.phoEleVeto[i] == True
					#and abs(ggIn.phoEta[i]) < 1.479 # isEB
					and abs(ggIn.phoEta[i]) < 1.442 # isEB
					):
				nPhotons += 1
				evtST += ggIn.phoEt[i]
		if nPhotons != nPhoCut_:
			continue 
		#vPho0.SetPtEtaPhiE(ggIn.phoEt[0],ggIn.phoEta[0],ggIn.phoPhi[0],ggIn.phoE[0])
		#vPho1.SetPtEtaPhiE(ggIn.phoEt[1],ggIn.phoEta[1],ggIn.phoPhi[1],ggIn.phoE[1])
		#if (vPho0+vPho1).M() < massCut_:
		#	continue
		
		# Jet selection
		if (not doPhoTrg) and (ggIn.HLTJet>>33)&1 == 0: # HLT_PFHT900
			continue
		nJets = 0
		evtHT = 0
		for i in range(ggIn.nJet):
			if (ggIn.jetPt[i] > 30.0
					and abs(ggIn.jetEta[i]) < 2.4
					and ggIn.jetPFLooseId[i] != False # for some reason == True doesnt work
					and ggIn.jetPUidFullDiscriminant[i] > 0.62
					# Medium or Tight:
					and	ggIn.jetCHF[i] > 0.
					and	ggIn.jetCEF[i] < 0.99
					and	ggIn.jetNCH[i] > 0.
					#and ggIn.jetNHF[i] < 0.99
					#and ggIn.jetNEF[i] < 0.99
					# Tight only:
					and ggIn.jetNHF[i] < 0.90
					and ggIn.jetNEF[i] < 0.90
					):
				nJets += 1
				evtST += ggIn.jetPt[i]
				evtHT += ggIn.jetPt[i]
		if nJets < 2 or evtHT < HTcut_:
		#if nJets < 2 or evtHT < HTcut_ or nJets > 3:
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
					and ggIn.muIDbit[i]>>2&1 == 1 # >>0:loose, >>1:med, >>2:tight, >>3:soft, >>4:highpT
					#and ggIn.muIsTightID[i] == 1
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

#_____ Call main() ______#
if __name__ == '__main__':
	  main()
