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
parser.add_argument('-r','--DeltaR',  default=0., help='DeltaR(pho,jet) cut.',type=float)
args = parser.parse_args()

# Initialize selection params
nPhoCut_  = 1
PhoEtLead = 20. 
PhoEtAll  = 20. 
if args.sel == 'A':
	nPhoCut_ = 2
	PhoEtLead = 35.
	PhoEtAll  = 25. 
	#PhoEtLead = 10.
	#PhoEtAll  = 10. 
if args.sel == 'C':
	nPhoCut_ = 0
HTcut_   = args.ht
runEra   = args.era
minDeltaR = args.DeltaR
print " >> Running STEALTH 2016 Data Selection:",args.sel
print " >> HT cut:",HTcut_
print " >> Era: 2016%s"%runEra
print " >> DeltaR(pho,jt): %.2f"%minDeltaR

## MAIN ##
def main():

	# Keep time
	sw = ROOT.TStopwatch()
	sw.Start()

	# For DeltaR cut
	phoEta0 = 0.
	phoPhi0 = 0.
	phoEta1 = 0.
	phoPhi1 = 0.

	# Load input TTrees into TChain
	eosDir = "~/eos/cms/store/user/mandrews"
	#eosDir = "/eos/uscms/store/user/mba2012"
	#ggInStr = "%s/MC/ggSKIMS/GJet_*_DoubleEMEnriched_SKIM_2Pho10Loose.root"%(eosDir)
	#ggInStr = "%s/MC/ggSKIMS/QCD_*_DoubleEMEnriched_SKIM_2Pho10Loose.root"%(eosDir)
	#ggInStr = "%s/DATA/ggSKIMS/JetHT_Run2016%s_SepRereco_HLTPFJet450HT900_SKIM.root"%(eosDir,runEra)
	#ggInStr = "%s/DATA/stNTUPLES/DoubleEG_SepRereco_Run2016%s_selA_HT60_DeltaR00.root"%(eosDir,runEra)
	ggInStr = "%s/DATA/stNTUPLES/JetHT_SepRereco_Run2016%s_selB_HT800_DeltaR00.root"%(eosDir,runEra)
	#ggInStr = "%s/DATA/ggSKIMS/DoubleEG_Run2016%s_SepRereco_HLTDiPho3018M90_SKIM.root"%(eosDir,runEra)
	ggIn = ROOT.TChain("ggNtuplizer/EventTree")
	ggIn.Add(ggInStr)
	nEvts = ggIn.GetEntries()
	print " >> Input file(s):",ggInStr
	print " >> nEvts:",nEvts

	# Initialize output file as empty clone
	outFileStr = "%s/DATA/stNTUPLES/JetHT_SepRereco_Run2016%s_sel%s_HT%d_DeltaR%02d.root"%(eosDir,runEra,args.sel,HTcut_,minDeltaR*10.)
	#outFileStr = "%s/DATA/stNTUPLES/DoubleEG_SepRereco_Run2016%s_sel%s_HT%d_DeltaR%02d.root"%(eosDir,runEra,args.sel,HTcut_,minDeltaR*10.)
	#outFileStr = "%s/MC/stNTUPLES/GJet_sel%s_HT%d_Pho10_DeltaR%02d.root"%(eosDir,args.sel,HTcut_,minDeltaR*10.)
	#outFileStr = "%s/MC/stNTUPLES/QCD_sel%s_HT%d_Pho10_DeltaR%02d.root"%(eosDir,args.sel,HTcut_,minDeltaR*10.)
	#outFileStr = "test.root"
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
	#massCut_ = 90.

	##### EVENT SELECTION START #####
	nAcc = 0
	iEvtStart = 0
	iEvtEnd   = nEvts
	#iEvtEnd   = 10000
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
		if ggIn.isData and args.sel == 'A' and ggIn.HLTPho>>14&1 == False: # HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90
		  continue
		nPhotons = 0
		for i in range(ggIn.nPho):
			if (ggIn.phoEt[0] < PhoEtLead):
				continue
			if not (ggIn.phoEt[i] > PhoEtAll
					and ggIn.phoIDbit[i]>>1&1 == 1 # >>0:loose, >>1:medium, >>2:tight
					and ggIn.phoEleVeto[i] == True
					#and abs(ggIn.phoEta[i]) < 1.479 # isEB
					and abs(ggIn.phoEta[i]) < 1.442 # isEB
					):
				continue
			nPhotons += 1
			evtST += ggIn.phoEt[i]
		if nPhotons != nPhoCut_:
			continue 
		phoEta0 = ggIn.phoEta[0]
		phoPhi0 = ggIn.phoPhi[0]
		if args.sel == 'A':
			phoEta1 = ggIn.phoEta[1]
			phoPhi1 = ggIn.phoPhi[1]
		#vPho0.SetPtEtaPhiE(ggIn.phoEt[0],ggIn.phoEta[0],ggIn.phoPhi[0],ggIn.phoE[0])
		#vPho1.SetPtEtaPhiE(ggIn.phoEt[1],ggIn.phoEta[1],ggIn.phoPhi[1],ggIn.phoE[1])
		#if (vPho0+vPho1).M() < massCut_:
		#	continue
		
		# Jet selection
		if ggIn.isData and args.sel == 'B' and ggIn.HLTJet>>33&1 == False and ggIn.HLTJet>>18&1 == False: # HLT_PFHT900,PFJet450, resp.
			continue
		nJets = 0
		evtHT = 0
		for i in range(ggIn.nJet):
			if not (ggIn.jetPt[i] > 30.0
					and abs(ggIn.jetEta[i]) < 2.4
					and ggIn.jetPFLooseId[i] != False # for some reason == True doesnt work
					and ggIn.jetPUidFullDiscriminant[i] > 0.62
					# Medium or Tight:
					and ggIn.jetCHF[i] > 0.
					and ggIn.jetCEF[i] < 0.99
					and ggIn.jetNCH[i] > 0.
					#and ggIn.jetNHF[i] < 0.99
					#and ggIn.jetNEF[i] < 0.99
					# Tight only:
					and ggIn.jetNHF[i] < 0.90
					and ggIn.jetNEF[i] < 0.90
					):
				continue
			if np.hypot(phoEta0-ggIn.jetEta[i],phoPhi0-ggIn.jetPhi[i]) < minDeltaR: #DeltaR(pho[0],jet[i])
				continue
			if args.sel == 'A' and np.hypot(phoEta1-ggIn.jetEta[i],phoPhi1-ggIn.jetPhi[i]) < minDeltaR: #DeltaR(pho[1],jet[i])
				continue
			nJets += 1
			evtST += ggIn.jetPt[i]
			evtHT += ggIn.jetPt[i]
		#if evtHT > HTcut_:
		#	print ggIn.event,nJets,evtHT
		if nJets < 2 or evtHT < HTcut_:
		#if nJets < 2 or evtHT < HTcut_ or nJets > 3:
			continue

		# Electron veto
		nEle = 0
		for i in range(ggIn.nEle):
			if not (ggIn.elePt[i] > 15.0
					and ggIn.eleIDbit[i]>>3&1 == True # >>0:veto, >>1:loose, >>2:medium, >>3:tight
					and abs(ggIn.eleEta[i]) < 2.5
					and abs(ggIn.eleDz[i]) < 0.1
					and ggIn.elePFPUIso[i] < 0.1
					):
				continue
			nEle += 1
		if nEle != 0:
			continue 

		# Muon veto
		nMu = 0
		for i in range(ggIn.nMu):
			if not (ggIn.muPt[i] > 15.0
					and ggIn.muIDbit[i]>>2&1 == True # >>0:loose, >>1:med, >>2:tight, >>3:soft, >>4:highpT
					#and ggIn.muIsTightID[i] == True
					#and ggIn.muIsLooseID[i] == True
					and ggIn.muPFPUIso[i] < 0.12
					):
				continue
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
