import os
import sys
import numpy as np
import ROOT

# Keep time
print " >> Running Loose Photon Skim..."
sw = ROOT.TStopwatch()
sw.Start()

# Load input TTrees into TChain
ggInStr = "/eos/uscms/store/user/mba2012/JetHT/JetHT_2016G_Pho30Loose_*.root"
ggIn = ROOT.TChain("ggNtuplizer/EventTree")
ggIn.Add(ggInStr)
nEvts = ggIn.GetEntries()
isMC = False
xsec = 108000000.
print " >> Input file(s):",ggInStr
print " >> nEvts:",nEvts
if isMC:
	print " >> isMC, xsec:",xsec

# Initialize output file as empty clone
outFileStr = "test.root"
outFile = ROOT.TFile(outFileStr, "RECREATE")
outDir = outFile.mkdir("ggNtuplizer")
outDir.cd()
ggOut = ggIn.CloneTree(0) 
print " >> Output file:",outFileStr

# Initialize output branches 
evtWgt_1_pb = np.zeros(1, dtype=float)
b_evtWgt_1_pb = ggOut.Branch('b_evtWgt_1_pb', evtWgt_1_pb, 'evtWgt_1_pb/D')

##### EVENT SELECTION START #####
nAcc = 0
iEvtStart = 0
#iEvtEnd   = nEvts
iEvtEnd   = 20000
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

	# Photon selection
	nPhotons = 0
	for i in range(ggIn.nPho):
		if (ggIn.phoEt[i] > 15.0 and ggIn.phoIDbit[i]>>0&1 != 0): # >>0:loose, >>1:medium, >>2:tight
			nPhotons += 1
	if nPhotons < 1:
		continue 

	# Write this evt to output tree
	if isMC:
		evtWgt_1_pb[0] = xsec / nEntries
	else:
		evtWgt_1_pb[0] = 1.
	ggOut.Fill()
	nAcc += 1

##### EVENT SELECTION END #####
outFile.Write()
outFile.Close()

sw.Stop()
print " >> nAccepted evts:",nAcc,"/",iEvtEnd-iEvtStart,"(",100.*nAcc/(iEvtEnd-iEvtStart),"% )"
print " >> Real time:",sw.RealTime()/60.,"minutes"
print " >> CPU time: ",sw.CpuTime() /60.,"minutes"
