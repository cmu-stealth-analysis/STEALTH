import os
import sys
import numpy as np
import argparse
import ROOT

# Register command line options
parser = argparse.ArgumentParser(description='Photon skimmer')
parser.add_argument('-e','--era', required=True, help='Run Era',type=str)
parser.add_argument('-n','--nPhoCut', default=1, required=False, help='Min number of photons',type=int)
args = parser.parse_args()

# Keep time
print " >> Running Loose Photon Skim..."
sw = ROOT.TStopwatch()
sw.Start()

runEra = args.era
nPhoCut = args.nPhoCut

# Load input TTrees into TChain
eosDir = "/eos/cms/store/group/phys_smp/ggNtuples/13TeV/data/V08_00_26_01"
ggInStr = "%s/job_DoubleEG_Run2016%s_FebReminiAOD/ggtree_data_*.root"%(eosDir,runEra)
#eosDir = "/eos/uscms/store/user/lpcsusystealth/DATA"
#ggInStr = "%s/JetHT_FebReminiAOD/crab_job_JetHT_Run2016%s_FebReminiAOD*/*/000*/ggtree_data_*.root"%(eosDir,runEra)
ggIn = ROOT.TChain("ggNtuplizer/EventTree")
ggIn.Add(ggInStr)
nEvts = ggIn.GetEntries()
print " >> Input file(s):",ggInStr
print " >> nEvts:",nEvts

# Initialize output file as empty clone
eosDir = "/eos/cms/store/user/mandrews/DATA/ggSKIMS"
#outFileStr = "%s/DoubleEG_Run2016%s_ReminiAOD_HLTDiPho3018M90_SKIM.root"%(eosDir,runEra)
#outFileStr = "%s/JetHT_Run2016%s_ReminiAOD_HLTPFJet450HT900_SKIM.root"%(eosDir,runEra)
outFileStr = 'test.root'
outFile = ROOT.TFile(outFileStr, "RECREATE")
outDir = outFile.mkdir("ggNtuplizer")
outDir.cd()
ggOut = ggIn.CloneTree(0) 
print " >> Output file:",outFileStr


##### EVENT SELECTION START #####

# Event range to process
iEvtStart = 0
iEvtEnd   = nEvts
#iEvtEnd   = 10000 
print " >> Processing entries: [",iEvtStart,"->",iEvtEnd,")"

nAcc = 0
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

    # Photon skim by trigger path
    if (ggIn.HLTPho>>14)&1 == 0: # HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90
        continue

    # Photon skim by number of photons
    #nPhotons = 0
    #for i in range(ggIn.nPho):
    #   if (ggIn.phoEt[i] > 20.0 and ggIn.phoIDbit[i]>>0&1 != 0): # >>0:loose, >>1:medium, >>2:tight
    #       nPhotons += 1
    #if nPhotons < nPhoCut:
    #   continue 

    # Write this evt to output tree
    ggOut.Fill()
    nAcc += 1

##### EVENT SELECTION END #####
outFile.Write()
outFile.Close()

sw.Stop()
print " >> nAccepted evts:",nAcc,"/",iEvtEnd-iEvtStart,"(",100.*nAcc/(iEvtEnd-iEvtStart),"% )"
print " >> Real time:",sw.RealTime()/60.,"minutes"
print " >> CPU time: ",sw.CpuTime() /60.,"minutes"
