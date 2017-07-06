import os
import sys
import numpy as np
import argparse
import ROOT

# Register command line options
parser = argparse.ArgumentParser(description='Photon skimmer')
parser.add_argument('-n','--nPhoCut', default=1, required=False, help='Min number of photons',type=int)
parser.add_argument('-x','--xsec', required=True, help='xsec of decay',type=int)
parser.add_argument('-g','--nGenEvts', required=True, help='Generated events in decay',type=int)
parser.add_argument('-b','--batch', required=True, help='Batch to process',type=int)
args = parser.parse_args()

# Keep time
print " >> Running Loose Photon Skim..."
sw = ROOT.TStopwatch()
sw.Start()

nPhoCut = args.nPhoCut
xsec = args.xsec
nGenEvts = args.nGenEvts
batch = args.batch

# Load input TTrees into TChain
eosDir = "/eos/uscms/store/user/lpcsusystealth/MC"
ggInStr = "%s/DiPhotonJetsBox_M40_80-Sherpa/crab_job_DiPhotonJetsBox_M40_80-Sherpa/170330_160556/0000/ggtree_mc_*.root"%(eosDir)
#ggInStr = "%s/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa/crab_job_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa/170330_160525/000%d/ggtree_mc_*.root"%(eosDir,batch)
#ggInStr = "%s/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8/crab_job_DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8/170323_142746/000%d/ggtree_mc_*.root"%(eosDir,batch)
ggIn = ROOT.TChain("ggNtuplizer/EventTree")
ggIn.Add(ggInStr)
nEvts = ggIn.GetEntries()
print " >> Input file(s):",ggInStr
print " >> nEvts:",nEvts

# Initialize output file as empty clone
eosDir = "/eos/cms/store/user/mandrews/MC/ggSKIMS"
#outFileStr = "%s/DiPhotonJetsBox_M40_80-Sherpa_SKIM_%d.root"%(eosDir,batch)
#outFileStr = "%s/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa_SKIM_%d.root"%(eosDir,batch)
#outFileStr = "%s/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8_SKIM_%d.root"%(eosDir,batch)
outFileStr = 'test.root'
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
iEvtEnd   = nEvts
#iEvtEnd   = 1000 
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

    # Photon skim by number of photons
    nPhotons = 0
    for i in range(ggIn.nPho):
       if (ggIn.phoEt[i] > 20.0 and ggIn.phoIDbit[i]>>0&1 != 0): # >>0:loose, >>1:medium, >>2:tight
           nPhotons += 1
    if nPhotons < nPhoCut:
       continue 

    # Write this evt to output tree
    evtWgt_1_pb[0] = xsec / nGenEvts
    ggOut.Fill()
    nAcc += 1

##### EVENT SELECTION END #####
outFile.Write()
outFile.Close()

sw.Stop()
print " >> nAccepted evts:",nAcc,"/",iEvtEnd-iEvtStart,"(",100.*nAcc/(iEvtEnd-iEvtStart),"% )"
print " >> Real time:",sw.RealTime()/60.,"minutes"
print " >> CPU time: ",sw.CpuTime() /60.,"minutes"
