import os
import sys
import numpy as np
import argparse
import ROOT

from tmProgressBar import tmProgressBar

# Register command line options
parser = argparse.ArgumentParser(description='Photon skimmer')
parser.add_argument('-e','--era', required=True, help='Run Era',type=str)
# parser.add_argument('-n','--nPhoCut', default=1, required=False, help='Min number of photons',type=int)
args = parser.parse_args()

# Keep time
print " >> Running Loose Photon Skim..."
sw = ROOT.TStopwatch()
sw.Start()

runEra = args.era
# runEraList = ['B', 'C', 'D', 'E', 'F']
# nPhoCut = args.nPhoCut

listOfInputFiles = []
inputFileNamesFileObject = open('inputFiles_2016%s.txt'%(runEra), 'r')
for inputFileName in inputFileNamesFileObject:
    listOfInputFiles.append(inputFileName.strip())
inputFileNamesFileObject.close()
  

# Load input TTrees into TChain
eosDir = "/eos/cms/store/group/phys_smp/ggNtuples/13TeV/data/V08_00_26_04"
# eosDir = "root://cmseos.fnal.gov//store/user/lpcsusystealth/DATA"
# ggInStr = "%s/job_DoubleEG_Run2016%s_FebReminiAOD/ggtree_data_*.root"%(eosDir,runEra)
# eosDir = "root://cmsxrootd.fnal.gov///store/user/lpcsusystealth/DATA"
# ggInStr = "%s/JetHT_FebReminiAOD/crab_job_JetHT_Run2016%s_FebReminiAOD*/*/000*/ggtree_data_*.root"%(eosDir,runEra)
ggIn = ROOT.TChain("ggNtuplizer/EventTree")
# ggIn.Add(ggInStr)
# for runEra in runEraList:
#     ggInStr = "%s/job_DoubleEG_Run2016%s_FebReminiAOD*/DoubleEG/crab_job_DoubleEG_Run2016%s_FebReminiAOD*/*/*/ggtree_data_*.root"%(eosDir,runEra,runEra)
#     # ggInStr = "%s/JetHT_FebReminiAOD/crab_job_JetHT_Run2016%s_FebReminiAOD*/*/000*/ggtree_data_*.root"%(eosDir,runEra)
#     ggIn.Add(ggInStr)
# ggIn.Add("%s/job_DoubleEG_Run2016*_FebReminiAOD*/DoubleEG/crab_job_DoubleEG_Run2016*_FebReminiAOD*/*/*/ggtree_data_*.root"%(eosDir))
for inputFile in listOfInputFiles:
    print "Adding...", inputFile
    ggIn.Add(inputFile)
nEvts = ggIn.GetEntries()
print " >> nEvts:",nEvts

# Initialize output file as empty clone
eosDir = "/eos/cms/store/user/tmudholk/stealth/ggSKIMS"
outFileStr = "%s/DoubleEG_Run2016%s_FebReMiniAOD_SKIM_DoubleFake.root"%(eosDir, runEra)
# outFileStr = "%s/JetHT_Runs2016_ReminiAOD_HLTPFJet450HT900_SKIM.root"%(eosDir)
# outFileStr = 'test.root'
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
progressBar = tmProgressBar(nEvts)
progressBar.initializeTimer()
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
    if jEvt % 1000 == 0:
        # print " .. Processing entry",jEvt
        progressBar.updateBar(1.0*jEvt/nEvts, jEvt)
    

    # Photon skim by trigger path
    # if (ggIn.HLTPho>>14)&1 == 0: # HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90
    #     continue

    # Photon skim criteria
    nFakePhotonsPassingHigherCut = 0
    nFakePhotonsPassingLowerCut = 0
    for i in range(ggIn.nPho):
        # if (ggIn.phoEt[i] > 20.0 and ggIn.phoIDbit[i]>>0&1 != 0): # >>0:loose, >>1:medium, >>2:tight
        #     nPhotons += 1
        if (ggIn.phoIDbit[i]>>1&1 == 1): # if photon passes medium ID...
            continue                     # then veto the event
        photonPT = ggIn.phoEt[i]
        if (ggIn.phoHoverE[i] > 0.0396 and ggIn.phoSigmaIEtaIEtaFull5x5[i] > 0.01022 and ggIn.phoPFChIso[i] < 1.295 and ggIn.phoPFNeuIso[i] < (5.931 + 0.0163*photonPT + 0.000014*photonPT*photonPT) and ggIn.phoPFPhoIso[i] < (3.630 + 0.0047*photonPT)):
            if (photonPT > 35.):
                nFakePhotonsPassingHigherCut += 1
            if (photonPT > 25.):
                nFakePhotonsPassingLowerCut += 1

    if not(nFakePhotonsPassingLowerCut == 2 and nFakePhotonsPassingHigherCut >= 1):
        continue

    # Write this evt to output tree
    ggOut.Fill()
    nAcc += 1

print " >> finished loop"
##### EVENT SELECTION END #####
outFile.Write()
outFile.Close()

sw.Stop()
print " >> nAccepted evts:",nAcc,"/",iEvtEnd-iEvtStart,"(",100.*nAcc/(iEvtEnd-iEvtStart),"% )"
print " >> Real time:",sw.RealTime()/60.,"minutes"
print " >> CPU time: ",sw.CpuTime() /60.,"minutes"
