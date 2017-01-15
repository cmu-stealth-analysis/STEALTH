import os
import sys
import glob
import ROOT
import re
import time

runEra = "C"
eosDir = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/DATA'

#INPATH = 'root://cms-xrd-global.cern.ch//store/user/mandrews/DATA/JetHT/crab_job_JetHT_Run2016%s_*/1701*/000*/ggtree_data_*.root'%(runEra)
INPATH = '%s/JetHT/crab_job_JetHT_Run2016%s_SepRereco/1701*/000*/ggtree_data_*.root'%(eosDir,runEra)
#INPATH = '%s/JetHT/crab_job_JetHT_Run2016C_SepRereco/170105_214838/000*/ggtree_data_*.root'
#INPATH = '%s/JetHT/crab_job_JetHT_Run2016D_SepRereco/170105_215031/000*/ggtree_data_*.root'
#INPATH = '%s/JetHT/crab_job_JetHT_Run2016%s_PRv*/1701*/000*/ggtree_data_*.root'%(eosDir,runEra)
#INPATH = '%s/JetHT/crab_job_JetHT_Run2016H_PRv2-notFinished/170107_121055/000*/ggtree_data_*.root'
#INPATH = '%s/JetHT/crab_job_JetHT_Run2016H_PRv3/170105_213805/000*/ggtree_data_*.root'
#INPATH = 'root://cms-xrd-global.cern.ch//store/user/mandrews/job_JetHT_Run2016B_SepRereco/JetHT/crab_job_JetHT_Run2016B_SepRereco/161218_020346/000*/ggtree_data*.root'

print " >> Merging files in:",INPATH

# Keep time
sw = ROOT.TStopwatch()
sw.Start()


fDir = []
count = 0
print " >> Input files:",len(glob.glob(INPATH))
for f in glob.glob(INPATH):
	fDir.append(f)
	count += 1
print " >> Input chain:",len(fDir)
if len(fDir) != len(glob.glob(INPATH)):
	print " !! Files in chain do not match files in directory !!"

ggIn = ROOT.TChain("ggNtuplizer/EventTree")
ggIn.SetMaxTreeSize(10000000000) # 100 GB
for _f in fDir:
	f = re.sub('/afs/cern.ch/user/m/mandrews/eos/cms', 'root://cms-xrd-global.cern.ch/', f)
	ggIn.Add(f)
print " >> Input evts:",ggIn.GetEntries()

# Initialize output file as empty clone
outFileStr = "%s/ggSKIMS/JetHT_Run2016%s_SepRereco_HLTPFJet450HT900_SKIM.root"%(eosDir,runEra)
print " >> Output file:",outFileStr
outFile = ROOT.TFile(outFileStr, "RECREATE")
outDir = outFile.mkdir("ggNtuplizer")
outDir.cd()
ggIn.Merge(outFile,0,"fast")
#time.sleep(60)
print " >> Merge done. Checking..."
#ggOut = ROOT.TChain("ggNtuplizer/EventTree")
#ggOut.Add(outFileStr)
#print " >> Output evts:",ggOut.GetEntries()

sw.Stop()
print "Real time: " + str(sw.RealTime() / 60.0) + " minutes"
print "CPU time:  " + str(sw.CpuTime() / 60.0) + " minutes"
