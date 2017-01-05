import os
import sys
import glob
import ROOT
import time

#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016E_SepRereco/JetHT/crab_job_JetHT_Run2016E_SepRereco/161215_172950/*/ggtree_data*.root'
#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016G_SepRereco/JetHT/crab_job_JetHT_Run2016G_SepRereco/161216_100555/000*/ggtree_data*.root'
#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016G_SepRereco/JetHT/crab_job_JetHT_Run2016G_SepRereco/161216_100555/0001/ggtree_data*.root'
#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016G_SepRereco/JetHT_Run2016G_SepRereco_HLTPFHT200250900_Merge*.root'
#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016D_SepRereco/JetHT/crab_job_JetHT_Run2016D_SepRereco/161219_170650/000*/ggtree_data*.root'
#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016F_SepRereco_1/JetHT/crab_job_JetHT_Run2016F_SepRereco_1/161221_115728/0000/ggtree_data*.root'
#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016F_SepRereco_1-Missing/JetHT/crab_job_JetHT_Run2016F_SepRereco_1-Missing/161226_104136/0000/ggtree_data*.root'
#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016F_SepRereco_2/JetHT/crab_job_JetHT_Run2016F_SepRereco_2/161221_115841/0000/ggtree_data*.root'
#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016H_PRv2/JetHT/crab_job_JetHT_Run2016H_PRv2/161220_145307/000*/ggtree_data*.root'
#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016H_PRv3/JetHT/crab_job_JetHT_Run2016H_PRv3/161220_145447/0000/ggtree_data*.root'
INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016B_SepRereco/JetHT/crab_job_JetHT_Run2016B_SepRereco/161218_020346/000*/ggtree_data*.root'
#INPATH = 'root://cms-xrd-global.cern.ch//store/user/mandrews/job_JetHT_Run2016B_SepRereco/JetHT/crab_job_JetHT_Run2016B_SepRereco/161218_020346/000*/ggtree_data*.root'
#INPATH = '/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016C_SepRereco/JetHT/crab_job_JetHT_Run2016C_SepRereco/161219_164839/0000/ggtree_data*.root'

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
for f in fDir:
	ggIn.Add(f)
print " >> Input evts:",ggIn.GetEntries()

# Initialize output file as empty clone
#outFileStr = "/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016E_SepRereco/JetHT_Run2016E_SepRereco_HLTPFHT200250900_Merge.root"
#outFileStr = "/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016G_SepRereco/JetHT_Run2016G_SepRereco_HLTPFHT200250900.root"
#outFileStr = "/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016D_SepRereco/JetHT_Run2016D_SepRereco_HLTPFHT200250900_Merge.root"
#outFileStr = "/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016F_SepRereco_1/JetHT_Run2016F_SepRereco_HLTPFHT200250900_Merge.root"
#outFileStr = "/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016H_PRv2/JetHT_Run2016H_SepRereco_HLTPFHT200250900_Merge.root"
outFileStr = "/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016B_SepRereco/JetHT_Run2016B_SepRereco_HLTPFHT200250900_Merge.root"
#outFileStr = "root://cms-xrd-global.cern.ch//store/user/mandrews/job_JetHT_Run2016B_SepRereco/JetHT_Run2016B_SepRereco_HLTPFHT200250900_Merge.root"
#outFileStr = "/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/job_JetHT_Run2016C_SepRereco/JetHT_Run2016C_SepRereco_HLTPFHT200250900_Merge.root"
print " >> Output file:",outFileStr
outFile = ROOT.TFile(outFileStr, "RECREATE")
outDir = outFile.mkdir("ggNtuplizer")
outDir.cd()
ggIn.Merge(outFile,0,"fast")
time.sleep(60)
print " >> Waiting..."
ggOut = ROOT.TChain("ggNtuplizer/EventTree")
ggOut.Add(outFileStr)
print " >> Output evts:",ggOut.GetEntries()

sw.Stop()
print "Real time: " + str(sw.RealTime() / 60.0) + " minutes"
print "CPU time:  " + str(sw.CpuTime() / 60.0) + " minutes"
