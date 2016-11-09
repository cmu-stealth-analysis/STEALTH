import os
import sys
import glob
import array
import numpy as np
import ROOT


sw = ROOT.TStopwatch()
sw.Start()

ggIn = ROOT.TChain("ggNtuplizer/EventTree")
#ggIn.Add("root://cmsxrootd.fnal.gov///store/group/phys_smp/ggNtuples/13TeV/data/V07_04_14_00/GoldenJSON/job_DoubleMu_Run2015C_Oct05_miniAOD.root")
#ggIn.Add("~mandrews/eos/cms/store/user/mandrews/2015ggNtuple-CMSSW76/DoubleEG-Run2015C_25ns-16Dec2015-v1-MINIAOD.root")
#ggIn.Add("~mandrews/eos/cms/store/user/mandrews/2015ggNtuple-CMSSW76/DoubleEG-Run2015D-16Dec2015-v1-MINIAOD-part1.root")
#for f in glob.glob("/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/SinglePhoton/Run2016B/160823_172053/0000/*.root"):
#for f in glob.glob("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160823_172053/0000/*.root"):
#	ggIn.Add(f)
#for f in glob.glob("/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/*.root"):
for i in range(201,333):
	ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_"+str(i)+".root")
	print "adding i="+str(i)
#for f in glob.glob("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/*.root"):
#	ggIn.Add(f)
#ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_1.root")
#ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_2.root")

#ggIn.SetBranchStatus("tau*", 0)
nEntries = ggIn.GetEntries()
print "entries="+str(nEntries)

#outFile = ROOT.TFile("SinglePhoton_2016B_evtSt.root", "RECREATE")
outFile = ROOT.TFile("/afs/cern.ch/user/m/mandrews/eos/cms/store/user/mandrews/SinglePhoton/Run2016B/SinglePhoton_2016B_4.root", "RECREATE")
#outFile = ROOT.TFile("DoubleEG_2015D_evtSt.root", "RECREATE")
#outFile = ROOT.TFile("SingleElectron_evtSt.root", "RECREATE")
outDir = outFile.mkdir("ggNtuplizer")
outDir.cd()
ggOut = ggIn.CloneTree()


outFile.Write()
outFile.Close()

sw.Stop()
print "Real time: " + str(sw.RealTime() / 60.0) + " minutes"
print "CPU time:  " + str(sw.CpuTime() / 60.0) + " minutes"
