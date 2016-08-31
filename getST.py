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
#for f in glob.glob("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/*.root"):
#	ggIn.Add(f)
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/SinglePhoton_2016B_1.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_1.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_10.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_100.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_101.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_102.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_103.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_104.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_105.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_106.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_107.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_108.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_109.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_11.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_110.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_111.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_112.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_113.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_114.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_115.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_116.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_117.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_118.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_119.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_12.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_120.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_121.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_122.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_123.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_124.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_125.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_126.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_127.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_128.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_129.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_13.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_130.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_131.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_132.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_133.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_134.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_135.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_136.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_137.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_138.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_139.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_14.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_140.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_141.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_142.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_143.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_144.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_145.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_146.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_147.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_148.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_149.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_15.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_150.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_151.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_152.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_153.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_154.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_155.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_156.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_157.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_158.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_159.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_16.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_160.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_161.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_162.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_163.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_164.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_165.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_166.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_167.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_168.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_169.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_17.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_170.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_171.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_172.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_173.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_174.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_175.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_176.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_177.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_178.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_179.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_18.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_180.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_181.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_182.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_183.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_184.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_185.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_186.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_187.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_188.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_189.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_19.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_190.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_191.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_192.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_193.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_194.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_195.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_196.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_197.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_198.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_199.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_2.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_20.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_200.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_201.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_202.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_203.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_204.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_205.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_206.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_207.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_208.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_209.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_21.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_210.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_211.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_212.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_213.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_214.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_215.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_216.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_217.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_218.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_219.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_22.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_220.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_221.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_222.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_223.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_224.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_225.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_226.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_227.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_228.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_229.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_23.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_230.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_231.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_232.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_233.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_234.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_235.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_236.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_237.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_238.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_239.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_24.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_240.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_241.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_242.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_243.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_244.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_245.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_246.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_247.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_248.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_249.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_25.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_250.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_251.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_252.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_253.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_254.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_255.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_256.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_257.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_258.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_259.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_26.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_260.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_261.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_262.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_263.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_264.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_265.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_266.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_267.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_268.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_269.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_27.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_270.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_271.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_272.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_273.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_274.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_275.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_276.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_277.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_278.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_279.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_28.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_280.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_281.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_282.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_283.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_284.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_285.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_286.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_287.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_288.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_289.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_29.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_290.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_291.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_292.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_293.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_294.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_295.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_296.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_297.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_298.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_299.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_3.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_30.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_300.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_301.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_302.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_303.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_304.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_305.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_306.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_307.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_308.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_309.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_31.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_310.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_311.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_312.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_313.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_314.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_315.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_316.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_317.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_318.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_319.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_32.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_320.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_321.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_322.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_323.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_324.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_325.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_326.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_327.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_328.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_329.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_33.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_330.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_331.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_332.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_34.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_35.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_36.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_37.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_38.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_39.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_4.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_40.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_41.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_42.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_43.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_44.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_45.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_46.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_47.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_48.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_49.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_5.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_50.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_51.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_52.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_53.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_54.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_55.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_56.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_57.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_58.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_59.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_6.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_60.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_61.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_62.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_63.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_64.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_65.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_66.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_67.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_68.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_69.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_7.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_70.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_71.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_72.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_73.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_74.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_75.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_76.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_77.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_78.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_79.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_8.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_80.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_81.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_82.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_83.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_84.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_85.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_86.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_87.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_88.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_89.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_9.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_90.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_91.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_92.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_93.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_94.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_95.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_96.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_97.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_98.root")
ggIn.Add("root://cms-xrd-global.cern.ch//store/user/mandrews/SinglePhoton/Run2016B/160827_150350/0000/ggtree_data_99.root")

#ggIn.SetBranchStatus("tau*", 0)
nEntries = ggIn.GetEntries()

outFile = ROOT.TFile("SinglePhoton_2016B_evtSt.root", "RECREATE")
#outFile = ROOT.TFile("DoubleEG_2015D_evtSt.root", "RECREATE")
#outFile = ROOT.TFile("SingleElectron_evtSt.root", "RECREATE")
outDir = outFile.mkdir("ggNtuplizer")
outDir.cd()
ggOut = ggIn.CloneTree(0)


# evtSt = array.array('f', [0])
nJets_  = np.zeros(1, dtype=int)
evtSt_  = np.zeros(1, dtype=float)
StDiff_ = np.zeros(1, dtype=float)
b_nJets = ggOut.Branch("b_nJets", nJets_, "b_nJets/I")
b_evtSt = ggOut.Branch("b_evtSt", evtSt_, "b_evtSt/D")
b_StDiff = ggOut.Branch("b_StDiff", StDiff_, "b_StDiff/D")

count = 0
print "nEntries: " + str(nEntries) 
#for jEvt in range(500000):
for jEvt in range(nEntries):

    # Initialize event
    iEvt = ggIn.LoadTree(jEvt)
    if iEvt < 0:
        break
    nb = ggIn.GetEntry(jEvt)
    if nb <= 0:
        continue
    if jEvt % 100000 == 0:
        print "Processing entry " + str(jEvt)

    evtSt = 0.

    # Photon selection
    nPhotons = 0
    for i in range(ggIn.nPho):
        if (ggIn.phoEt[i] > 100.0 and 
            abs(ggIn.phoEta[i]) < 1.479 and #isEB
            ggIn.phoIDbit[i]>>1&1 == 1): #medium photonID
            nPhotons += 1
            evtSt += ggIn.phoEt[i]
    if nPhotons != 1:
       continue 

    # Jet selection
    nJets = 0
    sumHt = 0
    for i in range(ggIn.nJet):
        if (ggIn.jetPt[i] > 30.0 and
            #ggIn.jetNHF[i] < 0.90 and
            #ggIn.jetNEF[i] < 0.90 and
            ggIn.jetNHF[i] < 0.99 and
            ggIn.jetNEF[i] < 0.99 and
            ggIn.jetCHF[i] > 0. and
            ggIn.jetCEF[i] < 0.99 and
            ggIn.jetNCH[i] > 0. and
            abs(ggIn.jetEta[i]) < 2.4 and
            ggIn.jetPFLooseId[i] == 1):
            nJets += 1
            evtSt += ggIn.jetPt[i]
            sumHt += ggIn.jetPt[i]
    if nJets < 2 or sumHt < 700:
        continue

    # Electron veto
    nEle = 0
    for i in range(ggIn.nEle):
        if (ggIn.elePt[i] > 15.0 and 
            abs(ggIn.eleEta[i]) < 2.5 and
            ggIn.eleIDbit[i]>>3&1 == 1 and # tight electronID
            abs(ggIn.eleDz[i]) < 0.1 and
            ggIn.elePFPUIso[i] < 0.1):
            nEle += 1
    if nEle != 0:
       continue 

    # Muon veto
    nMu = 0
    for i in range(ggIn.nMu):
        if (ggIn.muPt[i] > 15.0 and 
            ggIn.muIsTightID[i] == 1 and
            ggIn.muPFPUIso[i] < 0.12):
            nMu += 1
    if nMu != 0:
        continue

    # MET selection
    if ggIn.pfMET > 15.:
        evtSt += ggIn.pfMET

    # Fill output tree
    evtSt_[0] = evtSt
    nJets_[0] = nJets
    StDiff_[0] = evtSt - ggIn.pfMETsumEt
    ggOut.Fill()
    count += 1
    #if count % 10 == 0:
    #    print " >> S_T = " + str(evtSt_[0])


outFile.Write()
outFile.Close()

sw.Stop()
print "Real time: " + str(sw.RealTime() / 60.0) + " minutes"
print "CPU time:  " + str(sw.CpuTime() / 60.0) + " minutes"
