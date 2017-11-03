#!/bin/bash

WORKDIR="/afs/cern.ch/user/t/tmudholk/public/research/stealth/from_michael/STEALTH"
export X509_USER_PROXY=/afs/cern.ch/user/t/tmudholk/private/x509up_u83667
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/t/tmudholk/public/research/ecaldqm/cmssw/CMSSW_9_2_13/src && eval `scramv1 runtime -sh`
export PYTHONPATH=/afs/cern.ch/user/t/tmudholk/public/tmPyUtils:$PYTHONPATH
cd $WORKDIR

# ./submitJobs_skimKinematic_Data.py --inputFilePath=/eos/cms/store/user/tmudholk/stealth/ggSKIMS/DoubleEG_Run2016${era}_FebReMiniAOD_SKIM_DoubleFake.root --outputFolder=/eos/cms/store/user/tmudholk/stealth/ggKinematicSkims/ --outputFilePrefix=DoubleEG_Run2016${era}_FebReMiniAOD_KinematicSkim_DoubleFake --logDirectory=/afs/cern.ch/user/t/tmudholk/public/research/stealth/from_michael/STEALTH/logs/

./submitJobs_skimKinematic_Data.py --inputFromFile --inputFilePath=/afs/cern.ch/user/t/tmudholk/public/research/stealth/from_michael/STEALTH/inputFileList_JetHT_era${1}.txt --outputFolder=/eos/cms/store/user/tmudholk/stealth/ggKinematicSkims_JetHT/ --outputFilePrefix=HTJet_Run2016${1}_FebReMiniAOD_KinematicSkim_DoubleFake --logDirectory=/afs/cern.ch/user/t/tmudholk/public/research/stealth/from_michael/STEALTH/logs/ --nEvtsPerOutputFile=100000 --queue=8nh
