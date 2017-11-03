#!/bin/bash

WORKDIR="/afs/cern.ch/user/t/tmudholk/public/research/stealth/from_michael/STEALTH"
export X509_USER_PROXY=/afs/cern.ch/user/t/tmudholk/private/x509up_u83667
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/t/tmudholk/public/research/ecaldqm/cmssw/CMSSW_9_2_13/src && eval `scramv1 runtime -sh`
export PYTHONPATH=/afs/cern.ch/user/t/tmudholk/public/tmPyUtils:$PYTHONPATH

INPUT_FROM_FILE_FLAG=""
if [ "$6" == "inputFromFile" ]; then
    INPUT_FROM_FILE_FLAG=" --inputFromFile"
fi

cd ${WORKDIR} && \
    echo "PWD=${PWD}" && \
    echo "Starting kinematic skim" && \
    echo "Python version: " && python --version && \
    ./skimKinematic_Data.py${INPUT_FROM_FILE_FLAG} --inputFilePath=${1} --outputFilePath=${2} --counterStartInclusive=${3} --counterEndInclusive=${4} 2>&1 | tee ${5}
