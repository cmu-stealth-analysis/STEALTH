#!/bin/bash

WORKDIR="/afs/cern.ch/user/t/tmudholk/public/research/stealth/from_michael/STEALTH"
export X509_USER_PROXY=/afs/cern.ch/user/t/tmudholk/private/x509up_u83667
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/t/tmudholk/public/research/ecaldqm/cmssw/CMSSW_9_2_13/src && eval `scramv1 runtime -sh`
export PYTHONPATH=/afs/cern.ch/user/t/tmudholk/public/tmPyUtils:$PYTHONPATH

cd ${WORKDIR} && \
    echo "PWD=${PWD}" && \
    echo "Starting skim for era = ${1}" && \
    echo "Python version: " && python --version && \
    ./skimHLTPaths_forFakePhotons.py -e ${1} 2>&1 | tee ${WORKDIR}/skimLog_era${1}.txt
