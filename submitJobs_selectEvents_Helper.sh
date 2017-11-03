#!/bin/bash

WORKDIR="/afs/cern.ch/user/t/tmudholk/public/research/stealth/from_michael/STEALTH"
cd ${WORKDIR}
source setupEnv.sh

echo "PWD=${PWD}" && \
    echo "Starting event selection" && \
    echo "Python version: " && python --version && \
    ./selectEvents.py --inputFilePath=${1} --outputFilePath=${2} --counterStartInclusive=${3} --counterEndInclusive=${4} --photonSelectionType=${5} 2>&1 | tee ${6}
