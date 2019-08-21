#!/bin/bash

function xrdmv_with_check {
    if [ "${#}" != 2  ]; then
        echo "ERROR: number of arguments passed to \"${FUNCNAME}\": ${#}"
        exit 1
    fi
    xrdcp -f ${1} ${2} 2>&1
    XRDEXIT=${?}
    if [[ ${XRDEXIT} -ne 0 ]]; then
        rm -f *.root
        echo "exit code ${XRDEXIT}, failure in xrdcp"
        exit ${XRDEXIT}
    fi
    rm ${1}
}

cd ${_CONDOR_SCRATCH_DIR}

# Source CMSSW environment
echo "Sourcing CMSSW environment..."
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc630
eval `scramv1 project CMSSW CMSSW_9_2_13`
cd CMSSW_9_2_13/src/
eval `scramv1 runtime -sh` # cmsenv is not an alias on the workers
echo "CMSSW: "${CMSSW_BASE}
cd ${_CONDOR_SCRATCH_DIR}

echo "Extracting tmUtils tarball..."
./extract_tmUtilsTarball.sh && rm -f tmUtils.tar.gz

echo "Setting env variables..."
export PYTHONPATH=${_CONDOR_SCRATCH_DIR}/tmPyUtils:${PYTHONPATH}
export TMPYUTILS=${_CONDOR_SCRATCH_DIR}/tmPyUtils/
export TMCPPUTILS=${_CONDOR_SCRATCH_DIR}/tmCPPUtils/

echo "Extracting event selection tarball..."
cd ${_CONDOR_SCRATCH_DIR}
./extract_eventSelectionTarball.sh && rm -f eventSelection.tar.gz

cd ${_CONDOR_SCRATCH_DIR}
# echo "Current directory structure: (excluding CMSSW)"
# ls -I "CMSSW*" -R

set -x
echo "PWD=${PWD}" && echo "Starting event selection" && ./eventSelection/bin/runEventSelection inputFilesList=${1} isMC=${2} disableJetSelection=${3} counterStartInclusive=${4} counterEndInclusive=${5} year=${6}

ALLJETSPREFIX=""
if [ "${3}" == "true" ]; then
    ALLJETSPREFIX="_allJets"
fi

EOSPREFIX=root://cmseos.fnal.gov/
echo "Copying selections..."
xrdmv_with_check selection_signal.root ${EOSPREFIX}${7}/selection_${9}${ALLJETSPREFIX}_${6}_signal_begin_${4}_end_${5}.root
xrdmv_with_check selection_control_fakefake.root ${EOSPREFIX}${7}/selection_${9}${ALLJETSPREFIX}_${6}_control_fakefake_begin_${4}_end_${5}.root
xrdmv_with_check selection_control_mediumfake.root ${EOSPREFIX}${7}/selection_${9}${ALLJETSPREFIX}_${6}_control_mediumfake_begin_${4}_end_${5}.root
echo "Finished copying selections!"

echo "Copying statistics histograms..."
xrdmv_with_check statisticsHistograms.root ${EOSPREFIX}${8}/statistics_${9}${ALLJETSPREFIX}_${6}_begin_${4}_end_${5}.root
echo "Finished copying statistics!"

echo "All done!"
set +x
