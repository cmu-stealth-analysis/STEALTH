#!/bin/bash

function xrdmv_with_check {
    if [ "${#}" != 2  ]; then
        echo "ERROR: number of arguments passed to \"${FUNCNAME}\": ${#}"
        exit 1
    fi
    xrdcp -f ${1} ${2} 2>&1
    XRDEXIT=${?}
    if [[ ${XRDEXIT} -ne 0 ]]; then
        rm *.root
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
echo "PWD=${PWD}" && echo "Starting event selection" && ./eventSelection/bin/runEventSelection inputFilesList=${1} isMC=${2} counterStartInclusive=${3} counterEndInclusive=${4} year=${5}

MCDATAPREFIX="data"
if [ "${2}" == "true" ]; then
    MCDATAPREFIX="MC"
fi

EOSPREFIX=root://cmseos.fnal.gov/
echo "Copying selections..."
xrdmv_with_check selection_signal.root ${EOSPREFIX}${6}/selection_${MCDATAPREFIX}_${5}_signal_begin_${3}_end_${4}.root
xrdmv_with_check selection_control_fakefake.root ${EOSPREFIX}${6}/selection_${MCDATAPREFIX}_${5}_control_fakefake_begin_${3}_end_${4}.root
xrdmv_with_check selection_control_mediumfake.root ${EOSPREFIX}${6}/selection_${MCDATAPREFIX}_${5}_control_mediumfake_begin_${3}_end_${4}.root
echo "Finished copying selections!"

echo "Copying statistics histograms..."
xrdmv_with_check statisticsHistograms.root ${EOSPREFIX}${7}/statistics_${MCDATAPREFIX}_${5}_begin_${3}_end_${4}.root
echo "Finished copying statistics!"

echo "All done!"
set +x
