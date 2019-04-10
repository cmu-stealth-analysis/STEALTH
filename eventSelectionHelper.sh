#!/bin/bash

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
echo "PWD=${PWD}" && echo "Starting event selection" && ./eventSelection/bin/runEventSelection inputFilesList=${1} outputFilePath=${2} isMC=${3} counterStartInclusive=${4} counterEndInclusive=${5} photonSelectionType=${6} year=${7}

OUTDIR=root://cmseos.fnal.gov//store/user/lpcsusystealth/${8}
echo "Copying main output..."
xrdcp -f ${2} ${OUTDIR}/${FILE} 2>&1
XRDEXIT=$?
if [[ $XRDEXIT -ne 0 ]]; then
    rm *.root
    echo "exit code $XRDEXIT, failure in xrdcp"
    exit $XRDEXIT
fi
echo "Finished copying main output!"
rm ${2}

if [ "${3}" = "true" ]; then
    if [ "${6}" = "medium" ]; then
        echo "Copying MC statistics output..."
        OUTDIR=root://cmseos.fnal.gov//store/user/lpcsusystealth/MCSelectionStatistics
        xrdcp -f MCStatisticsDetails.root ${OUTDIR}/MCSelectionStatistics_optimized${7}_${6}_begin_${4}_end_${5}.root 2>&1
        XRDEXIT=$?
        if [[ $XRDEXIT -ne 0 ]]; then
            rm *.root
            echo "exit code $XRDEXIT, failure in xrdcp"
            exit $XRDEXIT
        fi
        echo "Finished copying MC statistics output!"
    fi
    rm MCStatisticsDetails.root
else
    echo "Copying trigger efficiency output..."
    OUTDIR=root://cmseos.fnal.gov//store/user/lpcsusystealth/triggerStatistics
    xrdcp -f triggerEfficiencyRawEventCounters.root ${OUTDIR}/triggerEfficiencyRawEventCounters_cuts${7}_${6}_begin_${4}_end_${5}.root 2>&1
    XRDEXIT=$?
    if [[ $XRDEXIT -ne 0 ]]; then
        rm *.root
        echo "exit code $XRDEXIT, failure in xrdcp"
        exit $XRDEXIT
    fi
    echo "Finished copying trigger efficiency output!"
    rm triggerEfficiencyRawEventCounters.root
fi
echo "All done!"
set +x
