#!/bin/bash

function xrdmv_with_check {
    if [ "${#}" != 2  ]; then
        echo "ERROR: number of arguments passed to \"${FUNCNAME}\": ${#}"
        exit 1
    fi
    xrdcp --verbose --force --path --streams 15 ${1} ${2} 2>&1
    XRDEXIT=${?}
    if [[ ${XRDEXIT} -ne 0 ]]; then
        rm -f *.root
        echo "exit code ${XRDEXIT}, failure in xrdcp"
        exit ${XRDEXIT}
    fi
    rm ${1}
}

cd ${_CONDOR_SCRATCH_DIR}

# Make sure there's exactly one file beginning with "x509"
x509Matches=`find . -maxdepth 2 -type f -name x509*`
NMatchesMinusOne=`echo ${x509Matches} | tr -cd " \t" | wc -c`
if [ "${NMatchesMinusOne}" != "0" ]; then
    echo "ERROR: More than one file found beginning with x509"
    exit 1
fi

# Move file to proxy
mkdir proxy
mv -v x509* proxy/
x509Matches=`find proxy/ -maxdepth 2 -type f -name x509*`
export X509_USER_PROXY=${_CONDOR_SCRATCH_DIR}/${x509Matches}
echo "Set X509_USER_PROXY=${X509_USER_PROXY}"
echo "voms output:"
voms-proxy-info --all

cd ${_CONDOR_SCRATCH_DIR}

# Source CMSSW environment
echo "Sourcing CMSSW environment..."
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700
eval `scramv1 project CMSSW CMSSW_10_2_10`
cd CMSSW_10_2_10/src/
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
echo "PWD=${PWD}" && echo "Starting event selection" && ./eventSelection/bin/runEventSelection inputPathsFile=${1} selectionType=${2} disablePhotonSelection=${3} disableJetSelection=${4} lineNumberStartInclusive=${5} lineNumberEndInclusive=${6} year=${7} invertElectronVeto=${8}
EVT_SELECTION_STATUS="${?}"
if [ "${EVT_SELECTION_STATUS}" != "0" ]; then
    echo "Error in event selection: exit with code ${EVT_SELECTION_STATUS}"
    exit ${EVT_SELECTION_STATUS}
fi

PHOTONSELECTIONPREFIX=""
if [ "${3}" == "true" ]; then
    PHOTONSELECTIONPREFIX="_noPhotonSelection"
fi

JETSELECTIONPREFIX=""
if [ "${4}" == "true" ]; then
    JETSELECTIONPREFIX="_noJetSelection"
fi

ELECTRONVETOPREFIX=""
if [ "${8}" == "true" ]; then
    ELECTRONVETOPREFIX="_invertElectronVeto"
fi

echo "Copying selections..."
if [[ "${3}" == "true" ]]; then
    echo "Copying selection..."
    xrdmv_with_check selection_unified.root ${9}/${10}/selection_${2}${PHOTONSELECTIONPREFIX}${JETSELECTIONPREFIX}${ELECTRONVETOPREFIX}_${7}_unified_begin_${5}_end_${6}.root
    echo "Finished copying selection!"
else
    if [[ "${2}" == "data_singlephoton" || "${2}" =~ ^MC_GJet16_singlephoton[0-9]*$ || "${2}" =~ ^MC_GJet17_singlephoton[0-9]*$ || "${2}" =~ ^MC_GJet18_singlephoton[0-9]*$ || "${2}" =~ ^MC_QCD_singlephoton[0-9]*$ ]]; then
        echo "Copying selections..."
        xrdmv_with_check selection_control_singlemedium.root ${9}/${10}/selection_${2}${PHOTONSELECTIONPREFIX}${JETSELECTIONPREFIX}${ELECTRONVETOPREFIX}_${7}_control_singlemedium_begin_${5}_end_${6}.root
        xrdmv_with_check selection_control_singleloose.root ${9}/${10}/selection_${2}${PHOTONSELECTIONPREFIX}${JETSELECTIONPREFIX}${ELECTRONVETOPREFIX}_${7}_control_singleloose_begin_${5}_end_${6}.root
        xrdmv_with_check selection_control_singlefake.root ${9}/${10}/selection_${2}${PHOTONSELECTIONPREFIX}${JETSELECTIONPREFIX}${ELECTRONVETOPREFIX}_${7}_control_singlefake_begin_${5}_end_${6}.root
        echo "Finished copying selections!"
    elif [ "${2}" == "data_jetHT" ]; then
        echo "Copying statistics histograms..."
        xrdmv_with_check statisticsHistograms.root ${9}/${11}/statistics_${2}${PHOTONSELECTIONPREFIX}${JETSELECTIONPREFIX}${ELECTRONVETOPREFIX}_${7}_begin_${5}_end_${6}.root
        echo "Finished copying statistics!"
    else
        echo "Copying selections..."
        xrdmv_with_check selection_signal.root ${9}/${10}/selection_${2}${PHOTONSELECTIONPREFIX}${JETSELECTIONPREFIX}${ELECTRONVETOPREFIX}_${7}_signal_begin_${5}_end_${6}.root
        xrdmv_with_check selection_signal_loose.root ${9}/${10}/selection_${2}${PHOTONSELECTIONPREFIX}${JETSELECTIONPREFIX}${ELECTRONVETOPREFIX}_${7}_signal_loose_begin_${5}_end_${6}.root
        xrdmv_with_check selection_control_fakefake.root ${9}/${10}/selection_${2}${PHOTONSELECTIONPREFIX}${JETSELECTIONPREFIX}${ELECTRONVETOPREFIX}_${7}_control_fakefake_begin_${5}_end_${6}.root
        echo "Finished copying selections!"
        echo "Copying statistics histograms..."
        xrdmv_with_check statisticsHistograms.root ${9}/${11}/statistics_${2}${PHOTONSELECTIONPREFIX}${JETSELECTIONPREFIX}${ELECTRONVETOPREFIX}_${7}_begin_${5}_end_${6}.root
        echo "Finished copying statistics!"
    fi
fi

echo "All done!"
set +x
