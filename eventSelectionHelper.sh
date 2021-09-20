#!/bin/bash

INPUTPATHSFILE=${1}
SELECTIONTYPE=${2}
DISABLEPHOTONSEL=${3}
DISABLEJETSEL=${4}
LSTART=${5}
LEND=${6}
YEAR=${7}
INV_ELE_VETO=${8}
MCWEIGHT=${9}
EOS_PREFIX=${10}
SEL_OUTPUT_FOLDER_PATH=${11}
STATS_OUTPUT_FOLDER_PATH=${12}

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
echo "PWD=${PWD}" && echo "Starting event selection" && ./eventSelection/bin/runEventSelection inputPathsFile=${INPUTPATHSFILE} selectionType=${SELECTIONTYPE} disablePhotonSelection=${DISABLEPHOTONSEL} disableJetSelection=${DISABLEJETSEL} lineNumberStartInclusive=${LSTART} lineNumberEndInclusive=${LEND} year=${YEAR} invertElectronVeto=${INV_ELE_VETO} MCBackgroundWeight=${MCWEIGHT}
EVT_SELECTION_STATUS="${?}"
if [ "${EVT_SELECTION_STATUS}" != "0" ]; then
    echo "Error in event selection: exit with code ${EVT_SELECTION_STATUS}"
    exit ${EVT_SELECTION_STATUS}
fi

PHOTONSELECTIONPREFIX=""
if [ "${DISABLEPHOTONSEL}" == "true" ]; then
    PHOTONSELECTIONPREFIX="_noPhotonSelection"
fi

JETSELECTIONPREFIX=""
if [ "${DISABLEJETSEL}" == "true" ]; then
    JETSELECTIONPREFIX="_noJetSelection"
fi

ELECTRONVETOPREFIX=""
if [ "${INV_ELE_VETO}" == "true" ]; then
    ELECTRONVETOPREFIX="_invertElectronVeto"
fi

OVERALL_PREFIX="${PHOTONSELECTIONPREFIX}${JETSELECTIONPREFIX}${ELECTRONVETOPREFIX}"

echo "Copying selections..."
if [[ "${DISABLEPHOTONSEL}" == "true" ]]; then
    echo "Copying selection..."
    xrdmv_with_check selection_unified.root ${EOS_PREFIX}/${SEL_OUTPUT_FOLDER_PATH}/selection_${SELECTIONTYPE}${OVERALL_PREFIX}_${YEAR}_unified_begin_${LSTART}_end_${LEND}.root
    echo "Finished copying selection!"
else
    if [[ "${SELECTIONTYPE}" == "data_singlephoton" || "${SELECTIONTYPE}" =~ ^MC_GJet16_singlephoton[0-9]*$ || "${SELECTIONTYPE}" =~ ^MC_GJet17_singlephoton[0-9]*$ || "${SELECTIONTYPE}" =~ ^MC_GJet18_singlephoton[0-9]*$ || "${SELECTIONTYPE}" =~ ^MC_QCD16_singlephoton[0-9]*$ || "${SELECTIONTYPE}" =~ ^MC_QCD17_singlephoton[0-9]*$ || "${SELECTIONTYPE}" =~ ^MC_QCD18_singlephoton[0-9]*$ ]]; then
        echo "Copying selections..."
        xrdmv_with_check selection_control_singlemedium.root ${EOS_PREFIX}/${SEL_OUTPUT_FOLDER_PATH}/selection_${SELECTIONTYPE}${OVERALL_PREFIX}_${YEAR}_control_singlemedium_begin_${LSTART}_end_${LEND}.root
        xrdmv_with_check selection_control_singleloose.root ${EOS_PREFIX}/${SEL_OUTPUT_FOLDER_PATH}/selection_${SELECTIONTYPE}${OVERALL_PREFIX}_${YEAR}_control_singleloose_begin_${LSTART}_end_${LEND}.root
        xrdmv_with_check selection_control_singlefake.root ${EOS_PREFIX}/${SEL_OUTPUT_FOLDER_PATH}/selection_${SELECTIONTYPE}${OVERALL_PREFIX}_${YEAR}_control_singlefake_begin_${LSTART}_end_${LEND}.root
        echo "Finished copying selections!"
    elif [ "${SELECTIONTYPE}" == "data_jetHT" ]; then
        echo "Copying statistics histograms..."
        xrdmv_with_check statisticsHistograms.root ${EOS_PREFIX}/${STATS_OUTPUT_FOLDER_PATH}/statistics_${SELECTIONTYPE}${OVERALL_PREFIX}_${YEAR}_begin_${LSTART}_end_${LEND}.root
        echo "Finished copying statistics!"
    else
        echo "Copying selections..."
        xrdmv_with_check selection_signal.root ${EOS_PREFIX}/${SEL_OUTPUT_FOLDER_PATH}/selection_${SELECTIONTYPE}${OVERALL_PREFIX}_${YEAR}_signal_begin_${LSTART}_end_${LEND}.root
        xrdmv_with_check selection_signal_loose.root ${EOS_PREFIX}/${SEL_OUTPUT_FOLDER_PATH}/selection_${SELECTIONTYPE}${OVERALL_PREFIX}_${YEAR}_signal_loose_begin_${LSTART}_end_${LEND}.root
        xrdmv_with_check selection_control_fakefake.root ${EOS_PREFIX}/${SEL_OUTPUT_FOLDER_PATH}/selection_${SELECTIONTYPE}${OVERALL_PREFIX}_${YEAR}_control_fakefake_begin_${LSTART}_end_${LEND}.root
        echo "Finished copying selections!"
        echo "Copying statistics histograms..."
        xrdmv_with_check statisticsHistograms.root ${EOS_PREFIX}/${STATS_OUTPUT_FOLDER_PATH}/statistics_${SELECTIONTYPE}${OVERALL_PREFIX}_${YEAR}_begin_${LSTART}_end_${LEND}.root
        echo "Finished copying statistics!"
    fi
fi

echo "All done!"
set +x
