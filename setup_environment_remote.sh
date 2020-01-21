#!/bin/bash

function xrdcp_with_check {
    if [ "${#}" != 2  ]; then
        echo "ERROR: number of arguments passed to \"${FUNCNAME}\": ${#}"
        exit 1
    fi
    xrdcp --verbose --force --path --streams 15 ${1} ${2} 2>&1
    XRDEXIT=${?}
    if [[ ${XRDEXIT} -ne 0 ]]; then
        echo "exit code ${XRDEXIT}, failure in xrdcp"
        exit ${XRDEXIT}
    fi
}

function xrdmv_with_check {
    if [ "${#}" != 2  ]; then
        echo "ERROR: number of arguments passed to \"${FUNCNAME}\": ${#}"
        exit 1
    fi
    xrdcp_with_check "${1}" "${2}"
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

echo "Starting job on: `date`" #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node

# Source CMSSW environment
echo "Sourcing CMSSW environment..."
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700
xrdcp --verbose --force --path --streams 15 root://cmseos.fnal.gov//store/user/lpcsusystealth/combineToolCMSSW/CMSSW10210.tar.gz .
tar -xzf CMSSW10210.tar.gz && rm CMSSW10210.tar.gz
cd CMSSW_10_2_10/src/ && scramv1 b ProjectRename && eval `scramv1 runtime -sh` && cd ../..
echo "CMSSW version: ${CMSSW_BASE}"

echo "Extracting tmUtils tarball..."
./extract_tmUtilsTarball.sh && rm -f tmUtils.tar.gz

echo "Setting env variables..."
export PYTHONPATH=${_CONDOR_SCRATCH_DIR}/tmPyUtils:${PYTHONPATH}
export TMPYUTILS=${_CONDOR_SCRATCH_DIR}/tmPyUtils/
export TMCPPUTILS=${_CONDOR_SCRATCH_DIR}/tmCPPUtils/
