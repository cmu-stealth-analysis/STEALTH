#!/bin/bash

INPUT_FROM_FILE_FLAG=""
if [ "${9}" == "inputFromFile" ]; then
    INPUT_FROM_FILE_FLAG=" --inputFromFile"
fi

cd ${_CONDOR_SCRATCH_DIR}
mkdir tmPyUtils
mv tmProgressBar.py tmPyUtils/.
mv tmGeneralUtils.py tmPyUtils/.
mv __init__.py tmPyUtils/.

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc630
eval `scramv1 project CMSSW CMSSW_9_2_13`
cd CMSSW_9_2_13/src/
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
echo "CMSSW: "$CMSSW_BASE
cd ${_CONDOR_SCRATCH_DIR}
export PYTHONPATH=${_CONDOR_SCRATCH_DIR}/tmPyUtils:${PYTHONPATH}

echo "PWD=${PWD}" && \
    echo "Starting event selection" && \
    echo "Python version: " && python --version && \
    ./skim_tree_original.py${INPUT_FROM_FILE_FLAG} --inputFilePath=${1} --outputFilePath=${2} --counterStartInclusive=${3} --counterEndInclusive=${4} --photonSelectionType=${5} --year=${6} --JECUncertainty=${7}

echo "*******************************************"
OUTDIR=root://cmseos.fnal.gov//store/user/lpcsusystealth/${8}
echo "xrdcp output for condor"
for FILE in *.root
do
  echo "xrdcp -f ${FILE} ${OUTDIR}/${FILE}"
  xrdcp -f ${FILE} ${OUTDIR}/${FILE} 2>&1
  XRDEXIT=$?
  if [[ $XRDEXIT -ne 0 ]]; then
    rm *.root
    echo "exit code $XRDEXIT, failure in xrdcp"
    exit $XRDEXIT
  fi
  rm ${FILE}
done
