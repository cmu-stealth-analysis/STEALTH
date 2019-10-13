export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source ${VO_CMS_SW_DIR}/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700
cd ${STEALTH_CMSSW_BASE}/src && eval `scramv1 runtime -sh`
cd ${STEALTH_ROOT}
