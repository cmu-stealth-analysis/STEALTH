WORKDIR="/afs/cern.ch/user/t/tmudholk/public/research/stealth/from_michael/STEALTH"
export X509_USER_PROXY=/afs/cern.ch/user/t/tmudholk/private/x509up_u83667
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/t/tmudholk/public/research/ecaldqm/cmssw/CMSSW_9_2_13/src && eval `scramv1 runtime -sh`
export PYTHONPATH=/afs/cern.ch/user/t/tmudholk/public/tmPyUtils:$PYTHONPATH
cd ${WORKDIR}
