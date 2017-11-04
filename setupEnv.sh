BASEHOSTNAME=$(echo $HOSTNAME | sed "s|\([^\.]*\)\.cern\.ch|\1|" | sed "s|\([^\.]*\)\.fnal\.gov|\1|")

if [[ "$BASEHOSTNAME" =~ ^lxplus[0-9]{3,4}$ ]]; then
    WORKDIR="/afs/cern.ch/user/t/tmudholk/public/research/stealth/from_michael/STEALTH"
    export X509_USER_PROXY=/afs/cern.ch/user/t/tmudholk/private/x509up_u83667
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    cd /afs/cern.ch/user/t/tmudholk/public/research/ecaldqm/cmssw/CMSSW_9_2_13/src && eval `scramv1 runtime -sh`
    export PYTHONPATH=/afs/cern.ch/user/t/tmudholk/public/tmPyUtils:$PYTHONPATH
    cd ${WORKDIR}
elif [[ "$BASEHOSTNAME" =~ ^cmslpc[0-9]{2}$ ]]; then
    WORKDIR="/uscms/home/tmudholk/private/stealth/STEALTH"
    export X509_USER_PROXY=/uscms/home/tmudholk/private/x509up_u50838
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    cd /uscms/home/tmudholk/private/stealth/cmssw/CMSSW_9_2_13/src && eval `scramv1 runtime -sh`
    export PYTHONPATH=/afs/cern.ch/user/t/tmudholk/public/tmPyUtils:$PYTHONPATH
    cd ${WORKDIR}
fi
