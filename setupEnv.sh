set -x

BASEHOSTNAME=$(echo $HOSTNAME | sed "s|\([^\.]*\)\.cern\.ch|\1|" | sed "s|\([^\.]*\)\.fnal\.gov|\1|")

if [[ "$BASEHOSTNAME" =~ ^lxplus[0-9]{3,4}$ ]]; then
    WORKDIR="/afs/cern.ch/user/t/tmudholk/public/research/stealth/from_michael/STEALTH"
    export X509_USER_PROXY=/afs/cern.ch/user/t/tmudholk/private/x509up_u83667
    set +x && source /cvmfs/cms.cern.ch/cmsset_default.sh && set -x
    export SCRAM_ARCH=slc6_amd64_gcc630
    cd /afs/cern.ch/user/t/tmudholk/public/research/stealth/from_michael/StealthProduction/CMSSW_9_4_9_cand2/src && set +x && eval `scramv1 runtime -sh` && export PYTHONPATH=/afs/cern.ch/user/t/tmudholk/public/tmPyUtils:$PYTHONPATH && set -x
    echo "Set default CMSSW_BASE = ${CMSSW_BASE}"
    cd ${WORKDIR}
elif [[ "$BASEHOSTNAME" =~ ^cmslpc[0-9]{2}$ ]]; then
    WORKDIR="/uscms/home/tmudholk/private/stealth/STEALTH"
    export X509_USER_PROXY=/uscms/home/tmudholk/private/x509up_u50838
    set +x && source /cvmfs/cms.cern.ch/cmsset_default.sh && set -x
    export SCRAM_ARCH=slc6_amd64_gcc630
    cd /uscms/home/tmudholk/private/stealth/cmssw/CMSSW_9_4_9_cand2/src && set +x && eval `scramv1 runtime -sh` && export PYTHONPATH=/afs/cern.ch/user/t/tmudholk/public/tmPyUtils:$PYTHONPATH && set -x
    echo "Set default CMSSW_BASE = ${CMSSW_BASE}"
    cd ${WORKDIR}
else
    echo "Sourcing env only supported for lxplus and FNAL nodes at the moment."
fi

set +x
