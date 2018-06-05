#!/bin/bash

INPUT_FROM_FILE_FLAG=""
if [ "${10}" == "inputFromFile" ]; then
    INPUT_FROM_FILE_FLAG=" inputFromFile"
fi

rm -f ${1} && touch ${1}

echo "universe = vanilla" >> ${1}
echo "Executable = submitJobs_selectEvents_Helper_Condor.sh" >> ${1}
echo "Should_Transfer_Files = YES" >> ${1}
echo "WhenToTransferOutput = ON_EXIT" >> ${1}
if [ "${10}" == "inputFromFile" ]; then
    echo "transfer_input_files = /uscms/homes/t/tmudholk/private/tmPyUtils/tmProgressBar.py, /uscms/homes/t/tmudholk/private/tmPyUtils/tmGeneralUtils.py, /uscms/homes/t/tmudholk/private/tmPyUtils/__init__.py, /uscms/homes/t/tmudholk/private/stealth/STEALTH/selectEvents.py, /uscms/homes/t/tmudholk/private/stealth/STEALTH/${2}" >> ${1}
else
    echo "transfer_input_files = /uscms/homes/t/tmudholk/private/tmPyUtils/tmProgressBar.py, /uscms/homes/t/tmudholk/private/tmPyUtils/tmGeneralUtils.py, /uscms/homes/t/tmudholk/private/tmPyUtils/__init__.py, /uscms/homes/t/tmudholk/private/stealth/STEALTH/selectEvents.py" >> ${1}
fi
echo "Output = log_${9}.stdout" >> ${1}
echo "Error = log_${9}.stderr" >> ${1}
echo "Log = log_${9}.log" >> ${1}
echo "notify_user = tmudholk@cern.ch" >> ${1}
echo "x509userproxy = \$ENV(X509_USER_PROXY)" >> ${1}
echo "Arguments = ${2} ${3} ${4} ${5} ${6} ${7} ${8}${INPUT_FROM_FILE_FLAG}" >> ${1}
echo "Queue 1" >> ${1}
