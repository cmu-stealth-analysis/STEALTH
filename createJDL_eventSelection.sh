#!/bin/bash

rm -f ${1} && touch ${1}

echo "universe = vanilla" >> ${1}
echo "Executable = submitJobs_selectEvents_Helper_Condor.sh" >> ${1}
echo "Should_Transfer_Files = YES" >> ${1}
echo "WhenToTransferOutput = ON_EXIT" >> ${1}
echo "transfer_input_files = /uscms/homes/t/tmudholk/private/tmPyUtils/tmProgressBar.py, /uscms/homes/t/tmudholk/private/tmPyUtils/tmGeneralUtils.py, /uscms/homes/t/tmudholk/private/tmPyUtils/__init__.py, /uscms/homes/t/tmudholk/private/stealth/STEALTH/selectEvents.py, /uscms/homes/t/tmudholk/private/stealth/STEALTH/${2}" >> ${1}
echo "Output = log_selectEvents_begin_${3}_end_${4}.stdout" >> ${1}
echo "Error = log_selectEvents_begin_${3}_end_${4}.stderr" >> ${1}
echo "Log = log_selectEvents_begin_${3}_end_${4}.log" >> ${1}
echo "notify_user = tmudholk@cern.ch" >> ${1}
echo "x509userproxy = \$ENV(X509_USER_PROXY)" >> ${1}
echo "Arguments = ${2} ${3} ${4} ${5} ${6}" >> ${1}
echo "Queue 1" >> ${1}
