#!/bin/bash

rm -f ${1} && touch ${1}

echo "universe = vanilla" >> ${1}
echo "Executable = submitJobs_skimKinematic_Data_Helper_Condor.sh" >> ${1}
echo "Should_Transfer_Files = YES" >> ${1}
echo "WhenToTransferOutput = ON_EXIT" >> ${1}
echo "transfer_input_files = /uscms/homes/t/tmudholk/private/tmPyUtils/tmProgressBar.py, /uscms/homes/t/tmudholk/private/tmPyUtils/__init__.py, /uscms/homes/t/tmudholk/private/stealth/STEALTH/skimKinematic_Data.py" >> ${1}
echo "Output = log_skimKinematic_begin_${3}_end_${4}.stdout" >> ${1}
echo "Error = log_skimKinematic_begin_${3}_end_${4}.stderr" >> ${1}
echo "Log = log_skimKinematic_begin_${3}_end_${4}.log" >> ${1}
echo "notify_user = tmudholk@cern.ch" >> ${1}
echo "x509userproxy = $ENV(X509_USER_PROXY)" >> ${1}
echo "Arguments = ${2} ${3} ${4} ${5}" >> ${1}
echo "Queue 1" >> ${1}
