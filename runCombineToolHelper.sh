#!/bin/bash

echo "Running combine tool for gluino mass bin ${3}, neutralino mass bin ${4}."
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc630
cd /uscms/home/tmudholk/private/stealth/cmssw/CMSSW_9_4_9_cand2/src && eval `scramv1 runtime -sh` && cd /uscms/home/tmudholk/private/stealth/STEALTH/${1}
echo "CMSSW version: "$CMSSW_BASE
echo "combine tool path:"
which combine | cat

combine -M AsymptoticLimits ${2}_dataCard_gluinoMassBin${3}_neutralinoMassBin${4}.txt -n "_${2}_gluinoMassBin${3}_neutralinoMassBin${4}"
rsync --progress -av higgsCombine_${2}_gluinoMassBin${3}_neutralinoMassBin${4}.* /uscms/home/tmudholk/private/stealth/STEALTH/${5}/ && rm higgsCombine_${2}_gluinoMassBin${3}_neutralinoMassBin${4}.*
cd /uscms/home/tmudholk/private/stealth/STEALTH
echo "combine tool ran successfully for gluino mass bin ${3}, neutralino mass bin ${4}."
