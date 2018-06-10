#!/bin/bash

echo "Running combine tool for gluino mass ${2}, neutralino mass ${3}."
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc630
cd /uscms/homes/t/tmudholk/private/stealth/cmssw/CMSSW_9_4_4/src && eval `scramv1 runtime -sh` && cd /uscms/home/tmudholk/private/stealth/STEALTH/${1}
echo "CMSSW version: "$CMSSW_BASE
echo "combine tool path:"
which combine | cat

combine -M AsymptoticLimits dataCard_gluinoMass_${2}_neutralinoMass_${3}.txt -n "_gluinoMass_${2}_neutralinoMass_${3}"
rsync --progress -av higgsCombine_gluinoMass_${2}_neutralinoMass_${3}.* /uscms/home/tmudholk/private/stealth/STEALTH/${4}/ && rm higgsCombine_gluinoMass_${2}_neutralinoMass_${3}.*
echo "Resourcing old environment..."
cd /uscms/home/tmudholk/private/stealth/STEALTH && source setupEnv.sh
echo "combine tool ran successfully for gluino mass ${2}, neutralino mass ${3}."
