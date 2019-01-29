#!/bin/bash

cd ${_CONDOR_SCRATCH_DIR}

echo "Output path set to ${1}, prefix set to ${2}."
echo "Running combine tool for gluino mass bin ${3}, neutralino mass bin ${4}."

echo "Starting job on: `date`" #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc630
xrdcp -s root://cmseos.fnal.gov//store/user/lpcsusystealth/combineToolCMSSW/CMSSW949cand2.tar.gz .
tar -xzf CMSSW949cand2.tar.gz && rm CMSSW949cand2.tar.gz
cd CMSSW_9_4_9_cand2/src/ && scramv1 b ProjectRename && eval `scramv1 runtime -sh` && cd ../..
echo "CMSSW version: ${CMSSW_BASE}"
echo "combine tool path:"
which combine | cat
echo "Running nominal datacard..."

combine -M AsymptoticLimits ${2}_dataCard_gluinoMassBin${3}_neutralinoMassBin${4}.txt -n "_${2}_gluinoMassBin${3}_neutralinoMassBin${4}" --rMax=10
xrdcp -f higgsCombine_${2}_gluinoMassBin${3}_neutralinoMassBin${4}.AsymptoticLimits.mH120.root ${1}/higgsCombine_${2}_gluinoMassBin${3}_neutralinoMassBin${4}.AsymptoticLimits.mH120.root && rm higgsCombine_${2}_gluinoMassBin${3}_neutralinoMassBin${4}.AsymptoticLimits.mH120.root
XRDEXIT=$?
if [[ $XRDEXIT -ne 0 ]]; then
    rm *.root
    echo "exit code $XRDEXIT, failure in xrdcp"
    exit $XRDEXIT
fi

echo "Running datacard with cross sections shifted down..."
combine -M AsymptoticLimits ${2}_dataCard_gluinoMassBin${3}_neutralinoMassBin${4}_crossSectionsDown.txt -n "_${2}_gluinoMassBin${3}_neutralinoMassBin${4}_crossSectionsDown" --rMax=10
xrdcp -f higgsCombine_${2}_gluinoMassBin${3}_neutralinoMassBin${4}_crossSectionsDown.AsymptoticLimits.mH120.root ${1}/higgsCombine_${2}_gluinoMassBin${3}_neutralinoMassBin${4}_crossSectionsDown.AsymptoticLimits.mH120.root && rm higgsCombine_${2}_gluinoMassBin${3}_neutralinoMassBin${4}_crossSectionsDown.AsymptoticLimits.mH120.root
XRDEXIT=$?
if [[ $XRDEXIT -ne 0 ]]; then
    rm *.root
    echo "exit code $XRDEXIT, failure in xrdcp"
    exit $XRDEXIT
fi

echo "Running datacard with cross sections shifted up..."
combine -M AsymptoticLimits ${2}_dataCard_gluinoMassBin${3}_neutralinoMassBin${4}_crossSectionsUp.txt -n "_${2}_gluinoMassBin${3}_neutralinoMassBin${4}_crossSectionsUp" --rMax=10
xrdcp -f higgsCombine_${2}_gluinoMassBin${3}_neutralinoMassBin${4}_crossSectionsUp.AsymptoticLimits.mH120.root ${1}/higgsCombine_${2}_gluinoMassBin${3}_neutralinoMassBin${4}_crossSectionsUp.AsymptoticLimits.mH120.root && rm higgsCombine_${2}_gluinoMassBin${3}_neutralinoMassBin${4}_crossSectionsUp.AsymptoticLimits.mH120.root
XRDEXIT=$?
if [[ $XRDEXIT -ne 0 ]]; then
    rm *.root
    echo "exit code $XRDEXIT, failure in xrdcp"
    exit $XRDEXIT
fi
cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_9_4_9_cand2

echo "combine tool ran successfully for gluino mass bin ${3}, neutralino mass bin ${4}."
