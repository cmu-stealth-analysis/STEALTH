#!/bin/bash

cd ${_CONDOR_SCRATCH_DIR}

echo "Output path set to ${1}, prefix set to ${2}, initial rmax set to ${5}."
echo "Running combine tool for eventProgenitor mass bin ${3}, neutralino mass bin ${4}."

echo "Starting job on: `date`" #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700
xrdcp --verbose --force --path --streams 15 root://cmseos.fnal.gov//store/user/lpcsusystealth/combineToolCMSSW/CMSSW10210.tar.gz .
tar -xzf CMSSW10210.tar.gz && rm CMSSW10210.tar.gz
cd CMSSW_10_2_10/src/ && scramv1 b ProjectRename && eval `scramv1 runtime -sh` && cd ../..
echo "CMSSW version: ${CMSSW_BASE}"
echo "combine tool path:"
which combine | cat
echo "Running combine tool for all three input datacards..."

set -x
for crossSectionsSuffix in "" "_crossSectionsDown" "_crossSectionsUp"; do
    echo "Trying crossSectionsSuffix=\"${crossSectionsSuffix}\""
    RUNNING_RMAX="${5}"
    IS_CONVERGENT="false"
    while [ "${IS_CONVERGENT}" = "false" ]; do
        echo "No convergent result found yet. Trying --rMax=${RUNNING_RMAX}..."
        combine -M AsymptoticLimits "${2}_dataCard_eventProgenitorMassBin${3}_neutralinoMassBin${4}${crossSectionsSuffix}.txt" -n "_${2}_eventProgenitorMassBin${3}_neutralinoMassBin${4}${crossSectionsSuffix}" --rMax="${RUNNING_RMAX}"
        ./checkLimitsConvergence.py --inputROOTFile "higgsCombine_${2}_eventProgenitorMassBin${3}_neutralinoMassBin${4}${crossSectionsSuffix}.AsymptoticLimits.mH120.root" > tmp_limitsCheck.txt
        IS_CONVERGENT=`cat tmp_limitsCheck.txt | tr -d '\n'` # tr -d '\n' deletes all newlines
        rm tmp_limitsCheck.txt
        RUNNING_RMAX_NEW=`python -c "print(${RUNNING_RMAX}/2.0)"`
        RUNNING_RMAX="${RUNNING_RMAX_NEW}"
    done
    xrdcp --verbose --force --path --streams 15 "higgsCombine_${2}_eventProgenitorMassBin${3}_neutralinoMassBin${4}${crossSectionsSuffix}.AsymptoticLimits.mH120.root" "${1}/higgsCombine_${2}_eventProgenitorMassBin${3}_neutralinoMassBin${4}${crossSectionsSuffix}.AsymptoticLimits.mH120.root" && rm "higgsCombine_${2}_eventProgenitorMassBin${3}_neutralinoMassBin${4}${crossSectionsSuffix}.AsymptoticLimits.mH120.root"
    XRDEXIT=$?
    if [[ $XRDEXIT -ne 0 ]]; then
        rm *.root
        echo "exit code $XRDEXIT, failure in xrdcp"
        exit $XRDEXIT
    fi
done

cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_10_2_10

echo "combine tool ran successfully for eventProgenitor mass bin ${3}, neutralino mass bin ${4}."
set +x
