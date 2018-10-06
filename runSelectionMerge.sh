#!/bin/bash

YEARS=("2016" "2017")
SELECTIONTYPES=("medium" "mediumfake" "fake")
N_SELECTION_TYPES=${#SELECTIONTYPES[@]}
MAX_SELECTION_INDEX=$((N_SELECTION_TYPES-1))
SELECTION_OUTPUT_FOLDERS=("combinedSignal" "combinedControl" "combinedControl")
SELECTION_NAMES=("DoubleMedium" "OneMediumOneFake" "DoubleFake")
# XROOTPREFIX="root://cmseos.fnal.gov/"

echo "EOSPREFIX=${EOSPREFIX}"
echo "CMSSW_BASE=${CMSSW_BASE}"

if [ "${1}" == "data" ]; then
    echo "Running selection merge script for data..."
    # for SELECTIONTYPE in ${SELECTIONTYPES[@]}; do
    for SELECTION_INDEX in `seq 0 ${MAX_SELECTION_INDEX}`; do
        SELECTIONTYPE=${SELECTIONTYPES[SELECTION_INDEX]}
        SELECTION_OUTPUT_FOLDER=${SELECTION_OUTPUT_FOLDERS[SELECTION_INDEX]}
        SELECTION_NAME=${SELECTION_NAMES[SELECTION_INDEX]}
        for YEAR in ${YEARS[@]}; do
            echo "Starting merge jobs for year: ${YEAR}, selection type: ${SELECTIONTYPE}"
            set -x && ./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/${SELECTIONTYPE}/DoubleEG_${YEAR}_${SELECTIONTYPE}_begin_*.root" --outputFilePath /uscms/home/tmudholk/nobackup/merged/data_DoubleEG_${YEAR}_${SELECTIONTYPE}.root && xrdcp -f /uscms/home/tmudholk/nobackup/merged/data_DoubleEG_${YEAR}_${SELECTIONTYPE}.root ${EOSPREFIX}/store/user/lpcsusystealth/selections/${SELECTION_OUTPUT_FOLDER}/data_DoubleEG_${YEAR}_${SELECTION_NAME}.root && rm -v /uscms/home/tmudholk/nobackup/merged/data_DoubleEG_${YEAR}_${SELECTIONTYPE}.root && set +x
        done
    done
elif [ "${1}" == "JECNominal" ]; then
    set -x && ./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/medium/MCProduction_2018_medium_begin_*.root" --outputFilePath /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium.root && xrdcp -f /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium.root ${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedSignal/MC_2018Production_DoubleMedium.root && rm -v /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium.root && set +x
elif [ "${1}" == "JECUp" ]; then
    set -x && ./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/medium/MCProduction_2018_medium_JECUp_begin_*.root" --outputFilePath /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium_JECUp.root && xrdcp -f /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium_JECUp.root ${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedSignal/MC_2018Production_DoubleMedium_JECUp.root && rm -v /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium_JECUp.root && set +x
elif [ "${1}" == "JECDown" ]; then
    set -x && ./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/medium/MCProduction_2018_medium_JECDown_begin_*.root" --outputFilePath /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium_JECDown.root && xrdcp -f /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium_JECDown.root ${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedSignal/MC_2018Production_DoubleMedium_JECDown.root && rm -v /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium_JECDown.root && set +x
else
    echo "ERROR:Unrecognized selection: ${1}"
fi
