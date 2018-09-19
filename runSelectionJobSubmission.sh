#!/bin/bash

YEARS=("2016" "2017")
SELECTIONTYPES=("medium" "mediumfake" "fake")

if [ "${1}" == "data" -o "${1}" == "" ]; then
    echo "Running selection submissions for data..."
    for SELECTIONTYPE in ${SELECTIONTYPES[@]}; do
        for YEAR in ${YEARS[@]}; do
            echo "Submitting selection jobs for year: ${YEAR}, selection type: ${SELECTIONTYPE}"
            set -x && ./submitEventSelectionJobs.py --inputFilesList inputFileList_data_DoubleEG_${YEAR}.txt --isMC false --photonSelectionType ${SELECTIONTYPE} --year ${YEAR} --JECUncertainty 0 --outputFilePrefix DoubleEG_${YEAR}_${SELECTIONTYPE}_oldHLT_noInvMassCut --outputDirectory selections/DoublePhoton/${SELECTIONTYPE}_oldHLT_noInvMassCut && set +x
        done
    done
fi

# if [ "${1}" == "MC" -o "${1}" == "" ]; then
#     # For JEC uncertainty 0, submit MC selections for all selection types
#     for SELECTIONTYPE in ${SELECTIONTYPES[@]}; do
#         set -x && ./submitEventSelectionJobs.py --inputFilesList inputFileList_MC_2018Production.txt --isMC true --photonSelectionType ${SELECTIONTYPE} --year -1 --JECUncertainty 0 --outputFilePrefix MCProduction_2018_${SELECTIONTYPE} --outputDirectory selections/DoublePhoton/${SELECTIONTYPE} && set +x
#     done
#     # For JEC uncertainties, submit MC selections only for the doublle medium selection
#     OTHER_JEC_VALUES=("-1" "1")
#     OTHER_JEC_VALUE_NAMES=("JECDown" "JECUp")
#     N_JEC_VALUES=${#OTHER_JEC_VALUES[@]}
#     MAX_JEC_INDEX=$((N_JEC_VALUES-1))
#     for JEC_INDEX in `seq 0 ${MAX_JEC_INDEX}`; do
#         JEC_VALUE=${OTHER_JEC_VALUES[JEC_INDEX]}
#         JEC_NAME=${OTHER_JEC_VALUE_NAMES[JEC_INDEX]}
#         set -x && ./submitEventSelectionJobs.py --inputFilesList inputFileList_MC_2018Production.txt --isMC true --photonSelectionType medium --year -1 --JECUncertainty ${JEC_VALUE} --outputFilePrefix MCProduction_2018_medium_${JEC_NAME} --outputDirectory selections/DoublePhoton/medium && set +x
#     done
# fi

if [ "${1}" != "" ]; then
    if [ "${1}" != "data" -a "${1}" != "MC" ]; then
        echo "Unrecognized argument: ${1}"
    fi
fi
