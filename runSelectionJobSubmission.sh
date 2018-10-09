#!/bin/bash

cd /uscms/home/tmudholk/private/stealth/STEALTH
source setupEnv.sh
set -x
YEARS=("2016" "2017")
SELECTIONTYPES=("medium" "mediumfake" "fake")
OPTIONAL_IDENTIFIER=""
JOBTYPE_TO_RUN=""
USECACHE="false"
set +x

for argValuePair in "$@"; do
    argName=`echo ${argValuePair} | cut --delimiter="=" --fields=1`
    argValue=`echo ${argValuePair} | cut --delimiter="=" --fields=2`
    EQUALSIGNTEST=`echo ${argValuePair} | grep "="`
    if [ "${EQUALSIGNTEST}" = "" ]; then
        echo "Argument not in format \"arg=val\": ${argValuePair}"
        exit
    elif [ "${argName}" = "" -o "${argValue}" = "" ]; then
        echo "Malformed argument: \"${argValuePair}\"; should be in the format arg=val"
        exit
    elif [ "${argName}" = "identifier" ]; then
        set -x
        OPTIONAL_IDENTIFIER="_${argValue}"
        set +x
    elif [ "${argName}" = "type" ]; then
        if [ "${argValue}" = "data" ]; then
            set -x
            JOBTYPE_TO_RUN="data"
            set +x
        elif [ "${argValue}" = "MC" ]; then
            set -x
            JOBTYPE_TO_RUN="MC"
            set +x
        else
            echo "Unrecognized selection job type: \"${argValue}\""
            exit
        fi
    elif [ "${argName}" = "use_cache" ]; then
        if [ "${argValue}" = "true" ]; then
            USECACHE="true"
        else
            echo "Unrecognized value for \"use_cache\": \"${argValue}\""
        fi
    else
        echo "Unrecognized argument: \"${argName}\""
        exit
    fi
done

for SELECTIONTYPE in ${SELECTIONTYPES[@]}; do
    OUTPUT_DIRECTORY="/store/user/lpcsusystealth/selections/DoublePhoton/${SELECTIONTYPE}${OPTIONAL_IDENTIFIER}"
    echo "Please check output..."
    set -x && eos root://cmseos.fnal.gov ls "${OUTPUT_DIRECTORY}" && set +x
    read -p "Is this OK? (enter \"y\" for confirmation): " IS_OK
    if [ "${IS_OK}" = "y" ]; then
        eos root://cmseos.fnal.gov rm -r ${OUTPUT_DIRECTORY}
        eos root://cmseos.fnal.gov mkdir -p ${OUTPUT_DIRECTORY}
    else
        echo "Aborting at your request..."
        exit
    fi
done

if [ "${JOBTYPE_TO_RUN}" = "data" -o "${JOBTYPE_TO_RUN}" = "" ]; then
    echo "Running selection submissions for data..."
    for YEAR in ${YEARS[@]}; do
        echo "Fetching number of events in input files list for ${YEAR} data..."
        if [ "${USECACHE}" = "false" ]; then
            echo "Removing cached number of events if it exists..."
            rm -f cached_nEvents_inputFileList_data_DoubleEG_${YEAR}.txt
        fi
        if [ ! -f cached_nEvents_inputFileList_data_DoubleEG_${YEAR}.txt ]; then
            echo "Caching number of events..."
            ./getNEvts.py --inputFilesList "inputFileList_data_DoubleEG_${YEAR}.txt" | tee cached_nEvents_inputFileList_data_DoubleEG_${YEAR}.txt
        fi
    done
    for SELECTIONTYPE in ${SELECTIONTYPES[@]}; do
        for YEAR in ${YEARS[@]}; do
            TOTAL_N_EVENTS=`cat cached_nEvents_inputFileList_data_DoubleEG_${YEAR}.txt | tr -d '\n'` # tr -d '\n' deletes all newlines
            echo "Submitting selection jobs for year: ${YEAR}, selection type: ${SELECTIONTYPE}"
            set -x && ./submitEventSelectionJobs.py --inputFilesList inputFileList_data_DoubleEG_${YEAR}.txt --nEventsInInputFilesList ${TOTAL_N_EVENTS} --isMC false --photonSelectionType ${SELECTIONTYPE} --year ${YEAR} --JECUncertainty 0 --outputFilePrefix DoubleEG_${YEAR}_${SELECTIONTYPE}${OPTIONAL_IDENTIFIER} --outputDirectory selections/DoublePhoton/${SELECTIONTYPE}${OPTIONAL_IDENTIFIER} && set +x
        done
    done
fi

if [ "${JOBTYPE_TO_RUN}" = "MC" -o "${JOBTYPE_TO_RUN}" = "" ]; then
    echo "Running selection submissions for MC..."
    echo "Fetching number of events in input files list for MC..."
    if [ "${USECACHE}" = "false" ]; then
        echo "Removing cached number of events if it exists..."
        rm -f cached_nEvents_inputFileList_MC_2018Production.txt
    fi
    if [ ! -f cached_nEvents_inputFileList_MC_2018Production.txt ]; then
        echo "Caching number of events..."
        ./getNEvts.py --inputFilesList "inputFileList_MC_2018Production.txt" | tee cached_nEvents_inputFileList_MC_2018Production.txt
    fi
    set -x && TOTAL_MC_EVENTS=`cat cached_nEvents_inputFileList_MC_2018Production.txt | tr -d '\n'` && set +x
    # For JEC uncertainty 0, submit MC selections for all selection types
    for SELECTIONTYPE in ${SELECTIONTYPES[@]}; do
        set -x && ./submitEventSelectionJobs.py --inputFilesList inputFileList_MC_2018Production.txt --nEventsInInputFilesList ${TOTAL_MC_EVENTS} --isMC true --photonSelectionType ${SELECTIONTYPE} --year -1 --JECUncertainty 0 --outputFilePrefix MCProduction_2018_${SELECTIONTYPE}${OPTIONAL_IDENTIFIER} --outputDirectory selections/DoublePhoton/${SELECTIONTYPE}${OPTIONAL_IDENTIFIER} && set +x
    done
    # For JEC uncertainties, submit MC selections only for the double medium selection
    OTHER_JEC_VALUES=("-1" "1")
    OTHER_JEC_VALUE_NAMES=("JECDown" "JECUp")
    N_JEC_VALUES=${#OTHER_JEC_VALUES[@]}
    MAX_JEC_INDEX=$((N_JEC_VALUES-1))
    for JEC_INDEX in `seq 0 ${MAX_JEC_INDEX}`; do
        JEC_VALUE=${OTHER_JEC_VALUES[JEC_INDEX]}
        JEC_NAME=${OTHER_JEC_VALUE_NAMES[JEC_INDEX]}
        set -x && ./submitEventSelectionJobs.py --inputFilesList inputFileList_MC_2018Production.txt --nEventsInInputFilesList ${TOTAL_MC_EVENTS} --isMC true --photonSelectionType medium --year -1 --JECUncertainty ${JEC_VALUE} --outputFilePrefix MCProduction_2018_medium_${JEC_NAME}${OPTIONAL_IDENTIFIER} --outputDirectory selections/DoublePhoton/medium${OPTIONAL_IDENTIFIER} && set +x
    done
fi
