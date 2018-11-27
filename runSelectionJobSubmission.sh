#!/bin/bash

cd /uscms/home/tmudholk/private/stealth/STEALTH
source setupEnv.sh

read -p "Move old logs to archives? (enter \"y\" for confirmation): " IS_OK
if [ "${IS_OK}" = "y" ]; then
    mkdir -p ~/nobackup/archives/logs && mv condor_working_directory/* ~/nobackup/archives/logs/.
    mkdir -p ~/nobackup/archives/logs/mergeLogs && mv mergeLogs/* ~/nobackup/archives/logs/mergeLogs/.
fi

set -x
YEARS=("2016" "2017")
N_YEARS=${#YEARS[@]}
INPUTFILELISTS_DATA=("inputFileList_data_DoubleEG_2016_ntuplized80X.txt" "inputFileList_data_DoubleEG_2017_ntuplized949cand2.txt")
INPUTFILELISTS_MC=("inputFileList_MC_2018Production_ntuplized80X.txt" "inputFileList_MC_2018Production_ntuplized949cand2.txt")
echo "N_YEARS=${N_YEARS}"
MAX_YEAR_INDEX=$((N_YEARS-1))
echo "MAX_YEAR_INDEX=${MAX_YEAR_INDEX}"
SELECTIONTYPES=("medium" "mediumfake" "fake")
OPTIONAL_IDENTIFIER=""
JOBTYPE_TO_RUN=""
USECACHE="false"
DRYRUNFLAG=""
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
    elif [ "${argName}" = "dry_run" ]; then
        if [ "${argValue}" = "true" ]; then
            DRYRUNFLAG=" --isDryRun"
        else
            echo "Unrecognized value for \"dry_run\": \"${argValue}\""
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
        # eos root://cmseos.fnal.gov rm -r ${OUTPUT_DIRECTORY}
        eos root://cmseos.fnal.gov mkdir -p ${OUTPUT_DIRECTORY}
    else
        echo "Aborting at your request..."
        exit
    fi
done

if [ "${JOBTYPE_TO_RUN}" = "data" -o "${JOBTYPE_TO_RUN}" = "" ]; then
    echo "Running selection submissions for data..."
    for YEAR_INDEX in `seq 0 ${MAX_YEAR_INDEX}`; do
        YEAR=${YEARS[YEAR_INDEX]}
        INPUTFILELIST_DATA=${INPUTFILELISTS_DATA[YEAR_INDEX]}
        echo "Fetching number of events in input files list for ${YEAR} data..."
        if [ "${USECACHE}" = "false" ]; then
            echo "Removing cached number of events if it exists..."
            rm -f cached_nEvents_${INPUTFILELIST_DATA}
        fi
        if [ ! -f cached_nEvents_${INPUTFILELIST_DATA} ]; then
            echo "Caching number of events..."
            ./getNEvts.py --inputFilesList "${INPUTFILELIST_DATA}" | tee cached_nEvents_${INPUTFILELIST_DATA}
        fi
    done
    for SELECTIONTYPE in ${SELECTIONTYPES[@]}; do
        for YEAR_INDEX in `seq 0 ${MAX_YEAR_INDEX}`; do
            YEAR=${YEARS[YEAR_INDEX]}
            INPUTFILELIST_DATA=${INPUTFILELISTS_DATA[YEAR_INDEX]}
            TOTAL_N_EVENTS=`cat cached_nEvents_${INPUTFILELIST_DATA} | tr -d '\n'` # tr -d '\n' deletes all newlines
            echo "Submitting selection jobs for year: ${YEAR}, selection type: ${SELECTIONTYPE}"
            set -x && ./submitEventSelectionJobs.py${DRYRUNFLAG} --inputFilesList ${INPUTFILELIST_DATA} --nEventsInInputFilesList ${TOTAL_N_EVENTS} --isMC false --photonSelectionType ${SELECTIONTYPE} --year ${YEAR} --JECUncertainty 0 --outputFilePrefix DoubleEG_${YEAR}_${SELECTIONTYPE}${OPTIONAL_IDENTIFIER} --outputDirectory selections/DoublePhoton/${SELECTIONTYPE}${OPTIONAL_IDENTIFIER} && set +x
        done
    done
fi

if [ "${JOBTYPE_TO_RUN}" = "MC" -o "${JOBTYPE_TO_RUN}" = "" ]; then
    echo "Running selection submissions for MC..."
    for YEAR_INDEX in `seq 0 ${MAX_YEAR_INDEX}`; do
        YEAR=${YEARS[YEAR_INDEX]}
        INPUTFILELIST_MC=${INPUTFILELISTS_MC[YEAR_INDEX]}
        echo "Fetching number of events in input files list for ${YEAR}-optimized MC..."
        if [ "${USECACHE}" = "false" ]; then
            echo "Removing cached number of events if it exists..."
            rm -f cached_nEvents_${INPUTFILELIST_MC}
        fi
        if [ ! -f cached_nEvents_${INPUTFILELIST_MC} ]; then
            echo "Caching number of events..."
            ./getNEvts.py --inputFilesList "${INPUTFILELIST_MC}" | tee cached_nEvents_${INPUTFILELIST_MC}
        fi
    done

    for YEAR_INDEX in `seq 0 ${MAX_YEAR_INDEX}`; do
        YEAR=${YEARS[YEAR_INDEX]}
        INPUTFILELIST_MC=${INPUTFILELISTS_MC[YEAR_INDEX]}
        TOTAL_MC_EVENTS=`cat cached_nEvents_${INPUTFILELIST_MC} | tr -d '\n'` # tr -d '\n' deletes all newlines
        # For JEC uncertainty 0, submit MC selections for all selection types
        for SELECTIONTYPE in ${SELECTIONTYPES[@]}; do
            set -x && ./submitEventSelectionJobs.py${DRYRUNFLAG} --inputFilesList ${INPUTFILELIST_MC} --nEventsInInputFilesList ${TOTAL_MC_EVENTS} --isMC true --photonSelectionType ${SELECTIONTYPE} --year ${YEAR} --JECUncertainty 0 --outputFilePrefix MCProduction_2018_${SELECTIONTYPE}${OPTIONAL_IDENTIFIER}_optimized${YEAR} --outputDirectory selections/DoublePhoton/${SELECTIONTYPE}${OPTIONAL_IDENTIFIER} && set +x
        done
        # For JEC uncertainties, submit MC selections only for the double medium selection
        OTHER_JEC_VALUES=("-1" "1")
        OTHER_JEC_VALUE_NAMES=("JECDown" "JECUp")
        N_JEC_VALUES=${#OTHER_JEC_VALUES[@]}
        MAX_JEC_INDEX=$((N_JEC_VALUES-1))
        for JEC_INDEX in `seq 0 ${MAX_JEC_INDEX}`; do
            JEC_VALUE=${OTHER_JEC_VALUES[JEC_INDEX]}
            JEC_NAME=${OTHER_JEC_VALUE_NAMES[JEC_INDEX]}
            set -x && ./submitEventSelectionJobs.py${DRYRUNFLAG} --inputFilesList ${INPUTFILELIST_MC} --nEventsInInputFilesList ${TOTAL_MC_EVENTS} --isMC true --photonSelectionType medium --year ${YEAR} --JECUncertainty ${JEC_VALUE} --outputFilePrefix MCProduction_2018_medium_${JEC_NAME}${OPTIONAL_IDENTIFIER}_optimized${YEAR} --outputDirectory selections/DoublePhoton/medium${OPTIONAL_IDENTIFIER} && set +x
        done
    done
fi
