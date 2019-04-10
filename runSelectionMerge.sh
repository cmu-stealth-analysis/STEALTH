#!/bin/bash

cd /uscms/home/tmudholk/private/stealth/STEALTH
source setupEnv.sh
mkdir -p /uscmst1b_scratch/lpc1/3DayLifetime/tmudholk/merged/

YEARS=("2016" "2017")
SELECTIONTYPES=("medium" "mediumfake" "fake")
N_SELECTION_TYPES=${#SELECTIONTYPES[@]}
MAX_SELECTION_INDEX=$((N_SELECTION_TYPES-1))
SELECTION_OUTPUT_FOLDERS=("combinedSignal" "combinedControl" "combinedControl")
SELECTION_NAMES=("DoubleMedium" "OneMediumOneFake" "DoubleFake")
OPTIONAL_IDENTIFIER=""
OPTIONAL_IDENTIFIER_WITHOUT_UNDERSCORE=""
MERGETYPE_TO_RUN=""

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
        OPTIONAL_IDENTIFIER_WITHOUT_UNDERSCORE="${argValue}"
        set +x
    elif [ "${argName}" = "type" ]; then
        if [ "${argValue}" = "data" ]; then
            set -x
            MERGETYPE_TO_RUN="data"
            set +x
        elif [ "${argValue}" = "MC" ]; then
            set -x
            MERGETYPE_TO_RUN="MC"
            set +x
        else
            echo "Unrecognized selection job type: \"${argValue}\""
            exit
        fi
    else
        echo "Unrecognized argument: \"${argName}\""
        exit
    fi
done

echo "EOSPREFIX=${EOSPREFIX}"
echo "CMSSW_BASE=${CMSSW_BASE}"

EOSFolder="/store/user/lpcsusystealth/selections"
if [ "${OPTIONAL_IDENTIFIER}" != "" ]; then
    EOSFolder="/store/user/lpcsusystealth/selections/testsMerged/${OPTIONAL_IDENTIFIER_WITHOUT_UNDERSCORE}"
fi
echo "Please check output..."
set -x
eos root://cmseos.fnal.gov ls "${EOSFolder}/combinedControl"
CONTROL_DIR_EXISTS_TEST="$?"
eos root://cmseos.fnal.gov ls "${EOSFolder}/combinedSignal"
SIGNAL_DIR_EXISTS_TEST="$?"
set +x
read -p "Is this OK? (enter \"y\" for confirmation): " IS_OK
if [ "${IS_OK}" != "y" ]; then
    echo "Aborting at your request..."
    exit
fi
if [ "${CONTROL_DIR_EXISTS_TEST}" -ne 0 ]; then
    eos root://cmseos.fnal.gov mkdir -p ${EOSFolder}/combinedControl
fi
if [ "${SIGNAL_DIR_EXISTS_TEST}" -ne 0 ]; then
    eos root://cmseos.fnal.gov mkdir -p ${EOSFolder}/combinedSignal
fi

if [ "${MERGETYPE_TO_RUN}" = "data" -o "${JOBTYPE_TO_RUN}" = "" ]; then
    echo "Starting screen sessions for merging data n-tuples..."
    for SELECTION_INDEX in `seq 0 ${MAX_SELECTION_INDEX}`; do
        SELECTIONTYPE=${SELECTIONTYPES[SELECTION_INDEX]}
        SELECTION_OUTPUT_FOLDER=${SELECTION_OUTPUT_FOLDERS[SELECTION_INDEX]}
        SELECTION_NAME=${SELECTION_NAMES[SELECTION_INDEX]}
        for YEAR in ${YEARS[@]}; do
            echo "Starting merge jobs for year: ${YEAR}, selection type: ${SELECTIONTYPE}"
            set -x
            LOGFILE="mergeLogs/merger${OPTIONAL_IDENTIFIER}_data_year${YEAR}_type${SELECTIONTYPE}.log"
            screen -S "merger_data_year${YEAR}_type${SELECTIONTYPE}" -d -m bash -c "cd /uscms/home/tmudholk/private/stealth/STEALTH && source setupEnv.sh && ./mergeFiles.py --inputFilePath \"${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/${SELECTIONTYPE}${OPTIONAL_IDENTIFIER}/DoubleEG_${YEAR}_${SELECTIONTYPE}${OPTIONAL_IDENTIFIER}_begin_*.root\" --outputFilePath /uscmst1b_scratch/lpc1/3DayLifetime/tmudholk/merged/data_DoubleEG_${YEAR}_${SELECTIONTYPE}${OPTIONAL_IDENTIFIER}.root > ${LOGFILE} 2>&1 && xrdcp -f /uscmst1b_scratch/lpc1/3DayLifetime/tmudholk/merged/data_DoubleEG_${YEAR}_${SELECTIONTYPE}${OPTIONAL_IDENTIFIER}.root ${EOSPREFIX}/${EOSFolder}/${SELECTION_OUTPUT_FOLDER}/data_DoubleEG_${YEAR}_${SELECTION_NAME}${OPTIONAL_IDENTIFIER}.root >> ${LOGFILE} 2>&1 && rm -v /uscmst1b_scratch/lpc1/3DayLifetime/tmudholk/merged/data_DoubleEG_${YEAR}_${SELECTIONTYPE}${OPTIONAL_IDENTIFIER}.root >> ${LOGFILE} 2>&1"
            set +x
        done
    done
fi

if [ "${MERGETYPE_TO_RUN}" = "MC" -o "${JOBTYPE_TO_RUN}" = "" ]; then
    echo "Starting screen sessions for merging MC n-tuples..."
    for YEAR in ${YEARS[@]}; do
        for SELECTION_INDEX in `seq 0 ${MAX_SELECTION_INDEX}`; do
            SELECTIONTYPE=${SELECTIONTYPES[SELECTION_INDEX]}
            SELECTION_OUTPUT_FOLDER=${SELECTION_OUTPUT_FOLDERS[SELECTION_INDEX]}
            SELECTION_NAME=${SELECTION_NAMES[SELECTION_INDEX]}
            set -x
            LOGFILE="mergeLogs/merger${OPTIONAL_IDENTIFIER}_MC_${SELECTIONTYPE}_year${YEAR}_type.log"
            screen -S "merger_MC_${SELECTIONTYPE}_year${YEAR}" -d -m bash -c "cd /uscms/home/tmudholk/private/stealth/STEALTH && source setupEnv.sh && ./mergeFiles.py --inputFilePath \"${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/${SELECTIONTYPE}${OPTIONAL_IDENTIFIER}/MCProduction_2018_${SELECTIONTYPE}${OPTIONAL_IDENTIFIER}_optimized${YEAR}_begin_*.root\" --outputFilePath /uscmst1b_scratch/lpc1/3DayLifetime/tmudholk/merged/MC_2018Production_${SELECTION_NAME}${OPTIONAL_IDENTIFIER}_optimized${YEAR}.root > ${LOGFILE} 2>&1 && xrdcp -f /uscmst1b_scratch/lpc1/3DayLifetime/tmudholk/merged/MC_2018Production_${SELECTION_NAME}${OPTIONAL_IDENTIFIER}_optimized${YEAR}.root ${EOSPREFIX}/${EOSFolder}/${SELECTION_OUTPUT_FOLDER}/MC_2018Production_${SELECTION_NAME}${OPTIONAL_IDENTIFIER}_optimized${YEAR}.root >> ${LOGFILE} 2>&1 && rm -v /uscmst1b_scratch/lpc1/3DayLifetime/tmudholk/merged/MC_2018Production_${SELECTION_NAME}${OPTIONAL_IDENTIFIER}_optimized${YEAR}.root >> ${LOGFILE} 2>&1"
            set +x
        done
    done
fi

sleep 10
while true; do
    clear
    for outputFile in mergeLogs/merger${OPTIONAL_IDENTIFIER}_*; do
        echo "Output of ${outputFile}:"
        tail -3 ${outputFile}
        echo ""
    done
    sleep 10
done
