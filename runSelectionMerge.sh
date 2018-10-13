#!/bin/bash

cd /uscms/home/tmudholk/private/stealth/STEALTH
source setupEnv.sh
mkdir -p /uscms/home/tmudholk/nobackup/merged/

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
        elif [ "${argValue}" = "JECNominal" ]; then
            set -x
            MERGETYPE_TO_RUN="JECNominal"
            set +x
        elif [ "${argValue}" = "JECUp" ]; then
            set -x
            MERGETYPE_TO_RUN="JECUp"
            set +x
        elif [ "${argValue}" = "JECDown" ]; then
            set -x
            MERGETYPE_TO_RUN="JECDown"
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

if [ "${MERGETYPE_TO_RUN}" = "data" ]; then
    echo "Running selection merge script for data..."
    for SELECTION_INDEX in `seq 0 ${MAX_SELECTION_INDEX}`; do
        SELECTIONTYPE=${SELECTIONTYPES[SELECTION_INDEX]}
        SELECTION_OUTPUT_FOLDER=${SELECTION_OUTPUT_FOLDERS[SELECTION_INDEX]}
        SELECTION_NAME=${SELECTION_NAMES[SELECTION_INDEX]}
        for YEAR in ${YEARS[@]}; do
            echo "Starting merge jobs for year: ${YEAR}, selection type: ${SELECTIONTYPE}"
            set -x && ./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/${SELECTIONTYPE}${OPTIONAL_IDENTIFIER}/DoubleEG_${YEAR}_${SELECTIONTYPE}${OPTIONAL_IDENTIFIER}_begin_*.root" --outputFilePath /uscms/home/tmudholk/nobackup/merged/data_DoubleEG_${YEAR}_${SELECTIONTYPE}${OPTIONAL_IDENTIFIER}.root && xrdcp -f /uscms/home/tmudholk/nobackup/merged/data_DoubleEG_${YEAR}_${SELECTIONTYPE}${OPTIONAL_IDENTIFIER}.root ${EOSPREFIX}/${EOSFolder}/${SELECTION_OUTPUT_FOLDER}/data_DoubleEG_${YEAR}_${SELECTION_NAME}${OPTIONAL_IDENTIFIER}.root && rm -v /uscms/home/tmudholk/nobackup/merged/data_DoubleEG_${YEAR}_${SELECTIONTYPE}${OPTIONAL_IDENTIFIER}.root && set +x
        done
    done
elif [ "${MERGETYPE_TO_RUN}" = "JECNominal" ]; then
    set -x && ./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/medium${OPTIONAL_IDENTIFIER}/MCProduction_2018_medium${OPTIONAL_IDENTIFIER}_begin_*.root" --outputFilePath /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium${OPTIONAL_IDENTIFIER}.root && xrdcp -f /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium${OPTIONAL_IDENTIFIER}.root ${EOSPREFIX}/${EOSFolder}/combinedSignal/MC_2018Production_DoubleMedium${OPTIONAL_IDENTIFIER}.root && rm -v /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium${OPTIONAL_IDENTIFIER}.root && set +x
elif [ "${MERGETYPE_TO_RUN}" = "JECUp" ]; then
    set -x && ./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/medium${OPTIONAL_IDENTIFIER}/MCProduction_2018_medium_JECUp${OPTIONAL_IDENTIFIER}_begin_*.root" --outputFilePath /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium${OPTIONAL_IDENTIFIER}_JECUp.root && xrdcp -f /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium${OPTIONAL_IDENTIFIER}_JECUp.root ${EOSPREFIX}/${EOSFolder}/combinedSignal/MC_2018Production_DoubleMedium${OPTIONAL_IDENTIFIER}_JECUp.root && rm -v /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium${OPTIONAL_IDENTIFIER}_JECUp.root && set +x
elif [ "${MERGETYPE_TO_RUN}" = "JECDown" ]; then
    set -x && ./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/medium${OPTIONAL_IDENTIFIER}/MCProduction_2018_medium_JECDown${OPTIONAL_IDENTIFIER}_begin_*.root" --outputFilePath /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium${OPTIONAL_IDENTIFIER}_JECDown.root && xrdcp -f /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium${OPTIONAL_IDENTIFIER}_JECDown.root ${EOSPREFIX}/${EOSFolder}/combinedSignal/MC_2018Production_DoubleMedium${OPTIONAL_IDENTIFIER}_JECDown.root && rm -v /uscms/home/tmudholk/nobackup/merged/MC_2018Production_DoubleMedium${OPTIONAL_IDENTIFIER}_JECDown.root && set +x
else
    echo "ERROR:Unrecognized selection: ${MERGETYPE_TO_RUN}"
fi
