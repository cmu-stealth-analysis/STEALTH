#!/bin/bash

cd /uscms/home/tmudholk/private/stealth/STEALTH && source setupEnv.sh && set -x

if [ -f analysis_running.lock ]; then
    echo "Found lock file -- only one script can run at a time..."
    exit
fi

touch analysis_running.lock

function removeLockAndExit() {
    rm -f analysis_running.lock
    set +x
    exit
}

COMMON_XROOT_PREFIX="root://cmseos.fnal.gov/"
INPUTDATADIR_CONTROL="${COMMON_XROOT_PREFIX}/store/user/lpcsusystealth/selections/combinedControl"
INPUTDATADIR_SIGNAL="${COMMON_XROOT_PREFIX}/store/user/lpcsusystealth/selections/combinedSignal"
INPUTDATADIR_MCSIGNAL="${COMMON_XROOT_PREFIX}/store/user/lpcsusystealth/selections/combinedSignal"

# L_total = L_2016 + L_2017 => deltaL_total/L_total = (deltaL_2016/L_2016)*(L_2016/L_total) + (deltaL_2017/L_2017)*(L_2017/L_total)
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-17-004/index.html, the 2017 uncertainty is 2.3 percent
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-17-001/index.html, the 2016 uncertainty is 2.5 percent
INTEGLUMI="83780.0"
INTEGLUMI_FRACTIONALERROR="0.024"
DATAPATTERNCONTROL="data_DoubleEG_201*.root"
DATAPATTERNSIGNAL="data_DoubleEG_201*_DoubleMedium.root"
MCPATTERNSIGNAL="MC_2018Production_DoubleMedium.root"
MCPATTERNSIGNAL_JECUP="MC_2018Production_DoubleMedium_JECUp.root"
MCPATTERNSIGNAL_JECDOWN="MC_2018Production_DoubleMedium_JECDown.root"
MCLUMIYEARIDENTIFIER=""
YEARIDENTIFIER=""
OPTIONAL_IDENTIFIER=""
OPTIONAL_IDENTIFIER_WITHOUT_UNDERSCORE=""
USE_STANDARD_MC_SELECTION="true"
SPECIFIC_STEP_INDEX=""
STEPS_TO_RUN=(1 2 3 4 5 6 7)

for argValuePair in "$@"; do
    argName=`echo ${argValuePair} | cut --delimiter="=" --fields=1`
    argValue=`echo ${argValuePair} | cut --delimiter="=" --fields=2`
    EQUALSIGNTEST=`echo ${argValuePair} | grep "="`
    if [ "${EQUALSIGNTEST}" = "" ]; then
        echo "Argument not in format \"arg=val\": ${argValuePair}"
        removeLockAndExit
    elif [ "${argName}" = "" -o "${argValue}" = "" ]; then
        echo "Malformed argument: \"${argValuePair}\"; should be in the format arg=val"
        removeLockAndExit
    elif [ "${argName}" = "identifier" ]; then
        OPTIONAL_IDENTIFIER="_${argValue}"
        OPTIONAL_IDENTIFIER_WITHOUT_UNDERSCORE="${argValue}"
    elif [ "${argName}" = "year" ]; then
        if [ "${argValue}" = "2016" ]; then
            INTEGLUMI="37760.0"
            INTEGLUMI_FRACTIONALERROR="0.025"
            DATAPATTERNCONTROL="data_DoubleEG_2016_*.root"
            DATAPATTERNSIGNAL="data_DoubleEG_2016_DoubleMedium.root"
            MCLUMIYEARIDENTIFIER="_lumi2016"
            YEARIDENTIFIER="_2016"
        elif [ "${argValue}" = "2017" ]; then
            INTEGLUMI="46020.0"
            INTEGLUMI_FRACTIONALERROR="0.023"
            DATAPATTERNCONTROL="data_DoubleEG_2017_*.root"
            DATAPATTERNSIGNAL="data_DoubleEG_2017_DoubleMedium.root"
            MCLUMIYEARIDENTIFIER="_lumi2017"
            YEARIDENTIFIER="_2017"
        else
            echo "Unrecognized year: \"${argValue}\""
        fi
    elif [ "${argName}" = "useStandardMCSelection" ]; then
        if [ "${argValue}" = "false" ]; then
            USE_STANDARD_MC_SELECTION="false"
        else
            echo "Invalid value for \"useStandardMCSelection\": ${argValue}"
        fi
    elif [ "${argName}" = "specificStep" ]; then
        SPECIFIC_STEP_INDEX="${argValue}"
    elif [ "${argName}" = "chain" ]; then
        if [ "${argValue}" = "data" ]; then
            STEPS_TO_RUN=(1 2)
        elif [ "${argValue}" = "MC" ]; then
            STEPS_TO_RUN=(3 4)
        elif [ "${argValue}" = "combine" ]; then
            STEPS_TO_RUN=(5 6)
        elif [ "${argValue}" = "limits" ]; then
            STEPS_TO_RUN=(7)
        else
            echo "Unrecognized chain: \"${argValue}\""
            removeLockAndExit
        fi
    else
        echo "Unrecognized argument: \"${argName}\""
        removeLockAndExit
    fi
done

if [ "${OPTIONAL_IDENTIFIER}" != "" ]; then
    INPUTDATADIR_CONTROL="${COMMON_XROOT_PREFIX}/store/user/lpcsusystealth/selections/testsMerged/${OPTIONAL_IDENTIFIER_WITHOUT_UNDERSCORE}/combinedControl"
    INPUTDATADIR_SIGNAL="${COMMON_XROOT_PREFIX}/store/user/lpcsusystealth/selections/testsMerged/${OPTIONAL_IDENTIFIER_WITHOUT_UNDERSCORE}/combinedSignal"
    if [ "${YEARIDENTIFIER}" = "_2016" -o "${YEARIDENTIFIER}" = "_2017" ]; then
        DATAPATTERNCONTROL="data_DoubleEG${YEARIDENTIFIER}_*.root"
        DATAPATTERNSIGNAL="data_DoubleEG${YEARIDENTIFIER}_DoubleMedium${OPTIONAL_IDENTIFIER}.root"
    elif [ "${YEARIDENTIFIER}" = "" ]; then
        DATAPATTERNCONTROL="data_DoubleEG_201*.root"
        DATAPATTERNSIGNAL="data_DoubleEG_201*_DoubleMedium${OPTIONAL_IDENTIFIER}.root"
    else
        echo "Unknown YEARIDENTIFIER: ${YEARIDENTIFIER}"
        removeLockAndExit
    fi
    
    if [ ${USE_STANDARD_MC_SELECTION} = "false" ]; then
        INPUTDATADIR_MCSIGNAL="${COMMON_XROOT_PREFIX}/store/user/lpcsusystealth/selections/testsMerged/${OPTIONAL_IDENTIFIER_WITHOUT_UNDERSCORE}/combinedSignal"
        MCPATTERNSIGNAL="MC_2018Production_DoubleMedium${OPTIONAL_IDENTIFIER}.root"
        MCPATTERNSIGNAL_JECUP="MC_2018Production_DoubleMedium${OPTIONAL_IDENTIFIER}_JECUp$.root"
        MCPATTERNSIGNAL_JECDOWN="MC_2018Production_DoubleMedium${OPTIONAL_IDENTIFIER}_JECDown.root"
    fi
fi

function runStep() {
    case ${1} in
        1)
            ./getDataEventHistogramsAndSystematics.py --inputFilePath "${INPUTDATADIR_CONTROL}/${DATAPATTERNCONTROL}" --outputPrefix control${YEARIDENTIFIER} --allowHigherNJets
            ;;
        2)
            ./getDataEventHistogramsAndSystematics.py --inputFilePath "${INPUTDATADIR_SIGNAL}/${DATAPATTERNSIGNAL}" --outputPrefix signal${YEARIDENTIFIER} --isSignal
            ;;
        3)
            ./getMCSystematics/bin/getEventHistograms inputMCPath=${INPUTDATADIR_MCSIGNAL}/${MCPATTERNSIGNAL} inputMCPath_JECUp=${INPUTDATADIR_MCSIGNAL}/${MCPATTERNSIGNAL_JECUP} inputMCPath_JECDown=${INPUTDATADIR_MCSIGNAL}/${MCPATTERNSIGNAL_JECDOWN} outputPrefix=MC_2018${MCLUMIYEARIDENTIFIER} integratedLuminosity=${INTEGLUMI}
            ;;
        4)
            ./getMCSystematics/bin/getMCUncertainties inputPath=analysis/MCEventHistograms/MC_2018${MCLUMIYEARIDENTIFIER}_savedObjects.root outputPrefix=MC_2018${MCLUMIYEARIDENTIFIER}
            ;;
        5)
            ./createDataCards.py --outputPrefix "fullChain${YEARIDENTIFIER}" --inputFile_MCEventHistograms "analysis/MCEventHistograms/MC_2018${MCLUMIYEARIDENTIFIER}_savedObjects.root" --inputFile_MCUncertainties "analysis/MCSystematics/MC_2018${MCLUMIYEARIDENTIFIER}_MCUncertainties_savedObjects.root" --inputFile_dataSystematics "analysis/dataSystematics/signal${YEARIDENTIFIER}_dataSystematics.dat" --inputFile_dataSystematics_sTScaling "analysis/dataSystematics/control${YEARIDENTIFIER}_dataSystematics_sTScaling.dat" --inputFile_dataSystematics_eventCounters "analysis/dataSystematics/signal${YEARIDENTIFIER}_eventCounters.dat" --luminosityUncertainty ${INTEGLUMI_FRACTIONALERROR}
            ;;
        6)
            ./runCombineTool.py --dataCardsPrefix fullChain${YEARIDENTIFIER} --minGluinoMass 975.0
            ;;
        7)
            ./plotLimits.py --combineOutputPrefix fullChain${YEARIDENTIFIER} --outputSuffix fullChain${YEARIDENTIFIER} --minGluinoMass 1000.0 --maxGluinoMass 1750.0
            ;;
        *)
            echo "Unrecognized or empty step index: ${1}"
            removeLockAndExit
            ;;
    esac
}

if [ -d "analysis" ]; then
    if [ -d "analysis_original" ]; then
        echo "ERROR: directory \"analysis_original\" exists."
        removeLockAndExit
    fi
    mkdir -p analysis_original && rsync --progress -av analysis/ analysis_original/ && rm -rf analysis/
fi
if [ -d "analysis${OPTIONAL_IDENTIFIER}" ] && [ "${OPTIONAL_IDENTIFIER}" != "" ]; then
    mkdir -p analysis && rsync --progress -av analysis${OPTIONAL_IDENTIFIER}/ analysis/ && rm -rf analysis${OPTIONAL_IDENTIFIER}/
fi
mkdir -p analysis/{MCEventHistograms,MCSystematics,combineToolOutputs,dataCards,dataEventHistograms,dataSystematics,limitPlots,signalContamination}

if [ "${SPECIFIC_STEP_INDEX}" != "" ]; then
    runStep ${SPECIFIC_STEP_INDEX}
else
    for stepIndex in ${STEPS_TO_RUN[@]}; do
        runStep ${stepIndex}
    done
fi
if [ "${OPTIONAL_IDENTIFIER}" != "" ]; then
    mkdir -p analysis${OPTIONAL_IDENTIFIER} && rsync --progress -av analysis/ analysis${OPTIONAL_IDENTIFIER}/ && rm -rf analysis
fi
if [ -d "analysis_original" ]; then
    mkdir -p analysis && rsync --progress -av analysis_original/ analysis/ && rm -rf analysis_original
fi

removeLockAndExit
