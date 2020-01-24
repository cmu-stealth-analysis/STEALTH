#!/bin/bash

cd ${_CONDOR_SCRATCH_DIR}

OUTPUTPATH=${1}
OUTPUTPREFIX=${2}
EVENTPROGENITORMASSBIN=${3}
NEUTRALINOMASSBIN=${4}
INITIALRMAX=${5}
CROSSSECTIONSFILENAME=${6}
MCTEMPLATEPATH=${7}
MCHISTOGRAMS_SIGNAL=${8}
MCHISTOGRAMS_SIGNAL_LOOSE=${9}
MCHISTOGRAMS_CONTROL=${10}
MCUNCERTAINTIES_SIGNAL=${11}
MCUNCERTAINTIES_SIGNAL_LOOSE=${12}
MCUNCERTAINTIES_CONTROL=${13}
LUMINOSITY_UNCERTAINTY=${14}
EOSANALYSISAREA=${15}
RUNUNBLINDEDSTRING=${16}
ADDLOOSESIGNALSTRING=${17}

source setup_environment_remote.sh

echo "combine tool path:"
which combine | cat

set -x
echo "Starting to run combine chain..."

mv -v MCTemplateReader.py commonFunctions.py ${TMPYUTILS}/
mkdir "data"
mv -v *.dat data/
mv -v data/STRegionBoundaries.dat ./

REGIONSTOUSE="signal,control"
if [ "${ADDLOOSESIGNALSTRING}" == "true" ]; then
    REGIONSTOUSE="signal,signal_loose,control"
else
    if [ ! "${ADDLOOSESIGNALSTRING}" == "false" ]; then
	echo "ERROR: ADDLOOSESIGNALSTRING can only take values \"true\" or \"false\". Currently, ADDLOOSESIGNALSTRING: ${ADDLOOSESIGNALSTRING}"
	exit 1
    fi
fi

UNBLINDED_RUN_FLAG=""
if [ "${RUNUNBLINDEDSTRING}" == "true" ]; then
    UNBLINDED_RUN_FLAG=" --runUnblinded"
else
    if [ ! "${RUNUNBLINDEDSTRING}" == "false" ]; then
	echo "ERROR: RUNUNBLINDEDSTRING can only take values \"true\" or \"false\". Currently, RUNUNBLINDEDSTRING: ${RUNUNBLINDEDSTRING}"
	exit 1
    fi
fi

crossSectionsScales=( "nominal" "down" "up" )
declare -A crossSectionSuffixes=( ["nominal"]="" ["down"]="_crossSectionsDown" ["up"]="_crossSectionsUp" )
declare -A crossSectionsScaleValues=( ["nominal"]="0" ["down"]="-1" ["up"]="1" )
# declare -A bestFitSignalStrengths
for crossSectionsScale in "${crossSectionsScales[@]}"; do
    crossSectionSuffix=${crossSectionSuffixes["${crossSectionsScale}"]}
    crossSectionsScaleValue=${crossSectionsScaleValues["${crossSectionsScale}"]}
    # Step 1: create the datacard without the ST scaling systematics
    ./createDataCard.py --outputPrefix "forFit_${OUTPUTPREFIX}${crossSectionSuffix}" --outputDirectory "." --eventProgenitorMassBin ${EVENTPROGENITORMASSBIN} --neutralinoMassBin ${NEUTRALINOMASSBIN} --crossSectionsFile ${CROSSSECTIONSFILENAME} --crossSectionsScale ${crossSectionsScaleValue} --MCTemplatePath ${MCTEMPLATEPATH} --inputFile_MCEventHistograms_signal ${MCHISTOGRAMS_SIGNAL} --inputFile_MCEventHistograms_signal_loose ${MCHISTOGRAMS_SIGNAL_LOOSE} --inputFile_MCEventHistograms_control ${MCHISTOGRAMS_CONTROL} --inputFile_MCUncertainties_signal ${MCUNCERTAINTIES_SIGNAL} --inputFile_MCUncertainties_signal_loose ${MCUNCERTAINTIES_SIGNAL_LOOSE} --inputFile_MCUncertainties_control ${MCUNCERTAINTIES_CONTROL} --inputFile_dataSystematics_signal "data/signal_dataSystematics.dat" --inputFile_dataSystematics_signal_loose "data/signal_loose_dataSystematics.dat" --inputFile_dataSystematics_control "data/control_dataSystematics.dat" --inputFile_dataSystematics_expectedEventCounters_signal "data/signal_eventCounters.dat" --inputFile_dataSystematics_expectedEventCounters_signal_loose "data/signal_loose_eventCounters.dat" --inputFile_dataSystematics_expectedEventCounters_control "data/control_eventCounters.dat" --inputFile_dataSystematics_observedEventCounters_signal "data/signal_observedEventCounters.dat" --inputFile_dataSystematics_observedEventCounters_signal_loose "data/signal_loose_observedEventCounters.dat" --inputFile_dataSystematics_observedEventCounters_control "data/control_observedEventCounters.dat" --luminosityUncertainty ${LUMINOSITY_UNCERTAINTY} --regionsToUse ${REGIONSTOUSE}${UNBLINDED_RUN_FLAG}
    xrdcp_with_check "forFit_${OUTPUTPREFIX}${crossSectionSuffix}_dataCard_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt" "${EOSANALYSISAREA}/dataCards/noSTScalingUncertainties/forFit_${OUTPUTPREFIX}${crossSectionSuffix}_dataCard_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt"
    # list files
    echo "After step 1, list of files:"
    ls -alh

    # Step 2: get the best fit value
    echo "Starting to get best fit signal strength for crossSectionSuffix=\"${crossSectionSuffix}\"".
    RUNNING_RMAX="${INITIALRMAX}"
    IS_CONVERGENT="false"
    while [ "${IS_CONVERGENT}" = "false" ]; do
        echo "No convergent result found yet. Trying --rMax=${RUNNING_RMAX}..."
	combine -M MultiDimFit "forFit_${OUTPUTPREFIX}${crossSectionSuffix}_dataCard_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt" -n "_forFit_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}" --rMax="${RUNNING_RMAX}" --robustFit=1
	./checkBestFitConvergence.py --inputROOTFile "higgsCombine_forFit_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.MultiDimFit.mH120.root" --rMaxPassed "${RUNNING_RMAX}" > tmp_bestFitCheck.txt 2>&1
        IS_CONVERGENT=`cat tmp_bestFitCheck.txt | tr -d '\n'` # tr -d '\n' deletes all newlines
        rm -v -r -f tmp_bestFitCheck.txt
        RUNNING_RMAX_NEW=`python -c "print(${RUNNING_RMAX}/1.9)"`
        RUNNING_RMAX="${RUNNING_RMAX_NEW}"
    done
    xrdcp_with_check "higgsCombine_forFit_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.MultiDimFit.mH120.root" "${EOSANALYSISAREA}/multiDimOutputs/higgsCombine_forFit_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.MultiDimFit.mH120.root"
    ./readBestFitFromMultiDimOutput.py --inputROOTFile "higgsCombine_forFit_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.MultiDimFit.mH120.root" > multiDimOutput.txt 2>&1
    BEST_FIT=`cat multiDimOutput.txt | tr -d '\n'` # tr -d '\n' deletes all newlines
    xrdcp_with_check "multiDimOutput.txt" "${EOSANALYSISAREA}/bestFitSignalStrengths/multiDimOutput_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt"
    # list files
    echo "After step 2, list of files:"
    ls -alh

    # Step 3: Get the ST scaling systematic from the best fit signal value
    echo "Fetching ST scaling systematics with crossSectionSuffix=\"${crossSectionSuffix}\""
    ./getSTScalingSystematics.py --inputFile_nEvtsExpected "data/control_eventCounters.dat" --inputFile_shapeSystematics "data/control_dataSystematics.dat" --inputFile_nSignalEvents ${MCHISTOGRAMS_CONTROL} --eventProgenitorMassBin ${EVENTPROGENITORMASSBIN} --neutralinoMassBin ${NEUTRALINOMASSBIN} --bestFitSignalStrength ${BEST_FIT}
    xrdcp_with_check "dataSystematics_scaling.dat" "${EOSANALYSISAREA}/dataSystematics_scaling/dataSystematics_scaling_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt"
    # list files
    echo "After step 3, list of files:"
    ls -alh
    mv -v dataSystematics_scaling.dat data/

    # Step 4: Create data card with ST scaling uncertainties calculated in the previous step
    ./createDataCard.py --outputPrefix "${OUTPUTPREFIX}${crossSectionSuffix}" --outputDirectory "." --eventProgenitorMassBin ${EVENTPROGENITORMASSBIN} --neutralinoMassBin ${NEUTRALINOMASSBIN} --crossSectionsFile ${CROSSSECTIONSFILENAME} --crossSectionsScale ${crossSectionsScaleValue} --MCTemplatePath ${MCTEMPLATEPATH} --inputFile_MCEventHistograms_signal ${MCHISTOGRAMS_SIGNAL} --inputFile_MCEventHistograms_signal_loose ${MCHISTOGRAMS_SIGNAL_LOOSE} --inputFile_MCEventHistograms_control ${MCHISTOGRAMS_CONTROL} --inputFile_MCUncertainties_signal ${MCUNCERTAINTIES_SIGNAL} --inputFile_MCUncertainties_signal_loose ${MCUNCERTAINTIES_SIGNAL_LOOSE} --inputFile_MCUncertainties_control ${MCUNCERTAINTIES_CONTROL} --inputFile_dataSystematics_signal "data/signal_dataSystematics.dat" --inputFile_dataSystematics_signal_loose "data/signal_loose_dataSystematics.dat" --inputFile_dataSystematics_control "data/control_dataSystematics.dat" --inputFile_dataSystematics_sTScaling "data/dataSystematics_scaling.dat" --useSTScalingSystematics --inputFile_dataSystematics_expectedEventCounters_signal "data/signal_eventCounters.dat" --inputFile_dataSystematics_expectedEventCounters_signal_loose "data/signal_loose_eventCounters.dat" --inputFile_dataSystematics_expectedEventCounters_control "data/control_eventCounters.dat" --inputFile_dataSystematics_observedEventCounters_signal "data/signal_observedEventCounters.dat" --inputFile_dataSystematics_observedEventCounters_signal_loose "data/signal_loose_observedEventCounters.dat" --inputFile_dataSystematics_observedEventCounters_control "data/control_observedEventCounters.dat" --luminosityUncertainty ${LUMINOSITY_UNCERTAINTY} --regionsToUse ${REGIONSTOUSE}${UNBLINDED_RUN_FLAG}
    xrdcp_with_check "${OUTPUTPREFIX}${crossSectionSuffix}_dataCard_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt" "${EOSANALYSISAREA}/dataCards/withSTScalingUncertainties/${OUTPUTPREFIX}${crossSectionSuffix}_dataCard_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt"
    # list files
    echo "After step 4, list of files:"
    ls -alh

    # Step 5: Run the data card
    RUNNING_RMAX="${INITIALRMAX}"
    IS_CONVERGENT="false"
    while [ "${IS_CONVERGENT}" = "false" ]; do
        echo "No convergent result found yet. Trying --rMax=${RUNNING_RMAX}..."
        combine -M AsymptoticLimits "${OUTPUTPREFIX}${crossSectionSuffix}_dataCard_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt" -n "_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}" --rMax="${RUNNING_RMAX}"
        ./checkLimitsConvergence.py --inputROOTFile "higgsCombine_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.AsymptoticLimits.mH120.root" > tmp_limitsCheck.txt
        IS_CONVERGENT=`cat tmp_limitsCheck.txt | tr -d '\n'` # tr -d '\n' deletes all newlines
        rm -v -r -f tmp_limitsCheck.txt
        RUNNING_RMAX_NEW=`python -c "print(${RUNNING_RMAX}/1.9)"`
        RUNNING_RMAX="${RUNNING_RMAX_NEW}"
    done
    # list files
    echo "After step 5, list of files:"
    ls -alh

    # Step 6: Copy important output file to EOS
    xrdcp_with_check "higgsCombine_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.AsymptoticLimits.mH120.root" "${OUTPUTPATH}/higgsCombine_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.AsymptoticLimits.mH120.root" && rm "higgsCombine_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.AsymptoticLimits.mH120.root"
    # list files
    echo "After step 6, list of files:"
    ls -alh
done

cd ${_CONDOR_SCRATCH_DIR}
echo "combine tool ran successfully for eventProgenitor mass bin ${EVENTPROGENITORMASSBIN}, neutralino mass bin ${NEUTRALINOMASSBIN}."
echo "Removing everything else..."
rm -v -r -f *_dataCard_*.txt
rm -v -r -f data
rm -v -r -f STRegionBoundaries.dat
rm -v -r -f *_dataCard_*.root
rm -r -f multiDimOutput.txt

cleanup

set +x
