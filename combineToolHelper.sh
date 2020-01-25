#!/bin/bash

cd ${_CONDOR_SCRATCH_DIR}

OUTPUTPATH=${1}
OUTPUTPREFIX=${2}
EVENTPROGENITORMASSBIN=${3}
NEUTRALINOMASSBIN=${4}
CROSSSECTIONSFILENAME=${5}
MCTEMPLATEPATH=${6}
MCHISTOGRAMS_SIGNAL=${7}
MCHISTOGRAMS_SIGNAL_LOOSE=${8}
MCHISTOGRAMS_CONTROL=${9}
MCUNCERTAINTIES_SIGNAL=${10}
MCUNCERTAINTIES_SIGNAL_LOOSE=${11}
MCUNCERTAINTIES_CONTROL=${12}
LUMINOSITY_UNCERTAINTY=${13}
EOSANALYSISAREA=${14}
RUNUNBLINDEDSTRING=${15}
ADDLOOSESIGNALSTRING=${16}

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
for crossSectionsScale in "${crossSectionsScales[@]}"; do
    crossSectionSuffix=${crossSectionSuffixes["${crossSectionsScale}"]}
    crossSectionsScaleValue=${crossSectionsScaleValues["${crossSectionsScale}"]}

    # Step 1: create the datacard
    ./createDataCard.py --outputPrefix "${OUTPUTPREFIX}${crossSectionSuffix}" --outputDirectory "." --eventProgenitorMassBin ${EVENTPROGENITORMASSBIN} --neutralinoMassBin ${NEUTRALINOMASSBIN} --crossSectionsFile ${CROSSSECTIONSFILENAME} --crossSectionsScale ${crossSectionsScaleValue} --MCTemplatePath ${MCTEMPLATEPATH} --inputFile_MCEventHistograms_signal ${MCHISTOGRAMS_SIGNAL} --inputFile_MCEventHistograms_signal_loose ${MCHISTOGRAMS_SIGNAL_LOOSE} --inputFile_MCEventHistograms_control ${MCHISTOGRAMS_CONTROL} --inputFile_MCUncertainties_signal ${MCUNCERTAINTIES_SIGNAL} --inputFile_MCUncertainties_signal_loose ${MCUNCERTAINTIES_SIGNAL_LOOSE} --inputFile_MCUncertainties_control ${MCUNCERTAINTIES_CONTROL} --inputFile_dataSystematics_signal "data/signal_dataSystematics.dat" --inputFile_dataSystematics_signal_loose "data/signal_loose_dataSystematics.dat" --inputFile_dataSystematics_control "data/control_dataSystematics.dat" --inputFile_dataSystematics_expectedEventCounters_signal "data/signal_eventCounters.dat" --inputFile_dataSystematics_expectedEventCounters_signal_loose "data/signal_loose_eventCounters.dat" --inputFile_dataSystematics_expectedEventCounters_control "data/control_eventCounters.dat" --inputFile_dataSystematics_observedEventCounters_signal "data/signal_observedEventCounters.dat" --inputFile_dataSystematics_observedEventCounters_signal_loose "data/signal_loose_observedEventCounters.dat" --inputFile_dataSystematics_observedEventCounters_control "data/control_observedEventCounters.dat" --luminosityUncertainty ${LUMINOSITY_UNCERTAINTY} --regionsToUse ${REGIONSTOUSE}${UNBLINDED_RUN_FLAG}
    xrdcp_with_check "${OUTPUTPREFIX}${crossSectionSuffix}_dataCard_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt" "${EOSANALYSISAREA}/dataCards/combinedFit/${OUTPUTPREFIX}${crossSectionSuffix}_dataCard_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt"
    echo "After step 1, list of files:"
    ls -alh

    # Step 2: Run the data card
    IS_CONVERGENT="false"
    RUNNING_RMAX="20.0"
    while [ "${IS_CONVERGENT}" = "false" ]; do
        echo "Trying --rMax=${RUNNING_RMAX}..."
	combine -M AsymptoticLimits -d "${OUTPUTPREFIX}${crossSectionSuffix}_dataCard_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt" -n "_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}" -v 1 -V --rMax="${RUNNING_RMAX}"
	./checkLimitsConvergence.py --inputROOTFile "higgsCombine_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.AsymptoticLimits.mH120.root" > tmp_bestFitCheck.txt 2>&1
        IS_CONVERGENT=`cat tmp_bestFitCheck.txt | tr -d '\n'` # tr -d '\n' deletes all newlines
        rm -v -r -f tmp_bestFitCheck.txt
        RUNNING_RMAX_NEW=`python -c "print(${RUNNING_RMAX}/10.0)"`
        RUNNING_RMAX="${RUNNING_RMAX_NEW}"
    done
    echo "After step 2, list of files:"
    ls -alh

    # Step 3: Copy important output file to EOS
    xrdcp_with_check "higgsCombine_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.AsymptoticLimits.mH120.root" "${OUTPUTPATH}/higgsCombine_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.AsymptoticLimits.mH120.root" && rm "higgsCombine_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.AsymptoticLimits.mH120.root"
    echo "After step 3, list of files:"
    ls -alh
done

cd ${_CONDOR_SCRATCH_DIR}
echo "combine tool ran successfully for eventProgenitor mass bin ${EVENTPROGENITORMASSBIN}, neutralino mass bin ${NEUTRALINOMASSBIN}."
echo "Removing everything else..."
rm -v -r -f *_dataCard_*.txt
rm -v -r -f data
rm -v -r -f STRegionBoundaries.dat
rm -v -r -f *_dataCard_*.root
cleanup

set +x
