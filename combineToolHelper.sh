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

REGIONSTOUSE="signal"
if [ "${ADDLOOSESIGNALSTRING}" == "true" ]; then
    REGIONSTOUSE="signal,signal_loose"
else
    if [ ! "${ADDLOOSESIGNALSTRING}" == "false" ]; then
        echo "ERROR: ADDLOOSESIGNALSTRING can only take values \"true\" or \"false\". Currently, ADDLOOSESIGNALSTRING: ${ADDLOOSESIGNALSTRING}"
        exit 1
    fi
fi

# UNBLINDED_RUN_FLAG=" --usePoissonForAsimov"
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
RUNNING_RMAX="20.0"
MIN_RMAX="0.001"
for crossSectionsScale in "${crossSectionsScales[@]}"; do
    crossSectionSuffix=${crossSectionSuffixes["${crossSectionsScale}"]}
    crossSectionsScaleValue=${crossSectionsScaleValues["${crossSectionsScale}"]}

    # Step 1: create the datacard
    ./createDataCard.py --outputPrefix "${OUTPUTPREFIX}${crossSectionSuffix}" --outputDirectory "." --eventProgenitorMassBin ${EVENTPROGENITORMASSBIN} --neutralinoMassBin ${NEUTRALINOMASSBIN} --crossSectionsFile ${CROSSSECTIONSFILENAME} --crossSectionsScale ${crossSectionsScaleValue} --MCTemplatePath ${MCTEMPLATEPATH} --inputFile_MCEventHistograms_signal ${MCHISTOGRAMS_SIGNAL} --inputFile_MCEventHistograms_signal_loose ${MCHISTOGRAMS_SIGNAL_LOOSE} --inputFile_MCEventHistograms_control ${MCHISTOGRAMS_CONTROL} --inputFile_MCUncertainties_signal ${MCUNCERTAINTIES_SIGNAL} --inputFile_MCUncertainties_signal_loose ${MCUNCERTAINTIES_SIGNAL_LOOSE} --inputFile_MCUncertainties_control ${MCUNCERTAINTIES_CONTROL} --inputFile_dataSystematics_signal "data/signal_dataSystematics.dat" --inputFile_dataSystematics_signal_loose "data/signal_loose_dataSystematics.dat" --inputFile_dataSystematics_control "data/control_dataSystematics.dat" --inputFile_dataSystematics_scaling_signal "data/signal_GJet17_STComparisons_scalingSystematics.dat" --inputFile_dataSystematics_scaling_signal_loose "data/signal_loose_GJet17_STComparisons_scalingSystematics.dat" --inputFile_dataSystematics_scaling_control "data/control_GJet17_STComparisons_scalingSystematics.dat" --inputFile_dataSystematics_scalingQuality "data/scalingQualitySystematics_combined.dat" --inputFile_dataSystematics_MCShapeAdjustment_signal "data/adjustments_all_MC_GJet_signal.dat" --inputFile_dataSystematics_MCShapeAdjustment_signal_loose "data/adjustments_all_MC_GJet_signal_loose.dat" --inputFile_dataSystematics_MCShapeAdjustment_control "data/adjustments_all_MC_GJet_control.dat" --inputFile_dataSystematics_dataMCRatioAdjustment_signal "data/ratio_adjustment_2017_data_singlemedium.dat" --inputFile_dataSystematics_dataMCRatioAdjustment_signal_loose "data/ratio_adjustment_2017_data_singleloose.dat" --inputFile_dataSystematics_expectedEventCounters_signal "data/signal_eventCounters.dat" --inputFile_dataSystematics_expectedEventCounters_signal_loose "data/signal_loose_eventCounters.dat" --inputFile_dataSystematics_expectedEventCounters_control "data/control_eventCounters.dat" --inputFile_dataSystematics_observedEventCounters_signal "data/signal_observedEventCounters.dat" --inputFile_dataSystematics_observedEventCounters_signal_loose "data/signal_loose_observedEventCounters.dat" --inputFile_dataSystematics_observedEventCounters_control "data/control_observedEventCounters.dat" --luminosityUncertainty ${LUMINOSITY_UNCERTAINTY} --regionsToUse ${REGIONSTOUSE}${UNBLINDED_RUN_FLAG}
    xrdcp_with_check "${OUTPUTPREFIX}${crossSectionSuffix}_dataCard_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt" "${EOSANALYSISAREA}/dataCards/combinedFit/${OUTPUTPREFIX}${crossSectionSuffix}_dataCard_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt"
    echo "Step 1 done."

    # Step 2: Run the data card
    IS_CONVERGENT="false"
    RUNNING_RMAX="20.0"
    while [ "${IS_CONVERGENT}" != "true" ]; do
        echo "Trying --rMax=${RUNNING_RMAX}..."
        if (( $(echo "${RUNNING_RMAX} < ${MIN_RMAX}" | bc -l) )); then
            echo "Hit lower limit on RUNNING_RMAX"
            break
        fi
        combine -M AsymptoticLimits -d "${OUTPUTPREFIX}${crossSectionSuffix}_dataCard_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt" -n "_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}" -v 1 -V --expectSignal 0 --rMax="${RUNNING_RMAX}"
        ./checkLimitsConvergence.py --inputROOTFile "higgsCombine_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.AsymptoticLimits.mH120.root" --checkObservedLimit > tmp_bestFitCheck.txt 2>&1
        IS_CONVERGENT=`cat tmp_bestFitCheck.txt | tr -d '\n'` # tr -d '\n' deletes all newlines
        rm -v -r -f tmp_bestFitCheck.txt
        RUNNING_RMAX_NEW=`python -c "print(${RUNNING_RMAX}/10.0)"`
        RUNNING_RMAX="${RUNNING_RMAX_NEW}"
    done
    if [[ "${IS_CONVERGENT}" == "false" ]]; then
        echo "Combine output does not converge at crossSectionScale: ${crossSectionScale}, event progenitor mass bin: ${EVENTPROGENITORMASSBIN}, neutralino mass bin: ${NEUTRALINOMASSBIN}"
        continue
    fi
    echo "Step 2 done."

    # Step 3: Copy important output file to EOS
    xrdcp_with_check "higgsCombine_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.AsymptoticLimits.mH120.root" "${OUTPUTPATH}/higgsCombine_${OUTPUTPREFIX}${crossSectionSuffix}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.AsymptoticLimits.mH120.root"
    echo "Step 3 done."
done

# Step 4: get best-fit values for parameters of interest (obsolete now, but was used to get best-fit values of the scaling parameters used in an earlier version of the analysis, preserving here for reference)
# Run multiDimFit on nominal datacard and transfer multi dim fit output to EOS
# First get the expected upper limit from the nominal cross-section output and write it to file "tmp_rmax.txt"
python -c "import commonFunctions; commonFunctions.write_ten_times_expected_upper_limit_from_combine_output_to_file(combineOutputFilePath=\"higgsCombine_${OUTPUTPREFIX}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.AsymptoticLimits.mH120.root\", outputFilePath=\"tmp_rmax.txt\")"
RMAX_TO_USE=`cat tmp_rmax.txt | tr -d '\n'` # tr -d '\n' deletes all newlines
rm -v -r -f tmp_rmax.txt
if [ ${RMAX_TO_USE} = "unavailable" ]; then
    echo "rmax not available for MultiDimFit, not calculating best fit values..."
else
    # Run combine tool
    combine -M MultiDimFit --saveFitResult -d "${OUTPUTPREFIX}_dataCard_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt" -n "_${OUTPUTPREFIX}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}" --expectSignal 0 -v 1 -V --rMax=${RMAX_TO_USE}
    # Copy multidimfit to EOS
    xrdcp_with_check "multidimfit_${OUTPUTPREFIX}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.root" "${OUTPUTPATH}/multidimfit_${OUTPUTPREFIX}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.root"
    echo "Step 4 done."
fi

# Step 5: Create datacard with added signal and copy it over
./createDataCard.py --outputPrefix "${OUTPUTPREFIX}" --outputDirectory "." --eventProgenitorMassBin ${EVENTPROGENITORMASSBIN} --neutralinoMassBin ${NEUTRALINOMASSBIN} --crossSectionsFile ${CROSSSECTIONSFILENAME} --crossSectionsScale 0 --MCTemplatePath ${MCTEMPLATEPATH} --inputFile_MCEventHistograms_signal ${MCHISTOGRAMS_SIGNAL} --inputFile_MCEventHistograms_signal_loose ${MCHISTOGRAMS_SIGNAL_LOOSE} --inputFile_MCEventHistograms_control ${MCHISTOGRAMS_CONTROL} --inputFile_MCUncertainties_signal ${MCUNCERTAINTIES_SIGNAL} --inputFile_MCUncertainties_signal_loose ${MCUNCERTAINTIES_SIGNAL_LOOSE} --inputFile_MCUncertainties_control ${MCUNCERTAINTIES_CONTROL} --inputFile_dataSystematics_signal "data/signal_dataSystematics.dat" --inputFile_dataSystematics_signal_loose "data/signal_loose_dataSystematics.dat" --inputFile_dataSystematics_control "data/control_dataSystematics.dat" --inputFile_dataSystematics_scaling_signal "data/signal_GJet17_STComparisons_scalingSystematics.dat" --inputFile_dataSystematics_scaling_signal_loose "data/signal_loose_GJet17_STComparisons_scalingSystematics.dat" --inputFile_dataSystematics_scaling_control "data/control_GJet17_STComparisons_scalingSystematics.dat" --inputFile_dataSystematics_scalingQuality "data/scalingQualitySystematics_combined.dat" --inputFile_dataSystematics_MCShapeAdjustment_signal "data/adjustments_all_MC_GJet_signal.dat" --inputFile_dataSystematics_MCShapeAdjustment_signal_loose "data/adjustments_all_MC_GJet_signal_loose.dat" --inputFile_dataSystematics_MCShapeAdjustment_control "data/adjustments_all_MC_GJet_control.dat" --inputFile_dataSystematics_dataMCRatioAdjustment_signal "data/ratio_adjustment_2017_data_singlemedium.dat" --inputFile_dataSystematics_dataMCRatioAdjustment_signal_loose "data/ratio_adjustment_2017_data_singleloose.dat" --inputFile_dataSystematics_expectedEventCounters_signal "data/signal_eventCounters.dat" --inputFile_dataSystematics_expectedEventCounters_signal_loose "data/signal_loose_eventCounters.dat" --inputFile_dataSystematics_expectedEventCounters_control "data/control_eventCounters.dat" --inputFile_dataSystematics_observedEventCounters_signal "data/signal_observedEventCounters.dat" --inputFile_dataSystematics_observedEventCounters_signal_loose "data/signal_loose_observedEventCounters.dat" --inputFile_dataSystematics_observedEventCounters_control "data/control_observedEventCounters.dat" --luminosityUncertainty ${LUMINOSITY_UNCERTAINTY} --regionsToUse ${REGIONSTOUSE}${UNBLINDED_RUN_FLAG} --addSignalToBackground
xrdcp_with_check "WITH_ADDED_SIGNAL_${OUTPUTPREFIX}_dataCard_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt" "${EOSANALYSISAREA}/dataCards/combinedFit/WITH_ADDED_SIGNAL_${OUTPUTPREFIX}_dataCard_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt"
echo "Step 5 done."

# Step 6: Run MultiDimFit on datacard with added signal
# Run combine tool
combine -M MultiDimFit --saveFitResult -d "WITH_ADDED_SIGNAL_${OUTPUTPREFIX}_dataCard_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.txt" -n "_WITH_ADDED_SIGNAL_${OUTPUTPREFIX}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}" --expectSignal 0 -v 1 -V # no rMax needed because we're already adding a signal with strength 1
# Copy multidimfit to EOS
xrdcp_with_check "multidimfit_WITH_ADDED_SIGNAL_${OUTPUTPREFIX}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.root" "${OUTPUTPATH}/multidimfit_WITH_ADDED_SIGNAL_${OUTPUTPREFIX}_eventProgenitorMassBin${EVENTPROGENITORMASSBIN}_neutralinoMassBin${NEUTRALINOMASSBIN}.root"
echo "Step 6 done."

cd ${_CONDOR_SCRATCH_DIR}
echo "combine tool ran successfully for eventProgenitor mass bin ${EVENTPROGENITORMASSBIN}, neutralino mass bin ${NEUTRALINOMASSBIN}."
echo "Removing everything else..."
rm -r -f *_dataCard_*.txt
rm -r -f data
rm -r -f STRegionBoundaries.dat
rm -r -f *_dataCard_*.root
rm -r -f higgsCombine_*.root
rm -r -f multidimfit_*.root
cleanup

set +x
