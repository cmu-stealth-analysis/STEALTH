#!/bin/bash

ANALYSIS_SOURCE="/uscms/home/tmudholk/nobackup/analysisAreas/analysis"
AN_DESTINATION="/uscms/home/tmudholk/private/stealth/analysis_note_git/img"
AN_DESTINATION_TABLES="/uscms/home/tmudholk/private/stealth/analysis_note_git/tex/tables"
AN_DESTINATION_STATS="/uscms/home/tmudholk/private/stealth/analysis_note_git/tex/stats"

COPY_COMMAND="rsync --quiet --checksum --archive"

echo "Copying rho-optimization plots..."
for SIGNALTYPE in "signal" "signal_loose" "control"; do
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/dataSystematics/${SIGNALTYPE}_rhoNLL.pdf ${AN_DESTINATION}/optimizingRho/
done

echo "Copying GJet MC single and double photon fits..."
${COPY_COMMAND} ${ANALYSIS_SOURCE}/fits_doublephoton/*.pdf ${AN_DESTINATION}/fits_doublephoton/
${COPY_COMMAND} ${ANALYSIS_SOURCE}/fits_singlephoton/*.pdf ${AN_DESTINATION}/fits_singlephoton/

echo "Copying toy MC data and kernel estimate plots..."
for SIGNALTYPE in "signal" "signal_loose" "control"; do
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/dataSystematics/${SIGNALTYPE}_toyMCDataAndKernelEstimates.pdf ${AN_DESTINATION}/systematics/
done

echo "Copying shape systematics plots..."
for SIGNALTYPE in "signal" "signal_loose" "control"; do
    for STRegionIndex in `seq 2 7`; do
        ${COPY_COMMAND} ${ANALYSIS_SOURCE}/dataSystematics/${SIGNALTYPE}_shapeSystematics_STRegion${STRegionIndex}.pdf ${AN_DESTINATION}/systematics/
    done
done

echo "Copying rho systematics plots..."
for SIGNALTYPE in "signal" "signal_loose" "control"; do
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/dataSystematics/${SIGNALTYPE}_kernelPDF_rhoValues.pdf ${AN_DESTINATION}/systematics/
done

echo "Copying expected ST shapes from control selection..."
for NJETSBIN in `seq 4 6`; do
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/publicationPlots/STDistributions_control_${NJETSBIN}Jets.pdf ${AN_DESTINATION}/signalExpected/
done

echo "Copying observed and expected ST shapes and tables..."
for SIGNALTYPE in "signal" "signal_loose"; do
    for NJETSBIN in `seq 4 6`; do
        for BKGTYPE in "blinded" "preFit" "postFit"; do
	    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/publicationPlots/STDistributions_${BKGTYPE}_${SIGNALTYPE}_${NJETSBIN}Jets.pdf ${AN_DESTINATION}/signalExpected/
	done
	${COPY_COMMAND} ${ANALYSIS_SOURCE}/publicationPlots/STDistributions_postFit_${SIGNALTYPE}_${NJETSBIN}Jets_table.tex ${AN_DESTINATION_TABLES}/
    done
done

echo "Copying limits plots..."
for PRODUCTIONTYPE in "squark" "gluino"; do
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/publicationPlots/${PRODUCTIONTYPE}_*Limits.pdf ${AN_DESTINATION}/results/
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/publicationPlots/${PRODUCTIONTYPE}_injectedSignalModel_bestFitSignalStrength.pdf ${AN_DESTINATION}/results/
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/publicationPlots/${PRODUCTIONTYPE}_METUncCorrelationStudy.pdf ${AN_DESTINATION}/results/
done

echo "Copying some remaining control selection plots..."
${COPY_COMMAND} ${ANALYSIS_SOURCE}/dataEventHistograms/control_kernelPDF_normJetsBin.pdf ${AN_DESTINATION}/STShapes/
for NJETSBIN in `seq 3 6`; do
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/dataEventHistograms/control_kernelPDF_${NJETSBIN}Jets.pdf ${AN_DESTINATION}/STShapes/
done

echo "Copying signal contamination plots..."
for PRODUCTIONTYPE in "squark" "gluino"; do
    for SIGNALTYPE in "signal" "signal_loose"; do
        for STRegionIndex in `seq 1 7`; do # 2 jets: all ST regions
            ${COPY_COMMAND} ${ANALYSIS_SOURCE}/limits/signalContamination_cleaned_${SIGNALTYPE}_${PRODUCTIONTYPE}_STRegion${STRegionIndex}_2Jets.pdf ${AN_DESTINATION}/signalContamination/
        done
        for NJETSBIN in `seq 4 6`; do
            ${COPY_COMMAND} ${ANALYSIS_SOURCE}/limits/signalContamination_cleaned_${SIGNALTYPE}_${PRODUCTIONTYPE}_STRegion1_${NJETSBIN}Jets.pdf ${AN_DESTINATION}/signalContamination/  # 4-6 jets: only ST region 1
        done
	for STRegionIndex in `seq 2 7`; do
	    for NJETSBIN in `seq 4 6`; do
		for MONITORED_QUANTITY_LABEL in "fractionalSignalCorrection"; do
		    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/limits/${MONITORED_QUANTITY_LABEL}_${SIGNALTYPE}_${PRODUCTIONTYPE}_STRegion${STRegionIndex}_${NJETSBIN}Jets.pdf ${AN_DESTINATION}/signalContamination/
		done
	    done
	done
    done
done

echo "Copying MC weights..."
for PRODUCTIONTYPE in "squark" "gluino"; do
    for SIGNALTYPE in "signal" "signal_loose"; do
        for WEIGHTTYPE in "prefiring" "scale_factor"; do
            ${COPY_COMMAND} ${ANALYSIS_SOURCE}/MCEventHistograms/MC_stealth_${PRODUCTIONTYPE}_all_${SIGNALTYPE}_${WEIGHTTYPE}_weights.pdf ${AN_DESTINATION}/MCWeights/
        done
    done
done

echo "Copying missing HEM uncertainties..."
for SIGNALTYPE in "signal" "signal_loose"; do
    for UNCERTAINTY_TYPE in "JEC" "missingHEM"; do
        ${COPY_COMMAND} ${ANALYSIS_SOURCE}/MCSystematics/MC_stealth_gluino_all_${SIGNALTYPE}_${UNCERTAINTY_TYPE}UncertaintyDown_STRegion7_6Jets.pdf ${AN_DESTINATION}/missingHEM/
    done
done

echo "Copying f-statistic data from analysis logs..."
for SIGNALTYPE in "signal" "signal_loose"; do
    cat ${ANALYSIS_SOURCE}/analysisLogs/step_GJetMC_doublephoton_${SIGNALTYPE}.log | grep -A 8 "p-values for binned fit comparisons using f-statistic" | tail -n 8 > ${AN_DESTINATION_TABLES}/f_statistic_pvalues_${SIGNALTYPE}.tex
done

echo "Copying best fit values for sqrt fit from analysis logs..."
for SIGNALTYPE in "signal" "signal_loose"; do
    cat ${ANALYSIS_SOURCE}/analysisLogs/step_GJetMC_doublephoton_${SIGNALTYPE}.log | grep -A 8 "Best fit values for sqrt fit" | tail -n 8 > ${AN_DESTINATION_TABLES}/best_fit_values_sqrt_fit_${SIGNALTYPE}.tex
done

echo "Copying Asimov statistics checks..."
for IDENTIFIER in "gluinoA" "gluinoB" "gluinoC" "squarkA" "squarkB" "squarkC"; do
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/${IDENTIFIER}/asimov_signal0/r_best_fit_${IDENTIFIER}_signal_zero.tex ${AN_DESTINATION_STATS}/
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/${IDENTIFIER}/asimov_signal0/diffNuisances_${IDENTIFIER}_signal_zero.tex ${AN_DESTINATION_STATS}/
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/${IDENTIFIER}/asimov_signal0/impacts_${IDENTIFIER}_signal_zero.pdf ${AN_DESTINATION}/stats/
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/${IDENTIFIER}/asimov_signal0/covariance_fit_b_${IDENTIFIER}_signal_zero_high_res.pdf ${AN_DESTINATION}/stats/
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/${IDENTIFIER}/asimov_signal0/covariance_fit_s_${IDENTIFIER}_signal_zero_high_res.pdf ${AN_DESTINATION}/stats/
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/${IDENTIFIER}/asimov_signal0/interesting_correlations_covariance_fit_b_${IDENTIFIER}_signal_zero.tex ${AN_DESTINATION_STATS}/
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/${IDENTIFIER}/asimov_signal0/interesting_correlations_covariance_fit_s_${IDENTIFIER}_signal_zero.tex ${AN_DESTINATION_STATS}/
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/${IDENTIFIER}/asimov_signal1/r_best_fit_${IDENTIFIER}_signal_injected.tex ${AN_DESTINATION_STATS}/
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/${IDENTIFIER}/asimov_signal1/diffNuisances_${IDENTIFIER}_signal_injected.tex ${AN_DESTINATION_STATS}/
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/${IDENTIFIER}/asimov_signal1/impacts_${IDENTIFIER}_signal_injected.pdf ${AN_DESTINATION}/stats/
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/${IDENTIFIER}/asimov_signal1/covariance_fit_b_${IDENTIFIER}_signal_injected_high_res.pdf ${AN_DESTINATION}/stats/
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/${IDENTIFIER}/asimov_signal1/covariance_fit_s_${IDENTIFIER}_signal_injected_high_res.pdf ${AN_DESTINATION}/stats/
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/${IDENTIFIER}/asimov_signal1/interesting_correlations_covariance_fit_b_${IDENTIFIER}_signal_injected.tex ${AN_DESTINATION_STATS}/
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/${IDENTIFIER}/asimov_signal1/interesting_correlations_covariance_fit_s_${IDENTIFIER}_signal_injected.tex ${AN_DESTINATION_STATS}/
done

echo "Copying statistics checks on the data..."
${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/data/r_best_fit_data.tex ${AN_DESTINATION_STATS}/
${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/data/diffNuisances_data.tex ${AN_DESTINATION_STATS}/
${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/data/impacts_data.pdf ${AN_DESTINATION}/stats/
${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/data/covariance_fit_b_data_high_res.pdf ${AN_DESTINATION}/stats/
${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/data/covariance_fit_s_data_high_res.pdf ${AN_DESTINATION}/stats/
${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/data/interesting_correlations_covariance_fit_b_data.tex ${AN_DESTINATION_STATS}/
${COPY_COMMAND} ${ANALYSIS_SOURCE}/statisticsChecks/data/interesting_correlations_covariance_fit_s_data.tex ${AN_DESTINATION_STATS}/
for SIGNALTYPE in "signal" "signal_loose"; do
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/dataEventHistograms/normalization_vs_observation_Poisson_errors_${SIGNALTYPE}.tex ${AN_DESTINATION_TABLES}/
done

unset COPY_COMMAND
unset ANALYSIS_SOURCE
unset AN_DESTINATION
unset AN_DESTINATION_TABLES
unset AN_DESTINATION_STATS
