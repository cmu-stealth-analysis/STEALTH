#!/bin/bash

ANALYSIS_SOURCE="/uscms/home/tmudholk/nobackup/analysisAreas/analysis"
AN_DESTINATION="/uscms/home/tmudholk/private/stealth/analysis_note_git/img"

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

echo "Copying expected ST shapes..."
for SIGNALTYPE in "signal" "signal_loose" "control"; do
    for NJETSBIN in `seq 4 6`; do
        for PRODUCTIONTYPE in "squark" "gluino"; do
            ${COPY_COMMAND} ${ANALYSIS_SOURCE}/publicationPlots/STDistributions_${PRODUCTIONTYPE}_${SIGNALTYPE}_${NJETSBIN}Jets.pdf ${AN_DESTINATION}/signalExpected/
        done
    done
done

echo "Copying limits plots..."
for PRODUCTIONTYPE in "squark" "gluino"; do
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/publicationPlots/${PRODUCTIONTYPE}_expectedLimits.pdf ${AN_DESTINATION}/results/
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
            ${COPY_COMMAND} ${ANALYSIS_SOURCE}/signalContamination/MC_stealth_${PRODUCTIONTYPE}_all_${SIGNALTYPE}_signalContamination_STRegion${STRegionIndex}_2Jets.pdf ${AN_DESTINATION}/signalContamination/
        done
        for NJETSBIN in `seq 4 6`; do
            ${COPY_COMMAND} ${ANALYSIS_SOURCE}/signalContamination/MC_stealth_${PRODUCTIONTYPE}_all_${SIGNALTYPE}_signalContamination_STRegion1_${NJETSBIN}Jets.pdf ${AN_DESTINATION}/signalContamination/  # 4-6 jets: only ST region 1
        done
    done
done

echo "Copying limits plots..."
for PRODUCTIONTYPE in "squark" "gluino"; do
    ${COPY_COMMAND} ${ANALYSIS_SOURCE}/publicationPlots/${PRODUCTIONTYPE}_injectedSignalModel_bestFitSignalStrength.pdf ${AN_DESTINATION}/results/
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

unset COPY_COMMAND
unset ANALYSIS_SOURCE
unset AN_DESTINATION
