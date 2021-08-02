#!/bin/bash

OUTPUT_FOLDER_TMP="${1}"
DATACARD_TEMPLATE_PARENT_FOLDER="${2}"
DATACARD_TEMPLATE_FILE_NAME="${3}"
OUTPUT_FOLDER_WITH_EOS_PREFIX="${4}"

echo "Creating fit diagnostics file from template: ${DATACARD_TEMPLATE_PARENT_FOLDER}/${DATACARD_TEMPLATE_FILE_NAME}, in temporary output folder: ${OUTPUT_FOLDER_TMP}, then copying file EOS folder: ${OUTPUT_FOLDER_WITH_EOS_PREFIX}"

mkdir -p ${OUTPUT_FOLDER_TMP}
cd ${OUTPUT_FOLDER_TMP}
xrdcp --nopbar --force ${DATACARD_TEMPLATE_PARENT_FOLDER}/${DATACARD_TEMPLATE_FILE_NAME} ./

# Need to add the line "shapes * * FAKE" to "fake" a shape analysis, because that is required by the FitDiagnostics tool
sed -i '/number of nuisance parameters/a shapes * * FAKE' ${DATACARD_TEMPLATE_FILE_NAME}

combine -M FitDiagnostics -d ${DATACARD_TEMPLATE_FILE_NAME} --forceRecreateNLL --expectSignal 0 --saveShapes --saveWithUncertainties

xrdcp --nopbar --force fitDiagnostics.root ${OUTPUT_FOLDER_WITH_EOS_PREFIX}/fitDiagnostics.root
rm -v -f combine_logger.out ${DATACARD_TEMPLATE_FILE_NAME} higgsCombineTest.FitDiagnostics.mH120.root fitDiagnostics.root
