# Meant to be sourced from this folder

cd miscScripts
make
MAKERESULT=$?
if [ ! ${MAKERESULT} -eq 0 ]; then
    echo "ERROR: make failed with status ${MAKERESULT}"
    cd ..
    return
fi
cd ..

# SELECTION="singlemedium"
# IDENTIFIER="MC_GJet17"
# YEARSTRING="2017"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/normBinOptimization_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries_normOptimization.dat"
# PDFNSTBINS="25"
# PDFSTMIN="1000.0"
# PDFSTMAX="3500.0"
# PARAMETERSINPUTFILE=""

# SELECTION="singlemedium"
# IDENTIFIER="data"
# YEARSTRING="2017"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/normBinOptimization_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries_normOptimization.dat"
# PDFNSTBINS="25"
# PDFSTMIN="1000.0"
# PDFSTMAX="3500.0"
# PARAMETERSINPUTFILE="${ANALYSISROOT}/normBinOptimization_singlephoton/fitParameters_unbinnedFit_${YEARSTRING}_MC_GJet17_${SELECTION}.dat"

# SELECTION="singleloose"
# IDENTIFIER="MC_GJet17"
# YEARSTRING="2017"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/normBinOptimization_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries_normOptimization.dat"
# PDFNSTBINS="25"
# PDFSTMIN="1000.0"
# PDFSTMAX="3500.0"
# PARAMETERSINPUTFILE=""

# SELECTION="singleloose"
# IDENTIFIER="data"
# YEARSTRING="2017"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/normBinOptimization_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries_normOptimization.dat"
# PDFNSTBINS="25"
# PDFSTMIN="1000.0"
# PDFSTMAX="3500.0"
# PARAMETERSINPUTFILE="${ANALYSISROOT}/normBinOptimization_singlephoton/fitParameters_unbinnedFit_${YEARSTRING}_MC_GJet17_${SELECTION}.dat"

# SELECTION="signal"
# IDENTIFIER="MC_GJet"
# YEARSTRING="all"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/normBinOptimization_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries_normOptimization.dat"
# PDFNSTBINS="25"
# PDFSTMIN="1000.0"
# PDFSTMAX="3500.0"
# PARAMETERSINPUTFILE=""

SELECTION="signal_loose"
IDENTIFIER="MC_GJet"
YEARSTRING="all"
SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
OUTPUTFOLDER="${ANALYSISROOT}/normBinOptimization_doublephoton"
STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries_normOptimization.dat"
PDFNSTBINS="25"
PDFSTMIN="1000.0"
PDFSTMAX="3500.0"
PARAMETERSINPUTFILE=""

READPARAMETERSFROMFILE=""
if [[ -n "${PARAMETERSINPUTFILE}" ]]; then
    READPARAMETERSFROMFILE=" readParametersFromFile=${PARAMETERSINPUTFILE}"
fi

echo running ./miscScripts/bin/optimizeNormBin "sourceFilePath=${SOURCEFILEPATH}" "outputFolder=${OUTPUTFOLDER}" "selection=${SELECTION}" "identifier=${IDENTIFIER}" "yearString=${YEARSTRING}" "STBoundariesSourceFile=${STBOUNDARIESSOURCEFILE}" "PDF_nSTBins=${PDFNSTBINS}" "PDF_STMin=${PDFSTMIN}" "PDF_STMax=${PDFSTMAX}" ${READPARAMETERSFROMFILE}
./miscScripts/bin/optimizeNormBin "sourceFilePath=${SOURCEFILEPATH}" "outputFolder=${OUTPUTFOLDER}" "selection=${SELECTION}" "identifier=${IDENTIFIER}" "yearString=${YEARSTRING}" "STBoundariesSourceFile=${STBOUNDARIESSOURCEFILE}" "PDF_nSTBins=${PDFNSTBINS}" "PDF_STMin=${PDFSTMIN}" "PDF_STMax=${PDFSTMAX}" ${READPARAMETERSFROMFILE}
