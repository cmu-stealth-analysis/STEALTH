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

SELECTION="singleloose"
IDENTIFIER="MC_GJet17"
YEARSTRING="2017"
SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
OUTPUTFOLDER="${ANALYSISROOT}/normBinOptimization_singlephoton"
STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries_normOptimization.dat"
PDFNSTBINS="25"
PDFSTMIN="1000.0"
PDFSTMAX="3500.0"

# SELECTION="singlemedium"
# IDENTIFIER="data"
# YEARSTRING="2017"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/normBinOptimization_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries_normOptimization.dat"
# PDFNSTBINS="25"
# PDFSTMIN="1000.0"
# PDFSTMAX="3500.0"

# SELECTION="singleloose"
# IDENTIFIER="data"
# YEARSTRING="2017"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/normBinOptimization_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries_normOptimization.dat"
# PDFNSTBINS="13"
# PDFSTMIN="900.0"
# PDFSTMAX="3500.0"

echo running ./miscScripts/bin/optimizeNormBin "sourceFilePath=${SOURCEFILEPATH}" "outputFolder=${OUTPUTFOLDER}" "selection=${SELECTION}" "identifier=${IDENTIFIER}" "yearString=${YEARSTRING}" "STBoundariesSourceFile=${STBOUNDARIESSOURCEFILE}"
./miscScripts/bin/optimizeNormBin "sourceFilePath=${SOURCEFILEPATH}" "outputFolder=${OUTPUTFOLDER}" "selection=${SELECTION}" "identifier=${IDENTIFIER}" "yearString=${YEARSTRING}" "STBoundariesSourceFile=${STBOUNDARIESSOURCEFILE}" "PDF_nSTBins=${PDFNSTBINS}" "PDF_STMin=${PDFSTMIN}" "PDF_STMax=${PDFSTMAX}"
