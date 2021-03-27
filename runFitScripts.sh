# Meant to be sourced from this folder

cd fitScripts
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
# OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
# PDFNSTBINS="50"
# PARAMETERSINPUTFILES=""

# SELECTION="singlemedium"
# IDENTIFIER="data"
# YEARSTRING="2017"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
# PDFNSTBINS="50"
# PARAMETERSINPUTFILES="${ANALYSISROOT}/fits_singlephoton/unbinned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat,${ANALYSISROOT}/fits_singlephoton/binned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat"

# SELECTION="singleloose"
# IDENTIFIER="MC_GJet17"
# YEARSTRING="2017"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
# PDFNSTBINS="50"
# PARAMETERSINPUTFILES=""

SELECTION="singleloose"
IDENTIFIER="data"
YEARSTRING="2017"
SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton"
STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
PDFNSTBINS="50"
PARAMETERSINPUTFILES="${ANALYSISROOT}/fits_singlephoton/unbinned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat,${ANALYSISROOT}/fits_singlephoton/binned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat"

# SELECTION="signal"
# IDENTIFIER="MC_GJet"
# YEARSTRING="all"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat"
# PDFNSTBINS="25"
# PARAMETERSINPUTFILES=""

# SELECTION="signal_loose"
# IDENTIFIER="MC_GJet"
# YEARSTRING="all"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat"
# PDFNSTBINS="25"
# PARAMETERSINPUTFILES=""

READPARAMETERSFROMFILES=""
if [[ -n "${PARAMETERSINPUTFILES}" ]]; then
    READPARAMETERSFROMFILES=" readParametersFromFiles=${PARAMETERSINPUTFILES}"
fi

echo running ./fitScripts/bin/runFits "sourceFilePath=${SOURCEFILEPATH}" "outputFolder=${OUTPUTFOLDER}" "selection=${SELECTION}" "identifier=${IDENTIFIER}" "yearString=${YEARSTRING}" "STBoundariesSourceFile=${STBOUNDARIESSOURCEFILE}" "PDF_nSTBins=${PDFNSTBINS}" ${READPARAMETERSFROMFILES}
./fitScripts/bin/runFits "sourceFilePath=${SOURCEFILEPATH}" "outputFolder=${OUTPUTFOLDER}" "selection=${SELECTION}" "identifier=${IDENTIFIER}" "yearString=${YEARSTRING}" "STBoundariesSourceFile=${STBOUNDARIESSOURCEFILE}" "PDF_nSTBins=${PDFNSTBINS}" ${READPARAMETERSFROMFILES}
