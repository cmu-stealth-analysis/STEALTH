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
# PLOTCONCISE=""

# SELECTION="singlemedium"
# IDENTIFIER="data"
# YEARSTRING="2017"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
# PDFNSTBINS="50"
# PARAMETERSINPUTFILES="${ANALYSISROOT}/fits_singlephoton/unbinned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat,${ANALYSISROOT}/fits_singlephoton/binned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat"
# PLOTCONCISE=""

SELECTION="singleloose"
IDENTIFIER="MC_GJet17"
YEARSTRING="2017"
SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton"
STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
PDFNSTBINS="50"
PARAMETERSINPUTFILES=""
PLOTCONCISE=""

# SELECTION="singleloose"
# IDENTIFIER="data"
# YEARSTRING="2017"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
# PDFNSTBINS="50"
# PARAMETERSINPUTFILES="${ANALYSISROOT}/fits_singlephoton/unbinned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat,${ANALYSISROOT}/fits_singlephoton/binned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat"
# PLOTCONCISE=""

# SELECTION="signal"
# IDENTIFIER="MC_GJet"
# YEARSTRING="all"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat" # important to get proper adjustments to be given as input to the combine tool
# PDFNSTBINS="25"
# PARAMETERSINPUTFILES=""
# PLOTCONCISE="true"

# SELECTION="signal_loose"
# IDENTIFIER="MC_GJet"
# YEARSTRING="all"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat" # important to get proper adjustments to be given as input to the combine tool
# PDFNSTBINS="25"
# PARAMETERSINPUTFILES=""
# PLOTCONCISE="true"

# SELECTION="control"
# IDENTIFIER="MC_GJet"
# YEARSTRING="all"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat" # important to get proper adjustments to be given as input to the combine tool
# PDFNSTBINS="25"
# PARAMETERSINPUTFILES=""
# PLOTCONCISE="true"

READPARAMETERSFROMFILES=""
if [[ -n "${PARAMETERSINPUTFILES}" ]]; then
    READPARAMETERSFROMFILES=" readParametersFromFiles=${PARAMETERSINPUTFILES}"
fi

ARG_PLOTCONCISE=""
if [[ -n "${PLOTCONCISE}" ]]; then
    ARG_PLOTCONCISE=" plotConcise=true"
fi

echo running ./fitScripts/bin/runFits "sourceFilePath=${SOURCEFILEPATH}" "outputFolder=${OUTPUTFOLDER}" "selection=${SELECTION}" "identifier=${IDENTIFIER}" "yearString=${YEARSTRING}" "STBoundariesSourceFile=${STBOUNDARIESSOURCEFILE}" "PDF_nSTBins=${PDFNSTBINS}"${READPARAMETERSFROMFILES}${ARG_PLOTCONCISE}
./fitScripts/bin/runFits "sourceFilePath=${SOURCEFILEPATH}" "outputFolder=${OUTPUTFOLDER}" "selection=${SELECTION}" "identifier=${IDENTIFIER}" "yearString=${YEARSTRING}" "STBoundariesSourceFile=${STBOUNDARIESSOURCEFILE}" "PDF_nSTBins=${PDFNSTBINS}"${READPARAMETERSFROMFILES}${ARG_PLOTCONCISE}

unset SELECTION
unset IDENTIFIER
unset YEARSTRING
unset SOURCEFILEPATH
unset OUTPUTFOLDER
unset STBOUNDARIESSOURCEFILE
unset PDFNSTBINS
unset PARAMETERSINPUTFILES
unset PLOTCONCISE
