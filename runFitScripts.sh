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

DEBUG="false"
if [ $# -eq 1 ]; then
    if [ "${1}" == "d" ]; then
        DEBUG="true"
    else
        echo "ERROR: unrecognized arguments: "
        echo "$@"
        return
    fi
else
    if [ $# -gt 1 ]; then
        echo "ERROR: at most one argument expected. Current list of arguments: "
        echo "$@"
        return
    fi
fi

# SELECTION="singlemedium"
# IDENTIFIER="MC_GJet17"
# YEARSTRING="2017"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
# PDFNSTBINS="50"
# PARAMETERSINPUTFILES=""
# PLOTCONCISE="true"

# SELECTION="singlemedium"
# IDENTIFIER="data"
# YEARSTRING="2017"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
# PDFNSTBINS="50"
# PARAMETERSINPUTFILES="${ANALYSISROOT}/fits_singlephoton/unbinned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat,${ANALYSISROOT}/fits_singlephoton/binned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat,${STEALTH_ROOT}/STRegionBoundaries.dat"
# PLOTCONCISE="true"

SELECTION="singleloose"
IDENTIFIER="MC_GJet17"
YEARSTRING="2017"
SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton"
STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
PDFNSTBINS="50"
PARAMETERSINPUTFILES=""
PLOTCONCISE="true"

# SELECTION="singleloose"
# IDENTIFIER="data"
# YEARSTRING="2017"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
# PDFNSTBINS="50"
# PARAMETERSINPUTFILES="${ANALYSISROOT}/fits_singlephoton/unbinned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat,${ANALYSISROOT}/fits_singlephoton/binned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat,${STEALTH_ROOT}/STRegionBoundaries.dat"
# PLOTCONCISE="true"

# SELECTION="signal"
# IDENTIFIER="MC_GJet"
# YEARSTRING="all"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat" # important to get proper adjustments to be given as input to the combine tool
# PDFNSTBINS="25"
# PARAMETERSINPUTFILES=""
# PLOTCONCISE=""

# SELECTION="signal_loose"
# IDENTIFIER="MC_GJet"
# YEARSTRING="all"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat" # important to get proper adjustments to be given as input to the combine tool
# PDFNSTBINS="25"
# PARAMETERSINPUTFILES=""
# PLOTCONCISE=""

# SELECTION="control"
# IDENTIFIER="MC_GJet"
# YEARSTRING="all"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat" # important to get proper adjustments to be given as input to the combine tool
# PDFNSTBINS="25"
# PARAMETERSINPUTFILES=""
# PLOTCONCISE="true"

# SELECTION="signal"
# IDENTIFIER="MC_QCD"
# YEARSTRING="all"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat" # important to get proper adjustments to be given as input to the combine tool
# PDFNSTBINS="25"
# PARAMETERSINPUTFILES="${ANALYSISROOT}/fits_doublephoton/unbinned_fitParameters_${YEARSTRING}_MC_GJet_${SELECTION}.dat,${ANALYSISROOT}/fits_doublephoton/binned_fitParameters_${YEARSTRING}_MC_GJet_${SELECTION}.dat,${STEALTH_ROOT}/STRegionBoundaries.dat"
# PLOTCONCISE="true"

# SELECTION="signal_loose"
# IDENTIFIER="MC_QCD"
# YEARSTRING="all"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat" # important to get proper adjustments to be given as input to the combine tool
# PDFNSTBINS="25"
# PARAMETERSINPUTFILES="${ANALYSISROOT}/fits_doublephoton/unbinned_fitParameters_${YEARSTRING}_MC_GJet_${SELECTION}.dat,${ANALYSISROOT}/fits_doublephoton/binned_fitParameters_${YEARSTRING}_MC_GJet_${SELECTION}.dat,${STEALTH_ROOT}/STRegionBoundaries.dat"
# PLOTCONCISE="true"

# SELECTION="control"
# IDENTIFIER="MC_QCD"
# YEARSTRING="all"
# SOURCEFILEPATH="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat" # important to get proper adjustments to be given as input to the combine tool
# PDFNSTBINS="25"
# PARAMETERSINPUTFILES="${ANALYSISROOT}/fits_doublephoton/unbinned_fitParameters_${YEARSTRING}_MC_GJet_${SELECTION}.dat,${ANALYSISROOT}/fits_doublephoton/binned_fitParameters_${YEARSTRING}_MC_GJet_${SELECTION}.dat,${STEALTH_ROOT}/STRegionBoundaries.dat"
# PLOTCONCISE="true"

READPARAMETERSFROMFILES=""
if [[ -n "${PARAMETERSINPUTFILES}" ]]; then
    READPARAMETERSFROMFILES=" readParametersFromFiles=${PARAMETERSINPUTFILES}"
fi

ARG_PLOTCONCISE=""
if [[ -n "${PLOTCONCISE}" ]]; then
    ARG_PLOTCONCISE=" plotConcise=true"
fi

if [ "${DEBUG}" == "false" ]; then
    echo running ./fitScripts/bin/runFits "sourceFilePath=${SOURCEFILEPATH}" "outputFolder=${OUTPUTFOLDER}" "selection=${SELECTION}" "identifier=${IDENTIFIER}" "yearString=${YEARSTRING}" "STBoundariesSourceFile=${STBOUNDARIESSOURCEFILE}" "PDF_nSTBins=${PDFNSTBINS}"${READPARAMETERSFROMFILES}${ARG_PLOTCONCISE}
    ./fitScripts/bin/runFits "sourceFilePath=${SOURCEFILEPATH}" "outputFolder=${OUTPUTFOLDER}" "selection=${SELECTION}" "identifier=${IDENTIFIER}" "yearString=${YEARSTRING}" "STBoundariesSourceFile=${STBOUNDARIESSOURCEFILE}" "PDF_nSTBins=${PDFNSTBINS}"${READPARAMETERSFROMFILES}${ARG_PLOTCONCISE}
elif [ "${DEBUG}" == "true" ]; then
    echo running gdb --args fitScripts/bin/runFits "sourceFilePath=${SOURCEFILEPATH}" "outputFolder=${OUTPUTFOLDER}" "selection=${SELECTION}" "identifier=${IDENTIFIER}" "yearString=${YEARSTRING}" "STBoundariesSourceFile=${STBOUNDARIESSOURCEFILE}" "PDF_nSTBins=${PDFNSTBINS}"${READPARAMETERSFROMFILES}${ARG_PLOTCONCISE}
    gdb --args fitScripts/bin/runFits "sourceFilePath=${SOURCEFILEPATH}" "outputFolder=${OUTPUTFOLDER}" "selection=${SELECTION}" "identifier=${IDENTIFIER}" "yearString=${YEARSTRING}" "STBoundariesSourceFile=${STBOUNDARIESSOURCEFILE}" "PDF_nSTBins=${PDFNSTBINS}"${READPARAMETERSFROMFILES}${ARG_PLOTCONCISE}
fi

unset SELECTION
unset IDENTIFIER
unset YEARSTRING
unset SOURCEFILEPATH
unset OUTPUTFOLDER
unset STBOUNDARIESSOURCEFILE
unset PDFNSTBINS
unset PARAMETERSINPUTFILES
unset PLOTCONCISE
