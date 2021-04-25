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
# SOURCEFILEPATHS="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
# PDFNSTBINS="50"
# PARAMETERSINPUTFILES=""
# PLOTCONCISE="true"

# SELECTION="singlemedium"
# IDENTIFIER="data"
# YEARSTRING="2017"
# SOURCEFILEPATHS="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
# PDFNSTBINS="50"
# PARAMETERSINPUTFILES="${ANALYSISROOT}/fits_singlephoton/unbinned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat,${ANALYSISROOT}/fits_singlephoton/binned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat,${STEALTH_ROOT}/STRegionBoundaries.dat"
# PLOTCONCISE="true"

# SELECTION="singleloose"
# IDENTIFIER="MC_GJet17"
# YEARSTRING="2017"
# SOURCEFILEPATHS="${EOSPREFIX}/store/user/lpcsusystealth/selections/combined_DoublePhoton_lowerSTThreshold/merged_selection_${IDENTIFIER}_singlephoton_${YEARSTRING}_control_${SELECTION}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton_temp"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
# PDFNSTBINS="50"
# PARAMETERSINPUTFILES=""
# PLOTCONCISE="true"
# FETCHMCWEIGHTS="true"
# MINALLOWEDEMST="200.0"

# SELECTION="singleloose"
# IDENTIFIER="data"
# YEARSTRING="2017"
# SOURCEFILEPATHS="${EOSPREFIX}/store/user/lpcsusystealth/selections/combined_DoublePhoton_lowerSTThreshold/merged_selection_${IDENTIFIER}_singlephoton_${YEARSTRING}_control_${SELECTION}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton_temp"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
# PDFNSTBINS="50"
# PARAMETERSINPUTFILES="${ANALYSISROOT}/fits_singlephoton_temp/unbinned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat,${ANALYSISROOT}/fits_singlephoton_temp/binned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat,${STEALTH_ROOT}/STRegionBoundaries.dat"
# PLOTCONCISE="true"
# FETCHMCWEIGHTS="false"
# MINALLOWEDEMST="200.0"

# SELECTION="singlefake"
# IDENTIFIER="MC_GJet17"
# YEARSTRING="2017"
# SOURCEFILEPATHS="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
# PDFNSTBINS="50"
# PARAMETERSINPUTFILES=""
# PLOTCONCISE="true"

# SELECTION="singlefake"
# IDENTIFIER="data"
# YEARSTRING="2017"
# SOURCEFILEPATHS="${ANALYSISROOT}/STDistributions_singlephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_singlephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundariesFineBinned.dat"
# PDFNSTBINS="50"
# PARAMETERSINPUTFILES="${ANALYSISROOT}/fits_singlephoton/unbinned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat,${ANALYSISROOT}/fits_singlephoton/binned_fitParameters_${YEARSTRING}_MC_GJet17_${SELECTION}.dat,${STEALTH_ROOT}/STRegionBoundaries.dat"
# PLOTCONCISE="true"

SELECTION="signal"
IDENTIFIER="MC_GJet"
YEARSTRING="all"
SOURCEFILEPATHS="${EOSPREFIX}/store/user/lpcsusystealth/selections/combined_DoublePhoton_lowerSTThreshold/merged_selection_MC_GJet16_2016_${SELECTION}.root,${EOSPREFIX}/store/user/lpcsusystealth/selections/combined_DoublePhoton_lowerSTThreshold/merged_selection_MC_GJet17_2017_${SELECTION}.root,${EOSPREFIX}/store/user/lpcsusystealth/selections/combined_DoublePhoton_lowerSTThreshold/merged_selection_MC_GJet18_2018_${SELECTION}.root"
OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton_temp"
STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat" # important to get proper adjustments to be given as input to the combine tool
PDFNSTBINS="25"
PARAMETERSINPUTFILES=""
PLOTCONCISE=""
FETCHMCWEIGHTS="true"
MINALLOWEDEMST="-1.0"

# SELECTION="signal_loose"
# IDENTIFIER="MC_GJet"
# YEARSTRING="all"
# SOURCEFILEPATHS="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat" # important to get proper adjustments to be given as input to the combine tool
# PDFNSTBINS="25"
# PARAMETERSINPUTFILES=""
# PLOTCONCISE=""

# SELECTION="control"
# IDENTIFIER="MC_GJet"
# YEARSTRING="all"
# SOURCEFILEPATHS="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat" # important to get proper adjustments to be given as input to the combine tool
# PDFNSTBINS="25"
# PARAMETERSINPUTFILES=""
# PLOTCONCISE="false"

# SELECTION="control"
# IDENTIFIER="data"
# YEARSTRING="all"
# SOURCEFILEPATHS="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat" # important to get proper adjustments to be given as input to the combine tool
# PDFNSTBINS="25"
# PARAMETERSINPUTFILES="${ANALYSISROOT}/fits_doublephoton/unbinned_fitParameters_${YEARSTRING}_MC_GJet_${SELECTION}.dat,${ANALYSISROOT}/fits_doublephoton/binned_fitParameters_${YEARSTRING}_MC_GJet_${SELECTION}.dat,${STEALTH_ROOT}/STRegionBoundaries.dat"
# PLOTCONCISE="true"

# SELECTION="signal"
# IDENTIFIER="MC_QCD"
# YEARSTRING="all"
# SOURCEFILEPATHS="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat" # important to get proper adjustments to be given as input to the combine tool
# PDFNSTBINS="25"
# PARAMETERSINPUTFILES="${ANALYSISROOT}/fits_doublephoton/unbinned_fitParameters_${YEARSTRING}_MC_GJet_${SELECTION}.dat,${ANALYSISROOT}/fits_doublephoton/binned_fitParameters_${YEARSTRING}_MC_GJet_${SELECTION}.dat,${STEALTH_ROOT}/STRegionBoundaries.dat"
# PLOTCONCISE="true"

# SELECTION="signal_loose"
# IDENTIFIER="MC_QCD"
# YEARSTRING="all"
# SOURCEFILEPATHS="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
# OUTPUTFOLDER="${ANALYSISROOT}/fits_doublephoton"
# STBOUNDARIESSOURCEFILE="${STEALTH_ROOT}/STRegionBoundaries.dat" # important to get proper adjustments to be given as input to the combine tool
# PDFNSTBINS="25"
# PARAMETERSINPUTFILES="${ANALYSISROOT}/fits_doublephoton/unbinned_fitParameters_${YEARSTRING}_MC_GJet_${SELECTION}.dat,${ANALYSISROOT}/fits_doublephoton/binned_fitParameters_${YEARSTRING}_MC_GJet_${SELECTION}.dat,${STEALTH_ROOT}/STRegionBoundaries.dat"
# PLOTCONCISE="true"

# SELECTION="control"
# IDENTIFIER="MC_QCD"
# YEARSTRING="all"
# SOURCEFILEPATHS="${ANALYSISROOT}/STDistributions_doublephoton/distributions_${YEARSTRING}_${SELECTION}_${IDENTIFIER}.root"
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

ARG_FETCHMCWEIGHTS=""
if [ "${FETCHMCWEIGHTS}" == "false" ]; then
    ARG_FETCHMCWEIGHTS=" fetchMCWeights=false"
else
    ARG_FETCHMCWEIGHTS=" fetchMCWeights=true"
fi

if [ "${DEBUG}" == "false" ]; then
    echo running ./fitScripts/bin/runFits "sourceFilePaths=${SOURCEFILEPATHS}" "outputFolder=${OUTPUTFOLDER}" "selection=${SELECTION}" "identifier=${IDENTIFIER}" "yearString=${YEARSTRING}" "STBoundariesSourceFile=${STBOUNDARIESSOURCEFILE}" "PDF_nSTBins=${PDFNSTBINS}" "minAllowedEMST=${MINALLOWEDEMST}"${READPARAMETERSFROMFILES}${ARG_PLOTCONCISE}${ARG_FETCHMCWEIGHTS}
    ./fitScripts/bin/runFits "sourceFilePaths=${SOURCEFILEPATHS}" "outputFolder=${OUTPUTFOLDER}" "selection=${SELECTION}" "identifier=${IDENTIFIER}" "yearString=${YEARSTRING}" "STBoundariesSourceFile=${STBOUNDARIESSOURCEFILE}" "PDF_nSTBins=${PDFNSTBINS}" "minAllowedEMST=${MINALLOWEDEMST}"${READPARAMETERSFROMFILES}${ARG_PLOTCONCISE}${ARG_FETCHMCWEIGHTS}
elif [ "${DEBUG}" == "true" ]; then
    echo running gdb --args fitScripts/bin/runFits "sourceFilePaths=${SOURCEFILEPATHS}" "outputFolder=${OUTPUTFOLDER}" "selection=${SELECTION}" "identifier=${IDENTIFIER}" "yearString=${YEARSTRING}" "STBoundariesSourceFile=${STBOUNDARIESSOURCEFILE}" "PDF_nSTBins=${PDFNSTBINS}" "minAllowedEMST=${MINALLOWEDEMST}"${READPARAMETERSFROMFILES}${ARG_PLOTCONCISE}${ARG_FETCHMCWEIGHTS}
    gdb --args fitScripts/bin/runFits "sourceFilePaths=${SOURCEFILEPATHS}" "outputFolder=${OUTPUTFOLDER}" "selection=${SELECTION}" "identifier=${IDENTIFIER}" "yearString=${YEARSTRING}" "STBoundariesSourceFile=${STBOUNDARIESSOURCEFILE}" "PDF_nSTBins=${PDFNSTBINS}" "minAllowedEMST=${MINALLOWEDEMST}"${READPARAMETERSFROMFILES}${ARG_PLOTCONCISE}${ARG_FETCHMCWEIGHTS}
fi

unset SELECTION
unset IDENTIFIER
unset YEARSTRING
unset SOURCEFILEPATHS
unset OUTPUTFOLDER
unset STBOUNDARIESSOURCEFILE
unset PDFNSTBINS
unset PARAMETERSINPUTFILES
unset PLOTCONCISE
