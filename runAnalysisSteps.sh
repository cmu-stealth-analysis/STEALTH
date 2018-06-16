#!/bin/bash

check_for_one_decimal_point() {
    case "$@" in  
        *.*.*)
            echo "At least two decimal points: $@"
            exit
            ;;
        *.*)
            echo "Found in expected format: $@"
            ;;
        *)
            echo "No decimal point found: $@"
            exit
            ;;
    esac
}

checkNEntries() {
    echo "NArgs:"
    echo "$#"
    case "$#" in
        1)
            echo "Found in expected format."
            ;;
        0)
            echo "No arguments passed."
            exit
            ;;
        *)
            echo "More than one argument passed"
            exit
            ;;
    esac
}

COMMON_XROOT_PREFIX="root://cmseos.fnal.gov/"
INPUTDATADIR="${COMMON_XROOT_PREFIX}/store/user/lpcsusystealth/selections/combined"

cd /uscms/home/tmudholk/private/stealth/STEALTH && source setupEnv.sh && set -x

INTEGLUMI=""
DATAPATTERNDOUBLEFAKE=""
DATAPATTERNMEDIUMFAKE=""
DATAPATTERNDOUBLEMEDIUM=""
MCPATTERNDOUBLEFAKE="MC_2018Production_DoubleFake.root"
MCPATTERNMEDIUMFAKE="MC_2018Production_OneMediumOneFake.root"
MCPATTERNDOUBLEMEDIUM="MC_2018Production_DoubleMedium.root"

case ${2} in
    2016)
        INTEGLUMI="37760.0"
        DATAPATTERNDOUBLEFAKE="data_DoubleEG_2016_DoubleFake.root"
        DATAPATTERNMEDIUMFAKE="data_DoubleEG_2016_OneMediumOneFake.root"
        DATAPATTERNDOUBLEMEDIUM="data_DoubleEG_2016_DoubleMedium.root"
        ;;
    2017)
        INTEGLUMI="46020.0"
        DATAPATTERNDOUBLEFAKE="data_DoubleEG_2017_DoubleFake.root"
        DATAPATTERNMEDIUMFAKE="data_DoubleEG_2017_OneMediumOneFake.root"
        DATAPATTERNDOUBLEMEDIUM="data_DoubleEG_2017_DoubleMedium.root"
        ;;
    2016Plus2017)
        INTEGLUMI="83780.0"
        DATAPATTERNDOUBLEFAKE="data_DoubleEG_201*_DoubleFake.root"
        DATAPATTERNMEDIUMFAKE="data_DoubleEG_201*_OneMediumOneFake.root"
        DATAPATTERNDOUBLEMEDIUM="data_DoubleEG_201*_DoubleMedium.root"
        ;;
    *)
        echo "Unrecognized or empty year: ${2}"
        set +x && exit
        ;;
esac

case ${1} in
    1)
        ./getSystematicsAndCreateDataCardTemplate.py --inputFilePath "${INPUTDATADIR}/${DATAPATTERNDOUBLEFAKE}" --outputDirectoryPrefix DoubleEG_${2}_processedJun18_DoubleFake --allowHigherNJets && ./getSystematicsAndCreateDataCardTemplate.py --inputFilePath "${INPUTDATADIR}/${DATAPATTERNMEDIUMFAKE}" --outputDirectoryPrefix DoubleEG_${2}_processedJun18_OneMediumOneFake --allowHigherNJets
        ;;
    2)
        FOURJETSSCALE_DOUBLEFAKE=`cat analysis/DoubleEG_${2}_processedJun18_DoubleFake_*/sTScalingSystematics.dat | grep "4   " | sed "s|4 [ ]*\(.*$\)|\1|"`
        check_for_one_decimal_point "${FOURJETSSCALE_DOUBLEFAKE}"
        echo "fake+fake: 4 jets scale = ${FOURJETSSCALE_DOUBLEFAKE}"
        FOURJETSSCALE_MEDIUMFAKE=`cat analysis/DoubleEG_${2}_processedJun18_OneMediumOneFake_*/sTScalingSystematics.dat | grep "4   " | sed "s|4 [ ]*\(.*$\)|\1|"`
        check_for_one_decimal_point "${FOURJETSSCALE_MEDIUMFAKE}"
        echo "medium+fake: 4 jets scale = ${FOURJETSSCALE_MEDIUMFAKE}"
        FOURJETSSCALE=`python -c "print max(abs(1.0-${FOURJETSSCALE_DOUBLEFAKE}), abs(1.0-${FOURJETSSCALE_MEDIUMFAKE}))"`
        echo "Overall 4 jets scale = ${FOURJETSSCALE}"
        FIVEJETSSCALE_DOUBLEFAKE=`cat analysis/DoubleEG_${2}_processedJun18_DoubleFake_*/sTScalingSystematics.dat | grep "5   " | sed "s|5 [ ]*\(.*$\)|\1|"`
        check_for_one_decimal_point "${FIVEJETSSCALE_DOUBLEFAKE}"
        echo "fake+fake: 5 jets scale = ${FIVEJETSSCALE_DOUBLEFAKE}"
        FIVEJETSSCALE_MEDIUMFAKE=`cat analysis/DoubleEG_${2}_processedJun18_OneMediumOneFake_*/sTScalingSystematics.dat | grep "5   " | sed "s|5 [ ]*\(.*$\)|\1|"`
        check_for_one_decimal_point "${FIVEJETSSCALE_MEDIUMFAKE}"
        echo "medium+fake: 5 jets scale = ${FIVEJETSSCALE_MEDIUMFAKE}"
        FIVEJETSSCALE=`python -c "print max(abs(1.0-${FIVEJETSSCALE_DOUBLEFAKE}), abs(1.0-${FIVEJETSSCALE_MEDIUMFAKE}))"`
        echo "Overall 5 jets scale = ${FIVEJETSSCALE}"
        SIXJETSSCALE_DOUBLEFAKE=`cat analysis/DoubleEG_${2}_processedJun18_DoubleFake_*/sTScalingSystematics.dat | grep "6   " | sed "s|6 [ ]*\(.*$\)|\1|"`
        check_for_one_decimal_point "${SIXJETSSCALE_DOUBLEFAKE}"
        echo "fake+fake: >=6 jets scale = ${SIXJETSSCALE_DOUBLEFAKE}"
        SIXJETSSCALE_MEDIUMFAKE=`cat analysis/DoubleEG_${2}_processedJun18_OneMediumOneFake_*/sTScalingSystematics.dat | grep "6   " | sed "s|6 [ ]*\(.*$\)|\1|"`
        check_for_one_decimal_point "${SIXJETSSCALE_MEDIUMFAKE}"
        echo "medium+fake: >=6 jets scale = ${SIXJETSSCALE_MEDIUMFAKE}"
        SIXJETSSCALE=`python -c "print max(abs(1.0-${SIXJETSSCALE_DOUBLEFAKE}), abs(1.0-${SIXJETSSCALE_MEDIUMFAKE}))"`
        echo "Overall 6 jets scale = ${SIXJETSSCALE}"
        ./getSystematicsAndCreateDataCardTemplate.py --inputFilePath "${INPUTDATADIR}/${DATAPATTERNDOUBLEMEDIUM}" --outputDirectoryPrefix DoubleEG_${2}_processedJun18_DoubleMedium --sTScalingFractionalUncertainty_4Jets ${FOURJETSSCALE} --sTScalingFractionalUncertainty_5Jets ${FIVEJETSSCALE} --sTScalingFractionalUncertainty_geq6Jets ${SIXJETSSCALE} --generateDataCardTemplate
        ;;
    3)
        DOUBLEMEDIUMDIR=`find analysis/ -type d -name "*DoubleEG_${2}_*DoubleMedium*"`
        checkNEntries "${DOUBLEMEDIUMDIR}"
        ./getSignalContaminationAndCreateDataCards.py --inputMCPath ${INPUTDATADIR}/${MCPATTERNDOUBLEMEDIUM} --inputDataPath ${INPUTDATADIR}/${DATAPATTERNDOUBLEMEDIUM} --dataCardTemplate ${DOUBLEMEDIUMDIR}/dataCardTemplate.txt --outputPrefix data_${2}_processedJun18_DoubleMedium --analyze_nJetsBin 2 --analyze_nJetsBin 3 --totalIntegratedLuminosity ${INTEGLUMI}
        ./getSignalContaminationAndCreateDataCards.py --inputMCPath ${INPUTDATADIR}/${MCPATTERNDOUBLEFAKE} --inputDataPath ${INPUTDATADIR}/${DATAPATTERNDOUBLEFAKE} --dataCardTemplate "/dev/null" --outputPrefix data_${2}_processedJun18_DoubleFake --analyze_nJetsBin 2 --analyze_nJetsBin 3 --analyze_nJetsBin 4 --analyze_nJetsBin 5 --analyze_nJetsBin 6 --totalIntegratedLuminosity ${INTEGLUMI} --skipDataCardCreation
        ./getSignalContaminationAndCreateDataCards.py --inputMCPath ${INPUTDATADIR}/${MCPATTERNMEDIUMFAKE} --inputDataPath ${INPUTDATADIR}/${DATAPATTERNMEDIUMFAKE} --dataCardTemplate "/dev/null" --outputPrefix data_${2}_processedJun18_OneMediumOneFake --analyze_nJetsBin 2 --analyze_nJetsBin 3 --analyze_nJetsBin 4 --analyze_nJetsBin 5 --analyze_nJetsBin 6 --totalIntegratedLuminosity ${INTEGLUMI} --skipDataCardCreation
        ;;
    4)
        ./runCombineTool.py --dataCardsPrefix data_${2}_processedJun18_DoubleMedium --minGluinoMass 975.0
        ;;
    5)
        ./plotLimits.py --dataCardsPrefix data_${2}_processedJun18_DoubleMedium --outputSuffix data_${2}_processedJun18 --minGluinoMass 1000.0 --maxGluinoMass 1750.0
        ;;
    *)
        echo "Unrecognized or empty step index: ${1}"
        set +x && exit
        ;;
esac

set +x
