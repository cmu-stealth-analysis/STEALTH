#!/bin/bash

SOURCE_WITH_XRD_PREFIX="root://cmseos.fnal.gov//store/user/lpcsusystealth/analysisEOSAreas/analysis/dataCards/combinedFit"
TARGET_FOLDER="/uscms/home/tmudholk/nobackup/datacards"

function xrdcp_with_check {
    if [ "${#}" != 2  ]; then
        echo "ERROR: number of arguments passed to \"${FUNCNAME}\": ${#}"
        exit 1
    fi
    xrdcp --nopbar --force --path --streams 15 ${1} ${2} 2>&1
    XRDEXIT=${?}
    if [[ ${XRDEXIT} -ne 0 ]]; then
        echo "exit code ${XRDEXIT}, failure in xrdcp"
        exit ${XRDEXIT}
    fi
}

MASS_BINS=("gluino_23_73" "gluino_20_141" "gluino_18_9" "squark_21_57" "squark_19_125" "squark_14_9")
MASS_BIN_TITLES=("gluinoMass2100_neutralinoMass1000" "gluinoMass1950_neutralinoMass1850" "gluinoMass1850_neutralinoMass200" "squarkMass1850_neutralinoMass800" "squarkMass1750_neutralinoMass1650" "squarkMass1500_neutralinoMass200")

for INDEX in ${!MASS_BINS[@]}; do
    MASS_BIN=${MASS_BINS[${INDEX}]}
    GLUINO_OR_SQUARK=`echo ${MASS_BIN} | cut --delimiter="_" --fields=1`
    GS_MASS_BIN_INDEX=`echo ${MASS_BIN} | cut --delimiter="_" --fields=2`
    NEUT_MASS_BIN_INDEX=`echo ${MASS_BIN} | cut --delimiter="_" --fields=3`
    MASS_BIN_TITLE=${MASS_BIN_TITLES[${INDEX}]}
    xrdcp_with_check "${SOURCE_WITH_XRD_PREFIX}/${GLUINO_OR_SQUARK}_dataCard_eventProgenitorMassBin${GS_MASS_BIN_INDEX}_neutralinoMassBin${NEUT_MASS_BIN_INDEX}.txt" "${TARGET_FOLDER}/dataCard_${MASS_BIN_TITLE}.txt"
done

echo "Cards copied."

# To compare asymptotic with full limits:

# (expected limits)
# combine -M AsymptoticLimits -d "dataCard_gluinoMass2100_neutralinoMass1000.txt" -n "_gluinoMass2100_neutralinoMass1000" -v 1 -V --expectSignal 0 && combine -M AsymptoticLimits -d "dataCard_gluinoMass1850_neutralinoMass200.txt" -n "_gluinoMass1850_neutralinoMass200" -v 1 -V --expectSignal 0 && combine -M AsymptoticLimits -d "dataCard_gluinoMass1950_neutralinoMass1850.txt" -n "_gluinoMass1950_neutralinoMass1850" -v 1 -V --expectSignal 0 && combine -M AsymptoticLimits -d "dataCard_squarkMass1500_neutralinoMass200.txt" -n "_squarkMass1500_neutralinoMass200" -v 1 -V --expectSignal 0 && combine -M AsymptoticLimits -d "dataCard_squarkMass1850_neutralinoMass800.txt" -n "_squarkMass1850_neutralinoMass800" -v 1 -V --expectSignal 0 && combine -M AsymptoticLimits -d "dataCard_squarkMass1750_neutralinoMass1650.txt" -n "_squarkMass1750_neutralinoMass1650" -v 1 -V --expectSignal 0 && combine -M HybridNew --LHCmode LHC-limits -d "dataCard_gluinoMass2100_neutralinoMass1000.txt" -n "_gluinoMass2100_neutralinoMass1000" -v 1 -V --expectSignal 0 --expectedFromGrid=0.5 && combine -M HybridNew --LHCmode LHC-limits -d "dataCard_gluinoMass1850_neutralinoMass200.txt" -n "_gluinoMass1850_neutralinoMass200" -v 1 -V --expectSignal 0 --expectedFromGrid=0.5 && combine -M HybridNew --LHCmode LHC-limits -d "dataCard_gluinoMass1950_neutralinoMass1850.txt" -n "_gluinoMass1950_neutralinoMass1850" -v 1 -V --expectSignal 0 --expectedFromGrid=0.5 && combine -M HybridNew --LHCmode LHC-limits -d "dataCard_squarkMass1500_neutralinoMass200.txt" -n "_squarkMass1500_neutralinoMass200" -v 1 -V --expectSignal 0 --expectedFromGrid=0.5 && combine -M HybridNew --LHCmode LHC-limits -d "dataCard_squarkMass1850_neutralinoMass800.txt" -n "_squarkMass1850_neutralinoMass800" -v 1 -V --expectSignal 0 --expectedFromGrid=0.5 && combine -M HybridNew --LHCmode LHC-limits -d "dataCard_squarkMass1750_neutralinoMass1650.txt" -n "_squarkMass1750_neutralinoMass1650" -v 1 -V --expectSignal 0 --expectedFromGrid=0.5

# (observed limits)
# combine -M HybridNew --LHCmode LHC-limits -d "dataCard_gluinoMass2100_neutralinoMass1000.txt" -n "_gluinoMass2100_neutralinoMass1000" -v 1 -V --expectSignal 0 && combine -M HybridNew --LHCmode LHC-limits -d "dataCard_gluinoMass1850_neutralinoMass200.txt" -n "_gluinoMass1850_neutralinoMass200" -v 1 -V --expectSignal 0 && combine -M HybridNew --LHCmode LHC-limits -d "dataCard_gluinoMass1950_neutralinoMass1850.txt" -n "_gluinoMass1950_neutralinoMass1850" -v 1 -V --expectSignal 0 && combine -M HybridNew --LHCmode LHC-limits -d "dataCard_squarkMass1500_neutralinoMass200.txt" -n "_squarkMass1500_neutralinoMass200" -v 1 -V --expectSignal 0 && combine -M HybridNew --LHCmode LHC-limits -d "dataCard_squarkMass1850_neutralinoMass800.txt" -n "_squarkMass1850_neutralinoMass800" -v 1 -V --expectSignal 0 && combine -M HybridNew --LHCmode LHC-limits -d "dataCard_squarkMass1750_neutralinoMass1650.txt" -n "_squarkMass1750_neutralinoMass1650" -v 1 -V --expectSignal 0
