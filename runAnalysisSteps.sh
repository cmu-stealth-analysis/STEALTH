#!/bin/bash

cd /uscms/home/tmudholk/private/stealth/STEALTH && source setupEnv.sh && set -x

COMMON_XROOT_PREFIX="root://cmseos.fnal.gov/"
INPUTDATADIR_CONTROL="${COMMON_XROOT_PREFIX}/store/user/lpcsusystealth/selections/testsMerged/oldHLT_noInvMassCut/combinedControl"
INPUTDATADIR_SIGNAL="${COMMON_XROOT_PREFIX}/store/user/lpcsusystealth/selections/testsMerged/oldHLT_noInvMassCut/combinedSignal"
INPUTDATADIR_MCSIGNAL="${COMMON_XROOT_PREFIX}/store/user/lpcsusystealth/selections/combinedSignal"

INTEGLUMI=""
INTEGLUMI_FRACTIONALERROR=""
DATAPATTERNCONTROL=""
DATAPATTERNSIGNAL=""
MCPATTERNSIGNAL="MC_2018Production_DoubleMedium.root"
MCPATTERNSIGNAL_JECUP="MC_2018Production_DoubleMedium_JECUp.root"
MCPATTERNSIGNAL_JECDOWN="MC_2018Production_DoubleMedium_JECDown.root"
MCLUMIYEARIDENTIFIER=""
YEARIDENTIFIER=""

# L_total = L_2016 + L_2017 => deltaL_total/L_total = (deltaL_2016/L_2016)*(L_2016/L_total) + (deltaL_2017/L_2017)*(L_2017/L_total)
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-17-004/index.html, the 2017 uncertainty is 2.3 percent
# From http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/LUM-17-001/index.html, the 2016 uncertainty is 2.5 percent

case ${2} in
    2016)
        INTEGLUMI="37760.0"
        INTEGLUMI_FRACTIONALERROR="0.025"
        DATAPATTERNCONTROL="data_DoubleEG_2016_*_oldHLT_noInvMassCut.root"
        DATAPATTERNSIGNAL="data_DoubleEG_2016_DoubleMedium_oldHLT_noInvMassCut.root"
        MCLUMIYEARIDENTIFIER="_lumi2016"
        YEARIDENTIFIER="_2016"
        ;;
    2017)
        INTEGLUMI="46020.0"
        INTEGLUMI_FRACTIONALERROR="0.023"
        DATAPATTERNCONTROL="data_DoubleEG_2017_*_oldHLT_noInvMassCut.root"
        DATAPATTERNSIGNAL="data_DoubleEG_2017_DoubleMedium_oldHLT_noInvMassCut.root"
        MCLUMIYEARIDENTIFIER="_lumi2017"
        YEARIDENTIFIER="_2017"
        ;;
    2016Plus2017)
        INTEGLUMI="83780.0"
        INTEGLUMI_FRACTIONALERROR="0.024"
        DATAPATTERNCONTROL="data_DoubleEG_201*_oldHLT_noInvMassCut.root"
        DATAPATTERNSIGNAL="data_DoubleEG_201*_DoubleMedium_oldHLT_noInvMassCut.root"
        ;;
    *)
        echo "Unrecognized or empty year: ${2}"
        set +x && exit
        ;;
esac

case ${1} in
    1)
        ./getDataEventHistogramsAndSystematics.py --inputFilePath "${INPUTDATADIR_CONTROL}/${DATAPATTERNCONTROL}" --outputPrefix control${YEARIDENTIFIER} --allowHigherNJets
        ;;
    2)
        ./getDataEventHistogramsAndSystematics.py --inputFilePath "${INPUTDATADIR_SIGNAL}/${DATAPATTERNSIGNAL}" --outputPrefix signal${YEARIDENTIFIER} --isSignal
        ;;
    3)
        ./getMCSystematics/bin/getEventHistograms inputMCPath=${INPUTDATADIR_MCSIGNAL}/${MCPATTERNSIGNAL} inputMCPath_JECUp=${INPUTDATADIR_MCSIGNAL}/${MCPATTERNSIGNAL_JECUP} inputMCPath_JECDown=${INPUTDATADIR_MCSIGNAL}/${MCPATTERNSIGNAL_JECDOWN} outputPrefix=MC_2018${MCLUMIYEARIDENTIFIER} integratedLuminosity=${INTEGLUMI}
        ;;
    4)
        ./getMCSystematics/bin/getMCUncertainties inputPath=analysis/MCEventHistograms/MC_2018${MCLUMIYEARIDENTIFIER}_savedObjects.root outputPrefix=MC_2018${MCLUMIYEARIDENTIFIER}
        ;;
    5)
        ./createDataCards.py --outputPrefix "fullChain${YEARIDENTIFIER}" --inputFile_MCEventHistograms "analysis/MCEventHistograms/MC_2018${MCLUMIYEARIDENTIFIER}_savedObjects.root" --inputFile_MCUncertainties "analysis/MCSystematics/MC_2018${MCLUMIYEARIDENTIFIER}_MCUncertainties_savedObjects.root" --inputFile_dataSystematics "analysis/dataSystematics/signal${YEARIDENTIFIER}_dataSystematics.dat" --inputFile_dataSystematics_sTScaling "analysis/dataSystematics/control${YEARIDENTIFIER}_dataSystematics_sTScaling.dat" --inputFile_dataSystematics_eventCounters "analysis/dataSystematics/signal${YEARIDENTIFIER}_eventCounters.dat" --luminosityUncertainty ${INTEGLUMI_FRACTIONALERROR}
        ;;
    6)
        ./runCombineTool.py --dataCardsPrefix fullChain${YEARIDENTIFIER} --minGluinoMass 975.0
        ;;
    7)
        ./plotLimits.py --combineOutputPrefix fullChain${YEARIDENTIFIER} --outputSuffix fullChain${YEARIDENTIFIER} --minGluinoMass 1000.0 --maxGluinoMass 1750.0
        ;;
    *)
        echo "Unrecognized or empty step index: ${1}"
        set +x && exit
        ;;
esac

set +x
