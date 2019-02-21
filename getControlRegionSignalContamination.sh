#!/bin/bash

./getDataEventHistogramsAndSystematics.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/data_DoubleEG_201*_DoubleFake.root" --outputPrefix control_fakefake && ./getDataEventHistogramsAndSystematics.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/data_DoubleEG_201*_OneMediumOneFake.root" --outputPrefix control_mediumfake

./getMCSystematics/bin/getEventHistograms inputMCPathMain=${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/MC_2018Production_DoubleFake_optimized2017.root integratedLuminosityMain=41900.0 inputMCPathAux=${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/MC_2018Production_DoubleFake_optimized2016.root integratedLuminosityAux=35920.0 outputPrefix=fakefake && ./getMCSystematics/bin/getEventHistograms inputMCPathMain=${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/MC_2018Production_OneMediumOneFake_optimized2017.root integratedLuminosityMain=41900.0 inputMCPathAux=${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/MC_2018Production_OneMediumOneFake_optimized2016.root integratedLuminosityAux=35920.0 outputPrefix=mediumfake

./getMCSystematics/bin/getMCUncertainties inputNEventsFile=analysis/dataSystematics/control_fakefake_observedEventCounters.dat inputPath=analysis/MCEventHistograms/fakefake_savedObjects.root outputPrefix=fakefake unrestrictedSignalContamination=true && ./getMCSystematics/bin/getMCUncertainties inputNEventsFile=analysis/dataSystematics/control_mediumfake_observedEventCounters.dat inputPath=analysis/MCEventHistograms/mediumfake_savedObjects.root outputPrefix=mediumfake unrestrictedSignalContamination=true

rm -r analysis/MCSystematics/mediumfake_* && rm -r analysis/MCSystematics/fakefake_*
rm -r analysis/MCEventHistograms/mediumfake_* && rm -r analysis/MCEventHistograms/fakefake_*
