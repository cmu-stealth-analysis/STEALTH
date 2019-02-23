#!/bin/bash

# ./getDataEventHistogramsAndSystematics.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/data_DoubleEG_201*_DoubleFake.root" --outputPrefix control_fakefake && ./getDataEventHistogramsAndSystematics.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/data_DoubleEG_201*_OneMediumOneFake.root" --outputPrefix control_mediumfake

# ./getMCSystematics/bin/getEventHistograms inputMCPathMain=${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/MC_2018Production_DoubleFake_optimized2017.root integratedLuminosityMain=41900.0 inputMCPathAux=${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/MC_2018Production_DoubleFake_optimized2016.root integratedLuminosityAux=35920.0 outputPrefix=fakefake && ./getMCSystematics/bin/getEventHistograms inputMCPathMain=${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/MC_2018Production_OneMediumOneFake_optimized2017.root integratedLuminosityMain=41900.0 inputMCPathAux=${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/MC_2018Production_OneMediumOneFake_optimized2016.root integratedLuminosityAux=35920.0 outputPrefix=mediumfake

# ./getMCSystematics/bin/getMCUncertainties inputNEventsFile=analysis/dataSystematics/control_fakefake_observedEventCounters.dat inputPath=analysis/MCEventHistograms/fakefake_savedObjects.root outputPrefix=fakefake unrestrictedSignalContamination=true minGluinoMass=975.0 nGluinoMassBins=16 && ./getMCSystematics/bin/getMCUncertainties inputNEventsFile=analysis/dataSystematics/control_mediumfake_observedEventCounters.dat inputPath=analysis/MCEventHistograms/mediumfake_savedObjects.root outputPrefix=mediumfake unrestrictedSignalContamination=true minGluinoMass=975.0 nGluinoMassBins=16

# rm -r analysis/MCEventHistograms/mediumfake_* && rm -r analysis/MCEventHistograms/fakefake_*
# rm -r analysis/MCSystematics/mediumfake_* && rm -r analysis/MCSystematics/fakefake_*

./getMCSystematics/bin/getEventHistograms inputMCPathMain=${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/MC_2018Production_control_optimized2017.root integratedLuminosityMain=41900.0 inputMCPathAux=${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/MC_2018Production_control_optimized2016.root integratedLuminosityAux=35920.0 outputPrefix=control && ./getMCSystematics/bin/getMCUncertainties inputNEventsFile=analysis/dataSystematics/control_observedEventCounters.dat inputPath=analysis/MCEventHistograms/control_savedObjects.root outputPrefix=control unrestrictedSignalContamination=true minGluinoMass=975.0 nGluinoMassBins=16

./getMCSystematics/bin/getMCUncertainties inputPath=analysis/MCEventHistograms/MC_2018_savedObjects.root outputPrefix=MC_2018NewMinimum minGluinoMass=975.0 nGluinoMassBins=16 # Reproduce signal contamination plots with new minGluinoMass
rm -r analysis/MCSystematics/MC_2018NewMinimum*
for file in analysis/signalContamination/MC_2018NewMinimum_signalContamination_*; do
    suffix=`echo ${file} | cut --delimiter="/" --fields=3 | cut --delimiter="_" --fields=3-` && mv -v analysis/signalContamination/MC_2018NewMinimum_${suffix} analysis/signalContamination/MC_2018_${suffix}
done
