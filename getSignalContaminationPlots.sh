#!/bin/bash

./getMCSystematics/bin/getEventHistograms inputMCPathMain=${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/MC_2018Production_SingleMedium_optimized2017.root integratedLuminosityMain=41900.0 inputMCPathAux=${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/MC_2018Production_SingleMedium_optimized2016.root integratedLuminosityAux=35920.0 outputPrefix=control && ./getMCSystematics/bin/getMCUncertainties inputNEventsFile=analysis/dataSystematics/control_observedEventCounters.dat inputPath=analysis/MCEventHistograms/control_savedObjects.root outputPrefix=control unrestrictedSignalContamination=true minGluinoMass=975.0 nGluinoMassBins=16

./getMCSystematics/bin/getMCUncertainties inputPath=analysis/MCEventHistograms/MC_2018_savedObjects.root outputPrefix=MC_2018NewMinimum minGluinoMass=975.0 nGluinoMassBins=16 # Reproduce signal contamination plots with new minGluinoMass
rm -r analysis/MCSystematics/MC_2018NewMinimum*
for file in analysis/signalContamination/MC_2018NewMinimum_signalContamination_*; do
    suffix=`echo ${file} | cut --delimiter="/" --fields=3 | cut --delimiter="_" --fields=3-` && mv -v analysis/signalContamination/MC_2018NewMinimum_${suffix} analysis/signalContamination/MC_2018_${suffix}
done
