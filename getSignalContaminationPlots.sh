#!/bin/bash

./getMCSystematics/bin/getEventHistograms inputMCPathMain="${EOSPREFIX}/store/user/lpcsusystealth/selections/combined/merged_selection_MC_stealth_t5_2017_control_*.root" integratedLuminosityMain=41900.0 outputPrefix=control && ./getMCSystematics/bin/getMCUncertainties inputNEventsFile=analysis/dataSystematics/control_observedEventCounters.dat inputPath=analysis/MCEventHistograms/control_savedObjects.root outputPrefix=control unrestrictedSignalContamination=true minGluinoMass=975.0 nGluinoMassBins=16
