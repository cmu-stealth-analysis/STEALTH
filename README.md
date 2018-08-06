# Stealth SUSY Analysis Instructions
## Carnegie Mellon Group
### (M. Andrews, N. Bower, T. Mudholkar, M. Paulini, M. Sun, M. Weinberg)

### This documentation will be improved. For full details about the compulsory and optional arguments to pass to any executable (including the C++ executables!), please do: `executableName -h`; for example, `./mergeFiles.py -h`

## Creating N-tuples
### Event selection
More details to follow. Commands to run complete selection:

```
./submitJobs_selectEvents_Condor.py --inputFromFile --inputFilePath "inputFileList_data_DoubleEG_2016.txt" --outputDirectory "selections/DoublePhoton/DoubleFake" --outputFilePrefix "data_DoubleEG_2016_DoubleFake" --photonSelectionType fake --year 2016
./submitJobs_selectEvents_Condor.py --inputFromFile --inputFilePath "inputFileList_data_DoubleEG_2016.txt" --outputDirectory "selections/DoublePhoton/DoubleMedium" --outputFilePrefix "data_DoubleEG_2016_DoubleMedium" --photonSelectionType medium --year 2016
./submitJobs_selectEvents_Condor.py --inputFromFile --inputFilePath "inputFileList_data_DoubleEG_2016.txt" --outputDirectory "selections/DoublePhoton/OneMediumOneFake" --outputFilePrefix "data_DoubleEG_2016_OneMediumOneFake" --photonSelectionType mediumfake --year 2016
./submitJobs_selectEvents_Condor.py --inputFromFile --inputFilePath "inputFileList_data_DoubleEG_2017.txt" --outputDirectory "selections/DoublePhoton/DoubleFake" --outputFilePrefix "data_DoubleEG_2017_DoubleFake" --photonSelectionType fake --year 2017
./submitJobs_selectEvents_Condor.py --inputFromFile --inputFilePath "inputFileList_data_DoubleEG_2017.txt" --outputDirectory "selections/DoublePhoton/DoubleMedium" --outputFilePrefix "data_DoubleEG_2017_DoubleMedium" --photonSelectionType medium --year 2017
./submitJobs_selectEvents_Condor.py --inputFromFile --inputFilePath "inputFileList_data_DoubleEG_2017.txt" --outputDirectory "selections/DoublePhoton/OneMediumOneFake" --outputFilePrefix "data_DoubleEG_2017_OneMediumOneFake" --photonSelectionType mediumfake --year 2017
./submitJobs_selectEvents_Condor.py --inputFromFile --inputFilePath "inputFileList_MC_2018Production.txt" --outputDirectory "selections/DoublePhoton/DoubleMedium" --outputFilePrefix "MC_2018Production" --photonSelectionType mediumMC --JECUncertainty=0
./submitJobs_selectEvents_Condor.py --inputFromFile --inputFilePath "inputFileList_MC_2018Production.txt" --outputDirectory "selections/DoublePhoton/DoubleMedium" --outputFilePrefix "MC_2018Production_JECDown" --photonSelectionType mediumMC --JECUncertainty=-1
./submitJobs_selectEvents_Condor.py --inputFromFile --inputFilePath "inputFileList_MC_2018Production.txt" --outputDirectory "selections/DoublePhoton/DoubleMedium" --outputFilePrefix "MC_2018Production_JECUp" --photonSelectionType mediumMC --JECUncertainty=1
```
### Merging output files
```
./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/DoubleFake/data_DoubleEG_2016_DoubleFake_begin_*.root" --outputFilePath "${HOME}/nobackup/merged/data_DoubleEG_2016_DoubleFake.root"
./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/DoubleFake/data_DoubleEG_2017_DoubleFake_begin_*.root" --outputFilePath "${HOME}/nobackup/merged/data_DoubleEG_2017_DoubleFake.root"
./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/OneMediumOneFake/data_DoubleEG_2016_OneMediumOneFake_begin_*.root" --outputFilePath "${HOME}/nobackup/merged/data_DoubleEG_2016_OneMediumOneFake.root"
./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/OneMediumOneFake/data_DoubleEG_2017_OneMediumOneFake_begin_*.root" --outputFilePath "${HOME}/nobackup/merged/data_DoubleEG_2017_OneMediumOneFake.root"
./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/DoubleMedium/data_DoubleEG_2016_DoubleMedium_begin_*.root" --outputFilePath "${HOME}/nobackup/merged/data_DoubleEG_2016_DoubleMedium.root"
./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/DoubleMedium/data_DoubleEG_2017_DoubleMedium_begin_*.root" --outputFilePath "${HOME}/nobackup/merged/data_DoubleEG_2017_DoubleMedium.root"
./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/DoubleMedium/MC_2018Production_DoubleMedium_begin_*.root" --outputFilePath "${HOME}/nobackup/merged/MC_2018Production_DoubleMedium.root"
./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/DoubleMedium/MC_2018Production_DoubleMedium_JECUp_begin_*.root" --outputFilePath "${HOME}/nobackup/merged/MC_2018Production_DoubleMedium_JECUp.root"
./mergeFiles.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/DoublePhoton/DoubleMedium/MC_2018Production_DoubleMedium_JECDown_begin_*.root" --outputFilePath "${HOME}/nobackup/merged/MC_2018Production_DoubleMedium_JECDown.root"
```

## Full analysis chain
```
for step in `seq 1 7`; do ./runAnalysisSteps.sh ${step} 2016Plus2017; done
```
