#!/usr/bin/env python

from __future__ import print_function, division
import subprocess, argparse
import tmJDLInterface
import stealthEnv

inputArgumentsParser = argparse.ArgumentParser(description='Submit jobs that save sum of MC weights and PU weights.')
inputArgumentsParser.add_argument('--isProductionRun', action='store_true', help="By default, this script does not submit the actual jobs and instead only prints the shell command that would have been called. Passing this switch will execute the commands.")
inputArguments = inputArgumentsParser.parse_args()

n_subsamples = {
    "MC_HighHTQCD16": 7,
    "MC_HighHTQCD17": 8,
    "MC_HighHTQCD18": 8,
    "MC_GJetHT16": 5,
    "MC_GJetHT17": 4,
    "MC_GJetHT18": 4
}

input_file_list_and_output_details = []
for year_last_two_digits in [16, 17, 18]:
    year = 2000 + year_last_two_digits
    input_file_list_and_output_details.append(tuple(["DiPhotonJets_{y}".format(y=year), ["fileLists/inputFileList_MC_DiPhotonJets_{y}.txt".format(y=year)], "{eP}/{sER}/MCWeights/dataPU_{y4}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, y4=year), "{eP}/{sER}/MCWeights".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot), "sumMCWeights_DiPhotonJets_{y}.json".format(y=year), "PUWeights_DiPhotonJets_{y}.root".format(y=year)]))
    input_file_list_and_output_details.append(tuple(["DiPhotonJetsBox_{y}".format(y=year), ["fileLists/inputFileList_MC_DiPhotonJetsBox_{y}.txt".format(y=year)], "{eP}/{sER}/MCWeights/dataPU_{y4}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, y4=year), "{eP}/{sER}/MCWeights".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot), "sumMCWeights_DiPhotonJetsBox_{y}.json".format(y=year), "PUWeights_DiPhotonJetsBox_{y}.root".format(y=year)]))
    for MCBKGDatasetID in ["HighHTQCD", "GJetHT"]:
        for index_subsample in range(1, 1+n_subsamples["MC_{did}{y2}".format(did=MCBKGDatasetID, y2=year_last_two_digits)]):
            input_file_list_and_output_details.append(tuple(["{did}{i}_{y}".format(did=MCBKGDatasetID, i=index_subsample, y=year), ["fileLists/inputFileList_MC_{did}{i}_{y}.txt".format(did=MCBKGDatasetID, i=index_subsample, y=year)], "{eP}/{sER}/MCWeights/dataPU_{y4}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, y4=year), "{eP}/{sER}/MCWeights".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot), "sumMCWeights_{did}{i}_{y}.json".format(did=MCBKGDatasetID, i=index_subsample, y=year), "PUWeights_{did}{i}_{y}.root".format(did=MCBKGDatasetID, i=index_subsample, y=year)]))

for year in [2016, 2017, 2018]:
    for tDesignationIndex in [5, 6]:
        input_file_list_and_output_details.append(tuple(["stealth_t{td}_{y}".format(td=tDesignationIndex, y=year), ["fileLists/inputFileList_MC_Fall17_stealth_t{td}Wg.txt".format(td=tDesignationIndex)], "{eP}/{sER}/MCWeights/dataPU_{y4}.root".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot, y4=year), "{eP}/{sER}/MCWeights".format(eP=stealthEnv.EOSPrefix, sER=stealthEnv.stealthEOSRoot), "sumMCWeights_stealth_t{td}_{y}.json".format(td=tDesignationIndex, y=year), "PUWeights_stealth_t{td}_{y}.root".format(td=tDesignationIndex, y=year)]))

condor_output_folder = "{c}/stealth_MCWeights".format(c=stealthEnv.condorWorkAreaRoot)
subprocess.check_call("mkdir -p {o}".format(o=condor_output_folder), shell=True, executable="/bin/bash")
subprocess.check_call("rsync --checksum -a {r}/miscUtils/getSumMCWeightsHelper.sh {o}/getSumMCWeightsHelper.sh".format(r=stealthEnv.stealthRoot, o=condor_output_folder), shell=True, executable="/bin/bash")

# Make sure the tarballs to transfer are up to date
subprocess.check_call("cd {r}/miscUtils && ./update_MCReweightingScriptsTarball.sh && cd {r}".format(r=stealthEnv.stealthRoot), shell=True, executable="/bin/bash")
subprocess.check_call("cd {u} && ./update_tmUtilsTarball.sh && cd {r}".format(u=stealthEnv.tmUtilsParent, r=stealthEnv.stealthRoot), shell=True, executable="/bin/bash")

for identifier, input_paths_files_list, dataPUSourceWithXRDPrefix, outputEOSFolderWithXRDPrefix, outputFileName_sumWeights, outputFileName_PUWeights in input_file_list_and_output_details:
    print("Running for identifier: {ident}, input_paths_files_list: {i}, outputEOSFolderWithXRDPrefix: {oeos}, outputFileName_sumWeights: {op}, outputFileName_PUWeights: {opu}".format(ident=identifier, i=str(input_paths_files_list), oeos=outputEOSFolderWithXRDPrefix, op=outputFileName_sumWeights, opu=outputFileName_PUWeights))

    filesToTransfer = ["{xP}".format(xP=stealthEnv.x509Proxy), "{u}/tmUtils.tar.gz".format(u=stealthEnv.tmUtilsParent), "{u}/extract_tmUtilsTarball.sh".format(u=stealthEnv.tmUtilsParent), "{r}/setup_environment_remote.sh".format(r=stealthEnv.stealthRoot), "{r}/miscUtils/MCReweightingScripts.tar.gz".format(r=stealthEnv.stealthRoot), "{r}/miscUtils/extract_MCReweightingScriptsTarball.sh".format(r=stealthEnv.stealthRoot)]
    input_paths_files_string = ""
    for dataset_file in input_paths_files_list:
        filesToTransfer.append("{r}/{i}".format(r=stealthEnv.stealthRoot, i=dataset_file))
        input_paths_files_string += (dataset_file.replace("fileLists/", "") + ",")
    input_paths_files_string = input_paths_files_string[:-1] # To remove the last comma

    jdlInterface = tmJDLInterface.tmJDLInterface(processName=identifier, scriptPath="getSumMCWeightsHelper.sh", outputDirectoryRelativePath=condor_output_folder) # works even if "outputDirectoryRelativePath" is an absolute path
    jdlInterface.addFilesToTransferFromList(filesToTransfer)
    # Arguments for script:
    jdlInterface.addScriptArgument(input_paths_files_string) # Argument 1: comma-separated list of paths to files with input paths
    jdlInterface.addScriptArgument(dataPUSourceWithXRDPrefix) # Argument 2: PU source with XRD prefix
    jdlInterface.addScriptArgument(outputEOSFolderWithXRDPrefix) # Argument 3: output EOS folder with XRD prefix
    jdlInterface.addScriptArgument(outputFileName_sumWeights) # Argument 4: name for output json that stores sum of weights
    jdlInterface.addScriptArgument(outputFileName_PUWeights) # Argument 5: name for output root file that stores PU weights
    jdlInterface.writeToFile()
    submissionCommand = "cd {o} && condor_submit {ident}.jdl && cd {r}".format(o=condor_output_folder, ident=identifier, r=stealthEnv.stealthRoot)
    print ("Generated command: {sC}".format(sC=submissionCommand))
    if not(inputArguments.isProductionRun):
        print("Not submitting because isProductionRun flag is not set explicitly.")
    else:
        subprocess.check_call(submissionCommand, shell=True, executable="/bin/bash")
        print ("Submitted.")
